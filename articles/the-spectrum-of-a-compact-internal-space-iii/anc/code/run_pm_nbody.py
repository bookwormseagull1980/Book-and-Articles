"""Run the internal particle-mesh trajectory solver."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import time

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from cosmic_cg.intermediate import write_output_inventory
from cosmic_cg.interaction_adapters import INTERACTION_TABLES, pm_interaction_outputs
from cosmic_cg.io.pm_audit import audit_pm_run, write_pm_audit_markdown, write_pm_audit_report
from cosmic_cg.io.provenance import create_run_paths, write_provenance
from cosmic_cg.io.tables import write_rows
from cosmic_cg.params import load_model_params, write_param_snapshot
from cosmic_cg.pm_nbody import (
    PM_COST_COLUMNS,
    PM_HALO_COLUMNS,
    PM_POWER_COLUMNS,
    PM_SNAPSHOT_COLUMNS,
    PM_TRANSFER_COLUMNS,
    PM_TREE_COLUMNS,
    PMNBodyConfig,
    available_physical_memory_bytes,
    config_to_dict as pm_config_to_dict,
    estimate_pm_resources,
    generate_pm_history,
    pm_cost_row,
)
from cosmic_cg.result_tables import directory_size_bytes
from cosmic_cg.run_plan import generation_redshifts


def main() -> int:
    args = parse_args()
    start = time.perf_counter()
    redshifts = _redshifts(args)
    config = PMNBodyConfig(
        dynamics=args.dynamics,
        particles_per_dim=args.particles_per_dim,
        mesh_size=args.mesh_size,
        box_size_mpch=args.box_size_mpch,
        seed=args.seed,
        z_initial=args.z_initial,
        steps_per_interval=args.steps_per_interval,
        initial_displacement_fraction=args.initial_displacement_fraction,
        initial_velocity_factor=args.initial_velocity_factor,
        force_strength=args.force_strength,
        mode_bins=args.mode_bins,
        density_threshold_sigma=args.density_threshold_sigma,
        density_threshold_delta=args.density_threshold_delta,
        fof_min_particles=args.fof_min_particles,
        max_halos_per_snapshot=args.max_halos_per_snapshot,
        neighbor_links=args.neighbor_links,
    )
    if args.dry_run:
        _print_dry_run(args=args, config=config, redshifts=redshifts)
        return 0
    _check_memory_budget(args=args, config=config, redshift_count=len(redshifts))

    run_paths = create_run_paths(args.results_root, args.run_id)
    params = load_model_params(args.params_file or None)
    write_param_snapshot(params, run_paths.root / "params_snapshot.json")

    resolved = _resolved_config(args, config, redshifts)
    (run_paths.root / "config_resolved.yaml").write_text(_render_config(resolved), encoding="utf-8")

    outputs = generate_pm_history(
        run_id=args.run_id,
        redshifts=redshifts,
        params=params,
        config=config,
        progress=_progress if args.progress else None,
    )
    write_rows(run_paths.tables / "pm_snapshot_table.csv", outputs["snapshot"], PM_SNAPSHOT_COLUMNS)
    write_rows(run_paths.tables / "pm_power_table.csv", outputs["power"], PM_POWER_COLUMNS)
    write_rows(run_paths.tables / "pm_halo_table.csv", outputs["halo"], PM_HALO_COLUMNS)
    write_rows(run_paths.tables / "pm_tree_table.csv", outputs["tree"], PM_TREE_COLUMNS)
    write_rows(run_paths.tables / "pm_transfer_table.csv", outputs["transfer"], PM_TRANSFER_COLUMNS)
    interaction_outputs = pm_interaction_outputs(
        snapshot_rows=outputs["snapshot"],
        power_rows=outputs["power"],
        halo_rows=outputs["halo"],
        transfer_rows=outputs["transfer"],
    )
    interaction_counts = {}
    for table_key, (filename, columns, output_key) in INTERACTION_TABLES.items():
        rows = interaction_outputs[output_key]
        write_rows(run_paths.tables / filename, rows, columns)
        interaction_counts[table_key] = len(rows)

    relationship_rows = sum(len(rows) for rows in outputs.values()) + sum(interaction_counts.values())
    relationship_bytes = directory_size_bytes(run_paths.tables)
    cost = pm_cost_row(
        run_id=args.run_id,
        config=config,
        redshift_count=len(redshifts),
        relationship_rows=relationship_rows,
        relationship_bytes_written=relationship_bytes,
    )
    write_rows(run_paths.tables / "pm_cost_table.csv", [cost], PM_COST_COLUMNS)

    manifest = {
        "model": "finite_window_internal_pm_nbody",
        "run_id": args.run_id,
        "redshifts": redshifts,
        "pm_config": pm_config_to_dict(config),
        "storage_policy": "relationship_tables_only",
        "external_data": [],
        "state_decomposition": {
            "dynamical_state": [
                "particle_phase_space",
                "mesh_density",
                "resolved_force",
                "halo_catalog",
                "transfer_graph",
            ],
            "physical_coarse_graining": "particle_mesh_assignment",
            "activity_selector": "resolved_mesh_modes",
            "external_data": [],
        },
        "closure": {
            "map": "mesh_resolved_poisson_force",
            "residual_norm": "unresolved_modes_above_mesh_nyquist",
            "audit": "mass_conservation_power_bins_particle_id_overlap",
        },
        "dynamics": {
            "solver": "particle_mesh",
            "force": "spectral_poisson",
            "integrator": "kick_drift_kick",
            "tracker": "particle_id_overlap",
            "time_variable": "scale_factor",
            "clock": args.dynamics,
            "force_prefactor": "3*Omega_m0/(2*a^2*E(a))" if args.dynamics == "cosmological" else "Q(a)=force_strength/max(a,1e-4)",
            "momentum": "p=a^3*E(a)*dx/da" if args.dynamics == "cosmological" else "u=dx/da",
            "drift_prefactor": "1/(a^3*E(a))" if args.dynamics == "cosmological" else "1",
            "initial_velocity": "f(a_initial)*displacement/a_initial" if args.dynamics == "cosmological" else "declared multiplier",
            "force_strength": None if args.dynamics == "cosmological" else args.force_strength,
        },
        "output_tables": {
            "pm_snapshot_table": len(outputs["snapshot"]),
            "pm_power_table": len(outputs["power"]),
            "pm_halo_table": len(outputs["halo"]),
            "pm_tree_table": len(outputs["tree"]),
            "pm_transfer_table": len(outputs["transfer"]),
            "pm_cost_table": 1,
            **interaction_counts,
        },
    }
    (run_paths.root / "generation_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    write_provenance(run_paths.root / "provenance.json", run_id=args.run_id, config=resolved)
    write_output_inventory(run_paths.root)

    report = audit_pm_run(
        run_paths.root,
        expected_particles_per_dim=args.particles_per_dim,
        expected_mesh_size=args.mesh_size,
    )
    write_pm_audit_report(report, run_paths.quality_checks / "pm_audit_report.json")
    write_pm_audit_markdown(report, run_paths.quality_checks / "pm_audit_report.md")
    write_output_inventory(run_paths.root)
    if report.status == "fail":
        raise SystemExit(f"PM audit failed: {run_paths.quality_checks / 'pm_audit_report.md'}")

    print(
        json.dumps(
            {
                "run_id": args.run_id,
                "result_dir": str(run_paths.root),
                "particles_per_dim": args.particles_per_dim,
                "mesh_size": args.mesh_size,
                "redshift_count": len(redshifts),
                "relationship_rows": relationship_rows,
                "relationship_bytes": directory_size_bytes(run_paths.tables),
                "elapsed_seconds": time.perf_counter() - start,
                "audit_status": report.status,
            },
            indent=2,
        )
    )
    return 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-id", default="pm_dynamics_001")
    parser.add_argument("--results-root", default=str(PROJECT_ROOT / "results"))
    parser.add_argument("--params-file", default="")
    parser.add_argument("--dynamics", choices=["cosmological", "comparison"], default="cosmological")
    parser.add_argument("--redshifts", default="")
    parser.add_argument("--particles-per-dim", type=int, default=32)
    parser.add_argument("--mesh-size", type=int, default=64)
    parser.add_argument("--box-size-mpch", type=float, default=128.0)
    parser.add_argument("--seed", type=int, default=20260904)
    parser.add_argument("--z-initial", type=float, default=120.0)
    parser.add_argument("--steps-per-interval", type=int, default=4)
    parser.add_argument("--initial-displacement-fraction", type=float, default=0.22)
    parser.add_argument("--initial-velocity-factor", type=float, default=0.35)
    parser.add_argument("--force-strength", type=float, default=0.35)
    parser.add_argument("--mode-bins", type=int, default=24)
    parser.add_argument("--density-threshold-sigma", type=float, default=1.5)
    parser.add_argument("--density-threshold-delta", type=float, default=0.8)
    parser.add_argument("--fof-min-particles", type=int, default=12)
    parser.add_argument("--max-halos-per-snapshot", type=int, default=96)
    parser.add_argument("--neighbor-links", type=int, default=3)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--skip-memory-check", action="store_true")
    parser.add_argument("--memory-safety-factor", type=float, default=1.25)
    parser.add_argument("--progress", action="store_true")
    return parser.parse_args()


def _redshifts(args: argparse.Namespace) -> list[float]:
    if not args.redshifts:
        return list(generation_redshifts())
    return [float(item.strip()) for item in args.redshifts.split(",") if item.strip()]


def _resolved_config(args: argparse.Namespace, config: PMNBodyConfig, redshifts: list[float]) -> dict[str, object]:
    resolved = vars(args).copy()
    resolved["redshift_values"] = redshifts
    resolved["pm_config"] = pm_config_to_dict(config)
    resolved["storage_policy"] = "relationship_tables_only"
    resolved["external_data"] = []
    resolved["state_decomposition"] = {
        "dynamical_state": [
            "particle_phase_space",
            "mesh_density",
            "resolved_force",
            "halo_catalog",
            "transfer_graph",
        ],
        "physical_coarse_graining": "particle_mesh_assignment",
        "activity_selector": "resolved_mesh_modes",
        "external_data": [],
    }
    return resolved


def _render_config(config: dict[str, object]) -> str:
    lines = ["# Resolved internal particle-mesh configuration"]
    for key in sorted(config):
        lines.append(f"{key}: {json.dumps(config[key], ensure_ascii=False)}")
    return "\n".join(lines) + "\n"


def _print_dry_run(*, args: argparse.Namespace, config: PMNBodyConfig, redshifts: list[float]) -> None:
    estimate = estimate_pm_resources(redshift_count=len(redshifts), config=config)
    available = available_physical_memory_bytes()
    required = int(float(estimate["estimated_peak_memory_gb"]) * 1024**3 * args.memory_safety_factor)
    payload = {
        "run_id": args.run_id,
        "particles_per_dim": args.particles_per_dim,
        "mesh_size": args.mesh_size,
        "redshift_count": len(redshifts),
        "redshifts": redshifts,
        "estimated_peak_memory_gb": estimate["estimated_peak_memory_gb"],
        "memory_safety_factor": args.memory_safety_factor,
        "required_memory_gb": required / 1024**3,
        "available_memory_gb": None if available is None else available / 1024**3,
        "raw_bytes_all_snapshots": estimate["raw_bytes_all_snapshots"],
        "storage_policy": "relationship_tables_only",
        "external_data": [],
        "state_decomposition": {
            "dynamical_state": [
                "particle_phase_space",
                "mesh_density",
                "resolved_force",
                "halo_catalog",
                "transfer_graph",
            ],
            "physical_coarse_graining": "particle_mesh_assignment",
            "activity_selector": "resolved_mesh_modes",
            "external_data": [],
        },
        "can_start": True if available is None else available >= required,
    }
    print(json.dumps(payload, indent=2))


def _check_memory_budget(*, args: argparse.Namespace, config: PMNBodyConfig, redshift_count: int) -> None:
    if args.skip_memory_check:
        return
    estimate = estimate_pm_resources(redshift_count=redshift_count, config=config)
    available = available_physical_memory_bytes()
    if available is None:
        return
    required = int(float(estimate["estimated_peak_memory_gb"]) * 1024**3 * args.memory_safety_factor)
    if available < required:
        message = {
            "particles_per_dim": args.particles_per_dim,
            "mesh_size": args.mesh_size,
            "estimated_peak_memory_gb": estimate["estimated_peak_memory_gb"],
            "memory_safety_factor": args.memory_safety_factor,
            "required_memory_gb": required / 1024**3,
            "available_memory_gb": available / 1024**3,
            "action": "reduce particle or mesh resolution, close other applications, or rerun with --skip-memory-check",
        }
        raise SystemExit("memory budget check failed: " + json.dumps(message, indent=2))


def _progress(message: str) -> None:
    print(message, flush=True)


if __name__ == "__main__":
    raise SystemExit(main())
