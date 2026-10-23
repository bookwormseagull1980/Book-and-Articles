"""Run the three-dimensional finite-window field generator."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys
import time

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from cosmic_cg.field3d import (
    FIELD_CLOSURE_COLUMNS,
    FIELD_COST_COLUMNS,
    FIELD_EPOCH_RELATION_COLUMNS,
    FIELD_MODE_COLUMNS,
    FIELD_OBJECT_COLUMNS,
    FIELD_RELATION_SUMMARY_COLUMNS,
    FIELD_TRANSFER_COLUMNS,
    FIELD_WINDOW_COLUMNS,
    Field3DConfig,
    available_physical_memory_bytes,
    config_to_dict as field_config_to_dict,
    estimate_field_resources,
    field_cost_row,
    generate_field_history,
)
from cosmic_cg.generation import InternalSpectrumConfig, config_to_dict as spectrum_config_to_dict
from cosmic_cg.interaction_adapters import field_interaction_outputs
from cosmic_cg.interaction_adapters import INTERACTION_TABLES
from cosmic_cg.intermediate import write_output_inventory
from cosmic_cg.io.field_audit import audit_field_run, write_field_audit_markdown, write_field_audit_report
from cosmic_cg.io.figure_data import write_field_figure_data
from cosmic_cg.io.provenance import create_run_paths, write_provenance
from cosmic_cg.io.tables import write_rows
from cosmic_cg.params import load_model_params, write_param_snapshot
from cosmic_cg.result_tables import directory_size_bytes
from cosmic_cg.run_plan import generation_redshifts


def main() -> int:
    args = parse_args()
    start = time.perf_counter()
    spectrum_config = InternalSpectrumConfig(
        c_nl=args.c_nl,
        k_nl0_hmpc=args.k_nl0_hmpc,
        n_s=args.n_s,
        k_pivot_hmpc=args.k_pivot_hmpc,
        k_min_hmpc=args.k_min_hmpc,
        k_max_hmpc=args.k_max_hmpc,
        k_points=args.k_points,
        endpoint_radius_factor=args.endpoint_radius_factor,
        spectrum_shape=args.spectrum_shape,
        uv_envelope_k_hmpc=args.uv_envelope_k_hmpc,
        uv_envelope_alpha=args.uv_envelope_alpha,
        uv_envelope_beta=args.uv_envelope_beta,
    )
    field_config = Field3DConfig(
        grid_size=args.grid_size,
        field_box_size_mpch=args.field_box_size_mpch,
        global_box_size_mpch=args.global_box_size_mpch,
        k_cut_hmpc=args.field_k_cut_hmpc,
        gaussian_sigma_mpch=args.gaussian_sigma_mpch,
        closure_floor=args.closure_floor,
        seed=args.seed,
        mode_bins=args.mode_bins,
        tail_width=args.tail_width,
        peak_threshold_sigma=args.peak_threshold_sigma,
        peak_exclusion_cells=args.peak_exclusion_cells,
        peak_patch_radius_cells=args.peak_patch_radius_cells,
        peak_slab_cells=args.peak_slab_cells,
        max_peak_candidates=args.max_peak_candidates,
        max_peaks_per_snapshot=args.max_peaks_per_snapshot,
        neighbor_links=args.neighbor_links,
        object_step_gyr=args.object_step_gyr,
    )
    redshifts = _redshifts(args)
    if args.dry_run:
        _print_dry_run(args=args, field_config=field_config, redshifts=redshifts)
        return 0
    _check_memory_budget(args=args, field_config=field_config, redshift_count=len(redshifts))

    run_paths = create_run_paths(args.results_root, args.run_id)
    params = load_model_params(args.params_file or None)
    write_param_snapshot(params, run_paths.root / "params_snapshot.json")

    config = _resolved_config(args, spectrum_config, field_config, redshifts)
    (run_paths.root / "config_resolved.yaml").write_text(_render_config(config), encoding="utf-8")

    outputs = generate_field_history(
        run_id=args.run_id,
        redshifts=redshifts,
        params=params,
        spectrum_config=spectrum_config,
        field_config=field_config,
        progress=_progress if args.progress else None,
    )

    write_rows(run_paths.tables / "field_window_table.csv", outputs["window"], FIELD_WINDOW_COLUMNS)
    write_rows(run_paths.tables / "field_mode_table.csv", outputs["mode"], FIELD_MODE_COLUMNS)
    write_rows(run_paths.tables / "field_object_table.csv", outputs["object"], FIELD_OBJECT_COLUMNS)
    write_rows(run_paths.tables / "field_transfer_table.csv", outputs["transfer"], FIELD_TRANSFER_COLUMNS)
    write_rows(
        run_paths.tables / "field_relation_summary_table.csv",
        outputs["relation_summary"],
        FIELD_RELATION_SUMMARY_COLUMNS,
    )
    write_rows(
        run_paths.tables / "field_epoch_relation_table.csv",
        outputs["epoch_relation"],
        FIELD_EPOCH_RELATION_COLUMNS,
    )
    write_rows(run_paths.tables / "field_closure_table.csv", outputs["closure"], FIELD_CLOSURE_COLUMNS)
    interaction_outputs = field_interaction_outputs(
        window_rows=outputs["window"],
        mode_rows=outputs["mode"],
        object_rows=outputs["object"],
        transfer_rows=outputs["transfer"],
        epoch_rows=outputs["epoch_relation"],
    )
    interaction_counts = {}
    for table_key, (filename, columns, output_key) in INTERACTION_TABLES.items():
        rows = interaction_outputs[output_key]
        write_rows(run_paths.tables / filename, rows, columns)
        interaction_counts[table_key] = len(rows)

    relationship_rows = sum(len(rows) for rows in outputs.values()) + sum(interaction_counts.values())
    relationship_bytes = directory_size_bytes(run_paths.tables)
    cost = field_cost_row(
        run_id=args.run_id,
        field_config=field_config,
        redshift_count=len(redshifts),
        relationship_rows=relationship_rows,
        relationship_bytes_written=relationship_bytes,
    )
    write_rows(run_paths.tables / "field_cost_table.csv", [cost], FIELD_COST_COLUMNS)
    figure_paths = write_field_figure_data(run_paths.root, run_paths.figure_data)
    write_output_inventory(run_paths.root)

    manifest = {
        "model": "finite_window_3d_field_generation",
        "run_id": args.run_id,
        "redshifts": redshifts,
        "field_config": field_config_to_dict(field_config),
        "spectrum": spectrum_config_to_dict(spectrum_config),
        "storage_policy": "relationship_tables_only",
        "external_data": [],
        "state_decomposition": {
            "physical_coarse_graining": "G_sigma",
            "activity_selector": "A_Lambda",
            "field_state": ["X_sigma", "X_A", "X_T"],
            "external_data": [],
        },
        "closure": {
            "map": "spectral_tail_moment_proxy",
            "residual_norm": "unresolved_gaussian_coarse_power",
            "epsilon0": args.closure_floor,
            "norm": "log_power_integral",
        },
        "output_tables": {
            "field_window_table": len(outputs["window"]),
            "field_mode_table": len(outputs["mode"]),
            "field_object_table": len(outputs["object"]),
            "field_transfer_table": len(outputs["transfer"]),
            "field_relation_summary_table": len(outputs["relation_summary"]),
            "field_epoch_relation_table": len(outputs["epoch_relation"]),
            "field_closure_table": len(outputs["closure"]),
            **interaction_counts,
        },
        "figure_data": [path.name for path in figure_paths],
    }
    (run_paths.root / "generation_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    write_provenance(run_paths.root / "provenance.json", run_id=args.run_id, config=config)
    write_output_inventory(run_paths.root)

    report = audit_field_run(run_paths.root, expected_grid_size=args.grid_size)
    write_field_audit_report(report, run_paths.quality_checks / "field_audit_report.json")
    write_field_audit_markdown(report, run_paths.quality_checks / "field_audit_report.md")
    if report.status == "fail":
        raise SystemExit(f"3D field audit failed: {run_paths.quality_checks / 'field_audit_report.md'}")

    print(
        json.dumps(
            {
                "run_id": args.run_id,
                "result_dir": str(run_paths.root),
                "grid_size": args.grid_size,
                "redshift_count": len(redshifts),
                "relationship_rows": relationship_rows,
                "relationship_bytes": relationship_bytes,
                "elapsed_seconds": time.perf_counter() - start,
                "audit_status": report.status,
            },
            indent=2,
        )
    )
    return 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-id", default="generated_universe_field_001")
    parser.add_argument("--results-root", default=str(PROJECT_ROOT / "results"))
    parser.add_argument("--params-file", default="")
    parser.add_argument("--redshifts", default="")
    parser.add_argument("--grid-size", type=int, default=512)
    parser.add_argument("--field-box-size-mpch", type=float, default=512.0)
    parser.add_argument("--global-box-size-mpch", type=float, default=14000.0)
    parser.add_argument("--field-k-cut-hmpc", type=float, default=1.0)
    parser.add_argument("--gaussian-sigma-mpch", type=float, default=1.0)
    parser.add_argument("--closure-floor", type=float, default=1.0e-12)
    parser.add_argument("--seed", type=int, default=20260903)
    parser.add_argument("--mode-bins", type=int, default=48)
    parser.add_argument("--tail-width", type=float, default=0.75)
    parser.add_argument("--peak-threshold-sigma", type=float, default=2.5)
    parser.add_argument("--peak-exclusion-cells", type=int, default=5)
    parser.add_argument("--peak-patch-radius-cells", type=int, default=3)
    parser.add_argument("--peak-slab-cells", type=int, default=16)
    parser.add_argument("--max-peak-candidates", type=int, default=8192)
    parser.add_argument("--max-peaks-per-snapshot", type=int, default=384)
    parser.add_argument("--neighbor-links", type=int, default=3)
    parser.add_argument("--object-step-gyr", type=float, default=0.25)
    parser.add_argument("--c-nl", type=float, default=1.0)
    parser.add_argument("--k-nl0-hmpc", type=float, default=0.2)
    parser.add_argument("--n-s", type=float, default=0.965)
    parser.add_argument("--k-pivot-hmpc", type=float, default=0.2)
    parser.add_argument("--k-min-hmpc", type=float, default=1.0e-4)
    parser.add_argument("--k-max-hmpc", type=float, default=200.0)
    parser.add_argument("--k-points", type=int, default=768)
    parser.add_argument("--endpoint-radius-factor", type=float, default=1.0)
    parser.add_argument("--spectrum-shape", choices=("finite_window_envelope", "power_law"), default="finite_window_envelope")
    parser.add_argument("--uv-envelope-k-hmpc", type=float, default=20.0)
    parser.add_argument("--uv-envelope-alpha", type=float, default=2.0)
    parser.add_argument("--uv-envelope-beta", type=float, default=2.0)
    parser.add_argument("--dry-run", action="store_true")
    parser.add_argument("--skip-memory-check", action="store_true")
    parser.add_argument("--memory-safety-factor", type=float, default=1.25)
    parser.add_argument("--progress", action="store_true")
    return parser.parse_args()


def _redshifts(args: argparse.Namespace) -> list[float]:
    if not args.redshifts:
        return list(generation_redshifts())
    return [float(item.strip()) for item in args.redshifts.split(",") if item.strip()]


def _resolved_config(
    args: argparse.Namespace,
    spectrum_config: InternalSpectrumConfig,
    field_config: Field3DConfig,
    redshifts: list[float],
) -> dict[str, object]:
    config = vars(args).copy()
    config["redshift_values"] = redshifts
    config["spectrum_config"] = spectrum_config_to_dict(spectrum_config)
    config["field_config"] = field_config_to_dict(field_config)
    config["external_data"] = []
    config["storage_policy"] = "relationship_tables_only"
    return config


def _render_config(config: dict[str, object]) -> str:
    lines = ["# Resolved 3D field-generation configuration"]
    for key in sorted(config):
        lines.append(f"{key}: {json.dumps(config[key], ensure_ascii=False)}")
    return "\n".join(lines) + "\n"


def _print_dry_run(*, args: argparse.Namespace, field_config: Field3DConfig, redshifts: list[float]) -> None:
    estimate = estimate_field_resources(redshift_count=len(redshifts), field_config=field_config)
    available = available_physical_memory_bytes()
    required = int(float(estimate["estimated_peak_memory_gb"]) * 1024**3 * args.memory_safety_factor)
    payload = {
        "run_id": args.run_id,
        "grid_size": args.grid_size,
        "redshift_count": len(redshifts),
        "redshifts": redshifts,
        "estimated_peak_memory_gb": estimate["estimated_peak_memory_gb"],
        "memory_safety_factor": args.memory_safety_factor,
        "required_memory_gb": required / 1024**3,
        "available_memory_gb": None if available is None else available / 1024**3,
        "raw_density_bytes_all_snapshots": estimate["raw_density_bytes_all_snapshots"],
        "gaussian_sigma_mpch": field_config.gaussian_sigma_mpch,
        "closure_floor": field_config.closure_floor,
        "storage_policy": "relationship_tables_only",
        "external_data": [],
        "can_start": True if available is None else available >= required,
    }
    print(json.dumps(payload, indent=2))


def _check_memory_budget(*, args: argparse.Namespace, field_config: Field3DConfig, redshift_count: int) -> None:
    if args.skip_memory_check:
        return
    estimate = estimate_field_resources(redshift_count=redshift_count, field_config=field_config)
    available = available_physical_memory_bytes()
    if available is None:
        return
    required = int(float(estimate["estimated_peak_memory_gb"]) * 1024**3 * args.memory_safety_factor)
    if available < required:
        message = {
            "grid_size": args.grid_size,
            "estimated_peak_memory_gb": estimate["estimated_peak_memory_gb"],
            "memory_safety_factor": args.memory_safety_factor,
            "required_memory_gb": required / 1024**3,
            "available_memory_gb": available / 1024**3,
            "action": "reduce grid size, close other applications, or rerun with --skip-memory-check",
        }
        raise SystemExit("memory budget check failed: " + json.dumps(message, indent=2))


def _progress(message: str) -> None:
    print(message, flush=True)


if __name__ == "__main__":
    raise SystemExit(main())
