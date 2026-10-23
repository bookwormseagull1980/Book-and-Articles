"""Generate a finite-window universe history from the internal equations."""

from __future__ import annotations

import argparse
import json
import math
import time
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from cosmic_cg.entropy_audit import audit_entropy
from cosmic_cg.generation import (
    InternalSpectrumConfig,
    config_to_dict,
    generate_epoch,
    generate_objects,
)
from cosmic_cg.intermediate import (
    ENDPOINT_COLUMNS,
    HALO_SELECTION_COLUMNS,
    INPUT_VALIDATION_COLUMNS,
    MERGER_GROWTH_COLUMNS,
    POWER_SAMPLE_COLUMNS,
    write_output_inventory,
)
from cosmic_cg.interactions import (
    CARRIER_HISTORY_COLUMNS,
    EMERGENCE_EVENT_COLUMNS,
    INTERACTION_BUDGET_COLUMNS,
    INTERACTION_CROSSING_COLUMNS,
    INTERACTION_DOMINANCE_COLUMNS,
    INTERACTION_STRENGTH_COLUMNS,
    generate_interaction_outputs,
)
from cosmic_cg.io.audit import audit_run, write_audit_markdown, write_audit_report
from cosmic_cg.io.figure_data import write_standard_figure_data
from cosmic_cg.io.provenance import create_run_paths, write_provenance
from cosmic_cg.io.schema_catalog import SCHEMA_CATALOG_COLUMNS, schema_catalog_rows
from cosmic_cg.io.tables import write_rows
from cosmic_cg.knl import dimensionless_power, residual_floor_fraction, tail_power_fraction, tail_power_moment_fraction
from cosmic_cg.params import load_model_params, write_param_snapshot
from cosmic_cg.result_tables import (
    BAND_TABLE_COLUMNS,
    CLOSURE_TABLE_COLUMNS,
    COST_TABLE_COLUMNS,
    EVENT_TABLE_COLUMNS,
    NORMALIZATION_TABLE_COLUMNS,
    OBJECT_SUMMARY_TABLE_COLUMNS,
    OBJECT_TABLE_COLUMNS,
    TAIL_TABLE_COLUMNS,
    TRANSFER_SUMMARY_TABLE_COLUMNS,
    TRANSFER_TABLE_COLUMNS,
    WINDOW_TABLE_COLUMNS,
    band_table_rows,
    closure_table_row,
    cost_table_row,
    event_table_rows,
    normalization_table_rows,
    object_summary_rows,
    object_table_rows,
    tail_table_row,
    transfer_summary_rows,
    window_table_row,
)
from cosmic_cg.run_plan import generation_redshifts
from cosmic_cg.tail_ratio import halo_tail_ratio


def main() -> int:
    args = parse_args()
    start = time.perf_counter()
    run_paths = create_run_paths(args.results_root, args.run_id)
    config = _resolved_config(args)
    (run_paths.root / "config_resolved.yaml").write_text(_render_config(config), encoding="utf-8")

    params = load_model_params(args.params_file or None)
    write_param_snapshot(params, run_paths.root / "params_snapshot.json")
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

    global_rows = []
    window_rows = []
    band_rows_all = []
    tail_rows = []
    object_rows = []
    merger_rows = []
    transfer_rows = []
    event_rows = []
    normalization_rows = []
    closure_rows = []
    power_sample_rows = []
    endpoint_rows = []

    for z in _redshifts(args):
        epoch = generate_epoch(z=z, params=params, config=spectrum_config)
        k_box = 2.0 * 3.141592653589793 / args.box_size_mpch
        entropy = audit_entropy(
            s_abs=params.s_abs,
            k_nl_hmpc=epoch.knl.k_nl_hmpc,
            k_box_hmpc=k_box,
            k_cut_hmpc=args.k_cut_hmpc,
        )
        generated_objects, generated_links, generated_transfers = generate_objects(
            run_id=args.run_id,
            epoch=epoch,
            params=params,
            objects_per_band=args.objects_per_band,
            object_step_gyr=args.object_step_gyr,
        )
        tail = halo_tail_ratio(
            generated_objects,
            m_nl_msun_h=epoch.m_nl_msun_h,
            r_nl_kpc_h=epoch.r_nl_kpc_h,
            max_mass_factor=args.max_mass_factor,
        )
        omega_tail = tail_power_fraction(epoch.k_hmpc, epoch.p_mpc3_h3, epoch.knl.k_nl_hmpc)
        tail_moment = tail_power_moment_fraction(epoch.k_hmpc, epoch.p_mpc3_h3, epoch.knl.k_nl_hmpc)
        omega_floor = residual_floor_fraction(
            epoch.k_hmpc,
            epoch.p_mpc3_h3,
            epoch.knl.k_nl_hmpc,
            q_floor=args.tail_floor_q,
            sigma_g=args.tail_floor_sigma,
        )
        closure = _closure_audit(
            epoch=epoch,
            gaussian_sigma_mpch=args.gaussian_sigma_mpch,
            k_resolved_hmpc=args.k_cut_hmpc,
            activity_tail_width=args.activity_tail_width,
            epsilon0=args.closure_floor,
        )
        transfer_budget = _transfer_budget(generated_transfers, generated_objects, omega_tail, tail_moment, omega_floor)
        global_row = {
            "run_id": args.run_id,
            "source": "internal_generation",
            "snapshot": f"a{_tag(epoch.a)}",
            "z": epoch.z,
            "a": epoch.a,
            "box_size_mpch": args.box_size_mpch,
            "k_box_hmpc": k_box,
            "gaussian_sigma_mpch": args.gaussian_sigma_mpch,
            "gaussian_kernel": "exp[-0.5*(sigma*k)^2]",
            "activity_endpoint_hmpc": epoch.knl.k_nl_hmpc,
            "activity_selector": "log_gaussian_tail_with_cut",
            "C_NL": args.c_nl,
            "k_NL_hmpc": epoch.knl.k_nl_hmpc,
            "M_NL_msun_h": epoch.m_nl_msun_h,
            "R_NL_kpc_h": epoch.r_nl_kpc_h,
            "v_vir_MNL_kms": epoch.v_vir_mnl_kms,
            "S_abs": params.s_abs,
            "S_AW": entropy.s_aw,
            "k_cut_hmpc": args.k_cut_hmpc,
            "D_S": entropy.d_s,
            "Omega_tail": omega_tail,
            "Omega_floor": omega_floor,
            "Omega_excess": max(0.0, omega_tail - omega_floor),
            "R_h_mean": tail.r_h,
            "R_h_median": tail.median_per_halo,
            "R_h_p16": tail.p16_per_halo,
            "R_h_p84": tail.p84_per_halo,
            "R_v_optional": "",
            "v_lambda": "",
            "tail_moment_s1": tail_moment,
            "tail_injection_s0": transfer_budget["i_t_s0"],
            "tail_injection_s1": transfer_budget["i_t_s1"],
            "tail_drain_s0": transfer_budget["d_t_s0"],
            "tail_drain_s1": transfer_budget["d_t_s1"],
            "closure_map": "spectral_tail_moment_proxy",
            "closure_active_norm": closure["active_norm"],
            "closure_model_norm": closure["tail_model_norm"],
            "closure_residual_norm": closure["residual_norm"],
            "closure_error": closure["closure_error"],
            "closure_epsilon0": args.closure_floor,
            "closure_norm": "log_power_integral",
            "closure_status": closure["quality_flag"],
            "spectrum_source": f"internal_{args.spectrum_shape}",
            "quality_flag": "internal_generation",
            "notes": "generated from the finite-window internal equations",
        }
        bands = band_table_rows(global_row, epoch.k_hmpc, epoch.p_mpc3_h3)
        global_rows.append(global_row)
        band_rows_all.extend(bands)
        tail_rows.append(tail_table_row(global_row))
        closure_rows.append(closure_table_row(global_row))
        object_rows.extend(generated_objects)
        merger_rows.extend(generated_links)
        transfer_rows.extend(generated_transfers)
        event_rows.extend(event_table_rows(global_row, bands))
        normalization_rows.extend(normalization_table_rows(global_row, len(generated_objects)))
        power_sample_rows.extend(_power_rows(epoch))
        endpoint_rows.append(_endpoint_row(epoch))

    _attach_endpoint_velocity(global_rows)
    window_rows = [window_table_row(row) for row in global_rows]
    table_start = time.perf_counter()
    write_rows(run_paths.intermediate / "input_validation.csv", _input_rows(args, spectrum_config), INPUT_VALIDATION_COLUMNS)
    write_rows(run_paths.intermediate / "power_spectrum_samples.csv", power_sample_rows, POWER_SAMPLE_COLUMNS)
    write_rows(run_paths.intermediate / "endpoint_solution.csv", endpoint_rows, ENDPOINT_COLUMNS)
    write_rows(run_paths.intermediate / "halo_selection.csv", _halo_selection_rows(object_rows, global_rows), HALO_SELECTION_COLUMNS)
    write_rows(run_paths.intermediate / "merger_growth_links.csv", _merger_rows(merger_rows), MERGER_GROWTH_COLUMNS)

    write_rows(run_paths.tables / "global_window.csv", global_rows)
    write_rows(run_paths.tables / "halo_objects.csv", object_rows)
    write_rows(run_paths.tables / "window_table.csv", window_rows, WINDOW_TABLE_COLUMNS)
    write_rows(run_paths.tables / "band_table.csv", band_rows_all, BAND_TABLE_COLUMNS)
    write_rows(run_paths.tables / "tail_table.csv", tail_rows, TAIL_TABLE_COLUMNS)
    write_rows(run_paths.tables / "closure_table.csv", closure_rows, CLOSURE_TABLE_COLUMNS)
    write_rows(run_paths.tables / "object_table.csv", object_table_rows(object_rows), OBJECT_TABLE_COLUMNS)
    write_rows(run_paths.tables / "transfer_table.csv", transfer_rows, TRANSFER_TABLE_COLUMNS)
    write_rows(run_paths.tables / "event_table.csv", event_rows, EVENT_TABLE_COLUMNS)
    write_rows(run_paths.tables / "normalization_table.csv", normalization_rows, NORMALIZATION_TABLE_COLUMNS)
    write_rows(run_paths.tables / "object_summary_table.csv", object_summary_rows(object_rows), OBJECT_SUMMARY_TABLE_COLUMNS)
    write_rows(run_paths.tables / "transfer_summary_table.csv", transfer_summary_rows(transfer_rows), TRANSFER_SUMMARY_TABLE_COLUMNS)
    interaction_outputs = generate_interaction_outputs(
        run_id=args.run_id,
        global_rows=global_rows,
        band_rows=band_rows_all,
        object_rows=object_table_rows(object_rows),
        transfer_rows=transfer_rows,
    )
    write_rows(run_paths.tables / "interaction_strength_table.csv", interaction_outputs["strength"], INTERACTION_STRENGTH_COLUMNS)
    write_rows(run_paths.tables / "interaction_budget_table.csv", interaction_outputs["budget"], INTERACTION_BUDGET_COLUMNS)
    write_rows(run_paths.tables / "interaction_dominance_table.csv", interaction_outputs["dominance"], INTERACTION_DOMINANCE_COLUMNS)
    write_rows(run_paths.tables / "interaction_crossing_table.csv", interaction_outputs["crossing"], INTERACTION_CROSSING_COLUMNS)
    write_rows(run_paths.tables / "carrier_history_table.csv", interaction_outputs["carrier"], CARRIER_HISTORY_COLUMNS)
    write_rows(run_paths.tables / "emergence_event_table.csv", interaction_outputs["emergence"], EMERGENCE_EVENT_COLUMNS)
    write_rows(run_paths.tables / "schema_catalog.csv", schema_catalog_rows(), SCHEMA_CATALOG_COLUMNS)

    table_seconds = time.perf_counter() - table_start
    cost_rows = [
        cost_table_row(
            row,
            run_paths.root,
            k_ref_hmpc=args.k_ref_hmpc,
            time_exponent=args.time_exponent,
            access_seconds=0.0,
            audit_seconds=table_seconds,
        )
        for row in global_rows
    ]
    write_rows(run_paths.tables / "cost_table.csv", cost_rows, COST_TABLE_COLUMNS)

    manifest = {
        "model": "finite_window_internal_generation",
        "redshifts": [row["z"] for row in global_rows],
        "object_count": len(object_rows),
        "transfer_count": len(transfer_rows),
        "interaction_strength_count": len(interaction_outputs["strength"]),
        "interaction_channel_count": len({row["channel_id"] for row in interaction_outputs["strength"]}),
        "interaction_band_count": len({row["band_id"] for row in interaction_outputs["strength"]}),
        "spectrum": config_to_dict(spectrum_config),
        "inputs": ["parameter_snapshot", "run_configuration"],
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
    }
    (run_paths.root / "generation_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    write_provenance(run_paths.root / "provenance.json", run_id=args.run_id, config=config)
    write_output_inventory(run_paths.root)
    if args.write_figure_data:
        write_standard_figure_data(run_paths.root, run_paths.figure_data)
        write_output_inventory(run_paths.root)

    report = audit_run(run_paths.root, "generated_universe")
    write_audit_report(report, run_paths.quality_checks / "audit_report.json")
    write_audit_markdown(report, run_paths.quality_checks / "audit_report.md")
    if report.status == "fail":
        raise SystemExit(f"run audit failed: {run_paths.quality_checks / 'audit_report.md'}")

    print(
        json.dumps(
            {
                "run_id": args.run_id,
                "result_dir": str(run_paths.root),
                "redshift_count": len(global_rows),
                "object_rows": len(object_rows),
                "transfer_rows": len(transfer_rows),
                "figure_data_written": args.write_figure_data,
                "elapsed_seconds": time.perf_counter() - start,
                "audit_status": report.status,
            },
            indent=2,
        )
    )
    return 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-id", default="generated_universe_001")
    parser.add_argument("--results-root", default=str(PROJECT_ROOT / "results"))
    parser.add_argument("--params-file", default="")
    parser.add_argument("--redshifts", default="")
    parser.add_argument("--objects-per-band", type=int, default=24)
    parser.add_argument("--object-step-gyr", type=float, default=0.2)
    parser.add_argument("--box-size-mpch", type=float, default=14000.0)
    parser.add_argument("--k-cut-hmpc", type=float, default=20.0)
    parser.add_argument("--k-ref-hmpc", type=float, default=200.0)
    parser.add_argument("--time-exponent", type=float, default=2.0)
    parser.add_argument("--gaussian-sigma-mpch", type=float, default=1.0)
    parser.add_argument("--activity-tail-width", type=float, default=0.75)
    parser.add_argument("--closure-floor", type=float, default=1.0e-12)
    parser.add_argument("--max-mass-factor", type=float, default=1.0)
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
    parser.add_argument("--tail-floor-q", type=float, default=2.0)
    parser.add_argument("--tail-floor-sigma", type=float, default=0.5)
    parser.add_argument("--write-figure-data", action="store_true")
    return parser.parse_args()


def _redshifts(args: argparse.Namespace) -> list[float]:
    if not args.redshifts:
        return list(generation_redshifts())
    return [float(item.strip()) for item in args.redshifts.split(",") if item.strip()]


def _resolved_config(args: argparse.Namespace) -> dict[str, object]:
    return vars(args).copy()


def _render_config(config: dict[str, object]) -> str:
    lines = ["# Resolved generated-universe configuration"]
    for key in sorted(config):
        lines.append(f"{key}: {json.dumps(config[key], ensure_ascii=False)}")
    return "\n".join(lines) + "\n"


def _input_rows(args: argparse.Namespace, spectrum_config: InternalSpectrumConfig) -> list[dict[str, object]]:
    return [
        {
            "schema": "internal_spectrum",
            "path": "internal",
            "ok": True,
            "row_count": spectrum_config.k_points,
            "missing_columns": "",
        },
        {
            "schema": "internal_object_rule",
            "path": "internal",
            "ok": True,
            "row_count": args.objects_per_band,
            "missing_columns": "",
        },
    ]


def _power_rows(epoch) -> list[dict[str, object]]:
    delta2 = dimensionless_power(epoch.k_hmpc, epoch.p_mpc3_h3)
    rows = []
    spacing = _median_log_spacing(epoch.k_hmpc)
    for k_value, p_value, d2 in zip(epoch.k_hmpc, epoch.p_mpc3_h3, delta2):
        ratio = k_value / epoch.knl.k_nl_hmpc
        distance = max(0.0, math.log(ratio) / spacing)
        rows.append(
            {
                "z": epoch.z,
                "a": epoch.a,
                "k_hmpc": k_value,
                "p_mpc3_h3": p_value,
                "delta2": d2,
                "k_over_knl": ratio,
                "sector": "active" if k_value <= epoch.knl.k_nl_hmpc else "tail",
                "tail_weight_s1": math.exp(-distance),
            }
        )
    return rows


def _endpoint_row(epoch) -> dict[str, object]:
    delta2 = dimensionless_power(epoch.k_hmpc, epoch.p_mpc3_h3)
    i = max(0, min(epoch.knl.crossing_index, len(epoch.k_hmpc) - 2))
    return {
        "z": epoch.z,
        "c_nl": epoch.knl.c_nl,
        "k_nl_hmpc": epoch.knl.k_nl_hmpc,
        "crossing_index": epoch.knl.crossing_index,
        "method": epoch.knl.method,
        "k_left_hmpc": epoch.k_hmpc[i],
        "k_right_hmpc": epoch.k_hmpc[i + 1],
        "delta2_left": delta2[i],
        "delta2_right": delta2[i + 1],
    }


def _attach_endpoint_velocity(rows: list[dict[str, object]]) -> None:
    ordered = sorted(range(len(rows)), key=lambda index: float(rows[index]["a"]))
    if len(ordered) < 2:
        if ordered:
            rows[ordered[0]]["v_lambda"] = 0.0
        return
    for position, index in enumerate(ordered):
        if position == 0:
            left = ordered[position]
            right = ordered[position + 1]
        elif position == len(ordered) - 1:
            left = ordered[position - 1]
            right = ordered[position]
        else:
            left = ordered[position - 1]
            right = ordered[position + 1]
        k_left = float(rows[left]["k_NL_hmpc"])
        k_right = float(rows[right]["k_NL_hmpc"])
        a_left = float(rows[left]["a"])
        a_right = float(rows[right]["a"])
        rows[index]["v_lambda"] = (math.log(k_right) - math.log(k_left)) / (math.log(a_right) - math.log(a_left))


def _transfer_budget(
    transfer_rows: list[dict[str, object]],
    object_rows: list[dict[str, object]],
    omega_tail: float,
    tail_moment: float,
    omega_floor: float,
) -> dict[str, float]:
    strengths = [float(row["a_pq"]) for row in transfer_rows if row.get("a_pq") not in ("", None)]
    weights = [float(row.get("branch_weight", 1.0)) for row in object_rows if row.get("branch_weight") not in ("", None)]
    if strengths:
        injection_s0 = sum(strengths) / len(strengths)
    else:
        injection_s0 = 0.0
    if strengths and weights:
        weight_mean = sum(weights) / len(weights)
        injection_s1 = injection_s0 * weight_mean
    else:
        injection_s1 = 0.0
    drain_s0 = max(0.0, omega_tail - tail_moment)
    drain_s1 = max(0.0, tail_moment - omega_floor)
    return {
        "i_t_s0": injection_s0,
        "i_t_s1": injection_s1,
        "d_t_s0": drain_s0,
        "d_t_s1": drain_s1,
    }


def _closure_audit(
    *,
    epoch,
    gaussian_sigma_mpch: float,
    k_resolved_hmpc: float,
    activity_tail_width: float,
    epsilon0: float,
) -> dict[str, float | str]:
    delta2 = dimensionless_power(epoch.k_hmpc, epoch.p_mpc3_h3)
    coarse = [
        value * _gaussian_weight(k_value, gaussian_sigma_mpch) ** 2
        for k_value, value in zip(epoch.k_hmpc, delta2)
    ]
    activity = [
        _activity_weight(k_value, epoch.knl.k_nl_hmpc, k_resolved_hmpc, activity_tail_width)
        for k_value in epoch.k_hmpc
    ]
    active_norm = _masked_log_integral(
        epoch.k_hmpc,
        coarse,
        lambda k_value: 0.0 < k_value <= min(epoch.knl.k_nl_hmpc, k_resolved_hmpc),
    )
    tail_model_norm = _masked_log_integral(
        epoch.k_hmpc,
        [value * weight**2 for value, weight in zip(coarse, activity)],
        lambda k_value: epoch.knl.k_nl_hmpc < k_value <= k_resolved_hmpc,
    )
    total_norm = _masked_log_integral(epoch.k_hmpc, coarse, lambda k_value: k_value > 0.0)
    residual_norm = max(0.0, total_norm - active_norm - tail_model_norm)
    denom = active_norm + tail_model_norm + epsilon0
    closure_error = residual_norm / denom if denom > 0.0 else 0.0
    return {
        "active_norm": active_norm,
        "tail_model_norm": tail_model_norm,
        "total_norm": total_norm,
        "residual_norm": residual_norm,
        "closure_error": closure_error,
        "quality_flag": "closure_audited" if closure_error <= 0.1 else "residual_tracked",
    }


def _gaussian_weight(k_hmpc: float, sigma_mpch: float) -> float:
    if sigma_mpch <= 0.0:
        return 1.0
    return math.exp(-0.5 * (sigma_mpch * k_hmpc) ** 2)


def _activity_weight(k_hmpc: float, k_nl_hmpc: float, k_resolved_hmpc: float, tail_width: float) -> float:
    if k_hmpc <= 0.0 or k_hmpc > k_resolved_hmpc:
        return 0.0
    if k_hmpc <= k_nl_hmpc:
        return 1.0
    return math.exp(-0.5 * (math.log(k_hmpc / max(k_nl_hmpc, 1.0e-30)) / tail_width) ** 2)


def _masked_log_integral(k_values, values, predicate) -> float:
    pairs = [(math.log(k), value) for k, value in zip(k_values, values) if predicate(k)]
    if len(pairs) < 2:
        return 0.0
    return sum(0.5 * (right[1] + left[1]) * (right[0] - left[0]) for left, right in zip(pairs, pairs[1:]))


def _halo_selection_rows(object_rows: list[dict[str, object]], global_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    v_by_z = {float(row["z"]): float(row["v_vir_MNL_kms"]) for row in global_rows}
    m_by_z = {float(row["z"]): float(row["M_NL_msun_h"]) for row in global_rows}
    rows = []
    for row in object_rows:
        z = float(row["z"])
        rows.append(
            {
                "run_id": row["run_id"],
                "source": row["source"],
                "snapshot": row["snapshot"],
                "z": z,
                "halo_id": row["halo_id"],
                "mass_msun_h": row["mass_msun_h"],
                "radius_kpc_h": row["radius_kpc_h"],
                "sigma_v_kms": row["sigma_v_kms"],
                "vmax_kms": "",
                "sigma_source": row["sigma_source"],
                "sigma_used_kms": row["sigma_v_kms"],
                "mass_over_m_nl": float(row["mass_msun_h"]) / m_by_z[z],
                "r_h_object": float(row["sigma_v_kms"]) / v_by_z[z],
                "selected_for_r_h": float(row["mass_msun_h"]) <= m_by_z[z],
                "quality_flag": row["quality_flag"],
            }
        )
    return rows


def _merger_rows(merger_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    rows = []
    for row in merger_rows:
        mass_from = float(row["mass_from_msun_h"])
        mass_to = float(row["mass_to_msun_h"])
        time_from = float(row["time_from_gyr"])
        time_to = float(row["time_to_gyr"])
        delta_time = time_to - time_from
        growth_rate = "" if delta_time <= 0 or mass_from <= 0 or mass_to <= 0 else math.log(mass_to / mass_from) / delta_time
        rows.append(
            {
                "source": row["source"],
                "halo_id_from": row["halo_id_from"],
                "halo_id_to": row["halo_id_to"],
                "snapshot_from": row["snapshot_from"],
                "snapshot_to": row["snapshot_to"],
                "z_from": row["z_from"],
                "z_to": row["z_to"],
                "mass_from_msun_h": mass_from,
                "mass_to_msun_h": mass_to,
                "delta_mass_msun_h": mass_to - mass_from,
                "time_from_gyr": time_from,
                "time_to_gyr": time_to,
                "delta_time_gyr": delta_time,
                "growth_rate_log_per_gyr": growth_rate,
                "quality_flag": row.get("quality_flag", "internal_generation"),
            }
        )
    return rows


def _median_log_spacing(k_hmpc) -> float:
    values = sorted(math.log(k1 / k0) for k0, k1 in zip(k_hmpc, k_hmpc[1:]) if k1 > k0)
    mid = len(values) // 2
    if len(values) % 2:
        return values[mid]
    return 0.5 * (values[mid - 1] + values[mid])


def _tag(value: float) -> str:
    return f"{value:g}".replace("-", "m").replace(".", "p")


if __name__ == "__main__":
    raise SystemExit(main())
