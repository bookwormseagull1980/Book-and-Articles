"""Completeness checks for run output directories."""

from __future__ import annotations

import csv
from dataclasses import asdict, dataclass
import json
from pathlib import Path
from typing import Any

from cosmic_cg.io.schema_catalog import SCHEMA_CATALOG_COLUMNS
from cosmic_cg.intermediate import (
    ENDPOINT_COLUMNS,
    HALO_SELECTION_COLUMNS,
    INPUT_VALIDATION_COLUMNS,
    MERGER_GROWTH_COLUMNS,
    OUTPUT_INVENTORY_COLUMNS,
    POWER_SAMPLE_COLUMNS,
)
from cosmic_cg.interactions import (
    CARRIER_HISTORY_COLUMNS,
    EMERGENCE_EVENT_COLUMNS,
    INTERACTION_BUDGET_COLUMNS,
    INTERACTION_CROSSING_COLUMNS,
    INTERACTION_DOMINANCE_COLUMNS,
    INTERACTION_STRENGTH_COLUMNS,
)
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
)


RUN_TYPES = {
    "generated_universe",
}

REQUIRED_FILES: dict[str, tuple[str, ...]] = {
    "generated_universe": (
        "config_resolved.yaml",
        "params_snapshot.json",
        "generation_manifest.json",
        "provenance.json",
        "intermediate/input_validation.csv",
        "intermediate/power_spectrum_samples.csv",
        "intermediate/endpoint_solution.csv",
        "intermediate/halo_selection.csv",
        "intermediate/merger_growth_links.csv",
        "intermediate/output_inventory.csv",
        "tables/global_window.csv",
        "tables/halo_objects.csv",
        "tables/window_table.csv",
        "tables/band_table.csv",
        "tables/tail_table.csv",
        "tables/closure_table.csv",
        "tables/object_table.csv",
        "tables/transfer_table.csv",
        "tables/cost_table.csv",
        "tables/event_table.csv",
        "tables/normalization_table.csv",
        "tables/object_summary_table.csv",
        "tables/transfer_summary_table.csv",
        "tables/interaction_strength_table.csv",
        "tables/interaction_budget_table.csv",
        "tables/interaction_dominance_table.csv",
        "tables/interaction_crossing_table.csv",
        "tables/carrier_history_table.csv",
        "tables/emergence_event_table.csv",
        "tables/schema_catalog.csv",
    ),
}

REQUIRED_COLUMNS: dict[str, tuple[str, ...]] = {
    "tables/global_window.csv": (
        "run_id",
        "source",
        "snapshot",
        "z",
        "a",
        "box_size_mpch",
        "k_box_hmpc",
        "gaussian_sigma_mpch",
        "activity_endpoint_hmpc",
        "C_NL",
        "k_NL_hmpc",
        "M_NL_msun_h",
        "R_NL_kpc_h",
        "v_vir_MNL_kms",
        "S_abs",
        "S_AW",
        "k_cut_hmpc",
        "D_S",
        "v_lambda",
        "Omega_tail",
        "Omega_floor",
        "Omega_excess",
        "R_h_mean",
        "R_h_median",
        "R_h_p16",
        "R_h_p84",
        "tail_moment_s1",
        "tail_injection_s0",
        "tail_injection_s1",
        "tail_drain_s0",
        "tail_drain_s1",
        "closure_active_norm",
        "closure_model_norm",
        "closure_residual_norm",
        "closure_error",
        "quality_flag",
    ),
    "tables/halo_objects.csv": (
        "run_id",
        "source",
        "snapshot",
        "z",
        "halo_id",
        "mass_msun_h",
        "radius_kpc_h",
        "v_vir_kms",
        "sigma_source",
        "x",
        "y",
        "z_pos",
        "vx",
        "vy",
        "vz",
        "t_dyn_gyr",
        "B_h",
        "B_class",
        "quality_flag",
    ),
    "tables/window_table.csv": WINDOW_TABLE_COLUMNS,
    "tables/band_table.csv": BAND_TABLE_COLUMNS,
    "tables/tail_table.csv": TAIL_TABLE_COLUMNS,
    "tables/closure_table.csv": CLOSURE_TABLE_COLUMNS,
    "tables/object_table.csv": OBJECT_TABLE_COLUMNS,
    "tables/transfer_table.csv": TRANSFER_TABLE_COLUMNS,
    "tables/interaction_strength_table.csv": INTERACTION_STRENGTH_COLUMNS,
    "tables/interaction_budget_table.csv": INTERACTION_BUDGET_COLUMNS,
    "tables/interaction_dominance_table.csv": INTERACTION_DOMINANCE_COLUMNS,
    "tables/interaction_crossing_table.csv": INTERACTION_CROSSING_COLUMNS,
    "tables/carrier_history_table.csv": CARRIER_HISTORY_COLUMNS,
    "tables/emergence_event_table.csv": EMERGENCE_EVENT_COLUMNS,
    "tables/cost_table.csv": COST_TABLE_COLUMNS,
    "tables/event_table.csv": EVENT_TABLE_COLUMNS,
    "tables/normalization_table.csv": NORMALIZATION_TABLE_COLUMNS,
    "tables/object_summary_table.csv": OBJECT_SUMMARY_TABLE_COLUMNS,
    "tables/transfer_summary_table.csv": TRANSFER_SUMMARY_TABLE_COLUMNS,
    "tables/schema_catalog.csv": SCHEMA_CATALOG_COLUMNS,
    "intermediate/input_validation.csv": INPUT_VALIDATION_COLUMNS,
    "intermediate/power_spectrum_samples.csv": POWER_SAMPLE_COLUMNS,
    "intermediate/endpoint_solution.csv": ENDPOINT_COLUMNS,
    "intermediate/halo_selection.csv": HALO_SELECTION_COLUMNS,
    "intermediate/merger_growth_links.csv": MERGER_GROWTH_COLUMNS,
    "intermediate/output_inventory.csv": OUTPUT_INVENTORY_COLUMNS,
    "figure_data/figure2_moving_endpoint.csv": (
        "run_id",
        "source",
        "snapshot",
        "z",
        "C_NL",
        "k_box_hmpc",
        "k_NL_hmpc",
        "S_abs",
        "S_AW",
        "k_cut_hmpc",
        "D_S",
        "quality_flag",
    ),
    "figure_data/figure3_tail_activity.csv": (
        "run_id",
        "source",
        "snapshot",
        "z",
        "Omega_tail",
        "R_h_mean",
        "R_h_median",
        "R_h_p16",
        "R_h_p84",
        "quality_flag",
    ),
    "figure_data/figure4_object_fate.csv": (
        "run_id",
        "source",
        "snapshot",
        "z",
        "halo_id",
        "mass_msun_h",
        "t_dyn_gyr",
        "B_h",
        "B_class",
        "quality_flag",
    ),
}

NUMERIC_COLUMNS: dict[str, tuple[str, ...]] = {
    "tables/global_window.csv": (
        "z",
        "a",
        "box_size_mpch",
        "k_box_hmpc",
        "gaussian_sigma_mpch",
        "activity_endpoint_hmpc",
        "C_NL",
        "k_NL_hmpc",
        "M_NL_msun_h",
        "R_NL_kpc_h",
        "v_vir_MNL_kms",
        "S_abs",
        "S_AW",
        "k_cut_hmpc",
        "D_S",
        "Omega_tail",
        "Omega_floor",
        "Omega_excess",
        "R_h_mean",
        "R_h_median",
        "R_h_p16",
        "R_h_p84",
        "v_lambda",
        "tail_moment_s1",
        "tail_injection_s0",
        "tail_injection_s1",
        "tail_drain_s0",
        "tail_drain_s1",
        "closure_active_norm",
        "closure_model_norm",
        "closure_residual_norm",
        "closure_error",
    ),
    "tables/halo_objects.csv": (
        "z",
        "mass_msun_h",
        "radius_kpc_h",
        "v_vir_kms",
        "x",
        "y",
        "z_pos",
        "vx",
        "vy",
        "vz",
        "t_dyn_gyr",
    ),
    "tables/window_table.csv": (
        "z",
        "a",
        "box_size_mpch",
        "k_box_hmpc",
        "gaussian_sigma_mpch",
        "activity_endpoint_hmpc",
        "c_nl",
        "k_nl_hmpc",
        "m_nl_msunh",
        "r_nl_mpch",
        "v_vir_mnl_kms",
        "s_abs",
        "s_aw",
        "k_cut_hmpc",
        "d_s",
        "v_lambda",
        "closure_active_norm",
        "closure_model_norm",
        "closure_residual_norm",
        "closure_error",
    ),
    "tables/band_table.csv": (
        "z",
        "a",
        "k_min_hmpc",
        "k_max_hmpc",
        "r_min_mpch",
        "r_max_mpch",
        "m_min_msunh",
        "m_max_msunh",
        "omega_band",
    ),
    "tables/tail_table.csv": (
        "z",
        "a",
        "omega_tail_p",
        "r_h",
        "e_t_s0",
        "e_t_s1",
        "i_t_s0",
        "i_t_s1",
        "d_t_s0",
        "d_t_s1",
        "omega_floor",
        "omega_excess",
        "closure_error",
    ),
    "tables/closure_table.csv": (
        "z",
        "a",
        "gaussian_sigma_mpch",
        "activity_endpoint_hmpc",
        "active_norm",
        "tail_model_norm",
        "residual_norm",
        "closure_error",
        "epsilon0",
    ),
    "tables/object_table.csv": (
        "z",
        "x",
        "y",
        "zpos",
        "vx",
        "vy",
        "vz",
        "mass_msunh",
        "r_vir_kpch",
        "b_h",
    ),
    "tables/transfer_table.csv": (
        "z",
        "a_pq",
        "g_pq",
        "tau_x",
        "tau_y",
        "tau_z",
        "j_x",
        "j_y",
        "j_z",
    ),
    "tables/interaction_strength_table.csv": (
        "z",
        "a",
        "source_term",
        "sink_term",
        "net_term",
        "raw_strength",
        "chi",
        "state_factor",
        "scale_factor",
        "carrier_factor",
    ),
    "tables/interaction_budget_table.csv": (
        "z",
        "a",
        "source_total",
        "sink_total",
        "net_total",
        "chi_total",
        "n_channels",
    ),
    "tables/interaction_dominance_table.csv": (
        "z",
        "a",
        "dominant_chi",
        "second_chi",
        "dominance_gap",
        "chi_sum",
    ),
    "tables/interaction_crossing_table.csv": (
        "z_left",
        "z_right",
        "a_left",
        "a_right",
    ),
    "tables/carrier_history_table.csv": (
        "z",
        "a",
        "presence",
        "source_total",
        "sink_total",
        "net_total",
    ),
    "tables/emergence_event_table.csv": (
        "z",
        "a",
        "value",
    ),
    "tables/cost_table.csv": (
        "z",
        "a",
        "lambda_active_hmpc",
        "k_ref_hmpc",
        "delta_s",
        "rho_n_pred",
        "rho_c_pred",
        "active_bytes",
        "residual_bytes",
        "audit_bytes",
        "total_bytes",
        "access_seconds",
        "audit_seconds",
    ),
    "tables/event_table.csv": ("z", "a", "k_nl_hmpc", "m_nl_msunh", "s_aw"),
    "tables/normalization_table.csv": ("z", "condition_value"),
    "tables/object_summary_table.csv": (
        "n_objects",
        "forming_fraction",
        "balanced_fraction",
        "stripping_fraction",
    ),
    "tables/transfer_summary_table.csv": ("n_pairs",),
    "intermediate/power_spectrum_samples.csv": ("z", "a", "k_hmpc", "p_mpc3_h3", "delta2", "k_over_knl", "tail_weight_s1"),
    "intermediate/endpoint_solution.csv": (
        "z",
        "c_nl",
        "k_nl_hmpc",
        "crossing_index",
        "k_left_hmpc",
        "k_right_hmpc",
        "delta2_left",
        "delta2_right",
    ),
    "intermediate/halo_selection.csv": (
        "z",
        "mass_msun_h",
        "radius_kpc_h",
        "sigma_used_kms",
        "mass_over_m_nl",
        "r_h_object",
    ),
    "intermediate/merger_growth_links.csv": (
        "z_from",
        "z_to",
        "mass_from_msun_h",
        "mass_to_msun_h",
        "delta_mass_msun_h",
        "time_from_gyr",
        "time_to_gyr",
        "delta_time_gyr",
    ),
    "intermediate/output_inventory.csv": ("bytes",),
    "figure_data/figure2_moving_endpoint.csv": (
        "z",
        "C_NL",
        "k_box_hmpc",
        "k_NL_hmpc",
        "S_abs",
        "S_AW",
        "k_cut_hmpc",
        "D_S",
    ),
    "figure_data/figure3_tail_activity.csv": (
        "z",
        "Omega_tail",
        "R_h_mean",
        "R_h_median",
        "R_h_p16",
        "R_h_p84",
    ),
    "figure_data/figure4_object_fate.csv": (
        "z",
        "mass_msun_h",
        "t_dyn_gyr",
    ),
}

JSON_FILES = {
    "params_snapshot.json",
    "generation_manifest.json",
    "provenance.json",
}


@dataclass(frozen=True)
class AuditIssue:
    severity: str
    code: str
    target: str
    message: str


@dataclass(frozen=True)
class AuditReport:
    run_dir: str
    run_type: str
    status: str
    checked_files: tuple[str, ...]
    checked_tables: tuple[str, ...]
    issues: tuple[AuditIssue, ...]

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def detect_run_type(run_dir: str | Path) -> str:
    root = Path(run_dir)
    if (root / "generation_manifest.json").is_file():
        return "generated_universe"
    return "unknown"


def audit_run(run_dir: str | Path, run_type: str | None = None) -> AuditReport:
    root = Path(run_dir)
    selected_type = detect_run_type(root) if run_type in (None, "auto") else str(run_type)
    issues: list[AuditIssue] = []
    checked_files: list[str] = []
    checked_tables: list[str] = []

    if selected_type not in RUN_TYPES:
        issues.append(_issue("error", "unknown_run_type", str(root), f"unknown run type: {selected_type}"))
        return _report(root, selected_type, checked_files, checked_tables, issues)

    for rel_path in REQUIRED_FILES[selected_type]:
        checked_files.append(rel_path)
        path = root / rel_path
        if not path.is_file():
            issues.append(_issue("error", "missing_file", rel_path, "required output file is missing"))
            continue
        if path.stat().st_size == 0:
            issues.append(_issue("error", "empty_file", rel_path, "required output file is empty"))
        if rel_path in JSON_FILES:
            _check_json(path, rel_path, root, issues)

    for rel_path in _required_tables(selected_type):
        checked_tables.append(rel_path)
        path = root / rel_path
        if path.is_file():
            _check_csv(path, rel_path, issues)

    _check_result_invariants(root, selected_type, issues)
    return _report(root, selected_type, checked_files, checked_tables, issues)


def write_audit_report(report: AuditReport, path: str | Path) -> None:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report.to_dict(), indent=2, sort_keys=True), encoding="utf-8")


def write_audit_markdown(report: AuditReport, path: str | Path) -> None:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    lines = [
        "# Run Audit Report",
        "",
        f"- Run directory: `{report.run_dir}`",
        f"- Run type: `{report.run_type}`",
        f"- Status: `{report.status}`",
        "",
        "## Checked Files",
        "",
    ]
    for rel_path in report.checked_files:
        lines.append(f"- `{rel_path}`")
    lines.extend(["", "## Checked Tables", ""])
    for rel_path in report.checked_tables:
        lines.append(f"- `{rel_path}`")
    lines.extend(["", "## Issues", ""])
    if report.issues:
        lines.extend(["| Severity | Code | Target | Message |", "|---|---|---|---|"])
        for issue in report.issues:
            lines.append(
                f"| {issue.severity} | {issue.code} | `{issue.target}` | {_escape_markdown_cell(issue.message)} |"
            )
    else:
        lines.append("No issues found.")
    lines.append("")
    output.write_text("\n".join(lines), encoding="utf-8")


def _required_tables(run_type: str) -> tuple[str, ...]:
    return tuple(rel_path for rel_path in REQUIRED_FILES[run_type] if rel_path.endswith(".csv"))


def _check_json(path: Path, rel_path: str, root: Path, issues: list[AuditIssue]) -> None:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        issues.append(_issue("error", "invalid_json", rel_path, f"cannot parse JSON: {exc}"))
        return

    if rel_path == "provenance.json":
        if not isinstance(payload, dict):
            issues.append(_issue("error", "invalid_provenance", rel_path, "provenance payload is not an object"))
            return
        if not payload.get("run_id"):
            issues.append(_issue("error", "missing_run_id", rel_path, "provenance does not record run_id"))
        elif str(payload["run_id"]) != root.name:
            issues.append(
                _issue("warning", "run_id_mismatch", rel_path, f"provenance run_id is {payload['run_id']!r}")
            )
        if not isinstance(payload.get("config"), dict):
            issues.append(_issue("warning", "missing_config", rel_path, "provenance does not include config object"))

    if rel_path == "generation_manifest.json":
        if not isinstance(payload, dict):
            issues.append(_issue("error", "invalid_generation_manifest", rel_path, "generation manifest is not an object"))
        else:
            for key in ("model", "redshifts", "object_count", "transfer_count", "spectrum", "inputs"):
                if key not in payload:
                    issues.append(_issue("error", "generation_manifest_missing_key", rel_path, f"missing key: {key}"))
            if payload.get("external_data") != []:
                issues.append(_issue("error", "external_data_not_empty", rel_path, "external_data must be an empty list"))
            state = payload.get("state_decomposition")
            if not isinstance(state, dict):
                issues.append(_issue("error", "missing_state_decomposition", rel_path, "state_decomposition is missing"))
            elif state.get("physical_coarse_graining") != "G_sigma" or state.get("activity_selector") != "A_Lambda":
                issues.append(_issue("error", "invalid_state_decomposition", rel_path, str(state)))
            closure = payload.get("closure")
            if not isinstance(closure, dict) or not closure.get("map"):
                issues.append(_issue("error", "missing_closure_manifest", rel_path, "closure map is missing"))


def _check_csv(path: Path, rel_path: str, issues: list[AuditIssue]) -> None:
    try:
        with path.open("r", newline="", encoding="utf-8-sig") as fh:
            reader = csv.DictReader(fh)
            header = reader.fieldnames or []
            rows = list(reader)
    except csv.Error as exc:
        issues.append(_issue("error", "invalid_csv", rel_path, f"cannot parse CSV: {exc}"))
        return

    missing = sorted(set(REQUIRED_COLUMNS.get(rel_path, ())) - set(header))
    if missing:
        issues.append(_issue("error", "missing_columns", rel_path, "missing columns: " + ", ".join(missing)))
    if not rows:
        issues.append(_issue("error", "empty_table", rel_path, "table contains no data rows"))

    for column in NUMERIC_COLUMNS.get(rel_path, ()):
        for index, row in enumerate(rows, start=1):
            value = row.get(column, "")
            if value == "":
                issues.append(_issue("warning", "missing_numeric_value", f"{rel_path}:{index}:{column}", "value is blank"))
                continue
            try:
                float(value)
            except ValueError:
                issues.append(_issue("error", "invalid_numeric_value", f"{rel_path}:{index}:{column}", value))


def _check_result_invariants(root: Path, run_type: str, issues: list[AuditIssue]) -> None:
    if run_type == "generated_universe":
        _check_generation_manifest_counts(root, issues)

    global_path = root / "tables" / "global_window.csv"
    if global_path.is_file():
        rows = _safe_read_rows(global_path)
        for index, row in enumerate(rows, start=1):
            _require_positive(row, "C_NL", f"tables/global_window.csv:{index}", issues)
            _require_positive(row, "k_NL_hmpc", f"tables/global_window.csv:{index}", issues)
            _require_nonnegative(row, "gaussian_sigma_mpch", f"tables/global_window.csv:{index}", issues)
            _require_positive(row, "activity_endpoint_hmpc", f"tables/global_window.csv:{index}", issues)
            _require_positive(row, "S_abs", f"tables/global_window.csv:{index}", issues)
            _require_nonnegative(row, "S_AW", f"tables/global_window.csv:{index}", issues)
            _require_nonnegative(row, "D_S", f"tables/global_window.csv:{index}", issues)
            _require_nonnegative(row, "closure_error", f"tables/global_window.csv:{index}", issues)

    closure_path = root / "tables" / "closure_table.csv"
    if closure_path.is_file():
        for index, row in enumerate(_safe_read_rows(closure_path), start=1):
            target = f"tables/closure_table.csv:{index}"
            _require_nonnegative(row, "gaussian_sigma_mpch", target, issues)
            _require_positive(row, "activity_endpoint_hmpc", target, issues)
            _require_nonnegative(row, "active_norm", target, issues)
            _require_nonnegative(row, "tail_model_norm", target, issues)
            _require_nonnegative(row, "residual_norm", target, issues)
            _require_nonnegative(row, "closure_error", target, issues)

    halo_path = root / "tables" / "halo_objects.csv"
    if run_type == "generated_universe" and halo_path.is_file():
        allowed_classes = {"forming", "balanced", "stripping"}
        rows = _safe_read_rows(halo_path)
        for index, row in enumerate(rows, start=1):
            if row.get("B_class") not in allowed_classes:
                issues.append(
                    _issue("error", "invalid_b_class", f"tables/halo_objects.csv:{index}:B_class", row.get("B_class", ""))
                )
            if row.get("B_class") != "stripping":
                _require_numeric(row, "B_h", f"tables/halo_objects.csv:{index}", issues)

    band_path = root / "tables" / "band_table.csv"
    if band_path.is_file():
        for index, row in enumerate(_safe_read_rows(band_path), start=1):
            target = f"tables/band_table.csv:{index}"
            k_min = _float_or_issue(row, "k_min_hmpc", target, issues)
            k_max = _float_or_issue(row, "k_max_hmpc", target, issues)
            r_min = _float_or_issue(row, "r_min_mpch", target, issues)
            r_max = _float_or_issue(row, "r_max_mpch", target, issues)
            m_min = _float_or_issue(row, "m_min_msunh", target, issues)
            m_max = _float_or_issue(row, "m_max_msunh", target, issues)
            omega = _float_or_issue(row, "omega_band", target, issues)
            if k_min is not None and k_max is not None and k_min >= k_max:
                issues.append(_issue("error", "invalid_band_order", f"{target}:k", "k_min_hmpc must be below k_max_hmpc"))
            if r_min is not None and r_max is not None and r_min >= r_max:
                issues.append(_issue("error", "invalid_band_order", f"{target}:r", "r_min_mpch must be below r_max_mpch"))
            if m_min is not None and m_max is not None and m_min >= m_max:
                issues.append(_issue("error", "invalid_band_order", f"{target}:m", "m_min_msunh must be below m_max_msunh"))
            if omega is not None and not (0.0 <= omega <= 1.0):
                issues.append(_issue("error", "fraction_out_of_range", f"{target}:omega_band", str(omega)))

    cost_path = root / "tables" / "cost_table.csv"
    if cost_path.is_file():
        for index, row in enumerate(_safe_read_rows(cost_path), start=1):
            target = f"tables/cost_table.csv:{index}"
            for column in ("rho_n_pred", "rho_c_pred"):
                value = _float_or_issue(row, column, target, issues)
                if value is not None and not (0.0 <= value <= 1.0):
                    issues.append(_issue("error", "fraction_out_of_range", f"{target}:{column}", str(value)))

    transfer_path = root / "tables" / "transfer_table.csv"
    if transfer_path.is_file():
        for index, row in enumerate(_safe_read_rows(transfer_path), start=1):
            target = f"tables/transfer_table.csv:{index}"
            a_pq = _float_or_issue(row, "a_pq", target, issues)
            g_pq = _float_or_issue(row, "g_pq", target, issues)
            if a_pq is not None and a_pq < 0.0:
                issues.append(_issue("error", "negative_value", f"{target}:a_pq", str(a_pq)))
            if g_pq is not None and not (-1.0 <= g_pq <= 1.0):
                issues.append(_issue("error", "alignment_out_of_range", f"{target}:g_pq", str(g_pq)))

    _check_interaction_tables(root, issues)

    summary_path = root / "tables" / "object_summary_table.csv"
    if summary_path.is_file():
        for index, row in enumerate(_safe_read_rows(summary_path), start=1):
            target = f"tables/object_summary_table.csv:{index}"
            total = 0.0
            complete = True
            for column in ("forming_fraction", "balanced_fraction", "stripping_fraction"):
                value = _float_or_issue(row, column, target, issues)
                if value is None:
                    complete = False
                    continue
                if not (0.0 <= value <= 1.0):
                    issues.append(_issue("error", "fraction_out_of_range", f"{target}:{column}", str(value)))
                total += value
            if complete and abs(total - 1.0) > 1.0e-6:
                issues.append(_issue("warning", "fraction_sum", target, f"object fractions sum to {total}"))


def _check_interaction_tables(root: Path, issues: list[AuditIssue]) -> None:
    strength_path = root / "tables" / "interaction_strength_table.csv"
    dominance_path = root / "tables" / "interaction_dominance_table.csv"
    if not strength_path.is_file() or not dominance_path.is_file():
        return
    strength_rows = _safe_read_rows(strength_path)
    dominance_rows = _safe_read_rows(dominance_path)
    groups = sorted({_interaction_key(row) for row in strength_rows})
    for key in groups:
        selected = [row for row in strength_rows if _interaction_key(row) == key]
        chi_sum = sum(float(row["chi"]) for row in selected)
        target = f"tables/interaction_strength_table.csv:{key}"
        if abs(chi_sum - 1.0) > 1.0e-8:
            issues.append(_issue("error", "interaction_chi_sum", target, f"chi values sum to {chi_sum}"))
        for row in selected:
            chi = float(row["chi"])
            raw = float(row["raw_strength"])
            if chi < 0.0 or chi > 1.0:
                issues.append(_issue("error", "interaction_chi_range", target, str(chi)))
            if raw < 0.0:
                issues.append(_issue("error", "interaction_raw_strength", target, str(raw)))
    dominance_by_key = {_interaction_key(row): row for row in dominance_rows}
    for key in groups:
        selected = [row for row in strength_rows if _interaction_key(row) == key]
        if key not in dominance_by_key:
            issues.append(_issue("error", "interaction_missing_dominance", key, "dominance row is missing"))
            continue
        strongest = max(selected, key=lambda row: float(row["chi"]))
        row = dominance_by_key[key]
        if row["dominant_channel"] != strongest["channel_id"]:
            issues.append(
                _issue(
                    "error",
                    "interaction_dominance_mismatch",
                    f"tables/interaction_dominance_table.csv:{key}",
                    f"dominant channel is {row['dominant_channel']} but strongest is {strongest['channel_id']}",
                )
            )
        chi_sum = float(row["chi_sum"])
        if abs(chi_sum - 1.0) > 1.0e-8:
            issues.append(_issue("error", "interaction_dominance_sum", key, f"chi_sum is {chi_sum}"))
        gap = float(row["dominance_gap"])
        if gap < -1.0e-12:
            issues.append(_issue("error", "negative_dominance_gap", key, str(gap)))


def _check_generation_manifest_counts(root: Path, issues: list[AuditIssue]) -> None:
    manifest_path = root / "generation_manifest.json"
    if not manifest_path.is_file():
        return
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return
    if not isinstance(manifest, dict):
        return
    checks = (
        ("object_count", root / "tables" / "object_table.csv"),
        ("transfer_count", root / "tables" / "transfer_table.csv"),
    )
    for key, table_path in checks:
        if key not in manifest or not table_path.is_file():
            continue
        expected = int(manifest[key])
        actual = len(_safe_read_rows(table_path))
        if expected != actual:
            issues.append(
                _issue(
                    "error",
                    "manifest_count_mismatch",
                    f"generation_manifest.json:{key}",
                    f"manifest records {expected} rows but {table_path.name} contains {actual}",
                )
            )
    redshifts = manifest.get("redshifts")
    window_path = root / "tables" / "window_table.csv"
    if isinstance(redshifts, list) and window_path.is_file():
        actual = len(_safe_read_rows(window_path))
        if len(redshifts) != actual:
            issues.append(
                _issue(
                    "error",
                    "manifest_redshift_mismatch",
                    "generation_manifest.json:redshifts",
                    f"manifest records {len(redshifts)} epochs but window_table.csv contains {actual}",
                )
            )


def _safe_read_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8-sig") as fh:
        return list(csv.DictReader(fh))


def _interaction_key(row: dict[str, str]) -> str:
    return f"{float(row['z']):.12g}|{float(row['a']):.12g}|{row['band_id']}"


def _require_numeric(row: dict[str, str], column: str, target: str, issues: list[AuditIssue]) -> None:
    value = row.get(column, "")
    if value == "":
        issues.append(_issue("error", "missing_value", f"{target}:{column}", "value is blank"))
        return
    try:
        float(value)
    except ValueError:
        issues.append(_issue("error", "invalid_numeric_value", f"{target}:{column}", value))


def _require_positive(row: dict[str, str], column: str, target: str, issues: list[AuditIssue]) -> None:
    value = _float_or_issue(row, column, target, issues)
    if value is not None and value <= 0.0:
        issues.append(_issue("error", "nonpositive_value", f"{target}:{column}", str(value)))


def _require_nonnegative(row: dict[str, str], column: str, target: str, issues: list[AuditIssue]) -> None:
    value = _float_or_issue(row, column, target, issues)
    if value is not None and value < 0.0:
        issues.append(_issue("error", "negative_value", f"{target}:{column}", str(value)))


def _float_or_issue(row: dict[str, str], column: str, target: str, issues: list[AuditIssue]) -> float | None:
    value = row.get(column, "")
    if value == "":
        issues.append(_issue("error", "missing_value", f"{target}:{column}", "value is blank"))
        return None
    try:
        return float(value)
    except ValueError:
        issues.append(_issue("error", "invalid_numeric_value", f"{target}:{column}", value))
        return None


def _report(
    root: Path,
    run_type: str,
    checked_files: list[str],
    checked_tables: list[str],
    issues: list[AuditIssue],
) -> AuditReport:
    if any(issue.severity == "error" for issue in issues):
        status = "fail"
    elif issues:
        status = "warning"
    else:
        status = "pass"
    return AuditReport(
        run_dir=str(root),
        run_type=run_type,
        status=status,
        checked_files=tuple(checked_files),
        checked_tables=tuple(checked_tables),
        issues=tuple(issues),
    )


def _issue(severity: str, code: str, target: str, message: str) -> AuditIssue:
    return AuditIssue(severity=severity, code=code, target=target, message=message)


def _escape_markdown_cell(value: str) -> str:
    return value.replace("|", "\\|")
