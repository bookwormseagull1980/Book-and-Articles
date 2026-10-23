"""Figure-source table exports derived from stored result tables."""

from __future__ import annotations

from pathlib import Path

from cosmic_cg.io.tables import read_rows, write_rows


FIGURE1_COLUMNS = (
    "run_id",
    "component",
    "table_name",
    "row_count",
    "status",
)

FIGURE2_COLUMNS = (
    "run_id",
    "source",
    "snapshot",
    "z",
    "a",
    "gaussian_sigma_mpch",
    "activity_endpoint_hmpc",
    "c_nl",
    "k_box_hmpc",
    "k_nl_hmpc",
    "m_nl_msunh",
    "r_nl_mpch",
    "v_lambda",
    "s_abs",
    "s_aw",
    "k_cut_hmpc",
    "d_s",
    "closure_error",
    "endpoint_status",
)

FIGURE3_COLUMNS = (
    "run_id",
    "source",
    "snapshot",
    "z",
    "a",
    "band_id",
    "k_min_hmpc",
    "k_max_hmpc",
    "r_min_mpch",
    "r_max_mpch",
    "m_min_msunh",
    "m_max_msunh",
    "representative_structure",
    "omega_band",
    "activity_status",
    "channel_quality",
)

FIGURE4_COLUMNS = (
    "run_id",
    "source",
    "snapshot",
    "z",
    "a",
    "omega_tail_p",
    "r_v",
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
    "tail_channel",
    "q_t",
)

FIGURE5_COLUMNS = (
    "run_id",
    "source",
    "snapshot",
    "z_bin",
    "mass_bin",
    "n_objects",
    "b_h_median",
    "b_h_p16",
    "b_h_p84",
    "forming_fraction",
    "balanced_fraction",
    "stripping_fraction",
    "r_h_median",
    "selection",
    "quality_flag",
)

FIGURE6_COLUMNS = (
    "run_id",
    "source",
    "snapshot",
    "z_bin",
    "source_channel",
    "mass_ratio_bin",
    "environment_bin",
    "n_pairs",
    "a_pq_median",
    "a_pq_p16",
    "a_pq_p84",
    "g_pq_median",
    "g_pq_p16",
    "g_pq_p84",
    "quality_flag",
)

FIGURE7_COLUMNS = (
    "run_id",
    "source",
    "snapshot",
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
    "rho_total_measured",
)

FIGURE8_FIELD_COLUMNS = (
    "panel",
    "run_id",
    "source",
    "snapshot",
    "z",
    "a",
    "scale_band",
    "k_mid_hmpc",
    "gaussian_sigma_mpch",
    "activity_endpoint_hmpc",
    "power_fraction",
    "active_power_fraction",
    "tail_power_fraction",
    "cumulative_power_fraction",
    "closure_error",
    "object_count",
    "mass_median_msunh",
    "b_h_median",
    "transfer_count",
    "a_pq_median",
    "abs_g_pq_median",
    "endpoint_velocity",
    "relation_transition",
    "quality_flag",
)

LEGACY_FIGURE2_COLUMNS = (
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
)

LEGACY_FIGURE3_COLUMNS = (
    "run_id",
    "source",
    "snapshot",
    "z",
    "Omega_tail",
    "R_h_mean",
    "R_h_median",
    "R_h_p16",
    "R_h_p84",
    "R_v_optional",
    "Omega_floor",
    "Omega_excess",
    "tail_moment_s1",
    "tail_injection_s0",
    "tail_injection_s1",
    "tail_drain_s0",
    "tail_drain_s1",
    "quality_flag",
)

LEGACY_FIGURE4_COLUMNS = (
    "run_id",
    "source",
    "snapshot",
    "z",
    "halo_id",
    "mass_msun_h",
    "t_dyn_gyr",
    "t_growth_gyr",
    "B_h",
    "B_class",
    "quality_flag",
)


def write_standard_figure_data(run_dir: str | Path, output_dir: str | Path | None = None) -> list[Path]:
    root = Path(run_dir)
    out = Path(output_dir) if output_dir is not None else root / "figure_data"
    out.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []

    state_rows = _state_rows(root)
    if state_rows:
        path1 = out / "fig01_state_architecture.csv"
        write_rows(path1, state_rows, FIGURE1_COLUMNS)
        written.append(path1)

    window_path = root / "tables" / "window_table.csv"
    if window_path.is_file():
        path2 = out / "fig02_cosmic_active_window_history.csv"
        write_rows(path2, [_select(row, FIGURE2_COLUMNS) for row in read_rows(window_path)], FIGURE2_COLUMNS)
        written.append(path2)

    band_path = root / "tables" / "band_table.csv"
    if band_path.is_file():
        path3 = out / "fig03_scale_mass_bands.csv"
        write_rows(path3, [_select(row, FIGURE3_COLUMNS) for row in read_rows(band_path)], FIGURE3_COLUMNS)
        written.append(path3)

    tail_path = root / "tables" / "tail_table.csv"
    if tail_path.is_file():
        path4 = out / "fig04_residual_activity.csv"
        write_rows(path4, [_select(row, FIGURE4_COLUMNS) for row in read_rows(tail_path)], FIGURE4_COLUMNS)
        written.append(path4)

    object_summary_path = root / "tables" / "object_summary_table.csv"
    if object_summary_path.is_file():
        path5 = out / "fig05_object_fate.csv"
        write_rows(path5, [_select(row, FIGURE5_COLUMNS) for row in read_rows(object_summary_path)], FIGURE5_COLUMNS)
        written.append(path5)

    transfer_summary_path = root / "tables" / "transfer_summary_table.csv"
    if transfer_summary_path.is_file():
        path6 = out / "fig06_transfer_channels.csv"
        write_rows(path6, [_select(row, FIGURE6_COLUMNS) for row in read_rows(transfer_summary_path)], FIGURE6_COLUMNS)
        written.append(path6)

    cost_path = root / "tables" / "cost_table.csv"
    if cost_path.is_file():
        path7 = out / "fig07_cost_audit.csv"
        write_rows(path7, [_select(row, FIGURE7_COLUMNS) for row in read_rows(cost_path)], FIGURE7_COLUMNS)
        written.append(path7)

    global_path = root / "tables" / "global_window.csv"
    if global_path.is_file():
        rows = read_rows(global_path)
        path2 = out / "figure2_moving_endpoint.csv"
        path3 = out / "figure3_tail_activity.csv"
        write_rows(path2, [_select(row, LEGACY_FIGURE2_COLUMNS) for row in rows], LEGACY_FIGURE2_COLUMNS)
        write_rows(path3, [_select(row, LEGACY_FIGURE3_COLUMNS) for row in rows], LEGACY_FIGURE3_COLUMNS)
        written.extend([path2, path3])

    halo_path = root / "tables" / "halo_objects.csv"
    if halo_path.is_file():
        rows = read_rows(halo_path)
        path4 = out / "figure4_object_fate.csv"
        write_rows(path4, [_select(row, LEGACY_FIGURE4_COLUMNS) for row in rows], LEGACY_FIGURE4_COLUMNS)
        written.append(path4)

    return written


def write_field_figure_data(run_dir: str | Path, output_dir: str | Path | None = None) -> list[Path]:
    root = Path(run_dir)
    out = Path(output_dir) if output_dir is not None else root / "figure_data"
    out.mkdir(parents=True, exist_ok=True)
    rows = _field_figure_rows(root)
    if not rows:
        return []
    path = out / "fig08_field_relations.csv"
    write_rows(path, rows, FIGURE8_FIELD_COLUMNS)
    return [path]


def _select(row: dict[str, str], columns: tuple[str, ...]) -> dict[str, str]:
    return {column: row.get(column, "") for column in columns}


def _state_rows(root: Path) -> list[dict[str, object]]:
    components = (
        ("active field", "window_table.csv"),
        ("scale bands", "band_table.csv"),
        ("residual tail", "tail_table.csv"),
        ("generated objects", "object_table.csv"),
        ("transfer record", "transfer_table.csv"),
        ("cost audit", "cost_table.csv"),
    )
    rows = []
    run_id = root.name
    for component, table_name in components:
        path = root / "tables" / table_name
        if path.is_file():
            rows.append(
                {
                    "run_id": run_id,
                    "component": component,
                    "table_name": table_name,
                    "row_count": len(read_rows(path)),
                    "status": "recorded",
                }
            )
    return rows


def _field_figure_rows(root: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    mode_path = root / "tables" / "field_mode_table.csv"
    if mode_path.is_file():
        for row in read_rows(mode_path):
            rows.append(
                _field_figure_row(
                    panel="shell",
                    row=row,
                    scale_band=row.get("scale_band", ""),
                    k_mid_hmpc=row.get("k_mid_hmpc", ""),
                    gaussian_sigma_mpch=row.get("gaussian_sigma_mpch", ""),
                    activity_endpoint_hmpc=row.get("activity_endpoint_hmpc", ""),
                    power_fraction=row.get("power_fraction", ""),
                    active_power_fraction=row.get("active_power_fraction", ""),
                    tail_power_fraction=row.get("tail_power_fraction", ""),
                    cumulative_power_fraction=row.get("cumulative_power_fraction", ""),
                    quality_flag=row.get("activity_sector", ""),
                )
            )
    summary_path = root / "tables" / "field_relation_summary_table.csv"
    if summary_path.is_file():
        for row in read_rows(summary_path):
            rows.append(
                _field_figure_row(
                    panel="object_summary",
                    row=row,
                    scale_band=row.get("scale_band", ""),
                    object_count=row.get("object_count", ""),
                    mass_median_msunh=row.get("mass_median_msunh", ""),
                    b_h_median=row.get("b_h_median", ""),
                    transfer_count=row.get("transfer_count", ""),
                    a_pq_median=row.get("a_pq_median", ""),
                    abs_g_pq_median=row.get("abs_g_pq_median", ""),
                    quality_flag=row.get("relation_status", ""),
                )
            )
    closure_path = root / "tables" / "field_closure_table.csv"
    if closure_path.is_file():
        for row in read_rows(closure_path):
            rows.append(
                _field_figure_row(
                    panel="closure",
                    row=row,
                    gaussian_sigma_mpch=row.get("gaussian_sigma_mpch", ""),
                    activity_endpoint_hmpc=row.get("activity_endpoint_hmpc", ""),
                    closure_error=row.get("closure_error", ""),
                    quality_flag=row.get("quality_flag", ""),
                )
            )
    epoch_path = root / "tables" / "field_epoch_relation_table.csv"
    if epoch_path.is_file():
        for row in read_rows(epoch_path):
            rows.append(
                _field_figure_row(
                    panel="epoch_transition",
                    row=row,
                    endpoint_velocity=row.get("endpoint_velocity", ""),
                    relation_transition=row.get("relation_transition", ""),
                    quality_flag=row.get("quality_flag", ""),
                )
            )
    return rows


def _field_figure_row(panel: str, row: dict[str, str], **values: object) -> dict[str, object]:
    output = {column: "" for column in FIGURE8_FIELD_COLUMNS}
    output.update(
        {
            "panel": panel,
            "run_id": row.get("run_id", ""),
            "source": row.get("source", ""),
            "snapshot": row.get("snapshot", row.get("snapshot_to", "")),
            "z": row.get("z", row.get("z_to", "")),
            "a": row.get("a", row.get("a_to", "")),
        }
    )
    output.update(values)
    return output
