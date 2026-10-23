"""Adapters from generated run tables to interaction ledgers."""

from __future__ import annotations

from collections import defaultdict
import json
import math
from pathlib import Path
from typing import Mapping, Sequence

from cosmic_cg.interactions import (
    CARRIER_HISTORY_COLUMNS,
    EMERGENCE_EVENT_COLUMNS,
    INTERACTION_BUDGET_COLUMNS,
    INTERACTION_CROSSING_COLUMNS,
    INTERACTION_DOMINANCE_COLUMNS,
    INTERACTION_STRENGTH_COLUMNS,
    generate_interaction_outputs,
)
from cosmic_cg.io.tables import read_rows, write_rows
from cosmic_cg.result_tables import DEFAULT_MASS_BANDS


INTERACTION_TABLES = {
    "interaction_strength_table": ("interaction_strength_table.csv", INTERACTION_STRENGTH_COLUMNS, "strength"),
    "interaction_budget_table": ("interaction_budget_table.csv", INTERACTION_BUDGET_COLUMNS, "budget"),
    "interaction_dominance_table": ("interaction_dominance_table.csv", INTERACTION_DOMINANCE_COLUMNS, "dominance"),
    "interaction_crossing_table": ("interaction_crossing_table.csv", INTERACTION_CROSSING_COLUMNS, "crossing"),
    "carrier_history_table": ("carrier_history_table.csv", CARRIER_HISTORY_COLUMNS, "carrier"),
    "emergence_event_table": ("emergence_event_table.csv", EMERGENCE_EVENT_COLUMNS, "emergence"),
}

FIXED_BAND_RANGES_HMPC = {
    "dwarf_halo": (1.0, 20.0),
    "galaxy_halo": (0.1, 1.0),
    "group_halo": (0.03, 0.1),
    "cluster_halo": (0.01, 0.03),
    "large_environment": (0.001, 0.01),
}


def write_field_interaction_tables(run_root: str | Path) -> dict[str, int]:
    root = Path(run_root)
    tables = root / "tables"
    outputs = field_interaction_outputs(
        window_rows=read_rows(tables / "field_window_table.csv"),
        mode_rows=read_rows(tables / "field_mode_table.csv"),
        object_rows=read_rows(tables / "field_object_table.csv"),
        transfer_rows=read_rows(tables / "field_transfer_table.csv"),
        epoch_rows=read_rows(tables / "field_epoch_relation_table.csv"),
    )
    counts = _write_interaction_tables(tables, outputs)
    _update_manifest_counts(root, counts)
    _refresh_cost_table(root)
    return counts


def write_pm_interaction_tables(run_root: str | Path) -> dict[str, int]:
    root = Path(run_root)
    tables = root / "tables"
    outputs = pm_interaction_outputs(
        snapshot_rows=read_rows(tables / "pm_snapshot_table.csv"),
        power_rows=read_rows(tables / "pm_power_table.csv"),
        halo_rows=read_rows(tables / "pm_halo_table.csv"),
        transfer_rows=read_rows(tables / "pm_transfer_table.csv"),
    )
    counts = _write_interaction_tables(tables, outputs)
    _update_manifest_counts(root, counts)
    _refresh_cost_table(root)
    return counts


def field_interaction_outputs(
    *,
    window_rows: Sequence[Mapping[str, object]],
    mode_rows: Sequence[Mapping[str, object]],
    object_rows: Sequence[Mapping[str, object]],
    transfer_rows: Sequence[Mapping[str, object]],
    epoch_rows: Sequence[Mapping[str, object]],
) -> dict[str, list[dict[str, object]]]:
    run_id = str(window_rows[0]["run_id"]) if window_rows else ""
    global_rows = _field_global_rows(window_rows, epoch_rows)
    band_rows = _band_rows_from_modes(window_rows, mode_rows, power_key="power_fraction")
    return generate_interaction_outputs(
        run_id=run_id,
        global_rows=global_rows,
        band_rows=band_rows,
        object_rows=object_rows,
        transfer_rows=transfer_rows,
    )


def pm_interaction_outputs(
    *,
    snapshot_rows: Sequence[Mapping[str, object]],
    power_rows: Sequence[Mapping[str, object]],
    halo_rows: Sequence[Mapping[str, object]],
    transfer_rows: Sequence[Mapping[str, object]],
) -> dict[str, list[dict[str, object]]]:
    run_id = str(snapshot_rows[0]["run_id"]) if snapshot_rows else ""
    global_rows = _pm_global_rows(snapshot_rows, power_rows, halo_rows)
    band_rows = _band_rows_from_modes(snapshot_rows, power_rows, power_key="power_fraction")
    object_rows = [
        {
            "z": row["z"],
            "a": row["a"],
            "mass_msunh": row["mass_msunh"],
            "b_h": row["b_h"],
        }
        for row in halo_rows
    ]
    return generate_interaction_outputs(
        run_id=run_id,
        global_rows=global_rows,
        band_rows=band_rows,
        object_rows=object_rows,
        transfer_rows=transfer_rows,
    )


def write_interaction_tables_for_run(run_root: str | Path) -> dict[str, int]:
    root = Path(run_root)
    tables = root / "tables"
    if (tables / "field_window_table.csv").is_file():
        return write_field_interaction_tables(root)
    if (tables / "pm_snapshot_table.csv").is_file():
        return write_pm_interaction_tables(root)
    raise ValueError(f"unsupported run table layout: {root}")


def _field_global_rows(
    window_rows: Sequence[Mapping[str, object]],
    epoch_rows: Sequence[Mapping[str, object]],
) -> list[dict[str, object]]:
    velocity_by_snapshot = {}
    for row in epoch_rows:
        velocity = abs(_float(row, "endpoint_velocity", 0.0))
        velocity_by_snapshot[str(row["snapshot_from"])] = velocity
        velocity_by_snapshot[str(row["snapshot_to"])] = velocity
    output = []
    for row in window_rows:
        omega_tail = _float(row, "tail_power_fraction", 0.0)
        closure_error = _float(row, "closure_error", 0.0)
        output.append(
            {
                "run_id": row["run_id"],
                "source": row["source"],
                "snapshot": row["snapshot"],
                "z": row["z"],
                "a": row["a"],
                "k_NL_hmpc": row["k_nl_hmpc"],
                "M_NL_msun_h": row["m_nl_msunh"],
                "Omega_tail": omega_tail,
                "Omega_excess": closure_error / (1.0 + closure_error) if closure_error >= 0.0 else 0.0,
                "R_h_mean": _float(row, "measured_sigma_delta", 0.0),
                "v_lambda": velocity_by_snapshot.get(str(row["snapshot"]), 0.0),
            }
        )
    return output


def _pm_global_rows(
    snapshot_rows: Sequence[Mapping[str, object]],
    power_rows: Sequence[Mapping[str, object]],
    halo_rows: Sequence[Mapping[str, object]],
) -> list[dict[str, object]]:
    powers_by_snapshot = _group_by(power_rows, "snapshot")
    halos_by_snapshot = _group_by(halo_rows, "snapshot")
    endpoint_by_snapshot = {}
    for snapshot, selected in powers_by_snapshot.items():
        endpoint_by_snapshot[snapshot] = _power_endpoint(selected)
    velocities = _endpoint_velocities(snapshot_rows, endpoint_by_snapshot)
    output = []
    for row in snapshot_rows:
        snapshot = str(row["snapshot"])
        selected_power = powers_by_snapshot.get(snapshot, [])
        selected_halos = halos_by_snapshot.get(snapshot, [])
        endpoint = endpoint_by_snapshot.get(snapshot, 0.2)
        tail = _tail_power_fraction(selected_power, endpoint)
        output.append(
            {
                "run_id": row["run_id"],
                "source": row["source"],
                "snapshot": row["snapshot"],
                "z": row["z"],
                "a": row["a"],
                "k_NL_hmpc": endpoint,
                "M_NL_msun_h": _median([_float(item, "mass_msunh", 0.0) for item in selected_halos])
                or 12.0 * _float(row, "particle_mass_msunh", 1.0),
                "Omega_tail": tail,
                "Omega_excess": _high_k_fraction(selected_power),
                "R_h_mean": _median([abs(_float(item, "b_h", 0.0)) for item in selected_halos]),
                "v_lambda": velocities.get(snapshot, 0.0),
            }
        )
    return output


def _band_rows_from_modes(
    state_rows: Sequence[Mapping[str, object]],
    mode_rows: Sequence[Mapping[str, object]],
    *,
    power_key: str,
) -> list[dict[str, object]]:
    modes_by_snapshot = _group_by(mode_rows, "snapshot")
    output = []
    for row in state_rows:
        snapshot = str(row["snapshot"])
        selected_modes = modes_by_snapshot.get(snapshot, [])
        endpoint = _float(row, "k_nl_hmpc", _float(row, "k_NL_hmpc", 0.2))
        for band_id, mass_min, mass_max, label in DEFAULT_MASS_BANDS:
            k_min, k_max = FIXED_BAND_RANGES_HMPC[band_id]
            output.append(
                {
                    "run_id": row["run_id"],
                    "source": row["source"],
                    "snapshot": row["snapshot"],
                    "z": row["z"],
                    "a": row["a"],
                    "band_id": band_id,
                    "k_min_hmpc": k_min,
                    "k_max_hmpc": k_max,
                    "m_min_msunh": mass_min,
                    "m_max_msunh": mass_max,
                    "omega_band": _band_power_fraction(selected_modes, k_min, k_max, power_key),
                    "activity_status": _activity_status(k_min, k_max, endpoint),
                    "representative_structure": label,
                }
            )
    return output


def _write_interaction_tables(tables_dir: Path, outputs: Mapping[str, Sequence[Mapping[str, object]]]) -> dict[str, int]:
    counts = {}
    for table_key, (filename, columns, output_key) in INTERACTION_TABLES.items():
        rows = list(outputs[output_key])
        write_rows(tables_dir / filename, rows, columns)
        counts[table_key] = len(rows)
    return counts


def _update_manifest_counts(root: Path, counts: Mapping[str, int]) -> None:
    path = root / "generation_manifest.json"
    if not path.is_file():
        return
    manifest = json.loads(path.read_text(encoding="utf-8"))
    output_tables = dict(manifest.get("output_tables", {}))
    output_tables.update(counts)
    manifest["output_tables"] = output_tables
    path.write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")


def _refresh_cost_table(root: Path) -> None:
    manifest_path = root / "generation_manifest.json"
    if not manifest_path.is_file():
        return
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    output_tables = manifest.get("output_tables", {})
    relationship_rows = sum(
        int(value)
        for key, value in output_tables.items()
        if not str(key).endswith("_cost_table")
    )
    relationship_bytes = _directory_size_bytes(root / "tables")
    model = manifest.get("model")
    if model == "finite_window_3d_field_generation":
        from cosmic_cg.field3d import FIELD_COST_COLUMNS, Field3DConfig, field_cost_row

        config = Field3DConfig(**manifest["field_config"])
        row = field_cost_row(
            run_id=str(manifest["run_id"]),
            field_config=config,
            redshift_count=len(manifest.get("redshifts", [])),
            relationship_rows=relationship_rows,
            relationship_bytes_written=relationship_bytes,
        )
        write_rows(root / "tables" / "field_cost_table.csv", [row], FIELD_COST_COLUMNS)
    elif model == "finite_window_internal_pm_nbody":
        from cosmic_cg.pm_nbody import PM_COST_COLUMNS, PMNBodyConfig, pm_cost_row

        config = PMNBodyConfig(**manifest["pm_config"])
        row = pm_cost_row(
            run_id=str(manifest["run_id"]),
            config=config,
            redshift_count=len(manifest.get("redshifts", [])),
            relationship_rows=relationship_rows,
            relationship_bytes_written=relationship_bytes,
        )
        write_rows(root / "tables" / "pm_cost_table.csv", [row], PM_COST_COLUMNS)


def _directory_size_bytes(root: Path) -> int:
    return sum(path.stat().st_size for path in root.rglob("*") if path.is_file())


def _group_by(rows: Sequence[Mapping[str, object]], key: str) -> dict[str, list[Mapping[str, object]]]:
    grouped: dict[str, list[Mapping[str, object]]] = defaultdict(list)
    for row in rows:
        grouped[str(row[key])].append(row)
    return dict(grouped)


def _power_endpoint(rows: Sequence[Mapping[str, object]]) -> float:
    ordered = sorted(rows, key=lambda row: _float(row, "k_mid_hmpc", 0.0))
    if not ordered:
        return 0.2
    for row in ordered:
        if _float(row, "cumulative_power_fraction", 0.0) >= 0.5:
            return _float(row, "k_mid_hmpc", 0.2)
    return _float(ordered[-1], "k_mid_hmpc", 0.2)


def _endpoint_velocities(
    snapshot_rows: Sequence[Mapping[str, object]],
    endpoint_by_snapshot: Mapping[str, float],
) -> dict[str, float]:
    ordered = sorted(snapshot_rows, key=lambda row: _float(row, "a", 0.0))
    output: dict[str, float] = {}
    for left, right in zip(ordered, ordered[1:]):
        left_k = max(1.0e-12, endpoint_by_snapshot.get(str(left["snapshot"]), 0.2))
        right_k = max(1.0e-12, endpoint_by_snapshot.get(str(right["snapshot"]), 0.2))
        delta_ln_a = math.log(max(1.0e-12, _float(right, "a", 0.0))) - math.log(max(1.0e-12, _float(left, "a", 0.0)))
        value = abs(math.log(right_k / left_k) / delta_ln_a) if abs(delta_ln_a) > 1.0e-12 else 0.0
        output[str(left["snapshot"])] = value
        output[str(right["snapshot"])] = value
    return output


def _band_power_fraction(
    rows: Sequence[Mapping[str, object]],
    k_min: float,
    k_max: float,
    power_key: str,
) -> float:
    total = 0.0
    for row in rows:
        k_mid = _float(row, "k_mid_hmpc", 0.0)
        if k_min <= k_mid < k_max:
            total += max(0.0, _float(row, power_key, 0.0))
    return total


def _tail_power_fraction(rows: Sequence[Mapping[str, object]], endpoint: float) -> float:
    return sum(max(0.0, _float(row, "power_fraction", 0.0)) for row in rows if _float(row, "k_mid_hmpc", 0.0) > endpoint)


def _high_k_fraction(rows: Sequence[Mapping[str, object]]) -> float:
    if not rows:
        return 0.0
    k_max = max(_float(row, "k_mid_hmpc", 0.0) for row in rows)
    threshold = 0.5 * k_max
    return sum(max(0.0, _float(row, "power_fraction", 0.0)) for row in rows if _float(row, "k_mid_hmpc", 0.0) >= threshold)


def _activity_status(k_min: float, k_max: float, endpoint: float) -> str:
    if k_min <= endpoint <= k_max:
        return "frontier"
    if k_max < endpoint:
        return "active"
    return "tail"


def _median(values: Sequence[float]) -> float:
    selected = sorted(value for value in values if math.isfinite(value) and value > 0.0)
    if not selected:
        return 0.0
    mid = len(selected) // 2
    if len(selected) % 2:
        return selected[mid]
    return 0.5 * (selected[mid - 1] + selected[mid])


def _float(row: Mapping[str, object], key: str, default: float) -> float:
    try:
        value = row.get(key, default)
        if value in ("", None):
            return default
        return float(value)
    except (TypeError, ValueError):
        return default
