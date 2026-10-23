"""Summarize a completed three-dimensional field-generation run."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import sys
from typing import Iterable, Mapping, Sequence

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from cosmic_cg.intermediate import write_output_inventory
from cosmic_cg.io.tables import read_rows, write_rows


KEY_EPOCH_COLUMNS = (
    "run_id",
    "z",
    "a",
    "k_nl_hmpc",
    "lambda_nl_mpch",
    "m_nl_msunh",
    "s_aw_global",
    "measured_sigma_delta",
    "active_power_fraction",
    "tail_power_fraction",
    "peak_count",
    "transfer_count",
    "object_mass_median_msunh",
    "object_radius_median_kpch",
    "b_h_median",
    "forming_fraction",
    "balanced_fraction",
    "stripping_fraction",
    "a_pq_median",
    "abs_g_pq_median",
    "dominant_scale_band",
)

STAGE_COLUMNS = (
    "run_id",
    "stage_marker",
    "z",
    "a",
    "selected_by",
    "value",
    "description",
)


def main() -> int:
    args = parse_args()
    root = Path(args.run_dir)
    analysis_dir = root / "analysis"
    analysis_dir.mkdir(parents=True, exist_ok=True)

    tables = _load_tables(root)
    key_epochs = _key_epoch_rows(root.name, tables, _requested_redshifts(args.key_redshifts))
    stage_rows = _stage_rows(root.name, tables)
    summary = _summary_payload(root, tables, key_epochs, stage_rows)

    write_rows(analysis_dir / "field_key_epoch_table.csv", key_epochs, KEY_EPOCH_COLUMNS)
    write_rows(analysis_dir / "field_stage_markers.csv", stage_rows, STAGE_COLUMNS)
    (analysis_dir / "field_result_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False, sort_keys=True),
        encoding="utf-8",
    )
    (analysis_dir / "field_result_interpretation.md").write_text(
        _interpretation_markdown(root.name, summary, key_epochs, stage_rows),
        encoding="utf-8",
    )
    write_output_inventory(root)
    print(json.dumps({"run_dir": str(root), "analysis_dir": str(analysis_dir), "status": "written"}, indent=2))
    return 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True)
    parser.add_argument("--key-redshifts", default="100,10,3,1,0,-0.666667")
    return parser.parse_args()


def _load_tables(root: Path) -> dict[str, list[dict[str, str]]]:
    names = {
        "window": "field_window_table.csv",
        "mode": "field_mode_table.csv",
        "object": "field_object_table.csv",
        "transfer": "field_transfer_table.csv",
        "summary": "field_relation_summary_table.csv",
        "epoch": "field_epoch_relation_table.csv",
        "cost": "field_cost_table.csv",
    }
    return {key: read_rows(root / "tables" / name) for key, name in names.items()}


def _requested_redshifts(raw: str) -> list[float]:
    return [float(item.strip()) for item in raw.split(",") if item.strip()]


def _key_epoch_rows(
    run_id: str,
    tables: Mapping[str, Sequence[Mapping[str, str]]],
    requested: Sequence[float],
) -> list[dict[str, object]]:
    window = tables["window"]
    objects_by_z = _group_by_z(tables["object"])
    transfers_by_z = _group_by_z(tables["transfer"])
    summary_by_z = _group_by_z(tables["summary"])
    rows = []
    for target_z in requested:
        base = _nearest_z_row(window, target_z)
        z_key = _z_key(base)
        object_rows = objects_by_z.get(z_key, [])
        transfer_rows = transfers_by_z.get(z_key, [])
        summary_rows = summary_by_z.get(z_key, [])
        fate = _aggregate_fate(summary_rows)
        row = {
            "run_id": run_id,
            "z": _float(base, "z"),
            "a": _float(base, "a"),
            "k_nl_hmpc": _float(base, "k_nl_hmpc"),
            "lambda_nl_mpch": 2.0 * math.pi / _float(base, "k_nl_hmpc"),
            "m_nl_msunh": _float(base, "m_nl_msunh"),
            "s_aw_global": _float(base, "s_aw_global"),
            "measured_sigma_delta": _float(base, "measured_sigma_delta"),
            "active_power_fraction": _float(base, "active_power_fraction"),
            "tail_power_fraction": _float(base, "tail_power_fraction"),
            "peak_count": int(float(base["peak_count"])),
            "transfer_count": int(float(base["transfer_count"])),
            "object_mass_median_msunh": _median(_values(object_rows, "mass_msunh")),
            "object_radius_median_kpch": _median(_values(object_rows, "radius_kpch")),
            "b_h_median": _median(_values(object_rows, "b_h")),
            "forming_fraction": fate["forming_fraction"],
            "balanced_fraction": fate["balanced_fraction"],
            "stripping_fraction": fate["stripping_fraction"],
            "a_pq_median": _median(_values(transfer_rows, "a_pq")),
            "abs_g_pq_median": _median(abs(value) for value in _values(transfer_rows, "g_pq")),
            "dominant_scale_band": base["dominant_scale_band"],
        }
        rows.append(row)
    return rows


def _stage_rows(run_id: str, tables: Mapping[str, Sequence[Mapping[str, str]]]) -> list[dict[str, object]]:
    window = list(tables["window"])
    summary_by_z = _group_by_z(tables["summary"])
    stages = []
    selectors = (
        (
            "tail_onset",
            "tail_power_fraction",
            lambda row: _float(row, "tail_power_fraction") > 0.0,
            lambda row: _float(row, "tail_power_fraction"),
            "first stored epoch with nonzero residual-tail power",
        ),
        (
            "tail_majority",
            "tail_power_fraction",
            lambda row: _float(row, "tail_power_fraction") >= _float(row, "active_power_fraction"),
            lambda row: _float(row, "tail_power_fraction"),
            "first stored epoch where residual-tail power equals or exceeds active power",
        ),
        (
            "unit_field_variance",
            "measured_sigma_delta",
            lambda row: _float(row, "measured_sigma_delta") >= 1.0,
            lambda row: _float(row, "measured_sigma_delta"),
            "first stored epoch with density variance above unity",
        ),
        (
            "object_balance_dominance",
            "balanced_fraction",
            lambda row: _aggregate_fate(summary_by_z.get(_z_key(row), []))["balanced_fraction"] >= 0.5,
            lambda row: _aggregate_fate(summary_by_z.get(_z_key(row), []))["balanced_fraction"],
            "first stored epoch where balanced objects are at least half of the selected peaks",
        ),
        (
            "future_branch",
            "z",
            lambda row: _float(row, "z") < 0.0,
            lambda row: _float(row, "z"),
            "first stored epoch in the future branch of the generated history",
        ),
    )
    for name, column, predicate, value_of, description in selectors:
        row = _first(window, predicate)
        if row is None:
            continue
        stages.append(
            {
                "run_id": run_id,
                "stage_marker": name,
                "z": _float(row, "z"),
                "a": _float(row, "a"),
                "selected_by": column,
                "value": value_of(row),
                "description": description,
            }
        )
    return stages


def _summary_payload(
    root: Path,
    tables: Mapping[str, Sequence[Mapping[str, str]]],
    key_epochs: Sequence[Mapping[str, object]],
    stage_rows: Sequence[Mapping[str, object]],
) -> dict[str, object]:
    window = tables["window"]
    objects = tables["object"]
    transfers = tables["transfer"]
    cost = tables["cost"][0] if tables["cost"] else {}
    return {
        "run_id": root.name,
        "redshift_count": len(window),
        "redshift_range": [_float(window[-1], "z"), _float(window[0], "z")] if window else [],
        "key_epochs": list(key_epochs),
        "stage_markers": list(stage_rows),
        "window_ranges": {
            key: _range(_values(window, key))
            for key in (
                "k_nl_hmpc",
                "m_nl_msunh",
                "s_aw_global",
                "measured_sigma_delta",
                "active_power_fraction",
                "tail_power_fraction",
            )
        },
        "object_ranges": {
            "mass_msunh": _range(_values(objects, "mass_msunh")),
            "radius_kpch": _range(_values(objects, "radius_kpch")),
            "b_h": _range(_values(objects, "b_h")),
        },
        "transfer_ranges": {
            "a_pq": _range(_values(transfers, "a_pq")),
            "abs_g_pq": _range(abs(value) for value in _values(transfers, "g_pq")),
            "separation_mpch": _range(_values(transfers, "separation_mpch")),
        },
        "counts": {
            "field_window_table": len(tables["window"]),
            "field_mode_table": len(tables["mode"]),
            "field_object_table": len(tables["object"]),
            "field_transfer_table": len(tables["transfer"]),
            "field_relation_summary_table": len(tables["summary"]),
            "field_epoch_relation_table": len(tables["epoch"]),
        },
        "cost": dict(cost),
    }


def _interpretation_markdown(
    run_id: str,
    summary: Mapping[str, object],
    key_epochs: Sequence[Mapping[str, object]],
    stage_rows: Sequence[Mapping[str, object]],
) -> str:
    counts = summary["counts"]
    cost = summary["cost"]
    ranges = summary["window_ranges"]
    lines = [
        f"# Field Result Interpretation for {run_id}",
        "",
        "## Numerical record",
        "",
        (
            f"The run stores {counts['field_window_table']} window rows, "
            f"{counts['field_mode_table']} shell rows, "
            f"{counts['field_object_table']} object rows, "
            f"{counts['field_transfer_table']} transfer rows, and "
            f"{counts['field_epoch_relation_table']} adjacent-epoch rows."
        ),
        (
            f"The cost table records grid size {cost.get('grid_size')}, "
            f"{cost.get('voxel_count')} voxels, {cost.get('relationship_bytes_written')} "
            f"saved relationship bytes, and saved-to-raw ratio {cost.get('saved_to_raw_ratio')}."
        ),
        "",
        "## Main physical trends",
        "",
        (
            "The nonlinear endpoint moves from "
            f"{_fmt_range(ranges['k_nl_hmpc'])} h/Mpc across the stored history. "
            "The associated nonlinear mass range is "
            f"{_fmt_range(ranges['m_nl_msunh'])} Msun/h."
        ),
        (
            "The active-window span covers "
            f"{_fmt_range(ranges['s_aw_global'])}, while the generated field variance covers "
            f"{_fmt_range(ranges['measured_sigma_delta'])}."
        ),
        (
            "The active and residual-tail power fractions cover "
            f"{_fmt_range(ranges['active_power_fraction'])} and "
            f"{_fmt_range(ranges['tail_power_fraction'])}, respectively."
        ),
        "",
        "## Automatically selected stage markers",
        "",
    ]
    for row in stage_rows:
        lines.append(
            f"- {row['stage_marker']}: z={_fmt(row['z'])}, a={_fmt(row['a'])}, "
            f"{row['selected_by']}={_fmt(row['value'])}."
        )
    lines.extend(["", "## Key epochs", ""])
    for row in key_epochs:
        lines.append(
            "- "
            f"z={_fmt(row['z'])}, k_NL={_fmt(row['k_nl_hmpc'])} h/Mpc, "
            f"M_NL={_fmt(row['m_nl_msunh'])} Msun/h, "
            f"sigma={_fmt(row['measured_sigma_delta'])}, "
            f"tail={_fmt(row['tail_power_fraction'])}, "
            f"B_h,50={_fmt(row['b_h_median'])}, "
            f"|G|,50={_fmt(row['abs_g_pq_median'])}."
        )
    return "\n".join(lines) + "\n"


def _group_by_z(rows: Sequence[Mapping[str, str]]) -> dict[str, list[Mapping[str, str]]]:
    grouped: dict[str, list[Mapping[str, str]]] = {}
    for row in rows:
        grouped.setdefault(_z_key(row), []).append(row)
    return grouped


def _z_key(row: Mapping[str, str]) -> str:
    return f"{float(row['z']):.8g}"


def _nearest_z_row(rows: Sequence[Mapping[str, str]], target_z: float) -> Mapping[str, str]:
    return min(rows, key=lambda row: abs(_float(row, "z") - target_z))


def _aggregate_fate(rows: Sequence[Mapping[str, str]]) -> dict[str, float]:
    total = sum(int(float(row["object_count"])) for row in rows)
    if total <= 0:
        return {"forming_fraction": 0.0, "balanced_fraction": 0.0, "stripping_fraction": 0.0}
    return {
        key: sum(int(float(row["object_count"])) * _float(row, key) for row in rows) / total
        for key in ("forming_fraction", "balanced_fraction", "stripping_fraction")
    }


def _values(rows: Iterable[Mapping[str, str]], key: str) -> list[float]:
    values = []
    for row in rows:
        try:
            value = float(row[key])
        except (KeyError, TypeError, ValueError):
            continue
        if math.isfinite(value):
            values.append(value)
    return values


def _median(values: Iterable[float]) -> float:
    sorted_values = sorted(value for value in values if math.isfinite(value))
    if not sorted_values:
        return 0.0
    mid = len(sorted_values) // 2
    if len(sorted_values) % 2:
        return sorted_values[mid]
    return 0.5 * (sorted_values[mid - 1] + sorted_values[mid])


def _range(values: Iterable[float]) -> dict[str, float]:
    selected = [value for value in values if math.isfinite(value)]
    if not selected:
        return {"min": 0.0, "max": 0.0}
    return {"min": min(selected), "max": max(selected)}


def _first(rows: Sequence[Mapping[str, str]], predicate) -> Mapping[str, str] | None:
    for row in rows:
        if predicate(row):
            return row
    return None


def _float(row: Mapping[str, str], key: str) -> float:
    return float(row[key])


def _fmt(value: object) -> str:
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return str(value)
    if numeric == 0.0:
        return "0"
    if abs(numeric) >= 1.0e4 or abs(numeric) < 1.0e-3:
        return f"{numeric:.4e}"
    return f"{numeric:.4g}"


def _fmt_range(row: Mapping[str, float]) -> str:
    return f"{_fmt(row['min'])}-{_fmt(row['max'])}"


if __name__ == "__main__":
    raise SystemExit(main())
