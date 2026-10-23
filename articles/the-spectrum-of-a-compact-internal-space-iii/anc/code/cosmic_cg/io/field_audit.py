"""Audit checks for three-dimensional field-generation runs."""

from __future__ import annotations

import csv
from dataclasses import dataclass
import json
import math
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from cosmic_cg.io.figure_data import FIGURE8_FIELD_COLUMNS
from cosmic_cg.field3d import (
    FIELD_CLOSURE_COLUMNS,
    FIELD_EPOCH_RELATION_COLUMNS,
    FIELD_MODE_COLUMNS,
    FIELD_OBJECT_COLUMNS,
    FIELD_RELATION_SUMMARY_COLUMNS,
    FIELD_TRANSFER_COLUMNS,
    FIELD_WINDOW_COLUMNS,
)


@dataclass(frozen=True)
class FieldAuditReport:
    status: str
    checks: tuple[dict[str, object], ...]


REQUIRED_TABLES: tuple[tuple[str, Sequence[str]], ...] = (
    ("tables/field_window_table.csv", FIELD_WINDOW_COLUMNS),
    ("tables/field_mode_table.csv", FIELD_MODE_COLUMNS),
    ("tables/field_object_table.csv", FIELD_OBJECT_COLUMNS),
    ("tables/field_transfer_table.csv", FIELD_TRANSFER_COLUMNS),
    ("tables/field_relation_summary_table.csv", FIELD_RELATION_SUMMARY_COLUMNS),
    ("tables/field_epoch_relation_table.csv", FIELD_EPOCH_RELATION_COLUMNS),
    ("tables/field_closure_table.csv", FIELD_CLOSURE_COLUMNS),
    ("figure_data/fig08_field_relations.csv", FIGURE8_FIELD_COLUMNS),
)


def audit_field_run(run_dir: str | Path, *, expected_grid_size: int | None = None) -> FieldAuditReport:
    root = Path(run_dir)
    checks = []
    checks.extend(_check_required_files(root))
    manifest = _manifest(root)
    checks.extend(_check_tables(root, manifest))
    checks.extend(_check_manifest(root, expected_grid_size))
    checks.extend(_check_epoch_relation_count(root, manifest))
    checks.extend(_check_numeric_ranges(root))
    checks.append(_check_no_raw_field_arrays(root))
    status = "pass" if all(check["ok"] for check in checks) else "fail"
    return FieldAuditReport(status=status, checks=tuple(checks))


def write_field_audit_report(report: FieldAuditReport, path: str | Path) -> None:
    payload = {"status": report.status, "checks": list(report.checks)}
    Path(path).write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def write_field_audit_markdown(report: FieldAuditReport, path: str | Path) -> None:
    lines = [
        "# 3D Field Run Audit",
        "",
        f"Status: {report.status}",
        "",
        "| Check | Status | Detail |",
        "|---|---|---|",
    ]
    for check in report.checks:
        marker = "pass" if check["ok"] else "fail"
        lines.append(f"| {check['name']} | {marker} | {check.get('detail', '')} |")
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def _check_required_files(root: Path) -> list[dict[str, object]]:
    required = [
        "config_resolved.yaml",
        "params_snapshot.json",
        "generation_manifest.json",
        "provenance.json",
    ]
    checks = []
    for relative in required:
        path = root / relative
        checks.append(
            {
                "name": f"required_file:{relative}",
                "ok": path.is_file(),
                "detail": "present" if path.is_file() else "missing",
            }
        )
    return checks


def _check_tables(root: Path, manifest: Mapping[str, object]) -> list[dict[str, object]]:
    checks = []
    for relative, columns in REQUIRED_TABLES:
        path = root / relative
        if not path.is_file():
            checks.append({"name": f"required_table:{relative}", "ok": False, "detail": "missing"})
            continue
        rows, headers = _read_csv(path)
        missing = [column for column in columns if column not in headers]
        min_rows = _minimum_rows(relative, manifest)
        checks.append(
            {
                "name": f"required_table:{relative}",
                "ok": not missing and len(rows) >= min_rows,
                "detail": f"rows={len(rows)} missing={','.join(missing)}",
            }
        )
    return checks


def _check_manifest(root: Path, expected_grid_size: int | None) -> list[dict[str, object]]:
    payload = _manifest(root)
    if not payload:
        return [{"name": "manifest", "ok": False, "detail": "missing"}]
    checks = [
        {
            "name": "manifest_model",
            "ok": payload.get("model") == "finite_window_3d_field_generation",
            "detail": str(payload.get("model")),
        },
        {
            "name": "manifest_storage_policy",
            "ok": payload.get("storage_policy") == "relationship_tables_only",
            "detail": str(payload.get("storage_policy")),
        },
        {
            "name": "manifest_external_data",
            "ok": payload.get("external_data") == [],
            "detail": json.dumps(payload.get("external_data", None)),
        },
        {
            "name": "manifest_state_decomposition",
            "ok": _has_state_decomposition(payload),
            "detail": str(payload.get("state_decomposition", {})),
        },
        {
            "name": "manifest_closure",
            "ok": isinstance(payload.get("closure"), dict) and bool(payload.get("closure", {}).get("map")),
            "detail": str(payload.get("closure", {})),
        },
    ]
    if expected_grid_size is not None:
        checks.append(
            {
                "name": "manifest_grid_size",
                "ok": int(payload.get("field_config", {}).get("grid_size", -1)) == expected_grid_size,
                "detail": str(payload.get("field_config", {}).get("grid_size", "")),
            }
        )
    return checks


def _check_epoch_relation_count(root: Path, manifest: Mapping[str, object]) -> list[dict[str, object]]:
    redshifts = manifest.get("redshifts", [])
    expected = max(0, len(redshifts) - 1) if isinstance(redshifts, list) else None
    if expected is None:
        return [{"name": "epoch_relation_count", "ok": False, "detail": "manifest redshifts missing"}]
    rows, _headers = _read_csv(root / "tables" / "field_epoch_relation_table.csv")
    return [
        {
            "name": "epoch_relation_count",
            "ok": len(rows) == expected,
            "detail": f"rows={len(rows)} expected={expected}",
        }
    ]


def _check_numeric_ranges(root: Path) -> list[dict[str, object]]:
    checks = []
    window_rows, _headers = _read_csv(root / "tables" / "field_window_table.csv")
    if window_rows:
        checks.append(_all_finite_positive("window_sigma", window_rows, "measured_sigma_delta"))
        checks.append(_all_finite_nonnegative("window_gaussian_sigma", window_rows, "gaussian_sigma_mpch"))
        checks.append(_all_finite_positive("window_k_nl", window_rows, "k_nl_hmpc"))
        checks.append(_fraction_range("window_power_fractions", window_rows, ("active_power_fraction", "tail_power_fraction")))
        checks.append(_all_nonnegative("window_closure_error", window_rows, "closure_error"))
    closure_rows, _headers = _read_csv(root / "tables" / "field_closure_table.csv")
    if closure_rows:
        checks.append(_all_finite_nonnegative("closure_gaussian_sigma", closure_rows, "gaussian_sigma_mpch"))
        checks.append(_all_nonnegative("closure_residual_norm", closure_rows, "residual_norm"))
        checks.append(_all_nonnegative("closure_error", closure_rows, "closure_error"))
    object_rows, _headers = _read_csv(root / "tables" / "field_object_table.csv")
    if object_rows:
        checks.append(_all_finite_positive("object_mass", object_rows, "mass_msunh"))
        checks.append(_all_finite_positive("object_radius", object_rows, "radius_kpch"))
        checks.append(_class_set("object_fate_class", object_rows, "fate_class", {"forming", "balanced", "stripping"}))
    transfer_rows, _headers = _read_csv(root / "tables" / "field_transfer_table.csv")
    if transfer_rows:
        checks.append(_all_nonnegative("transfer_strength", transfer_rows, "a_pq"))
        checks.append(_signed_unit_range("transfer_alignment", transfer_rows, "g_pq"))
    return checks


def _check_no_raw_field_arrays(root: Path) -> dict[str, object]:
    raw_suffixes = {".npy", ".npz", ".h5", ".hdf5"}
    raw = [path for path in root.rglob("*") if path.is_file() and path.suffix.lower() in raw_suffixes]
    return {
        "name": "raw_field_arrays",
        "ok": not raw,
        "detail": "none" if not raw else ",".join(str(path.relative_to(root)) for path in raw[:8]),
    }


def _manifest(root: Path) -> dict[str, object]:
    path = root / "generation_manifest.json"
    if not path.is_file():
        return {}
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return {}
    return payload if isinstance(payload, dict) else {}


def _minimum_rows(relative: str, manifest: Mapping[str, object]) -> int:
    if relative == "tables/field_epoch_relation_table.csv":
        redshifts = manifest.get("redshifts", [])
        if isinstance(redshifts, list) and len(redshifts) <= 1:
            return 0
    return 1


def _all_finite_positive(name: str, rows: Sequence[Mapping[str, str]], key: str) -> dict[str, object]:
    values = [_float(row.get(key, "")) for row in rows]
    ok = bool(values) and all(math.isfinite(value) and value > 0 for value in values)
    return {"name": name, "ok": ok, "detail": _range_detail(values)}


def _all_finite_nonnegative(name: str, rows: Sequence[Mapping[str, str]], key: str) -> dict[str, object]:
    values = [_float(row.get(key, "")) for row in rows]
    ok = bool(values) and all(math.isfinite(value) and value >= 0 for value in values)
    return {"name": name, "ok": ok, "detail": _range_detail(values)}


def _all_nonnegative(name: str, rows: Sequence[Mapping[str, str]], key: str) -> dict[str, object]:
    values = [_float(row.get(key, "")) for row in rows]
    ok = bool(values) and all(math.isfinite(value) and value >= 0 for value in values)
    return {"name": name, "ok": ok, "detail": _range_detail(values)}


def _signed_unit_range(name: str, rows: Sequence[Mapping[str, str]], key: str) -> dict[str, object]:
    values = [_float(row.get(key, "")) for row in rows]
    ok = bool(values) and all(math.isfinite(value) and -1.0000001 <= value <= 1.0000001 for value in values)
    return {"name": name, "ok": ok, "detail": _range_detail(values)}


def _fraction_range(name: str, rows: Sequence[Mapping[str, str]], keys: Iterable[str]) -> dict[str, object]:
    values = [_float(row.get(key, "")) for row in rows for key in keys]
    ok = bool(values) and all(math.isfinite(value) and -1.0e-10 <= value <= 1.0000001 for value in values)
    return {"name": name, "ok": ok, "detail": _range_detail(values)}


def _class_set(name: str, rows: Sequence[Mapping[str, str]], key: str, allowed: set[str]) -> dict[str, object]:
    found = {str(row.get(key, "")) for row in rows}
    return {"name": name, "ok": bool(found) and found <= allowed, "detail": ",".join(sorted(found))}


def _read_csv(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    if not path.is_file():
        return [], []
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        return rows, list(reader.fieldnames or [])


def _has_state_decomposition(payload: Mapping[str, object]) -> bool:
    value = payload.get("state_decomposition")
    if not isinstance(value, dict):
        return False
    return (
        value.get("physical_coarse_graining") == "G_sigma"
        and value.get("activity_selector") == "A_Lambda"
        and value.get("external_data") == []
    )


def _float(value: object) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def _range_detail(values: Sequence[float]) -> str:
    finite = [value for value in values if math.isfinite(value)]
    if not finite:
        return "empty"
    return f"n={len(finite)} min={min(finite):.6g} max={max(finite):.6g}"
