"""Audit checks for internal particle-mesh trajectory runs."""

from __future__ import annotations

import csv
from dataclasses import dataclass
import json
import math
from pathlib import Path
from typing import Mapping, Sequence

from cosmic_cg.pm_nbody import (
    PM_COST_COLUMNS,
    PM_HALO_COLUMNS,
    PM_POWER_COLUMNS,
    PM_SNAPSHOT_COLUMNS,
    PM_TRANSFER_COLUMNS,
    PM_TREE_COLUMNS,
)


@dataclass(frozen=True)
class PMAuditReport:
    status: str
    checks: tuple[dict[str, object], ...]


REQUIRED_TABLES: tuple[tuple[str, Sequence[str]], ...] = (
    ("tables/pm_snapshot_table.csv", PM_SNAPSHOT_COLUMNS),
    ("tables/pm_power_table.csv", PM_POWER_COLUMNS),
    ("tables/pm_halo_table.csv", PM_HALO_COLUMNS),
    ("tables/pm_tree_table.csv", PM_TREE_COLUMNS),
    ("tables/pm_transfer_table.csv", PM_TRANSFER_COLUMNS),
    ("tables/pm_cost_table.csv", PM_COST_COLUMNS),
)


def audit_pm_run(
    run_dir: str | Path,
    *,
    expected_particles_per_dim: int | None = None,
    expected_mesh_size: int | None = None,
) -> PMAuditReport:
    root = Path(run_dir)
    checks = []
    checks.extend(_check_required_files(root))
    manifest = _manifest(root)
    checks.extend(_check_tables(root, manifest))
    checks.extend(_check_manifest(root, expected_particles_per_dim, expected_mesh_size))
    checks.extend(_check_snapshot_count(root, manifest))
    checks.extend(_check_numeric_ranges(root))
    checks.append(_check_no_raw_arrays(root))
    status = "pass" if all(check["ok"] for check in checks) else "fail"
    return PMAuditReport(status=status, checks=tuple(checks))


def write_pm_audit_report(report: PMAuditReport, path: str | Path) -> None:
    payload = {"status": report.status, "checks": list(report.checks)}
    Path(path).write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


def write_pm_audit_markdown(report: PMAuditReport, path: str | Path) -> None:
    lines = [
        "# Particle-Mesh Run Audit",
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
        "intermediate/output_inventory.csv",
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


def _check_manifest(
    root: Path,
    expected_particles_per_dim: int | None,
    expected_mesh_size: int | None,
) -> list[dict[str, object]]:
    payload = _manifest(root)
    if not payload:
        return [{"name": "manifest", "ok": False, "detail": "missing"}]
    pm_config = payload.get("pm_config", {})
    external_data = payload.get("external_data", None)
    state = payload.get("state_decomposition", None)
    closure = payload.get("closure", None)
    dynamics = payload.get("dynamics", None)
    checks = [
        {
            "name": "manifest_model",
            "ok": payload.get("model") == "finite_window_internal_pm_nbody",
            "detail": str(payload.get("model")),
        },
        {
            "name": "manifest_storage_policy",
            "ok": payload.get("storage_policy") == "relationship_tables_only",
            "detail": str(payload.get("storage_policy")),
        },
        {
            "name": "manifest_external_data",
            "ok": external_data == [],
            "detail": str(external_data),
        },
        {
            "name": "manifest_state_decomposition",
            "ok": _has_pm_state_decomposition(state),
            "detail": str(state),
        },
        {
            "name": "manifest_closure",
            "ok": isinstance(closure, Mapping) and bool(closure.get("map")),
            "detail": str(closure),
        },
        {
            "name": "manifest_dynamics",
            "ok": _has_pm_dynamics(dynamics),
            "detail": str(dynamics),
        },
    ]
    if expected_particles_per_dim is not None:
        checks.append(
            {
                "name": "manifest_particles_per_dim",
                "ok": int(pm_config.get("particles_per_dim", -1)) == expected_particles_per_dim,
                "detail": str(pm_config.get("particles_per_dim", "")),
            }
        )
    if expected_mesh_size is not None:
        checks.append(
            {
                "name": "manifest_mesh_size",
                "ok": int(pm_config.get("mesh_size", -1)) == expected_mesh_size,
                "detail": str(pm_config.get("mesh_size", "")),
            }
        )
    return checks


def _check_snapshot_count(root: Path, manifest: Mapping[str, object]) -> list[dict[str, object]]:
    redshifts = manifest.get("redshifts", [])
    expected = len(redshifts) if isinstance(redshifts, list) else None
    rows, _headers = _read_csv(root / "tables" / "pm_snapshot_table.csv")
    return [
        {
            "name": "snapshot_count",
            "ok": expected is not None and len(rows) == expected,
            "detail": f"rows={len(rows)} expected={expected}",
        }
    ]


def _check_numeric_ranges(root: Path) -> list[dict[str, object]]:
    checks = []
    snapshot_rows, _headers = _read_csv(root / "tables" / "pm_snapshot_table.csv")
    if snapshot_rows:
        checks.append(_absolute_bound("snapshot_mass_conservation", snapshot_rows, "mass_conservation_error", 1.0e-5))
        checks.append(_all_finite_positive("snapshot_rms_delta", snapshot_rows, "rms_delta"))
        checks.append(_all_finite_positive("snapshot_particle_count", snapshot_rows, "particle_count"))
    power_rows, _headers = _read_csv(root / "tables" / "pm_power_table.csv")
    if power_rows:
        checks.append(_fraction_range("power_fraction", power_rows, "power_fraction"))
        checks.append(_fraction_range("cumulative_power_fraction", power_rows, "cumulative_power_fraction"))
        checks.append(_all_nonnegative("power_proxy", power_rows, "power_proxy"))
    halo_rows, _headers = _read_csv(root / "tables" / "pm_halo_table.csv")
    if halo_rows:
        checks.append(_all_finite_positive("halo_mass", halo_rows, "mass_msunh"))
        checks.append(_all_finite_positive("halo_radius", halo_rows, "radius_kpch"))
        checks.append(_class_set("halo_fate_class", halo_rows, "fate_class", {"forming", "balanced", "stripping"}))
    tree_rows, _headers = _read_csv(root / "tables" / "pm_tree_table.csv")
    if tree_rows:
        checks.append(_fraction_range("tree_overlap", tree_rows, "overlap_fraction"))
        checks.append(_all_finite_positive("tree_delta_time", tree_rows, "delta_time_gyr"))
    transfer_rows, _headers = _read_csv(root / "tables" / "pm_transfer_table.csv")
    if transfer_rows:
        checks.append(_all_nonnegative("transfer_strength", transfer_rows, "a_pq"))
        checks.append(_signed_unit_range("transfer_alignment", transfer_rows, "g_pq"))
        checks.append(_all_nonnegative("transfer_separation", transfer_rows, "separation_mpch"))
    return checks


def _check_no_raw_arrays(root: Path) -> dict[str, object]:
    raw_suffixes = {".npy", ".npz", ".h5", ".hdf5"}
    raw = [path for path in root.rglob("*") if path.is_file() and path.suffix.lower() in raw_suffixes]
    return {
        "name": "raw_arrays",
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
    redshifts = manifest.get("redshifts", [])
    redshift_count = len(redshifts) if isinstance(redshifts, list) else 0
    if relative in {"tables/pm_tree_table.csv", "tables/pm_transfer_table.csv"} and redshift_count <= 1:
        return 0
    return 1


def _absolute_bound(name: str, rows: Sequence[Mapping[str, str]], key: str, limit: float) -> dict[str, object]:
    values = [_float(row.get(key, "")) for row in rows]
    ok = bool(values) and all(math.isfinite(value) and abs(value) <= limit for value in values)
    return {"name": name, "ok": ok, "detail": _range_detail(values)}


def _all_finite_positive(name: str, rows: Sequence[Mapping[str, str]], key: str) -> dict[str, object]:
    values = [_float(row.get(key, "")) for row in rows]
    ok = bool(values) and all(math.isfinite(value) and value > 0 for value in values)
    return {"name": name, "ok": ok, "detail": _range_detail(values)}


def _all_nonnegative(name: str, rows: Sequence[Mapping[str, str]], key: str) -> dict[str, object]:
    values = [_float(row.get(key, "")) for row in rows]
    ok = bool(values) and all(math.isfinite(value) and value >= 0 for value in values)
    return {"name": name, "ok": ok, "detail": _range_detail(values)}


def _fraction_range(name: str, rows: Sequence[Mapping[str, str]], key: str) -> dict[str, object]:
    values = [_float(row.get(key, "")) for row in rows]
    ok = bool(values) and all(math.isfinite(value) and -1.0e-10 <= value <= 1.0000001 for value in values)
    return {"name": name, "ok": ok, "detail": _range_detail(values)}


def _signed_unit_range(name: str, rows: Sequence[Mapping[str, str]], key: str) -> dict[str, object]:
    values = [_float(row.get(key, "")) for row in rows]
    ok = bool(values) and all(math.isfinite(value) and -1.0000001 <= value <= 1.0000001 for value in values)
    return {"name": name, "ok": ok, "detail": _range_detail(values)}


def _class_set(name: str, rows: Sequence[Mapping[str, str]], key: str, allowed: set[str]) -> dict[str, object]:
    found = {str(row.get(key, "")) for row in rows}
    return {"name": name, "ok": bool(found) and found <= allowed, "detail": ",".join(sorted(found))}


def _has_pm_state_decomposition(value: object) -> bool:
    if not isinstance(value, Mapping):
        return False
    dynamical_state = value.get("dynamical_state")
    required = {"particle_phase_space", "mesh_density", "resolved_force", "halo_catalog", "transfer_graph"}
    return (
        isinstance(dynamical_state, list)
        and required <= set(str(item) for item in dynamical_state)
        and value.get("physical_coarse_graining") == "particle_mesh_assignment"
        and value.get("activity_selector") == "resolved_mesh_modes"
        and value.get("external_data") == []
    )


def _has_pm_dynamics(value: object) -> bool:
    if not isinstance(value, Mapping):
        return False
    return (
        value.get("solver") == "particle_mesh"
        and value.get("force") == "spectral_poisson"
        and value.get("integrator") == "kick_drift_kick"
        and value.get("tracker") == "particle_id_overlap"
        and value.get("time_variable") == "scale_factor"
        and (value.get("force_prefactor") == "Q(a)=force_strength/max(a,1e-4)"
             or (value.get("clock") == "cosmological"
                 and value.get("force_prefactor") == "3*Omega_m0/(2*a^2*E(a))"
                 and value.get("momentum") == "p=a^3*E(a)*dx/da"
                 and value.get("drift_prefactor") == "1/(a^3*E(a))"))
    )


def _read_csv(path: Path) -> tuple[list[dict[str, str]], list[str]]:
    if not path.is_file():
        return [], []
    with path.open("r", newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
        return rows, list(reader.fieldnames or [])


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
