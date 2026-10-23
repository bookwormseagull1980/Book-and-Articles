"""Intermediate tables retained by each analysis run."""

from __future__ import annotations

from dataclasses import asdict
import math
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from cosmic_cg.io.tables import write_rows
from cosmic_cg.knl import KNLResult, dimensionless_power
from cosmic_cg.tail_ratio import virial_velocity


INPUT_VALIDATION_COLUMNS = ("schema", "path", "ok", "row_count", "missing_columns")
POWER_SAMPLE_COLUMNS = (
    "z",
    "a",
    "k_hmpc",
    "p_mpc3_h3",
    "delta2",
    "k_over_knl",
    "sector",
    "tail_weight_s1",
)
ENDPOINT_COLUMNS = (
    "z",
    "c_nl",
    "k_nl_hmpc",
    "crossing_index",
    "method",
    "k_left_hmpc",
    "k_right_hmpc",
    "delta2_left",
    "delta2_right",
)
HALO_SELECTION_COLUMNS = (
    "run_id",
    "source",
    "snapshot",
    "z",
    "halo_id",
    "mass_msun_h",
    "radius_kpc_h",
    "sigma_v_kms",
    "vmax_kms",
    "sigma_source",
    "sigma_used_kms",
    "mass_over_m_nl",
    "r_h_object",
    "selected_for_r_h",
    "quality_flag",
)
MERGER_GROWTH_COLUMNS = (
    "source",
    "halo_id_from",
    "halo_id_to",
    "snapshot_from",
    "snapshot_to",
    "z_from",
    "z_to",
    "mass_from_msun_h",
    "mass_to_msun_h",
    "delta_mass_msun_h",
    "time_from_gyr",
    "time_to_gyr",
    "delta_time_gyr",
    "growth_rate_log_per_gyr",
    "quality_flag",
)
OUTPUT_INVENTORY_COLUMNS = ("relative_path", "kind", "bytes", "rows")


def write_intermediate_tables(
    output_dir: str | Path,
    *,
    input_reports: Iterable[object],
    k_hmpc: Sequence[float],
    p_mpc3_h3: Sequence[float],
    knl: KNLResult,
    halos: Iterable[Mapping[str, object]],
    merger_rows: Iterable[Mapping[str, object]],
    m_nl_msun_h: float,
    v_vir_mnl_kms: float,
    max_mass_factor: float,
) -> None:
    root = Path(output_dir)
    root.mkdir(parents=True, exist_ok=True)
    write_rows(root / "input_validation.csv", _input_validation_rows(input_reports), INPUT_VALIDATION_COLUMNS)
    write_rows(root / "power_spectrum_samples.csv", _power_sample_rows(k_hmpc, p_mpc3_h3, knl), POWER_SAMPLE_COLUMNS)
    write_rows(root / "endpoint_solution.csv", [_endpoint_row(k_hmpc, p_mpc3_h3, knl)], ENDPOINT_COLUMNS)
    write_rows(
        root / "halo_selection.csv",
        _halo_selection_rows(halos, m_nl_msun_h=m_nl_msun_h, v_vir_mnl_kms=v_vir_mnl_kms, max_mass_factor=max_mass_factor),
        HALO_SELECTION_COLUMNS,
    )
    write_rows(root / "merger_growth_links.csv", _merger_growth_rows(merger_rows), MERGER_GROWTH_COLUMNS)


def write_output_inventory(run_dir: str | Path) -> Path:
    root = Path(run_dir)
    rows = []
    for path in sorted(item for item in root.rglob("*") if item.is_file()):
        rel_path = path.relative_to(root).as_posix()
        rows.append(
            {
                "relative_path": rel_path,
                "kind": rel_path.split("/", 1)[0] if "/" in rel_path else "root",
                "bytes": path.stat().st_size,
                "rows": _csv_row_count(path) if path.suffix.lower() == ".csv" else "",
            }
        )
    output = root / "intermediate" / "output_inventory.csv"
    write_rows(output, rows, OUTPUT_INVENTORY_COLUMNS)
    return output


def _input_validation_rows(input_reports: Iterable[object]) -> list[dict[str, object]]:
    rows = []
    for report in input_reports:
        if isinstance(report, Mapping):
            row = dict(report)
            missing = row.get("missing_columns", "")
        else:
            row = asdict(report)
            missing = getattr(report, "missing_columns", "")
        row["missing_columns"] = "|".join(missing) if isinstance(missing, tuple) else missing
        rows.append(row)
    return rows


def _power_sample_rows(k_hmpc: Sequence[float], p_mpc3_h3: Sequence[float], knl: KNLResult) -> list[dict[str, object]]:
    delta2 = dimensionless_power(k_hmpc, p_mpc3_h3)
    rows = []
    spacing = _median_log_spacing(k_hmpc)
    for k_value, p_value, d2 in zip(k_hmpc, p_mpc3_h3, delta2):
        ratio = k_value / knl.k_nl_hmpc
        sector = "active" if k_value <= knl.k_nl_hmpc else "tail"
        distance = max(0.0, math.log(ratio) / spacing)
        rows.append(
            {
                "z": knl.z,
                "a": 1.0 / (1.0 + knl.z),
                "k_hmpc": k_value,
                "p_mpc3_h3": p_value,
                "delta2": d2,
                "k_over_knl": ratio,
                "sector": sector,
                "tail_weight_s1": math.exp(-distance),
            }
        )
    return rows


def _endpoint_row(k_hmpc: Sequence[float], p_mpc3_h3: Sequence[float], knl: KNLResult) -> dict[str, object]:
    delta2 = dimensionless_power(k_hmpc, p_mpc3_h3)
    i = max(0, min(knl.crossing_index, len(k_hmpc) - 2))
    return {
        "z": knl.z,
        "c_nl": knl.c_nl,
        "k_nl_hmpc": knl.k_nl_hmpc,
        "crossing_index": knl.crossing_index,
        "method": knl.method,
        "k_left_hmpc": k_hmpc[i],
        "k_right_hmpc": k_hmpc[i + 1],
        "delta2_left": delta2[i],
        "delta2_right": delta2[i + 1],
    }


def _halo_selection_rows(
    halos: Iterable[Mapping[str, object]],
    *,
    m_nl_msun_h: float,
    v_vir_mnl_kms: float,
    max_mass_factor: float,
) -> list[dict[str, object]]:
    rows = []
    for halo in halos:
        mass = float(halo["mass_msun_h"])
        radius = float(halo["radius_kpc_h"])
        sigma, source = _sigma_value(halo, mass=mass, radius=radius)
        selected = 0.0 < mass <= max_mass_factor * m_nl_msun_h
        rows.append(
            {
                "run_id": halo.get("run_id", ""),
                "source": halo.get("source", ""),
                "snapshot": halo.get("snapshot", ""),
                "z": halo.get("z", ""),
                "halo_id": halo["halo_id"],
                "mass_msun_h": mass,
                "radius_kpc_h": radius,
                "sigma_v_kms": _blank_none(halo.get("sigma_v_kms")),
                "vmax_kms": _blank_none(halo.get("vmax_kms")),
                "sigma_source": source,
                "sigma_used_kms": sigma,
                "mass_over_m_nl": mass / m_nl_msun_h,
                "r_h_object": sigma / v_vir_mnl_kms,
                "selected_for_r_h": selected,
                "quality_flag": halo.get("quality_flag", "ok"),
            }
        )
    return rows


def _merger_growth_rows(merger_rows: Iterable[Mapping[str, object]]) -> list[dict[str, object]]:
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
                "source": row.get("source", ""),
                "halo_id_from": row["halo_id_from"],
                "halo_id_to": row["halo_id_to"],
                "snapshot_from": row.get("snapshot_from", ""),
                "snapshot_to": row.get("snapshot_to", ""),
                "z_from": row.get("z_from", ""),
                "z_to": row.get("z_to", ""),
                "mass_from_msun_h": mass_from,
                "mass_to_msun_h": mass_to,
                "delta_mass_msun_h": mass_to - mass_from,
                "time_from_gyr": time_from,
                "time_to_gyr": time_to,
                "delta_time_gyr": delta_time,
                "growth_rate_log_per_gyr": growth_rate,
                "quality_flag": row.get("quality_flag", "ok"),
            }
        )
    return rows


def _sigma_value(halo: Mapping[str, object], *, mass: float, radius: float) -> tuple[float, str]:
    if halo.get("sigma_v_kms") not in ("", None):
        return float(halo["sigma_v_kms"]), "sigma_v_kms"
    if halo.get("vmax_kms") not in ("", None):
        return float(halo["vmax_kms"]) / math.sqrt(2.0), "vmax_kms/sqrt2"
    return virial_velocity(mass, radius) / math.sqrt(2.0), "virial/sqrt2"


def _blank_none(value: object) -> object:
    return "" if value is None else value


def _median_log_spacing(k_hmpc: Sequence[float]) -> float:
    values = sorted(math.log(k1 / k0) for k0, k1 in zip(k_hmpc, k_hmpc[1:]) if k1 > k0)
    if len(values) == 0:
        return 1.0
    mid = len(values) // 2
    if len(values) % 2:
        return values[mid]
    return 0.5 * (values[mid - 1] + values[mid])


def _csv_row_count(path: Path) -> int:
    try:
        with path.open("r", encoding="utf-8-sig") as fh:
            return max(0, sum(1 for _ in fh) - 1)
    except UnicodeDecodeError:
        return 0
