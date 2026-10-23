"""Derived result tables for the manuscript figures."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
import math
from pathlib import Path
from typing import Iterable, Mapping, Sequence

from cosmic_cg.knl import band_power_fraction


DEFAULT_MASS_BANDS: tuple[tuple[str, float, float, str], ...] = (
    ("dwarf_halo", 1.0e9, 1.0e11, "dwarf halo scale"),
    ("galaxy_halo", 1.0e11, 1.0e13, "galaxy halo scale"),
    ("group_halo", 1.0e13, 1.0e14, "group scale"),
    ("cluster_halo", 1.0e14, 1.0e15, "cluster scale"),
    ("large_environment", 1.0e15, 1.0e17, "large scale environment"),
)

WINDOW_TABLE_COLUMNS = (
    "run_id",
    "source",
    "snapshot",
    "z",
    "a",
    "box_size_mpch",
    "k_box_hmpc",
    "gaussian_sigma_mpch",
    "gaussian_kernel",
    "activity_endpoint_hmpc",
    "activity_selector",
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
    "closure_status",
    "spectrum_source",
    "endpoint_status",
)

BAND_TABLE_COLUMNS = (
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

TAIL_TABLE_COLUMNS = (
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

OBJECT_TABLE_COLUMNS = (
    "run_id",
    "source",
    "snapshot",
    "z",
    "halo_id",
    "x",
    "y",
    "zpos",
    "vx",
    "vy",
    "vz",
    "mass_msunh",
    "r_vir_kpch",
    "sigma_v_kms",
    "spin_x",
    "spin_y",
    "spin_z",
    "b_h",
    "fate_class",
    "parent_id",
    "descendant_id",
    "branch_weight",
    "link_quality",
    "q_h",
)

TRANSFER_TABLE_COLUMNS = (
    "run_id",
    "source",
    "snapshot",
    "z",
    "source_id",
    "target_id",
    "source_channel",
    "a_pq",
    "g_pq",
    "tau_x",
    "tau_y",
    "tau_z",
    "j_x",
    "j_y",
    "j_z",
    "tensor_source",
    "spin_source",
    "vector_status",
)

COST_TABLE_COLUMNS = (
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

EVENT_TABLE_COLUMNS = (
    "run_id",
    "source",
    "snapshot",
    "event_name",
    "z",
    "a",
    "criterion",
    "k_nl_hmpc",
    "m_nl_msunh",
    "s_aw",
    "physical_reading",
)

NORMALIZATION_TABLE_COLUMNS = (
    "run_id",
    "source",
    "snapshot",
    "z",
    "condition_name",
    "condition_value",
    "unit",
    "role",
    "status",
)

OBJECT_SUMMARY_TABLE_COLUMNS = (
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

TRANSFER_SUMMARY_TABLE_COLUMNS = (
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

CLOSURE_TABLE_COLUMNS = (
    "run_id",
    "source",
    "snapshot",
    "z",
    "a",
    "gaussian_sigma_mpch",
    "activity_endpoint_hmpc",
    "closure_map",
    "active_norm",
    "tail_model_norm",
    "residual_norm",
    "closure_error",
    "epsilon0",
    "norm",
    "quality_flag",
)


@dataclass(frozen=True)
class EndpointMassMap:
    c_m: float
    rho_m0_msunh_mpch3: float


def endpoint_mass_map(m_nl_msunh: float, r_nl_kpch: float, k_nl_hmpc: float) -> EndpointMassMap:
    if m_nl_msunh <= 0 or r_nl_kpch <= 0 or k_nl_hmpc <= 0:
        raise ValueError("endpoint inputs must be positive")
    r_nl_mpch = r_nl_kpch / 1000.0
    rho = m_nl_msunh / ((4.0 * math.pi / 3.0) * r_nl_mpch**3)
    return EndpointMassMap(c_m=r_nl_mpch * k_nl_hmpc, rho_m0_msunh_mpch3=rho)


def window_table_row(global_row: Mapping[str, object]) -> dict[str, object]:
    return {
        "run_id": global_row["run_id"],
        "source": global_row["source"],
        "snapshot": global_row["snapshot"],
        "z": global_row["z"],
        "a": global_row["a"],
        "box_size_mpch": global_row["box_size_mpch"],
        "k_box_hmpc": global_row["k_box_hmpc"],
        "gaussian_sigma_mpch": global_row.get("gaussian_sigma_mpch", ""),
        "gaussian_kernel": global_row.get("gaussian_kernel", ""),
        "activity_endpoint_hmpc": global_row.get("activity_endpoint_hmpc", global_row["k_NL_hmpc"]),
        "activity_selector": global_row.get("activity_selector", ""),
        "c_nl": global_row["C_NL"],
        "k_nl_hmpc": global_row["k_NL_hmpc"],
        "m_nl_msunh": global_row["M_NL_msun_h"],
        "r_nl_mpch": float(global_row["R_NL_kpc_h"]) / 1000.0,
        "v_vir_mnl_kms": global_row["v_vir_MNL_kms"],
        "s_abs": global_row["S_abs"],
        "s_aw": global_row["S_AW"],
        "k_cut_hmpc": global_row["k_cut_hmpc"],
        "d_s": global_row["D_S"],
        "v_lambda": global_row.get("v_lambda", ""),
        "closure_active_norm": global_row.get("closure_active_norm", ""),
        "closure_model_norm": global_row.get("closure_model_norm", ""),
        "closure_residual_norm": global_row.get("closure_residual_norm", ""),
        "closure_error": global_row.get("closure_error", ""),
        "closure_status": global_row.get("closure_status", ""),
        "spectrum_source": global_row.get("spectrum_source", global_row["source"]),
        "endpoint_status": global_row.get("quality_flag", "ok"),
    }


def band_table_rows(
    global_row: Mapping[str, object],
    k_hmpc: Sequence[float],
    p_mpc3_h3: Sequence[float],
    mass_bands: Sequence[tuple[str, float, float, str]] = DEFAULT_MASS_BANDS,
) -> list[dict[str, object]]:
    mapper = endpoint_mass_map(
        float(global_row["M_NL_msun_h"]),
        float(global_row["R_NL_kpc_h"]),
        float(global_row["k_NL_hmpc"]),
    )
    rows = []
    for band_id, mass_min, mass_max, label in mass_bands:
        if mass_min <= 0 or mass_max <= mass_min:
            raise ValueError("mass bands must be positive and increasing")
        r_min = mass_to_radius_mpch(mass_min, mapper.rho_m0_msunh_mpch3)
        r_max = mass_to_radius_mpch(mass_max, mapper.rho_m0_msunh_mpch3)
        k_min = mapper.c_m / r_max
        k_max = mapper.c_m / r_min
        rows.append(
            {
                "run_id": global_row["run_id"],
                "source": global_row["source"],
                "snapshot": global_row["snapshot"],
                "z": global_row["z"],
                "a": global_row["a"],
                "band_id": band_id,
                "k_min_hmpc": k_min,
                "k_max_hmpc": k_max,
                "r_min_mpch": r_min,
                "r_max_mpch": r_max,
                "m_min_msunh": mass_min,
                "m_max_msunh": mass_max,
                "representative_structure": label,
                "omega_band": band_power_fraction(k_hmpc, p_mpc3_h3, k_min, k_max),
                "activity_status": activity_status(
                    k_min,
                    k_max,
                    float(global_row["k_box_hmpc"]),
                    float(global_row["k_NL_hmpc"]),
                ),
                "channel_quality": band_quality(k_min, k_max, k_hmpc[0], k_hmpc[-1]),
            }
        )
    return rows


def mass_to_radius_mpch(mass_msunh: float, rho_m0_msunh_mpch3: float) -> float:
    if mass_msunh <= 0 or rho_m0_msunh_mpch3 <= 0:
        raise ValueError("mass and density must be positive")
    return (3.0 * mass_msunh / (4.0 * math.pi * rho_m0_msunh_mpch3)) ** (1.0 / 3.0)


def activity_status(k_min_hmpc: float, k_max_hmpc: float, k_box_hmpc: float, k_nl_hmpc: float) -> str:
    if k_max_hmpc < k_box_hmpc:
        return "infrared"
    if k_min_hmpc <= k_nl_hmpc <= k_max_hmpc:
        return "frontier"
    if k_box_hmpc <= k_min_hmpc and k_max_hmpc < k_nl_hmpc:
        return "active"
    if k_min_hmpc < k_box_hmpc <= k_max_hmpc:
        return "partly_active"
    return "tail"


def band_quality(k_min_hmpc: float, k_max_hmpc: float, data_k_min_hmpc: float, data_k_max_hmpc: float) -> str:
    if k_max_hmpc < data_k_min_hmpc or k_min_hmpc > data_k_max_hmpc:
        return "outside_power_range"
    if k_min_hmpc < data_k_min_hmpc or k_max_hmpc > data_k_max_hmpc:
        return "edge_limited"
    return "ok"


def tail_table_row(global_row: Mapping[str, object]) -> dict[str, object]:
    omega_tail = global_row.get("Omega_tail", "")
    return {
        "run_id": global_row["run_id"],
        "source": global_row["source"],
        "snapshot": global_row["snapshot"],
        "z": global_row["z"],
        "a": global_row["a"],
        "omega_tail_p": omega_tail,
        "r_v": global_row.get("R_v_optional", ""),
        "r_h": global_row.get("R_h_median", ""),
        "e_t_s0": omega_tail,
        "e_t_s1": global_row.get("tail_moment_s1", ""),
        "i_t_s0": global_row.get("tail_injection_s0", ""),
        "i_t_s1": global_row.get("tail_injection_s1", ""),
        "d_t_s0": global_row.get("tail_drain_s0", ""),
        "d_t_s1": global_row.get("tail_drain_s1", ""),
        "omega_floor": global_row.get("Omega_floor", ""),
        "omega_excess": global_row.get("Omega_excess", omega_tail),
        "closure_error": global_row.get("closure_error", ""),
        "tail_channel": "halo_velocity_and_power_tail",
        "q_t": global_row.get("quality_flag", "ok"),
    }


def closure_table_row(global_row: Mapping[str, object]) -> dict[str, object]:
    return {
        "run_id": global_row["run_id"],
        "source": global_row["source"],
        "snapshot": global_row["snapshot"],
        "z": global_row["z"],
        "a": global_row["a"],
        "gaussian_sigma_mpch": global_row.get("gaussian_sigma_mpch", ""),
        "activity_endpoint_hmpc": global_row.get("activity_endpoint_hmpc", global_row["k_NL_hmpc"]),
        "closure_map": global_row.get("closure_map", "spectral_tail_moment_proxy"),
        "active_norm": global_row.get("closure_active_norm", ""),
        "tail_model_norm": global_row.get("closure_model_norm", ""),
        "residual_norm": global_row.get("closure_residual_norm", ""),
        "closure_error": global_row.get("closure_error", ""),
        "epsilon0": global_row.get("closure_epsilon0", ""),
        "norm": global_row.get("closure_norm", "log_power_integral"),
        "quality_flag": global_row.get("closure_status", global_row.get("quality_flag", "ok")),
    }


def object_table_rows(object_rows: Iterable[Mapping[str, object]]) -> list[dict[str, object]]:
    rows = []
    for row in object_rows:
        rows.append(
            {
                "run_id": row["run_id"],
                "source": row["source"],
                "snapshot": row["snapshot"],
                "z": row["z"],
                "halo_id": row["halo_id"],
                "x": row["x"],
                "y": row["y"],
                "zpos": row["z_pos"],
                "vx": row["vx"],
                "vy": row["vy"],
                "vz": row["vz"],
                "mass_msunh": row["mass_msun_h"],
                "r_vir_kpch": row["radius_kpc_h"],
                "sigma_v_kms": row.get("sigma_v_kms", ""),
                "spin_x": row.get("spin_x", ""),
                "spin_y": row.get("spin_y", ""),
                "spin_z": row.get("spin_z", ""),
                "b_h": row["B_h"],
                "fate_class": row["B_class"],
                "parent_id": row.get("parent_id", ""),
                "descendant_id": row.get("descendant_id", ""),
                "branch_weight": row.get("branch_weight", 1.0),
                "link_quality": row.get("link_quality", row.get("quality_flag", "ok")),
                "q_h": row.get("quality_flag", "ok"),
            }
        )
    return rows


def event_table_rows(global_row: Mapping[str, object], band_rows: Iterable[Mapping[str, object]]) -> list[dict[str, object]]:
    rows = [
        {
            "run_id": global_row["run_id"],
            "source": global_row["source"],
            "snapshot": global_row["snapshot"],
            "event_name": "active_frontier",
            "z": global_row["z"],
            "a": global_row["a"],
            "criterion": "Delta2(k_NL)=C_NL",
            "k_nl_hmpc": global_row["k_NL_hmpc"],
            "m_nl_msunh": global_row["M_NL_msun_h"],
            "s_aw": global_row["S_AW"],
            "physical_reading": "current nonlinear frontier",
        }
    ]
    for row in band_rows:
        if row.get("activity_status") != "frontier":
            continue
        rows.append(
            {
                "run_id": global_row["run_id"],
                "source": global_row["source"],
                "snapshot": global_row["snapshot"],
                "event_name": f"{row['band_id']}_frontier",
                "z": global_row["z"],
                "a": global_row["a"],
                "criterion": "band contains k_NL",
                "k_nl_hmpc": global_row["k_NL_hmpc"],
                "m_nl_msunh": global_row["M_NL_msun_h"],
                "s_aw": global_row["S_AW"],
                "physical_reading": row["representative_structure"],
            }
        )
    return rows


def normalization_table_rows(global_row: Mapping[str, object], object_count: int) -> list[dict[str, object]]:
    base = {
        "run_id": global_row["run_id"],
        "source": global_row["source"],
        "snapshot": global_row["snapshot"],
        "z": global_row["z"],
        "status": global_row.get("quality_flag", "ok"),
    }
    rows = []
    for condition, value, unit, role in (
        ("delta2_endpoint", global_row["C_NL"], "dimensionless", "endpoint definition"),
        ("k_nl_hmpc", global_row["k_NL_hmpc"], "h Mpc^-1", "generated endpoint"),
        ("m_nl_msunh", global_row["M_NL_msun_h"], "Msun h^-1", "generated nonlinear mass"),
        ("s_aw", global_row["S_AW"], "dimensionless", "active window entropy"),
        ("omega_tail_p", global_row["Omega_tail"], "dimensionless", "residual tail activity"),
        ("r_h_median", global_row["R_h_median"], "dimensionless", "generated object activity"),
        ("n_objects", object_count, "count", "generated object ensemble"),
    ):
        item = dict(base)
        item.update(
            {
                "condition_name": condition,
                "condition_value": value,
                "unit": unit,
                "role": role,
            }
        )
        rows.append(item)
    return rows


def object_summary_rows(
    object_rows: Iterable[Mapping[str, object]],
    mass_bands: Sequence[tuple[str, float, float, str]] = DEFAULT_MASS_BANDS,
) -> list[dict[str, object]]:
    objects = list(object_rows)
    if not objects:
        return []
    rows = []
    z_values = sorted({_z_bin(row["z"]) for row in objects}, key=_z_sort_key)
    for z_value in z_values:
        objects_at_z = [row for row in objects if _z_bin(row["z"]) == z_value]
        if not objects_at_z:
            continue
        first = objects_at_z[0]
        for band_id, mass_min, mass_max, label in mass_bands:
            selected = [row for row in objects_at_z if mass_min <= float(row["mass_msun_h"]) < mass_max]
            if not selected:
                continue
            classes = Counter(str(row.get("B_class", "")) for row in selected)
            b_values = [_as_float(row.get("B_h")) for row in selected]
            b_values = [value for value in b_values if math.isfinite(value)]
            r_values = [_as_float(row.get("r_h_object")) for row in selected]
            r_values = [value for value in r_values if math.isfinite(value)]
            n = len(selected)
            rows.append(
                {
                    "run_id": first["run_id"],
                    "source": first["source"],
                    "snapshot": first["snapshot"],
                    "z_bin": z_value,
                    "mass_bin": band_id,
                    "n_objects": n,
                    "b_h_median": _quantile_or_blank(b_values, 0.50),
                    "b_h_p16": _quantile_or_blank(b_values, 0.16),
                    "b_h_p84": _quantile_or_blank(b_values, 0.84),
                    "forming_fraction": classes.get("forming", 0) / n,
                    "balanced_fraction": classes.get("balanced", 0) / n,
                    "stripping_fraction": classes.get("stripping", 0) / n,
                    "r_h_median": _quantile_or_blank(r_values, 0.50),
                    "selection": label,
                    "quality_flag": _merge_quality_flags(selected),
                }
            )
    return rows


def transfer_table_rows(
    merger_rows: Iterable[Mapping[str, object]],
    object_by_id: Mapping[str, Mapping[str, object]],
    *,
    run_id: str,
    source: str,
    snapshot: str,
    z: float,
) -> list[dict[str, object]]:
    rows = []
    for link in merger_rows:
        source_id = str(link["halo_id_from"])
        target_id = str(link["halo_id_to"])
        obj = object_by_id.get(source_id, {})
        spin = (obj.get("spin_x", ""), obj.get("spin_y", ""), obj.get("spin_z", ""))
        spin_recorded = all(value not in ("", None) for value in spin)
        rows.append(
            {
                "run_id": run_id,
                "source": source,
                "snapshot": snapshot,
                "z": z,
                "source_id": source_id,
                "target_id": target_id,
                "source_channel": "main_branch_mass_link",
                "a_pq": "",
                "g_pq": "",
                "tau_x": "",
                "tau_y": "",
                "tau_z": "",
                "j_x": spin[0],
                "j_y": spin[1],
                "j_z": spin[2],
                "tensor_source": "",
                "spin_source": "recorded_spin" if spin_recorded else "",
                "vector_status": "spin_recorded" if spin_recorded else "tree_link_recorded",
            }
        )
    return rows


def transfer_summary_rows(transfer_rows: Iterable[Mapping[str, object]]) -> list[dict[str, object]]:
    all_rows = list(transfer_rows)
    if not all_rows:
        return []
    output = []
    groups = sorted(
        {
            (_z_bin(row["z"]), str(row.get("source_channel", "all")))
            for row in all_rows
        },
        key=lambda item: (_z_sort_key(item[0]), item[1]),
    )
    for z_value, source_channel in groups:
        rows = [
            row
            for row in all_rows
            if _z_bin(row["z"]) == z_value
            and str(row.get("source_channel", "all")) == source_channel
        ]
        first = rows[0]
        a_values = [_as_float(row.get("a_pq")) for row in rows]
        g_values = [_as_float(row.get("g_pq")) for row in rows]
        output.append(
            {
                "run_id": first["run_id"],
                "source": first["source"],
                "snapshot": first["snapshot"],
                "z_bin": z_value,
                "source_channel": source_channel,
                "mass_ratio_bin": "all",
                "environment_bin": "all",
                "n_pairs": len(rows),
                "a_pq_median": _quantile_or_blank(a_values, 0.50),
                "a_pq_p16": _quantile_or_blank(a_values, 0.16),
                "a_pq_p84": _quantile_or_blank(a_values, 0.84),
                "g_pq_median": _quantile_or_blank(g_values, 0.50),
                "g_pq_p16": _quantile_or_blank(g_values, 0.16),
                "g_pq_p84": _quantile_or_blank(g_values, 0.84),
                "quality_flag": _merge_quality_flags(rows, key="vector_status"),
            }
        )
    return output


def cost_table_row(
    global_row: Mapping[str, object],
    run_dir: str | Path,
    *,
    k_ref_hmpc: float,
    time_exponent: float,
    access_seconds: float,
    audit_seconds: float,
) -> dict[str, object]:
    lambda_active = float(global_row["k_NL_hmpc"])
    if k_ref_hmpc <= 0:
        raise ValueError("k_ref_hmpc must be positive")
    delta_s = max(0.0, math.log(k_ref_hmpc / lambda_active))
    ratio = min(1.0, lambda_active / k_ref_hmpc)
    root = Path(run_dir)
    table_bytes = directory_size_bytes(root / "tables")
    figure_bytes = directory_size_bytes(root / "figure_data")
    audit_bytes = (
        _file_size(root / "config_resolved.yaml")
        + _file_size(root / "params_snapshot.json")
        + _file_size(root / "generation_manifest.json")
        + _file_size(root / "provenance.json")
        + directory_size_bytes(root / "quality_checks")
    )
    total = table_bytes + figure_bytes + audit_bytes
    ref_bytes = table_bytes / max(ratio**3, 1.0e-300)
    return {
        "run_id": global_row["run_id"],
        "source": global_row["source"],
        "snapshot": global_row["snapshot"],
        "z": global_row["z"],
        "a": global_row["a"],
        "lambda_active_hmpc": lambda_active,
        "k_ref_hmpc": k_ref_hmpc,
        "delta_s": delta_s,
        "rho_n_pred": ratio**3,
        "rho_c_pred": ratio ** (3.0 + time_exponent),
        "active_bytes": table_bytes,
        "residual_bytes": figure_bytes,
        "audit_bytes": audit_bytes,
        "total_bytes": total,
        "access_seconds": access_seconds,
        "audit_seconds": audit_seconds,
        "rho_total_measured": total / ref_bytes if ref_bytes > 0 else "",
    }


def _interval(lo: object, hi: object) -> str:
    if lo in ("", None) or hi in ("", None):
        return ""
    return f"{lo},{hi}"


def _as_float(value: object) -> float:
    if value in ("", None):
        return math.nan
    try:
        return float(value)
    except (TypeError, ValueError):
        return math.nan


def _quantile_or_blank(values: Sequence[float], q: float) -> float | str:
    finite = sorted(value for value in values if math.isfinite(value))
    if not finite:
        return ""
    if q <= 0:
        return finite[0]
    if q >= 1:
        return finite[-1]
    pos = q * (len(finite) - 1)
    lo = math.floor(pos)
    hi = math.ceil(pos)
    if lo == hi:
        return finite[lo]
    frac = pos - lo
    return finite[lo] * (1.0 - frac) + finite[hi] * frac


def _merge_quality_flags(rows: Iterable[Mapping[str, object]], key: str = "quality_flag") -> str:
    flags = sorted({str(row.get(key, "")) for row in rows if row.get(key, "") not in ("", None)})
    return "|".join(flags) if flags else "ok"


def _z_bin(value: object) -> str:
    try:
        return f"{float(value):.4g}"
    except (TypeError, ValueError):
        return str(value)


def _file_size(path: Path) -> int:
    return path.stat().st_size if path.is_file() else 0


def directory_size_bytes(path: str | Path) -> int:
    root = Path(path)
    if not root.exists():
        return 0
    return sum(item.stat().st_size for item in root.rglob("*") if item.is_file())


def _z_sort_key(value: str) -> float:
    try:
        return -float(value)
    except ValueError:
        return math.inf
