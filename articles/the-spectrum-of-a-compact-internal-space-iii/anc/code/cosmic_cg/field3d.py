"""Three-dimensional finite-window field generation."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import ctypes
import math
from typing import Callable, Mapping, Sequence

from cosmic_cg.generation import (
    InternalSpectrumConfig,
    cosmic_age_gyr,
    generate_epoch,
    growth_factor,
    omega_m0,
    spectrum_normalization_amplitude,
)
from cosmic_cg.knl import dimensionless_power
from cosmic_cg.object_level import alignment_cosine, dynamical_time_gyr
from cosmic_cg.result_tables import DEFAULT_MASS_BANDS, mass_to_radius_mpch

try:
    import numpy as np
except ModuleNotFoundError:  # pragma: no cover
    np = None  # type: ignore[assignment]


RHO_CRIT_MSUNH_MPCH3 = 2.77536627e11


@dataclass(frozen=True)
class Field3DConfig:
    grid_size: int = 512
    field_box_size_mpch: float = 512.0
    global_box_size_mpch: float = 14000.0
    k_cut_hmpc: float = 1.0
    gaussian_sigma_mpch: float = 1.0
    closure_floor: float = 1.0e-12
    seed: int = 20260903
    mode_bins: int = 48
    tail_width: float = 0.75
    peak_threshold_sigma: float = 2.5
    peak_exclusion_cells: int = 5
    peak_patch_radius_cells: int = 3
    peak_slab_cells: int = 16
    max_peak_candidates: int = 8192
    max_peaks_per_snapshot: int = 384
    neighbor_links: int = 3
    object_step_gyr: float = 0.25


FIELD_WINDOW_COLUMNS = (
    "run_id",
    "source",
    "snapshot",
    "z",
    "a",
    "grid_size",
    "voxel_count",
    "field_box_size_mpch",
    "global_box_size_mpch",
    "dx_mpch",
    "gaussian_sigma_mpch",
    "gaussian_kernel",
    "activity_endpoint_hmpc",
    "activity_selector",
    "k_box_field_hmpc",
    "k_box_global_hmpc",
    "k_nyquist_hmpc",
    "k_cut_hmpc",
    "k_resolved_hmpc",
    "c_nl",
    "k_nl_hmpc",
    "m_nl_msunh",
    "s_aw_global",
    "target_sigma_delta",
    "measured_sigma_delta",
    "active_mode_count",
    "tail_mode_count",
    "active_power_fraction",
    "tail_power_fraction",
    "closure_active_norm",
    "closure_model_norm",
    "closure_residual_norm",
    "closure_error",
    "closure_status",
    "field_mean_delta",
    "field_min_delta",
    "field_max_delta",
    "field_skewness",
    "field_kurtosis",
    "peak_count",
    "transfer_count",
    "dominant_scale_band",
    "dominant_relation",
    "quality_flag",
)

FIELD_MODE_COLUMNS = (
    "run_id",
    "source",
    "snapshot",
    "z",
    "a",
    "mode_bin",
    "gaussian_sigma_mpch",
    "activity_endpoint_hmpc",
    "k_min_hmpc",
    "k_max_hmpc",
    "k_mid_hmpc",
    "mode_count",
    "active_mode_count",
    "tail_mode_count",
    "power_proxy",
    "power_fraction",
    "active_power_fraction",
    "tail_power_fraction",
    "cumulative_power_fraction",
    "delta2_mean",
    "coarse_delta2_mean",
    "gaussian_weight_mean",
    "activity_weight_mean",
    "window_weight_mean",
    "scale_band",
    "activity_sector",
)

FIELD_OBJECT_COLUMNS = (
    "run_id",
    "source",
    "snapshot",
    "z",
    "object_id",
    "rank",
    "ix",
    "iy",
    "iz",
    "x_mpch",
    "y_mpch",
    "z_mpch",
    "peak_delta",
    "mean_delta_patch",
    "delta_contrast_to_snapshot",
    "support_voxels",
    "mass_msunh",
    "radius_kpch",
    "vx_proxy",
    "vy_proxy",
    "vz_proxy",
    "spin_x",
    "spin_y",
    "spin_z",
    "b_h",
    "growth_link_mass_ratio",
    "fate_class",
    "scale_band",
    "quality_flag",
)

FIELD_TRANSFER_COLUMNS = (
    "run_id",
    "source",
    "snapshot",
    "z",
    "source_id",
    "target_id",
    "rank_source",
    "rank_target",
    "separation_mpch",
    "mass_ratio",
    "a_pq",
    "g_pq",
    "tau_x",
    "tau_y",
    "tau_z",
    "j_x",
    "j_y",
    "j_z",
    "source_channel",
    "quality_flag",
)

FIELD_RELATION_SUMMARY_COLUMNS = (
    "run_id",
    "source",
    "snapshot",
    "z",
    "a",
    "scale_band",
    "object_count",
    "mass_total_msunh",
    "mass_median_msunh",
    "mass_p16_msunh",
    "mass_p84_msunh",
    "b_h_median",
    "b_h_p16",
    "b_h_p84",
    "forming_fraction",
    "balanced_fraction",
    "stripping_fraction",
    "transfer_count",
    "a_pq_median",
    "a_pq_p16",
    "a_pq_p84",
    "abs_g_pq_median",
    "abs_g_pq_p16",
    "abs_g_pq_p84",
    "relation_status",
)

FIELD_EPOCH_RELATION_COLUMNS = (
    "run_id",
    "source",
    "snapshot_from",
    "snapshot_to",
    "z_from",
    "z_to",
    "a_from",
    "a_to",
    "delta_ln_a",
    "delta_ln_k_nl",
    "endpoint_velocity",
    "delta_sigma",
    "delta_peak_count",
    "delta_transfer_count",
    "dominant_scale_from",
    "dominant_scale_to",
    "relation_transition",
    "quality_flag",
)

FIELD_CLOSURE_COLUMNS = (
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

FIELD_COST_COLUMNS = (
    "run_id",
    "source",
    "grid_size",
    "redshift_count",
    "voxel_count",
    "half_spectrum_modes",
    "field_box_size_mpch",
    "estimated_peak_memory_gb",
    "raw_density_bytes_per_snapshot",
    "raw_density_bytes_all_snapshots",
    "relationship_rows",
    "relationship_bytes_written",
    "saved_to_raw_ratio",
    "storage_policy",
)


def require_numpy():
    if np is None:
        raise RuntimeError("NumPy is required for the 3D field generator")
    return np


def config_to_dict(config: Field3DConfig) -> dict[str, object]:
    return asdict(config)


def generate_field_history(
    *,
    run_id: str,
    redshifts: Sequence[float],
    params,
    spectrum_config: InternalSpectrumConfig,
    field_config: Field3DConfig,
    progress: Callable[[str], None] | None = None,
) -> dict[str, list[dict[str, object]]]:
    _validate_config(field_config)
    npx = require_numpy()
    setup = _prepare_spectral_grid(field_config)
    base_modes = _base_modes(field_config, npx)
    outputs: dict[str, list[dict[str, object]]] = {
        "window": [],
        "mode": [],
        "object": [],
        "transfer": [],
        "relation_summary": [],
        "epoch_relation": [],
        "closure": [],
    }
    previous_objects: dict[str, list[dict[str, object]]] = {}
    previous_time_gyr: float | None = None
    for z in redshifts:
        if progress is not None:
            progress(f"field snapshot z={z:g}")
        epoch = generate_epoch(z=z, params=params, config=spectrum_config)
        time_gyr = cosmic_age_gyr(epoch.z, params=params)
        realization = _realize_snapshot(
            epoch=epoch,
            base_modes=base_modes,
            setup=setup,
            params=params,
            spectrum_config=spectrum_config,
            field_config=field_config,
        )
        objects = _extract_objects(
            run_id=run_id,
            epoch=epoch,
            field=realization["field"],
            params=params,
            field_config=field_config,
        )
        delta_time_gyr = (
            max(time_gyr - previous_time_gyr, 1.0e-12)
            if previous_time_gyr is not None
            else field_config.object_step_gyr
        )
        _attach_growth_from_previous(objects, previous_objects, delta_time_gyr)
        previous_objects = _objects_by_band(objects)
        previous_time_gyr = time_gyr
        transfers = _transfer_rows(
            run_id=run_id,
            epoch=epoch,
            objects=objects,
            field_config=field_config,
        )
        relation_rows = _relation_summary_rows(run_id, epoch, objects, transfers)
        mode_rows = _mode_rows(
            run_id=run_id,
            epoch=epoch,
            setup=setup,
            realization=realization,
            field_config=field_config,
            objects=objects,
            transfers=transfers,
        )
        window_row = _window_row(
            run_id=run_id,
            epoch=epoch,
            setup=setup,
            realization=realization,
            field_config=field_config,
            objects=objects,
            transfers=transfers,
            mode_rows=mode_rows,
        )
        outputs["window"].append(window_row)
        outputs["mode"].extend(mode_rows)
        outputs["object"].extend(objects)
        outputs["transfer"].extend(transfers)
        outputs["relation_summary"].extend(relation_rows)
        outputs["closure"].append(_closure_row(run_id, epoch, realization, field_config))
        if len(outputs["window"]) >= 2:
            outputs["epoch_relation"].append(_epoch_relation_row(run_id, outputs["window"][-2], outputs["window"][-1]))
        del realization
    return outputs


def estimate_field_resources(
    *,
    redshift_count: int,
    field_config: Field3DConfig,
    relationship_rows: int = 0,
    relationship_bytes_written: int = 0,
) -> dict[str, object]:
    _validate_config(field_config)
    n = field_config.grid_size
    voxel_count = n**3
    half_modes = n * n * (n // 2 + 1)
    real_bytes = voxel_count * 4
    complex_bytes = half_modes * 8
    half_real_bytes = half_modes * 4
    peak_memory = (
        2 * real_bytes
        + 2 * complex_bytes
        + 4 * half_real_bytes
        + voxel_count
        + half_modes * 2
    )
    raw_all = real_bytes * redshift_count
    ratio = relationship_bytes_written / raw_all if raw_all > 0 else 0.0
    return {
        "grid_size": n,
        "redshift_count": redshift_count,
        "voxel_count": voxel_count,
        "half_spectrum_modes": half_modes,
        "estimated_peak_memory_gb": peak_memory / 1024**3,
        "raw_density_bytes_per_snapshot": real_bytes,
        "raw_density_bytes_all_snapshots": raw_all,
        "relationship_rows": relationship_rows,
        "relationship_bytes_written": relationship_bytes_written,
        "saved_to_raw_ratio": ratio,
    }


def available_physical_memory_bytes() -> int | None:
    class MemoryStatusEx(ctypes.Structure):
        _fields_ = [
            ("dwLength", ctypes.c_ulong),
            ("dwMemoryLoad", ctypes.c_ulong),
            ("ullTotalPhys", ctypes.c_ulonglong),
            ("ullAvailPhys", ctypes.c_ulonglong),
            ("ullTotalPageFile", ctypes.c_ulonglong),
            ("ullAvailPageFile", ctypes.c_ulonglong),
            ("ullTotalVirtual", ctypes.c_ulonglong),
            ("ullAvailVirtual", ctypes.c_ulonglong),
            ("ullAvailExtendedVirtual", ctypes.c_ulonglong),
        ]

    status = MemoryStatusEx()
    status.dwLength = ctypes.sizeof(MemoryStatusEx)
    try:
        ok = ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(status))
    except AttributeError:
        return None
    if not ok:
        return None
    return int(status.ullAvailPhys)


def field_cost_row(
    *,
    run_id: str,
    field_config: Field3DConfig,
    redshift_count: int,
    relationship_rows: int,
    relationship_bytes_written: int,
) -> dict[str, object]:
    estimate = estimate_field_resources(
        redshift_count=redshift_count,
        field_config=field_config,
        relationship_rows=relationship_rows,
        relationship_bytes_written=relationship_bytes_written,
    )
    return {
        "run_id": run_id,
        "source": "internal_3d_field",
        "grid_size": estimate["grid_size"],
        "redshift_count": estimate["redshift_count"],
        "voxel_count": estimate["voxel_count"],
        "half_spectrum_modes": estimate["half_spectrum_modes"],
        "field_box_size_mpch": field_config.field_box_size_mpch,
        "estimated_peak_memory_gb": estimate["estimated_peak_memory_gb"],
        "raw_density_bytes_per_snapshot": estimate["raw_density_bytes_per_snapshot"],
        "raw_density_bytes_all_snapshots": estimate["raw_density_bytes_all_snapshots"],
        "relationship_rows": estimate["relationship_rows"],
        "relationship_bytes_written": estimate["relationship_bytes_written"],
        "saved_to_raw_ratio": estimate["saved_to_raw_ratio"],
        "storage_policy": "relationship_tables_only",
    }


def _validate_config(config: Field3DConfig) -> None:
    if config.grid_size < 8 or config.grid_size % 2:
        raise ValueError("grid_size must be an even integer greater than or equal to 8")
    if config.field_box_size_mpch <= 0 or config.global_box_size_mpch <= 0:
        raise ValueError("box sizes must be positive")
    if config.k_cut_hmpc <= 0 or config.mode_bins <= 0:
        raise ValueError("k_cut_hmpc and mode_bins must be positive")
    if config.gaussian_sigma_mpch < 0 or config.closure_floor <= 0:
        raise ValueError("gaussian_sigma_mpch must be nonnegative and closure_floor must be positive")
    if config.tail_width <= 0:
        raise ValueError("tail_width must be positive")
    if config.peak_threshold_sigma <= 0:
        raise ValueError("peak_threshold_sigma must be positive")
    if config.peak_exclusion_cells <= 0 or config.peak_patch_radius_cells <= 0 or config.peak_slab_cells <= 0:
        raise ValueError("peak radii must be positive")
    if config.max_peaks_per_snapshot <= 0 or config.neighbor_links <= 0:
        raise ValueError("object and link counts must be positive")


def _prepare_spectral_grid(config: Field3DConfig) -> dict[str, object]:
    npx = require_numpy()
    n = config.grid_size
    dx = config.field_box_size_mpch / n
    kx = (2.0 * math.pi * npx.fft.fftfreq(n, d=dx)).astype(npx.float32)
    ky = (2.0 * math.pi * npx.fft.fftfreq(n, d=dx)).astype(npx.float32)
    kz = (2.0 * math.pi * npx.fft.rfftfreq(n, d=dx)).astype(npx.float32)
    k2 = (
        kx[:, None, None] * kx[:, None, None]
        + ky[None, :, None] * ky[None, :, None]
        + kz[None, None, :] * kz[None, None, :]
    ).astype(npx.float32)
    k_mag = npx.sqrt(k2, dtype=npx.float32)
    del k2
    k_mag[0, 0, 0] = 0.0
    k_nyquist = math.pi / dx
    k_resolved = min(config.k_cut_hmpc, k_nyquist)
    bin_edges = npx.geomspace(
        max(2.0 * math.pi / config.field_box_size_mpch, 1.0e-8),
        k_resolved,
        config.mode_bins + 1,
        dtype=npx.float32,
    )
    bin_index = _mode_bin_index(k_mag, bin_edges, k_resolved)
    bin_counts = npx.bincount(bin_index.ravel(), minlength=config.mode_bins + 1).astype(npx.int64)
    return {
        "dx": dx,
        "k_mag": k_mag,
        "bin_edges": bin_edges,
        "bin_index": bin_index,
        "bin_counts": bin_counts,
        "k_box_field": 2.0 * math.pi / config.field_box_size_mpch,
        "k_box_global": 2.0 * math.pi / config.global_box_size_mpch,
        "k_nyquist": k_nyquist,
        "k_resolved": k_resolved,
    }


def _mode_bin_index(k_mag, bin_edges, k_resolved: float):
    npx = require_numpy()
    valid = (k_mag > 0.0) & (k_mag <= k_resolved)
    scaled = npx.zeros(k_mag.shape, dtype=npx.uint16)
    if npx.any(valid):
        log_min = float(npx.log(bin_edges[0]))
        log_step = float(npx.log(bin_edges[-1] / bin_edges[0])) / (len(bin_edges) - 1)
        raw = npx.floor((npx.log(npx.maximum(k_mag, bin_edges[0])) - log_min) / log_step).astype(npx.int32) + 1
        raw = npx.clip(raw, 1, len(bin_edges) - 1)
        scaled[valid] = raw[valid].astype(npx.uint16)
    return scaled


def _base_modes(config: Field3DConfig, npx):
    rng = npx.random.default_rng(config.seed)
    n = config.grid_size
    white = rng.normal(0.0, 1.0, size=(n, n, n)).astype(npx.float32)
    white -= npx.mean(white, dtype=npx.float64).astype(npx.float32)
    modes = npx.fft.rfftn(white, axes=(0, 1, 2))
    modes[0, 0, 0] = 0.0
    del white
    return modes


def _realize_snapshot(
    *,
    epoch,
    base_modes,
    setup: Mapping[str, object],
    params,
    spectrum_config: InternalSpectrumConfig,
    field_config: Field3DConfig,
) -> dict[str, object]:
    npx = require_numpy()
    k_mag = setup["k_mag"]
    factor = _spectrum_factor(k_mag, z=epoch.z, params=params, spectrum_config=spectrum_config)
    gaussian_weight = _gaussian_weight(k_mag, field_config.gaussian_sigma_mpch)
    activity_weight = _window_weight(
        k_mag,
        k_nl=epoch.knl.k_nl_hmpc,
        k_cut=float(setup["k_resolved"]),
        tail_width=field_config.tail_width,
    )
    factor *= gaussian_weight
    factor *= activity_weight
    modes = base_modes * factor
    del factor, gaussian_weight, activity_weight
    field = npx.fft.irfftn(
        modes,
        s=(field_config.grid_size,) * 3,
        axes=(0, 1, 2),
    ).astype(npx.float32, copy=False)
    field -= npx.mean(field, dtype=npx.float64).astype(npx.float32)
    target_sigma = _target_sigma(epoch=epoch, field_config=field_config, k_resolved=float(setup["k_resolved"]))
    measured = float(npx.std(field, dtype=npx.float64))
    scale = target_sigma / measured if measured > 0 else 0.0
    if scale != 0.0:
        field *= npx.float32(scale)
        modes *= npx.float32(scale)
    sector = _mode_power_counts(
        modes,
        setup["bin_index"],
        k_mag,
        field_config.mode_bins,
        k_nl=epoch.knl.k_nl_hmpc,
        k_cut=float(setup["k_resolved"]),
    )
    power_by_bin = sector["power"]
    active_counts = sector["active_counts"]
    tail_counts = sector["tail_counts"]
    active_power_by_bin = sector["active_power"]
    tail_power_by_bin = sector["tail_power"]
    total_power = float(npx.sum(power_by_bin[1:], dtype=npx.float64))
    active_power = float(npx.sum(active_power_by_bin[1:], dtype=npx.float64))
    tail_power = float(npx.sum(tail_power_by_bin[1:], dtype=npx.float64))
    closure = _closure_audit(epoch=epoch, field_config=field_config, k_resolved=float(setup["k_resolved"]))
    return {
        "field": field,
        "modes": modes,
        "power_by_bin": power_by_bin,
        "active_power_by_bin": active_power_by_bin,
        "tail_power_by_bin": tail_power_by_bin,
        "active_counts": active_counts,
        "tail_counts": tail_counts,
        "target_sigma": target_sigma,
        "measured_sigma": float(npx.std(field, dtype=npx.float64)),
        "field_mean": float(npx.mean(field, dtype=npx.float64)),
        "field_min": float(npx.min(field)),
        "field_max": float(npx.max(field)),
        "field_skewness": _field_moment_ratio(field, 3),
        "field_kurtosis": _field_moment_ratio(field, 4),
        "active_power_fraction": active_power / total_power if total_power > 0 else 0.0,
        "tail_power_fraction": tail_power / total_power if total_power > 0 else 0.0,
        "closure": closure,
    }


def _spectrum_factor(k_mag, *, z: float, params, spectrum_config: InternalSpectrumConfig):
    npx = require_numpy()
    k_safe = npx.maximum(k_mag, npx.float32(1.0e-30))
    amplitude = spectrum_normalization_amplitude(spectrum_config)
    growth = growth_factor(z, params=params)
    shape = npx.power(k_safe / spectrum_config.k_pivot_hmpc, spectrum_config.n_s, dtype=npx.float32)
    if spectrum_config.spectrum_shape == "finite_window_envelope":
        envelope = npx.power(
            1.0 + npx.power(k_safe / spectrum_config.uv_envelope_k_hmpc, spectrum_config.uv_envelope_alpha),
            -spectrum_config.uv_envelope_beta,
        ).astype(npx.float32)
        shape *= envelope
    elif spectrum_config.spectrum_shape != "power_law":
        raise ValueError(f"unknown spectrum shape: {spectrum_config.spectrum_shape}")
    factor = npx.sqrt(amplitude * growth**2 * shape, dtype=npx.float32)
    factor[k_mag <= 0.0] = 0.0
    return factor


def _window_weight(k_mag, *, k_nl: float, k_cut: float, tail_width: float):
    npx = require_numpy()
    weight = npx.zeros(k_mag.shape, dtype=npx.float32)
    active = (k_mag > 0.0) & (k_mag <= min(k_nl, k_cut))
    weight[active] = 1.0
    tail = (k_mag > k_nl) & (k_mag <= k_cut)
    if npx.any(tail):
        ratio = npx.maximum(k_mag[tail] / max(k_nl, 1.0e-30), 1.0)
        weight[tail] = npx.exp(-0.5 * npx.square(npx.log(ratio) / tail_width)).astype(npx.float32)
    return weight


def _gaussian_weight(k_mag, gaussian_sigma_mpch: float):
    npx = require_numpy()
    if gaussian_sigma_mpch <= 0:
        return npx.ones(k_mag.shape, dtype=npx.float32)
    return npx.exp(-0.5 * npx.square(k_mag * npx.float32(gaussian_sigma_mpch))).astype(npx.float32)


def _target_sigma(
    *,
    epoch,
    field_config: Field3DConfig,
    k_resolved: float,
) -> float:
    npx = require_numpy()
    k = npx.asarray(epoch.k_hmpc, dtype=npx.float64)
    delta2 = npx.asarray(dimensionless_power(epoch.k_hmpc, epoch.p_mpc3_h3), dtype=npx.float64)
    gaussian = _gaussian_weight(k.astype(npx.float32), field_config.gaussian_sigma_mpch).astype(npx.float64)
    activity = _window_weight(k.astype(npx.float32), k_nl=epoch.knl.k_nl_hmpc, k_cut=k_resolved, tail_width=field_config.tail_width).astype(npx.float64)
    valid = (k > 0.0) & (k <= k_resolved)
    if npx.count_nonzero(valid) < 2:
        return 0.0
    x = npx.log(k[valid])
    y = delta2[valid] * gaussian[valid] ** 2 * activity[valid] ** 2
    variance = _trapz(y, x)
    return math.sqrt(max(float(variance), 0.0))


def _closure_audit(*, epoch, field_config: Field3DConfig, k_resolved: float) -> dict[str, float | str]:
    npx = require_numpy()
    k = npx.asarray(epoch.k_hmpc, dtype=npx.float64)
    delta2 = npx.asarray(dimensionless_power(epoch.k_hmpc, epoch.p_mpc3_h3), dtype=npx.float64)
    gaussian = _gaussian_weight(k.astype(npx.float32), field_config.gaussian_sigma_mpch).astype(npx.float64)
    activity = _window_weight(
        k.astype(npx.float32),
        k_nl=epoch.knl.k_nl_hmpc,
        k_cut=k_resolved,
        tail_width=field_config.tail_width,
    ).astype(npx.float64)
    coarse = delta2 * gaussian**2
    active_mask = (k > 0.0) & (k <= min(epoch.knl.k_nl_hmpc, k_resolved))
    tail_mask = (k > epoch.knl.k_nl_hmpc) & (k <= k_resolved)
    full_mask = k > 0.0
    active_norm = _masked_log_integral(coarse, k, active_mask)
    tail_model_norm = _masked_log_integral(coarse * activity**2, k, tail_mask)
    total_norm = _masked_log_integral(coarse, k, full_mask)
    residual_norm = max(0.0, total_norm - active_norm - tail_model_norm)
    denom = active_norm + tail_model_norm + field_config.closure_floor
    closure_error = residual_norm / denom if denom > 0 else 0.0
    return {
        "active_norm": active_norm,
        "tail_model_norm": tail_model_norm,
        "total_norm": total_norm,
        "residual_norm": residual_norm,
        "closure_error": closure_error,
        "quality_flag": "closure_audited" if closure_error <= 0.1 else "residual_tracked",
    }
def _mode_power_counts(modes, bin_index, k_mag, mode_bins: int, *, k_nl: float, k_cut: float) -> dict[str, object]:
    npx = require_numpy()
    mode_power = (modes.real * modes.real + modes.imag * modes.imag).astype(npx.float32)
    active = (k_mag > 0.0) & (k_mag <= min(k_nl, k_cut))
    tail = (k_mag > k_nl) & (k_mag <= k_cut)
    power = npx.bincount(bin_index.ravel(), weights=mode_power.ravel(), minlength=mode_bins + 1)
    active_power = npx.bincount(
        bin_index.ravel(),
        weights=(mode_power * active).ravel(),
        minlength=mode_bins + 1,
    )
    tail_power = npx.bincount(
        bin_index.ravel(),
        weights=(mode_power * tail).ravel(),
        minlength=mode_bins + 1,
    )
    active_counts = npx.bincount(bin_index.ravel(), weights=active.ravel(), minlength=mode_bins + 1).astype(npx.int64)
    tail_counts = npx.bincount(bin_index.ravel(), weights=tail.ravel(), minlength=mode_bins + 1).astype(npx.int64)
    del mode_power, active, tail
    return {
        "power": power,
        "active_power": active_power,
        "tail_power": tail_power,
        "active_counts": active_counts,
        "tail_counts": tail_counts,
    }


def _extract_objects(
    *,
    run_id: str,
    epoch,
    field,
    params,
    field_config: Field3DConfig,
) -> list[dict[str, object]]:
    npx = require_numpy()
    sigma = float(npx.std(field, dtype=npx.float64))
    if sigma <= 0:
        return []
    threshold = field_config.peak_threshold_sigma * sigma
    candidates = _peak_candidates(field, threshold=threshold, field_config=field_config)
    rows = []
    dx = field_config.field_box_size_mpch / field_config.grid_size
    rho_m = RHO_CRIT_MSUNH_MPCH3 * omega_m0(params)
    for rank, flat_index in enumerate(candidates, start=1):
        ix, iy, iz = _unravel_index(int(flat_index), field_config.grid_size)
        patch = _patch(field, ix, iy, iz, field_config.peak_patch_radius_cells)
        support = patch > threshold
        support_voxels = int(npx.count_nonzero(support))
        if support_voxels == 0:
            continue
        patch_values = patch[support]
        mean_delta = float(npx.mean(patch_values, dtype=npx.float64))
        mass = rho_m * dx**3 * float(npx.sum(npx.maximum(1.0 + patch_values, 0.0), dtype=npx.float64))
        if mass <= 0:
            continue
        radius_kpch = (3.0 * support_voxels * dx**3 / (4.0 * math.pi)) ** (1.0 / 3.0) * 1000.0
        velocity = _local_gradient_proxy(field, ix, iy, iz, dx)
        position = ((ix + 0.5) * dx, (iy + 0.5) * dx, (iz + 0.5) * dx)
        spin = _spin_proxy(position, velocity, field_config.field_box_size_mpch)
        b_h, fate = _field_fate(
            mass_msunh=mass,
            radius_kpch=radius_kpch,
            mean_delta=mean_delta,
            epoch=epoch,
            field_config=field_config,
        )
        rows.append(
            {
                "run_id": run_id,
                "source": "internal_3d_field",
                "snapshot": f"a{_tag(epoch.a)}",
                "z": epoch.z,
                "object_id": f"fld_z{_tag(epoch.z)}_{rank:04d}",
                "rank": rank,
                "ix": ix,
                "iy": iy,
                "iz": iz,
                "x_mpch": position[0],
                "y_mpch": position[1],
                "z_mpch": position[2],
                "peak_delta": float(field[ix, iy, iz]),
                "mean_delta_patch": mean_delta,
                "delta_contrast_to_snapshot": mean_delta / sigma,
                "support_voxels": support_voxels,
                "mass_msunh": mass,
                "radius_kpch": radius_kpch,
                "vx_proxy": velocity[0],
                "vy_proxy": velocity[1],
                "vz_proxy": velocity[2],
                "spin_x": spin[0],
                "spin_y": spin[1],
                "spin_z": spin[2],
                "b_h": b_h,
                "growth_link_mass_ratio": "",
                "fate_class": fate,
                "scale_band": _mass_band(mass),
                "quality_flag": "field_peak_patch",
            }
        )
    return rows


def _peak_candidates(field, *, threshold: float, field_config: Field3DConfig) -> list[int]:
    npx = require_numpy()
    n = field_config.grid_size
    candidate_parts = []
    value_parts = []
    per_slab_limit = max(field_config.max_peak_candidates // max(math.ceil(n / field_config.peak_slab_cells), 1), 128)
    for x0 in range(0, n, field_config.peak_slab_cells):
        x1 = min(n, x0 + field_config.peak_slab_cells)
        slab = field[x0:x1]
        values = slab.reshape(-1)
        if values.size == 0:
            continue
        limit = min(per_slab_limit, values.size)
        local = npx.argpartition(values, -limit)[-limit:]
        local_values = values[local]
        keep = local_values > threshold
        if not npx.any(keep):
            continue
        candidate_parts.append(local[keep] + x0 * n * n)
        value_parts.append(local_values[keep])
    if not candidate_parts:
        return []
    candidate_flat = npx.concatenate(candidate_parts)
    values = npx.concatenate(value_parts)
    limit = min(field_config.max_peak_candidates, len(candidate_flat))
    if len(candidate_flat) > limit:
        keep = npx.argpartition(values, -limit)[-limit:]
        candidate_flat = candidate_flat[keep]
        values = values[keep]
    values = field.ravel()[candidate_flat]
    order = npx.argsort(values)[::-1]
    selected: list[int] = []
    selected_coords: list[tuple[int, int, int]] = []
    radius2 = field_config.peak_exclusion_cells**2
    for item in candidate_flat[order]:
        coord = _unravel_index(int(item), field_config.grid_size)
        if all(_periodic_distance2(coord, other, field_config.grid_size) >= radius2 for other in selected_coords):
            selected.append(int(item))
            selected_coords.append(coord)
            if len(selected) >= field_config.max_peaks_per_snapshot:
                break
    return selected


def _transfer_rows(
    *,
    run_id: str,
    epoch,
    objects: Sequence[Mapping[str, object]],
    field_config: Field3DConfig,
) -> list[dict[str, object]]:
    if len(objects) < 2:
        return []
    rows: list[dict[str, object]] = []
    for index, target in enumerate(objects):
        neighbours = _nearest_objects(index, objects, field_config.neighbor_links, field_config.field_box_size_mpch)
        for source_index, distance in neighbours:
            source = objects[source_index]
            rows.append(_transfer_row(run_id, epoch, source, target, distance))
    return rows


def _transfer_row(run_id: str, epoch, source: Mapping[str, object], target: Mapping[str, object], distance: float) -> dict[str, object]:
    source_mass = float(source["mass_msunh"])
    target_mass = float(target["mass_msunh"])
    mass_ratio = min(source_mass, target_mass) / max(source_mass, target_mass)
    separation = (
        float(target["x_mpch"]) - float(source["x_mpch"]),
        float(target["y_mpch"]) - float(source["y_mpch"]),
        float(target["z_mpch"]) - float(source["z_mpch"]),
    )
    relative_velocity = (
        float(target["vx_proxy"]) - float(source["vx_proxy"]),
        float(target["vy_proxy"]) - float(source["vy_proxy"]),
        float(target["vz_proxy"]) - float(source["vz_proxy"]),
    )
    torque = _cross(separation, relative_velocity)
    spin = (float(target["spin_x"]), float(target["spin_y"]), float(target["spin_z"]))
    torque_norm = _norm(torque)
    spin_norm = _norm(spin)
    a_pq = mass_ratio / max(distance, 1.0e-12) ** 3
    g_pq = alignment_cosine(spin, torque) if torque_norm > 0 and spin_norm > 0 else 0.0
    return {
        "run_id": run_id,
        "source": "internal_3d_field",
        "snapshot": f"a{_tag(epoch.a)}",
        "z": epoch.z,
        "source_id": source["object_id"],
        "target_id": target["object_id"],
        "rank_source": source["rank"],
        "rank_target": target["rank"],
        "separation_mpch": distance,
        "mass_ratio": mass_ratio,
        "a_pq": a_pq,
        "g_pq": g_pq,
        "tau_x": torque[0] * a_pq,
        "tau_y": torque[1] * a_pq,
        "tau_z": torque[2] * a_pq,
        "j_x": spin[0],
        "j_y": spin[1],
        "j_z": spin[2],
        "source_channel": "field_peak_neighbour_transfer",
        "quality_flag": "nearest_peak_relation",
    }


def _mode_rows(
    *,
    run_id: str,
    epoch,
    setup: Mapping[str, object],
    realization: Mapping[str, object],
    field_config: Field3DConfig,
    objects: Sequence[Mapping[str, object]],
    transfers: Sequence[Mapping[str, object]],
) -> list[dict[str, object]]:
    npx = require_numpy()
    edges = setup["bin_edges"]
    counts = setup["bin_counts"]
    power = realization["power_by_bin"]
    active_power = realization["active_power_by_bin"]
    tail_power = realization["tail_power_by_bin"]
    active_counts = realization["active_counts"]
    tail_counts = realization["tail_counts"]
    total_power = float(npx.sum(power[1:], dtype=npx.float64))
    cumulative_power = 0.0
    delta2_values = dimensionless_power(epoch.k_hmpc, epoch.p_mpc3_h3)
    rows = []
    for index in range(1, field_config.mode_bins + 1):
        k_min = float(edges[index - 1])
        k_max = float(edges[index])
        k_mid = math.sqrt(k_min * k_max)
        sector = _activity_sector(k_min, k_max, epoch.knl.k_nl_hmpc)
        gaussian_weight = _scalar_gaussian_weight(k_mid, field_config.gaussian_sigma_mpch)
        activity_weight = _scalar_window_weight(k_mid, epoch.knl.k_nl_hmpc, float(setup["k_resolved"]), field_config.tail_width)
        fraction = float(power[index] / total_power) if total_power > 0 else 0.0
        active_fraction = float(active_power[index] / total_power) if total_power > 0 else 0.0
        tail_fraction = float(tail_power[index] / total_power) if total_power > 0 else 0.0
        delta2_mean = _interpolate_log(epoch.k_hmpc, delta2_values, k_mid)
        cumulative_power += fraction
        rows.append(
            {
                "run_id": run_id,
                "source": "internal_3d_field",
                "snapshot": f"a{_tag(epoch.a)}",
                "z": epoch.z,
                "a": epoch.a,
                "mode_bin": index,
                "gaussian_sigma_mpch": field_config.gaussian_sigma_mpch,
                "activity_endpoint_hmpc": epoch.knl.k_nl_hmpc,
                "k_min_hmpc": k_min,
                "k_max_hmpc": k_max,
                "k_mid_hmpc": k_mid,
                "mode_count": int(counts[index]),
                "active_mode_count": int(active_counts[index]),
                "tail_mode_count": int(tail_counts[index]),
                "power_proxy": float(power[index]),
                "power_fraction": fraction,
                "active_power_fraction": active_fraction,
                "tail_power_fraction": tail_fraction,
                "cumulative_power_fraction": cumulative_power,
                "delta2_mean": delta2_mean,
                "coarse_delta2_mean": delta2_mean * gaussian_weight**2,
                "gaussian_weight_mean": gaussian_weight,
                "activity_weight_mean": activity_weight,
                "window_weight_mean": activity_weight,
                "scale_band": _k_band(k_mid, epoch),
                "activity_sector": sector,
            }
        )
    return rows


def _window_row(
    *,
    run_id: str,
    epoch,
    setup: Mapping[str, object],
    realization: Mapping[str, object],
    field_config: Field3DConfig,
    objects: Sequence[Mapping[str, object]],
    transfers: Sequence[Mapping[str, object]],
    mode_rows: Sequence[Mapping[str, object]],
) -> dict[str, object]:
    npx = require_numpy()
    active_count = int(npx.sum(realization["active_counts"][1:]))
    tail_count = int(npx.sum(realization["tail_counts"][1:]))
    dominant_scale = _dominant_scale_band(mode_rows)
    dominant_relation = _dominant_relation(objects, transfers)
    s_aw = math.log(epoch.knl.k_nl_hmpc / float(setup["k_box_global"]))
    closure = realization["closure"]
    return {
        "run_id": run_id,
        "source": "internal_3d_field",
        "snapshot": f"a{_tag(epoch.a)}",
        "z": epoch.z,
        "a": epoch.a,
        "grid_size": field_config.grid_size,
        "voxel_count": field_config.grid_size**3,
        "field_box_size_mpch": field_config.field_box_size_mpch,
        "global_box_size_mpch": field_config.global_box_size_mpch,
        "dx_mpch": setup["dx"],
        "gaussian_sigma_mpch": field_config.gaussian_sigma_mpch,
        "gaussian_kernel": "exp[-0.5*(sigma*k)^2]",
        "activity_endpoint_hmpc": epoch.knl.k_nl_hmpc,
        "activity_selector": "log_gaussian_tail_with_cut",
        "k_box_field_hmpc": setup["k_box_field"],
        "k_box_global_hmpc": setup["k_box_global"],
        "k_nyquist_hmpc": setup["k_nyquist"],
        "k_cut_hmpc": field_config.k_cut_hmpc,
        "k_resolved_hmpc": setup["k_resolved"],
        "c_nl": epoch.knl.c_nl,
        "k_nl_hmpc": epoch.knl.k_nl_hmpc,
        "m_nl_msunh": epoch.m_nl_msun_h,
        "s_aw_global": s_aw,
        "target_sigma_delta": realization["target_sigma"],
        "measured_sigma_delta": realization["measured_sigma"],
        "active_mode_count": active_count,
        "tail_mode_count": tail_count,
        "active_power_fraction": realization["active_power_fraction"],
        "tail_power_fraction": realization["tail_power_fraction"],
        "closure_active_norm": closure["active_norm"],
        "closure_model_norm": closure["tail_model_norm"],
        "closure_residual_norm": closure["residual_norm"],
        "closure_error": closure["closure_error"],
        "closure_status": closure["quality_flag"],
        "field_mean_delta": realization["field_mean"],
        "field_min_delta": realization["field_min"],
        "field_max_delta": realization["field_max"],
        "field_skewness": realization["field_skewness"],
        "field_kurtosis": realization["field_kurtosis"],
        "peak_count": len(objects),
        "transfer_count": len(transfers),
        "dominant_scale_band": dominant_scale,
        "dominant_relation": dominant_relation,
        "quality_flag": "relationship_tables_only",
    }


def _relation_summary_rows(
    run_id: str,
    epoch,
    objects: Sequence[Mapping[str, object]],
    transfers: Sequence[Mapping[str, object]],
) -> list[dict[str, object]]:
    rows = []
    bands = sorted({_mass_band(float(row["mass_msunh"])) for row in objects})
    for band in bands:
        selected = [row for row in objects if _mass_band(float(row["mass_msunh"])) == band]
        linked = [
            row
            for row in transfers
            if _object_band(row["target_id"], objects) == band
        ]
        masses = [float(row["mass_msunh"]) for row in selected]
        b_values = [float(row["b_h"]) for row in selected if math.isfinite(float(row["b_h"]))]
        a_values = [float(row["a_pq"]) for row in linked]
        g_values = [abs(float(row["g_pq"])) for row in linked]
        classes = {name: sum(1 for row in selected if row["fate_class"] == name) for name in ("forming", "balanced", "stripping")}
        n = len(selected)
        rows.append(
            {
                "run_id": run_id,
                "source": "internal_3d_field",
                "snapshot": f"a{_tag(epoch.a)}",
                "z": epoch.z,
                "a": epoch.a,
                "scale_band": band,
                "object_count": n,
                "mass_total_msunh": sum(masses),
                "mass_median_msunh": _quantile(masses, 0.50),
                "mass_p16_msunh": _quantile(masses, 0.16),
                "mass_p84_msunh": _quantile(masses, 0.84),
                "b_h_median": _quantile(b_values, 0.50),
                "b_h_p16": _quantile(b_values, 0.16),
                "b_h_p84": _quantile(b_values, 0.84),
                "forming_fraction": classes["forming"] / n if n else 0.0,
                "balanced_fraction": classes["balanced"] / n if n else 0.0,
                "stripping_fraction": classes["stripping"] / n if n else 0.0,
                "transfer_count": len(linked),
                "a_pq_median": _quantile(a_values, 0.50),
                "a_pq_p16": _quantile(a_values, 0.16),
                "a_pq_p84": _quantile(a_values, 0.84),
                "abs_g_pq_median": _quantile(g_values, 0.50),
                "abs_g_pq_p16": _quantile(g_values, 0.16),
                "abs_g_pq_p84": _quantile(g_values, 0.84),
                "relation_status": "field_peak_relation",
            }
        )
    return rows


def _attach_growth_from_previous(
    objects: list[dict[str, object]],
    previous_by_band: Mapping[str, Sequence[Mapping[str, object]]],
    step_gyr: float,
) -> None:
    for row in objects:
        band = str(row["scale_band"])
        previous = previous_by_band.get(band, ())
        if not previous:
            row["growth_link_mass_ratio"] = ""
            continue
        rank = min(int(row["rank"]) - 1, len(previous) - 1)
        previous_mass = float(previous[rank]["mass_msunh"])
        current_mass = float(row["mass_msunh"])
        row["growth_link_mass_ratio"] = current_mass / previous_mass if previous_mass > 0 else ""
        if step_gyr > 0 and previous_mass > 0:
            delta = current_mass - previous_mass
            if delta != 0.0:
                t_growth = current_mass * step_gyr / delta
                t_dyn = dynamical_time_gyr(current_mass, float(row["radius_kpch"]))
                row["b_h"] = t_dyn / t_growth
                row["fate_class"] = _class_from_bh(float(row["b_h"]), delta)


def _objects_by_band(objects: Sequence[Mapping[str, object]]) -> dict[str, list[Mapping[str, object]]]:
    output: dict[str, list[Mapping[str, object]]] = {}
    for row in objects:
        output.setdefault(str(row["scale_band"]), []).append(row)
    for rows in output.values():
        rows.sort(key=lambda item: int(item["rank"]))
    return output


def _epoch_relation_row(run_id: str, previous: Mapping[str, object], current: Mapping[str, object]) -> dict[str, object]:
    a_from = float(previous["a"])
    a_to = float(current["a"])
    k_from = float(previous["k_nl_hmpc"])
    k_to = float(current["k_nl_hmpc"])
    delta_ln_a = math.log(a_to / a_from) if a_from > 0 and a_to > 0 else ""
    delta_ln_k = math.log(k_to / k_from) if k_from > 0 and k_to > 0 else ""
    endpoint_velocity = delta_ln_k / delta_ln_a if isinstance(delta_ln_a, float) and delta_ln_a != 0.0 else ""
    scale_from = str(previous["dominant_scale_band"])
    scale_to = str(current["dominant_scale_band"])
    return {
        "run_id": run_id,
        "source": "internal_3d_field",
        "snapshot_from": previous["snapshot"],
        "snapshot_to": current["snapshot"],
        "z_from": previous["z"],
        "z_to": current["z"],
        "a_from": a_from,
        "a_to": a_to,
        "delta_ln_a": delta_ln_a,
        "delta_ln_k_nl": delta_ln_k,
        "endpoint_velocity": endpoint_velocity,
        "delta_sigma": float(current["measured_sigma_delta"]) - float(previous["measured_sigma_delta"]),
        "delta_peak_count": int(current["peak_count"]) - int(previous["peak_count"]),
        "delta_transfer_count": int(current["transfer_count"]) - int(previous["transfer_count"]),
        "dominant_scale_from": scale_from,
        "dominant_scale_to": scale_to,
        "relation_transition": "same_scale" if scale_from == scale_to else f"{scale_from}_to_{scale_to}",
        "quality_flag": "adjacent_epoch_relation",
    }


def _closure_row(run_id: str, epoch, realization: Mapping[str, object], field_config: Field3DConfig) -> dict[str, object]:
    closure = realization["closure"]
    return {
        "run_id": run_id,
        "source": "internal_3d_field",
        "snapshot": f"a{_tag(epoch.a)}",
        "z": epoch.z,
        "a": epoch.a,
        "gaussian_sigma_mpch": field_config.gaussian_sigma_mpch,
        "activity_endpoint_hmpc": epoch.knl.k_nl_hmpc,
        "closure_map": "spectral_tail_moment_proxy",
        "active_norm": closure["active_norm"],
        "tail_model_norm": closure["tail_model_norm"],
        "residual_norm": closure["residual_norm"],
        "closure_error": closure["closure_error"],
        "epsilon0": field_config.closure_floor,
        "norm": "log_power_integral",
        "quality_flag": closure["quality_flag"],
    }


def _field_fate(
    *,
    mass_msunh: float,
    radius_kpch: float,
    mean_delta: float,
    epoch,
    field_config: Field3DConfig,
) -> tuple[float, str]:
    t_dyn = dynamical_time_gyr(mass_msunh, radius_kpch)
    local_activity = max(mean_delta, 1.0e-6)
    t_growth = field_config.object_step_gyr / local_activity
    b_h = t_dyn / t_growth
    return b_h, _class_from_bh(b_h, mean_delta)


def _class_from_bh(b_h: float, signed_growth: float) -> str:
    if signed_growth < 0:
        return "stripping"
    if not math.isfinite(b_h):
        return "balanced"
    if b_h > 1.0:
        return "forming"
    return "balanced"


def _nearest_objects(
    index: int,
    objects: Sequence[Mapping[str, object]],
    count: int,
    box_size: float,
) -> list[tuple[int, float]]:
    target = objects[index]
    distances = []
    for other_index, other in enumerate(objects):
        if other_index == index:
            continue
        distance = _periodic_physical_distance(target, other, box_size)
        distances.append((other_index, distance))
    distances.sort(key=lambda item: item[1])
    return distances[:count]


def _periodic_physical_distance(a: Mapping[str, object], b: Mapping[str, object], box_size: float) -> float:
    diffs = []
    for key in ("x_mpch", "y_mpch", "z_mpch"):
        raw = abs(float(a[key]) - float(b[key]))
        diffs.append(min(raw, box_size - raw))
    return math.sqrt(sum(value * value for value in diffs))


def _patch(field, ix: int, iy: int, iz: int, radius: int):
    npx = require_numpy()
    n = field.shape[0]
    axes = [npx.arange(index - radius, index + radius + 1) % n for index in (ix, iy, iz)]
    return field[npx.ix_(axes[0], axes[1], axes[2])]


def _local_gradient_proxy(field, ix: int, iy: int, iz: int, dx: float) -> tuple[float, float, float]:
    n = field.shape[0]
    values = []
    for axis, index in enumerate((ix, iy, iz)):
        left = [ix, iy, iz]
        right = [ix, iy, iz]
        left[axis] = (index - 1) % n
        right[axis] = (index + 1) % n
        values.append(-(float(field[tuple(right)]) - float(field[tuple(left)])) / (2.0 * dx))
    return (values[0], values[1], values[2])


def _spin_proxy(
    position: Sequence[float],
    velocity: Sequence[float],
    box_size: float,
) -> tuple[float, float, float]:
    centered = tuple(float(value) - 0.5 * box_size for value in position)
    spin = _cross(centered, velocity)
    norm = _norm(spin)
    if norm <= 0:
        return (0.0, 0.0, 1.0)
    return (spin[0] / norm, spin[1] / norm, spin[2] / norm)


def _dominant_scale_band(mode_rows: Sequence[Mapping[str, object]]) -> str:
    if not mode_rows:
        return ""
    row = max(mode_rows, key=lambda item: float(item["power_fraction"]))
    return str(row["scale_band"])


def _dominant_relation(objects: Sequence[Mapping[str, object]], transfers: Sequence[Mapping[str, object]]) -> str:
    if not objects:
        return "field_window_only"
    if not transfers:
        return "field_objects"
    median_a = _quantile([float(row["a_pq"]) for row in transfers], 0.50)
    median_g = _quantile([abs(float(row["g_pq"])) for row in transfers], 0.50)
    if median_a > 0.0 and median_g >= 0.5:
        return "coherent_transfer"
    if median_a > 0.0:
        return "local_transfer"
    return "object_population"


def _mass_band(mass_msunh: float) -> str:
    for band_id, mass_min, mass_max, _label in DEFAULT_MASS_BANDS:
        if mass_min <= mass_msunh < mass_max:
            return band_id
    if mass_msunh < DEFAULT_MASS_BANDS[0][1]:
        return "sub_dwarf_field_peak"
    return "super_cluster_field_peak"


def _k_band(k_hmpc: float, epoch) -> str:
    mapper_density = epoch.m_nl_msun_h / ((4.0 * math.pi / 3.0) * (epoch.r_nl_kpc_h / 1000.0) ** 3)
    for band_id, mass_min, mass_max, _label in DEFAULT_MASS_BANDS:
        r_min = mass_to_radius_mpch(mass_min, mapper_density)
        r_max = mass_to_radius_mpch(mass_max, mapper_density)
        k_min = epoch.r_nl_kpc_h / 1000.0 * epoch.knl.k_nl_hmpc / r_max
        k_max = epoch.r_nl_kpc_h / 1000.0 * epoch.knl.k_nl_hmpc / r_min
        if k_min <= k_hmpc <= k_max:
            return band_id
    return "outside_named_mass_bands"


def _activity_sector(k_min: float, k_max: float, k_nl: float) -> str:
    if k_min <= k_nl <= k_max:
        return "frontier"
    if k_max < k_nl:
        return "active"
    return "tail"


def _object_band(object_id: object, objects: Sequence[Mapping[str, object]]) -> str:
    for row in objects:
        if row["object_id"] == object_id:
            return str(row["scale_band"])
    return ""


def _scalar_window_weight(k_value: float, k_nl: float, k_cut: float, tail_width: float) -> float:
    if k_value <= 0.0 or k_value > k_cut:
        return 0.0
    if k_value <= k_nl:
        return 1.0
    return math.exp(-0.5 * (math.log(k_value / max(k_nl, 1.0e-30)) / tail_width) ** 2)


def _scalar_gaussian_weight(k_value: float, gaussian_sigma_mpch: float) -> float:
    if gaussian_sigma_mpch <= 0:
        return 1.0
    return math.exp(-0.5 * (gaussian_sigma_mpch * k_value) ** 2)


def _interpolate_log(x_values: Sequence[float], y_values: Sequence[float], x: float) -> float:
    if x <= x_values[0]:
        return float(y_values[0])
    if x >= x_values[-1]:
        return float(y_values[-1])
    lo = 0
    hi = len(x_values) - 1
    log_x = math.log(x)
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if x_values[mid] <= x:
            lo = mid
        else:
            hi = mid
    x0 = math.log(x_values[lo])
    x1 = math.log(x_values[hi])
    y0 = float(y_values[lo])
    y1 = float(y_values[hi])
    t = (log_x - x0) / (x1 - x0)
    return y0 * (1.0 - t) + y1 * t


def _trapz(y, x) -> float:
    npx = require_numpy()
    if len(y) < 2:
        return 0.0
    return float(npx.sum(0.5 * (y[1:] + y[:-1]) * (x[1:] - x[:-1]), dtype=npx.float64))


def _masked_log_integral(values, k_values, mask) -> float:
    npx = require_numpy()
    selected = npx.asarray(mask, dtype=bool)
    if npx.count_nonzero(selected) < 2:
        return 0.0
    return _trapz(values[selected], npx.log(k_values[selected]))


def _field_moment_ratio(field, order: int) -> float:
    npx = require_numpy()
    sigma = float(npx.std(field, dtype=npx.float64))
    if sigma <= 0:
        return 0.0
    total = 0.0
    count = 0
    slab_cells = 16
    for start in range(0, field.shape[0], slab_cells):
        slab = field[start : start + slab_cells].astype(npx.float64, copy=False) / sigma
        total += float(npx.sum(npx.power(slab, order), dtype=npx.float64))
        count += int(slab.size)
    return total / count if count else 0.0


def _quantile(values: Sequence[float], q: float) -> float | str:
    finite = sorted(value for value in values if math.isfinite(value))
    if not finite:
        return ""
    if q <= 0:
        return finite[0]
    if q >= 1:
        return finite[-1]
    position = q * (len(finite) - 1)
    lo = math.floor(position)
    hi = math.ceil(position)
    if lo == hi:
        return finite[lo]
    frac = position - lo
    return finite[lo] * (1.0 - frac) + finite[hi] * frac


def _periodic_distance2(a: Sequence[int], b: Sequence[int], n: int) -> int:
    total = 0
    for left, right in zip(a, b):
        distance = abs(left - right)
        distance = min(distance, n - distance)
        total += distance * distance
    return total


def _unravel_index(flat_index: int, n: int) -> tuple[int, int, int]:
    plane = n * n
    ix = flat_index // plane
    rem = flat_index - ix * plane
    iy = rem // n
    iz = rem - iy * n
    return ix, iy, iz


def _cross(a: Sequence[float], b: Sequence[float]) -> tuple[float, float, float]:
    return (
        float(a[1]) * float(b[2]) - float(a[2]) * float(b[1]),
        float(a[2]) * float(b[0]) - float(a[0]) * float(b[2]),
        float(a[0]) * float(b[1]) - float(a[1]) * float(b[0]),
    )


def _norm(values: Sequence[float]) -> float:
    return math.sqrt(sum(float(value) * float(value) for value in values))


def _tag(value: float) -> str:
    return f"{value:g}".replace("-", "m").replace(".", "p")
