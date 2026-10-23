"""Internal finite-window generation of cosmic structure histories."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import math
from typing import Iterable, Sequence

from cosmic_cg.halo_model import endpoint_virial_velocity, virial_radius_from_mass
from cosmic_cg.knl import KNLResult, find_knl
from cosmic_cg.object_level import alignment_cosine, halo_fate
from cosmic_cg.params import CosmicCGParams
from cosmic_cg.result_tables import DEFAULT_MASS_BANDS
from cosmic_cg.tail_ratio import virial_velocity


RHO_CRIT_MSUNH_MPCH3 = 2.77536627e11


@dataclass(frozen=True)
class InternalSpectrumConfig:
    c_nl: float = 1.0
    k_nl0_hmpc: float = 0.2
    n_s: float = 0.965
    k_pivot_hmpc: float = 0.2
    k_min_hmpc: float = 1.0e-4
    k_max_hmpc: float = 200.0
    k_points: int = 768
    endpoint_radius_factor: float = 1.0
    spectrum_shape: str = "finite_window_envelope"
    uv_envelope_k_hmpc: float = 20.0
    uv_envelope_alpha: float = 2.0
    uv_envelope_beta: float = 2.0


@dataclass(frozen=True)
class GeneratedEpoch:
    z: float
    a: float
    k_hmpc: tuple[float, ...]
    p_mpc3_h3: tuple[float, ...]
    knl: KNLResult
    m_nl_msun_h: float
    r_nl_kpc_h: float
    v_vir_mnl_kms: float


def generate_epoch(
    *,
    z: float,
    params: CosmicCGParams,
    config: InternalSpectrumConfig,
) -> GeneratedEpoch:
    k_values = tuple(logspace(config.k_min_hmpc, config.k_max_hmpc, config.k_points))
    p_values = tuple(internal_power_spectrum(k_values, z=z, params=params, config=config))
    knl = find_knl(k_values, p_values, z=z, c_nl=config.c_nl)
    omega_m = omega_m0(params)
    r_nl_mpch = config.endpoint_radius_factor / knl.k_nl_hmpc
    m_nl = nonlinear_mass_from_radius(r_nl_mpch, omega_m=omega_m)
    r_nl_kpc_h, v_nl = endpoint_virial_velocity(m_nl, omega_m=omega_m, h=hubble_h(params))
    return GeneratedEpoch(
        z=z,
        a=1.0 / (1.0 + z),
        k_hmpc=k_values,
        p_mpc3_h3=p_values,
        knl=knl,
        m_nl_msun_h=m_nl,
        r_nl_kpc_h=r_nl_kpc_h,
        v_vir_mnl_kms=v_nl,
    )


def internal_power_spectrum(
    k_hmpc: Sequence[float],
    *,
    z: float,
    params: CosmicCGParams,
    config: InternalSpectrumConfig,
) -> list[float]:
    amplitude = _normalization_amplitude(config)
    growth = growth_factor(z, params=params)
    return [amplitude * growth**2 * internal_spectrum_shape(k, config=config) for k in k_hmpc]


def generate_objects(
    *,
    run_id: str,
    epoch: GeneratedEpoch,
    params: CosmicCGParams,
    objects_per_band: int,
    object_step_gyr: float,
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    objects: list[dict[str, object]] = []
    growth_links: list[dict[str, object]] = []
    transfers: list[dict[str, object]] = []
    omega_m = omega_m0(params)
    h = hubble_h(params)
    time_to = cosmic_age_gyr(epoch.z, params=params)
    time_from = max(0.0, time_to - object_step_gyr)
    for band_id, mass_min, mass_max, _label in DEFAULT_MASS_BANDS:
        for index, mass in enumerate(_sample_masses(mass_min, mass_max, objects_per_band)):
            halo_id = f"z{_tag(epoch.z)}_{band_id}_{index:03d}"
            radius = virial_radius_from_mass(mass, omega_m=omega_m, h=h)
            v_vir = virial_velocity(mass, radius)
            activity = _activity_weight(mass, epoch.m_nl_msun_h)
            signed_growth = _signed_growth_rate(epoch.z, activity)
            mass_from = mass * math.exp(-signed_growth * object_step_gyr)
            fate = halo_fate(
                halo_id=halo_id,
                mass_msun_h=mass,
                radius_kpc_h=radius,
                delta_mass_msun_h=mass - mass_from,
                delta_time_gyr=object_step_gyr,
            )
            phase = _phase_value(epoch.z, band_id, index)
            sigma = max(1.0e-12, v_vir / math.sqrt(2.0) * (0.8 + 0.25 * activity + 0.05 * math.sin(phase)))
            position = _position(epoch.a, band_id, index)
            velocity = _velocity(v_vir, phase)
            spin = _spin(phase)
            r_h_object = sigma / epoch.v_vir_mnl_kms
            object_row = {
                "run_id": run_id,
                "source": "internal_generation",
                "snapshot": f"a{_tag(epoch.a)}",
                "z": epoch.z,
                "halo_id": halo_id,
                "mass_msun_h": mass,
                "radius_kpc_h": radius,
                "v_vir_kms": v_vir,
                "sigma_v_kms": sigma,
                "sigma_source": "internal_virial_dispersion",
                "x": position[0],
                "y": position[1],
                "z_pos": position[2],
                "vx": velocity[0],
                "vy": velocity[1],
                "vz": velocity[2],
                "spin_x": spin[0],
                "spin_y": spin[1],
                "spin_z": spin[2],
                "t_dyn_gyr": fate.t_dyn_gyr,
                "t_growth_gyr": "" if fate.t_growth_gyr is None else fate.t_growth_gyr,
                "B_h": fate.b_h,
                "B_class": fate.class_name,
                "parent_id": "",
                "descendant_id": halo_id,
                "r_h_object": r_h_object,
                "branch_weight": activity,
                "link_quality": "internal_generation",
                "quality_flag": "internal_generation",
            }
            objects.append(object_row)
            growth_link = {
                "source": "internal_generation",
                "snapshot_from": f"a{_tag(scale_factor_from_age(time_from, params=params))}",
                "snapshot_to": f"a{_tag(epoch.a)}",
                "z_from": redshift_from_age(time_from, params=params),
                "z_to": epoch.z,
                "time_from_gyr": time_from,
                "time_to_gyr": time_to,
                "halo_id_from": halo_id,
                "halo_id_to": halo_id,
                "mass_from_msun_h": mass_from,
                "mass_to_msun_h": mass,
                "quality_flag": "internal_generation",
            }
            growth_links.append(growth_link)
            transfers.append(_transfer_row(run_id, epoch, halo_id, spin, position, velocity, activity, mass, mass_from))
    transfers.extend(_pair_transfer_rows(run_id, epoch, objects))
    return objects, growth_links, transfers


def nonlinear_mass_from_radius(r_mpch: float, *, omega_m: float) -> float:
    if r_mpch <= 0:
        raise ValueError("radius must be positive")
    return (4.0 * math.pi / 3.0) * RHO_CRIT_MSUNH_MPCH3 * omega_m * r_mpch**3


def growth_factor(z: float, *, params: CosmicCGParams) -> float:
    if z <= -1.0:
        raise ValueError("redshift must be greater than -1")
    a = 1.0 / (1.0 + z)
    return _linear_growth_unnormalized(a, params=params) / _linear_growth_unnormalized(1.0, params=params)


def cosmic_age_gyr(z: float, *, params: CosmicCGParams) -> float:
    if z <= -1.0:
        raise ValueError("redshift must be greater than -1")
    omega_m = omega_m0(params)
    omega_l = max(params.omega_lambda, 1.0e-12)
    a = 1.0 / (1.0 + z)
    from cosmic_cg.units import HBAR_GEV_S, SECONDS_PER_GYR
    h0_inv_gyr = HBAR_GEV_S / (params.h0_gev * SECONDS_PER_GYR)
    return 2.0 / (3.0 * math.sqrt(omega_l)) * math.asinh(math.sqrt(omega_l / omega_m) * a**1.5) * h0_inv_gyr


def redshift_from_age(age_gyr: float, *, params: CosmicCGParams) -> float:
    if age_gyr <= 0:
        return 1.0e6
    lo, hi = -0.95, 1.0e6
    for _ in range(160):
        mid = 0.5 * (lo + hi)
        if cosmic_age_gyr(mid, params=params) > age_gyr:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def scale_factor_from_age(age_gyr: float, *, params: CosmicCGParams) -> float:
    z = redshift_from_age(age_gyr, params=params)
    return 1.0 / (1.0 + z)


def omega_m0(params: CosmicCGParams) -> float:
    return params.omega_b + params.omega_sigma


def hubble_h(params: CosmicCGParams) -> float:
    from cosmic_cg.units import HBAR_GEV_S, MPC_IN_KM
    if not math.isfinite(params.h0_gev) or params.h0_gev <= 0:
        raise ValueError("The Hubble scale must be finite and positive")
    return params.h0_gev * MPC_IN_KM / (100.0 * HBAR_GEV_S)


def config_to_dict(config: InternalSpectrumConfig) -> dict[str, object]:
    return asdict(config)


def logspace(start: float, stop: float, count: int) -> list[float]:
    if start <= 0 or stop <= start:
        raise ValueError("logspace endpoints must be positive and increasing")
    if count < 2:
        raise ValueError("logspace requires at least two samples")
    log_start = math.log(start)
    step = (math.log(stop) - log_start) / (count - 1)
    return [math.exp(log_start + i * step) for i in range(count)]


def internal_spectrum_shape(k_hmpc: float, *, config: InternalSpectrumConfig) -> float:
    """Return the dimensionless internal spectral shape at a wavenumber."""

    if k_hmpc <= 0 or config.k_pivot_hmpc <= 0:
        raise ValueError("wavenumbers must be positive")
    base = (k_hmpc / config.k_pivot_hmpc) ** config.n_s
    if config.spectrum_shape == "power_law":
        return base
    if config.spectrum_shape == "finite_window_envelope":
        return base * finite_window_envelope(
            k_hmpc,
            k_uv_hmpc=config.uv_envelope_k_hmpc,
            alpha=config.uv_envelope_alpha,
            beta=config.uv_envelope_beta,
        )
    raise ValueError(f"unknown spectrum shape: {config.spectrum_shape}")


def finite_window_envelope(
    k_hmpc: float,
    *,
    k_uv_hmpc: float,
    alpha: float,
    beta: float,
) -> float:
    """Return the analytic high-k envelope used by the internal generator."""

    if k_hmpc <= 0 or k_uv_hmpc <= 0:
        raise ValueError("wavenumbers must be positive")
    if alpha <= 0 or beta <= 0:
        raise ValueError("envelope parameters must be positive")
    return (1.0 + (k_hmpc / k_uv_hmpc) ** alpha) ** (-beta)


def spectrum_normalization_amplitude(config: InternalSpectrumConfig) -> float:
    return _normalization_amplitude(config)


def _normalization_amplitude(config: InternalSpectrumConfig) -> float:
    target_p = config.c_nl * 2.0 * math.pi**2 / config.k_nl0_hmpc**3
    shape = internal_spectrum_shape(config.k_nl0_hmpc, config=config)
    return target_p / shape


def _linear_growth_unnormalized(a: float, *, params: CosmicCGParams) -> float:
    if a <= 0:
        raise ValueError("scale factor must be positive")
    omega_m = omega_m0(params)
    omega_l = params.omega_lambda
    return 2.5 * omega_m * _e_of_a(a, omega_m=omega_m, omega_l=omega_l) * _growth_integral(a, omega_m=omega_m, omega_l=omega_l)


def _growth_integral(a: float, *, omega_m: float, omega_l: float) -> float:
    lower = 1.0e-8
    if a <= lower:
        return 2.0 * a**2.5 / (5.0 * omega_m**1.5)
    n = 512
    x0 = math.log(lower)
    x1 = math.log(a)
    step = (x1 - x0) / n
    total = 0.0
    for index in range(n + 1):
        x = x0 + index * step
        aa = math.exp(x)
        weight = 0.5 if index in (0, n) else 1.0
        total += weight / (aa**2 * _e_of_a(aa, omega_m=omega_m, omega_l=omega_l) ** 3)
    return total * step


def _e_of_a(a: float, *, omega_m: float, omega_l: float) -> float:
    return math.sqrt(omega_m / a**3 + omega_l)


def _sample_masses(mass_min: float, mass_max: float, count: int) -> Iterable[float]:
    if count <= 0:
        return []
    log_min = math.log(mass_min)
    step = (math.log(mass_max) - log_min) / count
    return [math.exp(log_min + (index + 0.5) * step) for index in range(count)]


def _activity_weight(mass: float, m_nl: float) -> float:
    distance = abs(math.log(mass / m_nl))
    return math.exp(-0.5 * distance)


def _signed_growth_rate(z: float, activity: float) -> float:
    if z < 0:
        return -0.04 * abs(z) / (1.0 + z) * (0.4 + activity)
    return 0.08 * (1.0 + z) ** 0.6 * activity


def _phase_value(z: float, band_id: str, index: int) -> float:
    return 0.37 * index + 0.13 * len(band_id) + 0.07 * z


def _position(a: float, band_id: str, index: int) -> tuple[float, float, float]:
    phase = _phase_value(a, band_id, index)
    radius = 1000.0 * (1.0 + 0.15 * len(band_id))
    return (
        radius * math.sin(phase),
        radius * math.cos(phase * 1.7),
        radius * math.sin(phase * 0.7),
    )


def _velocity(v_vir: float, phase: float) -> tuple[float, float, float]:
    scale = v_vir / math.sqrt(3.0)
    return (
        scale * math.cos(phase),
        scale * math.sin(phase * 1.3),
        scale * math.cos(phase * 0.6),
    )


def _spin(phase: float) -> tuple[float, float, float]:
    return (
        math.cos(phase),
        math.sin(phase),
        math.cos(0.5 * phase),
    )


def _transfer_row(
    run_id: str,
    epoch: GeneratedEpoch,
    halo_id: str,
    spin: Sequence[float],
    position: Sequence[float],
    velocity: Sequence[float],
    activity: float,
    mass: float,
    mass_from: float,
) -> dict[str, object]:
    torque = _cross(position, velocity)
    torque_norm = math.sqrt(sum(value * value for value in torque))
    mass_fraction = abs(mass - mass_from) / mass
    a_pq = activity * mass_fraction
    g_pq = alignment_cosine(spin, torque) if torque_norm > 0 else 0.0
    return {
        "run_id": run_id,
        "source": "internal_generation",
        "snapshot": f"a{_tag(epoch.a)}",
        "z": epoch.z,
        "source_id": f"{halo_id}_prev",
        "target_id": halo_id,
        "source_channel": "generated_branch_transfer",
        "a_pq": a_pq,
        "g_pq": g_pq,
        "tau_x": torque[0] * a_pq,
        "tau_y": torque[1] * a_pq,
        "tau_z": torque[2] * a_pq,
        "j_x": spin[0],
        "j_y": spin[1],
        "j_z": spin[2],
        "tensor_source": "internal_window_tensor",
        "spin_source": "internal_spin",
        "vector_status": "internal_generation",
    }


def _pair_transfer_rows(
    run_id: str,
    epoch: GeneratedEpoch,
    objects: Sequence[dict[str, object]],
) -> list[dict[str, object]]:
    if len(objects) < 2:
        return []
    ordered = sorted(objects, key=lambda row: (float(row["mass_msun_h"]), str(row["halo_id"])))
    rows = []
    for index, target in enumerate(ordered):
        source = ordered[index - 1]
        rows.append(_pair_transfer_row(run_id, epoch, source, target))
    return rows


def _pair_transfer_row(
    run_id: str,
    epoch: GeneratedEpoch,
    source: dict[str, object],
    target: dict[str, object],
) -> dict[str, object]:
    source_position = (float(source["x"]), float(source["y"]), float(source["z_pos"]))
    target_position = (float(target["x"]), float(target["y"]), float(target["z_pos"]))
    source_velocity = (float(source["vx"]), float(source["vy"]), float(source["vz"]))
    target_velocity = (float(target["vx"]), float(target["vy"]), float(target["vz"]))
    separation = tuple(target_position[i] - source_position[i] for i in range(3))
    relative_velocity = tuple(target_velocity[i] - source_velocity[i] for i in range(3))
    torque = _cross(separation, relative_velocity)
    spin = (float(target["spin_x"]), float(target["spin_y"]), float(target["spin_z"]))
    source_mass = float(source["mass_msun_h"])
    target_mass = float(target["mass_msun_h"])
    mass_ratio = min(source_mass, target_mass) / max(source_mass, target_mass)
    activity = math.sqrt(float(source["branch_weight"]) * float(target["branch_weight"]))
    distance = max(_norm(separation), 1.0e-12)
    speed = _norm(relative_velocity)
    sigma_sum = max(float(source["sigma_v_kms"]) + float(target["sigma_v_kms"]), 1.0e-12)
    coupling = activity * mass_ratio / (1.0 + distance / 2000.0)
    a_pq = coupling * speed / sigma_sum
    g_pq = alignment_cosine(spin, torque) if _norm(torque) > 0 else 0.0
    return {
        "run_id": run_id,
        "source": "internal_generation",
        "snapshot": f"a{_tag(epoch.a)}",
        "z": epoch.z,
        "source_id": source["halo_id"],
        "target_id": target["halo_id"],
        "source_channel": "generated_pair_tidal_transfer",
        "a_pq": a_pq,
        "g_pq": g_pq,
        "tau_x": torque[0] * a_pq,
        "tau_y": torque[1] * a_pq,
        "tau_z": torque[2] * a_pq,
        "j_x": spin[0],
        "j_y": spin[1],
        "j_z": spin[2],
        "tensor_source": "internal_pair_tidal_tensor",
        "spin_source": "internal_spin",
        "vector_status": "internal_generation",
    }


def _cross(a: Sequence[float], b: Sequence[float]) -> tuple[float, float, float]:
    return (
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    )


def _norm(values: Sequence[float]) -> float:
    return math.sqrt(sum(float(value) * float(value) for value in values))


def _tag(value: float) -> str:
    return f"{value:g}".replace("-", "m").replace(".", "p")
