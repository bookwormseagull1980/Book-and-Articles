"""Tail-activity diagnostics based on generated object ensembles."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Iterable, Mapping

from .units import G_KPC_KM2_S2_MSUN


@dataclass(frozen=True)
class TailRatioResult:
    r_h: float
    v_rms_tail_kms: float
    v_vir_mnl_kms: float
    n_halos: int
    median_per_halo: float
    p16_per_halo: float
    p84_per_halo: float
    sigma_source: str


def virial_velocity(mass_msun_h: float, radius_kpc_h: float) -> float:
    """Return virial velocity in km/s for matching h-scaled mass and radius."""

    if mass_msun_h <= 0:
        raise ValueError("mass must be positive")
    if radius_kpc_h <= 0:
        raise ValueError("radius must be positive")
    return math.sqrt(G_KPC_KM2_S2_MSUN * mass_msun_h / radius_kpc_h)


def halo_tail_ratio(
    halos: Iterable[Mapping[str, float]],
    *,
    m_nl_msun_h: float,
    r_nl_kpc_h: float,
    mass_key: str = "mass_msun_h",
    sigma_key: str = "sigma_v_kms",
    vmax_key: str = "vmax_kms",
    radius_key: str = "radius_kpc_h",
    max_mass_factor: float = 1.0,
) -> TailRatioResult:
    """Compute the halo tail-activity ratio R_h."""

    if max_mass_factor <= 0:
        raise ValueError("max_mass_factor must be positive")
    selected = []
    sigma_source = sigma_key
    for halo in halos:
        mass = float(halo[mass_key])
        if mass <= 0 or mass > max_mass_factor * m_nl_msun_h:
            continue
        sigma = _velocity_dispersion(halo, sigma_key, vmax_key, mass_key, radius_key)
        source = _velocity_source(halo, sigma_key, vmax_key)
        if source != sigma_key:
            sigma_source = source
        selected.append((mass, sigma))

    if not selected:
        raise ValueError("no halos selected for R_h calculation")

    numerator = sum(mass * sigma * sigma for mass, sigma in selected)
    denominator = sum(mass for mass, _ in selected)
    v_rms = math.sqrt(numerator / denominator)
    v_ref = virial_velocity(m_nl_msun_h, r_nl_kpc_h)
    ratios = sorted(sigma / v_ref for _, sigma in selected)
    return TailRatioResult(
        r_h=v_rms / v_ref,
        v_rms_tail_kms=v_rms,
        v_vir_mnl_kms=v_ref,
        n_halos=len(selected),
        median_per_halo=_quantile(ratios, 0.50),
        p16_per_halo=_quantile(ratios, 0.16),
        p84_per_halo=_quantile(ratios, 0.84),
        sigma_source=sigma_source,
    )


def _velocity_dispersion(
    halo: Mapping[str, float],
    sigma_key: str,
    vmax_key: str,
    mass_key: str,
    radius_key: str,
) -> float:
    if sigma_key in halo and halo[sigma_key] is not None:
        value = float(halo[sigma_key])
    elif vmax_key in halo and halo[vmax_key] is not None:
        value = float(halo[vmax_key]) / math.sqrt(2.0)
    else:
        value = virial_velocity(float(halo[mass_key]), float(halo[radius_key])) / math.sqrt(2.0)
    if value <= 0:
        raise ValueError("velocity dispersion proxy must be positive")
    return value


def _velocity_source(halo: Mapping[str, float], sigma_key: str, vmax_key: str) -> str:
    if sigma_key in halo and halo[sigma_key] is not None:
        return sigma_key
    if vmax_key in halo and halo[vmax_key] is not None:
        return f"{vmax_key}/sqrt2"
    return "virial/sqrt2"


def _quantile(values: list[float], q: float) -> float:
    if not values:
        raise ValueError("quantile requires at least one value")
    if q <= 0:
        return values[0]
    if q >= 1:
        return values[-1]
    pos = q * (len(values) - 1)
    lo = math.floor(pos)
    hi = math.ceil(pos)
    if lo == hi:
        return values[lo]
    frac = pos - lo
    return values[lo] * (1.0 - frac) + values[hi] * frac
