"""Minimal halo-model utilities used by the baseline diagnostics."""

from __future__ import annotations

import math

from .tail_ratio import virial_velocity


def virial_radius_from_mass(
    mass_msun_h: float,
    *,
    omega_m: float = 0.3149557,
    h: float = 0.6745,
    overdensity: float = 200.0,
) -> float:
    """Return R_200m in kpc/h for a mass in Msun/h."""

    if mass_msun_h <= 0:
        raise ValueError("mass must be positive")
    if omega_m <= 0 or h <= 0 or overdensity <= 0:
        raise ValueError("cosmological parameters must be positive")
    rho_crit_msun_h2_mpc3 = 2.77536627e11
    rho_m_msun_h2_mpc3 = omega_m * rho_crit_msun_h2_mpc3
    radius_mpc_h = (3.0 * mass_msun_h / (4.0 * math.pi * overdensity * rho_m_msun_h2_mpc3)) ** (1.0 / 3.0)
    return radius_mpc_h * 1000.0


def endpoint_virial_velocity(
    mass_msun_h: float,
    *,
    omega_m: float = 0.3149557,
    h: float = 0.6745,
    overdensity: float = 200.0,
) -> tuple[float, float]:
    """Return radius in kpc/h and virial velocity in km/s for M_NL."""

    radius = virial_radius_from_mass(mass_msun_h, omega_m=omega_m, h=h, overdensity=overdensity)
    return radius, virial_velocity(mass_msun_h, radius)

