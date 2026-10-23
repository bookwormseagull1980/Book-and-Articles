"""Physical constants and unit conversions used by the diagnostics."""

from __future__ import annotations

G_KPC_KM2_S2_MSUN = 4.30091727003628e-6
KPC_PER_KM_S_TO_GYR = 0.977792221673
KM_S_PER_MPC_TO_GYR_INV = 1.022712165045695e-3
HBAR_GEV_S = 6.582119569e-25
MPC_IN_KM = 3.0856775814913673e19
SECONDS_PER_GYR = 31557600.0e9


def h_scaled_mass_radius_to_physical_ratio(mass_msun_h: float, radius_kpc_h: float) -> float:
    """Return M/R in Msun/kpc for matching h-scaled mass and radius inputs."""

    if radius_kpc_h <= 0:
        raise ValueError("radius must be positive")
    return mass_msun_h / radius_kpc_h

