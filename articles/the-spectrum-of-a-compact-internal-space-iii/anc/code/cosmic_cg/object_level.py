"""Halo-level fate and tidal-spin diagnostics."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Literal, Sequence

from .units import G_KPC_KM2_S2_MSUN, KPC_PER_KM_S_TO_GYR


FateClass = Literal["forming", "balanced", "stripping"]


@dataclass(frozen=True)
class HaloFate:
    halo_id: str
    b_h: float
    t_dyn_gyr: float
    t_growth_gyr: float | None
    class_name: FateClass


def dynamical_time_gyr(mass_msun_h: float, radius_kpc_h: float) -> float:
    """Return sqrt(R^3/GM) in Gyr for matching h-scaled mass and radius."""

    if mass_msun_h <= 0:
        raise ValueError("mass must be positive")
    if radius_kpc_h <= 0:
        raise ValueError("radius must be positive")
    time_kpc_per_kms = math.sqrt(radius_kpc_h**3 / (G_KPC_KM2_S2_MSUN * mass_msun_h))
    return time_kpc_per_kms * KPC_PER_KM_S_TO_GYR


def growth_time_gyr(mass_msun_h: float, delta_mass_msun_h: float, delta_time_gyr: float) -> float | None:
    """Return M/dMdt in Gyr."""

    if mass_msun_h <= 0:
        raise ValueError("mass must be positive")
    if delta_time_gyr <= 0:
        raise ValueError("delta_time_gyr must be positive")
    if delta_mass_msun_h == 0:
        return None
    return mass_msun_h * delta_time_gyr / delta_mass_msun_h


def halo_fate(
    *,
    halo_id: str,
    mass_msun_h: float,
    radius_kpc_h: float,
    delta_mass_msun_h: float,
    delta_time_gyr: float,
    forming_threshold: float = 1.0,
    balanced_threshold: float = 0.0,
) -> HaloFate:
    """Classify a halo using B_h=t_dyn/t_growth."""

    t_dyn = dynamical_time_gyr(mass_msun_h, radius_kpc_h)
    t_growth = growth_time_gyr(mass_msun_h, delta_mass_msun_h, delta_time_gyr)
    b_h = math.nan if t_growth is None else t_dyn / t_growth
    if delta_mass_msun_h < 0:
        class_name: FateClass = "stripping"
    elif math.isnan(b_h):
        class_name = "balanced"
    elif b_h > forming_threshold:
        class_name: FateClass = "forming"
    elif b_h >= balanced_threshold:
        class_name = "balanced"
    else:
        class_name = "stripping"
    return HaloFate(halo_id=halo_id, b_h=b_h, t_dyn_gyr=t_dyn, t_growth_gyr=t_growth, class_name=class_name)


def tidal_torque(inertia: Sequence[Sequence[float]], tidal: Sequence[Sequence[float]]) -> tuple[float, float, float]:
    """Return tau_i=epsilon_ijk I_jl T_lk for 3x3 inertia and tidal tensors."""

    _check_3x3(inertia, "inertia")
    _check_3x3(tidal, "tidal")
    tau = []
    for i in range(3):
        total = 0.0
        for j in range(3):
            for k in range(3):
                eps = _levi_civita(i, j, k)
                if eps == 0:
                    continue
                for ell in range(3):
                    total += eps * inertia[j][ell] * tidal[ell][k]
        tau.append(total)
    return (tau[0], tau[1], tau[2])


def projected_torque(torque: Sequence[float], spin: Sequence[float]) -> float:
    """Project a torque vector along the halo spin direction."""

    if len(torque) != 3 or len(spin) != 3:
        raise ValueError("torque and spin must be three-vectors")
    spin_norm = math.sqrt(sum(float(x) * float(x) for x in spin))
    if spin_norm <= 0:
        raise ValueError("spin vector norm must be positive")
    return sum(float(torque[i]) * float(spin[i]) for i in range(3)) / spin_norm


def alignment_cosine(a: Sequence[float], b: Sequence[float]) -> float:
    """Return the cosine of the angle between two three-vectors."""

    if len(a) != 3 or len(b) != 3:
        raise ValueError("inputs must be three-vectors")
    norm_a = math.sqrt(sum(float(x) * float(x) for x in a))
    norm_b = math.sqrt(sum(float(x) * float(x) for x in b))
    if norm_a <= 0 or norm_b <= 0:
        raise ValueError("vector norms must be positive")
    return sum(float(a[i]) * float(b[i]) for i in range(3)) / (norm_a * norm_b)


def _check_3x3(matrix: Sequence[Sequence[float]], name: str) -> None:
    if len(matrix) != 3 or any(len(row) != 3 for row in matrix):
        raise ValueError(f"{name} must be a 3x3 matrix")


def _levi_civita(i: int, j: int, k: int) -> int:
    if len({i, j, k}) < 3:
        return 0
    if (i, j, k) in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
        return 1
    return -1
