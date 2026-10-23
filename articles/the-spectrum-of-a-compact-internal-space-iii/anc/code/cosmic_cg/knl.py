"""Nonlinear endpoint and tail-power calculations."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Iterable, Sequence


@dataclass(frozen=True)
class KNLResult:
    z: float
    c_nl: float
    k_nl_hmpc: float
    crossing_index: int
    method: str = "log_linear_interpolation"


def dimensionless_power(k_hmpc: Sequence[float], p_mpc3_h3: Sequence[float]) -> list[float]:
    """Return Delta^2(k)=k^3 P(k)/(2*pi^2)."""

    _check_pair(k_hmpc, p_mpc3_h3)
    return [(k**3) * p / (2.0 * math.pi**2) for k, p in zip(k_hmpc, p_mpc3_h3)]


def find_knl(
    k_hmpc: Sequence[float],
    p_mpc3_h3: Sequence[float],
    *,
    z: float = 0.0,
    c_nl: float = 1.0,
) -> KNLResult:
    """Find the first crossing Delta^2(k,z)=C_NL in ascending k."""

    if c_nl <= 0:
        raise ValueError("c_nl must be positive")
    delta2 = dimensionless_power(k_hmpc, p_mpc3_h3)
    for i in range(1, len(k_hmpc)):
        y0 = delta2[i - 1] - c_nl
        y1 = delta2[i] - c_nl
        if y0 == 0:
            return KNLResult(z=z, c_nl=c_nl, k_nl_hmpc=k_hmpc[i - 1], crossing_index=i - 1)
        if y0 * y1 <= 0:
            k_cross = _log_linear_crossing(k_hmpc[i - 1], delta2[i - 1], k_hmpc[i], delta2[i], c_nl)
            return KNLResult(z=z, c_nl=c_nl, k_nl_hmpc=k_cross, crossing_index=i - 1)
    raise ValueError("no Delta^2 crossing found in supplied k range")


def tail_power_fraction(
    k_hmpc: Sequence[float],
    p_mpc3_h3: Sequence[float],
    k_nl_hmpc: float,
) -> float:
    """Return the fraction of integrated Delta^2 dlnk above k_NL."""

    _check_pair(k_hmpc, p_mpc3_h3)
    if k_nl_hmpc <= 0:
        raise ValueError("k_nl_hmpc must be positive")
    delta2 = dimensionless_power(k_hmpc, p_mpc3_h3)
    total = _integrate_log_k(k_hmpc, delta2, k_hmpc[0], k_hmpc[-1])
    if total <= 0:
        raise ValueError("total integrated power must be positive")
    tail = _integrate_log_k(k_hmpc, delta2, k_nl_hmpc, k_hmpc[-1])
    return tail / total


def tail_power_moment_fraction(
    k_hmpc: Sequence[float],
    p_mpc3_h3: Sequence[float],
    k_nl_hmpc: float,
    *,
    s: float = 1.0,
) -> float:
    """Return an endpoint-weighted tail fraction."""

    _check_pair(k_hmpc, p_mpc3_h3)
    if k_nl_hmpc <= 0:
        raise ValueError("k_nl_hmpc must be positive")
    if s < 0:
        raise ValueError("s must be non-negative")
    delta2 = dimensionless_power(k_hmpc, p_mpc3_h3)
    total = _integrate_log_k(k_hmpc, delta2, k_hmpc[0], k_hmpc[-1])
    if total <= 0:
        raise ValueError("total integrated power must be positive")
    points = _sample_with_boundaries(k_hmpc, delta2, k_nl_hmpc, k_hmpc[-1])
    if len(points) < 2:
        return 0.0
    spacing = _median_log_spacing(k_hmpc)
    moment = 0.0
    for (k0, y0), (k1, y1) in zip(points, points[1:]):
        midpoint = math.sqrt(k0 * k1)
        shell_distance = max(0.0, math.log(midpoint / k_nl_hmpc) / spacing)
        weight = math.exp(-s * shell_distance)
        moment += 0.5 * (y0 + y1) * (math.log(k1) - math.log(k0)) * weight
    return moment / total


def residual_floor_fraction(
    k_hmpc: Sequence[float],
    p_mpc3_h3: Sequence[float],
    k_nl_hmpc: float,
    *,
    q_floor: float = 2.0,
    sigma_g: float = 0.5,
) -> float:
    """Return the Gaussian floor fraction retained by the tail audit."""

    _check_pair(k_hmpc, p_mpc3_h3)
    if k_nl_hmpc <= 0:
        raise ValueError("k_nl_hmpc must be positive")
    if q_floor <= 0 or sigma_g <= 0:
        raise ValueError("floor parameters must be positive")
    delta2 = dimensionless_power(k_hmpc, p_mpc3_h3)
    total = _integrate_log_k(k_hmpc, delta2, k_hmpc[0], k_hmpc[-1])
    if total <= 0:
        raise ValueError("total integrated power must be positive")
    points = _sample_with_boundaries(k_hmpc, delta2, k_nl_hmpc, k_hmpc[-1])
    if len(points) < 2:
        return 0.0
    floor = 0.0
    for (k0, y0), (k1, y1) in zip(points, points[1:]):
        midpoint = math.sqrt(k0 * k1)
        q_value = midpoint / k_nl_hmpc
        weight = 0.5 * math.erfc(math.log(q_value / q_floor) / (math.sqrt(2.0) * sigma_g))
        floor += 0.5 * (y0 + y1) * (math.log(k1) - math.log(k0)) * weight
    return max(0.0, min(floor / total, 1.0))


def band_power_fraction(
    k_hmpc: Sequence[float],
    p_mpc3_h3: Sequence[float],
    k_min_hmpc: float,
    k_max_hmpc: float,
) -> float:
    """Return the fraction of integrated Delta^2 dlnk in a finite band."""

    _check_pair(k_hmpc, p_mpc3_h3)
    if k_min_hmpc <= 0 or k_max_hmpc <= 0:
        raise ValueError("band endpoints must be positive")
    if k_min_hmpc >= k_max_hmpc:
        raise ValueError("k_min_hmpc must be smaller than k_max_hmpc")
    delta2 = dimensionless_power(k_hmpc, p_mpc3_h3)
    total = _integrate_log_k(k_hmpc, delta2, k_hmpc[0], k_hmpc[-1])
    if total <= 0:
        raise ValueError("total integrated power must be positive")
    lower = max(k_min_hmpc, k_hmpc[0])
    upper = min(k_max_hmpc, k_hmpc[-1])
    if lower >= upper:
        return 0.0
    return _integrate_log_k(k_hmpc, delta2, lower, upper) / total


def _log_linear_crossing(k0: float, y0: float, k1: float, y1: float, target: float) -> float:
    if k0 <= 0 or k1 <= 0 or y0 <= 0 or y1 <= 0:
        frac = (target - y0) / (y1 - y0)
        return k0 + frac * (k1 - k0)
    x0 = math.log(k0)
    x1 = math.log(k1)
    ly0 = math.log(y0)
    ly1 = math.log(y1)
    frac = (math.log(target) - ly0) / (ly1 - ly0)
    return math.exp(x0 + frac * (x1 - x0))


def _integrate_log_k(k: Sequence[float], y: Sequence[float], lower: float, upper: float) -> float:
    if lower >= upper:
        return 0.0
    points = _sample_with_boundaries(k, y, lower, upper)
    area = 0.0
    for (k0, y0), (k1, y1) in zip(points, points[1:]):
        area += 0.5 * (y0 + y1) * (math.log(k1) - math.log(k0))
    return area


def _sample_with_boundaries(
    k: Sequence[float],
    y: Sequence[float],
    lower: float,
    upper: float,
) -> list[tuple[float, float]]:
    lower = max(lower, k[0])
    upper = min(upper, k[-1])
    points: list[tuple[float, float]] = [(lower, _interp_log_linear(k, y, lower))]
    for kk, yy in zip(k, y):
        if lower < kk < upper:
            points.append((kk, yy))
    points.append((upper, _interp_log_linear(k, y, upper)))
    return points


def _interp_log_linear(k: Sequence[float], y: Sequence[float], x: float) -> float:
    if x <= k[0]:
        return y[0]
    if x >= k[-1]:
        return y[-1]
    for i in range(1, len(k)):
        if k[i] >= x:
            if y[i - 1] <= 0 or y[i] <= 0:
                frac = (x - k[i - 1]) / (k[i] - k[i - 1])
                return y[i - 1] + frac * (y[i] - y[i - 1])
            lx0, lx1 = math.log(k[i - 1]), math.log(k[i])
            ly0, ly1 = math.log(y[i - 1]), math.log(y[i])
            frac = (math.log(x) - lx0) / (lx1 - lx0)
            return math.exp(ly0 + frac * (ly1 - ly0))
    return y[-1]


def _median_log_spacing(k: Sequence[float]) -> float:
    values = sorted(math.log(k1 / k0) for k0, k1 in zip(k, k[1:]) if k1 > k0)
    if not values:
        raise ValueError("at least two increasing k samples are required")
    mid = len(values) // 2
    if len(values) % 2:
        return values[mid]
    return 0.5 * (values[mid - 1] + values[mid])


def _check_pair(k: Sequence[float], p: Sequence[float]) -> None:
    if len(k) != len(p):
        raise ValueError("k and p arrays must have the same length")
    if len(k) < 2:
        raise ValueError("at least two samples are required")
    previous = -math.inf
    for value in k:
        if value <= 0:
            raise ValueError("k values must be positive")
        if value <= previous:
            raise ValueError("k values must be strictly increasing")
        previous = value
    if any(value < 0 for value in p):
        raise ValueError("power values must be non-negative")
