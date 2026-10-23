"""Active-window cost scalings."""

from __future__ import annotations

import math


def active_mode_count(k_nl_hmpc: float, k_box_hmpc: float) -> float:
    if k_nl_hmpc <= 0 or k_box_hmpc <= 0:
        raise ValueError("wavenumbers must be positive")
    return (k_nl_hmpc / k_box_hmpc) ** 3


def active_spacetime_cost(k_nl_hmpc: float, k_box_hmpc: float, *, time_exponent: float = 2.0) -> float:
    if time_exponent < 0:
        raise ValueError("time_exponent must be non-negative")
    ratio = k_nl_hmpc / k_box_hmpc
    if ratio <= 0:
        raise ValueError("wavenumber ratio must be positive")
    return ratio ** (3.0 + time_exponent)

