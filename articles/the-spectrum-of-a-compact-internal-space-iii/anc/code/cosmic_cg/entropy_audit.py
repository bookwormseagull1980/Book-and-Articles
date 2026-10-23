"""Window entropy and truncation-deficit diagnostics."""

from __future__ import annotations

from dataclasses import dataclass
import math


@dataclass(frozen=True)
class EntropyAudit:
    s_abs: float
    s_aw: float
    d_s: float
    k_box_hmpc: float
    k_nl_hmpc: float
    k_cut_hmpc: float | None


def window_entropy(k_nl_hmpc: float, k_box_hmpc: float) -> float:
    """Return S_AW=log(k_NL/k_box)."""

    if k_nl_hmpc <= 0 or k_box_hmpc <= 0:
        raise ValueError("wavenumbers must be positive")
    if k_nl_hmpc <= k_box_hmpc:
        raise ValueError("k_nl_hmpc must exceed k_box_hmpc")
    return math.log(k_nl_hmpc / k_box_hmpc)


def entropy_deficit(k_nl_hmpc: float, k_cut_hmpc: float) -> float:
    """Return the unresolved span max(0, log(k_NL/k_cut))."""

    if k_nl_hmpc <= 0 or k_cut_hmpc <= 0:
        raise ValueError("wavenumbers must be positive")
    return max(0.0, math.log(k_nl_hmpc / k_cut_hmpc))


def audit_entropy(
    *,
    s_abs: float,
    k_nl_hmpc: float,
    k_box_hmpc: float,
    k_cut_hmpc: float | None = None,
) -> EntropyAudit:
    d_s = 0.0 if k_cut_hmpc is None else entropy_deficit(k_nl_hmpc, k_cut_hmpc)
    return EntropyAudit(
        s_abs=s_abs,
        s_aw=window_entropy(k_nl_hmpc, k_box_hmpc),
        d_s=d_s,
        k_box_hmpc=k_box_hmpc,
        k_nl_hmpc=k_nl_hmpc,
        k_cut_hmpc=k_cut_hmpc,
    )
