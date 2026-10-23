# -*- coding: utf-8 -*-
"""Finite Fourier matter trajectory diagnostic.

The stored clock is a.  rhs first forms the conformal-time equations
delta_eta = -theta-div(delta*v),
theta_eta = -Hcal*theta - laplacian(|v|^2)/2 - 3*Omega_m0*delta/(2*a),
then divides both by a*Hcal.  Here theta=div(v) and
v_hat=-i*k*theta_hat/k^2.  Angular wavenumbers, Hermitian initial data and
strict two-thirds dealiasing are used.  The initial amplitude convention
is diagnostic, not a calibrated realization of the production PM state.

Snapshot secants are diagnostics, not interval defect bounds.
Density positivity does not certify absence of shell crossing.
No imported viscous constant or power-decay ratio certifies this flow.
The earlier result files are historical; corrected outputs use a new name.
"""
import csv
import math
import os
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from cosmic_cg.params import load_model_params                    # noqa: E402
from cosmic_cg.generation import (                                # noqa: E402
    generate_epoch, InternalSpectrumConfig, internal_power_spectrum,
    growth_factor, omega_m0)
from cosmic_cg.knl import residual_floor_fraction                  # noqa: E402

RUN_DIR = str(Path(__file__).resolve().parents[1] / "results" / "pm_dynamics_64_protocol_001")
GRID = 128
BOX_MPCH = 128.0
SEED = 20260904
S = 2
K_RESCALE = BOX_MPCH / (2.0 * math.pi)
EPOCH_A = [0.009900990099009901, 0.01639344262295082, 0.024390243902439025,
           0.03225806451612903, 0.047619047619047616, 0.0625,
           0.07692307692307693, 0.09090909090909091, 0.1111111111111111,
           0.14285714285714285, 0.16666666666666666, 0.2, 0.25,
           0.3333333333333333, 0.4, 0.5, 0.5882352941176471,
           0.6666666666666666, 0.7692307692307693, 1.0, 1.25, 1.5, 2.0, 3.0]
A_REF = EPOCH_A[0]


def e_of_a(a, om_m0, om_l):
    return math.sqrt(om_m0 / a ** 3 + om_l)


def build_grid():
    kx = 2.0 * math.pi * np.fft.fftfreq(GRID, d=BOX_MPCH / GRID)
    ky = 2.0 * math.pi * np.fft.fftfreq(GRID, d=BOX_MPCH / GRID)
    kz = 2.0 * math.pi * np.fft.rfftfreq(GRID, d=BOX_MPCH / GRID)
    KX, KY, KZ = np.meshgrid(kx, ky, kz, indexing="ij")
    ksq = KX * KX + KY * KY + KZ * KZ
    inv = np.zeros_like(ksq)
    np.divide(1.0, ksq, out=inv, where=ksq > 0)
    return KX, KY, KZ, ksq, inv


def advective_terms(delta_hat, theta_hat, G):
    KX, KY, KZ, ksq, inv = G
    v0 = -1j * KX * inv * theta_hat
    v1 = -1j * KY * inv * theta_hat
    v2 = -1j * KZ * inv * theta_hat
    axes = (0, 1, 2)
    delta = np.fft.irfftn(delta_hat, s=(GRID,) * 3, axes=axes)
    vx = np.fft.irfftn(v0, s=(GRID,) * 3, axes=axes)
    vy = np.fft.irfftn(v1, s=(GRID,) * 3, axes=axes)
    vz = np.fft.irfftn(v2, s=(GRID,) * 3, axes=axes)
    div_delta_v = (np.fft.rfftn(delta * vx) * (1j * KX)
                   + np.fft.rfftn(delta * vy) * (1j * KY)
                   + np.fft.rfftn(delta * vz) * (1j * KZ))
    half_lap_v2 = -0.5 * ksq * np.fft.rfftn(vx * vx + vy * vy + vz * vz)
    return div_delta_v, half_lap_v2


def rhs(delta_hat, theta_hat, a, om_m0, om_l, G):
    cutoff = (2.0 * math.pi / BOX_MPCH) * GRID / 3.0
    mask = (np.abs(G[0]) < cutoff) & (np.abs(G[1]) < cutoff) & (np.abs(G[2]) < cutoff)
    delta_hat, theta_hat = delta_hat * mask, theta_hat * mask
    div_dv, half_lap_v2 = advective_terms(delta_hat, theta_hat, G)
    hh = a * e_of_a(a, om_m0, om_l)
    d_delta = -theta_hat - div_dv
    d_theta = -hh * theta_hat - half_lap_v2 - 1.5 * om_m0 * delta_hat / a
    # RK stages advance a, not conformal time.
    return mask * d_delta / (a * hh), mask * d_theta / (a * hh)


def rk4_step(delta_hat, theta_hat, a, da, om_m0, om_l, G):
    k1d, k1t = rhs(delta_hat, theta_hat, a, om_m0, om_l, G)
    k2d, k2t = rhs(delta_hat + 0.5 * da * k1d, theta_hat + 0.5 * da * k1t,
                   a + 0.5 * da, om_m0, om_l, G)
    k3d, k3t = rhs(delta_hat + 0.5 * da * k2d, theta_hat + 0.5 * da * k2t,
                   a + 0.5 * da, om_m0, om_l, G)
    k4d, k4t = rhs(delta_hat + da * k3d, theta_hat + da * k3t,
                   a + da, om_m0, om_l, G)
    return (delta_hat + (da / 6.0) * (k1d + 2 * k2d + 2 * k3d + k4d),
            theta_hat + (da / 6.0) * (k1t + 2 * k2t + 2 * k3t + k4t))


def bin_edges():
    return np.geomspace(0.049, 64.0 * math.pi / BOX_MPCH, 25)


def active_norm_s2(delta_hat, ksq, k_nl, edges):
    total = 0.0
    for k_lo, k_hi in zip(edges[:-1], edges[1:]):
        if k_hi <= k_nl:
            k_mid = math.sqrt(k_lo * k_hi)
            k_norm = k_mid * K_RESCALE
            sel = (ksq >= k_lo * k_lo) & (ksq < k_hi * k_hi)
            total += ((1.0 + k_norm * k_norm) ** S
                      * float(np.sum(np.abs(delta_hat[sel]) ** 2)))
    return total


def spectrum_on_grid(params, config, z, ksq):
    """Declared internal power spectrum interpolated onto the 3D grid."""
    from cosmic_cg.generation import logspace, _normalization_amplitude
    kg = logspace(config.k_min_hmpc, config.k_max_hmpc, config.k_points)
    pg = np.array(internal_power_spectrum(kg, z=z, params=params,
                                          config=config))
    km = np.sqrt(ksq).ravel()
    p = np.interp(km, kg, pg)
    return p.reshape(ksq.shape)


def main():
    steps = int(sys.argv[1]) if len(sys.argv) > 1 else 4
    params = load_model_params()
    config = InternalSpectrumConfig()
    om_m0 = float(omega_m0(params))
    om_l = float(params.omega_lambda)
    G = build_grid()
    KX, KY, KZ, ksq, inv = G
    edges = bin_edges()

    rng = np.random.default_rng(SEED)
    z0 = 1.0 / A_REF - 1.0
    p0 = spectrum_on_grid(params, config, z0, ksq)
    vol_factor = (2.0 * math.pi / BOX_MPCH) ** 3
    amp = np.sqrt(np.maximum(p0 * vol_factor, 0.0))
    g = np.fft.rfftn(rng.standard_normal((GRID,) * 3)) / math.sqrt(GRID**3)
    cutoff = (2.0 * math.pi / BOX_MPCH) * GRID / 3.0
    mask = (np.abs(KX) < cutoff) & (np.abs(KY) < cutoff) & (np.abs(KZ) < cutoff)
    delta_hat = amp * g * mask
    initial_delta = delta_hat.copy()
    delta_hat[ksq == 0] = 0.0
    D0 = growth_factor(z0, params=params)
    a0p = A_REF * 1.05
    Dp = growth_factor(1.0 / a0p - 1.0, params=params)
    f0 = (math.log(Dp) - math.log(D0)) / (math.log(a0p) - math.log(A_REF))
    theta_hat = -A_REF * e_of_a(A_REF, om_m0, om_l) * f0 * delta_hat

    rows = []
    a_cur = A_REF
    for idx, a in enumerate(EPOCH_A):
        if a > A_REF:
            da_total = a - a_cur
            for _ in range(steps):
                da = da_total / steps
                delta_hat, theta_hat = rk4_step(
                    delta_hat, theta_hat, a_cur, da, om_m0, om_l, G)
                a_cur += da
        z = 1.0 / a - 1.0
        D_ref = growth_factor(z, params=params)
        delta_ref = (D_ref / D0) * initial_delta
        delta_ref[ksq == 0] = 0.0
        epoch = generate_epoch(z=z, params=params, config=config)
        k_nl = float(epoch.knl.k_nl_hmpc)
        # single-source certificate: tail share computed from the evolved
        # field's own power above the endpoint, not from the declared spectrum
        ksq_np = ksq
        tail_sel = ksq_np > k_nl * k_nl
        tail_pow = float(np.sum(np.abs(delta_hat[tail_sel]) ** 2))
        tot_pow = float(np.sum(np.abs(delta_hat) ** 2))
        theta_rf = tail_pow / tot_pow if tot_pow > 0 else 0.0
        norm_s2 = active_norm_s2(delta_hat, ksq, k_nl, edges)
        resid = delta_hat - delta_ref
        resid_norm_s2 = active_norm_s2(resid, ksq, k_nl, edges)
        # energy audit (Parseval)
        ek = 0.5 * float(np.sum(np.abs(theta_hat) ** 2 * inv))
        ep = -0.75 * (om_m0 / a) * float(np.sum(np.abs(delta_hat) ** 2 * inv))
        axes = (0, 1, 2)
        delta_r = np.fft.irfftn(delta_hat, s=(GRID,) * 3, axes=axes)
        min_1d = float(np.min(1.0 + delta_r))
        # Endpoint secant diagnostic; it is not an interval defect bound.
        if a > A_REF and idx > 0:
            da_prev = (a - rows[-1]['a'])
            d_delta_rhs, d_theta_rhs = rhs(delta_hat, theta_hat, a,
                                            om_m0, om_l, G)
            resid_eq = (delta_hat - delta_prev) / da_prev - d_delta_rhs
            resid_eq_norm = float(np.sum(np.abs(resid_eq) ** 2 * inv))
        else:
            resid_eq_norm = 0.0
        delta_prev = delta_hat.copy()
        rows.append({"a": a, "z": z, "k_nl": k_nl, "norm_s2": norm_s2,
                     "resid_norm_s2": resid_norm_s2, "theta_rf": theta_rf,
                     "ek": ek, "ep": ep, "min_1d": min_1d,
                     "resid_eq": resid_eq_norm})

    ref = rows[0]["norm_s2"]
    for r in rows:
        r["norm_s"] = math.sqrt(r["norm_s2"] / ref)
        r["resid_share"] = math.sqrt(r["resid_norm_s2"] / ref)
        r["e_tot"] = r["ek"] + r["ep"]
        r["valid"] = r["min_1d"] > 0.0

    out = os.path.join(RUN_DIR, f"resolved_spectral_corrected_s{steps}_table.csv")
    with open(out, "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=[
            "a", "z", "k_nl", "norm_s", "resid_share", "min_1d",
            "density_positive", "secant_defect_squared", "status"])
        writer.writeheader()
        for row in rows:
            if not all(math.isfinite(row[key]) for key in ("norm_s", "resid_share", "resid_eq")):
                raise RuntimeError("Nonfinite trajectory; refine steps and restrict the interval")
            writer.writerow(dict(a=row["a"], z=row["z"], k_nl=row["k_nl"],
                                 norm_s=row["norm_s"], resid_share=row["resid_share"],
                                 min_1d=row["min_1d"], density_positive=row["valid"],
                                 secant_defect_squared=row["resid_eq"],
                                 status="trajectory_diagnostic_not_interval_certificate"))
    print("Corrected a-time trajectory diagnostic:", out)
    print("No C_flow certificate is inferred from snapshot secants or tail powers.")
    print("Use verify_finite_window_closure.py for the proved local finite-model certificate.")


if __name__ == "__main__":
    main()
