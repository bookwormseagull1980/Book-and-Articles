# -*- coding: utf-8 -*-
"""Read PM density-power diagnostics without claiming a flow certificate.

Binned density powers do not determine the joint density/velocity Jacobian
or the interpolation defect.  This retained entry point writes only its
actual power readings.  Appendix B's independent closed finite-model
certificate is evaluated by verify_finite_window_closure.py.
"""
import csv
import math
import os
import sys

RUN_DIR = os.path.join(
    r"D:\论文撰写\论文9-粗粒化", "results", "pm_dynamics_64_protocol_001")

S = 2                # spectral Sobolev index
BOX_MPCH = 128.0     # PM box size (config_resolved.yaml)
K_RESCALE = BOX_MPCH / (2.0 * math.pi)   # angular wavenumber to lattice index

sys.path.insert(0, r"D:\论文撰写\论文9-粗粒化")
from cosmic_cg.params import load_model_params          # noqa: E402
from cosmic_cg.generation import (                       # noqa: E402
    generate_epoch, InternalSpectrumConfig)
from cosmic_cg.knl import tail_power_fraction, residual_floor_fraction  # noqa: E402


def load_pm_power():
    """Per-epoch list of (a, [(k_mid, power), ...]) from the PM power table."""
    rows = []
    path = os.path.join(RUN_DIR, "tables", "pm_power_table.csv")
    with open(path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh):
            rows.append(row)
    epochs = {}
    for row in rows:
        a = float(row["a"])
        epochs.setdefault(a, []).append(
            (float(row["k_mid_hmpc"]), float(row["k_max_hmpc"]),
             float(row["power_proxy"])))
    return epochs


def active_norm_s2(bins, k_nl):
    """||X_A||_s^2 over bins with k_max <= k_NL, s-weight <k_norm>^{2s}."""
    total = 0.0
    for k_mid, k_max, power in bins:
        if k_max <= k_nl:
            k_norm = k_mid * K_RESCALE
            japanese = 1.0 + k_norm * k_norm
            total += (japanese ** S) * power
    return total


def main():
    params = load_model_params()
    config = InternalSpectrumConfig()
    powers = load_pm_power()
    out = os.path.join(RUN_DIR, "pm_power_diagnostics.csv")
    with open(out, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["a", "k_nl", "active_weighted_power",
                         "sampled_tail_fraction", "status"])
        for a, bins in sorted(powers.items()):
            epoch = generate_epoch(z=1/a-1, params=params, config=config)
            endpoint = float(epoch.knl.k_nl_hmpc)
            total = sum(power for _, _, power in bins)
            tail = sum(power for mid, _, power in bins if mid >= endpoint)
            writer.writerow([a, endpoint, active_norm_s2(bins, endpoint),
                             tail/total if total > 0 else "",
                             "binned_density_diagnostic_not_flow_bound"])
    print("PM density-power diagnostics:", out)
    print("No propagation constant can be inferred from these columns alone.")


if __name__ == "__main__":
    main()
