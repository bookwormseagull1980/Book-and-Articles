"""Exact-arithmetic certificate for Appendix B's finite spectral advance."""

import json
from fractions import Fraction as Q
from pathlib import Path


def certificate(step=Q(1, 10000), duration=Q(1, 100)):
    if step <= 0 or duration <= 0:
        raise ValueError("Step and duration must be positive")
    count = duration / step
    if count.denominator != 1:
        raise ValueError("Positive duration must contain an integer number of steps")
    radius = Q(3, 100)
    linear, bilinear = Q(7, 2), Q(3)
    lip = linear + 2 * bilinear * radius
    speed = linear * radius + bilinear * radius**2
    if lip * duration >= 1:
        raise ValueError("Rational exponential bound not available")
    amplitude = Q(1, 100)
    trajectory = [amplitude]
    for _ in range(int(count)):
        amplitude += step * amplitude
        trajectory.append(amplitude)
    initial_error = Q(1, 75000)
    compression_defect = bilinear * radius**2 * duration
    integration_defect = lip * speed * step * duration / 2
    amplification = 1 / (1 - lip * duration)
    error = amplification * (initial_error + compression_defect + integration_defect)
    interpolant_norm = 2 * max(trajectory)
    tolerance = Q(1, 10000)
    return {
        "steps": int(count),
        "initial_error": initial_error,
        "compression_defect_integral": compression_defect,
        "integration_defect_integral": integration_defect,
        "same_norm_lipschitz": lip,
        "amplification_upper": amplification,
        "uniform_error_upper": error,
        "interpolant_norm_upper": interpolant_norm,
        "density_lower": 1 - radius,
        "tube_pass": interpolant_norm + error < radius,
        "half_tolerance_pass": error < tolerance / 2,
        "reference": "three-mode Galerkin matter system",
        "retained": "one-mode system, exact rational Euler updates",
        "scope": "local finite-model certificate, not the production PM history",
    }


def main():
    result = certificate()
    if not result["tube_pass"] or not result["half_tolerance_pass"]:
        raise RuntimeError("Certificate did not close")
    path = Path(__file__).resolve().parents[1] / "tables" / "finite_window_certificate.json"
    path.write_text(json.dumps(result, default=str, indent=2) + "\n", encoding="utf-8")
    print("Exact error bound:", result["uniform_error_upper"])
    print("Decimal upper bound:", float(result["uniform_error_upper"]))
    print("Tube and half-tolerance: PASS")
    print(path)


if __name__ == "__main__":
    main()
