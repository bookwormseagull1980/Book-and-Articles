"""Cosmological clock for a flat pressureless matter--vacuum realization.

With u=dx/da, p=a^3 E(a)u eliminates the exact Hubble damping.
The force g solves div(g)=-delta; these coefficients do not derive that
Poisson closure or a microscopic cold-source distribution.
"""
from functools import lru_cache
import math
import numpy as np

from cosmic_cg.params import validate_model_params


def expansion(a, params):
    if not math.isfinite(a) or a <= 0:
        raise ValueError("Scale factor must be finite and positive")
    return math.sqrt((params.omega_b + params.omega_sigma) / a**3 + params.omega_lambda)


def coefficients(a, params):
    e = expansion(a, params)
    omega = (params.omega_b + params.omega_sigma) / (a**3 * e**2)
    return (3.0 - 1.5 * omega) / a, 1.5 * omega / a**2


@lru_cache(maxsize=1)
def _quadrature():
    return np.polynomial.legendre.leggauss(8)


def drift_kick_integrals(a0, a1, params):
    """Eighth-order Gaussian quadrature of the two physical coefficients."""
    if not a1 >= a0 > 0:
        raise ValueError("Require 0 < a0 <= a1")
    nodes, weights = _quadrature()
    a = (a0 + a1) / 2 + (a1 - a0) * nodes / 2
    omega0 = params.omega_b + params.omega_sigma
    e = np.sqrt(omega0 / a**3 + params.omega_lambda)
    factor = (a1 - a0) / 2
    drift = factor * np.dot(weights, 1 / (a**3 * e))
    kick = factor * np.dot(weights, 1.5 * omega0 / (a**2 * e))
    return float(drift), float(kick)


def growing_rate(a, params):
    """f=d log(D)/d log(a), from the same matter--vacuum growth solution."""
    from scipy.integrate import quad
    validate_model_params(params)
    expansion(a, params)
    omega0 = params.omega_b + params.omega_sigma
    integral = quad(lambda x: x**1.5 / (omega0 + params.omega_lambda * x**3)**1.5,
                    0.0, a, epsabs=1e-13, epsrel=1e-12)[0]
    if not math.isfinite(integral) or integral <= 0:
        raise ValueError("The matter growing-mode integral must be finite and positive")
    e = expansion(a, params)
    return -1.5 * omega0 / (a**3 * e**2) + 1 / (a**2 * e**3 * integral)


def linear_growth_advance(a0, a1, steps, params, *, displacement=1.0, velocity=None):
    """Independent one-mode clock check: its exact linear force is g=D."""
    validate_model_params(params)
    if not isinstance(steps, int) or steps <= 0 or not a1 > a0 > 0:
        raise ValueError("Invalid integration interval or step count")
    if velocity is None:
        velocity = growing_rate(a0, params) * displacement / a0
    p = a0**3 * expansion(a0, params) * velocity
    d = float(displacement)
    step = (a1-a0)/steps
    for n in range(steps):
        a, end = a0+n*step, a0+(n+1)*step
        middle = (a+end)/2
        _, k1 = drift_kick_integrals(a, middle, params)
        drift, _ = drift_kick_integrals(a, end, params)
        _, k2 = drift_kick_integrals(middle, end, params)
        p += k1*d
        d += drift*p
        p += k2*d
    return d, p/(a1**3 * expansion(a1, params))
