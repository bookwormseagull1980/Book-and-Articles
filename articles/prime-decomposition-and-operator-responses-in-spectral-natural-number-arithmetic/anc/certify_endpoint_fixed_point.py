"""Outward interval certificate for the diagonal endpoint fixed point.

Python standard library only. Every arithmetic endpoint is
rounded outwards. Pi, square root, and exponential use explicit enclosures.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import sys
from dataclasses import dataclass
from decimal import Context, Decimal, ROUND_CEILING, ROUND_FLOOR
from functools import lru_cache
from pathlib import Path


VERSION = "1.1.0"
PRECISION = 220
DOWN = Context(prec=PRECISION, rounding=ROUND_FLOOR)
UP = Context(prec=PRECISION, rounding=ROUND_CEILING)


def exact_decimal(value):
    """Accept exact decimal input, never an implicit binary-float conversion."""
    if type(value) not in (int, str, Decimal):
        raise TypeError("Use an integer, decimal string, or Decimal")
    result = Decimal(value)
    if not result.is_finite():
        raise ValueError("A finite decimal is required")
    return result


def configure_precision(digits):
    """Reset the serial calculation context and all precision-dependent caches."""
    global PRECISION, DOWN, UP
    if type(digits) is not int or digits < 220:
        raise ValueError("Certification requires at least 220 decimal digits")
    PRECISION = digits
    DOWN = Context(prec=digits, rounding=ROUND_FLOOR)
    UP = Context(prec=digits, rounding=ROUND_CEILING)
    exp_point.cache_clear()


@dataclass(frozen=True)
class Interval:
    lo: Decimal
    hi: Decimal

    def __post_init__(self):
        if not isinstance(self.lo, Decimal) or not isinstance(self.hi, Decimal):
            raise TypeError("Interval endpoints must be Decimal values")
        if not self.lo.is_finite() or not self.hi.is_finite() or self.lo > self.hi:
            raise ValueError("Invalid finite interval")

    @staticmethod
    def point(value):
        value = exact_decimal(value)
        return Interval(value, value)

    @staticmethod
    def coerce(value):
        return value if isinstance(value, Interval) else Interval.point(value)

    def __add__(self, other):
        other = self.coerce(other)
        return Interval(DOWN.add(self.lo, other.lo), UP.add(self.hi, other.hi))

    __radd__ = __add__

    def __neg__(self):
        return Interval(self.hi.copy_negate(), self.lo.copy_negate())

    def __sub__(self, other):
        return self + (-self.coerce(other))

    def __rsub__(self, other):
        return self.coerce(other) - self

    def __mul__(self, other):
        other = self.coerce(other)
        products = [(a, b) for a in (self.lo, self.hi) for b in (other.lo, other.hi)]
        return Interval(min(DOWN.multiply(a, b) for a, b in products),
                        max(UP.multiply(a, b) for a, b in products))

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = self.coerce(other)
        if other.lo <= 0 <= other.hi:
            raise ZeroDivisionError("Divisor interval contains zero")
        reciprocal = Interval(DOWN.divide(Decimal(1), other.hi),
                              UP.divide(Decimal(1), other.lo))
        return self * reciprocal

    def __rtruediv__(self, other):
        return self.coerce(other) / self

    def __pow__(self, exponent):
        if not isinstance(exponent, int) or exponent < 0:
            raise ValueError("Use a nonnegative integer power")
        result = Interval.point(1)
        base = self
        while exponent:
            if exponent % 2:
                result = result * base
            base = base * base
            exponent //= 2
        return result

    def abs_upper(self):
        return max(self.lo.copy_abs(), self.hi.copy_abs())

    def json(self):
        return {"lower": str(self.lo), "upper": str(self.hi)}


def sqrt_point(value):
    value = exact_decimal(value)
    if value < 0:
        raise ValueError("Negative radicand")
    if value == 0:
        return Interval.point(0)
    numerator, denominator = value.as_integer_ratio()
    places = PRECISION - 1 - value.adjusted() // 2
    if places >= 0:
        numerator *= 10 ** (2 * places)
    else:
        denominator *= 10 ** (-2 * places)
    root = math.isqrt(numerator // denominator)

    def scaled(integer):
        return Decimal((0, tuple(map(int, str(integer))), -places))

    upper = root if root * root * denominator == numerator else root + 1
    return Interval(scaled(root), scaled(upper))


def sqrt_interval(value):
    value = Interval.coerce(value)
    return Interval(sqrt_point(value.lo).lo, sqrt_point(value.hi).hi)


def atan_inverse(denominator):
    if type(denominator) is not int or denominator < 2:
        raise ValueError("The inverse-arctangent denominator must be an integer >= 2")
    # The omitted alternating tail has the sign of the next term and is
    # bounded by that term. All partial-sum arithmetic is outward rounded.
    total = Interval.point(0)
    power = denominator
    cutoff = Decimal("1e-" + str(PRECISION + 15))
    k = 0
    while True:
        term = Interval.point(1) / ((2 * k + 1) * power)
        if term.hi < cutoff:
            tail = Interval(Decimal(0), term.hi)
            return total + (tail if k % 2 == 0 else -tail)
        total = total + (term if k % 2 == 0 else -term)
        power *= denominator * denominator
        k += 1


def pi_interval():
    return 16 * atan_inverse(5) - 4 * atan_inverse(239)


@lru_cache(maxsize=None, typed=True)
def exp_point(value):
    value = exact_decimal(value)
    if value == 0:
        return Interval.point(1)
    if value < 0:
        return 1 / exp_point(value.copy_negate())
    half_steps = max(0, value.adjusted() * 4)
    reduced = Interval.point(value) / (2 ** half_steps)
    while reduced.hi > Decimal("0.125"):
        half_steps += 1
        reduced = Interval.point(value) / (2 ** half_steps)
    total = term = Interval.point(1)
    cutoff = Decimal("1e-" + str(PRECISION + 15))
    n = 0
    while True:
        next_term = term * reduced / (n + 1)
        remainder = next_term / (1 - reduced / (n + 2))
        if remainder.hi < cutoff:
            total = total + Interval(Decimal(0), remainder.hi)
            break
        total = total + next_term
        term = next_term
        n += 1
    for _ in range(half_steps):
        total = total * total
    return total


def exp_interval(value):
    value = Interval.coerce(value)
    return Interval(exp_point(value.lo).lo, exp_point(value.hi).hi)


def weight_sum(k):
    return Interval.point(k * (k + 1) * (k - 2)) / 3


def weight(k, m):
    return (2 * m * (k - m) * (1 - Interval.point(2) / (k - 1))
            + Interval.point(k * (k * k - 1)) / (3 * (k - 1) ** 2))


def decomposition_values(n_max):
    divisors = [0] * (n_max + 1)
    for d in range(1, n_max + 1):
        for multiple in range(d, n_max + 1, d):
            divisors[multiple] += 1
    return [max(0, count - 2) for count in divisors]


def gaussian_tail(beta, n_max, moment=0, shift=0):
    # The factor 4 is proved only for W(k) and k*W(k), with k >= 4.
    if (type(n_max) is not int or n_max < 4
            or type(moment) is not int or moment not in (0, 1)):
        raise ValueError("The tail bound requires cutoff >= 4 and moment 0 or 1")
    beta = Interval.coerce(beta)
    if beta.lo <= 0:
        raise ValueError("The Gaussian parameter must be strictly positive")
    k = n_max + 1
    first = (k ** moment) * weight_sum(k) * exp_interval(-beta * (k * k - shift))
    ratio = 4 * exp_interval(-beta * (2 * n_max + 3))
    if ratio.hi >= 1:
        raise ArithmeticError("Tail ratio is not contractive")
    return first / (1 - ratio)


def response_majorant(beta, n_max):
    finite = sum((k * weight_sum(k) * exp_interval(-beta * (k * k - 9)) / 4
                  for k in range(4, n_max + 1)), Interval.point(0))
    tail = gaussian_tail(beta, n_max, moment=1, shift=9) / 4
    return finite + Interval(Decimal(0), tail.hi)


def make_retention(beta_uv, beta_ir, n_max, r_values):
    data = {}
    for k in range(3, n_max + 1):
        data[k] = [(weight(k, m), r_values[m] + r_values[k - m])
                   for m in range(1, k)]
    decays = {beta: {k: exp_interval(-beta * (k * k - 9))
                      for k in range(3, n_max + 1)}
              for beta in (beta_uv, beta_ir)}
    leading = exp_interval(-9 * (beta_ir - beta_uv))

    def response(beta, tau):
        result = Interval.point(0)
        for k in range(3, n_max + 1):
            inner = sum((a * exp_interval(-tau * h) for a, h in data[k]),
                        Interval.point(0))
            result = result + decays[beta][k] * inner
        return result

    def retention(tau):
        return leading * response(beta_ir, tau) / response(beta_uv, tau)

    return retention, leading


def certify(n_max):
    if type(n_max) is not int or n_max < 30:
        raise ValueError("The published certificate requires cutoff >= 30")
    pi = pi_interval()
    beta_uv = pi ** 2 / (2 * (sqrt_interval(2 * pi) - sqrt_interval(2)) ** 2)
    beta_alt = pi ** 2 / (4 * (sqrt_interval(pi) - 1) ** 2)
    if max(beta_uv.lo, beta_alt.lo) > min(beta_uv.hi, beta_alt.hi):
        raise ArithmeticError("Endpoint formulas do not overlap")
    beta_ir = pi * beta_uv
    b_uv = response_majorant(beta_uv, n_max)
    b_ir = response_majorant(beta_ir, n_max)
    retention, leading = make_retention(beta_uv, beta_ir, n_max,
                                       decomposition_values(n_max))
    contraction = leading * (b_uv + b_ir)
    if not 0 <= contraction.lo <= contraction.hi < 1:
        raise ArithmeticError("A contraction bound strictly below one is required")
    first = retention(Interval.point(0))
    second = retention(first)
    step_bound = (second - first).abs_upper()
    q_upper = Interval.point(contraction.hi)
    iteration_error = q_upper / (1 - q_upper) * Interval.point(step_bound)
    tail_uv = gaussian_tail(beta_uv, n_max)
    denominator_floor = 4 * exp_interval(-9 * beta_uv)
    map_error = 2 * tail_uv / denominator_floor
    truncation_error = map_error / (1 - q_upper)
    total_error = iteration_error + truncation_error
    full_root = second + Interval(total_error.hi.copy_negate(), total_error.hi)
    one_step_error = q_upper * leading + map_error
    one_step_root = first + Interval(one_step_error.hi.copy_negate(), one_step_error.hi)
    if max(full_root.lo, one_step_root.lo) > min(full_root.hi, one_step_root.hi):
        raise ArithmeticError("One-step and two-step certificates do not overlap")

    bracket = Interval(Decimal("2.429273782017158e-35"),
                       Decimal("2.429273782017160e-35"))
    rounded_value = Interval.point("2.429273782017159e-35")
    absolute_error = (full_root - rounded_value).abs_upper()
    claims = {
        "B_UV_lt_3.579e-12": b_uv.hi < Decimal("3.579e-12"),
        "B_IR_lt_4.279e-39": b_ir.hi < Decimal("4.279e-39"),
        "contraction_lt_8.70e-47": contraction.hi < Decimal("8.70e-47"),
        "iterate_difference_lt_4e-98": step_bound < Decimal("4e-98"),
        "iteration_error_lt_4e-144": iteration_error.hi < Decimal("4e-144"),
        "second_iterate_interval_width_lt_1e-230":
            UP.subtract(second.hi, second.lo) < Decimal("1e-230"),
        "tail_lt_1e-1721": tail_uv.hi < Decimal("1e-1721"),
        "truncation_error_lt_1e-1705": truncation_error.hi < Decimal("1e-1705"),
        "root_in_decimal_bracket": bracket.lo < full_root.lo and full_root.hi < bracket.hi,
        "absolute_error_lt_1e-50": absolute_error < Decimal("1e-50"),
    }
    intervals = {
        "pi": pi, "beta_uv": beta_uv, "beta_ir": beta_ir,
        "B_uv": b_uv, "B_ir": b_ir, "contraction_majorant": contraction,
        "leading_bound": leading, "kappa_1": first, "kappa_2": second,
        "iteration_error_bound": iteration_error,
        "gaussian_tail_bound": tail_uv, "map_error_bound": map_error,
        "truncation_error_bound": truncation_error,
        "fixed_point": full_root, "one_step_fixed_point": one_step_root,
    }
    return {
        "certificate_version": VERSION,
        "python_version": platform.python_version(),
        "decimal_precision": PRECISION,
        "truncation_N": n_max,
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "intervals": {name: value.json() for name, value in intervals.items()},
        "absolute_error_from_displayed_value_upper": str(absolute_error),
        "second_iterate_interval_width_upper": str(UP.subtract(second.hi, second.lo)),
        "claims": claims,
        "all_claims_verified": all(claims.values()),
    }


def read_certificate(path):
    """Validate provenance before recomputing the archived settings."""
    reference = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(reference, dict):
        raise ValueError("Certificate must be a JSON object")
    if reference.get("certificate_version") != VERSION:
        raise ValueError("Certificate version does not match this source")
    if reference.get("script_sha256") != hashlib.sha256(Path(__file__).read_bytes()).hexdigest():
        raise ValueError("Certificate source checksum does not match this file")
    if (type(reference.get("decimal_precision")) is not int
            or reference["decimal_precision"] < 220
            or type(reference.get("truncation_N")) is not int
            or reference["truncation_N"] < 30):
        raise ValueError("Invalid certificate precision or cutoff")
    if reference.get("all_claims_verified") is not True:
        raise ValueError("Archived certificate is not successful")
    if not isinstance(reference.get("python_version"), str):
        raise ValueError("Missing archived Python version")
    return reference


def verify_recomputed(reference, computed):
    """Compare every numerical field; the runtime version is provenance only."""
    if set(reference) != set(computed):
        raise ValueError("Certificate fields differ from the current schema")
    for key in computed:
        if key != "python_version" and reference[key] != computed[key]:
            raise ValueError("Recomputed certificate differs in field: " + key)
    if computed["all_claims_verified"] is not True:
        raise ArithmeticError("Recomputed certificate failed")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version", action="version", version=VERSION)
    parser.add_argument("--digits", type=int, help="decimal precision (default: 220)")
    parser.add_argument("--cutoff", type=int, help="last retained shell (default: 30)")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--output", type=Path, help="write a new JSON file; never overwrite")
    mode.add_argument("--verify", type=Path, help="recompute and verify an archived JSON file")
    args = parser.parse_args()
    if args.verify and (args.digits is not None or args.cutoff is not None):
        parser.error("--verify uses the archived precision and cutoff")
    try:
        if args.output and args.output.exists():
            raise FileExistsError("Output already exists; choose a new filename")
        reference = read_certificate(args.verify) if args.verify else None
        digits = reference["decimal_precision"] if reference else (
            args.digits if args.digits is not None else 220)
        cutoff = reference["truncation_N"] if reference else (
            args.cutoff if args.cutoff is not None else 30)
        configure_precision(digits)
        certificate = certify(cutoff)
        if reference is not None:
            verify_recomputed(reference, certificate)
        if args.output:
            # Exclusive creation also protects the archive against a concurrent write.
            with args.output.open("x", encoding="utf-8", newline="\n") as stream:
                stream.write(json.dumps(certificate, indent=2) + "\n")
    except (ValueError, TypeError, ArithmeticError, OSError) as error:
        print("ERROR: " + str(error), file=sys.stderr)
        return 2
    print(json.dumps({key: certificate[key] for key in
                      ("certificate_version", "python_version", "decimal_precision",
                       "truncation_N", "all_claims_verified", "claims")}, indent=2))
    print("fixed_point =", certificate["intervals"]["fixed_point"])
    print("absolute_error_upper =", certificate["absolute_error_from_displayed_value_upper"])
    if not certificate["all_claims_verified"]:
        return 1
    if reference is not None:
        print("VERIFIED: source checksum, numerical intervals, and all claims match")
        print("Python archived/current:", reference["python_version"], platform.python_version())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
