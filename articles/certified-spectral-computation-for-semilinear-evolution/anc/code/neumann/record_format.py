"""Input normalization, integer parameters, and bounded record parsing.

This module contains no trajectory construction or acceptance estimate.
"""

from fractions import Fraction
import json
from pathlib import Path
import re
import sys

if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)

EQUATION = "u_t=u_xx+u^2; Neumann [0,pi]"
T = Fraction(1, 16)
DEFAULT_INITIAL = [Fraction(1, 4), Fraction(1, 16), Fraction(1, 64)]
INTEGER = re.compile(r"-?(0|[1-9][0-9]*)\Z")
RATIONAL = re.compile(r"-?(0|[1-9][0-9]*)(/[1-9][0-9]*)?\Z")


class Rejected(ValueError):
    pass


def require(condition, message):
    if not condition:
        raise Rejected(message)


def ceil_log2_rational(value):
    """Smallest integer k for which 2**k >= value, using integers only."""
    value = Fraction(value)
    require(value > 0, "logarithm requires a positive rational")
    k = value.numerator.bit_length() - value.denominator.bit_length()
    if k >= 0:
        if (value.denominator << k) < value.numerator:
            k += 1
    elif value.denominator < (value.numerator << -k):
        k += 1
    return k


def guaranteed_parameters(m):
    require(type(m) is int and m >= 1, "precision must be a positive integer")
    cutoff = m + 3
    dimension = cutoff + 1
    lipschitz = 16 * dimension * dimension
    g = (32 * lipschitz - 1).bit_length()
    steps = 1 << (g - 4)
    amplification = 1 << (2 * dimension * dimension)
    xi = Fraction(1, (1 << m) * 8 * amplification) / (1 + T)
    degree = ceil_log2_rational(64 * (1 + (1 << g)) / xi)
    weights = ((1 << (2 * (cutoff + 1))) - 1) // 3
    bits = ceil_log2_rational(64 * (degree + 2)**2 * (1 + (1 << g)) * weights / xi)
    return {"cutoff": cutoff, "steps": steps, "degree": degree,
            "fractional_bits": bits, "g": g, "xi": xi}


def format_limits(m, input_bytes=0):
    p = guaranteed_parameters(m)
    limits = {"cutoff": 2 * (m + 4), "steps": p["steps"],
              "degree": max(m + 8, p["degree"]),
              "fractional_bits": max(64, p["fractional_bits"])}
    limits["integer_characters"] = limits["fractional_bits"] + 3
    limits["record_bytes"] = (8192 + 2 * input_bytes +
        limits["steps"] * (limits["degree"] + 1) *
        (limits["cutoff"] + 1) * (limits["integer_characters"] + 4))
    return limits


def normalize_input(data):
    require(isinstance(data, dict), "input must be an object")
    require(data.get("format") == "neumann-input-v1", "unknown input format")
    require(data.get("equation_id") == EQUATION, "wrong input equation")
    require(data.get("T") == "1/16", "wrong input time interval")
    raw = data.get("initial")
    require(isinstance(raw, list) and len(raw) >= 1, "initial coefficients must be dense")
    require(all(isinstance(x, str) and RATIONAL.fullmatch(x) for x in raw),
            "initial coefficients must be rational strings")
    values = [Fraction(x) for x in raw]
    require(sum((1 << (2*j)) * abs(x) for j, x in enumerate(values)) <= 1,
            "initial A4 budget fails")
    canonical = {"format": "neumann-input-v1", "equation_id": EQUATION,
                 "T": "1/16", "initial": [str(x) for x in values]}
    encoded_bytes = len(json.dumps(canonical, separators=(",", ":")).encode("utf-8"))
    return canonical, values, encoded_bytes


def input_from_coefficients(values):
    return normalize_input({"format": "neumann-input-v1", "equation_id": EQUATION,
                            "T": "1/16", "initial": [str(x) for x in values]})[0]


def load_input(path):
    return normalize_input(json.loads(Path(path).read_text(encoding="utf-8")))[0]


def load_record(path, expected_input, expected_m):
    _, _, input_bytes = normalize_input(expected_input)
    limits = format_limits(expected_m, input_bytes)
    path = Path(path)
    require(path.stat().st_size <= limits["record_bytes"], "record byte bound exceeded")
    largest = max(expected_m, *(limits[k] for k in
                  ("cutoff", "steps", "degree", "fractional_bits")))
    header_digits = len(str(largest))

    def bounded_integer(text):
        require(len(text.lstrip("-")) <= header_digits, "oversized JSON integer")
        return int(text)

    with path.open("r", encoding="utf-8") as stream:
        return json.load(stream, parse_int=bounded_integer,
                         parse_float=lambda _: (_ for _ in ()).throw(Rejected("floating JSON field")))
