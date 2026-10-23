"""Construct dyadic Taylor records from densely encoded rational data."""

import argparse
from fractions import Fraction
import json
from math import comb
from pathlib import Path
import time
from record_format import (DEFAULT_INITIAL, EQUATION, format_limits,
                           guaranteed_parameters, input_from_coefficients,
                           load_input, normalize_input, require)

INITIAL = DEFAULT_INITIAL


def round_numerator(value, denominator):
    scaled = value.numerator * denominator
    magnitude = abs(scaled) // value.denominator
    return -magnitude if scaled < 0 else magnitude


def integer_product(left, right, cutoff):
    """Return the integer coefficients of twice the projected product."""
    result = [0] * (cutoff + 1)
    for i, x in enumerate(left):
        if not x:
            continue
        for j, y in enumerate(right):
            if not y:
                continue
            term = x*y
            result[abs(i-j)] += term
            if i+j <= cutoff:
                result[i+j] += term
    return result


def taylor_piece(node, cutoff, steps, degree, bits, guaranteed=False):
    g = (16 * steps).bit_length() - 1
    q = max(bits, g + 1)
    Q = 1 << q
    linear_factor = 1 << (2*q-g)
    product_factor = 1 << (q-g-1)
    coefficients = [[x << (q-bits) for x in node]]
    rounded = [list(node)]
    common_denominator = Q
    # The proved Cauchy bound makes the omitted exact Taylor coefficients
    # round to zero; no floating stopping test is used.
    active_degree = min(degree, (bits + 4)//3) if guaranteed else degree
    for ell in range(active_degree):
        nonlinear = [0] * (cutoff + 1)
        for a in range(ell//2 + 1):
            factor = comb(ell, a) * (1 if 2*a == ell else 2)
            term = integer_product(coefficients[a], coefficients[ell-a], cutoff)
            for j in range(cutoff + 1):
                nonlinear[j] += factor * term[j]
        row = [product_factor * nonlinear[j] - linear_factor*j*j*coefficients[ell][j]
               for j in range(cutoff + 1)]
        coefficients.append(row)
        common_denominator *= (ell + 1) * Q * Q
        rounded.append([(-1 if x < 0 else 1) * ((abs(x) << bits)//common_denominator)
                        for x in row])
    rounded.extend([[0]*(cutoff+1) for _ in range(degree-active_degree)])
    return rounded


def make_record(m, cutoff=None, steps=4, degree=None, bits=64,
                initial=None, mode="practical", progress=None):
    external = input_from_coefficients(INITIAL) if initial is None else initial
    canonical, values, input_bytes = normalize_input(external)
    if mode == "guaranteed":
        p = guaranteed_parameters(m)
        cutoff, steps, degree, bits = (p[k] for k in
            ("cutoff", "steps", "degree", "fractional_bits"))
    else:
        require(mode == "practical", "unknown parameter mode")
        cutoff = m+3 if cutoff is None else cutoff
        degree = m+4 if degree is None else degree
    limits = format_limits(m, input_bytes)
    for key, value in (("cutoff", cutoff), ("steps", steps), ("degree", degree),
                       ("fractional_bits", bits)):
        require(type(value) is int and 1 <= value <= limits[key], "parameter outside format: " + key)
    require(steps & (steps-1) == 0, "step count must be a power of two")
    denominator = 1 << bits
    projected = (values + [Fraction(0)] * (cutoff+1))[:cutoff+1]
    node = [round_numerator(x, denominator) for x in projected]
    pieces = []
    for i in range(steps):
        if mode == "guaranteed":
            require(sum((1 << (2*j))*abs(x) for j,x in enumerate(node))*4 <= 9*denominator,
                    "computed node exceeds the proved Taylor disk")
        rounded = taylor_piece(node, cutoff, steps, degree, bits, mode == "guaranteed")
        pieces.append([[str(x) for x in row] for row in rounded])
        node = [sum(row[j] for row in rounded) for j in range(cutoff+1)]
        if progress is not None and (i+1) % max(1, steps//16) == 0:
            progress(i+1, steps)
    return {"format": "neumann-dyadic-v2", "equation_id": EQUATION,
            "T": "1/16", "input_reference": canonical["initial"],
            "m": m, "mode": mode, "cutoff": cutoff, "steps": steps,
            "degree": degree, "fractional_bits": bits, "pieces": pieces,
            "observables": ["mean", "value_at_zero"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--m", type=int, required=True)
    parser.add_argument("--mode", choices=("practical", "guaranteed"), default="practical")
    parser.add_argument("--cutoff", type=int)
    parser.add_argument("--steps", type=int, default=4)
    parser.add_argument("--degree", type=int)
    parser.add_argument("--bits", type=int, default=64)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    started = time.perf_counter()
    record = make_record(args.m, args.cutoff, args.steps, args.degree, args.bits,
                         load_input(args.input), args.mode)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(record, separators=(",", ":"))+"\n", encoding="utf-8")
    print(json.dumps({"record": str(args.output), "producer_seconds": time.perf_counter()-started,
                      "record_bytes": args.output.stat().st_size}))


if __name__ == "__main__":
    main()
