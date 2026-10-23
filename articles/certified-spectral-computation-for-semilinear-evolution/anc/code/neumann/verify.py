"""Independently certify a polynomial trajectory against external input.

The verifier imports no constructor or Taylor recurrence.  It checks finite
arrays and expands their differential residual by exact integer arithmetic.
"""

import argparse
from fractions import Fraction
import json
from pathlib import Path
import time
from record_format import (DEFAULT_INITIAL, EQUATION, INTEGER, Rejected, T,
                           format_limits, guaranteed_parameters,
                           input_from_coefficients, load_input, load_record,
                           normalize_input, require)


def fraction_report(value):
    return {"exact": str(value), "decimal": float(value)}


def verify(record, expected_m, expected_input=None):
    external = input_from_coefficients(DEFAULT_INITIAL) if expected_input is None else expected_input
    canonical, initial, input_bytes = normalize_input(external)
    limits = format_limits(expected_m, input_bytes)
    require(isinstance(record, dict), "record must be an object")
    require(record.get("format") == "neumann-dyadic-v2", "unknown record format")
    require(record.get("equation_id") == EQUATION, "wrong equation")
    require(record.get("T") == "1/16", "wrong time interval")
    require(record.get("input_reference") == canonical["initial"], "initial input mismatch")
    require(type(record.get("m")) is int and record["m"] == expected_m,
            "precision request mismatch")
    require(record.get("observables") == ["mean", "value_at_zero"], "unknown observations")
    mode = record.get("mode")
    require(mode in ("practical", "guaranteed"), "unknown parameter mode")
    parameters = {}
    for key in ("cutoff", "steps", "degree", "fractional_bits"):
        value = record.get(key)
        require(type(value) is int and 1 <= value <= limits[key], "parameter bound: " + key)
        parameters[key] = value
    cutoff, steps, degree, bits = (parameters[k] for k in
        ("cutoff", "steps", "degree", "fractional_bits"))
    require(steps & (steps-1) == 0, "step count must be a power of two")
    if mode == "guaranteed":
        p = guaranteed_parameters(expected_m)
        require(all(parameters[k] == p[k] for k in parameters), "theorem parameter mismatch")

    # Validate dimensions and string lengths before coefficient integers,
    # dyadic denominators, or residual arrays are constructed.
    pieces = record.get("pieces")
    require(isinstance(pieces, list) and len(pieces) == steps, "piece count mismatch")
    for piece in pieces:
        require(isinstance(piece, list) and len(piece) == degree+1, "degree mismatch")
        for row in piece:
            require(isinstance(row, list) and len(row) == cutoff+1, "mode count mismatch")
            for value in row:
                require(isinstance(value, str) and len(value) <= bits+3,
                        "coefficient string length bound")
                require(INTEGER.fullmatch(value) is not None, "coefficient is not an integer string")
    denominator = 1 << bits
    table = [[[int(v) for v in row] for row in piece] for piece in pieces]
    require(all(abs(v) <= 3*denominator for piece in table for row in piece for v in row),
            "coefficient magnitude bound")
    projected = (initial + [Fraction(0)]*(cutoff+1))[:cutoff+1]
    initial_error = sum((1+j)*abs(Fraction(table[0][0][j], denominator)-projected[j])
                        for j in range(cutoff+1))
    initial_error4 = sum((1 << (2*j))*abs(Fraction(table[0][0][j], denominator)-projected[j])
                         for j in range(cutoff+1))
    norm_max = Fraction(0)
    residual_max = Fraction(0)
    residual4_max = Fraction(0)
    piece_residuals = []
    previous_endpoint = None
    weights4 = [1 << (2*j) for j in range(cutoff+1)]
    for piece in table:
        norm = Fraction(sum(weights4[j]*abs(v) for row in piece for j,v in enumerate(row)),
                        denominator)
        require(norm <= 3, "polynomial A4 budget fails")
        norm_max = max(norm_max, norm)
        if previous_endpoint is not None:
            require(previous_endpoint == piece[0], "trajectory is discontinuous")
        previous_endpoint = [sum(row[j] for row in piece) for j in range(cutoff+1)]
        active = [(ell, [(j,v) for j,v in enumerate(row) if v])
                  for ell,row in enumerate(piece) if any(row)]
        residual = [[0]*(cutoff+1) for _ in range(2*degree+1)]
        for ell,row in active:
            for j,v in row:
                residual[ell][j] += 2*denominator*j*j*v
                if ell:
                    residual[ell-1][j] += 2*denominator*16*steps*ell*v
        # Paired time indices account for both ordered factors; spatial
        # coefficients are expanded directly by the cosine product law.
        for a,(ell,left) in enumerate(active):
            for other,right in active[a:]:
                multiplicity = 1 if ell == other else 2
                target = residual[ell+other]
                for j,x in left:
                    for k,y in right:
                        value = multiplicity*x*y
                        target[abs(j-k)] -= value
                        if j+k <= cutoff:
                            target[j+k] -= value
        residual_denominator = 2*denominator*denominator
        bound = Fraction(sum((1+j)*abs(v) for row in residual for j,v in enumerate(row)),
                         residual_denominator)
        bound4 = Fraction(sum(weights4[j]*abs(v) for row in residual for j,v in enumerate(row)),
                          residual_denominator)
        piece_residuals.append(bound)
        residual_max = max(residual_max, bound)
        residual4_max = max(residual4_max, bound4)

    spatial_error = Fraction(5, 3*(1 << (cutoff+1)))
    time_error = (initial_error + T*residual_max)/(1-5*T)
    total_error = spatial_error + time_error
    target = Fraction(1, 1 << (expected_m+2))
    require(total_error < target, "certified total exceeds target")
    endpoint = [Fraction(v, denominator) for v in previous_endpoint]
    report = {"status": "ACCEPT", "m": expected_m, "mode": mode, **parameters,
        "input_reference": canonical["initial"],
        "max_polynomial_A4": fraction_report(norm_max),
        "initial_error_Alin": fraction_report(initial_error),
        "initial_error_A4": fraction_report(initial_error4),
        "piece_residuals_Alin": [fraction_report(x) for x in piece_residuals],
        "max_residual_Alin": fraction_report(residual_max),
        "max_residual_A4": fraction_report(residual4_max),
        "spatial_error": fraction_report(spatial_error),
        "time_and_rounding_error": fraction_report(time_error),
        "total_error": fraction_report(total_error), "target": fraction_report(target),
        "endpoint_observations": {"mean": fraction_report(endpoint[0]),
            "value_at_zero": fraction_report(sum(endpoint)),
            "common_absolute_error": fraction_report(total_error)}}
    if mode == "guaranteed":
        xi = guaranteed_parameters(expected_m)["xi"]
        report["theorem_strong_checks"] = {
            "xi": str(xi), "initial_below_xi": initial_error4 <= xi,
            "residual_below_xi": residual4_max <= xi}
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("record", type=Path)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--expected-m", type=int, required=True)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    started = time.perf_counter()
    try:
        external = load_input(args.input)
        result = verify(load_record(args.record, external, args.expected_m), args.expected_m, external)
    except (Rejected, ValueError, TypeError, KeyError, RecursionError) as error:
        result = {"status": "REJECT", "reason": str(error)}
    result["verifier_seconds"] = time.perf_counter()-started
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    raise SystemExit(0 if result["status"] == "ACCEPT" else 1)


if __name__ == "__main__":
    main()
