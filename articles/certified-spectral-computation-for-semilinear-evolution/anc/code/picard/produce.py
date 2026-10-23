"""Execute fixed Picard iterations and retain every finite path record."""

import argparse
from fractions import Fraction
import json
from pathlib import Path
import time
from record_format import EQUATION, load_input, normalize_input, parameters


def round_towards_zero(value, denominator):
    numerator = value.numerator*denominator
    return (-1 if numerator < 0 else 1)*(abs(numerator)//value.denominator)


def construct(external_input, m):
    canonical, matrix, polynomial, _ = normalize_input(external_input)
    dimension = len(matrix)
    p = parameters(dimension,len(polynomial)-1,m)
    degree, beta = p["degree"],p["fractional_bits"]
    denominator = 1 << beta
    current = [[0]*dimension for _ in range(degree+1)]
    iterates = [[[str(v) for v in row] for row in current]]
    for _ in range(p["iterations"]):
        following = []
        for ell in range(degree+1):
            exact = list(polynomial[ell]) if ell < len(polynomial) else [Fraction(0)]*dimension
            if ell:
                for i in range(dimension):
                    exact[i] += sum(matrix[i][j]*current[ell-1][j]
                                    for j in range(dimension))/(ell*denominator)
            following.append([round_towards_zero(v,denominator) for v in exact])
        current = following
        iterates.append([[str(v) for v in row] for row in current])
    endpoint = [sum(row[j] for row in current) for j in range(dimension)]
    return {"format":"picard-execution-v1", "equation_id":EQUATION,
            "input_reference":canonical, **p, "iterates":iterates,
            "terminal":[str(v) for v in endpoint]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input",type=Path,required=True)
    parser.add_argument("--m",type=int,required=True)
    parser.add_argument("--output",type=Path,required=True)
    args = parser.parse_args()
    started = time.perf_counter()
    record = construct(load_input(args.input),args.m)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(record,separators=(",",":"))+"\n",encoding="utf-8")
    print(json.dumps({"producer_seconds":time.perf_counter()-started,
                      "record_bytes":args.output.stat().st_size}))


if __name__ == "__main__":
    main()
