"""Replay a finite Picard execution by independent rational checks.

The constructor is not imported.  Each update is checked against its exact
coefficient formula and the defining interval for truncation towards zero.
"""

import argparse
from fractions import Fraction
import json
from pathlib import Path
import time
from record_format import (EQUATION, INTEGER, Rejected, load_input, load_record,
                           normalize_input, parameters, require, vector_norm,
                           weighted_norm)


def reported(value):
    return {"exact":str(value),"decimal":float(value)}


def verify(record, expected_input, expected_m):
    canonical, matrix, polynomial, _ = normalize_input(expected_input)
    p = parameters(len(matrix),len(polynomial)-1,expected_m)
    require(isinstance(record,dict),"record must be an object")
    require(record.get("format") == "picard-execution-v1","unknown record format")
    require(record.get("equation_id") == EQUATION,"wrong equation")
    require(record.get("input_reference") == canonical,"external input mismatch")
    for key,value in p.items():
        require(type(record.get(key)) is int and record[key] == value,"parameter mismatch: "+key)
    dimension, degree, beta = p["dimension"],p["degree"],p["fractional_bits"]
    tables = record.get("iterates")
    require(isinstance(tables,list) and len(tables) == p["iterations"]+1,"iterate count mismatch")
    for table in tables:
        require(isinstance(table,list) and len(table) == degree+1,"degree mismatch")
        for row in table:
            require(isinstance(row,list) and len(row) == dimension,"vector dimension mismatch")
    terminal = record.get("terminal")
    require(isinstance(terminal,list) and len(terminal) == dimension,"terminal dimension mismatch")
    for table in tables+[ [terminal] ]:
        for row in table:
            for value in row:
                require(isinstance(value,str) and len(value) <= beta+3,"coefficient string length bound")
                require(INTEGER.fullmatch(value) is not None,"coefficient is not an integer string")
    denominator = 1 << beta
    converted = [[[int(v) for v in row] for row in table] for table in tables]
    terminal = [int(v) for v in terminal]
    require(all(abs(v) <= denominator for table in converted for row in table for v in row)
            and all(abs(v) <= denominator for v in terminal),"coefficient magnitude bound")
    require(all(v == 0 for row in converted[0] for v in row),"initial iterate is not zero")
    strong_bounds = []
    for table in converted:
        bound = Fraction(sum((1 << ell)*max(abs(v) for v in row)
                             for ell,row in enumerate(table)),denominator)
        require(bound <= 1,"strong invariant ball fails")
        strong_bounds.append(bound)
    step_errors = []
    strong_step_errors = []
    exact_image_bounds = []
    error = Fraction(1)
    for iteration in range(p["iterations"]):
        previous = converted[iteration]
        following = converted[iteration+1]
        image = []
        for ell in range(degree+2):
            row = [polynomial[ell][j] if ell < len(polynomial) else Fraction(0)
                   for j in range(dimension)]
            if ell:
                for j in range(dimension):
                    for k in range(dimension):
                        row[j] += matrix[j][k]*Fraction(previous[ell-1][k],ell*denominator)
            image.append(row)
        image_bound = weighted_norm(image)
        require(image_bound <= Fraction(3,4),"exact image bound fails")
        exact_image_bounds.append(image_bound)
        differences = []
        for ell in range(degree+1):
            for j, exact in enumerate(image[ell]):
                stored = following[ell][j]
                scaled = exact*denominator
                if scaled >= 0:
                    valid = 0 <= stored <= scaled < stored+1
                else:
                    valid = stored-1 < scaled <= stored <= 0
                require(valid,"update does not match truncation towards zero")
            differences.append([Fraction(following[ell][j],denominator)-image[ell][j]
                                for j in range(dimension)])
        differences.append([-v for v in image[-1]])
        step_error = sum(vector_norm(row) for row in differences)
        strong_step_error = weighted_norm(differences)
        require(step_error <= Fraction(1,1 << p["working_bits"]),"weak update error exceeds budget")
        step_errors.append(step_error)
        strong_step_errors.append(strong_step_error)
        error = error/8+step_error
    expected_endpoint = [sum(row[j] for row in converted[-1]) for j in range(dimension)]
    require(terminal == expected_endpoint,"terminal value mismatch")
    theorem_error = Fraction(1,8**p["iterations"])+Fraction(8,7*(1 << p["working_bits"]))
    target = Fraction(1,1 << (expected_m+2))
    require(error <= theorem_error < target,"global approximation budget fails")
    centers = [Fraction(v,denominator) for v in terminal]
    return {"status":"ACCEPT",**p,"max_strong_norm":reported(max(strong_bounds)),
        "max_exact_image_norm":reported(max(exact_image_bounds)),
        "max_weak_step_error":reported(max(step_errors)),
        "max_strong_step_difference":reported(max(strong_step_errors)),
        "step_errors":[str(v) for v in step_errors],
        "certified_path_error":reported(error),"uniform_path_error":reported(theorem_error),
        "target":reported(target),"terminal_centers":[str(v) for v in centers],
        "terminal_intervals":[[str(v-error),str(v+error)] for v in centers]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("record",type=Path)
    parser.add_argument("--input",type=Path,required=True)
    parser.add_argument("--expected-m",type=int,required=True)
    parser.add_argument("--report",type=Path)
    args = parser.parse_args()
    started = time.perf_counter()
    try:
        external = load_input(args.input)
        report = verify(load_record(args.record,external,args.expected_m),external,args.expected_m)
    except (Rejected,ValueError,TypeError,KeyError,RecursionError) as error:
        report = {"status":"REJECT","reason":str(error)}
    report["verifier_seconds"] = time.perf_counter()-started
    if args.report:
        args.report.parent.mkdir(parents=True,exist_ok=True)
        args.report.write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))
    raise SystemExit(0 if report["status"] == "ACCEPT" else 1)


if __name__ == "__main__":
    main()
