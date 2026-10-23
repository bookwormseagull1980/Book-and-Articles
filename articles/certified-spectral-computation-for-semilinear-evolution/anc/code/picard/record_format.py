"""Finite formats for a coupled linear Volterra calculation."""

from fractions import Fraction
import json
from pathlib import Path
import re
import sys

if hasattr(sys, "set_int_max_str_digits"):
    sys.set_int_max_str_digits(0)

EQUATION = "y'=B*y+p'; y(0)=p(0); t in [0,1]"
INTEGER = re.compile(r"-?(0|[1-9][0-9]*)\Z")
RATIONAL = re.compile(r"-?(0|[1-9][0-9]*)(/[1-9][0-9]*)?\Z")


class Rejected(ValueError):
    pass


def require(condition, message):
    if not condition:
        raise Rejected(message)


def vector_norm(row):
    return max(abs(x) for x in row)


def weighted_norm(coefficients):
    return sum((1 << ell)*vector_norm(row) for ell, row in enumerate(coefficients))


def normalize_input(data):
    require(isinstance(data, dict), "input must be an object")
    require(data.get("format") == "picard-linear-input-v1", "unknown input format")
    require(data.get("equation_id") == EQUATION, "wrong equation")
    matrix = data.get("matrix_B")
    polynomial = data.get("polynomial_p")
    require(isinstance(matrix, list) and len(matrix) >= 1, "matrix must be dense and nonempty")
    dimension = len(matrix)
    require(all(isinstance(row, list) and len(row) == dimension for row in matrix),
            "matrix must be square")
    require(isinstance(polynomial, list) and len(polynomial) >= 1,
            "polynomial coefficient table must be dense and nonempty")
    require(all(isinstance(row, list) and len(row) == dimension for row in polynomial),
            "polynomial vector dimension mismatch")
    require(all(isinstance(v, str) and RATIONAL.fullmatch(v)
                for table in (matrix, polynomial) for row in table for v in row),
            "input coefficients must be rational strings")
    matrix = [[Fraction(v) for v in row] for row in matrix]
    polynomial = [[Fraction(v) for v in row] for row in polynomial]
    require(max(sum(abs(v) for v in row) for row in matrix) <= Fraction(1,8),
            "matrix row-sum budget fails")
    require(weighted_norm(polynomial) <= Fraction(1,2), "polynomial source budget fails")
    canonical = {"format":"picard-linear-input-v1", "equation_id":EQUATION,
                 "matrix_B":[[str(v) for v in row] for row in matrix],
                 "polynomial_p":[[str(v) for v in row] for row in polynomial]}
    encoded_bytes = len(json.dumps(canonical,separators=(",",":")).encode("utf-8"))
    return canonical, matrix, polynomial, encoded_bytes


def make_input(matrix, polynomial):
    return normalize_input({"format":"picard-linear-input-v1", "equation_id":EQUATION,
        "matrix_B":[[str(v) for v in row] for row in matrix],
        "polynomial_p":[[str(v) for v in row] for row in polynomial]})[0]


def load_input(path):
    def reject_number(_):
        raise Rejected("input coefficients must be encoded as rational strings")
    with Path(path).open("r", encoding="utf-8") as stream:
        data = json.load(stream,parse_int=reject_number,parse_float=reject_number,
                         parse_constant=reject_number)
    return normalize_input(data)[0]


def parameters(dimension, source_degree, m):
    require(type(m) is int and m >= 1, "precision must be a positive integer")
    require(type(dimension) is int and dimension >= 1, "invalid external dimension")
    require(type(source_degree) is int and source_degree >= 0, "invalid external source degree")
    b = m + 4
    degree = max(source_degree, b+3)
    return {"dimension":dimension, "source_degree":source_degree, "m":m,
            "working_bits":b, "degree":degree, "fractional_bits":b+degree+4,
            "iterations":2*m+8}


def limits(specification, input_bytes):
    p = specification
    entries = ((p["iterations"]+1)*(p["degree"]+1)+1)*p["dimension"]
    return {"integer_characters":p["fractional_bits"]+3,
            "record_bytes":8192+2*input_bytes+entries*(p["fractional_bits"]+7)}


def load_record(path, expected_input, expected_m):
    _, matrix, polynomial, input_bytes = normalize_input(expected_input)
    p = parameters(len(matrix),len(polynomial)-1,expected_m)
    bounds = limits(p,input_bytes)
    path = Path(path)
    require(path.stat().st_size <= bounds["record_bytes"], "record byte bound exceeded")
    digits = len(str(max(p.values())))

    def bounded_integer(text):
        require(len(text.lstrip("-")) <= digits, "oversized JSON header integer")
        return int(text)

    def reject_number(_):
        raise Rejected("noninteger JSON numeric field")

    with path.open("r",encoding="utf-8") as stream:
        return json.load(stream,parse_int=bounded_integer,parse_float=reject_number,
                         parse_constant=reject_number)
