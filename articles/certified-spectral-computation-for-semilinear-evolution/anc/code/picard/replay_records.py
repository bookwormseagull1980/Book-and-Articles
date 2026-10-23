"""Verify every published Picard record without importing its constructor."""

from fractions import Fraction
import hashlib
import json
from pathlib import Path
from benchmarks import coefficient_error, exact_nilpotent_solution
from record_format import load_input, load_record, normalize_input
from verify import verify


def main():
    results = Path(__file__).resolve().parents[2]/"results"/"picard"
    summary = json.loads((results/"summary.json").read_text(encoding="utf-8"))
    count = exact_checks = 0
    for row in summary["rows"]:
        external = load_input(results/row["input_file"])
        path = results/row["record_file"]
        record = load_record(path,external,row["m"])
        report = verify(record,external,row["m"])
        assert report["certified_path_error"]["exact"] == row["certified_path_error"]["exact"]
        assert hashlib.sha256(path.read_bytes()).hexdigest() == row["record_sha256"]
        if "exact_polynomial_error_bound" in row:
            _,matrix,polynomial,_ = normalize_input(external)
            exact = exact_nilpotent_solution(matrix,polynomial)
            assert coefficient_error(record,exact) <= Fraction(report["certified_path_error"]["exact"])
            exact_checks += 1
        count += 1
    print(json.dumps({"status":"ACCEPT","records":count,"exact_solution_checks":exact_checks,
                      "constructor_imported":False}))


if __name__ == "__main__":
    main()
