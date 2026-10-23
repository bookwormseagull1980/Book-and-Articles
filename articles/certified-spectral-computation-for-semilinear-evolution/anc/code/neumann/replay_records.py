"""Check every published trajectory without importing its constructor."""

from fractions import Fraction
import hashlib
import json
from pathlib import Path
from record_format import load_input, load_record
from verify import verify


def main():
    results = Path(__file__).resolve().parents[2]/"results"/"neumann"
    count = 0
    baseline = json.loads((results/"summary.json").read_text(encoding="utf-8"))
    for row in baseline["rows"]:
        path = results/f'm{row["m"]}.record.json'
        external = load_input(results/"baseline.input.json")
        report = verify(load_record(path,external,row["m"]),row["m"],external)
        assert report["total_error"]["exact"] == row["total_error"]["exact"]
        assert hashlib.sha256(path.read_bytes()).hexdigest() == row["record_sha256"]
        count += 1
    matrix = json.loads((results/"matrix_summary.json").read_text(encoding="utf-8"))
    for row in matrix["rows"]:
        path = results/"matrix"/row["record_file"]
        external = load_input(results/"matrix"/row["input_file"])
        report = verify(load_record(path,external,row["m"]),row["m"],external)
        assert report["total_error"]["exact"] == row["total_error"]["exact"]
        assert hashlib.sha256(path.read_bytes()).hexdigest() == row["record_sha256"]
        if row["case"] == "E2-constant":
            center = Fraction(report["endpoint_observations"]["mean"]["exact"])
            assert abs(center-Fraction(16,31)) <= Fraction(report["time_and_rounding_error"]["exact"])
        count += 1
    print(json.dumps({"status":"ACCEPT","records":count,
                      "construction_imported":False,"arithmetic":"exact"}))


if __name__ == "__main__":
    main()
