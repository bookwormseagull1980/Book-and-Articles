"""Rebuild the three records and invoke the standalone exact verifier."""

import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys

from format_table import write_rows
from record_format import DEFAULT_INITIAL, input_from_coefficients


def main():
    source = Path(__file__).resolve().parent
    results = source.parents[1] / "results" / "neumann"
    results.mkdir(parents=True, exist_ok=True)
    input_file = results / "baseline.input.json"
    input_file.write_text(json.dumps(input_from_coefficients(DEFAULT_INITIAL), indent=2)+"\n", encoding="utf-8")
    rows = []
    for m in (4, 6, 8):
        record = results / f"m{m}.record.json"
        report = results / f"m{m}.report.json"
        produced = subprocess.run([
            sys.executable, str(source / "produce.py"), "--input", str(input_file), "--m", str(m),
            "--output", str(record)], check=True, text=True, capture_output=True)
        production = json.loads(produced.stdout)
        checked = subprocess.run([
            sys.executable, str(source / "verify.py"), str(record),
            "--input", str(input_file), "--expected-m", str(m), "--report", str(report)],
            check=True, text=True, capture_output=True)
        verification = json.loads(checked.stdout)
        rows.append({
            "m": m,
            "cutoff": verification["cutoff"],
            "steps": verification["steps"],
            "degree": verification["degree"],
            "fractional_bits": verification["fractional_bits"],
            "max_polynomial_A4": verification["max_polynomial_A4"],
            "max_residual_Alin": verification["max_residual_Alin"],
            "spatial_error": verification["spatial_error"],
            "time_and_rounding_error": verification["time_and_rounding_error"],
            "total_error": verification["total_error"],
            "target": verification["target"],
            "producer_seconds": production["producer_seconds"],
            "verifier_seconds": verification["verifier_seconds"],
            "record_bytes": record.stat().st_size,
            "record_sha256": hashlib.sha256(record.read_bytes()).hexdigest(),
        })
    summary = {
        "python": platform.python_version(),
        "platform": platform.platform(),
        "arithmetic": "exact rational and integer; floating output is descriptive only",
        "timings": "single local runs; platform dependent; excludes interpreter startup",
        "rows": rows,
    }
    (results / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    write_rows(summary, results / "table_rows.tex")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
