"""Audit a completed three-dimensional field-generation result directory."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from cosmic_cg.io.field_audit import audit_field_run, write_field_audit_markdown, write_field_audit_report


def main() -> int:
    args = parse_args()
    report = audit_field_run(args.run_dir, expected_grid_size=args.expected_grid_size)
    qc_dir = Path(args.run_dir) / "quality_checks"
    qc_dir.mkdir(parents=True, exist_ok=True)
    write_field_audit_report(report, qc_dir / "field_audit_report.json")
    write_field_audit_markdown(report, qc_dir / "field_audit_report.md")
    print(json.dumps({"run_dir": args.run_dir, "status": report.status}, indent=2))
    if args.strict and report.status == "fail":
        return 1
    return 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True)
    parser.add_argument("--expected-grid-size", type=int, default=None)
    parser.add_argument("--strict", action="store_true")
    return parser.parse_args()


if __name__ == "__main__":
    raise SystemExit(main())
