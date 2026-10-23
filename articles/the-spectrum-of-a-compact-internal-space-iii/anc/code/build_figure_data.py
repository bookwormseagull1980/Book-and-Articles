"""Build standard figure-source tables for an existing run."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from cosmic_cg.io.audit import audit_run, write_audit_markdown, write_audit_report
from cosmic_cg.io.figure_data import write_standard_figure_data
from cosmic_cg.intermediate import write_output_inventory


def main() -> int:
    args = parse_args()
    run_dir = Path(args.run_dir)
    written = write_standard_figure_data(run_dir)
    write_output_inventory(run_dir)
    report = audit_run(run_dir, None if args.run_type == "auto" else args.run_type)
    qc_dir = run_dir / "quality_checks"
    write_audit_report(report, qc_dir / "audit_report.json")
    write_audit_markdown(report, qc_dir / "audit_report.md")
    print(
        json.dumps(
            {
                "run_dir": str(run_dir),
                "written": [str(path) for path in written],
                "audit_status": report.status,
            },
            indent=2,
        )
    )
    if args.strict and report.status == "fail":
        return 1
    return 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True)
    parser.add_argument(
        "--run-type",
        default="auto",
        choices=("auto", "generated_universe"),
    )
    parser.add_argument("--strict", action="store_true")
    return parser.parse_args()


if __name__ == "__main__":
    raise SystemExit(main())
