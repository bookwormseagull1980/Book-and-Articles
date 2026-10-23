"""Build a compact index of local result directories."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from cosmic_cg.io.audit import audit_run
from cosmic_cg.io.tables import write_rows


INDEX_COLUMNS = (
    "run_id",
    "run_type",
    "audit_status",
    "issue_count",
    "has_intermediate",
    "has_global_window",
    "has_halo_objects",
    "has_window_table",
    "has_band_table",
    "has_tail_table",
    "has_object_table",
    "has_transfer_table",
    "has_cost_table",
    "has_schema_catalog",
    "has_generation_manifest",
    "has_figure_data",
    "result_dir",
)


def main() -> int:
    args = parse_args()
    results_root = Path(args.results_root)
    rows = []
    for run_dir in sorted(path for path in results_root.iterdir() if path.is_dir()):
        report = audit_run(run_dir)
        if report.run_type != "generated_universe" and not args.include_legacy:
            continue
        rows.append(
            {
                "run_id": run_dir.name,
                "run_type": report.run_type,
                "audit_status": report.status,
                "issue_count": len(report.issues),
                "has_intermediate": has_files(run_dir / "intermediate", "*.csv"),
                "has_global_window": (run_dir / "tables" / "global_window.csv").is_file(),
                "has_halo_objects": (run_dir / "tables" / "halo_objects.csv").is_file(),
                "has_window_table": (run_dir / "tables" / "window_table.csv").is_file(),
                "has_band_table": (run_dir / "tables" / "band_table.csv").is_file(),
                "has_tail_table": (run_dir / "tables" / "tail_table.csv").is_file(),
                "has_object_table": (run_dir / "tables" / "object_table.csv").is_file(),
                "has_transfer_table": (run_dir / "tables" / "transfer_table.csv").is_file(),
                "has_cost_table": (run_dir / "tables" / "cost_table.csv").is_file(),
                "has_schema_catalog": (run_dir / "tables" / "schema_catalog.csv").is_file(),
                "has_generation_manifest": (run_dir / "generation_manifest.json").is_file(),
                "has_figure_data": has_files(run_dir / "figure_data", "*.csv"),
                "result_dir": str(run_dir),
            }
        )
    output = Path(args.output)
    if not output.is_absolute():
        output = PROJECT_ROOT / output
    write_rows(output, rows, INDEX_COLUMNS)
    print(output)
    return 0


def has_files(directory: Path, pattern: str) -> bool:
    return directory.is_dir() and any(directory.glob(pattern))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results-root", default=str(PROJECT_ROOT / "results"))
    parser.add_argument("--output", default=str(PROJECT_ROOT / "tables" / "results_index.csv"))
    parser.add_argument("--include-legacy", action="store_true")
    return parser.parse_args()


if __name__ == "__main__":
    raise SystemExit(main())
