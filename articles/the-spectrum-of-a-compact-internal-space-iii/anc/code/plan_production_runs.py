"""Write the production run grid used by the finite-window analysis."""

from __future__ import annotations

import argparse
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from cosmic_cg.io.tables import write_rows
from cosmic_cg.run_plan import RUN_GRID_COLUMNS, production_run_grid, render_run_grid_markdown


def main() -> int:
    args = parse_args()
    rows = production_run_grid(objects_per_band=args.objects_per_band)
    output_csv = Path(args.output_csv)
    output_md = Path(args.output_md)
    if not output_csv.is_absolute():
        output_csv = PROJECT_ROOT / output_csv
    if not output_md.is_absolute():
        output_md = PROJECT_ROOT / output_md
    write_rows(output_csv, rows, RUN_GRID_COLUMNS)
    output_md.parent.mkdir(parents=True, exist_ok=True)
    output_md.write_text(render_run_grid_markdown(rows), encoding="utf-8")
    print(f"wrote {output_csv}")
    print(f"wrote {output_md}")
    return 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--objects-per-band", type=int, default=24)
    parser.add_argument("--output-csv", default=str(PROJECT_ROOT / "tables" / "production_run_grid.csv"))
    parser.add_argument("--output-md", default=str(PROJECT_ROOT / "PRODUCTION_RUN_GRID_2026-09-03.md"))
    return parser.parse_args()


if __name__ == "__main__":
    raise SystemExit(main())
