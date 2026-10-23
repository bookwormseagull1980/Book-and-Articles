"""Build interaction-ledger tables from an existing generated run."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from cosmic_cg.interaction_adapters import write_interaction_tables_for_run
from cosmic_cg.intermediate import write_output_inventory


def main() -> int:
    args = parse_args()
    run_dir = Path(args.run_dir)
    counts = write_interaction_tables_for_run(run_dir)
    write_output_inventory(run_dir)
    print(json.dumps({"run_dir": str(run_dir), "written": counts}, indent=2, sort_keys=True))
    return 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True, help="existing run directory")
    return parser.parse_args()


if __name__ == "__main__":
    raise SystemExit(main())
