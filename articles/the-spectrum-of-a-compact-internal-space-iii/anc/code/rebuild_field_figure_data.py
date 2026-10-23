"""Rebuild field figure-source tables from stored relationship tables."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from cosmic_cg.intermediate import write_output_inventory
from cosmic_cg.io.figure_data import write_field_figure_data


def main() -> int:
    args = parse_args()
    run_dir = Path(args.run_dir)
    figure_dir = run_dir / "figure_data"
    written = write_field_figure_data(run_dir, figure_dir)
    _update_manifest(run_dir, written)
    write_output_inventory(run_dir)
    print(json.dumps({"run_dir": str(run_dir), "figure_data": [path.name for path in written]}, indent=2))
    return 0


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-dir", required=True)
    return parser.parse_args()


def _update_manifest(run_dir: Path, written: list[Path]) -> None:
    manifest_path = run_dir / "generation_manifest.json"
    if not manifest_path.is_file():
        return
    payload = json.loads(manifest_path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict):
        return
    payload["figure_data"] = [path.name for path in written]
    manifest_path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")


if __name__ == "__main__":
    raise SystemExit(main())
