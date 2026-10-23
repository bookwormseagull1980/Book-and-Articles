"""Run provenance and versioned output directories."""

from __future__ import annotations

from dataclasses import dataclass
import json
from pathlib import Path
import platform
import sys
from typing import Any

from cosmic_cg import __version__


@dataclass(frozen=True)
class RunPaths:
    run_id: str
    root: Path
    intermediate: Path
    tables: Path
    figure_data: Path
    figures: Path
    logs: Path
    quality_checks: Path


def create_run_paths(results_root: str | Path, run_id: str) -> RunPaths:
    root = Path(results_root) / run_id
    if root.exists() and any(root.iterdir()):
        raise FileExistsError(f"run directory already exists: {root}")
    paths = RunPaths(
        run_id=run_id,
        root=root,
        intermediate=root / "intermediate",
        tables=root / "tables",
        figure_data=root / "figure_data",
        figures=root / "figures",
        logs=root / "logs",
        quality_checks=root / "quality_checks",
    )
    for path in paths.__dict__.values():
        if isinstance(path, Path):
            path.mkdir(parents=True, exist_ok=True)
    return paths


def write_provenance(path: str | Path, *, run_id: str, config: dict[str, Any]) -> None:
    payload = {
        "run_id": run_id,
        "package_version": __version__,
        "python": sys.version,
        "platform": platform.platform(),
        "config": config,
    }
    Path(path).write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
