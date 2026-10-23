"""Run-grid definitions for internally generated universe histories."""

from __future__ import annotations

from typing import Iterable

from cosmic_cg.result_tables import DEFAULT_MASS_BANDS


PAST_REDSHIFTS: tuple[float, ...] = (
    100.0,
    60.0,
    40.0,
    30.0,
    20.0,
    15.0,
    12.0,
    10.0,
    8.0,
    6.0,
    5.0,
    4.0,
    3.0,
    2.0,
    1.5,
    1.0,
    0.7,
    0.5,
    0.3,
    0.0,
)

FUTURE_SCALE_FACTORS: tuple[float, ...] = (1.25, 1.5, 2.0, 3.0)

RUN_GRID_COLUMNS = (
    "phase",
    "run_id",
    "z",
    "a",
    "epoch_role",
    "mass_band_id",
    "mass_min_msun_h",
    "mass_max_msun_h",
    "object_target",
    "spectrum_mode",
    "object_mode",
    "transfer_mode",
    "expected_local_mb",
    "expected_wall_time",
    "output_scope",
)


def generation_redshifts() -> tuple[float, ...]:
    future = tuple(1.0 / a - 1.0 for a in FUTURE_SCALE_FACTORS)
    return PAST_REDSHIFTS + future


def production_run_grid(
    *,
    objects_per_band: int = 24,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for z in generation_redshifts():
        rows.append(
            _row(
                phase="window_history",
                run_id=f"window_z{_tag(z)}",
                z=z,
                a=1.0 / (1.0 + z),
                epoch_role=_epoch_role(z),
                mass_band_id="all",
                mass_min_msun_h="",
                mass_max_msun_h="",
                object_target=0,
                spectrum_mode="internal_finite_window_envelope",
                object_mode="none",
                transfer_mode="none",
                expected_local_mb=1.0,
                expected_wall_time="seconds",
                output_scope="window, band, tail",
            )
        )
        for band_id, mass_min, mass_max, _label in DEFAULT_MASS_BANDS:
            rows.append(
                _row(
                    phase="object_history",
                    run_id=f"object_z{_tag(z)}_{band_id}",
                    z=z,
                    a=1.0 / (1.0 + z),
                    epoch_role=_epoch_role(z),
                    mass_band_id=band_id,
                    mass_min_msun_h=mass_min,
                    mass_max_msun_h=mass_max,
                    object_target=objects_per_band,
                    spectrum_mode="internal_finite_window_envelope",
                    object_mode="generated_representative_ensemble",
                    transfer_mode="generated_branch_transfer,generated_pair_tidal_transfer",
                    expected_local_mb=10.0,
                    expected_wall_time="seconds to minutes",
                    output_scope="object fate, transfer channel, summaries",
                )
            )
    return rows


def render_run_grid_markdown(rows: Iterable[dict[str, object]]) -> str:
    rows = list(rows)
    phases = _count_by(rows, "phase")
    lines = [
        "# Generated Universe Run Grid",
        "",
        "This file records the internal finite-window generation grid.",
        "",
        "## Phase Counts",
        "",
        "| Phase | Rows |",
        "|---|---:|",
    ]
    for phase, count in phases:
        lines.append(f"| {phase} | {count} |")
    lines.extend(
        [
            "",
            "## Epochs",
            "",
            "| Redshift | Scale factor | Role |",
            "|---:|---:|---|",
        ]
    )
    for z in generation_redshifts():
        lines.append(f"| {z:.6g} | {1.0 / (1.0 + z):.6g} | {_epoch_role(z)} |")
    lines.extend(
        [
            "",
            "## Output Policy",
            "",
            "The generated-universe run writes canonical CSV tables, intermediate CSV checkpoints, JSON provenance, and an audit report. Figure-source tables are written after table verification.",
            "",
        ]
    )
    return "\n".join(lines)


def _row(**items: object) -> dict[str, object]:
    return {column: items.get(column, "") for column in RUN_GRID_COLUMNS}


def _epoch_role(z: float) -> str:
    if z > 0:
        return "past"
    if z == 0:
        return "present"
    return "future"


def _tag(value: float) -> str:
    return f"{value:g}".replace("-", "m").replace(".", "p")


def _count_by(rows: list[dict[str, object]], key: str) -> list[tuple[str, int]]:
    counts: dict[str, int] = {}
    for row in rows:
        value = str(row[key])
        counts[value] = counts.get(value, 0) + 1
    return sorted(counts.items())
