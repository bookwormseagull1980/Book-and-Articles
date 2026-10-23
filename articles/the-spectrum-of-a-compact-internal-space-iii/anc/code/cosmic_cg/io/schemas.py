"""Shared schemas for generated spectra, objects, branch links, and run outputs."""

from __future__ import annotations

from dataclasses import asdict, dataclass, fields
from typing import Any, Iterable


SCHEMA_VERSION = "0.1"


@dataclass(frozen=True)
class PowerSpectrum:
    source: str
    snapshot: str
    z: float
    k_hmpc: list[float]
    p_mpc3_h3: list[float]
    kind: str
    box_size_mpch: float


@dataclass(frozen=True)
class HaloRecord:
    source: str
    snapshot: str
    z: float
    halo_id: str
    mass_msun_h: float
    radius_kpc_h: float
    position_x: float
    position_y: float
    position_z: float
    velocity_x: float
    velocity_y: float
    velocity_z: float
    n_particles: int
    sigma_v_kms: float | None = None
    vmax_kms: float | None = None
    spin_x: float | None = None
    spin_y: float | None = None
    spin_z: float | None = None
    parent_id: str | None = None
    descendant_id: str | None = None
    quality_flag: str = "ok"


@dataclass(frozen=True)
class MergerLink:
    source: str
    snapshot_from: str
    snapshot_to: str
    z_from: float
    z_to: float
    time_from_gyr: float
    time_to_gyr: float
    halo_id_from: str
    halo_id_to: str
    mass_from_msun_h: float
    mass_to_msun_h: float


def dataclass_to_row(record: Any) -> dict[str, Any]:
    return asdict(record)


def rows_to_dataclasses(cls: type, rows: Iterable[dict[str, str]]) -> list[Any]:
    names = {field.name: field for field in fields(cls)}
    parsed = []
    for row in rows:
        values: dict[str, Any] = {}
        for name, field in names.items():
            if name not in row:
                continue
            values[name] = _parse_value(row[name], field.type)
        parsed.append(cls(**values))
    return parsed


def _parse_value(value: str, target_type: Any) -> Any:
    if value == "":
        return None
    text = str(target_type)
    if "int" in text and "float" not in text:
        return int(float(value))
    if "float" in text:
        return float(value)
    return value
