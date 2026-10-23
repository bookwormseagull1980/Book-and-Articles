"""Internal particle-mesh dynamics for trajectory audits."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import ctypes
import math
from typing import Callable, Mapping, Sequence

from cosmic_cg.generation import cosmic_age_gyr, omega_m0
from cosmic_cg.object_level import alignment_cosine, dynamical_time_gyr, halo_fate
from cosmic_cg.params import CosmicCGParams
from cosmic_cg.cosmological_dynamics import drift_kick_integrals, expansion, growing_rate

try:
    import numpy as np
except ModuleNotFoundError:  # pragma: no cover
    np = None  # type: ignore[assignment]


RHO_CRIT_MSUNH_MPCH3 = 2.77536627e11


@dataclass(frozen=True)
class PMNBodyConfig:
    dynamics: str = "cosmological"
    particles_per_dim: int = 32
    mesh_size: int = 64
    box_size_mpch: float = 128.0
    seed: int = 20260904
    z_initial: float = 120.0
    steps_per_interval: int = 4
    initial_displacement_fraction: float = 0.22
    initial_velocity_factor: float = 0.35
    force_strength: float = 0.35
    mode_bins: int = 24
    density_threshold_sigma: float = 1.5
    density_threshold_delta: float = 0.8
    fof_min_particles: int = 12
    max_halos_per_snapshot: int = 96
    neighbor_links: int = 3


PM_SNAPSHOT_COLUMNS = (
    "run_id",
    "source",
    "snapshot",
    "step",
    "z",
    "a",
    "particle_count",
    "mesh_size",
    "mesh_cell_count",
    "box_size_mpch",
    "dx_mesh_mpch",
    "particle_mass_msunh",
    "mean_delta",
    "rms_delta",
    "min_delta",
    "max_delta",
    "kinetic_proxy",
    "potential_proxy",
    "total_energy_proxy",
    "mass_conservation_error",
    "momentum_norm_proxy",
    "halo_count",
    "tree_link_count",
    "transfer_count",
    "quality_flag",
)

PM_POWER_COLUMNS = (
    "run_id",
    "source",
    "snapshot",
    "z",
    "a",
    "mode_bin",
    "k_min_hmpc",
    "k_max_hmpc",
    "k_mid_hmpc",
    "mode_count",
    "power_proxy",
    "power_fraction",
    "cumulative_power_fraction",
    "quality_flag",
)

PM_HALO_COLUMNS = (
    "run_id",
    "source",
    "snapshot",
    "z",
    "a",
    "halo_id",
    "rank",
    "particle_count",
    "mass_msunh",
    "radius_kpch",
    "x_mpch",
    "y_mpch",
    "z_mpch",
    "vx_proxy",
    "vy_proxy",
    "vz_proxy",
    "spin_x",
    "spin_y",
    "spin_z",
    "parent_id",
    "descendant_id",
    "match_particle_fraction",
    "b_h",
    "t_dyn_gyr",
    "t_growth_gyr",
    "fate_class",
    "quality_flag",
)

PM_TREE_COLUMNS = (
    "run_id",
    "source",
    "snapshot_from",
    "snapshot_to",
    "z_from",
    "z_to",
    "a_from",
    "a_to",
    "parent_id",
    "descendant_id",
    "shared_particle_count",
    "parent_particle_count",
    "descendant_particle_count",
    "overlap_fraction",
    "mass_from_msunh",
    "mass_to_msunh",
    "delta_mass_msunh",
    "delta_time_gyr",
    "b_h",
    "fate_class",
    "quality_flag",
)

PM_TRANSFER_COLUMNS = (
    "run_id",
    "source",
    "snapshot",
    "z",
    "a",
    "source_id",
    "target_id",
    "rank_source",
    "rank_target",
    "separation_mpch",
    "mass_ratio",
    "a_pq",
    "g_pq",
    "tau_x",
    "tau_y",
    "tau_z",
    "j_x",
    "j_y",
    "j_z",
    "source_channel",
    "quality_flag",
)

PM_COST_COLUMNS = (
    "run_id",
    "source",
    "particles_per_dim",
    "mesh_size",
    "redshift_count",
    "particle_count",
    "mesh_cell_count",
    "box_size_mpch",
    "estimated_peak_memory_gb",
    "raw_particle_bytes_per_snapshot",
    "raw_mesh_bytes_per_snapshot",
    "raw_bytes_all_snapshots",
    "relationship_rows",
    "relationship_bytes_written",
    "saved_to_raw_ratio",
    "storage_policy",
)


def require_numpy():
    if np is None:
        raise RuntimeError("NumPy is required for the particle-mesh solver")
    return np


def config_to_dict(config: PMNBodyConfig) -> dict[str, object]:
    return asdict(config)


def available_physical_memory_bytes() -> int | None:
    class MemoryStatusEx(ctypes.Structure):
        _fields_ = [
            ("dwLength", ctypes.c_ulong),
            ("dwMemoryLoad", ctypes.c_ulong),
            ("ullTotalPhys", ctypes.c_ulonglong),
            ("ullAvailPhys", ctypes.c_ulonglong),
            ("ullTotalPageFile", ctypes.c_ulonglong),
            ("ullAvailPageFile", ctypes.c_ulonglong),
            ("ullTotalVirtual", ctypes.c_ulonglong),
            ("ullAvailVirtual", ctypes.c_ulonglong),
            ("ullAvailExtendedVirtual", ctypes.c_ulonglong),
        ]

    status = MemoryStatusEx()
    status.dwLength = ctypes.sizeof(MemoryStatusEx)
    try:
        ok = ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(status))
    except AttributeError:
        return None
    if not ok:
        return None
    return int(status.ullAvailPhys)


def estimate_pm_resources(
    *,
    redshift_count: int,
    config: PMNBodyConfig,
    relationship_rows: int = 0,
    relationship_bytes_written: int = 0,
) -> dict[str, object]:
    _validate_config(config)
    particle_count = config.particles_per_dim**3
    mesh_cell_count = config.mesh_size**3
    half_modes = config.mesh_size * config.mesh_size * (config.mesh_size // 2 + 1)
    particle_bytes = particle_count * (6 * 4 + 8)
    mesh_bytes = mesh_cell_count * 4
    complex_bytes = half_modes * 16
    peak_memory = particle_bytes * 2 + mesh_bytes * 8 + complex_bytes * 6
    raw_all = (particle_bytes + 4 * mesh_bytes) * redshift_count
    ratio = relationship_bytes_written / raw_all if raw_all > 0 else 0.0
    return {
        "particles_per_dim": config.particles_per_dim,
        "mesh_size": config.mesh_size,
        "redshift_count": redshift_count,
        "particle_count": particle_count,
        "mesh_cell_count": mesh_cell_count,
        "box_size_mpch": config.box_size_mpch,
        "estimated_peak_memory_gb": peak_memory / 1024**3,
        "raw_particle_bytes_per_snapshot": particle_bytes,
        "raw_mesh_bytes_per_snapshot": mesh_bytes,
        "raw_bytes_all_snapshots": raw_all,
        "relationship_rows": relationship_rows,
        "relationship_bytes_written": relationship_bytes_written,
        "saved_to_raw_ratio": ratio,
    }


def pm_cost_row(
    *,
    run_id: str,
    config: PMNBodyConfig,
    redshift_count: int,
    relationship_rows: int,
    relationship_bytes_written: int,
) -> dict[str, object]:
    estimate = estimate_pm_resources(
        redshift_count=redshift_count,
        config=config,
        relationship_rows=relationship_rows,
        relationship_bytes_written=relationship_bytes_written,
    )
    return {
        "run_id": run_id,
        "source": "internal_pm_nbody",
        "particles_per_dim": estimate["particles_per_dim"],
        "mesh_size": estimate["mesh_size"],
        "redshift_count": estimate["redshift_count"],
        "particle_count": estimate["particle_count"],
        "mesh_cell_count": estimate["mesh_cell_count"],
        "box_size_mpch": estimate["box_size_mpch"],
        "estimated_peak_memory_gb": estimate["estimated_peak_memory_gb"],
        "raw_particle_bytes_per_snapshot": estimate["raw_particle_bytes_per_snapshot"],
        "raw_mesh_bytes_per_snapshot": estimate["raw_mesh_bytes_per_snapshot"],
        "raw_bytes_all_snapshots": estimate["raw_bytes_all_snapshots"],
        "relationship_rows": estimate["relationship_rows"],
        "relationship_bytes_written": estimate["relationship_bytes_written"],
        "saved_to_raw_ratio": estimate["saved_to_raw_ratio"],
        "storage_policy": "relationship_tables_only",
    }


def generate_pm_history(
    *,
    run_id: str,
    redshifts: Sequence[float],
    params: CosmicCGParams,
    config: PMNBodyConfig,
    progress: Callable[[str], None] | None = None,
) -> dict[str, list[dict[str, object]]]:
    _validate_config(config)
    npx = require_numpy()
    targets = sorted([_target_state(z) for z in redshifts], key=lambda item: item["a"])
    initial = _target_state(config.z_initial)
    if targets and targets[0]["a"] < initial["a"]:
        raise ValueError("all output redshifts must be later than z_initial")
    setup = _spectral_setup(config)
    setup["params"] = params
    particle_mass = _particle_mass(config, params)
    state = _initial_state(config, setup, npx)
    outputs: dict[str, list[dict[str, object]]] = {"snapshot": [], "power": [], "halo": [], "tree": [], "transfer": []}
    previous_halos: list[dict[str, object]] = []
    previous_meta: dict[str, float] | None = None
    current_a = float(initial["a"])
    for step, target in enumerate(targets):
        if progress is not None:
            progress(f"pm snapshot z={target['z']}")
        current_a = _advance_to(state, current_a, float(target["a"]), config, setup)
        snapshot = _snapshot_outputs(
            run_id=run_id,
            step=step,
            target=target,
            state=state,
            setup=setup,
            params=params,
            particle_mass=particle_mass,
            config=config,
            previous_halos=previous_halos,
            previous_meta=previous_meta,
        )
        outputs["snapshot"].append(snapshot["snapshot_row"])
        outputs["power"].extend(snapshot["power_rows"])
        outputs["halo"].extend(snapshot["halo_rows"])
        outputs["tree"].extend(snapshot["tree_rows"])
        outputs["transfer"].extend(snapshot["transfer_rows"])
        previous_halos = snapshot["internal_halos"]
        previous_meta = {"z": float(target["z"]), "a": float(target["a"])}
    return outputs


def _validate_config(config: PMNBodyConfig) -> None:
    if config.dynamics not in {"cosmological", "comparison"}:
        raise ValueError("Unknown trajectory dynamics")
    if config.particles_per_dim < 4:
        raise ValueError("particles_per_dim must be at least 4")
    if config.mesh_size < config.particles_per_dim:
        raise ValueError("mesh_size must be at least particles_per_dim")
    if config.mesh_size % 2:
        raise ValueError("mesh_size must be even")
    if config.box_size_mpch <= 0:
        raise ValueError("box_size_mpch must be positive")
    if config.z_initial <= -1.0:
        raise ValueError("z_initial must be greater than -1")
    if config.steps_per_interval <= 0 or config.mode_bins <= 0:
        raise ValueError("step and mode counts must be positive")
    if config.fof_min_particles <= 0 or config.max_halos_per_snapshot <= 0:
        raise ValueError("halo limits must be positive")


def _target_state(z: float) -> dict[str, float]:
    if z <= -1.0:
        raise ValueError("redshift must be greater than -1")
    return {"z": float(z), "a": 1.0 / (1.0 + float(z))}


def _particle_mass(config: PMNBodyConfig, params: CosmicCGParams) -> float:
    rho_m = RHO_CRIT_MSUNH_MPCH3 * omega_m0(params)
    return rho_m * config.box_size_mpch**3 / config.particles_per_dim**3


def _spectral_setup(config: PMNBodyConfig) -> dict[str, object]:
    npx = require_numpy()
    mesh = config.mesh_size
    dx = config.box_size_mpch / mesh
    kx = (2.0 * math.pi * npx.fft.fftfreq(mesh, d=dx)).astype(npx.float32)
    ky = (2.0 * math.pi * npx.fft.fftfreq(mesh, d=dx)).astype(npx.float32)
    kz = (2.0 * math.pi * npx.fft.rfftfreq(mesh, d=dx)).astype(npx.float32)
    k2 = (
        kx[:, None, None] * kx[:, None, None]
        + ky[None, :, None] * ky[None, :, None]
        + kz[None, None, :] * kz[None, None, :]
    ).astype(npx.float32)
    k2[0, 0, 0] = npx.inf
    k_mag = npx.sqrt(k2, dtype=npx.float32)
    k_mag[0, 0, 0] = 0.0
    k_resolved = math.pi / dx
    edges = npx.geomspace(2.0 * math.pi / config.box_size_mpch, k_resolved, config.mode_bins + 1, dtype=npx.float32)
    bin_index = _mode_bin_index(k_mag, edges)
    return {
        "dx": dx,
        "kx": kx,
        "ky": ky,
        "kz": kz,
        "k2": k2,
        "k_mag": k_mag,
        "bin_edges": edges,
        "bin_index": bin_index,
        "bin_counts": npx.bincount(bin_index.ravel(), minlength=config.mode_bins + 1).astype(npx.int64),
    }


def _mode_bin_index(k_mag, edges):
    npx = require_numpy()
    valid = k_mag > 0.0
    bins = npx.zeros(k_mag.shape, dtype=npx.uint16)
    if npx.any(valid):
        raw = npx.searchsorted(edges, k_mag, side="right").astype(npx.int32)
        raw = npx.clip(raw, 1, len(edges) - 1)
        bins[valid] = raw[valid].astype(npx.uint16)
    return bins


def _initial_state(config: PMNBodyConfig, setup: Mapping[str, object], npx) -> dict[str, object]:
    rng = npx.random.default_rng(config.seed)
    n = config.particles_per_dim
    dxp = config.box_size_mpch / n
    grid = npx.indices((n, n, n), dtype=npx.float32).reshape(3, -1).T
    positions = ((grid + 0.5) * dxp).astype(npx.float32)
    displacement_grid = _initial_displacement_grid(config, setup, rng)
    displacement = _sample_vector_field(displacement_grid, positions, config)
    target_rms = config.initial_displacement_fraction * dxp
    rms = float(npx.sqrt(npx.mean(npx.sum(displacement * displacement, axis=1), dtype=npx.float64)))
    displacement_scale = npx.float32(1.0)
    if rms > 0:
        displacement_scale = npx.float32(target_rms / rms)
        displacement *= displacement_scale
    positions = (positions + displacement) % npx.float32(config.box_size_mpch)
    velocities = (config.initial_velocity_factor * displacement).astype(npx.float32)
    if config.dynamics == "cosmological":
        a0 = 1.0 / (1.0 + config.z_initial)
        # u=dx/da.  A growing displacement has u=f D/a; the comparison
        # model's independently assigned velocity multiplier is not used.
        growth_rate = growing_rate(a0, setup["params"])
        velocities = (growth_rate / a0 * displacement).astype(npx.float64)
    ids = npx.arange(positions.shape[0], dtype=npx.int64)
    result={"positions": positions.astype(npx.float32), "velocities": velocities, "ids": ids}
    if setup.get("capture_initial_sheet",False):
        result["_initial_sheet"]={"displacement_grid":displacement_grid,
                                   "displacement_scale":float(displacement_scale),
                                   "growth_rate":growth_rate if config.dynamics=="cosmological" else None}
    return result


def _initial_displacement_grid(config: PMNBodyConfig, setup: Mapping[str, object], rng):
    npx = require_numpy()
    mesh = config.mesh_size
    white = rng.normal(0.0, 1.0, size=(mesh, mesh, mesh)).astype(npx.float32)
    white -= npx.mean(white, dtype=npx.float64).astype(npx.float32)
    modes = npx.fft.rfftn(white)
    k_mag = setup["k_mag"]
    smooth_k = 0.35 * math.pi / float(setup["dx"])
    modes *= npx.exp(-0.5 * npx.square(k_mag / smooth_k)).astype(modes.dtype)
    k2 = setup["k2"]
    fields = []
    for axis in ("kx", "ky", "kz"):
        component_k = -1j * _axis_array(setup, axis) * modes / k2
        component_k[0, 0, 0] = 0.0
        fields.append(npx.fft.irfftn(component_k, s=(mesh, mesh, mesh), axes=(0, 1, 2)).astype(npx.float32))
    del white, modes
    return tuple(fields)


def _axis_array(setup: Mapping[str, object], axis: str):
    if axis == "kx":
        return setup["kx"][:, None, None]
    if axis == "ky":
        return setup["ky"][None, :, None]
    return setup["kz"][None, None, :]


def _advance_to(
    state: dict[str, object],
    current_a: float,
    target_a: float,
    config: PMNBodyConfig,
    setup: Mapping[str, object],
) -> float:
    if target_a < current_a:
        raise ValueError("target scale factor must not decrease")
    if target_a == current_a:
        return current_a
    if config.dynamics == "cosmological":
        return _advance_cosmological(state, current_a, target_a, config, setup)
    steps = max(1, config.steps_per_interval)
    da = (target_a - current_a) / steps
    a = current_a
    for _ in range(steps):
        midpoint = a + 0.5 * da
        acceleration = _particle_acceleration(state["positions"], config, setup)
        state["velocities"] += _pm_scale(midpoint, config) * _npx_float(0.5 * da) * acceleration
        state["positions"] = (state["positions"] + _npx_float(da) * state["velocities"]) % _npx_float(config.box_size_mpch)
        acceleration = _particle_acceleration(state["positions"], config, setup)
        state["velocities"] += _pm_scale(midpoint, config) * _npx_float(0.5 * da) * acceleration
        a += da
    return target_a


def _advance_cosmological(state, current_a, target_a, config, setup):
    """Physical a-clock, p=a^3 E u; both drag and force are retained."""
    npx = require_numpy()
    params = setup["params"]
    step = (target_a - current_a) / config.steps_per_interval
    state["positions"] = npx.asarray(state["positions"], dtype=npx.float64)
    p = current_a**3 * expansion(current_a, params) * state["velocities"].astype(npx.float64)
    acceleration = _particle_acceleration(state["positions"], config, setup)
    for n in range(config.steps_per_interval):
        a, end = current_a + n*step, current_a + (n+1)*step
        middle = (a+end)/2
        _, k1 = drift_kick_integrals(a, middle, params)
        drift, _ = drift_kick_integrals(a, end, params)
        _, k2 = drift_kick_integrals(middle, end, params)
        p += k1 * acceleration
        state["positions"] = (state["positions"] + drift*p) % config.box_size_mpch
        acceleration = _particle_acceleration(state["positions"], config, setup)
        p += k2 * acceleration
    state["velocities"] = p / (target_a**3 * expansion(target_a, params))
    return target_a


def _npx_float(value: float):
    return require_numpy().float32(value)


def _pm_scale(a: float, config: PMNBodyConfig):
    return _npx_float(config.force_strength / max(a, 1.0e-4))


def _particle_acceleration(positions, config: PMNBodyConfig, setup: Mapping[str, object]):
    delta = _cic_delta(positions, config)
    acceleration_grid, _phi = _grid_acceleration(delta, config, setup, return_potential=False)
    return _sample_vector_field(acceleration_grid, positions, config)


def _cic_delta(positions, config: PMNBodyConfig):
    npx = require_numpy()
    mesh = config.mesh_size
    dx = config.box_size_mpch / mesh
    coord = positions / npx.float32(dx)
    base = npx.floor(coord).astype(npx.int32)
    frac = coord - base.astype(npx.float32)
    density = npx.zeros((mesh, mesh, mesh), dtype=npx.float32)
    for sx in (0, 1):
        wx = frac[:, 0] if sx else 1.0 - frac[:, 0]
        ix = (base[:, 0] + sx) % mesh
        for sy in (0, 1):
            wy = frac[:, 1] if sy else 1.0 - frac[:, 1]
            iy = (base[:, 1] + sy) % mesh
            for sz in (0, 1):
                wz = frac[:, 2] if sz else 1.0 - frac[:, 2]
                iz = (base[:, 2] + sz) % mesh
                npx.add.at(density, (ix, iy, iz), (wx * wy * wz).astype(npx.float32))
    mean = positions.shape[0] / float(mesh**3)
    density /= npx.float32(mean)
    density -= npx.float32(1.0)
    return density


def _grid_acceleration(delta, config: PMNBodyConfig, setup: Mapping[str, object], *, return_potential=True):
    npx = require_numpy()
    mesh = config.mesh_size
    delta_k = npx.fft.rfftn(delta)
    phi_k = -delta_k / setup["k2"]
    phi_k[0, 0, 0] = 0.0
    components = []
    for axis in ("kx", "ky", "kz"):
        acc_k = -1j * _axis_array(setup, axis) * phi_k
        acc_k[0, 0, 0] = 0.0
        components.append(npx.fft.irfftn(acc_k, s=(mesh, mesh, mesh), axes=(0, 1, 2)).astype(npx.float32))
    phi = (npx.fft.irfftn(phi_k, s=(mesh, mesh, mesh), axes=(0, 1, 2)).astype(npx.float32)
           if return_potential else None)
    return tuple(components), phi


def _sample_vector_field(vector_grid, positions, config: PMNBodyConfig):
    npx = require_numpy()
    return npx.stack([_sample_scalar_field(component, positions, config) for component in vector_grid], axis=1).astype(npx.float32)


def _sample_scalar_field(field, positions, config: PMNBodyConfig):
    npx = require_numpy()
    mesh = config.mesh_size
    dx = config.box_size_mpch / mesh
    coord = positions / npx.float32(dx)
    base = npx.floor(coord).astype(npx.int32)
    frac = coord - base.astype(npx.float32)
    values = npx.zeros(positions.shape[0], dtype=npx.float32)
    for sx in (0, 1):
        wx = frac[:, 0] if sx else 1.0 - frac[:, 0]
        ix = (base[:, 0] + sx) % mesh
        for sy in (0, 1):
            wy = frac[:, 1] if sy else 1.0 - frac[:, 1]
            iy = (base[:, 1] + sy) % mesh
            for sz in (0, 1):
                wz = frac[:, 2] if sz else 1.0 - frac[:, 2]
                iz = (base[:, 2] + sz) % mesh
                values += (wx * wy * wz).astype(npx.float32) * field[ix, iy, iz]
    return values


def _snapshot_outputs(
    *,
    run_id: str,
    step: int,
    target: Mapping[str, float],
    state: Mapping[str, object],
    setup: Mapping[str, object],
    params: CosmicCGParams,
    particle_mass: float,
    config: PMNBodyConfig,
    previous_halos: Sequence[Mapping[str, object]],
    previous_meta: Mapping[str, float] | None,
) -> dict[str, object]:
    npx = require_numpy()
    positions = state["positions"]
    velocities = state["velocities"]
    delta = _cic_delta(positions, config)
    _acceleration_grid, phi = _grid_acceleration(delta, config, setup)
    internal_halos, halo_rows, tree_rows = _halo_rows(
        run_id=run_id,
        step=step,
        target=target,
        positions=positions,
        velocities=velocities,
        ids=state["ids"],
        delta=delta,
        params=params,
        particle_mass=particle_mass,
        config=config,
        previous_halos=previous_halos,
        previous_meta=previous_meta,
    )
    transfer_rows = _transfer_rows(run_id, target, internal_halos, config)
    power_rows = _power_rows(run_id, target, delta, setup, config)
    snapshot_row = {
        "run_id": run_id,
        "source": "internal_pm_nbody",
        "snapshot": f"a{_tag(target['a'])}",
        "step": step,
        "z": target["z"],
        "a": target["a"],
        "particle_count": positions.shape[0],
        "mesh_size": config.mesh_size,
        "mesh_cell_count": config.mesh_size**3,
        "box_size_mpch": config.box_size_mpch,
        "dx_mesh_mpch": setup["dx"],
        "particle_mass_msunh": particle_mass,
        "mean_delta": float(npx.mean(delta, dtype=npx.float64)),
        "rms_delta": float(npx.std(delta, dtype=npx.float64)),
        "min_delta": float(npx.min(delta)),
        "max_delta": float(npx.max(delta)),
        "kinetic_proxy": float(0.5 * npx.mean(npx.sum(velocities * velocities, axis=1), dtype=npx.float64)),
        "potential_proxy": float(0.5 * npx.mean(delta * phi, dtype=npx.float64)),
        "total_energy_proxy": float(0.5 * npx.mean(npx.sum(velocities * velocities, axis=1), dtype=npx.float64) + 0.5 * npx.mean(delta * phi, dtype=npx.float64)),
        "mass_conservation_error": abs(float(npx.sum(delta, dtype=npx.float64))) / max(1.0, config.mesh_size**3),
        "momentum_norm_proxy": float(npx.sqrt(npx.sum(npx.square(npx.mean(velocities, axis=0, dtype=npx.float64))))),
        "halo_count": len(halo_rows),
        "tree_link_count": len(tree_rows),
        "transfer_count": len(transfer_rows),
        "quality_flag": "particle_mesh_trajectory",
    }
    return {
        "snapshot_row": snapshot_row,
        "power_rows": power_rows,
        "halo_rows": halo_rows,
        "tree_rows": tree_rows,
        "transfer_rows": transfer_rows,
        "internal_halos": internal_halos,
    }


def _halo_rows(
    *,
    run_id: str,
    step: int,
    target: Mapping[str, float],
    positions,
    velocities,
    ids,
    delta,
    params: CosmicCGParams,
    particle_mass: float,
    config: PMNBodyConfig,
    previous_halos: Sequence[Mapping[str, object]],
    previous_meta: Mapping[str, float] | None,
) -> tuple[list[dict[str, object]], list[dict[str, object]], list[dict[str, object]]]:
    components = _density_components(delta, positions, config)
    internal_halos = []
    halo_rows = []
    tree_rows = []
    age = cosmic_age_gyr(float(target["z"]), params=params)
    previous_age = cosmic_age_gyr(float(previous_meta["z"]), params=params) if previous_meta is not None else None
    for rank, component in enumerate(components):
        particle_indices = component["particle_indices"]
        if len(particle_indices) < config.fof_min_particles:
            continue
        particle_ids = ids[particle_indices]
        halo_positions = positions[particle_indices]
        halo_velocities = velocities[particle_indices]
        center = _periodic_center(halo_positions, config.box_size_mpch)
        velocity = _mean_vector(halo_velocities)
        radius_mpch = _radius_mpch(halo_positions, center, config.box_size_mpch)
        spin = _spin_proxy(halo_positions, halo_velocities, center, velocity, config.box_size_mpch)
        mass = len(particle_indices) * particle_mass
        parent = _best_parent(particle_ids, previous_halos)
        delta_time = max(1.0e-4, age - previous_age) if previous_age is not None else max(age, 1.0e-4)
        mass_from = float(parent["mass_msunh"]) if parent is not None else max(0.0, mass - particle_mass * len(particle_ids))
        delta_mass = mass - mass_from
        fate = halo_fate(
            halo_id=f"pm_{step:03d}_{rank:03d}",
            mass_msun_h=max(mass, 1.0),
            radius_kpc_h=max(radius_mpch * 1000.0, 1.0e-6),
            delta_mass_msun_h=delta_mass if delta_mass != 0.0 else 1.0e-8 * mass,
            delta_time_gyr=delta_time,
        )
        b_h = 0.0 if not math.isfinite(fate.b_h) else fate.b_h
        halo_id = f"pm_{step:03d}_{rank:03d}"
        parent_id = "" if parent is None else str(parent["halo_id"])
        match_fraction = 0.0 if parent is None else float(parent["overlap_fraction"])
        internal = {
            "halo_id": halo_id,
            "rank": rank,
            "particle_ids": particle_ids,
            "particle_set": set(int(value) for value in particle_ids.tolist()),
            "particle_count": int(len(particle_ids)),
            "mass_msunh": mass,
            "radius_kpch": radius_mpch * 1000.0,
            "position": center,
            "velocity": velocity,
            "spin": spin,
            "b_h": b_h,
            "fate_class": fate.class_name,
        }
        internal_halos.append(internal)
        halo_rows.append(
            {
                "run_id": run_id,
                "source": "internal_pm_nbody",
                "snapshot": f"a{_tag(target['a'])}",
                "z": target["z"],
                "a": target["a"],
                "halo_id": halo_id,
                "rank": rank,
                "particle_count": int(len(particle_ids)),
                "mass_msunh": mass,
                "radius_kpch": radius_mpch * 1000.0,
                "x_mpch": center[0],
                "y_mpch": center[1],
                "z_mpch": center[2],
                "vx_proxy": velocity[0],
                "vy_proxy": velocity[1],
                "vz_proxy": velocity[2],
                "spin_x": spin[0],
                "spin_y": spin[1],
                "spin_z": spin[2],
                "parent_id": parent_id,
                "descendant_id": halo_id,
                "match_particle_fraction": match_fraction,
                "b_h": b_h,
                "t_dyn_gyr": fate.t_dyn_gyr,
                "t_growth_gyr": "" if fate.t_growth_gyr is None else fate.t_growth_gyr,
                "fate_class": fate.class_name,
                "quality_flag": "particle_id_tracked",
            }
        )
        if parent is not None and previous_meta is not None:
            tree_rows.append(
                {
                    "run_id": run_id,
                    "source": "internal_pm_nbody",
                    "snapshot_from": f"a{_tag(previous_meta['a'])}",
                    "snapshot_to": f"a{_tag(target['a'])}",
                    "z_from": previous_meta["z"],
                    "z_to": target["z"],
                    "a_from": previous_meta["a"],
                    "a_to": target["a"],
                    "parent_id": parent_id,
                    "descendant_id": halo_id,
                    "shared_particle_count": int(parent["shared_particle_count"]),
                    "parent_particle_count": int(parent["parent_particle_count"]),
                    "descendant_particle_count": int(len(particle_ids)),
                    "overlap_fraction": match_fraction,
                    "mass_from_msunh": mass_from,
                    "mass_to_msunh": mass,
                    "delta_mass_msunh": delta_mass,
                    "delta_time_gyr": delta_time,
                    "b_h": b_h,
                    "fate_class": fate.class_name,
                    "quality_flag": "particle_id_overlap",
                }
            )
    return internal_halos, halo_rows, tree_rows


def _density_components(delta, positions, config: PMNBodyConfig) -> list[dict[str, object]]:
    npx = require_numpy()
    threshold = max(config.density_threshold_delta, float(npx.mean(delta, dtype=npx.float64) + config.density_threshold_sigma * npx.std(delta, dtype=npx.float64)))
    active = npx.flatnonzero(delta > npx.float32(threshold))
    if active.size == 0:
        return []
    values = delta.ravel()[active]
    active = active[npx.argsort(values)[::-1]]
    active_set = set(int(value) for value in active.tolist())
    visited: set[int] = set()
    cell_particles = _cell_particle_map(positions, config)
    components = []
    for seed in active.tolist():
        seed = int(seed)
        if seed in visited:
            continue
        queue = [seed]
        visited.add(seed)
        cells = []
        particle_indices = []
        while queue:
            cell = queue.pop()
            cells.append(cell)
            particle_indices.extend(cell_particles.get(cell, ()))
            for neighbor in _cell_neighbors(cell, config.mesh_size):
                if neighbor in active_set and neighbor not in visited:
                    visited.add(neighbor)
                    queue.append(neighbor)
        if len(particle_indices) >= config.fof_min_particles:
            components.append({"cells": cells, "particle_indices": sorted(set(particle_indices))})
        if len(components) >= config.max_halos_per_snapshot:
            break
    components.sort(key=lambda item: len(item["particle_indices"]), reverse=True)
    return components[: config.max_halos_per_snapshot]


def _cell_particle_map(positions, config: PMNBodyConfig) -> dict[int, list[int]]:
    npx = require_numpy()
    mesh = config.mesh_size
    dx = config.box_size_mpch / mesh
    cell = npx.floor(positions / npx.float32(dx)).astype(npx.int64) % mesh
    flat = (cell[:, 0] * mesh + cell[:, 1]) * mesh + cell[:, 2]
    mapping: dict[int, list[int]] = {}
    for index, value in enumerate(flat.tolist()):
        mapping.setdefault(int(value), []).append(index)
    return mapping


def _cell_neighbors(flat: int, mesh: int) -> tuple[int, ...]:
    ix = flat // (mesh * mesh)
    rem = flat - ix * mesh * mesh
    iy = rem // mesh
    iz = rem - iy * mesh
    values = []
    for dx, dy, dz in ((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)):
        nx = (ix + dx) % mesh
        ny = (iy + dy) % mesh
        nz = (iz + dz) % mesh
        values.append((nx * mesh + ny) * mesh + nz)
    return tuple(values)


def _periodic_center(points, box_size: float) -> tuple[float, float, float]:
    npx = require_numpy()
    center = []
    for axis in range(3):
        angle = 2.0 * math.pi * points[:, axis] / box_size
        mean_sin = float(npx.mean(npx.sin(angle), dtype=npx.float64))
        mean_cos = float(npx.mean(npx.cos(angle), dtype=npx.float64))
        value = math.atan2(mean_sin, mean_cos) * box_size / (2.0 * math.pi)
        center.append(value % box_size)
    return (center[0], center[1], center[2])


def _mean_vector(values) -> tuple[float, float, float]:
    npx = require_numpy()
    mean = npx.mean(values, axis=0, dtype=npx.float64)
    return (float(mean[0]), float(mean[1]), float(mean[2]))


def _periodic_offsets(points, center: Sequence[float], box_size: float):
    npx = require_numpy()
    offsets = points - npx.asarray(center, dtype=npx.float32)
    offsets -= npx.rint(offsets / npx.float32(box_size)) * npx.float32(box_size)
    return offsets


def _radius_mpch(points, center: Sequence[float], box_size: float) -> float:
    npx = require_numpy()
    offsets = _periodic_offsets(points, center, box_size)
    radius2 = npx.sum(offsets * offsets, axis=1)
    return float(npx.sqrt(npx.max(radius2))) if radius2.size else 0.0


def _spin_proxy(points, velocities, center: Sequence[float], mean_velocity: Sequence[float], box_size: float) -> tuple[float, float, float]:
    npx = require_numpy()
    offsets = _periodic_offsets(points, center, box_size)
    relative_velocity = velocities - npx.asarray(mean_velocity, dtype=npx.float32)
    spin = npx.sum(npx.cross(offsets, relative_velocity), axis=0, dtype=npx.float64)
    return (float(spin[0]), float(spin[1]), float(spin[2]))


def _best_parent(particle_ids, previous_halos: Sequence[Mapping[str, object]]) -> dict[str, object] | None:
    if not previous_halos:
        return None
    current = set(int(value) for value in particle_ids.tolist())
    best: dict[str, object] | None = None
    for parent in previous_halos:
        shared = len(current & parent["particle_set"])
        if shared <= 0:
            continue
        fraction = shared / max(1, len(current))
        candidate = {
            "halo_id": parent["halo_id"],
            "mass_msunh": parent["mass_msunh"],
            "shared_particle_count": shared,
            "parent_particle_count": parent["particle_count"],
            "overlap_fraction": fraction,
        }
        if best is None or fraction > float(best["overlap_fraction"]):
            best = candidate
    return best


def _transfer_rows(run_id: str, target: Mapping[str, float], halos: Sequence[Mapping[str, object]], config: PMNBodyConfig) -> list[dict[str, object]]:
    rows = []
    for source in halos:
        distances = []
        for target_halo in halos:
            if source["halo_id"] == target_halo["halo_id"]:
                continue
            distance = _periodic_distance(source["position"], target_halo["position"], config.box_size_mpch)
            distances.append((distance, target_halo))
        distances.sort(key=lambda item: item[0])
        for distance, target_halo in distances[: config.neighbor_links]:
            rows.append(_transfer_row(run_id, target, source, target_halo, distance, config.box_size_mpch))
    return rows


def _transfer_row(
    run_id: str,
    target: Mapping[str, float],
    source: Mapping[str, object],
    target_halo: Mapping[str, object],
    distance: float,
    box_size: float,
) -> dict[str, object]:
    mass_ratio = float(source["mass_msunh"]) / max(float(target_halo["mass_msunh"]), 1.0)
    separation = _separation_vector(source["position"], target_halo["position"], distance, box_size)
    relative_velocity = tuple(float(target_halo["velocity"][i]) - float(source["velocity"][i]) for i in range(3))
    torque = _cross(separation, relative_velocity)
    spin = target_halo["spin"]
    a_pq = mass_ratio / max(distance, 1.0e-12) ** 3
    g_pq = _safe_alignment(spin, torque)
    return {
        "run_id": run_id,
        "source": "internal_pm_nbody",
        "snapshot": f"a{_tag(target['a'])}",
        "z": target["z"],
        "a": target["a"],
        "source_id": source["halo_id"],
        "target_id": target_halo["halo_id"],
        "rank_source": source["rank"],
        "rank_target": target_halo["rank"],
        "separation_mpch": distance,
        "mass_ratio": mass_ratio,
        "a_pq": a_pq,
        "g_pq": g_pq,
        "tau_x": torque[0] * a_pq,
        "tau_y": torque[1] * a_pq,
        "tau_z": torque[2] * a_pq,
        "j_x": spin[0],
        "j_y": spin[1],
        "j_z": spin[2],
        "source_channel": "pm_particle_neighbour_transfer",
        "quality_flag": "particle_mesh_relation",
    }


def _separation_vector(
    source: Sequence[float],
    target: Sequence[float],
    distance: float,
    box_size: float,
) -> tuple[float, float, float]:
    if distance <= 0:
        return (0.0, 0.0, 0.0)
    values = []
    for index in range(3):
        delta = float(target[index]) - float(source[index])
        delta -= round(delta / box_size) * box_size
        values.append(delta)
    return (values[0], values[1], values[2])


def _periodic_distance(source: Sequence[float], target: Sequence[float], box_size: float) -> float:
    total = 0.0
    for left, right in zip(source, target):
        delta = abs(float(left) - float(right))
        delta = min(delta, box_size - delta)
        total += delta * delta
    return math.sqrt(total)


def _power_rows(run_id: str, target: Mapping[str, float], delta, setup: Mapping[str, object], config: PMNBodyConfig) -> list[dict[str, object]]:
    npx = require_numpy()
    modes = npx.fft.rfftn(delta)
    power = npx.abs(modes) ** 2 / float(config.mesh_size**6)
    bin_index = setup["bin_index"]
    power_by_bin = npx.bincount(bin_index.ravel(), weights=power.ravel(), minlength=config.mode_bins + 1)
    total_power = float(npx.sum(power_by_bin[1:], dtype=npx.float64))
    edges = setup["bin_edges"]
    counts = setup["bin_counts"]
    cumulative = 0.0
    rows = []
    for index in range(1, config.mode_bins + 1):
        fraction = float(power_by_bin[index] / total_power) if total_power > 0 else 0.0
        cumulative += fraction
        k_min = float(edges[index - 1])
        k_max = float(edges[index])
        rows.append(
            {
                "run_id": run_id,
                "source": "internal_pm_nbody",
                "snapshot": f"a{_tag(target['a'])}",
                "z": target["z"],
                "a": target["a"],
                "mode_bin": index,
                "k_min_hmpc": k_min,
                "k_max_hmpc": k_max,
                "k_mid_hmpc": math.sqrt(k_min * k_max),
                "mode_count": int(counts[index]),
                "power_proxy": float(power_by_bin[index]),
                "power_fraction": fraction,
                "cumulative_power_fraction": cumulative,
                "quality_flag": "pm_density_power",
            }
        )
    return rows


def _cross(left: Sequence[float], right: Sequence[float]) -> tuple[float, float, float]:
    return (
        float(left[1]) * float(right[2]) - float(left[2]) * float(right[1]),
        float(left[2]) * float(right[0]) - float(left[0]) * float(right[2]),
        float(left[0]) * float(right[1]) - float(left[1]) * float(right[0]),
    )


def _safe_alignment(left: Sequence[float], right: Sequence[float]) -> float:
    try:
        return alignment_cosine(left, right)
    except ValueError:
        return 0.0


def _tag(value: float) -> str:
    text = f"{value:.6g}".replace("-", "m").replace(".", "p")
    return text
