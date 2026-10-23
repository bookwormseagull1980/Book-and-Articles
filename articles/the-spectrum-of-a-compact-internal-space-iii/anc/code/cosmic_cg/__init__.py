"""Finite-window diagnostics for cosmological structure formation."""

from .entropy_audit import EntropyAudit, entropy_deficit, window_entropy
from .knl import (
    KNLResult,
    band_power_fraction,
    dimensionless_power,
    find_knl,
    tail_power_fraction,
    tail_power_moment_fraction,
)
from .params import CosmicCGParams, load_model_params
from .tail_ratio import TailRatioResult, halo_tail_ratio, virial_velocity
from .field3d import Field3DConfig, available_physical_memory_bytes, estimate_field_resources, generate_field_history
from .pm_nbody import (
    PMNBodyConfig,
    available_physical_memory_bytes as available_pm_physical_memory_bytes,
    estimate_pm_resources,
    generate_pm_history,
)

__all__ = [
    "EntropyAudit",
    "Field3DConfig",
    "KNLResult",
    "PMNBodyConfig",
    "CosmicCGParams",
    "TailRatioResult",
    "available_physical_memory_bytes",
    "available_pm_physical_memory_bytes",
    "band_power_fraction",
    "dimensionless_power",
    "entropy_deficit",
    "estimate_field_resources",
    "estimate_pm_resources",
    "find_knl",
    "generate_field_history",
    "generate_pm_history",
    "halo_tail_ratio",
    "load_model_params",
    "tail_power_fraction",
    "tail_power_moment_fraction",
    "virial_velocity",
    "window_entropy",
]

__version__ = "0.1.0"
