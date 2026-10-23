"""Declared constants for the generated finite-window history."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import hashlib
import json
import math
import os
from pathlib import Path
from typing import Any, Mapping


DEFAULT_PARAM_PATH_ENV = "COSMIC_CG_PARAMS"


@dataclass(frozen=True)
class CosmicCGParams:
    """Parameter snapshot required by the finite-window calculation."""

    m_p_gev: float
    gamma_m_integral: float
    h0_gev: float
    omega_lambda: float
    omega_b: float
    omega_sigma: float
    rho_lambda_gev4: float
    m_nu1_ev: float | None = None
    kl: float | None = None
    source_path: str | None = None
    source_hash: str | None = None
    source_run_id: str | None = None
    source_reference: str | None = None
    snapshot_hash: str | None = None

    @property
    def s_abs(self) -> float:
        return self.gamma_m_integral

    @property
    def sigma0_over_sigmac(self) -> float:
        return math.sqrt(math.pi) * math.exp(-self.gamma_m_integral)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


DEFAULT_PARAMS = CosmicCGParams(
    m_p_gev=2.435323595526305e18,
    gamma_m_integral=139.2537061093592,
    h0_gev=1.4388498142964775e-42,
    omega_lambda=0.6850442965140158,
    omega_b=0.049408192038365814,
    omega_sigma=0.2655475114476184,
    rho_lambda_gev4=2.5233889268202705e-47,
    m_nu1_ev=0.002607415413287687,
    kl=2.4935343325226915,
    source_path=None,
    source_hash=None,
)


KEY_ALIASES = {
    "m_p_gev": ["M_P", "M_P_GeV", "M_P_GEV"],
    "gamma_m_integral": ["entropy_integral", "ir_flow_int_gamma", "int_gamma_M", "gamma_M_integral", "gamma_M_int"],
    "h0_gev": ["H0_GEV", "H0_GeV"],
    "omega_lambda": ["Omega_Lambda", "Omega_L"],
    "omega_b": ["Omega_b"],
    "omega_sigma": ["Omega_Sigma", "Omega_DM"],
    "rho_lambda_gev4": ["rho_Lambda", "rho_Lambda_GeV4"],
    "m_nu1_ev": ["m_nu1", "m_nu1_eV"],
    "kl": ["L_Gg", "kL", "kappa_L"],
}


def load_model_params(path: str | os.PathLike[str] | None = None) -> CosmicCGParams:
    """Load the declared parameter subset from a JSON file.

    Accept a V4 parameter store or this package's flat snapshot.  Explicit
    input is complete: a missing required quantity is never silently replaced.
    The packaged current snapshot binds default production to a V4 release.
    """

    raw_path = path if path is not None else os.environ.get(DEFAULT_PARAM_PATH_ENV)
    if not raw_path:
        packaged = Path(__file__).resolve().parents[2] / "tables" / "v4_parameter_snapshot.json"
        if packaged.exists():
            raw_path = packaged
        else:
            raise FileNotFoundError("The formal V4 parameter snapshot is missing")
    resolved = Path(raw_path).expanduser()
    if not resolved.exists():
        raise FileNotFoundError(f"parameter file not found: {resolved}")

    raw = json.loads(resolved.read_text(encoding="utf-8"))
    values = _extract_values(raw)
    data: dict[str, Any] = {}
    for field, aliases in KEY_ALIASES.items():
        matches = [values[key] for key in [field, *aliases] if key in values and values[key] is not None]
        if any(isinstance(value, bool) for value in matches):
            raise ValueError(f"A boolean is not a physical parameter: {field}")
        if matches and any(not math.isclose(float(v), float(matches[0]), rel_tol=1e-12, abs_tol=0.0) for v in matches[1:]):
            raise ValueError(f"Conflicting definitions for {field}")
        value = matches[0] if matches else None
        if value is None:
            if field in {"m_nu1_ev", "kl"}:
                data[field] = None
            else:
                raise ValueError(f"Required parameter {field} is missing")
        else:
            data[field] = float(value)
    nested = isinstance(raw.get("parameters"), Mapping)
    data["source_path"] = None
    data["snapshot_hash"] = sha256_file(resolved)
    data["source_hash"] = data["snapshot_hash"] if nested else raw.get("source_hash")
    records = raw.get("parameters", {})
    used = {key for field, aliases in KEY_ALIASES.items() for key in [field, *aliases] if key in records}
    runs = {records[key].get("run_id") for key in used if isinstance(records[key], Mapping)} - {None}
    if len(runs) > 1:
        raise ValueError("Parameters combine more than one V4 release")
    data["source_run_id"] = next(iter(runs), raw.get("source_run_id"))
    data["source_reference"] = raw.get("source_reference", "Coarse-Graining Genesis Framework V4.0")
    result = CosmicCGParams(**data)
    validate_model_params(result)
    return result


def validate_model_params(params: CosmicCGParams) -> None:
    for field in KEY_ALIASES:
        value = getattr(params, field)
        if value is None:
            continue
        if not math.isfinite(value) or value < 0 or (field not in {"omega_b", "omega_sigma", "omega_lambda", "rho_lambda_gev4", "m_nu1_ev"} and value == 0):
            raise ValueError(f"Invalid physical value for {field}")
    if not math.isclose(params.omega_b + params.omega_sigma + params.omega_lambda, 1.0, abs_tol=1e-10):
        raise ValueError("This flat matter--vacuum realization requires density fractions summing to one")
    if params.omega_b + params.omega_sigma <= 0:
        raise ValueError("The matter growing-mode realization requires positive total matter density")


def write_param_snapshot(params: CosmicCGParams, output_path: str | os.PathLike[str]) -> None:
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(params.to_dict(), indent=2, sort_keys=True), encoding="utf-8")


def sha256_file(path: str | os.PathLike[str]) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def _extract_values(obj: Any, out: dict[str, Any] | None = None) -> dict[str, Any]:
    if out is None:
        out = {}
    if isinstance(obj, Mapping):
        records = obj.get("parameters", obj)
        for key, value in records.items():
            if isinstance(value, Mapping) and "value" in value:
                out[key] = value["value"]
            elif isinstance(value, (int, float)) and not isinstance(value, bool):
                out[key] = value
    return out


def _first_present(values: Mapping[str, Any], aliases: list[str]) -> Any | None:
    for key in aliases:
        if key in values:
            return values[key]
    return None
