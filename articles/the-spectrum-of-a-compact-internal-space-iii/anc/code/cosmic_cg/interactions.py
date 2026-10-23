"""Interaction-strength spectra for generated cosmic histories."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Iterable, Mapping, Sequence

from cosmic_cg.result_tables import DEFAULT_MASS_BANDS


INTERACTION_STRENGTH_COLUMNS = (
    "run_id",
    "source",
    "snapshot",
    "z",
    "a",
    "band_id",
    "channel_id",
    "family",
    "carrier",
    "target",
    "mediator",
    "source_term",
    "sink_term",
    "net_term",
    "raw_strength",
    "chi",
    "state_factor",
    "scale_factor",
    "carrier_factor",
    "output_table",
    "physical_reading",
)

INTERACTION_BUDGET_COLUMNS = (
    "run_id",
    "source",
    "snapshot",
    "z",
    "a",
    "band_id",
    "carrier",
    "source_total",
    "sink_total",
    "net_total",
    "chi_total",
    "n_channels",
    "quality_flag",
)

INTERACTION_DOMINANCE_COLUMNS = (
    "run_id",
    "source",
    "snapshot",
    "z",
    "a",
    "band_id",
    "dominant_channel",
    "dominant_family",
    "dominant_chi",
    "second_channel",
    "second_family",
    "second_chi",
    "dominance_gap",
    "chi_sum",
    "quality_flag",
)

INTERACTION_CROSSING_COLUMNS = (
    "run_id",
    "source",
    "band_id",
    "left_snapshot",
    "right_snapshot",
    "z_left",
    "z_right",
    "a_left",
    "a_right",
    "from_channel",
    "to_channel",
    "from_family",
    "to_family",
    "crossing_kind",
    "quality_flag",
)

CARRIER_HISTORY_COLUMNS = (
    "run_id",
    "source",
    "snapshot",
    "z",
    "a",
    "carrier",
    "presence",
    "source_total",
    "sink_total",
    "net_total",
    "dominant_channel",
    "quality_flag",
)

EMERGENCE_EVENT_COLUMNS = (
    "run_id",
    "source",
    "snapshot",
    "z",
    "a",
    "event_name",
    "band_id",
    "criterion",
    "value",
    "dominant_channel",
    "physical_reading",
    "quality_flag",
)


@dataclass(frozen=True)
class InteractionChannel:
    channel_id: str
    family: str
    carrier: str
    target: str
    mediator: str
    source_weight: float
    sink_weight: float
    state_exponents: Mapping[str, float]
    scale_exponents: Mapping[str, float]
    carrier_exponents: Mapping[str, float]
    output_table: str
    physical_reading: str


INTERACTION_CHANNELS = (
    InteractionChannel(
        "metric_seed",
        "primordial",
        "metric_spectrum",
        "density_field",
        "curvature",
        1.0,
        0.02,
        {"early_weight": 1.2, "radiation_weight": 0.4},
        {"horizon_weight": 1.0},
        {"field_presence": 1.0},
        "interaction_strength_table.csv",
        "coarse spectral seed carried by the density field",
    ),
    InteractionChannel(
        "weak_rate",
        "weak",
        "thermal_species",
        "neutrino_component",
        "weak_current",
        0.8,
        0.3,
        {"temperature_weight": 1.1, "expansion_rate": -0.5},
        {"horizon_weight": 0.6},
        {"radiation_presence": 0.7},
        "interaction_strength_table.csv",
        "weak exchange strength relative to the generated state",
    ),
    InteractionChannel(
        "nuclear_binding",
        "nuclear",
        "baryon_component",
        "light_element_record",
        "nuclear_reaction",
        0.7,
        0.35,
        {"temperature_weight": 0.8, "baryon_weight": 1.0, "expansion_rate": -0.4},
        {"horizon_weight": 0.5},
        {"matter_presence": 0.5},
        "interaction_strength_table.csv",
        "binding source carried by the baryonic component",
    ),
    InteractionChannel(
        "radiation_drag",
        "radiation",
        "photon_baryon_field",
        "baryon_velocity",
        "thomson_drag",
        0.9,
        0.6,
        {"radiation_weight": 1.1, "ionization_weight": 0.8},
        {"sound_weight": 0.9},
        {"radiation_presence": 1.0, "matter_presence": 0.4},
        "interaction_strength_table.csv",
        "radiation coupling in the field carrier",
    ),
    InteractionChannel(
        "atomic_binding",
        "electromagnetic",
        "baryon_component",
        "neutral_component",
        "atomic_transition",
        0.75,
        0.4,
        {"cooling_weight": 1.0, "ionization_inverse": 0.8},
        {"sound_weight": 0.4},
        {"matter_presence": 0.7},
        "interaction_strength_table.csv",
        "atomic binding contribution to carrier conversion",
    ),
    InteractionChannel(
        "free_streaming",
        "free_streaming",
        "neutrino_component",
        "spectral_tail",
        "phase_mixing",
        0.55,
        0.15,
        {"early_weight": 0.4, "expansion_rate": 0.3},
        {"small_scale_weight": 1.0},
        {"tail_presence": 0.8},
        "tail_table.csv",
        "phase mixing source for the residual tail",
    ),
    InteractionChannel(
        "gravitational_growth",
        "gravity",
        "matter_field",
        "active_window",
        "Newtonian_potential",
        1.0,
        0.1,
        {"growth_weight": 1.3, "matter_weight": 0.8},
        {"frontier_weight": 1.2},
        {"field_presence": 0.7, "tail_presence": 0.2},
        "global_window.csv",
        "gravitational growth of the nonlinear endpoint",
    ),
    InteractionChannel(
        "gas_cooling",
        "gas",
        "baryon_gas",
        "compact_baryon_component",
        "radiative_cooling",
        0.65,
        0.25,
        {"cooling_weight": 1.1, "matter_weight": 0.6},
        {"halo_scale_weight": 0.8},
        {"object_presence": 0.6, "matter_presence": 0.4},
        "object_table.csv",
        "cooling source for compact baryonic carriers",
    ),
    InteractionChannel(
        "object_virialization",
        "object",
        "active_window",
        "object_layer",
        "self_gravity",
        0.85,
        0.2,
        {"growth_weight": 0.8, "object_weight": 1.0},
        {"halo_scale_weight": 1.0},
        {"object_presence": 1.0},
        "object_table.csv",
        "conversion of nonlinear field content into object records",
    ),
    InteractionChannel(
        "tidal_transfer",
        "object",
        "object_pair",
        "spin_record",
        "tidal_tensor",
        0.6,
        0.15,
        {"object_weight": 1.0, "tail_weight": 0.3},
        {"halo_scale_weight": 0.9},
        {"object_presence": 1.0},
        "transfer_table.csv",
        "directional transfer between object carriers",
    ),
    InteractionChannel(
        "feedback_return",
        "feedback",
        "compact_object",
        "gas_and_tail",
        "coarse_feedback",
        0.45,
        0.5,
        {"object_weight": 1.0, "tail_weight": 0.6},
        {"halo_scale_weight": 0.7},
        {"object_presence": 0.8, "tail_presence": 0.5},
        "tail_table.csv",
        "return of compact information to gas and tail carriers",
    ),
    InteractionChannel(
        "residual_feedback",
        "residual",
        "spectral_tail",
        "active_window",
        "coarse_source",
        0.5,
        0.35,
        {"tail_weight": 1.1, "growth_weight": 0.3},
        {"frontier_weight": 1.0},
        {"tail_presence": 1.0},
        "tail_table.csv",
        "feedback of residual tail activity into the active window",
    ),
    InteractionChannel(
        "dark_energy_suppression",
        "dark_energy",
        "background",
        "endpoint_motion",
        "expansion",
        0.3,
        0.85,
        {"future_weight": 1.2, "expansion_rate": 0.8},
        {"horizon_weight": 0.8},
        {"field_presence": 0.5},
        "global_window.csv",
        "background suppression of future endpoint motion",
    ),
)


def generate_interaction_outputs(
    *,
    run_id: str,
    global_rows: Sequence[Mapping[str, object]],
    band_rows: Sequence[Mapping[str, object]],
    object_rows: Sequence[Mapping[str, object]],
    transfer_rows: Sequence[Mapping[str, object]],
    channels: Sequence[InteractionChannel] = INTERACTION_CHANNELS,
) -> dict[str, list[dict[str, object]]]:
    object_stats = _object_stats(object_rows)
    transfer_stats = _transfer_stats(transfer_rows)
    strength_rows = interaction_strength_rows(
        run_id=run_id,
        global_rows=global_rows,
        band_rows=band_rows,
        object_stats=object_stats,
        transfer_stats=transfer_stats,
        channels=channels,
    )
    budget_rows = interaction_budget_rows(strength_rows)
    dominance_rows = interaction_dominance_rows(strength_rows)
    crossing_rows = interaction_crossing_rows(dominance_rows)
    carrier_rows = carrier_history_rows(strength_rows)
    event_rows = emergence_event_rows(dominance_rows, carrier_rows)
    return {
        "strength": strength_rows,
        "budget": budget_rows,
        "dominance": dominance_rows,
        "crossing": crossing_rows,
        "carrier": carrier_rows,
        "emergence": event_rows,
    }


def interaction_strength_rows(
    *,
    run_id: str,
    global_rows: Sequence[Mapping[str, object]],
    band_rows: Sequence[Mapping[str, object]],
    object_stats: Mapping[tuple[str, str], Mapping[str, float]],
    transfer_stats: Mapping[tuple[str, str], Mapping[str, float]],
    channels: Sequence[InteractionChannel] = INTERACTION_CHANNELS,
) -> list[dict[str, object]]:
    state_rows = _interaction_state_rows(global_rows)
    global_by_key = {_epoch_key(row): row for row in state_rows}
    grouped_bands: dict[str, list[Mapping[str, object]]] = {}
    for row in band_rows:
        grouped_bands.setdefault(_epoch_key(row), []).append(row)
    output: list[dict[str, object]] = []
    for epoch_key in sorted(global_by_key, key=lambda value: -float(value.split("|")[0])):
        global_row = global_by_key[epoch_key]
        bands = grouped_bands.get(epoch_key, _default_interaction_bands(global_row))
        state = state_coordinates(global_row)
        raw_rows = []
        for band in sorted(bands, key=lambda row: str(row["band_id"])):
            band_state = band_coordinates(band)
            z_key = _z_key(global_row["z"])
            obj_state = object_stats.get((z_key, str(band["band_id"])), {})
            transfer_state = transfer_stats.get((z_key, str(band["band_id"])), {})
            for channel in channels:
                source, sink, raw, factors = evaluate_channel(channel, state, band_state, obj_state, transfer_state)
                raw_rows.append((band, channel, source, sink, raw, factors))
        totals_by_band: dict[str, float] = {}
        for band, _channel, _source, _sink, raw, _factors in raw_rows:
            band_id = str(band["band_id"])
            totals_by_band[band_id] = totals_by_band.get(band_id, 0.0) + raw
        for band, channel, source, sink, raw, factors in raw_rows:
            band_id = str(band["band_id"])
            total = totals_by_band[band_id]
            chi = raw / total if total > 0 else 0.0
            output.append(
                {
                    "run_id": run_id,
                    "source": "internal_generation",
                    "snapshot": global_row["snapshot"],
                    "z": global_row["z"],
                    "a": global_row["a"],
                    "band_id": band_id,
                    "channel_id": channel.channel_id,
                    "family": channel.family,
                    "carrier": channel.carrier,
                    "target": channel.target,
                    "mediator": channel.mediator,
                    "source_term": source,
                    "sink_term": sink,
                    "net_term": source - sink,
                    "raw_strength": raw,
                    "chi": chi,
                    "state_factor": factors["state"],
                    "scale_factor": factors["scale"],
                    "carrier_factor": factors["carrier"],
                    "output_table": channel.output_table,
                    "physical_reading": channel.physical_reading,
                }
            )
    return output


def evaluate_channel(
    channel: InteractionChannel,
    state: Mapping[str, float],
    band: Mapping[str, float],
    objects: Mapping[str, float],
    transfers: Mapping[str, float],
) -> tuple[float, float, float, dict[str, float]]:
    variables = dict(state)
    variables.update(band)
    variables.update(objects)
    variables.update(transfers)
    state_factor = _product(channel.state_exponents, variables)
    scale_factor = _product(channel.scale_exponents, variables)
    carrier_factor = _product(channel.carrier_exponents, variables)
    source = channel.source_weight * state_factor * scale_factor * carrier_factor
    sink = channel.sink_weight * (1.0 + 0.5 * state_factor) * (0.5 + carrier_factor)
    raw = max(0.0, source + 0.25 * abs(sink))
    return source, sink, raw, {"state": state_factor, "scale": scale_factor, "carrier": carrier_factor}


def state_coordinates(row: Mapping[str, object]) -> dict[str, float]:
    z = float(row["z"])
    a = float(row["a"])
    k_nl = _float(row, "k_NL_hmpc", 0.0)
    m_nl = _float(row, "M_NL_msun_h", 0.0)
    omega_tail = _float(row, "Omega_tail", 0.0)
    omega_excess = _float(row, "Omega_excess", 0.0)
    r_h = _float(row, "R_h_mean", 0.0)
    v_lambda = abs(_float(row, "v_lambda", 0.0))
    return {
        "early_weight": 1.0 / (1.0 + a),
        "future_weight": max(0.0, a - 1.0) / a,
        "temperature_weight": max(1.0, 1.0 + z),
        "expansion_rate": math.sqrt(max(1.0e-30, (1.0 + z) ** 3 + 1.0)),
        "growth_weight": 1.0 / (1.0 + v_lambda),
        "matter_weight": m_nl / (m_nl + 1.0e12),
        "radiation_weight": (1.0 + z) / (1.0 + z + 3400.0) if z > -1.0 else 0.0,
        "baryon_weight": 0.16,
        "cooling_weight": 1.0 / (1.0 + max(0.0, z) / 20.0),
        "ionization_weight": max(0.0, min(1.0, (1.0 + z) / 1200.0)),
        "ionization_inverse": 1.0 - max(0.0, min(1.0, (1.0 + z) / 1200.0)),
        "tail_weight": omega_excess,
        "tail_presence": omega_tail,
        "object_weight": r_h,
        "field_presence": max(1.0e-12, min(1.0, k_nl / 0.2)),
    }


def band_coordinates(row: Mapping[str, object]) -> dict[str, float]:
    band_id = str(row["band_id"])
    omega = float(row["omega_band"])
    k_min = float(row["k_min_hmpc"])
    k_max = float(row["k_max_hmpc"])
    k_mid = math.sqrt(k_min * k_max)
    return {
        "omega_band": max(1.0e-12, omega),
        "frontier_weight": 1.0 if row.get("activity_status") == "frontier" else 0.35,
        "small_scale_weight": max(1.0e-12, min(20.0, k_mid / 0.2)),
        "halo_scale_weight": _halo_scale_weight(band_id),
        "sound_weight": 1.0 / (1.0 + k_mid),
        "horizon_weight": 1.0 / (1.0 + k_mid / 1.0e-3),
    }


def interaction_budget_rows(rows: Iterable[Mapping[str, object]]) -> list[dict[str, object]]:
    rows = list(rows)
    groups = sorted(
        {(_epoch_key(row), str(row["band_id"]), str(row["carrier"])) for row in rows},
        key=lambda item: (-float(item[0].split("|")[0]), item[1], item[2]),
    )
    output = []
    for epoch_key, band_id, carrier in groups:
        selected = [row for row in rows if _epoch_key(row) == epoch_key and str(row["band_id"]) == band_id and row["carrier"] == carrier]
        first = selected[0]
        output.append(
            {
                "run_id": first["run_id"],
                "source": first["source"],
                "snapshot": first["snapshot"],
                "z": first["z"],
                "a": first["a"],
                "band_id": band_id,
                "carrier": carrier,
                "source_total": _sum(selected, "source_term"),
                "sink_total": _sum(selected, "sink_term"),
                "net_total": _sum(selected, "net_term"),
                "chi_total": _sum(selected, "chi"),
                "n_channels": len(selected),
                "quality_flag": "internal_generation",
            }
        )
    return output


def interaction_dominance_rows(rows: Iterable[Mapping[str, object]]) -> list[dict[str, object]]:
    rows = list(rows)
    groups = sorted(
        {(_epoch_key(row), str(row["band_id"])) for row in rows},
        key=lambda item: (-float(item[0].split("|")[0]), item[1]),
    )
    output = []
    for epoch_key, band_id in groups:
        selected = [row for row in rows if _epoch_key(row) == epoch_key and str(row["band_id"]) == band_id]
        ordered = sorted(selected, key=lambda row: float(row["chi"]), reverse=True)
        first = ordered[0]
        second = ordered[1] if len(ordered) > 1 else ordered[0]
        chi_sum = _sum(selected, "chi")
        output.append(
            {
                "run_id": first["run_id"],
                "source": first["source"],
                "snapshot": first["snapshot"],
                "z": first["z"],
                "a": first["a"],
                "band_id": band_id,
                "dominant_channel": first["channel_id"],
                "dominant_family": first["family"],
                "dominant_chi": first["chi"],
                "second_channel": second["channel_id"],
                "second_family": second["family"],
                "second_chi": second["chi"],
                "dominance_gap": float(first["chi"]) - float(second["chi"]),
                "chi_sum": chi_sum,
                "quality_flag": "internal_generation",
            }
        )
    return output


def interaction_crossing_rows(rows: Iterable[Mapping[str, object]]) -> list[dict[str, object]]:
    rows_by_band: dict[str, list[Mapping[str, object]]] = {}
    for row in rows:
        rows_by_band.setdefault(str(row["band_id"]), []).append(row)
    output = []
    for band_id, selected in sorted(rows_by_band.items()):
        ordered = sorted(selected, key=lambda row: float(row["a"]))
        for left, right in zip(ordered, ordered[1:]):
            if left["dominant_channel"] == right["dominant_channel"]:
                continue
            output.append(
                {
                    "run_id": right["run_id"],
                    "source": right["source"],
                    "band_id": band_id,
                    "left_snapshot": left["snapshot"],
                    "right_snapshot": right["snapshot"],
                    "z_left": left["z"],
                    "z_right": right["z"],
                    "a_left": left["a"],
                    "a_right": right["a"],
                    "from_channel": left["dominant_channel"],
                    "to_channel": right["dominant_channel"],
                    "from_family": left["dominant_family"],
                    "to_family": right["dominant_family"],
                    "crossing_kind": "dominance_change",
                    "quality_flag": "internal_generation",
                }
            )
    return output


def carrier_history_rows(rows: Iterable[Mapping[str, object]]) -> list[dict[str, object]]:
    budget = interaction_budget_rows(rows)
    output = []
    for row in budget:
        presence = min(1.0, max(0.0, float(row["chi_total"])))
        dominant = _dominant_for_carrier(rows, row)
        output.append(
            {
                "run_id": row["run_id"],
                "source": row["source"],
                "snapshot": row["snapshot"],
                "z": row["z"],
                "a": row["a"],
                "carrier": row["carrier"],
                "presence": presence,
                "source_total": row["source_total"],
                "sink_total": row["sink_total"],
                "net_total": row["net_total"],
                "dominant_channel": dominant,
                "quality_flag": row["quality_flag"],
            }
        )
    return output


def emergence_event_rows(
    dominance_rows: Sequence[Mapping[str, object]],
    carrier_rows: Sequence[Mapping[str, object]],
) -> list[dict[str, object]]:
    _ = carrier_rows
    return _dominance_change_events(dominance_rows)


def _dominance_change_events(rows: Sequence[Mapping[str, object]]) -> list[dict[str, object]]:
    events = []
    grouped: dict[str, list[Mapping[str, object]]] = {}
    for row in rows:
        grouped.setdefault(str(row["band_id"]), []).append(row)
    for band_id, selected in sorted(grouped.items()):
        ordered = sorted(selected, key=lambda row: float(row["a"]))
        for left, right in zip(ordered, ordered[1:]):
            if left["dominant_channel"] == right["dominant_channel"]:
                continue
            events.append(
                {
                    "run_id": right["run_id"],
                    "source": right["source"],
                    "snapshot": right["snapshot"],
                    "z": right["z"],
                    "a": right["a"],
                    "event_name": "dominant_interaction_change",
                    "band_id": band_id,
                    "criterion": "dominant channel changes between adjacent stored states",
                    "value": right["dominant_chi"],
                    "dominant_channel": right["dominant_channel"],
                    "physical_reading": (
                        f"{left['dominant_channel']} gives way to {right['dominant_channel']} "
                        f"in the stored {band_id} trajectory"
                    ),
                    "quality_flag": right["quality_flag"],
                }
            )
    return events


def _object_stats(rows: Iterable[Mapping[str, object]]) -> dict[tuple[str, str], dict[str, float]]:
    grouped: dict[tuple[str, str], list[Mapping[str, object]]] = {}
    for row in rows:
        grouped.setdefault((_z_key(row["z"]), _band_from_mass(float(row["mass_msunh"]))), []).append(row)
    output = {}
    for key, selected in grouped.items():
        b_values = [abs(float(row["b_h"])) for row in selected]
        output[key] = {
            "object_presence": min(1.0, len(selected) / 24.0),
            "object_weight": _mean(b_values) if b_values else 0.0,
        }
    return output


def _transfer_stats(rows: Iterable[Mapping[str, object]]) -> dict[tuple[str, str], dict[str, float]]:
    grouped: dict[tuple[str, str], list[Mapping[str, object]]] = {}
    for row in rows:
        grouped.setdefault((_z_key(row["z"]), "all"), []).append(row)
    output: dict[tuple[str, str], dict[str, float]] = {}
    for (epoch_key, _band), selected in grouped.items():
        strength = _mean([float(row["a_pq"]) for row in selected])
        output[(epoch_key, "all")] = {"transfer_weight": strength}
        for band_id in ("dwarf_halo", "galaxy_halo", "group_halo", "cluster_halo", "large_environment"):
            output[(epoch_key, band_id)] = {"transfer_weight": strength}
    return output


def _product(exponents: Mapping[str, float], values: Mapping[str, float]) -> float:
    total = 1.0
    for key, exponent in exponents.items():
        value = max(1.0e-12, float(values.get(key, 0.0)))
        total *= value**exponent
    return total


def _sum(rows: Iterable[Mapping[str, object]], key: str) -> float:
    return sum(float(row[key]) for row in rows)


def _mean(values: Sequence[float]) -> float:
    return sum(values) / len(values) if values else 0.0


def _dominant_for_carrier(rows: Iterable[Mapping[str, object]], budget_row: Mapping[str, object]) -> str:
    selected = [
        row
        for row in rows
        if _epoch_key(row) == _epoch_key(budget_row)
        and row["carrier"] == budget_row["carrier"]
        and row["band_id"] == budget_row["band_id"]
    ]
    if not selected:
        return ""
    return str(max(selected, key=lambda row: float(row["chi"]))["channel_id"])


def _epoch_key(row: Mapping[str, object]) -> str:
    return f"{float(row['z']):.12g}|{float(row['a']):.12g}"


def _z_key(value: object) -> str:
    return f"{float(value):.12g}"


def _interaction_state_rows(global_rows: Sequence[Mapping[str, object]]) -> list[dict[str, object]]:
    rows = {(_z_key(row["z"])): dict(row) for row in global_rows}
    for z in _log_redshift_grid(max_z=1.0e12, min_z=100.0, count=16):
        key = _z_key(z)
        if key in rows:
            continue
        a = 1.0 / (1.0 + z)
        rows[key] = {
            "run_id": global_rows[0]["run_id"] if global_rows else "",
            "source": "internal_generation",
            "snapshot": f"a{_tag(a)}",
            "z": z,
            "a": a,
            "k_NL_hmpc": 0.0,
            "M_NL_msun_h": 0.0,
            "Omega_tail": 0.0,
            "Omega_excess": 0.0,
            "R_h_mean": 0.0,
            "v_lambda": 0.0,
        }
    return [rows[key] for key in sorted(rows, key=lambda value: -float(value))]


def _default_interaction_bands(row: Mapping[str, object]) -> list[dict[str, object]]:
    fixed_ranges = {
        "dwarf_halo": (1.0, 20.0),
        "galaxy_halo": (0.1, 1.0),
        "group_halo": (0.03, 0.1),
        "cluster_halo": (0.01, 0.03),
        "large_environment": (0.001, 0.01),
    }
    rows = []
    for band_id, mass_min, mass_max, label in DEFAULT_MASS_BANDS:
        k_min, k_max = fixed_ranges[band_id]
        rows.append(
            {
                "z": row["z"],
                "a": row["a"],
                "band_id": band_id,
                "k_min_hmpc": k_min,
                "k_max_hmpc": k_max,
                "m_min_msunh": mass_min,
                "m_max_msunh": mass_max,
                "omega_band": 0.0,
                "activity_status": "not_resolved",
                "representative_structure": label,
            }
        )
    return rows


def _log_redshift_grid(*, max_z: float, min_z: float, count: int) -> list[float]:
    hi = math.log1p(max_z)
    lo = math.log1p(min_z)
    step = (hi - lo) / (count - 1)
    return [math.expm1(hi - index * step) for index in range(count)]


def _float(row: Mapping[str, object], key: str, default: float) -> float:
    value = row.get(key, default)
    if value in ("", None):
        return default
    return float(value)


def _tag(value: float) -> str:
    return f"{value:g}".replace("-", "m").replace(".", "p")


def _band_from_mass(mass_msunh: float) -> str:
    if mass_msunh < 1.0e11:
        return "dwarf_halo"
    if mass_msunh < 1.0e13:
        return "galaxy_halo"
    if mass_msunh < 1.0e14:
        return "group_halo"
    if mass_msunh < 1.0e15:
        return "cluster_halo"
    return "large_environment"


def _halo_scale_weight(band_id: str) -> float:
    weights = {
        "dwarf_halo": 0.8,
        "galaxy_halo": 1.0,
        "group_halo": 0.9,
        "cluster_halo": 0.75,
        "large_environment": 0.45,
    }
    return weights.get(band_id, 0.5)
