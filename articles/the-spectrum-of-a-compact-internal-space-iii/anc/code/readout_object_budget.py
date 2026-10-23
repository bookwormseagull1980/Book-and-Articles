# -*- coding: utf-8 -*-
"""Recompute particle-count conservation and growth from existing PM tables.

Disjointness follows from the generator's component assignment.  This
readout checks its persisted count consequences, not unavailable raw IDs.
Missing flux or event records never become zero-valued certificates.
"""
import csv
import math
from collections import defaultdict
from pathlib import Path

RUN = Path(__file__).resolve().parents[1] / "results" / "pm_dynamics_64_protocol_001"


def rows(path, required):
    with Path(path).open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        missing = set(required) - set(reader.fieldnames or ())
        if missing:
            raise ValueError(f"{path}: missing columns {sorted(missing)}")
        return list(reader)


def assemble(snap, halo, tree):
    if not snap:
        raise ValueError("Empty snapshot table")
    if len({row["snapshot"] for row in snap}) != len(snap):
        raise ValueError("Duplicate snapshot ID")
    if not all(math.isfinite(float(row["a"])) and float(row["a"]) > 0 for row in snap):
        raise ValueError("Invalid scale factor")
    snap = sorted(snap, key=lambda row: float(row["a"]))
    if any(float(b["a"]) <= float(a["a"]) for a, b in zip(snap, snap[1:])):
        raise ValueError("Scale factors are not strictly increasing")
    n = int(snap[0]["particle_count"])
    mass = float(snap[0]["particle_mass_msunh"])
    if n <= 0 or not math.isfinite(mass) or mass <= 0:
        raise ValueError("Invalid particle mass or count")
    by_snap = defaultdict(list)
    by_id = {}
    for obj in halo:
        key = (obj["snapshot"], obj["halo_id"])
        if key in by_id:
            raise ValueError("Duplicate object record")
        count = int(obj["particle_count"])
        measured = float(obj["mass_msunh"])
        if count <= 0 or not math.isclose(measured, count * mass, rel_tol=1e-12):
            raise ValueError("Object mass does not match particle count")
        by_id[key] = obj
        by_snap[obj["snapshot"]].append(obj)
    known = defaultdict(int)
    parent_shared = defaultdict(int)
    seen_links = set()
    growth = []
    for link in tree:
        key = (link["snapshot_from"], link["snapshot_to"], link["descendant_id"])
        if key in seen_links:
            raise ValueError("More than one best-parent link for a descendant")
        seen_links.add(key)
        parent = by_id[(link["snapshot_from"], link["parent_id"])]
        child = by_id[(link["snapshot_to"], link["descendant_id"])]
        shared = int(link["shared_particle_count"])
        if not 0 <= shared <= min(int(parent["particle_count"]), int(child["particle_count"])):
            raise ValueError("Impossible shared-particle count")
        parent_key = (*key[:2], link["parent_id"])
        parent_shared[parent_key] += shared
        if parent_shared[parent_key] > int(parent["particle_count"]):
            raise ValueError("Parent particles counted more than once")
        known[key[:2]] += shared
        mp, md = float(link["mass_from_msunh"]), float(link["mass_to_msunh"])
        if not math.isclose(mp, float(parent["mass_msunh"]), rel_tol=1e-12):
            raise ValueError("Parent link mass mismatch")
        if not math.isclose(md, float(child["mass_msunh"]), rel_tol=1e-12):
            raise ValueError("Descendant link mass mismatch")
        dt, td = float(link["delta_time_gyr"]), float(child["t_dyn_gyr"])
        if min(mp, md, dt, td) <= 0 or not all(map(math.isfinite, (mp, md, dt, td))):
            raise ValueError("Invalid logarithmic-growth domain")
        growth.append(dict(snapshot=child["snapshot"], halo_id=child["halo_id"],
                           b_log=td * math.log(md / mp) / dt,
                           b_saved=float(child["b_h"]), status="linked_positive_masses"))
    budget = []
    previous = None
    for row in snap:
        key = row["snapshot"]
        if int(row["particle_count"]) != n or float(row["particle_mass_msunh"]) != mass:
            raise ValueError("Particle population changed")
        objects = by_snap.pop(key, [])
        if len(objects) != int(row["halo_count"]):
            raise ValueError("Snapshot/catalog count mismatch")
        no = sum(int(obj["particle_count"]) for obj in objects)
        if not 0 <= no <= n:
            raise ValueError("Object allocation exceeds total population")
        momentum = float(row["momentum_norm_proxy"])
        if not math.isfinite(momentum) or momentum < 0:
            raise ValueError("Invalid momentum norm")
        record = dict(snapshot=key, a=float(row["a"]), object_particles=no,
                      diffuse_particles=n-no, total_particles=n, net_particle_error=0,
                      shared_lower="", shared_upper="", diffuse_in_lower="",
                      diffuse_in_upper="", object_out_lower="", object_out_upper="",
                      momentum_drift_per_mass_upper="", reentry_status="not_recorded",
                      angular_flux_status="not_recorded")
        if previous is not None:
            pair = (previous["snapshot"], key)
            np = previous["object_particles"]
            lo = max(known.pop(pair, 0), np + no - n)
            hi = min(np, no)
            if lo > hi:
                raise ValueError("Inconsistent parent intersections")
            record.update(shared_lower=lo, shared_upper=hi,
                          diffuse_in_lower=no-hi, diffuse_in_upper=no-lo,
                          object_out_lower=np-hi, object_out_upper=np-lo,
                          net_particle_error=(no-np)+((n-no)-(n-np)),
                          momentum_drift_per_mass_upper=momentum+previous["momentum"])
        budget.append(record)
        previous = dict(snapshot=key, object_particles=no, momentum=momentum)
    if by_snap or known:
        raise ValueError("Object/tree epoch not on the snapshot chain")
    return budget, growth


def write_table(path, records):
    if not records:
        raise ValueError("Cannot certify an empty result")
    with Path(path).open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def main():
    tables = RUN / "tables"
    snap = rows(tables / "pm_snapshot_table.csv",
                ["snapshot", "a", "particle_count", "particle_mass_msunh", "halo_count",
                 "momentum_norm_proxy", "mass_conservation_error"])
    halo = rows(tables / "pm_halo_table.csv",
                ["snapshot", "halo_id", "particle_count", "mass_msunh", "t_dyn_gyr", "b_h"])
    tree = rows(tables / "pm_tree_table.csv",
                ["snapshot_from", "snapshot_to", "parent_id", "descendant_id",
                 "shared_particle_count", "mass_from_msunh", "mass_to_msunh", "delta_time_gyr"])
    budget, growth = assemble(snap, halo, tree)
    write_table(RUN / "pm_object_budget.csv", budget)
    write_table(RUN / "pm_object_growth_log.csv", growth)
    print(f"Certified integer mass partitions: {len(budget)} epochs")
    print(f"Linked logarithmic growth readings: {len(growth)}")
    print(f"Final object particles: {budget[-1]['object_particles']}")
    print("Max momentum drift per mass bound:",
          max(float(r["momentum_drift_per_mass_upper"] or 0) for r in budget))
    print("Re-entry events and angular transfer: not certified by the saved tables")


if __name__ == "__main__":
    main()
