"""Publish the current finite cosmological-clock and active-evolution checks.

Use --rebind-only with --params-file to bind identical numerical inputs
without repeating the calculation. --publish-input also publishes the
current parameter snapshot. A run without --rebind-only evaluates the full
numerical verification.
"""
import argparse
from fractions import Fraction
from dataclasses import replace
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.integrate import solve_ivp

from cosmic_cg.active_spectral import verify_spectral_pair
from cosmic_cg.cosmological_dynamics import coefficients, expansion, growing_rate, linear_growth_advance
from cosmic_cg.params import load_model_params, write_param_snapshot
from cosmic_cg.pm_nbody import PMNBodyConfig, _spectral_setup, _initial_state, _advance_to, _cic_delta
from cosmic_cg.pm_certificate import certify_pm_trajectory


ROOT = Path(__file__).resolve().parents[1]


def evaluate(params):
    eds = replace(params, omega_b=0.0, omega_sigma=1.0, omega_lambda=0.0)
    linear = []
    for steps in (32, 64, 128):
        displacement, velocity = linear_growth_advance(0.1, 1.0, steps, eds, displacement=0.1, velocity=1.0)
        linear.append({"steps": steps, "displacement_error": abs(displacement-1.0), "velocity_error": abs(velocity-1.0)})
    a0, a1 = 0.2, 0.4
    f0 = growing_rate(a0, params)
    def rhs(a, state):
        drag, force = coefficients(a, params)
        return [state[1], -drag*state[1]+force*state[0]]
    reference = solve_ivp(rhs, (a0, a1), [1.0, f0/a0], method="DOP853", rtol=1e-12, atol=1e-13)
    if not reference.success:
        raise RuntimeError(reference.message)
    displacement, velocity = linear_growth_advance(a0, a1, 128, params)
    clock_error = max(abs(displacement-reference.y[0, -1]), abs(velocity-reference.y[1, -1]))
    positions = []
    pm_rows = []
    for steps in (16, 32, 64):
        config = PMNBodyConfig(dynamics="cosmological", particles_per_dim=16, mesh_size=16,
                               box_size_mpch=32.0, z_initial=1/a0-1,
                               steps_per_interval=steps, initial_displacement_fraction=0.002,
                               mode_bins=8, max_halos_per_snapshot=16)
        if steps == 64:
            pm_certificate=certify_pm_trajectory(config,params,a0=a0,a1=a1,return_state=True)
            state=pm_certificate.pop("_final_state")
        else:
            setup = _spectral_setup(config)
            setup["params"] = params
            state = _initial_state(config, setup, np)
            _advance_to(state, a0, a1, config, setup)
        density = _cic_delta(state["positions"], config)
        positions.append(state["positions"].copy())
        pm_rows.append({"steps": steps,
                        "density_mean_error": abs(float(np.mean(density, dtype=np.float64))),
                        "minimum_density": float(np.min(1+density)),
                        "rms_density": float(np.std(density)),
                        "particle_count": int(len(state["ids"])),
                        "force_evaluations": steps+1})
    diffs = []
    for first, second in zip(positions[:-1], positions[1:]):
        delta = (first-second+16.0) % 32.0-16.0
        diffs.append(float(np.sqrt(np.mean(delta**2))))
    active = verify_spectral_pair()
    active_scaling = verify_spectral_pair(full_modes=32, duration=Fraction(1,20), step=Fraction(1,1000))
    checks = {"EdS_second_order": all(linear[i]["displacement_error"] > 3.0*linear[i+1]["displacement_error"] for i in range(2)),
              "current_background_clock": clock_error < 2e-5,
              "PM_step_refinement": diffs[1] < 0.6*diffs[0],
              "PM_mass_conservation": max(r["density_mean_error"] for r in pm_rows) < 2e-7,
              "PM_density_positive": min(r["minimum_density"] for r in pm_rows) >= -2e-7,
              "PM_trajectory_enclosure":pm_certificate["all_checks_passed"],
              "PM_initial_Newtonian_force":pm_certificate["initial_Newtonian_force_observables"]["all_checks_passed"],
              "active_evolution": active["all_checks_passed"],
              "active_scaling": active_scaling["all_checks_passed"]}
    checks = {key: bool(value) for key, value in checks.items()}
    if not all(checks.values()):
        raise RuntimeError(f"Cosmological verification failed: {checks}; clock={clock_error}; PM={diffs}")
    return {"scope": "specified cold-source Vlasov--Poisson realization; finite numerical validation",
            "upstream_run_id": params.source_run_id, "upstream_parameter_sha256": params.source_hash,
            "upstream_reference": params.source_reference,
            "clock": {"time_variable": "a", "momentum": "p=a^3 E(a) dx/da", "EdS": linear,
                      "current_background_interval": [a0, a1], "current_background_max_error": clock_error},
            "PM": {"particles": 4096, "mesh": [16, 16, 16], "box_Mpc_over_h": 32.0,
                   "seed": 20260904, "initial_displacement_fraction": 0.002,
                   "initial_velocity": "growing mode f(a0) displacement/a0",
                   "scope": "finite CIC particle trajectory enclosure, mass and time-refinement checks",
                   "rows": pm_rows, "successive_position_RMS_Mpc_over_h": diffs,
                   "trajectory_certificate":pm_certificate},
            "active_spectral": active, "active_scaling": active_scaling,
            "checks": checks, "all_checks_passed": all(checks.values())}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--params-file")
    parser.add_argument("--output", default=str(ROOT / "results" / "cosmological_active_current" / "results.json"))
    parser.add_argument("--publish-input", action="store_true")
    parser.add_argument("--rebind-only", action="store_true",
                        help="Verify identical numerical inputs and bind the saved evidence; no evolution or certificate evaluation")
    args = parser.parse_args()
    params = load_model_params(args.params_file)
    params = replace(params, source_reference="Coarse-Graining Genesis Framework 4.0.0.dev20261005 (unpublished revision)")
    if args.publish_input and (not args.params_file or not params.source_run_id or not params.source_hash):
        raise ValueError("Publication requires a fingerprinted V4 release input")
    if args.rebind_only:
        from verify_pm_certificate_identity import verify, parameter_values, MANIFEST
        previous = load_model_params(ROOT / "tables" / "v4_parameter_snapshot.json")
        if parameter_values(previous) != parameter_values(params):
            raise ValueError("Identity-only binding requires bitwise-equal consumed parameters")
        checks = verify(params)
        result = json.loads((ROOT / "results/cosmological_active_current/results.json").read_text(encoding="utf-8"))
        result.update(upstream_run_id=params.source_run_id,
                      upstream_parameter_sha256=params.source_hash,
                      upstream_reference=params.source_reference)
        result["identity_binding"] = checks
    else:
        result = evaluate(params)
        files = sorted((ROOT / "code").rglob("*.py"))
        result["source_sha256"] = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest() for p in files}
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    write_param_snapshot(params, output.parent / "params_snapshot.json")
    if args.publish_input:
        write_param_snapshot(params, ROOT / "tables" / "v4_parameter_snapshot.json")
        if args.rebind_only:
            deployment = json.loads(MANIFEST.read_text(encoding="utf-8"))
            deployment["current_input_binding"] = {
                "source_run_id": params.source_run_id,
                "source_hash": params.source_hash,
                "source_reference": params.source_reference,
                "consumed_values_identical": True,
                "new_force_evaluations": 0, "certificate_recomputed": False}
            MANIFEST.write_text(json.dumps(deployment, indent=2, sort_keys=True)+"\n", encoding="utf-8")
    print(json.dumps({"all_checks_passed": True, "upstream_run_id": params.source_run_id,
                      "active_error_upper": result["active_spectral"]["uniform_error_upper"],
                      "active_advanced_coordinates": result["active_spectral"]["work"]["advancement_float"]["rhs_coordinates"],
                      "output": str(output)}, indent=2))


if __name__ == "__main__":
    main()
