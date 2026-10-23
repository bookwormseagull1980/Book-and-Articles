import json
import math
from dataclasses import replace
from fractions import Fraction as Q
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import numpy as np

from cosmic_cg.params import load_model_params, write_param_snapshot
from cosmic_cg.cosmological_dynamics import coefficients, expansion, growing_rate, linear_growth_advance
from cosmic_cg.active_spectral import matter_rhs, sustained_active_run, verify_spectral_pair, ReserveBudget
from cosmic_cg.pm_nbody import PMNBodyConfig, _advance_to


class ParameterBindingTests(unittest.TestCase):
    def test_flat_snapshot_changes_are_consumed_and_round_trip(self):
        current = load_model_params()
        changed = replace(current, omega_b=current.omega_b+0.001, omega_sigma=current.omega_sigma-0.001,
                          h0_gev=1.1*current.h0_gev)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/"input.json"
            write_param_snapshot(changed, path)
            loaded = load_model_params(path)
        self.assertEqual(loaded.h0_gev, changed.h0_gev)
        self.assertEqual(loaded.omega_b, changed.omega_b)
        self.assertEqual(loaded.source_run_id, current.source_run_id)

    def test_missing_required_input_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/"input.json"
            path.write_text(json.dumps({"omega_b": 0.05}), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "missing"):
                load_model_params(path)

    def test_nested_current_entropy_key_and_conflicts(self):
        current = load_model_params()
        names = {"m_p_gev": "M_P", "gamma_m_integral": "entropy_integral", "h0_gev": "H0_GEV",
                 "omega_lambda": "Omega_Lambda", "omega_b": "Omega_b", "omega_sigma": "Omega_Sigma",
                 "rho_lambda_gev4": "rho_Lambda"}
        records = {name: {"value": getattr(current, field), "run_id": "one"} for field, name in names.items()}
        records["entropy_integral"]["value"] += 0.1
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/"input.json"
            path.write_text(json.dumps({"parameters": records}), encoding="utf-8")
            self.assertAlmostEqual(load_model_params(path).s_abs, current.s_abs+0.1)
            records["ir_flow_int_gamma"] = {"value": current.s_abs, "run_id": "one"}
            path.write_text(json.dumps({"parameters": records}), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Conflicting"):
                load_model_params(path)

    def test_boolean_nested_quantity_rejected_and_massless_neutrino_allowed(self):
        current = replace(load_model_params(), m_nu1_ev=0)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/"input.json"
            write_param_snapshot(current, path)
            self.assertEqual(load_model_params(path).m_nu1_ev, 0)
            data = {key: {"value": value} for key, value in current.to_dict().items() if isinstance(value, (int,float))}
            data["m_p_gev"] = {"value": True}
            path.write_text(json.dumps({"parameters": data}), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "boolean"):
                load_model_params(path)


class CosmologicalClockTests(unittest.TestCase):
    def test_absent_matter_growing_branch_is_rejected(self):
        vacuum = replace(load_model_params(), omega_b=0, omega_sigma=0, omega_lambda=1)
        with self.assertRaisesRegex(ValueError, "positive total matter"):
            growing_rate(0.2, vacuum)

    def test_eds_coefficients_and_growth_refinement(self):
        params = replace(load_model_params(), omega_b=0, omega_sigma=1, omega_lambda=0)
        drag, force = coefficients(0.2, params)
        self.assertAlmostEqual(drag, 7.5)
        self.assertAlmostEqual(force, 37.5)
        errors = [abs(linear_growth_advance(0.1, 1, n, params, displacement=0.1, velocity=1)[0]-1) for n in (32,64)]
        self.assertGreater(errors[0]/errors[1], 3.0)

    def test_force_free_momentum_and_hubble_damping(self):
        params = load_model_params()
        state = {"positions": np.ones((1,3)), "velocities": np.ones((1,3))}
        config = PMNBodyConfig(dynamics="cosmological", steps_per_interval=16)
        with patch("cosmic_cg.pm_nbody._particle_acceleration", lambda *args: np.zeros((1,3))):
            _advance_to(state, 0.2, 0.4, config, {"params": params})
        expected = 0.2**3*expansion(0.2, params)/(0.4**3*expansion(0.4, params))
        np.testing.assert_allclose(state["velocities"], expected, rtol=1e-14)
        self.assertLess(expected, 1)


class ActiveEvolutionTests(unittest.TestCase):
    def test_nonrational_certificate_clock_is_rejected(self):
        with self.assertRaises(TypeError):
            sustained_active_run(step=0.0001)
        with self.assertRaises(TypeError):
            sustained_active_run(full_modes=True)

    def test_actual_quadratic_tail_from_a_single_mode(self):
        rhs, _ = matter_rhs([Q(1,100), -Q(1,100)], 1, 3)
        self.assertEqual(rhs[1], Q(1,10000))
        self.assertEqual(rhs[4], -Q(1,20000))
        self.assertEqual(rhs[2], 0)

    def test_production_certifies_without_reference_or_full_history(self):
        with patch("scipy.integrate.solve_ivp", side_effect=AssertionError("reference used in production")):
            result = sustained_active_run()
        self.assertTrue(result["all_checks_passed"])
        self.assertEqual(result["steps"], 128)
        self.assertEqual(result["final_active_modes"], 2)
        self.assertEqual(result["events"][0]["time"], 0.0)
        self.assertEqual(result["work"]["rejected_attempts"], 1)
        self.assertGreater(result["work"]["rejected_advancement_float"]["multiplications"], 0)
        self.assertEqual(result["work"]["records_and_storage"]["validation_state_scalars_stored"], 0)
        self.assertNotIn("_validation_history", result)
        self.assertNotIn("full_certified_reference", result)
        self.assertGreater(result["work"]["certificate"]["integer_quadratic_products"], 0)

    def test_shared_reference_checks_both_strict_certificates(self):
        result = verify_spectral_pair()
        full = result["full_certified_reference"]
        self.assertTrue(result["all_checks_passed"] and full["all_checks_passed"])
        self.assertEqual(result["steps"], full["steps"])
        self.assertEqual(result["tolerance"], full["tolerance"])
        self.assertLessEqual(result["max_sampled_error_fraction_of_bound"], 1.00000001)
        self.assertLessEqual(full["max_sampled_error_fraction_of_bound"], 1.00000001)
        self.assertEqual(result["verification_work"]["dense_reference_evaluations"], 1)
        self.assertLess(result["work"]["advancement_float"]["rhs_coordinates"],
                        full["work"]["advancement_float"]["rhs_coordinates"])
        self.assertNotIn("_validation_history", full)

    def test_budget_recurrence_is_a_selection_policy(self):
        budget = ReserveBudget(Q(100,99), 1e-8, 1e-4, 128)
        allowances = [budget.advance() for _ in range(128)]
        self.assertTrue(all(0 < value <= 1e-4 for value in allowances))
        self.assertEqual(allowances[-1], 1e-4)
        self.assertEqual(budget.work()["float_power_evaluations"], 1)


if __name__ == "__main__":
    unittest.main()
