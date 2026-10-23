import csv
import importlib.util
import tempfile
import unittest
from fractions import Fraction as Q
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def module(name):
    spec = importlib.util.spec_from_file_location(name, ROOT / "code" / f"{name}.py")
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


CERT = module("verify_finite_window_closure")
BUDGET = module("readout_object_budget")


class ProofClosureTests(unittest.TestCase):
    def test_exact_certificate(self):
        result = CERT.certificate()
        self.assertEqual(result["uniform_error_upper"], Q(15199313, 361200000000))
        self.assertTrue(result["tube_pass"])
        self.assertTrue(result["half_tolerance_pass"])

    def test_time_refinement_reduces_defect(self):
        coarse, fine = CERT.certificate(), CERT.certificate(step=Q(1, 20000))
        self.assertEqual(coarse["integration_defect_integral"],
                         2 * fine["integration_defect_integral"])
        self.assertLess(fine["uniform_error_upper"], coarse["uniform_error_upper"])

    def test_zero_step_rejected(self):
        with self.assertRaises(ValueError):
            CERT.certificate(step=Q(0))

    def test_laplace_boundary_identity(self):
        energy = [Q(2), Q(3), Q(5), Q(7)]
        z = Q(2, 3)
        ej = sum(e * z**(i+1) for i, e in enumerate(energy))
        enext = sum(e * z**(i+1) for i, e in enumerate(energy[1:]))
        self.assertEqual(enext, ej / z - energy[0])

    def test_soft_linear_tail_is_not_zero(self):
        x, selector = Q(2), Q(1, 2)
        self.assertEqual(selector * (x-selector*x), Q(1, 2))

    def test_moving_coordinate_uses_right_inverse(self):
        a, adot, full_error = Q(1, 2), Q(3), Q(2)
        active_error = a * full_error
        self.assertEqual(adot/a * active_error, adot * full_error)
        self.assertNotEqual(adot * active_error, adot * full_error)

    def test_missing_flux_column_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "missing.csv"
            with path.open("w", newline="", encoding="utf-8") as handle:
                csv.writer(handle).writerows([["snapshot", "a_pq"], ["s0", "1"]])
            with self.assertRaises(ValueError):
                BUDGET.rows(path, ["net_total"])

    def test_actual_particle_budget(self):
        tables = BUDGET.RUN / "tables"
        read = lambda name: BUDGET.rows(tables / name, [])
        result, growth = BUDGET.assemble(
            read("pm_snapshot_table.csv"), read("pm_halo_table.csv"), read("pm_tree_table.csv"))
        self.assertEqual(len(result), 24)
        self.assertEqual(result[-1]["object_particles"], 17050)
        self.assertEqual(len(growth), 401)
        for row in result:
            self.assertEqual(row["net_particle_error"], 0)
            self.assertEqual(row["object_particles"] + row["diffuse_particles"], 262144)
            self.assertEqual(row["reentry_status"], "not_recorded")
        for row in result[1:]:
            self.assertLessEqual(row["shared_lower"], row["shared_upper"])

    def test_growth_is_not_saved_fractional_proxy(self):
        import math
        self.assertNotEqual(math.log(2), Q(1, 2))

    def test_corrected_fourier_clock_and_gravity(self):
        import numpy as np
        import math
        integrator = module("run_resolved_spectral_integrator")
        integrator.GRID = 8
        integrator.BOX_MPCH = 2 * math.pi
        grid = integrator.build_grid()
        self.assertAlmostEqual(grid[0][1, 0, 0], 1)
        delta = np.zeros_like(grid[3], dtype=complex)
        theta = np.zeros_like(delta)
        delta[1, 0, 0] = delta[-1, 0, 0] = 1e-6
        dd, dt = integrator.rhs(delta, theta, 0.25, 1.0, 0.0, grid)
        # Hcal=2, d eta/d a=2, so gravitational coefficient is -12.
        self.assertAlmostEqual(dt[1, 0, 0].real, -12e-6)
        self.assertAlmostEqual(abs(dd[1, 0, 0]), 0)


if __name__ == "__main__":
    unittest.main()
