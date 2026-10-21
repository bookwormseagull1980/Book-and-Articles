"""Primitive, end-to-end, and command-line checks for the certificate."""

import copy
import json
import subprocess
import sys
import tempfile
import unittest
from decimal import Context, Decimal
from fractions import Fraction
from pathlib import Path

import certify_endpoint_fixed_point as cert
from certify_endpoint_fixed_point import Interval, exp_point, pi_interval, sqrt_point, weight, weight_sum

HERE = Path(__file__).resolve().parent


class IntervalChecks(unittest.TestCase):
    def setUp(self):
        cert.configure_precision(220)

    def tearDown(self):
        cert.configure_precision(220)

    def assert_contains(self, interval, exact):
        self.assertLessEqual(Fraction(interval.lo), exact)
        self.assertGreaterEqual(Fraction(interval.hi), exact)

    def test_arithmetic_contains_exact_endpoint_extrema(self):
        samples = [Interval(Decimal(a), Decimal(b)) for a, b in
                   [("-3.17", "-1.01"), ("-0.001", "0.02"),
                    ("0.0000001", "3.7"), ("7.3333333333333333333", "8.2")]]
        for x in samples:
            for y in samples:
                for a in (Fraction(x.lo), Fraction(x.hi)):
                    for b in (Fraction(y.lo), Fraction(y.hi)):
                        self.assert_contains(x + y, a + b)
                        self.assert_contains(x - y, a - b)
                        self.assert_contains(x * y, a * b)
                        if not y.lo <= 0 <= y.hi:
                            self.assert_contains(x / y, a / b)
        self.assert_contains(Interval.point(1) / 7, Fraction(1, 7))
        with self.assertRaises(ZeroDivisionError):
            Interval.point(1) / samples[1]

    def test_square_root_brackets_by_exact_squaring(self):
        for text in ("0", "0.02", "2", "4", "123456789.12345", "1e-1001", "1e1001"):
            value = Decimal(text)
            root = sqrt_point(value)
            self.assertLessEqual(Fraction(root.lo) ** 2, Fraction(value))
            self.assertGreaterEqual(Fraction(root.hi) ** 2, Fraction(value))

    def test_pi_against_exact_rational_machin_sum(self):
        def rational_atan(q):
            total = sum((Fraction((-1) ** k, (2 * k + 1) * q ** (2 * k + 1))
                         for k in range(400)), Fraction(0))
            return total, total + Fraction(1, 801 * q ** 801)
        a_lo, a_hi = rational_atan(5)
        b_lo, b_hi = rational_atan(239)
        pi = pi_interval()
        self.assert_contains(pi, 16 * a_lo - 4 * b_hi)
        self.assert_contains(pi, 16 * a_hi - 4 * b_lo)

    def test_exp_against_separate_decimal_implementation(self):
        reference = Context(prec=300)
        for text in ("-12768", "-3978", "-1", "-1e-35", "0", "0.125", "4.14", "100"):
            value = Decimal(text)
            computed = reference.exp(value)
            interval = exp_point(value)
            self.assertLessEqual(interval.lo, computed)
            self.assertGreaterEqual(interval.hi, computed)

    def test_exp_against_exact_rational_taylor_bounds(self):
        for x in (Fraction(1, 8), Fraction(1), Fraction(4)):
            term = Fraction(1)
            total = term
            for n in range(1, 501):
                term *= x / n
                total += term
            remainder = term * x / 501 / (1 - x / 502)
            value = Decimal(x.numerator) / Decimal(x.denominator)
            interval = exp_point(value)
            self.assert_contains(interval, total)
            self.assert_contains(interval, total + remainder)

    def test_precision_change_discards_cached_intervals(self):
        old = exp_point(Decimal(1))
        self.assertGreater(exp_point.cache_info().currsize, 0)
        cert.configure_precision(260)
        self.assertEqual(exp_point.cache_info().currsize, 0)
        new = exp_point(Decimal(1))
        self.assertGreaterEqual(new.lo, old.lo)
        self.assertLessEqual(new.hi, old.hi)
        self.assertLess(Fraction(new.hi) - Fraction(new.lo),
                        Fraction(old.hi) - Fraction(old.lo))

    def test_exact_input_and_interval_validation(self):
        exp_point(Decimal(1))
        for value in (1.0, True, "NaN", "Infinity", "-Infinity"):
            with self.subTest(value=value):
                with self.assertRaises((ValueError, TypeError)):
                    Interval.point(value)
                with self.assertRaises((ValueError, TypeError)):
                    exp_point(value)
        with self.assertRaises(ValueError):
            Interval(Decimal(2), Decimal(1))
        with self.assertRaises(TypeError):
            Interval(0, 1)

    def test_precision_and_function_domains(self):
        for digits in (0, 219, 220.0, True):
            with self.assertRaises((ValueError, TypeError)):
                cert.configure_precision(digits)
        with self.assertRaises(ValueError):
            sqrt_point(-1)
        for denominator in (0, 1, -2, 5.0):
            with self.assertRaises(ValueError):
                cert.atan_inverse(denominator)
        self.assertEqual(exp_point(Decimal(0)), Interval.point(1))

    def test_tail_domain_guards(self):
        for cutoff, moment, beta in ((3, 0, 1), (4, 2, 1), (4, True, 1),
                                     (4, 0, 0), (4, 0, -1)):
            with self.assertRaises(ValueError):
                cert.gaussian_tail(beta, cutoff, moment=moment)
        with self.assertRaises(ArithmeticError):
            cert.gaussian_tail("0.001", 4)
        with self.assertRaises(ValueError):
            cert.certify(29)

    def test_tail_majorizes_finite_partial_sums(self):
        for moment in (0, 1):
            for shift in (0, 9):
                tail = cert.gaussian_tail(1, 4, moment=moment, shift=shift)
                partial = sum((k ** moment * weight_sum(k)
                               * exp_point(Decimal(-(k * k - shift)))
                               for k in range(5, 31)), Interval.point(0))
                self.assertGreater(tail.lo, 0)
                self.assertLessEqual(partial.hi, tail.hi)

    def test_divisor_sieve_against_direct_count(self):
        values = cert.decomposition_values(100)
        for n in range(1, 101):
            expected = sum(n % d == 0 for d in range(2, n))
            self.assertEqual(values[n], expected)

    def test_fibre_coefficients_against_rational_sum(self):
        for k in range(3, 36):
            exact_sum = Fraction(0)
            enclosed_sum = Interval.point(0)
            for m in range(1, k):
                exact = (2 * m * (k - m) * (1 - Fraction(2, k - 1))
                         + Fraction(k * (k * k - 1), 3 * (k - 1) ** 2))
                self.assertGreaterEqual(exact, 0)
                self.assert_contains(weight(k, m), exact)
                exact_sum += exact
                enclosed_sum = enclosed_sum + weight(k, m)
            expected = Fraction(k * (k + 1) * (k - 2), 3)
            self.assertEqual(exact_sum, expected)
            self.assert_contains(enclosed_sum, expected)
            self.assert_contains(weight_sum(k), expected)


class CertificateChecks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cert.configure_precision(220)
        cls.primary = cert.certify(30)
        cert.configure_precision(260)
        cls.check = cert.certify(32)
        cert.configure_precision(220)

    def test_full_certificates_and_cross_precision_overlap(self):
        for result in (self.primary, self.check):
            self.assertTrue(result["all_claims_verified"])
            self.assertTrue(all(result["claims"].values()))
        a = self.primary["intervals"]["fixed_point"]
        b = self.check["intervals"]["fixed_point"]
        self.assertLessEqual(max(Decimal(a["lower"]), Decimal(b["lower"])),
                             min(Decimal(a["upper"]), Decimal(b["upper"])))

    def test_archived_certificates_reproduce_every_field(self):
        for filename, computed in (("endpoint_fixed_point_v1.json", self.primary),
                                   ("endpoint_fixed_point_check_v1.json", self.check)):
            reference = cert.read_certificate(HERE / filename)
            cert.verify_recomputed(reference, computed)

    def test_changed_numerical_result_is_rejected(self):
        altered = copy.deepcopy(self.primary)
        altered["intervals"]["fixed_point"]["lower"] = "0"
        with self.assertRaisesRegex(ValueError, "intervals"):
            cert.verify_recomputed(altered, self.primary)

    def test_schema_and_claim_changes_are_rejected(self):
        altered = copy.deepcopy(self.primary)
        del altered["claims"]
        with self.assertRaisesRegex(ValueError, "schema"):
            cert.verify_recomputed(altered, self.primary)
        altered = copy.deepcopy(self.primary)
        altered["claims"]["absolute_error_lt_1e-50"] = False
        with self.assertRaisesRegex(ValueError, "claims"):
            cert.verify_recomputed(altered, self.primary)

    def test_runtime_metadata_can_differ(self):
        altered = copy.deepcopy(self.primary)
        altered["python_version"] = "another runtime"
        cert.verify_recomputed(altered, self.primary)

    def test_invalid_archives_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "invalid.json"
            for field, value in (("script_sha256", "0" * 64),
                                 ("certificate_version", "0.0.0"),
                                 ("all_claims_verified", False),
                                 ("decimal_precision", 20), ("truncation_N", 3)):
                altered = copy.deepcopy(self.primary)
                altered[field] = value
                path.write_text(json.dumps(altered), encoding="utf-8")
                with self.subTest(field=field), self.assertRaises(ValueError):
                    cert.read_certificate(path)
            path.write_text("[]", encoding="utf-8")
            with self.assertRaises(ValueError):
                cert.read_certificate(path)


class CommandLineChecks(unittest.TestCase):
    def run_cli(self, *args):
        return subprocess.run([sys.executable, "-B", str(HERE / "certify_endpoint_fixed_point.py"),
                               *args], text=True, capture_output=True, timeout=120)

    def test_help_and_version(self):
        self.assertEqual(self.run_cli("--help").returncode, 0)
        result = self.run_cli("--version")
        self.assertEqual(result.returncode, 0)
        self.assertEqual(result.stdout.strip(), cert.VERSION)

    def test_invalid_options_fail(self):
        for args in (("--digits", "219"), ("--cutoff", "29"),
                     ("--verify", "missing.json", "--digits", "220"),
                     ("--output", "a.json", "--verify", "b.json")):
            with self.subTest(args=args):
                self.assertEqual(self.run_cli(*args).returncode, 2)

    def test_output_never_overwrites_an_existing_file(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "keep.json"
            path.write_text("original contents", encoding="utf-8")
            result = self.run_cli("--output", str(path))
            self.assertEqual(result.returncode, 2)
            self.assertEqual(path.read_text(encoding="utf-8"), "original contents")

    def test_malformed_json_fails_cleanly(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "invalid.json"
            path.write_text("{invalid", encoding="utf-8")
            result = self.run_cli("--verify", str(path))
            self.assertEqual(result.returncode, 2)
            self.assertIn("ERROR:", result.stderr)
            self.assertNotIn("Traceback", result.stderr)


if __name__ == "__main__":
    unittest.main()
