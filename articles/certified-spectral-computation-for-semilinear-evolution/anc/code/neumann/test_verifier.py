"""Acceptance, malformed-record, and mathematical residual checks."""

from copy import deepcopy
from fractions import Fraction
import json
from pathlib import Path
import re
import tempfile
import unittest
from unittest.mock import patch

from format_table import power_of_ten, render_rows, upper_scientific
from verify import Rejected, verify
from produce import make_record
from record_format import (DEFAULT_INITIAL, format_limits, guaranteed_parameters,
                           input_from_coefficients, load_record, normalize_input)


RESULTS = Path(__file__).resolve().parents[2] / "results" / "neumann"


class RecordTests(unittest.TestCase):
    def setUp(self):
        self.record = json.loads((RESULTS / "m4.record.json").read_text(encoding="utf-8"))

    def constant_record(self):
        changed = deepcopy(self.record)
        initial = changed["pieces"][0][0]
        for piece in changed["pieces"]:
            piece[0] = initial.copy()
            for row in piece[1:]:
                row[:] = ["0"] * len(row)
        return changed

    def test_published_records(self):
        for m in (4, 6, 8):
            record = json.loads((RESULTS / f"m{m}.record.json").read_text(encoding="utf-8"))
            self.assertEqual(verify(record, m)["status"], "ACCEPT")

    def test_manual_constant_residual(self):
        changed = self.constant_record()
        changed["m"] = 1
        result = verify(changed, 1)
        # For the stated initial polynomial, the Alin norm of A u0-u0^2
        # is 529/8192 + 2(31/1024) + 3(27/512) + 4/1024 + 5/8192.
        self.assertEqual(Fraction(result["max_residual_Alin"]["exact"]),
                         Fraction(1179, 4096))

    def test_no_recurrence_requirement(self):
        changed = deepcopy(self.record)
        old = int(changed["pieces"][-1][-1][0])
        changed["pieces"][-1][-1][0] = str(old + 1)
        self.assertEqual(verify(changed, 4)["status"], "ACCEPT")

    def test_large_residual_rejected(self):
        changed = self.constant_record()
        changed["claimed_total_error"] = "0"
        with self.assertRaisesRegex(Rejected, "total exceeds"):
            verify(changed, 4)

    def test_discontinuity_rejected(self):
        changed = deepcopy(self.record)
        changed["pieces"][1][0][0] = str(int(changed["pieces"][1][0][0]) + 1)
        with self.assertRaisesRegex(Rejected, "discontinuous"):
            verify(changed, 4)

    def test_norm_rejected(self):
        changed = deepcopy(self.record)
        changed["pieces"][0][0][0] = str(4 << changed["fractional_bits"])
        with self.assertRaisesRegex(Rejected, "magnitude bound|A4 budget"):
            verify(changed, 4)

    def test_precision_rejected(self):
        with self.assertRaisesRegex(Rejected, "precision request mismatch"):
            verify(self.record, 8)

    def test_initial_data_rejected(self):
        changed = deepcopy(self.record)
        changed["input_reference"][0] = "1/8"
        with self.assertRaisesRegex(Rejected, "initial input mismatch"):
            verify(changed, 4)

    def test_equation_rejected(self):
        changed = deepcopy(self.record)
        changed["equation_id"] = "u_t=u_xx-u^2; Neumann [0,pi]"
        with self.assertRaisesRegex(Rejected, "wrong equation"):
            verify(changed, 4)

    def test_malformed_coefficient_rejected(self):
        changed = deepcopy(self.record)
        changed["pieces"][0][0][0] = "1/2"
        with self.assertRaisesRegex(Rejected, "integer string"):
            verify(changed, 4)

    def test_missing_mode_rejected(self):
        changed = deepcopy(self.record)
        changed["pieces"][0][0].pop()
        with self.assertRaisesRegex(Rejected, "mode count"):
            verify(changed, 4)

    def test_external_input_binding(self):
        other = input_from_coefficients([Fraction(1,8)])
        with self.assertRaisesRegex(Rejected, "initial input mismatch"):
            verify(self.record, 4, other)

    def test_dense_initial_tail_and_quantization(self):
        coefficients = [Fraction(1,3)] + [Fraction(0)]*11 + [Fraction(1, 4**13)]
        external = input_from_coefficients(coefficients)
        record = make_record(4, initial=external)
        result = verify(record, 4, external)
        self.assertGreater(len(coefficients)-1, record["cutoff"])
        self.assertGreater(Fraction(result["initial_error_Alin"]["exact"]), 0)
        self.assertLess(Fraction(result["total_error"]["exact"]), Fraction(1,64))

    def test_guaranteed_parameter_hand_values(self):
        p = guaranteed_parameters(1)
        self.assertEqual(tuple(p[k] for k in ("cutoff","steps","degree","fractional_bits")),
                         (4,1024,75,96))
        self.assertEqual(p["xi"], Fraction(1,19140298416324608))

    def test_guaranteed_records_and_exact_constant_solution(self):
        summary = json.loads((RESULTS/"matrix_summary.json").read_text(encoding="utf-8"))
        rows = [row for row in summary["rows"] if row["mode"] == "guaranteed"]
        self.assertEqual(len(rows),3)
        for row in rows:
            external = json.loads((RESULTS/"matrix"/row["input_file"]).read_text(encoding="utf-8"))
            path = RESULTS/"matrix"/row["record_file"]
            result = verify(load_record(path,external,row["m"]),row["m"],external)
            checks = result["theorem_strong_checks"]
            self.assertTrue(checks["initial_below_xi"])
            self.assertTrue(checks["residual_below_xi"])
            if row["case"] == "E2-constant":
                center = Fraction(result["endpoint_observations"]["mean"]["exact"])
                self.assertLessEqual(abs(center-Fraction(16,31)),
                                     Fraction(result["time_and_rounding_error"]["exact"]))

    def test_false_guaranteed_header_rejected(self):
        changed = deepcopy(self.record)
        changed["mode"] = "guaranteed"
        with self.assertRaisesRegex(Rejected, "theorem parameter mismatch"):
            verify(changed,4)

    def test_large_headers_rejected_before_coefficients(self):
        for key in ("cutoff", "steps", "degree", "fractional_bits"):
            changed = deepcopy(self.record)
            changed[key] = 10**100
            changed["pieces"] = None
            with self.assertRaisesRegex(Rejected, "parameter bound"):
                verify(changed, 4)

    def test_large_coefficient_rejected_before_integer_conversion(self):
        changed = deepcopy(self.record)
        changed["pieces"][0][0][0] = "9"*10000
        with self.assertRaisesRegex(Rejected, "string length"):
            verify(changed, 4)

    def test_file_preflight_and_json_header_limit(self):
        external = input_from_coefficients(DEFAULT_INITIAL)
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/"record.json"
            path.write_text('{"fractional_bits":'+"9"*1000+'}', encoding="utf-8")
            with self.assertRaisesRegex(Rejected, "oversized JSON integer"):
                load_record(path, external, 4)
            # The size rejection happens before opening a record body.
            with patch.object(Path, "stat", return_value=type("S",(),{"st_size":10**30})()):
                with self.assertRaisesRegex(Rejected, "byte bound"):
                    load_record(path, external, 4)

    def test_input_budget_rejected(self):
        with self.assertRaisesRegex(Rejected, "initial A4 budget"):
            input_from_coefficients([Fraction(1),Fraction(1,16)])

    def test_constant_solution_exact_enclosure(self):
        external = input_from_coefficients([Fraction(1,2)])
        result = verify(make_record(4,initial=external),4,external)
        exact = Fraction(16,31)
        center = Fraction(result["endpoint_observations"]["mean"]["exact"])
        self.assertLessEqual(abs(center-exact), Fraction(result["total_error"]["exact"]))

    def test_upward_table_format(self):
        summary = json.loads((RESULTS / "summary.json").read_text(encoding="utf-8"))
        values = [Fraction(0), Fraction(1, 3), Fraction(99999, 10000),
                  Fraction(1, 10**30), Fraction(100001, 10000)]
        values += [Fraction(row[key]["exact"]) for row in summary["rows"]
                   for key in ("max_residual_Alin", "total_error")]
        for value in values:
            text, rounded = upper_scientific(value)
            self.assertGreaterEqual(rounded, value)
            if value:
                match = re.fullmatch(r"\$([0-9]\.[0-9]{3})\\times10\^\{(-?[0-9]+)\}\$", text)
                self.assertIsNotNone(match)
                reconstructed = Fraction(match.group(1)) * power_of_ten(int(match.group(2)))
                self.assertEqual(reconstructed, rounded)
                self.assertLess(rounded - value, power_of_ten(int(match.group(2)) - 3))
        self.assertEqual(upper_scientific(Fraction(99999, 10000))[0],
                         r"$1.000\times10^{1}$")
        lines = render_rows(summary).splitlines()
        self.assertEqual(lines[1], r"\newcommand{\NeumannResultRows}{%")
        self.assertEqual(lines[-1], "}")
        rows = [line for line in lines if " & " in line]
        self.assertEqual(len(rows), 3)
        self.assertTrue(all(len(row.split(" & ")) == 7 and row.endswith(r"\\") for row in rows))


if __name__ == "__main__":
    unittest.main(verbosity=2)
