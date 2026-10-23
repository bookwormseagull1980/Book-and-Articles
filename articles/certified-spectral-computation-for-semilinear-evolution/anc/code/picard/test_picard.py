"""Mathematical, replay, input-binding, and preflight checks."""

from copy import deepcopy
from fractions import Fraction as F
import json
from pathlib import Path
import re
import tempfile
import unittest
from unittest.mock import patch
from benchmarks import coefficient_error, exact_nilpotent_solution
from produce import construct
from record_format import (Rejected, load_record, make_input, normalize_input)
from run_examples import cases, upward
from verify import verify

RESULTS = Path(__file__).resolve().parents[2]/"results"/"picard"


class PicardTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        matrix, polynomial, _ = cases()["P3-nilpotent2"]
        cls.external = make_input(matrix,polynomial)
        cls.original = construct(cls.external,6)

    def setUp(self):
        self.record = deepcopy(self.original)

    def test_nonzero_coupled_exact_solution(self):
        report = verify(self.record,self.external,6)
        self.assertEqual(report["terminal_centers"],["153/512","-3/32"])
        _,matrix,polynomial,_ = normalize_input(self.external)
        exact = exact_nilpotent_solution(matrix,polynomial)
        self.assertEqual(coefficient_error(self.record,exact),0)

    def test_zero_matrix_time_dependent_solution(self):
        matrix, polynomial, _ = cases()["P2-forcing"]
        external = make_input(matrix,polynomial)
        record = construct(external,1)
        report = verify(record,external,1)
        _,matrix,polynomial,_ = normalize_input(external)
        self.assertEqual(coefficient_error(record,polynomial),0)
        self.assertEqual(report["status"],"ACCEPT")

    def test_high_degree_uses_weak_step_error(self):
        matrix, polynomial, _ = cases()["P8-highdegree"]
        external = make_input(matrix,polynomial)
        report = verify(construct(external,20),external,20)
        budget = F(1,1 << report["working_bits"])
        self.assertLessEqual(F(report["max_weak_step_error"]["exact"]),budget)
        self.assertGreater(F(report["max_strong_step_difference"]["exact"]),budget)
        self.assertLessEqual(F(report["max_strong_norm"]["exact"]),1)

    def test_every_published_family(self):
        summary = json.loads((RESULTS/"summary.json").read_text(encoding="utf-8"))
        rows = [row for row in summary["rows"] if row["m"] == 6]
        self.assertEqual(len(rows),8)
        for row in rows:
            external = json.loads((RESULTS/row["input_file"]).read_text(encoding="utf-8"))
            record = load_record(RESULTS/row["record_file"],external,6)
            self.assertEqual(verify(record,external,6)["status"],"ACCEPT")

    def test_tampered_update_rejected(self):
        self.record["iterates"][3][1][0] = str(int(self.record["iterates"][3][1][0])+1)
        with self.assertRaisesRegex(Rejected,"update does not match"):
            verify(self.record,self.external,6)

    def test_initial_iterate_must_be_zero(self):
        self.record["iterates"][0][0][0] = "1"
        with self.assertRaisesRegex(Rejected,"initial iterate"):
            verify(self.record,self.external,6)

    def test_terminal_must_be_exact_sum(self):
        self.record["terminal"][0] = str(int(self.record["terminal"][0])+1)
        with self.assertRaisesRegex(Rejected,"terminal value mismatch"):
            verify(self.record,self.external,6)

    def test_external_input_binding(self):
        changed = deepcopy(self.external)
        changed["polynomial_p"][0][0] = "1/8"
        with self.assertRaisesRegex(Rejected,"external input mismatch"):
            verify(self.record,changed,6)

    def test_precision_and_equation_binding(self):
        with self.assertRaisesRegex(Rejected,"parameter mismatch: m"):
            verify(self.record,self.external,7)
        self.record["equation_id"] = "different equation"
        with self.assertRaisesRegex(Rejected,"wrong equation"):
            verify(self.record,self.external,6)

    def test_boolean_and_huge_headers(self):
        for key in ("m","dimension","degree","fractional_bits","iterations"):
            for value in (True,10**100):
                changed = deepcopy(self.record)
                changed[key] = value
                changed["iterates"] = None
                with self.assertRaisesRegex(Rejected,"parameter mismatch"):
                    verify(changed,self.external,6)

    def test_overlong_coefficient_is_rejected_before_conversion(self):
        self.record["iterates"][1][0][0] = "9"*10000
        with self.assertRaisesRegex(Rejected,"string length"):
            verify(self.record,self.external,6)

    def test_bad_coefficient_and_ragged_array(self):
        self.record["iterates"][1][0][0] = "1e50"
        with self.assertRaisesRegex(Rejected,"integer string"):
            verify(self.record,self.external,6)
        self.record = deepcopy(self.original)
        self.record["iterates"][1][0].pop()
        with self.assertRaisesRegex(Rejected,"dimension mismatch"):
            verify(self.record,self.external,6)

    def test_input_budgets_and_dimensions(self):
        with self.assertRaisesRegex(Rejected,"row-sum budget"):
            make_input([[F(1,4)]],[[F(1,4)]])
        with self.assertRaisesRegex(Rejected,"source budget"):
            make_input([[0]],[[F(1,2)],[F(1,8)]])
        bad = deepcopy(self.external)
        bad["matrix_B"][1].pop()
        with self.assertRaisesRegex(Rejected,"square"):
            normalize_input(bad)

    def test_file_and_header_preflight(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory)/"record.json"
            path.write_text('{"degree":'+"9"*1000+'}',encoding="utf-8")
            with self.assertRaisesRegex(Rejected,"oversized JSON header"):
                load_record(path,self.external,6)
            with patch.object(Path,"stat",return_value=type("S",(),{"st_size":10**30})()):
                with self.assertRaisesRegex(Rejected,"record byte bound"):
                    load_record(path,self.external,6)

    def test_outward_table_rounding(self):
        summary = json.loads((RESULTS/"summary.json").read_text(encoding="utf-8"))
        values = [F(1,3),F(99999,10000),F(1,10**40)]
        values += [F(row["certified_path_error"]["exact"]) for row in summary["rows"]]
        for value in values:
            match = re.fullmatch(r"\$([0-9]\.[0-9]{3})\\times10\^\{(-?[0-9]+)\}\$",upward(value))
            self.assertIsNotNone(match)
            rounded = F(match.group(1))*F(10)**int(match.group(2))
            self.assertGreaterEqual(rounded,value)


if __name__ == "__main__":
    unittest.main(verbosity=2)
