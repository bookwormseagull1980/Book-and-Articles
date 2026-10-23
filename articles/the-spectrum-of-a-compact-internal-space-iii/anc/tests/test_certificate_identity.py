"""Identity checks must never silently repeat a cosmological calculation."""
import math
import unittest
from dataclasses import replace
from pathlib import Path
from unittest.mock import patch

import verify_pm_certificate_identity as identity
import p9_cic_replacement_validation as cache
from cosmic_cg.params import load_model_params


class CertificateIdentityTests(unittest.TestCase):
    def test_saved_package(self):
        report = identity.verify()
        self.assertTrue(report["all_checks_passed"])
        self.assertEqual(report["new_force_evaluations"], 0)
        self.assertFalse(report["certificate_recomputed"])

    def test_one_ulp_input_change_rejected(self):
        params = load_model_params()
        with self.assertRaisesRegex(ValueError, "consumed numerical parameter"):
            identity.verify(replace(params, omega_b=math.nextafter(params.omega_b, math.inf)))

    def test_cache_corruption_rejected(self):
        original = identity.digest
        def changed(path):
            return "0" * 64 if Path(path).suffix == ".npz" else original(path)
        with patch.object(identity, "digest", side_effect=changed):
            with self.assertRaisesRegex(ValueError, "Published evidence changed"):
                identity.verify()

    def test_missing_cache_does_not_advance(self):
        original = Path.exists
        def exists(path):
            return False if path.name == "p9_actual_trajectory_20261004.npz" else original(path)
        with patch.object(Path, "exists", exists), patch.object(cache, "_initial_state") as advance:
            with self.assertRaisesRegex(FileNotFoundError, "requires the saved real trajectory"):
                cache.original_trajectory()
            advance.assert_not_called()

    def test_changed_source_identity_rejected(self):
        with patch.object(identity, "digest", return_value="0" * 64):
            with self.assertRaisesRegex(ValueError, "Deployed source differs"):
                identity.source_record_matches("0" * 64, identity.ANC / "code/p9_continuous_residual_certificate.py")

    def test_changed_formal_result_rejected(self):
        with patch.object(identity, "numerical_payload_digest", return_value="0" * 64):
            with self.assertRaisesRegex(ValueError, "saved formal numerical result changed"):
                identity.verify()


if __name__ == "__main__":
    unittest.main()
