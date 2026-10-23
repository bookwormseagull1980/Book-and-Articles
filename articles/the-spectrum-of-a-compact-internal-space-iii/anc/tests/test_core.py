from pathlib import Path
import math
import unittest

from cosmic_cg.entropy_audit import audit_entropy, entropy_deficit, window_entropy
from cosmic_cg.generation import (
    InternalSpectrumConfig,
    finite_window_envelope,
    generate_epoch,
    generate_objects,
    growth_factor,
    internal_spectrum_shape,
)
from cosmic_cg.field3d import Field3DConfig, estimate_field_resources, generate_field_history
from cosmic_cg.interactions import INTERACTION_CHANNELS, generate_interaction_outputs
from cosmic_cg.io.schema_catalog import schema_catalog_rows
from cosmic_cg.knl import (
    band_power_fraction,
    dimensionless_power,
    find_knl,
    residual_floor_fraction,
    tail_power_fraction,
    tail_power_moment_fraction,
)
from cosmic_cg.object_level import alignment_cosine, halo_fate, projected_torque, tidal_torque
from cosmic_cg.params import DEFAULT_PARAMS, load_model_params
from cosmic_cg.pm_nbody import PMNBodyConfig, estimate_pm_resources, generate_pm_history
from cosmic_cg.result_tables import band_table_rows, endpoint_mass_map, object_summary_rows
from cosmic_cg.run_plan import generation_redshifts, production_run_grid
from cosmic_cg.tail_ratio import halo_tail_ratio, virial_velocity


class KNLTests(unittest.TestCase):
    def test_dimensionless_power_and_crossing(self):
        k = [0.1, 0.2, 0.4]
        p = [3947.841760435743, 2467.4011002723395, 925.275412602518]
        delta2 = dimensionless_power(k, p)
        self.assertAlmostEqual(delta2[1], 1.0, places=12)
        result = find_knl(k, p, c_nl=1.0)
        self.assertAlmostEqual(result.k_nl_hmpc, 0.2, places=12)

    def test_tail_power_fraction(self):
        k = [0.1, 0.2, 0.4, 1.0]
        p = [3947.841760435743, 2467.4011002723395, 925.275412602518, 118.4352528130723]
        value = tail_power_fraction(k, p, 0.2)
        self.assertGreater(value, 0.0)
        self.assertLess(value, 1.0)
        band = band_power_fraction(k, p, 0.2, 0.4)
        moment = tail_power_moment_fraction(k, p, 0.2, s=1.0)
        floor = residual_floor_fraction(k, p, 0.2)
        self.assertGreater(band, 0.0)
        self.assertLessEqual(moment, value)
        self.assertGreaterEqual(floor, 0.0)
        self.assertLessEqual(floor, value)


class EntropyTests(unittest.TestCase):
    def test_entropy_audit(self):
        self.assertAlmostEqual(window_entropy(0.2, 0.05), math.log(4.0))
        self.assertAlmostEqual(entropy_deficit(0.2, 0.1), math.log(2.0))
        audit = audit_entropy(s_abs=139.253706, k_nl_hmpc=0.2, k_box_hmpc=0.05, k_cut_hmpc=0.3)
        self.assertEqual(audit.d_s, 0.0)


class TailRatioTests(unittest.TestCase):
    def test_halo_tail_ratio(self):
        halos = [
            {"mass_msun_h": 1.0e12, "radius_kpc_h": 220.0, "sigma_v_kms": 95.0},
            {"mass_msun_h": 3.0e12, "radius_kpc_h": 310.0, "sigma_v_kms": 150.0},
        ]
        result = halo_tail_ratio(halos, m_nl_msun_h=1.0e13, r_nl_kpc_h=500.0)
        self.assertGreater(result.r_h, 0.0)
        self.assertEqual(result.n_halos, 2)
        self.assertGreater(virial_velocity(1.0e12, 220.0), 0.0)


class ObjectLevelTests(unittest.TestCase):
    def test_halo_fate_classes(self):
        forming = halo_fate(halo_id="a", mass_msun_h=1.0e12, radius_kpc_h=200.0, delta_mass_msun_h=5.0e11, delta_time_gyr=0.1)
        stripping = halo_fate(halo_id="b", mass_msun_h=1.0e12, radius_kpc_h=200.0, delta_mass_msun_h=-1.0e11, delta_time_gyr=0.1)
        self.assertEqual(forming.class_name, "forming")
        self.assertEqual(stripping.class_name, "stripping")

    def test_tidal_torque_projection(self):
        inertia = [[2.0, 0.1, 0.0], [0.1, 1.0, 0.0], [0.0, 0.0, 0.5]]
        tidal = [[0.0, 1.0, 0.0], [1.0, 0.0, 0.2], [0.0, 0.2, 0.0]]
        torque = tidal_torque(inertia, tidal)
        self.assertEqual(len(torque), 3)
        self.assertAlmostEqual(alignment_cosine([1, 0, 0], [1, 0, 0]), 1.0)
        self.assertIsInstance(projected_torque(torque, [1, 0, 0]), float)


class ResultTableTests(unittest.TestCase):
    def test_band_rows_and_schema_catalog(self):
        k = [0.1, 0.2, 0.4, 1.0]
        p = [3947.841760435743, 2467.4011002723395, 925.275412602518, 118.4352528130723]
        global_row = {
            "run_id": "test",
            "source": "fixture",
            "snapshot": "000",
            "z": 0.0,
            "a": 1.0,
            "k_box_hmpc": 0.05,
            "k_NL_hmpc": 0.2,
            "M_NL_msun_h": 1.0e13,
            "R_NL_kpc_h": 500.0,
        }
        mapper = endpoint_mass_map(1.0e13, 500.0, 0.2)
        rows = band_table_rows(global_row, k, p)
        self.assertGreater(mapper.rho_m0_msunh_mpch3, 0.0)
        self.assertEqual(len(rows), 5)
        self.assertTrue(any(row["activity_status"] == "frontier" for row in rows))
        self.assertTrue(any(row["table_name"] == "band_table.csv" for row in schema_catalog_rows()))

    def test_object_summary_rows_save_sample_counts(self):
        rows = object_summary_rows(
            [
                {
                    "run_id": "test",
                    "source": "fixture",
                    "snapshot": "000",
                    "z": 0.0,
                    "mass_msun_h": 1.0e12,
                    "B_h": 1.2,
                    "B_class": "forming",
                    "r_h_object": 0.4,
                    "quality_flag": "ok",
                },
                {
                    "run_id": "test",
                    "source": "fixture",
                    "snapshot": "000",
                    "z": 0.0,
                    "mass_msun_h": 3.0e12,
                    "B_h": -0.1,
                    "B_class": "stripping",
                    "r_h_object": 0.5,
                    "quality_flag": "ok",
                },
            ]
        )
        self.assertEqual(rows[0]["n_objects"], 2)
        self.assertAlmostEqual(rows[0]["forming_fraction"] + rows[0]["stripping_fraction"], 1.0)


class ParamsTests(unittest.TestCase):
    def test_defaults(self):
        params = load_model_params(None)
        self.assertAlmostEqual(params.s_abs, DEFAULT_PARAMS.gamma_m_integral)
        self.assertGreater(params.sigma0_over_sigmac, 0.0)


class RunPlanTests(unittest.TestCase):
    def test_generation_grid_has_past_present_future(self):
        redshifts = generation_redshifts()
        self.assertTrue(any(z > 0 for z in redshifts))
        self.assertIn(0.0, redshifts)
        self.assertTrue(any(z < 0 for z in redshifts))

    def test_production_grid_covers_required_phases(self):
        rows = production_run_grid(objects_per_band=10)
        phases = {row["phase"] for row in rows}
        self.assertIn("window_history", phases)
        self.assertIn("object_history", phases)
        self.assertEqual(sum(1 for row in rows if row["phase"] == "window_history"), len(generation_redshifts()))


class GenerationTests(unittest.TestCase):
    def test_growth_factor_extends_to_future(self):
        self.assertLess(growth_factor(1.0, params=DEFAULT_PARAMS), 1.0)
        self.assertAlmostEqual(growth_factor(0.0, params=DEFAULT_PARAMS), 1.0)
        self.assertGreater(growth_factor(-0.5, params=DEFAULT_PARAMS), 1.0)
        self.assertGreater(growth_factor(-2.0 / 3.0, params=DEFAULT_PARAMS), growth_factor(-0.5, params=DEFAULT_PARAMS))

    def test_internal_spectrum_shape_has_high_k_envelope(self):
        config = InternalSpectrumConfig()
        low = internal_spectrum_shape(0.2, config=config)
        high = internal_spectrum_shape(200.0, config=config)
        power_only = InternalSpectrumConfig(spectrum_shape="power_law")
        self.assertGreater(finite_window_envelope(0.2, k_uv_hmpc=20.0, alpha=2.0, beta=2.0), 0.0)
        self.assertLess(high / low, internal_spectrum_shape(200.0, config=power_only) / internal_spectrum_shape(0.2, config=power_only))

    def test_generated_epoch_crosses_endpoint(self):
        epoch = generate_epoch(z=0.0, params=DEFAULT_PARAMS, config=InternalSpectrumConfig(k_points=128))
        self.assertAlmostEqual(epoch.knl.k_nl_hmpc, 0.2, places=2)
        self.assertGreater(epoch.m_nl_msun_h, 0.0)

    def test_generated_objects_include_branch_and_pair_transfers(self):
        epoch = generate_epoch(z=0.0, params=DEFAULT_PARAMS, config=InternalSpectrumConfig(k_points=128))
        objects, links, transfers = generate_objects(
            run_id="test",
            epoch=epoch,
            params=DEFAULT_PARAMS,
            objects_per_band=3,
            object_step_gyr=0.2,
        )
        channels = {row["source_channel"] for row in transfers}
        self.assertEqual(len(objects), 15)
        self.assertEqual(len(links), 15)
        self.assertIn("generated_branch_transfer", channels)
        self.assertIn("generated_pair_tidal_transfer", channels)
        self.assertTrue(all(0.0 <= float(row["a_pq"]) for row in transfers))
        self.assertTrue(all(-1.0 <= float(row["g_pq"]) <= 1.0 for row in transfers))


class Field3DTests(unittest.TestCase):
    def test_field_history_writes_relationship_layers(self):
        outputs = generate_field_history(
            run_id="field_test",
            redshifts=[0.0],
            params=DEFAULT_PARAMS,
            spectrum_config=InternalSpectrumConfig(k_points=128),
            field_config=Field3DConfig(
                grid_size=16,
                mode_bins=8,
                max_peaks_per_snapshot=6,
                max_peak_candidates=64,
                peak_slab_cells=4,
                neighbor_links=2,
            ),
        )
        self.assertEqual(len(outputs["window"]), 1)
        self.assertEqual(len(outputs["mode"]), 8)
        self.assertGreater(len(outputs["object"]), 0)
        self.assertGreater(len(outputs["transfer"]), 0)
        self.assertEqual(len(outputs["closure"]), 1)
        window = outputs["window"][0]
        self.assertEqual(window["grid_size"], 16)
        self.assertEqual(window["gaussian_kernel"], "exp[-0.5*(sigma*k)^2]")
        self.assertAlmostEqual(window["activity_endpoint_hmpc"], window["k_nl_hmpc"])
        self.assertGreaterEqual(window["closure_error"], 0.0)
        self.assertGreater(window["measured_sigma_delta"], 0.0)
        self.assertEqual(window["quality_flag"], "relationship_tables_only")
        mode = outputs["mode"][0]
        self.assertIn("coarse_delta2_mean", mode)
        self.assertIn("gaussian_weight_mean", mode)
        self.assertIn("activity_weight_mean", mode)
        self.assertTrue(all(float(row["a_pq"]) >= 0.0 for row in outputs["transfer"]))
        self.assertTrue(all(-1.0 <= float(row["g_pq"]) <= 1.0 for row in outputs["transfer"]))

    def test_field_resource_estimate_tracks_raw_arrays(self):
        estimate = estimate_field_resources(redshift_count=24, field_config=Field3DConfig(grid_size=512))
        self.assertEqual(estimate["voxel_count"], 512**3)
        self.assertGreater(estimate["raw_density_bytes_all_snapshots"], estimate["raw_density_bytes_per_snapshot"])
        self.assertGreater(estimate["estimated_peak_memory_gb"], 0.0)


class PMNBodyTests(unittest.TestCase):
    def test_pm_history_tracks_particles_and_object_links(self):
        config = PMNBodyConfig(
            particles_per_dim=8,
            mesh_size=16,
            box_size_mpch=32.0,
            z_initial=5.0,
            steps_per_interval=1,
            mode_bins=8,
            density_threshold_sigma=0.25,
            density_threshold_delta=0.0,
            fof_min_particles=2,
            max_halos_per_snapshot=8,
            neighbor_links=1,
        )
        outputs = generate_pm_history(
            run_id="pm_test",
            redshifts=[3.0, 1.0, 0.0],
            params=DEFAULT_PARAMS,
            config=config,
        )
        self.assertEqual(len(outputs["snapshot"]), 3)
        self.assertEqual(len(outputs["power"]), 24)
        self.assertGreater(len(outputs["halo"]), 0)
        self.assertGreater(len(outputs["tree"]), 0)
        self.assertGreater(len(outputs["transfer"]), 0)
        self.assertTrue(all(abs(float(row["mass_conservation_error"])) < 1.0e-5 for row in outputs["snapshot"]))
        self.assertTrue(all(0.0 <= float(row["overlap_fraction"]) <= 1.0 for row in outputs["tree"]))
        self.assertTrue(all(float(row["a_pq"]) >= 0.0 for row in outputs["transfer"]))
        self.assertTrue(all(-1.0 <= float(row["g_pq"]) <= 1.0 for row in outputs["transfer"]))

    def test_pm_resource_estimate_tracks_particle_and_mesh_arrays(self):
        config = PMNBodyConfig(particles_per_dim=16, mesh_size=32)
        estimate = estimate_pm_resources(redshift_count=4, config=config)
        self.assertEqual(estimate["particle_count"], 16**3)
        self.assertEqual(estimate["mesh_cell_count"], 32**3)
        self.assertGreater(estimate["raw_bytes_all_snapshots"], estimate["raw_particle_bytes_per_snapshot"])
        self.assertGreater(estimate["estimated_peak_memory_gb"], 0.0)


class InteractionTests(unittest.TestCase):
    def test_interaction_strengths_are_normalized_by_state_and_band(self):
        global_rows = [
            {
                "run_id": "test",
                "source": "internal_generation",
                "snapshot": "a1",
                "z": 0.0,
                "a": 1.0,
                "k_NL_hmpc": 0.2,
                "M_NL_msun_h": 1.0e13,
                "Omega_tail": 0.8,
                "Omega_excess": 0.4,
                "R_h_mean": 0.5,
                "v_lambda": -0.2,
            }
        ]
        band_rows = [
            {
                "z": 0.0,
                "a": 1.0,
                "band_id": "galaxy_halo",
                "k_min_hmpc": 0.1,
                "k_max_hmpc": 0.3,
                "omega_band": 0.2,
                "activity_status": "frontier",
            }
        ]
        object_rows = [
            {
                "z": 0.0,
                "a": 1.0,
                "mass_msunh": 1.0e12,
                "b_h": 0.4,
            }
        ]
        transfer_rows = [
            {
                "z": 0.0,
                "a": 1.0,
                "a_pq": 0.1,
            }
        ]
        outputs = generate_interaction_outputs(
            run_id="test",
            global_rows=global_rows,
            band_rows=band_rows,
            object_rows=object_rows,
            transfer_rows=transfer_rows,
        )
        strengths = [
            row
            for row in outputs["strength"]
            if float(row["z"]) == 0.0 and row["band_id"] == "galaxy_halo"
        ]
        dominance = [
            row
            for row in outputs["dominance"]
            if float(row["z"]) == 0.0 and row["band_id"] == "galaxy_halo"
        ][0]
        self.assertEqual(len(strengths), len(INTERACTION_CHANNELS))
        self.assertAlmostEqual(sum(float(row["chi"]) for row in strengths), 1.0)
        strongest = max(strengths, key=lambda row: float(row["chi"]))
        self.assertEqual(dominance["dominant_channel"], strongest["channel_id"])


if __name__ == "__main__":
    unittest.main()
