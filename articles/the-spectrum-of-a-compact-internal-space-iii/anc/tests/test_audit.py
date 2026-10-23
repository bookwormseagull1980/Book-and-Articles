import json
from pathlib import Path
import tempfile
import unittest

from cosmic_cg.field3d import (
    FIELD_CLOSURE_COLUMNS,
    FIELD_COST_COLUMNS,
    FIELD_EPOCH_RELATION_COLUMNS,
    FIELD_MODE_COLUMNS,
    FIELD_OBJECT_COLUMNS,
    FIELD_RELATION_SUMMARY_COLUMNS,
    FIELD_TRANSFER_COLUMNS,
    FIELD_WINDOW_COLUMNS,
    Field3DConfig,
    field_cost_row,
    generate_field_history,
)
from cosmic_cg.generation import InternalSpectrumConfig
from cosmic_cg.io.audit import audit_run, detect_run_type
from cosmic_cg.io.field_audit import audit_field_run
from cosmic_cg.io.figure_data import write_field_figure_data
from cosmic_cg.io.pm_audit import audit_pm_run
from cosmic_cg.io.tables import write_rows
from cosmic_cg.params import DEFAULT_PARAMS
from cosmic_cg.pm_nbody import (
    PM_COST_COLUMNS,
    PM_HALO_COLUMNS,
    PM_POWER_COLUMNS,
    PM_SNAPSHOT_COLUMNS,
    PM_TRANSFER_COLUMNS,
    PM_TREE_COLUMNS,
    PMNBodyConfig,
    generate_pm_history,
    pm_cost_row,
)


class AuditTests(unittest.TestCase):
    def test_detects_generated_universe(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "generated"
            _write_generated_run(root)
            self.assertEqual(detect_run_type(root), "generated_universe")

    def test_generated_universe_audit_passes(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "generated"
            _write_generated_run(root)
            report = audit_run(root, "generated_universe")
            self.assertEqual(report.status, "pass")

    def test_missing_table_fails(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "generated"
            _write_generated_run(root)
            (root / "tables" / "object_table.csv").unlink()
            report = audit_run(root, "generated_universe")
            self.assertEqual(report.status, "fail")

    def test_legacy_directory_is_not_a_generated_run(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "legacy"
            (root / "tables").mkdir(parents=True)
            (root / "tables" / "global_window.csv").write_text("run_id\nlegacy\n", encoding="utf-8")
            self.assertEqual(detect_run_type(root), "unknown")

    def test_field_run_audit_passes(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "field"
            _write_field_run(root)
            report = audit_field_run(root, expected_grid_size=16)
            self.assertEqual(report.status, "pass")

    def test_pm_run_audit_passes(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp) / "pm"
            _write_pm_run(root)
            report = audit_pm_run(root, expected_particles_per_dim=8, expected_mesh_size=16)
            self.assertEqual(report.status, "pass")


def _write_generated_run(root: Path) -> None:
    (root / "tables").mkdir(parents=True)
    (root / "intermediate").mkdir(parents=True)
    (root / "quality_checks").mkdir(parents=True)
    (root / "figure_data").mkdir(parents=True)
    (root / "config_resolved.yaml").write_text("run_id: generated\n", encoding="utf-8")
    (root / "params_snapshot.json").write_text(json.dumps({"S_abs": 139.253706}), encoding="utf-8")
    (root / "generation_manifest.json").write_text(
        json.dumps(
            {
                "model": "finite_window_internal_generation",
                "redshifts": [0.0],
                "object_count": 1,
                "transfer_count": 1,
                "spectrum": {"spectrum_shape": "finite_window_envelope"},
                "inputs": ["parameter_snapshot", "run_configuration"],
                "storage_policy": "relationship_tables_only",
                "external_data": [],
                "state_decomposition": {
                    "physical_coarse_graining": "G_sigma",
                    "activity_selector": "A_Lambda",
                    "field_state": ["X_sigma", "X_A", "X_T"],
                    "external_data": [],
                },
                "closure": {
                    "map": "spectral_tail_moment_proxy",
                    "residual_norm": "unresolved_gaussian_coarse_power",
                    "epsilon0": 1e-12,
                    "norm": "log_power_integral",
                },
            }
        ),
        encoding="utf-8",
    )
    (root / "provenance.json").write_text(json.dumps({"run_id": "generated", "config": {}}), encoding="utf-8")
    _write_csv(
        root / "intermediate" / "input_validation.csv",
        "schema,path,ok,row_count,missing_columns\ninternal,internal,True,8,\n",
    )
    _write_csv(
        root / "intermediate" / "power_spectrum_samples.csv",
        "z,a,k_hmpc,p_mpc3_h3,delta2,k_over_knl,sector,tail_weight_s1\n0,1,0.2,1,0.001,1,active,1\n",
    )
    _write_csv(
        root / "intermediate" / "endpoint_solution.csv",
        "z,c_nl,k_nl_hmpc,crossing_index,method,k_left_hmpc,k_right_hmpc,delta2_left,delta2_right\n0,1,0.2,0,interp,0.1,0.3,0.5,2\n",
    )
    _write_csv(
        root / "intermediate" / "halo_selection.csv",
        "run_id,source,snapshot,z,halo_id,mass_msun_h,radius_kpc_h,sigma_v_kms,vmax_kms,sigma_source,sigma_used_kms,mass_over_m_nl,r_h_object,selected_for_r_h,quality_flag\n"
        "generated,internal_generation,a1,0,h1,1000000000000,200,100,,internal,100,0.1,0.2,True,internal_generation\n",
    )
    _write_csv(
        root / "intermediate" / "merger_growth_links.csv",
        "source,halo_id_from,halo_id_to,snapshot_from,snapshot_to,z_from,z_to,mass_from_msun_h,mass_to_msun_h,delta_mass_msun_h,time_from_gyr,time_to_gyr,delta_time_gyr,growth_rate_log_per_gyr,quality_flag\n"
        "internal_generation,h0,h1,a0,a1,0.1,0,900000000000,1000000000000,100000000000,13.6,13.8,0.2,0.5,internal_generation\n",
    )
    _write_csv(root / "intermediate" / "output_inventory.csv", "relative_path,kind,bytes,rows\nx,root,1,\n")
    _write_csv(
        root / "tables" / "global_window.csv",
        "run_id,source,snapshot,z,a,box_size_mpch,k_box_hmpc,gaussian_sigma_mpch,gaussian_kernel,activity_endpoint_hmpc,activity_selector,C_NL,k_NL_hmpc,M_NL_msun_h,R_NL_kpc_h,v_vir_MNL_kms,S_abs,S_AW,k_cut_hmpc,D_S,Omega_tail,Omega_floor,Omega_excess,R_h_mean,R_h_median,R_h_p16,R_h_p84,R_v_optional,v_lambda,tail_moment_s1,tail_injection_s0,tail_injection_s1,tail_drain_s0,tail_drain_s1,closure_active_norm,closure_model_norm,closure_residual_norm,closure_error,spectrum_source,quality_flag,notes\n"
        "generated,internal_generation,a1,0,1,14000,0.00045,1,gaussian,0.2,activity,1,0.2,10000000000000,500,300,139.253706,6,20,0,0.6,0.2,0.4,0.3,0.2,0.1,0.4,,-0.4,0.1,0.03,0.01,0.5,0.02,1,0.5,0.01,0.006,internal_finite_window_envelope,internal_generation,ok\n",
    )
    _write_csv(root / "tables" / "halo_objects.csv", _halo_header() + _halo_row())
    _write_csv(
        root / "tables" / "window_table.csv",
        "run_id,source,snapshot,z,a,box_size_mpch,k_box_hmpc,gaussian_sigma_mpch,gaussian_kernel,activity_endpoint_hmpc,activity_selector,c_nl,k_nl_hmpc,m_nl_msunh,r_nl_mpch,v_vir_mnl_kms,s_abs,s_aw,k_cut_hmpc,d_s,v_lambda,closure_active_norm,closure_model_norm,closure_residual_norm,closure_error,closure_status,spectrum_source,endpoint_status\n"
        "generated,internal_generation,a1,0,1,14000,0.00045,1,gaussian,0.2,activity,1,0.2,10000000000000,0.5,300,139.253706,6,20,0,-0.4,1,0.5,0.01,0.006,closure_audited,internal_finite_window_envelope,internal_generation\n",
    )
    _write_csv(
        root / "tables" / "band_table.csv",
        "run_id,source,snapshot,z,a,band_id,k_min_hmpc,k_max_hmpc,r_min_mpch,r_max_mpch,m_min_msunh,m_max_msunh,representative_structure,omega_band,activity_status,channel_quality\n"
        "generated,internal_generation,a1,0,1,galaxy,0.1,0.3,0.3,1.0,100000000000,10000000000000,galaxy,0.2,frontier,ok\n",
    )
    _write_csv(
        root / "tables" / "tail_table.csv",
        "run_id,source,snapshot,z,a,omega_tail_p,r_v,r_h,e_t_s0,e_t_s1,i_t_s0,i_t_s1,d_t_s0,d_t_s1,omega_floor,omega_excess,closure_error,tail_channel,q_t\n"
        "generated,internal_generation,a1,0,1,0.6,,0.2,0.6,0.1,0.03,0.01,0.5,0.02,0.2,0.4,0.006,halo_velocity_and_power_tail,internal_generation\n",
    )
    _write_csv(
        root / "tables" / "closure_table.csv",
        "run_id,source,snapshot,z,a,gaussian_sigma_mpch,activity_endpoint_hmpc,closure_map,active_norm,tail_model_norm,residual_norm,closure_error,epsilon0,norm,quality_flag\n"
        "generated,internal_generation,a1,0,1,1,0.2,spectral_tail_moment_proxy,1,0.5,0.01,0.006,1e-12,log_power_integral,closure_audited\n",
    )
    _write_csv(root / "tables" / "object_table.csv", _object_header() + _object_row())
    _write_csv(
        root / "tables" / "transfer_table.csv",
        "run_id,source,snapshot,z,source_id,target_id,source_channel,a_pq,g_pq,tau_x,tau_y,tau_z,j_x,j_y,j_z,tensor_source,spin_source,vector_status\n"
        "generated,internal_generation,a1,0,h0,h1,generated_pair_tidal_transfer,0.1,0.2,1,2,3,1,0,0,internal,internal,internal_generation\n",
    )
    _write_csv(
        root / "tables" / "cost_table.csv",
        "run_id,source,snapshot,z,a,lambda_active_hmpc,k_ref_hmpc,delta_s,rho_n_pred,rho_c_pred,active_bytes,residual_bytes,audit_bytes,total_bytes,access_seconds,audit_seconds,rho_total_measured\n"
        "generated,internal_generation,a1,0,1,0.2,200,6.9,1e-9,1e-15,100,0,20,120,0,0.01,1.2e-9\n",
    )
    _write_csv(
        root / "tables" / "event_table.csv",
        "run_id,source,snapshot,event_name,z,a,criterion,k_nl_hmpc,m_nl_msunh,s_aw,physical_reading\n"
        "generated,internal_generation,a1,active_frontier,0,1,Delta2,0.2,10000000000000,6,frontier\n",
    )
    _write_csv(
        root / "tables" / "normalization_table.csv",
        "run_id,source,snapshot,z,condition_name,condition_value,unit,role,status\n"
        "generated,internal_generation,a1,0,k_nl_hmpc,0.2,h Mpc^-1,endpoint,internal_generation\n",
    )
    _write_csv(
        root / "tables" / "object_summary_table.csv",
        "run_id,source,snapshot,z_bin,mass_bin,n_objects,b_h_median,b_h_p16,b_h_p84,forming_fraction,balanced_fraction,stripping_fraction,r_h_median,selection,quality_flag\n"
        "generated,internal_generation,a1,0,galaxy,1,0.5,0.5,0.5,0,1,0,0.2,galaxy,internal_generation\n",
    )
    _write_csv(
        root / "tables" / "transfer_summary_table.csv",
        "run_id,source,snapshot,z_bin,source_channel,mass_ratio_bin,environment_bin,n_pairs,a_pq_median,a_pq_p16,a_pq_p84,g_pq_median,g_pq_p16,g_pq_p84,quality_flag\n"
        "generated,internal_generation,a1,0,generated_pair_tidal_transfer,all,all,1,0.1,0.1,0.1,0.2,0.2,0.2,internal_generation\n",
    )
    _write_csv(
        root / "tables" / "interaction_strength_table.csv",
        "run_id,source,snapshot,z,a,band_id,channel_id,family,carrier,target,mediator,source_term,sink_term,net_term,raw_strength,chi,state_factor,scale_factor,carrier_factor,output_table,physical_reading\n"
        "generated,internal_generation,a1,0,1,galaxy,gravitational_growth,gravity,matter_field,active_window,Newtonian_potential,1,0.1,0.9,1,1,1,1,1,global_window.csv,gravity channel\n",
    )
    _write_csv(
        root / "tables" / "interaction_budget_table.csv",
        "run_id,source,snapshot,z,a,band_id,carrier,source_total,sink_total,net_total,chi_total,n_channels,quality_flag\n"
        "generated,internal_generation,a1,0,1,galaxy,matter_field,1,0.1,0.9,1,1,internal_generation\n",
    )
    _write_csv(
        root / "tables" / "interaction_dominance_table.csv",
        "run_id,source,snapshot,z,a,band_id,dominant_channel,dominant_family,dominant_chi,second_channel,second_family,second_chi,dominance_gap,chi_sum,quality_flag\n"
        "generated,internal_generation,a1,0,1,galaxy,gravitational_growth,gravity,1,gravitational_growth,gravity,1,0,1,internal_generation\n",
    )
    _write_csv(
        root / "tables" / "interaction_crossing_table.csv",
        "run_id,source,band_id,left_snapshot,right_snapshot,z_left,z_right,a_left,a_right,from_channel,to_channel,from_family,to_family,crossing_kind,quality_flag\n"
        "generated,internal_generation,galaxy,a0,a1,1,0,0.5,1,metric_seed,gravitational_growth,primordial,gravity,dominance_change,internal_generation\n",
    )
    _write_csv(
        root / "tables" / "carrier_history_table.csv",
        "run_id,source,snapshot,z,a,carrier,presence,source_total,sink_total,net_total,dominant_channel,quality_flag\n"
        "generated,internal_generation,a1,0,1,matter_field,1,1,0.1,0.9,gravitational_growth,internal_generation\n",
    )
    _write_csv(
        root / "tables" / "emergence_event_table.csv",
        "run_id,source,snapshot,z,a,event_name,band_id,criterion,value,dominant_channel,physical_reading,quality_flag\n"
        "generated,internal_generation,a1,0,1,dominance_change,galaxy,computed dominance change,1,gravitational_growth,interaction event,internal_generation\n",
    )
    _write_csv(
        root / "tables" / "schema_catalog.csv",
        "table_name,column_name,dtype,unit,required,allowed_range,description\nwindow_table.csv,z,float64,redshift,True,z > -1,epoch\n",
    )


def _halo_header() -> str:
    return (
        "run_id,source,snapshot,z,halo_id,mass_msun_h,radius_kpc_h,v_vir_kms,"
        "sigma_v_kms,sigma_source,x,y,z_pos,vx,vy,vz,spin_x,spin_y,spin_z,"
        "t_dyn_gyr,t_growth_gyr,B_h,B_class,parent_id,descendant_id,quality_flag\n"
    )


def _halo_row() -> str:
    return (
        "generated,internal_generation,a1,0,h1,1000000000000,200,140,100,internal,"
        "1,2,3,10,20,30,1,0,0,1.5,3,0.5,balanced,h0,h1,internal_generation\n"
    )


def _object_header() -> str:
    return (
        "run_id,source,snapshot,z,halo_id,x,y,zpos,vx,vy,vz,mass_msunh,"
        "r_vir_kpch,sigma_v_kms,spin_x,spin_y,spin_z,b_h,fate_class,"
        "parent_id,descendant_id,branch_weight,link_quality,q_h\n"
    )


def _object_row() -> str:
    return (
        "generated,internal_generation,a1,0,h1,1,2,3,10,20,30,1000000000000,"
        "200,100,1,0,0,0.5,balanced,h0,h1,1,internal_generation,internal_generation\n"
    )


def _write_csv(path: Path, text: str) -> None:
    path.write_text(text, encoding="utf-8")


def _write_field_run(root: Path) -> None:
    (root / "tables").mkdir(parents=True)
    (root / "figure_data").mkdir(parents=True)
    (root / "quality_checks").mkdir(parents=True)
    field_config = Field3DConfig(
        grid_size=16,
        mode_bins=6,
        max_peaks_per_snapshot=4,
        max_peak_candidates=32,
        peak_slab_cells=4,
        neighbor_links=1,
    )
    outputs = generate_field_history(
        run_id="field",
        redshifts=[1.0, 0.0],
        params=DEFAULT_PARAMS,
        spectrum_config=InternalSpectrumConfig(k_points=96),
        field_config=field_config,
    )
    write_rows(root / "tables" / "field_window_table.csv", outputs["window"], FIELD_WINDOW_COLUMNS)
    write_rows(root / "tables" / "field_mode_table.csv", outputs["mode"], FIELD_MODE_COLUMNS)
    write_rows(root / "tables" / "field_object_table.csv", outputs["object"], FIELD_OBJECT_COLUMNS)
    write_rows(root / "tables" / "field_transfer_table.csv", outputs["transfer"], FIELD_TRANSFER_COLUMNS)
    write_rows(
        root / "tables" / "field_relation_summary_table.csv",
        outputs["relation_summary"],
        FIELD_RELATION_SUMMARY_COLUMNS,
    )
    write_rows(
        root / "tables" / "field_epoch_relation_table.csv",
        outputs["epoch_relation"],
        FIELD_EPOCH_RELATION_COLUMNS,
    )
    write_rows(root / "tables" / "field_closure_table.csv", outputs["closure"], FIELD_CLOSURE_COLUMNS)
    write_rows(
        root / "tables" / "field_cost_table.csv",
        [
            field_cost_row(
                run_id="field",
                field_config=field_config,
                redshift_count=2,
                relationship_rows=sum(len(rows) for rows in outputs.values()),
                relationship_bytes_written=1,
            )
        ],
        FIELD_COST_COLUMNS,
    )
    write_field_figure_data(root, root / "figure_data")
    (root / "config_resolved.yaml").write_text("run_id: field\n", encoding="utf-8")
    (root / "params_snapshot.json").write_text(json.dumps({"S_abs": 139.253706}), encoding="utf-8")
    (root / "generation_manifest.json").write_text(
        json.dumps(
            {
                "model": "finite_window_3d_field_generation",
                "field_config": {"grid_size": 16},
                "storage_policy": "relationship_tables_only",
                "external_data": [],
                "state_decomposition": {
                    "physical_coarse_graining": "G_sigma",
                    "activity_selector": "A_Lambda",
                    "field_state": ["X_sigma", "X_A", "X_T"],
                    "external_data": [],
                },
                "closure": {
                    "map": "spectral_tail_moment_proxy",
                    "residual_norm": "unresolved_gaussian_coarse_power",
                    "epsilon0": field_config.closure_floor,
                    "norm": "log_power_integral",
                },
                "redshifts": [1.0, 0.0],
            }
        ),
        encoding="utf-8",
    )
    (root / "provenance.json").write_text(json.dumps({"run_id": "field", "config": {}}), encoding="utf-8")


def _write_pm_run(root: Path) -> None:
    (root / "tables").mkdir(parents=True)
    (root / "intermediate").mkdir(parents=True)
    (root / "quality_checks").mkdir(parents=True)
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
    redshifts = [3.0, 1.0, 0.0]
    outputs = generate_pm_history(
        run_id="pm",
        redshifts=redshifts,
        params=DEFAULT_PARAMS,
        config=config,
    )
    write_rows(root / "tables" / "pm_snapshot_table.csv", outputs["snapshot"], PM_SNAPSHOT_COLUMNS)
    write_rows(root / "tables" / "pm_power_table.csv", outputs["power"], PM_POWER_COLUMNS)
    write_rows(root / "tables" / "pm_halo_table.csv", outputs["halo"], PM_HALO_COLUMNS)
    write_rows(root / "tables" / "pm_tree_table.csv", outputs["tree"], PM_TREE_COLUMNS)
    write_rows(root / "tables" / "pm_transfer_table.csv", outputs["transfer"], PM_TRANSFER_COLUMNS)
    write_rows(
        root / "tables" / "pm_cost_table.csv",
        [
            pm_cost_row(
                run_id="pm",
                config=config,
                redshift_count=len(redshifts),
                relationship_rows=sum(len(rows) for rows in outputs.values()),
                relationship_bytes_written=1,
            )
        ],
        PM_COST_COLUMNS,
    )
    (root / "config_resolved.yaml").write_text("run_id: pm\n", encoding="utf-8")
    (root / "params_snapshot.json").write_text(json.dumps({"S_abs": 139.253706}), encoding="utf-8")
    (root / "generation_manifest.json").write_text(
        json.dumps(
            {
                "model": "finite_window_internal_pm_nbody",
                "pm_config": {"particles_per_dim": 8, "mesh_size": 16},
                "storage_policy": "relationship_tables_only",
                "external_data": [],
                "state_decomposition": {
                    "dynamical_state": [
                        "particle_phase_space",
                        "mesh_density",
                        "resolved_force",
                        "halo_catalog",
                        "transfer_graph",
                    ],
                    "physical_coarse_graining": "particle_mesh_assignment",
                    "activity_selector": "resolved_mesh_modes",
                    "external_data": [],
                },
                "closure": {
                    "map": "mesh_resolved_poisson_force",
                    "residual_norm": "unresolved_modes_above_mesh_nyquist",
                    "audit": "mass_conservation_power_bins_particle_id_overlap",
                },
                "dynamics": {
                    "solver": "particle_mesh",
                    "force": "spectral_poisson",
                    "integrator": "kick_drift_kick",
                    "tracker": "particle_id_overlap",
                    "time_variable": "scale_factor",
                    "force_prefactor": "Q(a)=force_strength/max(a,1e-4)",
                    "force_strength": 0.35,
                },
                "redshifts": redshifts,
            }
        ),
        encoding="utf-8",
    )
    (root / "provenance.json").write_text(json.dumps({"run_id": "pm", "config": {}}), encoding="utf-8")
    _write_csv(root / "intermediate" / "output_inventory.csv", "relative_path,kind,bytes,rows\nx,root,1,\n")


if __name__ == "__main__":
    unittest.main()
