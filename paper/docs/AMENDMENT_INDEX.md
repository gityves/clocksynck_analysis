# Amendment index: archive amendments -> migrated scripts

Amendments A1-A108c, V and D notes were registered in the archive (`$SCCC_ARCHIVE/docs/PLAN_internal_desync.md`,
results in `$SCCC_ARCHIVE/docs/RESULTS_preliminary.md`; `$SCCC_ARCHIVE` =
/home/maxine/Documents/paychere/scCircadianCancer). This table maps each one to the script that computes it here (paths
relative to `paper/`), derived from the archive's generated index (2026-10-22) and the migration map of
`common/migration/migrate_from_archive.py`. "archive only": not migrated (not cited by the paper plan); use the archive
script if it is ever needed. Amendments missing here have no computation or were not ported in the archive either: the
legacy code (`$SCCC_ARCHIVE/legacy/`, catalogue in `legacy/catalogue/`) is then the reference. New amendments (A109+)
are registered in `docs/registered_experiments.md` and added here when their scripts are validated.

| amendment | scripts |
|---|---|
| A1 | `results/extended_data/scripts/calibration/compute_orpheus_cells.py`, `results/extended_data/scripts/simulations/compute_celltype_truth.py`, `results/figure_2/scripts/atlas_fits/compute_atlas_fits.py`, `results/figure_2/scripts/organs/compute_organs_timed.py`, `results/figure_4/scripts/resilience/compute_contrasts_a17.py` |
| A3 | `results/figure_2/scripts/organs/compute_chen_age.py` |
| A4 | `results/figure_1/scripts/windows/compute_window_check.py`; archive only: `figure_6_lung_breaks/compute_a6_tstates.py` |
| A5 | archive only: `figure_6_lung_breaks/compute_a5_ild.py`, `figure_6_lung_breaks/compute_a6_tstates.py`, `figure_S04_host_identity/compute_axis_reading.py`, `figure_S04_host_identity/compute_common_axis.py` |
| A6 | archive only: `figure_6_lung_breaks/compute_a6_stress.py`, `figure_6_lung_breaks/compute_a6_tstates.py` |
| A7 | archive only: `figure_6_lung_breaks/compute_a7_activation.py` |
| A8 | `results/extended_data/scripts/calibration/compute_small_n.py` |
| A9 | `results/extended_data/scripts/host_identity/compute_cortex_modules.py`; archive only: `figure_6_lung_breaks/compute_a13.py`, `figure_S04_host_identity/compute_axis_reading.py`, `figure_S04_host_identity/compute_common_axis.py` |
| A10 | `results/extended_data/scripts/calibration/compute_small_n.py`, `results/extended_data/scripts/gate/compute_gate_power.py` |
| A12 | archive only: `figure_6_lung_breaks/compute_a13.py`, `figure_S04_host_identity/compute_axis_reading.py`, `figure_S04_host_identity/compute_common_axis.py` |
| A13 | archive only: `figure_6_lung_breaks/compute_a13.py` |
| A16 | `results/extended_data/scripts/simulations/compute_celltype_truth.py`, `results/figure_2/scripts/atlas_fits/compute_atlas_fits.py`, `results/figure_4/scripts/resilience/compute_contrasts_a17.py` |
| A17 | `results/extended_data/scripts/calibration/compute_anchored_simulations.py`, `results/extended_data/scripts/calibration/compute_orpheus_cells.py`, `results/extended_data/scripts/gate/compute_gate_power.py`, `results/extended_data/scripts/gate/compute_waveform_sensitivity.py`, `results/extended_data/scripts/host_identity/compute_anchors.py`, `results/extended_data/scripts/host_identity/compute_coherence.py`, `results/extended_data/scripts/host_identity/compute_cortex_modules.py`, `results/extended_data/scripts/host_identity/compute_host_identity.py`, `results/extended_data/scripts/simulations/compute_celltype_truth.py`, `results/figure_2/scripts/atlas_fits/compute_atlas_fits.py`, `results/figure_2/scripts/atlas_fits/compute_step7_confirmations.py`, `results/figure_2/scripts/backbone/compute_family_summary.py`, `results/figure_2/scripts/backbone/compute_followups_a18.py`, `results/figure_2/scripts/backbone/compute_per_cell.py`, `results/figure_2/scripts/gene_rhythms/compute_atlas_gene_rhythms.py`, `results/figure_2/scripts/gene_rhythms/compute_donor_curve.py`, `results/figure_2/scripts/gene_rhythms/compute_gene_map.py`, `results/figure_3/scripts/gtex_cues/compute_gc_atlas.py`, `results/figure_4/scripts/resilience/compute_a45_learned.py`, `results/figure_4/scripts/resilience/compute_age.py`, `results/figure_4/scripts/resilience/compute_contrasts_a17.py`, `results/figure_4/scripts/resilience/compute_procurement.py`, `results/figure_4/scripts/spread_multiplicity/compute_multiplicity.py`, `results/figure_4/scripts/spread_multiplicity/compute_phase_cancer.py`, `results/figure_5/scripts/cancer/compute_anchored_cancer.py`, `results/figure_5/scripts/cancer/compute_anchors_cancer.py`, `results/figure_5/scripts/cancer/compute_cycling_a24.py`, `results/figure_5/scripts/lung_breaks/compute_a18_followups.py`, `results/figure_5/scripts/lung_breaks/compute_decomposition.py`, `results/figure_5/scripts/lung_breaks/compute_modules.py`; archive only: `figure_6_lung_breaks/compute_a13.py`, `figure_6_lung_breaks/compute_a5_ild.py`, `figure_6_lung_breaks/compute_a6_stress.py`, `figure_6_lung_breaks/compute_a6_tstates.py`, `figure_S04_host_identity/compute_axis_reading.py`, `figure_S04_host_identity/compute_common_axis.py`, `figure_S08_gene_rhythms/compute_bulk_clock.py` |
| A18 | `results/extended_data/scripts/calibration/compute_orpheus_cells.py`, `results/extended_data/scripts/host_identity/compute_cortex_modules.py`, `results/figure_1/scripts/windows/compute_window_check.py`, `results/figure_2/scripts/backbone/compute_cancer_split.py`, `results/figure_2/scripts/backbone/compute_followups_a18.py`, `results/figure_2/scripts/backbone/compute_per_cell.py`, `results/figure_2/scripts/backbone/compute_positive_control.py`, `results/figure_4/scripts/resilience/compute_contrasts_a17.py`, `results/figure_4/scripts/resilience/compute_procurement.py`, `results/figure_4/scripts/spread_multiplicity/compute_spread_inputs.py`, `results/figure_5/scripts/cancer/compute_anchored_cancer.py`, `results/figure_5/scripts/lung_breaks/compute_a18_followups.py`, `results/figure_5/scripts/lung_breaks/compute_modules.py` |
| A19 | `results/extended_data/scripts/calibration/compute_age_acceptance.py`, `results/extended_data/scripts/calibration/compute_orpheus_bias.py`, `results/extended_data/scripts/calibration/compute_orpheus_cells.py`, `results/extended_data/scripts/simulations/compute_celltype_truth.py`, `results/figure_4/scripts/resilience/compute_age.py`, `results/figure_4/scripts/resilience/compute_procurement.py` |
| A20 | `results/extended_data/scripts/host_identity/compute_group_anchors.py`, `results/extended_data/scripts/simulations/compute_hierarchy_validation.py`, `results/figure_2/scripts/atlas_fits/compute_atlas_fits.py`, `results/figure_4/scripts/spread_multiplicity/compute_hierarchy.py`, `results/figure_4/scripts/spread_multiplicity/compute_multiplicity.py`, `results/figure_4/scripts/spread_multiplicity/compute_spread_inputs.py` |
| A21 | `results/extended_data/scripts/host_identity/compute_group_anchors.py`, `results/extended_data/scripts/simulations/compute_a22_colon.py`, `results/figure_1/scripts/species/compute_colon_injury.py` |
| A22 | `results/extended_data/scripts/host_identity/compute_group_anchors.py`, `results/extended_data/scripts/simulations/compute_a22_colon.py`, `results/figure_4/scripts/resilience/compute_a45_learned.py`, `results/figure_4/scripts/resilience/compute_masld_gut.py`, `results/figure_4/scripts/spread_multiplicity/compute_spread_inputs.py` |
| A23 | `results/extended_data/scripts/simulations/compute_hierarchy_validation.py`, `results/figure_4/scripts/spread_multiplicity/compute_hierarchy.py`, `results/figure_4/scripts/spread_multiplicity/compute_phase_cancer.py`, `results/figure_4/scripts/spread_multiplicity/compute_spread_inputs.py` |
| A24 | `results/figure_5/scripts/cancer/compute_cycling_a24.py` |
| A25 | `results/extended_data/scripts/gate/compute_extended_clock.py`, `results/extended_data/scripts/gate/compute_extended_e4.py`, `results/extended_data/scripts/simulations/compute_r26_cohorts.py`, `results/figure_1/scripts/species/compute_baboon.py`, `results/figure_1/scripts/windows/compute_baboon_windows.py`, `results/figure_1/scripts/windows/compute_breast_timed.py`, `results/figure_1/scripts/windows/compute_chen_windows.py`, `results/figure_1/scripts/windows/compute_internal_time.py`, `results/figure_1/scripts/windows/compute_nontemporal.py`, `results/figure_3/scripts/gtex_cues/compute_gc_li_breast.py` |
| A26 | `results/extended_data/scripts/simulations/compute_r26_cohorts.py`, `results/figure_1/scripts/windows/compute_nontemporal.py` |
| A27 | `results/figure_2/scripts/gene_rhythms/compute_atlas_gene_rhythms.py`, `results/figure_2/scripts/gene_rhythms/compute_atlas_gene_rhythms_modes.py`, `results/figure_2/scripts/gene_rhythms/compute_donor_curve.py`, `results/figure_2/scripts/gene_rhythms/compute_gtex_truth.py`, `results/figure_3/scripts/cues_lost/compute_walltime_hosts.py`, `results/figure_3/scripts/gtex_cues/compute_dip_arm.py`, `results/figure_3/scripts/gtex_cues/compute_gc_atlas.py`, `results/figure_3/scripts/gtex_cues/compute_gc_li_breast.py`, `results/figure_3/scripts/gtex_cues/compute_glucocorticoid.py`, `results/figure_3/scripts/gtex_cues/compute_gonogo_diagnostics.py`, `results/figure_3/scripts/gtex_cues/compute_internal_external.py`; archive only: `figure_S08_gene_rhythms/compute_antiphase_modes.py` |
| A28 | `results/extended_data/scripts/gate/compute_extended_clock.py`, `results/extended_data/scripts/gate/compute_extended_e4.py`, `results/extended_data/scripts/gate/compute_gene_screen.py`, `results/extended_data/scripts/gate/compute_sheep_reference_primary.py` |
| A29 | archive only: `figure_S08_gene_rhythms/compute_bulk_clock.py` |
| A30 | archive only: `figure_S08_gene_rhythms/compute_antiphase.py`, `figure_S08_gene_rhythms/compute_antiphase_modes.py` |
| A31 | `results/figure_1/scripts/method/compute_celltype_truth.py` (the `gene_rhythms` copy was removed 2026-10-25: identical estimates) |
| A31b | `results/figure_2/scripts/gene_rhythms/compute_atlas_gene_rhythms.py`, `results/figure_2/scripts/gene_rhythms/compute_atlas_gene_rhythms_modes.py`, `results/figure_2/scripts/gene_rhythms/compute_donor_curve.py`, `results/figure_2/scripts/gene_rhythms/compute_gene_map.py`, `results/figure_2/scripts/gene_rhythms/compute_gtex_truth.py`; archive only: `figure_S08_gene_rhythms/compute_antiphase.py`, `figure_S08_gene_rhythms/compute_antiphase_modes.py` |
| A32 | `results/figure_2/scripts/gene_rhythms/compute_blood_clock.py`; archive only: `figure_S08_gene_rhythms/compute_blood_gtex_arms.py` |
| A33 | `results/figure_3/scripts/cues_lost/compute_ventilator_gc.py`, `results/figure_3/scripts/cues_lost/compute_ventilator_tissues.py` |
| A34 | `results/extended_data/scripts/gate/compute_baboon_reference.py`, `results/extended_data/scripts/gate/compute_sheep_reference_primary.py`, `results/figure_1/scripts/species/compute_sheep.py`, `results/figure_1/scripts/species/compute_sheep_study_gates.py` |
| A35 | `results/figure_2/scripts/gene_rhythms/compute_blood_clock.py` |
| A36 | `results/extended_data/scripts/gate/compute_baboon_reference.py`, `results/figure_1/scripts/method/compute_baboon_truth.py`, `results/figure_1/scripts/species/compute_baboon.py`, `results/figure_1/scripts/species/compute_baboon_stability.py`, `results/figure_1/scripts/windows/compute_baboon_windows.py`, `results/figure_4/scripts/resilience/compute_mouse_ageing.py` |
| A37 | `results/figure_1/scripts/method/compute_skin_truth.py` |
| A39 | `results/figure_2/scripts/gene_rhythms/compute_donor_curve.py`, `results/figure_2/scripts/gene_rhythms/compute_gene_map.py` |
| A40 | `results/extended_data/scripts/gate/compute_baboon_reference.py`, `results/figure_1/scripts/species/compute_baboon.py`, `results/figure_1/scripts/windows/compute_baboon_windows.py` |
| A41 | `results/extended_data/scripts/gate/compute_baboon_reference.py`, `results/extended_data/scripts/gate/compute_sheep_reference_primary.py`, `results/figure_1/scripts/species/compute_baboon.py` |
| A42 | `results/figure_1/scripts/species/compute_baboon.py`, `results/figure_2/scripts/organs/compute_brain_periphery_species.py`, `results/figure_2/scripts/organs/compute_species_rstar_full.py` |
| A43 | `results/figure_4/scripts/resilience/compute_masld_gut.py` |
| A44 | `results/figure_4/scripts/resilience/compute_masld_gut.py` |
| A45 | `results/figure_4/scripts/resilience/compute_a45_learned.py`, `results/figure_4/scripts/resilience/compute_masld_gut.py`, `results/figure_5/scripts/psoriasis_details/compute_taurus_amplitude.py` |
| A46 | `results/figure_3/scripts/cues_lost/compute_walltime_hosts.py`, `results/figure_3/scripts/gtex_cues/compute_ventilation_dose.py` |
| A47 | `results/figure_3/scripts/gtex_cues/compute_metabolic.py` |
| A48 | `results/figure_4/scripts/resilience/compute_mouse_ageing.py` |
| A50 | `results/figure_2/scripts/gene_rhythms/compute_gene_map.py`, `results/figure_2/scripts/gene_rhythms/compute_gtex_truth.py` |
| A52 | `results/figure_5/scripts/cancer/compute_matched_cues_tcga.py`, `results/figure_5/scripts/cancer/compute_tcga_paired.py`, `results/figure_5/scripts/psoriasis_details/compute_cue_context.py`, `results/figure_5/scripts/psoriasis_details/compute_cues_diseases.py`, `results/figure_5/scripts/psoriasis_details/compute_cues_matched.py` |
| A53 | `results/figure_5/scripts/psoriasis_details/compute_cue_context.py`, `results/figure_5/scripts/psoriasis_details/compute_cues_companion.py`, `results/figure_5/scripts/psoriasis_details/compute_cues_diseases.py`, `results/figure_5/scripts/psoriasis_details/compute_other_diseases.py`, `results/figure_5/scripts/psoriasis_details/compute_taurus_amplitude.py` |
| A54 | `results/figure_5/scripts/psoriasis/compute_cues.py`, `results/figure_5/scripts/psoriasis/compute_cues_gtex.py`, `results/figure_5/scripts/psoriasis/compute_dose.py`, `results/figure_5/scripts/psoriasis/compute_replication.py`, `results/figure_5/scripts/psoriasis/compute_treatment.py`, `results/figure_5/scripts/psoriasis_details/compute_bin_mix.py`, `results/figure_5/scripts/psoriasis_details/compute_cue_context.py`, `results/figure_5/scripts/psoriasis_details/compute_cues_diseases.py`, `results/figure_5/scripts/psoriasis_details/compute_layers.py`, `results/figure_5/scripts/psoriasis_details/compute_other_diseases.py`, `results/figure_5/scripts/psoriasis_details/compute_sc_attribution.py` |
| A55 | `results/figure_5/scripts/psoriasis_details/compute_cue_context.py`, `results/figure_5/scripts/psoriasis_details/compute_cues_diseases.py`, `results/figure_5/scripts/psoriasis_details/compute_other_diseases.py` |
| A56 | `results/figure_5/scripts/psoriasis/compute_dose.py`, `results/figure_5/scripts/psoriasis/compute_treatment.py` |
| A57 | `results/figure_5/scripts/psoriasis/compute_response.py` |
| A58 | `results/figure_5/scripts/psoriasis/compute_sc_status.py`, `results/figure_5/scripts/psoriasis_details/compute_gtex_skin.py`, `results/figure_5/scripts/psoriasis_details/compute_sc_attribution.py` |
| A59 | `results/figure_5/scripts/psoriasis_details/compute_cues_companion.py`, `results/figure_5/scripts/psoriasis_details/compute_cues_diseases.py`, `results/figure_5/scripts/psoriasis_details/compute_other_diseases.py` |
| A60 | `results/figure_5/scripts/psoriasis_details/compute_bin_mix.py`, `results/figure_5/scripts/psoriasis_details/compute_sc_attribution.py` |
| A61 | `results/figure_5/scripts/psoriasis_details/compute_gtex_skin.py`, `results/figure_5/scripts/psoriasis_details/compute_sc_attribution.py` |
| A62 | `results/figure_5/scripts/psoriasis_details/compute_layers.py` |
| A63 | `results/figure_5/scripts/psoriasis_details/compute_bin_mix.py`, `results/figure_5/scripts/psoriasis_details/compute_differentiation.py` |
| A64 | `results/figure_5/scripts/psoriasis_details/compute_differentiation.py` |
| A65 | `results/figure_5/scripts/psoriasis_details/compute_bin_mix.py`, `results/figure_5/scripts/psoriasis_details/compute_differentiation.py` |
| A66 | `results/figure_5/scripts/psoriasis_details/compute_differentiation.py` |
| A67 | `results/figure_5/scripts/psoriasis/compute_cues.py`, `results/figure_5/scripts/psoriasis_details/compute_timed_layers.py` |
| A68 | `results/figure_5/scripts/psoriasis_details/compute_gtex_skin.py` |
| A69 | `results/figure_5/scripts/psoriasis_details/compute_bin_mix.py` |
| A71 | `results/figure_1/scripts/species/compute_mouse_skin_timed.py` |
| A72 | `results/figure_1/scripts/species/compute_mouse_skin_timed.py` |
| A73 | `results/figure_5/scripts/psoriasis/compute_dose.py`, `results/figure_5/scripts/psoriasis/compute_response.py` |
| A74 | `results/figure_5/scripts/psoriasis/compute_response.py` |
| A76 | `results/figure_5/scripts/psoriasis/compute_cues_gtex.py`, `results/figure_5/scripts/psoriasis_details/compute_timed_layers.py` |
| A77 | `results/figure_5/scripts/psoriasis_details/compute_layers.py` |
| A78 | `results/figure_5/scripts/psoriasis/compute_cues.py`, `results/figure_5/scripts/psoriasis/compute_cues_gtex.py`, `results/figure_5/scripts/psoriasis_details/compute_cues_companion.py`, `results/figure_5/scripts/psoriasis_details/compute_cues_diseases.py` |
| A79 | `results/figure_5/scripts/psoriasis_details/compute_cues_companion.py` |
| A81 | `results/figure_3/scripts/gtex_cues/compute_cue_reference.py`, `results/figure_5/scripts/psoriasis/compute_cues_gtex.py`, `results/figure_5/scripts/psoriasis_details/compute_cue_context.py` |
| A81b | `results/figure_3/scripts/gtex_cues/compute_cue_reference.py`, `results/figure_3/scripts/gtex_cues/compute_internal_external.py`, `results/figure_4/scripts/resilience/compute_gtex_sex_bmi.py`, `results/figure_5/scripts/cancer/compute_cue_following.py`, `results/figure_5/scripts/cancer/compute_matched_cues_tcga.py`, `results/figure_5/scripts/psoriasis/compute_cues_gtex.py`, `results/figure_5/scripts/psoriasis_details/compute_cue_context.py`, `results/figure_5/scripts/psoriasis_details/compute_cues_diseases.py`, `results/figure_5/scripts/psoriasis_details/compute_cues_matched.py` |
| A82 | `results/figure_5/scripts/psoriasis_details/compute_cues_diseases.py` |
| A83 | `results/figure_5/scripts/cancer/compute_malignant_gate.py`, `results/figure_5/scripts/cancer/compute_tumour_gate.py` |
| A84 | `results/figure_3/scripts/gtex_cues/compute_internal_external.py`, `results/figure_5/scripts/cancer/compute_matched_cues_tcga.py`, `results/figure_5/scripts/psoriasis_details/compute_cues_matched.py` |
| A85 | `results/figure_5/scripts/cancer/compute_cue_following.py` |
| A86 | `results/figure_5/scripts/cancer/compute_malignant_gate.py` |
| A87 | `results/figure_5/scripts/cancer/compute_decomposition.py`, `results/figure_5/scripts/cancer/compute_proliferation.py`, `results/figure_5/scripts/lung_breaks/compute_decomposition.py` |
| A88 | `results/figure_5/scripts/psoriasis_details/compute_cue_context.py` |
| A89 | `results/figure_3/scripts/gtex_cues/compute_cue_hierarchy.py`, `results/figure_3/scripts/gtex_cues/compute_drug_targets.py`, `results/figure_3/scripts/gtex_cues/compute_internal_external.py` |
| A90 | `results/figure_3/scripts/gtex_cues/compute_cue_hierarchy.py` |
| A91 | `results/figure_3/scripts/gtex_cues/compute_chiral_age.py`, `results/figure_3/scripts/gtex_cues/compute_covariate_sensitivity.py`, `results/figure_3/scripts/gtex_cues/compute_drug_targets.py`, `results/figure_3/scripts/gtex_cues/compute_internal_external.py`, `results/figure_4/scripts/resilience/compute_gtex_chiral_age.py` |
| A92 | `results/figure_3/scripts/cues_lost/compute_walltime_hosts.py`, `results/figure_3/scripts/gtex_cues/compute_coupling_sources.py` |
| A93 | `results/figure_3/scripts/cues_lost/compute_walltime_hosts.py`, `results/figure_3/scripts/gtex_cues/compute_coupling_sources.py` |
| A94 | `results/figure_4/scripts/resilience/compute_gtex_age.py`, `results/figure_4/scripts/resilience/compute_gtex_age_continuous.py`, `results/figure_4/scripts/resilience/compute_gtex_age_tiebreak.py`, `results/figure_4/scripts/resilience/compute_gtex_sex_bmi.py` |
| A94b | `results/figure_4/scripts/resilience/compute_gtex_age_continuous.py` |
| A96 | `results/figure_4/scripts/resilience/compute_gtex_sex_bmi.py` |
| A97 | `results/figure_5/scripts/cancer/compute_delay.py`, `results/figure_5/scripts/cancer/compute_proliferation.py` |
| A98 | `results/figure_5/scripts/lung_breaks/compute_decomposition.py` |
| A99 | `results/figure_4/scripts/resilience/compute_gtex_age.py`, `results/figure_4/scripts/resilience/compute_gtex_sex_bmi.py` |
| A100 | `results/figure_3/scripts/gtex_cues/compute_covariate_sensitivity.py`, `results/figure_3/scripts/gtex_cues/compute_drug_targets.py`, `results/figure_3/scripts/gtex_cues/compute_internal_external.py`, `results/figure_4/scripts/resilience/compute_gtex_chiral_age.py` |
| A101 | `results/figure_3/scripts/gtex_cues/compute_chiral_age.py`, `results/figure_3/scripts/gtex_cues/compute_covariate_sensitivity.py`, `results/figure_4/scripts/resilience/compute_gtex_chiral_age.py` |
| A101b | `results/figure_3/scripts/gtex_cues/compute_chiral_age.py`, `results/figure_3/scripts/gtex_cues/compute_chiral_prep_counts.py`, `results/figure_3/scripts/gtex_cues/compute_covariate_sensitivity.py` |
| A102 | `results/figure_3/scripts/gtex_cues/compute_chiral_prep_counts.py`, `results/figure_3/scripts/gtex_cues/compute_covariate_sensitivity.py` |
| A102b | `results/figure_3/scripts/gtex_cues/compute_chiral_prep_counts.py`, `results/figure_3/scripts/gtex_cues/compute_covariate_sensitivity.py`, `results/figure_3/scripts/gtex_cues/compute_internal_external.py` |
| A103 | `results/figure_2/scripts/organs/compute_brain_artefact.py`, `results/figure_2/scripts/organs/compute_brain_artefact_d4.py`, `results/figure_2/scripts/organs/compute_brain_periphery_species.py`, `results/figure_2/scripts/organs/compute_species_rstar_full.py` |
| A104 | `results/figure_2/scripts/backbone/compute_oscillator_gtex.py` |
| A105 | `results/figure_4/scripts/natural_experiments/compute_covid_couplings.py`, `results/figure_4/scripts/natural_experiments/compute_covid_gates.py`, `results/figure_4/scripts/natural_experiments/compute_covid_ventilation.py`, `results/figure_4/scripts/natural_experiments/compute_gtex_reference.py` |
| A105b | `results/figure_4/scripts/natural_experiments/compute_covid_gates.py` |
| A105c | `results/figure_4/scripts/natural_experiments/compute_covid_ventilation.py` |
| A106 | `results/figure_2/scripts/embryo/compute_placenta.py`, `results/figure_4/scripts/natural_experiments/compute_placenta_sensitivity.py` |
| A107 | `results/figure_2/scripts/embryo/compute_fetal_atlases.py`, `results/figure_4/scripts/natural_experiments/compute_fetal_secondary.py` |
| A107b | `results/figure_2/scripts/embryo/compute_fetal_atlases.py` |
| A108 | `results/figure_2/scripts/embryo/compute_embryo_emergence.py`, `results/figure_2/scripts/embryo/compute_embryo_lineages.py`, `results/figure_2/scripts/embryo/compute_embryo_timed.py` |
| A108b | `results/figure_2/scripts/embryo/compute_embryo_timed.py` |
| A108c | `results/figure_2/scripts/embryo/compute_embryo_emergence.py`, `results/figure_2/scripts/embryo/compute_embryo_lineages.py` |
