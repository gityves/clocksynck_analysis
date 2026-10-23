# Subsample arguments that reproduce the archive's subsample outputs (M2, 2026-10-23)

Recovered by the M2 checks (the archive READMEs rarely state them). With these arguments and default seeds the
migrated scripts write files identical to `$SCCC_ARCHIVE/outputs/paper/<old_folder>/subsample/` (exceptions in
`paper/docs/provisory_plan_paper.md`, "M2 result per group"). `N` = `--subsample N`.

| group | script: N (other options) |
|---|---|
| method | gtex_benchmark 24 (`--ncores 2`); celltype_truth 16, `--window` 24; baboon_truth 24; skin_truth 24; tcga_cofe 24 (`--ncores 2`); cyclops_anchoring: not run (> 30 min at 30) |
| windows | internal_time 120 (cohorts: 150); baboon_windows 12; chen_windows 60; breast_timed 20; nontemporal 120; window_check 120 |
| species | colon_injury 120 (also `--counts raw`); sheep 120; baboon 12; baboon_stability 12 (`--n-perm 20`); sheep_study_gates 60 (`--n-perm 20`); mouse_skin_timed 120 |
| organs | organs_timed 60 (fast), 80 (`--cohorts ventilator`); brain_artefact 60; brain_artefact_d4 60; brain_periphery_species 200; species_rstar_full 200; chen_age 60 |
| backbone | oscillator_gtex 30; family_summary 30; per_cell 30; cancer_split <= 12 (floored at 12); positive_control 30; followups_a18 40 |
| embryo | all scripts 60 |
| atlas_fits | atlas_fits 30; step7_confirmations 30 |
| gene_rhythms | gtex_truth 60; atlas_gene_rhythms 60 (`--n-null 5`); atlas_gene_rhythms_modes 60 (`--atlases breast --n-null 5`); donor_curve 60; gene_map 100; celltype_truth 60; blood_clock 150 |
| cues_lost | all scripts 60 |
| gtex_cues | cue_reference, cue_hierarchy, internal_external, drug_targets, chiral_age, covariate_sensitivity, ventilation_dose 150; metabolic, coupling_sources 200; glucocorticoid, dip_arm 250; gonogo_diagnostics 150 (diag1) / 200 (diag2, diag4); gc_atlas 30; chiral_prep_counts 500 |
| resilience | contrasts_a17 40; procurement 60; age 60; masld_gut 50; mouse_ageing 20; gtex_age 500; gtex_age_continuous 600 (`--draws 3`); gtex_chiral_age 600; gtex_sex_bmi 250; a45_learned >= 70 |
| spread_multiplicity | spread_inputs 40; phase_cancer 30; hierarchy (`--main-only` in the archive) |
| natural_experiments | covid_gates 20; covid_couplings 24; covid_ventilation 25; fetal_secondary 24; placenta_sensitivity 20; gtex_reference 4 |
| calibration | small_n 30; anchored_simulations 4; orpheus_bias 100; orpheus_cells 30; vonmises_calibration 20 |
| gate | gate_power 60; extended_clock 150 (fails at 60: empty donor set); baboon_reference 12; waveform_sensitivity 40 (`--only lung_hlca heart islets`); extended_e4 60; gene_screen 60; sheep_reference_primary: N not recovered (30 gives 13 tissues, archive 16) |
| host_identity | all scripts 30, cortex_modules 40 |
| simulations | hierarchy_validation 3; a22_colon 10; celltype_truth 3 (skin; colon > 25 min); r26_cohorts 40 |
| cancer | anchors_cancer 30; cue_following 150 (12 crashes: no classes); others: README N |
| lung_breaks | modules 40; a18_followups 80; decomposition 80 |
| psoriasis | README N 12; sc_status 20 (`--n-perm 8 --n-boot 8`). Run cues, and psoriasis_details cues_matched and cue_context, with `--jobs 1` (> 30 GB with 2 jobs) |
| psoriasis_details | README N 12; gtex_skin 60 (`--n-perm 100`); taurus_amplitude 30 (`--n-perm 6 --n-pseudo 4`); bin_mix 24 (`--n-boot 4`); differentiation 24 (`--n-boot 6 --n-perm 6`); sc_attribution 24 (`--n-boot 6`) |
| cancer_pipeline | 12 |
