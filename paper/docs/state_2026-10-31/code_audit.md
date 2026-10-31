> Appendix to `../STATE_2026-10-31.md`. A Sonnet subagent compiled it read-only on 2026-10-31; the lead spot-checked it
> and made the two corrections marked "Lead check".

# Code-to-result tracing audit, sccc (audit date 2026-10-31)

Read only. Evidence: `git log` of the outer repo (`paper/`) and of the nested `clocksync/` repo (the outer repo ignores `/clocksync/`, so `git log -- clocksync` in the outer repo returns nothing), `paper/results/logs/jobs.tsv` (204 lines), output mtimes, `paper/results/logs/{m3_todo.txt,m3_rerun_2026-10-29_compare.md,phase3_compare.md,class_fix_tables.md,rebuild_and_stamp.md}`, plan, STATUS, HANDOFF, REVIEW. Times are local, from `git log`/`jobs.tsv` (jobs.tsv column 2 = end time; times below are end times unless "start").

Caveats up front:
- `clocksync/CHANGELOG.md` has NOT been updated since the first version (`git diff 1c01018 HEAD -- CHANGELOG.md` is empty). Every change below is classified from commit message plus `git show` diff, not from the changelog.
- jobs.tsv exists only from 2026-10-24 (commit 57a6f96). Runs before that (archive-era, M1/M2 subsample checks) have no ledger line. Runs launched outside `run.sh` (a bare python) are not in it either.
- Outputs are git-ignored; mtimes can be touched by later partial reruns, so "last run" is the newest file mtime in the group's output folder cross-checked with jobs.tsv.
- No script was run. "CURRENT" means "no relevant default-changing change found by commit message and diff reading", not a re-verification.

## 1. Code changes since 2026-10-23

### 1a. clocksync (nested repo, `/home/maxine/Documents/paychere/scCircadianOverview/clocksync`)

Commit before the window: 1c01018 (2026-10-22 15:33, API additions) and f7f3c86 (2026-10-23 20:46, pickling).

| commit | date (time) | what | class | results touched |
|---|---|---|---|---|
| f7f3c86 | 10-23 20:46 | callable modules pickle by name (joblib) | bit-identical (infrastructure) | none |
| 0a288ac | 10-24 12:45 | `gene_columns(missing=ignore/warn/raise)`, all-absent guard | opt-in; "default behaviour unchanged" | none (paper loaders made to `raise` in 387a214, see 1b) |
| 09d0142 | 10-25 00:53 | factor EM rewritten as 2x2 algebra by component | numerically identical on converged fits (max diff 5e-11 / 2e-12 stated in message); NOT identical on fits that hit `max_iter` (m3_todo 10-25 ~22h: stop at different points, z diff 0.19, skin endothelial and breast_hbca endothelial) | default-changing only for non-converged fits |
| 25ec04b | 10-25 20:43 | `gene_rhythms` reports non-converged hosts (RuntimeWarning, `.attrs['hosts_not_converged']`) | opt-in diagnostic | none by itself |
| 13b57a3 | 10-25 20:48 | `gene_rhythms`: refit non-converged host with 100,000 EM maps; still failing = not estimable (rows dropped) | DEFAULT-CHANGING for non-converged hosts | gene_rhythms, any gene_rhythms caller |
| e1b6658 | 10-26 10:26 | `gene_rhythms`: compartment with n <= design columns not estimable (previously se 0, p 0) | DEFAULT-CHANGING for saturated designs | A39 donor subsamples (rerun 10-26) |
| 613e224 | 10-26 13:14 | `covariate`: within-stratum permutation null is the DEFAULT (`null='simulated'` kept); gamma not estimable at own-noise bound | DEFAULT-CHANGING | every `Dataset.covariate` caller (see 2) |
| e89defd | 10-28 17:02 | Phase 1: bounded/flagged ccorr, responsive pseudo sets, `drop_genes`, `z2/amplitude/coupling_debiased`, `boundary_type`, per-arm host gate in `contrast`, `n_boot_failed`, Hartung-Knapp meta, below-paper-settings warnings, lymphatic excluded from `suggest_anchors` | "defaults bit-identical" (message); diff confirms additions are new columns / opt-in args (`require_host_gate_each_arm=False`, `pseudo_responsive=False`, `drop_genes=None`, `nan_on_flag=False`). EXCEPTIONS found in the diff: `suggest_anchors` now excludes lymphatic endothelium by default; new warnings. The HANDOFF/results_preliminary regression check ("defaults bit-identical on gate, fit, contrast, coherence, paired, ccorr, meta, pseudo_waveforms; 47 tests") is the only equivalence evidence; I did not rerun it | opt-in / bit-identical, except `suggest_anchors` | `suggest_anchors` is called only in `extended_data/scripts/host_identity/anchor_rule.py` |
| bc6c419 | 10-28 20:50 | `anchor_readout='loo'` opt-in | opt-in; "defaults bit-identical" | none unless a script passes it (A155, Phase 3 readouts) |
| 6868d5d | 10-29 07:25 | opt-in `LEGACY_ALIASES` (`canonical(legacy=True)`, `gene_columns(legacy=True)`) | opt-in; default unchanged | none yet; Phase 8 item would make it default (see 4) |
| 9db4ed8 | 10-30 18:25 | `p_log_z2`, `p_log_psi` in `Contrast`; `host.reliability`, `Dataset.fit(attenuation=)` | opt-in / bit-identical ("existing columns bit-identical") | A163/A165 only |

Pre-existing uncommitted files inside clocksync (not ours, per HANDOFF): `.gitignore`, `pyproject.toml`, `tests/test_additions.py`, `src/clocksync/utils.py`, `tests/__init__.py`. They are still uncommitted today, so the working tree that produced results may differ from HEAD (`utils.py` is untracked; I did not diff it).

### 1b. paper/common (outer repo)

| commit | date (time) | what | class | results touched |
|---|---|---|---|---|
| cb4a452 | 10-23 | builder `tabula_sapiens_organs` | new data | A109b |
| 57a6f96, abefb52, 06ed03a, bfaa0b0, d0bc4d7, b46a6bd | 10-24/25 | run.sh, jobs.sh, watch_jobs.sh, m3_rerun.sh, cli.py (M3_JOBS/M3_AFTER) | infrastructure | none |
| 09cfe2c, b1831bb, 9e25333 | 10-24 | new builders: greenwell2019, jouffe2022_bmal1ko, manella2021, vipko, xin2021, adipose/muscle timed, breast_timed (+ age) | new data | A113, A129, A130 only |
| f784ddf | 10-24 | `anchor_robustness.py` | new analysis | WP0 item 8 |
| 387a214 | 10-24 | `gtex.py`, `natexp_common`: `gene_columns` raises on missing genes | DEFAULT-CHANGING (error instead of silent -1 index) | GTEx and natural-experiment loaders; outputs of figure_3/gtex_cues (10-26 00:16, after) and natural_experiments (10-26 12:08, after) already post-date it |
| 281bb26 | 10-28 14:39 | class rules: case-sensitive whole-word "T cell"/"B cell" (mast, goblet, tuft, myofibroblast, club were misclassed); heart lymphocyte no longer takes mast | DEFAULT-CHANGING (pseudobulk tables) | 23 A17 tables CHANGED (`class_fix_tables.md`): organ_adipose, organ_breast, organ_colon, organ_eye_anterior, organ_ileum, organ_optic_nerve, kong_colon, kong_ileum, ia_spleen, ia_jejunum, ia_lung, heart_global, heart_dcm, heart_dcm_healthy, heart_myocarditis, kong_cd, periodontitis, ia_crohn, *_scored twins |
| 588d37d | 10-28 14:50 | heart mast/adipocyte/mesothelial classes; COARSE takes back colonocyte, mucous neck, foveolar, T-helper/T follicular, melanocyte, trabecular meshwork, astrocyte/oligodendrocyte (PI) | DEFAULT-CHANGING | same tables plus eye/optic nerve (melanocyte, glia); `cue_units/build.py` also edited |
| 0dcbea6, a0f8a36 | 10-28 16:13 / 16:55 | A17 anchors relearned on class-fixed tables; `anchors_a17.csv` rewritten, pre-fix copy `anchors_a17_pre_classfix.csv`, `Target(new_classes)` | DEFAULT-CHANGING (anchor table) | every A17 fit/contrast/gate |
| bb8e4c8 | 10-28 17:02 | rule-version stamp in parquet metadata, `StaleTableError`/legacy warnings, kidney_sc mast out of T (103 cells, 0.9 %), COARSE `basophil` word boundary | stamp = infrastructure; kidney_sc fix DEFAULT-CHANGING (small) | kidney_sc tables rebuilt, no compartment change (`rebuild_and_stamp.md`, 48 tables, 0 further changed) |
| 0eb3676 | 10-28 17:04 | rebuild_and_stamp.py; all 48 A17 tables rebuilt and stamped; "no compartment change, anchors stand" | data refresh | A17 data |
| 674224a | 10-29 02:50 | `anchors_a17.csv` changed in one target (organ_bone_marrow, gate fails, no fit uses it); `anchors_a17_pre_phase3.csv` kept | DEFAULT-CHANGING (one unused row) | none |

No `paper/common` code change after 10-29 02:50. Note: `e89defd` named in the brief exists only in the nested clocksync repo (the outer repo has no such hash).

## 2. Per-result staleness

Legend: CURRENT / STALE-POSSIBLE / RERUN-CONFIRMED as defined in the brief. "post-fix" = after 281bb26 (10-28 14:39), anchors relearned a0f8a36 (10-28 16:55), Phase 1 clocksync e89defd (10-28 17:02), phase 3 commit 674224a (10-29 02:50).

### 2a. Figure groups (`paper/results/figure_*`, `extended_data`)

M3 first-pass full runs were 10-24 to 10-27 (pre-class-fix). Last ok run per group below.

| group | last ok run (jobs.tsv end; newest output mtime) | later default-changing change touching it | flag | evidence / note |
|---|---|---|---|---|
| figure_1/method | 10-25 21:11 (76,550 s) | 09d0142/13b57a3 (EM, non-converged hosts) landed 10-25 00:53-20:48, during or just before the run; compute_gtex_benchmark and compute_celltype_truth call `gene_rhythms`. Anchors not used | CURRENT (minor risk: non-converged hosts only) | M3 compared to archive: A31/A36/A37 as archive, TCGA COFE identical, benchmark as legacy except CYCLOPS draw 1 (m3_todo 10-25 ~22:40). Not rerun after 13b57a3 as a group |
| figure_1/species | 10-25 22:35 | none | CURRENT | M3 compared 10-26 ~00:45 |
| figure_1/windows | 10-29 18:09 (rerun, start 10-29 ~13:11; 17,879 s); earlier 10-25 12:32, cohorts part 10-26 22:24 | class-fix, anchors | RERUN-CONFIRMED | `m3_rerun_2026-10-29_compare.md`: no verdict change; R25 T3 p 0.176 to 0.068 (driver ia_spleen deceased r 0.199 to 0.788), still not supported; spleen anchors learned on 9 |
| figure_2/atlas_fits | 10-29 00:26 (run.sh ok 10-28 21:11 after 3 interrupted runs 10-28 14:48-16:12) | none after | RERUN-CONFIRMED | `phase3_compare.md`: 24 of 34 gates pass, same set; 0 flips |
| figure_2/backbone | `run.sh` last 10-26 01:19 (pre-fix); `compute_per_cell` 10-28 22:19; `compute_family_summary` 10-29 00:37 | class-fix | RERUN-CONFIRMED for family medians and per-cell (phase3_compare: learned lymphocyte 0.089 to 0.125, endothelial 0.549 to 0.526; per-cell A18 identical). Other backbone scripts (e.g. `a18_pd_braak`, `a18_thinning`, hierarchy) not rerun after 10-26 | STALE-POSSIBLE for scripts other than the two above | not enumerated script by script (evidence missing) |
| figure_2/gene_rhythms | 10-29 13:11 (rerun); earlier 10-25 18:10, 10-26 09:23 | 13b57a3, e1b6658, class-fix | RERUN-CONFIRMED | compare md: colon R_err up to 0.338; A31b primary fails before and after; heart AA frame genes R > null95 to not; heart LV and adipose v3 p cross 0.05 (post hoc, no criterion) |
| figure_2/embryo | 10-25 20:56 | none | CURRENT | M3 identical |
| figure_2/organs | 10-26 13:30 (43,855 s) | `organ_*` single-cell atlases were not in this group's data (GTEx timed) | CURRENT | M3 compared 10-26 14:31, no verdict change. Not rerun after Phase 1 |
| figure_3/cues_lost | 10-25 23:14 | none found | CURRENT | M3 compared 10-26 00:45 |
| figure_3/gtex_cues | 10-26 00:16 | 613e224 (10-26 13:14): `compute_metabolic.py`, `compute_ventilation_dose.py` call `Dataset.covariate` with the old simulated null | STALE-POSSIBLE | docs record the permutation null only for A19 age (`results_preliminary.md` 2026-10-26, Q10); no later rerun of these two scripts found in jobs.tsv. A164 (10-30, exploratory) is a new analysis |
| figure_4/natural_experiments | 10-26 12:08 | 613e224 (1 h later): `compute_covid_ventilation.py` uses `ds.covariate` | STALE-POSSIBLE | no rerun after 13:14 in jobs.tsv |
| figure_4/resilience | `run.sh` 10-26 10:53 (a17 contrasts, `a17_contrast.csv` mtime 10-26 11:39); `compute_age.py` 10-30 20:36 (rerun on permutation null, old copy `*_2026-10-26.csv`) | class-fix, anchors, Phase 1 for `a17_contrast*.csv`, `a17_gates.csv` | STALE-POSSIBLE for A17 contrasts (the Figure 4 disease contrasts). RERUN-CONFIRMED only for A19 age | `compute_contrasts_a17.py` was interrupted 10-28 16:13 (jobs.tsv, exit 130) and never rerun; `results_preliminary.md` "Not rerun: the Figure 4 contrasts (A153)"; "Figure 4 waits for A153/Q15". The review's pre-fix examples (dcm_healthy heart:lymphocyte d -0.292 on the mast-contaminated class) remain in this file |
| figure_4/spread_multiplicity | 10-26 13:19 (`compute_phase_cancer`); `run.sh` failed 10-26 11:39 (needed a17_contrast.csv) | reads `resilience/a17_contrast.csv` (pre-fix) | STALE-POSSIBLE | inherits resilience staleness |
| figure_5/cancer | 10-29 12:09 (rerun 10-29, 9,375 s); also `compute_cue_following` 10-28 22:21 | class-fix (cue_units), anchors | RERUN-CONFIRMED | A83-A87, A97 byte-identical; point estimates move <= 1.2e-4; interval ends up to 0.030 |
| figure_5/cancer_pipeline | 10-25 12:32 (23,743 s, 72 cohorts) | cancer builders independent of A17 rules; no clocksync default change relevant | CURRENT | M3 only; README states legacy comparison |
| figure_5/lung_breaks | 10-29 18:30 (rerun) | class-fix, anchors | RERUN-CONFIRMED | thinning: dcm_healthy endothelial p_thinned 0.025 to 0.0498, one permutation from 0.05 (borderline, no verdict change) |
| figure_5/psoriasis | 10-25 02:32 (6,776 s) | skin tables not in the changed list; 09d0142 only for non-converged | CURRENT | M3 point estimates equal archive; README A54 swap p at 10,000 swaps |
| figure_5/psoriasis_details | `run.sh` killed_memory 10-25 05:56; remainder `compute_cues_matched.py` 10-26 12:34 and `compute_cue_context.py` 10-26 12:44 (jobs.tsv: ok, `--jobs 1`) | none | CURRENT, but group `run.sh` never exited ok (see 3). Note `compute_cue_context.py` reads cue units; `cue_units/build.py` edited 10-28 | A59 no-MF sensitivity PASS to FAIL as legacy; `psoriasis_details_m3_compare.md` |
| extended_data/calibration | 10-26 11:13 (A10, A19 acceptance); `compute_age_acceptance.py` rerun 10-26 14:57 after 613e224 | none after | CURRENT | A19 acceptance holds within MC error under permutation null |
| extended_data/gate | 10-29 18:30 (rerun, 4,066 s) | class-fix, anchors | RERUN-CONFIRMED | 34 healthy gates, 25 pass, same set; organ_colon acrophase p 0.006 to 0.043 (near 0.05) |
| extended_data/host_identity | `run.sh` 10-26 12:16; `compute_coherence` 10-28 21:10, `compute_anchors` 10-28 21:10, `compute_host_identity` 10-28 22:17, `compute_group_anchors` 10-28 23:40 | class-fix | RERUN-CONFIRMED (the scripts that were rerun); anchor frozen table reproduced 70/70 on old tables (REVIEW C4) | `phase3_compare.md`: colon flag unstable (0.795) to stable (0.915), V17.1/V17.3 shift slightly |
| extended_data/simulations | `run.sh` interrupted 10-26 15:51; r26 per-read tasks done 10-27 01:02; windows cohorts rerun 10-26 22:24 | none | CURRENT | R26 and A25 cohorts rescored after CNV fix |

### 2b. Exploratory amendment folders (last output mtime; `jobs.tsv` end)

Boundaries: class-fix 10-28 14:39, anchors 10-28 16:55, Phase 1 10-28 17:02. A folder whose last run precedes 10-28 17:02 ran on the pre-Phase-1 code, but Phase 1 defaults are bit-identical (only new columns), so only the class and anchor changes matter, and only for scripts that read A17/atlas tables.

| folder | last run | uses class-fixed A17/atlas tables or anchors? | flag |
|---|---|---|---|
| A109_coupling_nulls | 10-28 18:18 | TS, GTEx | CURRENT (post-Phase-1 `compute_h109_responsive`) |
| A110, A111, A113, A132, A133, A135, WP3_scout, WP4_scout | 10-24 | GTEx or timed bulk; no class tables | CURRENT |
| A115_A116_psoriasis | 10-27 | skin | CURRENT |
| A121_embryo_timing | 10-27 22:50 | embryo | CURRENT |
| A122_cue_controls | 10-28 23:33 | cue_units, anchors | RERUN-CONFIRMED (H122.1/2 FAIL/FAIL before and after; Holm 0.547 to 0.645 / 0.507) |
| A123_working_hours | 10-28 16:57 (a123b to 10-28 12:28; truth iii 16:46) | GTEx | CURRENT (docs: "A123 needs no rerun"). Uses `gene_rhythms` (13b57a3 earlier than runs) |
| A124_night_shift | 10-28 17:11 | not verified | CURRENT (probable) |
| A126_menstrual, A127_scanb, A161 | 10-29 | bulk | CURRENT (after all changes) |
| A134_cancer_coherence | 10-29 03:43 (bounded) | cancer | RERUN-CONFIRMED (bounded readout, post hoc; permutation p 0.017 vs bootstrap CI [-0.16, 0.27]) |
| A136_glucocorticoid | no outputs | n/a | not run |
| A138_cancer_estimability | score 10-27 18:19; responsive 10-28 18:17 (after Phase 1); H138.5 bounded 10-29 03:54 | cancer, unbounded ccorr in original `score.py` | STALE-POSSIBLE for the original unbounded `tests.csv`; bounded rescoring (post hoc) exists, H138.5 not estimable on bounded |
| A139_waveform_refit | 10-28 05:17 (run interrupted, exit 130 in jobs.tsv, 26,163 s) | no | outputs from an interrupted run; status of A139 verdict not traced here (evidence missing) |
| A140_dip_cyclops | 10-29 06:42 | no class tables | CURRENT |
| A141_brain_stability | 10-27 04:18 | brain GTEx | CURRENT |
| A142_organ_blocks | script 10-27 14:57; bounded 10-29 03:45 | GTEx | RERUN-CONFIRMED (bounded: r* 0.639 vs q05 0.936 stands) |
| A143_brain_artefact | 10-27 15:05; bounded 10-29 03:52 | GTEx | RERUN-CONFIRMED |
| A144_strong_clocks | split-half 10-27 15:11; cell decomposition 10-28 22:20 | A17 tables | decomposition RERUN-CONFIRMED (D144.2 identical to 3 decimals); split-half STALE-POSSIBLE (pre-fix, not in phase3_compare) |
| A145_two_step | 10-28 10:00 (jobs end 10-28 09:24, 40,111 s) | simulation | CURRENT |
| A146_tme_time | 10-28 12:01 | cancer bulk | CURRENT |
| A147_pull_order | 10-28 16:56 | simulation | CURRENT |
| A148_psoriasis_decomposition | 10-28 10:39 | skin | CURRENT |
| A149_paired_tumour_normal | script 10-28 12:36; bounded 10-29 04:11; per-cancer 10-29 04:11 | cancer | original `tests.csv` unbounded ccorr: STALE-POSSIBLE (REVIEW C1); bounded outputs RERUN-CONFIRMED; H149.1 0.490 and H149.3 0.707 hold, pairing not clock-specific |
| A150_atlas_specificity | 10-28 13:11 (before class-fix at 14:39) | atlas A17 tables, family medians | STALE-POSSIBLE: results_preliminary "Not rerun: ... A150" (phase3 section 4) |
| A151_fibroblast_pairing | 10-28 13:14, 13:27; bounded 10-29 04:01 | cancer | bounded RERUN-CONFIRMED (post hoc); registered statistic broke down (C1) |
| A152_lymphocyte_residency | 10-28 22:47 | heart class fix | RERUN-CONFIRMED (heart lymphocyte 0.556 to 0.327; H152.1 identical) |
| A153_corrected_contrasts | 10-30 20:08 (parked on Q15; eligibility screen script) | A17 tables, current anchors | CURRENT (ran after all changes), but A153 has no primaries yet |
| A154_gate_calibration | 10-30 17:05 (N1 200 reps, 57,643 s) | simulation | CURRENT |
| A155, A156, A157, A158 (sim_2026-10-29) | 10-28 17:21 / 17:22 / 20:32 / 10-29 04:33 | simulation; A155 uses LOO anchor_readout (bc6c419, 20:50, after A155 17:21; opt-in) | CURRENT |
| A159_specific_coupling | 10-28 21:36 | atlases (Phase 1 pseudo sets) | CURRENT (post-fix) |
| A160 | 10-29 03:12 | skin sc | CURRENT |
| A163, A164, A165, A166, A167, A168, A169, A170, A171, A172, A173 | 10-30 to 10-31 | after all changes | CURRENT |
| A19_acceptance_check | 10-26 11:49 | uses simulated null by design (null calibration check) | CURRENT |
| A25_prior_baseline | 10-28 20:31 | n/a | CURRENT |
| psoriasis_robustness_2026-10-28 | 10-28 18:11 | skin | CURRENT |
| WP0_integrity | `anchor_robustness` 10-29 02:47; `q3_anchored` 10-28 23:30; `a21_shift_without_lymphatic` 10-27 14:21 | anchors | anchor_robustness and q3 RERUN-CONFIRMED (265 rows: sign_change 60 to 61; only ia_crohn rows flip); `a21_shift_without_lymphatic` STALE-POSSIBLE (pre-fix) |

### 2c. Archive numbers never rerun in this repo
| item | evidence |
|---|---|
| A5, A6, A7, A12, A13, A29, A30, A32, cycling checks of the 2026-10-04 cancer plan; builders `hlca_tstates_a7`, `natri_tstates_a7` | not migrated (`provisory_plan_paper.md` lines 33-34); stay in the archive (`AMENDMENT_INDEX.md` "archive only") |
| Archive amendments A1-A108c that no figure table cites | migration covers only cited amendments (`migrate_from_archive.py`); no repo-side rerun |
| Datasets built through `from_legacy()` (legacy pickles, "27 legacy-loaded A17 datasets"; `rebuild_and_stamp` rebuilt them under current rules so their A17 compartments are current, but the cell-level inputs, e.g. archive `outputs/` tables of the cancer cell builders, still come from the archive) | plan line 36; `results_preliminary.md` Legacy audit |
| `method compute_cyclops_anchoring.py`, simulations `compute_celltype_truth.py` colon arm (V13) | plan line 66: not run on a subsample; M3 full run only (jobs.tsv has `figure_1/method/run.sh` ok, so presumably run; per-script status not verified) |
| Archive A17 contrast numbers in figure 4 text (the 24 + 20 disease contrasts) | superseded by figure_4 M3 (pre-fix) but not rerun post-fix; A153 parked |

## 3. M3 status

M3 = full rerun of all migrated groups, compared with the archive snapshot (`archive_snapshots_m3/ARCHIVE_SNAPSHOT_<group>.tsv`: 17 snapshot hash lists; snapshots removed after each compare, backups in `$SCCC_DATA/m3_snapshot_backup/`). Plan line 42 (draft 2026-10-27) says "M3 done, every group rerun and compared".

| group | M3 run (jobs.tsv) | compared vs archive | result | pending |
|---|---|---|---|---|
| method | ok 10-25 21:11 (21.3 h) | yes, m3_todo 10-25 | A31/A36/A37 as archive; benchmark as legacy except CYCLOPS draw 1 (run-to-run) | no |
| species, embryo, cues_lost, gtex_cues, backbone, natural_experiments | ok 10-25 to 10-26 | yes (m3_todo, READMEs) | unchanged | backbone partially rerun post-fix (2a) |
| lung_breaks | ok 10-25 20:56 | yes, identical to 1e-9 | matched | rerun 10-29 post-fix |
| gene_rhythms | ok 10-25; M3b 10-26 09:23 | yes | only retried-host rows move; A31b fails, A39 n45 fails as archive; A39 rerun 10-26 12:14 | rerun 10-29 |
| calibration | ok 10-26 11:13 (14 h) | yes | V12 = legacy, A10 rule 15 stands, A19 3b as archive; A19 acceptance FAILS the null criterion under the simulated null, then holds under the permutation null (erratum, Q10) | no |
| gate | ok 10-26 11:13; rerun 10-29 | yes | unchanged | no |
| host_identity | ok 10-26 12:16 | yes (agent, verified) | Natri T p corrected to 0.001 | script reruns 10-28 |
| organs | ok 10-26 13:30 (12.2 h) | yes (`organs_m3_compare.md`) | 25-tissue medians 0.300 vs 1.538 match (earlier F5 mismatch flag was wrong) | no |
| atlas_fits | ok 10-26 12:45 (2.8 h); rerun 10-28 21:11 | yes (`atlas_fits_m3_compare.md`) | 0 gate flips; step7 GTEx join bug fixed | rerun post-fix |
| resilience | ok 10-26 10:53 | yes (`resilience_m3_compare.md`) | Natri T p corrected | A17 contrasts stale post-fix; `compute_contrasts_a17.py` rerun interrupted 10-28 |
| spread_multiplicity | run.sh FAILED 10-26 11:39 (missing a17_contrast.csv); `compute_multiplicity.py` 13:09 and `compute_phase_cancer.py` 13:19 ok | yes (`spread_multiplicity_m3_compare.md`) | erratum R20.4/R23.1 (boundary intervals) | group `run.sh` never exited ok, run piecewise; depends on pre-fix a17_contrast.csv |
| simulations | run.sh interrupted 10-26 15:51; r26 split and rerun ok 10-27 01:02 | yes (`simulations_m3_compare.md`) | R26 and A25 cohorts had used CNV tables (8 cohorts failed); fixed with `cnv=False`, rerun | group `run.sh` never exited ok, run piecewise |
| windows | ok 10-25 12:32; cohorts 10-26 22:24 (50/50 ok); rerun 10-29 18:09 | yes (README) | R25 rescored | rerun post-fix |
| cancer | ok 10-24 22:25 (9,005 s) and 10-29 | yes | A25/R26 cohort tables fixed during M3, A17 step 7 GTEx join fixed | no |
| cancer_pipeline | ok 10-25 12:32 (6.6 h serial) | yes (README) | matched | no |
| psoriasis | first start interrupted 10-25 00:37 (130, stray lane), run ok 10-25 02:32 (34.1 GB peak) | yes (README table of rows marked M3) | point estimates equal the archive; swap p at registered sizes | no |
| psoriasis_details | run.sh killed_memory 10-25 05:56 (peak 30.6 GB; `m3_psoriasis_details.log` line 1449 "run: PSS 30 GB > 30 GB, killing"); first stray start interrupted 10-25 00:39 | yes (`psoriasis_details_m3_compare.md`, README) | no registered verdict changes; A59 no-MF sensitivity p 0.0480 (PASS) to 0.0539 (FAIL), as legacy | `run.sh` never ended ok. The remainder (`compute_cues_matched.py` ok 10-26 12:34, `compute_cue_context.py` ok 10-26 12:44, `--jobs 1`, peak 6.9 / 6.1 GB) completed OK. Answer to the example question: the group was completed piecewise at `--jobs 1`; a full `run.sh` pass was not repeated |
| figure_4 spread_multiplicity phase_cancer | ok | see above | | |

No group is "pending" or "killed" in the M3 sense apart from the two groups whose `run.sh` row shows non-ok in the ledger (psoriasis_details killed_memory; spread_multiplicity failed; simulations interrupted), each completed by split reruns and compared. The post-class-fix rerun covered five groups (cancer, gene_rhythms, windows, lung_breaks, gate) plus atlas_fits, family summary, per-cell, host_identity scripts. NOT rerun post-fix: Figure 4 (resilience A17 contrasts, spread_multiplicity), method, species, organs, embryo, cues_lost, gtex_cues, natural_experiments, calibration, simulations, cancer_pipeline, psoriasis, psoriasis_details, and the rest of backbone (`results_preliminary.md` "Not rerun"; HANDOFF section 2). Plan 8.3 intends a final full regeneration (`m3_rerun.sh` all groups) on frozen code, not yet done.

Missing evidence: per-group compare logs exist only for 7 groups (`*_m3_compare.md`); the others are documented in group READMEs and m3_todo.txt. The plan statement "Snapshots removed (backups in the session scratchpad)" cannot be checked from this repo (scratchpad gone).

## 4. Pending code that could change results

| item | source | could modify results if added? |
|---|---|---|
| Acrophase permutation null in `waveform.py` too stringent / structure-preserving (two-block gene phases) | `clocksync/to_add.txt` (2026-10-21) | YES: would change every gate p-value and hence gate passes; no matching change in code or A154 registration found beyond A154 calibration |
| Default worker cap of 4 for a public repo | `to_add.txt` | No (speed only) |
| Reference choice in disease-only data; make brain analysis a general test | `to_add.txt`; Phase 8; A135/A142 | Possible if the default anchor rule changes |
| Rename shortened variables; remove A-ID references; unused arguments | `to_add.txt` (`_donor_values` unused arg already removed in 1c01018) | No if done mechanically; run the regression check |
| Anchor-rule move into clocksync (`anchor_rule`, identity), docstrings without amendment IDs, type hints, three tutorials, calibration tests as slow CI | ROADMAP 8.4 (3-4 weeks), plan WP8 | Anchor-rule move could change defaults; tests/docs do not |
| Bounded ccorr: clip at 1 or not (`nan_on_flag` does not clip) | HANDOFF section 5 item 6 | YES for any readout built on ccorr (A138/A149/A151 family) |
| `LEGACY_ALIASES` as default in `genes.canonical` (4 A81b cue genes lost to aliases in A160) | HANDOFF section 5 (Phase 8), STATUS 2026-10-29 | YES: psoriasis A81b, cue panels, any gene-set lookup |
| Host-reliability attenuation default (`fit(attenuation=)`, A163 corrections) | 9db4ed8; A163 scored, "df reliability correction fails validation, c stays without intervals" (commit c0d94cb) | Not added as default; none pending |
| API defaults below paper settings (n_boot 200, n_perm) | REVIEW 2.4, `config.py` warnings | Public defaults only; paper scripts pass explicit sizes |
| Phase 8.3 full figure regeneration | ROADMAP 8.3, plan | Not code, but will supersede every row marked STALE-POSSIBLE |
| A153 (corrected Figure 4 contrasts: per-arm gate, untruncated, LOO) | Q15 pending | Will replace Figure 4 A17 contrasts (the stale ones) |
| `to_add.txt` TODO/FIXME in code | `grep TODO|FIXME` over `clocksync/src` and `paper/common` | none found |

## 5. Uncommitted and pre-existing modified files

| location | status |
|---|---|
| outer repo `git status` | one untracked path: `paper/exploratory/A158_bmal1_bias/sim_2026-10-29/output/run1_without_noise_hj0/` (simulation output, git-ignored siblings). No modified tracked files |
| outer repo ignored (not git-tracked) | `clocksync/` (whole nested repo), `.claude/`, `.pytest_cache/`, `paper/common/.pytest_cache/`, all `output/` folders |
| nested `clocksync` repo, modified | `.gitignore`, `pyproject.toml`, `tests/test_additions.py` |
| nested `clocksync` repo, untracked | `src/clocksync/utils.py`, `tests/__init__.py` |
| ownership | HANDOFF section 1: these five are "not ours: never commit them" |
| risk | `utils.py` is untracked inside the package; Lead check (2026-10-31): nothing in `clocksync/src` or `paper/` imports it |
| known untracked outputs | exploratory output folders are mostly git-ignored; results in `paper/results/*/output/` and `m3_snapshot_backup` live outside git. Comparison md files are git-ignored (m3_rerun_2026-10-29_compare.md) |
| unrelated | `scratchpad/compute_a163_n100.py` ran via run.sh (jobs.tsv line, 10-30 19:59). Lead check: it is a symlink to `paper/exploratory/A163_reliability_correction/scripts/compute_a163.py`, so the run is reproducible from git |
