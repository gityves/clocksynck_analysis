# Agent instructions: analyses to strengthen the internal-time paper, figures and narrative

## Purpose, scope and rules

The agent runs nine work packages (WP0-WP8) that close the paper's main weaknesses: a fragile Figure 3 control, no mechanism, no living-human validation of internal time beyond one cohort, and an unexplained cancer break; WP8 adds tissue proxies for night work, breast cancer and social jetlag. Everything uses public data or data already held by the project; no new experiments and no new access requests.

**Positives.** The paper's pre-specified disease positives are the cancer break (A86, A87, A97) and the psoriatic local-clock delay (A54, A78, A81b); WP5 and WP6 strengthen them. The perturbation controls of WP3 (feeding, A111, A112, A129) and WP8 (simulated night shift, A124) are interesting additions, not required for the paper: run them, but no main claim waits on them, except that the A126-A127 cue results are read as claims only if A124 passes.

**Confirmatory family.** Five hypotheses carry the headline claims and are tested with Holm's correction over the five (family-wise alpha 0.05, each one-sided as registered): H109.1 and H109.2 (Figure 3: organs coupled without cues, and the coupling is clock rotation, not a shared displacement), H114.1 (cancer: MYC beyond proliferation, purity-adjusted), H115.1 and H116.1 (psoriasis: trait component, and reversal specific to treatment). Every other hypothesis is secondary (own alpha 0.05, labelled secondary in the text) or descriptive. Report both the raw and the Holm-adjusted p of the five. The registration of each of the five states its membership in the family.

**Out of scope (do not attempt):** brain-bank times of death (HBCC, AMP-PD), the circadiomics left-ventricle request, ICU temperature records, Synapse or dbGaP datasets not already held (Ruzicka 2024, AMP RA, ImmuNexUT), and any wet-lab or animal experiment. GTEx protected phenotypes already held under phs000424 may be used, aggregates only, as in A46-A99.

## Repository and migration (read before any work package)

The project moved on 2026-10-23 from the archive `$SCCC_ARCHIVE` = `/home/maxine/Documents/paychere/scCircadianCancer` (read only) to this repository. Amendment ids A1-A108c, V and D notes refer to the archive's `docs/PLAN_internal_desync.md` and `docs/RESULTS_preliminary.md`; new amendments start at A109.

### Layout

| path | content |
| --- | --- |
| `clocksync/` | the method package (own git repository; public API, tests, tutorials) |
| `paper/common/` | `paths.py` (all paths; read it first), `cli.py` (`--subsample`, `--seed`, `--jobs`), `data.py` (`load(name)` of standard pseudobulk tables), `gtex.py` (GTEx v10, controlled access), `anchors.py` + `anchors_a17.csv` (frozen A17 anchors), `benchmark/` (COFE, CYCLOPS, CHIRAL, TimeTeller), `run.sh` (RAM-capped runner, default 30 GB) |
| `paper/common/datasets/<name>/build.py` | dataset builders: raw public files -> standard tables in `$SCCC_DATA/pseudobulk/`; indexes `INDEX_{atlases,bulk,natexp_cancer}.md` |
| `paper/results/<figure>/scripts/<group>/` | validated analyses, one group (folder of `compute_*.py`, `_shared.py`, `README.md`, `run.sh`) per topic, in the figure whose claim it supports; outputs in `paper/results/<figure>/output/<group>/` (`common.paths.out_dir(group)`) |
| `paper/exploratory/<A1xx_slug>/` | one folder per new registered hypothesis: `scripts/`, `output/` (`common.paths.exploratory_dir`), `README.md`. Moved to a `results/` group once the hypothesis is scored and kept for a figure |
| `paper/data/` | small data (GTEx waveform, CHIRAL DIP, GSE233242 ids; `Gtex_phenotypes/`: controlled access, git-ignored) |
| `paper/docs/` | this plan; `registered_experiments.md` (new amendments, append only); `results_preliminary.md` (new results and errata to archive results, append only); `METHOD.md`; `AMENDMENT_INDEX.md` (archive amendment -> migrated script) |

Groups per figure (`common.paths.GROUPS`): figure_1 `method`, `windows`, `species`; figure_2 `organs`, `backbone`, `embryo`, `atlas_fits`, `gene_rhythms`; figure_3 `cues_lost`, `gtex_cues`; figure_4 `resilience`, `spread_multiplicity`, `natural_experiments`; figure_5 `cancer`, `cancer_pipeline`, `lung_breaks`, `psoriasis`, `psoriasis_details`; extended_data `calibration`, `gate`, `host_identity`, `simulations`. A group stays in one figure folder even when some of its outputs go to Extended Data; the figure tables below say where each panel goes.

### What was migrated (M1, done 2026-10-23, commit cb8b896)

- The archive's clean analyses (`paper/code/`, all on the public clocksync API, each checked against the legacy outputs in its README) for every amendment the figure tables below cite, by `paper/common/migration/migrate_from_archive.py`. Path expressions were rewritten to the same files (report: `paper/common/migration/migration_report.md`). Archive folder -> group: `FOLDERS` in that script.
- Not migrated (not cited): A5, A6, A7, A12, A13, A29, A30, A32 and the cycling checks of the cancer analysis plan (2026-10-04); builders `hlca_tstates_a7`, `natri_tstates_a7`. They stay in the archive (`docs/AMENDMENT_INDEX.md`, "archive only").
- Archive full-run outputs were copied to `paper/results/<figure>/output/<group>/` with their sha256 in `ARCHIVE_SNAPSHOT.tsv` (`snapshot_outputs.py`), so scripts that read another group's tables run before any rerun. A file whose hash still matches is an archive number. Tables with a per-donor column were not copied.
- Still read from the archive (through `common.paths.ARCHIVE` / `LEGACY`): legacy tables that builders convert with `from_legacy()` when `$SCCC_DATA/pseudobulk/<name>.parquet` does not exist; legacy outputs that READMEs compare against (`legacy/desync/**/output/`, `legacy/desync/common_output/`); `outputs/` tables of the cancer cell builders; the legacy code of amendments listed as not ported.

### Migration status and remaining steps

- **M2, subsample equivalence (2026-10-23).** Every migrated compute script run with its README's `--subsample N` and compared with the archive's subsample outputs (same seeds). Result per group: see the table at the end of this section. A script listed as failing or differing must be fixed (or its difference explained in its group README) before it is used in a work package.
- **M3, full reruns, per figure, on demand.** Before a figure's numbers go into text, rerun its groups in full (`run.sh` without `--subsample`, registered permutation and bootstrap sizes, `paper/common/run.sh` for RAM caps) and compare point estimates with the archive snapshot (`ARCHIVE_SNAPSHOT.tsv`; READMEs give the expected values). Record agreement or explained differences in the group README; a full rerun replaces the snapshot files (delete the snapshot rows of rerun files). Symlink the input tables the figure reads (`$SCCC_DATA/pseudobulk/*.parquet`, GTEx aggregates) into `results/<figure>/data/`. Long runs (> 20 min) are documented scripts, launched in the background, at most 4 parallel jobs.
- **M4, cut the data dependency on the archive.** For each dataset a work package uses, save the standard table once (`common.data.save(name, *build_module.from_legacy())`, or `build()` from raw files) so later runs read `$SCCC_DATA/pseudobulk/` only. When every used dataset is saved, `from_legacy()` paths are dead code for this repository.
- **M5, before public release (user's call).** Remove amendment ids from user-facing docstrings, delete `paper/common/migration/`, delete per-donor scratch (none migrated), choose licence.

### Rules for the work packages in this repository

- **Registration.** Write each amendment in `paper/docs/registered_experiments.md` (dated, hypothesis ids, estimand, null, pass rule, data rule, power or MDE) and commit it before scoring. The work-package texts below are drafts of those amendments: copy them, complete the data-dependent parts after the scout (labelled "written after the scout, before scoring"), and commit. Results and run notes go to `paper/docs/results_preliminary.md` with a `src:` pointer to the output table; corrections of archive results (WP0) go there under "Errata to the archive", never into the archive.
- **Code.** New analyses in `paper/exploratory/<A1xx_slug>/scripts/`, following the conventions of the migrated groups (docstring with amendment, hypotheses, estimand, null, outputs; `common.cli.parser`; `--subsample` check first; seeds fixed). Use the public clocksync API; extend clocksync (with a test) when a tool is generally useful, otherwise a helper in `paper/common/`. Do not edit a migrated group for a new hypothesis; import its helpers or copy what is needed.
- **Base code.** Each work package names the migrated scripts it builds on ("Base code"). `paper/docs/AMENDMENT_INDEX.md` maps any other archive amendment to its script.

### M2 result per group

All 22 groups run after the migration and reproduce the archive's subsample outputs file for file (identical to 1e-9 relative, most byte-identical), with the arguments in `paper/common/migration/subsample_args.md`. No migration breakage was found. Differences that remain are not caused by the migration:

| group, script | status | action before use |
| --- | --- | --- |
| backbone `compute_cancer_split.py` (A18 per-cell, cancer) | **resolved (WP0 item 2, 2026-10-23)**: the archive subsample file predates a change of the random draw order (interleaved -> batch); the migrated script equals the current archive script (byte-identical subsample outputs), whose full run matches the legacy groups; not clocksync | M3 full rerun as planned; see backbone README |
| resilience `compute_gtex_sex_bmi.py` (A96, A99, A94 cue-core) | same donors, values differ slightly (log radius 1.4e-2); the archive file predates clocksync edits of 2026-10-22 | explain in M3 (full rerun vs archive full-run snapshot) before A118 or WP0 item 6 use it |
| calibration `compute_age_acceptance.py` (A19 acceptance) | archive not reproducible: seed used Python's salted `hash()`. **Fixed 2026-10-23** (`zlib.crc32`); numbers will differ from the archive | rerun in M3; erratum if the A19 acceptance verdict changes |
| calibration `compute_small_n.py` (A10) | `main()` was never called (archive too), so `run.sh` skipped it. **Fixed** (`__main__` guard); reproduces the archive | none |
| windows `compute_internal_time.py --jobs > 1` | failed: callable clocksync modules could not be pickled. **Fixed in clocksync** (f7f3c86, with a test); `--jobs 4` reproduces the archive | none |
| gate `compute_sheep_reference_primary.py` (A41) | runs; archive subsample N not recovered (13 vs 16 tissues) | compare in M3 on the full run |
| method `compute_cyclops_anchoring.py`; simulations `compute_celltype_truth.py` (colon arm, V13) | not run: > 25-30 min on a subsample | M3 full run only |
| gtex_cues a101_tissues, a46 CIs; species a21 (cellbender log z2), a34_amplitude; organs phase_resid_corr_fast; backbone a18_pd_braak, a18_thinning; cancer a84, a24_cancer; lung_breaks decomposition CIs; psoriasis a81b_h2 (n_boot) | archive subsample files are stale (written before the archive script's last edit) or used other resample sizes; the newer archive runs match | none (M3 compares with the full-run snapshot) |
| psoriasis `compute_cues.py`; psoriasis_details `compute_cues_matched.py`, `compute_cue_context.py` | exceed 30 GB with `--jobs 2` | run with `--jobs 1` (about 6 GB) |

`paper/common/run.sh` forwards script arguments (`RUN_CAP_GB=30 paper/common/run.sh <script> [args...]`).


## WP0. Integrity fixes (no new data; do first)

These remove inconsistencies a reviewer would find in the supplement. None needs a pre-registration; each gets a dated run note in `results_preliminary.md` (corrections of archive numbers under "Errata to the archive"). Code and outputs: `paper/exploratory/WP0_integrity/`.

1. **Missing tables.** Recompute the fibrotic-lung T-cell Fisher combined p (HLCA p 0.027, Natri p 0.009; quoted as 0.002) from the A17 contrast outputs and write it to `paper/exploratory/WP0_integrity/output/` (inputs: the A17 contrast tables `$SCCC_ARCHIVE/legacy/desync/common_output/a17_unified_*.csv`; script `results/figure_4/scripts/resilience/compute_contrasts_a17.py`). Search the archive `RESULTS_preliminary.md` for every number cited by the figure tables below without a `src` pointer; produce the table or list the number as unsupported in the errata.
2. **Malignant cells: per-cell synchrony vs coherence.** The 2026-10-08 copy-number check reads malignant cells as no less synchronous than normal epithelium (dR -0.09); A86/A87 read a loss of clock-shaped share and gene coherence. On the same 79 patients (A87 set), compute side by side: A86 share, A87 coherence and alignment, the anchored coupling, and the per-cell D, lambda and R of `x5_vm_human.py --cnv` (migrated as `results/figure_2/scripts/backbone/compute_cancer_split.py`, A18 part 2; legacy `$SCCC_ARCHIVE/legacy/scripts/x5_vm_human.py`). Simulate a compartment with incoherent gene phases (gene phases drawn uniformly, amplitude unchanged) at the real depths and show what each readout returns. Deliverable: a one-paragraph reconciliation for Methods and a decision on which cancer per-cell numbers stay in the paper.
   *Base code:* `results/figure_5/scripts/cancer/{compute_malignant_gate,compute_decomposition}.py` (A86, A87), `results/extended_data/scripts/calibration/compute_vonmises_calibration.py` (V11 simulations at human depths), `clocksync.cells`, `clocksync.simulate`.
3. **Pre-registration scorecard.** Script that parses every amendment and its run note (archive `docs/PLAN_internal_desync.md` and `docs/RESULTS_preliminary.md` for A1-A108c; `paper/docs/registered_experiments.md` and `results_preliminary.md` from A109; commits from both git histories) into one table: amendment, date, commit, hypothesis id, label (primary, secondary, exploratory, post hoc), estimate, interval, p, MDE, verdict (pass, fail, uninformative, infeasible), output file. This becomes Extended Data Table 1 and the backbone of the cover letter.
4. **Estimator consistency.** Confirm that every main-text coupling uses the A17+ anchored estimator with A22 confirmation. List any figure panel still built on free-fit or pre-A17 values (for example the early near-1 immune couplings, A20 own-noise "damping") and flag it for removal.
5. **CHIRAL reanalysis robustness (A101b, A102b).** Rerun the death-type-matched count comparison with TMM normalisation (as CHIRAL) instead of library-size CPM, and add CHIRAL's sex comparison under the same matching. Report beside the existing numbers; no new hypothesis. *Base code:* `results/figure_3/scripts/gtex_cues/{compute_chiral_prep_counts,compute_chiral_age,compute_covariate_sensitivity}.py`, `results/figure_4/scripts/resilience/compute_gtex_chiral_age.py`.

**6. Gene-symbol audit (after the A85 correction of 2026-10-22).** The legacy A85 script matched gene symbols against the Ensembl `var_names` of `atlas_wg/*/pseudobulk.h5ad`, so every gene resolved to the same column. List every script that reads an Ensembl-indexed h5ad and matches symbols with `gene_columns(a.var_names, ...)` or similar; check first those behind A31b, A39, A50 and A96. Search `paper/` (migrated code, `common/gtex.py`, `common/datasets/`) and the legacy scripts those results came from. For each, rerun with matching on `feature_name`, record which results change, and add a dated erratum to `results_preliminary.md`. *Base code:* A31b `results/figure_2/scripts/gene_rhythms/compute_atlas_gene_rhythms.py`, A39 `.../gene_rhythms/compute_donor_curve.py`, A50 `.../gene_rhythms/compute_gene_map.py`, A96 `results/figure_4/scripts/resilience/compute_gtex_sex_bmi.py` with builders `common/datasets/{islets_a96,masld_pooled_a96}`. All new builds (A118, A122) match on `feature_name` from the start, and every script reports the number of requested genes found.

**8. Anchor robustness as a standard readout.** In datasets with disease donors only, or where the anchors were chosen by rule rather than learned, the host depends on the anchor choice. For every contrast cited in Figures 4-5, report the block agreement r\* between the anchor block and the remaining compartments (`Dataset.agreement`, as the brain-periphery test of A103) and the contrast recomputed on the host of all compartments (leave one out). Add this as a general function to clocksync (with a test) if `agreement` does not already cover it, and as a column of the contrast tables. A contrast whose sign changes with the host is flagged in the figure table.

**7. A22 V17.2 note.** The legacy coherence values came from an empty stratum read as NaN (pairs not centred). Record the centred values from the port (e\_control 0.45-0.55; `results/extended_data/scripts/simulations/compute_a22_colon.py`, README) as an erratum to the archive run note; the verdict is unchanged.

## WP1. Is the coupling of organs without cues a clock? (A109, GTEx; highest priority)

Figure 3 rests on inter-organ coherence C staying at 0.41-0.46 in ventilated and slow deaths (A92, A93). C has no clock-specific null, and A25 found a shared non-temporal post-mortem state lying inside the clock plane. A109 tests whether C is clock-shaped and whether it reflects clocks that keep their phase relations, not a shared push. Two confounds are fixed in advance. (1) Procurement: ventilated donors are beating-heart organ donors, fast deaths are post-mortem collections, so they differ in collection as well as in cues; Hardy 3 and Hardy 4 deaths are post-mortem collections without daily cues, and the fast-vs-Hardy 4 contrast is reported beside every fast-vs-ventilated one. (2) Amplitude: C falls when the clock-plane variance falls relative to noise, with no change in coupling; every C is reported with the per-tissue clock-plane variance (z2 V + 2 psi, A104) of the same group.

*Base code:* `results/figure_3/scripts/gtex_cues/compute_coupling_sources.py` (A92, A93: C, groups, equal-n subsampling, GC score), `results/figure_3/scripts/cues_lost/compute_walltime_hosts.py`, `.../gtex_cues/compute_glucocorticoid.py` (A27 pull direction), `results/figure_1/scripts/windows/compute_nontemporal.py` (A26 non-temporal direction), `common/gtex.py`. Nulls: `clocksync.paired.pseudo_waveforms` (expression-matched pseudo-clock sets with the clock's loadings) and `Waveform.permuted` (permuted acrophases).

**Data.** GTEx v10, the 12 A92 tissues; groups fast (Hardy 1-2), Hardy 3, Hardy 4, ventilated (Hardy 0). Clock states as A92 (GTEx waveform projection, Poisson variance, centred within tissue and group). Equal n by subsampling as A93.

**Nulls (fixed before scoring).**

1. *Pseudo-clock genes.* 200 expression-matched 14-gene sets (`clocksync.paired.pseudo_waveforms`, matched on mean log CPM per tissue, clock and D6 clock-output genes excluded) carrying the GTEx waveform loadings; C recomputed for each.
2. *Permuted-acrophase planes (secondary null).* 200 waveforms with gene acrophases permuted; C recomputed. The clock genes peak in two antiphase blocks, so many permutations keep a clock-like plane and this null is conservative; it is reported, not required.
3. *Shuffled pairing.* Donor labels shuffled within tissue (expected C about 0; sanity check).

**Primary (H109.1, confirmatory family).** In ventilated donors, C\_clock exceeds the pseudo-clock null. Statistic: C\_clock minus the pseudo-clock null median, one-sided empirical p (Holm over the family). Secondary: the same against the acrophase null.

**Secondary.**

- The same in Hardy 3, Hardy 4 and fast deaths (fast is the positive reference).
- *Excess coherence*, C\_clock minus pseudo-clock median, per group with donor bootstrap: does the clock-specific part fall from fast to ventilated, and by how much?
- *Non-temporal component.* In fast deaths, regress each tissue's state on (cos, sin) of time of death; the residual cross-tissue covariance defines the non-temporal shared direction. Report the share of ventilated C along that direction vs orthogonal to it.
- *Processing check.* C restricted to tissue pairs collected by different procedures (exclude pairs dissected together if GTEx annotations allow).
- *Procurement check.* H109.1 and H109.2 repeated in Hardy 4 (post-mortem, no cues); a pass in ventilated donors with a failure in Hardy 4 is read as procurement-dependent.

**H109.2 Rotating clocks or a shared displacement? (primary, confirmatory family).** A shared perturbation (steroids, catecholamines, hypoxia, the agonal phase) that pushes every organ's clock genes the same way would create cross-organ covariance without phase coupling. Two direct signatures separate the readings: running clocks keep the organs' phase relations, and their shared states scatter around the plane; a displacement gives real-valued (0 or 12 h) relations along one axis.

- *Offset preservation (primary statistic).* For each tissue pair, the phase of the noise-corrected complex covariance, arg E[z\_a conj(z\_b)], in ventilated donors vs the same in fast deaths. Statistic: circular correlation (Fisher-Lee) over the 66 pairs; null: tissue labels permuted jointly in the ventilated group (10,000). Pass: correlation > 0, one-sided p (Holm over the family), and the median absolute phase of the ventilated pair relations, folded to [0, 12 h], is not closer to 0 or 12 h than in fast deaths (one-sided permutation p > 0.05 for "closer").
- *Circularity (secondary).* Pooled ratio |E[z\_a z\_b]| / |E[z\_a conj(z\_b)]| over pairs (pseudo-covariance over covariance): near 0 for rotating clocks with scattered phases, near 1 for a displacement. Compared with simulations at the ventilated design (same donors, tissues, depths, 200 each): (i) clocks rotating with uniform phases at the observed coupling; (ii) a pure displacement of the observed size along the A27 glucocorticoid direction (17.7 h on the GTEx clock); (iii) mixtures at 25 / 50 / 75 % displacement. Report the mixture fraction that best matches the data.
- *Descriptive.* Major-axis share of the host covariance; Rayleigh and axial R of donors' host phases; angle between the major axis and the glucocorticoid direction, modulo 12 h.
- Reading: offset preservation failing (relations collapse to 0 / 12 h, or no correlation with the fast-death relations) means the coupling is a shared state, not phase coupling; the coupling claim is removed whatever H109.1 shows.

**A109b, independent replication in Tabula Sapiens (secondary).** Tabula Sapiens donors are brain-dead, ventilated organ donors with many organs each, outside GTEx; the archive's multi-organ gate already passes on them (healthy-atlas table, p\_acrophase 0.017; check the number in the archive results before citing). Per donor x organ pseudobulks of stromal-vascular cells (fibroblasts, endothelium, mural cells; the A17 anchor families), organs with >= 5 donors. Statistics: C against the pseudo-clock null (as H109.1) and offset preservation against the GTEx fast-death pair relations for the organs both share (as H109.2). One-sided p < 0.05 each. Feasibility rule: >= 6 organs shared by >= 5 donors, otherwise descriptive. Builder: `common/datasets/cells_census` or a new `tabula_sapiens_organs` builder (Census, feature\_name matching).

**Reading.** H109.1 and H109.2 pass: organs share a clock-shaped state without cues and keep their phase relations; Figure 3 stands (with A109b as replication if it passes). H109.1 passes and H109.2 fails: organs share a donor state, not a running clock; Figure 3 is reframed as "organs share a donor state", and WP2 is dropped from the main text. H109.1 fails: remove the coupling claim. A pass in ventilated donors only (Hardy 4 failing) is reported as procurement-dependent.

## WP2. What keeps organs coupled? Memory, active coupling, and exclusions (A110)

*Base code:* A104 `results/figure_2/scripts/backbone/compute_oscillator_gtex.py`, A33 `results/figure_3/scripts/cues_lost/compute_ventilator_tissues.py`, Eraslan stromal-vascular fraction `results/figure_2/scripts/atlas_fits/compute_step7_confirmations.py`, A46 `results/figure_3/scripts/gtex_cues/compute_ventilation_dose.py`, WP1 code.

**Superseded (2026-10-23):** after A109/A109c and the PI's framing decision (Q1), WP2 is registered as A110 in `registered_experiments.md` on Hardy 4 / Hardy 3, with Hardy 0 organ donors as the SCN contrast. The draft below is kept for its base code and tests. Run only if H109.1 and H109.2 pass. Three explanations remain after WP1: memory of the last common phase, active coupling among peripheral organs that free-run as a group, or a residual external cue. A110 decides between the first two and tests the cues brain-dead ICU patients are unlikely to have. Data, groups and states as A92/A93; protected phenotypes aggregates only.

**H110.1 Does coupling decay with time without cues? (primary of A110, secondary for the paper).** Ventilation days go with case mix (cause of death, sepsis, steroids), few donors pass 10 days, and the date of brain death or of cue loss is not recorded, so a pre-specified contrast is preferred to model selection between curves.

- Statistic: C in the top vs the bottom quintile of ventilation days (equal n), difference with donor bootstrap; null: day labels permuted within cause-of-death class (Hardy 0 cause categories of the protected phenotypes; classes with < 10 donors pooled), 2,000. Prediction of memory alone: top < bottom. Report the MDE of the difference at this design.
- Reading: no decay with the MDE below the C drop that a memory model gives at the top-quintile median days for a period spread of 0.3 h/day (computed by simulation before scoring, as A92's power note) supports maintained coupling; a decay supports memory.
- Descriptive: C in sliding windows of 60 donors ordered by days, with the M0 memory curve, C(t) = C0 exp(-(2 pi sigma\_T t / 24)^2 / 2), and the plateau curve, C(t) = C\_inf + (C0 - C\_inf) exp(-t / tau), fitted for display; the upper bound on sigma\_T under M0. Context only: published spreads of human fibroblast periods (Brown et al. 2005) and peripheral desynchrony in SCN-lesioned mice (Yoo et al. 2004; Tahara et al. 2012). Verify each paper's numbers before citing.

**H110.2 Signatures of active peripheral coupling (secondary, Holm over a-c).**

- (a) *Strong clocks carry the network.* The A104 statistic, Spearman(own amplitude A, coupling c) over peripheral tissues, in ventilated donors. Prediction rho > 0, as in fast deaths; also compare slopes.
- (b) *Stroma.* Retained coupling per tissue (ventilated / fast, A33 machinery) against the Eraslan stromal-vascular fraction; Spearman > 0, one-sided permutation p < 0.05. Within-tissue check: donor-level stromal-vascular marker score (COL1A1, COL1A2, DCN, PECAM1, VWF, RGS5).
- (c) *Hubs.* Split ventilated donors in halves (seed fixed). In half 1, define hubs as the 4 peripheral tissues with the largest own amplitude A. In half 2, for every non-hub tissue, coupling to a host of all other tissues minus coupling to a host without the hubs; same in fast deaths. Statistic: mean drop, ventilated minus fast; null = removing 4 random non-hub tissues. Prediction > 0: hubs carry more of the coupling once external cues are gone.

**H110.3-H110.5 Exclusion tests (each reported as excluded or not).**

- *Temperature (H110.3).* Unlikely: brain-dead donors tend to lose hypothalamic thermoregulation and ICU care controls temperature. Scores: cold-inducible CIRBP, RBM3; heat shock HSPA1A, HSPA1B, HSPH1, DNAJB1; residualised on ischaemic time and its square per tissue. Excluded if the drop in C after regressing them out is below 0.05 (the glucocorticoid benchmark) with the upper bound below 0.05 and not above 200 random expression-matched 6-gene scores.
- *A functioning SCN (H110.4).* Already weakened: donors with a primary brain lesion are as coupled as the others (A92). Coupling of Brain - Hypothalamus vs the median other brain region to the peripheral host, in fast deaths and Hardy 3-4. SCN-driven body supported only if the hypothalamus exceeds that median in Hardy 3-4 with bootstrap lower bound > 0.
- *Feeding axis (H110.5).* Groups: digestive (colon transverse and sigmoid, stomach, oesophagus mucosa and muscularis, liver, pancreas) vs other A92 tissues. Within-digestive minus digestive-to-other coherence, ventilated minus fast; group-label permutation. Excluded if not > 0 at p < 0.05.

**Deliverable.** A Figure 3 panel showing the three explanations with their verdicts. If WP1 passes, H110.1 passes and at least one of H110.2a-c passes, the abstract sentence becomes: the human periphery behaves as a self-coupled oscillator network that holds together without the brain and without daily cues, and free-runs away from wall time. Related precedent to check before citing: hepatocytes stay synchronised without the SCN or other extrahepatic clocks (Sinturel et al. 2021).

## WP3. Feeding and genetic controls: organ desynchrony, the lock measure, a tissue without its clock (A111, A112, A129; public mouse data; interesting, not required)

A49 (Deota) could not serve as a positive control because night-time feeding produced almost no desynchrony. Daytime and inverted feeding do. Arrhythmic feeding gives a timed truth for the cue/locked distinction.

*Base code:* mouse orthologue convention (A48) `results/figure_4/scripts/resilience/compute_mouse_ageing.py`; Deota builder `common/datasets/mouse_ageing_deota` (A49); timed truth `clocksync.timed` (`Dataset.timed`, `timed_loadings`, `cosinor`). New builders go to `common/datasets/<name>/build.py` with a `PROVENANCE.md` under `$SCCC_DATA/<name>/`.

**Datasets (scout first; take accessions from the papers).**

| dataset | design | tissues | role |
| --- | --- | --- | --- |
| [Manella et al. 2021, Nat Metab](https://pmc.ncbi.nlm.nih.gov/articles/PMC7611072) | ad libitum vs daytime feeding, liver-clock mutants, around the clock | liver, WAT, lung, quadriceps, kidney, heart | primary positive control |
| [Xin et al. 2021, iScience](https://pmc.ncbi.nlm.nih.gov/articles/PMC8050734) | day- vs night-restricted feeding, 1 week, 24 h | several peripheral tissues (TPM only) | replication; pseudo-counts as V25.5 |
| [Greenwell et al. 2019, GSE118967](https://omicsdi.org/dataset/geo/GSE118967) | ad libitum, night-only, arrhythmic (1/8 of food every 3 h) | liver | lock validation |

**A111, positive control (Manella primary, Xin replication).**

- Units: one mouse per time point and condition; tissues as compartments; host anchored on all tissues except liver (all except the tissue read, leave-one-out). Human GTEx waveform with mouse orthologues (A48 convention); ZT folded to 24 h.
- Timed truth first: per tissue and condition, cosinor phase of the clock-plane state on ZT; truth change = phase(daytime feeding) minus phase(ad libitum) per tissue, relative to the tissue mean.
- H111.1 (offsets): the untimed change in tissue offsets relative to the host correlates with the timed change, circular correlation (Fisher-Lee) > 0, exact permutation p < 0.05 over tissues.
- H111.2 (coupling): liver coupling to the non-liver host falls under daytime feeding (untimed d < 0), in the direction of the timed liver-vs-rest R^2 change. Permutation of condition labels within ZT.
- Informativeness rule (as A49): the timed truth must show a liver phase shift >= 3 h relative to the other tissues; otherwise report uninformative.
- H111.3 (secondary, a positive control for a break): in the liver-clock mutant arm under ad libitum feeding, liver coupling to the non-liver host is lower than in wild type (untimed d < 0 at a common host spread, d\*), or the liver fails the cross gate while the other tissues pass; condition labels permuted within ZT, one-sided p < 0.05. Check the mutant's design in the paper first (which clock gene, which cells) and write it into the amendment.
- Descriptive: does the rest of the body stay together without a liver clock (C among non-liver tissues, mutant vs wild type)?

**A112, lock validation (Greenwell).**

- Single tissue, so the host is the liver's own clock-plane state per sample. For each gene, lock = |coefficient on the host| in arrhythmic feeding divided by the same in ad libitum (A81 definition, both timed and untimed versions).
- Gene sets fixed now: mouse orthologues of the GTEx A81b cue set and locked set (40 each), core clock genes as reference.
- H112.1: cue-set orthologues lose more lock under arrhythmic feeding than locked-set orthologues (median log lock ratio, cue minus locked < 0), gene-label permutation p < 0.05.
- Truth check: Greenwell's own rhythmicity calls (rhythmic in ad libitum, lost in arrhythmic) should be enriched among cue orthologues.
- Reading: a pass shows the GTEx-derived cue/locked split corresponds to feeding-driven vs clock-driven rhythms in a controlled experiment, across species.

**A129, a tissue without its clock (scout, then register).** No work package so far shows that the method reports a break when a tissue's own clock is removed.

- Scout public multi-tissue transcriptomes of tissue-specific clock knockout or reconstitution mice with >= 2 tissues per mouse: e.g. Koronowski et al. 2019 (liver-only Bmal1 reconstitution), Welz et al. 2019 (epidermis-only reconstitution), and tissue-specific Bmal1 knockouts. Verify accessions and designs; use only open matrices.
- If found (>= 6 mice per genotype with the tissues of interest): H129.1, the tissue lacking its clock has lower coupling to the host of the other tissues than in controls (d\* < 0), or fails the cross gate, while a companion tissue with an intact clock does not change (specificity on the same permutations). Timed truth where ZT is recorded. Register after the scout, before scoring.
- If none: state it as a limitation; A111 H111.3 is the only genetic control.

## WP4. Internal time in living humans (A113)

V25.5 (breast surgery, 42 patients) is the only living-hospital validation of per-patient internal time. A113 adds living timed cohorts, including one where the feeding rhythm is flattened.

*Base code:* A28 E1 projector `results/extended_data/scripts/gate/compute_extended_clock.py`; A35 within-subject arm `results/figure_2/scripts/gene_rhythms/compute_blood_clock.py`; A25 window estimator `results/figure_1/scripts/windows/{compute_internal_time,_internal_time}.py`; V25.5 `.../windows/compute_breast_timed.py`.

**Datasets (scout first).**

- GSE108539 (timed skeletal muscle) and GSE104674 (timed adipose): listed in the plan as timed human bulk sets, never used. Confirm subjects, sampling times, repeated biopsies per subject, platform and count type.
- Geneva human muscle with hourly isocaloric feeding (Perrin et al. 2018; Loizides-Mangold et al. 2017): serial biopsies every 4 h. Use only if a transcriptome matrix is public.
- Smith et al. 2024 (JCEM, hourly liquid meals, muscle every 4 h): use only if genome-wide or at least the 14 clock genes are public.

**Estimators.**

- *Single sample, single tissue:* clock-plane phase from the 14-gene GTEx projector (A28 E1); error after a common offset.
- *Within subject:* states centred on the subject mean (A35 arm c); ordering of each subject's biopsies.
- *Cohort window:* A25 axis and curvature on subsets mimicking hospital windows (08-18 h, 09-17 h) built from the timed samples, flip resolved by the hospital prior (V25.5 rule).

**Primary (H113.1), per dataset.** Circular correlation of estimated phase with recorded clock time > 0.3, permutation p < 0.05 (subject-level permutation when biopsies are repeated). Pass in >= 2 datasets.

**Secondary.**

- Median per-sample error after the common offset, against V25.5 (1.1 h) and GTEx single-tissue (2.6 h).
- Window subsets: Spearman of the tangent coordinate with clock time; direction correct under the hospital prior.
- Hourly-feeding cohort: does the muscle clock order samples as well as under normal meals? Descriptive, the human counterpart of ICU feeding.
- Common offset (internal frame minus wall clock) per dataset, next to V25.5 (-2.2 h).

**Reading.** Passing in two or more living cohorts turns internal-time estimation from biopsies into a validated readout and supports the translational paragraph of Figure 3.

**A123, the GTEx organ map from working-hours deaths only (written before scoring; the 8-tissue window results of the benchmark have been seen).** Reproduce the CHIRAL organ figure using only donors who died in working hours, with times hidden.

- *Base code:* `results/figure_1/scripts/method/compute_gtex_benchmark.py` (benchmark truth, CHIRAL runner, windows), `common/benchmark/chiral.py`, `results/figure_2/scripts/organs/compute_organs_timed.py`.
- Donors: GTEx v10 fast deaths (Hardy 1-2) who died 08-18 h; every tissue with >= 20 such donors.
- Estimates: untimed reference-free tissue phase offsets and relative amplitudes; per-donor internal time (A25 tangent coordinate, flip by the hospital prior).
- Truth (i): per-tissue cosinor on recorded time of all fast deaths, reduced by complex SVD (the benchmark truth). Truth (ii): CHIRAL's published per-tissue phases and amplitudes (naef-lab/CHIRAL, commit 1f78a9a).
- H123.1: offsets vs truth (i), circular correlation > 0.7, tissue-permutation p < 0.05. H123.2: relative amplitudes vs truth (i), Pearson > 0.7. H123.3: donor ordering, Spearman of the tangent coordinate with recorded time > 0.4 within the window, with the direction correct.
- Secondary: the same against truth (ii); CHIRAL (pychiral) run on the same 08-18 h donors as the comparator that assigns phases.
- Deliverable: a CHIRAL-style organ figure built from working-hours donors only, for Figure 1.

**A130, internal time against known chronotype epidemiology.** Chronotype advances with age after early adulthood and is earlier in women until about 50 (Roenneberg et al. 2004; verify before citing). In GTEx the age advance was already seen (A94: about 0.2-0.5 h per decade, peripheral tissues) and the sex phase difference was estimated (A99: +0.35 h, SE 0.31): both are reported as consistent with that epidemiology, labelled post hoc.

- H130.1 (secondary, pre-registered replication outside GTEx): per-subject internal minus wall time against age, pooled over the living cohorts of A113 and V25.5 that record age (cohort as stratum, A25 internal time with the hospital prior where sampling is in working hours, timed phase where times span the day). Slope < 0 (earlier with age), one-sided permutation of age within cohort, p < 0.05. Feasibility: >= 30 subjects with age and an age range >= 25 years, otherwise descriptive.
- Descriptive: women minus men, same estimand, by age above or below 50.

**A131, sun time within time zones (scout, exploratory).** Human chronotype follows sun time within a time zone (Roenneberg et al. 2007; verify). GTEx samples carry a collection-centre code (`SMCENTER`). Scout whether public GTEx documentation maps the codes to locations, and whether the donor's place of death is near the centre (OPO catchment). If yes: centre longitude offset from its time-zone meridian against the donor's host phase minus time of death (fast deaths; A93 H93.3 machinery), centre as the unit, descriptive with interval. Centres differ in procedures, so any trend is a candidate, not a finding. Aggregates per centre only.

## WP5. Why malignant cells lose clock coherence (A114)

**Superseded (2026-10-25, PI, Q9 e):** replaced by A138 (what can be estimated about the clock in tumours; CoFE /
CYCLOPS orderings tested by simulation), registered in `registered_experiments.md`; run before WP6.

A86, A87 and A97 show that malignant cells keep amplitude but lose gene coherence, more so in proliferative tumours. A114 asks whether oncogenic transcription, copy number or proliferation itself explains it, and whether proliferating normal cells lose coherence too. In bulk tumours, proliferation rises with purity, and a mixture of incoherent malignant cells with coherent stroma becomes less coherent as purity rises, so every bulk coherence contrast is adjusted for purity.

*Base code:* A87 and A97 `results/figure_5/scripts/cancer/{compute_decomposition,compute_proliferation,compute_delay}.py`, A24 `.../cancer/compute_cycling_a24.py`, A86 `.../cancer/compute_malignant_gate.py`; builders `common/datasets/{tcga_paired,cancer_sc,cells_cancer,cells_cancer_cycle,colon_chrono}`; A21 `results/figure_1/scripts/species/compute_colon_injury.py`.

**Data.** TCGA (recount3, the 11 local projects), plus open GDC gene-level copy number and published open aneuploidy and purity/ploidy scores (download only open-access tiers). A21 colon chrono-atlas (known time). A86 single-cell cancer studies. Census lymphoma (67b6b9ac), ovarian (b252b015) and prostate (881e0e6b) cancer sets.

**H114.0 Purity check of A97 (secondary, run first).** Purity per tumour from open estimates (ESTIMATE, or ABSOLUTE purity where open; record the source). The A97 proliferation contrast of coherence within purity tertiles (tertile as stratum, labels permuted within cancer x purity tertile, 1,000), and coherence against purity alone. Reading: the A97 association holds within purity tertiles (one-sided p < 0.05) = proliferation, not purity; otherwise the bulk proliferation claim is restated as purity-confounded and rests on the single-cell malignant arm (A86/A87, same contrast with per-cell cycling fraction).

**H114.1 MYC beyond proliferation (confirmatory family).** Mechanism to test, known in the literature: MYC represses the clock through BMAL1 (Altman et al. 2015; Shostak et al. 2016; verify before citing).

- Score: Hallmark MYC\_TARGETS\_V1 and V2, clock and A97 G2M genes removed; residualised on the A97 proliferation score and on purity within cancer.
- Statistic: mean over cancers of gene coherence (A87) in the top minus the bottom tertile of the residual MYC score; tertile labels permuted within cancer (1000).
- Pass: < 0, one-sided p, Holm over the confirmatory family. Companion: the A97 proliferation contrast within MYC halves. Single-cell check (secondary): the same contrast in malignant cells of the A86 studies, MYC score per patient.

**H114.2 Copy number at the clock loci.**

- Per tumour, copy-number deviation at the 14 clock-gene loci; and the aneuploidy score.
- (a) Per-gene dosage: for each of the 14 clock genes, regress the gene's deviation from the tumour's clock-plane fit (A87 residual) on its own copy-number deviation, pooled over cancers (cancer as stratum, purity as covariate); and recompute the tumour-normal coherence gap after removing the dosage-explained part of each gene. (Excluding every tumour with an arm-level change at any of 14 loci would leave too few tumours.)
- (b) Coherence by aneuploidy tertile within cancer, proliferation held by stratification.
- Pass for "copy number explains it": (a) the gap shrinks > 50 % after the dosage correction, with bootstrap support. Otherwise copy number is ruled out as the main cause.

**H114.3 Do proliferating normal cells lose coherence? (known time)**

- A21 colon, undamaged and regenerating arms: per arm, coherence (A87) of stem/TA, of colonocytes and of secretory cells, from timed and untimed fits; A24 author labels for proliferating vs non-proliferating.
- Also A86 adjacent-normal epithelium split into cycling and non-cycling cells (A24 top-quartile rule), human.
- Reading: if proliferating normal cells keep coherence while malignant cells lose it, the break is malignancy-specific; if they lose it too, it is a property of proliferation.

**Descriptive.**

- A86 share and A87 components in the three new Census tumour sets (no adjacent normal in lymphoma; compare with the matched normal tissue in Census where available).
- Per-gene view: which clock genes lose their phase relation first (rank of per-gene phase deviation in high- vs low-proliferation tumours).

## WP6. Psoriasis: trait or state, and is the reversal specific to treatment? (A115, A116)

The lesional local-clock delay is replicated (A54, A78, A81b), but two questions a reviewer will ask are open: is part of it already in the patient's uninvolved skin, and does treatment reverse it more than placebo (A56 post hoc p 0.28)?

*Base code:* A54 `results/figure_5/scripts/psoriasis/compute_replication.py`, A78 / A81b core-cue estimator `.../psoriasis/{_cues,compute_cues,compute_cues_gtex}.py`, A56 `.../psoriasis/{compute_treatment,compute_dose}.py`, A57 `.../psoriasis/compute_response.py`, A61-A69 `results/figure_5/scripts/psoriasis_details/`; builders `common/datasets/skin_paired_{discovery,heldout,treatment,other}`; GTEx cue sets `results/figure_3/output/gtex_cues/a81b_*.csv`.

**A115, trait vs state (data already on disk).**

- Data rule: the A54/A53 psoriasis and atopic dermatitis series that also contain healthy-control skin from different persons (check each series' metadata; for example GSE121212 is expected to include healthy controls). Same platform within a series.
- Estimand: the within-sample core-vs-cue offset (A78 estimator: arg E\[c\_core c\_cue\*\] over samples of one group), using the GTEx cue set (A81b) and, as sensitivity, the A78 skin systemic set. This needs no pairing, so unpaired healthy controls can be compared.
- Check first that the estimand reproduces the paired result: lesional minus uninvolved core-cue offset should be about -1.4 to -1.7 h.
- H115.1: uninvolved psoriatic minus healthy core-cue offset < 0, pooled over series (A54 weights), one-sided group-label permutation within series p < 0.05. A pass means part of the delay is a trait of psoriatic skin.
- Descriptive: the same for atopic dermatitis non-lesional vs healthy (specificity).

**A116, placebo-controlled pooling.**

- Scout GEO for psoriasis trials with a placebo arm and lesional plus non-lesional biopsies at baseline and at a follow-up visit (any systemic or biologic agent). Data rule: >= 10 placebo patients with both visits per series, log-scale or count matrix with the 14 genes. GSE117468 and GSE201827 are included.
- Estimand: change of the lesional offset (baseline minus follow-up, A56 estimator) in treated minus placebo patients, pooled over series.
- H116.1: treated minus placebo < 0, arm labels permuted within series (2000), one-sided p < 0.05.
- Secondary: the same with PASI75 response in place of arm, within treated patients only (separates response from arm, which A57 could not).

**Descriptive, keratinocyte route.** In series with healthy controls, does the keratinocyte differentiation score (A61 sets) differ between uninvolved psoriatic and healthy skin? It links A115 to the transit model.

## WP7. Timed human brain and the remaining fragile results (A117-A122)

*Base code:* A3 `results/figure_2/scripts/organs/compute_chen_age.py`, V25.3 `results/figure_1/scripts/windows/compute_chen_windows.py`; A96-style cue and locked builders `common/datasets/{islets_a96,masld_pooled_a96}` and `_atlas/heldout.py`; A19 procurement `results/figure_4/scripts/resilience/compute_procurement.py`; A108 `results/figure_2/scripts/embryo/`; A85 `results/figure_5/scripts/cancer/compute_cue_following.py`; A98 `results/figure_5/scripts/lung_breaks/compute_decomposition.py`.

**A117, timed brain regions (public only).** The brain claims rest on GTEx, sheep and one timed cortex study (Chen 2016).

- Scout: Li et al. 2013 (PNAS; several brain regions with time of death, controls and major depression) and Seney et al. 2019 (Nat Commun; prefrontal cortex with time of death, schizophrenia). Use a dataset only if its expression matrix and times are open on GEO; otherwise report infeasible and stop.
- H117.1 (validation): untimed between-region offsets vs timed region phases, circular correlation > 0, p < 0.05 (as A36 readout 2).
- H117.2 (age): log2 amplitude older vs younger donors, pooled over regions, < 0 (replicates A3 on independent donors).
- Descriptive: schizophrenia vs controls (timed amplitude and between-region coupling), as a timed check of the unreplicated HBCC drops; MDD vs controls.

**A118, amplitude behind the breaks.** A104 shows coupling tracks amplitude, so each coupling drop needs an amplitude readout that is not referenced to the host.

- Rebuild pseudobulks with the A81b cue and locked genes for Natri, HLCA fibrotic, PD DMV and AD lesional skin (A96 approach).
- Per break: log(core/cue) and log(core/locked) case minus control, with the A96 permutation and bootstrap.
- Reading: a drop in core/cue with locked/cue unchanged = the clock itself is damped in that cell type; no change = the coupling drop is relative only.

**A119, lung T-cell procurement (scout, then test or demote).**

- Scout open single-cell lung cohorts (CELLxGENE, GEO) with fibrotic ILD and controls procured the same way (both surgical, or both explant/biopsy), >= 15 per group with T cells.
- If found: A19 P1b design (d and d\*, permutation within study). If none: state in the paper that the T-cell result cannot be separated from procurement, and move it to Extended Data.

**A120, embryo replication (scout, then test or demote).**

- Scout for public mouse embryonic or fetal transcriptomes (bulk or single cell, any tissue, E15-P0) sampled at >= 4 recorded clock times.
- If found: the A108 H108.1 statistic (R of host phase vs clock time, one-sided permutation), pre-registered before scoring. If none: the embryo panel's place is decided by A121 alone.

**A121, the embryo clock against the known timing of fetal clock emergence (Qiu 2024 data already built; written before scoring).** The A108 result agrees in timing with prior work: no cell-autonomous oscillation in E10 hearts and clear oscillation at E18 (Umemura et al. 2017), and fetal peaks stabilising near dusk around E15.5 in in utero PER2::LUC imaging (a 2023 preprint). A121 tests two predictions of that picture; H121.1 must pass, and H121.2 must pass or be uninformative, for the Figure 2 embryo panel to stay in the main text.

- *H121.1, onset within the window.* Alignment with clock time is stronger in E17.0-E18.75 than in E15.0-E16.75. Per half: the H108.1 statistic R (8 embryos each, batch intercepts, host refitted per half). Statistic: R(late) minus R(early); null by permuting clock times within batch in both halves (10,000), one-sided p < 0.05. Secondary: the same per lineage (hepatocytes, mesenchyme), descriptive.
- *H121.2, phase against the known fetal timing.* Convert the common offset of the embryo host (host phase minus clock time) into mouse ZT with a human-to-mouse frame offset calibrated on every timed mouse dataset in hand that uses the same orthologue projection (A21 colon, A48 mouse ages, A72 mouse skin; V25.4 alone failed as specified): offset = inverse-variance weighted circular mean of the per-dataset offsets, with its SE by bootstrap over datasets and samples. The light schedule is not recorded: report the predicted ZT of the embryo clock peak for lights-on at 06:00, 07:00 and 08:00. Pass if, under the standard schedule (07:00), the published fetal peak lies within the 95 % interval of the predicted ZT (interval from the offset SE and the embryo host SE). If the calibrated offset's SE exceeds 1.5 h, or the datasets disagree by more than 3 h, H121.2 is uninformative and reported descriptively. Verify the published fetal peak in the preprint before fixing it in the amendment.
- *Heart, stated in advance.* Cardiomyocytes were not rhythmic in either half (A108c; p 0.22, 0.28). The paper says "consistent with" the heart literature, not "replicates".
- *Reading.* H121.1 passes and H121.2 passes (or is uninformative): the Figure 2 embryo panel, labelled exploratory, framed as consistent with and extending known late emergence. H121.1 fails, or H121.2 fails: the embryo goes back to Extended Data.

**A122, controls for the corrected A85 result (written after seeing it: immune cue-following index 0.182, structural 0.095, difference +0.087, one-sided p 0.002; tables in `paper/results/figure_5/output/cancer/`).** Immune classes have a small host loading, so the index |b| / (|a| + |b|) can rise mechanically. These controls decide whether the result carries a Figure 2 panel. Gene matching on `feature_name`.

- *H122.1, locked-set control (primary).* Rerun A85 with the A81b locked set in place of the cue set. Statistic: (immune minus structural, cue index) minus (immune minus structural, locked index), pooled over atlases as A85. Pass: > 0, one-sided p < 0.05 (atlas bootstrap).
- *H122.2, pseudo-cue sets.* 200 expression-matched random 40-gene sets (excluding clock, cue and locked genes). Pass: the observed immune-minus-structural difference exceeds their 95th percentile.
- *Consistency (descriptive).* Sign of the difference in each of the 9 atlases; within atlas, the class index against log pseudobulk depth, and the difference after adjusting for depth.
- *Composition (descriptive).* Variance of the cue state's own expression by class: do immune classes simply express the cue genes more?
- *Reading.* Both pass: Figure 2 panel, "immune cells take relatively more of their timing from body-wide signals". Otherwise Extended Data, described as a consequence of low host coupling.

**Optional, if time allows: the neural sign reversal (A104 post hoc).** Pre-register as a test: across neural compartments of independent atlases not used to see the pattern, Spearman(A, c) < 0. Report as a candidate if it passes.

## WP8. Circadian disruption, breast cancer and body mass: tissue proxies in public data (A124-A128)

No public tissue dataset records night work or sleep timing. WP8 uses the processes thought to link them to disease, which public data do record: imposed misalignment, the menstrual cycle, proliferation, body mass, and the weekly social clock. The readout is the core-vs-cue offset (A78, A81b): a local clock lagging or leading the body-wide cues in the same samples. It cannot show causation; it provides the tissue intermediate that epidemiology lacks. Run A124 first: the other cue-based tests are read as claims only if it passes.

*Base code:* core-cue offset estimator `results/figure_5/scripts/psoriasis/_cues.py` (A78/A115 estimator, array convention of A53), GTEx cue and locked sets `results/figure_3/output/gtex_cues/a81b_*.csv`; breast and endometrium builders `common/datasets/{endometrium,endometrium_scored,breast_hbca_genome,organ_breast,breast_timed_gse233242}`; A97 `results/figure_5/scripts/cancer/compute_proliferation.py`; A93 H93.3 machinery `results/figure_3/scripts/gtex_cues/compute_coupling_sources.py`.

### A124, positive control: misalignment imposed by a simulated night shift

- **Data.** [Wefers et al. 2018, GSE106800](https://omicsdi.org/dataset/geo/GSE106800): 14 men, skeletal muscle after a control protocol and after a 12 h shift of the behavioural cycle (microarray). Scout first: biopsies per subject, sampling times, platform. Secondary: [Kervezee et al. 2018, GSE107537](https://omicsdi.org/dataset/geo/GSE107537), PBMC every 4 h before and after a simulated shift (blood readout weak, A35; descriptive).
- **Estimator.** Arrays: log2 intensity, equal gene weights (A53 convention). Core state (GTEx 14-gene waveform), cue state and locked state (A81b sets). Core-cue offset arg E\[c\_core c\_cue\*\] over subjects, per condition (the A115 estimator); within-subject change if biopsies are paired. If several time points per condition exist, add the timed version (cosinor per state).
- **H124.1.** The core-cue offset changes between control and misalignment, two-sided, condition labels swapped within subject (exact or 10,000), p < 0.05.
- **Specificity.** The core-locked offset changes less than the core-cue offset (difference of absolute changes > 0, same permutations).
- **Reading.** Pass: the readout detects imposed misalignment in human tissue, and A126-A127 cue results become claims. Fail: they stay descriptive.

### A125, mouse chronic jet lag (scout, then test)

- Scout: Hadadi et al. 2020 (chronic jet lag in a spontaneous mammary tumour model; tumour and mammary transcriptomes) and Kettner et al. 2016 (chronic jet lag and liver cancer; liver time series). Use only public matrices.
- If timed multi-tissue or timed tumour data exist: does the untimed method separate internal desynchrony (organs or compartments drifting apart, A111 logic) from misalignment with the light cycle (whole body shifted, A93 logic)? Pre-register the specific contrast once the design is known.

### A126, breast: the menstrual cycle as a natural experiment

Rationale: tumour incoherence grows with proliferation (A97); each luteal phase drives a wave of breast epithelial proliferation; lifetime number of cycles is a breast cancer risk factor. If physiological proliferation perturbs the epithelial clock like tumours do, the clock offers a mechanistic reading of that risk factor.

- **Data.**
  - [Murrow et al. 2022, GSE198732](https://pmc.ncbi.nlm.nih.gov/articles/PMC9590200): premenopausal breast scRNA-seq, per-donor hormone state, parity, obesity (primary).
  - Breast atlases with menopausal status, BMI and BRCA carriers: Kumar et al. 2023 (HBCA), Reed et al. 2024, Pal et al. 2021. Verify accessions and metadata at scout.
  - GTEx breast, women under 45 (timed; BMI held).
  - Endometrium across the cycle: GSE266265 (already built, menstrual phase recorded), Wang et al. 2020, Garcia-Alonso et al. 2021.
- **Scores (fixed now).** Luteal/progesterone score in hormone-responsive luminal cells: mean z of PGR, TNFSF11, WNT4, AREG; use the authors' hormone-state annotation when given. Proliferation: the A97 G2M set.
- **Proliferation control first (H126.0).** Progesterone regulates clock genes directly, so a menstrual-phase contrast mixes hormone with proliferation. The control is within tissue: in endometrial epithelium (single-cell series) and in breast luminal cells, A87 coherence of cycling minus non-cycling cells of the same donors (A24 top-quartile rule), donor-level label swaps, one-sided p < 0.05. If it fails, breast results are descriptive. Descriptive: endometrial coherence in proliferative minus secretory phase.
- **H126.1 (B1, primary).** Breast epithelial coherence (A87; luminal hormone-responsive, luminal progenitor, basal) is lower in donors in the top vs bottom tertile of the luteal score, pooled over Murrow and the atlases (study as stratum), one-sided permutation p < 0.05. Secondary: epithelial phase relative to the stromal host (A23 y4) by luteal tertile.
- **H126.2 (B2).** Core-cue offset in hormone-responsive cells varies with the luteal score (circular-linear slope), two-sided p < 0.05.
- **H126.3 (B3, obesity).** Epithelial coupling and coherence against BMI (A19 covariate model; Murrow obese vs not), two-sided.
- **Descriptive.** GTEx breast: alignment with time of death and coupling to other organs by luteal score. BRCA1/2 carriers vs non-carriers: coherence of normal epithelium (precancer). Parity (Murrow).
- **Feasibility.** >= 10 donors per tertile, or the test is reported infeasible.

### A127, breast tumours: scale and cue misalignment

- **SCAN-B (GSE96058; about 3,000 tumours with follow-up; verify access).** H127.1: A87 coherence lower in the top vs bottom proliferation tertile within PAM50 subtype and purity tertile (replicates A97 at scale, purity-adjusted as H114.0), one-sided. H127.2: ER-negative vs ER-positive coherence (replicates the A97 breast observation), one-sided. Exploratory: survival (Cox with age, stage, subtype; off-plane score as A97).
- **Single-cell tumours (Wu et al. 2021; Pal et al. 2021; verify accessions).** H127.3: within patients, core-cue offset of malignant cells minus that of stromal cells (fibroblast, endothelial) differs from 0, two-sided, patient-level label swaps; locked set as control. The psoriasis analogue: is the tumour's local clock misaligned from body-wide cues while its stroma is not?
- **Li 2024 breast (V25.5), descriptive.** Per-patient misalignment (internal minus wall time) against menopausal status, BMI and tumour proliferation, if recorded.

### A128, social jetlag: the weekly clock in GTEx (descriptive)

- **Data.** GTEx fast deaths with time of death; weekday from the held date of death (aggregates only). Ventilated donors as negative control (no weekly social schedule).
- **Estimand.** Per donor, host phase minus clock time (A93 H93.3 machinery); negative = internal clock delayed.
- **D128.1 (descriptive, not a hypothesis test).** Deaths on Monday-Tuesday vs Thursday-Friday, circular mean difference with its interval and the weekday-permutation p (2,000). The MDE (about 1 h with about 240 donors over seven days) exceeds any plausible tissue effect of social jetlag, so a null is uninformative and is reported with the MDE.
- **Negative control.** No weekday effect in ventilated donors.
- **Descriptive.** Interaction with BMI and age (younger people have more social jetlag); weekend deaths; the week after vs before daylight-saving transitions (low power).
- **Reading.** Descriptive only: a clear delay would be a candidate for a registered test in another cohort, not a claim.

## Order of work and decision rules

Run in this order; packages on the same line can run in parallel. Each outcome changes the figures as listed.

1. **WP0** (integrity; item 2 first resolves the `compute_cancer_split.py` null difference of M2) and **WP1** (A109, coupling nulls; A109b Tabula Sapiens).
2. **WP2** (A110, mechanism), only if H109.1 and H109.2 pass; **WP3** (A111-A112 feeding, A129 scout) and **WP4** (A113 living humans; A123 GTEx working-hours donors; A130 chronotype; A131 scout) in parallel.
3. **WP5** (A114, cancer) and **WP6** (A115-A116, psoriasis).
4. **WP7** (A117-A122): scouts first for A117, A119 and A120; A121 and A122 run on data already built.
5. **WP8** (A124-A128): A124 first, in parallel with WP5-WP6; A126-A127 cue results are claims only if A124 passes; A128 any time (GTEx data in hand).

| outcome | consequence for the paper |
| --- | --- |
| H109.1 and H109.2 pass (Holm) | Figure 3 as planned; WP2 results become its mechanism panel; A109b as replication panel if it passes |
| H109.1 passes, H109.2 fails (no offset preservation) | Figure 3 reframed as "organs share a donor state"; drop WP2 from the main text |
| H109.1 fails | Remove the coupling-without-cues claim; Figure 3 keeps loss of wall time and drug-target rhythms |
| Pass in ventilated donors, fail in Hardy 4 | Coupling reported as procurement-dependent; Figure 3 claim restricted accordingly |
| A111 passes | Feeding-control panel in Figure 1 (an addition, not a required positive); replaces the failed A49 control |
| A111 H111.3 or A129 passes | Genetic control panel (a tissue without its clock is detected as a break) in Figure 5 next to the taxonomy; otherwise a limitation sentence |
| A112 passes | Cue/locked split validated across species; shown in Figure 3 |
| A113 passes in >= 2 cohorts | Internal-time readout in living people is a headline translational result |
| H114.0 fails (A97 not seen within purity tertiles) | Bulk proliferation claim restated as purity-confounded; the proliferation panel rests on single-cell malignant cells |
| A114: MYC or copy number explains incoherence | Mechanism sentence in Figure 5; otherwise state that proliferation tracks it, cause open |
| A130 H130.1 passes | Internal time matches chronotype epidemiology outside GTEx: Figure 2 panel with the GTEx age and sex readouts; otherwise GTEx readouts in ED, labelled post hoc |
| A115 passes | Psoriasis delay partly a trait of patients' skin; changes the psoriasis claim in Figure 5 |
| A116 passes | Reversal specific to treatment; otherwise "follows clinical response, not shown to be specific" |
| A118: core damped in a break | That break described as damping; otherwise as a relative loss of coupling |
| A119 infeasible or failing; A121 failing | Lung T cells (A119) or the embryo panel (A121) move to Extended Data; the embryo panel stays in Figure 2 if A121 passes, even when A120 finds no replication data |
| A122 passes (locked-set and pseudo-cue controls) | Immune cue-following panel in Figure 2; otherwise Extended Data, described relative to host coupling |
| A123 passes | CHIRAL-style organ map from working-hours donors in Figure 1; otherwise the 8-tissue benchmark window result stays as the only window panel |
| A124 passes (imposed misalignment detected) | Misalignment readout validated: panel in Figure 1 next to the feeding control; A126-A127 cue results become claims. Otherwise they are descriptive |
| A126: endometrial control and H126.1 pass | Figure 5 panel: physiological proliferation perturbs the breast clock as tumours do; Discussion links it to lifetime menstrual cycles as a risk factor |
| A127 H127.1-H127.2 pass | SCAN-B replication of the proliferation and ER findings in Figure 5; H127.3 as a panel if A124 passes |
| A128 (descriptive) | A sentence with its interval and MDE in the Discussion; a clear delay is named as a candidate for a registered test elsewhere |

**Reporting.** After each package: amendment in `paper/docs/registered_experiments.md` (committed before scoring), run note and result in `paper/docs/results_preliminary.md` with `src` pointers, a row in `paper/docs/METHOD.md` when it validates the method, scorecard row, a line in `paper/docs/AMENDMENT_INDEX.md` once a script is promoted to `results/`, and a one-line message to the PI with the verdict. Subagent numbers are checked against the output files before they enter the docs.

## Narrative arc

**Positives.** The disease positives are the cancer and psoriasis breaks (Figure 5); the feeding and night-shift controls, if they pass, add perturbation evidence but are not required.

**One sentence:** human organs keep internal time with each other rather than with the clock on the wall; this coupling survives the loss of daily cues, chronic disease and age, and breaks only locally, in three distinct ways.

**Working title:** *The human body keeps internal time: organs stay coupled when the day is lost, and clocks break locally in disease.* Target: Cell, with a presubmission enquiry.

1. **The problem.** Human clocks are read from time-stamped samples, and almost no human tissue has them. Existing untimed methods order samples and fail when sampling is confined to working hours (benchmark).
2. **The idea, validated.** Read relations rather than times: how each organ or cell type follows the shared clock state of the same person. These relations are identified without times; they match known time across species and perturbations, and they survive working-hours sampling: the GTEx organ map is recovered from deaths between 08:00 and 18:00, and patients operated in the morning or afternoon are placed correctly (Figure 1).
3. **Organisation in health.** Each person carries an internal time that departs from clock time by hours. Organs and cell types follow it on a stromal-vascular backbone; the brain is partly separate, in humans and sheep; immune cells take relatively more of their timing from body-wide signals; strong clocks are the coupled ones; in the mouse embryo the shared clock appears late, mesenchyme first (Figure 2).
4. **When daily cues are lost, the body keeps its own time.** In ventilated and slowly dying patients organs lose wall time but stay coupled; if WP1-WP2 pass, they hold as a self-sustained peripheral network that free-runs away from wall time. Rhythmic genes and drug targets follow internal time, so clock-time dosing stops tracking them in hospital, while a biopsy can estimate internal time (Figure 3).
5. **Chronic conditions rarely move it.** Across about 20 disease contrasts, ageing to 70, sex and body mass, coupling holds within stated detection limits; same-time paired designs exclude hidden tissue-wide damping; the reported age dampening of rhythms reflects how donors died (Figure 4).
6. **Where clocks break, they break locally.** Malignant cells keep amplitude but lose coherence, growing with proliferation; psoriatic plaques run their core clock behind on-time systemic cues, a delay that shrinks with clinical response and is partly explained by keratinocyte differentiation; fibrotic-lung T cells lose relative amplitude, if the procurement check allows (Figure 5).
7. **Discussion.** Internal time as a measurable clinical variable; a tissue readout of circadian disruption to test night work, breast cancer and social jetlag where epidemiology places them (WP8); a peripheral network that holds without cues; limits (observational, GTEx frame, post-mortem states, ages above 70, donor numbers per cell type).

## Figures

Five main figures plus a graphical abstract. Each figure states one claim; its table lists the evidence, its source, how it supports the claim, and whether it goes in the main figure or in Extended Data (ED). **New** marks evidence that depends on a work package above; "if passes" means its place is decided by that test.

### Figure 1. Clock relations and internal time can be read from samples taken once, at unknown times, including in working hours

| content | source | how it supports the claim | main or ED |
| --- | --- | --- | --- |
| Concept: one sample per person, compartments, shared host state, coupling, relative amplitude, offset | schematic | Defines the readouts and why they need no times | main |
| Benchmark vs COFE, CYCLOPS, CHIRAL, TimeTeller: organ offsets r 0.934 / 0.900 (n 80), whole day and 08-18 h; donor-count sweep <!-- src: archive outputs/paper/fig1_method/gtex/sweep_summary.csv --> | benchmark (Fig. 3 code) | The only untimed method that recovers organ offsets, also in working hours | main |
| GTEx organ map (CHIRAL-style offsets, amplitudes, donor order) from 08-18 h deaths only | **new** A123 | The published organ map is recovered without any donor outside working hours | main if passes |
| Heart surgery patients, morning vs afternoon | Jia et al. 2024 heart samples | Living patients sampled in working hours are placed correctly on the internal day | main |
| Breast surgery: per-patient internal time vs recorded resection time (1.1 h) | V25.5 | Per-patient ordering in a living hospital cohort | main |
| Timed truths: baboon (offsets: median error 0.455 h; peripheral 0.383, neural 0.497) <!-- src: figure_1/output/species/a36_readouts.csv -->, regenerating mouse colon (9 h shift recovered blind), mouse at three ages, living skin layers | A36, A21, A48, A37 | Relations match known time across species, ages and a perturbation | main (baboon, colon); ED (rest) |
| Feeding control: daytime feeding desynchronises the liver, detected untimed | **new** A111 | The method sees a known desynchrony between organs (addition, not required) | main if passes |
| Living timed muscle and adipose cohorts | **new** A113 | Internal time estimated in living people beyond one cohort | main if passes, else ED |
| Working-hours windows in baboon; human vs baboon reference | A40, A41 | Robust to window width and to the species of the reference | ED |

### Figure 2. In health, organs and cell types follow a shared internal time built on a stromal-vascular backbone, with the brain partly separate

| content | source | how it supports the claim | main or ED |
| --- | --- | --- | --- |
| Donors' shared departure from time of death (hours); organs follow it with unequal gains | 01\_gtex\_organs | Each person carries an internal time distinct from clock time | main |
| Brain vs periphery r\* 0.64, robust to collection delay, death type and PER1/PER2 | A103 | The brain is partly separate, not by artefact | main |
| Sheep r\* 0.63; brain couples less in sheep and baboon | A42 | The brain-periphery split is conserved | main |
| Anchored coupling by cell family across 22 healthy atlases (17 studies, 14 tissues) <!-- src: legacy/desync/common_output/a17_unified_healthy.csv --> | A17 unified | Stromal-vascular backbone; parenchyma lower; lymphocytes barely coupled | main |
| Healthy mouse colon against known time (stroma high, T and B near 0) | A21 R21.1 | The ranking holds against a timed truth | main |
| Immune cells lean more on body-wide cues (index 0.18 vs 0.10, p 0.002) | A85 corrected; **new** A122 | Explains why immune cells barely follow the tissue clock | main if passes, else ED |
| Clock strength and coupling: pre-registered H104.1 (rho < 0) not supported; observed GTEx rho +0.73 (post hoc: stronger clocks are more coupled); H104.2 not replicated (sheep, baboon) | A104 | Post hoc reading, labelled as such | main |
| Mouse embryo: shared clock only late in gestation (E15-E18.75: H108.1 R 0.52, p 0.016, n 16; clock-shape gate borderline, p 0.08); earlier windows not clock-shaped; mesenchyme first, aligned to the mother's day | A108, A108b, A108c; **new** A121 (A120 if data exist) | The backbone is built late in gestation, consistent with known late emergence | main (exploratory label) if passes, else ED |
| Neurons carry the brain host, glia less; retina | 02, 04, A18 3.5 | The brain's backbone is neuronal | ED |
| Age: small phase advance in the periphery; prefrontal amplitude loss above 60 | A94, A3; **new** A117 if feasible | Ageing changes timing, not organisation, in the periphery | ED |
| Internal time follows chronotype epidemiology: earlier with age, women earlier | A94, A99 (post hoc); **new** A130 replication | A known-answer test of internal time in health | main if A130 passes, else ED |
| Gene-rhythm maps per cell type, lung and prefrontal cortex | A50 | Resource | ED |

### Figure 3. When daily cues are lost, organs lose wall time but stay coupled to each other

| content | source | how it supports the claim | main or ED |
| --- | --- | --- | --- |
| Fast, intermediate, slow and ventilated deaths: wall-time alignment falls, inter-organ coherence holds, against clock-specific nulls; organ phase relations preserved | A93; **new** A109 (H109.1, H109.2), A109b Tabula Sapiens | The core observation, with its controls and an independent cohort | main |
| A step, not a drift: no decay with ventilation days; no reset at intubation | A46, A92, A95 | Rules out a slow free-run and a common reset | main |
| Liver, muscle and adrenal loosen; the brain leaves the body | A33, A93 | Which organs depend on behavioural cues | main |
| Circularity and displacement simulations; decay with ventilation days; strong-clock, stroma and hub signatures; temperature, SCN and feeding excluded | **new** A109 H109.2, A110 | Mechanism: a self-sustained peripheral network | main if passes |
| Rhythmic genes follow internal time (liver 2,691 vs 0); cue/locked split validated in mouse | A91, A102; **new** A112 | Gene rhythms follow the body's clock, not the clock on the wall | main |
| Drug targets lose their wall-time rhythm in hospital; internal time estimable from a biopsy | A100, A102, V25.5; **new** A113 | Translational consequence | main |
| Glucocorticoid axis in ventilated donors | A27 | Part of the shared state, small part of the coupling | ED |

### Figure 4. Chronic disease, age, sex and body mass rarely change the clock network

| content | source | how it supports the claim | main or ED |
| --- | --- | --- | --- |
| Forest plot of all disease contrasts with MDEs (coupling d\*) | A17-A22 | Nulls are quantified, not just non-significant | main |
| Hierarchical pooling by disease class and cell family; neurodegeneration margins | A20, A22 | Resilience holds when cohorts are pooled | main |
| Paired same-time designs: no tissue-wide damping in tumours, inflamed gut, skin | V52, A52, A53 | Rules out losses invisible to a host-referenced readout | main |
| Ageing to 70, sex, BMI, T2D, ventilation dose | A46, A47, A94, A99 | The main covariates leave the network intact | main |
| CHIRAL's age dampening reproduced, then explained by death type | A101b, A102b; WP0 robustness | Corrects a published reading with the same data | main |
| COVID-19 decedents keep organ clocks and endothelial sharing | A105, A105b | Resilience under severe acute inflammation | main |
| Fatty liver: liver-wide loss carried by procurement; inflamed gut stroma intact | A43, A44, A96, A45 | Apparent breaks traced to design | ED |

### Figure 5. Where clocks break, they break locally, in distinct ways

| content | source | how it supports the claim | main or ED |
| --- | --- | --- | --- |
| Malignant cells: clock-shaped share near null while the stroma of the same tumours is intact | A86 | The break is in the cancer cells, not the tumour | main |
| Decomposition: amplitude kept, coherence lost, partly misaligned | A87 | Defines the cancer break as incoherence | main |
| Incoherence grows with proliferation (purity-adjusted); MYC, copy-number and proliferating-normal controls | A97; **new** A114 (H114.0-H114.3) | Links the break to tumour biology | main |
| A tissue without its own clock reads as a break | **new** A111 H111.3, A129 | Genetic control of the break readout | main if passes, else limitation |
| Lesional clocks run late: -1.50 h in discovery (H54.1-3 pass) and -0.93 h in held-out series (H54.2, p 0.004); held-out H54.1 (variance) and H54.3 (genes) not met | A54 | The delay replicates | main |
| Core clock delayed while systemic cues are on time (skin set and independent GTEx set) | A78, A81b | The delay is in the local clock, not in its inputs | main |
| The delay shrinks with treatment and clinical response; placebo pooling; trait vs state | A56, A57, A74; **new** A115, A116 | The delay follows disease activity | main |
| Keratinocyte clock advances with differentiation; the lesional state mix explains part of the delay | A64, A67, A72, A69 | A cellular explanation for part of the delay | main (one panel) |
| Specific to psoriasis (atopic dermatitis, other skin diseases, gut) | A59, A79, A82 | Not a generic inflammation effect | ED |
| Fibrotic-lung T cells lose relative amplitude, two cohorts | A98; **new** A118, A119 | A third mode, if procurement is resolved | main if A119 passes, else ED |
| Taxonomy: incoherent, damped, desynchronised, delayed | schematic | Summarises the modes of breakage | main |
| Per-cell synchrony vs coherence in malignant cells | WP0 item 2 | Reconciles the two cancer readouts | ED |

### Extended Data

1. Pre-registration scorecard (all amendments, verdicts, MDEs).
2. Estimator validation: A17 acceptance simulations, small-n calibration, anchors, A22 group confirmation.
3. A25 internal collection window, its failures (V25.3, A26) and the non-temporal post-mortem direction.
4. Gate power and diagnosis of failed fits (blood, lymphoid organs).
5. Waveform sensitivity, human vs baboon reference.
6. Glucocorticoid axis in ventilated donors.
7. Gene-level atlas go/no-go and donor-count curve (A27, A31, A39).
8. Per-cell amplitude and synchrony, with the WP0 reconciliation.
9. Development: human fetal and placental attempts (A106, A107, A107b) and embryo robustness (A108b leave-one-day-out, lineage halves of A108c).
10. Covariate sensitivity (A102) and per-disease contrast tables.
11. Circadian disruption proxies (WP8): simulated night shift, menstrual cycle, SCAN-B, weekday analysis, including negative results.
