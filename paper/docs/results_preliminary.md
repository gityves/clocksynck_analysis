# Preliminary results (append only)

Results of the amendments registered in `registered_experiments.md` (A109 onwards), each with its verdict and a
`<!-- src: ... -->` pointer to the output table it was read from. Results A1-A108c are in the archive
(`$SCCC_ARCHIVE/docs/RESULTS_preliminary.md`, frozen); corrections to them go in the next section, dated, with the
archive section they correct.

## Errata to the archive

### 2026-10-23, WP0 item 7: A22 V17.2 coherence values (corrects archive RESULTS "V17.2, A21 colon", line 1262)

The legacy V17.2 coherences were computed on a pseudobulk table whose empty `stratum` column was read as NaN, so pair
coherence was not centred within stratum (legacy e_control 0.16-0.32; delta -0.11, SE 0.07; case excess 0.044). The
port (`results/extended_data/scripts/simulations/compute_a22_colon.py`, colon as one centred stratum) gives e_control
0.45-0.55 over the four anchors (lymphatic EC 0.448, endothelial 0.517, smooth muscle 0.544, fibroblast 0.546) and
delta -0.17 to -0.40; lymphatic EC: delta -0.399 (SE 0.164), e_case 0.048. The verdict is unchanged: no anchor is
dropped, lymphatic EC is not dropped, V17.2 fails as in the archive. Other legacy coherence runs on tables with an empty
stratum may carry the same artefact (to be checked where cited).
<!-- src: paper/results/extended_data/output/simulations/a22_colon_confirm.csv (archive-port snapshot, sha256 e00e28cc...) -->

### 2026-10-23, WP0 item 2 (first part): `compute_cancer_split.py` D_null difference of M2 (no archive number changes)

The M2 difference in `D_null` was not a migration or clocksync change: the archive's stored subsample file was written
by an earlier version of the archive script that drew each permuted plane and its splits interleaved; the archive's
current script, its full run and the migrated script draw all planes first. A copy of the current archive script and the
migrated script give byte-identical outputs at `--subsample 12 --seed 0`; the archive full run matches the legacy
`x5_vm_human.py --cnv` groups (2,249 groups, D_null to 7e-13). Both orders are valid draws of the A18 null; D is
deterministic. The per-group `D_null` (mean over 10 permuted planes) is Monte-Carlo noisy (the two orders differ by up to
144 per group at equal mean), so per-cell null-corrected numbers are read pooled over groups. The reconciliation of
per-cell synchrony vs A86/A87 coherence (rest of item 2) follows.
<!-- src: paper/exploratory/WP0_integrity/scratch/cancer_split/{arch,mig}/a18_cancer_groups.csv (identical, cmp) -->

### 2026-10-23, WP0 item 1 (first part): fibrotic-lung T-cell Fisher combined p (supports archive A17 quote "0.002")

Fisher's combination of the two A17 learned-anchor permutation p values for lung T cells, HLCA fibrotic (0.0270) and
Natri ILD (0.0090), both two-sided (|null| >= |diff|): chi2 16.65, 4 df, combined p 0.0023. The archive's "0.002" is
this value rounded; the table now exists. The trace of the other figure numbers follows.
<!-- src: paper/exploratory/WP0_integrity/output/fisher_lung_tcell.csv (inputs: $SCCC_ARCHIVE/legacy/desync/common_output/a17_unified_contrast.csv, rows hlca_fibrotic and ild_natri x lung:T, column p_perm_learned) -->

### 2026-10-23, WP0 items 1 and 3: sources of the figure-table numbers; pre-registration scorecard (draft)

**Figure numbers (item 1).** All 19 numbers quoted in the plan's figure tables were traced: 8 supported by a migrated
output, 8 by a full-run archive table only (benchmark `outputs/paper/fig1_method/gtex/sweep_summary.csv`, A21 colon
`legacy/desync/21_colon_injury/output/r214_phase.csv`, A42 `r1_sheep_readout.csv`, A104
`legacy/desync/validation/oscillator/output/a104_gtex.csv`, A108 `outputs/paper/figure_3_embryo/a108_timed.csv`; these
enter M3), 3 design parameters; none unsupported. Value mismatches with the quoted text: benchmark offset r 0.934 /
0.900 (quoted 0.91 / 0.88); baboon offsets 0.455 h (quoted 0.5 h); healthy atlases 22 targets / 17 studies (quoted
"about 30"). Uses that contradict the archive verdict (for the PI; figure tables not changed): A104 rho +0.734 is quoted
as support but H104.1 predicted rho < 0 (one-sided p 1.00, "not supported"; H104.2 not replicated); A54 held-out
replicates the delay (H54.2, -0.93 h vs discovery -1.50 h) but H54.1 and H54.3 fail; mouse embryo "from E15": H108.1
passes in the late window (R 0.52, p 0.016) but the window gate does not (p_acrophase 0.08).
<!-- src: paper/exploratory/WP0_integrity/output/figure_number_sources.csv (columns status, value_in_file, archive_verdict, verdict_contradiction) -->

**Scorecard (item 3, draft).** `scripts/scorecard.py` parses 120 amendments (175 hypothesis rows) from both
repositories; a verdict is parsed for 79 amendments (structured fields 62 rows, keyword pass 56 rows, erratum 1: A85
now pass after the 2026-10-22 correction). 41 amendments need manual curation before the table becomes Extended Data
Table 1. Keyword verdicts are flagged in `parse_note`.
<!-- src: paper/exploratory/WP0_integrity/output/scorecard.csv -->

## A109 (scored 2026-10-23): is the coupling of organs without cues a clock?

Registered: `registered_experiments.md` A109 (commits 67d3962, scout note 9171b55). Full run, registered sizes.
Holm over the confirmatory family: H114.1, H115.1, H116.1 not yet scored; the bounds below hold for any outcome of
those three (Holm-adjusted p <= 5 x raw).

**H109.1 (confirmatory): PASS.** Ventilated donors (377, equal-n draws of 103): C_clock 0.463 vs pseudo-clock null median
0.154 (sd 0.010); excess 0.310 (donor bootstrap 95 % CI 0.225-0.348); one-sided p 0.005 (1/201, the floor at 200
sets); Holm-adjusted p <= 0.025. MDE was 0.018. Secondary: also above the acrophase-permuted null (median 0.282,
excess 0.181, p 0.005) and the shuffled-pairing floor (0.123). Fast, Hardy 3, Hardy 4 all exceed the pseudo-clock null
(excess 0.469, 0.246, 0.255; each p 0.005). Excess falls from fast to ventilated by 0.159 (CI 0.079-0.230) and from fast
to Hardy 4 by 0.214 (0.148-0.300); ventilated minus Hardy 4 +0.055 (-0.010 to 0.158). Processing check: on the 11
pairs with whole blood (different collection kit) ventilated C 0.302 vs null 0.153 (p 0.005); 55 same-kit pairs 0.499.
<!-- src: paper/exploratory/A109_coupling_nulls/output/h109_1_null.csv; excess_bootstrap.csv; procedure_check.csv -->

**H109.2 (confirmatory): FAIL in ventilated donors.** Fisher-Lee correlation of the 66 pair offsets, ventilated vs
fast: rho 0.261, one-sided p 0.165 (10,000 joint tissue-label permutations); Holm-adjusted p >= 0.165. Second
criterion also fails: the ventilated pair phases are closer to 0 / 12 h than in fast deaths (median folded distance
0.60 h vs 1.02 h, p_closer 0.001). Equal-n sensitivity: rho 0.342, p 0.078; p_closer 0.141. Power at registration: 80 %
at rho about 0.29.
**Procurement check (secondary, registered): Hardy 4 passes both criteria**: rho 0.784, p 0.0002; folded distance
1.06 h vs 1.02 h, p_closer 0.57 (equal-n rho 0.633, p 0.003). Hardy 3: rho 0.505, p 0.015, p_closer 0.82.
**Circularity (secondary):** pooled pseudo-covariance ratio in ventilated donors 0.456; simulations at the ventilated
design: rotating 0.141 (sd 0.017), 25 % displacement 0.280, 50 % 0.509, pure displacement 1.001; best-matching
displacement fraction 0.44. Descriptive per-pair medians of the ratio: fast 0.19, Hardy 4 0.33, Hardy 3 0.46,
ventilated 0.48.
**Non-temporal component (secondary):** the non-temporal shared axis learned in fast deaths (3.9 h) carries 0.70 of the
ventilated shared covariance (fast 0.56, Hardy 4 0.54); share along the glucocorticoid direction: ventilated 0.68, fast
0.54, Hardy 4 0.58. Host phases: Rayleigh R ventilated 0.014, axial R 0.26 (fast 0.08 / 0.03).
<!-- src: paper/exploratory/A109_coupling_nulls/output/h109_2_offsets.csv; circularity_sims.csv; pair_relations.csv (column pseudo_cov_ratio); nontemporal_share.csv -->

**Reading (registered rules).** H109.1 passes, H109.2 fails in ventilated donors: by the registered reading, organs of
ventilated donors share a donor state rather than running clocks with preserved phase relations; about half of their
shared covariance behaves as a displacement along one axis (circularity, axial R, non-temporal and GC shares). In
Hardy 4 deaths (post-mortem, no daily cues), the phase relations of fast deaths are preserved. The decision-table
outcome "H109.1 passes, H109.2 fails" applies to ventilated donors (Figure 3 reframed, WP2 out of the main text); the
Hardy 4 pass was registered as a secondary procurement check. Decision on the Figure 3 claim referred to the PI.

### 2026-10-23, WP0 item 6: gene-symbol audit (no erratum)

Every script and builder behind the cited results that reads an Ensembl-indexed input (atlas_wg, Census builders,
`_atlas/census.py`, scbc_full) matches on `feature_name` or a symbol column; only the legacy A85 script was affected,
already corrected on 2026-10-22. No result changes. Latent risk: `clocksync.genes.gene_columns` returns -1 silently for
a missing symbol (and -1 indexes the last column); a `missing=` option with an all-absent guard is proposed for
clocksync (diff and test not yet applied).
<!-- src: paper/exploratory/WP0_integrity/output/gene_audit/ (proposed_fix.diff); scripts/gene_audit/count_*.py -->

### 2026-10-23, WP0 item 4: estimator consistency of the figure-table couplings (for the PI)

Of 30 figure-table rows, 8 use the A17+ anchored estimator with A22 confirmation, 9 cite no coupling, 13 are flagged:
Fig 1 A21 colon (an anchor decouples, V17.2 fails); Fig 2 GTEx organ departures (pre-A17 timed), A103/A42 r* and A104
(block or all-tissue hosts, not the A17 rule), A85 index (leave-class-out mean host), neurons/retina panels (pre-A17
free fit; A18 3.5 is anchored); Fig 3 A92/A93 coherence C and the A33 tissue host (not the anchored c); Fig 4 A20
own-noise "damping" (withdrawn in A22) and COVID A105 (2 anchors, no A22 confirmation); Fig 5 fibrotic-lung T cells
(A5 coupling 0.45 was a free fit; whether `a17_unified_contrast.csv` reflects the A22 refit without the ild_natri
capillary EC anchor is unresolved). Several flags are by design (r*, C and block hosts are different readouts, not
couplings); the figure text must name the readout. Scorecard: 108/120 amendments with a verdict after curation.
<!-- src: paper/exploratory/WP0_integrity/output/estimator_consistency.csv (column flag); output/scorecard.csv -->

## A109b (scored 2026-10-23): Tabula Sapiens replication (secondary)

**Not replicated; underpowered.** C of stromal-vascular pseudobulks across 9 organs (12 donors, 36 pairs with >= 3
shared donors): C_clock 0.683 vs pseudo-clock null median 0.675, excess +0.009, one-sided p 0.438 (MDE 0.104).
Sensitivity (7 pairs with >= 5 shared donors): 0.619 vs 0.531, excess +0.087, p 0.204 (MDE 0.187). At 3-5 shared
donors the null floor of |correlation| is 0.53-0.68, so only a large excess could be detected. Offset test (4 organs
shared with GTEx, 6 pairs, 24 permutations): rho 0.091, p 0.417; registered as uninformative (power 0.06-0.13). The
archive's Tabula Sapiens multi-organ gate (p_acrophase 0.017, all compartments, 14 donors) is a different readout and
is not contradicted. Reading: Tabula Sapiens adds no evidence for or against H109.1 at this design; no replication
panel.
<!-- src: paper/exploratory/A109_coupling_nulls/output/a109b_c_null.csv; a109b_offsets.csv; a109b_power.csv -->

### 2026-10-23, WP0 item 5: CHIRAL reanalysis with TMM normalisation; sex comparison under death-type matching (no verdict changes)

CPM rows reproduce `figure_3/output/gtex_cues/a102b_a101b_tissues.csv` exactly; TMM (numpy implementation of edgeR
calcNormFactors, hand-tested; no edgeR cross-check available) on the same draws. Total log ratio of rhythmic-gene
counts, tissue bootstrap (not comparable to A101b's donor bootstrap):
- Age, old vs young (37 tissues; matched 29): CHIRAL subsample R1 CPM -0.80 (p 0.024), TMM -0.67 (p 0.013); pooled CPM
  -0.73 (0.008), TMM -0.58 (0.034); Hardy-matched CPM -0.13 (0.78), TMM -0.20 (0.65); Hardy x RIN-matched -0.02 / +0.005.
  Matched minus pooled: CPM +0.57 [0.11, 1.01] (p 0.024), TMM +0.39 [-0.10, 0.84] (p 0.092).
- Sex, female vs male (38 tissues, new): R1 CPM +0.79 (p < 0.001), TMM +0.61 (0.001); pooled +0.58 / +0.51;
  Hardy-matched +0.27 (0.15) / +0.23 (0.16); Hardy x RIN-matched +0.16 / +0.17. Matched minus pooled -0.32 (CPM, p < 0.001).
Reading: under TMM, CHIRAL's age dampening reproduces and vanishes with death-type matching, as with CPM (the
matched-minus-pooled difference is weaker, p 0.092; H101b was already a fail). CHIRAL's female excess of rhythmic genes
shrinks by one half to two thirds with death-type matching and is no longer significant.
<!-- src: paper/exploratory/WP0_integrity/output/chiral_tmm/count_comparison_tmm_summary.csv; sex_comparison_matched_summary.csv (columns estimate, lo, hi, p) -->

## A109c (scored 2026-10-23): Hardy 4 and wall time (H109c.1, H109c.2; secondary)

**H109c.1: no detectable wall-time alignment in Hardy 4, at limited power.** R 0.097 vs time-of-death permutation null
median 0.101 (sd 0.024), one-sided p 0.56 (103 donors, 1,000 permutations). Null q95 0.145; R for 80 % power 0.208.
By the registered reading the MDE condition (< 0.15) is met by the null q95 but not by the 80 %-power definition:
entrainment as strong as in ventilated donors (0.15) cannot be excluded, weaker entrainment than in fast deaths (0.53)
can.
**H109c.2: PASS.** Hardy 4 offset preservation survives removal of the wall-time component: residual rho 0.750, p 0.0002,
p_closer 0.39 (raw 0.784, 0.0002, 0.57). Without the 4 ventilator cases (99 donors): 0.745, p 0.0004, p_closer 0.37;
hospital inpatients only (85): 0.805, p 0.0001, p_closer 0.62.
Reading: the preserved phase relations of Hardy 4 are not inherited from (first-harmonic) wall-time entrainment, and do
not depend on the ventilator cases.
<!-- src: paper/exploratory/A109_coupling_nulls/output/a109c_walltime.csv; a109c_walltime_power.csv; a109c_offsets_residual.csv -->

## A109c D109c.3-4 (scored 2026-10-23, descriptive): the 0 / 12 h push in fatal COVID

**Delorey (4 organs, 23 donors, 8-14 shared donors per pair, C 0.485): uninformative.** Pooled circularity ratio 0.833
lies inside the simulated range of every model at this design (rotating median 0.64 [0.35, 1.08], 50 % mixture 0.81,
displacement 1.01); median folded pair distance 3.17 h (sims: rotating 3.1, mixture 2.2, displacement 1.4 h; GTEx fast
1.02, ventilated 0.60 h). Only heart-lung maps to an A92 GTEx tissue, so no Fisher-Lee correlation.
**Melms lung (8 cell classes, 19 COVID vs 7 surgical controls; C at the 0.99 simulation cap in both arms): no evidence
of a push in COVID.** Ratio COVID 0.677 vs controls 0.984 (permutation p 0.008); sims, rotating / 50 % / displacement:
COVID 0.454 / 0.650 / 1.004, controls 0.648 / 0.776 / 0.999. Folded distance 1.19 vs 1.65 h (p 0.56).
**Note on the circularity readout (written after these results).** For rotating clocks, E[z_a z_b] vanishes only if the
donors' host phases are uniform; any concentration of host phases (sampling window, synchronised schedule) raises the
ratio without a displacement. The surgical controls, sampled in operating hours, are expected near 1 for that reason,
and the ventilated GTEx ratio (0.456; axial R 0.26) can come from either a displacement or host phases concentrated on
one axis. The ratio is therefore not diagnostic on its own; offset preservation (H109.2, independent of the host-phase
distribution) is the decisive readout. Reading: the fatal COVID data do not show the GTEx-ventilated 0 / 12 h collapse,
but at these sizes they could not have shown it reliably.
<!-- src: paper/exploratory/A109_coupling_nulls/output/a109c_covid_delorey.csv; a109c_covid_melms.csv; a109c_covid_sims.csv (column dataset) -->

## A132 (2026-10-23): brain death or ventilation? (secondary)

**H132.1 not supported (no separation at this power).** Within Hardy 0, organ donors (brain dead) are not closer to
0 / 12 h than ventilated post-mortem cases: folded-distance difference +0.115 h at equal n (79 vs 79, 65 pairs; p 0.72;
MDE 0.54 h) and -0.057 h at full n (298 vs 79; p 0.68; MDE 0.44 h). Both groups are closer to 0 / 12 h than fast deaths
(median folded distance 0.75 and 0.77 h vs 1.01-1.02 h; p_closer 0.013 and 0.035; full-n organ donors 0.72 h, p_closer
0.012). Under the registered reading the push to 0 / 12 h goes with ventilation / ICU care rather than with brain
death, but the MDE exceeds the A109 ventilated-minus-fast gap, so a partial difference is not excluded.
<!-- src: paper/exploratory/A132_brain_death/output/a132_groups.csv; a132_contrast.csv -->

**Descriptive (registered as a per-group readout, no between-group test registered).** Offset preservation differs:
ventilated post-mortem cases keep the fast-death pair relations (Fisher-Lee rho 0.690, p 0.0006), organ donors do not
(rho 0.161, p 0.19 at equal n; 0.093, p 0.36 at full n). A test of this difference would be post hoc.
<!-- src: paper/exploratory/A132_brain_death/output/a132_groups.csv -->

**H132.2 not supported.** Top minus bottom ventilation-duration tertile (organ donors, cuts 61 / 101 h, 100 vs 100):
-0.192 h in the predicted direction, p 0.155 (MDE 0.45 h). Per tertile: folded 0.62 h (p_closer 0.009; rho -0.089) vs
0.79 h (p_closer 0.044; rho 0.247, p 0.16).
<!-- src: paper/exploratory/A132_brain_death/output/a132_groups.csv; a132_contrast.csv -->

D132.3 dropped (scout note). Scoring MDEs differ slightly from the scout values (different random draws; the scoring
draw is fixed by seed in `compute_brain_death.py`).

### 2026-10-24, WP0 item 2 (second part): per-cell D vs A86/A87 in malignant cells (inconclusive; simulation not calibrated)
Same 79 patients (A87 set). Real readouts, malignant / normal epithelium / fibroblast: A86 share 0.084 / 0.212 / 0.495
(permuted-acrophase null 0.058 / 0.044 / 0.064); A87 coherence 0.228 / 0.782 / 0.840, alignment 0.273 / 0.642 / 0.773
(d_coherence -0.554 [-0.725, -0.091], d_align -0.370 [-0.644, -0.076]); anchored coupling 0.284 / 0.325 / 0.461 (cohort
medians); per-cell D 4.01 (SE 0.68) / 7.11 (SE 8.05) / 17.1 (SE 4.9); R 0.462 / 0.462 / 0.487; malignant minus normal
d_R -0.178 (SE 0.20, 52 paired patients; A18: -0.249). Share and A87 components equal the stored A86/A87 outputs.
<!-- src: paper/exploratory/WP0_integrity/output/malignant_readouts.csv -->

Simulation (50 reps; malignant gene phases coherent vs uniform, amplitude unchanged): per-cell D 4.19 (SD 0.87) vs
-0.44 (SD 1.23) at the simulated depth; share 0.034 vs 0.030 and coherence 0.31 vs 0.27 (they separate only at x100 UMI:
share 0.35 vs 0.05, coherence 0.81 vs 0.28). **Not usable for the reconciliation yet:** the simulated normal epithelium
does not reproduce the real one (share 0.036 vs 0.212; coupling about 0; alignment degenerate), so the simulated depth or
between-patient variance is unrealistic, and "pseudobulk readouts sit at null at real depth" is a property of the
simulation, not shown for the data (the real normal epithelium is well above null at its depth).
<!-- src: paper/exploratory/WP0_integrity/output/incoherence_sim.csv -->

Lead's reading (no decision yet): real malignant D (4.0) sits at the coherent-scenario value and far from the
incoherent one, while A87 coherence is low. Either the A87 coherence gap is partly a depth / observability effect, or
per-cell D is insensitive in this regime; this bears on the A86/A87 main claim ("coherence lost"). Next (item 2, part 3):
calibrate the simulation to the real normal epithelium, and downsample real normal-epithelium counts to malignant depth
and recompute share, coherence and D. Cancer lambda, R and d_R stay out of the text meanwhile (d_R is not significant).

### 2026-10-24, A132b (POST HOC, PI request): offset preservation differs between the two ventilated groups
Over the 65 pairs with >= 20 shared donors in fast deaths and both groups (79 donors each), Fisher-Lee rho vs fast:
ventilated post-mortem 0.678, brain-dead organ donors 0.161; difference 0.517, two-sided donor-label permutation p 0.008
(null SD 0.204). Post hoc: chosen after A132's per-group readouts. Both groups are pushed toward 0 / 12 h (A132: folded
0.77 and 0.75 h vs 1.01-1.02 h in fast deaths). Description: ventilated donors not declared brain dead keep the
fast-death order of organ phase relations while these relations are pulled toward 0 / 12 h; brain-dead organ donors
show the pull without the order. A132 rho 0.690 for post-mortem cases was over 66 pairs (mask with one group).
<!-- src: paper/exploratory/A132_brain_death/output/a132b_rho_difference.csv; a132_groups.csv -->

**A132b descriptive (POST HOC, PI request): pull with the order kept, or reordering.** Signed deviation of each pair
phase from its nearest pole (0 or 12 h), group vs fast deaths (pairs with >= 20 shared donors; ratio = median
group / fast deviation over same-pole pairs with |fast deviation| >= 1 h; 1 = no pull, 0 = full pull):

| group | n | Spearman (p, tissue-label perm.) | ratio | same pole | Rayleigh R of time of death |
|---|---|---|---|---|---|
| ventilated post-mortem | 79 | 0.649 (0.001) | 0.58 | 1.00 | 0.067 |
| organ donors (equal n) | 79 | 0.257 (0.066) | 0.40 | 0.97 | 0.138 |
| organ donors (all) | 298 | 0.095 (0.37) | -0.01 | 1.00 | 0.110 |
| Hardy 4 | 103 | 0.788 (0.001) | 0.77 | 1.00 | 0.072 |
| fast (reference) | 244 | | | | 0.099 |

Reading: in ventilated post-mortem cases the deviations keep their order and shrink to about 0.6 of the fast-death
values, as expected for running clocks plus a shared non-temporal displacement; in organ donors the order is mostly
lost. Same pole is uninformative (most fast-death pairs lie within 6 h of 0). Times of death cover the day: R is near
its value under uniform times (about 0.10 at n = 79, 0.05 at n = 298); a concentration of host phases alone would not
pull offsets toward 0 / 12 h (it rescales E[z_a conj z_b] without changing its phase).
<!-- src: paper/exploratory/A132_brain_death/output/a132b_shrinkage.csv -->

### 2026-10-24, WP0 item 8: anchor robustness (block agreement r*, all-LOO host contrasts)
265 target x contrast rows (Figures 4-5), A17 learned anchors (pre-A22 snapshot); the anchored contrasts recompute the
A17 table exactly. All-LOO host: every compartment, the target left out (`Dataset(anchors=all)`). 200 permutations
(registered 1,000), 100 bootstraps. r* undefined with < 4 anchors (raw r reported as a lower bound). Helpers in
`paper/common/anchor_robustness.py` (tests pass); no clocksync change needed.
- **Primary targets: no sign change in the headline contrasts.** Anchored -> all-LOO: ILD (Natri) lung T -0.452 ->
  -0.458 (r* 0.92); HLCA fibrotic lung T -0.140 -> -0.173 (r* 1.00); T2D islet beta -0.397 -> -0.318 (r* 1.00); A44
  liver hepatocyte -0.290 -> -0.272 (r 0.86); AD skin T -0.302 -> -0.294 (r 0.75); A45 gut T -0.062 -> -0.040; DCM heart
  fibroblast -0.035 -> -0.017; CKD kidney T 0.766 -> 0.787.
- **Sign changes: 60 of 265, mostly near zero.** Endometriosis stromal +0.024 -> -0.073 (SE 0.15; 3 anchors).
  Brain (Figure 4): r* 0.42-0.67 and 1-7 of 9 compartments change sign per contrast (PD_mssm 7, ASCVD_hbcc 6, VaD_mssm
  5; `pd` 15 of 34, r* 0.53): **the brain host is anchor-dependent; brain disease contrasts are not anchor-robust.**
  Cancer per-cohort tumour-minus-normal epithelial coupling (9 of 21 cohorts estimable): 3 sign changes (ERP149257,
  GSE181294 near zero; GSE270680 +0.062 -> -0.538, r 0.20); the per-cohort anchored host is weak in several cohorts
  (r 0.05-0.49).
- Not covered: skin_ssc skin T and endometriosis lymphocyte (boundary), A87 malignant (no donor contrast), Figure 4
  age-slope forest (A46/A47), A22-confirmed anchors. Readout per row: `sign_change`, `r`, `r_star`.
<!-- src: paper/exploratory/WP0_integrity/output/anchor_robustness/anchor_robustness.csv -->

### 2026-10-24, A133 V133.1: simulation gate FAILS; A133 stops, the method is unchanged
Circular-part estimator (rank-1 Takagi fit of the pseudo-covariance subtracted from the Hermitian cross-covariance),
100 reps per cell, A109 generative model (fast-death offsets, ventilated amplitudes and noise). Gate (n 79):
(i) rho with truth at f 0.5 = 0.646 (< 0.8), and raw is not worse for a common-direction displacement (raw 0.663);
(ii) agreement with raw at f 0 = 0.941 (< 0.95); (iii) concentrated host phases, f 0: 0.782-0.808 (< 0.9; raw 0.856+).
At n 79 the raw estimator itself reaches only 0.88 with no displacement; the axial fit adds variance (about -0.05).
The decomposition helps only for per-tissue displacement directions (f 0.44: raw 0.36 vs circular 0.67 at n 79; 0.37 vs
0.90 at n 298). Registered consequence: no circular-part estimator in clocksync, H133.1 not run.
<!-- src: paper/exploratory/A133_axial_decomposition/output/v133_gate.csv; v133_grid.csv -->

**Post hoc reading of A132/A132b against the same simulations (raw estimator, R 0).** A shared displacement along one
common direction keeps the order of the offsets and pulls them toward 0 / 12 h: at f 0.44, n 79, rho 0.70 (SD 0.12) and
folded distance 0.71 h; observed in ventilated post-mortem cases (n 79): rho 0.68, 0.77 h, i.e. running clocks with the
fast-death relations plus a shared displacement of about this size. No cell reproduces brain-dead organ donors: at
n 298 a folded distance near the observed 0.72 h (f about 0.2-0.44: 0.78-0.58 h) goes with rho 0.90-0.95, and even f
0.8 keeps rho 0.55 (folded 0.32 h); observed rho 0.095. Their loss of order is therefore not a larger displacement of
the same kind: the relations between organ clocks themselves differ (or a displacement with tissue-specific directions
plus a common pull, not simulated jointly).
<!-- src: paper/exploratory/A133_axial_decomposition/output/v133_grid.csv; paper/exploratory/A132_brain_death/output/a132_groups.csv; a132b_rho_difference.csv -->

### 2026-10-24, WP0 item 2 (third part): the malignant coherence gap is not depth; per-cell D sees no loss
Same 79 patients. **Depth:** malignant pseudobulks are deeper than normal epithelium, not shallower (clock UMI per cell
4.31 vs 2.29; cells per patient 732 vs 751), so depth cannot explain their lower readouts. Thinning the normal epithelium
(20 seeds) lowered it further (391 cells, 2.10 UMI per cell; about a quarter of the malignant pseudobulk counts): share
0.102 (SD 0.028; full depth 0.212, malignant 0.084), coherence 0.680 (0.082; malignant 0.228), alignment 0.676 (0.029;
malignant 0.273). The thinning was not depth-matched (the agent's reading "depth explains the share gap" is withdrawn:
the thinned normal is far shallower than the malignant compartment). Coherence and alignment survive a 4-fold
shallower normal pseudobulk: **the A87 coherence and alignment gaps are not a depth effect.**
**Per-cell D:** malignant 3.89, normal 3.90 on the same patients (thinned normal 4.54, SD 2.75).
<!-- src: paper/exploratory/WP0_integrity/output/malignant_downsample.csv -->

Simulation calibrated to the normal epithelium (B 0.7, H 2, sd_off 4 h; share 0.208, coherence 0.772, alignment 0.617,
D 4.0 within 1 SE of the real values; coupling 0.114 vs 0.325 and lambda 2.07 vs 3.18 not reached). Malignant
compartment at real depth, coherent / uniformly random gene phases / 3 h gene-phase jitter: share 0.40 / 0.054 / 0.21;
coherence 0.84 / 0.30 / 0.58; D 4.70 / 0.42 / 2.85. Real malignant: share 0.084, coherence 0.228, D 3.89.
<!-- src: paper/exploratory/WP0_integrity/output/sim_calibration.csv; incoherence_sim_calibrated.csv -->

**Reading (lead).** The pseudobulk readouts sit at the scrambled-gene-phase scenario, per-cell D at the coherent one;
no simulated scenario fits both. Scrambled gene phases within malignant cells are rejected by D. A87 itself shows the
other component: off-plane variance of the clock genes across tumours is about 3-fold higher (d_log_off +1.16
[0.89, 1.79]), clock-plane variance is not lower (d_log_plane +0.30 [-0.27, 0.94]). Candidate reconciliation, not yet
tested: malignant cells keep a coherent clock within each cell, while across tumours clock-gene expression carries large
non-clock variation (copy number, oncogenic programmes) that hides the shared phase structure. This would restate the
A86/A87 claim from "gene coherence lost" to "clock intact in cells, between-tumour clock-gene variation dominated by
non-clock sources". Question to the PI (Q8). Cancer lambda, R and d_R stay out of the text.

## A110 (2026-10-24): what keeps organs coupled in Hardy 4? (H110.1 infeasible, not scored)
- **H110.2a PASS (Holm over a-b).** Over 23 peripheral tissues, Hardy 4: Spearman(own amplitude, coupling) 0.498,
  one-sided p 0.007, Holm 0.014 (fast deaths 0.745). Strong clocks still carry the coupling.
- **H110.2b not supported, opposite sign.** Retained coupling (Hardy 4 / fast at d*) vs Eraslan stromal-vascular
  fraction over 7 tissues: rho -0.857, p 0.994 (smooth-muscle-inclusive fraction -0.821). Seven points; not read
  further.
- **H110.2c (descriptive).** Hubs (Hardy 4 half 1): oesophagus GE junction, visceral adipose, lung, oesophagus
  muscularis; mean drop on removal 0.010 (Hardy 4) vs 0.003 (fast), difference 0.007, p 0.32 against 100 random sets.
- **H110.3: temperature excluded.** C 0.413 -> 0.395 after regressing out the cold and heat scores; drop 0.017
  [-0.000, 0.032] (< 0.05), random 6-gene scores median 0.016 (q95 0.038), p 0.44.
- **H110.4 (descriptive fallback; the organ-donor contrast is infeasible).** Hardy 3-4 (81) vs fast (145), coupling to the
  peripheral host at d*: hypothalamus +0.087, median other brain region -0.121; difference 0.208, permutation p 0.002
  (MDE 0.18; on raw d 0.180, p 0.012). Other brain regions lose coupling to the body in Hardy 3-4, the hypothalamus
  does not. Not an SCN test (dissections need not contain the SCN; illness and drugs confounded).
- **H110.5: feeding not detected, weak test.** Within-digestive minus digestive-to-other coherence, Hardy 4 minus fast:
  -0.0005, p 0.49 (MDE 0.20), on 2 within-digestive pairs (43 cross) at this fast draw. Excluded by the registered
  rule; with 2 pairs this exclusion carries little weight.
<!-- src: paper/exploratory/A110_mechanism/output/a110_h2.csv; a110_h2b_tissues.csv; a110_h2c_hubs.csv; a110_h3.csv; a110_h4.csv; a110_h5.csv -->

## A113 / A130 (2026-10-24): internal time in living humans
**H113.1: adipose passes; muscle passes only with an estimator chosen after the registered default failed (post hoc).**
The registration says "circular correlation"; the A35 default (Jammalamadaka-SenGupta) is documented in clocksync as
ill-defined when one sample is evenly spread (`clocksync/src/clocksync/stats.py`), which is the muscle design (4-hourly
grid over 24 h). Fisher-Lee was adopted after the JS value was seen, so the muscle pass is labelled post hoc.
- Adipose (GSE104674, 48 biopsies / 12 men, vs ZT): JS 0.888, Fisher-Lee 0.822, within-subject-permutation p 0.001 for
  both; median error after the common offset 1.11 h.
- Muscle (GSE108539, 57 / 10, 13-gene refit of the GTEx waveform): JS -0.854 (p 1.00); Fisher-Lee 0.576 (p 0.001);
  resultant 0.775; median error 1.30 h; common offset -0.63 h; per-subject Fisher-Lee median 0.80, positive in 10 of 10.
  Sensitivity (ARNTL row dropped from the published waveform): same conclusions (FL 0.585, error 1.25 h).
- Hospital-window subsets (muscle): Spearman of the tangent coordinate with clock time 0.885 (08-18 h, n 29, 3 times)
  and 0.798 (09-17 h, n 20, 2 times), direction correct; estimated window width unreliable (14.7 h vs 9.6 h).
- Reference: V25.5 error 1.1 h, offset -2.2 h; GTEx single tissue 2.6 h.
<!-- src: paper/exploratory/A113_living_time/output/a113_summary.csv; a113_windows.csv -->

**H130.1 not supported (p 0.054).** Internal minus wall time against age, cohort as stratum, 54 subjects (V25.5 42 women,
GSE104674 12 men; 33-93 y): -0.339 h per decade, one-sided permutation p 0.054 (null SD 0.20 h per decade). Per cohort
(descriptive): V25.5 -0.31 h per decade (Spearman -0.20), adipose -0.67 (-0.28). Direction as the epidemiology and as
A94 (0.2-0.5 h per decade), not significant.
<!-- src: paper/exploratory/A113_living_time/output/a130_summary.csv -->

## WP3 (2026-10-24): feeding and genetic controls (A111, A112, A129)
- **Informativeness.** Timed liver shift relative to the other tissues: Manella Alb-Cre ad libitum -> daytime feeding
  +5.46 h [4.78, 6.16] (liver +11.7 h, WAT +10.6 h, lung +1.9 h: the "rest" is heterogeneous); Xin night -> day
  restricted +10.11 h [9.21, 10.88]. Both informative.
- **H111.1 (untimed offset change tracks the timed change).** Xin PASS: Fisher-Lee 0.988, exact p 0.042 (the minimum
  with 4 tissues). Manella FAIL by design: 0.986 but exact p 0.167 is the minimum with 3 tissues (a registration flaw
  noticed only at scoring; the plan's 6 tissues were partly qPCR).
- **H111.2 (liver coupling falls under daytime feeding) FAIL.** Manella liver d -0.556 (p 0.010), d* -0.370 (p 0.060),
  but the timed liver-vs-rest R^2 change is +0.164 [-0.093, 0.347] (registered direction condition not met) and WAT
  and lung drop as much (d -0.753, -0.736). Xin d -0.172 (p 0.21), d* -0.011.
- **H111.3 / A129 PASS, specific: a liver without its clock decouples from the body.** Hepatocyte Bmal1 KO vs Alb-Cre,
  ad libitum, 24 vs 24 mice: liver coupling 0.967 -> 0.676, d* -0.287 [-0.411, -0.118], p 0.001; companions on the same
  permutations WAT d* -0.009 [-0.089, 0.055], lung +0.019 [-0.045, 0.087]; liver minus WAT -0.278 (p 0.001), minus
  lung -0.306 (p 0.001). Same with CIART kept (d* -0.325, p 0.001). The liver gate still passes in the knockout (the
  alternative criterion is not met). Daytime-feeding arm descriptive: d* -0.138, p 0.20.
- **H112.1 FAIL (wrong sign).** Greenwell arrhythmic vs ad libitum, untimed lock ratio, cue minus locked (24 vs 26
  mapped genes): +0.230 [-0.096, 0.819], p (less) 0.91. Truth check dropped (supplement not accessible).
- **Clock-loss readouts (descriptive).** Knockout liver offset to the WAT + lung host -1.31 h untimed (timed -0.98 h),
  p 0.001; WAT-lung coherence unchanged without a liver clock (0.958 vs 0.960). SCN-specific Bmal1 KO (VIP-Cre, CT0 /
  CT12, 6 vs 6 mice, ad libitum): liver-heart coherence 0.947 vs 0.912 (dC -0.035 [-0.178, 0.107]), no coupling loss
  (liver d* +0.10); small sample and 2 time points.
<!-- src: paper/exploratory/A111_feeding_controls/output/a111_informativeness.csv; a111_h111_1.csv; a111_h111_2.csv; a111_h111_3.csv; a112_h112_1.csv; a111_clockloss_liver_offset.csv; a111_pair_coherence.csv; a111_clockloss_coupling.csv -->
Deviations: genes failing QC (CIART, and CLOCK absent from Manella lung knockout files) were dropped from all tissues of
a dataset (`Dataset` takes one gene set); two-tissue companion contrasts use a local engine on `clocksync._core.host`,
verified identical to the public contrast for the liver. Builder docstring of `greenwell2019` cites the wrong journal
(Cell Reports 2019, PMID 30995463).

## A135 (2026-10-24): brain reference compartment, control donors only: no single reference
Five independent control sets (PsychAD MSSM 254 donors, HBCC 92, RADC 26; SEA-AD DLPFC 41; Parkinson 22). Rank sum of
clock strength, gate share and centrality per compartment; first-ranked family per set: excitatory neurons in 2 of 5
(MSSM, HBCC; gate 1.0), astrocytes 1 (RADC; gate share 0.43), inhibitory neurons 1 (SEA-AD), vascular 1 (Parkinson).
No family reaches the registered 2/3 (4 of 5). Post hoc, neurons merged: 3 of 5, still below. The subsample check gave
the same decision with different winners: first rank is not stable. Registered consequence: no brain reference is
justified; brain disease panels go to Extended Data with their anchor dependence (WP0 item 8); part 2 not run.
<!-- src: paper/exploratory/A135_brain_reference/output/a135_decision.csv; a135_ranks.csv; a135_criteria.csv -->
Choices beyond the registration: >= 15 control donors per compartment and >= 15 shared donors per pair (20 left the
Parkinson gate without pairs); nested control sets collapsed to the largest per cohort.

## A134 (2026-10-24): malignant coherence after removing inferred copy number and oncogenic programmes
- **H134.1 passes by p, small effect.** Regressing each clock gene on its own inferred copy-number score raises malignant
  coherence 0.228 -> 0.291 (gain +0.063) and normal epithelium 0.782 -> 0.768 (-0.014); malignant minus normal gain
  +0.077, one-sided permutation p 0.017 (null SD 0.039, MDE 0.097); patient-bootstrap interval [-0.160, 0.266] includes
  0. Copy number accounts for about a tenth of the 0.55 gap. Per-gene correlation with own-locus score 0.09-0.50
  (CRY2 0.50, NR1D2 0.44).
- **H134.2 not supported.** MYC targets, hypoxia and cycling jointly: -0.092, p 0.96 (each alone: MYC -0.055, hypoxia
  -0.088, cycling +0.003).
- **D134.3 uninformative.** Per-cell D with patient-specific vs common gene phases (64 patients): intervals span about
  +-10 (malignant minus normal cross-fit gain -0.56 [-14.8, 13.2]).
- Reading: copy-number variation of the clock loci explains a small, detectable part of the lower malignant coherence;
  oncogenic programmes do not; most of the gap remains unexplained. The registered restatement ("tumour-to-tumour
  variation breaks the shared phase structure") is supported only for that small part (PI, Q8 follow-up).
<!-- src: paper/exploratory/A134_cancer_coherence/output/a134_h1_h2.csv; a134_h1_h2_gene_effects.csv; a134_d134_3.csv -->
