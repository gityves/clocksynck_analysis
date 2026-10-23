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
