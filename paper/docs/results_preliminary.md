# Preliminary results (append only)

Results of the amendments registered in `registered_experiments.md` (A109 onwards), each with its verdict and a
`<!-- src: ... -->` pointer to the output table it was read from. Results A1-A108c are in the archive
(`$SCCC_ARCHIVE/docs/RESULTS_preliminary.md`, frozen); corrections to them go in the next section, dated, with the
archive section they correct.

## Errata to the archive

### 2026-10-26, PI Q10: A19 part 2 re-read with a within-stratum age permutation null (corrects archive RESULTS "A19 Part 2, Secondary")

The archive validated `Dataset.covariate` with a simulated null. At the real count tables it is miscalibrated: null
rejection reached 0.220, with 14 of 24 cells above the registered 0.07 (PI_QUESTIONS Q10). With the permutation null
(clocksync 613e224, the default), null rejection is 0.015-0.080, with one cell above 0.07 (expected about 1.7 by Monte
Carlo error at 200 draws). A19 part 2 re-read on the real ages:
- **Primaries unchanged.** H19.1 gamma -0.035 per decade (SE 0.072, p 0.630; archive -0.03, 0.07, 0.64). H19.2
  +0.004 (SE 0.040, p 0.922; archive -0.004, 0.036, 0.92). Both are not supported.
- **Secondary, neuron own noise.**
  - Pooled neuron gamma -0.116 (SE 0.054, p 0.031; 5 studies, 22 compartments; archive -0.12, SE 0.05, p 0.013,
    4 studies). The port has no design check, so six compartments the legacy check excluded enter.
  - The archive's "each per-compartment p <= 0.01" no longer holds. Of the six PsychAD compartments it names, MSSM
    excitatory (0.010) and MSSM CGE (0.005) remain. HBCC excitatory (0.020), HBCC MGE (0.015), MSSM MGE (0.025) and HBCC
    other (0.090) do not. Bias-corrected gamma keeps its sign and size (-0.16 to -0.25).
  - Read: the decrease of neuron own noise with age is nominal (p 0.031 uncorrected), a secondary, in the same direction
    in both PsychAD cohorts.
- **Calibration caveat.** Planted locked-amplitude effects are not recovered in the lung atlas, where fits sit at the
  coupling boundary (`calibration/README.md`). Lung age readouts stay descriptive.
<!-- src: paper/results/figure_4/output/resilience/a19_age_{pooled,compartments}.csv; paper/results/extended_data/output/calibration/age_acceptance_summary.csv (permutation null, 2026-10-26) -->

### 2026-10-26, M3: A20 R20.4 and A23 R23.1 intervals at the boundary (corrects archive RESULTS "R20.4" and "R23.1")

The full M3 rerun of `figure_4/scripts/spread_multiplicity` reproduces the A20 / A23 inputs exactly (y1-y3 of 230 rows to
5e-11; y4 of 194 rows) but refits the host in every bootstrap replicate (SEs 4-5 % larger for y1-y3, 11 % for y4). Three
archive statements rest on intervals that now reach 0:
- R20.4 "neurons are lower in every brain disease class": neurodegeneration x neuron d* -0.069 (-0.138 to +0.001) and
  metabolic x neuron -0.087 (-0.186 to +0.003) include 0; psychiatric x neuron -0.104 (-0.204 to -0.018) still excludes it.
  Point estimates within 0.004 of the archive. Read: neurons lower in psychiatric contrasts; in neurodegeneration and
  metabolic disease the same direction at the boundary.
- R23.1: neurodegeneration myeloid -0.57 h (-1.12 to +0.001) and glia +0.40 h (-0.002 to 0.88) now touch 0 (archive
  -0.64 and +0.43 h with intervals excluding 0); metabolic x neuron +1.08 h (+0.004 to 2.63) is the only interval
  excluding 0.
No pre-registered primary changes (R20.1-R20.3 resilience calls unchanged; R23.2 epithelium unchanged).
<!-- src: paper/results/figure_4/output/spread_multiplicity/{a20_cells,a20_readouts,a23_cells}.csv (M3, 2026-10-26) -->

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

## A123 (2026-10-24): the GTEx organ map from working-hours deaths only: H123.1-3 pass
Fast deaths (Hardy 1-2) who died 08-18 h: 111 donors, 40 tissues with >= 20 window donors; untimed estimates against
truth (i) (cosinor on recorded time of all fast deaths, complex SVD) and truth (ii) (re-derived from CHIRAL's per-sample
times); 10,000 tissue permutations, 200 bootstraps.
- **H123.1 PASS:** tissue offsets vs truth (i): circular correlation 0.945 [0.907, 0.972], permutation p < 0.001, median
  error 0.35 h (truth (ii): 0.831 [0.683, 0.924], 39 tissues). CHIRAL run on the same window donors: -0.165 (p 0.84).
- **H123.2 PASS:** relative amplitudes vs truth (i): Pearson (log2) 0.976 [0.956, 0.987] (truth (ii) 0.938; CHIRAL
  0.928).
- **H123.3 PASS:** donor ordering within the window, tangent coordinate with the hospital prior vs recorded time:
  Spearman 0.613 [0.476, 0.723], 111 donors, direction correct (the curvature-resolved flip without prior gives -0.613);
  estimated window width 15.3 h vs true 10.2 h (overestimated, as in A113). CHIRAL consensus: circular correlation
  0.30, median error 3.7 h.
- Sensitivity, the all-donor GTEx waveform instead of the cross-fitted one: unchanged (offsets 0.945, amplitudes 0.970).
Reading: the CHIRAL-style organ map (offsets, amplitudes, donor order) is recovered from working-hours deaths alone,
where CHIRAL itself fails on offsets. Figure 1 panel qualifies ("main if passes").
<!-- src: paper/exploratory/A123_working_hours/output/tests.csv; donor_time.csv; organ_map.csv -->

### 2026-10-25, M3 windows: A25 at full size; V25.1 acceptance now passes (archive: fails at 6 h windows)

First full run of the migrated A25 (`results/figure_1/scripts/windows/compute_internal_time.py`; the archive port ran
V25.1 and R25 on subsamples only).
- **V25.1 (simulation) PASS** at the registered acceptance (n = 40, coupling 0.8, w <= 12 h): axis error 0.34 / 0.27 /
  0.40 h and flip correct 0.91 / 1.00 / 1.00 at w = 6 / 9 / 12 h (archive: 0.53 / 0.38 / 0.40 h, flip 0.69 / 0.85 /
  0.99, failing at 6 h). Not a like-for-like replication: the port simulates counts with `clocksync.simulate` (archive:
  simulated states) and rebuilds Sigma_f from posterior host means plus covariance; the archive's estimated widths were
  biased upwards (w_est 15-16 h at true 6-12 h, coupling 0.8, n = 40; here 5.5-11.8 h). Which change removes the bias
  is not isolated. Quote V25.1 from this run with that note, or not at all.
- **V25.2 PASS** as in the archive: median |axis error| 1.26 h, flip correct 0.78 over 18 windows with r < 0.6.
- **R25 T1 PASS, T2 PASS, T3 not supported**, as in the archive. Corrected 2026-10-26: the first M3 read used the
  copy-number-split cancer tables and lost 8 cohorts; rerun with legacy's main tables. All 50 cohorts fit, with 17
  interpretable living cohorts (archive 17).
  - T1: mean axis 11.33 h, R 0.804, V test p 0.0014 (archive 11.01 h, p 0.010).
  - T2: 14 of 17 in 07-19 h, p 0.006 (archive 15 of 17, p 0.001).
  - T3: r living 0.275 (32 cohorts) vs deceased 0.332 (11), p 0.176 (archive p 0.214).
  - Rerun on the relearned anchors (2026-10-29): T1 R 0.816, p 0.0011 (PASS); T2 unchanged; T3 deceased r 0.469,
    p 0.068 (still not supported; the move comes from one cohort, `ia_spleen` deceased). See "M3 rerun of the
    anchor-dependent groups" below.
  - Anchor sets equal legacy in all 50 cohorts. Reported n counts every donor of the table; legacy counted donors with at
    least one anchor compartment (25 of 50 differ). Windows README, Deviations.
- Window check (GTEx 08-18 h vs others) at 500 bootstraps: observed / predicted 1.37 [1.05, 1.80] (archive [1.05, 1.80]).
<!-- src: paper/results/figure_1/output/windows/internal_time_{sim_summary,tests,cohorts}.csv (cohorts: rerun 2026-10-26 22:24); window_check.csv; logs/m3_windows.log (V25.2 line) -->

## A140 Q2 (i), (iii) (scored 2026-10-26): GSE266919 immune compartments do not pass the gate; P140.1-3 not estimable

Scored as fixed in the dated note (`registered_experiments.md`, A140 Q2 readouts; commit ff64bdd). Five compartments
(B, CD4_T, CD8_T, NK, myeloid), pseudobulks kept at >= 20 cells. Primary set: 43 of the 44 Pre-treatment samples (one has
no compartment with >= 20 cells). Sensitivity set: 77 of 78 samples.
- **Gate fails in both sets.**
  - The shared covariance between compartments of a biopsy is real (pairing p 0.001 in both).
  - It is not clock-shaped (acrophase p 0.459 primary, 0.540 all-78). The clock share is 0.025 against a null mean of
    0.028 (all-78: 0.025 vs 0.031).
  - Per CLAUDE.md, the failed gate is uninformative, not "decoupled": couplings are withheld, and P140.1-P140.3
    (composition and antigen presentation along the host phase) are not estimable.
- **Anchors.** The learned rule kept all five compartments (status ok); the fixed rule gives none. Pooled family scores
  at tissue scope were myeloid 0.199 (lower 0.113) and lymphocyte 0.123 (lower 0.040). The note quoted the organ_breast
  all-tissue table (0.202 / 0.091), so the scores differ, but the anchor outcome is the same.
- **Host.** Stability 1.000 (median similarity 0.980). The host is well identified but is not a clock: the dominant
  shared state of the immune compartments within a biopsy is non-clock.
- **Orientation check** (simulated timed cohort, 5 seeds): a larger host phase is later. R(phase - t) is 0.984-0.989,
  against 0.025-0.188 for R(phase + t).
- **Reading (descriptive, not a test of H140.1).** This is the situation H140.1 describes: strong shared non-clock
  variation across immune compartments of the same sample. On these data, cross-type concordance of an ordering that
  is not gated for clock shape can arise from it.
- **Process note.** During development the readout code path ran once on the primary set forced past the gate (40
  bootstraps, 100 permutations), as a code check. Its numbers were not kept or read into any decision.
- **Next.** Part (ii) (the published pipeline) waits for the Q1 simulator gate.
<!-- src: paper/exploratory/A140_dip_cyclops/output/q2/{part_i_summary,part_iii_readouts,orientation,primary_families}.csv -->

## A140 Q1 simulator gate (2026-10-27): FAILS; Q1 stops, H140.1-3 not scored

Arm D (one clock shared by three compartments, coupling 0.8, whole-day sampling), K = 1 non-clock factor, seed 1,500 genes
without core clock genes, 50 replicates. CYCLOPS 2.0 (port `common/benchmark/cyclops.py`, cell-type covariate, 10 models per
fit) phases against the true time: median circular correlation 0.067 (min 0.004), 0 of 50 at >= 0.7 (registered gate). The
three compartment phases agree with each other: median signed circular r 0.973, some replicates reflected at -0.99. The
DIP covers a median 0.29 of the circle. The clocksync gate on the same data passes in 43 of 50 replicates. The clock is
observable, but CYCLOPS orders the samples along the non-clock factor.
- **Deviations that may cause the failure:** 10 instead of 80 CYCLOPS models per fit (cost), and the factor calibrated on
  GSE266919 Pre samples (not gate-passing; eigenvalues 78.7, 50.6, 47.8 include Poisson noise).
- **Post hoc diagnostic** (8 replicates per scale, factor sd scaled): 0.1 gives median truth r 0.30 (1 of 8 >= 0.7); 0.25
  gives 0.16; 0.5 gives 0.16. Without the factor, CYCLOPS failed in its Julia code (8 of 8).
- **Reading:** this is not a test of H140.1-3. At this factor size the published pipeline does not recover a coupled
  clock that clocksync detects, and its compartments agree because they share the factor. Whether that holds at the
  template's 80 models and a cleaner calibration is PI_QUESTIONS Q11.
<!-- src: paper/exploratory/A140_dip_cyclops/output/q1/{simulator_gate,diag_summary}.csv, gate/D_K1_r*.csv -->

## A141 (scored 2026-10-27): brain contrast stability against noise; 4 of 12 contrasts are anchor-dependent beyond noise

Descriptive, as registered (PI, Q6). Brain contrast sets of WP0 item 8 (12 contrasts, 5 cohorts, A17 learned anchors).
The observed sign changes between the A17 and all-LOO hosts, recomputed here, equal the WP0 table.
- **Step 1, coverage:** the host phases cover the full day in most sets. Rayleigh R is within the uniform-sampling
  simulation range for MSSM (controls 0.027), HBCC (0.052), SEA-AD (0.063) and PD-atlas controls (0.071). It is above the
  range for RADC controls (0.195, n 26) and for PD-atlas all donors (0.149). GTEx fast brain: R 0.048 (13 tissues). No
  Rayleigh p is below 0.12.
- **Step 2, amplitude:** the median compartment amplitude is 0.558 (69 compartments) against the GTEx brain-region
  range 0.438-0.668 (median 0.576). Only 13 of 69 fall inside the range (29 above, 27 below), and cohort medians run from
  0.325 (RADC) to 0.760 (HBCC). Same scale on average, broad across compartments.
- **Step 3, decision:**
  - 8 of 12 contrasts are noise-driven, so not estimable at this design: AD, DLBD and SCZ (MSSM), BD and SCZ (HBCC),
    AD and DM (RADC), AD (SEA-AD).
  - 4 of 12 are not: PD_mssm (7 flips vs simulated 95th percentile 3), VaD_mssm (5 vs 4), ASCVD_hbcc (6 vs 3), and the
    PD atlas (15 of 34 vs 5).
- **Step 4, map** (agreement = same sign and |delta d| <= pooled SE): PD_mssm 0 of 9, VaD_mssm 3 of 9, ASCVD_hbcc 1 of 9,
  PD atlas 13 of 34. Expectations hold only partly. PsychAD neuron classes agree in 0-2 of 4. PD-atlas astrocytes and
  oligodendrocytes agree in 5 of 10, neurons in 2 of 8.
- **Reading:**
  - The simulation draws only sampling noise around one fitted host with uniform phases. "Not noise-driven" therefore
    means the anchored and all-LOO hosts differ by more than sampling noise: a structural anchor dependence (consistent
    with A135, no unique brain reference), not by itself a disease effect.
  - The noise-driven contrasts cannot be read. The other four depend on which compartments define the brain host.
  - Brain disease panels stay in Extended Data with this stated.
  - Per the registration, a time-of-death request is drafted for the PI to decide (outward-facing; not sent).
- **Choices beyond the registration** (README): state-level simulator, because `simulate` takes a scalar depth;
  per-contrast simulations; controls = the contrast with the most controls.
<!-- src: paper/exploratory/A141_brain_stability/output/{step1_coverage,step1_gtex,step2_summary,step4_summary}.csv, step3/step3_summary.csv -->

## WP0 / PI Q3 (2026-10-27): figure panels recomputed on the anchored estimator

Re-estimation with the current estimator, not a new hypothesis (`paper/exploratory/WP0_integrity/scripts/q3_anchored.py`;
targets read on their host, an anchor on the leave-one-out host of the other anchors). Gates as in
`a17_healthy_gates.csv` (all six pass).
- **GTEx organ departures (Figure 2, main): not recomputed.** The panel is timed: residuals against the recorded time of
  death, with no host, so it has no anchored version. PI (2026-10-27): keep it, with the readout named "timed, against
  recorded time of death, not the host-referenced estimator".
- **Neurons vs glia (Figure 2, ED), neuron-minus-glia mean coupling, 500 donor bootstraps, both frozen A17 rules:**

  | atlas (donors) | learned rule: anchors | learned | fixed (neuronal) rule |
  |---|---|---|---|
  | PsychAD HBCC (145) | OPC, astro, oligo, microglia | -0.111 [-0.250, 0.048] | 0.634 [0.528, 0.706] |
  | PsychAD MSSM (257) | OPC, astro, oligo, microglia | -0.155 [-0.231, 0.063] | 0.681 [0.610, 0.737] |
  | SEA-AD DLPFC (44) | OPC, astro, oligo, microglia, vascular | 0.102 [-0.017, 0.281] | 0.466 [0.311, 0.604] |
  | SEA-AD MTG (47) | OPC, astro, oligo, microglia, vascular | -0.003 [-0.129, 0.283] | 0.460 [0.303, 0.586] |
  | PD atlas controls (24) | 6 neuron + 5 microglia | 0.684 [-0.029, 0.727] | 0.675 [0.331, 0.726] |
  | retina (95) | 7 neuron + Mueller, astrocyte | 0.286 [0.167, 0.394] | 0.299 [0.176, 0.403] |

  - **Reading:** in the four cortical atlases the difference follows the anchors: -0.16 to +0.10 (every interval
    includes 0) under glial anchors, +0.46 to +0.68 under neuronal ones. "Neurons carry the brain host" holds only under
    neuronal anchors (consistent with A135, A141: no principled brain reference).
  - In the retina, neurons exceed glia under both rules, including when Mueller cells and astrocytes are anchors (archive
    A2, free fit: +0.17 [0.04, 0.35]).
  - PI (2026-10-27): the panel is restated and kept in ED (retina under both rules; cortex anchor-defined, shown under
    both rules). The message "the brain's backbone is neuronal" is dropped.
- **Fibrotic-lung T cells (Figure 5), A98 decomposition on the A22-confirmed anchors** (Natri: lung:EC_capillary
  dropped; HLCA: no drop, so identical to the M3 A98 table), 200 bootstraps, case minus control:
  - **HLCA fibrotic:** relative amplitude |z| 0.659 -> 0.000, diff -0.659 [-0.965, -0.126] (contributes). Coupling
    -0.147 [-0.276, 0.014]; own noise psi -0.225 [-0.329, 0.074].
  - **Natri ILD:** |z| 0.344 -> 0.000, diff -0.344 [-1.194, 0.197] (no change). Coupling -0.434 [-0.715, 0.174]. The
    control T fit is at the boundary after the anchor drop (M3, with the capillary anchor: |z| 0.360, diff -0.360
    [-1.259, 0.202]).
  - **Fibroblast specificity companion (same draws):** no change in either cohort (|z| diff 0.003 Natri, 0.163 HLCA;
    every interval includes 0).
  - **Reading:** in both cohorts the case T-cell loading on the host plane goes to the zero boundary. The interval
    excludes 0 in HLCA only. The panel stands as "loss of relative amplitude, interval excluding 0 in one of two
    cohorts", still conditional on A119 (procurement) for the main text. The A5 free-fit 0.45 is replaced.
<!-- src: paper/exploratory/WP0_integrity/output/q3_anchored/{panel2_neuron_glia,panel2_fits,panel3_differences,panel3_groups,panel3_anchors}.csv; paper/results/figure_5/output/lung_breaks/decomposition_differences.csv (M3); paper/results/logs/q3_anchored.log -->

## A140 Q1 revised simulator gate (2026-10-27): fails at full factor; CYCLOPS cannot start without the factor

Revision as registered (Q11 (b), 80 CYCLOPS models, factor calibrated on gate-passing tumour immune pseudobulks) and the
post hoc stop note (2b680e7). Arm D (one clock, three compartments, coupling 0.8, whole day), K = 1, seed 1,500 genes
without core clock genes.
- **Factor scale 1 (50 replicates): FAILS.** CYCLOPS phase vs truth, median circular r 0.060 (q90 0.178, max 0.273);
  0 of 50 >= 0.7. DIP vs truth 0.072. Cross-compartment concordance median 0.967; median rhythmic genes called per
  compartment 536. clocksync gate passes 24 of 50 (below the 35 of 50 H140.3 asks for in arm A; descriptive).
- **Scale 0.5 (14 replicates, interrupted queue; descriptive):** median 0.079 (max 0.204), concordance 0.979; gate 11 of 14.
- **Scale 0 (20 replicates):** every CYCLOPS fit fails (20 of 20, status error). Cause, reproduced on replicate 0:
  `SVDReduceDimensions` finds no eigengene above `eigen_contr_var` (0.05) among the seed genes (`findlast` returns
  nothing; Julia `MethodError`), i.e. without the factor no seed-gene principal component carries 5 % of the variance.
  clocksync gate passes 20 of 20.
- **Consequence (registered):** scale 0.25 runs only if scale 0 passes, so it does not run; H140.1-3 are not scored.
  Whether to try the registered core-clock seed is a design change after a result (PI).
- **Reading:** at this design CYCLOPS 2.0 orders the samples along the shared non-clock factor when there is one (high
  cross-compartment concordance, no clock recovery) and cannot start when there is none. clocksync's gate keeps
  sensitivity when the factor is small (20 / 20) and loses it at full factor size (24 / 50).
<!-- src: paper/exploratory/A140_dip_cyclops/output/q1rev/gate/D_K1_r*.csv, D_K1_s0.5_r*.csv, D_K1_s0_r*.csv (columns truth_r_median, dip_truth_r, cc_median, n_rhy_min, gate_pass, status); paper/results/logs/q1_run.log, q1_run_rev_gate_s1.log -->

### 2026-10-27, PI Q3 follow-up: the A21 colon shift without the decoupling anchor (Figure 1 panel restored)

Re-estimation, not a new hypothesis. A22 V17.2 failed because the rule did not drop lymphatic EC, known to decouple after
DSS (coherence 0.448 -> 0.048), so the R-DC host included a decoupled anchor. R21.4 recomputed as in M3 with anchors
fibroblast, smooth muscle, pericyte and endothelial only (500 animal bootstraps):
- Epithelial mean shift R-DC minus C-DC: untimed -8.95 h [-10.87, -7.04] vs known time -9.03 h [-11.03, -7.09];
  difference 0.08 h. With lymphatic EC as anchor (M3): -8.70 h vs -9.40 h. The R21.4 acceptance (within 3 h, both
  intervals excluding 0) is met either way.
- Per epithelial class, untimed vs known time: colonocyte -7.44 vs -7.72 h; secretory -9.12 vs -8.74 h; stem -10.27 vs
  -10.66 h.
- Lymphatic EC as a target: -1.02 h [-8.10, 4.38] (known time +1.48 h); no shift detectable.
- Reading: the blind recovery of the epithelial shift does not depend on the decoupled anchor. The Figure 1 panel is
  restored (PI, 2026-10-27) with the anchors stated: set by hand from prior knowledge, because the automatic rule
  failed V17.2. It is a known-answer check of the readout, not of the anchor rule.
<!-- src: paper/exploratory/WP0_integrity/output/a21_shift_without_lymphatic.csv; paper/results/figure_1/output/species/a21_r214_phase.csv (M3) -->

## A142 (scored 2026-10-27): the brain split is brain-specific; digestive organs carry a smaller departure of their own

GTEx v10 fast deaths, 14-gene plane, A103 statistic generalised to any tissue block; 1,000 random size-matched peripheral
blocks per comparison; 200 donor bootstraps.
- **H142.1 PASS (Holm p 0.002).** Brain vs periphery r* 0.639 [0.530, 0.728], 205 donors. Random 13-tissue peripheral
  blocks: median 0.978, 5th percentile 0.936, SD 0.025. One-sided p 0.001 (the floor). Same without PER1/PER2 (0.585)
  and after ischaemic-time and Hardy adjustment (0.640), each p 0.001. The brain separates far more than any group of
  peripheral organs.
- **H142.2 PASS (Holm p 0.028).** Digestive (7 tissues) vs the other 20 peripheral tissues: r* 0.883 [0.819, 0.971],
  160 donors. Random 7-tissue blocks: median 0.970, 5th percentile 0.900. One-sided p 0.028.
  - Sensitivities: adjusted 0.884 (p 0.031); without PER1/PER2 0.896 (p 0.056).
  - The digestive departure is real but about a third the size of the brain's (1 - r*: 0.12 vs 0.36).
- **Descriptive blocks:**
  - Food-exposed (stomach, transverse colon, oesophageal mucosa, liver; 107 donors): r* 0.922, p 0.31. Too few tissues
    and donors to read: its reference SD is 0.11, and its bootstrap upper bound is undefined.
  - Cardiovascular 0.975 (p 0.57) and endocrine/adipose 0.964 (p 0.47): at the reference.
- **Deviation (data rule, not registered):** reference blocks with an undefined r* (low-donor tissues) were dropped:
  0 of 1,000 at k = 7 and 13, 15 at k = 4, 0 at k = 5 (base plane).
- **Reading:**
  - At tissue level, the brain's partly separate internal time is specific to the brain. The A143 processing checks
    are running.
  - Digestive organs depart from the rest of the body as a group, more than an arbitrary organ group does. This is
    consistent with meal timing (PI hypothesis) but not shown to be meal timing: GTEx has no meal times, and the block
    includes smooth-muscle layers.
<!-- src: paper/exploratory/A142_organ_blocks/output/{tests,targets,reference_summary}.csv -->

## A144 part 2 (D144.2, descriptive; 2026-10-27): per-cell amplitude vs within-sample synchrony as carriers of coupling: inconclusive

33 eligible healthy or control compartments (AD_mssm 7, AD_radc 5, SEA-AD DLPFC 6, HLCA fibrotic controls 5, HLCA lung
5, PD controls 5; none dropped). Coupling c from each group's anchored refit, reproducing the A18 rho to 6.5e-16.
logit(c) ~ log lambda + log R + target, 1,000 bootstraps of compartments within target:

| coefficient | estimate [95 %] |
|---|---|
| beta_lambda (per-cell amplitude) | 1.15 [0.46, 2.08] |
| beta_R (within-sample synchrony) | 1.88 [0.75, 2.83] |
| beta_R - beta_lambda | 0.73 [-0.48, 1.90] |
| beta_rho (mechanical reference, logit(c) ~ log rho + target) | 1.69 [0.89, 2.47] |

- Both coefficients are positive, as the shared log rho term forces.
- Synchrony's coefficient exceeds per-cell amplitude's, but the interval of the difference includes 0. So the
  population view is not shown.
- The theory's signature, a per-cell amplitude coefficient below the mechanical reference, is not shown either:
  beta_lambda 1.15 vs beta_rho 1.69, overlapping intervals.
- At 33 rows, mostly brain, this cannot separate the explanations.
- Caveat: HLCA fibrotic controls and HLCA lung may share donors (not checked).
<!-- src: paper/exploratory/A144_strong_clocks/output/{d144_2_coefficients,d144_2_rows}.csv -->

## A143 (scored 2026-10-27): the brain split is not explained by sample quality or isolation batch; no hypothalamus exception

GTEx v10 fast deaths, A103 base statistic; v8 sample attributes for 80 % of brain and 89 % of peripheral samples (note 1:
no autolysis term).
- **H143.1 PASS.** r* 0.627 [0.550, 0.712] after adjusting every tissue for RIN, ischaemic time, Hardy and centre. On the
  same samples without adjustment: 0.621 [0.514, 0.723]. Quality does not move the split. The autolysis sensitivity is
  not estimable: only 2 brain regions are scored, so the brain reliability is undefined.
- **H143.2 FAIL (no brain-specific batch effect detected; Holm p 0.38).**
  - Brain-region pairs from the same isolation batch correlate 0.039 more than pairs from different batches (30 pairs,
    p 0.19; null SD 0.044, q95 0.070).
  - Peripheral pairs show the same excess, detected with more pairs: 0.040 (155 pairs, p 0.0035).
  - Reading: isolation batch carries a small shared departure in every tissue, of the same size in brain and periphery.
    It is not brain-specific, and an excess above about 0.07 in brain pairs is excluded at this power.
- **H143.3 FAIL.** Hypothalamus alignment ratio 0.582 vs the median of the other regions (about 0.66); difference
  -0.082 [-0.248, 0.152], p 0.74.
  - Other regions span 0.37 (frontal cortex BA9) to 0.81 (hippocampus); the spinal cord is at 1.18 and follows the body.
    Whether this spread exceeds noise was not tested.
  - Hardy 3-4 ratios (descriptive) are noisy (-0.83 to 1.59).
  - The A110 slow-death observation is a different estimand and stands as descriptive.
- **Reading (registered rule):** H143.1 passes and H143.2 fails, so the split is not explained by sample quality or by
  isolation batch. The brain claim keeps "partly separate" in Figure 2, with A142 (brain-specific) and A143 cited. No
  hypothalamus panel. Not excluded: brain-only processing steps without a recorded variable (dissection protocol).
<!-- src: paper/exploratory/A143_brain_artefact/output/{coverage,h143_1,h143_2,h143_2_pairs,h143_3_regions,tests}.csv -->

## A144 part 1 (scored 2026-10-27): the strong-clock association is not shared noise; it lies in the shared part

GTEx fast deaths, 27 peripheral tissues, 100 random donor halves.
- **H144.1 PASS.** Amplitude from one half against coupling from the other: cross-half rho 0.690 (5-95 % over splits
  0.645-0.732). Within the same half: 0.718. Permutation null mean -0.001, SD 0.184; p 1e-4 (the floor). Noise shared
  within donors accounts for 0.028 of the association (5-95 %: -0.015 to 0.088).
- **D144.1.**
  - Across tissues, log locked variance varies about 3-fold more than log unlocked variance (SD 1.02 vs 0.38).
  - Coupling follows the locked part (rho 0.86), not the unlocked part (-0.08).
- **Reading:** strong peripheral clocks are strong by the part they share with the body. This is consistent with cell
  synchrony maintained by body-wide signals, or with stronger cue input, rather than single-oscillator rigidity, which
  these data do not test (part 2 was inconclusive). The Figure 2 row stays, worded this way.
<!-- src: paper/exploratory/A144_strong_clocks/output/{h144_1,d144_1_components,splits}.csv -->

## A138 power step (2026-10-27, before any observed A138 statistic): H138.4 powered, H138.5 weak, H138.1 test anti-conservative
- **H138.4** (TCGA, 517 tumour-normal pairs, 10 gated cancers): null 0.122 ± 0.021 (q95 0.157), MDE 0.174, smallest p
  1e-4. Powered.
- **H138.5** (single cell, 52 paired patients, 6 cohorts): null 0.443 ± 0.092 (q95 0.600), MDE 0.677, smallest p 1e-4.
  Detectable only near the null's upper range. Fibroblast companion (ERP149257 only): null 0.526 ± 0.219, MDE 1.05,
  so uninformative.
- **H138.1** (`clocksync.simulate` at the observed designs of the Q1 cohorts; 30 replicates per c, 100 bootstraps per
  fit, note 4): share of replicates whose pooled lower bound > 0 is 0.133 at c = 0, 0.667 at 0.1, 1.000 at 0.2, 0.3
  and 0.5. MDE 0.2 (registered rule: >= 24 of 30).
  - **The registered test is anti-conservative:** 4 of 30 null replicates pass (false-positive rate 0.133; exact 95 %
    interval about 0.04-0.31) against a nominal 0.05. The pooled estimate at c = 0 is biased upward (mean 0.030, q95
    0.083), as expected for a coupling bounded at 0; 100 rather than 1,000 bootstraps may add to it.
  - Fit failures: 20 of 2,550 cohort fits (2, 2, 4, 4 and 8 from c = 0 to 0.5), all a resampled SRP273482 without
    the normal fibroblast anchor; those cohorts leave that replicate's pool. 16 cohorts per pool.
<!-- src: paper/exploratory/A138_cancer_estimability/output/power/{power_summary,h138_1_sim,h138_1_sim_cohorts}.csv -->

## A115 (scored 2026-10-27, descriptive): uninvolved psoriatic skin already runs its core clock late against healthy skin
Healthy rows added from the raw files with the original builders; the rebuilt lesional and uninvolved rows equal the
processed ones exactly in all 4 series (acceptance passed, nothing dropped). 185 uninvolved vs 138 healthy samples.
- **D115.1 (uninvolved minus healthy core-cue offset, pooled):** -1.806 h [-3.811, -0.689] (2,000 bootstraps);
  one-sided permutation p 0.0084 (10,000). Pooled offsets: uninvolved -0.099 h, healthy +1.706 h. Coherence 0.389
  (uninvolved) and 0.585 (healthy). Permutation null SD 0.768 h (MDE 1.91 h; registered blind estimate 2.19-2.89 h).
  Registered reading rule: "a trait component of at least about 0.7 h".
- **Per series:** all 4 negative, none alone below p 0.05: GSE13355 -1.52 h (p 0.058), GSE121212 -4.73 (0.076),
  GSE201827 -0.64 (0.31), GSE268958 -1.92 (0.087).
- **Leave one series out (post hoc):** -1.33 to -2.53 h, p 0.008-0.048 (2,000 permutations).
- **Atopic dermatitis (GSE121212, 21 vs 38):** -0.43 h [-9.36, 8.70], p 0.42, MDE 6.9 h. Uninformative.
- **Keratinocyte differentiation score** (DIFF minus BASAL; Hedges g, uninvolved minus healthy): psoriasis pooled
  +0.47 [0.24, 0.70], all 4 series positive; atopic dermatitis -0.08 [-0.61, 0.44].
- **Caveats:**
  - Healthy controls are different persons, possibly sampled at other body sites or in other batches. The
    permutation assumes the two groups are exchangeable within a series, which such differences would break; site and
    batch cannot be checked: none of the 4 series records body site, batch or
    processing date in GEO (only GSE268958 records sex).
  - The trait estimate (1.8 h) exceeds the within-group lesional-minus-uninvolved estimate on 15 series (1.05 h; paired
    estimator 1.70 h). Different series and estimator, but the size is surprising and needs the site / batch check
    before it is read as biology.
  - Uninvolved skin is also more differentiated than healthy skin (g 0.47), so the clock delay and the keratinocyte
    state co-occur; the transit model would predict that, but it is not tested here.
<!-- src: paper/exploratory/A115_A116_psoriasis/output/{a115_acceptance,a115_d115_1,a115_atopic,a115_keratinocyte,a115_loso}.csv -->

## A138 Q3 and Q1 (scored 2026-10-27): tumour clocks follow their own patient's normal tissue (H138.4 pass); malignant host coupling not shown (H138.1 fail)
Regenerated nulls equal the committed power nulls (mean, SD, q95, estimable flags).
- **H138.4 PASS** (TCGA, 517 tumour-normal pairs, 10 gated cancers): pooled Poisson-corrected |correlation| of tumour
  and own-patient normal clock states 0.365; pairing null 0.122 ± 0.021 (q95 0.157); p 1e-4 (the floor). Pseudo-clock
  companion (200 expression-matched gene sets): median 0.194, q95 0.261, none reaches 0.365, so the pairing signal is
  clock-specific. Per cancer: BLCA 0.52, LIHC 0.49, LUAD 0.49, PRAD 0.43, KIRC 0.43, BRCA 0.32, THCA 0.30, KIRP 0.29,
  COAD 0.21, UCEC 0.08; p < 0.01 in 7 of 10 (COAD 0.19, KIRP 0.072, UCEC 0.87).
- **H138.5 (descriptive):** malignant vs adjacent-normal epithelium after removing the stromal host, 52 patients, 6
  cohorts: 0.583, null 0.443 ± 0.092, p 0.069; below the MDE (0.677). Not detected at this power. Fibroblast companion
  (1 cohort, 6 patients): 0.240, p 0.87, uninformative. GSE131907 excluded by note 4 (ii).
- **H138.1 FAIL** (note 5 rule): pooled anchored coupling of malignant cells 0.076 [0.014, 0.139] (DerSimonian-Laird,
  17 cohorts, 225 patients, tau2 0), below the calibrated threshold 0.083; p 3/31 = 0.097. The registered lower-bound
  rule would pass (0.014 > 0): the case note 5 was written for. 15 of 17 cohorts have a lower bound of 0. The power
  simulation used the design count (270 patients); fewer patients were fitted (225), which widens the true null, so
  the threshold is if anything lenient. Reading: malignant host coupling, if present, is below about 0.2 (MDE).
- **Holm (provisional; Q2 unscored):** H138.4 adjusted 2e-4 (over the scored two) or 6e-4 (family of 6, unscored taken
  as not smaller); H138.1 0.097 or 0.48.
- **D138.1 (descriptive):** median shares of the clock-plane variance (host-coupled z2V, patient-specific 2psi,
  Poisson): malignant 0.18 / 0.79 / 0.08 (17 cohorts); normal epithelium 0.34 / 0.49 / 0.20 (9); stroma 0.53 / 0.27 /
  0.18. Patient-specific share malignant minus normal epithelium, 7 cohorts with both: mean +0.195, positive in 6 of 7.
  The between-cell part has no registered estimator (the A18 per-group D was reused, 11 of 17 cohorts, very noisy) and
  is not read.
- **Reading (Q8 tension):** the clock state of a tumour keeps a patient-specific relation to the same patient's normal
  tissue, yet malignant cells carry little variance coupled to the stromal host and most of their clock variance is
  patient-specific. This matches suggestion 5 of the status note: malignant clocks are not gone; their between-tumour
  variation is dominated by patient-specific sources. It contradicts "no fixed timing relationship with healthy tissue"
  (CoFE, Li et al. 2024) at bulk level.
<!-- src: paper/exploratory/A138_cancer_estimability/output/score/{tests,per_unit,holm,d138_1,d138_1_contrast}.csv -->

## A146 part 1 (scored 2026-10-27): breast tumour programs follow recorded surgery time in the predicted direction, not internal time
GSE233242, 42 patients, 10 programs (AURKB absent; P1 on 9 genes), 10,000 time permutations, seed 146.
- **H146.1 FAIL (internal time):** omnibus S -0.056, p 0.55 (null SD 0.437, MDE 1.09).
- **H146.2 PASS (recorded resection time):** S 1.339, p 0.0006, Holm 0.0012 (null SD 0.381, MDE 0.94). Strongest:
  T-cell activation residual (P4) t 3.48 (p 0.0011), PD-1 residual (P10) 3.22 (0.0014), MHC-I (P2) 2.15 (0.018).
- **Internal minus recorded (D146.1):** -1.39 [-2.40, -0.47] (2,000 patient bootstraps): recorded time explains the
  programs better, the opposite of the registered clock reading.
- **Mapping offset (shown, not tested):** recorded-time S peaks at the registered mapping (offset 0: 1.34; -2 h: 1.28;
  +-6 h: +-0.44); internal-time S stays between -0.29 and 0.29 at every offset.
- **Adjacent normal (descriptive, P2-P10):** internal time S 1.03 (p 0.007; interferon t 2.90, phagocytosis 2.14,
  MHC-I 1.93); recorded time 0.64 (p 0.028; MHC-I 3.53).
- **Reading:** in tumours, the time-of-day pattern predicted from the mouse atlas is present against the hospital
  clock but not against the tumour's internal (clock-gene) time; the H146.1 MDE was below the H146.2 effect, so this
  is not only power. Candidate origins (not tested): procedural or scheduling covariates that track the hour of surgery
  (pre-operative fasting duration, anaesthesia, case order), or wall-time systemic cues acting on immune programs
  independently of the tumour clock (A85: immune cells lean on body-wide cues). In adjacent normal tissue the programs
  follow internal time (descriptive). Caveat: `est_day` is a noisy per-patient estimate (Spearman 0.56 with recorded
  time, MAE 1.1 h), which attenuates any internal-time association.
<!-- src: paper/exploratory/A146_tme_time/output/{part1_programs,part1_tests,part1_per_program,part1_offset_sweep,part1_internal_vs_recorded}.csv -->

## A146 part 2 (scored 2026-10-27): pan-cancer single cell, programs against internal time: not detected
16 gate-passing cohorts, 213 scored patients (226 timed), 10,000 within-cohort time permutations.
- **H146.3 FAIL:** omnibus S 0.528, p 0.092, Holm (H146.1-3) 0.184; null SD 0.395, MDE 1.00. Per program all |t| < 2
  (largest: T-cell activation in T/NK 1.83, PD-L1 in malignant 1.31, T fraction 1.04; proliferation -0.94).
- A146 Holm over the three: H146.1 0.55, **H146.2 0.0018 (pass)**, H146.3 0.18.
- **Descriptive:** clock-keeping cancers (PRAD, LIHC; 25 patients) S -0.33; non-keeping (COAD, 13) -0.08; breast
  (ERP130387, 11) -0.48; all far below their MDEs (1.0-1.15). Offset sweep: S peaks at 0 to +2 h (0.53), no offset
  significant.
- **Reading:** across cancers the predicted organisation is not detected against internal time; an effect as large as
  the recorded-time effect in breast (1.34) would have been detected. Together with part 1, the time-of-day signal of
  the tumour microenvironment, where seen, follows the hospital clock, not the tumour's internal time. Without recorded
  times outside GSE233242, part 2 cannot test the recorded-time reading.
<!-- src: paper/exploratory/A146_tme_time/output/{part2_tests,part2_per_program,part2_per_cohort,part2_subsets,part2_offset_sweep}.csv -->

## A121 (scored 2026-10-27): embryo clock onset not shown within the late window (H121.1 fail); phase uninformative (H121.2)
Qiu 2024 mouse embryo atlas, host as H108.1 (late-window statistic reproduces H108.1: R 0.516).
- **H121.1 FAIL:** R(E15.0-16.75) 0.392 (p 0.41), R(E17.0-18.75) 0.594 (p 0.081); difference 0.202, one-sided p 0.226
  (10,000 within-batch permutations), Holm 0.45. 8 embryos per half: underpowered for a moderate rise.
- **Lineages (secondary, descriptive, uncorrected):** hepatocytes 0.375 -> 0.766 (late p 0.008; difference p 0.048);
  mesoderm 0.157 -> 0.655 (difference p 0.022).
- **H121.2 UNINFORMATIVE (registered rule):** calibrated frame offset 21.54 h (SE 0.30 h) from 11 timed mouse datasets
  (A21, A48, A72), but the datasets disagree by 3.07 h (> 3 h; A48 aged, the noisiest, against A72 GSE190032). The
  embryo offset itself has SE 3.29 h (16 embryos), so a 95 % interval spans about 13 h. Deviation: the registration
  did not define "the embryo clock peak"; three definitions were computed (`h121.csv`, column `definition`) and they
  disagree on overlap with ZT12-14, which does not matter under the uninformative verdict but must be fixed before any
  rerun.
- **Registered reading:** H121.1 fails, so the embryo panel goes to Extended Data (PI to confirm, below).
<!-- src: paper/exploratory/A121_embryo_timing/output/{h121,frame_offsets,lineages,embryo_offset_boot}.csv -->

## A139 (scored 2026-10-28, primary cell 08-18 h, n = 80, 100 replicates): gene-timing refit and flags fail; the E4 likelihood-ratio flag nearly calibrated
Deviation (allowed by the registration's ordering note): the full grid (about 90 h) was not run; the primary cell
decides V139.1 and V139.2, the rest of D139.3 is not computed (run stopped after the 400 primary tasks; cache kept).
Implementation choices listed in `scripts/` docstrings (E-b simultaneous update with gauge fixing, bootstrap capped
at 3 rounds; target coupling 0.8, offset +3 h).
- **V139.1 FAIL.** Median absolute error, shifted / unshifted genes (limits 1.5 / 1.0 h): S1 (one gene +6 h, E-a)
  1.05 / 0.89 h (passes); S3 (three genes, E-a) 1.88 / 1.12 h (fails both); SW (waveform wrong in every compartment,
  E-b) 1.71 / 0.71 h (fails shifted). Against the gene-centred truth E4 can read: 0.77-1.00 h shifted.
- **V139.2 FAIL.** False flags in S0 (no shift): bootstrap flag up to 0.15 per gene (E-a; 13 of 14 genes above the
  bound) and 0.54 (E-b); power to flag a 6 h shift: S1 1.00, S3 0.67 (minimum gene), SW 0.99.
- **E4 likelihood-ratio flag (descriptive):** false-flag rate 0.06 in S0 (95 % upper bound 0.115), power 1.00 in S1,
  S3 and SW.
- **Reading (registered):** both fail, so the conserved waveform stays an assumption, stated as such, and no per-gene
  timing result is quoted as gene-specific. The likelihood-ratio flag is close to calibrated and fully powered here;
  a future registered check could adopt it (not done now). No failed fits.
<!-- src: paper/exploratory/A139_waveform_refit/output/{v139,d139_3,gene_level,failures}.csv -->

## A147 (scored 2026-10-28): pull-versus-order readout: order validated, shrink ratio biased at n 79 (overall FAIL)
A133 generative model, 200 replicates per cell, seed 147, 0 failures.
- **V147.1 PASS (order separates loss from pull, n 79):** AUC 0.977 between pull f 0.44 and loss q 0.5 (n 298: 1.000);
  at f 0.8, 0.805 of replicates keep order above the loss median (n 298: 0.960). Order null calibrated: false-positive
  rate 0.050-0.065 in the full-loss arm.
- **V147.2 FAIL:** the shrink ratio falls monotonically with the pull (Spearman of cell medians with f -1.0 at both n;
  medians 0.81, 0.65, 0.45, 0.14 at f 0, 0.2, 0.44, 0.8, n 79), but with no pull it reads 0.809 at n 79 (limit
  [0.85, 1.15]; n 298: 0.929): biased low by noise in the reference deviations.
- **Registered consequence:** the readout stays descriptive (Extended Data); not added to clocksync.
- **Implementation choice that matters (not specified in the registration):** the simulated reference had the group's
  n (79). The real fast-death reference is larger (244 in A132b), which would reduce the bias. A rerun with the real
  reference size would be post hoc; parked for the PI. A draft `pull_order` function and tests are in
  `paper/exploratory/A147_pull_order/scripts/`.
<!-- src: paper/exploratory/A147_pull_order/output/{v147,d147,failures}.csv -->

## A122 (scored 2026-10-28): the A85 immune cue-following excess is not specific to cue genes (H122.1, H122.2 fail)
A85 reproduced exactly (immune 0.182, structural 0.095, difference +0.087, p 0.002; 70 classes, 9 atlases).
- **H122.1 FAIL:** with the locked set the immune-minus-structural difference is +0.080 (immune 0.217, structural
  0.137); cue minus locked +0.007 [-0.090, 0.098], p 0.44 (atlas bootstrap), Holm 0.55.
- **H122.2 FAIL:** 200 expression-matched pseudo-cue sets give immune-minus-structural differences with median 0.072,
  95th percentile 0.114 (observed 0.087), p 0.27, Holm 0.55.
- **Descriptive:** the cue difference is positive in all 8 atlases with both classes; cue minus locked is positive in
  5, negative in 3. The index falls with pseudobulk depth across classes (Spearman -0.31 for both sets); adjusting for
  depth reduces the immune coefficient from 0.128 to 0.106 (cue) and 0.098 to 0.072 (locked).
- **Registered consequence:** Extended Data, described as a consequence of low host coupling (the index rises for any
  gene set in classes with a small host loading). The Figure 2 immune panel is dropped. A146's reading no longer
  cites A85 as support for "immune cells lean on body-wide cues".
<!-- src: paper/exploratory/A122_cue_controls/output/{reproduction,h122,per_atlas,depth,composition,pseudo}.csv -->

## A145 (scored 2026-10-28): two-step orderer pipelines recover coupling nearly as well as clocksync; no registered advantage survives Holm
Simulation (14 clock genes, 100 donors, 4 anchors, 5 targets with c 0-1), 30 replicates per arm, 0 orderer failures;
CYCLOPS with the 80-model template. Coupling MAE per method (clocksync / CHIRAL / CYCLOPS / CoFE / timed oracle):
| arm | clocksync | CHIRAL | CYCLOPS | CoFE | oracle |
|---|---|---|---|---|---|
| F full day | 0.044 | 0.046 | 0.055 | 0.047 | 0.034 |
| W 08-18 h | 0.060 | 0.097 | 0.074 | 0.206 | 0.048 |
| S_small (shared factor 0.066 of plane variance) | 0.117 | 0.115 | 0.122 | 0.117 | 0.086 |
| S_large (0.22) | 0.239 | 0.234 | 0.278 | 0.258 | 0.174 |
- **H145.1 (working hours) FAIL after Holm:** clocksync 0.060 vs best two-step (CYCLOPS) 0.074; clocksync better in
  21 of 30 replicates, sign-test p 0.021, Holm 0.064. Offsets in W: median error 0.14 h (clocksync) vs 0.25-0.34 h.
- **H145.2 (shared factor) FAIL:** 0.239 vs CHIRAL 0.234 (13 of 30), p 0.82.
- **H145.3 (no clock) FAIL:** clocksync's gate passed in 0.133 of no-clock replicates (4 of 30), and its false-claim
  share (gate pass and c > 0.2) was 0.127 of targets, against CoFE two-step 0.06 (CHIRAL 0.10, CYCLOPS 0.15); p 0.99.
- **D145.1 (full day):** all pipelines equivalent (Spearman with true c 0.96-0.98; offset circular r 0.998).
- **D145.2 (real known answer, Bmal1 liver):** every pipeline sees the coupling drop; clocksync -0.293 (A111 reference
  -0.291), timed oracle -0.186, CYCLOPS -0.137, CHIRAL -0.101, CoFE -0.092.
- **Reading (registered):** coupling, offsets and amplitudes are recoverable by a two-step orderer pipeline when
  sampling covers the day, and nearly so under working hours in this clock-gene-only simulation; clocksync's
  specificity claims H145.1-3 are not supported (H145.1 nominal only). Two cautions for the method text: (i) the gate
  passed 4 of 30 times with a shared non-clock factor and no clock (nominal 0.05; binomial p about 0.06), the same
  sensitivity A140 showed in the other direction; (ii) here orderers see only clock genes, while A123/A140 show that
  genome-scale orderers fail under working hours or with strong non-clock structure, so the advantage on real data
  rests on A123 (CHIRAL r -0.17 vs 0.945) and A140, not on A145.
<!-- src: paper/exploratory/A145_two_step/output/{h145,d145_1,d145_2,failures}.csv -->

## A138 H138.8 (scored 2026-10-28, secondary): the tumour-normal clock match does not depend on tumour purity
PanCanAtlas ABSOLUTE purity (475 of 517 pairs called); H138.4 reproduced (0.3647) before splitting.
- **H138.8 FAIL:** low- minus high-purity match D +0.025 (pooled low 0.427, high 0.402; 8 cancers, 221 / 213 pairs),
  p 0.35; null SD 0.065, MDE 0.16 (less than half the H138.4 statistic). BLCA and UCEC dropped (< 10 pairs per half).
- **Per cancer:** mixed signs (PRAD +0.33, KIRC +0.19; LIHC -0.19, THCA -0.15, LUAD -0.11; others within +-0.05).
  Pair-level Spearman of purity with the pair contribution -0.32 to +0.15, all p >= 0.19.
- **ESTIMATE sensitivity (descriptive):** D -0.029, p 0.67 (ESTIMATE vs ABSOLUTE purity Spearman 0.38-0.80).
- **Reading (registered, "fail with D near 0 and the MDE small"):** the bulk match is not carried mainly by the
  non-malignant fraction; malignant cells contribute to it. With H138.1 (weak coupling to the stromal host) this points
  to a patient-specific component shared by malignant and normal tissue that does not run through the tumour stroma
  (e.g. the systemic state or time of surgery shared by both samples), consistent in direction with H138.5 (0.58 vs
  null 0.44, p 0.069). Caveat: the purity range within a cancer is moderate (median split), and ESTIMATE was
  re-implemented (R package unavailable).
<!-- src: paper/exploratory/A138_cancer_estimability/output/h138_8/{tests,per_cancer,pair_level,estimate_sensitivity,h138_4_reproduction}.csv -->

## A148 (2026-10-28): stopped at the registered acceptance check (bulk deconvolution not valid)
The deconvolved lesional-minus-uninvolved keratinocyte shares did not move as in single cell: cycling -0.020
[-0.050, 0.011], most differentiated -0.019 [-0.045, 0.005] (3 series; GSE13355 is gene-centred log2, so CPM is
undefined and it was excluded; 4-series sensitivity with imputed levels also fails). Deconvolution quality was poor:
immune share 0 everywhere, 46-82 % of the squared signal unexplained, because the registered gene set (most variable
across the five keratinocyte bins) carries almost no non-keratinocyte markers. No predicted offset was computed;
H148.1-2 not scored.
- **Post hoc diagnostic (not registered):** with genes variable across all nine references, the direction is right
  but weak (cycling +0.019 [-0.012, 0.050], most differentiated -0.020 [-0.042, 0.002]). Not used.
- **Reading:** bulk microarray / RNA-seq deconvolution into keratinocyte states is too fragile here to decompose the
  delay. The differentiation argument rests on the single-cell results (A61, A63, A69) and on the A115 descriptive:
  uninvolved skin is more differentiated than healthy (g +0.47), which the transit model predicts would run ahead,
  not 1.8 h behind.
<!-- src: paper/exploratory/A148_psoriasis_decomposition/output/{acceptance,shares,genes_selected,acceptance_posthoc}.csv -->

## Erratum to the A138 H138.8 reading (2026-10-28)
The H138.8 reading said the shared tumour-normal component "does not run through the tumour stroma". That reads a
failed test (H138.1: malignant coupling to the in-sample stromal host 0.076 [0.014, 0.139], MDE about 0.2) as a
negative, which the statistics rule excludes. Corrected reading: the bulk match is carried in part by malignant cells
(H138.8); whether it runs through the tumour stroma is not resolved at the single-cell designs. The sizes are not in
conflict: a Poisson-corrected state correlation of 0.365 corresponds to about 0.13 of shared clock variance (0.365^2,
back-of-envelope, different models and samples), inside the H138.1 interval and near the descriptive malignant
host-coupled share (0.18, D138.1); H138.5 points the same way (0.58 vs null 0.44, p 0.069). The H138.1 and H138.5
designs are underpowered for an effect of that size. Follow-up: pooled paired tumour/normal single-cell cohorts (A149).

## A146 note 2 checks (scored 2026-10-28, secondary family "A146 checks"): the breast recorded-time pattern is not scheduling, handling or a generic trend; the glucocorticoid control is in direction but not detected
GSE233242, 42 tumours, part 1 permutation stream.
- **H146.4 FAIL (glucocorticoid score falls with resection hour):** t -1.61, one-sided p 0.054, Holm 0.054 (null SD
  1.02 in t units; Spearman -0.22). In the predicted direction, not detected. Normal tissue (descriptive): t +0.98.
- **H146.5 PASS (scheduling):** S with programs residualised on age and subtype 1.234 (unadjusted 1.339), p 0.0005,
  Holm 0.001. Resection hour is not associated with age (Spearman 0.18, p 0.25), subtype (Kruskal-Wallis p 0.27;
  Luminal A 28, TNBC 9, Luminal B 3, HER2 2), T-cell abundance (0.14, p 0.37) or myeloid abundance (0.19, p 0.22).
- **D146.2 (handling):** immediate-early score t -0.29, hypoxia score t +1.04 (CA9 absent); S after residualising on
  both 1.403, p 0.0004.
- **D146.3 (trend vs predicted phases):** with the 10 predicted hours shuffled among programs, S averages 0.342 (q95
  1.126); 0.97 % of shuffles reach the registered 1.339, so the predicted phase structure, not a common trend, carries
  the score. Immune-cell programs S 1.354 (p 0.0027); the other programs S 1.318 (p 0.025; proliferation t 1.03,
  MHC-I 2.15, interferon 1.15, PD-L1 0.94): not confined to immune cells.
- **Reading (registered):** scheduling by age or subtype, the handling markers and a generic monotone trend do not
  explain the pattern. The glucocorticoid positive control is in the predicted direction but not detected (p 0.054),
  so the shared-wake-time systemic-cue reading is neither confirmed nor excluded; fasting duration remains confounded
  with the hour. Programs of tumour-intrinsic as well as immune origin follow the hour. Replication in a second cohort
  with recorded resection times remains the decisive step.
<!-- src: paper/exploratory/A146_tme_time/output/{checks_tests,checks_descriptive}.csv -->

## A123b (scored 2026-10-28): on the 40-tissue organ map, untimed orderers fail in working hours and recover with the whole day; clocksync holds in both
GTEx v10 fast deaths, 40 tissues; W = 111 donors who died 08-18 h, F = all fast deaths (up to 244 donors). No method
failed except CoFE clock genes in 1 W tissue.
- **H123.4 PASS (W, truth (i)), all three:** clocksync offset circular r 0.945 vs CoFE genome-wide -0.283 (difference
  1.23 [0.88, 1.48]), CoFE clock genes -0.417 (1.36 [0.95, 1.69], 39 tissues), CYCLOPS 2.0 -0.296 (1.24 [0.81,
  1.67]); each one-sided p 0.0005 (floor), Holm 0.0015. CHIRAL (known): -0.165.
- **Offsets, circular r vs truth (i), W / F; median absolute error W:** clocksync 0.945 / 0.966, 0.35 h; TimeTeller
  (supervised) 0.754 / 0.950, 0.50 h; timed oracle 0.738 / (F is truth (i) itself), 0.61 h; CHIRAL -0.165 / 0.877,
  1.18 h; CoFE clock genes -0.417 / 0.856, 1.36 h; CoFE genome-wide -0.283 / 0.139, 1.41 h; CYCLOPS -0.296 / 0.115,
  1.83 h. Truth (ii) gives the same ranking (clocksync W 0.831).
- **Donor order (consensus |circular r| with recorded time, best rotation and reflection), W / F:** CYCLOPS 0.47 /
  0.50, TimeTeller 0.46 / 0.25, CoFE clock 0.32 / 0.22, CHIRAL 0.30 / 0.16, CoFE genome-wide 0.06 / 0.09 (clocksync
  donor order in A123: Spearman 0.613 within the window).
- **Reading:** the real-data counterpart of A145. With whole-day sampling, two-step orderers built on clock genes
  (CHIRAL, CoFE clock genes) recover the organ map nearly as well as clocksync (0.86-0.88 vs 0.97), as A145 found on
  idealised data; genome-wide CoFE and CYCLOPS do not. In working hours every orderer loses the map (r <= -0.17),
  while clocksync keeps it (0.945) and beats the supervised TimeTeller and the timed oracle on the window donors. The
  distinctive property of the method on real data is partial-day coverage, not coupling estimation as such.
<!-- src: paper/exploratory/A123_working_hours/output/a123b/{h123_4,tests,donor_order,runtime}.csv -->
- **Post hoc check (2026-11-02, `scripts/posthoc_reflect.py`, `posthoc_reflect.csv`): reflection and stretch allowed.**
  The orderer maps were compared without a global reflection, and orderers spread a 10 h window over the whole circle.
  Taking the better of the two directions and a least-squares stretch, design W offsets vs truth (i): CHIRAL |r| 0.165
  (Spearman 0.127), CoFE clock 0.417 (0.400), CoFE genome-wide 0.283 (0.141), CYCLOPS 0.296 (0.406); clocksync 0.945
  (0.911). Design F: CHIRAL 0.877, CoFE clock 0.856. The working-hours loss of the orderers is not a reflection or
  scaling artefact.
<!-- src: paper/exploratory/A123_working_hours/output/a123b/posthoc_reflect.csv -->

## A149 (scored 2026-10-28, secondary): single-cell tumour / normal pairing exists but is not clock-specific; malignant cells covary with the tumour stroma, specificity not tested
27 cohorts (18 existing, 9 new), 190 paired patients for S1; S2 61 (9 cohorts), S3 51 (7). Nulls equal the power step.
- **H149.1 nominal pass, not clock-specific:** malignant vs own adjacent-normal epithelium 0.534, null 0.446 +- 0.042
  (q95 0.517), p 0.022, Holm 0.044. Registered companion: 200 expression-matched pseudo-clock sets give median 0.579,
  q95 0.767; the clock statistic is below their median. By the registered rule the pairing is not read as
  clock-specific: tumour and normal epithelium of a patient share structure in any matched gene set at least as much as
  in clock genes (shared processing, ambient RNA, genotype or patient state). Contrast with TCGA bulk, where the clock
  exceeded every pseudo-clock set (H138.4).
- **H149.2 FAIL:** after removing the tumour stromal host, 0.550 vs null 0.539 (MDE 0.768), p 0.45.
- **H149.3 PASS (no specificity control registered):** malignant vs tumour stromal host 0.607, null 0.410, p 0.0077,
  Holm 0.023; companion normal epithelium vs the same host 0.666 (null 0.421, p 0.0018). The malignant and normal
  epithelial clocks covary with the tumour stroma of the same patient about equally. Without a pseudo-clock control,
  and with H149.1 showing strong non-clock sharing within patients, this is not evidence of clock coupling.
- **D149.1 / D149.2 (descriptive):** retained share 0.07; per cancer S1 excess lung +0.22 (42 patients, p 0.004),
  colorectal +0.09 (82, p 0.087), gastric -0.02 (21), liver +0.01 (14); S3 lung +0.23 (p 0.045), colorectal +0.13
  (p 0.17); I^2 0 for S1-S3 (4 or 2 cancer types); keeping minus at-null S1 +0.05. **Method flaw:** the patient
  bootstrap intervals of the excess do not contain their point estimates (e.g. S1 0.088 [0.133, 0.540]):
  resampling with replacement duplicates patients and inflates |correlation|. These intervals and the retained-share
  interval are invalid and not reported further.
- **Reading:** the single-cell pairing between a patient's malignant and normal epithelium is real but carried by
  non-clock structure; the clock-specific tumour-normal match remains a bulk TCGA result (H138.4). Whether it runs
  through the stroma cannot be settled with these designs: H149.3 needs a pseudo-clock control before reading.
  Cancer types do not differ detectably; lung shows the largest pairing.
<!-- src: paper/exploratory/A149_paired_tumour_normal/output/score/{tests,pooled,per_cohort,pseudo_companion,d149_1,d149_2_heterogeneity}.csv; output/units.csv; output/power/nulls.csv -->

## A150 (scored 2026-10-28, secondary): Figure 2 atlas coupling is clock-specific overall and for the stromal-vascular and epithelial families; myeloid, lymphocyte and mural are not distinguishable from pseudo-clocks; ambient over-correction lowers lung epithelium
5 gate-passing atlases with genome-wide twins (lung_hlca, heart, organ_breast, organ_muscle, organ_colon), 52
compartments, A17 learned anchors. Reproduction of A17 couplings on the genome tables: exact in 4 atlases; lung_hlca
Spearman 0.998, max difference 0.009 (467 rows differ in depth / cell counts between the genome and clock builds).
- **H150.1 PASS:** mean coupling 0.404 vs 200 pseudo-clock sets mean 0.249, q95 0.316, max 0.381; p 0.005 (floor).
  Per atlas: heart 0.596, muscle 0.606, breast 0.459, lung 0.356, colon 0.203. 45 of 52 compartments have positive
  excess, 13 exceed their own pseudo q95.
- **D150.1 (family median coupling, clock vs pseudo median / q95):** fibroblast-stromal 0.670 vs 0.305 / 0.494;
  endothelial 0.511 vs 0.176 / 0.342; epithelial 0.517 vs 0.267 / 0.420; "other" 0.727 vs 0.410 / 0.653 (all above
  q95). Below q95 ("not clock-specific at this design"): mural 0.613 vs 0.412 / 0.647; myeloid 0.271 vs 0.168 / 0.299;
  lymphocyte 0.094 vs 0.087 / 0.237; parenchyma and neuron (one compartment each, wide q95).
- **D150.2 (ambient over-correction):** median contamination 0.030 (muscle 0.141). Couplings before vs after Spearman
  0.943 (lung 0.849, others >= 0.94). Family median change: epithelial -0.195 (driven by lung: AT1 0.680 -> 0.485,
  basal 0.376 -> 0.123, multiciliated 0.403 -> 0.121, secretory 0.571 -> 0.309; breast and colon epithelium unchanged),
  mural -0.061, neuron -0.059, myeloid -0.056, endothelial -0.045, fibroblast-stromal -0.033, lymphocyte -0.014.
  Family order: epithelial falls from 4th to 7th; the rest keeps its order. By the registered reading ambient RNA is
  not excluded for the epithelial family (and marginally mural / myeloid). Caveat: the correction is an upper bound,
  and lung basal and secretory cells have median contamination 0, so their drop comes through the host (anchors
  corrected), not their own counts.
- **Reading:** the stromal-vascular backbone (fibroblasts, endothelium) is clock-specific and robust to ambient
  over-correction. The low end (myeloid, lymphocytes) is not clock-specific: their small couplings are at the level
  matched non-clock genes reach, which fits "barely coupled" but means their coupling values carry no clock information.
  Lung epithelial coupling is sensitive to ambient correction. Coverage: 5 of 24 gate-passing atlases. Brought to the PI
  before any Figure 2 wording change (registered).
<!-- src: paper/exploratory/A150_atlas_specificity/output/{h150_1,per_family,per_compartment,d150_2,reproduction}.csv -->

## A151 (scored 2026-10-28, secondary): tumour fibroblasts pair strongly with the patient's normal fibroblasts; clock specificity not determinable (correction breakdown in the pseudo sets)
- **H151.1 PASS:** tumour vs own adjacent-normal fibroblasts, 98 patients / 17 cohorts: 0.860 vs null 0.504 +- 0.058
  (MDE 0.650), p 1e-4, Holm 2e-4. Colorectal 0.864 (42, p 2e-4), lung 0.929 (24, p 1e-4). Endothelial companion: 0.984
  vs 0.521 (103 / 15, p 1e-4); liver 0.837, lung 0.842.
- **H151.2 FAIL:** 200 pseudo-clock sets pooled: median 1.013, q95 1.693; p 0.77. **Not interpretable as "non-clock"**:
  the Poisson-corrected |correlation| exceeds 1 in 25 % of per-cohort pseudo values (up to 46 % in a cohort) and in
  several clock values (e.g. 1.22, 1.06), i.e. the variance correction breaks down at fibroblast pseudobulk depths for
  expression-matched genes, and the pooled pseudo distribution is driven by these values. The A149 / A138 estimability
  rule (null mean <= 1) does not catch this for the pseudo sets.
- **Reading:** a strong patient-specific pairing of tumour and normal-tissue stroma exists (fibroblasts and
  endothelium); whether it is clock-specific cannot be decided with this statistic. A bounded comparison (e.g.
  uncorrected |correlation| for clock and pseudo sets alike, or pseudo sets matched on depth-normalised counts and
  Poisson-corrected variance) would be post hoc and needs PI approval. The same breakdown may have inflated the A149
  pseudo-clock values (median 0.579); this is checked before A149's reading is used in text.
<!-- src: paper/exploratory/A151_fibroblast_pairing/output/{tests,pooled,per_cohort,pseudo,power_nulls,units}.csv -->

## A149 note 5 / A151 note 2 (scored 2026-10-28, POST HOC): with the bounded statistic, tumour stroma follows the patient's normal stroma clock-specifically; malignant-normal epithelium pairing stays non-specific
Uncorrected |complex correlation| (cannot exceed 1), clock vs the same 200 pseudo-clock sets per cohort, n-weighted.
- **Fibroblasts (98 patients, 17 cohorts):** clock 0.699 vs pseudo median 0.458, q95 0.542; no pseudo set reaches the
  clock value, p 0.005 (floor).
- **Endothelium (103, 15):** clock 0.701 vs 0.432 / 0.527, p 0.005 (floor).
- **Malignant vs adjacent-normal epithelium (190, 27):** clock 0.433 vs 0.394 / 0.461, p 0.17: not clock-specific,
  as the registered companion (the corrected pseudo values exceeded 1 in 14 % of per-cohort sets, 25-31 % for stroma).
- **Reading (post hoc):** tumour fibroblasts and endothelium keep a clock-specific, patient-specific relation to the
  same patient's normal-tissue stroma, the single-cell stromal counterpart of TCGA H138.4; malignant epithelium shows
  patient-level sharing with normal epithelium that matched non-clock genes reproduce. Methods note: the
  Poisson-corrected |correlation| needs a bound or a breakdown check whenever it is compared with pseudo-clock sets at
  low depth; the registered A151 H151.2 verdict (fail) stands as scored and is superseded in reading by this check.
<!-- src: paper/exploratory/A151_fibroblast_pairing/output/{posthoc_pooled,posthoc_per_cohort}.csv -->

## A152 (scored 2026-10-28, secondary): residency does not carry lymphocyte coupling; the A17 lymphocyte classes include mast cells
A6 residency rule; A17 setup reproduced exactly (max |difference| 8e-17); per-cell Census sums reproduce the atlas
lymphocyte pseudobulks exactly. 1,000 within-donor label swaps; minimum p 0.001.
- **H152.1 FAIL:** heart + breast mean resident-like minus circulating-like coupling -0.092; null 0.004 +- 0.085, q95
  0.148, MDE 0.220; p 0.86.
- **Per tissue (c_res / c_circ, difference, p):** lung 0.091 / 0.077, +0.013, p 0.32 (95 donors with both, MDE 0.085);
  heart 0.043 / 0.106, -0.063, p 0.66 (23 donors, LV nuclei, about 40 cells and 50k UMI per pseudobulk, both at the
  boundary, MDE 0.317); breast 0.427 / 0.549, -0.122, p 0.90 (111 donors). CD8 only: lung 0.029 / 0.000 (p 0.15),
  breast 0.251 / 0.453 (p 0.74). Resident-like share per donor 0.48-0.52. All split tables pass the gate. Sensitivities
  (both-state donors only; tertiles on analysed cells) keep the signs.
- **D152.0 (lung reconciliation):** A17 lung:T is 0.153 on the clock table (the 0.159 quoted in the registration is the
  A150 genome-twin value); T as one class from the state builder 0.142-0.145; the split 0.091 / 0.077. Under the current
  estimator lung T cells are weakly coupled whichever way they are cut; the archive A6 values (0.61 / 1.00) do not
  reproduce.
- **Finding (class mapping):** the A17 class regex assigns Census "mast cell" to `heart:lymphocyte` and
  `breast:T_other` ("t cell" inside "mast cell"). Heart atlas: 2,927 mast cells against 16,636 lymphocytes and 4,302
  NK cells. The A17 heart lymphocyte coupling (0.556) therefore mixes mast cells (tissue-resident) with lymphocytes; the
  residency split, which excludes mast cells, gives 0.04-0.11 in heart. Mast-cell diagnostic running (post hoc).
- **Reading:** residency does not confer coupling (consistent with archive A6 in healthy lung); in breast,
  circulating-like T cells couple at least as much as resident-like ones. Heart is uninformative at its depth.
<!-- src: paper/exploratory/A152_lymphocyte_residency/output/{h1521,perm_summary,states_summary,fits_observed,lung_repro_check,part2_repro_check,repro_a17,gates}.csv -->

**A152 mast diagnostic (2026-10-28, POST HOC, descriptive)** (`A152_lymphocyte_residency/output/mast_diagnostic.csv`):
heart:lymphocyte A17 0.556 (38 donors; mast share median 0.28) -> without mast cells 0.328 (36, boundary), heart:mast
alone 0.509 (26 donors). breast:T_other 0.089 (mast share 0.09) -> 0.070; breast:mast 9 donors, not estimable.
**Class-mapping audit (2026-10-28, lead):** `_atlas/organ.py` COARSE is matched case-insensitively (`census.classify`,
default) with `T_other = "T cell|thymocyte"` before the epithelial, mast and fibroblast patterns, so Census labels
containing "...t cell" are classed as T_other: mast cells, goblet cells, tuft cells, myofibroblast cells and "alveolar
type 1 fibroblast cell". `heart.py` has the same issue (`lymphocyte|T cell|natural killer`, mast cells). Affected
builders: organs.py (organ_* atlases), kong.py, genome.py, immune_aging.py, heldout.py, heart.py, and possibly
cue_units (own COARSE copy; to check). Not affected: lung_hlca / Natri (author annotations), Tabula Sapiens
(case-sensitive). Examples: colon_crc goblet 16,142 + tuft 2,845 + mast 3,834 cells in colon:T_other; Kong colon and
ileum epithelial goblet / tuft and stromal myofibroblasts in T_other. The archive A17 run has the same mapping (A17
reproduces it exactly). Fix and reruns brought to the PI.

## A123c (scored 2026-10-28): against a truth from out-of-window donors only, H123.5 fails (offsets 0.679), H123.6 passes (amplitudes 0.965)

Truth (iii): cosinor on recorded time of the fast deaths who died outside 08-18 h, complex SVD as truth (i).
- 133 truth donors and 111 test donors, disjoint (asserted).
- 39 of 40 A123 tissues have >= 20 truth donors; Kidney - Cortex (19) is dropped.
- The A123 clocksync estimate is regenerated in memory and matches `organ_map.csv` to 4 decimals.

**Hypotheses:**
- **H123.5 FAIL.** Offsets vs truth (iii): circular correlation 0.679 [0.467, 0.839], tissue-permutation p 1e-4 (the
  floor), median error 0.58 h. The point estimate is below the registered 0.7; the p criterion alone passes.
- **H123.6 PASS.** log2 relative amplitudes vs truth (iii): Pearson 0.965 [0.948, 0.980], p 1e-4.

**Descriptive** (offset r [interval]; amplitude Pearson):
- **The truths agree imperfectly:** truth (i) vs truth (iii) 0.821 [0.657, 0.925]; amplitude 0.978. Truth (iii) rests
  on night deaths only (a 14 h arc, 133 donors).
- **Other methods vs truth (iii):**
  - the timed oracle on the window donors' recorded times: 0.491 [0.168, 0.736] (p 8e-4);
  - CHIRAL 0.023;
  - CoFE genome-wide -0.109;
  - CoFE clock -0.326 (38 tissues);
  - CYCLOPS -0.227.
  clocksync stays well above every untimed orderer, and above the timed oracle of the same donors.
- **Split sensitivity:** donors split once, seed 123, 65 / 68. The waveform comes from half 1, clocksync is refitted, and
  the truth comes from half 2.
  - At >= 20 half-2 donors (35 tissues): 0.776 [0.632, 0.879]. The stored estimate gives 0.775.
  - At >= 10 (the same 39 tissues as H123.5): 0.804 [0.684, 0.896]. The stored estimate gives 0.807.
  - Unexplained, not yet checked: the same estimate agrees better with the half-2 truth (68 donors) than with the full
    truth (iii) (133 donors).
- **Prior-only baseline, per-donor hour on the 111 window donors** (median absolute error):
  - clocksync 1.73 h with the error-minimising common offset, 2.09 h with the circular mean set to 13:00 (uses no truth);
  - always 13:00: 2.62 h;
  - best constant: 2.37 h (it uses the truth).

**Choices beyond the registration** (README):
- The truth (iii) DESeq2 size factors come from out-of-window samples only.
- `compare` uses its hard-coded 1,000 tissue bootstraps.
- The split is by donor, greedy within tissue.
- The >= 10 half-2 variant and the best-constant baseline were added.

**Reading (PI to decide):**
- Amplitudes hold against a disjoint truth.
- Offsets fall just below the registered threshold, with a wide interval.
- The ceiling set by truth (iii)'s own noise is unknown, since truth (i) vs truth (iii) is only 0.821.
- This removes or reframes a main-figure claim (Figure 1 organ map), so it is raised with the PI before any wording
  changes.
<!-- src: paper/exploratory/A123_working_hours/output/a123c/{tests.csv,truth_iii.csv,README.md}; logs/compute_truth_iii.log -->
- **PI decision (2026-10-28):** (a) the Figure 1 panel is kept with the verdict as scored. It shows amplitudes confirmed
  against the disjoint truth, offsets 0.68 (below the registered 0.7) but above every untimed orderer and the timed
  oracle, and both truths. (b) A post hoc descriptive diagnostic is added: per-tissue residuals, the bootstrap noise
  ceiling of truth (iii), and the split-half distribution (`output/a123c/diag/`). Class-fix reruns are held after the
  anchor relearn until roadmap phases 1-2 are done.

### A123c post hoc (2026-11-02, descriptive): the estimate is far more reproducible than the disjoint truth

`scripts/posthoc_estimate_split.py`: the 111 window donors split into random halves balanced within tissue (40 splits),
clocksync refitted on each half with the A123 cross-fitted waveform.
- **Half map vs half map:** offset circular r median 0.873 [0.735, 0.942]. The truth (iii) halves agree at 0.544
  [0.317, 0.714] (`diag/split_half.csv`, 200 splits).
- **Each half (about 55 donors) vs truth (i):** 0.900 and 0.915; vs truth (iii): 0.637 and 0.669.
- Spearman-Brown reliability at full n (linear formula applied to circular r, descriptive): estimate 0.93, truth (iii)
  0.70; r(estimate, truth iii) corrected for both 0.681 / sqrt(0.93 x 0.70) = 0.84.
- The full map refitted here gives 0.945 vs truth (i) and 0.681 vs truth (iii) (stored H123.5 0.679; the difference
  comes from truth (iii) read back from its 4-decimal csv).
- **Reading:** the working-hours map reproduces across disjoint donor halves (0.87 at half the donors), so the H123.5
  shortfall comes from the night-only truth (reliability 0.54 at half size), not from the estimate.
<!-- src: paper/exploratory/A123_working_hours/output/a123c/diag/{estimate_split_half.csv,split_half.csv} -->

## A17 note (2026-10-28): anchors relearned on the class-fixed tables (provisional until roadmap phase 2)

`compute_coherence.py` then `compute_anchors.py` were rerun on the class-fixed tables, with the unchanged rule.
- **New compartments in the learning:** the 28 dataset/compartment pairs the fix created
  (`common/class_fix_new_compartments.csv`), through `Target(new_classes=True)`. A first attempt with every table
  compartment was discarded: it added compartments the original A17 excluded on purpose in 22 targets.
- **Check:** on the pre-fix coherence outputs, the code reproduces the frozen table 70/70 for both rules.
- **Changes:** 12 of 70 targets changed, and 3 of 301 target x family rows changed learned status. Learned flag changes,
  verified against the two tables:
  - `ia_lung`: the lymphocyte family is no longer learned (B, CD4_T, CD8_T, ILC, NK and T_other removed), and ia:mast is
    added (myeloid family).
  - `organ_colon`: colon:endothelial is removed and colon:macrophage_DC is added. The endothelial family goes from
    learned to not learned, and the myeloid family from not learned to learned.
  - New learned anchors: heart:adipocyte in dcm and dcm_healthy.
  - kong_colon / kong_ileum T_other keep their learned flag after the relearn.
- **Fixed rule:** unchanged everywhere.
- **Tabula Sapiens:** rebuilt under the current rules, with no change in rows or cells.
- **Cancer anchor rows** (214) are carried over unchanged: their tables do not use the atlas class rules.
- **Status:** no fit has been rerun (PI hold). The relearn is repeated after the roadmap phase 2 legacy-table audit if
  any table changes.
<!-- src: paper/common/anchors_a17.csv vs anchors_a17_pre_classfix.csv; paper/results/logs/anchor_relearn_diff.md; paper/results/extended_data/output/host_identity/a17_anchor_families.csv -->

### 2026-10-28, A123c post hoc diagnostic (descriptive; the H123.5 verdict stands)

Truth (iii) is noisy. The estimate's gap to it is mostly truth noise plus one tissue.
- **Noise ceiling:**
  - Split-half reliability of truth (iii) (200 donor splits): r between half truths, median 0.544 [0.317, 0.714].
    Spearman-Brown full-truth reliability 0.695, which gives a ceiling of about 0.83 for a perfect estimate. A
    bootstrap proxy gives 0.93.
  - The estimate's r disattenuated for truth noise only is 0.73-0.81. Truth (i) vs truth (iii) disattenuated is
    0.88-0.99.
  - These are descriptive estimates (linear attenuation applied to circular r).
- **The favourable split:** the registered split result (0.807) is a favourable draw. Over 400 random half truths,
  r(estimate, half) has median 0.596 and only 5 % are >= 0.807. Construction checks found no bug: the reference tissue,
  centring set and size factors are immaterial.
- **Whole Blood drives the gap.** Estimate minus truth (iii) is +3.78 h, while estimate minus truth (i) is +0.52 h.
  Whole Blood's truth (iii) cosinor R2 is the weakest of 39 (0.028, median 0.156), and leaving it out gives r 0.754.
  Next: Testis (+2.73 h), Stomach (-2.46 h, 21 donors) and Breast (-2.24 h).
- **Reading:** against the disjoint truth, the offsets sit near what truth (iii) noise allows. The amplitudes (0.965)
  are robust. For the Figure 1 legend: show both truths, the noise ceiling and the Whole Blood outlier.
<!-- src: paper/exploratory/A123_working_hours/output/a123c/diag/{bootstrap_truth,split_half,residuals,construction_checks}.csv -->

### 2026-10-28, roadmap phases 1-2: code and data integrity (no result changes)

- **clocksync (commit e89defd in clocksync/):** every addition is opt-in.
  - flagged and uncorrected ccorr;
  - responsive-gene pseudo sets and `drop_genes`;
  - debiased z2 / coupling with `boundary_type`;
  - per-arm host gate in `contrast`;
  - bootstrap failure counts;
  - Hartung-Knapp `meta`;
  - lymphatic excluded from `suggest_anchors`;
  - below-paper-setting warnings.
  A regression check shows the defaults bit-identical on gate, fit, contrast, coherence, paired, ccorr, meta and
  pseudo_waveforms. 47 tests pass.
- **Rule-version stamp:** stored tables carry it (`data.save`). `load` refuses stale tables and warns on legacy
  fallbacks.
- **Legacy audit:** 27 legacy-loaded A17 datasets. Only kidney_sc carried the "T cell" trap: 103 mast cells went to
  kidney:T (0.9 %). It is rebuilt.
- **Other rule fixes:** a COARSE substring trap was fixed: "basophil" matched "basophilic erythroblast". No compartment
  changes, because those cells do not reach the table thresholds.
- **Rebuild and stamp:** all 48 A17 tables were rebuilt under the current rules and stamped. No compartment changed
  beyond the fixes above (`results/logs/rebuild_and_stamp.md`), so the relearned anchors stand.
- **Tabula Sapiens:** unchanged under the new rules (the mast / goblet / tuft labels never reached 5 donors).
- **A22 stratum artefact:** besides V17.2, which was already recomputed, no cited coherence was uncentred.
- **Open, PI:** COARSE ordering choices remain: "granulocyte monocyte progenitor cell" goes to monocyte, and
  "mesenchymal stem cell" and "satellite stem cell" go to progenitor (Q14).

### 2026-10-28, A147 post hoc: reference group at the real fast-death size (244)

Same simulator and seeds, with the reference group drawn at n 244 instead of the group's n.
- **The shrink ratio at f = 0 is no longer biased:** n 79 gives 0.892, inside [0.85, 1.15] (registered run 0.809). n 298
  gives 0.911.
- **The order criterion V147.1b now misses at n 79:** 0.755 vs >= 0.80 (registered run 0.805, a pass). AUC 0.985 still
  passes.
- **Overall:** still FAIL (V147.1 fails at n 79).
- **Reading:** both n-79 criteria sit near their thresholds, so the run-to-run Monte Carlo noise at 200 replicates
  decides them. The readout works at n 298 under either reference size. Extended Data, as already decided, with this
  sensitivity stated.
<!-- src: paper/exploratory/A147_pull_order/output/posthoc_nref244/v147.csv; output/v147.csv -->

## A124 (scored 2026-10-28): simulated night shift in human muscle, misalignment not detected (H124.1, H124.2 fail; MDE 7.3 h)

GSE106800: 11 complete subjects, 44 arrays, 2,048 exact within-subject label swaps (smallest p 4.9e-4). Gene coverage
with `clocksync.genes.canonical`: core 14/14, cue 35/40, locked 36/40.
- **Data-rule error (lead):** the data rule expected 40 / 39, but the canonical table lacks 8 old platform names
  (FAM46B, MTERFD3, HIST2H2BE, C12orf23, SEPT4, C10orf116, CENPC1, FAM134B). The primary follows the rule as written
  (canonical); the full-coverage run is a post hoc sensitivity.

**Hypotheses:**
- **H124.1 FAIL.** |Δ core-cue| = 3.35 h (misaligned minus control -3.35 h, bootstrap [-7.45, 2.01]), p 0.211, Holm
  0.422.
- **H124.2 FAIL.** |Δ core-cue| − |Δ core-locked| = -0.43 h, p 0.771. Δ core-locked -3.78 h [-8.57, -0.60].
- **MDE:** 7.3 h (null circular SD 2.61 h; null q95 of |Δ| 4.76 h).

**Descriptive:**
- Coherence falls under misalignment: core-cue 0.55 -> 0.40, core-locked 0.71 -> 0.34.
- **The muscle core clock does not follow the shifted behaviour within the protocol.** The core state moves between
  07:00 and 19:00 in the same direction in both conditions (cosine 0.91), with about half the amplitude under
  misalignment (3.20 -> 1.67). The cue state is too noisy for a direction.
- With all 12 subjects: Δ core-cue -3.37 h.
- **Post hoc, full gene coverage** (8 platform aliases added; core 14, cue 40, locked 39): H124.1 |Δ| 3.12 h, p 0.140.
  H124.2 3.01 h, p 0.054, Holm 0.107. Δ core-locked -0.11 h. MDE 5.9 h. Neither passes. H124.2 depends strongly on the
  mapping of 8 genes.

**Reading:** not detected at n = 11 with two biopsy times per condition (MDE 7.3 h). The design cannot resolve shifts of
the size expected. Under the registration, A126-A127 stay descriptive. The halved evening-morning core amplitude is
consistent with a damped, not shifted, muscle clock under acute misalignment (descriptive).

**Tool gap:** `clocksync.genes.canonical` lacks common old platform names. Adding them is roadmap 8.4; it changes gene
coverage for array datasets and needs a regression check.
<!-- src: paper/exploratory/A124_night_shift/output/{tests,mde,descriptive_deltas,descriptive_time_state_difference,sensitivity_12_subjects}.csv; output/alias_sensitivity/tests.csv -->

## A156 (scored 2026-10-28): forced desynchrony in human blood, core-cue change not detected (H156.1, H156.2 fail; MDE 13.1 h)

GSE48113: 19 subjects with >= 5 arrays per day, 259 arrays, 10,000 within-subject day swaps. Genes mapped with
canonical plus the 8 A124 aliases (core 14/14, cue 39/40, locked 38/40).

**Hypotheses:**
- **H156.1 FAIL.** |circular mean Δ core-cue| 4.80 h (D4 minus D1 -4.80 h, bootstrap [-13.5, 4.4]), p 0.365, Holm
  0.731.
- **H156.2 FAIL.** |mean Δ core-cue| − |mean Δ core-locked| = 3.04 h, p 0.468.
- **MDE:** 13.1 h (null circular SD 4.69 h).

**Sensitivities** (none pass):
- 7 + 7 subjects (exact): 4.96 h, p 0.70.
- All 21 subjects: 2.40 h, p 0.49.
- Canonical names only: 0.37 h, p 0.97. The untimed offset moves by 4.4 h with the mapping of 5 cue genes.

**Timed descriptive (cosinor on clock time, D4 minus D1):**
- The blood core state stays near the central clock: -1.36 h [-4.7, 5.7], resultant 0.27.
- The cue state shifts with the displaced behaviour: -6.85 h [-11.1, -2.8].
- Per-subject coherence is low (|mean z| core-cue 0.19 on D1, 0.21 on D4).
- The untimed Δ core-cue does not track the timed core-minus-cue change: circular r 0.04 [-0.57, 0.64].

**Imputation:** 1,708 gene-array cells were imputed within the subject-day. In 40 subject-day-gene cases the gene was
missing in every array; these fell back to the same-day mean over subjects, a step not in the data rule (stated).

**Reading:** not detected. At 7 blood samples per subject-day the per-subject core-cue offset is too noisy (MDE 13 h),
and it depends on a handful of cue genes. The timed analysis shows the expected biology: cue genes follow behaviour and
the core clock stays. The untimed per-subject core-cue readout cannot resolve it in blood arrays. With A124, two
registered human perturbation tests are now uninformative for power. The core-cue estimator is reliable only where
clock and cue signals are strong (psoriasis skin, GTEx), which is stated as a limitation.
<!-- src: paper/exploratory/A156_forced_desynchrony/output/{tests,mde,sensitivities,descriptive_deltas,descriptive_timed,coherence_per_day}.csv -->

### 2026-10-28, roadmap 6.3 (GTEx long-term steroid users): infeasible, not registered

Aggregate counts of `MHSTRDLT` = yes by Hardy: H0 8, H1 0, H2 6, H3 4, H4 6. That is 24 users in total and 6 among fast
deaths. Below any matched design (>= 15 per arm, as in A132 D132.3). Dropped; no clock readout was done.

### 2026-10-28, psoriasis delay: gene-set robustness (post hoc, descriptive; motivated by A156)

The headlines reproduce exactly (15 gated series): core -1.426 h, cue +0.263 h, core minus cue -1.689 h; A54 held-out
core -0.932 h.
- **Leave-one-gene-out:** every estimate stays negative.
  - Core: range -1.53 to -0.89 h. Only NR1D1 shifts it by more than 0.5 h (+0.54 h; +0.31 h on the held-out series).
  - Cue genes: core minus cue ranges -1.76 to -1.56 h (largest shift 0.13 h).
- **Random 80 % subsets** (500 each): core minus cue [-1.89, -1.43] (cue subsets) and [-1.86, -0.98] (core subsets),
  100 % negative.
- **Pseudo-gene controls** (200 sets each): pseudo offsets are centred on 0. Responsive-gene sets on the core waveform
  give [-0.58, 0.40], matched sets [-0.61, 0.59]. 0 of 200 sets in any family reach the core offset. The cue offset
  (+0.26) sits inside the pseudo range.
- **Gene names:** the A124 alias names are not core or cue genes in these series (coverage 14 / 14 core, 36-37 / 40
  cue).
- **Reading:** unlike in blood (A156), the psoriasis core delay is robust to the gene sets and specific against
  responsive and matched non-clock genes. The delay sits in the core clock, not in the cue plane. One gene (NR1D1)
  carries about a third of the core shift; this is stated.
<!-- src: paper/exploratory/psoriasis_robustness_2026-10-28/output/{reproduction,leave_one_gene_out*,random_subsets_summary,pseudo_sets_summary,platform_coverage}.csv -->

## A109 / A138 note scored (2026-10-28): H109.1-R and H138.4-R pass; the shared states are not a GC or immediate-early response

The registered statistics reproduce to 1e-9 (H109.1) and exactly (H138.4). Responsive-gene sets: 200 per test, 28 of 29
pool genes present.
- **H109.1-R PASS** (Holm 0.0099). Ventilated C_clock 0.463 vs responsive median 0.214, q95 0.257, p 0.005 (the floor).
  Fast deaths 0.609 vs 0.197; Hardy 3 0.449 vs 0.275; Hardy 4 0.413 vs 0.198, each p 0.005.
- **H138.4-R PASS** (Holm 0.0099). Pairing 0.365 vs responsive median 0.252, q95 0.310, p 0.005. Per cancer, the
  responsive median exceeds the pairing in THCA and UCEC, and the responsive q95 exceeds it in COAD, LUAD and PRAD.
- **PER1/PER2-free:** unchanged. Ventilated C 0.463 (vs matched q95 0.171); H138.4 pairing 0.367 (vs matched q95
  0.261; permutation p 1e-4).
- **Immediate-early control (H138.4):** the pooled IEG-set pairing is 0.232, not distinguishable from the matched null
  (p 0.19). Shared surgery stress pairs the two samples only weakly.
- **Stratified pairing nulls** (10,000):
  - within cancer x TSS: null mean 0.129, q95 0.166;
  - within cancer x tumour plate: mean 0.122;
  - within cancer x normal plate: mean 0.127;
  - observed 0.365, p 1e-4 in each.
  Batch contributes at most +0.008 to the null mean. Per cancer, COAD (p 0.11-0.30), KIRP (0.08-0.10) and UCEC
  (0.82-0.94) do not pass under stratification.
- **Caveat:** the responsive sets overlap heavily (Jaccard 0.79 in GTEx, 0.91 in TCGA; KIRC and LIHC 1.0). The responsive
  null is narrow in gene content.
- **Reading:** the review's main alternative for Figure 3 and Figure 5 (a shared GC / acute-stress pulse, or batch) is
  tested and does not explain the shared clock-gene state or the tumour-normal pairing. H109.1 can be worded "a
  clock-gene-shaped shared state, not a GC or immediate-early response"; "running clocks" still rests on H109.2.
<!-- src: paper/exploratory/A109_coupling_nulls/output/responsive/{h109_1_responsive,reproduction_check,jaccard_responsive_sets}.csv; paper/exploratory/A138_cancer_estimability/output/responsive/{h138_4_responsive,holm_responsive,stratified_null,per_cancer}.csv -->

## A157, A158 and the V25.5 prior baseline (scored 2026-10-28)

**A157: V157.1 FAILS at the boundary.** Errors were calibrated first: median error 1.31 h for muscle, 1.11 h for adipose.
500 replicates per design.
- **Muscle:** JS median 0.33 [5-95 %: -0.83, 0.87], JS < 0 in 37.2 %. Fisher-Lee 0.71 [0.59, 0.80], > 0.5 in 99.6 %.
  The muscle criteria pass.
- **Adipose:** JS < 0 in exactly 5.0 % (25 of 500), against the registered "< 5 %". This part fails, so V157.1 fails.
- **Reading, as registered:** muscle goes to Extended Data as a registered fail. The simulation shows the mechanism:
  JS is unstable on the evenly spread muscle grid, and the observed -0.854 lies inside its range.
- **Choices beyond the registration:** sample noise SD 1 h; the muscle waveform is GTEx without ARNTL, since the 13-gene
  refit needs time of death.

**A158: the gap is not estimator bias.**
- **Registration correction (lead):** A129 / D145.2 are on Manella 2021 (hepatocyte Bmal1 KO, timed). The registration
  wrongly named `jouffe2022_bmal1ko`, a whole-body KO sampled at one ZT, so the Manella design was simulated.
- **Bias:** both estimators shrink the true change by about 9 % (slopes 0.91 and 0.92, bias +0.02 to +0.05). In
  simulation clocksync's change is <= the oracle's (gap <= 0.014).
- **The observed values** correspond to true changes of -0.32 to -0.33 (clocksync) and -0.20 (oracle).
- **Anchor and truncation:** the liver as anchor has no effect (coupling is read on the leave-one-out host), and no KO
  fit reaches the z2 floor.
- **Reading:** the larger real-data clocksync drop (-0.293 vs -0.186) comes from the data, not from estimator bias. For
  example, the control coupling differs: clocksync 0.969 vs oracle 0.934. Methods sentence: "both estimators are
  conservative by about 9 % in simulation".

**V25.5 prior-only baseline (42 patients, median absolute error):**

| predictor | error (h) |
|---|---|
| method, error-minimising offset (as published) | 1.07 |
| method, circular mean set to 13:00 (no truth) | 1.42 |
| always 13:00 | 2.10 |
| best constant (uses truth) | 1.67 |

The recorded-hour SD is 2.17 h. The method beats even the oracle-favoured constant.
<!-- src: paper/exploratory/A157_js_grid/output/a157_summary.csv; A158_bmal1_bias/output/{a158_bias,a158_reproduce}.csv; A25_prior_baseline/output/prior_baseline.csv -->

## A155 (scored 2026-10-28): reading anchors on a host that contains them inflates anchor amplitude at 3-4 anchors (decision rule met); offsets and the A123 map unaffected

200 replicates per cell. Readout (a) is the current full host; readout (b) is the leave-one-out host for anchors.
- **Amplitude:** the anchor amplitude bias under (a) (mean absolute, % of |z|) depends on the number of anchors:

  | anchors | bulk | single cell |
  |---|---|---|
  | 3 | 13.0-13.6 % | 18.2-19.6 % |
  | 4 | 9.0-9.1 % | 12.1-14.5 % |
  | 6 | 3.6-3.9 % | 4.6-5.8 % |
  | 10 | 1.5-1.9 % | 1.7-2.5 % |

- **Offset:** bias <= 0.06 h everywhere.
- **A123-like design** (40 anchors, n 111): offset r 0.9824 under both readouts (change -3e-6); amplitude r 0.9927 (change
  -1e-6).
- **Decision rule (registered):** met, since amplitude bias exceeds 5 % at typical designs (m 3-6). Leave-one-out becomes
  the anchor amplitude and offset readout in clocksync, pending PI confirmation (Q13, 0.4).
  - Coupling is already read on the leave-one-out host, so coupling results are unaffected.
  - Affected: anchor relative amplitudes (Figure 2 amplitude panels, anchor-vs-target amplitude comparisons).
  - A123 needs no rerun.
<!-- src: paper/exploratory/A155_anchor_self_inclusion/output/{decision,summary_a123like,summary_bias}.csv -->

## Phase 3 reruns on the relearned anchors, and A159 (scored 2026-10-29): no registered verdict changes; A159 H159.1-2 not scorable, D159.1 descriptive

**Reruns** (17 jobs, all ok; comparison in `results/logs/phase3_compare.md`). Relearned anchors, `anchor_readout="loo"`,
family-summary gates seeded by target name.
- **Anchors:** the confirmation relearn differs from the 2026-10-28 table only in organ_bone_marrow (gate fails; no fit
  uses it). The basophil COARSE fix moved basophilic erythroblasts out of `bone_marrow:granulocyte`, a change
  `rebuild_and_stamp.md` missed. Table: 736 relearned plus 214 cancer rows; previous table in
  `common/anchors_a17_pre_phase3.csv`; diff in `results/logs/anchor_relearn_diff.md`.
- **Registered or pre-specified verdicts:** none changed.
  - A85 H85 still passes: 0.087 to 0.097 [0.067, 0.143], 71 classes.
  - A122 H122.1-2 still fail (Holm 0.645 / 0.507).
  - A152 H152.1 identical (−0.092, p 0.856). A144 D144.2 identical.
  - A22 V17.1 dropped share 0.0159 (pass); A22 anchor drops 5.
  - A17 gate-passing healthy atlases: 24 of 34, the same set.
  - The A17 colon host is now stable (0.795 to 0.915).
  - WP0 anchor robustness gains sign flips only in ia_crohn (new compartments).
- **Family medians, learned host** (`figure_2/output/backbone/family_medians.csv`, variant all):
  - lymphocyte 0.089 to 0.125;
  - other 0.410 to 0.488;
  - endothelial 0.549 to 0.526;
  - parenchyma 0.545 to 0.565.
  Fixed and free medians moved by at most 0.04. Stroma, vascular and neurons remain above myeloid and lymphocyte under
  every host; lymphocyte remains last. Ten targets are dropped for failed gates (`family_dropped_targets.csv`).
- **Anchor readout loo vs full** (A17 fits): median anchor amplitude ratio 0.967 (189 anchor rows). Median offset change
  0.061 h for anchors and 0.006 h for others. organ_colon moves by up to 10 h because its anchor set changed.
- **A152:** heart lymphocyte 0.556 (38 donors) to 0.327 (36); `heart:mast` 0.509 (26 donors). The mast diagnostic is now
  vacuous (no mast cells left in the lymphocyte rows).
- **Not rerun** (read the anchor table; long or out of scope; listed in `phase3_compare.md` section 4): the Figure 4
  contrasts (A153), A150, gene rhythms, the cancer anchored scripts and others.

**A159:**
- **Coverage** (`exploratory/A159_specific_coupling/output/row_check.csv`):
  - The exact row check passes for heart, organ_breast and organ_colon.
  - It fails for lung_hlca (467 of 1,321 rows with count differences) and organ_muscle (1 row, 5 counts).
  - organ_adipose fails its gate.
  - No Figure 4 control arm has a genome-wide table on disk.
- **H159.1-2:** not scorable, because all 3 covered atlases are A150 atlases (`h159.csv`). No Holm was run.
- **D159.1 (descriptive; 31 compartments, 200 matched and 200 responsive sets each, 0 failed fits):**
  - Mean c_clock 0.387 vs mean matched median 0.227; median Δc 0.132.
  - Clock-specific: 3 of 31, all in breast — endothelial 0.763, epithelial 0.786, plasma 0.810.
  - Family median Δc (`d159_families.csv`, covered): endothelial 0.306, fibroblast_stromal 0.233, epithelial 0.192,
    mural 0.125, myeloid 0.066, lymphocyte 0.015.
  - The per-compartment pseudo q95 is 0.37-0.92 (`d159_compartments.csv`). The registered flag therefore admits only
    compartments with c above about 0.75: heart fibroblast fails at 0.718 (matched q95 0.912).
- **Consequence:** A153 has no primaries as registered: no control-arm tables, and a flag that is underpowered per
  compartment. A153 is parked until the PI decides (PI_QUESTIONS Q15).

## A160 (scored 2026-10-29): psoriatic delay within keratinocytes — not testable at this design (both MDEs > 1.7 h), descriptive only

Data: GSE173706 (on disk). All-gene pseudobulks were rebuilt from the per-library count matrices using the archived cell
annotations; the clock-gene counts match the archive. Code and outputs: `exploratory/A160_psoriasis_keratinocyte/`.
- **Feasibility** (`output/a160_feasibility.csv`):
  - 10 PN-PP pairs; patient 31170 drops at 132 PN keratinocytes, and 30696 V1/V3 count as two pairs.
  - 5 of 8 healthy libraries pass.
  - The keratinocyte variant uses 34 of 40 A81b cue genes (4 absent as old aliases).
  - Eligible bins: cycling, d2, d3, d4 (8-9 pairs); d1 has 7 pairs and is excluded.
- **Blind power** (PN plus healthy only, 2,000 draws; `a160_blind_power.csv`):
  - Keratinocyte null SD 2.17 h (boot) / 2.00 h (split), MDE 5.39 h.
  - Bins: 1.40 / 1.56 h, MDE 3.88 h.
  - Neither variant meets the 1.7 h gate, so H160.1-2 have no verdict.
- **Descriptive readout** (`a160_tests.csv`, `a160_descriptives.csv`, exact swaps 2^10 = 1,024, one-sided p):
  - Keratinocyte PP − PN core-cue −2.82 h (p 0.014). Without NR1D1 −1.57 h (p 0.22); core alone −2.61 h (p 0.32).
  - Bins −2.83 h (p 0.10). Without NR1D1 +0.91 h. Per bin (`a160_per_bin.csv`): −7.3, +1.0, −1.8, −5.0 h, bootstrap SE
    4.9-6.0.
  - Whole skin on the same 10 pairs: −1.35 h (p 0.058). With all 11 pairs: +0.73 h.
  - PN − healthy: keratinocyte +3.36 h [−2.72, 8.16]; whole skin −0.68 h [−6.91, 4.96].
- **Deviations from the registration:**
  - 10 pairs instead of 11, by the data rule.
  - The pseudo-pair null resamples, because the pool (15 libraries) is smaller than 2n. The gate uses the larger MDE of
    the bootstrap and split schemes (the A115 precedent).
  - The estimand is read as A115's within-group offset (one projector per set over both groups). A first blind run with
    a between-sample paired correlation (null SD 7.6 h, PN/healthy only) was discarded before any PP statistic was
    computed.
- **Reading:** the single-cell series is too small to separate keratinocyte-intrinsic delay from composition. The
  point estimate is negative and leans on NR1D1, as in bulk. Even the whole-skin pseudobulk of these libraries does not
  stably reproduce the bulk delay: −1.35 h, and +0.73 h with one more pair. Roadmap 6.9 composition control stays open;
  it needs a larger paired single-cell cohort (the 2026 AD/psoriasis atlases in the scout, accession unverified) or
  bulk deconvolution.

## Phase 5 bounded re-scoring (2026-10-29; roadmap 5.3-5.5; descriptive / post hoc, registered verdicts unchanged)

Outputs are in `exploratory/<amendment>/bounded_2026-10-29/output/`. A149, A151 and H138.5 use the bounded
`stats.ccorr(nan_on_flag)`. A142/A143 r* is a Spearman-Brown disattenuated correlation, not a ccorr, so it gets a
reliability analogue of the flag: NaN when reliability ≤ 2 SE.

- **A142** (`targets_bounded.csv`):
  - Brain r* 0.639, reference q05 0.936, p 0.001, reproduced and not flagged.
  - Breakdown share: 0 of 1000 reference blocks at k = 7 and 13; 7.2 % at k = 4.
  - On the bounded uncorrected r: brain 0.575 vs reference q05 0.802, p 0.001. Digestive 0.782, p 0.207, so the
    digestive departure (r* p 0.029) does not hold on r.
  - The review's r* > 1 in 19 % of k = 13 reference blocks is sampling scatter around a reference median of 0.978.
  - "Brain departs most" stands.
- **A143** (`h143_1_bounded.csv`, `h143_3_*_bounded.csv`):
  - H143.1 adjusted r* 0.627, unchanged.
  - H143.3 on bounded alternatives (share a_P/(a_P+a_B), r_P, uncorrected ratio): the hypothalamus is 10th of 13 in fast
    deaths (11th on r_P). Hypothalamus minus the median of the other regions: −0.030 on share (p 0.74). The review's
    "6th of 13" is the descriptive Hardy 3-4 column.
  - The hypothalamus stays Discussion-only, with no "follows the body" wording.
- **A134 H134.1** (`h134_1_bootstrap.csv`):
  - The null permutes whole-patient copy-number covariate rows, malignant and normal independently, with expression
    fixed. It tests association against unrelated covariates, not the precision of the gain.
  - Bootstrap SD 0.102 vs null SD 0.039; a patient-resampled permutation null has SD 0.043 (p 0.035). The wide bootstrap
    therefore reflects uneven patient leverage of the real association.
  - Bootstrap p: share ≤ 0 is 0.297, bootstrap-t 0.226, centred 0.118; CI [−0.160, 0.266].
  - Wording: report the CI with the p. H134.1 stays weak.
- **A149** (`tests_bounded.csv`, `pseudo_companion_bounded.csv`; flagged units S1 10/27, S2 6/9, S3 4/9):
  - H149.1 0.534 → 0.490 (17 cohorts, 131 patients), p 0.011, Holm 0.021: still a pass.
  - H149.2 → 0.213 on 3 cohorts and 28 patients: not informative.
  - H149.3 0.607 → 0.707 (5 cohorts), p 0.0004, Holm 0.0012: still a pass and stronger.
  - Pseudo companion: clock 0.490 vs pseudo median 0.445, q95 0.602, p 0.338. Still not clock-specific; the old
    unbounded pseudo q95 was 4.5e4.
  - D149.1 retained share is meaningless (−2.97).
- **A151** (`tests_bounded.csv`; flagged fibroblast 5/17, endothelial 5/15):
  - H151.1 0.856 (12 cohorts, 65 patients), p 1e-4.
  - H151.2 bounded pseudo median 0.664, q95 0.906, p 0.141. Not readable: 79 % of fibroblast cohort × pseudo-set values
    are flagged, so the surviving null is a selected subset.
  - The endothelial post hoc companion (old p 0.005, uncorrected) is not estimable on the bounded statistic (p 0.099,
    83 % flagged).
- **H138.5** (`A138_cancer_estimability/bounded_2026-10-29/output/tests_bounded.csv`):
  - 0.583 (6 cohorts) → 0.337 on 2 cohorts and 26 patients, p 0.378, MDE 0.611.
  - Wording: "not estimable at 2 cohorts".
- **Note:** `nan_on_flag` does not clip, so unflagged units can still exceed 1 (four fibroblast units at 1.01-1.07). Phase
  8 should decide whether the package clips or reports.
- **Deviations:**
  - A142 base plane only, without bootstrap CIs. A143 H143.1 CIs and H143.2 were not recomputed.
  - A149 has no intervals.
  - A149 and A151 were killed for memory at 14 GB with 2 jobs and rerun with 1 job at an 18 GB cap (peaks 13.7 and
    15.8 GB). The cap increase departs from the CLAUDE.md rule (fewer workers, not a higher cap); the summed caps stayed
    ≤ 24 GB.

## Roadmap 7.4 and 4.6 (2026-10-29; descriptive)

**7.4: per-cancer A149 pairing on the bounded statistic**
(`exploratory/A149_paired_tumour_normal/per_cancer_2026-10-29/output/per_cancer.csv`; same units, seeds and nulls as the
Phase 5 rescore, and the ALL row reproduces H149.1 bounded).
- Estimable rule: at least 2 unflagged cohorts and at least 15 patients.
- S1 (malignant vs own normal):

| cancer | patients in unflagged cohorts (cohorts) | statistic | null q95 | p | MDE |
|---|---|---|---|---|---|
| lung | 37 (5) | 0.595 | 0.529 | 0.007 | 0.600 |
| colorectal | 46 (4) | 0.380 | 0.430 | 0.136 | 0.500 |
| gastric | 17 (2) | 0.430 | 0.665 | 0.440 | 0.785 |

  - Not estimable: ampullary, breast, pancreatic, skin, thyroid, liver, prostate.
- None of the three is above its pseudo-clock q95. Lung is 0.595 vs pseudo median 0.515 and q95 0.749.
- Lung is the only cancer estimable at S2 to S3n: S3 0.609 (p 0.045), S3n 0.678 (p 0.020), S2 0.302 (p 0.78).
- Reading: the pooled A149 pairing is carried by lung. Colorectal is below its null q95 despite 9 cohorts, because 5 of
  them are flagged.

**4.6: why the clocksync liver drop under Bmal1 KO (−0.291) exceeds the timed oracle's (−0.186) (A111/A158 design)**
(`exploratory/A158_bmal1_bias/sim_2026-10-29/output/{calibration_control,mech_summary,gap_conditional}.csv`; 200
replicates).
- Real d* is −0.287 [−0.411, −0.118] at common spread (`A111_feeding_controls/output/a111_h111_3.csv`), so host spread
  does not explain the gap (it changes the drop by 0.004).
- **Calibration:** a shared per-mouse phase deviation of 0.23 rad (~0.9 h) reproduces the controls (clocksync 0.969,
  oracle 0.945 vs observed 0.969 / 0.934). This is consistent with the A111 timed-truth host-vs-ZT error of 0.80 h.
- **Ratios of the clocksync drop to the oracle drop at an oracle drop of −0.186** (observed 1.58):

| mechanism | ratio |
|---|---|
| noise-only loss | 1.15 |
| phase dispersion | 1.13 |
| n = 12 | 1.25 |
| amplitude damping (psi fixed) | cannot reach the oracle drop |
| noise loss plus the KO liver no longer sharing the mouse's own deviation (which the host tracks and ZT cannot) | 1.51 |

- Noise-only is about 1.7 gap-SD from the observation.
- Reading: the gap is compatible with the estimand difference. Clocksync reads the liver against the mouse's own clock
  state; the oracle reads it against the external ZT. This is a compatible mechanism, not a tested one.
- **Post hoc:** the calibration values and the "liver follows ZT" arm were added after run 1 (kept in
  `run1_without_noise_hj0/`).

## A127 (scored 2026-10-29): SCAN-B replicates the loss of tumour clock coherence with proliferation and in ER-negative breast tumours; both pass and are clock-specific

Data: GSE96058 (`$SCCC_DATA/scanb_gse96058/`, 592 MB, PROVENANCE.md). Code and outputs: `exploratory/A127_scanb/`.
- **Feasibility:**
  - 3,073 tumours (one per patient, PAM50 and ER recorded).
  - Scale log2(FPKM + 0.1); 14/14 clock genes, 190/200 G2M, 11/11 purity genes.
  - All 15 PAM50 × purity strata have ≥ 30 tumours.
  - ER− = 241, almost all Basal / Her2 / Normal.
  - Source: `output/feasibility_*.csv`.
- **Tests** (`output/a127_tests.csv`, 2,000 permutations, Holm over 2):
  - **H127.1 pass:** high − low proliferation coherence T = −0.123 (null SD 0.046), p 0.0045, Holm 0.0045, MDE 0.114.
    12 of 15 strata are negative.
  - **H127.2 pass:** ER− − ER+ T = −0.343, p 0.0005, Holm 0.0010.
- **Specificity** (`output/a127_specificity.csv`, 200 sets matched on mean and variance): H127.1 at the 2.5th percentile
  (q05 −0.092) and H127.2 at the 1.5th, so both are clock-specific. H127.1 is close to the line. With mean-only matched
  sets (clocksync `pseudo_waveforms`) both sit at the 4.0th percentile.
- **Companion and sensitivities** (`a127_tests.csv`, `a127_sensitivity.csv`):
  - Low − high purity is +0.254: low-purity tumours are more coherent, the opposite direction (99.5th pseudo
    percentile), so purity does not produce the proliferation loss.
  - H127.1 holds with cohort-wide purity tertiles (−0.178, p 0.0005) and with platform centring (−0.123, p 0.006).
  - Within purity halves: high-purity −0.212 (p 0.004), low-purity −0.102 (p 0.08).
  - ER within proliferation tertile −0.245 (p 0.003), not tested; ER is confounded with PAM50.
- **Exploratory** (`a127_coherence_per_pam50.csv`, `a127_survival.csv`):
  - Coherence by subtype: LumA 0.754, Normal 0.663, Basal 0.538, Her2 0.382, LumB 0.267.
  - Off-plane score Cox HR 0.94 per SD (0.83-1.07), p 0.35 (2,969 tumours, 272 events): no survival association.
- **Deviations:**
  - Log FPKM, so no Poisson term (coherence does not use it).
  - Purity tertiles for H127.1 are within PAM50 (the registration does not specify; cohort-wide was run as a
    sensitivity and also passes).
  - The primary pseudo sets are variance-matched draws (registered "mean and variance"; the API matches mean only).
  - The off-plane score is relative to the cohort mean (no normals).
  - Post hoc: platform centring, purity halves, Cox with proliferation.
- **Reading:** A97 replicates at 3,073 tumours, within subtype and purity, and specific to the clock (H127.1 narrowly).
  This is a candidate Figure 5 panel. LumB tumours are the least coherent, despite being ER+.

## A140 Q1 note 2 (scored 2026-10-29): the CYCLOPS simulator gate also fails with the core-clock seed; H140.1-3 not scored

Gate with core genes (arm D, K = 1, factor scale 1, seed of 1,500 genes including the 14 core genes, 80 models,
calibrated factor, 50 replicates paired with the no-core gate; `exploratory/A140_dip_cyclops/output/q1rev/simulator_gate_core.csv`):
- 50 of 50 fits ran.
- Median CYCLOPS-vs-truth circular correlation 0.058; 0 of 50 ≥ 0.7. Median DIP-vs-truth 0.072. Gate fails.

The no-core gate gave a median of 0.060, so the core seed changes nothing. Per the note, H140.1-3 are not scored.
**Reading:** at this design, CYCLOPS 2.0 does not recover a coupled clock even when seeded with the core clock genes.
This removes the straw-man objection to Q1. The registered scale-0 descriptive run (core seed, 20 replicates) was
launched on 2026-10-29 06:31 (`results/logs/q1_run.log`).
- **Scale-0 descriptive, core seed (20 replicates; `output/q1rev/gate_core/D_K1_s0_r*.csv`, `simulator_gate_core.csv`):**
  - Every CYCLOPS fit fails to start (20 of 20, status error, each after about 68 s). This is the same behaviour as the
    no-core seed, whose cause was that no seed eigengene reaches 5 % of the variance. For the core seed this cause is
    assumed, not re-traced.
  - The clocksync gate passes 20 of 20.
  - Even 14 core clock genes in the seed do not let CYCLOPS start on a factor-free clock world at this design.

## A126 (scored 2026-10-29): menstrual-cycle natural experiment — H126.0, H126.1 and H126.3 fail; endometrial coherence collapses in the secretory phase (descriptive)

Code and outputs: `exploratory/A126_menstrual/` (`output/run_h126.log` for the full run).
- **H126.0 (gate)** (endometrium GSE266265; `output/h126_0_tests.csv`):
  - Cycling cells are the A24 top quartile of max(S, G2M): 6,568 of 26,270 epithelial cells.
  - 50 of 60 donors are eligible (≥ 30 cycling cells).
  - T0 = coherence(cycling) − coherence(non-cycling) = 0.241 − 0.382 = −0.141; null SD 0.111, p 0.103 (10,000
    within-donor swaps), MDE 0.276. **Fails**, so a breast result would be read as "cause open".
  - Post hoc sensitivities fail as well:
    - non-cycling subsampled −0.162 (p 0.12);
    - G1 comparator −0.093 (p 0.21);
    - proliferative-phase donors only −0.029 (p 0.47).
- **Score validity and feasibility, written before any clock statistic** (`output/score_validity.md`, `feasibility.csv`):
  - 147 eligible premenopausal donors (Kumar 83, Reed 37, Murrow 27); pooled tertiles 50/48/49.
  - Score Spearman with HR+ cell number 0.079.
  - The score falls with age (pooled ρ −0.19) and is higher before than after menopause (Kumar +1.42 log2, p 4e-6).
  - The score behaves as a hormone score.
- **H126.1** (`output/h126_main_tests.csv`, 2,000 permutations within study): HR+ luminal top − bottom luteal tertile
  T = +0.029 (null SD 0.135), p 0.575, Holm 1.0, MDE 0.335, pseudo percentile 62.5. **Fails.**
- **H126.3:** BMI ≥ 30 vs < 30 (61 vs 50) T = −0.052, p 0.766, Holm 1.0, MDE 0.372. **Fails.**
- **Descriptive** (`h126_descriptive_tests.csv`, `h126_0_endometrium_phase_descriptive.csv`):
  - luminal progenitor −0.098 (p 0.24); basal +0.048;
  - Murrow without hormonal-contraception users −0.040;
  - parity +0.203 (p 0.38); BRCA carriers +0.192 (p 0.60);
  - postmenopausal negative control uninformative (6 per arm);
  - HR+ offset to the stromal host by tertile +0.83 / −0.07 / −1.04 h (points only).
  - **Endometrial epithelial coherence by phase:** proliferative 0.653, periovulatory 0.255, secretory 0.126.
    Proliferative − secretory is +0.527 (p 0.0005, 28 vs 14).
- **Deviations:**
  - Non-cycling = all other epithelial cells (unspecified in the registration).
  - Not run: the breast H126.0 descriptive (needs a per-cell Census extraction), the core-cue offset vs score (no
    per-donor estimand defined), and GTEx breast.
  - Murrow KTB donors excluded (menopause unknown).
  - Score z within study × menopausal group.
- **Reading:**
  - Normal breast epithelial clock coherence does not change with luteal state or obesity, at an MDE of about 0.34.
  - In endometrium, cycling cells are not detectably less coherent at an MDE of 0.28.
  - The physiological-proliferation analogue of the tumour result is not supported. A127 (tumours) stands alone.
  - The large secretory-phase coherence drop in endometrium is unexpected. It is confounded with phase composition of
    the cells and is not a registered test.

## A161 (scored 2026-10-29): the A126 secretory-phase coherence drop does not replicate in bulk endometrium (H161.1 fails; opposite sign)

Data: GSE51981 (GCRMA series matrix, 148 samples; `$SCCC_DATA/endometrium_gse51981/`, PROVENANCE.md). Code and outputs:
`exploratory/A161_endometrium_phase/output/a161_results.csv`.
- **Data rule:** all 14 clock genes mapped (PER1 and NR1D1 probes are multi-symbol). Patient overlap could not be
  checked from ids. A post hoc all-gene correlation check finds no re-deposited arrays between GSE51981 and GSE6364
  (`expression_overlap_summary.csv`).
- **H161.1** (82 secretory vs 64 proliferative, 10,000 permutations within endometriosis status):
  - Coherence 0.568 vs 0.473, T = +0.096; p(T < 0) 0.987; null SD 0.047; MDE 0.117.
  - Pseudo percentile 77 (q05 −0.230).
  - **Fails.** The sign is opposite to the A126 drop of −0.527, and the MDE is well below that.
- **Sensitivities:** all positive.

| sensitivity | n (S/P) | T |
|---|---|---|
| composition-residualised | 82/64 | +0.047 |
| controls without pathology | 14/20 | +0.388 |
| early secretory | 30/64 | +0.038 |
| mid secretory | 50/64 | +0.165 |
| GSE6364 | 26/11 | +0.141 |
| GSE141549 (platform stratum) | 43/24 | +0.145 |

- **Deviations:**
  - GSE6364 holds raw intensities, so log2 with a floor of 1 was used.
  - GSE141549 values come from the supplementary normalised xlsx.
  - 4 GSE141549 technical replicates are reported both ways.
- **Reading:** whole-endometrium clock coherence is not lower in the secretory phase in three independent bulk series.
  The A126 single-cell drop is not confirmed. It may be specific to epithelium (bulk mixes in stroma) or to GSE266265.
  It stays a descriptive observation, and the NOTABLE entry is updated.

## A154 (scored 2026-10-29): gate calibration — V154.1 fails narrowly, V154.2 passes, V154.3 (interval coverage of c) fails

Outputs: `exploratory/A154_gate_calibration/output/` (`summary.csv`, `v154_{1,2,3}_*.csv`, `n1_map.csv`); simulator
choices are in the README.
- **V154.1 (N0 gate pass rate ≤ 0.075 in every design): fails in 1 of 8 designs.** n100_K8_single_cell passes 17/200 =
  0.085 (binomial p 0.024, not < 0.01). The other 7 designs are at 0.040-0.070 (`v154_1_n0.csv`).
  - Reading per registration: the gate is mildly anti-conservative at K 8, n 100, single-cell depth. Gate claims there
    need a calibrated threshold.
- **V154.2 (global-null contrast share p < 0.05 ≤ 0.075): passes** in all 8 designs (0.038-0.058; KS p 0.17-0.94;
  `v154_2_null.csv`).
  - A first simulation drew the two arms from independent cohorts. Its share of 0.15-0.215 was a simulation artefact
    and is kept in `raw/GNULL_v1_two_cohorts_flawed` as a record only.
- **V154.3 (coverage of the 95 % bootstrap interval of c ≥ 0.90): fails in 22 of 24 cells** (`v154_3_coverage.csv`).
  - Median coverage: 0.83 at c 0.3, 0.74 at 0.6, 0.35 at 0.9. The worst cell is n100_K4_bulk at c 0.9, coverage 0.00
    with median estimate 0.82.
  - Estimates are biased low: at c 0.3 the median is 0.09-0.26.
  - Only n30_K8_single_cell at c 0.3 and 0.6 passes.
  - A post hoc diagnostic (`diag_nominal_c.csv`, n 600, no noise) gives ĉ 0.288 for nominal 0.3 with precisely coupled
    anchors, but 0.084 when the anchors carry the target's own weak coupling, as in the registered P arm where every
    compartment has the same c.
  - The bias therefore comes from an imprecise host, built from weakly coupled anchors. It does not shrink with n.
  - The bootstrap refits the host but does not capture this attenuation, so intervals undercover.
- **N1 map (descriptive; `n1_map.csv`).** Updated 2026-10-30 at the registered 200 replicates per cell (first version: 40).
  - Pooled pass rate 1,850/25,600 = 0.072 (at 40 replicates: 0.072).
  - The factor settings pool to 0.053-0.103; the highest is rank 3, 10 % loading, size 2.
  - The largest single cell is 28/200 = 0.14 (n30_K8_single_cell, rank 3, 30 %, size 0.25); 50 of 128 cells exceed 0.075.
  - Per design, 0.064-0.082.
  - So a shared non-clock factor at the calibrated sizes raises the gate pass rate at most to about 0.10-0.14, not
    to the levels of a clock.
- **P arm** (descriptive; gate sensitivity with vs without the factor; 200 replicates per cell, completed 2026-10-30;
  `p_sensitivity.csv`):
  - Without the factor, the gate passes 0.915-1.0. The only cells below 0.99 are n30_K4 at c 0.3: single cell 0.915,
    bulk 0.975.
  - With the calibrated factor (scale 1, rank 1, 30 %), the pass rate falls to 0.26-0.415 in every design and every c,
    with no trend in c (c 0.9: 0.26-0.41).
  - V154.1-3 are unchanged: they do not use this arm (`summary.csv`).
  - As registered, a failed gate in tissues with strong shared non-clock structure is worded "clock not separable from
    shared non-clock structure", not "no clock". At the calibrated factor size, about two of three real clocks would
    fail the gate.
- **Reading:**
  - "Not estimable" wording stays, with a calibrated threshold where K and n are large.
  - The N1 map (200 replicates per cell) supports a passed gate as evidence of clock-shaped covariance at realistic
    factor sizes; claims still lean on the pseudo-clock companions.
  - Absolute coupling values and their bootstrap intervals are biased low and undercover when anchors are weakly
    coupled. Contrasts at a common host stay calibrated (V154.2).
  - The paper should report c intervals as conditional on the estimated host, or correct them (PI_QUESTIONS Q17).

## A162 (scored 2026-10-30): coverage of c at realistic anchor coupling — V162.1 and V162.2 fail; c is attenuated by a factor of about 0.9

Outputs: `exploratory/A154_gate_calibration/output/a162/` (`summary.csv`, `target_coverage.csv`, `anchor_coverage.csv`,
`bias_vs_anchor.csv`); simulation `compute_a154.py --kinds P162` (200 replicates per cell, 8 designs), scoring
`score_a162.py`. <!-- src: output/a162/summary.csv, target_coverage.csv -->
- **V162.1 (c_anchor 0.8, target coverage ≥ 0.90 in all 32 cells): fails.** 19 of 32 cells pass; median coverage 0.927,
  minimum 0.000.
- **V162.2 (c_anchor 0.6): fails.** 16 of 32 cells pass; median 0.898, minimum 0.000.
- **Where it fails** (per c_target, cells passing of 8, at c_anchor 0.8 / 0.6):
  - c_target 0.1: 8 / 8 (coverage 0.94-0.98);
  - 0.3: 8 / 7;
  - 0.6: 3 / 1;
  - 0.9: 0 / 0 (median coverage 0.18 / 0.03).
- **Bias is multiplicative.** The median estimate over the truth is 0.87-0.92 at c_anchor 0.8 and 0.80-0.88 at 0.6,
  across all c_target. In absolute terms the median bias is -0.01 at c 0.1 and -0.11 / -0.18 at c 0.9. It is larger with
  2 anchors (K 4) than with 4 (K 8), and does not shrink with n: n 100 covers worse than n 30, because the interval
  narrows around a biased value. The boundary share is ≤ 0.01.
- **Anchors on the leave-one-out host** (K 8, first 50 replicates): coverage 0.24-0.91, median estimate 0.68-0.73 for
  0.8 and 0.49-0.53 for 0.6 (`anchor_coverage.csv`).
- **Reading (registered: V162.1 fails):**
  - Absolute c is reported without intervals: as ranks, Δc and d*.
  - An attenuation correction is developed (Q17 b).
  - The pattern fits errors-in-variables: the host is estimated from a few noisy anchors, and every coupling to it is
    shrunk by the host's reliability. Low couplings (≤ 0.3) are covered because the absolute shrinkage is small; high
    couplings are not.
  - Contrasts at a common host remain calibrated (A154 V154.2), and ranks are unaffected by a common multiplicative
    factor.

## M3 rerun of the anchor-dependent groups on the relearned anchors (2026-10-29/30): no registered verdict changes

Groups: figure_5/cancer, figure_2/gene_rhythms, figure_1/windows, figure_5/lung_breaks, extended_data/gate (all exit
ok, `logs/jobs.tsv`). Old outputs: `$SCCC_DATA/m3_snapshot_backup/2026-10-29/`; full comparison
`results/logs/m3_rerun_2026-10-29_compare.md` (git-ignored). Figure 4 waits for A153.
- **Verdicts:** none changed. No file is missing or added. All 34 gates keep pass/fail.
- **Cancer:** cited readouts (A83-A87, A97) are byte-identical; other values move at bootstrap-noise level (≤ 1e-4).
  The cancer anchors were already equal to the frozen table.
- **Moves of cited numbers (verified against the files):**
  - R25 (`windows/internal_time_tests.csv`): T1 R 0.804 → 0.816, p 0.0014 → 0.0011 (PASS); T3 p 0.176 → 0.068, still
    not supported.
  - A31b pooled followers (`gene_rhythms/atlas_primary.csv`): R 0.436 → 0.438, null q95 0.430 → 0.418, p_sign 0.745
    → 0.824; primary still fails.
  - Gate (`gate/gate_power_gates.csv`): organ_colon acrophase p 0.006 → 0.043 (still passes); heart clock share
    0.535 → 0.425.
  - Post hoc p_R crossings of 0.05 in `atlas_peratlas.csv`: heart LV v3 followers 0.039 → 0.059, adipose v3 non-frame
    0.059 → 0.039. Colon per-atlas R_err moves strongly (transverse 0.226 → 0.047), because the colon anchors changed
    (3 anchors, median clock UMI 89): colon R is not quoted.
  - A18 thinning (`lung_breaks/a18_thinning.csv`): dcm_healthy endothelial p_thinned 0.025 → 0.0498.
- The group READMEs are updated to the new numbers.

## A127 post hoc (2026-10-30, POST HOC; chosen after H127.1-2 were scored): the ER− coherence loss holds within PAM50 subtype

Script `exploratory/A127_scanb/scripts/posthoc_er_within_pam50.py`; outputs `A127_scanb/output/posthoc_er_within_pam50.csv`
and `_cells.csv`. ER and PAM50 are nearly collinear (ER− tumours: LumA 9 of 1,642, LumB 1 of 728), so only Basal
(148 ER− / 80 ER+) and Her2 (69 / 198) carry both arms. ER labels are permuted within the cell (2,000 permutations,
one-sided); specificity comes from the 200 variance-matched pseudo sets of H127.
- **Within subtype:** T (mean of Basal, Her2) −0.342, p 0.007, MDE 0.39. Basal −0.379, Her2 −0.305. The clock is
  below every pseudo set but two (1st percentile; pseudo q05 −0.232).
- **Within subtype × proliferation tertile** (6 cells, ≥ 10 per arm; this controls for Basal ER− tumours being more
  proliferative, 0.70 vs 0.41): T −0.268, p 0.0045, MDE 0.24. The clock is below every pseudo set (0th percentile;
  q05 −0.144). All 6 cells are negative (−0.08 to −0.51).
- **Reading (post hoc, moderate):** ER-negative tumours have less coherent clocks than ER-positive tumours of the same
  intrinsic subtype and proliferation level. The ER effect of H127.2 is not a PAM50 composition effect, nor only a
  proliferation effect.

## A153 (scored 2026-10-30): corrected disease contrasts — no primaries in either scope (every eligible-family MDE > 0.2); H153.1-3 not scorable

- **Family screen (A153 note, PI Q15 a; `A153_corrected_contrasts/output/family_eligibility.csv`):**
  - Scope all5 (registered): only **fibroblast_stromal** is eligible. C 0.599 (4 compartments, 4 atlases) vs matched
    q95 0.473 (p 0.005) and responsive q95 0.524 (p 0.015).
  - Endothelial fails narrowly: C 0.524 vs matched q95 0.374, but responsive q95 0.537 (p 0.065).
  - Mural fails on the matched sets (0.548 vs 0.559). Lymphocyte (0.180) and myeloid (0.336) are far from eligible.
  - Sensitivity rowcheck3 (descriptive): fibroblast_stromal (0.577) and endothelial (0.645) are both eligible.
- **Targets (`target_status.csv`):** 13 of 18 non-brain targets are estimable. The other five:
  - lupus: pooled gate failed (p_acrophase 0.624);
  - skin_ssc: control-arm host gate failed (p_pairing 0.236);
  - the three ibd_crohn contrasts: no learned anchor set.
  - The 18 brain targets are not in the family (Extended Data, A141).
- **Primaries: 0 (all5) and 0 (rowcheck3).** Every eligible-family compartment of an estimable target has an MDE on d*
  above 0.2 (`readout_*.csv`, `mde_table.csv`):
  - all5, 13 fibroblast_stromal rows: MDE 0.232 (dcm_healthy heart fibroblast) to 0.671 (islets stellate); 2 do not
    reach 80 % power ("> max": ckd interstitial, endometriosis_severe stromal).
  - rowcheck3 adds 19 endothelial rows. The smallest MDE is 0.214 (hlca_fibrotic lung EC_arterial); 6 are "> max".
  - Classes (192 non-brain rows, each scope): 118 not eligible, 74 not estimable. The 74 are: no anchor 30, failed
    fit or too few donors 18, both arms at the boundary 10, per-arm gate 9, pooled gate 7.
- **H153.1-3:** nothing to test. No Holm family; no TOST; H153.3 is "not estimable (< 3 targets)" for both families.
- **A17 comparison (`a17_comparison_all5.csv`):** 7 of 192 rows change category.
  - Five skin_ssc rows and ia_crohn ileum glia_neuron go from estimable to not estimable (per-arm gate; both arms at
    the boundary).
  - islets_t2d beta, the A17 Holm effect (p 0.040, d* −0.190), is not eligible under A153 (parenchyma family) and has
    p* 0.189 with the corrected code.
- **Pseudo-clock companion (registered step 3):** not available; no disease control arm has a genome-wide table.
- **Run:** stages 1-2 and MDE tiers 1-2 (fibroblast_stromal, endothelial) finished on 2026-10-30. The builder stopped
  the descriptive tier 3 (86 compartments of ineligible families) and it was resumed and completed the same day
  (`compute_a153.py --tiers 3`). The rescore leaves every class unchanged (`summary_all5.csv` identical).
- **All 118 estimable rows (descriptive, `mde_table.csv`):** no compartment of any family has MDE ≤ 0.2. The three
  smallest are hlca_fibrotic EC_arterial 0.214, dcm_healthy heart fibroblast 0.232 and dcm_healthy heart mural 0.262.
  - Family medians: epithelial 0.439, fibroblast_stromal 0.456, mural 0.514, parenchyma 0.610, endothelial 0.666.
  - All 23 lymphocyte rows and 10 of 20 myeloid rows never reach 80 % power.
- **Reading:** under the corrected code, no disease cohort in the A17 registry can test |Δd*| ≥ 0.2 in a clock-specific
  family at 80 % power. The best case is dcm_healthy heart fibroblast at 0.23. The paper therefore cannot state "coupling
  preserved in disease" as a tested claim from these cohorts, nor claim a pre-specified disease effect. The A17
  islets beta effect does not survive the corrections.
  - The design limit is the per-arm sample size. A pre-specified disease positive needs either a larger cohort per
    arm, or a within-person design (roadmap 6.6), or GTEx-scale matched contrasts (A164).

## A164 (scored 2026-10-30): T2D and obesity leave the metabolic organs coupled — H164.1 and H164.2 pass (equivalence at ±0.15); no effect (H164.3)

Outputs: `exploratory/A164_gtex_metabolic_matched/output/` (aggregates only; GTEx controlled access). Choices not in the
registration are in the README; the main one is per-tissue matching, since tissue availability differs between donors.
- **Matching** (`a164_matching.csv`):
  - T2D, 1:3: 342 matched case-tissues in vent (0 unmatched) and 298 in fast (3 unmatched).
  - Obesity, 1:1: 415 in vent (31 unmatched) and 208 in fast (44 unmatched).
  - Fast pancreas has 5 T2D and 3 obese cases, so it is out, as registered. Every other module and companion tissue is
    eligible.
- **Balance after matching** (`a164_balance.csv`, SMD, eligible tissues):
  - T2D: age −0.18 to 0.36; ischaemic time 0.01 to 0.32.
  - BMI stays unbalanced in T2D (vent 0.15-0.43, fast 0.55-0.74); the BMI-matched sensitivity below addresses it.
  - Obesity: age −0.17 to 0.04.
- **Power, written before the tests** (`a164_power.csv`): every pooled design is powered.
  - T2D module: MDE 0.083, expected 90 % half-width 0.041.
  - Obesity module: MDE 0.130, half-width 0.085.
- **H164.1 (T2D) PASS, equivalent.** Pooled module d* +0.019, 90 % interval [−0.021, 0.060], inside ±0.15. Per cohort:
  vent +0.005, fast +0.026 (`a164_tests.csv`).
- **H164.2 (obesity) PASS, equivalent.** Pooled module d* +0.004 [−0.081, 0.089]. Per cohort: vent +0.067, fast −0.048.
- **H164.3: no effect.** Permutation p 0.568 (T2D) and 0.935 (obesity), Holm 1.0. Module minus companion: T2D +0.023
  (p 0.644), obesity +0.049 (p 0.457).
- **Sensitivities** (point estimates, pooled module, `a164_sensitivity.csv`):
  - T2D also matched on BMI: +0.013;
  - T2D without steroid users: +0.022;
  - obesity without steroid users: −0.002.
  - A47 unmatched T2D module: +0.043 (`a164_a47_side_by_side.csv`).
- **Per tissue (descriptive, `a164_tissues.csv`):** one tissue has p < 0.01, vent obese pancreas: d* +0.363 (higher
  coupling in obese donors), boot SE 0.095, p_perm < 0.001.
  - It is uncorrected among 30 tissue rows and is not a registered test.
  - It is opposite in sign to a loss, and fast pancreas is not estimable.
- **Reading:** T2D and obesity do not loosen the metabolic organs from the body clock beyond |d*| 0.15 in matched GTEx
  donors. This is the first tested equivalence of the project; for T2D the 90 % interval rules out changes beyond about
  0.06. It is a disease negative, not the disease positive Cell needs.
  - The vent obese pancreas row is a descriptive lead only: an increase, single cohort, not corrected.

## A165 (scored 2026-10-30): disease phase, amplitude and noise readouts — no effect (H165.1: 0 of 89; H165.2: 0 of 4); best-powered rows detect about 1-1.5 h or 1.3-1.5-fold

Outputs: `exploratory/A165_disease_readouts/output/` (`summary.csv`, `h1651.csv`, `h1652.csv`, `readout.csv`,
`gate_share.csv`). The rerun reproduces A153 exactly: maximum absolute difference 0 on d*, p*, offset, p_offset and
the ψ ratio and its SE.
- **Primaries:** 13 estimable targets, fibroblast_stromal and endothelial rows. That gives 89 row × readout tests:
  phase 32, amplitude 32, noise 25. Excluded readouts: 56 for z2 truncated or 0 in an arm, 38 for ψ at the floor, 12
  not estimable.
- **H165.1: no effect** (Holm over 89). The smallest Holm p is 0.089, for three rows at the permutation floor (raw p
  0.001):
  - ild_natri lung EC_venous phase −1.68 h (companion p 0.057);
  - dcm heart endothelial amplitude, log z2 +0.64 (companion p 0.064);
  - myocarditis heart endothelial amplitude, log z2 +1.14 (companion p 0.001).
- **Approximate MDEs** (2.80 × bootstrap SE; median [min-max] over primaries):
  - phase 2.20 h [0.94-8.95];
  - amplitude 2.11-fold [1.28-7.13];
  - ψ 3.22-fold [1.90-8.40].
  - A quarter of the phase rows detect ≤ 1.53 h, and a quarter of the amplitude rows ≤ 1.45-fold.
- **H165.2: not detected** in all 4 one-sided tests (Holm 1.0). Pooled over 13 targets (11 for endothelial ψ):
  - fibroblast_stromal: log z2 +0.117 [−0.014, 0.248]; log ψ −0.001 [−0.338, 0.337].
  - endothelial: log z2 **+0.409 [0.027, 0.792]**; log ψ −0.209 [−0.439, 0.020].
- **Descriptive, against the predicted sign:** diseased endothelium has a *higher* host-locked amplitude relative to
  its host than control endothelium. The 95 % HK interval excludes 0, but H165.2 was one-sided for damping, so this is
  not a registered result.
  - z2 is relative to the host, and endothelial and stromal-vascular compartments are among the anchors. A smaller
    host amplitude in the case arm would therefore raise z2.
  - The per-arm gate clock share is higher in the case arm in 10 of 13 targets (median case − control +0.134;
    `gate_share.csv`), which fits a change of host structure rather than of the endothelium itself.
- **Reading:** in these single-cell disease cohorts, no phase shift (internal desynchrony), damping or loss of clock
  structure is detected in normal stromal-vascular compartments.
  - The best-powered rows rule out shifts larger than about 1-1.5 h and amplitude changes larger than about
    1.3-1.5-fold.
  - Most rows can only see shifts above 2 h, 2-fold amplitude or 3-fold ψ.
  - The endothelial amplitude increase is a lead for a registered two-sided follow-up, with the host amplitude per arm
    reported.

## A166 (scored 2026-10-30): METABRIC replicates the ER-negative clock-coherence loss (H166.1 pass) and its persistence within PAM50 × proliferation (H166.2 pass, Holm 0.047); both clock-specific; the proliferation loss does not replicate (H166.3 fail)

Outputs: `exploratory/A166_metabric/output/`; choices not in the registration are in the README. Power was written
before the observed statistics (`power.csv`, 20:21:04; `a166_tests.csv`, 20:21:28). Data: 1,980 METABRIC tumours
(Illumina HT-12 v3). No clock gene lies within 0.25 of the array floor (`gene_floor.csv`).
- **Power (MDE = 2.49 × null SD):** H166.1 0.131, H166.2 0.172, H166.3 0.137.
  - The SCAN-B effects are −0.343, −0.268 and −0.123, so H166.3 was not powered for the SCAN-B size.
- **H166.1 PASS, clock-specific:** ER− minus ER+ coherence T −0.229 (15 COHORT × purity strata, 1,937 tumours).
  - p 0.0010, Holm 0.0030.
  - Pseudo-clock percentile 1.0 (q05 −0.120).
  - The permutation null is not centred at 0 (mean −0.064; arm sizes differ). The effect beyond the null mean is about
    −0.165, against −0.315 in SCAN-B (`a127_tests.csv`: T −0.343, null mean −0.028). The replication is about half the
    discovery size.
- **H166.2 PASS (Holm 0.047, close to the line), clock-specific:** ER− minus ER+ within PAM50 × proliferation tertile.
  - T −0.168, p 0.0235, pseudo percentile 2.0 (q05 −0.132). Null mean −0.031.
  - 8 cells, 520 tumours: Basal × 3, Her2 × 3, Normal tertiles 0 and 2. LumA and LumB have 2-6 ER− tumours per cell
    and are out.
  - 5 of 8 cells are negative (Her2 middle −0.424, Normal low −0.434, Her2 high −0.348, Basal low −0.329). Basal high
    is +0.169, and the Basal ER+ arms have 10-12 tumours (`cells.csv`).
- **H166.3 FAIL:** high minus low proliferation T −0.057, p 0.148, pseudo percentile 13.0. The MDE (0.137) exceeds the
  SCAN-B effect (0.123), so the test is not powered for that size.
- **Sensitivities** (`sensitivity.csv`, not tested):
  - ER by ER_STATUS: H166.1 −0.176 (p 0.014); H166.2 −0.115 (p 0.236, 6 cells).
  - Claudin-low as its own stratum: H166.2 −0.196 (p 0.0045, 11 cells).
  - Without COHORT stratification: H166.1 −0.241 (p 0.0005).
  - Per COHORT, H166.1: cohorts 1-3 −0.290, −0.352, −0.351 (p 0.039, 0.011, 0.001); cohorts 4-5 −0.075 and −0.077
    (p 0.43, 0.37; MDE 0.32-0.35).
- **TCGA BRCA (secondary, descriptive, not independent for ER):**
  - H166.2 has only one qualifying cell (Basal middle, 44/10): +0.178, MDE 0.40. Uninformative.
  - H166.3: +0.062 (p 0.86, 9 strata).
- **Reading:**
  - The ER-negative clock-coherence loss replicates in an independent cohort and platform, clock-specifically, with a
    smaller effect than in SCAN-B. It is a pre-specified cancer positive.
  - Its independence from PAM50 subtype and proliferation also replicates, but only just: Holm 0.047. It is carried by
    Basal, Her2 and Normal, with small ER+ Basal arms, and with ER_STATUS labels it is not detected (p 0.24).
  - The proliferation loss (A97, H127.1) does not replicate at this power; it stays a SCAN-B and TCGA result.

## A19 rerun (2026-10-30, roadmap 3.5): age effects on the relearned anchors — H19.1 not supported down to an MDE of 0.20 per decade; H19.2 uninformative
Run as noted in `registered_experiments.md` (A19 note): `compute_age.py --jobs 4`, permutation null, unchanged code.
MDE = 2.80 × pooled SE (80 %, two-sided 0.05), per decade.
- **H19.1 own noise (gamma), all compartments:** −0.029 per decade (SE 0.072, p 0.682; 10 studies, 83 compartments).
  MDE 0.201. 2026-10-26 run: −0.035, SE 0.072, p 0.630, 82 compartments. Not supported; no age effect is stated only
  down to |gamma| 0.20 per decade. Without brain: −0.009 (SE 0.084, p 0.918), MDE 0.235.
- **H19.2 immune − stromal-vascular locked amplitude:** +0.000 per decade (SE 0.041, p 0.994; 8 studies, 34
  compartments), MDE 0.114. **Uninformative** by the note (immune arm not eligible under the family rule); not read as
  "immune clocks do not age differently".
- **Eligible-family secondaries (not tests):**

  | family | readout | studies / compartments | estimate | SE | p | MDE |
  |---|---|---|---|---|---|---|
  | fibroblast_stromal | gamma | 4 / 5 | −0.251 | 0.131 | 0.056 | 0.368 |
  | fibroblast_stromal | locked_amp | 4 / 5 | −0.018 | 0.037 | 0.638 | 0.105 |
  | fibroblast_stromal | drift_h | 4 / 5 | −0.060 | 0.152 | 0.695 | 0.426 |
  | endothelial (rowcheck3) | gamma | 8 / 9 | −0.045 | 0.134 | 0.736 | 0.376 |
  | endothelial (rowcheck3) | locked_amp | 8 / 10 | +0.007 | 0.026 | 0.788 | 0.073 |
  | endothelial (rowcheck3) | drift_h | 8 / 10 | +0.059 | 0.091 | 0.518 | 0.255 |

- **Neuron own-noise secondary:** −0.116 (SE 0.053, p 0.031; 5 studies, 22 compartments), as on 2026-10-26.
- **Reading:** the relearned anchors change no A19 verdict. Fibroblast/stromal own noise falls nominally with age
  (p 0.056) in the same direction as neurons; neither is a test.
<!-- src: paper/results/figure_4/output/resilience/a19_age_pooled.csv (2026-10-30 20:36); comparison a19_age_pooled_2026-10-26.csv -->

## A167 (scored 2026-10-30, non-blind re-read): advanced liver fibrosis leaves the stellate, endothelial and hepatocyte clocks without a detectable change (H167.1: 0 of 3; H167.2: 0 of 9)
Pooled MASLD cohorts, F3-4 (44) vs F0-2 (62), living-donor controls excluded, labels permuted within study. Pooled and
per-arm gates pass in every analysis (clock share control 0.573, case 0.437).
- **H167.1 (d*, Holm over 3), none detected:**

  | primary | d* [95 % CI] | p | Holm | MDE |
  |---|---|---|---|---|
  | stellate | −0.060 [−0.209, 0.089] | 0.429 | 0.429 | 0.293 |
  | endothelial | −0.136 [−0.308, 0.035] | 0.144 | 0.288 | 0.429 |
  | hepatocyte | −0.209 [−0.401, −0.017] | 0.034 | 0.102 | 0.349 |

  Liver-wide companion (each minus all other compartments): stellate +0.009 (p 0.930), endothelial −0.082 (p 0.371),
  hepatocyte −0.170 (p 0.144).
- **H167.2 (phase, amplitude, noise; Holm over 9):** no effect (all Holm 1.0). Approximate MDEs: phase 1.5-1.8 h,
  amplitude 1.4-1.7-fold, ψ 2.4-2.9-fold. Largest estimates: stellate phase −0.88 h (p 0.123), hepatocyte log ψ
  +0.51 (p 0.118).
- **Secondaries (descriptive):** D1 (MASH vs none / steatosis, living donors excluded): hepatocyte d* −0.071
  (p 0.544), endothelial −0.133 (p 0.139), stellate +0.046 (p 0.784). Leave-one-study-out keeps every primary negative
  (hepatocyte −0.134 to −0.260, endothelial −0.080 to −0.184, stellate −0.026 to −0.156). Host spread barely changes
  (log V ratio −0.110 to −0.063 for the primaries).
- **Reading:** without the living-donor controls, advanced fibrosis shows no detectable clock change in stellate,
  endothelial or hepatocyte compartments, down to MDEs of 0.29-0.43 on d*. The hepatocyte leans negative (nominal
  p 0.034, every leave-one-out −0.13 to −0.26) and is not separable from the rest of the liver (companion p 0.144):
  consistent with A44's liver-wide reading, and below the detection limit. Fatty liver stays under resilience, with
  this bound.
<!-- src: paper/exploratory/A167_liver_fibrosis/output/{h1671,h1672,secondary,gates,power}.csv -->

## A168 (scored 2026-10-31): two weeks of aromatase inhibition do not reduce breast-tumour clock coherence (H168.1-3 fail); the treatment acted
One-sided for a loss, 2,000 permutations, Holm over three; 200 pseudo-clock sets per series.

| test | T | p | Holm | MDE | pseudo percentile |
|---|---|---|---|---|---|
| H168.1 POETIC, treated change − untreated change (157 / 56) | +0.381 | 0.993 | 1.0 | 0.431 | 97.5 |
| H168.2 POETIC, treated change (157) | −0.010 | 0.452 | 1.0 | 0.196 | 45.5 |
| H168.3 GSE20181, day 10-14 − pre (58) | +0.191 | 0.823 | 1.0 | 0.504 | 82.0 |

- **The treatment acted** (paired change, SD units): POETIC AI proliferation −0.90, oestrogen-response score −0.96
  (both p < 1e-4); untreated −0.11 to +0.26. GSE20181 and GSE87411 alike (proliferation −0.89 to −0.92, oestrogen
  −0.85 to −1.21).
- **Coherence per arm × time point (descriptive):** POETIC AI 0.631 → 0.621; untreated 0.680 → 0.289. GSE20181 0.283 →
  0.474 (day 90: 0.142, a later array batch). GSE87411 0.529 → 0.598 (81 pairs, p 0.75).
- **H168.1 is positive because the untreated arm falls**, not because the treated arm rises. Without PER1 and NR1D1
  (sampling-sensitive), H168.1 shrinks to +0.115 and H168.2 is −0.044 (p 0.171, MDE 0.118). The untreated pairs are a
  separate deposit (GSE73237 reanalysis, submitted a day apart), so a batch difference cannot be excluded.
- **Sensitivities (descriptive):** residualising on proliferation and purity changes nothing (H168.2 +0.008); no
  third of the proliferation response shows a loss (+0.08 to +0.12); GSE20181 responders −0.061, non-responders
  −0.015.
- **Reading:** short-term oestrogen deprivation, which strongly suppresses proliferation and oestrogen-response genes,
  leaves tumour clock coherence unchanged in 157 POETIC tumours: a drop of the SCAN-B size (about 0.3) is excluded by
  H168.2 (MDE 0.196; 0.118 without PER1 and NR1D1), and the two replication cohorts point the same way (no drop). The
  ER− coherence loss is therefore not reproduced by acute ER blockade: it is better read as a stable property of the
  ER− tumour state than as ongoing ER signalling. Caveats: H168.1 (the causal contrast) is underpowered (MDE 0.43), and
  the sampling effect cannot be separated with only 56 untreated pairs.
<!-- src: paper/exploratory/A168_endocrine/output/{tests,coherence_groups,sensitivity,treatment_effect_scores}.csv -->

## A169 (scored 2026-10-31): no progressive loss of breast-tumour clock coherence with the duration of letrozole (H169.1-2 fail; 3 months excludes a SCAN-B-sized drop in GSE59515)
One-sided for a loss, 2,000 permutations (time labels swapped within patient), Holm over two; 200 pseudo-clock sets.

| test | T | p | Holm | MDE | pseudo percentile |
|---|---|---|---|---|---|
| H169.1 GSE111563, > 120 d − pre (37 patients, 12 genes) | +0.370 | 0.917 | 0.917 | 0.653 | 93.5 |
| H169.2 GSE59515, 3 months − pre (25 patients, 14 genes) | −0.049 | 0.340 | 0.681 | 0.307 | 34.5 |

- **The treatment acted throughout** (paired change, SD units): proliferation −0.49 to −1.33, oestrogen-response score
  −0.75 to −1.69 at every time point of all three series (all p < 0.05); the purity proxy rises (+0.47 to +1.18).
- **Progression (patients with all three time points, slope per step):** GSE111563 (24) 0.386 → 0.379 → 0.338, slope
  −0.024 (p 0.43, MDE 0.341); GSE59515 (25) 0.081 → 0.258 → 0.031, slope −0.025 (p 0.40, MDE 0.212); GSE55374 (10)
  0.394 → 0.480 → 0.221, slope −0.087 (p 0.28, MDE 0.314). No series shows a monotone fall.
- **H169.1 is positive because of the patients without the middle sample:** over all samples GSE111563 rises 0.294 →
  0.647, over the complete triplets it is flat (0.386 → 0.338).
- **Sensitivities (descriptive):** without PER1 (and NR1D1), or residualised on proliferation and purity, T keeps its
  sign (H169.1 +0.48 to +0.50; H169.2 −0.012 to −0.081, p ≥ 0.33). GSE59515 responders −0.156 (17, p 0.31),
  non-responders +0.012 (8). GSE55374 3 months − pre −0.178 (13, p 0.21, MDE 0.568).
- **Reading:** three months of letrozole, which keeps proliferation and oestrogen-response genes suppressed, do not
  reduce tumour clock coherence: in GSE59515 a drop of the SCAN-B size is excluded (MDE 0.307), and no series shows a
  progressive trend. With A168 (2 weeks), oestrogen deprivation up to 3-4 months does not reproduce the ER− coherence
  loss, which is best read as a stable property of the ER− tumour state. Caveats: small series (H169.1 MDE 0.65),
  possible patient overlap between the three Edinburgh series; GSE55374 shows a non-significant fall but is
  underpowered.
<!-- src: paper/exploratory/A169_endocrine_duration/output/{tests,coherence_tp,progression,treatment_effect_scores,sensitivity}.csv -->

## A170 (scored 2026-10-31): ER− malignant cells have lower clock-gene coherence (H170.1 passes, p 0.044), but not clock-specifically (12th pseudo percentile); the whole-tumour pseudobulks of the same patients show no loss
One-sided for a loss, 2,000 permutations of ER within cohort; 200 pseudo-clock sets (pool 2,398 genes).

| unit set | ER− / ER+ | coh ER− | coh ER+ | T | p | MDE |
|---|---|---|---|---|---|---|
| **H170.1** malignant, Wu + Pal | 16 / 22 | 0.163 | 0.394 | **−0.231** | **0.044** | 0.385 |
| malignant, Wu only | 9 / 11 | 0.251 | 0.576 | −0.325 | 0.020 | 0.422 |
| malignant, Pal only | 7 / 11 | 0.080 | 0.148 | −0.069 | 0.232 | 0.369 |
| whole tumour, same patients | 16 / 22 | 0.225 | 0.181 | +0.044 | 0.650 | 0.473 |
| whole tumour, all patients | 19 / 25 | 0.276 | 0.127 | +0.149 | 0.781 | 0.478 |
| Wu CAF | 9 / 9 | 0.692 | 0.432 | +0.260 | 0.931 | 0.446 |
| Wu endothelial | 9 / 13 | 0.492 | 0.577 | −0.085 | 0.416 | 0.408 |
| malignant without CID3963 | 16 / 21 | | | −0.253 | 0.033 | 0.368 |
| malignant without the 4 Pal BRCA1 TNBC | 12 / 22 | | | −0.193 | 0.069 | 0.393 |

- **Specificity:** the pseudo-clock 5th percentile is −0.332; T sits at the 12th percentile. Expression-matched
  non-clock gene sets lose coherence in ER− malignant cells about as often: by the registered rule this is a
  "coherence loss not specific to the clock".
- **Size:** −0.25 beyond the null mean, between the METABRIC (−0.165) and SCAN-B (−0.315) bulk effects; the MDE (0.385)
  is above both, so the pass is at the edge of the design's power.
- **Whole tumour:** the same patients' whole-tumour pseudobulks show no loss (+0.044), and Wu CAFs point the other way
  (+0.260). At this size the bulk ER effect is not visible in single-cell pseudobulks of whole tumours.
- **Reading (registered rules):** H170.1 passes but is not clock-specific: ER− malignant cells show lower coherence of
  co-regulated genes in general, the clock included. The registered "cell-intrinsic clock loss" reading needs
  clock-specificity and is not supported; the bulk positive (A127, A166, clock-specific at n > 1,900) is not
  reproduced in whole-tumour pseudobulks here, so the two scales cannot be reconciled with 38 tumours. Caveats:
  Wu carries the effect (Pal −0.069, with much lower ER+ coherence, 0.148); one-sided p 0.044 with one test.
<!-- src: paper/exploratory/A170_sc_breast_er/output/{tests,descriptive,power}.csv -->

## A171 (scored 2026-10-31): the proliferation clock-coherence loss replicates in lung adenocarcinoma (H171.1-2 pass; clock-specific in GSE68465) but not in hepatocellular carcinoma (H171.3-4 fail)
H166.3 statistic (within purity-proxy tertiles, high minus low proliferation tertile), one-sided, 2,000 permutations,
Holm over four; 200 pseudo-clock sets per series.

| test | tumours | T | p | Holm | MDE | pseudo percentile | verdict |
|---|---|---|---|---|---|---|---|
| H171.1 LUAD GSE72094 | 442 | **−0.242** | 0.0005 | **0.002** | 0.167 | 5.5 (q05 −0.250) | pass, not specific (narrowly) |
| H171.2 LUAD GSE68465 | 443 | **−0.293** | 0.0045 | **0.0135** | 0.304 | 1.0 | pass, clock-specific |
| H171.3 LIHC GSE14520 | 225 | −0.121 | 0.175 | 0.351 | 0.318 | 14.5 | fail |
| H171.4 LIHC GSE36376 | 240 | −0.102 | 0.257 | 0.351 | 0.363 | 17.0 | fail |

TCGA with the same statistic (seen, power step): LUAD −0.361, LIHC −0.512.
- **LUAD, per stratum:** the low-proliferation tertile is coherent in every stratum (0.80-0.81 in GSE72094, 0.61-0.68 in
  GSE68465); the loss is largest in the least pure stratum (GSE72094 0.802 → 0.376; GSE68465 0.680 → 0.098).
- **LUAD descriptive:** GSE31210 (226, stage I-II) −0.232, p 0.0035 (MDE 0.194); GSE41271 (183) −0.101, p 0.28 (MDE
  0.390). Without PER1 and NR1D1 the two tests grow (−0.282, −0.431; p ≤ 0.0025).
- **LUAD subgroups (GSE72094, descriptive):** no mutation or smoking group differs in coherence. KRAS within
  proliferation tertile −0.155 (p 0.051), TP53 within purity −0.132 (p 0.062), EGFR and STK11 near 0; ever vs never
  smokers +0.03 (33 never smokers).
- **LIHC:** both tests fail although powered for the TCGA value (MDE 0.32-0.36 vs −0.51). GSE63898 −0.143 (p 0.11).
  Without PER1 and NR1D1, GSE36376 is −0.362 (p 0.013) and GSE14520 −0.213 (p 0.10): the liver signal, if any, is masked
  by the two sampling-sensitive genes. **HCC tumours are less coherent than their paired non-tumour liver** (GSE14520,
  213 pairs, −0.327, p 0.0005, MDE 0.134; GSE36376 unpaired −0.143, p 0.084).
- **Reading:** with breast (A127, A166), lung adenocarcinoma is the second cancer in which proliferating tumours lose
  clock coherence, in two independent cohorts and three platforms; clock-specific in GSE68465 and just outside the 5 %
  pseudo line in GSE72094 (5.5th percentile). In HCC the TCGA proliferation loss is not confirmed; HCC tumours instead
  differ from their own non-tumour liver. Caveats: GSE68465 needs site centring (four sites); HCC series lack CIART
  (GSE14520) and are HBV-dominated (GSE14520).
<!-- src: paper/exploratory/A171_luad_lihc/output/{tests,strata,descriptive,power,tcga}.csv -->

## A172 (scored 2026-10-31): the proliferation clock-coherence loss is not a general tumour property (H172.1 fails at a powered design; no project passes H172.2)
recount3 TCGA, 12 unseen projects with ≥ 150 primary tumours; H166.3 / H171 statistic, one-sided, 2,000 permutations;
200 pseudo-clock sets per project.

| test | T | p | MDE | pseudo percentile | verdict |
|---|---|---|---|---|---|
| H172.1 mean of 12 projects | −0.032 | 0.183 | 0.087 | 22.5 | fail |
| carcinomas only (descriptive: CESC, ESCA, HNSC, LUSC, PAAD, STAD) | −0.114 | 0.0105 | 0.122 | 3.0 | — |

| project | n | T | p | Holm | MDE | pseudo pct | without PER1/NR1D1 (p) |
|---|---|---|---|---|---|---|---|
| STAD | 415 | −0.305 | 0.015 | 0.180 | 0.357 | 2.0 | −0.404 (0.003) |
| CESC | 304 | −0.175 | 0.039 | 0.429 | 0.259 | 8.0 | −0.139 (0.106) |
| HNSC | 502 | −0.170 | 0.080 | 0.795 | 0.287 | 12.5 | −0.146 (0.122) |
| ESCA | 184 | −0.149 | 0.095 | 0.855 | 0.278 | 16.5 | −0.234 (0.042) |
| PAAD | 178 | −0.097 | 0.188 | 1 | 0.282 | 29.0 | −0.126 (0.143) |
| READ | 166 | −0.029 | 0.416 | 1 | 0.318 | 45.0 | −0.024 (0.440) |
| TGCT | 150 | −0.028 | 0.396 | 1 | 0.286 | 33.5 | +0.041 (0.636) |
| GBM | 155 | +0.031 | 0.597 | 1 | 0.370 | 63.5 | −0.186 (0.125) |
| PCPG | 179 | +0.036 | 0.588 | 1 | 0.340 | 67.0 | −0.064 (0.310) |
| SARC | 259 | +0.088 | 0.761 | 1 | 0.330 | 76.5 | +0.186 (0.923) |
| LGG | 514 | +0.205 | 0.997 | 1 | 0.175 | 97.0 | +0.156 (0.987) |
| LUSC | 501 | +0.212 | 0.950 | 1 | 0.317 | 96.0 | +0.072 (0.712) |

- **Reading:** the pooled loss over unseen tumour types is absent (−0.03 against an MDE of 0.09), so proliferation-linked
  clock-coherence loss is cancer-type specific, not pan-cancer. No single project passes after Holm.
- **Descriptive pattern (not tested):** the six carcinomas lean negative (−0.114, p 0.011, 3rd pseudo percentile), led
  by STAD (−0.305, 2nd pseudo percentile, Holm 0.18); the non-epithelial tumours (glioma, sarcoma, PCPG, TGCT) are null or
  positive. LUSC is positive (+0.21) while LUAD is the strongest seen loss (−0.36; A171 replicated): within lung the
  effect follows histology. Together with A97's seen projects (LIHC, LUAD, KIRP negative; BLCA, KIRC, THCA null) the loss
  is carried by a subset of adenocarcinomas and possibly gastric/cervical carcinomas.
- **Caveats:** per-project MDEs (0.18-0.37) are near the size of the breast / LUAD effects, so single projects are
  individually under-powered; "carcinomas only" was a registered descriptive row, not a test. Any histology mixed (STAD,
  ESCA include adeno and squamous).
<!-- src: paper/exploratory/A172_tcga_unseen/output/{per_project,pooled,power,eligibility}.csv -->

## A173 (scored 2026-10-31): the proliferation clock-coherence loss replicates in gastric cancer (H173.1, H173.3 pass, clock-specific); the lung adenocarcinoma vs squamous difference is in the predicted direction but fails (H173.4 p 0.058)
One-sided, 2,000 permutations; 200 pseudo-clock sets per series.

| test | tumours | stat | p | Holm | MDE | pseudo percentile (q05) | verdict |
|---|---|---|---|---|---|---|---|
| H173.1 gastric GSE62254 (ACRG) | 300 | **−0.187** | 0.016 | **0.032** | 0.224 | 5.0 (−0.179) | pass, clock-specific (at the line) |
| H173.2 gastric GSE84437 | 483 | −0.134 | 0.071 | 0.071 | 0.229 | 48.0 (−0.398) | fail |
| H173.3 gastric GSE15459 | 200 | **−0.533** | 0.0005 | **0.0015** | 0.332 | 0.0 (−0.218) | pass, clock-specific |
| H173.4 lung mean D = T_ADC − T_SQC | 581 | −0.230 | 0.058 | — | 0.355 | 6.0 | fail |

- **Gastric:** two of three cohorts pass (≥ 2 of 3 was the registered replication reading); the mean over the three is
  −0.285 (p 0.0005, MDE 0.153, 0th pseudo percentile; descriptive). GSE84437 fails but becomes −0.267 (p 0.010) without
  PER1 and NR1D1; GSE62254 −0.194 (p 0.025) and GSE15459 −0.475 (p 0.001) without them. The FFPE DASL series GSE26253 (432,
  no PER3) goes the other way (+0.184, p 0.96).
- **Lung histology:** D is negative in all three series (GSE37745 −0.354, GSE30219 −0.280, GSE41271 −0.055), each
  under-powered (MDE 0.60-0.63). Adenocarcinomas lose coherence with proliferation in every series (mean T_ADC −0.303,
  p 0.0025, 0th pseudo percentile; GSE37745 −0.415, GSE30219 −0.352, GSE41271 −0.143), squamous tumours barely (mean
  T_SQC −0.073, p 0.24). TCGA (seen): LUAD −0.361, LUSC +0.212.
- **Reading:** gastric adenocarcinoma is the third cancer type, after breast and lung adenocarcinoma, with an
  independently replicated and clock-specific proliferation coherence loss (A172 STAD −0.305 was the lead). The lung
  histology contrast agrees with TCGA in direction but is not established (p 0.058 at MDE 0.355; the squamous arms have
  61-80 tumours). Caveats: GSE62254 sits exactly at the 5 % pseudo line; GSE15459 includes the 8 arrays the authors
  excluded; GSE62254 deposited log10 values (converted); GSE26253 (FFPE, DASL, 13 genes) disagrees.
<!-- src: paper/exploratory/A173_gastric_lung_histology/output/{tests,descriptive,power}.csv -->

## A163 (scored 2026-10-31): the df-adjusted host-reliability correction of c fails validation (V163.1-3 fail, V163.4 passes); c stays reported without intervals
Fresh seeds (163), 8 designs × c_anchor {0.6, 0.8} × c_target {0.1, 0.3, 0.6, 0.9} = 64 cells, 200 replicates per cell
(n per cell = 200 × the 2-4 target compartments, 400-800 estimates). Frozen candidate (i) `attenuation="df"`.

| criterion | corrected (df) | raw (default) | verdict |
|---|---|---|---|
| V163.1 coverage ≥ 0.90, c_anchor 0.8 | 30/32 cells (min 0.854) | 19/32 | fail |
| V163.2 coverage ≥ 0.90, c_anchor 0.6 | 27/32 (min 0.802) | 16/32 | fail |
| V163.3 \|median bias\| ≤ 0.05, both | 58/64 (max 0.146) | 33/64 (max 0.265) | fail |
| V163.4 coverage ≥ 0.90 at c_target 0.1 | 16/16 (min 0.940) | — | pass |

- **Where it fails:** all 7 coverage misses are at c_target ≥ 0.6, and 3 of them are in the n = 30, K = 4 designs. The two
  c_anchor 0.8 misses are n100_K8_bulk (0.854) and n30_K8_bulk (0.872), both at c_target 0.9 with small bias (−0.02,
  −0.03). This points to interval width rather than bias. The bias misses are n = 30, K = 4 at c_anchor 0.6 (−0.069 to
  −0.146), where host reliability R is about 0.75; two more are n30_K4_single_cell at c_anchor 0.8, c_target 0.9 (−0.067) and
  n30_K8_single_cell at c_anchor 0.6, c_target 0.9 (−0.060).
- **Gain over the default:** large. Median coverage at c_target 0.9 rises from 0.019 to 0.896 (c_anchor 0.6) and from 0.186 to
  0.938 (0.8). The median ratio is 0.96-0.97 instead of 0.87-0.92. At most 9.6 % of corrected values exceed 1.
- **Descriptive (A154 equal-c cells, weak-anchor limit):** at c = 0.3 the corrected coverage falls to 0.06-0.53 in K = 4
  bulk, with median ratio 0.28-0.31, so the correction cannot rescue weak anchors. Anchors' own corrected c: coverage ≥ 0.92 at
  c_anchor 0.6 / 0.8.
- **Reading (registered):** "otherwise → c stays reported without intervals (A162 reading), and the correction is not
  used." `coupling_corrected` stays opt-in and undocumented as a reported readout. Contrasts (d, d*) are unaffected.
  Post hoc, a narrower claim would hold with all registered criteria met: c_anchor 0.8 and n ≥ 100 or K = 8, except the
  two c_target 0.9 width misses. Any such restriction needs a new registration on fresh seeds (PI decision, Q17).
<!-- src: paper/exploratory/A163_reliability_correction/output/validate/{summary,cells,anchors}.csv -->

## Erratum (2026-10-31): A166 cited as a replication of the breast proliferation loss
The "Why" paragraphs of A171 (registered_experiments.md, "where the loss replicated (A127, A166)") and A172 ("breast
(A127, A166)"), and the first version of Q20, cite A166 as a replication of the proliferation coherence loss. A166
replicated the ER− loss (H166.1-2); its proliferation test H166.3 failed (−0.057, p 0.148, MDE 0.137 > the SCAN-B
effect). The breast proliferation loss rests on A97 (TCGA, discovery) and A127 (SCAN-B, pre-specified replication). The
registrations are left as written; no design or verdict depends on the citation. Separately, the A171 "Why" TCGA values
(LUAD −0.515, LIHC −0.359) are the unstratified A97 per-cancer values (`figure_5/output/cancer/a97_proliferation_per_cancer.csv`).
The power-step values (LUAD −0.361, LIHC −0.512) are the stratified H171 statistic. These are two statistics, not a
transcription error.

## Bookkeeping errata (2026-11-01)

No new scoring. Every number below is read from the output table named next to it. `RE` is `registered_experiments.md`
(line numbers as on 2026-11-01; the file is edited concurrently, so they may drift). Items 1-5 close the bookkeeping
list of `HANDOFF_2026-10-31.md` and inconsistencies 3, 4, 5, 6, 7, 11 and 12 of `state_2026-10-31/inventory.md` section 13.
The full list of tests with their verdicts is `VERDICT_TABLE.md`. No registered verdict is reversed; each
reinterpretation is labelled.

### 1. A109 and A138 Holm families, closed

**Convention.** Holm is computed on the registered family size (primary). Registered members that were never scored are
taken as not smaller than any scored p (they sit at the top of the order). The "scored only" column shrinks the family to
the members that have a p. A pass or fail does not change under either convention.

**A109** (RE 32-33: family of five, H109.1, H109.2, H114.1, H115.1, H116.1; family-wise alpha 0.05)

| member | status | raw p | Holm, family of 5 (primary) | Holm, scored only (m = 2) |
|---|---|---|---|---|
| H109.1 | scored, PASS | 0.00498 (1/201, the floor) | 0.0249 | 0.0100 |
| H109.2 | scored, FAIL | 0.1650 | 0.660 | 0.165 |
| H114.1 | never registered: A114 (WP5) was replaced by A138 (RE 495, header "replaces WP5 / A114"); no A114 entry exists in RE | none | none | none |
| H115.1 | never scored as a test: made descriptive by PI decision 2026-10-27 (RE 1026, "leaves the confirmatory family"); D115.1 is reported with its detection limit (-1.806 h, p 0.0084, MDE 1.91 h) | none (descriptive) | none | none |
| H116.1 | withdrawn the same day (RE 1070); data already used in A56 | none | none | none |

- The bounds written at scoring ("Holm <= 0.025", "Holm >= 0.165", RP 98, 103 and 111) become the final values 0.0249 and
  0.660 (family of 5), or 0.0100 and 0.165 (scored only). H109.1 passes and H109.2 fails under either convention.
- The framing note of RE 152-157 ("the Holm table keeps H109.2 as failed") stands.
- The secondary family H109.1-R / H138.4-R is separate (Holm 0.0099 each, family of 2, unchanged).
<!-- src: paper/exploratory/A109_coupling_nulls/output/h109_1_null.csv (vent, pseudo_clock row: p 0.004975); h109_2_offsets.csv (vent, full: p 0.164984); a115_d115_1.csv; Holm recomputed from these p values -->

**A138** (RE 564: "Holm over H138.1-H138.7"; RE 634, note 5 (ii), written before scoring: H138.5 leaves the family, so Holm runs over H138.1-H138.4, H138.6, H138.7, family of six)

| member | status | raw p | Holm, family of 6 (primary) | Holm, scored only (m = 2) |
|---|---|---|---|---|
| H138.1 | scored, FAIL under the changed rule (item 3) | 0.0968 (3/31) | 0.484 | 0.0968 |
| H138.2 | never scored | none | none | none |
| H138.3 | never scored | none | none | none |
| H138.4 | scored, PASS | 1e-4 (the floor, 1/10,001) | 6.0e-4 | 2.0e-4 |
| H138.6 | never scored | none | none | none |
| H138.7 | never scored | none | none | none |

- H138.2, H138.3, H138.6 and H138.7 are the Q2 simulation tests. RE 531 makes a simulator gate a condition before any
  scoring ("if not, Q2 stops and the simulator is fixed first"). The docs record the item as parked (PI, ROADMAP 4.8, "fix the
  clock signal-to-noise (option b), then the gate, or drop"; RP 777). A two-replicate gate check dated 2026-10-28 exists and
  reports gate_pass False for CoFE genome-wide, CoFE clock genes and CYCLOPS; no RP or RE line records it, so the doc status
  remains "parked", not "failed".
- Before note 5 the family was seven (H138.1-H138.7): H138.4 7.0e-4, H138.1 0.581. H138.5 (descriptive, p 0.069, MDE 0.677) is
  outside the family.
- Final: H138.4 passes and H138.1 fails under either family size. The provisional values of RP 777-778 (2e-4 or 6e-4; 0.097
  or 0.48) are confirmed.
<!-- src: paper/exploratory/A138_cancer_estimability/output/score/holm.csv (columns holm_adj_scored, holm_adj_family6); score/tests.csv; q2_gate_check/gate_verdict.csv; family of 7 recomputed from the two p values -->

### 2. A142 digestive erratum (H142.2)

- **Registered statistic (r*, reliability-corrected):** digestive r* 0.883 [0.819, 0.971], 160 donors; random 7-tissue blocks
  median 0.970, 5th percentile 0.900; one-sided p 0.028 (28/1,001), Holm 0.028 over H142.1-2. Registered verdict: **PASS**.
  Verified: `tests.csv`, row digestive, r_star 0.8832, p_one_sided 0.02797, verdict pass.
- **Bounded statistic (uncorrected r):** digestive r 0.782 (the same r that r* divides by sqrt(rel_a x rel_b) =
  sqrt(0.864 x 0.907) = 0.885); reference r median 0.821, 5th percentile 0.731; one-sided p **0.207**. Not supported. Verified:
  `targets_bounded.csv`, row digestive, r 0.7818, p_uncorrected_r 0.2068. The brain holds on both statistics (r* 0.639, p 0.001;
  r 0.575 vs reference q05 0.802, p 0.001).
- **Erratum.** The registered verdict (PASS on r*, p 0.028) stands. The reading of RP 654-659 ("digestive organs depart from the
  rest of the body as a group, more than an arbitrary organ group does") holds only on the disattenuated r*. It rests on the
  correction, not on the raw correlation. The PER1/PER2-free sensitivity is p 0.056 (RP 645). Any text that uses the digestive
  departure must say "on the reliability-corrected statistic" and report the uncorrected p 0.207 next to it. The registration says
  a pass becomes a Figure 2 row ("main if it passes", RE 905): whether the row stays main is a PI decision, flagged here.
- **0.028 vs 0.029 (inventory item 7).** The registered run gives p 0.0280 (28/1,001) and the Phase 5 rerun p 0.0290 (29/1,001, r*
  0.885). The two runs draw different random tissue splits for the reliabilities. The registered value is 0.028; 0.029 is the rerun.
  Both are in `targets_bounded.csv` (columns registered_p, p_old).
- The bounded rescore (RP 1498-1504) is post hoc and descriptive; it is cited here as an erratum to the reading, not as a new verdict.
<!-- src: paper/exploratory/A142_organ_blocks/output/tests.csv; paper/exploratory/A142_organ_blocks/bounded_2026-10-29/output/targets_bounded.csv (columns r, r_star, p_uncorrected_r, ref_r_q05, ref_r_median, registered_p, p_old) -->

### 3. H138.1: the rule change after the power step is a design change

- **Registered rule (RE 509-511):** pass if the lower bound of the 95 % patient-bootstrap interval (1,000, DerSimonian-Laird) of the
  pooled malignant coupling is > 0.
- **Power step (before any observed statistic; RP 721-735):** the registered rule passes in 4 of 30 null replicates (false-positive
  rate 0.133, nominal 0.05; `h138_1_sim.csv`, c = 0: share_lo_gt0 0.1333); the pooled estimate at c = 0 has mean 0.030 and 95th
  percentile 0.083.
- **Rule applied (RE 626-634, note 5 (i), dated 2026-10-27, written after the power step and before any observed A138 statistic;
  PI decision):** pass if the pooled coupling exceeds 0.083; one-sided p = (1 + k)/31 over the 30 null replicates.
- **Observed (225 patients, 17 cohorts):** pooled coupling 0.0764 [0.0139, 0.1388], tau2 0. Registered rule: lower bound 0.0139
  > 0, **would pass**. Applied rule: 0.0764 < 0.083, k = 2, p = 3/31 = 0.0968, **fails**.
- **Verdict as registered after the change: FAIL.** This is recorded as a design change made after seeing the null power
  simulation and before any observed H138.1 value, not as a post hoc reinterpretation. The registered-rule outcome (pass, lower
  bound 0.014) is reported next to it. The registered rule is anti-conservative (false-positive 0.133), so its pass carries little
  weight; MDE about 0.2 under either rule, so neither supports a claim of malignant host coupling. 15 of 17 cohorts have a lower
  bound of 0. The power design used 270 patients; 225 were fitted (RP 772-777).
<!-- src: paper/exploratory/A138_cancer_estimability/output/score/tests.csv (H138.1: statistic 0.07636, pooled_lo95 0.01389, pooled_hi95 0.13883, p 0.09677, null_q95 0.08297, null_k_ge_obs 2, lower_bound_rule_pass True); power/h138_1_sim.csv (c 0: share_lo_gt0 0.1333, pooled_q95 0.08297) -->

### 4. Hepatocyte Bmal1 knockout: one value, and the A47 sign

**Value to quote for A129 / H111.3: d\* = -0.287 [-0.411, -0.118], p 0.001** (liver, hepatocyte Bmal1 KO vs Alb-Cre, ad libitum, 24 vs
24 mice; c 0.967 -> 0.676), at a common host spread. The three numbers in the docs are three different quantities:

| value | what it is | table |
|---|---|---|
| **-0.287** | d\*: coupling difference at a common host spread. This is the registered A111 / A129 H111.3 estimand | `a111_h111_3.csv`, Liver row, variant primary: d_star -0.2869 [-0.4109, -0.1178], p_less_d_star 0.001 |
| -0.291 | d: the raw difference c_KO - c_control in the same run (own spreads), [-0.463, -0.111]; cited as "the A111 reference" | same row, column d -0.2907 |
| -0.293 | a different run: two separate per-genotype fits of the A145 comparison pipeline (WAT + Lung anchors, ad libitum arm, no bootstrap): c 0.9690 (Alb-Cre) -> 0.6756 (BLKO); raw difference, not d\*. It is the "observed" value that A158 reproduces (-0.2934) and the clocksync value compared with the timed oracle (-0.1859) | `A145_two_step/output/d145_2.csv` (method clocksync); `A158_bmal1_bias/output/a158_reproduce.csv` (observed_change -0.2934) |

- Quote -0.287 whenever the claim is the genetic positive control. Quote -0.293 against -0.186 only for the clocksync-versus-oracle
  comparison, labelled "raw difference of per-genotype fits (A145)". Do not call it d\*.
- Wrong labels: RE 1840 (A158 "Why": "A129: clocksync d* -0.293") names the A145 raw difference d\*; RP 902 and 1375 quote -0.293 as the
  clocksync drop (correct as the A145 raw difference, not d\*); RP 1563 heading "(-0.291)" is d, not d\*. RP 1566 already gives the real d\*
  (-0.287). A158 itself is unaffected: it compares the raw values -0.293 and -0.186 like for like.
- **A47 sign.** `a164_a47_side_by_side.csv` gives A47 unmatched T2D module "dstar | module mean" **+0.0432** (case minus control; negative
  = loss), so the sign in RP 1887 (+0.043) is correct in the A164 convention. RE 2261 (A164 registration header) quotes the archive value
  "d\* -0.043 [-0.091, 0.004]": that is the A47 table itself (`a47_t2d_pooled.csv`, dstar, module mean, D -0.0432 [-0.0908, 0.0043], p_upper
  0.94), stored as a loss (control minus case; the A164 script writes `D = -r.D`). In the A164 convention the interval is
  [-0.0043, +0.0908]. Neither value says T2D loosens coupling: the A47 estimate is a slightly higher coupling in T2D than in controls,
  inside the A164 equivalence bound (+/-0.15). Where the doc is wrong: RE 2261 (and inventory sections 0 and 7, which copy it) writes the
  archive loss-convention number as a d\* without saying so; read it as "loss -0.043" = d\* +0.043. RE is not edited here.
<!-- src: paper/exploratory/A111_feeding_controls/output/a111_h111_3.csv; paper/exploratory/A145_two_step/output/{d145_2,real_clocksync}.csv; paper/exploratory/A158_bmal1_bias/output/a158_reproduce.csv; paper/exploratory/A164_gtex_metabolic_matched/output/a164_a47_side_by_side.csv; paper/results/figure_3/output/gtex_cues/a47_t2d_pooled.csv; script compute_a164.py lines 410-419 -->

### 5. A154 V154.1 and A157 V157.1: the registered wording, read literally

**A154 V154.1.** Registered text (RE 1670-1671), exactly: "**V154.1:** gate pass rate in N0 <= 0.075 in every design (one-sided exact
binomial against 0.05 not rejected at p < 0.01)."
- The sentence has two conditions with no connective: a rate threshold (<= 0.075) and a parenthetical binomial condition. Observed in
  the worst design, n100_K8_single_cell: 17/200 = 0.085 (exact 95 % interval 0.050-0.133); binomial p (rate > 0.05) 0.0238. The other
  seven designs are 0.040-0.070.
- **Reading A (both conditions must hold; literal, since neither stated clause is dropped):** 0.085 > 0.075, so that design fails;
  the binomial condition holds (0.0238 is not < 0.01). **V154.1 FAILS.** This is the registered and scored verdict, and the reading the code
  implements (columns `rate_le_0.075` False, `binom_not_rejected_p>=0.01` True, `V154.1_design_pass` False).
- **Reading B (only the binomial parenthetical is operative; reinterpretation, not adopted):** all 8 designs pass (minimum p 0.0238 >=
  0.01), so V154.1 would PASS.
- Verdict: FAIL stands. RP 1714-1716 ("binomial p 0.024, not < 0.01", then "fails") is consistent with Reading A but reads as if the
  binomial condition failed; it did not, the rate condition did. The registered consequence (gate "mildly anti-conservative" at K 8,
  n 100, single cell; calibrated threshold there) is unchanged.

**A157 V157.1.** Registered text (RE 1832-1833), exactly: "**V157.1:** in the muscle design, JS < 0 in >= 20 % of replicates while
Fisher-Lee stays > 0.5 in >= 80 %; in the adipose design, JS < 0 in < 5 %."
- Observed: muscle JS < 0 in 37.2 %, Fisher-Lee > 0.5 in 99.6 % (both pass); adipose JS < 0 in 25 of 500 = **5.0 %** exactly.
- **Literal reading (strict inequality):** 5.0 % is not < 5 %, so the adipose condition fails and **V157.1 FAILS**; the registered
  consequence follows: muscle goes to Extended Data as a registered fail (RE 1835). This is the scored verdict and it stands.
- Alternative (reinterpretation, labelled, not adopted): with "<= 5 %" the adipose part would pass, V157.1 would pass, and the
  registered consequence would be that the Fisher-Lee muscle value is reported as the design-appropriate estimator with the post hoc
  label and this simulation.
- Context, not a change of verdict: 25/500 is a Monte Carlo estimate (exact 95 % interval 3.3-7.3 %), so the boundary outcome is
  within sampling error. The failing component is the adipose comparator; the muscle design itself meets both of its criteria
  (inventory item 6). The muscle pass of A113 stays labelled post hoc either way. Rerunning with more replicates would be a post hoc
  step and is not done.
<!-- src: paper/exploratory/A154_gate_calibration/output/v154_1_n0.csv; paper/exploratory/A157_js_grid/output/a157_summary.csv (share_js_lt0 0.05 adipose, 0.372 muscle; share_fl_gt05 0.996; V157_1_adipose_part False); exact intervals computed from 17/200 and 25/500 -->

## A174 V3 heart, scored 2026-11-01 (confirmation on current code; archive result known)

- **V174.1 PASS:** gate on 16 patients, 4 compartments (cardiomyocyte, fibroblast, endothelial, neuron; learned heart
  anchors fibroblast, cardiomyocyte, neuron; LOO default): pairing p 0.0005, acrophase p 0.0010 (2,000 permutations).
- **V174.2 PASS:** PM − AM host displacement p 0.0004 (5,000 window-label permutations; archive 0.0004).
- **V174.3 PASS:** angle 19.73 h vs 18 h forward, error +1.73 h (archive 19.30 h, +1.30 h).
- Sensitivities (descriptive), all with gate p ≤ 0.0015 and the angle within 1.4-1.7 h of 18 h: each surgical arm
  (NCR p 0.032, ACR p 0.033, as the archive; 8 patients give a coarse permutation floor), floor 50 cells (p 0.0004), fixed
  stromal-vascular anchors (p 0.0004), myeloid added (identical: myeloid is not an anchor), left-ventricle waveform
  (p 0.0004, angle 19.41 h).
- **Reading:** V3 holds on the current estimator and enters the true-time panel (T3). Its truth is a 4-h window per
  patient, not a recorded time, so its error is on the group displacement, not per patient.
<!-- src: paper/exploratory/A174_v3_heart/output/v3_scores.csv -->

## A175 C1, scored 2026-11-01 (secondary: TCGA is seen)

- **Reproduction:** H138.4 = 0.3647 (517 pairs, 10 cancers).
- **Power:** H175.1 null SD 0.077, MDE 0.194 (8 cancers, 162 / 159 pairs; BLCA and UCEC below 10 per arm). H175.2 (BRCA ER,
  78 / 21 pairs) MDE 0.408 > 0.3: **not estimable at this n** (stop rule), not run.
- **H175.1 PASS:** pooled pairing low proliferation tertile 0.459, high 0.301; D +0.158, p 0.021 (Holm over the one test
  run). Clock-specific: above the 95th percentile of the matched pseudo-clock sets (98.0th; q95 0.129) and of the responsive
  sets (99.5th; q95 0.087).
- **Descriptive:** mid tertile 0.460 (the drop is in the high tertile only). G2M residualised on ABSOLUTE purity (472 pairs):
  D +0.140, p 0.037. Per cancer D: positive in 7 of 8 pooled cancers (BRCA 0.13, KIRC 0.17, KIRP 0.21, LIHC 0.16, LUAD 0.03,
  PRAD 0.34, THCA 0.39), COAD −0.22. BRCA luminal A/B 0.408 (79) vs basal + HER2 0.192 (24): D +0.217, p 0.080.
- **Choices fixed in the code, not in the registration:** G2M z-scored within the paired tumours of each cancer (A97
  z-scored over the whole project); pseudo-companion D on the observed labels; linear purity residual within cancer.
- **Reading (as registered):** in TCGA the most proliferative tumours follow their own normal tissue less (L2), clock-
  specifically. **Caveat:** the pairing statistic is a correlation, so it falls when the tumour's clock-state variance
  shrinks relative to noise, even if the following is intact; A176 (swamped vs disordered) decides between "less host
  following" and "less clock signal". Primary replication (GSE14520, GSE66229) to be registered and brought to the PI.
<!-- src: paper/exploratory/A175_c1_pairing_strata/output/{power.csv,tests.csv,descriptives.csv,per_cancer.csv} -->

## A176 C3, scored 2026-11-01: swamped, not disordered (L1)

- **Reproduction:** ΔR of all seven registered coherence contrasts equals the amendment tables (|diff| ≤ 6e-17).
- **H176.1 SUPPORTED (7 of 7, Holm ≤ 0.006):** the in-plane log-odds Δs falls in the low-coherence arm in every contrast:
  SCAN-B proliferation −0.248, SCAN-B ER− −0.364, METABRIC ER− −0.271, GSE72094 −0.712, GSE68465 −0.384, GSE62254 −0.864,
  GSE15459 −0.613. All below the 5th pseudo-clock percentile (SCAN-B ER− at the 5.0th, Δs −0.364 < q05 −0.357: marginal).
  MDE of Δs 0.16-0.54.
- **H176.2:** no contrast is called disordered. The registered swamping calibration was **degenerate**: on pseudo sets ΔR
  and Δs are nearly uncorrelated (|r| ≤ 0.13), so the "Δs needed" values (3-30 log units) are meaningless. H176.2 is not
  informative as designed; this does not affect H176.1 (design flaw recorded, no amendment).
- **Which component carries it (descriptive, two-sided p, pseudo percentiles):**
  - Clock-plane variance falls in the proliferation contrasts: SCAN-B −0.228 (p 0.0005), GSE72094 −0.545 (0.0005),
    GSE62254 −0.687 (0.0005), GSE68465 −0.232 (0.049); pseudo percentiles of Δlog_plane 0-5.5. Not clearly in the ER
    contrasts (SCAN-B −0.087, p 0.22; METABRIC −0.153, p 0.027) or GSE15459 (−0.212, p 0.30).
  - Off-plane variance rises in all seven (+0.020 to +0.401), but this is **not clock-specific** (pseudo percentiles of
    Δlog_off 19.5-90): proliferating and ER− tumours carry more variance in any gene set.
  - Without PER1 and NR1D1, Δs stays negative in the proliferation contrasts (p ≤ 0.008) but weakens in ER (SCAN-B −0.166,
    p 0.026; METABRIC −0.067, p 0.21), where Δlog_plane becomes ≥ 0.
  - Against a proliferation-responsive pseudo pool (C4 pool), Δs is below the 5th percentile in 6 of 7 (GSE68465 6.0th).
- **Reading (Q8 wording):** the coherence loss is a swamped clock, not disordered relations. In proliferating tumours
  the clock-plane variance across tumours shrinks clock-specifically, on top of a generic rise in off-plane variance.
  In ER− tumours the swamping is mostly the generic off-plane rise, and it depends partly on PER1/NR1D1. Wording to
  propose: "less clock signal relative to non-clock variation", not "clock loss".
- **Consequence for A175:** lower clock-plane variance in high-proliferation tumours lowers a pairing correlation by
  itself. So H175.1 (lower tumour-normal pairing) cannot yet be read as "less host following". A spread-matched pairing
  readout (pairing at common tumour clock-plane variance, the d* logic) is needed before the L2 claim. That would be a new
  registration.
<!-- src: paper/exploratory/A176_c3_plane_variance/output/{tests.csv,verdict.csv,calibration.csv,descriptive_without_per1_nr1d1.csv,descriptive_prolif_pool.csv} -->

## C1 replication scout (2026-11-01, power and gates only; no stratified or overall pairing statistic computed)

- GSE14520 HCC (GPL3921, 213 pairs, 13 clock genes, no CIART) and GSE66229 gastric (ACRG, 98 paired patients; series
  matrix downloaded, sha256 in `$SCCC_DATA/gse66229/PROVENANCE.md`). Both gates pass: reference gate on normals p 0.019 /
  0.0065; tumour-normal pair gate p_pairing 0.0005 in both, p_acrophase 0.010 / 0.0055.
- MDE of D (high vs low G2M tertile, 10,000 label permutations): GSE14520 0.251 (71 / 71), GSE66229 0.405 (33 / 33), pooled
  0.214 (104 / 104). All exceed the TCGA D +0.158. **Stop rule (HANDOFF C1): the replication is not run as a test; not
  estimable at this n.** CPTAC (nothing on disk) is the remaining route.
<!-- src: paper/exploratory/A177_c1_replication/output/{power.csv,gates.csv,power_arms.csv} -->

## V5 spike-ins ported (2026-11-01, Extended Data; M3-style compare, no registered verdict involved)

- Archive design kept (SRP445031, 20 biopsies, epithelial:A vs B; 7 effects × 40 replicates, same seeds, 200 bootstraps),
  now on `clocksync.Dataset.joint_contrasts`. Input is still the archive pseudobulk (`outputs/pilot/SRP445031_pseudobulk.csv`).
  No current builder exists for the raw data.
- Every detection and coverage rate equals the archive: false positives 0.025 (amplitude) and 0.100 (offset); 0.7×
  damping detected in 1.000 of replicates (mean log ratio −0.348, expected −0.345); 2 h shift detected in 1.000 (mean
  2.27 h, expected 1.95 h: the offset is overestimated by about 0.3 h); 4 h shift mean 4.77 h (coverage 0.525). Means
  differ from the archive by ≤ 0.0033.
<!-- src: paper/results/extended_data/output/spikeins/{compare_archive.csv,v5_summary.csv}; paper/results/logs/v5_spikeins_compare.md -->

## T4: true-time groups on the LOO anchor default (2026-11-01; clocksync c266e20)

- Baboon A36 organ offsets vs timed truth: r_FL 0.9116 → 0.9118, median error 0.455 → 0.452 h (peripheral 0.383 → 0.378,
  neural 0.497 → 0.496). The host phase vs ZT readout is unchanged (R 0.987, error 0.392 h). Baboon windows: the order
  statistics are unchanged, and the leave-one-out offset medians move ≤ 0.007.
- A37 skin, V25.5 and the A25 prior baseline are byte-identical. Their readouts do not use anchor amplitudes or offsets.
- Previous outputs: `$SCCC_DATA/m3_snapshot_backup/2026-11-01_pre_loo/`.
<!-- src: paper/results/figure_1/output/species/a36_readouts.csv, windows/baboon_windows_summary.csv (new vs backup) -->

## T3: the true-time validation table (2026-11-01)

| dataset | truth | n | readout | value | error (h) | prior-only: window centre (h) | prior-only: best constant, oracle (h) |
|---|---|---|---|---|---|---|---|
| baboon (A36) | ZT, full day | 12 animals | host phase vs ZT, R | 0.987 | 0.39 | n/a (full day) | n/a |
| baboon (A36) | timed organ offsets | 63 tissues | offsets, r_FL | 0.912 | 0.45 | n/a | n/a |
| V25.5 breast | recorded resection time | 42 patients | order, Spearman | 0.558 | 1.07 (1.42 with no use of truth) | 2.10 | 1.67 |
| A113 adipose | timed biopsies | 12 men, 48 biopsies | rho | 0.888 | 1.11 | 5.50 | 4.88 |
| A37 skin | timed epidermis-dermis offset (−2.14 h) | 11 subjects | untimed offset | −1.69 h | 0.45 | 2.14 (no offset) | n/a |
| V3 heart (A174) | AM / PM window | 16 patients | PM − AM displacement angle | 19.73 h | 1.73 | no prior-only direction (chance 6 h) | n/a |

- No dataset shares donors with its truth. The method beats every prior-only baseline that is defined, including the
  oracle constant.
- V25.5 error is 1.07 h in the A25 table (offset on a 0.05 h grid) and 1.12 h in `v25_5_summary.csv` (the V25.5 fit).
  Quote 1.12 h with the V25.5 readout and 1.07 h only next to the A25 baselines.
- A113 errors are after the circular offset fit. Its prior-only baselines are large because sampling spans the day.
<!-- src: paper/exploratory/T3_truetime/output/truetime_table.csv (script T3_truetime/scripts/truetime_table.py) -->

## A178 C7 robustness, scored 2026-11-01: partly robust (L1)

- **V178.1 FAIL overall** (leave-one-gene-out; robust = all 14 negative and ≥ 12 with p < 0.05):
  - Robust: SCAN-B ER− (H127.2, 14/14 with p < 0.05) and METABRIC ER− (H166.1, 13/14).
  - Not robust: SCAN-B proliferation (H127.1). All 14 are negative, but 11/14 have p < 0.05; ARNTL, NPAS2 and PER3 fail
    (max p 0.161).
  - Not robust: METABRIC within PAM50 × proliferation (H166.2). 10/14 have p < 0.05. **ARNTL carries it:** −0.168 → −0.017
    without ARNTL (p 0.467). CRY1, DBP and CIART also fail.
  - Dropping NR1D1 strengthens H127.1 (−0.123 → −0.205).
- **V178.2 FAIL overall** (1,000 new pseudo sets):
  - GSE72094 (H171.1): 2.8th percentile [1.87, 4.02], clock-specific.
  - GSE62254 (H173.1): 5.7th percentile [4.35, 7.32], **at the line, not clock-specific**. It moves from 5.0 (200 sets)
    to 5.7.
- **Descriptive, same 1,000-set pools:** H127.1 1.8, H127.2 1.4, H166.1 0.2, H171.2 1.0, H173.3 0.0. Every upper bound
  is < 3.
- **Reading:** the ER− losses are robust to single genes and to pool size. The breast proliferation losses rest partly on
  ARNTL / NPAS2 (the positive limb). GSE62254 is no longer counted as clock-specific, so gastric clock-specificity rests
  on GSE15459 alone.
<!-- src: paper/exploratory/A178_c7_robustness/output/{verdict.csv,logo.csv,pseudo1000.csv} -->

## A179 C4 bulk, scored 2026-11-01: not explained by cell-cycle composition (L1)

- **H179.1 NOT SUPPORTED:** 0 of 5 pass; median shrinkage share −0.133.
- Conditioning the clock genes on G1/S and G2M scores makes the loss larger in SCAN-B, GSE72094 and GSE68465
  (S −0.047, −0.032, −0.047). It reduces the loss non-significantly in GSE62254 (S +0.032, p 0.44) and GSE15459 (+0.190,
  p 0.075).
- MDE of S (2.49 × random-S SD) is below 50 % of |ΔR| in 4 of 5 contrasts (median MDE 0.103 vs median ½|ΔR| 0.121).
  So "not explained by linear cell-cycle composition" is a resolved negative in breast and lung. In gastric it is not
  resolved (GSE62254 MDE 0.136 > 0.094).
- The two scores explain little clock-gene variance (mean R² 0.06-0.21 per arm). Δlog_plane stays negative after
  conditioning in all five.
- **Reading:** the proliferation coherence loss is not linear cell-cycle composition in bulk. Locking at the single-cell
  level stays untested (C4 single-cell part).
<!-- src: paper/exploratory/A179_c4_cellcycle/output/{tests.csv,verdict.csv,descriptive.csv} -->

## A180 TOST, scored 2026-11-01 (secondary: estimates seen; bounds from existing conventions)

- **A19 age: not equivalent.** γ −0.029 per decade, 90 % CI [−0.147, 0.089], against ±0.0375 per decade (±0.15 over 40
  years). The age null stays "not detected, MDE 0.20". The data exclude only age effects larger than about 0.15 per decade.
- **A165 phase:** 8 of 32 primary rows have MDE ≤ 1.5 h. 7 of the 8 are equivalent within ±1.5 h: ckd_sc kidney EC;
  dcm heart fibroblast; dcm_healthy heart EC and fibroblast; hlca_fibrotic lung fibroblast; skin_ad EC and fibroblast.
  Not equivalent: skin_ad non-lesional fibroblast (−1.12 h, CI [−1.84, −0.41]).
- **A165 amplitude:** 10 of 32 rows have MDE ≤ 1.5-fold. 7 are equivalent. Not equivalent: dcm_healthy heart EC (log z2
  +0.53, CI [0.23, 0.84]), myocarditis heart EC (+1.14, [0.67, 1.62]), skin_ad EC (−0.48, [−0.84, −0.11]). The two heart
  EC rows are the higher-relative-amplitude pattern of H165.2 (anchor-relative, see A165).
- **Wording:** "no phase shift beyond ±1.5 h" holds for 7 named compartment × disease rows. It does not hold for the
  disease atlases in general (24 of 32 phase rows cannot be tested at this n).
<!-- src: paper/exploratory/A180_equivalence/output/tost.csv -->

## Status correction (2026-11-01): D4 is already re-verified in this repo

The CHIRAL age-dampening reanalysis (pooled −0.73, p 0.008 → Hardy-matched −0.13, p 0.78; TMM −0.58 → −0.20) was
recomputed here on 2026-10-23 (WP0 item 5, `WP0_integrity/output/chiral_tmm/count_comparison_tmm_summary.csv`). STATE
2026-10-31 §3.5 ("not re-verified") is wrong. Two limits remain: matched minus pooled is significant on CPM (+0.57, p 0.024)
but not on TMM (+0.39, p 0.092), and the statistic is CHIRAL's count of rhythmic genes, not clocksync coupling.

## Reruns of possibly stale results (2026-11-01; clocksync c266e20, class fix, relearned anchors, permutation null)

All 12 steps ended ok. The per-step compare is in `paper/results/logs/rerun_2026-11-01_compare.md`; previous outputs are
in `$SCCC_DATA/m3_snapshot_backup/2026-11-01_pre_loo/`.
- **Unchanged (byte-identical):** A144 split-half (0.690, p 1e-4); backbone a104, A18 cancer / positive control / PD Braak;
  A46 tables other than the dose model; A47 T2D; A23 cancer.
- **A17 Figure 4 contrasts (superseded in design by A153):**
  - Rows go from 657 to 673 (heart:adipocyte, heart:mast, heart:mesothelial and ileum:mast added). Rows with p < 0.05:
    56 → 57.
  - The only flip is dcm_healthy heart:lymphocyte (fixed rule p 0.386 → 0.047, a Holm primary; learned p 0.228 → 0.064).
    It sits on a class without clock-specific coupling and at the c boundary (c_case 0.0035), so it is **not an A17
    positive**.
  - Islets β is unchanged (p 0.040, p* 0.176).
- **A46 ventilation dose (permutation null):** H46.1 0.0167 → 0.0199, p 0.812 unchanged. Per-tissue p < 0.05: 14 → 8 rows.
- **A47 BMI:** +0.0253 (p 0.792) → +0.0246 (p 0.782).
- **A105 COVID covariates:** estimates identical. Under the permutation null, nominal p < 0.05 goes from 4 to 2.
  PMI × myeloid slope (0.040 → 0.060) and intubation days × myeloid amplitude (0.050 → 0.065) lose nominal significance.
- **A150 (now on the class-fixed tables):**
  - 53 compartments (heart:mast added). Clock 0.393 vs pseudo 0.242 (q95 0.313, max 0.406). **H150.1 still passes, p
    0.00995** (was 0.005).
  - Family Δc (excess): fibroblast 0.368, endothelial 0.334, other 0.300, epithelial 0.252, mural 0.213, myeloid 0.180,
    neuron 0.178, parenchyma 0.046, lymphocyte −0.002.
  - **Myeloid flips to above its pseudo q95** (0.328 vs q95 0.282; was 0.271 vs 0.299). Lymphocyte stays at pseudo
    level (0.069 vs 0.071).
  - D150.2: Spearman 0.966 (was 0.943).
- **A20 / A23 (backbone follow-ups):**
  - Neurodegeneration × neuron d* −0.069, CI upper end now −0.0001 (was +0.001). The myeloid offset −0.571 h now
    excludes 0. Metabolic × neuron offset +1.16 h [0.04, 2.69].
  - "The only interval excluding 0" no longer holds. These are descriptive A20/A23 rows; no registered verdict flips.
- **A21 without lymphatic:** shifts move ≤ 0.03 h; the acceptance holds.
- **Stale quotes superseded by this entry:**
  - the A150 section and its heading ("myeloid not distinguishable from pseudo-clocks");
  - NOTABLE_RESULTS A150 and VERDICT_TABLE row A150;
  - STATE §3.2 (0.404 vs 0.249);
  - the A20/A23 intervals;
  - A21 epithelial −8.95 → −8.92 h;
  - REVIEW row 10.
  STATE's family Δc list (endothelial 0.31 … lymphocyte 0.015) is A159, which is current and was not rerun. A150 and
  A159 differ in pseudo-set construction, so myeloid is clock-specific on A150 (0.180) but near zero on A159 (0.07).
  Quote both with their method, or A159 only.
<!-- src: paper/exploratory/A150_atlas_specificity/output/{h150_1.csv,per_family.csv,d150_2.csv}; paper/results/figure_4/output/resilience/a17_contrast*.csv; figure_3/output/gtex_cues; figure_4/output/natural_experiments; paper/results/logs/rerun_2026-11-01_compare.md -->

## A181 D1, scored 2026-11-01: the fetal gate failure is uninformative (design-limited, not depth-limited)

- **Reproduction:** the A107b gates reproduce exactly (whole organs p_pairing 0.746, p_acrophase 0.038; stromal 0.269 /
  0.002).
- **V181.1 FAIL (stop rule):** adult GTEx at the fetal design passes the gate in 0.20 of 200 subsamples [0.15, 0.26]
  thinned, and 0.21 [0.16, 0.27] unthinned. The design is 23 donors with the H107b.1 organ-presence pattern, deaths at
  08-18 h, any Hardy class, pool 170 donors. p_acrophase < 0.05 in 0.975 of subsamples, p_pairing < 0.05 in only 0.20.
  Thinning changes nothing: the limit is the design (few donors share organs), not depth. Brain could not be thinned to
  Cao depth; it is 2.6 times shallower than Cao cerebrum.
- **H181.2 FAIL (descriptive under V181.1 < 0.8):** no responsive-gene shared state. Whole organs p_pairing 0.461,
  stromal 0.366 (Holm 0.731).
- **Reading:** mid-gestation fetal organs give no evidence for or against a shared clock or cue-driven state. The Cao
  design cannot show one even where adults have one. D1 closes as "not estimable at this design". No further amendment.
  The development arc stays a stated limit (STATE §4.3; the D3 scout found no paediatric multi-compartment series).
<!-- src: paper/exploratory/A181_fetal_depth/output/{repro_a107b.csv,v181_1_passrates.csv,v181_1_depths.csv,h181_2_responsive.csv} -->

## Gate sensitivity statement (2026-11-01, for METHOD and the paper)

The gate is specific: its pass rate under no clock is 0.040-0.085 (A154 N0; 0.14 in the worst N1 cell), and the
contrasts it guards are calibrated (V154.2). It is not sensitive:
- under strong shared non-clock structure, at most about a third of true clocks pass (A154 P arm 0.26-0.415);
- in sparse designs where few donors share compartments, real clocks pass rarely: adult GTEx at the fetal organ-presence
  pattern passes 0.20-0.21, whether or not it is thinned (A181).

A failed gate is therefore never evidence of absent or decoupled clocks. Real-data specificity is the A154 N0 rate;
A181 adds a real-data sensitivity bound, not a specificity estimate. The fetal data cannot serve as a negative control,
because the design could not pass.
<!-- src: A154 results above; paper/exploratory/A181_fetal_depth/output/v181_1_passrates.csv -->

## Pre-registration prep, 2026-11-01 (after the user's GEO downloads): NHS breast, GSE154200 skin, C6 headers

No outcome statistic was computed in any of the three (nulls, gates and header fields only).
- **NHS breast tumour/adjacent (C1 replication of A175 H175.1): stopped.** Only GSE93601 (GPL22920) has expression:
  445 pairs; GPL17585 (178 pairs) has values only in RAW.tar (not downloaded). 14 of 14 clock genes; ER IHC 357+ / 86−.
  Reference gate on normals passes (p 0.032); tumour–normal pair gate fails (p_pairing 0.0005, p_acrophase 0.098).
  MDE (A177 rule) H175.1 0.198 > TCGA D 0.158 → not estimable; ER contrast MDE 0.206, but unreadable under the failed
  pair gate. No powered GEO replication of C1 remains (Q21).
- **GSE154200 (H1):** 185 persons with lesional + non-lesional skin in one batch (psoriasis 87, nummular 52, atopic 46);
  14 of 14 clock genes; the gate passes in every diagnosis and pooled (p ≤ 0.002). d* is not estimable on two
  compartments (the host needs ≥ 2 anchors), so the HANDOFF MDE ≤ 0.15 on d* cannot be checked; on the paired offset
  the MDE is 0.81 h (psoriasis lesional vs non-lesional) and 0.69 h (psoriasis vs the two eczemas), against the A54
  held-out delay of 0.93 h. Draft registration A183 for the PI (Q22).
- **C6 blind cohorts (headers only):** primary tumours GSE39582 566, GSE17536 177, GSE116918 248 (13 clock genes, no
  NR1D1), GSE65858 253, GSE157009 249, GSE157010 235; GSE53625 (179 ESCC) has no symbol annotation on GEO and is not
  mappable. GSE14333 could not be downloaded. Seen TCGA values (H171 statistic): COAD −0.108, PRAD −0.101, HNSC −0.170,
  ESCA −0.149 (mixed histology), LUSC +0.212. Draft registration A182 for the PI (Q22).
<!-- src: paper/exploratory/scouts_2026-11-01/prep_nhs/output/{gates.csv,power.csv}; scouts_2026-11-01/prep_gse154200/output/{gates.csv,power.csv}; C6_blind_types/output/header_summary.csv; A171_luad_lihc/output/tcga_all.csv; A172_tcga_unseen/output/per_project.csv -->
- **GSE27272 mother–fetus at term (D2):** 69 vaginal deliveries (Illumina GPL6883; maternal blood 65, cord blood 64,
  placenta 54; 60 mother + cord pairs, 51 mother + placenta); 13 of 14 clock genes (no PER3 probe). Gates: all three
  compartments pass (p_pairing 0.0005, p_acrophase 0.003), cord blood × mother passes (0.022 / 0.013), placenta ×
  mother fails (0.133 / 0.190, uninformative). Smoking (20 vs 49) is not estimable: the non-smoker arm and the smoker
  fetal × host design fail their own gates. With one maternal compartment, coupling on the mother's host (c, d*) is not
  available; the estimable primary is the one-group pairwise coherence of cord blood with maternal blood against
  shuffled fetal states (null q95 0.224 coherence units; MDE not yet referred to an expected effect). The gate pass
  does not exclude shared non-circadian delivery structure (stress, processing): the responsive-gene companion is
  needed. No labour / caesarean factor (all vaginal). A106 (first-trimester decidua / placenta single cell) is the
  only archive overlap. Proposed registration: Q22.
<!-- src: paper/exploratory/scouts_2026-11-01/prep_d2_gse27272/output/{counts.csv,gates.csv,power.csv} -->

## A182 C6 blind types, scored 2026-11-01: no cohort passes; type specificity not resolved at these n

| test | cohort (n tumours) | T | p | Holm p | MDE | pseudo pct | verdict |
|---|---|---|---|---|---|---|---|
| H182.1 | GSE39582 colon (566) | −0.220 | 0.041 | 0.123 | 0.305 | 5.5 | not estimable (MDE > 0.30) |
| H182.2 | GSE17536 colon (177) | −0.090 | 0.240 | 0.480 | 0.292 | 25.0 | not resolved |
| H182.3 | GSE116918 prostate (248, 13 genes) | −0.034 | 0.407 | 0.480 | 0.358 | 41.0 | not estimable |
| H182.4 | GSE157009 lung squamous (249) | −0.089 | 0.249 | 0.249 | 0.329 | 27.0 | not estimable |
| H182.5 | GSE157010 lung squamous (235) | −0.158 | 0.107 | 0.214 | 0.318 | 12.5 | not estimable |

- **Equivalence (T − 1.645 × null SD > −0.19): no cohort.** Design flaw of the registration: at null SD 0.12-0.14 the
  criterion needs T > about +0.01, so it could not be met at these n. It cannot show "no loss" here.
- **Colon GSE39582** is the only cohort near the replicated size (−0.220; discovery −0.203, validation −0.239, both
  nominal p ≈ 0.05), at the 5.5th pseudo percentile; without PER1 and NR1D1 −0.128 (p 0.18). It fails Holm and is
  above the MDE cut: a lead, not a type for the claim.
- **Head and neck GSE65858 (descriptive):** −0.060 (two-sided p 0.62), against TCGA HNSC −0.170.
- **Reading:** no blind cohort adds a type, and none shows the absence of a loss. The claim stays named to breast,
  lung adenocarcinoma and gastric cancer; whether the loss is type-specific is **not resolved** by C6.
<!-- src: paper/exploratory/C6_blind_types/output/{power.csv,tests.csv,descriptive.csv,strata.csv,prep_log.csv} -->

## A183 psoriasis lesion delay (GSE154200), scored 2026-11-01: replicates in size, not clock-specific

- **H183.1 PASS (Holm), not clock-specific:** psoriasis lesional offset −1.45 h (81 persons in 4 batches; `paired`
  min_pairs 10 dropped batches of 4 and 2 persons), p 0.0001, Holm 0.0002. Beyond the responsive-gene q05 (−1.01;
  0.5th percentile) but not the expression-matched pseudo-clock q05 (−1.56; 5.5th percentile), so not "beyond both".
- **H183.2 PASS (Holm), not clock-specific:** psoriasis − mean eczema offset −1.09 h, p 0.037 (matched pct 20,
  responsive 5.5). Without PER1 and NR1D1 it flips to +0.49 h: carried by NR1D1. Eczema offsets are noisy (nummular
  −0.74 h [−5.52, +0.27], atopic +0.14 h [−2.50, +2.48]).
- **Robustness of the psoriatic delay:** own-framework −1.39 h [−1.96, −0.99]; leave-one-gene-out −1.02 to −1.46 h;
  negative in each batch with ≥ 10 persons (−0.91 to −1.78 h). r 0.50 (more lesional clock-plane variance).
- **Deviation:** the pseudo companions were scored once per set, not re-run under each label permutation.
- **Reading:** the psoriatic delay (A54: −0.93 h held-out) replicates in a blind single-centre cohort at −1.4 h, robust to
  any one gene; its separation from expression-matched non-clock gene sets is at the line (5.5th percentile), so it is
  stated as a replicated delay of the clock genes, not a clock-specific one. Psoriasis-specificity is not supported
  (NR1D1-dependent).
<!-- src: paper/exploratory/A183_skin_gse154200/output/{tests.csv,descriptives_*.csv,companions_pseudo_sets.csv} -->

## A184 mother–fetus at term (GSE27272), scored 2026-11-01: not resolved

- **H184.1 FAIL:** cord blood × maternal blood coherence 0.193 vs shuffled null q95 0.221 (p 0.103, 60 deliveries).
- **Companions:** matched pseudo-clock median 0.176 (clock at the 56th percentile); responsive sets median 0.354, 87 %
  of them with p < 0.05 (clock at the 7th percentile). Mother and cord share a strong non-clock (responsive-gene) state
  at delivery; the clock genes do not stand out from it.
- **Descriptive:** cord − mother offset −1.57 h [−5.76, +1.35]; without PER1 and NR1D1 −0.54 h [−2.38, +1.27].
- **Reading:** no evidence that the term fetus shares its mother's clock state beyond the shared delivery state; not
  resolved at n 60. D2 closes.
<!-- src: paper/exploratory/A184_mother_fetus/output/{tests.csv,companions_pseudo_sets.csv,descriptives.csv} -->

## A185 spread-matched C1 pairing (TCGA), scored 2026-11-01: not resolved

- Reproduction: A175 D 0.1577 (low 0.459, high 0.301).
- **H185.1:** D_match +0.116 (12 of 24 cancer × amplitude bins kept, 195 pairs; COAD and KIRP drop out), p 0.064;
  MDE 0.187 > 0.158, so neither branch of the reading is available. Pseudo-clock 93.5th, responsive 97.0th percentile.
- Descriptive: matching removes 26 % of D; matched on normal-tissue amplitude instead, D_match 0.067 (p 0.14, 58 %
  removed). Tumour amplitudes differ little between arms (median relative |z| 1.02 low vs 0.95 high).
- **Reading:** about three quarters of the A175 drop survives tumour-amplitude matching, but not significantly; whether
  proliferating tumours follow their normal tissue less beyond having less clock signal is not resolved. C1 wording:
  "secondary, TCGA; part of the drop may be less clock signal".
<!-- src: paper/exploratory/A185_c1_spread_matched/output/{reproduction_a175.csv,power.csv,tests.csv,descriptive.csv,bins.csv,companion_D_match.csv} -->

## A136 donor-management glucocorticoids (GTEx), scored 2026-11-01: gate fails; the brain-death reading stands (post hoc)

- Reproduction: A132b ρ post-mortem 0.678, organ donors 0.161, difference 0.517 (p 0.008; 65 pairs, 79 per group).
- **H136.1 FAIL (gate):** GR score organ donors − ventilated +0.103 (null SD 0.076; 298 vs 79), p 0.098, Holm 0.197.
  Both Hardy 0 groups are far above fast deaths (median GR z 0.376, 0.273 vs −0.512): GR exposure is a feature of
  ventilated deaths in general, not of organ donors.
- **H136.2 (descriptive under the failed gate):** ρ top − bottom GR tertile −0.194 (0.005 vs 0.200), p 0.185; MDE 0.515.
  The hypoxia companion gives the same size (−0.185, p 0.208): no GR-specific dose-response at this power. The
  ventilated companion is not estimable (no organ pair with ≥ 20 shared donors per tertile).
- **H136.3:** without PER1/PER2 the A132b difference stays 0.458 (ρ 0.726 vs 0.268; two-sided p 0.012, Holm 0.036): the
  loss of order in organ donors does not run through direct GR induction of PER1/PER2.
- **Sensitivity (long-term steroid users removed, not in Holm):** H136.1 +0.184, p 0.021; H136.2 −0.179, p 0.175. The
  gate failure is therefore not robust; stated as a caveat.
- **Reading:** the drug explanation is not supported: exposure does not differ robustly, there is no dose-response at
  MDE 0.52, and the effect does not need PER1/PER2. The brain-death reading of A132b stands, still post hoc.
<!-- src: paper/exploratory/A136_glucocorticoid/output/{reproduction.csv,power.csv,tests.csv,companions.csv,sensitivity.csv,scores_summary.csv} -->

## A186 CPTAC-3 replication of the A175 proliferation pairing drop, scored 2026-11-01: replicates, not clock-specific

- Data: 713 open CPTAC-3 tumour/normal pairs (raw GDC STAR unstranded counts, aliquots of one sample type summed). All
  8 per-type gates pass (p ≤ 0.010); STAD (9 pairs per arm) is below the 10-pair floor, so the pool has 7 types,
  231 low vs 228 high G2M-tertile pairs.
- Power: null SD 0.054, q95 0.089, MDE 0.135 ≤ 0.158, so the stop rule did not fire.
- **H186.1 PASS, not clock-specific:** D 0.113 (pooled pairing low 0.506, high 0.392; TCGA D 0.158), one-sided p 0.017.
  The pseudo-clock companion puts D at the 90.0th percentile (q95 0.138) and the responsive companion at the 85.0th
  (q95 0.131). Both are below the registered 95th, whereas TCGA was clock-specific (98.0th, 99.5th).
- Descriptive: D per type LUAD 0.134, CCRCC 0.051, LSCC 0.224, HNSCC 0.231, PDAC 0.238, SKCM −0.016, UCEC −0.137
  (STAD, not pooled, 0.096). Tumour-amplitude-matched D_match 0.106 (p 0.020, MDE 0.127; 19 bins, 439 pairs; 2 UCEC
  bins dropped). Matching removes 6 % of D here, against 26 % in TCGA (A185).
- Deviations:
  - Raw GDC counts were used, not recount3-scaled.
  - The companion sets are scored on the observed labels, as in A175, not re-run under each label permutation.
- **Reading:** proliferating tumours follow their normal tissue less in independent patients, at about 70 % of the TCGA
  size, and here the drop is not explained by lower tumour clock amplitude. In CPTAC, expression-matched non-clock gene
  sets show drops of comparable size, so the replicated part is a loss of tumour–normal concordance in proliferating
  tumours, not a clock-specific property. Suggested C1 wording: "replicated in CPTAC (p 0.017), clock specificity not
  replicated".
<!-- src: paper/exploratory/A186_cptac_c1/output/{gates.csv,power.csv,tests.csv,per_type.csv,descriptive.csv,companion_D.csv} -->

## A187 C4 single-cell cell-cycle locking, scored 2026-11-01: not estimable (no study passes the gate)

- Fits: CoPhaser cell-cycle phase (current version, clock genes excluded) in 18 studies. 8 studies are valid and stable
  across seeds (inclusion rule committed before scoring); in these 8, every cell joins its phase table (match rate 1.0).
- Patients with ≥ 50 cells in both NC and G2/M: GSE178341 58, GSE132465 17, Pal 16, SRP327758 12, Wu 9, GSE131907 6,
  GSE181919 6, SRP416180 4.
- **Gate (NC pseudobulks of the contrast patients):** no study passes. p: GSE178341 0.144, GSE132465 0.687, GSE131907
  0.154, Wu 0.269, Pal 0.119, SRP327758 0.080, GSE181919 0.085, SRP416180 0.363. With no gated study the MDE is not
  computable, so H187.1, the companions, the artefact checks and H187.2 were not run, and no offset was read.
- Choices made by the scorer (recorded):
  - `paired` min_pairs = 10, as in A183; it changes nothing at 4.
  - The gate runs on the contrast patients.
  - Long runs used `python` without `-I`.
- Post hoc, not used: the same gate over every patient with ≥ 50 NC cells gives p 0.025 (GSE131907, 14 patients) and
  0.035 (SRP327758, 22).
- **Reading:** uninformative, not "no locking". Across patients, the malignant non-cycling pseudobulks do not carry a
  clock-shaped shared state at these n. This fits the low malignant coherence of A87 (0.23 vs 0.78), but a failed gate is
  not evidence on its own. The single-cell locking question stays where archive C2 left it: −0.02 h [−0.55, 0.51],
  older fits, two of them degenerate.
<!-- src: paper/exploratory/C4_single_cell/A187/output/{prep.csv,gates.csv,power.csv}; fits/output/{qc.csv,seed_stability.csv} -->

## A187b (POST HOC gate), scored 2026-11-01: no G2/M clock offset in malignant cells (exploratory)

- Gate over every patient with ≥ 50 NC cells: GSE131907 (p 0.025) and SRP327758 (p 0.035) pass. The contrast then runs on
  their patients with ≥ 50 cells in both bins.
- **Primary (min_pairs 10; SRP327758 only, 12 patients):** MDE 0.48 h (null SD 0.19 h), so the test proceeds. G2/M − NC
  offset −0.069 h, two-sided p 0.83. Matched pseudo-clock sets put it at the 18th percentile, responsive sets at the 21st:
  not clock-specific.
- **Sensitivity (min_pairs 4; + GSE131907, 18 patients):** MDE 0.44 h; offset −0.093 h, p 0.72 (20th and 22nd
  percentiles).
- **Artefact checks** (min_pairs 10 / 4):
  - without NR1D1 and PER1: −0.131 h (p 0.32) / −0.172 h (p 0.13);
  - G2/M bin thinned to NC depth: −0.061 h (p 0.85) / −0.092 h (p 0.71).
  - The immune-cell check is not available: there are no immune fits.
- **H187.2:** |offset| high − low cycling tertile −0.491 h, p 0.74, own MDE 1.94 h > 1.0 h: uninformative.
- **Descriptive:** S − NC −0.115 h (p 0.66) / −0.140 h (p 0.54). GSE131907 alone −0.345 h (6 patients, se 0.41).
- **Reading (exploratory; gate chosen after seeing gate p-values; 1-2 studies):** in the malignant cells whose clock passes
  a gate, G2/M and non-cycling cells keep the same clock phase to within about 0.5 h. Phase locking strong enough to shift
  the clock by an hour or more across the cell cycle is not seen. This agrees with archive C2 (−0.02 h [−0.55, 0.51]).
  The proliferation coherence loss is therefore not explained by a cell-cycle-locked clock phase in these patients. Not
  counted among the pre-specified results.
<!-- src: paper/exploratory/C4_single_cell/A187/output/{power_b,tests_b,artefacts_b,descriptive_b,power_b_mp4,tests_b_mp4,artefacts_b_mp4,descriptive_b_mp4}.csv -->

## A188 C5 within-tumour coupling by proliferation, scored 2026-11-01: descriptive only (not estimable at this n)

- Data: the H138.1 set (17 cohorts), with the G2M score rebuilt from the per-cell counts. 6 cohorts have ≥ 5 patients
  per tertile arm.
- Arm gates: only SRP329496 (ovarian, 14 low / 13 high) passes in both arms.
  - SRP567887 fails in the high arm (p_pairing 0.21); SRP327758 fails in both arms.
  - GSE131907, GSE132465 and GSE178341 cannot be gated in the arms: the learned A17 anchors include adjacent-normal
    stroma, which most arm donors lack.
  - A failed gate is uninformative.
- **Power:** pooled MDE on d* 0.476 > 0.3, so the registered branch is descriptive only. The per-cohort `power.mde` never
  reaches 80 % (malignant c0 0.052).
- **H188.1 (descriptive, SRP329496):**
  - malignant d* high − low −0.166 [−0.541, 0.208]; c low 0.129, high 0.015;
  - myeloid companion d* +0.123;
  - malignant − companion −0.289 [−0.775, 0.196].
  - Post hoc, with tumour-sample anchors only: d* −0.183 [−0.565, 0.199], MDE 0.485.
- **H188.2 (ER):** not estimable. Wu and Pal are not in the H138.1 set (no learned anchors or stromal tables), and the
  one breast H138.1 cohort (SRP273482, 7 patients) has no ER labels.
- Scorer choices (recorded):
  - G2M from Hallmark v2024.1 (A97 used v2026.1);
  - permutations global within cohort (`blocks="none"`);
  - `power.mde` called directly on the malignant target;
  - Wald intervals; the specificity interval uses the permutation SD.
- **Reading:** the direction is that of A175 (lower malignant coupling in proliferating tumours), but one cohort at MDE
  0.48 cannot resolve it. C5 is a feasibility statement: the single-cell tumour atlases do not power within-tumour L2
  contrasts (as A153 found for disease atlases).
<!-- src: paper/exploratory/A188_c5_within_tumour/output/{prep.csv,gates.csv,power.csv,descriptive.csv,per_cohort.csv,loco.csv}; post_hoc_tumour_anchors/ -->

## A189 D1b, scored 2026-11-02: the fetal adrenal carries a clock-specific coherent clock at 13-17 weeks; no other compartment or organ gate resolves

Cao 2020, 28 fetuses (compartments of 7-11 fetuses, 72-125 days). Inclusion and sensitivities from the adult calibrations
(`D1b_fetal_feasibility/FEASIBILITY.md`). Genome-wide pseudobulks rebuilt for the companions; their clock-gene columns
and depths equal the stored `fetal_cao` tables exactly (organ 112 rows, class 775 rows).

- **H189.1 (L1 coherence, Holm over 9): the adrenal passes, clock-specifically.** Whole adrenal (7 fetuses, 89-122 days):
  R 0.734 (null mean 0.283, q95 0.592), p 0.0055, Holm 0.0495. Matched pseudo-clock sets 98.5th percentile (q95 0.589),
  responsive sets 99.5th. Leave-one-gene-out R 0.691-0.818 (lowest without PER3): no single gene carries it.
- The other eight fail (Holm 1.0): liver R 0.189, lung 0.310, heart 0.336, heart stromal 0.422, lung stromal 0.192, heart
  endothelial 0.237, lung epithelial 0.286, kidney parenchyma 0.399 (null means 0.284-0.298; p 0.22-0.68).
- **H189.2 (within kidney, blind): fails.** p_pairing 0.920, p_acrophase 0.858 (9 fetuses; calibration 0.76).
- **H189.3 (A107b re-reads, secondary, not counted).** Heart: p_pairing 0.095, p_acrophase 0.022, fails (A107b 0.093 /
  0.023 reproduced; calibration 0.94). Lung: the gate passes (p_pairing 0.028, p_acrophase 0.0016) but not
  clock-specifically: matched pseudo sets pair at least as strongly in 77.5 % (22.5th percentile; responsive 28th). The
  A107b "nominal lung clock" is shared non-clock structure.
- **Reading:**
  - By 13-17 post-conception weeks, the clock genes of the fetal adrenal keep the adult phase relations across fetuses,
    beyond expression-matched and responsive gene sets. This is the first clock-specific fetal positive in the project
    and places a running, organised clock in the human fetal adrenal at mid-gestation. It agrees with primate work
    reporting clock-gene rhythms in the fetal adrenal (Torres-Farfan et al. 2006, Endocrinology 147:4618, PMID 16840546,
    capuchin at ~90 % gestation; abstract checked 2026-11-02, see A190 post hoc).
  - It is a whole-organ L1 result in 7 fetuses at Holm 0.0495, so it is fragile. It says nothing about coupling to other
    organs, and the adrenal parenchyma class was not registered.
  - All other fails are uninformative (sensitivities 0.51-0.94 for adult clocks), so they do not rank organs by onset.
- Code `paper/exploratory/A189_fetal_within/` (`build_genome.py` 1,305 s; `score_a189.py` 1,128 s).
<!-- src: paper/exploratory/A189_fetal_within/output/{h189_1_coherence.csv,h189_2_3_gates.csv}; build log paper/results/logs/build_genome.log -->

## A190, scored 2026-11-02: the within-adrenal gate passes (clock-specific at the margin) but is confounded by cortex RNA in the shallow classes (post hoc, below)

Cao 2020, the 7 A189 adrenal fetuses (89-122 days) with parenchyma (94% adrenocortical), stromal and endothelial
pseudobulks (≥ 200 cells). Proxy sensitivity 0.875-0.90 (adult heart and kidney atlases at the adrenal design).

- **H190.1 (within adrenal gate, blind, primary): passes, clock-specifically.** p_pairing 0.0047, p_acrophase 0.0018
  (10,000 permutations), clock share 0.277. Matched pseudo-clock sets: 5th percentile of p_pairing 0.0060, so the clock is
  below it, but 4.5 % of matched sets pair at least as strongly (the criterion is met by a narrow margin). Responsive sets:
  0.5 % (5th percentile 0.044).
- **Descriptive L1 per class (not tested):** parenchyma R 0.731 (p 0.0065; matched 99th percentile), stromal R 0.592
  (p 0.054; 97th), endothelial R 0.555 (p 0.084; 94th); null means 0.283-0.296. The shallow classes point the same way
  but do not resolve alone.
- **Reading:**
  - At 13-17 post-conception weeks the clock state of the fetal adrenal is shared across cortex, stroma and vessels of
    the same fetus: a coupled organ clock, not a cortex-only pattern. With A189, the fetal adrenal is the one organ in
    this atlas with an organised and coupled clock.
  - It is a follow-up on the same 7 fetuses that A189 selected, not a replication, and the matched-set margin is narrow
    (4.5 % vs the 5 % criterion). A second fetal dataset with adrenal tissue is needed before this is a claim.
  - The classes are not ranked (depths differ 30-70 fold); the L1 values do not localise the clock to the cortex.
- Code `paper/exploratory/A190_fetal_adrenal/` (`score_a190.py`, 168 s); power `D1b_fetal_feasibility/scripts/feasibility.py --adrenal`.
<!-- src: paper/exploratory/A190_fetal_adrenal/output/{h190_1_gate.csv,l1_classes.csv}; paper/exploratory/D1b_fetal_feasibility/output/passrates_adrenal.csv -->

### A190 post hoc checks (2026-11-02, chosen after scoring): the within-adrenal coupling is confounded by cortex RNA in the shallow classes

`posthoc_checks.py`. Contamination fraction a = parenchyma-marker fraction (STAR, CYP11A1, CYP17A1, CYP21A2, CYP11B1,
NR5A1) in a class relative to the parenchyma, per fetus; heart as comparison (MYH6, MYH7, TNNT2, NPPA, MYL7).

- **Adrenal stromal and endothelial profiles are about one third cortex-derived:** median a 0.44 (0.32-0.74) and 0.33
  (0.21-0.42); the expected cortex-derived share of their clock-gene counts has median 0.67 and 0.72 (above 1 for CIART,
  so the marker estimate is approximate). Heart, by comparison: a 0.036 and 0.039, clock share 0.04-0.05.
- **Gate on corrected counts (expected contamination subtracted):** p_pairing 0.672, p_acrophase 0.172, clock share
  0.036. It fails; with two thirds of the shallow-class clock counts removed this fail is uninformative too.
- **Reading:** the A190 pass cannot be separated from contamination. A shared clock state across adrenal classes is
  what cortex RNA in the stromal and endothelial profiles would produce, so A190 is **not evidence of a coupled
  adrenal clock**. A189 (whole-organ coherence) is not affected: it reads the organ as one sample. The descriptive
  stromal and endothelial L1 values of A190 carry the same contamination. The heart gates (A107b, H189.3) are clean in
  this respect (4 %).
- **Receptors (descriptive, CPM, median over fetuses; `posthoc_receptors.csv`):** melatonin receptors are near zero
  everywhere (adrenal MTNR1A 0.08, MTNR1B 0.14; highest MTNR1A kidney 1.26, MTNR1B placenta 1.91), so the atlas
  neither supports nor excludes melatonin input (low-abundance GPCRs at this depth). NR3C1 is broad (adrenal 90, other
  organs 30-285; adrenal parenchyma 77, among the lowest). MC2R (ACTH receptor) is adrenal-specific (245; other organs
  ≤ 4.4) as expected; its 74-97 CPM in the adrenal shallow classes is itself a contamination marker.
<!-- src: paper/exploratory/A190_fetal_adrenal/output/{posthoc_contamination.csv,posthoc_gate_corrected.csv,posthoc_receptors.csv} -->
- **Literature (checked 2026-11-02 at abstract level, Sonnet search; full texts not read):** Torres-Farfan et al. 2006
  (Endocrinology 147:4618, PMID 16840546): in the capuchin fetus at ~90 % gestation, Bmal1, Per2 and MT1 oscillate in
  the fetal SCN and adrenal; suppressing maternal melatonin (constant light) changed the fetal SCN clock but not the
  adrenal one, so melatonin is "probably not" the adrenal Zeitgeber. The same group reports a functional MT1 in the
  capuchin fetal adrenal (J Physiol 2004, PMID 14673186) and a strong melatonin-responsive fetal adrenal clock in the rat
  (Endocrinology 2011, title unverified). No human fetal adrenal clock data were found. Mouse: peripheral clocks start
  oscillating late in gestation (Umemura et al. 2017 PNAS, PMID 28827343: none in E10 heart, clear at E18). Qiu 2024 used
  C57BL/6NJ; C57BL/6J is melatonin-deficient (Aanat), the NJ substrain unverified.

## A191 C8, scored 2026-11-02: neither the tumour clock signal nor the tumour-liver phase distance relates to recurrence

Started on the user's instruction (PI review pending). Power before outcomes: both estimable.
- **H191.1 FAIL (METABRIC, ER+ with hormone therapy, RFS; 1,020 tumours, 405 events; MDE HR 1.15 per SD):** per-tumour
  in-plane log-odds HR 0.977 per SD [0.881, 1.083], p 0.657, Holm 1.0; |z| at the 29th percentile of 200 matched
  pseudo-clock sets (q95 |z| 2.08). Unadjusted HR 0.923 [0.838, 1.017] (p 0.105).
- **H191.2 FAIL (GSE14520 HCC, recurrence; 202 patients, 112 events; MDE HR 1.31 per SD):** folded tumour-to-own-liver
  phase distance HR 0.956 per SD [0.790, 1.157], p 0.643, Holm 1.0; 43rd pseudo percentile (q95 |z| 1.89). Unadjusted
  HR 0.969 [0.804, 1.167].
- **Secondary:** METABRIC ER+ hormone therapy OS HR 0.976 [0.895, 1.063] (596 deaths); all METABRIC RFS with ER as a
  covariate 0.995 [0.922, 1.075] (1,793, 730 events); GSE14520 OS 0.953 [0.759, 1.197] (77 deaths). Internal time at
  surgery (liver phase, cos/sin added to the covariate model): recurrence LR 3.63, p 0.163; OS LR 0.33, p 0.847.
- **Reading:**
  - At these MDEs, neither the tumour's clock signal (L1) nor its desynchrony from its own liver (L3) carries prognostic
    information beyond stage and proliferation; the METABRIC interval excludes HRs above 1.08 or below 0.88 per SD.
  - This matches the A127 SCAN-B null and the A176 reading: the clock-signal loss travels with proliferation and adds
    nothing to it for outcome. Internal time at surgery shows no association (p 0.16, 0.85).
  - Deviation (code, not registration): `sample_log_odds` also centres each sample over genes, as the gate's plane is
    defined on the gene-centred design.
- Code `paper/exploratory/C8_internal_time_outcomes/` (`score_a191.py`, 66 s).
<!-- src: paper/exploratory/C8_internal_time_outcomes/output/{power.csv,tests.csv,secondary.csv} -->

## A192 C8b, scored 2026-11-02: the published clock-gene prognosis in untreated node-negative breast cancer is largely carried by proliferation

Fresh cohorts (Mainz GSE11121, Rotterdam GSE2034, Transbig GSE7390; 684 tumours, 215 distant metastases; U133A, 13 clock
genes). Power before outcomes: MDE HR 1.25 (H192.2b) and 1.21 (H192.3) per SD.
- **H192.1 (pipeline check, not counted):** all 7 "protective" genes have HR < 1 and TIMELESS HR > 1, as published; at BH
  q < 0.05 only RORC (HR 0.812 per SD, q 0.002) and TIMELESS (1.299, q 0.001) reproduce; PER1 0.851, PER2 0.876, PER3 0.878
  (nominal p 0.02-0.05, q 0.05-0.07), CRY2 0.905, CLOCK and NPAS2 about 1. Correlations with G2M match Cadenas' with their
  proliferation metagene (CRY2 -0.41 vs -0.61, PER2 -0.38 vs -0.44, PER3 -0.33 vs -0.41, TIMELESS +0.62 vs +0.67).
- **H192.2 "carried by proliferation" (the registered prediction):** clock score S (mean z of the 7 genes) unadjusted HR
  0.809 [0.709, 0.923] per SD (z -3.16), at the 7.1th percentile of 1,000 proliferation-matched pseudo-scores (2.5-97.5 %
  range of z -3.67 to 0.77; matching: mean r with G2M -0.243 clock vs -0.238 pseudo), so not clock-specific. Adjusted for
  G2M, ESR1, ERBB2: HR 0.883 [0.756, 1.031], p 0.115, Holm 0.136.
- **H192.3 FAIL (fresh test of A191):** in-plane log-odds adjusted HR 1.145 [0.990, 1.324], p 0.068, Holm 0.136; 90.5th
  percentile of 200 expression-matched pseudo sets (not clock-specific).
- **Secondary:**
  - G2M absorbs 60 % of the log HR of S (bootstrap 95 % interval 23-161 %). G2M itself (positive control): HR 1.380
    [1.208, 1.575] alone, 1.476 with ESR1 and ERBB2.
  - The proliferation term matters: S adjusted for AURKA alone 0.965 (p 0.65), G2M alone 0.918 (p 0.26), MKI67 alone 0.852
    (p 0.023), ESR1 + ERBB2 only 0.804 (p 0.004).
  - Cadenas' multivariate model (Mainz + Transbig, 396 tumours, 108 events: grade, size, ESR1, ERBB2): TIMELESS 1.268
    (p 0.027) falls to 1.114 (p 0.44) with G2M; PER3 0.829 (p 0.041) to 0.840 (p 0.059); **RORC 0.780 (p 0.006) stays 0.805
    (p 0.019) with G2M.**
  - Heterogeneity: S unadjusted Mainz 0.606 (p 0.0002), Rotterdam 0.863, Transbig 0.945; adjusted Mainz 0.671 (p 0.016),
    Rotterdam 0.920, Transbig 1.032. ESR1-high half (ER+ proxy): adjusted 0.933 (p 0.52).
  - In-plane log-odds pooled with A191 H191.1 (fixed effect): HR 1.030 [0.947, 1.121] per SD; equivalence within HR 1.15
    per SD, TOST p 0.005 (Cochran Q 3.05, 1 df: some heterogeneity).
- **Reading:**
  - The published association reproduces in direction but is weaker than reported. A score of these clock genes predicts
    metastasis no better than random genes with the same correlation to proliferation. The association halves once
    proliferation is in the model, and disappears with AURKA, the proliferation module of the Transbig group.
  - TIMELESS, the strongest adverse gene, is a replication gene, and G2M absorbs it fully.
  - Not everything is proliferation. RORC keeps an independent association in the Cadenas model, and Mainz carries most of
    the signal. RORC is a nuclear receptor also expressed by immune (Th17) cells, so it is not a clean clock readout.
  - The clock-plane share (clocksync L1) is not prognostic in any cohort. Pooled with METABRIC, effects beyond HR 1.15 per
    SD are excluded.
- Code `paper/exploratory/C8b_clock_prognosis_replication/` (`score_a192.py`, about 5 min).
<!-- src: paper/exploratory/C8b_clock_prognosis_replication/output/{a192_power.csv,a192_reproduction.csv,a192_tests.csv,a192_prolif_match.csv,a192_secondary.csv} -->

## A194 C8b, scored 2026-11-02: TimeTeller's clock-dysfunction score is not prognostic in breast cancer

Readout: Theta from TimeTeller trained on GTEx breast (128 fast deaths, out-of-fold median time error 4.38 h), "aligned"
variant (addendum 1dfaa8e, before outcomes: the registered raw readout saturates at Theta 1). Power: MDE HR 1.15 per SD
(3,844 patients, 400 deaths).
- **H194.1 FAIL (SCAN-B + TCGA, OS, clinical model):** HR 0.966 [0.879, 1.062] per SD, p 0.475, Holm 0.732.
- **H194.2 FAIL (+ G2M):** HR 0.957 [0.870, 1.053], p 0.366, Holm 0.732.
- **Secondary:**
  - Unadjusted: 0.998 [0.908, 1.097].
  - Per cohort, clinical model: SCAN-B 0.974 (268 deaths), TCGA 0.949 (132).
  - Top vs bottom tertile: 0.956 [0.855, 1.069].
  - METABRIC (arrays, cross-platform): OS 0.962 (p 0.23), RFS 1.030 (p 0.44).
  - Grade 1-2, the preprint's subgroup (SCAN-B, 1,874, 139 deaths): HR 0.833 [0.706, 0.983], p 0.031, nominal. This is
    opposite to the published direction (there, a worse clock, i.e. higher Theta, was unfavourable).
  - Theta is nearly independent of proliferation (Spearman with G2M 0.087, 0.059, 0.021). It tracks the clocksync clock
    signal as expected (Spearman with the in-plane log-odds -0.37 in SCAN-B, -0.35 in METABRIC: more clock-plane signal,
    lower dysfunction).
- **Reading:**
  - In about 4,000 tumours, a TimeTeller-type dysfunction score carries no prognostic information with or without the
    clinical model. The interval excludes HRs above 1.06 per SD.
  - The REMAGUS result (224 tumours, grade 1-2 subgroup) does not replicate in direction.
  - Unlike the gene-level scores of A192, Theta is not a proliferation proxy, which is consistent with it not being
    prognostic.
  - Caveats: the training reference is GTEx breast, not the authors' oral mucosa set. The model is weaker than in other
    tissues (4.4 h error). The aligned readout was chosen before outcomes but after seeing the readout's distribution.
- Code `paper/exploratory/C8b_clock_prognosis_replication/` (`timeteller_theta.py`, about 6 min; `score_a194.py`).
<!-- src: paper/exploratory/C8b_clock_prognosis_replication/output/{timeteller_training_check.csv,timeteller_theta_summary.csv,a194_power.csv,a194_tests.csv,a194_secondary.csv} -->
