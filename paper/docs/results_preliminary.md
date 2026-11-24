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
- **Correction (2026-11-06):** "swamped, not disordered" is not supported. H176.2, the only test of disordered relations,
  was not informative, so H176.1 shows less in-plane clock signal and nothing more. It does not separate swamping by
  non-clock variation from a damped or intra-tumour dephased clock or from rewired gene relations, which all lower the
  in-plane share and are clock-specific by construction. Wording: "less clock signal (in-plane share) in proliferating
  and ER− tumours; cause not resolved". Whether bulk data can separate these is registered as A212.
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
    the signal.
  - **Post hoc RORC check** (`posthoc_rorc.py`, after seeing the result):
    - RORC is not an immune proxy here: r -0.08 with an immune score (PTPRC, CD3E, CD2, CD247). It is only weakly tied to
      the other clock genes (r 0.23), G2M (-0.23) and ESR1 (+0.24).
    - Its HR is 0.843 [0.742, 0.958] (p 0.009) with G2M, ESR1 and ERBB2, and stays 0.844 with the immune score added.
    - The other six genes together are null after the same adjustment: 0.934 [0.801, 1.089], p 0.38.
    - So the residual is a single-gene RORC association, not a clock-wide one.
  - The clock-plane share (clocksync L1) is not prognostic in any cohort. Pooled with METABRIC, effects beyond HR 1.15 per
    SD are excluded.
- Code `paper/exploratory/C8b_clock_prognosis_replication/` (`score_a192.py`, about 5 min).
<!-- src: paper/exploratory/C8b_clock_prognosis_replication/output/{a192_power.csv,a192_reproduction.csv,a192_tests.csv,a192_prolif_match.csv,a192_secondary.csv,a192_posthoc_rorc.csv} -->

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

## A193 C8b, scored 2026-11-03: CYCLOPS magnitude does not predict survival in luminal A breast cancer (Li et al. 2024 not replicated)

- **Fits:** CYCLOPS 2.0, 3 fits per cohort, all completed (5 eigengenes each; 20 / 19 / 26 seed genes in TCGA / METABRIC /
  SCAN-B). Seed concordance of the magnitude (Spearman, tumours): TCGA 0.998-0.999, METABRIC 0.976-0.991, SCAN-B
  0.986-0.998. Total run 325 min, peak 12 GB.
- **Power:** H193.1 MDE HR 1.31 (≤ 2.0; 2,342 patients, 493 deaths); H193.2 MDE 1.14 per SD (≤ 1.5; 2,230, 451).
- **H193.1 FAIL (top tertile, age + G2M, METABRIC + SCAN-B):** HR 0.920 [0.761, 1.113], p 0.392, Holm 0.784.
- **H193.2 FAIL (per SD, clinical + G2M):** HR 0.978 [0.892, 1.073], p 0.638, Holm 0.784.
- **Reproduction in TCGA (their cohort, not counted):** their 5-year logistic model (age + M1): OR 0.810, p 0.620 (499
  patients, 31 deaths within 5 years); + MKI67 0.846 (p 0.698); unadjusted 0.810 (p 0.606). Cox OS top tertile, age + T + N:
  HR 0.837 [0.461, 1.519] (52 deaths). The published association is not reproduced with our refit (see deviations).
- **Secondary:**
  - Per cohort, top tertile, age + G2M: METABRIC 1.018 (p 0.87); SCAN-B 0.665 [0.447, 0.990], p 0.044, nominal and
    opposite to the published direction. SCAN-B per SD, clinical + G2M 0.803 (p 0.028), same caveat; unadjusted 1.039.
  - Pooled top tertile with MKI67 instead of G2M 0.883 (p 0.21); age only 0.919 (p 0.39).
  - METABRIC RFS: top tertile 1.107 (p 0.46); per SD 1.086 (p 0.22).
  - Their 5-year logistic endpoint: METABRIC OR 0.958 (p 0.87), SCAN-B 0.671 (p 0.08).
  - Magnitude is not a proliferation proxy in luminal A: Spearman with G2M −0.004 (METABRIC), −0.011 (SCAN-B), −0.190
    (TCGA); with the in-plane log-odds 0.109 (METABRIC), −0.075 (SCAN-B).
- **Deviations from the paper (all registered, stated here for the write-up):** CYCLOPS human seed list instead of the
  authors' Zhang/Meng list (not public); TCGA eigenspace from tumour + adjacent normal jointly, not the authors' normal-only
  space with GTEx/Manchester training samples; METABRIC and SCAN-B fitted on tumours only (no normal tissue); METABRIC
  re-read on all genes; SCAN-B on the cached gene subset. Cohort labels in the scorer were mapped to those of
  `cyclops_magnitude.csv` (a key mismatch fixed before the power step; no analysis change).
- **Reading:** in 2,342 luminal A tumours with 493 deaths, a CYCLOPS clock-strength magnitude carries no prognostic
  information; the tertile interval excludes HR above 1.11. It does not reproduce in TCGA either, so the published result
  may depend on the authors' normal-tissue eigenspace and seed list, which we could not rebuild. With A192 and A194, none of
  the three published clock-survival claims holds beyond proliferation in independent cohorts.
- Code `paper/exploratory/C8b_clock_prognosis_replication/` (`cyclops_magnitude.py --threads 8`, 325 min; `score_a193.py`).
<!-- src: paper/exploratory/C8b_clock_prognosis_replication/output/{cyclops_runs.csv,cyclops_magnitude.csv,a193_power.csv,a193_tests.csv,a193_reproduction.csv,a193_secondary.csv}; paper/results/logs/cyclops_magnitude.log -->

## A198 C8b, scored 2026-11-03: the tumour clock signal does not modify chemotherapy or endocrine benefit beyond proliferation

Readout: per-tumour in-plane log-odds (A191), METABRIC + SCAN-B, Cox OS stratified by cohort with s x T, s, T, G2M x T, G2M
and clinical covariates. Ratio = ratio of HRs per SD of s (treated vs untreated).
- **Power:** H198.1 MDE 1.19 (4,612 patients, 1,294 deaths, 1,431 chemotherapy); H198.2 MDE 1.27 (3,988 ER+, 1,032
  deaths, 3,354 endocrine). Both ≤ 1.5.
- **H198.1 FAIL (chemotherapy):** ratio 1.088 [0.937, 1.264], p 0.268, Holm 0.536; pseudo |z| 65.5th percentile.
- **H198.2 FAIL (endocrine, ER+):** ratio 1.068 [0.913, 1.249], p 0.410, Holm 0.536; pseudo 55.5th percentile.
- **Secondary:**
  - Without the G2M x T term: 1.074 (p 0.35) and 1.068 (p 0.41).
  - MKI67 in place of G2M: 1.054 and 1.094.
  - Chemotherapy, per cohort: METABRIC 1.202 [1.008, 1.433], p 0.040; SCAN-B 0.758 [0.569, 1.010], p 0.058. Nominal and
    in opposite directions, so heterogeneity, not a shared effect. METABRIC RFS 1.118 (p 0.23). Within ER+ 0.951 (p 0.62).
  - Endocrine, per cohort: METABRIC 1.078, SCAN-B 0.905; METABRIC RFS 1.041.
- **Reading:** no evidence of a clock x treatment interaction larger than a ratio of about 1.26 (chemotherapy) or 1.25
  (endocrine) per SD, in either direction (upper bounds of the intervals), and nothing above matched pseudo-clock sets.
  Treatment is not randomised, so this is observational. Wording: "no evidence of an interaction larger than X", not "the
  clock does not affect treatment response". Dosing time is not recorded, so chronotherapy is not addressed.
- Code `paper/exploratory/C8b_clock_prognosis_replication/scripts/score_a198.py` (about 5 min).
<!-- src: paper/exploratory/C8b_clock_prognosis_replication/output/{a198_power.csv,a198_tests.csv,a198_secondary.csv} -->

## A199 C8b, scored 2026-11-03: CYCLOPS magnitude barely tracks clock strength in GTEx, with or without whole-day coverage

GTEx fast deaths, 12 tissues (breast + 11 with most window donors), W = 08-18 h deaths (58-98 donors, median 84) vs Fs = as
many random fast deaths from the whole day; 72 CYCLOPS fits, all ok (108 min, peak 11 GB). Power: MDE 0.120 (≤ 0.20).
- **H199.1 FAIL:** median rho_Fs − rho_W 0.075, Wilcoxon p 0.259, Holm 0.259. Coverage does not detectably change what
  magnitude tracks.
- **H199.2 PASS (formally), negligible size:** median rho_Fs 0.055 (range −0.095 to 0.145), p 0.021, Holm 0.043. With
  whole-day sampling, magnitude carries a trace of clock strength (about 0.3 % of rank variance). Median rho_W 0.020 (−0.193
  to 0.241).
- **Secondary:**
  - R² of magnitude on the hour of death: W 0.016, Fs 0.017 (medians). Magnitude does not follow the hour.
  - |Spearman| with the first non-clock principal component: 0.062 / 0.054; magnitude is not a dominant non-clock axis either.
  - Seed concordance (minimum Spearman of three sample-order fits): W median 0.447 (0.138-0.918), Fs 0.425 (0.295-0.787).
    At 58-98 samples per tissue, the magnitude ranking is not reproducible across fits (tumour cohorts in A193, 499-1,642
    samples: ≥ 0.976).
  - CYCLOPS phase vs recorded time, |circular r|: W 0.069, Fs 0.101 (descriptive; range restriction in W).
- **Post hoc (truth reliability, `posthoc_truth_reliability.csv`):** timed clock strength computed from two disjoint halves
  of the clock genes agrees at Spearman 0.612 (median; 0.416-0.720). The truth is reliable, so the low rho is not truth noise.
- **Reading:**
  - The user's hypothesis (magnitude fails specifically in working hours) is not supported: magnitude is a poor
    clock-strength readout in single GTEx tissues whatever the coverage.
  - Main caveat: 58-98 samples per fit is small for CYCLOPS, and the fits are unstable across seeds. The result says
    nothing direct about magnitude in 500-1,600-sample tumour cohorts, where the fits are stable.
  - Supportable wording: "In GTEx tissues with recorded time of death, CYCLOPS magnitude was nearly unrelated to a timed
    measure of clock strength (median Spearman 0.06 with whole-day sampling, 0.02 in working hours)." Not "CYCLOPS cannot
    be used in hospital hours".
- Code `paper/exploratory/A199_cyclops_working_hours/scripts/a199.py` (`fit --jobs 3 --threads 4`, 108 min; `score`).
<!-- src: paper/exploratory/A199_cyclops_working_hours/output/{power.csv,fits.csv,tests.csv,per_tissue.csv,posthoc_truth_reliability.csv} -->

## A200 registered rerun of V25.5, scored 2026-11-03: breast-surgery internal time is clock-specific (seen data)

GSE233242, 42 patients, tumour + adjacent normal, recorded resection hour 08-18 h; V25.5 estimator unchanged.
- **H200.1 PASS:** Spearman of the internal-time coordinate with the recorded hour 0.558; 200 expression-matched
  pseudo-clock sets through the same pipeline: median −0.053, q95 0.215, none reached it (p 0.005, the floor; Holm
  0.006). 0 of 200 fits failed.
- **H200.2 PASS:** median error with no recorded time used (estimates centred on 13:00) 1.42 h vs 2.28 h for shuffled
  hours (p 0.003, Holm 0.006); pseudo-clock sets 2.51 h (median).
- **Secondary:** without PER1 0.560 (error 1.38 h); without PER1 and NR1D1 0.233 (1.97 h); normal alone 0.528, tumour alone
  0.517; breast waveform 0.548.
- **Reading:** the per-patient ordering by surgery hour is carried by the clock genes, not by matched genes that could track
  the hour (fasting, anaesthesia, handling), and it does not rest on glucocorticoid-induced PER1. It leans on NR1D1 (0.233
  without PER1 and NR1D1). Seen data: a robustness test, not a replication; a fresh cohort with recorded resection times
  remains the decisive step.
- Code `paper/exploratory/A200_breast_surgery_rerun/scripts/score_a200.py` (about 4 min).
<!-- src: paper/exploratory/A200_breast_surgery_rerun/output/{tests.csv,secondary.csv} -->

## A201 aorta jet lag (Auerbach et al. 2026), scored 2026-11-03: smooth muscle and fibroblasts shifted together after a 6 h advance (H201.1 not estimable by 0.03 h)

Pools of 2 aortas; aligned 14 libraries, misaligned 8 (5 days after a 6 h advance), Bmal1 KO 10; 14 clock genes.
- **Power:** SD of D 0.474 h (library jackknife; cell bootstrap 0.118), 1.645 SD = 0.78 h > 0.75 h: **H201.1 not
  estimable** (registered stop rule). Pool-to-pool variation dominates.
- **H201.0 PASS:** the clock phase moved: mean shift −2.86 h, 90 % [−3.95, −1.77] (about half of the 6 h advance after
  5 days).
- **H201.1 not estimable; descriptive:** SMC −2.85 h, fibroblast −2.87 h, D +0.02 h, 90 % [−0.76, +0.80]. The interval
  lies inside the registered ±1.5 h bounds, but the power criterion failed narrowly, so it is not counted as equivalence.
- **Secondary:**
  - Endothelial cells −2.29 h (≥ 30 cells in every library); macrophages −3.98 h (fewer than 30 cells in some libraries).
  - Per sex: female SMC −3.28 / fibroblast −3.24 h; male −2.41 / −2.48 h. Both sexes show D ≤ 0.07 h.
  - Amplitude misaligned / aligned: SMC 0.875, fibroblast 0.888 (the reduced amplitude the paper reports).
  - Bmal1 KO / aligned amplitude: SMC 0.091, fibroblast 0.240 (clock collapse as expected).
  - clocksync gate (libraries as donors, SMC + fibroblast): aligned p_pairing 0.001, p_acrophase 0.001, clock share 0.77;
    misaligned 0.001 / 0.001, 0.83; Bmal1 KO 0.002 / 0.003, 0.36. The gate does not fail in the KO: the residual
    covariance passes the shape check although the amplitude has collapsed, so the gate is not a clean negative control
    here (clock share halves).
- **Reading:** five days after a 6 h advance, aortic SMC, fibroblasts and endothelial cells sit at the same partial shift
  (about −2.3 to −2.9 h): misaligned with the light cycle together, with no sign of desynchrony between them (difference
  0.02 h). Not counted as equivalence (power rule missed by 0.03 h). One tissue, one time after the shift, n = 1 pool per
  sex x ZT. The paper's intracellular finding (Dbp, Nr1d1 about 1.2 h closer to Arntl within SMC) is about gene order
  within the clock, not between cell types, and is not tested here.
- Code `paper/exploratory/A201_aorta_jetlag/scripts/score_a201.py` (about 2 min).
<!-- src: paper/exploratory/A201_aorta_jetlag/output/{power.csv,tests.csv,secondary.csv,gate.csv} -->
- **Post hoc (2026-11-03, user question: is a pre-existing offset lost?; `posthoc_baseline_offset.py`):** under aligned
  conditions SMC and fibroblast clocks are nearly in phase: SMC − fibroblast +0.18 h (jackknife SD 0.14; Tempo's 4 male
  libraries +0.27, the other 10 +0.08; females +0.17, males +0.17). Misaligned: +0.26 h (SD 0.43). There is no sizeable
  baseline offset to lose, so the aorta shows cell types that start in phase and stay in phase while shifting; it cannot
  show preservation of a non-zero offset. Per gene (no intervals, descriptive): most genes within 0.6 h in both conditions;
  Cry2 +2.51 → +5.04 h and Rora −7.58 → −10.17 h move most (both low-amplitude in fibroblasts, unreliable without
  intervals); Nr1d1 +0.96 → +1.10 h.
<!-- src: paper/exploratory/A201_aorta_jetlag/output/{posthoc_baseline_offset.csv,posthoc_gene_offsets.csv} -->

## A195 single-cell NC vs cycling clock-plane variance, scored 2026-11-02: not estimable (no clock-shaped signal in malignant NC pseudobulks)

- Data: A187 malignant pseudobulks, 8 studies. H195.1: 150 patients with ≥ 50 NC and ≥ 50 cycling cells (100,964 NC,
  71,785 cycling cells). H195.2: 165 patients in 7 studies, 58 low vs 55 high cycling-fraction tertile.
- **Signal check FAIL:** pooled gene coherence of the NC pseudobulks R 0.243 (170 patients), 42nd percentile of 200
  expression-matched pseudo-clock sets (q95 0.567). The NC malignant pseudobulks carry no clock-shaped signal above
  matched non-clock genes, so both tests are uninformative as registered.
- **Power:** H195.1 MDE 0.158 (≤ 0.40); H195.2 MDE 0.630 and the ALL reference 0.586 (> 0.40): the decisive tumour-level
  test would not have been estimable even with a signal.
- Tests, companions and artefact checks were not run (stop rules).
- **Reading:** the single-cell tumour atlases cannot say whether the proliferation clock-signal loss sits in dividing
  cells. This agrees with A187 (gates pass in 2 of 8 studies) and A188 (not estimable). Ordering clock loss and
  proliferation needs an intervention design (paired biopsies before and during an anti-proliferative treatment).
- Code `paper/exploratory/C4_single_cell/A195/scripts/score_a195.py --power-only` (about 4 min).
<!-- src: paper/exploratory/C4_single_cell/A195/output/{prep.csv,signal_check.csv,power.csv} -->

## A196 proliferation arrest, scored 2026-11-02: coherence does not recover when proliferation is stopped (H196.1 not estimable; post hoc calibration uninformative)

- **Treatment check (NeoPalAna, GSE93204, C1D1 → C1D15, palbociclib added to anastrozole, 23 pairs):** proliferation
  −0.54 SD (Wilcoxon p 0.005). The arrest acted.
- **H196.1 not estimable** (MDE 0.478 > 0.40). Descriptive: T −0.139 (p 0.48), 21.5th percentile of the matched pseudo
  sets; residualised on proliferation and purity −0.098. Coherence did not rise.
- Descriptive NeoPalAna: BL → C1D1 (anastrozole, proliferation −1.04 SD) +0.316 (p 0.037, 22 pairs); BL → C1D15
  (−1.59 SD) −0.050 (p 0.79); C1D15 → surgery −0.243 (15 pairs). Coherence per time point: 0.312, 0.676, 0.550, 0.415.
- **H196.2 (post hoc) uninformative:** within these ER+ cohorts the pre-treatment samples show no proliferation–coherence
  relation to reverse (pooled low − high tertile 0.026, one-sided p 0.41; GSE87411 alone 0.484, p 0.018; the other three
  0.127, −0.417, −0.395). Per series T − P: POETIC −0.063, GSE20181 +0.408, GSE87411 −0.122, GSE93204 +0.010; pooled
  −0.019 (z −0.23).
- **Pooled observed change (descriptive, 4 series, 319 pairs, proliferation down 0.54-0.92 SD in each):** +0.019
  (SD 0.062), so about [−0.10, +0.14].
- **Reading:** stopping proliferation for 2 weeks (aromatase inhibition) or adding a CDK4/6 arrest does not measurably
  restore clock coherence. A rise of the size of the SCAN-B proliferation contrast (about 0.3) is excluded, but a rise
  scaled to the 0.5-0.9 SD drops achieved (about 0.1-0.15 if the loss were fully and acutely reversible) is not. Within
  ER+ tumours the cross-sectional proliferation–coherence relation itself is not present in these cohorts, which limits
  what reversal can show. Nothing here supports the loss being an acute, reversible byproduct of ongoing division;
  it does not establish the opposite direction either.
- Code `paper/exploratory/A196_proliferation_arrest/` (about 40 s).
<!-- src: paper/exploratory/A196_proliferation_arrest/output/{power,tests,descriptive,calibration_posthoc}.csv -->

## Post hoc, descriptive (2026-11-03): GTEx Whole Blood has no core-clock time signal; composition explains little of its weak time signal

User observation stated before running: blood shows no coherent core clock, its time-of-day genes follow cell composition.
GTEx v10 fast deaths with recorded time of death; Whole Blood (212 donors) vs five solid tissues (179-220). Aggregates only.
- **Clock genes:** Whole Blood median cosinor R² 0.0096, the 52nd percentile of all expressed genes (solid tissues
  0.081-0.350, 99.5th-99.96th percentile). Phase concordance with the GTEx waveform 0.430 (solid 0.854-0.937). Adjusting
  for composition changes neither (0.0088; 0.482).
- **Time genes (top 100, chosen in one donor half, read in the other):** Whole Blood held-out R² 0.021 vs 0.014 for
  random genes, so blood carries little time-of-death signal at all. Composition removes 6 % of it. 20 % of these genes
  correlate |r| > 0.5 with a composition score (solid tissues 0-7.5 %).
- **Composition rhythms:** weak (neutrophil R² 0.025, peak 20.6 h; erythroid 0.025; T 0.003).
- **Reading:** the absence of a coherent core clock in blood is clear (no better than a random gene). The composition
  explanation is only weakly supported in post-mortem blood: blood time genes lean towards composition markers, but
  composition removes little, because there is little time signal to remove. Post-mortem blood may also differ from
  timed blood in living people.
- Code `paper/exploratory/blood_composition_posthoc/scripts/blood_composition.py` (seconds).
<!-- src: paper/exploratory/blood_composition_posthoc/output/blood_composition.csv -->

## A197 dedifferentiation, scored 2026-11-03: no coherence loss with the ES-like index beyond proliferation and subtype (H197.1-2 fail)

ES-like score: BENPORATH_ES_1 without proliferation, purity and clock genes (254 genes SCAN-B, 293 METABRIC). Strata:
PAM50 × proliferation tertile (15 per cohort); 2,000 within-stratum permutations; 200 matched pseudo-clock sets.

| test | T (high − low ES-like) | p one-sided | Holm | MDE | pseudo percentile |
|---|---|---|---|---|---|
| H197.1 SCAN-B (3,073) | +0.092 | 0.979 | 0.979 | 0.114 | 98.0 |
| H197.2 METABRIC (1,756) | −0.067 | 0.100 | 0.200 | 0.131 | 13.5 |
| H197.3 GSE72094 (442, secondary) | −0.090 | 0.100 | — | 0.169 (not estimable) | 10.0 |
| H197.4 GSE68465 (443, secondary) | +0.008 | 0.546 | — | 0.193 (not estimable) | — (13 clock genes) |

- **Both primaries fail at MDE ≤ 0.15.** SCAN-B points the opposite way. A drop of 0.114 (SCAN-B) or 0.131 (METABRIC)
  is excluded at 80 % power.
- **Descriptive:** without proliferation strata, SCAN-B T −0.137 and GSE72094 −0.165, so the ES-like association runs
  through proliferation. The index stays correlated with proliferation after removing proliferation genes in RNA-seq and
  LUAD arrays (Spearman 0.72 SCAN-B, 0.80 GSE72094) but not in METABRIC (0.11). ES-like vs ER: −0.19 (SCAN-B), −0.42
  (METABRIC). PRC2-target score: +0.191 (SCAN-B), +0.019 (METABRIC). Leave-one-gene-out ranges: SCAN-B +0.007 to
  +0.109, METABRIC −0.124 to +0.008.
- **Reading (as registered):** the ES-like index adds nothing beyond proliferation and subtype, so the "clock is a
  property of differentiated tissue" message rests on ER− (within PAM50 × proliferation, A166 H166.2) and on histology
  (adenocarcinoma vs squamous), and is stated as a hypothesis. Caveat: one transcriptional proxy of dedifferentiation.
- Code `paper/exploratory/A197_dedifferentiation/scripts/score_a197.py` (about 3 min).
<!-- src: paper/exploratory/A197_dedifferentiation/output/{power,tests,descriptive}.csv -->

## A202 (scored 2026-11-04): the A164 T2D interval is too narrow; the equivalences hold on the registered re-reading

Known changes of coupling injected into the real GTEx counts (`clocksync.power.damped`) at the A164 matched design, labels
exchangeable within the matched sets; module = 5 (vent) or 4 (fast) metabolic tissues, pooled with the fixed A164 weights.
Realised module truths at ±0.15: T2D ±0.150, obesity ±0.145 (vent obese liver, c0 0.090, damped to 0).
- **V202.2 (equivalence size, arm P, 1,000 replicates per cell): T2D FAILS, obesity passes.** Share of replicates declared
  equivalent with the registered A164 half-width at δ = −0.15 / +0.15: T2D 0.126 / 0.095 (limit 0.075), obesity 0.053 /
  0.030.
- **Cause: the T2D bootstrap SE is too small, not bias.** Median bias ≤ 0.006 (T2D) and ≤ 0.014 (obesity) in every cell.
  SD of the pooled estimate over replicates: T2D 0.0335-0.0344 against the A164 bootstrap SE 0.0247 (ratio 1.36-1.39)
  and null SD 0.0334; obesity 0.048-0.053 against 0.0518 and 0.0523. Coverage of the 90 % interval with the A164 SE:
  T2D 0.76-0.77, obesity 0.89-0.91; with SD_sim or the null SD, 0.88-0.90 in both designs.
- **V202.1 (coverage of the full procedure, arm B, 100 replicates with 200 bootstrap replicates each): T2D FAILS narrowly,
  obesity passes.** T2D 0.84 [0.76, 0.90] at δ 0 and 0.95 at −0.15; obesity 0.89 and 0.94 (limit 0.85). The replicate
  bootstrap SEs for T2D (median 0.032 at δ 0, 0.037 at −0.15) are larger than the SE A164 obtained on the observed
  labels (0.0247, below the 5th percentile of the δ-0 replicates, 0.0252): the observed SE was a small draw of an already
  slightly anti-conservative procedure.
- **Registered re-reading of A164** (a re-reading after A202, not a new test): half-width 1.645 × max(SD_sim over the arm P
  cells, A164 null SD), no bias correction (|bias| ≤ 0.02).
  - T2D: +0.019, 90 % [−0.037, +0.076] (A164: [−0.021, +0.060]). Inside ±0.15: **H164.1 holds under the corrected
    interval.**
  - Obesity: +0.004 [−0.083, +0.091] (A164: [−0.081, +0.089]). Inside ±0.15.
- **Attenuation sensitivity (descriptive, arithmetic on the A164 intervals):** divided by the A162 attenuation range 0.8 /
  0.9, T2D [−0.027, +0.075] / [−0.024, +0.067], obesity [−0.101, +0.112] / [−0.090, +0.099]: inside ±0.15. Thinning
  defines the truth on the estimator's scale, so this stays a stated limit.
- **Reading:** T2D and obesity leave metabolic-organ coupling within ±0.15, now with an interval validated at the design.
  Text: "+0.019 (90 % interval −0.037 to +0.076, design-based)". Methods: the matched-set bootstrap SE of a pooled d*
  can be too small at this design; report max(bootstrap SE, permutation-null SD) for equivalence intervals.
- Deviation: the first arm B run stopped in the obesity cells ("need at least 5 donors"); the replicate table was then
  built on the donors of every tissue's sets, as A164 (README). Runs: arm P 1,154 s, arm B 10,987 s + 6,331 s (16 workers).
- Code `paper/exploratory/A202_a164_coverage/scripts/{compute_a202,score_a202}.py`.
<!-- src: paper/exploratory/A202_a164_coverage/output/{arm_p,arm_b,summary,reread,attenuation,per_cohort,design}.csv -->

### 2026-11-04, integrity (found by the A203 acceptance check): ARNTL is zero in the A18-era per-cell files

- **Finding.** In the per-cell clock-gene files built by `common/datasets/cells_census/build.py` (legacy
  `desync/common/extract_cells.py`), the ARNTL column is 0 in every cell: hlca_full, islets_t2d, natri2024_ild,
  parkinson_regions, psychad_mssm, psychad_radc, seaad_dlpfc. The builder requested Census genes by symbol and the Census
  2025-11-08 names ARNTL "BMAL1", so the column was never filled and no error was raised. The A203 regrouping of the HLCA
  cells reproduces the A17 `lung_hlca` table exactly for the 13 other genes, the cell counts and the depths (1,321 of
  1,321 rows) and fails only on ARNTL.
- **Not affected:** every pseudobulk table (built from the Census groups, which map the symbol), hence every registered
  result from pseudobulks; the A152 per-cell files (their extractor maps BMAL1 → ARNTL; ARNTL totals 10,843 and 2,640);
  the cancer per-cell files (`cells/cancer`, ARNTL present).
- **Affected (per-cell readouts only, none in the claims ledger):** A18 part 2 per-cell amplitude and synchrony
  (`figure_2/scripts/backbone/compute_per_cell.py`, `figure_2/output/backbone/a18_cell_split*.csv`), A144 part 2
  (D144.2, descriptive, inconclusive, reads `a18_cell_split.csv`), and the Extended Data per-cell calibration
  (`extended_data/scripts/calibration/compute_orpheus_cells.py`, `orpheus_cells_*.csv`) where it reads these datasets.
  These readouts used 13 of the 14 clock genes; whether they change is not known until rerun.
- **Fix:** the builder maps Census symbols to the panel (BMAL1 → ARNTL, as A152) and fails if a panel gene is missing;
  `check_columns` compares cell ids and the 13 other genes with the stored file. HLCA is being rebuilt for A203
  (`A203_cells_vs_donors/scripts/rebuild_hlca_cells.sh`); the other six datasets and the reruns of the affected
  readouts wait for a decision (PI_QUESTIONS).

## A203 (scored 2026-11-04): cells do not saturate at 100 per pseudobulk in HLCA; not a Methods design rule (V203.1 and V203.2 each fail on one metric)

HLCA healthy lung, the A17 `lung_hlca` rows regrouped from per-cell counts (acceptance: 1,321 of 1,321 rows exact, after
the ARNTL rebuild), 198 donors, 14 readout compartments, 50 split halves per level. Median over compartments and splits
of the between-half difference:

| arm, level | offset (h) | c |
|---|---|---|
| cells 10 / 25 / 50 / 100 / 200 / all | 3.98 / 2.43 / 1.45 / 1.07 / 0.93 / 0.86 | 0.000 / 0.045 / 0.066 / 0.074 / 0.075 / 0.090 |
| donors 50 / 100 / 150 / 198 (all cells) | 4.02 / 2.07 / 1.30 / 0.90 | 0.033 / 0.057 / 0.077 / 0.091 |
| thinning to 25 / 100 cells' UMI | 2.31 / 1.08 | 0.049 / 0.070 |

- **V203.1 (cells saturate at 100): FAILS on the offset** (1.07 / 0.86 = 1.25 > 1.10); passes on c (0.82).
- **V203.2 (donors limit): passes on the offset** (2.07 / 0.90 = 2.29 ≥ 1.25); **fails on c** (0.63).
- **The c metric is not interpretable.** c is truncated at 0, and the share of half pairs with a zero c rises as
  information falls (cells 10: 0.88, 25: 0.51, 100: 0.15, all: 0.06; donors 50: 0.23, 198: 0.06; `boundary.csv`). Less
  information pulls both halves to 0 and shrinks |c₁ − c₂|, so err_c decreases with fewer cells and donors. Its pass
  (V203.1) and its fail (V203.2) both come from the truncation. The registered metric was not fit for c (a design flaw,
  stated, not re-scored).
- **Offsets (interpretable):** both cells and donors limit precision at HLCA depth. Error falls steadily up to 200 cells
  and is still 8 % above the all-cells value there (m₉₅ not reached among the tested levels). Halving donors more than
  doubles the error (2.29 for 198 → 100; √2 expected for pure donor sampling, the excess plausibly from less stable host
  fits in 50-donor halves). Per compartment the low-coupled lymphocytes (B, NK/ILC) stay at 2.6-5.8 h at every cell level.
- **Cell heterogeneity adds little beyond counting noise:** binomial thinning to the UMI of m cells gives nearly the same
  offset error as subsampling m real cells (25: 2.31 vs 2.43 h; 100: 1.08 vs 1.07 h).
- **Reading:** not a Methods design rule (registered; user, 2026-11-04). Descriptive answer to the question: at current
  single-cell depth, clock-gene counts per pseudobulk still limit offset precision below a few hundred cells, so deeper
  sequencing or more cells per donor would help as well as more donors. Per-cell clock signal behaves like Poisson
  counts, so pseudobulks lose little relative to a cell-level model. One atlas.
- Runs: HLCA per-cell rebuild 1,197 s; split fits 983 s (16 workers).
- Code `paper/exploratory/A203_cells_vs_donors/scripts/{compute_a203,score_a203}.py`.
<!-- src: paper/exploratory/A203_cells_vs_donors/output/{acceptance,errors,summary,m95,families,boundary,compartments}.csv -->
- **Decision (user, 2026-11-04, PI_QUESTIONS Q23):** the healthy-atlas per-cell readouts (A18 part 2, A144 part 2, the
  per-cell calibration; plan ED item 8) are dropped; the six other per-cell files are not rebuilt (HLCA was rebuilt for
  A203). The malignant per-cell D of WP0 item 2 stays: the cancer per-cell files have ARNTL and malignant cells carry
  about 4.3 clock UMI per cell, where D separates coherent from scrambled gene phases. A203 bears on donor-level
  precision only, not on within-sample readouts.

## A204 (run 2026-11-05/06): simulator gate fails on TimeTeller, so A204 is not scored; COFE descriptives are reported, labelled as such
Registered gate (arm F, written before any readout): COFE |circular correlation| of tumour labels with true times 0.859
(>= 0.7, passes); TimeTeller median circular error 5.78 h (<= 3 h, **fails**). As registered, A204 stops: H204.1-H204.4 are
not scored. TimeTeller put 99-100 % of predictions inside 07-20 h even with whole-day sampling (`tt_cells.csv`
share_7_20), so it does not work on this simulator (likely the healthy-training vs tumour mesor mismatch or our 3 h
binning; not investigated). Nothing is concluded about the REMAGUS r ≈ 0.
**Descriptive only (not verdicts; the would-be readouts are in `tests.csv`):** COFE on 14 clock genes, 300 tumours,
about 60 matched normals, 30 replicates per cell:
| window, truth | true times: R / largest gap (h) | COFE labels: R / largest gap (h) | Rayleigh p > 0.05 on differences | clocksync pairing p < 0.05 |
|---|---|---|---|---|
| 08-18 h, c 0.9 | 0.715 / 11.03 | 0.043 / 0.83 | 18 of 30 | 30 of 30 |
| 08-18 h, c 0.5 | 0.712 / 10.95 | 0.023 / 0.57 | 22 of 30 | 30 of 30 |
| 08-18 h, uncoupled | 0.713 / 10.76 | 0.018 / 0.60 | 27 of 30 | 1 of 30 |
| whole day, c 0.9 | 0.054 / 0.48 | 0.024 / 0.52 | 13 of 30 | 30 of 30 |
- COFE spreads 08-18 h labels over the whole cycle (largest empty arc 0.83 h vs 11.0 h true) and its CV minimum is
  interior in every run (its validity signal). Its tumour-normal differences separate coupled (c 0.9) from uncoupled
  tumours with AUC 0.666 [0.518, 0.804] under working hours (0.928 whole day); clocksync's pairing correlation 1.000 in
  both, offset 3.05 h (truth 3 h). Even with whole-day sampling, COFE differences pass as uniform in 13 of 30 coupled
  replicates at about 60 pairs.
- Reading (descriptive): consistent with the mechanism stated in the prompt (full-cycle labels and weakly informative
  tumour-normal differences under working hours), but registered only conditional on a gate that failed; quote as a
  descriptive simulation, not a test. Genome-wide COFE was not simulated.
<!-- src: paper/exploratory/A204_orderer_uncoupling/output/{gate,tests,cofe_cells,auc,tt_cells}.csv -->

## A205 (scored 2026-11-06): under Saini-type partial shifts the gate passes (V205.1 PASS) but the offset is not bounded by the gene shifts (V205.2 FAIL)
Target compartment with ARNTL/CLOCK/NPAS2/CRY1 +6 h and the other 10 genes +12 h; 4 intact anchors; n 80; 100 replicates.
- **V205.1 PASS:** gate passes in 100 of 100 (08-18 h; binomial lower bound 0.96); also 100 % in every cell.
- **V205.2 FAIL:** offset in [6, 12] h in 1 of 100; circular mean 12.86 h (08-18 h), 14.26 h (whole day), i.e. beyond the
  larger gene shift. Caveat on the criterion: a coherent 12 h shift (P12) lands in [6, 12] in only 46-47 % because
  estimates scatter around the 12 h edge; the whole-day value 14.26 h is nevertheless clearly outside.
- **Descriptive:** incoherence reads as a weaker clock: coupling 0.695 vs 0.836 (P0, W), log relative amplitude -0.39
  (P12: -0.01); clock share 0.931 vs 0.976. With +3 h / +12 h: coupling 0.262, log amplitude -1.35, offset 11.28 h (W).
- **Reading:** the gate is a clock test, not a waveform test: a compartment whose loops shift unequally passes it and
  returns an offset that is not an average of the gene shifts (here 1-2 h beyond the larger one) together with lower |z|
  and c. A measured offset or amplitude drop therefore cannot be attributed to one loop, and a lower |z| alone does not
  distinguish a damped clock from an internally incoherent one. To add to METHOD.md assumptions (conserved waveform).
<!-- src: paper/exploratory/A205_partial_shift_gate/output/{tests,cells,per_rep}.csv -->

## A206 (scored 2026-11-06): a clockless liver with a feeding rhythm does not pass the single-set gate (H206.1 PASS, borderline)
GSE13093 liver, GTEx human waveform, `paired.reference_gate` (200 permutations, seed 206), 12 samples per set:
| set | clock share (null mean) | p |
|---|---|---|
| WT restricted feeding, even CTs | 0.692 (0.109) | 0.005 |
| WT restricted feeding, odd CTs | 0.742 (0.118) | 0.005 |
| Cry1/2 DKO restricted feeding | 0.158 (0.072) | 0.0597 |
| Cry1/2 DKO ad libitum | 0.182 (0.101) | 0.110 |
- Power companion passes (both WT halves). **H206.1 PASS**, but p 0.0597 is close to the threshold and the DKO share
  is still above its null mean: feeding leaves a weak, partly clock-like co-variation (cosinor amplitudes log2 0.01-0.66
  in DKO RF vs 0.18-3.78 in WT).
- Reading: in a single tissue, a feeding-driven rhythm without a clock is not read as a clock at n = 12, unlike the
  cross-compartment gate in Bmal1 knockouts (A201, A129), which kept passing. With the A201 result: the clock share,
  not the pass/fail, is the informative readout of clock loss. Caveat: CRY1/CRY2 are the knocked-out genes, so 2 of
  14 genes carry no signal by construction.
<!-- src: paper/exploratory/A206_cry_dko_gate/output/{gate,cosinor}.csv -->

## A207 (scored 2026-11-06): co-expression and order "disruption" are present in the same TCGA tumours that follow their normal's clock (H207.1, H207.2 PASS)
517 pairs, 10 cancers, 1,000 within-patient swaps (seed 207), Holm over two. Implemented by a Sonnet agent; numbers
checked against the files.
- **H207.1 PASS:** pooled ΔCCD (tumour − normal, Shilts reference, 12 genes) 0.890 (null mean 0.001), p 0.001, Holm
  0.002. Per cancer -0.09 (BLCA) to 2.14 (UCEC); 9 of 10 positive. Replicates Shilts 2018 on our table.
- **H207.2 PASS:** Grabe procedure (normals' PCA phase on the 14 genes, cosinor R² tumour − normal, median over genes)
  pooled -0.413; swap p 0.001, Holm 0.002; companion: 0 of 200 expression-matched pseudo sets at or below it (pseudo
  5th percentile -0.360, mean -0.305). **Caveat:** the normals define the phase in-sample, so ΔR² is negative for any
  gene set (all 200 pseudo sets have swap p ≈ 0.001); the pass rests on the companion only, and part of the clock
  excess may come from the clock genes having more normal-sample structure to lose.
- **Descriptive (n = 10 cancers):** Spearman of ΔCCD with the A138 H138.4 pairing corr -0.782 (p 0.0075); ΔR² with
  corr 0.152 (p 0.68). Cancers whose tumours keep the strongest pairing (LUAD, LIHC, BLCA, PRAD, KIRC: corr 0.42-0.52,
  p <= 0.005) have positive ΔCCD in 4 of 5 (0.07-1.33).
- **Reading (registered):** with H138.4, co-expression and order disruption coexist with tumour-normal clock pairing in
  the same cohorts, so "disrupted" (ΔCCD, Grabe R²) must not be read as "uncoupled". Post hoc observation, not a test:
  across cancers ΔCCD tracks weaker pairing (rho -0.78), so CCD and pairing share information at the cohort level.
<!-- src: paper/exploratory/A207_ccd_grabe_tcga/output/{tests,per_cancer,pseudo}.csv -->

## A204b (scored 2026-11-06, post hoc to A204's failed gate): with matched tissue means TimeTeller works, and a coupled tumour would correlate with sampling time at REMAGUS size (H204b.4 PASS)
- **Gate passes:** whole-day median error 1.23 h (A204, with independent per-gene tissue offsets: 5.78 h).
- **H204b.4 PASS:** 08-18 h, c_T 0.9, n 108: slope p < 0.05 in 30 of 30 (median slope 0.81, r 0.57). Power by c_T
  (08-18 h): 0.07 (uncoupled), 0.70 (0.25), 0.87 (0.5), 0.97 (0.75), 1.00 (0.9); MDE c_T 0.5. Narrow window 10-16 h
  (sensitivity): 0.03 / 0.67 / 0.77 / 0.97 / 0.97, MDE 0.75. Whole-day linear slopes are not meaningful (circular
  time) and not read.
- **Reading (registered for H204.4, applied here):** if TimeTeller times are as accurate as in this simulator, a tumour
  coupled at c >= 0.5 correlates with sampling time at n = 108 over a 08-18 h window; the REMAGUS r ≈ 0 (slope -0.027,
  p 0.80) is then not explained by the narrow window alone. Combined with A204 (post hoc, not a test): when the
  training tissue's per-gene mean levels differ from the test tissue's (sd 0.5 in log), TimeTeller's times are biased
  into 07-20 h and lose the link to time. Healthy-tissue training applied to tumours is that situation, so a mean-level
  mismatch is an alternative to tumour uncoupling that the REMAGUS analysis does not exclude. The window range of
  REMAGUS is not reported; both windows are assumptions.
<!-- src: paper/exploratory/A204_orderer_uncoupling/output/a204b/{gate,tests,tt_cells}.csv; A204 output/tt_cells.csv -->

## A208 (scored 2026-11-06): women and men have the same organ coupling and offsets in GTEx (H208.1 PASS, equivalence; H208.2 FAIL, no offset difference)
A164 procedure with sex as the label (implemented by a Sonnet agent; numbers checked against the files). 55 tissue x
cohort units (vent 24, fast 31), 2,981 matched pairs, 1,000 permutations, 200 bootstraps with host refit; interval
± 1.645 x max(bootstrap SE, null SD).
| | d* (female − male), 90 % | p | offset (h), 90 % | p |
|---|---|---|---|---|
| ventilated | +0.021 [-0.011, +0.053] | 0.19 | -0.058 [-0.270, +0.153] | 0.52 |
| fast | +0.018 [-0.015, +0.051] | 0.36 | +0.154 [-0.024, +0.333] | 0.11 |
| pooled | +0.020 [-0.003, +0.043] | 0.11 | +0.043 [-0.104, +0.189] | 0.50 |
- **H208.1 PASS:** pooled d* inside ±0.15 (MDE 0.032). Companion (3 lowest-coupling tissues by male c per cohort):
  +0.010, p 0.75.
- **H208.2 FAIL:** no sex offset (p 0.50); equivalent within ±1 h (and within about ±0.2 h). The pseudo-clock companion
  (50 sets) has |offset| q95 1.29 h, so no clock-specific offset. Same sign as the archive A99 timed +0.35 h, but about
  zero.
- **Reading:** organ coupling and organ phase offsets do not differ by sex in GTEx (coupling within ±0.04, offsets within
  ±0.2 h). A sex x dosing-time interaction in chronotherapy (Lévi 2010) would therefore not come from a sex difference in
  how organs are timed relative to the body's clock; host phase relative to wall time (A99 type) is a separate question
  and is not addressed (untimed readout).
<!-- src: paper/exploratory/A208_sex_coupling/output/{tests,pooled,companions,power}.csv -->

## A211 (scored 2026-11-06): regenerating liver keeps a coherent clock; synchronous proliferation does not swamp it (H211.1 FAIL, powered)
Data GSE180012 (`$SCCC_DATA/gse180012_phx/`): 24 PHx and 24 sham livers, 3 per time point, 4-32 h after surgery at
08:00-09:00; log2(RPKM + 1).
- **Gate passes:** G2M PHx − sham over 20-32 h +1.05 z, p 0.0005 (−1.11 at 4 h to +1.39 at 32 h); sham clock share 0.673
  vs 0.077, p 0.005.
- **H211.1 FAIL (powered):** coherence PHx 0.926 vs sham 0.953, T −0.027, p 0.153 (2,000 swaps within time point), MDE
  0.064; pseudo-clock 41st percentile (q05 −0.371).
- **Descriptive:** 4-16 h T −0.017 (p 0.38, MDE 0.137); 20-32 h (replicative) T −0.024 (p 0.073, MDE 0.040); neither
  near the pseudo q05. Cosinor over the 28 h series (24 h period, descriptive): median clock-gene amplitude 0.68 (PHx)
  vs 0.89 log2 (sham), median R² 0.40 vs 0.46; median acrophase difference −0.7 h (PER1 −3.6 h, CLOCK −3.2 h).
- **Reading (registered):** a powered fail with both gates passed: dividing hepatocytes keep their clock co-expression
  within 0.064 of sham, against tumour losses of −0.12 to −0.34 (A127, A166, A171). With A126 (endometrium), proliferation
  alone does not reproduce the tumour coherence loss; a tumour-state explanation is favoured. Limits: one mouse dataset;
  the G2M rise is relative within this design and not on the tumour tertile scale; PHx is also an injury; amplitude falls
  moderately (descriptive), so a damped but coherent clock is not excluded.
<!-- src: paper/exploratory/A211_liver_regeneration/output/{gate,power,tests,cosinor}.csv -->

## A210 (scored 2026-11-06): MYC amplification at matched proliferation: no loss in METABRIC (H210.1 FAIL, powered), clock-specific loss in TCGA (H210.2 PASS); not carried by MYC expression
- **H210.1 FAIL (METABRIC, primary):** 14 PAM50 x proliferation cells, 1,294 tumours (450 amplified). T (amplified −
  neutral coherence) +0.038; null mean −0.020, SD 0.052; p 0.86; MDE 0.130; pseudo 70th percentile.
- **H210.2 PASS (TCGA, secondary):** 30 project x proliferation cells in 13 projects, 2,728 tumours (844 amplified).
  T −0.129; null mean −0.044, SD 0.034 (so about −0.085 beyond the group-size bias); p 0.005, Holm 0.010, MDE 0.086;
  pseudo 0.5th percentile (q05 −0.064): clock-specific. 25 of 30 cells and 11 of 13 projects negative (HNSC −0.290, BRCA
  −0.234, LUAD −0.234; OV +0.008, ESCA +0.024).
- **Descriptive:** residual proliferation imbalance within cells, amplified − neutral: +0.045 z (METABRIC), +0.101 z
  (TCGA). Gain − neutral −0.033 / −0.049. MYC expression top − bottom tertile within cells +0.003 / +0.047: no loss.
- **Reading:** the registered primary fails at a powered design, and TCGA BRCA (−0.234) disagrees with METABRIC (+0.038).
  The TCGA loss is clock-specific and broad across tumour types, but MYC expression does not show it, so it is not a MYC
  activity effect; residual proliferation (+0.10 z) or copy-number burden correlated with MYC amplification (cf. A134)
  are not excluded. MYC amplification is not established as a cause of coherence loss beyond proliferation.
- **Deviation:** the tests stage first stopped on a helper name (`wf_for` is in `a166_lib`, not `a168_lib`); fixed
  before any statistic was printed.
<!-- src: paper/exploratory/A210_myc_coherence/output/{power,tests,descriptive,descriptive_cells,counts}.csv -->

## A209 (scored 2026-11-06): gate plus matched companion makes no false claim, but the no-clock arms were not exchangeable (see A209b)
- **H209.1 PASS (as registered):** arm N, gate-only false claims 26 of 150 (0.173; A145 had 0.127), with the matched
  companion 0 of 150 (p 1.4e-9 against 0.127).
- **H209.2 PASS but vacuous:** in N_resp the gate failed in all 30 replicates (the 2x responsive factor breaks the
  clock-shape check), so the companion was never needed.
- **H209.3 PASS:** arm W, targets with c >= 0.5: claims 86 of 90 with the companion vs 90 of 90 without (ratio 0.956);
  by true c: 0.25 → 0.167 (gate-only 0.733), 0.5 → 0.867, 0.75 and 1.0 → 1.00; no false claim at c 0.
- **Caveat found after scoring:** with clock amplitude 0, clock genes still received compartment clock-plane noise
  through their waveform amplitudes; the amplitude-0 pseudo genes did not, so pseudo sets read more coupled (median
  clock percentile 0 in N and N_resp) and H209.1 is strict by construction. A209b repeats N and N_resp with exchangeable
  pseudo genes. H209.3 (power) is not affected.
<!-- src: paper/exploratory/A209_gate_companion/output/{tests,rates,per_target}.csv -->

## A209b (scored 2026-11-06): with clock-like non-clock genes, still no false claim; the exchangeability check fails for a structural reason
- **H209b.1 PASS:** N_b, gate-only false claims 13 of 150 (0.087), with the matched companion 0 of 150 (p 1.4e-9 vs
  0.127, Holm 1.4e-9).
- **H209b.2 PASS:** N_resp_b + N_resp1_b, 0 of 300 (p 2.0e-18, Holm 4.0e-18); the gate failed in all 60 replicates, also
  with the 1x responsive loading.
- **Exchangeability check (registered descriptive) fails:** median clock-gene percentile among matched sets 0.010 (N_b),
  0.020 (N_resp1_b), 0.115 (N_resp_b), not about 0.5.
- **Reading:** the A145 no-clock arm is not clock-free: every compartment keeps its own clock-shaped state (psi, set from
  the target c), aligned with the clock waveform by definition, so it lowers the clock genes' shared fraction and
  cannot be given to pseudo genes, whose phases do not match the clock projector. What the arm does show: false coupling
  produced by a non-clock factor shared by all genes is carried at least as strongly by matched pseudo genes, and the
  companion removes every such claim; a factor biased to responsive clock genes breaks the gate (0 of 60). Not tested,
  and not separable by any method: a shared non-clock factor with the clock genes' own phase pattern.
<!-- src: paper/exploratory/A209_gate_companion/output/a209b/{tests,rates,per_target}.csv -->

## A212 (scored 2026-11-06): bulk statistics reject swamping for the proliferation loss but cannot tell damping from rewiring (V212.1 FAIL, narrow; V212.2 PASS)
- **Design (registered, plasmode):** real reference-arm tumours (B1 SCAN-B H127.1, 1,030 tumours in 15 strata; B2 METABRIC
  ER+ H166.1, 1,498 in 15), random halves, one half given a mechanism at the strength whose mean dR equals the observed
  dR (B1 −0.123, B2 −0.229): D damping / intra-tumour dephasing, S added non-clock factors, W shared rewiring, H
  per-tumour rewiring. Matched strengths B1: rho 0.83, S 0.72, W 1.8 h, H 2.2 h; B2: rho 0.83, H 2.5 h; S and W never
  reach the B2 loss (best −0.197, −0.207; run at their strongest value).
- **V212.1 FAIL (narrow):** AUC on (dlog_plane, dlog_off), W vs D 0.894, H vs D 0.909; bar 0.9 for both. B2: 0.889 and
  0.819. Rewired relations and a damped or dephased clock are not reliably separable in bulk.
- **V212.2 PASS:** swamping vs damping AUC 1.000 on both bases (W, H vs S also 1.000). Added non-clock factors raise
  off-plane variance strongly (and in-plane variance, through random loadings).
- **Observed points (descriptive, seen before registration):** SCAN-B proliferation (dlog_plane −0.228, dlog_off
  +0.020) lies outside the swamping cloud (Mahalanobis 11.6, beyond all 200 replicates) and inside D (73rd percentile),
  W (22nd) and H (36th). METABRIC ER (−0.153, +0.118) lies outside D and S (100th) and closest to W (77.5th, but W cannot
  reach the loss) and H (93.5th): no single mechanism fits; mixtures were not fitted.
- **Reading:** for the proliferation loss, "swamped" is excluded. The clock-plane variance falls without added
  non-clock variance, which is what a damped or dephased clock (lower per-cell amplitude, or clones or cells out of
  phase) or mildly rewired gene relations produce; bulk statistics cannot decide between these. Clonal heterogeneity
  fits either reading, so A213 tests it directly. The ER loss needs a mixture (plane fall plus an off-plane rise that
  A176 found not clock-specific).
- **Property of the statistic (README):** dR is not monotone in D, W and H strength (coherence rises again at rho <= 0.4
  or shifts >= 4-5 h); matched strengths lie on the monotone branch.
<!-- src: paper/exploratory/A212_coherence_mechanisms/output/{matched,tests,auc_all,observed,grid}.csv -->

## A213 (scored 2026-11-06): subclonal heterogeneity does not lower coherence (H213.1 FAIL, powered); copy-number burden does in TCGA (H213.2 PASS, clock-specific), not in METABRIC
- **Power note:** `registered_experiments.md` (commit e8b2f2e). ABSOLUTE subclonal genome fraction (PanCanAtlas,
  sha256 in `$SCCC_DATA/tcga/absolute_pancanatlas/PROVENANCE.md`); FGA from the A210 GISTIC matrices.
- **H213.1 FAIL (primary, powered):** TCGA, 69 project x proliferation cells in 23 projects, 5,499 tumours; top − bottom
  subclonal-fraction tertile T +0.021 (null −0.000, SD 0.023), p 0.82, MDE 0.058; pseudo 84th percentile. 18.7 % of
  called tumours have subclonal fraction 0 (ties split by table order).
- **H213.2 PASS, clock-specific (TCGA, FGA):** 6,024 tumours; T −0.063 (null −0.000, SD 0.023), p 0.003, Holm 0.009,
  MDE 0.058; pseudo 2nd percentile (q05 −0.047). Negative in 15 of 23 projects (PAAD −0.311, BRCA −0.290, LIHC −0.228).
- **H213.3 FAIL (METABRIC, FGA):** 15 cells, 1,176 tumours; T −0.032 (null 0.000, SD 0.055), p 0.28, Holm 0.56, MDE 0.136.
- **S213.4:** FGA-matched A210 MYC contrast (453 pairs in 25 cells, mean |FGA gap| 0.008): excess −0.065 vs the A210
  excess −0.085 (ratio 0.77), p 0.073. By the registered rule FGA does not carry the A210 result.
- **Descriptive:** Spearman subclonal ~ FGA 0.52 within project; FGA ~ proliferation 0.21. Residual imbalance top −
  bottom within cells: proliferation +0.076 z (H213.1), +0.142 z (H213.2), +0.021 z (H213.3); purity −0.061 / −0.002.
  Components (top − bottom): H213.1 dlog_plane −0.109, dlog_off −0.047; H213.2 dlog_plane −0.010, dlog_off +0.125,
  ds −0.136; H213.3 −0.106, −0.021.
- **Reading:** heterogeneity measured as subclonal genome fraction does not lower coherence at matched proliferation
  (bounded at 0.058), so clonal heterogeneity is not supported as a cause of the coherence results. Copy-number burden
  goes with a clock-specific loss in TCGA, but not in METABRIC (lower power), and the FGA arms keep a residual
  proliferation imbalance (+0.14 z). The FGA loss has the off-plane-rise pattern (in-plane variance unchanged), unlike
  the proliferation loss of A176/A212 (in-plane fall, off-plane flat). The A210 MYC result is not explained by FGA.
  Both TCGA passes (A210, A213 H213.2) are secondary: no METABRIC replication.
<!-- src: paper/exploratory/A213_heterogeneity_coherence/output/{power,tests,s213_4,descriptive,descriptive_spearman,descriptive_per_project}.csv -->

## A214 scout (2026-11-06, step 1 of the pooled single-cell round, not registered): per-cell questions not estimable with public tumour data
Prompt `paper/docs/prompts/sc_pooled_round.md`. No malignant-vs-normal D difference, amplitude or clone phase was computed.
- **Inventory.** On disk, 168 patients have >= 100 malignant and >= 100 normal epithelial cells (108 with adjacent-normal
  epithelium in 24 studies, 115 with in-tumour non-malignant epithelium). Clock UMI per cell (median of per-patient
  medians): malignant 2.0, adjacent normal 1.0, in-tumour normal 0.5; with the METHOD.md rule (median >= 2) on both
  groups, 38 patients remain. 3CA unreachable; its Zenodo export lacks patient, malignant and gene labels.
- **Effect sizes (exact von Mises model).** Per-gene acrophase shifts of SD 2 h lower D to 0.74 of the intact value
  (1 / 3 / 4 h: 0.92 / 0.51 / 0.30; uniform 0.04), over lam 1.5-2.5 and R 0.25-0.7: the 2 h equivalence bound is
  −0.26 D (about −1.0 at D 3.9). Damping vs dephasing at a matched pseudobulk amplitude ratio 0.5 / 0.7 / 0.85 differ by
  0.52-0.83 / 0.39-0.60 / 0.22-0.34 D, and that ratio is not identified without the malignant host coupling (H138.1).
- **Noise and power.** Split-half D per group (mean 0 by construction) gives per-group SE (rms) malignant 17.8, adjacent
  normal 25.2, in-tumour normal 23.7. MDE of the pooled paired D difference: 5.1-5.5 (108 patients), 7.3-7.9 (115),
  5.2-5.4 (168), 4.7-5.9 (38 under the UMI rule), lower bounds (between-patient heterogeneity excluded). Even full
  scrambling (about −3.8) is below the MDE. SE scales as clock UMI^−0.61 (SE 0.06) x cells^−0.41 (SE 0.06): MDE 1.0 needs
  about 233 patients at 8 and 100 at 16 clock UMI per cell in both compartments (x 4 with the present tails).
- **Deep-study scout** (plate chemistries included, user 2026-11-06): no study reaches 8 clock UMI per cell in both
  compartments. Moorman 2025 (CELLxGENE 829a3cd1, 10x 3' v3, 14 patients, about 15,000 UMI per cell) has median 2 clock
  UMI per malignant cell (mean 3.28) and 1-2 per normal epithelial cell (1.0-1.5 per 10,000 UMI). No Smart-seq3 /
  FLASH-seq / VASA-seq / Parse tumour study with malignant and normal epithelium was found.
- **Reading.** Rewired vs intact gene relations, damping vs dephasing, and their dependence on proliferation are not
  estimable from public single-cell tumour data; the limit is the depth of the normal epithelium. Clones: no clone calls
  on disk. A powered answer needs a deep dataset (>= 8-16 clock UMI per cell in both compartments, 100+ patients).
<!-- src: paper/exploratory/A214_sc_pooled_scout/{inventory/inventory_studies.csv,inventory/inventory_patients.csv,inventory/inventory_deep.csv,output/power/model_rewiring.csv,output/power/model_damp_dephase.csv,output/power/noise_groups.csv,output/power/power.csv,output/power/depth_scaling.csv} -->

## Erratum (2026-11-06) to WP0 item 2 (third part): "per-cell D sees no loss" is not informative
RP line 325 reads the equal per-cell D (malignant 3.89, normal 3.90, 79 patients) as "no loss". At real depth the paired
D difference has an MDE of about 5 (A214 scout, 108-168 patients; larger at 79), so equality is uninformative: a 2 h
rewiring (−1.0) and even full scrambling (−3.8) would go undetected in the paired contrast. What stands: the one-sample
malignant D (pooled SE about 1.2 over 228 groups) lies well above the scrambled-phase value (0.42 in the calibrated
simulation), so fully scrambled gene phases within malignant cells remain rejected. The WP0 simulation SDs (0.88 for
coherent malignant D) understate the real noise of the paired difference; the reconciliation "clock intact within cells"
(REVIEW 2026-10-28 Q8, STATUS 2026-10-27) should read "not scrambled; finer changes not estimable".
<!-- src: paper/exploratory/A214_sc_pooled_scout/output/power/{power,noise_groups}.csv; paper/exploratory/WP0_integrity/output/{malignant_downsample,incoherence_sim_calibrated}.csv -->

## A215 / A216 (scored 2026-11-06): TP53-altered tumours carry less clock signal at matched proliferation (H215.1 PASS, clock-specific); proliferation does not lower coherence more where p53 is lost (H215.4 FAIL); A216 not scored (precondition)
Power note: `registered_experiments.md` (commit fec75b7). TCGA MC3 public MAF (sha256 in `$SCCC_DATA/tcga/mc3/PROVENANCE.md`);
METABRIC mutation panel from the stored cBioPortal tarball. TP53 altered in 41.4 % (TCGA) and 32.9 % (METABRIC) of
tumours with a status; RB pathway 23.1 % and 5.8 %.
- **H215.1 PASS, clock-specific (TCGA TP53, primary):** 46 project x proliferation cells in 18 projects, 6,449 tumours;
  T −0.085 (null −0.013, SD 0.026; excess −0.072); p 0.004, Holm 0.012; MDE 0.066; below all 200 pseudo-clock sets (q05
  −0.045). 14 of 18 projects negative (KIRC −0.436, CESC −0.380, STAD −0.247, LUSC −0.223, LIHC −0.214, BRCA −0.209; PAAD
  +0.280, ESCA +0.351). LUAD (−0.132, 3 cells) was not blind (A171).
- **S215.5 (proliferation-matched pairs, 3,652 tumours in 44 cells):** −0.065, swap-null p 0.0115. Within cells the altered
  arm is +0.13 z more proliferative; matching on proliferation leaves the loss (−0.065 vs excess −0.072).
- **H215.2 FAIL (METABRIC TP53):** 12 cells, 1,519 tumours; T −0.106 (null −0.041; excess −0.065), p 0.124, Holm 0.248,
  MDE 0.139; pseudo 7th percentile. Same direction and size as TCGA, below this cohort's MDE.
- **H215.3 FAIL (TCGA RB pathway):** 49 cells, 6,586 tumours; T −0.050 (null −0.028; excess −0.022), p 0.203, Holm 0.248,
  MDE 0.068.
- **H215.4 FAIL, opposite sign (interaction):** 15 projects, 6,272 tumours; I +0.078 (null −0.051; excess +0.129), p 0.98,
  MDE 0.155; pseudo 89.5th percentile. The proliferation contrast is not steeper in TP53-altered tumours.
- **H216.1-2 not scored (registered precondition):** on the eligible tumours, P0 (proliferation contrast within project
  only) is not below its pseudo 5th percentile: TCGA −0.077 vs q05 −0.080 (marginal), METABRIC (within PAM50) +0.036 vs
  −0.128. Descriptive only: P1 at matched project x FGA x TP53 −0.030 (TCGA, pseudo 13.5th) and +0.010 (METABRIC);
  attenuation A 0.38 [−0.08, 1.30] (TCGA; bootstrap with replacement, which can bias coherence as in A149) and
  uninformative in METABRIC.
- **Descriptive:** H215.1 components dlog_plane −0.107, dlog_off +0.002 (in-plane fall, off-plane flat: the
  proliferation pattern of A176 / A212, unlike the FGA off-plane rise of A213). Truncating vs missense: TCGA −0.070 /
  −0.059, METABRIC −0.160 / −0.195.
- **Reading (registered branch "H215.1 passes, H215.4 fails"):** TP53 loss is a tumour-state marker with its own
  clock-specific loss of in-plane signal at matched proliferation, not a modifier of the proliferation effect. The
  checkpoint hypothesis in its strong form (proliferation lowers coherence only where checkpoints are lost) is not
  supported; RB-pathway loss shows no loss within MDE 0.068. With A213 (FGA, TCGA only) and A210 (MYC, TCGA only), TP53
  is the third genomic state with a clock-specific TCGA loss, and the first whose METABRIC estimate has the same size
  (excess −0.065 vs −0.072) though below that cohort's MDE. A216 cannot say whether proliferation is a marker: on these
  tumours the within-project proliferation loss itself is not clock-specific at the 5th percentile (TCGA, marginal).
<!-- src: paper/exploratory/A215_checkpoint_coherence/output/{power,tests,attenuation,descriptive,descriptive_per_project,descriptive_tp53_class,status_fractions}.csv -->

## POST HOC (2026-11-06, user question): the TCGA-only genomic losses in breast are carried by molecular subtype; qualification of A210 H210.2, A213 H213.2 and A215 H215.1
Question: are METABRIC and TCGA discordant because of cancer type? Within breast alone, TCGA BRCA was the more negative
in all three (MYC −0.234 vs METABRIC +0.038; FGA −0.290 vs −0.032; TP53 −0.209 vs −0.106), so cancer type is not the
reason. The registered TCGA cells are project x proliferation tertile; METABRIC cells are PAM50 x proliferation tertile.
TCGA BRCA recomputed both ways (descriptive, 2,000 within-cell permutations, no pseudo companion):

| contrast | cells: proliferation only (excess, p) | cells: PAM50 x proliferation (excess, p) |
|---|---|---|
| MYC amplified vs neutral | −0.205, 0.011 (603 tumours) | −0.084, 0.154 (387) |
| FGA top vs bottom tertile | −0.291, 0.001 (718) | −0.045, 0.228 (580) |
| TP53 altered vs intact | −0.191, 0.018 (1,006) | −0.008, 0.470 (696) |

The arms differ strongly in subtype (pooled over tertiles): Basal 36 / 34 / 44 % in the altered or top arm vs 4 / 2 / 3 %
in the reference arm, LumA 22 / 23 / 16 % vs 74 / 80 / 70 % (MYC / FGA / TP53). At matched subtype the BRCA losses shrink
by 59-96 % and match METABRIC. Single cancer type, post hoc, per-arm sizes 6-9 cells: not a powered test of "no loss".
**Consequence for the readings:** the pan-cancer TCGA passes (A210 H210.2, A213 H213.2, A215 H215.1) are not adjusted for
molecular subtype, which in breast carries most of the effect. They are read as "lower coherence in tumours of the
subtypes where these lesions are common", not as an effect of the lesion, until a subtype-stratified rerun across TCGA
projects (PanCanAtlas subtypes) is done. The METABRIC-TCGA discordance is a design difference (subtype strata), not
cancer type or platform.
<!-- src: paper/exploratory/A215_checkpoint_coherence/post_hoc_brca_pam50/{brca_pam50,brca_pam50_composition}.csv -->

## A217 / A218 (scored 2026-11-06): at matched molecular subtype the TCGA genomic-state losses are not detected (H217.1-3 FAIL); the proliferation loss is not clock-specific in the subtyped projects (A218 not scored)
Power note: `registered_experiments.md` (commit 17442cc). PanCanAtlas subtypes from the cBioPortal API (sha256 in
`$SCCC_DATA/tcga/pancan_subtypes/api/PROVENANCE.md`). 11 projects carry only the project name as subtype; with the
>= 2-level rule (choice before the tests), 10 projects enter (CESC, COAD, ESCA, GBM, HNSC, LGG, READ, SARC, STAD, TGCT,
UCEC; BRCA separate).
| test | cells / tumours | excess at matched subtype | p (Holm) | MDE | pseudo pct | original design, same tumours (excess) |
|---|---|---|---|---|---|---|
| H217.1 TP53 | 17 / 1,189 (7 projects) | +0.035 | 0.76 (1.0) | 0.128 | 85.5 | −0.068 |
| H217.2 FGA | 39 / 1,579 (9) | −0.015 | 0.32 (0.96) | 0.082 | 38.0 | −0.044 |
| H217.3 MYC | 9 / 387 (4) | +0.033 | 0.70 (1.0) | 0.166 | 34.0 | −0.051 |
- **All three FAIL.** At matched subtype the point estimates move to 0 or reverse (TP53 −0.068 → +0.035, FGA −0.044 →
  −0.015, MYC −0.051 → +0.033). The original-design losses on these tumours are themselves below the matched-design MDEs,
  so the fails bound the lesion effects at about 0.08-0.17 rather than excluding the original sizes. BRCA (descriptive):
  TP53 −0.188 → −0.007, FGA −0.286 → −0.043, MYC −0.263 → −0.084 (original vs matched excess).
- **A218 not scored (precondition):** P_proj (proliferation, within project x purity) −0.043 is not below its pseudo 5th
  percentile (−0.071) in these projects. Descriptive: P_sub −0.061 (pseudo 5.0th percentile), A 1.41 (jackknife range
  1.16-2.99); per project the larger losses are STAD −0.274, COAD −0.182, HNSC −0.156, CESC −0.154. BRCA: P_sub +0.049 vs
  P_proj −0.078. The proliferation loss is not smaller within subtype in these projects, but neither value is
  clock-specific at the registered level.
- **Reading (registered):** the pan-cancer TCGA genomic-state passes (A210 H210.2, A213 H213.2, A215 H215.1) are read as
  subtype composition, not lesion effects, in the projects with subtypes and in breast; the 11 projects without subtypes
  cannot be checked this way. With METABRIC (PAM50 strata, no loss for MYC and FGA, TP53 below MDE), no genomic lesion
  tested lowers tumour clock coherence beyond molecular subtype and proliferation. Whether subtype carries the proliferation
  loss is not settled: in breast the loss holds within PAM50 in SCAN-B (A127 H127.1) and not in METABRIC (A166 H166.3);
  in the other subtyped projects it is not clock-specific either way.
<!-- src: paper/exploratory/A217_subtype_matched/output/{power,tests,attenuation,attenuation_jackknife,brca_descriptive,descriptive_per_project,arm_composition,projects}.csv -->

## A219 (scored 2026-11-06): luminal A breast tumours keep more clock signal than other subtypes at matched proliferation (H219.2-3 PASS, clock-specific); outside breast, subtypes do not differ (H219.1 FAIL, powered)
Power note: `registered_experiments.md` (commit 726d173). Cells = cohort x proliferation tertile; subtypes with >= 15
tumours per cell.
- **H219.1 FAIL (primary, TCGA non-BRCA, 9 projects, 23 cells, 2,388 tumours):** S 0.073 vs null 0.076 (SD 0.010); excess
  −0.004, p 0.62, Holm 0.62, MDE 0.025; pseudo 14th percentile. At matched proliferation, subtypes of CESC, COAD, ESCA,
  HNSC, LGG, SARC, STAD, TGCT and UCEC do not differ in clock coherence beyond sampling noise. Unmatched (descriptive):
  excess +0.014, p 0.13.
- **H219.2 PASS, clock-specific (METABRIC + TCGA BRCA, 6 cells, 2,526 tumours):** LumA − (LumB + Her2 + Basal) +0.276
  (null −0.001, SD 0.062), p 0.0005, Holm 0.0015, MDE 0.155; above all 200 pseudo sets (q95 0.141). METABRIC +0.277, TCGA
  BRCA +0.275. Weighted coherence LumA vs others: METABRIC 0.647 vs 0.369, TCGA BRCA 0.489 vs 0.295.
- **H219.3 PASS, clock-specific, narrowly (SCAN-B, 3 cells, 2,865 tumours; unmatched ordering seen in A127):** +0.222
  (excess +0.207), p 0.0005, Holm 0.0015, MDE 0.121; pseudo 97.5th percentile (q95 0.189). LumA 0.706 vs others 0.332.
- **Descriptive:** between-subtype proliferation imbalance within cells 0.01-0.23 SD. Unmatched breast contrasts +0.218
  (METABRIC + BRCA) and +0.405 (SCAN-B): matching on proliferation removes half of the SCAN-B gap and little of the
  METABRIC + BRCA gap. Within Her2-enriched tumours, ER+ vs ER− coherence 0.482 vs 0.275 (METABRIC) and 0.498 vs 0.192
  (SCAN-B).
- **Reading (registered):** in breast, luminal A identity goes with more clock signal at equal proliferation, by about 0.2-0.28
  in coherence, replicated in three cohorts, clock-specifically. With the ER− losses within PAM50 (A127 post hoc, A166
  H166.2) and the Her2 ER split above, the luminal / ER state is the breast tumour-state correlate of clock coherence,
  beyond proliferation and beyond any genomic lesion tested (A217). Outside breast no subtype effect is detected (bounded at
  0.025 on S): the lineage effect is not general across cancer types. Not shown: that ER signalling acts on the clock
  (observational; ER and PAM50 are collinear).
<!-- src: paper/exploratory/A219_subtype_coherence/output/{power,tests,descriptive_per_group,descriptive_per_cohort,descriptive_unmatched,descriptive_cells_stat_and_prolif_imbalance,descriptive_er_within_subtype}.csv -->

## A222 (scored 2026-11-06): clock genes are altered at background rates; PER1 and PER3 losses ride on 17p and 1p36 deletions (H222.1-3 FAIL, powered)
Power note: `registered_experiments.md`. MC3 PASS-class coding calls, 9,730 samples (448 hypermutators > 1,000 calls
excluded); GISTIC2 by gene, 26 projects with >= 100 samples (9,681 samples); 1,000 random 14-gene sets.
- **H222.1 FAIL (mutation selection):** log(non-silent / silent) over the 14 clock genes 1.152 vs null 1.136 (SD 0.145),
  59th percentile, p 0.41, Holm 1.0, MDE 0.360 (a 1.43-fold excess). Mutation rate per Mb of CDS is near the length- and
  expression-matched median for every clock gene (DBP and TEF lowest, CIART highest).
- **H222.2 FAIL (deep-deletion focality):** clock-gene deep-deletion frequency minus its 50-gene flank 0.00019 vs null
  0.00027 (SD 0.00043), 48th percentile, p 0.52, MDE 0.00107.
- **H222.3 FAIL (amplification focality, two-sided):** 0.00043 vs null 0.00034 (SD 0.00039), 71st percentile, p 0.58.
- **Descriptive (registered):** PER1 is lost (GISTIC <= -1) in 37.4 % of tumours and TP53 in 38.2 %; 98.9 % of PER1 losses
  carry a TP53 loss, and PER1 is lost in 0.66 % of TP53-intact tumours (deep: 91.3 % of PER1 deep deletions carry a TP53
  deep deletion). PER3 is lost in 26.8 %; 97.7 % of these carry a loss of CHD5 (1p36 tumour-suppressor candidate, same
  block; choice of reference gene made here) and 93.6 % of TP73. Truncating share of non-silent calls: clock genes 2-32 %
  vs background median 12.9 %; CLOCK 16 / 50 (98.4th percentile of 10,077 genes) and ARNTL 9 / 39 (92.5th) highest, CRY2
  1 / 44 (1.1th) lowest. With 14 genes, one at the 98th percentile is expected about one time in five; not read as
  selection.
- **Reading (registered):** clock genes are mutated and lost as expected for their length, expression and position.
  "Frequent clock-gene alteration" reports background rates, and PER1 / PER3 "loss" is 17p (TP53) and 1p36 arm-level
  deletion. Clock genes are passengers in TCGA tumours.
<!-- src: paper/exploratory/A222_clock_alterations/output/{power,tests,null_draws,descriptive_genes,descriptive_17p_1p36,descriptive_truncating}.csv -->

## A221 (scored 2026-11-06): clock genes are lower in tumours than in paired normal tissue, but not uniformly; GTEx designs exaggerate the loss and invert it for ARNTL, NPAS2 and NR1D1 (H221.1-3 PASS)
Power note: `registered_experiments.md` (commit d151634). 13 TCGA projects, 621 tumour-normal pairs, GTEx recount3 (same
pipeline). 200 expression-matched random 14-gene sets per project (p floor 0.005).
- **H221.1 PASS (paired, same surgery):** L −0.330 log2 (random sets +0.091 on average, i.e. clock mean −0.238), p 0.005,
  Holm 0.010, MDE 0.163. Leave-one-gene-out (post hoc, approximate: random baseline kept): −0.265 (without PER1) to −0.384.
  **Not a uniform loss.** Mean paired LFC over projects: PER1 −1.07, PER3 −0.74, TEF −0.74, CRY2 −0.67, RORA −0.67, PER2
  −0.38; CRY1 −0.22, NR1D2 −0.20, CLOCK −0.01, ARNTL +0.02; DBP +0.18, NPAS2 +0.30, NR1D1 +0.38, CIART +0.47.
- **H221.2 PASS (design):** Delta −0.229 log2 (two-sided p 0.005, Holm 0.010, MDE 0.121): against GTEx, clock genes look
  more downregulated than against the paired normal. The GTEx design moves ARNTL (+0.02 → −1.01), NPAS2 (+0.30 → −0.74),
  NR1D1 (+0.38 → −0.93) and PER1 (−1.07 → −2.07) down, and PER2, PER3, TEF, NR1D2 up. Share of clock genes "significantly
  down" per project 7-71 % paired vs 21-86 % against GTEx (random genes 2-56 % vs 23-53 %; `share_down.csv`).
- **H221.3 PASS (mechanism):** across 182 gene x project points, the design difference follows the GTEx daytime offset
  (time of death 08-18 h minus all donors), Spearman +0.203, p 0.0005, Holm 0.0015: part of the GTEx distortion is
  sampling time.
- **Descriptive:** immediate-early genes are far lower in tumours than paired normal (10 genes; −2.38 to +0.46 log2 per
  project, 11 of 13 projects below −0.2), so the PER1 / TEF / RORA part of the paired loss may be stress / glucocorticoid
  response of the adjacent tissue rather than clock. Clock level vs ABSOLUTE purity within tumours: Spearman −0.29 to
  +0.29, no consistent sign (no simple composition gradient). GTEx by Hardy class: `gtex_by_hardy.csv`.
- **Reading (registered):** at the same surgery time, clock genes are lower in tumours than matched genes (real level
  loss), but the loss sits on PER1-3, CRY2, TEF and RORA, while ARNTL, NPAS2, NR1D1, DBP and CIART are flat or higher: not
  "clock genes are downregulated". GTEx comparators add a clock-specific distortion (−0.23 log2), partly sampling time,
  which turns ARNTL / NPAS2 / NR1D1 into apparent losses: published BMAL1 / NPAS2 / REV-ERBα "downregulation" from
  tumour-vs-GTEx designs is not supported by paired samples. Open (not tested): whether the PER / CRY2 / TEF pattern is a
  damped clock read at daytime surgery (morning-peaking genes fall, evening-peaking genes rise) or adjacent-tissue stress.
<!-- src: paper/exploratory/A221_tumour_normal_levels/output/{power,results,per_gene_project,share_down,immediate_early,purity,gtex_by_hardy,eligibility}.csv -->

## A223 (2026-11-06): not estimable (power rule); no level-based damping readout at 14 genes
Rotation-null MDE of the mean slope 2.8-3.2 (CPTAC 8 types, TCGA 13 projects; limit 0.5, complete damping at lam = 1 is
1). Gene-specific tumour-normal level changes dwarf the phase-predicted pattern at 14 genes, so the A221 observation that
morning-peaking PER genes fall while ARNTL / NPAS2 rise stays descriptive: neither a damping readout nor evidence against
one. Not computed: B at the observed phases.
<!-- src: paper/exploratory/A223_damping_levels/output/power.csv -->

## A220 (scored 2026-11-06): circadian risk signatures do no better than random gene signatures (H220.1-2 FAIL, powered); the "validation" rate is set by the cohort, not the list
Power note: `registered_experiments.md` (commit 8df048e). Published-style pipeline (univariate Cox filter, L1-Cox with CV
penalty, 10 held-out 30 % splits) on GOBP_CIRCADIAN_RHYTHM (137-154 expressed genes) vs 200 size- and expression-matched
random lists per cohort; 18 TCGA projects with >= 60 OS events, METABRIC, SCAN-B.
- **H220.1 FAIL (unadjusted held-out C, TCGA):** mean percentile of the circadian list among random lists 0.518 (null 0.5,
  SD 0.071), p 0.35, Holm 0.70, MDE 0.177.
- **H220.2 FAIL (incremental C over age, stage, proliferation, subtype):** 0.479, p 0.61, Holm 0.70, MDE 0.177.
- **H220.3 not scored (MDE 0.535):** percentiles METABRIC 0.775, SCAN-B 0.42 (unadjusted); 0.63, 0.675 (incremental).
- **Descriptive:** the circadian signature "validates" (median held-out median-split log-rank p < 0.05) in 8 of 20 cohorts
  (KIRC, LGG, LIHC, LUAD, MESO, UCEC, METABRIC, SCAN-B), exactly where random lists validate often (KIRC, LGG, METABRIC 100 %,
  SCAN-B 99.5 %, MESO 77 %, LIHC 33 %, LUAD 24 %, UCEC 18 %); in the other 12 random lists validate 0-17 % and the circadian
  list does not. The 14 core clock genes are worse than random 14-gene lists: percentile 0 in 6 cohorts, <= 0.10 in 12 of
  20, mean 0.227 (unadjusted C).
- **Reading (registered):** circadian risk signatures are the Venet effect: they validate where any gene list does
  (cohorts in which outcome tracks broad transcriptional programs such as proliferation or grade) and add nothing over
  random lists or over clinical covariates and proliferation. The core clock genes carry less outcome information than a
  random gene set of the same size (descriptive, post hoc reading of a registered descriptive).
<!-- src: paper/exploratory/A220_signature_benchmark/output/{power,tests,per_cohort,eligibility}.csv -->

## A224 (scored 2026-11-06): most of the luminal A coherence advantage is cis copy-number noise on the clock loci (H224.1 not clock-specific after adjustment); TP53 mutation lowers PER1 / PER2 beyond 17p dosage (H224.3 PASS, not adjusted for tumour state)
Power note: `registered_experiments.md` (commit af6bcfc). Dosage fitted with leave-the-gene-out phase covariates (user
caution: clock levels depend on sampling time); adjusted expression y - beta_g CN_g; 200 pseudo sets adjusted alike.
- **H224.1 (METABRIC + TCGA BRCA, 6 cells, 2,526 tumours): LumA − others +0.276 → +0.101 after adjustment**; permutation p
  0.011 (Holm 0.022) but pseudo 92nd percentile (q95 0.128): **not clock-specific, so not a registered pass** (registration:
  p < 0.05 and above q95; the script's verdict column reads "pass" on p alone). Attenuation 63 % (jackknife 95 % 26-101 %).
  MDE 0.133.
- **H224.2 FAIL (TCGA TP53, 46 cells, 6,449 tumours):** −0.085 → −0.052, p 0.097, Holm 0.097, pseudo 6.5th percentile;
  attenuation 38 % (−6 to 82 %). MDE 0.066.
- **H224.3 PASS (14 projects, 2,362 tumours without PER1 / TP53 loss and PER2 CN 0):** TP53 mutant − wild-type,
  phase- and dosage-adjusted mean log2 PER1 / PER2, −0.156 (two-sided p 0.0005, Holm 0.0015, MDE 0.079). Naive contrast over
  all tumours −0.233 (8,530). 2 x 2 means (mutation x 17p loss): WT / no loss +0.115 (3,807), WT / loss −0.076 (1,320),
  mutant / no loss −0.024 (1,416), mutant / loss −0.152 (1,987). Not adjusted for proliferation or subtype (as registered).
- **Dosage (before scoring, power note):** all 14 betas > 0 in TCGA (0.24-0.46 log2 per GISTIC step); CN carries 2-15 % of
  between-tumour variance.
- **Post hoc breakdown (not registered; `posthoc_*.csv`, `scripts/posthoc_a224.py`):**
  - Per cohort: METABRIC 0.277 → 0.117 (58 %), TCGA BRCA 0.275 → 0.085 (69 %).
  - Which genes: ARNTL alone 39 % (METABRIC 54 %, BRCA 24 %), PER3 alone 21 % (BRCA 41 %, METABRIC 2 %), PER2 11 %, NR1D1
    9 %, others <= 5 %; all but ARNTL 37 %. Not additive.
  - CN variance is 2-4 x higher in non-LumA tumours for most clock genes (e.g. METABRIC CLOCK 0.070 vs 0.288, ARNTL 0.160 vs
    0.363); PER1 is the exception (ratio 0.87-1.05). LumA is the copy-number-quiet subtype.
  - Cis vs trans: adding proliferation, subtype and FGA keeps 88-110 % of each TCGA beta (PER1 0.388 → 0.342); in METABRIC
    62-129 % (ARNTL 0.201 → 0.179). The slopes are mostly cis dosage, so over-adjustment by tumour-state biology is small in
    TCGA and partial in METABRIC (FGA is itself derived from CN and may absorb some cis signal).
  - H224.2 per project: the adjustment moves projects both ways (GBM −0.026 → +0.163, HNSC −0.043 → +0.077, BRCA −0.209 →
    −0.120; LUAD −0.132 → −0.202).
- **Reading (registered + post hoc qualification):** copy-number variation at the clock loci, larger in non-luminal-A
  tumours, lowers their between-tumour clock coherence; removing it leaves a LumA advantage of +0.10 that matched non-clock
  genes reproduce. The A219 H219.2-3 passes therefore mostly reflect genome stability of luminal A tumours, not a
  lineage-specific clock (qualified below). The TCGA TP53 coherence loss was compositional (A217) and partly positional
  (38 %). PER1 / PER2 are lower in TP53-mutant tumours without 17p loss, so "p53 and the clock" is not only positional, but the
  contrast is not adjusted for proliferation or subtype and is not a mechanism. Methods consequence: between-tumour coherence
  contrasts need a cis copy-number companion (clock and pseudo genes adjusted alike); the pseudo-clock companion alone did
  not catch this because random genes rarely share the clock genes' copy-number exposure in these subtypes.
<!-- src: paper/exploratory/A224_cis_copy_number/output/{power,tests,dosage,descriptive_naive,descriptive_2x2,pseudo_stat,posthoc_cells,posthoc_leave_one_gene,posthoc_cn_by_subtype,posthoc_beta_trans}.csv -->

## A219 qualification (2026-11-06, after A224): the luminal A advantage is mostly cis copy number
The A219 H219.2 pass (+0.276, clock-specific) stands as scored. After cis copy-number adjustment of clock and pseudo genes
alike (A224 H224.1) it is +0.101 and no longer above the adjusted pseudo sets (92nd percentile). The reading "luminal A / ER
state is the breast tumour-state correlate of clock coherence beyond proliferation" is withdrawn as stated; supportable:
"luminal A tumours are more clock-coherent, mostly because their clock loci vary less in copy number". H219.3 (SCAN-B) has
no copy number and is not re-examined; the ER-within-subtype descriptive (A219) is likewise not CN-adjusted.

## Correction to the A224 reading and the A219 qualification (2026-11-06, after A225 part A)
The A225 simulation (real joint CN profiles, A224 dosage and CN x phase terms, identical clocks in all subtypes) gives a raw
LumA − others D of −0.003 (−1 % of the observed +0.276), and the adjusted estimator equals the raw one there. Additive cis
copy-number noise therefore does **not** explain the luminal A gap, and the A224 statement "most of the luminal A coherence
advantage is cis copy-number noise" is withdrawn. A224's 63 % attenuation on real data comes from something its dosage
model does not describe (candidates: copy number correlated with clock-related expression within arms, i.e. biology or
tumour state, or a property of the pooled TCGA betas); diagnosis pending (A225 post hoc). Until then the A219 H219.2-3
result stands as scored (clock-specific, three cohorts), with the open caveat that cis-CN adjustment attenuates it for an
unexplained reason; the A219 qualification above is superseded.

## A225 (scored 2026-11-06): the copy-number question is not resolvable with bulk data; cis dosage noise is excluded as the cause of the luminal A gap, the A224 attenuation stays unexplained
Gate and power: `registered_experiments.md` (commit 5f97893).
- **Part A (registered):** with identical clocks and the real joint copy-number profiles, the raw LumA − others D is
  −0.003 (−1 % of the observed gap); raw and adjusted estimators pass the gate, masked fails (13 % rejections under the null;
  estimable only in METABRIC).
- **H225.1** = A224 H224.1 (adjusted estimator; +0.101, not clock-specific). **H225.2 (ARNTL trans test) not estimable**
  (no cell with enough copy-neutral tumours at the other 13 genes). **H225.3 FAIL** (PER1 −0.081, p 0.29; TEF +0.066, p 0.90;
  NR1D2 +0.054, p 0.62; MDE 0.18-0.26; pseudo sets estimable in only 5 / 21 / 41 of 200 after masking).
- **Post hoc diagnosis (`posthoc_*.csv`, `scripts/posthoc_a225.py`):** adjustment leaves LumA unchanged (METABRIC 0.646 →
  0.637, TCGA BRCA 0.488 → 0.493) and raises the other subtypes (0.369 → 0.520, 0.213 → 0.408). In the simulation the same
  adjustment raises them by +0.016 / +0.036 only. The simulated copy-number offsets match the real ones in size and in the
  share lying in the clock plane; the real data differ in copy number correlating with the clock-fitted part of some genes in
  non-LumA TCGA BRCA tumours (NR1D1 0.32 real vs 0.08 simulated; PER1 0.19 vs 0.06). The attenuation does not depend on how
  beta is estimated (cohort or pooled betas, with or without phase covariates: D 0.085-0.141). Pseudo-gene sets adjusted
  alike barely move (q95 0.141 → 0.128), so the adjustment acts on the clock genes specifically.
- **Reading:** the luminal A gap is not an additive copy-number artefact (part A). In non-luminal-A tumours, copy number at
  the clock loci is tied to their clock-gene co-variation beyond independent dosage, and removing it restores much of their
  coherence. Bulk data cannot say whether this is copy number acting on the oscillator (biology) or copy number marking a
  tumour state that also disorganises clock genes (confounding); the masked and trans designs that would separate them are
  not estimable at TCGA / METABRIC arm sizes. A219 stands as scored, with this caveat attached; neither "lineage sets the
  clock" nor "copy-number noise" is supported as the mechanism.
<!-- src: paper/exploratory/A225_cn_artefact_or_biology/output/{sim_calibration,sim_summary,power,tests,eligibility,posthoc_arm_coherence_real,posthoc_arm_coherence_sim,posthoc_offset_plane,posthoc_cn_vs_clock_residual,posthoc_variants_real}.csv -->

## A226 (2026-11-06): stromal host pairing cannot be confirmed in independent cohorts (not scored; too few new patients)
Fibroblasts 11 patients, endothelium 13, in the new cohorts that pass the A149 rules (MDE 0.38 / 0.40 vs threshold 0.12 /
0.13). Nearly all public same-patient tumour / normal single-cell studies are already A149 cohorts. The tumour-stroma ↔
normal-stroma clock-specific pairing stays a post hoc result (A149 note 5), reported as such; its registered-statistic
pairing (A151 H151.1, 0.860 vs null 0.504) is confirmed only as patient-specific, not as clock-specific.
<!-- src: paper/exploratory/A226_stroma_pairing_confirm/output/{eligibility,power}.csv; scout/SCOUT.md -->

## A227 (scored 2026-11-07): tumour myeloid and T/NK cells keep the patient's normal-tissue clock state, clock-specifically (H227.1-2 PASS); B / plasma and mast cells pair but not clock-specifically
Power note: `registered_experiments.md` (commit 0763cf8). 34 cohorts (A149 + A226), uncorrected |complex correlation|,
within-cohort permutation null (10,000), 200 pseudo-clock sets. First run killed for memory (8 workers, 39 GB); rerun with 2
workers, unchanged code.
- **H227.1 PASS (myeloid, 169 patients / 24 cohorts):** 0.581 vs null 0.361 (SD 0.036), p 1e-4, Holm 2e-4; pseudo median
  0.476, q95 0.548, p_clock 0.020: clock-specific.
- **H227.2 PASS (T/NK, 227 / 27):** 0.539 vs 0.320 (SD 0.031), p 1e-4, Holm 2e-4; pseudo median 0.461, q95 0.510, p_clock
  0.030: clock-specific.
- **H227.3 (B / plasma, 172 / 22):** pairing 0.420 vs 0.338, p 0.010, Holm 0.020; pseudo q95 0.465, p_clock 0.30: **not
  clock-specific**.
- **H227.5 (mast, 50 / 8):** 0.514 vs 0.390, p 0.036, Holm 0.036; pseudo q95 0.549, p_clock 0.17: **not clock-specific**.
  Mural not scored (37 patients).
- **Anchors (descriptive; mostly the A151 patients, not confirmatory):** fibroblasts 0.686 (pseudo q95 0.539), endothelium
  0.677 (0.512), both above every pseudo set.
- **Ordering (clock − pseudo median):** endothelium +0.255, fibroblasts +0.220, myeloid +0.106, mast +0.082, T/NK +0.077,
  B / plasma +0.020.
- **Per cancer type (descriptive, >= 12 patients):** lung strongest for all compartments (myeloid +0.181, T/NK +0.194,
  fibroblasts +0.329, endothelium +0.262); colorectal immune near 0 (myeloid +0.031, T/NK +0.003, B +0.033) while colorectal
  fibroblasts keep +0.218; liver and prostate positive (small n).
- **Reading (registered):** in the tumour, stromal and myeloid / T/NK compartments keep the same patient's normal-tissue clock
  state beyond matched non-clock genes; lymphoid B / plasma and mast cells share patient-specific structure that is not
  clock-specific. Immune passes carry the registered caveat: tumour and normal immune cells share the circulation (diurnal
  trafficking, composition), so they show host time reaching the tumour, not necessarily a local clock. The gradient
  stroma > myeloid ≈ T/NK > B follows residency. Colorectal immune cells pair without clock specificity (descriptive).
- **Note on the caveat (2026-11-07, wording only, verdicts unchanged):** a shared circulation does not argue against an
  internal clock: the shared signal is clock-gene, clock-shaped state beyond matched non-clock genes, so blood would act as
  a route of entrainment. Pairing does not separate three sources: the cells follow body time, a patient-specific immune
  clock state not tied to body time, or subset composition. A228 (registered 2026-11-07) tests the first (tumour immune vs
  normal-tissue fibroblasts) and the third (pairing within CD8, CD8-negative T and NK). Context: in healthy atlases
  lymphocyte coupling is not clock-specific (A150, 0.094 vs pseudo 0.087), so H227.2 is the project's first registered
  clock-specific T/NK result.
<!-- src: paper/exploratory/A227_tme_compartment_pairing/output/{power,tests,pooled,descriptive,per_cohort,eligibility,eligibility_by_cancer}.csv -->

## A228 (scored 2026-11-07): tumour immune states pair with normal fibroblasts but not clock-specifically (H228.1-2 fail); within T/NK, CD8-negative T cells keep the clock-specific pairing (H228.4 PASS), CD8 T and NK do not
Power note: `registered_experiments.md` (commit 9fbe4a8). 34 A227 cohorts; statistic, permutation null (10,000), 200
pseudo-clock sets and pooling as A227. First run killed for memory (2 workers, 30 GB cap); rerun with 3 workers under a 70 GB
cap (user, machine idle), unchanged code, peak 29.6 GB.
- **H228.1 FAIL (tumour T_NK vs normal fibroblast, 118 patients / 17 cohorts):** 0.441 vs null 0.366 (SD 0.043), p 0.044,
  Holm 0.044; pseudo median 0.394, q95 0.458, p_clock 0.15: pairing without clock specificity.
- **H228.2 FAIL (tumour myeloid vs normal fibroblast, 116 / 17):** 0.467 vs 0.370, p 0.015, Holm 0.029; pseudo q95 0.487,
  p_clock 0.11: not clock-specific.
- **Descriptive, vs normal endothelium:** T_NK 0.481 (q95 0.470, p_clock 0.035), myeloid 0.510 (q95 0.465, p_clock 0.005),
  both above pseudo q95 (not registered tests).
- **Positive control (tumour endothelium vs normal fibroblast, 76 / 13):** 0.642 vs pseudo q95 0.544, p_clock 0.005: two
  body-time compartments pair clock-specifically across cell types, so the cross-compartment design can detect it.
- **Same-sample normal immune vs normal fibroblast (ambient not excluded):** T_NK 0.490 (q95 0.504, p_clock 0.090); myeloid
  0.590 (q95 0.555, p_clock 0.025).
- **H228.4 PASS (CD8-negative T, 214 / 26):** 0.532 vs null 0.322 (SD 0.033), p 1e-4, Holm 3e-4; pseudo median 0.438, q95
  0.501, p_clock 0.015: clock-specific.
- **H228.3 (CD8 T, 191 / 23):** pairing 0.443 vs 0.321, p 4e-4, Holm 8e-4; pseudo q95 0.491, p_clock 0.31: **not
  clock-specific**.
- **H228.5 (NK, 67 / 9):** 0.492 vs 0.342, p 0.011, Holm 0.011; pseudo q95 0.543, p_clock 0.19: **not clock-specific**.
- **Per cancer type (descriptive, >= 12 patients):** colorectal at or below its pseudo median in every cross test (T_NK vs
  fibroblast −0.022, myeloid −0.009) and for CD8 T (−0.012); lung positive throughout (T_NK vs fibroblast +0.114, myeloid
  +0.244, CD8-negative T +0.220, positive control +0.313).
- **Reading (registered rules):** body time in the tumour immune compartments is not detected against normal fibroblasts
  (MDE 0.108; excess over null 0.076 for T_NK). By the registered rule this is "not detected", not evidence for a
  patient-specific immune trait; the endothelial reference, descriptive, is clock-specific for both, and lung carries the
  cross-compartment signal while colorectal sits at the pseudo level. Part 2 is the case the registration did not foresee
  (one of H228.3 / H228.4 passes): the clock-specific T/NK pairing survives inside CD8-negative T cells, so it is not
  CD4 / CD8 / NK mixing there; CD8 T cells pair without clock specificity. Composition inside the CD8-negative class (CD4
  naive / memory / Treg, gamma-delta, MAIT) remains possible.
<!-- src: paper/exploratory/A228_immune_body_time/output/{power,tests,pooled,per_cohort,pseudo,descriptive,subset_fraction_summary}.csv -->

## A229 (scored 2026-11-07): tumour endothelial VCAM1 and the tumour T/NK share follow the tumour endothelial clock state (secondaries H229.3, H229.5 PASS); the primaries, ICAM1 and CD8 share, do not
Power note: `registered_experiments.md` (commit 233e18f). 21 cohorts with tumour endothelium (of 34); statistic r = |sum y
conj z| / sqrt(sum y^2 sum |z|^2), within-cohort permutation of y (10,000), 200 pseudo-clock sets for z. 3 workers, 70 GB
cap (user), peak 32.8 GB, 228 s.
- **H229.1 FAIL (ICAM1, 190 patients / 21 cohorts):** 0.385 vs null 0.306 (SD 0.035), p 0.015, Holm 0.030; pseudo median
  0.335, q95 0.404, p_clock 0.11: associated, not clock-specific. ICAM1 ranks at the 87.5th percentile of 200
  expression-matched endothelial genes used as y.
- **H229.2 FAIL (CD8 share of tumour T cells, 172 / 20):** 0.312 vs null 0.314, p 0.51: no association.
- **H229.3 PASS (VCAM1, 190 / 21):** 0.426 vs null 0.308 (SD 0.034), p 8e-4, Holm 0.0024; pseudo q95 0.402, p_clock 0.020:
  clock-specific.
- **H229.4 FAIL (SELE):** 0.335 vs 0.306, p 0.20; p_clock 0.49.
- **H229.5 PASS, borderline (T/NK share of all tumour cells, 164 / 18):** 0.378 vs null 0.302 (SD 0.037), p 0.024, Holm
  0.048; pseudo q95 0.371, p_clock 0.040.
- **Implied peaks (descriptive, GTEx clock-time axis, `implied_peak_h`):** ICAM1 18.5 h, VCAM1 20.5 h, SELE 19.0 h (the
  three adhesion genes agree); T/NK share 4.8 h.
- **Adjacent-normal endothelium (descriptive):** VCAM1 0.438 (p 0.004; pseudo q95 0.447, p_clock 0.085), ICAM1 0.414
  (p 0.027; q95 0.453); SELE, CD8 share and T/NK share at null.
- **Per cancer type (descriptive):** lung carries the adhesion signal (clock minus pseudo median: ICAM1 +0.166, VCAM1
  +0.196, SELE +0.160; CD8 share +0.140); colorectal is at or below the pseudo level for ICAM1 (−0.077), SELE (−0.065)
  and CD8 share (−0.089).
- **Reading (registered rules):** the primaries fail: ICAM1 associates with the endothelial clock state but not beyond
  matched non-clock gene sets, and the CD8 share does not associate. Two secondaries pass: tumour endothelial VCAM1
  follows the endothelial clock state clock-specifically, and so, at the margin, does the T/NK share of the tumour. This
  is a partial human correlate of the Wang et al. (Cell 2024) gate: an adhesion molecule and T-cell abundance, not ICAM1
  or the CD8 share specifically. Association, not causation; lung-weighted as in A227 / A228. The implied hours are on
  the GTEx axis of human clock time and are not compared with the mouse ZT of Wang et al.
<!-- src: paper/exploratory/A229_endothelial_gate/output/{power,tests,pooled,per_cohort,pseudo,matched_genes}.csv -->

## A230 (scored 2026-11-07): no epithelial phase shift against the stroma in inflamed human gut (TAURUS H230.3 FAIL, MDE 1.7 h); GSE266616 not scored
Power note and design change: `registered_experiments.md` (commits c5ffb1c, 2a8fdc7). GSE266616 (blind): host-only gate
fails in both arms with 8 eligible patients: **H230.1, H230.2 not scored**. TAURUS (non-blind; stroma seen in A45 / A53):
25 patients, 82 non-inflamed and 60 inflamed samples, anchors fibroblast + endothelial, gate passes in both arms
(p_pairing 0.001).
- **H230.3 FAIL:** epithelial offset change (inflamed − non-inflamed, against the stromal host) +0.58 h (bootstrap SE
  0.83 h), within-patient swap p 0.39 (10,000); offset MDE 1.71 h. Pseudo-clock check (200 sets): q95 of |offset change|
  9.22 h, p_clock 0.77 (the pseudo offsets are wide because pseudo sets carry little phase signal; the check is weak for
  offsets and moot here).
- **Coupling (descriptive):** epithelial c 0.293 → 0.391, d +0.098 [−0.059, 0.302], d* +0.050, p_star 0.66; pseudo q95
  of |d*| 0.301.
- **Reading (registered):** a mouse-sized shift (about 9-12 h, regenerating colon) is not present in inflamed human
  intestinal epithelium; shifts below about 1.7 h are not excluded. Limits: TAURUS is non-blind, inflamed vs
  non-inflamed is confounded with treatment timing, chronic inflammation differs from regeneration 17 days after DSS,
  and the blind cohort could not be scored. The local clause of the headline is not supported by a human gut
  counterpart; it rests on the mouse colon, psoriasis (A115 descriptive) and tumours.
<!-- src: paper/exploratory/A230_ibd_local_offset/output/{tests,contrast_taurus,tests_gate_taurus,pseudo,power_taurus,power_gse266616,gate_gse266616}.csv -->

## A232 (scored 2026-11-07): tissue offsets track proliferation in fast deaths, but not clock-specifically, and not in ventilated donors or within tissues (H232.1-3 FAIL)
Power note: `registered_experiments.md` (commit 01b7f0c). Gate passes in both cohorts (fast 38 tissues, 244 donors;
vent 33 compartments, 381 donors; p_pairing 0.001). G2M = Hallmark v2024.1 (199 genes in GTEx, clock genes removed).
| | n tissues | estimate | p (Holm) | pseudo-clock q95 | random-set q95 | specific |
|---|---|---|---|---|---|---|
| H232.1 fast, Spearman offset vs G2M | 38 | rho 0.425 | 0.0089 (0.027) | 0.499 | 0.355 | no (clock check fails) |
| H232.2 vent | 30 | rho -0.068 | 0.72 (0.88) | 0.415 | 0.105 | no |
| H232.3 within tissue, top − bottom G2M tertile | 68 units | +0.044 h | 0.44 (0.88) | 1.27 h | - | no |
- **H232.1:** the rank association passes Holm and beats random expression-matched gene sets (G2M 0.425 vs random q95
  0.355; random-set median 0.292), but offsets from expression-matched pseudo-clock genes correlate with G2M about as
  strongly: |pseudo rho| ≥ 0.425 in 10 % of 200 sets (q95 0.499). Weighted slope 0.70 h per unit score. By the
  registered rule this is not a clock-specific association. The random-set median of 0.29 shows that fast-death tissue
  offsets align with a general axis of tissue expression.
- **H232.2:** no association in ventilated donors (between-tissue offset SD 0.61 h, reliability 0.70; MDE 0.52).
- **H232.3:** within tissues, donors with more proliferation are not phase-shifted (+0.044 h, MDE 0.137 h; fast +0.056,
  vent +0.027). Same sign as H232.1. The pseudo-clock |q95| (1.27 h) is wide, so this companion is weak here.
- **Reading (registered rules):** the Droin et al. (2019) prediction that proliferating organs run at a shifted clock
  phase is not supported in human GTEx: a cross-tissue association exists in fast deaths but cannot be separated from
  a tissue-expression artefact, it is absent in ventilated donors, and within tissues proliferation does not shift the
  offset by more than about 0.14 h. Tissues are the unit (MDE rho 0.46-0.52), so a modest cross-tissue link is not
  excluded.
<!-- src: paper/exploratory/A232_offset_proliferation/output/{score_tests,score_gate,score_pseudo_clock,score_random_sets,score_tissues,power}.csv -->

## A231 (scored 2026-11-08): inter-organ coupling and organ phase relations do not change with age in GTEx (H231.1 PASS equivalence; H231.2-3 FAIL, not detected)
Power note: `registered_experiments.md` (commit 4d6708a). Old 60-70 y vs young 21-49 y, 1:1 matched within cohort;
54 tissue x cohort units (vent 25, fast 29), 1,861 matched pairs. Holm over H231.1-3.
| | estimate (old − young) | 90 % interval | MDE_80 | p (Holm) | verdict |
|---|---|---|---|---|---|
| H231.1 pooled d* | +0.003 | [-0.027, 0.033] | 0.037 | 0.84 (1.0) | PASS: equivalent within ±0.15 |
| vent / fast d* | -0.015 / +0.021 | [-0.058, 0.028] / [-0.016, 0.058] | 0.051 / 0.052 | 0.47 / 0.32 | equivalent in both |
| H231.2 pooled offset | -0.045 h | [-0.27, 0.18] h | 0.20 h | 0.61 (1.0) | FAIL; equivalent within ±1 h |
| vent / fast offset | +0.032 / -0.113 h | [-0.33, 0.40] / [-0.37, 0.14] | 0.28 / 0.27 h | 0.78 / 0.30 | |
| H231.3 Q (reordering) | 54.5 (null mean 54.0, q95 73.5) | - | shift SD 0.45 h | 0.44 (1.0) | FAIL |
- **Companions:** d* in the 3 lowest-coupling tissues per cohort +0.016 (p 0.71). Pseudo-clock |offset| q95 1.54 h, so
  the offset companion could not have certified a shift below that size; no test reached it.
- **D231.4 timed comparator (fast, same matched donors, 29 tissues):** cosinor-on-time-of-death pooled offset +0.39 h
  (bootstrap SE 0.36) vs clocksync -0.11 h (SE 0.15); timed SE / clocksync SE 2.32. Per-tissue agreement circular
  correlation 0.64 (Pearson 0.65). The two estimates differ by 0.50 h, within their joint noise; clocksync is 2.3 times
  more precise on the same donors although they are timed.
- **D231.5 working hours: not estimable.** Among fast donors who died 08-18 h, no tissue has >= 20 donors per age arm
  and >= 10 matched pairs (young fast donors: 48 in total; eligibility counts in `score_d231_5_eligibility.csv`). The
  registered narrow-window comparison of the two approaches cannot be made on the age contrast. Code fix after the
  first scoring run crashed on the empty set: the guard records the status and the counts; no estimate was seen.
- **D231.6 host spread (descriptive):** fast: median spread old 3.92 vs young 4.09 (median log ratio -0.03; old lower in
  62 % of 29 tissues); vent: old 0.42 vs young 0.34 (log ratio +0.23; old lower in 8 % of 25). No consistent
  amplitude loss with age in the shared clock state; confounded by time of death (fast) and ICU state (vent).
- **Caveat:** T2D is more frequent in the old arm (median |SMD| 0.46 fast / 0.79 vent) and not adjusted; T2D leaves
  coupling unchanged within ±0.15 (A164 / A202).
- **Reading (registered rules):** the "ageing" clause of the headline becomes an equivalence for coupling: pooled d*
  old − young is within ±0.033 (90 %), in both cohorts, including the ventilated donors that no timed method can read.
  A mean organ shift relative to the body clock larger than 0.20 h, or organ-specific shifts with SD above 0.45 h, would
  have been detected 80 % of the time; none was seen. This does not contradict reports of rhythm dampening with age
  (Talamanca et al. 2023): amplitude loss shared by organs leaves coupling and offsets unchanged by construction.
<!-- src: paper/exploratory/A231_age_coupling/output/{score_tests,score_pooled,score_companion_d_star,score_pseudo,score_d231_4_timed,score_d231_5_eligibility,score_d231_5_workinghours,score_d231_6_spread,power,balance_summary}.csv -->

## A233 (scored 2026-11-08): hypertension, ischaemic heart disease, COPD, depression and smoking leave body-wide coupling unchanged (H233.c.1 PASS x5); no body-wide reordering or local organ shift detected (H233.c.2-3 FAIL x5)
Power note: `registered_experiments.md` (commit 4d6708a). Up to 3 matched controls per case within cohort; vent and
fast pooled. Holm over the five conditions per statistic.
| condition | body-wide d* [90 %] (MDE) | Q p (Holm) | module d* (MDE), p (Holm) | module offset h (MDE), p (Holm) | module − companion d*, p |
|---|---|---|---|---|---|
| hypertension | -0.001 [-0.021, 0.019] (0.028) | 0.23 (1.0) | +0.015 (0.054), 0.48 (1.0) | +0.17 (0.28), 0.15 (0.61) | +0.048, 0.26 |
| ischaemic heart disease | +0.009 [-0.017, 0.035] (0.035) | 0.97 (1.0) | -0.014 (0.087), 0.68 (1.0) | +0.20 (0.51), 0.34 (1.0) | -0.037, 0.53 |
| COPD | -0.017 [-0.043, 0.009] (0.035) | 0.32 (1.0) | +0.143 (0.22), 0.11 (0.56) | +0.04 (1.45), 0.96 (1.0) | +0.139, 0.15 |
| depression | +0.010 [-0.029, 0.049] (0.048) | 0.20 (1.0) | +0.104 (0.43), 0.59 (1.0) | +1.99 (3.07), 0.11 (0.56) | +0.140, 0.43 |
| smoking | +0.008 [-0.013, 0.030] (0.030) | 0.89 (1.0) | -0.074 (0.17), 0.28 (1.0) | -0.34 (1.21), 0.49 (1.0) | -0.093, 0.22 |
- **H233.c.1:** body-wide d* is equivalent to zero within ±0.15 for all five conditions; all 90 % intervals lie within
  ±0.05. **H233.c.2:** no condition reorders organ offsets (Q below the null q95 in all five). **H233.c.3:** no module
  coupling or offset change; the lung (COPD, smoking) and brain (depression, one fast unit) readouts are weak (module
  offset MDE 1.2-3.1 h) and only exclude large local effects. Pseudo-clock |offset| q95 3.2-10.0 h.
- **Descriptive per-cohort values (not tests):** hypertension module d* in vent +0.109 (p 0.004) but -0.032 in fast
  (pooled +0.015); the vent hypertension arm keeps an age imbalance (median |SMD| 0.54). Body-wide pooled offset (not a
  registered test): hypertension +0.12 h (p 0.040), depression +0.22 h (p 0.048); neither is corrected for the
  comparisons made and both are below 0.25 h.
- **Caveats:** comorbidity not removed (hypertension with T2D, SMD 0.66 / 0.72 vent / fast; smoking with COPD). BMI
  imbalance in COPD fast (0.49) and hypertension fast (0.35).
- **Reading (registered rules):** the "chronic disease" clause of the headline now rests on seven conditions (T2D and
  obesity from A164 / A202, plus these five): none changes inter-organ coupling by more than about 0.05 in d*, and none
  reorders organ phases detectably (organ-specific shift SD MDE 0.31-0.58 h). Local effects in lung and brain are not
  excluded below about 1 h.
<!-- src: paper/exploratory/A233_chronic_conditions/output/{score_tests,score_pooled,score_companion_module_minus_companion,score_pseudo_*,power,balance_summary}.csv -->

## A234 (scored 2026-11-08): internal time advances relative to clock time with age in GTEx fast deaths, in body and brain alike (H234.1 PASS, H234.3 PASS; H234.2 post hoc, no brain − body drift)
Power note: `registered_experiments.md` (commit 03239e5). Δ = internal (untimed host phase) − recorded time of death;
β > 0 = advance with age; strata Hardy x sex; 10,000 age permutations; 200 donor bootstraps with hosts refitted (all
succeeded). Holm over H234.1 and H234.3.
| | n | β (h per decade) | 90 % interval | MDE_80 | p (Holm) | verdict |
|---|---|---|---|---|---|---|
| H234.1 body − clock | 233 | +0.64 | [0.17, 1.11] | 0.70 | 0.010 (0.010) | PASS: advance |
| H234.3 brain − clock | 196 | +0.99 | [0.49, 1.49] | 0.85 | 0.001 (0.002) | PASS: advance |
| H234.2 brain − body (post hoc) | 192 | -0.01 | [-0.44, 0.41] | 0.66 | 0.96 | no drift detected |
- **Sensitivities:** Hardy 2 only: body +0.65 (p 0.050, 202 donors), brain +1.04 (p 0.005, 175); bracket age (archive
  comparability): body +0.55 (p 0.022); brain host without the cerebellar regions: brain − body -0.02 (p 0.92), brain −
  clock +1.04 (p 0.002).
- **Pseudo-clock companion (registered for H234.2 only):** uninformative as anticipated in the power note; pseudo hosts
  carry no time-of-day signal (Fisher-Lee with time of death ≤ 0.022), and their β sit at the estimator bound (|β| q95
  3.05; at the bound in 22-48 % of sets). The timed tests' specificity rests on D234.0 instead (both hosts follow time of
  death, Fisher-Lee 0.38 / 0.37).
- **Post hoc check (time of death vs age; chosen after scoring):** host phases are shrunk toward the mean, so an
  age-related time-of-death distribution could mimic a drift. Age is not associated with time of death (circular-linear
  R 0.05, permutation p 0.74 (body donors) / 0.80 (brain donors)), and adding cos / sin of time of death to a linear
  version of the model leaves β unchanged (body 0.570 → 0.576; brain 0.830 → 0.814). Script in the scratchpad; not part
  of the registered outputs.
- **H234.2 vs A231:** counting each donor once and using a brain-only host, the brain does not drift relative to the
  body (-0.01 h per decade, MDE 0.66). The A231 per-region pattern (brain regions old − young -1.15 h vs the A164 host)
  is not confirmed; it came from 10 regions on the same 20-25 matched pairs, read against a host that includes cortex
  and cerebellum.
- **Caveats:** the estimate is larger than the archive prior (0.35 h per decade) at low power (about 0.29 at the prior),
  so its size is probably inflated (winner's curse); the interval includes 0.35. Fast deaths only; time-of-death
  recording errors that change with age are not excluded; young fast donors are few (19 under 40 y).
- **Reading (registered rules):** with death type and sex held fixed and exact age, the shared internal time of GTEx
  donors runs ahead of clock time by about 0.5-1 h per decade of age, in the body and in the brain, and the brain keeps
  its phase relation to the body. Together with A231 (coupling and organ offsets unchanged with age): ageing shifts the
  timing of the whole body clock, not its organisation. Direction as the chronotype epidemiology (Roenneberg et al.
  2004; to be verified before citing).
<!-- src: paper/exploratory/A234_internal_time_age/output/{score_tests,score_sensitivity,score_pseudo,power,direction,counts}.csv -->

### 2026-11-08 sign note on A130 (pending PI confirmation, Q26)
A130 computes Δ = internal − wall time (`compute_living_time.py`, `IT.cdiff(a, b) = a − b`). Internal phase increases
with clock time (H113.1 adipose Fisher-Lee +0.82 vs ZT; V25.5 error 1.1 h after the offset), so Δ > 0 means the internal
clock is ahead (advance). The A130 registration read "slope < 0" as "earlier with age"; that label is inverted: its
-0.339 h per decade (one-sided p 0.054 for a negative slope) is a **delay** with age in the living cohorts, not
significant, opposite in sign to A234 (+0.64 h per decade, same convention, D234.0 checked). The phrase "direction as
the epidemiology and as A94" (A113 / A130 entry above) and the ledger row citing it are wrong as worded.
<!-- src: paper/exploratory/A113_living_time/output/{a130_summary,a113_summary}.csv; code paper/exploratory/A113_living_time/scripts/compute_living_time.py -->

### 2026-11-08 scout: no independent timed living cohort can replicate A234 (no hypothesis; nothing downloaded)
Floor as A130 (>= 30 subjects with age, range >= 25 y, known sampling clock time). GEO metadata checked directly for the
best lead, Wu et al. 2018 skin: GSE112660 (298 samples) and GSE139301 (533) list sex, tissue and collection time but **no
age**; the time-ordered arm is 20 men aged 21-49 y, and many population samples have collection time "NA". Others fail
on n, age range or recorded time: GSE326530 (20 subjects), GSE113883 TimeSignature (11), BodyTime (no age), sleep and
forced-desynchrony blood studies (GSE39445, GSE48113, GSE56931: young adults), GSE85358 skin (young vs old, no time),
del Olmo 2022 skin (population arm n 154, 20-74 y, no recorded times; timed arm 20 subjects, 21-49 y). Large blood
cohorts (Framingham, BIOS, EGCUT, GENOA) record age but no draw-time variable was found; access is controlled
(outward-facing, not pursued). A130 (54 subjects; null SD 0.20 h per decade, as precise as A234's 0.25 because living
subjects carry repeated timed samples) stays the only living comparison, and it disagrees in sign.

## A235 (scored 2026-11-09): no metabolic-organ phase shift in T2D or obesity; the multi-hour advance of high-fat-fed mice is excluded in fast deaths (H235.1 FAIL x2, equivalent within ±1 h)
Fast-death module (liver, both adiposes, skeletal muscle; A164 host), case − control offset, > 0 = advanced:

| condition | fast offset (h), 90 % | one-sided p / Holm | MDE | module − companion (h), p | vent offset (h), 90 % |
|---|---|---|---|---|---|
| T2D | +0.02 [−0.25, 0.30] | 0.46 / 0.46 | 0.37 | +0.46, 0.069 | +0.14 [−0.43, 0.70] |
| obesity | +0.20 [−0.23, 0.64] | 0.22 / 0.44 | 0.64 | +0.35, 0.20 | +0.07 [−0.70, 0.83] |

- **H235.1 FAIL for both** (no advance detected); both 90 % intervals lie within ±1 h, so a shift of an hour or more,
  and the multi-hour liver advance of high-fat-fed mice, is excluded in fast deaths. The pseudo-clock q95 is wide
  (2.7 h, 4.9 h): the specificity companion would not have confirmed a small shift.
- **Module − companion** is positive in both (+0.46 h, p 0.069; +0.35 h, p 0.20), not significant; the companion
  tissues (nerve, skins) drift slightly the other way in T2D (nerve −0.71 h).
- **Sensitivity host without the module anchors:** fast T2D +0.12 [−0.18, 0.42], obesity +0.19 [−0.28, 0.65]; unchanged,
  so absorption of a shared metabolic shift by the host does not hide one.
- **H235.2 (vent, descriptive):** no shift either; with no shift in fast deaths the negative control is uninformative
  about the feeding route.
- **Descriptive:** liver in fast deaths −0.57 h (T2D, 25 sets) and −0.53 h (obesity, 18 sets), the opposite sign to the
  mouse advance, single tissue, no test. Module coupling unchanged (pooled d* T2D +0.016 [−0.026, 0.058], obesity +0.001
  [−0.078, 0.081], agreeing with A164). Module log z2 ratio (relative amplitude, case vs control): T2D −0.003 (fast),
  −0.047 (vent); obesity −0.101, −0.127 (no null; a possible mild amplitude loss in obesity, untested).
- **Reading:** human metabolic disease does not reproduce the mouse high-fat-diet phase phenotype at the organ level,
  within ±0.3-0.6 h. Not a contradiction of the mouse work: the mouse shift needs rest-phase feeding, which obese people
  do not systematically have and which terminal illness removes; T2D and obesity are not separated (BMI imbalance
  0.68 in the T2D fast design).
<!-- src: paper/exploratory/A235_metabolic_offset/output/{score_tests,score_tissues,score_pooled,score_z2,score_pseudo_t2d,score_pseudo_obese}.csv -->

## A236 (scored 2026-11-10): the time-of-day state of the late mouse embryo is not shown to be a running clock, and H108.1 is not clock-specific (H236.1 FAIL, informative; H236.2 FAIL, uninformative; H236.3 FAIL)

Qiu 2024, E15.0-E18.75, 16 embryos in 4 batches, one embryo per batch at 0, 6, 12 and 18 h. Timed model: batch
intercepts, linear within-batch drift, cos and sin. Three deviations fixed in the power note before any observed
statistic (`registered_experiments.md`, A236): plasmode null for S1 (the acrophase null rejects in 0.41 of data sets
with no time-locked signal), simulated critical value for the rotation test, primaries on the 77 cell types present in
all embryos.

| test | statistic | value | critical value | p / Holm | power at the observed signal size |
|---|---|---|---|---|---|
| H236.1 gene content | S1, share of the time-locked harmonics in the clock plane | 0.069 | 0.657 (null mean 0.368) | 0.975 / 0.975 | 0.98 (clock-plane axis) |
| H236.2 rotation | jackknife t of the polarisation rho (rho 0.015) | 0.31 | 4.76 | 0.43 / 0.85 | 0.33 (rho 1), 0.14 (rho 0.5) |
| H236.3 specificity | R of H108.1 against pseudo-clock sets | 0.516 | 0.617 (responsive q95; matched 0.533) | 0.19 / 0.57 | - |

- **H236.1 FAIL, informative.** The 24-h harmonics of the 14 clock genes do not lie in the clock plane: S1 0.069, in
  the lowest 2.5 % of what the embryos' untimed variation gives by chance. Without the drift term S1 is 0.278 against
  a null mean of 0.283 (p 0.47); on all 100 cell types 0.062; leaving one gene out 0.052-0.271. A clock-plane axis of
  the size seen in the plane would have passed in 0.98 of the simulations.
  **H236.1b:** without PER1 and PER2, S1 0.059 (p 0.99; no drift 0.304, p 0.53). H108.1 without the two genes: R 0.435,
  p 0.053 (with them 0.516, p 0.014), so the alignment is not carried by the direct glucocorticoid targets; leaving one
  gene out, R ranges from 0.338 (without NPAS2) to 0.589 (without ARNTL).
- **Per gene (descriptive).** Drift model: NPAS2 (amplitude 0.45, peak 4.1 h) and CIART (0.22, 2.3 h) dominate; PER1
  peaks at 2.3 h and the other 11 genes between 14.1 and 17.0 h (amplitude <= 0.18). Peaks near 3 h and 15 h are the
  direction in which the drift term moves the harmonics (cos and sin shift by the same amount), so these harmonics
  are dominated by the contrast of 0 and 12 h against 6 and 18 h: NPAS2 is +0.28 at 6 h and -0.15, -0.09, -0.04 at 0,
  12 and 18 h. The design cannot tell a non-linear age effect from a 12-h component there. The antiphase of the two
  arms (11.4 h; reference 12.9 h) is NPAS2 alone: ARNTL and CLOCK peak with the E-box arm (15.7, 15.8 h).
  Without the drift term the E-box and D-box outputs (PER1, PER2, NR1D1, NR1D2, DBP, TEF) peak at 14.4-17.9 h and
  ARNTL, CLOCK, CRY1 and CIART at 23.1-0.8 h, with small amplitudes (ARNTL 0.024 against 0.504 in the reference, PER3
  0.003 against 0.273). The order resembles the adult one, but the share in the plane is what untimed variation gives
  (post hoc reading, no claim).
- **H236.2 FAIL, uninformative (power 0.33).** rho 0.015 (noise-corrected 0.037; leave-one-batch-out -0.031 to 0.040;
  all cell types 0.010); the time-locked share of the clock-plane harmonics is 0.42. At this signal size the plasmode
  gives rho -0.03 (sd 0.25) for an axis, 0.00 (sd 0.07) with no signal and 0.37 (sd 0.22) for a perfect rotation: the
  observed value sits with the first two, 1.6 sd below the third (post hoc). Without the drift term rho is 0.267
  (noise-corrected 0.98; t 1.80, p 0.19): the configuration the drift term was registered against (an age trend plus
  an axis imitates a rotation), not evidence of rotation. Named lineages, single cell types and no test: hepatocytes
  0.67 (se 0.91), cardiomyocytes 0.52 (se 0.31); trajectories with >= 5 cell types -0.14 (white blood cells) to 0.19.
- **H236.3 FAIL.** On genome-wide pseudobulks (24,204 genes with mean count >= 1), R 0.516 is the 93.0th percentile of
  200 expression-matched sets (q95 0.533) and the 81.5th of 200 responsive sets (q95 0.617). The timed test reaches
  p < 0.05 in 18.5 % of the matched and 41 % of the responsive sets (nominal 5 %). After residualisation on embryonic
  day: R 0.533, above the matched q95 (0.499; 97.5th percentile) but not the responsive one (0.612; 87.0th); 13.5 %
  and 27 % of the sets reach p < 0.05. Limit: expression matching leaves the same 14 responsive genes in every set (27
  available), so the 200 responsive sets are assignments of one gene set to the waveform, not independent draws.
- **Reading (as registered).** H236.1 fails: the wording stays "time-of-day state". H236.3 fails: the state is shared
  with cue-responsive genes, and to a lesser degree with arbitrary genes of the same expression. The reason is in the
  design: within a batch clock time and age move together (18 h in 6-h steps), so any gene with a developmental trend
  is time-locked, and permuting clock times does not remove it. **Consequence for A108:** H108.1 (R 0.52, p 0.014)
  stands as scored but is not clock-specific; it does not show a functioning clock and cannot serve to place embryo
  organs on external time. The fetal adrenal positive (A189) passed both companions and is not affected.
- **What this does not show.** It is not evidence against an embryonic clock: the rotation test has little power, the
  cardiomyocyte and hepatocyte estimates are positive, and entrained and driven oscillators cannot be separated
  cross-sectionally. Showing a running clock here needs a design in which age and clock time are crossed, or more
  than one cycle per litter age. PI question Q30.
<!-- src: paper/exploratory/A236_embryo_running_clock/output/{h236_primary,h236_1_content,h236_1_genes,h236_1_genes_nodrift,h236_1_profiles,h236_1_arms,h236_1b_timed,h236_2_rotation,h236_3_companions,h236_3_sets,power,power_info}.csv; build $SCCC_DATA/natexp/qiu2024/genome_late/PROVENANCE.md -->

## A237-A239 (scored 2026-11-11): no adult-like clock organisation detected before birth in human organs, gonad cell classes or the brain; the postnatal L1 controls all pass (every prenatal primary FAIL; H239.1 uninformative)

Three open datasets, one common design (age and sex adjustment, signed coherence C for L1, stratified gate for L2,
expression-matched and responsive pseudo-clock companions; `registered_experiments.md`, A237), registered together
before any clock statistic (commit `11aac5c`), one scoring run each. "Adult pass rate": the feasibility on adult
material cut to the design (alpha 0.05, no companions). Percentiles: clock statistic among 200 matched / 200
responsive pseudo sets.

**A237, Cardoso-Moreira 2019 (bulk RNA-seq, 7 organs, 4-20 wpc; postnatal samples of the same organs as control).**

| readout | n | statistic | p / Holm | matched / responsive pct | adult pass rate | verdict |
|---|---|---|---|---|---|---|
| control, postnatal forebrain | 21 | C 0.539 | 0.0143 | 99.0 / 99.5 | 0.880 | PASS |
| control, postnatal hindbrain | 25 | C 0.828 | 0.0001 | 100 / 100 | 1.000 | PASS |
| control, postnatal heart | 10 | C 0.841 | 0.0001 | 100 / 100 | 0.955 | PASS |
| control, postnatal kidney | 9 | C 0.708 | 0.0009 | 100 / 100 | 0.895 | PASS |
| control, postnatal liver | 12 | C 0.714 | 0.0014 | 100 / 100 | 0.975 | PASS |
| H237.1 L2, organs of one fetus, 9-20 wpc | 26 | clock share 0.071 (null mean 0.027); pairing z -0.60 | p_pairing 0.669, p_acrophase 0.008 | 7.5 / 2.5 | 0.825 | FAIL |
| H237.2 heart, 4-20 wpc | 34 | C -0.094 | 0.444 / 1.0 | 54.0 / 58.5 | 1.000 | FAIL |
| H237.2 liver | 37 | C 0.387 | 0.042 / 0.253 | 94.0 / 92.0 | 1.000 | FAIL (suggestive, not established; not clock-specific) |
| H237.2 forebrain | 32 | C 0.339 | 0.058 / 0.289 | 93.5 / 89.5 | 0.965 | FAIL |
| H237.2 hindbrain | 33 | C 0.259 | 0.129 / 0.516 | 78.5 / 80.0 | 1.000 | FAIL |
| H237.2 kidney | 27 | C -0.350 | 0.768 / 1.0 | 20.0 / 8.5 | 1.000 | FAIL |
| H237.2 ovary | 18 | C -0.194 | 0.571 / 1.0 | 31.5 / 34.5 | 0.940 | FAIL |
| H237.3 4-8 wpc: heart, liver, hindbrain, kidney | 14, 16, 16, 9 | C 0.081, -0.485, 0.411, 0.121 | 0.293, 0.889, 0.040, 0.296 / Holm 1.0, 1.0, 0.278, 1.0 | hindbrain 91.0 / 97.0 | 0.990, 0.985, 0.985, 0.875 | FAIL (hindbrain: suggestive, not established; not clock-specific) |
| H237.3 13-20 wpc: heart, liver, kidney | 9, 10, 10 | C -0.155, 0.203, 0.034 | 0.614, 0.155, 0.379 / Holm 1.0, 0.932, 1.0 | - | 0.870, 0.960, 0.870 | FAIL |

- **The controls pass on unrecorded times at n 9-25.** The five postnatal controls are clock-specific (all sets
  below the clock in four organs); leave-one-gene-out C 0.295-0.918. Postnatal testis (descriptive, adult pass rate
  0.465): C 0.682, p 0.0026, clock-specific.
- **H237.1 FAIL, informative within the registered limits** (adult pass rate 0.825, CI 0.765-0.875; no L2 control in
  this dataset). The clock-plane share of the shared covariance is above its acrophase null (p 0.008) but the
  individuals do not share it more than mis-paired ones (p_pairing 0.669), and other gene sets share more: 19.5 % of
  matched and 19.0 % of responsive sets reach p_pairing < 0.05. Other ranges (descriptive, adult pass rates
  0.140-0.630): p_pairing 0.393 (4-20), 0.757 (13-20), 0.423 (4-8 wpc); with the unlinked columns 0.345-0.641.
- **H237.2 / H237.3 FAIL.** No organ survives Holm. Two unadjusted p < 0.05 (liver 4-20 wpc, hindbrain 4-8 wpc), neither
  above the q95 of both pseudo sets: "suggestive, not established" by the registered rule. The other negatives
  (unadjusted p >= 0.05) are informative: adult pass rates 0.870-1.000 and the control of the organ passes (ovary has
  no control).
- **Descriptive.** 9-20 wpc windows: C -0.270 to 0.314, p >= 0.090. Testis 4-20 wpc: C 0.157, p 0.181. With the RNA
  quality number as a covariate: forebrain C 0.605, p 0.0044, above both companions (100th / 99.0th percentile; 26
  of the 32 libraries have the number; leave-one-gene-out 0.358-0.715); liver C 0.418, p 0.035 (96.5 / 93.0); the
  other organs p >= 0.103; the L2 gate at 9-20 wpc p_pairing 0.077, p_acrophase 0.0007 (pairing z 1.61; 78.5th /
  95.0th percentile). A covariate chosen as a sensitivity: no claim, a lead for the fetal forebrain at most.
- **Unsigned R against C in real data.** Postnatal forebrain: R 0.123, C 0.539. Prenatal heart 13-20 wpc: R 0.493, C
  -0.155; liver 9-20 wpc: R 0.510, C -0.148. The two statistics disagree in both directions, as in the simulation.

**A238, Garcia-Alonso 2022 (fetal gonad, scRNA, 50 donors, 6-21 wpf, four cell classes).**

| readout | n | statistic | p / Holm | matched / responsive pct | adult pass rate | verdict |
|---|---|---|---|---|---|---|
| H238.1 mesenchymal-stromal | 47 | C 0.037 | 0.272 / 0.815 | 65.5 / 64.5 | 1.000 | FAIL |
| H238.1 supporting | 46 | C -0.269 | 0.667 / 0.815 | 24.0 / 22.0 | 0.995 | FAIL |
| H238.1 endothelial | 24 | C 0.523 | 0.014 / 0.056 | 97.5 / 98.5 | 0.930 | FAIL (suggestive, not established; clock-specific) |
| H238.1 epithelial | 33 | C -0.038 | 0.376 / 0.815 | 52.0 / 56.5 | 0.985 | FAIL |
| H238.2 L2, four classes of one gonad | 47 | clock share 0.049 (null mean 0.042); pairing z 1.11 | p_pairing 0.137, p_acrophase 0.312 | 19.0 / 1.5 | 1.000 | FAIL |

- **H238.1 FAIL.** Endothelium is the closest: above both companions, leave-one-gene-out C 0.325-0.601, Holm 0.056.
  Per sex and age (descriptive, adult pass rates 0.650-0.745): male C 0.447 (n 11, p 0.031), female 0.263 (n 13, p
  0.159), from 13 wpf -0.016 (n 13, p 0.445): not carried by the later gonads. Stroma, supporting cells and
  epithelium are informative negatives (p >= 0.272; no positive control in this dataset).
- **H238.2 FAIL, informative.** The cell classes of one gonad share non-clock states strongly (68.5 % of matched and
  95.0 % of responsive sets reach p_pairing < 0.05) and the clock set less than most of them. Descriptive: three
  classes p_pairing 0.815, stromal-supporting pair 0.151; females p_pairing 0.030 with p_acrophase 0.222 (not a gate
  pass; 39.5th / 10.0th percentile).

**A239, Kang 2011 (exon arrays, 16 brain regions, 20 prenatal and 32 postnatal brains, 12 clock genes).**

| readout | n | statistic | p / Holm | matched / responsive pct | adult pass rate | verdict |
|---|---|---|---|---|---|---|
| control L1, postnatal brain mean | 32 | C 0.756 | 0.0003 | 99.5 / 100 | 1.00 | PASS |
| control L2, regions of one postnatal brain | 32 | clock share 0.314 (null mean 0.112); pairing z 95.10 | p_pairing 0.0001, p_acrophase 0.0036 | 95.0 / 62.0 (q95 95.03 / 122.74) | 0.98 | gate passes, not clock-specific: control not met |
| H239.1 L2, regions of one prenatal brain | 20 | clock share 0.059 (null mean 0.070); pairing z 14.08 | p_pairing 0.0001, p_acrophase 0.535 | 30.5 / 2.5 | 0.96 | FAIL, uninformative (control L2 not met) |
| H239.2 L1, prenatal brain mean | 20 | C -0.166 | 0.582 | 37.0 / 42.5 | 1.00 | FAIL, informative |

- **Control L1 PASS, H239.2 FAIL (informative).** Trajectory of the brain mean (descriptive): 13-22 pcw C 0.047 (n
  14, p 0.412); infancy and childhood C 0.735 (n 13, p 0.0029; 99.5th / 99.0th percentile); from 13 years C 0.663 (n
  19, p 0.0074; 100th / 100th). With RIN as a covariate: prenatal C -0.500 (p 0.881), postnatal 0.763 (p 0.0001).
  Unsigned R of the prenatal brains is 0.699 with C -0.166.
- **Control L2 not met, so H239.1 is uninformative as registered.** The postnatal gate passes both permutation tests
  but its pairing z is not above the responsive sets (62.0th percentile; 92.0th with RIN). Every pseudo set of every
  group reaches p_pairing < 0.05 (share 1.00): the regions of a brain share any gene set (brain identity, post-mortem
  interval, RNA quality), so the pairing statistic cannot single out the clock in this design. Post hoc reading, no
  claim: what differs between the groups is the clock shape of the shared covariance, the clock-plane share against
  its acrophase null: prenatal 0.059 against 0.070 (p 0.535), 13-22 pcw 0.057 against 0.078 (p 0.584), infancy and
  childhood 0.232 against 0.082 (p 0.0078), from 13 years 0.253 against 0.104 (p 0.0079), postnatal 0.314 against
  0.112 (p 0.0036).

**Reading (as registered for A237-A239 together).** None of the three conditions that would strengthen the
developmental claim is met: no shared state across the organs of a fetus (H237.1), no organ system passing in two
datasets (brain: forebrain, hindbrain and the Kang brain mean all fail; gonad: ovary and the four classes fail), no
L1 pass surviving Holm with both companions. With the L1 control passing in all six postnatal sets (five organs by
RNA-seq, the brain by arrays), the result is: no adult-like clock organisation is detected in human prenatal organs
at designs where adult and postnatal material pass; in the brain, the coherence across individuals is absent before
birth and present from infancy (descriptive trajectory). The fetal adrenal of A189 stays a single observation. Not
established but worth a replication: gonad endothelium (Holm 0.056, clock-specific), forebrain with the RNA quality
covariate (descriptive).
**Limits.** A negative needs clock phases that differ between individuals as they do in adults sampled in hospital
hours: fetal clocks entrained tightly to maternal time and collected in a narrow window, or clocks of low amplitude,
would be missed and cannot be told from no clock. The companions' cost is not calibrated; a Holm step is not in the
feasibility (worst-rank rates in the registration). Bulk composition changes with age (smooth part removed only).
PI question Q32.

**Post hoc: A189 with the signed coherence (not a registered test).** Fetal adrenal (n 7): C 0.926 (R 0.734), p
0.0001 (floor of 10,000 permutations), Holm over the nine compartments 0.0009, above both companions (q95 0.497
matched, 0.478 responsive; 100th percentile of both), leave-one-gene-out minimum 0.904. The eight other
compartments: C -0.310 to 0.321, p >= 0.129. The A189 verdicts are unchanged by the statistic.
<!-- src: paper/exploratory/A237_fetal_organs/output/{verdicts,l1,l2,l2_pseudo,design}.csv; paper/exploratory/A238_fetal_gonad/output/{verdicts,l1,l2,l2_pseudo,design}.csv; paper/exploratory/A239_fetal_brain/output/{verdicts,l1,l2,l2_pseudo,design}.csv; paper/exploratory/scouts_2026-11-10/feasibility_fetal/output/passrates_{cm_day_sex,gonad_sextrend,kang_all}.csv; paper/exploratory/A189_fetal_within/output/h189_1_signed_posthoc.csv -->

## A241 (scored 2026-11-11): whether the organs of a fetus share the clock state of its adrenal cannot be read in Cao 2020 (H241.1 FAIL, uninformative: adult pass rate 0.14)

Cao 2020, the 7 fetuses with an adrenal pseudobulk (89-122 days), whole organs with >= 1,000 cells; host-restricted
gate (organ x adrenal cross-covariances only, `min_shared` 3), 10,000 permutations; registered with its adult
calibration before any organ x adrenal statistic (commit `fb27d58`), one scoring run.

- **Power first: the design is not powered.** Adult donors (GTEx, deaths at 08-18 h) cut to the design and thinned
  to the fetal depths pass in 0.140 of 200 subsamples (95 % CI 0.095-0.196; `min_shared` 4: 0.125; 5: 0.060). The
  clock shape is seen (p_acrophase < 0.05 in 0.98 of the subsamples), the pairing of 7 individuals is not (p_pairing
  < 0.05 in 0.14).
- **H241.1 FAIL, uninformative.** `organs6` (lung 6, heart 4, liver 4, cerebellum 3, cerebrum 3, muscle 3 fetuses;
  6 organ x adrenal pairs): clock share 0.162 (null mean 0.042), p_acrophase 0.0071, p_pairing 0.706, pairing z
  -0.61 (31.5th percentile of the matched, 39.0th of the responsive pseudo sets; 7.5 % and 2.5 % of them reach
  p_pairing < 0.05). After removal of the expected adrenal RNA: p_pairing 0.685, p_acrophase 0.0068. The fetal
  result has the profile of the adult calibration (shape seen, pairing not): it does not separate "the organs follow
  the adrenal" from "they do not".
- **Adrenal RNA in the other organs is small.** Adrenal cortex marker fraction relative to the adrenal of the same
  fetus, in the six organs: median 0.001-0.016 per organ (maximum 0.077, heart); expected share of the clock-gene counts: median
  0.001-0.025 per organ, at most 0.101 for one gene (ARNTL, heart). Unlike the classes within the adrenal (A190),
  whole organs are not confounded by adrenal RNA.
- **Descriptive.** With the placenta (5 fetuses; 7 pairs): p_pairing 0.602, p_acrophase 0.0012, pairing z -0.32
  (36.0th / 34.5th percentile). One organ at a time: p_pairing 0.165 (cerebrum, 3 fetuses) to 1.0; lung with a
  linear age adjustment 0.256.
- **Reading (as registered):** uninformative, never "the organs do not follow the adrenal". The question needs more
  fetuses with an adrenal and other organs than any open dataset holds (Cao: 7).
- Code `paper/exploratory/A241_fetal_adrenal_host/scripts/{feasibility_a241.py (637 s), score_a241.py (1,664 s)}`.
<!-- src: paper/exploratory/A241_fetal_adrenal_host/output/{verdicts,l2,l2_pseudo,contamination,design_fetal,passrates,design}.csv -->

## A240 (scored 2026-11-11): the fetal negatives exclude a clock of 0.75 to 1 times the adult one, not a smaller one; a 4 h window of collection does not explain them

Adult material at the fetal designs of A237-A239 (same subsamples and seeds as the registered feasibility), with the
adult clock scaled by a (`clocksync.simulate.scale_clock`) or the GTEx pool restricted to deaths in a narrower window;
no fetal expression value is read. A calibration without pass / fail, registered before any scaled or windowed pass
rate (commit `fb27d58`), one run. 200 subsamples per readout and setting, 500 permutations, alpha 0.05.

- **Checks met.** At a = 1 every subsample reproduces the registered pass / fail (0 differing of 200 in the 14
  Cardoso-Moreira and the 5 gonad readouts, 0 of 100 in the Kang RNA-seq readout). At a = 0 the pass rate is 0.000 in
  all 21 readouts.
- **a80 (smallest a detected in >= 0.8 of the cohorts at it and above): 0.75 in 11 readouts, 1.0 in 10.** Of the 18
  negatives: 0.75 in 9 (heart 4-20 and 4-8 wpc, liver 4-8 wpc, hindbrain 4-20 wpc, kidney 4-20 wpc, the L2 of the
  organs of one fetus, gonad stroma, the L2 of the gonad classes, prenatal brain mean on RNA-seq; pass rates at 0.75:
  0.830-0.995) and 1.0 in 9 (heart, liver and kidney 13-20 wpc, kidney 4-8 wpc, forebrain and ovary 4-20 wpc, gonad
  supporting and epithelial classes, brain mean on the Kang arrays; pass rates at 0.75: 0.365-0.715). a50 is 0.75 in
  17 readouts, 0.5 in 2 (heart 4-20 wpc, L2 organs) and 1.0 in 2 (forebrain, arrays).
- **Below 0.75 the designs are blind.** Pass rate at a = 0.5: 0.000-0.745 (median 0.200; below 0.5 in 19 of 21
  readouts, the exceptions being heart 4-20 wpc, 0.745, and the L2 of the organs, 0.675); at a = 0.35: at most 0.110;
  at a = 0.25: at most 0.005.
- **On the platform of the prenatal brain negative (Kang arrays, 20 of the 32 postnatal brains):** 0.950 at a = 1,
  0.365 at 0.75, 0.000 at 0.5.
- **Window arm: restricting the adult deaths to 08-14 h or 08-12 h costs nothing.** 08-14 h (177 donors, 12
  readouts): pass rates 0.865-1.000, difference to the registered 08-18 h window -0.055 to +0.080 (median 0.000);
  08-12 h (113 donors, 10 readouts): 0.895-1.000, difference 0.000 to +0.085 (median +0.020). The median coherence
  changes by -0.073 to +0.054 and -0.061 to +0.049.
- **Descriptive.** Share of the clock-gene counts kept by the thinning, a = 0.75 to 0: 0.90-0.69 (organs), 0.81-0.48
  (gonad stand-ins), 0.92-0.73 (GTEx brain). a_excl (smallest a at which <= 0.05 of the adult subsamples fall at or
  below the fetal C): 0.25 (liver 4-8 wpc) to 1.0 (forebrain 4-20 wpc, hindbrain 4-8 wpc). The three suggestive
  readouts sit inside the adult distribution at a reduced or full clock: share of adult subsamples at or below the
  fetal C, liver 4-20 wpc 0.765 at a = 0.5; hindbrain 4-8 wpc 0.200 at 0.75; gonad endothelium 0.660 at 0.75 and
  0.140 at 1 (no a_excl). No bound is claimed for them.
- **Reading (as registered).** The informative negatives of A237-A239 exclude, in 8 of 10 such cohorts, a clock
  organised as in adults whose between-individual variation is at least 0.75 (9 readouts) or 1.0 (9 readouts) times
  the adult one; a clock of half the adult size would mostly have been missed and is not excluded. a bounds
  (amplitude) x (phase spread between individuals), but the adult rates hold in a 4 h window of deaths, so a
  collection restricted to a morning list does not explain the negatives; a lower amplitude, or phases more alike
  between fetuses than between adults at one clock time, stay open. No fetal verdict changes.
- **Limits (registered).** Conservative by construction, by an unknown margin (the scaling removes the non-clock
  variation of the clock plane with the clock; thinning removes counts): the pass rate at a is a lower bound of the
  detection rate, a80 an upper bound of the smallest detectable clock. Adult stand-ins and waveform; companions not
  calibrated; GTEx time of death is not a collection time; the array subsamples overlap (CI not valid).
- Code `paper/exploratory/A240_fetal_amplitude/scripts/amplitude_bound.py` (6,091 s, 8 workers).
<!-- src: paper/exploratory/A240_fetal_amplitude/output/{bound,curve,check_cm,check_gonad,check_kang}.csv; verdicts of A237-A239: A237_fetal_organs, A238_fetal_gonad, A239_fetal_brain/output/verdicts.csv -->

## A242 (scored 2026-11-11): the clock of the human fetal adrenal replicates in open data from other laboratories, from 12 weeks and before (H242.1 `late` PASS, `schroder` PASS; H242.2 `early` PASS), at a size below the adult one

Nine open cohorts, 59 fetal adrenals (one whole-adrenal pseudobulk or bulk library per fetus), none from Cao 2020;
L1 signed coherence C across fetuses, pooled over cohorts: cohort indicator and a common age polynomial regressed
out gene by gene, 10,000 per-gene permutations within cohort, both companions (200 sets each) among the genes
present in every cohort. Registered with its adult power calculation before any fetal clock statistic (commit
`5520253`). The scoring run stopped at the verdict step after writing `l1.csv` (a table without L2 tests was not
handled by `development.verdicts`); after the fix the same command was rerun with the same seed and gave an
identical `l1.csv` (same sha256): one set of results.

- **Power first: every registered readout is powered.** GTEx adrenal glands (93 donors, deaths at 08-18 h) standing
  for the fetuses, thinned to the clock-gene counts of each cohort, 200 subsamples, test at alpha 0.05: `late` 0.995
  (95 % CI 0.972-1.000; 0.930 at 0.75 of the adult clock, 0.375 at 0.5; without a clock 0.040), `schroder` 0.990
  (0.895; 0.410; 0.075), `early` 1.000 (1.000; 0.560; 0.065).
- **H242.1 `late` PASS** (fetuses from 12 weeks, the age range of A189; 19 fetuses of 7 cohorts: ENCODE 4, Kameneva
  3, Wang 3, Sasaki 3, Dong 2, Han 2, Roost 2): C 0.552, p 0.022 (Holm 0.022; null q95 0.462), clock-specific (99.0th
  percentile of the matched, 98.0th of the responsive pseudo sets; q95 0.409 and 0.420); leave-one-gene-out C
  0.356 (without NR1D1) to 0.643.
- **H242.1 `schroder` PASS** (13 bulk adrenals without recorded age, 3 sequencing batches; batch and one surrogate
  component of the non-clock genes regressed out): C 0.651, p 0.0024 (Holm 0.0048; null q95 0.406), clock-specific
  (100th percentile of both; q95 0.372 and 0.396); leave-one-gene-out 0.380 (without ARNTL) to 0.854.
- **H242.2 `early` PASS** (before 12 weeks, 5 to 11 weeks as stated; 25 fetuses of 5 cohorts: Kameneva 8, Wang 7,
  Sasaki 4, Kildisiute 4, Dong 2): C 0.732, p 0.0007 (null q95 0.400), clock-specific (100th percentile of both; q95
  0.355 and 0.363); leave-one-gene-out 0.608 (without NR1D1) to 0.878.
- **Descriptive (as registered).** No single cohort carries `late`: without one cohort at a time C 0.495 to 0.607,
  p 0.0099 to 0.046, clock-specific every time (97.0th percentile or above); the weakest is without ENCODE (C 0.495,
  p 0.046), whose fetuses come from the same tissue programme as Cao 2020. Without the two shallow cohorts
  (`late:deep`, 15 fetuses): C 0.634, p 0.0071. Every fetus with an age (`all`, 46 fetuses, 8 cohorts): C 0.750, p
  0.0006 (100th percentile of both). Schroder without and with two surrogate components: C 0.893 (p 0.0001) and 0.626
  (p 0.0056). Unsigned R: 0.584 (`late`), 0.658 (`schroder`), 0.070 (`early`).
- **The two cohorts large enough to be scored alone do not pass:** Kameneva (11 fetuses, 6-14 weeks) C 0.158, p 0.246
  (76.0th / 80.0th percentile); Wang (10 fetuses, 7-14 weeks) C 0.331, p 0.071 (93.5th / 90.5th). Adult pass rates
  at these designs are 0.995 and 0.985 for the full adult clock, 0.890 and 0.860 at 0.75, 0.360 and 0.330 at 0.5: a
  clock of adult size is excluded in each of them; one of half the adult size is not. The positive is a pooled one.
- **Size against adults** (share of the 200 adult subsamples with C at or below the fetal value; adult median C):
  `late` 3.0 % at the full adult clock (0.821), 16.5 % at 0.75 (0.703), 83.5 % at 0.5 (0.358); `schroder` 18.5 %,
  49.0 %, 94.0 %; `early` 6.5 %, 45.0 %, 100 %; Kameneva alone 0.0 %, 2.0 %, 16.5 %; Wang alone 0.5 %, 5.0 %, 47.5 %.
  The fetal coherence lies between that of half and of three quarters of the adult between-individual clock
  variation in the pooled readouts, lower in the two single cohorts. Size here is amplitude and spread of phases
  between fetuses together.
- **Reading (as registered).** `late` PASS: the clock-shaped covariation of the fetal adrenal (A189, Cao 2020, 7
  fetuses) is replicated across laboratories at the same age, on platforms and fetuses independent of Cao 2020: a
  two-dataset positive, to the PI (Q32). `schroder` PASS supports it on a third design. `early` PASS: the signal is
  present before 12 weeks. Not shown: that the adrenal starts before other organs. The negatives of A237-A239
  exclude a clock of 0.75 to 1 times the adult one (A240), and the adrenal signal read here is itself below 0.75 in
  `late`: a clock of this size in another organ would have been missed by those designs.
- **Limits.** Collection times are unrecorded and the cohort indicator removes what the fetuses of a cohort share;
  ages follow three conventions (up to 2 weeks apart); no sex adjustment; whole-adrenal pseudobulks (cell
  composition varies between fetuses; the matched and responsive companions are the control); `late` rests on p
  0.022 with 11 residual degrees of freedom; no leave-one-cohort-out was registered for `early`, of which Kameneva
  and Wang make 15 of 25 fetuses; adult stand-ins for power and size (35 % of the fetus x gene targets of `late`
  capped by the adult depth).
- Code `paper/exploratory/A242_fetal_adrenal_replication/scripts/{prep_adrenal.py (4 min 53 s), feasibility_a242.py
  (538 s, 8 workers), score_a242.py (42 s, 8 workers, peak 15.6 GB), position_a242.py}`; pooled cohorts with
  within-cohort permutations and surrogate covariates in `paper/common/development.py` (tests added).
<!-- src: paper/exploratory/A242_fetal_adrenal_replication/output/{verdicts,l1,design,design_cohorts,feasibility_passrates,feasibility_design,position_in_adult}.csv; paper/results/logs/jobs.tsv -->

## A245 (scored 2026-11-12): adhesion genes and T-cell content of healthy GTEx tissues do not follow the hours implied by A229 (H245.1-4 FAIL, powered)
Power note: `registered_experiments.md` (commit 511753c). 233 fast-death donors, 26 tissues (41-216 donors), 10,000
permutations of time within Hardy x sex; S = mean over tissues of the t at the predicted hour, R = phase-free companion.
- **H245.1 FAIL (VCAM1, 20.5 h, recorded time):** S -0.242 (null SD 0.405, MDE 1.01), one-sided p 0.73; matched-gene
  q95 0.691.
- **H245.2 FAIL (T-cell content, 4.8 h, recorded time):** S -0.680 (null SD 0.312, MDE 0.78), p 0.99: opposite in sign
  to the prediction; unchanged after residualising on HBB (-0.686).
- **H245.3 FAIL (recorded time):** ICAM1 (18.5 h) S 0.204, p 0.34; SELE (19.0 h) S 0.259, p 0.29.
- **H245.4 FAIL (internal time):** VCAM1 S -0.646 (p 0.95; pseudo-clock q95 1.419), T-cell content S -0.405 (p 0.90;
  pseudo q95 0.674). D = S(internal) - S(recorded): VCAM1 -0.404 [-0.918, 0.135], T cells +0.275 [-0.179, 0.792]
  (90 %, 2,000 donor bootstraps, hosts not refitted): neither axis is favoured.
- **Phase-free (descriptive, uncorrected):** against recorded time VCAM1 R 1.054 (null q95 0.977, p 0.032), T-cell
  content 0.699 (q95 0.764, p 0.077), PTPRC 0.698 (p 0.076), ICAM1 0.452 (p 0.66), SELE 0.386 (p 0.69); against internal
  time none below p 0.10 (ICAM1 1.027, p 0.11). Peak of the mean complex t over tissues (recorded time): VCAM1 13.6 h,
  T-cell content 16.4 h. Median cosinor amplitude over tissues: VCAM1 0.143, T-cell content 0.095 log2 units. Lung
  (176 donors): VCAM1 peak 14.2 h (amplitude 0.120), T-cell content 12.6 h (0.131). Offset sweep: recorded-time S of
  VCAM1 is largest at -6 h (0.97) and of T-cell content never above 0.02.
- **Ventilated donors (descriptive, internal time, 354 donors, 29 tissues; no pseudo-clock companion):** R 3.544 for
  VCAM1, 1.808 for T-cell content, 1.719 for ICAM1 (p 0.0005-0.0025); S at the predicted hours -3.241 (VCAM1), +1.207
  (T cells). In ventilated donors the host state also carries illness and treatment, and the adhesion genes are
  inflammation-inducible: not read as a clock effect.
- **Reading (registered rules):** the hours implied by the tumour endothelial clock state in A229 do not transfer to
  healthy tissue sampled over the whole day: at the predicted hours no readout is detected, with MDEs of 0.05 (T cells)
  to 0.11 (VCAM1) log2 units of amplitude in the median tissue. The registered "R passes while S fails" case holds for
  VCAM1 only on an uncorrected descriptive p of 0.032 (a midday peak, about 7 h from the prediction), which is not a
  claim. Not a refutation of the mouse gate: bulk tissue, post-mortem donors, and the T-cell score includes blood.
  A229's implied peaks stay descriptive and should not be quoted as hours.
- **Note (2026-11-12, after A248):** the blood positive control of A248 fails, so H245.2 is uninformative about T-cell
  biology (see A248); the registered verdict above is unchanged.
<!-- src: paper/exploratory/A245_trafficking_gtex/output/{power,tissues,direction,score_tests,score_descriptive,score_per_tissue,score_sweep,score_matched,score_pseudo}.csv -->

## A246 (scored 2026-11-12): tumour-microenvironment programs do not follow the clock state of the tumour at any phase (H246.1 FAIL; H246.2 associated but not clock-specific)
Power note: `registered_experiments.md` (commit 0f251fd). A227 cohorts; A229 statistic, 10,000 shared permutations, 200
pseudo-clock sets; nine programs scored on both states. Script by a Sonnet agent; numbers re-read in `tests.csv`.
- **H246.1 FAIL (state = tumour endothelial clock; 73-169 patients per program):** omnibus T -0.121 (null SD 0.340, MDE
  0.85), p 0.63; pseudo median 0.230, q95 0.825. No program below unadjusted p 0.05 (smallest: myeloid share 0.360 vs
  null 0.307, p 0.078).
- **H246.2 FAIL, not clock-specific (state = the program's own compartment; 162-318 patients):** omnibus T 1.268, p
  0.0003, Holm 0.0006, but pseudo-clock median 1.378 and q95 2.107 (p_clock 0.58): programs covary with the state of
  matched non-clock gene sets of their compartment as much as with its clock state.
- **H246.3 (programs):** only proliferation in malignant cells against their own state reaches Holm (r 0.378 vs null
  0.287, p 0.001, Holm 0.013), not clock-specific (pseudo q95 0.393, p_clock 0.12). All others fail.
- **Descriptive:** per cancer type the endothelial omnibus is between -0.379 (liver) and +0.043 (gastric), prostate
  +0.902 (85 program-patients, null SD 0.461); breast cohorts -0.130. Circular correlation of the implied peaks with the
  hours A146 transferred from the mouse: 0.302 (endothelial state), 0.855 (own state; nine programs, no test, and the
  associations it rests on are not clock-specific). Adjacent-normal tissue: endothelial-state omnibus 1.101 (p 0.0015;
  pseudo q95 1.072, p_clock 0.055, at the line), with T-cell activation (r 0.430, p_clock 0.020) and PD-L1 (0.464,
  p_clock 0.040) leading; own-state omnibus 1.149, not clock-specific (pseudo q95 1.873).
- **Reading (registered rules):** with A146 part 2, no clock-state organisation of these programs is detected in
  tumours, now without the mouse-to-human phase transfer and without a per-patient time estimate. The recorded-time
  result of A146 part 1 stands alone, as a hospital-clock effect that needs a second timed cohort. The own-compartment
  association is a generic covariation of expression states. The adjacent-normal pattern repeats the descriptive
  finding of A146 part 1 (programs follow body time in normal tissue, not in the tumour); it is descriptive and at the
  pseudo-clock line (D7: 1,000 sets before any wording).
<!-- src: paper/exploratory/A246_tme_clock_state/output/{power,eligibility,tests,pooled,per_cohort,pseudo,matched_sets}.csv -->

## A248 (scored 2026-11-12): GTEx blood does not show the published night peak of T cells (positive control H248.1 FAIL); glucocorticoid score and tissue T-cell hour not detected
Power note: `registered_experiments.md` (commit 79aa1f4). 208 donors with blood, 26 tissues, 10,000 permutations.
- **H248.1 FAIL (blood T-cell content at 02:10 h):** t -1.347 (null SD 1.007, MDE 2.51 = 0.274 log2, the published
  amplitude), one-sided p 0.91; fitted peak 15.3 h (phase-free p 0.34); -1.467 after residualising on the neutrophil
  score. Descriptive: naive score fitted peak 2.9 h (amplitude 0.048 log2, p 0.30 at 02:25 h); effector score fitted
  peak 13.3 h (0.199 log2; p 0.080 at 15:34 h).
- **H248.2 FAIL (GC score of 26 tissues at 10:40 h):** S 0.821 (null SD 0.641, MDE 1.60), p 0.10, Holm 0.20; above the
  matched q95 (0.675); fitted peak 7.9 h. Internal time: S 0.112.
- **H248.3 FAIL (blood, T-cell content vs GC):** partial r -0.011, p 0.45 (MDE 0.175).
- **H248.4 FAIL (post hoc to A245; tissue T-cell content at 14:00 h):** S 0.526 (null SD 0.316), p 0.045, Holm 0.090;
  below the matched q95 (0.632); fitted peak 16.4 h.
- **Descriptive:** tissue T-cell content against GC in the same sample: partial r -0.048 (p 0.024), -0.046 given the
  hour (p 0.036). Blood GC: fitted peak 11.9 h, amplitude 0.037 log2 (p 0.31). Ventilated donors (354, recorded time):
  GC of 29 tissues S 1.895 at 10:40 h (p 0.0015, fitted peak 10.7 h); blood and tissue T-cell content null (p 0.60,
  0.80). In ventilated care morning drug rounds (steroids) are an alternative to endogenous cortisol.
- **Reading (registered rules):** the positive control is not met: the rhythm of blood T cells that is established
  in living subjects (Dimitrov et al. 2009) is not recovered from RNA shares in the blood of fast-death donors, at a
  power of about 0.8 for the full published amplitude. By the registered reading, no T-cell statement from GTEx
  bulk tissue, A245's H245.2 included, is to be read as biology: A245's T-cell failure says the A229 hour is not
  seen in these data, not that tissue T-cell content lacks a daily rhythm. The glucocorticoid score points to the
  morning (fitted 7.9 h, beyond matched genes) but is not detected at this size; the tissue T-cell peak near 16 h
  stays a post hoc hint.
<!-- src: paper/exploratory/A248_glucocorticoid_tcells/output/{power,compartments,score_tests,score_descriptive,score_per_compartment,score_matched}.csv -->

## A247 (scored 2026-11-12): in glioma, cell classes of one tumour sample pair, but no more on clock genes than on matched genes (precondition P247.0 FAIL; H247.1-2 uninformative); MGMT does not follow the myeloid clock state (H247.3 FAIL)
Power note and note 1: `registered_experiments.md` (commits b6521de, 88fc3e5; the calibration was corrected before scoring and
a partial leak of the pairings is disclosed there). Extended GBmap, droplet libraries, pseudobulks with >= 100 cells; 10,000
permutations within study x assay; 200 matched and 200 responsive pseudo-clock sets. Scripts by a Sonnet agent; numbers
re-read in `tests.csv`.
- **P247.0 FAIL (myeloid vs oligodendrocyte, 53 patients / 8 strata):** pairing 0.499 vs null 0.373 (SD 0.068), p 0.035;
  matched sets median 0.506, q95 0.645 (p 0.54); responsive q95 0.628 (p 0.71): not clock-specific; after the foreign-RNA
  correction 0.465 (p 0.088).
- **H247.1 (malignant vs myeloid, 105 / 11):** 0.502 vs null 0.297 (SD 0.047), p 0.0003, Holm 0.0006; matched median 0.461,
  q95 0.594 (p 0.30); responsive q95 0.597 (p 0.39): pairing, not clock-specific; corrected 0.383 (p 0.037).
- **H247.2 (malignant vs oligodendrocyte, 42 / 5):** 0.432 vs null 0.326, p 0.086, Holm 0.086; matched q95 0.632; corrected
  0.295 (p 0.67): fail.
- **H247.3 FAIL (MGMT in malignant cells vs myeloid state, 85 / 10):** r 0.282 vs null 0.314 (SD 0.056, MDE 0.140), p 0.70;
  MGMT at the 22.5th percentile of matched malignant genes. **H247.4 (MGMT vs the malignant cells' own state, 93 / 10):**
  0.242 vs 0.302, p 0.88. The methylation-call test was not scored (power rule).
- **Descriptive:** T/NK vs myeloid 0.570 vs null 0.288 (p 0.0001), matched median 0.499, q95 0.679: not clock-specific
  (myeloid RNA bound in T/NK 0.258). Host gate on myeloid + oligodendrocyte anchors (51 patients): p_pairing 0.002,
  p_acrophase 0.064, pairing z 4.12 below both companion q95 (13.96, 16.37): no coupling estimate.
- **Reading (registered rules):** the precondition is not met, so H247.1-2 say nothing about a glioma clock tied to
  the patient's brain cells. What the data do show is that classes of one dissociated tumour share an expression state
  that matched non-clock gene sets carry as strongly as clock genes (matched medians 0.46-0.51 against permutation
  nulls of 0.30-0.37): same-sample covariation (ambient RNA, handling, patient state), not a clock readout. A
  within-sample design cannot isolate the clock here; the carcinoma pairings that are clock-specific (A151, A227) compare
  tumour with a separate normal sample. Dong et al.'s culture result is neither supported nor contradicted. MGMT
  expression shows no association with the clock state of the tumour's myeloid cells or of the malignant cells at an
  MDE of 0.13-0.14; with the precondition failing, whether that state carries body time is itself unshown, so this is
  "not detected", not evidence against a timed MGMT.
<!-- src: paper/exploratory/A247_glioma/output/{power,eligibility,tests,pooled,per_study,pseudo,matched_genes,descriptive,coupling,foreign_rna,foreign_clock_share,calibration_null_P247.0_all,calibration_null_H247.1_all,calibration_null_H247.2_all}.csv -->

## A249 (scored 2026-11-12): the chronotherapy window cannot be read in TCGA pairs (precondition P249.0 FAIL; H249.1-2 not shown); the tumour's clock-controlled genes follow the body time read in the normal sample, its drug-handling genes do not (descriptive); in the timed breast cohort the hour is seen in normal tissue more than in the tumour, for every gene set (H249.5, p 0.030 for set C)
Power note: `registered_experiments.md` (commit fc069bb). TCGA: 517 pairs, 10 cancers; S_N = normal-sample genes against the
tumour-sample clock state, S_T = tumour-sample genes against the normal-sample state; 10,000 permutations, 200 pseudo-clock
sets, 200 matched gene sets; null mean of S 0.120.
- **P249.0 FAIL (set O, clock-controlled genes, S_N):** 0.134 (null SD 0.011), p 0.096; pseudo-clock q95 0.156 (p_clock 0.48);
  matched q95 0.136. The clock state of the tumour sample does not order clock-controlled genes of the normal sample:
  the normal-side readouts below are uninformative.
- **H249.1 FAIL (set C, S_N):** 0.120 (p 0.50; MDE 0.023), below the pseudo-clock median (0.136). **H249.2 not shown:** S_T of
  set C 0.117 (p 0.70; MDE 0.018); D +0.003 [-0.024, 0.024].
- **H249.3 FAIL:** K S_N 0.128 (p 0.20), S_T 0.116 (p 0.67); F (761 drug targets) S_N 0.122 (p 0.41), S_T 0.122 (p 0.36).
- **Descriptive, the readable direction (body time from the normal sample):** the tumour's clock-controlled genes follow
  it: set O S_T 0.162 (p 0.0001; pseudo-clock q95 0.147, p_clock 0.015; matched q95 0.136, p 0.005), led by NFIL3 (r 0.206),
  RORC (0.195) and ARNTL2 (0.160); D for set O -0.028 [-0.051, 0.007]. Lung, thyroid, kidney and liver carry it (S_T 0.16-0.22).
  With the same reader, the tumour's drug-handling genes do not (set C 0.117, set F 0.122, both at the null).
  Same-sample readouts are higher throughout (0.18-0.24) and confounded.
- **H249.4 (single cell, 61 patients in 13 cohorts of 4 to 7):** nothing detected (set O S_N 0.525 vs null 0.492, p 0.13;
  set C S_T 0.534, p 0.13; no companion passed); blocks too small to carry a claim.
- **H249.5 (GSE233242, 42 patients, recorded hour):** mean cosinor modulus of set C 2.416 in normal tissue and 1.458 in
  the tumour, difference 0.958 (null SD 0.460), one-sided p 0.030. The same direction in every set: K 0.785 (p 0.065), O
  0.742 (p 0.082), F 0.556 (475 genes; normal 1.871, tumour 1.314; p 0.047). No matched-gene companion was registered here.
- **Reading (registered rules):** the window is not shown and, in TCGA, not testable: the tumour sample is too poor a
  reader of body time for the normal tissue's genes (the asymmetry noted at registration; the precondition fails). Read
  the other way, the result is informative about the tumour: its clock-controlled output genes are on the patient's
  body time, clock-specifically and beyond matched genes, which extends the tumour-normal pairing of the 14 clock genes
  (H138.4) to genes outside the waveform; its drug-handling genes and drug targets are not, at an MDE of 0.011-0.018
  against an excess of 0.042 for the clock-controlled set. Whether the same genes are on body time in normal tissue is
  the missing half. The timed breast cohort points that way (the hour is visible in normal tissue and near the null
  in the tumour) but for all gene sets alike, without a specificity control, in one cohort of 42 already used for
  A146 and A200: a hint that tumours lose hour-of-day structure broadly, not evidence of a drug-specific window.
<!-- src: paper/exploratory/A249_chrono_window/output/{power,score_tests,score_genes,score_descriptive,score_pseudo,score_matched}.csv -->

## A250 (scored 2026-11-12): the A249 asymmetry does not replicate as a clock-specific result in CPTAC (H250.1 detected, not clock-specific; H250.2 FAIL); NFIL3 alone replicates (H250.4); the window stays unreadable (H250.5)

- Data: 686 CPTAC-3 tumour-normal pairs, 7 types (A186 pool); one scoring run (ok, 403 s). Null mean of S about 0.086.
  Power rules met (MDE of S_T: set O 0.021, set C 0.015; Delta(O, C) 0.023).
- **H250.1 detected, not clock-specific (registered reading: a patient-level state shared by the two samples, not shown
  to be the clock):** set O in the tumour against the normal-sample state, S_T 0.109 (null 0.086, SD 0.008), p 0.006;
  pseudo-clock sets median 0.107, q95 0.131 (p_clock 0.42); matched genes q95 0.108 (p_matched 0.055). The excess over
  the null is 0.023, against 0.042 in TCGA.
- **H250.2 FAIL:** set C is detected in the tumour (S_T 0.097, p 0.036; not clock-specific, p_clock 0.57; not beyond
  matched genes, p 0.46), and Delta(O, C) is 0.012 (p 0.099; 90 % interval -0.005 to 0.031; MDE 0.023).
- **H250.3 FAIL:** Delta(O, F) 0.012 (p 0.036, Holm 0.072); Delta(O, K) 0.012 (p 0.106). Set F S_T 0.097 (p 0.008,
  p_clock 0.58); set K 0.097 (p 0.066).
- **H250.4:** NFIL3 PASS (r 0.155, null 0.086, p 0.0007, Holm 0.002; implied peak 18.9 h, TCGA 19.0 h); RORC FAIL
  (0.091, p 0.39); ARNTL2 FAIL (0.111, p 0.092, Holm 0.18; implied peak 21.8 h, TCGA 20.7 h).
- **H250.5 (window, second attempt): precondition FAIL again:** set O S_N 0.115, p 0.004, not clock-specific (pseudo
  q95 0.125, p_clock 0.16; beyond matched genes, p 0.010). Set C S_N 0.095 (p 0.11); D -0.002 (-0.020 to 0.018).
  The normal half stays unreadable from tumour-sample states.
- Descriptive: per type S_T of set O from 0.064 (SKCM) to 0.153 (PDAC), LUAD (215 pairs) 0.089; STAD (27 pairs, not
  pooled) 0.212. Same-sample readouts (confounded): set O 0.233 normal, 0.218 tumour. By G2M tertile, S_T of set O 0.180
  (low) and 0.174 (high), set C 0.148 and 0.164 (smaller strata, higher null: compare within tertile only). TCGA and
  CPTAC combined by pairs: set O S_T 0.132, set C 0.106. No purity strata (no table on disk).
- **Reading (registered rules):** in CPTAC every gene set of the tumour, clock-controlled or not, follows the state
  read in the normal sample slightly above the permutation null, and states built from non-clock gene sets do as well
  as the clock state (pseudo median 0.107 against 0.109 observed). The cross-sample association is therefore a
  patient-level state shared by the two samples (as H186.1 found for the pairing itself in CPTAC), not shown to be the
  clock; the clock-specific excess of TCGA is not reproduced, and the asymmetry between clock output and drug handling
  shrinks from 0.045 to 0.012. The A249 descriptive result stays descriptive and is not to be worded as a finding.
  What does replicate is one gene: NFIL3 in the tumour follows the normal sample's clock state in both cohorts, at the
  same implied hour (about 19 h).
<!-- src: paper/exploratory/A250_cptac_window/output/{power,score_tests,score_contrasts,score_genes,score_genes3,score_descriptive,score_pseudo,score_matched}.csv -->

## A251 (scored 2026-11-12): in living women, skin at three body sites and in two layers shares one clock state, clock-specifically (H251.1 PASS, H251.2 PASS x2); the untimed dermis-epidermis offset matches the timed one (H251.4 PASS); site offsets are within about 1.5 h but one interval exceeds the bound (H251.3 FAIL, 5 of 6)

- Data: Wu 2020 (GSE139301 epidermis, GSE139300 dermis), 154 women, untimed clinic biopsies, six compartments each in
  its own array batch; 146-150 women per pair. All six compartments pass the reference gate (clock share 0.59-0.76,
  null mean 0.09-0.12). One scoring run. Null mean of the pooled r 0.072.
- **P251.0 PASS (timed men, seen in the code test, see the power note):** the state follows the collection time in
  epidermis (modulus 1.62, null q95 0.45) and dermis (0.99, q95 0.32), p 0.0001 each (20 men, 79 samples per layer).
- **H251.1 PASS, clock-specific (primary):** the same layer at different sites, six pairs pooled: r 0.628 (null SD
  0.017, p 0.0001); matched pseudo-clock sets median 0.198, q95 0.380; responsive sets median 0.368, q95 0.420 (0 of 200
  sets reach it in either companion). Without PER1 and NR1D1: 0.637.
- **H251.2 PASS x2, clock-specific:** layers at the same site 0.813 (matched q95 0.373, responsive median 0.657, q95
  0.770); layers at different sites 0.594 (matched q95 0.212, responsive q95 0.382); p 0.0001 each. Without PER1 and
  NR1D1: 0.759 and 0.585.
- Per pair r 0.516 to 0.877 (15 pairs); the lowest are pairs with the buttock epidermis.
- **H251.3 FAIL by the registered rule (5 of 6 within +/- 1.5 h):** same-layer site offsets: dermis forearm-buttock
  -0.14 h (-0.48 to 0.21), forearm-cheek -0.85 (-1.11 to -0.58), buttock-cheek -0.71 (-0.99 to -0.40); epidermis
  forearm-buttock +0.73 (0.20 to 1.32), forearm-cheek -0.78 (-1.11 to -0.43), buttock-cheek -1.54 (-2.04 to -1.08). The
  cheek is the latest site in both layers (first minus second compartment). *Erratum 2026-11-23 (wording only): offsets
  are first minus second in the clocksync convention (> 0: advanced, peaks earlier; `compute_a251.offsets`, arg za conj zb),
  so the cheek is the most advanced site, its clock genes peaking earliest, not latest.*
- **H251.4 PASS:** dermis minus epidermis at the forearm: women (untimed) -1.58 h (-1.81 to -1.34); timed men -1.98 h
  (-3.15 to -0.95). The untimed estimate lies inside the timed interval; difference 0.40 h. Across all nine
  layer pairs the dermis is 0.5 to 2.8 h behind the epidermis.
- **Reading (registered rules):** in living people sampled once in clinic hours, skin compartments a metre apart are
  on one internal time: the pairing is far above the permutation null, above every expression-matched set and above
  every stress-responsive set, and does not depend on PER1 or NR1D1. Because each compartment was hybridised in its own
  batch and the cross-site pairs come from separate biopsies, a shared array batch or a shared biopsy cannot produce
  it. The responsive companion is itself high (median 0.37, and 0.66 within one biopsy): women also share a
  stress-type state across sites, which the clock pairing exceeds. Offsets between sites are small (0.1 to 1.5 h) but
  not all within the registered bound: the cheek runs later than the other two sites [erratum 2026-11-23: earlier, i.e. advanced; see above]; whether that is sun exposure, site
  or batch cannot be told (each compartment is one batch, and a gene-specific batch shift cannot move an offset
  estimated from centred states, but tissue composition can). The untimed layer offset reproduces the timed one within
  0.4 h. This is the living-cohort support the robustness arm lacked; it is one cohort, women only, on arrays, and age
  (20-74 y) is not public, so no age contrast is possible.
<!-- src: paper/exploratory/A251_skin_sites/output/{gates,power,score_tests,score_pairs,score_offsets,score_timed,score_descriptive}.csv -->

## A252 (scored 2026-11-12): within a living person, the clock state of adipose tissue and muscle advances with the hours between two biopsies, in the predicted direction and clock-specifically (H252.1 PASS, H252.2 PASS, H252.3 adipose PASS, H252.4 PASS x2); the implied time of day tracks the protocol's within about an hour

- Data: eight repeat-biopsy studies (184 + 30 units, 334 within-person changes), none used before; one scoring run.
  rho is the amplitude of the rotation in units of the GTEx waveform (0 = no rotation).
- **H252.1 PASS, clock-specific (primary; GSE118280, adipose, exact clock times 08:40, 12:20, 19:00; 15 men x 2
  visits):** rho 1.265 (null SD 0.215, p 0.0001); matched sets q95 0.109, responsive sets median 0.232, q95 0.738 (no set
  of 200 reaches it in either). By interval: +3.67 h rho 1.44, +10.33 h rho 1.23.
- **H252.2 PASS, clock-specific (seven studies with assumed morning times):** rho 2.970 (null SD 0.237, p 0.0001);
  matched q95 0.205; responsive median 0.401, q95 1.201. With the assumed first-biopsy time moved by -2 / +2 h: 2.67 / 2.48.
- **H252.3:** adipose (2 studies) PASS, rho 2.221 (p 0.0001; matched q95 0.351; responsive q95 2.120, p_responsive
  0.045: clock-specific, narrowly). Muscle (5 studies): **not scored** (power rule, MDE 0.76); descriptive value 3.218
  (responsive q95 1.147).
- **H252.4 PASS x2 (without PER1 and NR1D1):** GSE118280 rho 1.058 (responsive q95 0.292); seven studies 3.566
  (responsive q95 1.516); Holm p 0.0002 each.
- Direction (descriptive): the clock time of the interval's midpoint implied by the mean change, against the
  protocol's, per study and time point (14 points): GSE118280 10.1 vs 10.5 h and 11.1 vs 13.8 h; GSE252357 9.3 / 11.4 /
  13.3 vs 9.9 / 11.1 / 13.6 h; GSE151066 9.4 / 10.9 vs 9.3 / 10.8 h; GSE273041 9.6 / 10.7 vs 10.5 / 11.5 h; GSE120862 9.3 /
  11.3 vs 10.0 / 11.5 h; GSE231509 10.7 vs 10.0; GSE142401 11.8 vs 11.0; GSE303983 10.4 vs 10.1 h. Thirteen of 14 within
  1 h; within a study the implied time moves later with the later biopsy in all six studies with two or more
  intervals.
- Size (descriptive): rho per study and interval 1.2 to 6.5; adipose 1.2-1.4 with exact times, 1.2 and 4.8 with
  assumed times; muscle 1.5 to 6.5. The clock-plane change is larger than a rotation with the GTEx amplitude.
- Rest against exercise (descriptive, no claim): GSE252357 rest arm (5 persons, no intervention but breakfast) rho
  2.886, exercise arm 3.095, difference p 0.50; GSE120862 resting leg 2.646, exercised leg 2.300, p 0.41.
- Surgery and transplantation (descriptive, no clock or elapsed time): modulus of the mean change in atrium during
  bypass surgery 1.68 (matched q95 1.78; responsive median 7.71), in muscle around knee surgery 0.79 (matched q95
  1.15; responsive median 1.48), in the liver graft across transplantation 1.02 (matched q95 1.90; responsive median
  2.41): no clock-plane change beyond matched gene sets, while stress-responsive sets move several-fold more.
- **Reading (registered rules):** the first within-person validation: between two biopsies of one day the state moves
  along the clock's circle, in the direction and at the time of day the protocol predicts, in two tissues and eight
  laboratories, beyond expression-matched and stress-responsive gene sets and without PER1 and NR1D1. The resting arm
  moves as the exercised one, so at these n exercise is not needed for it and does not detectably add to it. Limits as
  registered: every interval but one holds a meal or an intervention; clock times of seven studies are protocol
  guesses (the result is stable to +/- 2 h); rho mixes amplitude and rotation. rho well above 1 says the living
  tissues' clock amplitude exceeds the GTEx waveform's (post-mortem, attenuated: review M2), or that the morning hours
  carry a larger change than a sinusoid allows; it is not an estimate of elapsed time. The three surgical datasets
  show no clock-specific change, but without times they cannot say whether the clock stalls under surgery or cold
  storage.
<!-- src: paper/exploratory/A252_elapsed_time/output/{design,power,score_tests,score_datasets,score_descriptive}.csv -->

## A253 (2026-11-12): clock relations between brain regions in depression and psychosis are not readable in public post-mortem cohorts: the controls' region pairing is matched by stress-responsive genes (P253.0 FAIL; H253.1-4 NOT SCORED)

- Data: GSE80655, GSE102556, GSE53987, GSE202537; 68 depressed, 122 psychotic and 65-79 control donors, 3 to 6 regions each.
  The depression patients of Li 2013 (GSE45642) were never deposited (controls only).
- **P253.0 FAIL (controls only):** regions of one control donor share a clock-plane state: R 0.734 against a
  donor-permutation null of 0.182 (SD 0.026, p 0.0001) and above every expression-matched pseudo-clock set (median
  0.387, q95 0.549, p 0.005). Stress-responsive sets reach the same level (median 0.648, q95 0.758; p_responsive 0.139).
- Power was not the limit: MDE of the pairing difference 0.134 (depression) and 0.168 (psychosis); of the offset
  statistic 0.34 and 0.41 h.
- **H253.1-4 not scored** (registered rule). The scoring stage was not run. Leak disclosed in the power note: the
  code test printed the pairing per true group (not offsets); it is not reported here and supports no statement.
- **Reading:** between regions of a post-mortem brain, donors share a state that glucocorticoid and immediate-early
  genes carry as strongly as clock genes (agonal state, pH, death circumstances), so a case-control difference in
  region pairing could not be attributed to the clock. This is the brain counterpart of the tumour lesson (a strong
  non-clock axis between persons) and agrees with WP0 item 8 and A141. A human test of internal desynchrony in
  depression needs living tissue or a design in which the responsive companion is low; public post-mortem brain is
  not it. The skin cohort (A251) shows the contrast: there the clock pairing exceeds the responsive sets.
<!-- src: paper/exploratory/A253_psych_regions/output/{design,precondition,power}.csv -->

## A255 (scored 2026-11-12): fat and skin of a living woman share a clock state that her cultured cell line does not carry (H255.1 PASS, H255.2 PASS); two tissues taken at one operation share it in two surgical cohorts (H255.3 PASS); co-twins share it across persons beyond their cell lines (H255.4 PASS x2), identically in monozygotic and dizygotic pairs

- Data: MuTHER (E-TABM-1140), 680 women with fat and skin (one punch biopsy), 794 / 679 with a lymphoblastoid cell line;
  GSE73439 (subcutaneous and omental fat, 51 women) and GSE133979 (rectus muscle and subcutaneous fat, 32 patients).
  GSE24883 left the pool at the power note (gates fail). One scoring run.
- Gates: fat p 0.010, skin p 0.015; the cell line has no clock-shaped covariance (share 0.042 against 0.065, p 0.75).
- **H255.1 PASS, clock-specific (primary):** fat-skin r 0.533 (null 0.034, SD 0.018, p 0.0001); matched sets median 0.121,
  q95 0.304 (p 0.010); responsive sets median 0.200, q95 0.325 (p 0.005). Without PER1 and NR1D1: 0.387. Offset fat minus
  skin +0.11 h (-0.15 to 0.41).
- **H255.2 PASS (primary; beyond the genome):** fat-cell line r 0.042 (p 0.24), skin-cell line 0.066 (p 0.055), neither
  above its companions (matched medians 0.147 and 0.105). r(fat, skin) minus their mean: 0.479 (90 % 0.432 to 0.510).
- **H255.3 PASS, clock-specific (surgical cohorts pooled, 83 persons):** r 0.686 (null 0.139, SD 0.052, p 0.0001; matched
  q95 0.396, responsive q95 0.530). GSE73439 0.648 (responsive q95 0.570), GSE133979 0.747 (responsive q95 0.533).
- **H255.4 PASS x2 (co-twins):** co-twin states pair in fat 0.520 (357 pairs) and skin 0.589 (271), p 0.0001 each; in the
  cell line 0.384 (359 pairs, p 0.0001). Tissue minus cell line: fat +0.135 (0.058 to 0.211), skin +0.205 (0.130 to
  0.282). By zygosity: fat 0.522 (MZ) / 0.520 (DZ); skin 0.587 / 0.591; cell line 0.407 / 0.373.
- Descriptive, unexplained: without PER1 and NR1D1 the skin-cell line pairing is 0.248 (fat-cell line 0.074); no
  companion was computed for these.
- **Reading (registered rules):** two organs of a living person, sampled once without a recorded time, share a
  clock state, beyond expression-matched and stress-responsive genes, in a twin cohort and in two surgical cohorts.
  The woman's own cell line, which carries her genome but not her body time, does not share it: the pairing is not
  constitutive clock-gene expression. With A251 this is the living counterpart of the GTEx statement (A109, A150).
  Co-twins share a clock state across persons to the same degree whether monozygotic or dizygotic, so that sharing
  is not genetic; twins usually attend together, which fits shared visit time, but their cell lines also pair
  (0.38), so part of it is shared laboratory handling of a twin pair, and only the excess over the cell line (0.14
  and 0.20) is attributable to the visit. Limits: fat and skin are from one punch biopsy (the responsive companion,
  q95 0.325, bounds what shared handling of the biopsy gives); women only; arrays; the 0.248 above shows that the
  cell-line control is not flat under every gene subset.
<!-- src: paper/exploratory/A255_living_tissues/output/{gates,power,score_tests,score_twins,score_descriptive}.csv -->

## A254 (scored 2026-11-12): no psoriasis-sized clock delay of the lesion in periodontitis, osteoarthritis, carotid plaque, fibroid or hidradenitis (H254.1 FAIL, powered: pooled offset +0.01 h, MDE 0.69 h; H254.2 FAIL x4); lesions carry less clock-plane variance than the person's own tissue, clock-specific at the line when pooled (H254.3)

- Data: five paired datasets, 184 persons; every non-lesional tissue passes the clock gate (p 0.005-0.015). Lesion and own
  tissue share a clock state in all (coherence 0.45 to 0.77). One scoring run.
- **H254.1 FAIL (primary, powered):** pooled offset lesion minus non-lesion +0.01 h (null SD 0.275 h; one-sided delay p
  0.52; MDE 0.69 h); companions: matched q05 -0.75 h, responsive q05 -0.67 h. A delay of the psoriatic size (1.4 h) is
  excluded for these lesions taken together.
- **H254.2 FAIL x4 (each disease, two-sided, Holm over four):** periodontitis -0.15 h (90 % -0.66 to 0.47; p 0.66);
  osteoarthritis +0.04 (-1.65 to 0.94; p 0.96); carotid plaque +0.53 (-1.05 to 2.34; p 0.56); hidradenitis -1.07 (-1.79 to
  -0.14; p 0.085, Holm 0.34; one-sided delay p 0.042; not clock-specific: companions' |offset| q95 6.9 and 4.4 h).
  Uterine fibroid +1.19 h (-0.43 to 3.01), not scored (MDE 2.56 h).
- **H254.3 (clock-plane variance, lesion over non-lesion, log ratio r):** pooled r -0.417 (SE 0.088; p 0.0001), below the
  q05 of both companions (matched -0.404, responsive -0.272): **PASS, clock-specific at the line.** Per disease (Holm
  over four): periodontitis -0.453 (p 0.0001, Holm 0.0004; inside the companions' range: matched q05 -0.864, responsive
  q05 -0.495); carotid plaque -0.545 (p 0.003, Holm 0.010; matched q05 -0.670); osteoarthritis +0.089 (p 0.74);
  hidradenitis -0.227 (p 0.37): none clock-specific alone. Fibroid -0.791 (p 0.025; outside both companions; not
  scored).
- Descriptive: the second periodontitis series (GSE10334, 64 persons; overlap with the first undocumented) gives offset
  -0.13 h (-0.71 to 0.63) and r -0.481 (p 0.0003). Barrett mucosa against the same patient's squamous oesophagus
  0.00 h and r +0.151; against gastric cardia +0.64 h and r +0.124. Without PER1 and NR1D1 the offsets stay within
  their intervals.
- **Reading (registered rules):** lesions of four non-skin tissues stay on the time of the person's own tissue:
  no shift, with a pooled bound of 0.7 h. The psoriatic delay is therefore not a general property of chronic
  lesions; within this set only hidradenitis, the one inflammatory skin disease, points the psoriasis way (-1.07 h,
  20 persons, not significant after Holm, not clock-specific), which is consistent with the delay being a feature of
  inflamed skin (A183: psoriasis and eczema alike). What lesions do share is a lower clock-plane variance (pooled
  ratio exp(-0.417) = 0.66), the direction seen in tumours; it is clock-specific only when pooled and at the line,
  and matched gene sets also lose variance in single diseases, so it reads as less clock signal in lesional tissue,
  cause open (composition, damping). For the claim "where clocks break, they break locally": the break is a loss
  of clock signal in the lesion, not a general desynchrony of the lesion from the body; the phase shift is specific
  to inflamed skin.
<!-- src: paper/exploratory/A254_paired_lesions/output/{power,score_datasets,score_pooled,score_descriptive}.csv -->

## A256 (scored 2026-11-13): the within-person rotation of the clock state replicates in seven fresh studies (H256.2 PASS, clock-specific; H256.3 PASS: 17 of 20 cells within 2 h of the protocol's time of day); the arms in which nothing is done are too few to be scored (H256.1 NOT SCORED, 13 units) but point the same way; interventions do not detectably add to it (H256.4)

- Data: seven studies not opened before the design addendum, 166 units in 17 arms (muscle, subcutaneous fat); one
  scoring run. rho as A252 (amplitude of the rotation in units of the GTEx waveform).
- **H256.2 PASS, clock-specific (all fresh units):** rho 2.533 (null SD 0.239, p 0.0001, Holm 0.0002); matched sets q95
  0.275, responsive sets median -0.380, q95 0.779 (no set of 200 reaches it in either). Without PER1 and NR1D1: 2.967. Per
  study: 3.39 (GSE40589; log base of the deposited ratios unknown, so its size is not comparable), 0.86 (GSE9105, no
  PER1 on the array), 2.08, 3.03, 2.28, 1.81, 0.81.
- **H256.3 PASS (direction):** the time of day implied by the mean change lies within 2 h of the protocol's midpoint in
  17 of 20 dataset x arm x time-point cells (5.2 expected under permuted labels; p 0.0001). The three outside: GSE59088
  resistance (17.8 vs 14.3 h) and endurance (16.3 vs 14.3 h, 2.03 h off), GSE163356 session 9 (13.4 vs 11.0 h).
- **H256.1 NOT SCORED (power rule: 13 units, MDE 2.57 against a limit of 1.0).** Descriptive: rho 2.668 (permutation p
  0.0003; matched q95 0.560, responsive q95 0.512; without PER1 and NR1D1 1.873). Saline alone for 4 h in fasted, supine
  persons: rho 2.92, implied midpoint 11.1 h against 11.0 h (7 units); rest alone from 13:00 to 15:30: rho 1.95, implied
  12.5 h against 14.3 h (6 units).
- **H256.4 (interventions against nothing done; GSE40589 and GSE59088):** difference of rho +0.653 (null mean 0.206, SD
  0.472; two-sided p 0.35): not detected, on 13 nothing-done units. By class over all studies (descriptive): nothing
  done 2.67, recovery after exercise 2.26, intervention 2.57.
- **Reading (registered rules):** A252 replicates on studies fixed before their data were opened: the clock state
  moves along the circle between two biopsies of a day, in the predicted direction, beyond matched and
  stress-responsive gene sets and without PER1 and NR1D1, in 16 studies overall (A252 and A256). Whether elapsed time
  alone suffices could not be tested at the registered power; what there is (saline and rest arms, 13 units) shows
  the same rotation as the intervention arms, and insulin clamps, GIP, exercise and recovery do not detectably change
  it. Together with the rest arm of A252 (5 persons, rho 2.89) the evidence for "time alone" is three small arms that
  agree, not a registered pass. Limits as registered; fasting continues in every nothing-done arm.
<!-- src: paper/exploratory/A256_elapsed_replication/output/{design,power,score_tests,score_cells,score_descriptive}.csv -->

## A257 (scored 2026-11-13; not blind, A124 scored the same data): after a 12 h shift of behaviour, the human muscle clock moves by only 2.7 h and its daily excursion falls to about a quarter (H257.1 PASS, H257.2 partial shift, H257.3 PASS)

- Data: GSE106800 (Wefers 2018), vastus lateralis at 07:00 and 19:00 after an aligned and after a 12 h-misaligned
  protocol; 12 aligned units, 11 subjects with all four arrays. One scoring run.
- **H257.1 PASS, clock-specific (positive control):** in the aligned condition the state moves from 07:00 to 19:00
  along the predicted chord: rho 1.336 (exact null over 4,096 swaps, SD 0.393, p 0.0002); matched sets q95 0.137,
  responsive sets median 0.014, q95 0.458 (no set of 200 reaches it).
- **H257.2 partial shift (registered reading):** the change from 07:00 to 19:00 under misalignment is turned by s =
  +2.68 h against the aligned change (90 % interval 1.81 to 3.59 h): it excludes 0 and 12 h and is not inside +/- 3 h.
  Without PER1 and NR1D1: +3.51 h. Per subject: 10 of 11 positive (-1.39 to +6.57 h). The muscle clock is about 2.7 h
  ahead of its aligned position (a 12 h shift has no direction of its own, so "towards behaviour" cannot be read from
  the sign).
- **H257.3 PASS (damping):** projected on the aligned chord, the misaligned change is R = 0.233 of the aligned one (90 %
  0.060 to 0.414); rho 1.322 vs 0.309, exact one-sided p 0.001; the difference (1.01) exceeds the q95 of the matched
  (0.238) and responsive (0.509) companions. A rotation of 2.7 h alone would leave cos(40 deg) = 0.77 of it: most of the
  loss is amplitude.
- **Reading:** in human muscle, three days of 12 h-shifted sleep, light, activity and meals move the clock by only
  2.7 h, against 12 h of behaviour, and flatten its 07:00 to 19:00 excursion to about a quarter. The muscle clock
  therefore stays close to central (light-entrained) time and loses amplitude rather than following behaviour: tissue
  misalignment in this design is a damped peripheral clock that lags the shift, not a peripheral clock re-set by
  feeding and activity. This formalises the A124 descriptive (same direction, cosine 0.91; amplitude 3.20 -> 1.67) with
  intervals and companions; the angle here is larger (2.7 h against about 1.6 h implied by the cosine) because the
  gene table includes the eight platform aliases and the projector is fitted over all arrays. Limits: not blind; two
  clock times per condition; 11 men; day 3 of the protocol.
- **Scope (note 2026-11-13):** A257 uses only the clock-plane projection of the method (per-sample state on the GTEx
  waveform) with a timed, within-person reference; it has one compartment, so no host state, coupling, relative
  amplitude, inter-compartment offset, gate or d* is estimated. H257.1 validates the projection step (the state tracks
  clock time within a person, as A252 / A256); it is not a positive control for host-referenced inference, and as a
  not-blind analysis it is not a pre-specified control of anything. Its place is a biological result (peripheral
  response to misalignment), with H257.1 beside A252 / A256 in Fig. 1.
- **Correction (note 2026-11-19, PAPER_PLAN_2026-11-19 section 8; the scored text above is kept):** R = 0.233 is the
  misaligned change projected on the aligned chord, not the excursion ratio. The excursion ratio is R / cos s =
  0.2335 / cos(2.677 h x 15 deg/h = 40.2 deg) = 0.31 (derived from `score_tests.csv`: R 0.2335, s 2.677 h). "Flatten its
  excursion to about a quarter" should read "to about a third"; a 2.7 h rotation alone would leave 0.77.
<!-- src: paper/exploratory/A257_nightshift_muscle/output/{score_tests,score_descriptive}.csv -->

## A244 (scored 2026-11-13): the clock-shaped covariation of the fetal adrenal is carried by the cortex (H244.1 cortex PASS, also after foreign-RNA correction); the medulla shows none at full adult power (H244.1 medulla FAIL); the shared state of the two classes is not estimable (H244.2)

- Data: six single-cell cohorts of A242, marker-defined classes (design addendum: cortex 34 fetuses, medulla 30, in 6
  cohorts each); adult pass rate at the full clock 1.000 for both classes, 0.995 / 0.910 at 0.75 and 0.520 / 0.335 at
  0.5. One scoring run (54 s).
- **H244.1 cortex PASS:** C 0.608 (p 0.0038, Holm 0.0076; null q95 0.351), above every matched (q95 0.300) and responsive
  (q95 0.451) set; corrected for medullary RNA 0.611 (p 0.0038); only fetuses with a >= 200-cell reference 0.764 (p
  0.0002, descriptive). Leave-one-gene-out 0.438 (without ARNTL) to 0.808.
- **H244.1 medulla FAIL:** C -0.048 (p 0.41; matched 57.5th, responsive 62.5th percentile); corrected -0.049; never
  positive in any leave-one-cohort-out (-0.19 to -0.01) or window. At the adult power table this excludes a medullary
  clock of 0.75 times the adult size (pass rate 0.910) or more, not one of half the size (0.335).
- **H244.2 not estimable** (design addendum: adult pass rate 0.395 with muscle as the second tissue); not computed.
- Descriptive: stroma C 0.322 (p 0.079). Cortex by window: from 12 weeks 0.169 (13 fetuses, p 0.27; adult pass rate
  0.910), before 12 weeks 0.652 (19 fetuses, p 0.0028, clock-specific). Cortex without one cohort: 0.369 to 0.743;
  without kameneva p 0.054 (24 fetuses), without wang 0.477 (p 0.020, not above the responsive q95 0.480); clock-specific
  in the other four. Position among adult subsamples at the full clock: cortex at the 1st percentile (adult median
  0.848), early cortex at the 14.5th. Unsigned R cortex 0.630, medulla 0.152. The registered depth-matched comparison
  of the two classes was not run.
- **Reading (registered rules):** "Cortex PASS, medulla negative at an adult pass rate >= 0.8: the cortex carries it; a
  medullary clock of the size given by the power table is excluded." The whole-adrenal signal of A242 lies in the
  steroidogenic cells, the lineage that carries the adult adrenal clock and makes glucocorticoids, not in the
  neural-crest medulla; class proportions do not explain it (the medulla, at the same depth range, shows nothing). The
  cortical signal is below the adult size (1st adult percentile), consistent with A242, and is strongest before 12
  weeks; from 12 weeks it is not detected in 13 fetuses although the adult power there is 0.91, which qualifies "from
  12 weeks" of A242 for the cortex alone (descriptive, not a registered test). It leans on the two largest cohorts.
<!-- src: paper/exploratory/A244_fetal_adrenal_classes/output/{design_classes,foreign_cohort,feasibility_passrates,l1,verdicts,verdicts_a244,position}.csv -->

## A258 (scored 2026-11-13): group contrasts are biased by unequal host spread up to the size of the equivalence bound; the arm-specific prior removes most of the bias but is not accepted by the registered criteria (V258.2 FAIL); the label-permutation test is anti-conservative at unequal spread; 3 to 55 rows of the equivalence families are flagged

- Simulation (`simulate_a258.py`): 36 design cells x 100 seeds (bias), 12 cells x 30 seeds with 200 permutations (size
  and power, 80 donors per arm); three runs, 435 s, 794 s, 795 s. Unit test of the correction run on 6 seeds at one
  design before the registered simulation (code check, disclosed).
- **V258.1 (bias of the pooled estimator, target, 80 donors per arm):** none at equal spread (-0.006 to +0.001). Case
  window 8-16 h (observed log V ratio -0.33 to -0.48): -0.015 to -0.038. Window 10-14 h (-0.79 to -1.25): -0.050 (anchor
  noise 0.002), -0.113 (0.05), -0.159 (0.15). Window 11-13 h (-0.97 to -1.62): -0.125, -0.182, -0.269. The anchor itself
  is biased as much (to -0.28).
- **The arm prior** brings the same cells to -0.033 to +0.018 (80 donors).
- **V258.2 FAIL (registered criteria):** (i) bias: 8 of 51 cells with log V ratio >= -1.0 exceed 0.02 (4 at >= 80 donors:
  largest -0.069 for the target and -0.080 for the anchor, both at 40 donors in the 11-13 h window with anchor noise 0.15);
  (ii) RMSE at equal spread: at most 1.058 times the pooled one (passes); (iii) power: at equal spread none lost (0.567 /
  0.367 / 0.400 against 0.567 / 0.367 / 0.367), but at unequal spread the pooled estimator rejects more (up to 0.567
  more), because its bias points the same way as the simulated loss; the registered criterion, taken over all cells,
  fails on that. No corrected column is added (part 2 not run).
- **Size of the label-permutation test (descriptive; 30 seeds per cell, nominal 0.05):** at equal spread 0.000 to
  0.033 for both estimators; at unequal spread the pooled test rejects in 0.033-0.300 (target) and up to 0.833 (anchor)
  of null data sets, the arm prior in 0.000-0.300 and up to 0.267. Permuting the labels exchanges the arms' spreads,
  so the null "same coupling, different spread" is not the one the permutation generates.
- **V258.3 (rows flagged, expected bias > 0.05 at the row's observed log V ratio; anchor noise 0.05 / 0.15):** A231 (age)
  0 / 3 of 54 rows (median |log V ratio| 0.18; largest expected bias 0.031 / 0.061); A233 (five chronic conditions) 6 /
  33 of 302 (0.20; 0.074 / 0.149); A235 (metabolic offsets; its claims are on offsets, d* is secondary there) 5 / 55 of
  262 (0.19; 0.062 / 0.124). A row is flagged from |log V ratio| 0.38 upwards at anchor noise 0.15. A164 / A202 / A208
  have no spread columns and are not mapped. Per-row table: `output/v258_3_rows.csv`.
- **Reading (registered rules):** "A failed V258.2 leaves the flag table as the answer." The pooled equivalence
  estimates of age and chronic conditions are exposed by at most about 0.03 (typical rows) and 0.06-0.15 (the flagged
  rows), always in the direction "the arm with the narrower spread reads lower": the pooled verdicts stand, the flagged
  rows are marked in the tables, and no sentence may say that d* removes the sampling spread. The arm prior is an
  improvement for 80 or more donors per arm and moderate spread differences, but it does not meet the registered bar,
  and neither estimator gives a calibrated permutation p at unequal spread: contrasts with |log V ratio| above about
  0.4 (31 % of the exploratory rows, review M1) have anti-conservative p values. A spread-matched null (permutations
  within strata of host spread, or a parametric bootstrap at each arm's own spread) is what is missing; it is a new
  amendment.
<!-- src: paper/exploratory/A258_spread_bias/output/{summary_bias,summary_size,summary_power,v258_1_bias,v258_2_acceptance,v258_size,v258_power,v258_3_rows,v258_3_summary}.csv -->

## Roadmap 2026-11-15, item A5: A91 on current code reproduces the archive exactly (no erratum)

`results/figure_3/scripts/gtex_cues/compute_internal_external.py` rerun on `clocksync` of 2026-11-15 (subsample 150 check,
then full through `run.sh`, 244 s, peak 2.4 GB); the 26-tissue selection table was recomputed (the cache was removed).
H91: 12 tissues, median lock difference 1.232, p 0.0005, PASS, as archived. All 312,018 gene calls identical; class
counts and noise floors identical in all 26 tissues (skeletal muscle 1,480 clock-following / 113 time-following / 676
both, floors 37 / 31.5 / 8.5; liver 2,717 / 5 / 326); lock values within 1.4e-9. A27 D5 summaries and host-vs-time
(R 0.631 fast deaths, 0.142 ventilated) identical. Still missing for the main text: matched and stress-responsive
companions (roadmap T2). Archived tables kept in `gtex_cues/archive_2026-10-25/`.
Reading limit (added 2026-11-15, after review): in unshifted GTEx donors the body time and the person's behaviour move
together, so A91's "clock-following" genes (B) follow the body time *or* behaviour; they cannot be read as "controlled
by the local clock". That partition is identifiable only where clock and behaviour are misaligned (A261 and phase T,
G2 (iii)).
<!-- src: paper/results/figure_3/output/gtex_cues/{a91_verdict,a91_tissues,a91_genes,diag5_summary,diag5_host_vs_time}.csv; archive_2026-10-25/ -->

## Roadmap 2026-11-15, item A1 (a, b): the scored offset contrasts in hours, and their spread bias (simulation)

- **(a) Table** (`exploratory/roadmap_A1_offset_bound/output/offsets_hours.csv`; no new statistic): women minus men
  0.043 h (90 % -0.104 to 0.189; MDE 0.166 h); old minus young -0.045 h (-0.274 to 0.184; 0.199 h); T2D 0.023 h
  (-0.252 to 0.299; 0.366 h) and obesity 0.202 h (-0.232 to 0.635; 0.638 h) in fast deaths; all four scored as
  equivalent within 1 h; pseudo-clock |offset| q95 1.290, 1.544, 2.716, 4.948 h. A233 (five chronic conditions)
  stores the module offset and p only (0.036 to 1.994 h, p 0.112 to 0.959), no interval or MDE: "not detected", not an
  equivalence.
- **(b) Simulation** (A258 design and code, `offset_diff_h` of `Dataset.contrast`; not a registered hypothesis; reading
  fixed in the script before the run): mean bias of the pooled-prior offset at most 0.061 h (case window 6-18 h) and
  0.060 h (8-16 h; observed log V ratio about -0.4), 0.106 h at 10-14 h (log V -1.0) and 0.228 h at 11-13 h (log V
  -1.3); Monte Carlo SEM about 0.06 h. Permutation size of p_offset (nominal 0.05, 30 seeds per cell, both priors): 0.000-0.133 at
  6-18 h and 8-16 h, 0.033-0.333 at 10-14 h and 11-13 h. A known +0.5 h shift is recovered within 0.16 h except at
  11-13 h (up to 0.34 h). The observed log V ratios of the scored families lie between -0.76 and +0.61 (A258 V258.3
  rows; median |log V| 0.18-0.20). **Reading:** the A208, A231 and A235 offset equivalences stand (an expected bias of
  at most about 0.1 h against a 1 h bound); their permutation p values are trustworthy at |log V| <= 0.4 and
  anti-conservative only at spreads narrower than any observed family.
<!-- src: paper/exploratory/roadmap_A1_offset_bound/output/{offsets_hours,offsets_bias_summary,offsets_size_summary,offsets_shift_summary}.csv; A258_spread_bias/output/v258_3_rows.csv -->

## A260 (scored 2026-11-15): the untimed procedure recovers clock dependence and feeding shifts of mouse liver genes (V260.1, V260.2, V260.3 PASS); a coupled liver shows no readable deviation of its own

- Data: Manella 2021, four arms x 24 mice, liver / WAT / lung; clock genes CIART and CLOCK dropped (A111 QC). Timed
  truth (non-clock liver genes): local 584, systemic 40. One scoring run (3.9 min) after three crashes before any
  statistic was computed (a host needs two anchors; CLOCK absent from the knockout lung; complex states).
- **V260.1 PASS:** in the knockout liver, the untimed host loading (host WAT + lung) separates systemic from local genes:
  AUC 0.997 (90 % 0.994 to 1.000).
- **V260.2 PASS:** untimed vs timed log amplitude ratio knockout / control, Spearman 0.870 (0.858 to 0.880; 1,545 genes);
  AUC local vs systemic 0.937 (0.901 to 0.966).
- **V260.3 PASS:** untimed vs timed DF-minus-AL shift relative to the lung, circular correlation 0.869 (0.856 to 0.883;
  572 genes; lung clock shift -2.31 h); knockout, reported: 0.886 (0.795 to 0.937; 24 genes).
- **D260.4 (descriptive):** kappa of the control liver against WAT + lung 0.001 (90 % -0.030 to 0.035; upper bound
  0.72 h); local genes with p < 0.01 on the liver's deviation 0.014 (other genes 0.019). Knockout liver kappa 0.368
  (0.077 to 0.774).
- **Reading (registered):** the untimed body reader is calibrated (G2 (ii) met for the body / local split and for the
  shift reader). The tissue-deviation term reads nothing where the local clock is coupled: in the control liver the
  genes that need the local clock do not load on the liver's deviation, because there is no deviation to load on.
  Local control is identified by removing the clock (knockout) or by moving times apart, not by a deviation regressor.
- Limits: truth and procedure on the same mice (optimistic); mouse labels assumed identical across tissue files;
  40 systemic genes. A descriptive column of the first output (median absolute shift difference) was mis-computed and
  removed from `tests.csv`; the code line is corrected, not rerun.
- Code: `exploratory/A260_mouse_calibration/scripts/score_a260.py`; host of two anchors passed to `rhythms.gene_rhythms`;
  V260.3 regresses every gene on the lung's own state (a host needs two compartments; registered "lung alone").
<!-- src: paper/exploratory/A260_mouse_calibration/output/{tests,descriptive,untimed_genes,truth_genes}.csv -->

## A261 (scored 2026-11-15; not blind for the clock genes): after a 12 h shift of behaviour, the rhythmic transcriptome of human muscle follows behaviour (shifted or evoked) while its clock moves 2.7 h (H261.1 PASS, H261.3 behaviour); single Hallmark programs are too weak to be read (H261.2: 1 of 28 gated)

- Power (scout): at 20 genes per program the three times are told apart with power >= 0.915 at level 0.01, truth
  excluded in at most 0.092 of simulated programs.
- **H261.1 PASS (positive control):** the 1,587 reference-rhythmic genes clear the rhythm gate (aligned loading 0.222
  vs matched q95 0.098) with an aligned argument of -1.41 h.
- **H261.3 (genome-wide program): behaviour.** Shift 11.65 h (90 % 10.54 to 12.76), excursion ratio 0.854 (0.498 to
  1.467); central excluded (p 0.001), tissue clock excluded (p 0.003), behaviour not (p 0.392). All 200 matched
  non-rhythmic sets of the same size are "unresolved" (none "behaviour").
- **H261.2:** 1 of 28 Hallmark programs passes the rhythm gate (TNF-alpha via NF-kB, 68 genes: aligned 0.504 vs q95
  0.500; ratio 0.035, called "tissue clock" because a flattened program is not excluded from r 0.233). The partition
  over Hallmark programs is not readable: their aligned loadings are 0.05 to 1.01 of the matched q95 (median 0.31).
  Over all 7,000 matched-set draws (35 programs x 200): unresolved 0.976, behaviour 0.017.
- **Descriptive, phase programs** (reference-rhythmic genes by 4 h bins of acrophase): gated 00-04 h (unresolved,
  10.56 h, ratio 0.59), 04-08 h (behaviour, -11.34 h, 0.56), 12-16 h (behaviour, -10.44 h, 0.69); not gated 08-12 h,
  16-20 h, 20-24 h (shifts -9.79, 11.90, -11.29 h).
- **Reading (registered; wording restricted 2026-11-15, see the masking limit):** the genome-wide daily program of
  muscle follows behaviour (shifted or evoked): most of the muscle's rhythmic transcriptome goes with sleep, activity and
  meals, while its core clock moved by only 2.68 h and damped (A257). It is not shown to be "shifted by 12 h". With
  A257, the muscle's clock and its rhythmic outputs come apart under misalignment. This is the reading Wefers 2018
  already gave for these arrays (clock genes kept their control pattern, fatty acid / PPAR genes changed; see
  `literature_positioning.md` section 11); A261 adds the program-level test against central time and the muscle clock
  in hours, with a matched companion, not a new direction.
- **Masking limit (added 2026-11-15, docs only, no new statistic):** with two clock times per condition, a reversed
  19:00 minus 07:00 difference does not distinguish a 12 h shift of an oscillation from an acute response to the
  fed / active versus fasted / rested state. The protocol (Wefers 2018, methods, read via Europe PMC) does not help:
  meals at 08:00, 12:30, 15:00, 20:00 and sleep 23:00-07:00 in the aligned condition, all "shifted by 12 h" in the
  misaligned one, so in both conditions one biopsy is at the behavioural morning (about 1 h after waking, overnight
  fasted) and the other at the behavioural evening (about 13 h awake, before dinner, fed during the day). The
  behavioural state at each biopsy therefore inverts exactly with the clock-time contrast; a fasting / waking response
  and a shifted rhythm predict the same data. The biopsy order relative to breakfast and the test day are not stated
  in the text (biopsies given as 08:00 / 20:00 in the methods and 07:00 / 19:00 in the figures). Only a design that
  samples a full cycle under a constant routine after the shift could separate the two.
- Limits: not blind for the clock genes; 11 men, two clock times, third day; light, sleep and meals move together; the
  reference cohort ate hourly while awake, so the rhythmic set leans towards sleep / activity-driven genes, and a
  meal-driven set may differ; the class is for the pooled program, not for each gene. The scout's power assumed
  reference amplitudes carry over to programs of 20 genes; on the arrays they do not (H261.2), so the scout
  overstated program-level power.
- Code: `exploratory/A261_nightshift_programs/scripts/score_a261.py` (one completed run, 1.9 min; two crashes in the
  matched-set draw of the genome-wide program after the Hallmark lines were printed; the rerun reproduces them, the
  draw being deterministic).
<!-- src: paper/exploratory/A261_nightshift_programs/output/{tests,programs,matched}.csv -->

## A259 (scored 2026-11-15): a tissue's clock departs from its person's body time by an upper bound of 1 to 3 h, and this departure is not shown to be clock (H259.1 FAIL, H259.2 not above companions in any tissue, H259.3 FAIL)

- Power first (registered addendum): at clock amplitude 1 and 200 donors the estimator is unbiased, detects 1.0 h with
  0.85 and, with no deviation, has a median upper bound of 0.70 h; a non-clock program shared by a donor's tissues
  raises kappa (by 0.13 at shared_sd 0.1). The addendum was committed 9 s after the scoring run was launched (a failed
  `git add` of an ignored output folder stopped the first commit); no output existed then (the output folder was empty).
- Data: GTEx fast deaths, 11 body tissues, 122 to 220 donors; MuTHER fat and skin, 680 women. One scoring run (472 s).
- **H259.1 FAIL:** median over the 11 targets of the upper 90 % bound of sigma 2.20 h; no target below 1 h. Per target
  kappa (sigma, 90 % interval in h): thyroid 0.073 (1.03; 0.65 to 1.32), visceral adipose 0.097 (1.19; 0.83 to 1.51),
  atrial appendage 0.131 (1.38; 0.96 to 1.82), tibial artery 0.155 (1.50; 1.18 to 1.81), tibial nerve 0.168 (1.56;
  1.16 to 1.91), subcutaneous adipose 0.253 (1.92; 1.63 to 2.20), skeletal muscle 0.289 (2.05; 1.68 to 2.45),
  oesophagus mucosa 0.320 (2.16; 1.71 to 2.57), skin not sun-exposed 0.529 (2.78; 2.34 to 3.26), skin sun-exposed
  0.544 (2.82; 2.42 to 3.27), lung 0.709 (3.22; 2.60 to 3.85).
- **H259.2: no tissue passes.** kappa > 0 in all 11 (bootstrap p 0.001 to 0.003, Holm 0.011), but none above the q95 of
  both companion families. Matched pseudo-clock sets alone reach a median kappa_tech of 0.049 (thyroid) to 0.355 (sun-
  exposed skin), 0.30 to 0.74 of the clock's kappa; their q95 is above the clock's kappa in 7 of 11 tissues (the clock
  is above it in tibial artery 0.155 vs 0.150, oesophagus mucosa 0.320 vs 0.292, non-sun-exposed skin 0.529 vs 0.528 and
  lung 0.709 vs 0.469); stress-responsive sets (median 0.39 to 2.0) have q95 above the clock's kappa in all 11.
  Cross-tissue companion share (matched, median) 0.017 to 0.036: shared non-clock variation is not what inflates kappa.
- **H259.3 FAIL:** MuTHER fat-skin sigma_pair 2.20 h (kappa 0.333, 90 % 0.161 to 0.557; upper 2.85 h against 1.41 h);
  matched kappa_tech q50 0.737 on these arrays, above the clock's kappa. GTEx, same pair (descriptive): subcutaneous
  adipose vs sun-exposed skin 2.29 h (upper 2.92), vs non-sun-exposed skin 2.80 h (3.51).
- **Reading (registered):** the bound in hours is not a headline (G3 component fails): the estimator bounds a tissue's
  person-specific departure from the body time at 1 to 3 h, and most or all of it is sample-level variation that
  non-clock genes of the same expression carry too. The "deviation" is not shown to be clock (H259.2), so no
  person-specific internal desynchrony between tissues is claimed, and none is excluded below about 1 to 3 h. For T2,
  a tissue-clock term cannot be told from sample-level noise in GTEx: the T2 tissue term is not registered; the GTEx
  part of the partition stays two-way (body vs day, A91), to which companions are to be added.
- Limits: sigma counts amplitude as phase; the matched companion's within-sample moment is the technical floor, not
  subtracted by registration; arrays (MuTHER) carry a larger technical floor than counts.
- Code: `exploratory/A259_tissue_deviation/scripts/score_a259.py`; `clocksync.deviation` (repository commits 40ded16,
  360abce).
<!-- src: paper/exploratory/A259_tissue_deviation/output/score/{tests,targets,living,companions}.csv; output/power/summary.csv -->

## A262 (scored 2026-11-15): after forced desynchrony, the genome-wide rhythmic program of blood neither stays on central time nor follows the displaced sleep: it damps to about a third and is displaced about 4 h (D1: H262.1 PASS, H262.2 "lost", at the level's edge); the night-shift PBMC design is uninformative (D2: gate not passed)

- Reference (GSE39445, sleep extension, extended wakefulness): 8,234 reference-rhythmic non-clock genes of 15,619; clock
  time = circadian phase + 4.11 h.
- **D1 GSE48113 (21 subjects, 8,234 genes):**
  - **H262.1 PASS:** aligned loading 1.105 above the matched q95 0.992; aligned argument -0.47 h.
  - **H262.2: "lost".** Shift 3.96 h (90 % 2.23 to 5.69), excursion ratio 0.355 (0.203 to 0.621). Central excluded
    (p < 0.0001), behaviour (12 h) excluded (p < 0.0001), lost not excluded (p 0.0101, against the level 0.01).
  - All 200 matched non-rhythmic sets are "mixed" (none "lost"). Their aligned loadings reach nearly the program's
    (q95 0.992 vs 1.105): genes that are flat under the reference's constant wakefulness carry clock-time structure in
    D1, which contains sleep. The gate passes, but narrowly.
  - Phase programs (descriptive): 00-04 h lost (shift 4.58 h, ratio 0.248, gated), 12-16 h mixed (3.61 h, 0.417,
    gated); the four smaller bins are not gated (shift 2.69 to 6.75 h, ratio 0.365 to 0.594, all mixed).
- **D2 GSE122725 dim-light arm (8 subjects, 7,074 genes): uninformative.** H262.1 FAIL: aligned loading 0.245 below the
  matched q95 0.254, even though the aligned argument is 0.33 h. Descriptive only: shift 0.72 h (-1.88 to 3.32), ratio
  0.445 (0.206 to 0.963); with central at 1.5 h, p central 0.034, behaviour 0.015 and lost 0.280 ("unresolved"). The
  central sensitivities at 0 h and 2 h give the same class.
- **Reading (registered):** in D1 the blood program does not move with the displaced sleep and meals (12 h is
  excluded); it loses about two thirds of its excursion and lies about 4 h from its day-1 timing. This is Archer
  2014's loss of rhythmicity (6.4 % -> 1.0 %) read at the program level. The registered condition for the
  tissue-specific sentence (both blood designs central, or one central and one lost) is not met, because D2 is
  uninformative. Muscle (A261: behaviour, ratio 0.854) and blood (D1: not behaviour, ratio 0.355) differ, but the
  sentence "which time a tissue's outputs follow is tissue-specific" is not written.
- Limits:
  - The "lost" call sits at the level's edge (p 0.0101).
  - The matched sets are almost as strong as the program in D1.
  - D1's sampling day contains the displaced sleep, so sleep-locked cell-composition changes can enter. Blood rhythms
    include diurnal composition; the program is not a cellular clock.
  - Subjects may overlap between the reference and D1.
  - Not blind for the GSE48113 clock and cue genes (A156: core state -1.36 h, cue state -6.85 h).
- Process:
  - Run 1 returned NaN for every D1 program, because GSE48113 arrays have missing values. Missing cells (156,920,
    across 15,619 genes x the arrays of 21 subjects) were then set to the gene's mean within the subject x condition
    unit, the A156 rule. This was decided without any D1 number, and is recorded in the script.
  - D2 is identical in both runs (its draws do not depend on the values). Run 1's tables are kept as `output/run1_*`.
- Code: `exploratory/A262_blood_misalignment/scripts/score_a262.py` (3 min).
<!-- src: paper/exploratory/A262_blood_misalignment/output/{tests,programs,matched}.csv, score.log -->

## A264 (scored 2026-11-15; NOT BLIND): the fetal adrenal cortex covariation survives further age and cohort adjustment and the steroidogenic activity score (H264.1, H264.2 PASS), but without PER1 and NR1D1 it is not clock-specific (H264.3): by the registered reading it rests on the two signal-induced clock genes and moves to the supplement
- Power first (`feasibility_a264.py`, 200 adult subsamples): every readout has an adult pass rate of 1.000 at the full
  clock (0.990 to 1.000 at 0.75; 0.365 to 0.545 at 0.5; false-positive rates 0.025 to 0.080). The adult steroid score
  shares almost no variance with the hour of death (R^2 0.004; 20 genes, mean pairwise r 0.485), so the adult power does
  not test separability from steroidogenic activity (limit registered).
- `cortex` (A244 design reproduced): C 0.608, p 0.0042; matched q95 0.367, responsive q95 0.482: PASS (A244: 0.608,
  0.0038; the companion draws use this amendment's seed).
- **H264.1 PASS** `cortex:age` (cubic age, cohort slopes for kameneva, kildisiute, sasaki, wang; 13 parameters): C 0.659,
  p 0.0021; matched q95 0.342, responsive q95 0.465.
- **H264.2 PASS** `cortex:steroid` (Xing 2010 ACTH score): C 0.526, p 0.0121; matched q95 0.319, responsive q95 0.462.
- **H264.3 not clock-specific** `cortex:-PER1-NR1D1` (12 genes): C 0.489, p 0.0290; above the matched q95 0.397 (98th
  percentile) but below the responsive q95 0.634 (76th percentile).
- Descriptive: `cortex:all3` C 0.438, p 0.0401, matched q95 0.417, responsive q95 0.565 (not clock-specific);
  `cortex:acute` (NR4A1-3 score) C 0.546, p 0.0066, clock-specific. Leave-one-gene-out with 14 genes: 0.438 (without
  ARNTL) to 0.808; unsigned R of `cortex:steroid` 0.482 (A244 design 0.630).
- **Reading as registered:** the covariation is not a smoother age trend, a cohort-specific age slope or the ACTH-driven
  steroidogenic program. Without PER1 and NR1D1 it remains above expression-matched sets and nominally significant, but
  it is not above the glucocorticoid / immediate-early sets: it is not shown to be a clock-network state, and the fetal
  result moves to the supplement as a hypothesis (the fetal exception counts as **fallen** in the tally of the lead
  prompt).
- **Limits.** (1) Not blind. (2) The responsive family draws from a pool of about 30 genes; with 12 genes its q95 rises
  from 0.465-0.482 (14 genes) to 0.634, so this verdict turns on the noisier companion; the adult power of the
  clock-specificity condition was not computed (the feasibility tests p only). (3) If the fetal clock drives
  steroidogenesis, the steroid adjustment removes clock signal too (registered limit 4); H264.2 passed nonetheless.
  (4) A static axis of another kind is not excluded by any of the three.
- Code: `exploratory/A264_fetal_adrenal_covariates/scripts/{design_a264,feasibility_a264,score_a264}.py` (278 s and
  21 s).
<!-- src: paper/exploratory/A264_fetal_adrenal_covariates/output/{verdicts,l1,feasibility_passrates,adult_scores}.csv -->

## A266 (scored 2026-11-15; NOT BLIND): without the signal-induced clock genes (PER1, PER2, NR1D1) the fetal adrenal cortex covariation is no longer significant (H266.1 fail, at an adult pass rate of 1.000); the responsive-gene relation is not estimable (H266.2); the fetal cortex expresses the ACTH receptor at a high level and no melatonin receptor (H266.3, descriptive): a hormone-driven response is the simpler reading
- Power first: `cortex:core` (11 genes) adult pass rate 1.000 / 0.995 / 0.630 at 1 / 0.75 / 0.5, false-positive rate
  0.055. H266.2: 0.500 at 34 adult glands (false-positive rate 0.075): not estimable, fetal statistic not computed.
- **H266.1 fail** `cortex:core`: C 0.417, p 0.0572; below the matched q95 0.493 (87th percentile) and the responsive
  q95 0.522 (74th). Descriptive: `cortex:core:steroid` C 0.311, p 0.106; `medulla:core` C -0.339, p 0.706.
- **H266.3 (descriptive):** median log2 cpm in the fetal cortex: MC2R 7.44 and MRAP 8.44 (detected in 34 of 34),
  NR3C1 3.19, HSD11B2 3.26; MTNR1A and MTNR1B detected in 0 of 29 cortex pseudobulks (symbol absent in dong and han).
  Medulla: MC2R 2.88, NR3C1 5.44, MTNR1A in 3 of 24. Adult adrenal (GTEx, 229): MC2R 6.38, MRAP 3.06, NR3C1 4.77;
  MTNR1A/B absent from the matrix. The receptor count table was rerun once with NaN-aware medians and detection counts
  after the first run returned NaN for MTNR1A/B (source symbols missing); the clock statistics are identical in both runs.
- **Reading as registered:** with PER1, PER2 and NR1D1 removed, the covariation is not significant; a core loop of
  three quarters of the adult size is excluded (adult pass rate 0.995 at 0.75). The covariation rests on the
  signal-induced clock genes: a hormone-driven response is the simpler reading; an entrained loop smaller than three
  quarters of the adult one is not excluded (0.630 at 0.5). The fetal cortex carries the ACTH receptor complex at
  adult-like or higher levels and no detectable melatonin receptor in these 3' single-cell data, which favours fetal
  ACTH / cAMP (and its own cortisol) over a direct maternal melatonin input as the driver; maternal cues acting through
  the fetal pituitary are not excluded.
- **Limits.** Not blind; signal-induced genes chosen from the literature (from memory); no sampling times; melatonin
  receptors are low-abundance and may be missed by droplet 3' assays; limits of A244 and A264.
- Code: `exploratory/A266_fetal_adrenal_drive/scripts/{design_a266,feasibility_a266,score_a266}.py` (358 s; power
  stage about 3 min; score 43 s).
<!-- src: paper/exploratory/A266_fetal_adrenal_drive/output/{verdicts,l1,feasibility_passrates,power_h266_2,receptors}.csv -->

## A243 (scored 2026-11-16): when open cohorts are pooled, the fetal lung and liver carry a clock-shaped covariation between fetuses (H243.1 PASS in two of four organs) of between half and three quarters of the adult size, as the adrenal; heart (powered at half the adult size) and kidney are negative
- Design (addendum 42fc9ad, before scoring): heart 72 fetuses in 8 cohorts, kidney 53 in 7, lung 38 in 5, liver 35 in
  3; every registered readout had an adult pass rate >= 0.8 at the full clock (lowest `lung:early` 0.950), so all were
  scored. At half the adult size only `heart:all` reached 0.8 (0.800); `all` readouts at 0.75: >= 0.995.
- **H243.1** (Holm over 4 organs):

  | readout | fetuses | C | p | Holm p | matched q95 (pct) | responsive q95 (pct) | verdict |
  |---|---|---|---|---|---|---|---|
  | `lung:all` | 38 | 0.530 | 0.0087 | 0.026 | 0.396 (98.0) | 0.429 (99.5) | **PASS** |
  | `liver:all` | 35 | 0.583 | 0.0062 | 0.025 | 0.460 (97.0) | 0.398 (100) | **PASS** |
  | `heart:all` | 72 | 0.032 | 0.329 | 0.658 | 0.360 (66.5) | 0.269 (73.0) | fail |
  | `kidney:all` | 53 | 0.017 | 0.439 | 0.658 | 0.413 (62.5) | 0.246 (64.0) | fail |

- **H243.2** (Holm over 8 windows): `lung:early` (13 fetuses, < 12 weeks) C 0.679, p 0.0022, Holm 0.018, above both
  companions: **PASS**. Lost at Holm (suggestive, not established): `lung:late` C 0.398, p 0.040; `liver:late` 0.497,
  p 0.024. Fail: `liver:early` 0.169 (p 0.235), `heart:late` 0.278 (0.133), `heart:early` -0.262, `kidney:late` 0.365
  (0.103), `kidney:early` -0.412.
- **Descriptive.** Without Cao 2020 and Cardoso-Moreira (`all:new`): lung 0.574 (p 0.0045), heart 0.142, kidney
  -0.145 (liver has no such readout). Leave one cohort out: lung stays significant without each cohort (lowest without
  Sountoulidis: 0.437, p 0.039); liver without Cardoso-Moreira 0.249 (p 0.19, 9 fetuses), without Cao 0.521 (p 0.015).
  Cohorts alone: lung Sountoulidis 0.545 (p 0.011, 17 fetuses, above both companions), lung Cao 0.190 (p 0.19), liver
  Cardoso-Moreira 0.240 (p 0.14), heart Cardoso-Moreira -0.030, heart Cui 0.024, kidney Cardoso-Moreira -0.489.
  Leave one gene out (minimum): lung 0.439 (without NPAS2), liver 0.499 (CRY1). Unsigned R: lung 0.575, liver 0.429.
  Position among the adult subsamples (share of adult C at or below the fetal C at 1 / 0.75 / 0.5): lung 0.0 / 8.5 /
  89.0 %, liver 0.5 / 46.0 / 100 %, `lung:early` 37.5 / 76.5 / 99.0 %; heart 0 / 0 / 0 %, kidney 0 / 0 / 3.0 %
  (adrenal, A242 `all`: 1.0 / 42.0 / 100 %).
- **Reading as registered:** a clock-shaped covariation between fetuses is present before birth in the lung and the
  liver as in the adrenal, between half and three quarters of each organ's adult size; the developmental result becomes
  "present in several organs, below the adult size" and the earlier negatives (A189 Cao lung 0.116, liver 0.114) are
  read as bounds on size. Heart: negative with an adult pass rate of 0.800 at half the adult size: a clock of half the
  adult size is excluded and the heart is behind the adrenal for the ages covered (4-28 weeks). Kidney: three quarters
  excluded (0.995), the comparison with the adrenal stays open (0.250 at 0.5). Lung is strongest before 12 weeks.
- **Post-scoring (not registered here):** in the adrenal the covariation did not survive the removal of the
  signal-induced clock genes (A266). C measures a clock-shaped covariation, not a 24 h period; the same two checks
  (core loop without PER1, PER2, NR1D1; glucocorticoid-target covariate) are registered for lung and liver as A269.
- **Limits.** Liver is largely haematopoietic at these ages and rests on pooling (no cohort alone is significant); lung
  rests in part on the deepest cohort (Sountoulidis, 5-14 weeks); unrecorded times; age conventions differ between
  cohorts; fetus identity across UK HDBR sources not testable; limits of A242.
- Code: `exploratory/A243_fetal_organs_pooled/scripts/{prep_a243,feasibility_a243,companions_dryrun_a243,score_a243,position_a243}.py`
  (feasibility 78 min; score 460 s).
<!-- src: paper/exploratory/A243_fetal_organs_pooled/output/{verdicts,l1,feasibility_passrates,position_in_adult}.csv -->

## A265 (scored 2026-11-15/16): COFE's published "uncoupled" reading is not reproduced on the TCGA pairs (uniform differences in 3 of 10 cancers), COFE's own labels pair tumour with normal (H265.1 PASS), and in a simulation on the real design COFE's published reading says "uncoupled" in 7 of 10 coupled replicates while clocksync detects coupling in 10 of 10 and never without it (H265.2: sampling-window reading)
- Data: 464 tumour-normal pairs in ten cancers (BLCA 19, BRCA 91, COAD 39, KIRC 72, KIRP 32, LIHC 49, LUAD 43, PRAD 51,
  THCA 49, UCEC 19); COFE genome-wide on each cancer's tumours (1,138-1,321 genes), normals projected (assumption).
- **(i) COFE's readout:** tumour-minus-normal label differences uniform (Rayleigh p >= 0.05) in 3 of 10 cancers (BRCA
  diff R 0.085, COAD 0.054, UCEC 0.174); concentrated in BLCA 0.509, KIRC 0.386, KIRP 0.425, LIHC 0.379, LUAD 0.369,
  PRAD 0.571, THCA 0.257 (p 0.039 to < 0.0001). Tumour labels cover the cycle (largest empty arc 0.23-0.67 h); CV minimum
  interior in all ten. **Not reproduced** (criterion: uniform in >= 6 of 10).
- **H265.1 PASS:** pooled R_c 0.302, null mean 0.245, p 0.0013 (10,000 permutations). clocksync (A138 statistic) on the
  same pairs: 0.361, null mean 0.129, p 0.0001, all ten cancers gated (descriptive).
- **H265.2 (simulation on the real design, 10 replicates per arm, real pairing broken, 365-gene GTEx program):**

  | arm | clocksync pass | COFE pairing pass | COFE published reading "uncoupled" | median uniform cancers | clocksync stat (median) | COFE R_c (median) |
  |---|---|---|---|---|---|---|
  | coupled, 08-18 h | 10 / 10 | 1 / 10 | 7 / 10 | 6.0 | 0.223 | 0.224 |
  | coupled, whole day | 10 / 10 | 1 / 10 | 5 / 10 | 5.5 | 0.334 | 0.240 |
  | uncoupled, 08-18 h | 0 / 10 | 0 / 10 | 5 / 10 | 5.5 | 0.132 | 0.238 |

- **Reading as registered:** the first reading holds (COFE "uncoupled" in >= 5 of 10 coupled surgical replicates: 7;
  clocksync >= 8 of 10: 10; clocksync <= 2 of 10 uncoupled: 0): COFE's published negative is not evidence against
  coupling under surgical hours, and A138 stands. Real data: H265.1 PASS, COFE's own labels pair tumour with normal,
  contradicting the published reading with COFE's output.
- **Observed, beyond the readings:** COFE's published readout is uninformative even with whole-day sampling (5 of 10
  coupled replicates "uncoupled"; the same rate, 5 of 10, without coupling), and COFE's pairing statistic detects the
  simulated coupling in 1 of 10 replicates in both coupled arms. Post hoc: the real-data COFE pairing (p 0.0013) is
  therefore unlikely to come from the clock coupling the simulation inserts; A267 tests what carries it.
- **Limits.** Projection of the normals is an assumption; GTEx-derived program (attenuated amplitudes); 10 replicates
  per arm; COFE transform deviates from DESeq2; 464 of A138's 517 pairs.
- Code: `exploratory/A265_cofe_head_to_head/scripts/{lib_a265,scout_a265,score_a265,simulate_a265}.py` (real part 45
  min; simulation 6.8 h at 12 workers).
<!-- src: paper/exploratory/A265_cofe_head_to_head/output/{real_tests,real_per_cancer,sim_summary,sim_reps}.csv -->

## A267 (scored 2026-11-16): COFE's tumour-normal pairing is not clock-specific: on the clock genes it is matched by expression-matched gene sets (H267.1 not clock-specific), and genome-wide it is no larger than generic expression axes (H267.2 FAIL); it is patient similarity, the clock evidence remains clocksync's A138
- **H267.1, not clock-specific:** COFE on the 14 clock genes, pooled R_c 0.290, null mean 0.188, p 0.0001; matched
  14-gene sets median 0.265, q95 0.321 (clock at the 77th percentile, 100 sets); responsive median 0.234, q95 0.254
  (100th percentile). Per cancer R_c 0.175 (KIRP) to 0.445 (KIRC). Agreement of clock-gene COFE tumour phases with the
  clocksync tumour state phase: 0.272 (BRCA) to 0.749 (THCA).
- **H267.2 FAIL:** the genome-wide COFE pairing of A265 (0.302) lies at the 82nd percentile of the 45 generic axes
  (pairs of the top 10 principal components; median 0.206, q95 0.340, maximum 0.361).
- Descriptive: genome-wide COFE without the 365 clock and GTEx-rhythmic genes of the A265 program: 0.250, null mean
  0.220, p 0.055.
- **Reading as registered:** COFE's clock-gene pairing is matched by non-clock sets (shared patient effects or a
  generic axis); for the clock, the evidence is clocksync's A138 (clock-specific against both families). The
  genome-wide COFE pairing of A265 is patient similarity, not evidence of clock coupling: A265 H265.1 is read as "COFE's
  genome-wide labels carry patient-shared structure". With A265, neither COFE's published negative nor its own positive
  pairing bears on tumour clock coupling; the head-to-head favours clocksync on these data.
- **Limits.** Projection of the normals is an assumption; short sparsity grid on 14 features; 100 sets per family; the
  responsive pool (about 29 genes) gives overlapping sets; principal-component axes are a generic companion.
- Code: `exploratory/A265_cofe_head_to_head/scripts/score_a267.py` (6.0 h, 4 workers).
<!-- src: paper/exploratory/A265_cofe_head_to_head/output/{a267_tests,a267_per_cancer,a267_sets,a267_axes}.csv -->

## A268 (scored 2026-11-16): glucocorticoid, cold and heat readers each keep a consistent phase relation to the body's clock state across GTEx tissues, but none beyond matched and stress-responsive gene sets (all six tests "not specific"), and none at the GTEx phase in living MuTHER tissue: these readers do not name the systemic signal
- **Tests** (12 body tissues, K = |sum beta| / sum |beta|, 10,000 permutations; Holm over three readers per host):

  | reader | host | K | p | Holm p | matched q95 (pct) | responsive q95 (pct) | MuTHER p | MuTHER vs GTEx phase (h) | verdict |
  |---|---|---|---|---|---|---|---|---|---|
  | glucocorticoid (FKBP5 + TSC22D3) | all | 0.703 | 0.0048 | 0.0096 | 0.812 (79.5) | 0.903 (65.0) | 0.0001 | 4.09 | not specific |
  | glucocorticoid | without PER1 | 0.716 | 0.0041 | 0.0082 | 0.819 (85.0) | 0.916 (63.0) | 0.0001 | 4.71 | not specific |
  | cold (CIRBP) | all | 0.786 | 0.0010 | 0.0030 | 0.830 (88.0) | 0.907 (70.0) | 0.041 | -6.14 | not specific |
  | cold | without PER1 | 0.782 | 0.0005 | 0.0015 | 0.829 (90.0) | 0.913 (59.5) | 0.045 | -6.06 | not specific |
  | heat (HSP90AA1) | all | 0.696 | 0.0053 | 0.0096 | 0.873 (76.0) | 0.907 (59.5) | 0.0001 | 7.19 | not specific |
  | heat | without PER1 | 0.689 | 0.0050 | 0.0082 | 0.875 (75.5) | 0.919 (50.5) | 0.0001 | 7.25 | not specific |

- MuTHER: the pooled fat + skin loading is significant for every reader (p 0.0001 to 0.045) but its phase lies 4.1 to
  7.3 h from the GTEx body phase (registered tolerance +-3 h): no reader is confirmed in living tissue.
- Brain (descriptive): projection of the cerebellum and cortex loadings relative to the median body tissue 0.96 to 2.09,
  every bootstrap interval includes 1 (e.g. glucocorticoid cerebellum 1.90 [-0.05, 4.85], cold cortex 0.96 [-0.16,
  3.06]): looser coupling of these readers in brain is not shown.
- **Reading as registered:** "not specific": the relations are shared by expression-matched and stress-responsive
  gene sets (matched median K 0.49-0.57 over the readers, responsive 0.59-0.65): a general state that many genes share with the body's clock
  state, not a named signal. These public data do not name the systemic signal by glucocorticoid, cold or heat
  readers. With MuTHER, the living-tissue relation is significant but at a different phase, so the GTEx relation is
  not transferable either.
- **Limits.** Relaxed reader rule (user decision); agonal stress in GTEx; MuTHER has two tissues and no times; readers
  of one or two genes; matched sets drawn genome-wide carry part of the clock-state signal themselves (A259: matched sets
  carry 0.30 to 0.74 of a tissue's apparent deviation), which makes "not specific" a conservative verdict.
- Code: `exploratory/A268_systemic_signal/scripts/score_a268.py` (power stage 51 min; score 12 min).
<!-- src: paper/exploratory/A268_systemic_signal/output/{tests,companions,muther,brain,power_muther}.csv -->

## A269 (scored 2026-11-16; NOT BLIND): the fetal lung and liver covariation survives adjustment for glucocorticoid-target activity (H269.2 PASS in both), but on the core loop without PER1, PER2 and NR1D1 neither organ passes: lung is just above p 0.05 (fail), liver is significant but not above the responsive companion (not clock-specific); both keep more of their signal on the core loop than the adrenal cortex did (A266)
- Power first (registration d8b7df2): adult pass rate at the full clock >= 0.975 for every readout; core loop at half
  the adult size: lung 0.840, liver 0.910. Scoring: two runs were killed for memory before any statistic (pyarrow with
  OMP_NUM_THREADS=1 peaked at 32 GB when reading the 88,884-gene pseudobulk); fixed by reading before the thread limit
  (commit after d8b7df2), no statistic changed or seen.
- **H269.1** (core loop, 11 genes; Holm over 2):

  | readout | fetuses | C | p | Holm p | matched q95 (pct) | responsive q95 (pct) | verdict |
  |---|---|---|---|---|---|---|---|
  | `lung:core` | 38 | 0.473 | 0.060 | 0.060 | 0.514 (92.5) | 0.206 (100) | fail |
  | `liver:core` | 35 | 0.686 | 0.0086 | 0.017 | 0.516 (99.5) | 0.712 (94.0) | not clock-specific |
  | adrenal cortex (A266) | 34 | 0.417 | 0.057 | - | 0.493 (87) | 0.522 (74) | fail |

- **H269.2** (14 genes plus the GC score; Holm over 2): `lung:gc` C 0.572, p 0.0034, Holm 0.0068, above both
  companions (99.0, 100): **PASS**; `liver:gc` C 0.589, p 0.0065, Holm 0.0068 (100, 100): **PASS**.
- Descriptive: `lung:core:gc` 0.450 (p 0.078); `liver:core:gc` 0.700 (p 0.007; matched 98.5, responsive 95.5);
  `lung:early:core` (13 fetuses) 0.647, p 0.014, above both companions (98.5, 99.5). Adult GC score: R^2 on the hour of
  death 0.004 (lung) and 0.015 (liver); mean pairwise r of the four GC genes 0.60 and 0.58.
- **Reading as registered.** H269.2: the covariation of lung and liver is not explained by glucocorticoid-target
  activity (as the adrenal was not by steroidogenic activity, A264). H269.1, lung: fail; by the registered reading a
  core loop of half the adult size is excluded and the lung reads as the adrenal. H269.1, liver: "driven or entrained,
  not separable" (above the matched, below the responsive q95; the responsive q95 on 11 genes is 0.712, as anticipated
  in limit 5). Neither organ carries a core loop shown to be clock-specific; the developmental result narrows to "the
  clock genes covary before birth in the adrenal, lung and liver, beyond glucocorticoid activity; a core loop is not
  shown in any".
- **Observed, beyond the readings (post hoc):** on the core loop lung keeps 89 % of its 14-gene C (0.473 of 0.530) and
  the liver more than all of it (0.686 of 0.583), whereas the adrenal cortex kept 69 % (0.417 of the A264 cortex 0.608);
  the lung fail is at p 0.060 and above the responsive companion. The lung negative is borderline rather than a collapse,
  and the early lung (before 12 weeks) is significant on the core loop above both companions (descriptive, 13 fetuses).
- **Limits.** Not blind; CORE and GC from the literature; GC score is four genes and a readout of glucocorticoid
  action, not of hormone level; other inputs not adjusted; liver haematopoietic and pooled; limits of A243.
- Code: `exploratory/A269_fetal_organs_drive/scripts/{design_a269,feasibility_a269,score_a269}.py` (power 5 min; score
  22 s).
<!-- src: paper/exploratory/A269_fetal_organs_drive/output/{verdicts,l1,feasibility_passrates,adult_scores}.csv -->

## A263 Part A, human, baboon and sheep (scored 2026-11-16; NOT BLIND): brain's lower coupling to the body's clock state is not specific to the clock genes: expression-matched and stress-responsive gene sets show the same or a larger drop (H263.1, H263.2, H263.3 fail); the brain exception falls as a clock statement (erratum to A142)
- **R2, L = log(median brain coupling / median peripheral coupling)** (peripheral-anchor host; 200 matched and 200
  responsive sets; registered rule: clock-specific if L_clock < q05 of both families; sets with undefined L, i.e.
  median brain coupling 0, are counted and left out, which favours the clock):

  | species | brain / periphery tissues | median coupling brain vs periphery | L clock | matched: valid, q05, median (pct at or below clock) | responsive: valid, q05, median (pct) | verdict |
  |---|---|---|---|---|---|---|
  | human (GTEx v10 fast deaths) | 13 / 26 | 0.406 vs 0.679 | -0.514 | 174, -3.207, -1.077 (75.3) | 199, -3.209, -1.306 (87.9) | **H263.1 fail** |
  | baboon | 22 / 35 | 0.420 vs 0.712 | -0.528 | 172, -2.509, -0.706 (59.3) | 175, -2.522, -0.770 (62.3) | **H263.3 fail** |
  | sheep (mixed strata) | 8 / 35 | 0.205 vs 0.416 | -0.707 | 182, -3.007, -0.287 (32.4) | 200, -1.552, -0.394 (30.0) | **H263.2 fail** |

- **R1 (human, untimed block agreement r* brain vs periphery):** clock 0.803; matched median 0.417 (q05 0.091),
  responsive median 0.376 (q05 0.131): the clock genes agree between brain and periphery more than non-clock sets do,
  not less: not clock-specific as a drop. **Sheep:** clock 0.614; matched median 0.445 (q05 0.087), responsive median
  0.711 (q05 0.339): above the matched sets, below the responsive median.
- **Reading as registered:** H263.1 fails (L_clock is not below the q05 of either family): brain's lower coupling is
  shared by non-clock gene sets, for the responsive family "brain responds less to the body's shared stress /
  hormonal state"; **the brain exception falls as a clock statement**, and A142 gets the erratum below. Baboon and
  sheep fail too. In sheep the clock drop is larger than most companions' (L -0.71 against medians -0.29 and -0.39), but
  not beyond the q05 of either family.
- **Observed, beyond the readings:** the non-clock sets show a larger brain drop than the clock genes (median L -1.08
  matched and -1.31 responsive in human against -0.51; 26 matched sets with zero brain coupling), and the clock genes
  keep the higher brain-periphery agreement (r* 0.803 against medians 0.42 and 0.38). Read together: brain shares less
  of every between-person state with the body than peripheral tissues do, and of the gene sets tested the clock genes
  are among the most shared, not the least.
- **Limits.** Not blind; R2 is new in human (A142 used timed departures) and R1 untimed; the responsive pool holds
  28-29 genes (overlapping sets); baboon genes restricted to one-to-one orthologs; pseudo sets carry the GTEx waveform.
- Runs: baboon 49 min (6 workers; two earlier runs lost at the write step, no result seen, output folder fix);
  human 4.4 h (8 workers; the runner left no ledger line, the outputs are complete: 401 sets, `tests_human.csv`);
  sheep 13.8 h (3 workers; no ledger line either, outputs complete: 401 sets, no errors, `tests_sheep.csv`).
- Code: `exploratory/A263_companions/scripts/{build_genomewide,score_a263,score_a263_human,score_a263_baboon}.py`.
<!-- src: paper/exploratory/A263_companions/output/{tests_human,sets_human,tests_baboon,sets_baboon,tests_sheep,sets_sheep}.csv -->

### Erratum to A142 (2026-11-16, from A263 H263.1)
A142 H142.1 (brain block r* 0.639 against random 13-tissue blocks, PASS) stands as registered: brain is less coupled
to the body's shared state than other organ blocks. It is **not shown to be clock-specific**: on GTEx fast deaths the
drop of brain coupling is matched by expression-matched and stress-responsive gene sets (A263). Wording from now on:
"brain is less coupled to the body's shared state", not "brain's clock is the exception".

## A263 Part B (scored 2026-11-16; NOT BLIND): the atlas couplings of Figure 2 are not clock-specific against the responsive family overall (H263.4 FAIL); fibroblast / stromal and "other" clear both families
Run `compute_a150.py --kind responsive` (4 workers, 12 min, peak 1.7 GB); fits reproduce A17 (max |dc| 0.009 lung_hlca,
~0 elsewhere); 200 of 200 sets valid, 0 failed fits.
- **H263.4 FAIL:** pooled mean coupling of the 53 compartments 0.393; responsive pooled sets median 0.321, q95 0.423
  (clock at the 87.5th percentile), p 0.129. Against matched sets (A150 H150.1) the same coupling passes (q95 0.313).
- **Per family (descriptive), clock median vs median / q95 of each family:**

  | family (compartments) | clock | matched median / q95 | responsive median / q95 |
  |---|---|---|---|
  | other (2) | 0.754 | 0.454 / 0.729 | 0.455 / 0.648 |
  | fibroblast / stromal (4) | 0.670 | 0.302 / 0.484 | 0.368 / 0.546 |
  | mural (4) | 0.613 | 0.401 / 0.705 | 0.406 / 0.584 |
  | endothelial (8) | 0.520 | 0.187 / 0.345 | 0.371 / 0.557 |
  | epithelial (7) | 0.517 | 0.265 / 0.430 | 0.348 / 0.554 |
  | parenchyma (1) | 0.461 | 0.415 / 0.858 | 0.420 / 0.675 |
  | neuron (1) | 0.449 | 0.272 / 0.703 | 0.475 / 0.751 |
  | myeloid (10) | 0.328 | 0.148 / 0.282 | 0.278 / 0.398 |
  | lymphocyte (16) | 0.069 | 0.071 / 0.217 | 0.220 / 0.372 |

  Above both q95: other, fibroblast / stromal. Responsive q95 only: mural. Matched q95 only: endothelial, epithelial,
  myeloid. Neither: parenchyma, neuron, lymphocyte (lymphocytes at the median of both families; neuron below the
  responsive median).
- **Reading as registered:** the organ-map coupling claims of Figure 2 are restated as "shared state", not clock.
- **Post hoc reading (2026-11-16, user and agent, after A270; replaces the registered wording for the paper):**
  (1) Absolute couplings are not "how clock-coupled a cell type is": any expression-matched gene set reaches a
  cell-type-dependent baseline (matched medians 0.07 lymphocytes to 0.45 "other"; causes untested), so couplings are
  compared across cell types only as the excess over the cell type's own companions. (2) The excess is clear for
  fibroblast / stromal cells (both families), present against matched genes for endothelial, epithelial and myeloid
  cells and borderline against responsive genes, absent for lymphocytes. (3) The responsive family is a harsh control
  (glucocorticoid targets and immediate-early genes); in GTEx deaths its shared state does not read time of death (A270
  H270.2), in the living atlases what it reads is untested. (4) Clock identity is shown by phase-based readouts (organ
  offsets vs timed truth: A36 r_FL 0.91, A123; host phase vs time: A270 clock above all 400 sets in both species;
  r* A263), at the organ level; no timed truth exists for atlas cell types.
<!-- src: paper/exploratory/A150_atlas_specificity/output/responsive/{h150_1,per_family,per_compartment,pseudo_pooled,reproduction}.csv; A150 matched: output/{h150_1,per_family}.csv -->

## A270 (scored 2026-11-16; NOT BLIND on the clock): the responsive sets read time in baboon (H270.1 PASS) but so do most matched sets; in GTEx they do not read time of death (H270.2 FAIL); the clock host reads time above all 400 sets in both
Runs: baboon 3.7 min, human 5.3 min (6 workers each, peak 0.6 / 0.8 GB); 0 errors.

| | baboon (ZT; 12 animals, 63 tissues) | human (GTEx fast deaths, time of death; 244 donors, 15 tissues) |
|---|---|---|
| clock: R, p_R, r_FL | 0.987, 0.001, 0.980 | 0.631, 0.001, 0.388 |
| responsive: share p_R < 0.05 | **0.555 (H270.1 PASS)** | **0.170 (H270.2 FAIL)** |
| matched: share p_R < 0.05 | 0.735 | 0.100 |
| median R responsive / matched | 0.537 / 0.622 | 0.073 / 0.054 |
| q95 R responsive / matched | 0.856 / 0.754 | 0.149 / 0.137 |
| median max(R, R reflected) resp. / matched | 0.681 / 0.644 | 0.082 / 0.069 |
| clock R above every set | yes (100 / 100 pct) | yes (100 / 100 pct) |

- **Reading (registered, split by species):** baboon: responsive sets carry time, but no more than matched sets (12
  animals around the clock: most genes read ZT, as Mure 2018); this does not single out the responsive family. GTEx:
  the responsive sets do not read time of death; their shared state is not internal time there, and the failures
  against them (A263, A269, scored on human data) stand as warnings that a non-time shared state reaches the clock-gene
  coupling. The clock host reads recorded time far above every companion in both species: specificity of the clock is
  shown on the timed readout, not on coupling size.
- **Limits:** not blind on the clock; pseudo sets carry the clock waveform's phases; GTEx time of death is recorded wall
  time; human tested in deaths only; the 3-set `--check` run on baboon was seen before scoring.
<!-- src: paper/exploratory/A270_responsive_time/output/{tests_baboon,sets_baboon,tests_human,sets_human}.csv -->

## A271 Part A (scored 2026-11-16; NOT BLIND on the matched family): the tumour-normal immune pairings are not clock-specific against the responsive family (H271.1 FAIL for myeloid, T/NK and CD8-negative T); fibroblasts and endothelium clear it
Runs: A227 tests 24 min (2 workers, peak 28.9 GB), A228 tests 21 min (3 workers, peak 33.2 GB); registered statistics
reproduce exactly (max |difference| 0.0).

| test | patients | pairing | responsive median / q95 | p_resp | Holm (3) |
|---|---|---|---|---|---|
| H227.1 myeloid | 169 | 0.581 | 0.576 / 0.637 | 0.438 | 0.597 |
| H227.2 T/NK | 227 | 0.539 | 0.509 / 0.569 | 0.199 | 0.597 |
| H228.4 CD8-negative T | 214 | 0.532 | 0.513 / 0.578 | 0.299 | 0.597 |

- **Descriptive:** B / plasma 0.420 (responsive median 0.490), mast 0.514 (0.604), CD8 T 0.443 (0.492), NK 0.492
  (0.513): below the responsive median. Anchors: fibroblast 0.686 (q95 0.638, p 0.015), endothelium 0.677 (0.659, p
  0.025). Positive control tumour endothelium vs normal fibroblast 0.642 (q95 0.607, p 0.030). Tumour immune vs normal
  stroma: T/NK vs endothelium 0.481 (q95 0.526), myeloid vs endothelium 0.510 (0.555), T/NK vs fibroblast 0.441 (0.552),
  myeloid vs fibroblast 0.467 (0.590).
- **Per cancer (descriptive, >= 12 patients), clock minus responsive median:** T/NK lung +0.136, liver +0.147, prostate
  +0.111, gastric -0.049, colorectal -0.099; myeloid lung +0.110, prostate +0.088, liver -0.038, colorectal -0.113.
- **Reading (registered):** the patient's tumour and normal immune cells share a state, but not more on the clock genes
  than on glucocorticoid / immediate-early genes: "shared patient state", without the clock. The A227 T/NK line is no
  longer a clock-specific lymphocyte result. Stromal compartments (fibroblasts, endothelium) stay clock-specific against
  both families, as in the healthy atlases. Tumour and normal samples share one operation (anaesthesia, steroid
  premedication); what the responsive shared state is in these samples is not tested.
<!-- src: paper/exploratory/A227_tme_compartment_pairing/output/responsive/{tests,pooled}.csv; paper/exploratory/A228_immune_body_time/output/responsive/tests.csv -->

## A272 (scored 2026-11-16; not blind for the anchored readout): without a host, stromal-vascular cells share more clock state with the other cells of a tissue than lymphocytes do, beyond both companion families (H272.2 PASS); their shared state itself is above expression-matched genes but not above stress-responsive genes (H272.1 fail against the responsive family only)

- Design: 5 healthy atlases (lung_hlca, heart, organ_breast, organ_muscle, organ_colon), 249 pairs of compartments,
  no host and no anchors; pairs (stromal-vascular, X) 60 in 5 atlases, (lymphocyte, X) 60 in 4. Primary statistic
  without the Poisson correction (design addendum, from the placebo run: the corrected variance is non-positive in 144
  of 249 clock pairs). Placebo (donors permuted): T_SV 0.161 against companion medians 0.150 and 0.153; D 0.005.
- **H272.1 fail against the responsive family only.** T_SV (median |r| of stromal-vascular cells with the other
  non-lymphocyte families, mean over atlases) 0.499; matched sets median 0.304, q95 0.399 (p 0.010); responsive sets
  median 0.434, q95 0.533 (p 0.149).
- **H272.2 PASS.** D = T_SV - T_LY 0.277; matched median 0.049, q95 0.117 (p 0.005); responsive median 0.089, q95 0.177
  (p 0.005). Lymphocytes with the other families: clock median |r| 0.220, matched 0.219, responsive 0.349.
- **D272.1.** Per-compartment centrality, family medians (clock; excess over matched; excess over responsive):
  mural 0.574 (+0.214; +0.081), fibroblast-stromal 0.529 (+0.260; +0.145), endothelial 0.473 (+0.165; +0.086),
  epithelial 0.415 (+0.164; +0.041), myeloid 0.409 (+0.153; +0.011), lymphocyte 0.243 (+0.028; -0.085). One compartment
  of 52 is above the q95 of both families (a fibroblast). By atlas, (stromal-vascular, X): heart 0.550, breast 0.585,
  muscle 0.545, lung 0.528, colon 0.289 (colon matched 0.218, responsive 0.442). Spearman of the pairwise excess over
  matched sets with the A150 anchored excess: 0.652 (52 compartments).
- **D272.2 offsets.** 14 pairs are above the q95 of both families; their offsets have median |offset| 0.62 h and a
  median 90 % interval width of 1.54 h; the largest is lung AT2 minus fibroblast, +1.55 h (0.87 to 2.20). One closed
  triangle only (closure 0.02 h): closure is not informative at this number of pairs.
- **D272.3.** After the A150 ambient over-correction of the clock counts: T_SV 0.423, D 0.235 (companions not
  corrected, so not a test). Poisson-corrected T_SV 0.834 over the pairs where it is defined; D not defined.
- **Reading (registered).** The backbone holds without a host as a ranking: the same partner cell types share more
  clock-gene state with stromal-vascular cells than with lymphocytes, and no companion family shows a difference of
  that size (H272.2). The level of the stromal-vascular shared state is not separable from the daily or acute hormone
  response in these atlases (H272.1 against the responsive sets; A270: such sets read time in the baboon, hardly in
  GTEx): it is reported as a shared state above expression-matched genes, not as a clock-specific level. The anchored
  result was not an effect of the anchoring (0.652 with the anchor-free excess). Lymphocytes share no more clock-gene
  state with the tissue than expression-matched genes do, and less than responsive genes.
- **Limits.** Five atlases, lung contributes the most pairs; colon is weak. |r| is attenuated by counting noise, so
  absolute values are lower bounds and families are compared through the companions only. Compartments of one sample
  share ambient RNA. Offsets: intervals of about 1.5 h, no lead or lag is interpreted (separate amendment after a
  calibration against timed truth).
- Code: `exploratory/A272_backbone_pairwise/scripts/score_a272.py` (one scoring run, 30 s; placebo runs before it).
<!-- src: paper/exploratory/A272_backbone_pairwise/output/{tests,pairs,centrality,family_matrix,closure,design}.csv, output/placebo/tests.csv -->

## A273 (scored 2026-11-16; calibration): untimed pairwise offsets between cell types of one sample are unbiased against known time in the mouse colon (V273.1, V273.2 PASS)

- Design: colon chrono-atlas, 8 mice per arm (ZT 1 / 7 / 13 / 19), 12 compartments in C-DC and in R-DC (66 pairs each),
  10 in R-PC. Primary set: 43 pairs of C-DC and R-DC with both compartments timed-rhythmic (timed coupling >= 0.5).
- **V273.1 PASS:** median |untimed - timed| 0.31 h; 97.7 % of pairs within 1 h; median jackknife SD 0.40 h.
- **V273.2 PASS:** slope of untimed on timed offsets 1.005 (Pearson 0.994). No shrinkage towards 0.
- By arm (timed-rhythmic pairs): C-DC 28 pairs, median error 0.24 h, slope 1.195, Pearson 0.959, with true offsets
  up to 1.79 h only; R-DC 15 pairs, 0.57 h, slope 1.000 (true offsets up to 11.6 h after injury); R-PC 3 pairs, 0.07 h.
- All 132 pairs of the two arms, including compartments without a timed rhythm: median error 0.84 h, slope 0.628,
  median |r| 0.475: an offset is not readable where a compartment has no rhythm; the A272 rule (offsets only for pairs
  above both companions) addresses this.
- **Descriptive, the timed offsets themselves.** In the undamaged colon the timed-rhythmic cell classes differ in clock
  phase by a median of 0.63 h and at most 1.79 h. Timed coupling to ZT: smooth muscle 0.970, endothelial 0.946,
  fibroblast 0.931, colonocyte 0.878, lymphatic endothelium 0.815, secretory 0.805, stem 0.803, macrophage 0.514,
  plasma 0.495, DC 0.223, B 0.167, T 0.000: the stromal-vascular classes are the most time-locked and the lymphocytes
  the least, with time known and no host.
- **Reading (registered).** Pairwise offsets between cell types of one sample are unbiased at the 1 h scale; the human
  offsets of A272 are limited by their own intervals. Bias only: precision in a narrow sampling window is not tested
  here (8 mice over 24 h).
- Code: `exploratory/A273_pairwise_offset_truth/scripts/score_a273.py` (one run).
<!-- src: paper/exploratory/A273_pairwise_offset_truth/output/{tests,pairs,compartments,design}.csv -->

## A275 (scored 2026-11-16; post hoc, descriptive): glucocorticoid targets and immediate-early genes are both shared between all cell types of a donor, lymphocytes included; the clock genes are shared between stromal-vascular and other tissue cells and not by lymphocytes

- Per-gene sharing (correlation over donors of a gene's expression in two compartments; median over genes and pairs,
  mean over the five A272 atlases), by pair class:

  | gene class | (SV, X) | (LY, X) | (SV, SV) | (SV, LY) |
  |---|---|---|---|---|
  | clock genes (14) | 0.333 | 0.172 | 0.461 | 0.151 |
  | glucocorticoid targets | 0.524 | 0.390 | 0.665 | 0.389 |
  | immediate-early genes | 0.535 | 0.392 | 0.647 | 0.377 |
  | expression-matched genes | 0.237 | 0.167 | 0.350 | 0.150 |

- Most shared genes in (SV, X) pairs: FKBP5 0.807, ZBTB16 0.711, ZFP36 0.697, PER1 0.690, PDK4 0.640, FOSB 0.604,
  FOS 0.604. Clock genes: PER1 0.690, PER3 0.514, NR1D2 0.450, DBP 0.425, NR1D1 0.383, PER2 0.365, TEF 0.331,
  CRY1 0.308, ARNTL 0.263, CLOCK 0.259, CIART 0.257, RORA 0.252, CRY2 0.217, NPAS2 0.061.
- **Reading (rule fixed in advance): both.** The responsive companion's level in A272 is carried by a hormonal state
  (FKBP5, ZBTB16, PDK4) and by an immediate-early state (ZFP36, FOS, JUN) alike; it is a handling state in part and a
  systemic state in part, and neither is separated from internal time here. What is specific to the clock genes is the
  pattern over cell types: responsive genes are shared with lymphocytes (0.39) about as with other cells, clock genes
  are shared with lymphocytes at the level of expression-matched genes (0.172 against 0.167) and above it elsewhere
  (0.333 against 0.237). PER1, the most shared clock gene, is also a direct glucocorticoid and immediate-early target.
- Limits: single-gene correlations, not corrected for counting noise (compare within a column and with the matched row
  only); post hoc.
- Code: `exploratory/A275_shared_gene_classes/scripts/score_a275.py`.
<!-- src: paper/exploratory/A275_shared_gene_classes/output/{classes,genes}.csv -->

## A271 Part B (scored 2026-11-16; NOT BLIND on the clock): "stronger clocks are more coupled" is clock-specific in the GTEx periphery (H271.2 PASS)
Run `score_a271_a104.py` (6 workers, 47 min, peak 2.2 GB); 244 donors, 40 tissues, 29,478 genes for matching; 0 errors;
the clock rows reproduce A104 exactly (periphery rho 0.734, 27 tissues; brain -0.115, 13).

| block | clock rho | matched median / q95 (p) | responsive median / q95 (p) | verdict |
|---|---|---|---|---|
| periphery | 0.734 | -0.160 / 0.258 (0.005) | 0.109 / 0.453 (0.005) | **H271.2 PASS** |
| brain (descriptive) | -0.115 | 0.220 / 0.742 (0.80) | -0.063 / 0.765 (0.53) | - |

- **Reading (registered):** the positive amplitude-coupling relation across peripheral tissues is a property of the
  clock genes, not of the shared z2 V term: any matched or responsive gene set gives rho near 0. The A104 panel stays.
- **Descriptive (not registered):** median tissue coupling in the GTEx periphery, clock 0.675 against matched 0.107
  and responsive 0.267 (set medians); brain 0.667 against 0.257 and 0.481. In GTEx bulk tissues (every tissue an anchor)
  the clock-gene coupling lies far above both companion families, unlike the single-cell atlas compartments (A263
  Part B).
<!-- src: paper/exploratory/A271_responsive_companions/output/{a104_tests,a104_sets}.csv -->

### A271 Part B, post hoc (not registered; chosen after H271.2 was scored): not a count artefact; correction on replication
Run `posthoc_a104_counts.py` (seconds). Per GTEx tissue of the A104 table, medians over donors of the summed clock-gene
counts (24,133 to 132,043 in the periphery), the smallest clock-gene count of a row and the library depth.

| periphery, 27 tissues | Spearman rho with own amplitude A | with coupling c |
|---|---|---|
| shared term z2 V (range 0.12 to 10.3) | 0.960 | 0.863 |
| own term 2 psi (range 0.90 to 4.14) | 0.504 | -0.085 |
| clock-gene counts | -0.106 | 0.093 |
| smallest clock-gene count | -0.323 | -0.228 |
| library depth | 0.178 | 0.134 |
| partial rho(A, c given clock-gene counts) | 0.752 | |

- **Counts:** neither A nor c rises with the clock-gene counts (A falls slightly with the smallest gene count), and the
  relation is unchanged with counts held fixed (0.752 vs 0.734). Low counts do not produce it.
- **Decomposition:** A tracks the shared term (0.960); the own term varies about 4-fold against about 85-fold for the
  shared term, so c follows how much shared rhythm a tissue carries. This reproduces the archive decomposition (A104).
  Brain block (13 tissues, descriptive): rho(A, clock counts) -0.533, partial rho(A, c) 0.197.
- **Across species:** the positive relation replicates.
  "H104.2 not replicated" in the archive refers to the predicted rho < 0. The observed rho is positive in sheep
  (+0.62, 35 tissues) and baboon timed gain (+0.77, 35), and +0.22 over 145 atlas cell classes (archive
  `docs/RESULTS_preliminary.md`, A104). The companion control (H271.2) covers GTEx only.
- **Reading:** of the two sources of misalignment that scale oppositely with amplitude (METHOD.md section 3, note
  2026-11-16), period mismatch predicts rho < 0 and phase noise or a common systemic drive rho > 0. The pre-registered
  rigidity prediction fails in three species. Which of noise and drive dominates is not separated here. A
  re-alignment-after-shift design would separate them (noise model: strong clocks re-align more slowly).
<!-- src: paper/exploratory/A271_responsive_companions/output/{a104_counts_rho,a104_counts_tissues}.csv; archive scCircadianCancer/docs/RESULTS_preliminary.md (A104 table) -->

## A274 (scored 2026-11-16; blind for these atlases, one scoring run): the cell-type ranking replicates in five atlases not used by A272 (H274.1 PASS), the stromal-vascular level is now above both companion families (H274.2 PASS), and the cell classes of the human cortex share a clock-gene state beyond both families (H274.3 PASS)

- Design (addendum): eight eligible atlases, each reproducing every clock count of its clock-gene table. Body: optic
  nerve (83 donors), anterior eye (67), kidney (26), Kong ileum (30), Kong colon (21). Cortex: PsychAD HBCC (145),
  SEA-AD DLPFC (44) and MTG (47; same cohort as DLPFC). Not eligible by the registered rule: islets (3 of 118 rows
  differ), lung control (no exact twin). 282 body pairs, 92 cortex pairs.
- Placebo (donors permuted within compartment and stratum): D 0.053 against q95 0.064 / 0.056 (matched /
  responsive), T_SV 0.209 against 0.217 / 0.220, T_CX 0.128 against 0.134 / 0.127 (p 0.050 against the responsive
  family, the one placebo test of six at its q95).
- **H274.1 PASS (replication of H272.2):** D = T_SV - T_LY 0.152 over optic nerve, ileum and colon; matched q95 0.104
  (p 0.020), responsive q95 0.109 (p 0.010). Per atlas, (SV, X) against (LY, X): optic nerve 0.510 / 0.335, ileum
  0.307 / 0.142, colon 0.436 / 0.322.
- **H274.2 PASS against both families:** T_SV 0.478 over the five body atlases; matched median 0.269, q95 0.344
  (p 0.005); responsive median 0.337, q95 0.405 (p 0.005). In A272 the level did not exceed the responsive family
  (p 0.149). Per atlas, the (SV, X) median against the responsive q95: optic nerve 0.510 against 0.438 (above); colon
  0.436 against 0.441, eye 0.611 against 0.627, kidney 0.524 against 0.615, ileum 0.307 against 0.375 (not above):
  the pooled statistic passes, four of the five atlases do not on their own.
- **H274.3 PASS:** T_CX 0.497; matched median 0.330, q95 0.466 (p 0.020); responsive median 0.321, q95 0.388
  (p 0.005). Per atlas (pairs of different classes): PsychAD 0.423 (matched q95 0.465, responsive 0.392), DLPFC
  0.557 (0.612, 0.510), MTG 0.510 (0.607, 0.486): each cohort is above the responsive q95 and below the matched q95
  on its own; the margin against the matched family comes from pooling.
- **Family matrix (descriptive; clock / matched / responsive set medians).** Body: fibroblast-endothelial 0.687 /
  0.271 / 0.419; fibroblast-myeloid 0.559 / 0.290 / 0.360; endothelial-epithelial 0.558 / 0.245 / 0.302;
  endothelial-myeloid 0.380 / 0.264 / 0.324; lymphocyte with fibroblast 0.156, endothelial 0.165, epithelial 0.221,
  myeloid 0.207, lymphocyte 0.214, all at or below their matched medians (0.211 to 0.255). Cortex: neuron-neuron
  0.865 / 0.692 / 0.631; neuron-endothelial 0.493 / 0.253 / 0.243; neuron-myeloid 0.510 / 0.320 / 0.330;
  glia-myeloid 0.517 / 0.316 / 0.347; glia-endothelial 0.470 / 0.258 / 0.295; neuron-glia 0.444 / 0.429 / 0.321;
  glia-glia 0.397 / 0.420 / 0.396. Neuron-glia and glia-glia pairs are not above the matched genes.
- **Offsets (descriptive; pairs above both companions).** Body: 46 of 282 pairs, median |offset| 0.63 h, largest
  4.32 h, median 90 % interval width 1.70 h, 63 % of intervals include 0. Cortex: 11 of 92 pairs, median |offset|
  1.02 h, largest 1.68 h, no interval includes 0. In PsychAD the three neuron classes lie within 0.6 h of each
  other and the vascular class differs from them by 1.0 to 1.5 h; in SEA-AD DLPFC the one readable neuron-vascular
  pair has the opposite sign (1.33 h [0.31, 2.52] against -1.53 h [-2.17, -0.94] for the same pair in PsychAD). The
  cortex offsets therefore do not agree between cohorts; no order of cortical cell classes is claimed.
- **Reading (registered).** The ranking (stromal-vascular cells share the clock state with the tissue, lymphocytes
  do not) is replicated in independent atlases and is written as a result. A shared clock state links the cortical
  cell classes: the brain is internally coherent, whatever its coupling to the body (A263).
- **Limits.** The lymphocyte side of H274.1 rests on the two gut atlases and three optic-nerve pairs; eye, optic
  nerve and cortex are post-mortem; nuclei in cortex; the two SEA-AD regions share donors; the cortex margin
  against matched genes is small (0.497 against 0.466) and the placebo sat near its own q95.
- Code: `exploratory/A274_celltype_map/scripts/score_a274.py` (design, placebo, one scoring run); loaders
  `common/datasets/_atlas/genome.py`.
<!-- src: paper/exploratory/A274_celltype_map/output/{tests,per_atlas,family_matrix,pairs,twin,design}.csv, placebo/tests.csv -->

## A276 (scored 2026-11-16; calibration): the pairwise test is anti-conservative (V276.1, V276.2 FAIL); cause: under no pairing the clock's statistic sits above its companions' by a shift of location, from the shape of the clock-plane states; the verdicts of H272.1-H274.3 are suspended pending the corrected test (A277)

- Design as registered: the 13 atlases and the exact rows, pairs and companion sets of A272 and A274; 500 placebo
  draws (donor labels permuted within compartment x stratum). Draw 0 (no permutation) reproduces all ten scored
  statistic x family values and q95 of A272 and A274 exactly (`repro.csv`). A six-draw run of the code was seen before
  the 500-draw run (its draws are draws 1-6 of the 500).
- **V276.1 FAIL, V276.2 FAIL.** Rejection rate over 500 draws (band <= 0.075 each; mean <= 0.065):

| statistic | matched | responsive | mean percentile of the clock (matched / responsive) |
|---|---|---|---|
| H272.1 T_SV | 0.078 | 0.154 | 0.580 / 0.647 |
| H272.2 D | 0.060 | 0.092 | 0.559 / 0.579 |
| H274.1 D | 0.062 | 0.086 | 0.534 / 0.535 |
| H274.2 T_SV | 0.100 | 0.112 | 0.588 / 0.580 |
| H274.3 T_CX | 0.102 | 0.088 | 0.616 / 0.558 |

  Mean rate 0.093. Per atlas and pair class the rates are closer to 0.05 (median 0.056 matched, 0.068 responsive;
  largest organ_colon SV_X 0.138 responsive): averaging over atlases shrinks the spread but not the shift, so the
  pooled statistics are the most affected.
- **Cause (registered diagnostics, `diagnose_a276.py`).** Location, not spread: the spread of the companions within a
  draw equals the spread of the clock over draws (e.g. H272.1 T_SV 0.0196 against 0.0174; H274.3 0.0111 against
  0.0102), while the clock's null statistic lies 0.12 to 0.55 SD above the companion median (H272.1 T_SV 0.162 against
  0.156 matched, 0.152 responsive). The clock-plane states are more isotropic (median anisotropy 0.760 against 0.821
  matched, 0.808 responsive), lighter-tailed (kurtosis ratio 2.30 against 2.74 and 2.52) and carried by more genes
  (effective number 5.8 against 5.2 and 4.7). Null |r| of two independent complex vectors falls as they become
  one-dimensional (E|r|^2 = 1/df whatever the shape, E|r| not); per pair the clock is above the companion mean in 59 %
  (matched) and 55 % (responsive) of 623 pairs, a shift of a few percent of |r| that the class medians and the mean over
  atlases turn into the rates above. Not a property of the clock's biology: a property of the statistic.
- **Consequence (registered reading).** The scored verdicts H272.1, H272.2, H274.1, H274.2, H274.3 are suspended, not
  reversed. The margins of the scored clock values over the companion q95 are large next to the null shift (e.g.
  H272.2 D 0.277 against q95 0.177, null shift 0.006), but the verdicts are restated only by the corrected test (A277).
- Code: `exploratory/A276_pairwise_calibration/scripts/{calibrate_a276,diagnose_a276}.py` (3 min, 6 workers; 4 min).
<!-- src: paper/exploratory/A276_pairwise_calibration/output/{tests,repro,per_atlas,draws,pair_null,diag_summary,diag_compartments}.csv -->

## A277 calibration (2026-11-16): the per-set null correction removes most of the bias but not all (V277.1 FAIL: responsive family 0.080 and 0.078; V277.2 PASS, mean 0.058); restated tests not read (registered stop)

- mu_j from N1 (500 draws), rejection over N2 (500 draws), band <= 0.075: H272.1 T_SV 0.056 / 0.080 (matched /
  responsive), H272.2 D 0.060 / 0.074, H274.1 D 0.036 / 0.038, H274.2 T_SV 0.046 / 0.054, H274.3 T_CX 0.058 / 0.078;
  mean 0.058.
- **Residual cause (descriptive, from the saved draws).** The responsive sets come from a pool of about 28 genes and
  overlap: within a draw their excesses correlate (mean 0.045 to 0.074, matched sets 0.002 to 0.009), so the
  reference spread within a draw is narrower than the clock's spread over draws (H272.1 T_SV 0.0166 against 0.0179;
  H274.3 0.0101 against 0.0107); the clock's own null spread is also slightly wider than a typical set's (0.0179 against
  0.0170). A q95 rule on correlated companions cannot be calibrated by a location correction alone.
- The restated tests H277.1-5 were not computed (the script stops before them). By the A277 rule the line stops for a
  report; the next amendment (A278) calibrates the test by construction and its real-data reading waits for the user.
- Code: `exploratory/A276_pairwise_calibration/scripts/corrected_a277.py` (5 min, 6 workers).
<!-- src: paper/exploratory/A276_pairwise_calibration/output/a277/{calibration.csv,null_draws.npz} -->

## A278 calibration (2026-11-16): the percentile test calibrated against its placebo distribution misses its band in one of ten combinations (V278 FAIL: H272.2 D responsive 0.082; mean 0.043); the band was too narrow for a split-sample check; restated tests not read

- Null pi from N1 (500 draws), rejection (p_cal < 0.05) over N2 (500): H272.1 T_SV 0.026 / 0.034 (matched /
  responsive), H272.2 D 0.074 / 0.082, H274.1 D 0.030 / 0.032, H274.2 T_SV 0.034 / 0.018, H274.3 T_CX 0.044 / 0.060;
  mean 0.043. Null rate of pi = 1 (clock above all 200 companions): 0.005 to 0.028, so the smallest attainable p_cal
  is about 0.006 to 0.03.
- **Reading.** Exchangeable draws make this test exact up to sampling error, and the rates scatter around 0.05 (eight
  of ten at or below it). The registered band (0.075) ignored the error of the N1 threshold: with both halves at 500
  draws the SD of a rate is about 0.014, and the chance that one of ten combinations exceeds 0.075 at a true rate of
  0.05 is about 0.3. V278 is recorded as failed by its rule; the band, not the test, is the likely cause. Not decided
  here: the next step is the user's (A278 as is, with 1000 draws as the null; or a larger null).
- **Scope of every placebo check (A276-A278).** The placebo removes all pairing; it calibrates the test under "no
  shared state for any gene set". The registered null of H272-H274 is "the clock shares no more than its companions",
  with companions that do share (|r| 0.3 to 0.4). Shape effects on |r| shrink as true sharing grows, so the placebo
  bias is an upper bound for the real null, but this is not shown.
- Restated tests not computed. Code: `calibrated_a278.py` (null only; `--read` not run).
<!-- src: paper/exploratory/A276_pairwise_calibration/output/a278/calibration.csv -->

## A279 (scored 2026-11-17; not blind for A201's per-cell-type shifts): untimed pairwise offsets are unbiased in the mouse aorta too (V279.1 PASS: median error 0.30 h, slope 1.04 over a 4.5 h span); macrophages do not couple measurably less to known time than the stromal-vascular cells (H279.1 FAIL); in the Bmal1 knockout the stromal-vascular pairs keep |r| 0.63 to 0.75 without any clock

- Design: Auerbach 2026 aorta (pools of 2 aortas), library x cell type pseudobulks of the clock genes, >= 30 cells,
  compartments in >= 10 libraries. Aligned arm: 14 libraries, SMC, fibroblast, EC (14 each), macrophage (10); ZT 0, 6,
  12, 18. Plumbing fix before scoring: class labels as in the source ("Fibroblast", "Macrophage"); one scoring run.
- **Timed (aligned):** coupling to known time SMC 0.965, fibroblast 0.925, EC 0.864, macrophage 0.810; phases EC -2.02,
  SMC -3.17, fibroblast -3.45, macrophage -4.83 h (relative to the waveform reference).
- **V279.1 PASS.** Six aligned pairs, all timed-rhythmic: median |untimed - timed| 0.30 h, 6 of 6 within 1 h, median
  jackknife SD 0.46 h; slope of untimed on timed offsets 1.04 (timed span 4.47 h); |r| 0.83 to 0.98. Largest error
  fibroblast-macrophage +0.64 h. Misaligned arm (descriptive, 3 pairs): median error 0.25 h, slope 0.82.
- **H279.1 FAIL.** Macrophage minus the smallest stromal-vascular coupling (EC) -0.055, jackknife 90 % [-0.252, 0.143].
  Macrophages are the least coupled class, but not measurably below EC with 10 libraries.
- **Bmal1 knockout (descriptive).** Timed coupling 0.000 to 0.200 in every class (no clock). Untimed |r| between
  stromal-vascular classes stays high: EC-fibroblast 0.753, EC-SMC 0.629, fibroblast-SMC 0.661; with macrophages 0.11 to
  0.35; untimed offsets then do not track timed ones (median error 3.8 h). A pairwise |r| can be large without a clock:
  pooled libraries share non-clock variation of the clock genes between the stromal-vascular classes. This is the
  reason the human pairwise statistic is read against companion gene sets (A272, A274), never alone.
- **Reading (registered).** With A273 (colon), untimed pairwise offsets between cell types are unbiased against known
  time in two mouse tissues, now over offsets of several hours (slope 1.04). No myeloid-below-stromal order is shown in
  the aorta. The knockout shows that the size of |r| is not a clock measure on its own.
- Limits: pools of 2 aortas, 1-3 pools per sex x ZT; macrophages in 10 libraries; one arm primary.
- Code: `exploratory/A279_aorta_pairwise_truth/scripts/score_a279.py` (run directly, one run; runtime not recorded).
<!-- src: paper/exploratory/A279_aorta_pairwise_truth/output/{tests,pairs,compartments,design}.csv -->

## Data note (2026-11-17, cell-type program A3): genome-wide twins of IBDverse ileum and rectum
- `$SCCC_DATA/atlas_wg/ibdverse` built from the local IBDverse h5ad (1,852,681 cells, raw/X) by
  `exploratory/celltype_map_builds/build_ibd_genome.py` with the grouping of the panel build (A14); 32,363 genes after
  canonical symbols. Builders `ibd_ileum_genome`, `ibd_rectum_genome` (same rules as `_atlas/ibd.py`).
- Twin check over all rows: rectum 1,677 rows, ileum 2,433 rows; every clock count, n_cells and depth equal (maximum
  difference 0), strata and healthy donors equal. Both are eligible by the A274 twin rule.
- Healthy donors with a lymphocyte and a stromal-vascular compartment (builders, counts only): rectum 61, ileum 53.
- First run killed for memory (13 GB > 12 GB cap: a dense temporary per chunk); rerun with in-place accumulation, peak
  10.8 GB.

## Scouts of the cell-type program (2026-11-17; counts and metadata only, no statistic)
- **A3, lymphocyte atlases without a genome table** (healthy A17 targets; donors with a lymphocyte and a
  stromal-vascular compartment): IBDverse rectum 61 and ileum 53 (built, exact twins, A280); skin systemic sclerosis
  controls 27 (GSE195452, MARS-seq), endometrium 21 (GSE266265), atopic dermatitis controls 19 (GSE328048): GEO sources
  with marker-based classes, exact twins not assured (legacy builders would need re-running with all genes); Tabula
  Sapiens 9 (all-gene table on disk, organ-specific compartments). Blood, bone marrow and the immune-aging lymphoid
  organs have lymphocytes but no stromal-vascular class.
- **C3:** PsychAD MSSM (1,032 donors, 257 healthy) and RADC (152 donors, 26 normal; no A17 healthy fit) have no
  all-gene table; built from the Census by `celltype_map_builds/run_psychad.sh` (running). Cortex atlases have no
  lymphocyte or fibroblast class (vascular is the only stromal-vascular class).
- **A4, timed mouse single-cell data on disk:** besides the colon chrono-atlas, only the Auerbach aorta aligned arm
  meets the design rule (scored as A279); skin dermis GSE223109 has two classes only (fibroblast, immune), the other
  candidates have one class or too few libraries.
- **C1, disease arms** (compartment pairs with >= 15 shared donors in each arm; control / case): HLCA fibrosis 198 / 71
  donors, SV_X 52 / 15, LY_X 27 / 19; lung ILD (Natri) 40 / 61, SV_X 15 / 41, LY_X 4 / 23; kidney CKD 26 / 37, SV_X
  10 / 10, LY_X 4 / 10; IBDverse Crohn's ileum 242 / 117, SV_X 4 / 4, LY_X 12 / 14; MASLD pooled (D1) 58 / 72, SV_X
  6 / 6, LY_X 6 / 9; cortex (PsychAD AD, SCZ, BD; SEA-AD AD) only SV = vascular: 7-8 SV_X pairs per arm; Kong Crohn's 6
  cases (no pair). Earlier anchored readouts on these arms: A153 (13 targets, every MDE > 0.2 on d*), A167 (liver
  fibrosis, 0 of 12), A230 (inflamed gut, MDE 1.7 h), A141 (brain contrasts anchor-dependent). No disease contrast is
  registered before Q36 and a power calculation on the calibrated test.
- **C2, the tumour in the pairwise language:** the question (malignant cells with the tumour's own stromal-vascular cells,
  against normal epithelium with the normal stromal-vascular cells) is A149 S3 / S3n without a host. A149 had 51
  patients in 7 cohorts at S3 (lung the only cancer estimable at S2 to S3n) and showed that non-clock sharing within a
  patient dominates single-cell tumour pairings (H149.1 below its pseudo-clock median). A pairwise version would rest on
  the same patients; its power calculation needs the calibrated null (Q36) and is left until then.
<!-- src: scratchpad scout reports (not in the repo); counts re-checked: IBDverse rectum 61, HLCA fibrosis 198 / 71 -->

## A281 (scored 2026-11-17; not blind for lymphocyte couplings and pairwise |r|): tissue lymphocytes carry no readable clock of their own: their clock genes are not clock-coherent between donors where depth-matched fibroblasts are (H281.1 FAIL: 2 of 9; H281.2 PASS: 5 of 7 powered atlases); every stromal-vascular and epithelial compartment is coherent
- Design: 12 healthy body atlases (A272, A274, IBDverse); per compartment with >= 15 donors, signed clock-gene
  coherence C across donors (`decompose.signed_coherence`, per-sample offset fitted out), 1000 per-gene permutations
  within strata. Nine atlases have a lymphocyte and a fibroblast-stromal compartment on >= 15 shared donors (18 to 117);
  muscle, anterior eye and kidney have none. One plumbing failure before any statistic (workers could not unpickle the
  file-loaded module; fixed as A274), then one scoring run.
- **H281.1 FAIL** (own time): lymphocyte compartments with p < 0.05: 2 of 9 (heart lymphocyte C 0.502, p 0.012; breast
  CD4 T 0.384, p 0.013). The others: lung T 0.237 (p 0.112), Kong colon CD4 T 0.344 (0.072), Kong ileum CD8 T 0.055,
  colon T 0.048, optic nerve T 0.055, IBDverse ileum T -0.426 (0.915), rectum T -0.452 (0.863).
- **H281.2 PASS** (weak clock): fibroblast-stromal compartments thinned to the lymphocytes' clock-gene counts pass in
  >= 80 % of 20 replicates in 7 atlases (lung, heart, breast, Kong ileum and colon, IBDverse ileum and rectum: 0.90 to
  1.00); the lymphocytes fail in 5 of these 7 (0.714 >= 2/3). Depth does not explain it: the IBDverse T cells have more
  clock-gene counts than the stromal cells (median 1,397 against 292 in ileum; 513 against 248 in rectum) and the most
  negative C; their C lies at the 12th and 6th percentile of 200 expression-matched gene sets.
- **All compartments (descriptive), median C and share with p < 0.05:** fibroblast-stromal 0.829 (13 of 13 = 1.000),
  mural 0.806 (1.000), endothelial 0.771 (1.000), epithelial 0.752 (1.000), parenchyma 0.746 (1.000), other 0.688
  (0.769), myeloid 0.396 (0.500), lymphocyte 0.025 (38 compartments, 0.132). Within myeloid (post hoc look):
  macrophages are coherent in most atlases (lung 0.675, heart 0.883, breast 0.626, optic nerve 0.691, eye 0.782 and
  0.689, Kong ileum 0.594, IBDverse 0.579 and 0.565, p <= 0.004), monocytes, mast cells and dendritic cells are not
  (lung monocyte 0.003, DC 0.054; mast cells about 0). Lymphocyte subsets: no subset is coherent in more than one
  atlas (plasma 1 of 6, B 1 of 7, CD8 T 0 of 4).
- **Reading (registered).** At depths where a stromal clock is read, the lymphocyte clock is not: tissue lymphocytes
  carry a weak or absent clock signal between donors, not a clock of their own that the tissue fails to share. With
  A272/A274 (no sharing beyond non-clock genes) and A275 (lymphocytes do share glucocorticoid-target and immediate-early
  states with the tissue), the simplest reading is that tissue lymphocytes respond to systemic signals but do not keep
  a readable clock between donors. Heart and breast lymphocytes are the exceptions.
- **Limits.** Between-donor clock variation depends on the sampling window, which is shared by all classes of an atlas
  (the thinned fibroblasts carry it); signed coherence reads the waveform's gene relations, so a lymphocyte clock with
  other phase relations would read low; lymphocyte classes pool resident and circulating cells (B2). The macrophage
  versus monocyte split is a post hoc observation, the motivation for B2, not a result.
- Code: `exploratory/A281_lymphocyte_clock/scripts/score_a281.py` (3 workers, 96 s, peak 5.2 GB; job ledger).
<!-- src: paper/exploratory/A281_lymphocyte_clock/output/{tests,pairs,compartments,design}.csv -->

## A283 (scored 2026-11-17; user's answer to Q36 (b); not blind for the A272/A274 statistics): the cell-type pairwise test is calibrated with a 5000-draw placebo null (V283 PASS), and every A272/A274 verdict stands on it: the stromal-vascular-over-lymphocyte ranking (H272.2, H274.1), stromal-vascular sharing beyond the matched family (H274.2) and cortical sharing (H274.3) pass; stromal-vascular sharing in the A272 atlases is still not above the responsive companions (H272.1)

- **Calibration (null only).** 5000 placebo draws (seeds [283, d]); one half as the null, rejection (p_cal < 0.05) over
  the other, both directions: 0.031 to 0.063 over the ten combinations (H274.3 T_CX the widest: 0.063 / 0.031 matched,
  0.061 / 0.040 responsive), all <= 0.0685. Null rate of pi = 1: 0.004 to 0.021 (smallest attainable p_cal).
- **Restated tests (null = all 5000 draws; p_cal matched / responsive):**
  - H283.1 = H272.1, T_SV 0.499 (companion medians 0.304 / 0.434): 0.010 / 0.299, **FAIL** (responsive), as before.
  - H283.2 = H272.2, D 0.277 (0.049 / 0.089): 0.007 / 0.014, **PASS**.
  - H283.3 = H274.1, D 0.152 (0.029 / -0.003): 0.020 / 0.014, **PASS**.
  - H283.4 = H274.2, T_SV 0.478 (0.269 / 0.337): 0.009 (decides) / 0.013, **PASS**.
  - H283.5 = H274.3, T_CX 0.497 (0.330 / 0.321): 0.038 / 0.011, **PASS**.
- **Per atlas (descriptive, same draws).** SV_X beats the matched family (p_cal < 0.05) in heart, HLCA, breast,
  kidney, anterior eye and optic nerve, and the responsive family only in optic nerve (0.014); colon, muscle and Kong
  ileum are not above either family. LY_X is above no family in any atlas (lowest p_cal: Kong colon 0.081 / 0.061).
  Cortex CX_diff is above the responsive family in HBCC (0.036) and DLPFC (0.029), not in MTG (0.061), and above the
  matched family in none (0.135 to 0.179); the pooled T_CX passes both.
- **Reading.** The anti-conservative q95 rule (A276) did not produce any of the cell-type verdicts: with a test
  calibrated by construction (rates around 0.05 in 5000 draws) all five keep their original outcome. What is clock-
  specific is the contrast (stromal-vascular cells share clock state more than lymphocytes do, beyond both companion
  families, in two independent atlas sets) and the cortical sharing pooled over three cohorts; the level of
  stromal-vascular sharing is matched by responsive genes in the A272 atlases and beaten in single atlases mostly by
  the matched family only. The A283 form is now the readout of the line (A280, A282 and later tests).
- Limits: the placebo removes all pairing (A278 scope note); pooled statistics weight atlases equally.
- Code: `exploratory/A276_pairwise_calibration/scripts/calibrated_a283.py` (1357 s, peak 3.2 GB, job ledger).
<!-- src: paper/exploratory/A276_pairwise_calibration/output/a283/calibration.csv, tests.csv, per_atlas.csv -->

## A280 (scored 2026-11-17 with the A283 form; blind for these atlases): the stromal-vascular-over-lymphocyte ranking replicates in IBDverse (H280.1 PASS, D 0.413); cortical sharing does not replicate in PsychAD MSSM, the largest cortex cohort (H280.2 FAIL); the offset-sign tie-break is not read (H280.3 gate failed)

- Design (addendum): IBDverse rectum (274 donors) and ileum (242), each SV_X 4 and LY_X 12 pairs; MSSM 257 healthy
  donors, CX_diff 27 pairs; all three exact genome twins. Null: 5000 placebo draws (seeds [280, 9, d]).
- **H280.1 PASS.** D = 0.413 (companion medians 0.016 matched, 0.048 responsive; pi 1.00 in both), p_cal 0.0036 /
  0.0164. The contrast is larger than in A272 (0.277) and A274 (0.152).
- **H280.2 FAIL.** MSSM T_CX 0.279 (medians 0.237 / 0.230; pi 0.645 / 0.775), p_cal 0.455 / 0.271. Not replicated at
  this power; MSSM's clock |r| between cortical families is lower than in HBCC, DLPFC and MTG (0.423 to 0.557, A283).
- **H280.3 not read.** Gate: inhibitory_CGE-vascular |r| 0.231 (95 donors), pi 0.68 / 0.74, p_cal 0.332 / 0.279. The
  offset (-2.16 h, 90 % [-4.35, +0.54]) decides nothing either way and is reported for completeness only.
- **Reading.** The lymphocyte gap now holds in three independent body sources, the third with 53-61 lymphocyte donors
  per tissue, all from living biopsies. The cortex statement weakens: pooled over HBCC, DLPFC and MTG it passes (A283
  H283.5), but per cohort it beats the matched family in none, and the largest cohort shows no clock-specific sharing.
  MSSM shares consortium and protocols with HBCC, so the difference is not one of processing.
- Code: `exploratory/A280_celltype_replication2/scripts/calibrated_a280.py --q36` (1132 s, peak 10.4 GB).
<!-- src: paper/exploratory/A280_celltype_replication2/output/calibrated/tests.csv, tiebreak.csv; output/design.csv -->

## A282 (scored 2026-11-17 with the A283 form; descriptive; not blind for the all-donor statistics): the ranking holds in living tissue alone: D = 0.225 over five living atlases, above both companion families (p_cal 0.011 / 0.011), with stromal-vascular sharing above both companion medians

- Arms (per donor, from atlas metadata; recall-only classes left out): living 7 atlases (D over HLCA, breast, colon,
  Kong ileum, Kong colon), post-mortem 7 (D over HLCA, heart, optic nerve; T_CX over the three cortex cohorts).
- **Question fixed in advance: yes.** Living arm: T_SV 0.450 (medians 0.279 / 0.359; p_cal 0.014 / 0.059), D 0.225 (0.028
  / -0.010; p_cal 0.011 / 0.011), T_LY 0.203 (0.243 / 0.331: below both companion medians).
- Post-mortem arm (descriptive): T_SV 0.547 (p_cal 0.010 / 0.012), D 0.261 (0.006 / 0.007), T_LY 0.264, T_CX 0.497
  (0.034 / 0.009).
- **The one atlas in both arms (HLCA):** SV_X 0.523 living vs 0.516 post-mortem; LY_X 0.195 vs 0.261; D 0.328 vs 0.255.
- **Per atlas, living:** SV_X above the matched family in breast and kidney (p_cal 0.007 each), above the responsive
  family only in HLCA (0.017); LY_X above no family (Kong colon closest, 0.076 / 0.060).
- **Post hoc (seen in the family matrix, different atlases per arm):** stromal-vascular pairs with myeloid cells share
  less in living tissue (fibroblast-myeloid 0.259 over 7 pairs) than post-mortem (0.643 over 13 pairs); lymphocyte pairs
  are low in both. Not tested; arms confound organ, donor type and the split.
- **Reading.** No agonal state shared by a donor's cells is needed for the stromal-vascular-over-lymphocyte ranking:
  it is as large in surgical and biopsy tissue as post-mortem, and within HLCA it is if anything larger in living
  donors. In living tissue lymphocytes share responsive genes with the other cells more than clock genes (T_LY below the
  responsive median).
- Limits: living tissue is surgical or biopsy (ischaemia, anaesthesia); pooled means weight atlases equally.
- Code: `exploratory/A282_pooled_living/scripts/score_a282.py` (1776 s, peak 2.9 GB).
<!-- src: paper/exploratory/A282_pooled_living/output/tests.csv, per_atlas.csv, family_matrix.csv, arms.csv -->

## A285 (scored 2026-11-17; calibration, not blind): in the bulk and tumour results the expression-matched q95 rule is calibrated, the stress-responsive one is not (placebo rejection 0.076 to 0.204); on the calibrated test the clock-specific verdicts stand for A138, A227, A255 surgical, and against matched genes everywhere, but A251 (skin sites) and A255 MuTHER fat-skin are no longer shown to be clock-specific against the responsive genes

- Placebo: pairing unit permuted within the scored strata, 2000 draws per test, every companion set re-paired
  identically; draw 0 reproduces each scored statistic and both companion q95 exactly (`repro.csv`).
- **V285, matched family (calibrated in every test):** rejection 0.0465 to 0.0595 (A251 0.047-0.057, A255 0.052-0.054,
  A138 0.050, A227 0.052-0.060). The user's expectation for large bulk cohorts holds for this family.
- **V285, responsive family (biased in every test but A227 H227.1):** A251 0.184 to 0.204; A255 MuTHER 0.131, surgical
  0.088; A138 0.076; A227 0.062 / 0.065 (not in its rule). Cause: the responsive sets are drawn from a pool of about 30
  genes, overlap, and so move together; with one to six compartment pairs, a clock above all of them under no pairing
  is not rare (placebo rate of pi = 1: A251 0.063-0.084, MuTHER 0.056).
- **H285, calibrated restatement (p_cal matched / responsive):**
  - A251 H251.1 0.004 / 0.075; H251.2a 0.007 / 0.084; H251.2b 0.008 / 0.064: the clock is above all 200 responsive
    sets in each, but that is attainable under no pairing in 6-8 % of draws. **Clock-specific against responsive
    genes: not shown (reversed from A251's "clock-specific");** against matched genes: stands.
  - A255 H255.1 MuTHER fat-skin 0.011 / 0.057 (floor 0.056): **not shown against responsive genes (reversed)**;
    H255.3 surgical pooled 0.005 / 0.017: **stands**.
  - A138 H138.4 0.004 / 0.011: **stands** (both families).
  - A227 H227.1 0.018 (matched) and H227.2 0.035: **stand**; responsive 0.439 / 0.200 (not clock-specific, as A271).
- **Reading.** The bias the user suspected is real but confined to the responsive family, and it matters where the
  statistic rests on few pairs (one pair in MuTHER, three to six in A251). There the responsive comparison cannot
  reach p < 0.05 at all: the result is "not decidable against responsive genes", not "matched by them". The pairing
  itself (each test's own permutation p) and the matched-family specificity are unchanged everywhere.
- **Not run: A263 Part A** (registered reason: per-draw refits, L undefined under no pairing); its verdicts are failures
  by wide margins.
- Code: `exploratory/A285_companion_calibration/scripts/placebo_a285.py` (A227: 997 s, peak 31.2 GB; others < 1 min).
<!-- src: paper/exploratory/A285_companion_calibration/output/{a251,a255,a138,a227}/{repro,tests}.csv -->

## A286 power note (2026-11-17; permuted labels only, no observed contrast computed): lung fibrosis is not powered for the pairwise readout in either cohort; no primary, A286 closed at the power stage

- **HLCA** (`hlca_fibrotic_genome`; healthy rows identical to `lung_hlca_genome`): 200 controls, 71 cases (65 lung
  explants). Only 4 study|assay strata hold both arms (97 donors): elsewhere arm = study, which the registered null
  cannot remove (S2 null mean 0.080 with all controls). Within study (48 controls, 49 cases; control SV_X 5 pairs,
  endothelial-X 0): MDE S1 0.238, S2 0.274, S3 0.387. Living controls only: 1 stratum with both arms (13 donors).
- **Natri ILD** (`lung_natri_genome`, exact twin of `lung_natri`): 40 controls, 61 cases, 2 processing sites with both
  arms (101 donors); control arm 23 pairs (SV_X 9, LY_X 2, fibroblast-X 0): MDE S1 0.178, S2 0.274, S3 0.228 (S3
  undefined in the control arm).
- **Rule (registered):** primary = S3, else S2, if MDE <= 0.10: neither cohort reaches it.
- **Reading.** Lung fibrosis cannot be read with the pairwise readout on public single-cell data: the control arms
  are small where they share a study with the cases (healthy lung atlases separate cases and controls by study), and
  every case is an explant taken at surgery while controls are donor lungs (sampling-time spread). Not a negative
  result. The tables (`hlca_fibrotic_genome`, `lung_natri_genome`) stay for later designs.
- Code: `exploratory/A286_lung_fibrosis_pairwise/scripts/score_a286.py --stage power [--cohort natri]` (HLCA 350 s,
  Natri 801 s).
<!-- src: paper/exploratory/A286_lung_fibrosis_pairwise/output/{power,design}.csv, output/natri/{power,design}.csv -->

## Data note (2026-11-17, from A284): Kong ileum "epithelial" contained the intraepithelial lymphocytes
- The coarse class rule (`organ.COARSE`, first match on Census `cell_type`) matches "intraepithelial lymphocyte" as
  epithelial. In `kong_ileum` / `kong_ileum_genome` the `ileum:epithelial` compartment therefore includes the IEL T
  cells: a median 9.6 % (max 23.2 %) of its cells per donor (5,004 of 55,991 cells). Kong colon has no such label.
- Affected readouts (Kong ileum epithelial rows): A274 (replication atlas), A281 (coherence per compartment), A282,
  A283 (per-atlas rows). Not rerun here; a T admixture would make the epithelial class look more lymphocyte-like, so
  any bias works against the stromal-vascular-over-lymphocyte contrast at the epithelial pairs. A284 uses the
  corrected rows. The class rule is left unchanged (other tables would move); to be fixed with a rebuild if the PI
  wants the earlier rows corrected.
<!-- src: paper/exploratory/A284_residency/output/table_checks.csv; check in the A284 build session (IEL share) -->

## A284 (scored 2026-11-17; B2 residency; labels fixed before scoring, scored once): sharing of the tissue's clock state does not increase with residency: resident T cells (H284.1 FAIL), plasma cells (H284.2 FAIL, wrong sign in every scorable atlas) and macrophages (H284.3 FAIL after Holm) do not share more than their circulating counterparts

- **Design (registered):** per atlas, donors with both classes; the deeper immune row thinned to the other's clock-gene
  count (10 replicates); Delta = median clock |r| of the resident class with the tissue classes minus the same for the
  circulating class; pooled = mean over atlases; null = the two rows of each donor swapped (2000 draws); Holm over three.
  Each class meets 4 tissue classes per atlas (lung 9), so each per-atlas Delta rests on 4 (9) pairs per class.
- **H284.1 T resident vs circulating (4 atlases): FAIL.** Pooled Delta 0.010, p_swap 0.337; p_cal matched 0.323,
  responsive 0.573. Per atlas: IBDverse ileum +0.084, breast +0.035, IBDverse rectum -0.050, Kong ileum -0.031.
- **H284.2 plasma vs B (6 registered atlases, 5 scorable): FAIL.** Delta is negative in every scorable atlas: breast
  -0.022, colon (Pelka) -0.069, Kong ileum -0.057, IBDverse ileum -0.010, rectum -0.020 (mean -0.036, post hoc: the
  registered pooled value is undefined, below). Plasma cells share less of the tissue's clock state than B cells.
- **H284.3 macrophage vs monocyte (2 atlases): FAIL.** Pooled Delta 0.092 (lung +0.038, breast +0.146), p_swap 0.033,
  Holm 0.067 as stored (0.067 to 0.100 with a valid H284.2 p); p_cal matched 0.017, responsive 0.073.
- **Scoring fault (found after scoring; no verdict depends on it).** Kong colon (16 donors with plasma and B, eligible by
  the registered >= 15 donors) yields no pair under the shared-donor rule (pooled df >= 15), so its Delta is NaN and
  the registered pooled mean over atlases is NaN. The stored H284.2 p_swap (0.0005) is an artefact of NaN comparisons
  and is not a result; it also lowered the stored Holm value of H284.3. H284.2 fails on sign (Delta < 0 in all five
  scorable atlases, so no pooling of them can be > 0); H284.3 fails Holm whatever the H284.2 p. Not rerun.
- **Reading against the stated outcomes: "none passes".** The lymphocyte gap (A272-A283) is not explained by residency:
  the chapter keeps "lymphocytes are apart" without a mechanism. Against the mouse colon expectation (plasma 0.495 vs B
  0.167, A21/A31) the B-lineage sign is reversed. Macrophage over monocyte is the only contrast in the predicted
  direction in both atlases and above the matched companions (p_cal 0.017), but not significant after correction nor
  above the responsive companions: a lead for a larger myeloid design, not a finding.
- **Limits.** Few pairs per class; author labels (circulating T in tissue may be recirculating memory; breast Tem may
  be tissue effectors); two myeloid atlases only.
- Code: `exploratory/A284_residency/scripts/score_a284.py --draws 2000 --jobs 6` (740 s, peak 6.2 GB).
<!-- src: paper/exploratory/A284_residency/output/{tests,per_atlas,design,table_checks}.csv; paper/results/logs/jobs.tsv -->

## A287 (scored 2026-11-17; not blind for A281): tissue lymphocytes have no clock-gene covariance beyond expression-matched genes, whatever its phase pattern, where depth-matched fibroblasts show it (H287.1 FAIL: 0 of 3 powered atlases; H287.2 PASS: 3 of 3 fail); "no readable clock", not "a different clock"

- Statistic (registered): U = share of the clock genes' correlation matrix in its two leading eigenvalues, after the
  per-sample shift is removed (no waveform); a compartment carries structure when U is above 95 % of 200 matched sets
  and above per-gene permutations (p < 0.05). Nine A281 pairs; fibroblast-stromal rows thinned to the lymphocytes'
  clock-gene count, 20 replicates.
- **Power: 3 atlases powered** (thinned fibroblasts carry structure in >= 80 % of replicates): Kong ileum 1.00, IBDverse
  ileum 1.00, rectum 1.00 (IBDverse fibroblasts are shallower than the T cells, so they are not thinned). Not powered:
  heart 0.50, breast 0.40, lung, Pelka colon, optic nerve, Kong colon 0.00. The registered minimum is 3.
- **H287.1 FAIL / H287.2 PASS:** lymphocytes carry structure in 0 of the 3 powered atlases; matched-set pi: Kong ileum
  CD8 T 0.335, IBDverse ileum T 0.740, rectum T 0.110 (U 0.427, 0.551, 0.350 against matched medians 0.449, 0.496,
  0.424). In IBDverse the T cells have more clock-gene counts than the fibroblasts (median 1,397 against 292 in ileum;
  513 against 248 in rectum) and still none.
- **Descriptive, all nine pairs unthinned:** no lymphocyte compartment carries structure (matched pi 0.110 to 0.845;
  every p_perm <= 0.003: the clock genes covary in lymphocytes, but no more than expression-matched genes do).
  Fibroblast-stromal: 5 of 9 carry structure (heart, optic nerve, Kong ileum, IBDverse ileum and rectum, pi 0.985 to
  1.000), and their leading-plane gene phases match the GTEx waveform without being told it (circular correlation 0.70
  to 0.88). Lung (pi 0.675), breast (0.785), Pelka colon and Kong colon fibroblasts do not carry structure on U: the
  statistic is less sensitive than the signed C of A281, which passes them.
- **Reading (registered): "no readable clock" holds for any phase pattern.** In the three atlases where a fibroblast
  clock is read at the lymphocytes' depth, the lymphocytes' clock genes covary between donors at the level of random
  expression-matched genes, so A281's failure is not a clock with other gene relations. With A275, tissue lymphocytes
  carry systemic (glucocorticoid, immediate-early) states but no clock between donors.
- **Limits.** Three powered atlases (the minimum), all gut; U misses the fibroblast clock in four atlases, so it is the
  weaker of the two statistics; the between-donor clock spread depends on each atlas's sampling window.
- Code: `exploratory/A287_lymphocyte_covariance/scripts/score_a287.py --jobs 6` (56 s, peak 3.6 GB).
<!-- src: paper/exploratory/A287_lymphocyte_covariance/output/{tests,compartments,thinned,phases,design}.csv; paper/results/logs/jobs.tsv -->

## A289 (scored 2026-11-17; NOT BLIND for six tissues and the group medians): the tissues outside the blood-brain barrier couple less, not more (H289.1, H289.2 FAIL, opposite direction); in the baboon the brain sits on the periphery's amplitude-coupling line (H289.3 FAIL, not extrapolated): the brain's lower coupling is what its weaker clocks predict
Run `score_a289.py` (seconds; reads the A42 readout-2 table, no refit); one scoring run.
- **H289.1 FAIL (baboon).** Pineal and pituitary median coupling 0.123 against 22 brain regions 0.420: Delta -0.297,
  permutation p 0.855 (276 assignments). Pineal 0.000, pituitary 0.246.
- **H289.2 FAIL (sheep).** Pineal and pituitary 0.082 against 8 brain regions 0.205: Delta -0.123, p 0.933 (45).
  Pineal 0.000 (43 animals), pituitary 0.164 (123).
- **H289.3 FAIL (baboon).** Peripheral line coupling = 0.265 + 0.350 log(own amplitude); brain residuals median -0.001,
  11 of 22 above (sign test p 0.584). Not extrapolated: brain amplitudes 0.81 to 2.63, 7 peripheral tissues in that
  range (periphery 1.32 to 5.28). Sheep (reported): median residual +0.014, 4 of 8 above (p 0.637), extrapolated (brain
  0.75 to 0.96, no peripheral tissue in range; periphery 0.58 to 2.35, slope 0.344).
- **Descriptive.** Baboon ocular tissues (behind the blood-retina barrier) 0.689: cornea 0.806, RPE 0.851, optic nerve
  head 0.573, retina 0.010. Medians of own amplitude: baboon brain 1.539, glands 1.537, periphery 3.524; sheep brain
  0.846, glands 1.100, periphery 1.334.
- **Reading (registered).** The barrier reading is not supported: the two tissues outside the barrier are the least
  coupled in both species, the pineal not at all. In the baboon, where brain and peripheral amplitudes overlap, the
  brain's lower coupling is what the periphery's amplitude-coupling relation predicts for its weaker clocks; the GTEx
  post hoc observation (13 of 13 brain regions above the line, extrapolated) does not replicate. Pineal and pituitary
  are glands, so H289.1-2 do not isolate the barrier (registered limit).
- **Post hoc (not registered).** A pineal coupling of 0 in both species fits a gland whose clock-gene rhythm is set by
  the SCN through its sympathetic input at night rather than by the body's peripheral state; GTEx has no pineal.

## A288 (scored 2026-11-17; not blind at the family level): macrophages share more clock state with the tissue than lymphocytes do, beyond both companion families (H288.2 PASS, D_MAC 0.148); their sharing level is above expression-matched genes but not above stress-responsive genes (H288.1 FAIL); alveolar over monocyte-derived macrophages misses (H288.3 FAIL, p 0.051)

- Design (registered, A283 form, `clocksync.pairwise`): ten healthy atlases with a macrophage class (IBDverse relabelled
  from the author labels); X = non-immune classes; pooled = mean over atlases with a pair; 200 matched + 200 responsive
  sets per atlas; 5000 placebo draws; Holm over H288.1-2 on the larger of the two family p.
- **H288.1 FAIL:** pooled T_MAC 0.409 (10 atlases); matched median 0.253, p_cal 0.0036; responsive median 0.367, pi
  0.910, p_cal 0.133. Holm 0.133.
- **H288.2 PASS:** pooled D_MAC = T_MAC - T_LY 0.148 (8 atlases); matched median 0.027, p_cal 0.0034; responsive median
  0.056, p_cal 0.0064; Holm 0.0128.
- **Per atlas (descriptive, same draws; p_cal matched / responsive).** T_MAC: lung 0.530 (0.010 / 0.327), breast 0.559
  (0.003 / 0.229), anterior eye 0.648 (0.012 / 0.073), optic nerve 0.519 (0.032 / 0.071), Kong colon 0.480 (0.021 /
  0.008), IBDverse rectum 0.286 (0.018 / 0.422), muscle 0.454 (0.289 / 0.658), IBDverse ileum 0.195 (0.101 / 0.534), Kong
  ileum 0.253 (0.489 / 0.459), Pelka colon 0.163 (0.856 / 0.982). D_MAC: lung 0.294 (0.003 / 0.006), breast 0.234
  (0.019 / 0.168), IBDverse rectum 0.189 (0.017 / 0.034), optic nerve 0.184, Kong colon 0.128, IBDverse ileum 0.101,
  Kong ileum 0.101, Pelka colon -0.043 (the one negative atlas).
- **H288.3 FAIL (HLCA decides):** alveolar minus monocyte-derived macrophages, 92 donors, 8 pairs each: Delta 0.165,
  p_swap 0.051; p_cal 0.077 / 0.123. Descriptive, IBDverse ileum CD163hi minus CD163lo (60 donors): -0.071, p_swap
  0.780. Rectum not eligible (5 donors with both).
- **Reading against the stated outcomes.** The registered combinations were "both pass", "H288.1 passes, H288.2 fails"
  and "H288.1 fails: macrophages are not shown to keep the shared state beyond non-clock genes; the section names
  stromal-vascular cells only". The last applies to the level. **Post hoc (not counted):** this is the pattern of
  stromal-vascular cells in A272 (H283.1 level FAIL against the responsive family, H283.2 contrast PASS); the
  stromal-vascular level test that passed (H274.2 = H283.4) was decided by the matched family alone, and on that rule
  T_MAC passes (0.0036). Macrophages therefore stand where stromal-vascular cells stood in A272: the contrast with
  lymphocytes is clock-specific, the level is matched by stress-responsive genes. The residency split within
  macrophages points the predicted way in lung (Delta 0.165) but misses the threshold, and the gut marker proxy
  points the other way.
- **Limits.** Most macrophage classes include dendritic cells; the matched-family level was known at the family level
  before registration; H288.3 rests on one atlas.
- Code: `exploratory/A288_macrophage_pairwise/scripts/build_a288.py`, `score_a288.py --jobs 6` (998 s, peak 5.9 GB).
<!-- src: paper/exploratory/A288_macrophage_pairwise/output/{tests,per_atlas,design,design_h3,h288_3}.csv; paper/results/logs/jobs.tsv -->

## A291 (scored 2026-11-17; NOT BLIND for the A104 table): measurement does not produce "stronger clocks are more coupled": a perfectly locked tissue would read a coupling of 0.97 at every GTEx tissue's real design (H291.1 PASS)

- Rows and fit as A104 (GTEx fast deaths, 27 peripheral tissues and 13 brain regions, every tissue an anchor,
  leave-one-out host); reproduction exact (periphery rho 0.734, brain -0.115). Per tissue, 100 simulations at true
  coupling 1 and 100 at its estimate with `clocksync.simulate.like` (real host means and posterior covariances,
  measurement covariances, strata, intercepts). One run (4 workers, 31 s).
- **Ceiling.** Perfect locking reads 0.966 to 0.977 in the periphery (median 0.972) and 0.937 to 0.956 in brain
  (0.946); no simulation reached the boundary. Observed couplings span 0.067 to 0.853 (periphery) and 0.308 to 0.706
  (brain). Design bias at the estimate: median -0.013 (periphery), -0.022 (brain).
- **H291.1 PASS.** rho(own amplitude, coupling / ceiling) 0.717 over 27 tissues, one-sided permutation p 0.0001
  (observed rho 0.734). The ceiling itself rises with amplitude (rho 0.408, p 0.020) but over a range of 0.011, two
  orders of magnitude below the spread of the couplings.
- **Brain (descriptive).** No relation, as in A104: rho -0.170 after correction (-0.115 observed).
- **Reading (registered).** Counting noise and host reliability do not produce the relation: a weak peripheral clock
  that was perfectly locked would be read as locked. The relation reflects own clock-plane variance that does not
  scale with the rhythm (A291 S2) or departures toward recent external cues (S3); A292 addresses the second.
- **Limits.** Non-count technical noise (dissection, composition, RNA quality) is part of psi and is not separated
  from biological own variation by this simulation; A259 found that matched non-clock genes carry much of a tissue's
  departure from body time, so part of psi is likely technical.
- Code: `exploratory/A291_attenuation/scripts/score_a291.py --jobs 4`.
<!-- src: paper/exploratory/A291_attenuation/output/{tests,tissues}.csv -->

## A292 (scored 2026-11-17; NOT BLIND for the A104 table, A259, A91, A270): weak clocks are not shown to track the hour of death beyond their host (H292.1 FAIL, opposite sign; uninformative by the registered power rule)

- Design as registered: GTEx fast deaths, 27 peripheral tissues, 244 donors; per tissue the states regressed on the
  leave-one-out host and on the hour of death rotated into the host frame; tracking weight w; H292.1 = rho(own
  amplitude, w) < 0 against 2000 permutations of the hour over donors. Aggregates only.
- **V292 (run and committed before scoring):** no tracking (S2): rejection 0.024 over 500 simulations (calibrated,
  conservative); tracking at w0 0.1: power 0.146; at w0 0.25: 0.422 (500 permutations inside each simulation).
  Power below 0.8: a failure is uninformative (registered).
- **H292.1 FAIL:** rho +0.425, one-sided p 0.971: the sign is the opposite of the prediction, tissues with stronger
  clocks carry more external-time weight beyond their host. w: median 0.014, range -0.256 to 0.189.
- **Companions (descriptive):** the clock's rho lies above 94 % of 100 matched sets (median -0.038) and 98 % of 100
  responsive sets (median 0.004).
- **Descriptive:** brain block (fast deaths, brain-only host): rho -0.066 (p 0.434), median w 0.070 (periphery
  0.014); ventilated donors, periphery (30 tissues): rho -0.081 (p 0.346), median w 0.013; ventilated brain not
  estimable (fewer than 2 brain tissues with >= 40 donors).
- **Reading (registered).** Failure at low power: uninformative for transient tracking (S3). With A291 PASS, own
  clock-plane variance that does not scale with the rhythm (S2) remains the reading of the amplitude-coupling line,
  by elimination and not shown; part of that variance is technical (A259).
- **Post hoc (not tested).** The registration assumed that host-error leakage into b is the same fraction for every
  tissue. With every tissue an anchor, leaving out a strong clock degrades its host more than leaving out a weak one,
  so leakage plausibly grows with amplitude; that would produce the positive sign, and the S2 simulation (host
  error drawn from each tissue's own posterior) may not capture it fully. The brain's higher w is on a different host
  and is not compared with the periphery.
- Code: `exploratory/A292_external_tracking/scripts/score_a292.py --stage calibrate | score --jobs 4` (calibration
  131 s; scoring 3059 s, mostly the 200 companion refits; run directly, not through `run.sh`).
<!-- src: paper/exploratory/A292_external_tracking/output/{calibration,tests,descriptive,tissues,sets}.csv -->

### A291 post hoc (2026-11-17; not registered, chosen after A291 and A292 were scored; descriptive, from the A104 table): tissues differ in the shared rhythm they carry, not in their own variation
- Coupling is c = S / (S + 2 psi) with S the shared clock-plane variance and 2 psi the tissue's own. Over the 27
  peripheral tissues, S spans 0.12 to 10.33 (85-fold; sd of log 1.02) and 2 psi spans 0.90 to 4.14 (5-fold; sd of log
  0.38). Own variance does not scale with the shared rhythm: Theil-Sen slope of log 2 psi on log S 0.07 [-0.08, 0.23];
  rho(own amplitude, S) 0.960, rho(own amplitude, 2 psi) 0.504.
- Reading (post hoc): in the periphery the amplitude-coupling relation is a difference in how much of the body's
  rhythm a tissue expresses above an own-variation floor of similar size in every tissue, not a difference in how
  much a tissue departs from the body's time. With A291 (the floor is not counting noise or host reliability), what
  the floor is (sample-level technical variation, or person-specific tissue clock variation of similar size
  everywhere) is not decided: A259 found matched non-clock genes carry 0.30 to 0.74 of a tissue's departure.
- Brain block (13 regions, own host; units not comparable with the periphery): S spans 2-fold, 2 psi 4-fold, and the
  two are negatively related (rho -0.434): between brain regions, own variation, not the shared rhythm, is what
  varies.
<!-- src: paper/results/figure_2/output/backbone/a104_gtex_tissues.csv (S = z2 V, own = 2 psi) -->

## A293 (scored 2026-11-17; NOT BLIND, see registration): slow deaths and ventilated living-brain donors keep the organ arrangement clock-specifically (H293.1, H293.2 PASS); cell classes keep their gain relative to endothelium across organs within a factor of about 1.4 (H293.4 PASS); a conserved stromal-vascular phase order is not shown by the registered sign test (H293.3 FAIL), though the pooled offset has the mouse sign

**Part 1 (level B; one run, 46 s; a smoke test on a 400-donor subsample with 20 sets was run and seen before it).**
Reproduction exact: Hardy 4 0.784, ventilated post-mortem 0.690, organ donors 0.161. Fisher-Lee rho of the pair phases
with fast deaths (200 matched and 200 responsive sets per tissue, same donors and pairs):

| group | n | pairs | clock rho (label p) | matched median / q95 | responsive median / q95 | verdict |
|---|---|---|---|---|---|---|
| slow deaths (Hardy 4) | 103 | 66 | 0.784 (0.0001) | 0.038 / 0.129 | 0.234 / 0.440 | **H293.1 PASS** |
| ventilated post-mortem (living brain) | 79 | 66 | 0.690 (0.001) | 0.025 / 0.108 | 0.145 / 0.318 | **H293.2 PASS** |
| Hardy 3 (descriptive) | 50 | 48 | 0.624 (0.002) | 0.033 / 0.183 | 0.327 / 0.573 | - |
| brain-dead organ donors (descriptive) | 79 | 65 | 0.161 (0.183) | 0.011 / 0.071 | 0.100 / 0.254 | - |

In both tested groups the clock is above all 200 matched and all 200 responsive sets. Non-clock genes carry little of
the fast-death arrangement (matched medians 0.01 to 0.04; responsive 0.10 to 0.33).
- **Reading (registered).** In slow deaths, whose internal time is not aligned to the hour of death (A109c), and in
  ventilated donors without meals or activity but with a living brain, the organs keep the fast-death arrangement on
  the clock genes beyond non-clock genes: the organ arrangement does not depend on alignment to the day or on
  behavioural cues. Brain-dead donors (descriptive, by decision): clock rho 0.161, not different from tissue-label
  permutation, at the 99th percentile of matched and 75th of responsive sets.

**Part 2 (level C; one run, seconds).** Mouse untimed offsets (timed in brackets), endothelial minus fibroblast: colon
+1.33 h (+1.30), aorta +1.07 h (+1.44); endothelial minus smooth muscle +0.96 (+0.65) and +0.93 (+1.16); smooth muscle
minus fibroblast +0.55 (+0.65) and +0.24 (+0.28): the same order in both tissues (positive = first class earlier, as the
untimed pairwise convention validated in A273).
- **H293.3 FAIL.** Human endothelial minus fibroblast, per atlas (readable pairs): 4 of 6 atlases with the mouse sign
  (kidney +1.50, lung +0.74, optic nerve +0.69, breast +0.24; anterior eye -0.39, Kong ileum -0.00), sign test p 0.344.
- **Descriptive.** Circular mean over the 6 atlases +0.46 h (atlas bootstrap 90 % 0.09 to 0.90); endothelial minus mural
  3 of 4 atlases, mean +0.96 h (0.13 to 1.84); mural minus fibroblast 1 of 2, +0.11 h.
- **Reading (registered).** No conserved cell-class phase order is shown by the registered test; at this resolution,
  timing is an organ property. The pooled human offset has the mouse sign and an interval excluding 0, about a third of
  the mouse size: a lead for a larger design, not a result.

**Part 3 (level C; one run, seconds).** Log gain relative to endothelium (natural log of the shared-rhythm amplitude),
random-effects over body atlases:

| family | atlases | pooled (SE) | tau | Q p | range over atlases |
|---|---|---|---|---|---|
| fibroblast_stromal | 13 | +0.228 (0.084) | 0.181 | 0.044 | -0.634 to +0.633 |
| mural | 7 | -0.079 (0.075) | 0.097 | 0.240 | -0.254 to +0.207 |
| parenchyma | 4 | -0.126 (0.123) | 0.145 | 0.205 | -0.438 to +0.310 |
| epithelial | 11 | -0.220 (0.126) | 0.305 | 0.001 | -1.234 to +0.345 |
| myeloid | 12 | -0.399 (0.117) | 0.302 | < 0.001 | -1.380 to +0.183 |

- **H293.4 PASS:** tau <= 0.35 in all five families (epithelial 0.305 and myeloid 0.302 close to the bound; Q shows
  heterogeneity beyond the bootstrap SEs in three families). Scale: GTEx tissues dominated by one class differ by log
  amplitude 0.76 to 0.98 (squamous: oesophagus mucosa against skin), 0.46 to 0.53 (adipocytes: visceral against
  subcutaneous and breast) and 0.25 to 0.32 (visceral smooth muscle: oesophagus against sigmoid colon) (A291 post hoc
  table; S ratios halved in log).
- **Reading (registered).** Relative gain is conserved across organs: fibroblasts about 1.26 times endothelium, mural
  cells and parenchyma close to it, epithelium about 0.8 and myeloid cells about 0.67, with between-organ spread
  (factor 1.1 to 1.36) below the organ differences of the same class in GTEx (factor 1.3 to 2.7). The organ
  differences of gain therefore come from composition and an organ-wide drive, not from cell classes responding
  differently in different organs. The organ-wide drive itself is not tested here (relative gains cancel it).
- **Limits.** Atlases have one or two compartments per family; the SE ignores the covariance with endothelium;
  "myeloid" mixes macrophages with monocytes and dendritic cells (A281: only macrophages carry a clock), which may
  explain part of its heterogeneity.
- Code: `exploratory/A293_organisation_levels/scripts/part1_death_types.py`, `part2_phase_order.py`,
  `part3_relative_gain.py`.
<!-- src: paper/exploratory/A293_organisation_levels/output/{part1/{groups,reproduction,sets},part2/{tests,atlases},part3/{families,atlases}}.csv -->

## A294 (scored 2026-11-17; NOT BLIND for A281 and A287): breast CD4 T coherence is not a clock (H294.1, H294.2 FAIL: "hormone-driven"); heart lymphocytes keep a coherence without the hormone-induced genes and above stress-responsive sets (H294.1, H294.2 PASS: "a weak lymphocyte clock signal"); post hoc, the core loop raises lymphocyte coherence in several atlases
Run `score_a294.py --jobs 4` (one scoring run, ledger ok, peak 4 GB). Reproduction: the 14-gene C and A281's matched
median and percentile reproduce in all nine pairs.
- **Power control:** thinned fibroblast-stromal rows pass the core loop in 100 % (breast) and 95 % (heart) of
  replicates: both informative.
- **Breast CD4 T (117 donors): FAIL / FAIL.** Core loop C 0.323, p 0.081; 14-gene C 0.384 below the responsive q95
  0.399 (92nd percentile). Registered reading: hormone-driven, not a lymphocyte clock.
- **Heart lymphocytes (36 donors): PASS / PASS.** Core loop C 0.597, p 0.014; 14-gene C 0.502 above the responsive
  q95 0.443 (96.5th percentile). Registered reading: a weak lymphocyte clock signal in the heart. Descriptive: the same
  C is inside the matched range (q95 0.584, 85th percentile), so not above random expression-matched genes.
- **Matched range for the slide (descriptive, 5th / 95th percentile of 200 sets):** the lymphocyte C lies inside it in
  seven of nine atlases, above it in breast (0.384 > 0.303), at the edge in heart; below the 5th percentile nowhere.
- **Post hoc (not registered; seen in the descriptive columns):** removing PER1, PER2 and NR1D1 raises the
  lymphocytes' coherence in most atlases: lung 0.24 -> 0.41 (p 0.075), IBDverse ileum -0.43 -> 0.40 (p 0.050), Kong
  colon 0.34 -> 0.54 (p 0.055), Kong ileum 0.06 -> 0.26 (p 0.24), optic nerve 0.06 -> 0.41 (p 0.096, power control
  0.20); IBDverse rectum stays negative (-0.40, p 0.85), Pelka colon ~0 (p 0.46). Six of nine lymphocyte compartments
  have core-loop p < 0.10. Reading, provisional: the hormone-induced clock genes move against the clock pattern in
  lymphocytes and may mask a weak core-loop coherence; A281 (14 genes) and A287 (14 genes, waveform-free) did not
  separate them. "No readable clock" holds on the full clock-gene set; on the core loop it is open.
- **Reading.** The two A281 passes are not one thing: breast is a hormonal pulse, heart a weak clock-pattern signal
  beyond a hormonal pulse but not beyond random genes. With the post hoc core-loop pattern, lymphocytes are better
  described as "no clock read on the full gene set; a weak core loop is possible", pending a blind test on fresh
  atlases (core loop, all lymphocyte compartments, power-controlled).
<!-- src: paper/exploratory/A294_lymphocyte_hormonal/output/{pairs,tests}.csv -->

## A295 (scored 2026-11-17; NOT BLIND for the organ map, A291, A293): a sample's composition does not move the phase of its clock state (H295.2 PASS, within +-0.13 h per SD); the lymphocyte dilution of gain is not detected (H295.1 FAIL); stromal-vascular content raises gain (descriptive)

- Design as registered: GTEx fast deaths, 27 peripheral tissues; per tissue, the clock state on its leave-one-out host
  with an interaction with a marker score (lymphocyte: CD3D, CD3E, CD2, CD247, MS4A1, CD79A, CD19; stromal-vascular:
  DCN, LUM, COL1A2, PECAM1, VWF, CDH5; standardised within tissue); donor bootstrap SE (200); random-effects pooling.
  One run (17 s).
- **H295.1 FAIL:** lymphocyte score, pooled gain change -0.010 per SD (SE 0.012), one-sided p 0.201; negative in 16
  of 27 tissues.
- **H295.2 PASS:** pooled phase change per SD +0.05 h (90 % -0.03 to +0.13) for the lymphocyte score and -0.05 h
  (-0.12 to +0.02) for the stromal-vascular score, both within +-0.5 h; between-tissue tau 0.16 h and 0.12 h.
- **Descriptive:** the stromal-vascular score raises gain, +0.041 per SD (SE 0.014; positive in 21 of 27 tissues).
- **Reading (registered).** H295.2 PASS with H295.1 FAIL: composition varies too little within tissues to move gain
  detectably through lymphocytes; the phase bound holds. Read with the stromal-vascular gain effect (not registered as
  a test): composition changes how loudly a sample expresses the body's time, not when, to within about 0.1 h per SD.
- Code: `exploratory/A295_composition/scripts/score_a295.py`.
<!-- src: paper/exploratory/A295_composition/output/{tests,tissues}.csv -->

## A296 (scored 2026-11-17; NOT BLIND for the organ offsets): the adrenal offset is a whole-clock property in GTEx (H296.1 PASS: +4.42 h with all clock genes, +4.61 h without PER1, PER2, NR1D1); in the baboon the outlier is the adrenal medulla (+5.9 h), not the cortex (+0.9 h)

- Design as registered: GTEx fast deaths, 27 peripheral tissues, timed on the hour of death; offsets relative to the
  mean over tissues with the full and the reduced waveform; donor bootstrap 200. One run (449 s; a 3-replicate smoke
  test was run first and its last lines, the baboon skin row, were seen).
- **H296.1 PASS.** Adrenal +4.42 h full, +4.61 h reduced, difference +0.19 h (90 % -0.58 to +1.19).
- **Secondary (all within 1 h):** lung +1.67 / +1.93 h (difference +0.26); skin sun-exposed +0.46 / +0.05 (-0.41);
  skin not sun-exposed +1.55 / +1.33 (-0.22). Tissues where dropping the three genes moves the offset by more than
  0.5 h: skeletal muscle (-1.36 h, 90 % -1.68 to -1.06), transverse colon (+0.79), pituitary (-0.61), sigmoid colon
  (+0.56).
- **Baboon (descriptive; timed on ZT, 63 tissues, median |offset| 1.47 h):** adrenal cortex +0.91 h (reduced +0.33),
  adrenal medulla +5.90 h (+4.94), lung +0.82 (+1.30), skin +1.97 (+1.94).
- **Reading (registered).** The adrenal sits about 4.5 h from the other organs on its core clock, not only on
  hormone-induced genes. Post hoc, from the baboon: the outlier there is the medulla, the part of the adrenal that is a
  sympathetic ganglion receiving direct preganglionic input, while the cortex sits with the other organs. GTEx adrenal
  samples mix cortex and medulla, so the human offset may come from the medulla; this fits the direct-neural-route
  hypothesis (as the pineal, A289) and is not tested.
- Code: `exploratory/A296_genewise_offsets/scripts/score_a296.py --n-boot 200`.
<!-- src: paper/exploratory/A296_genewise_offsets/output/{tests,offsets,genes,baboon}.csv -->

## A297 (scored 2026-11-17; NOT BLIND for the A111 shifts): removing the feeding rhythm keeps the mouse liver clock's amplitude far above what additive cues predict (H297.1 FAIL, above the interval): inputs set the liver clock's phase, the oscillator holds its amplitude

- Design as registered: one projector per dataset over its arms; per arm the timed complex amplitude of the liver
  state; within-ZT library bootstrap 2000. One run.
- **Calibration (inversion).** Manella daytime feeding against ad libitum: amplitude ratio 0.712 (90 % 0.659 to 0.772),
  shift 11.69 h; Xin day- against night-restricted: 0.438 (0.414 to 0.467), 11.29 h.
- **Prediction for removal, (1 + rho) / 2:** amplitude ratio 0.148 (0.117 to 0.181), shift +0.76 h (Manella); 0.288
  (0.274 to 0.301), +0.54 h (Xin).
- **Observed (Greenwell, arrhythmic against ad libitum):** amplitude ratio 0.724 (0.657 to 0.792), shift +1.57 h (0.94
  to 2.26). **H297.1 FAIL, above the interval** (also above the Xin prediction). Control, night-restricted against ad
  libitum: 1.037 (0.930 to 1.142), +0.74 h.
- **Reading (registered).** When the feeding rhythm is removed, the liver clock keeps about three quarters of its
  amplitude and moves 1.6 h, while additive cues predict a collapse to 15 to 30 %. An inverted cue moves it 11 to 12 h
  at 44 to 71 % of its amplitude. The liver clock behaves as a self-sustained oscillator: inputs set its phase, the
  oscillator sets its amplitude. Removing one cue (as in intensive care, A46) leaves it running near its phase;
  conflicting cues move it.
- **Limits.** Different laboratories, strains and protocols; liver only; amplitudes read in each dataset's own projector.
- Code: `exploratory/A297_feeding_model/scripts/score_a297.py`.
<!-- src: paper/exploratory/A297_feeding_model/output/{tests,arms}.csv -->

## A294b (scored 2026-11-17; NOT BLIND for the A104 table): the own variation of peripheral clocks is clock-specific, about 2.6 times that of expression-matched genes in every tissue (H294b.1 PASS); it grows only weakly with the shared rhythm (H294b.2: slope 0.22, neither the additive nor the individual-period signature)

- Design as registered: A104 rows and fit, 27 peripheral tissues; clock and 100 matched pseudo-clock sets. Reproduction
  exact (rho 0.734). One run through `run.sh` (1515 s, peak 2.8 GB, 4 workers).
- **H294b.1 PASS.** Median over tissues of log(2 psi_clock / median 2 psi_matched) 0.969 (tissue bootstrap 90 % 0.889 to
  1.077): the clock's own variance exceeds the matched median in all 27 tissues and the matched q95 in 10. Medians:
  clock 1.89, matched 0.71. The matched sets share little between tissues (median S 0.016 to 0.193, against the
  clock's 0.12 to 10.33).
- **H294b.2 (read because H294b.1 passed).** Theil-Sen slope of log excess own variance on log S_clock 0.223 (90 % 0.044 to
  0.332, 27 tissues): excludes both 0 and 1, so neither registered signature applies. Descriptive: slope of the clock's
  own variance on S 0.072; of the matched floor on S -0.025.
- **Reading (registered).** The floor is not generic: clock genes carry person-specific own variation beyond what
  random genes of their expression carry, of similar size in every peripheral tissue. It does not scale with the
  shared rhythm as individual clock periods would (slope 1); it grows weakly (excess proportional to about S^0.2, i.e.
  amplitude^0.45). Peripheral tissues are therefore not locked up to measurement noise: each carries a clock-specific
  departure from the body's time, of nearly constant size, and the amplitude-coupling relation is the gain over that
  departure.
- **Limits.** The matched floor is read through each set's own shared state, which is small; non-count technical
  variation specific to clock genes (none known) cannot be excluded. A259 bounded the departure at 1 to 3 h and found
  matched sets carry much of it on a different statistic (within-sample split-gene moment).
- Code: `exploratory/A294_floor/scripts/score_a294.py --jobs 4`.
<!-- src: paper/exploratory/A294_floor/output/{tests,tissues,sets}.csv -->

## A299 (scored 2026-11-17; NOT BLIND for the adrenal offset and the baboon medulla): the medulla content of a GTEx adrenal sample does not move its phase (H299.1 FAIL); cortex content does, towards the body's time (descriptive, 1.15 h per SD)

- Design as registered; one run (43 s). Primary donors Hardy 1-4: 75 adrenal samples with a leave-one-out host of 31
  tissues.
- **H299.1 FAIL.** Medulla score (TH, DBH, PNMT, CHGA, CHGB): phase change -0.25 h per SD (90 % -0.81 to +0.45),
  one-sided p 0.753; gain change -0.18 (-0.43 to +0.01).
- **Cortex score (specificity, descriptive):** phase change -1.15 h per SD (90 % -1.87 to -0.56), gain change +0.30
  (+0.15 to +0.49): samples with more cortex RNA sit closer to the body's time and express it more strongly. Fast
  deaths alone (60 samples): medulla -0.21 h (-0.82 to +0.65); cortex -0.86 h (-1.45 to -0.20), gain +0.29.
- **Reading (registered).** The human adrenal phase is not explained by its medulla content at this power; the adrenal
  offset stays an organ-level observation. Post hoc: the cortex effect is about nine times the largest composition
  effect on phase in the other tissues (A295, within 0.13 h per SD): in the adrenal, the part of a sample that is not
  cortex is phase-advanced, whatever it is (medulla, capsule, periadrenal fat); the medulla markers do not identify it.
- Code: `exploratory/A299_adrenal_medulla/scripts/score_a299.py`.
<!-- src: paper/exploratory/A299_adrenal_medulla/output/tests.csv -->

## A300 (scored 2026-11-17; NOT BLIND for the frontier pattern): frontier tissues do not change their clock relation to the body with the death season (H300.1 FAIL)

- Design as registered; one run (335 s). 19 tissues with >= 20 donors per arm (4 frontier: lung, both skin sites,
  oesophageal mucosa; stomach and transverse colon below the arm minimum); 66 summer and 53 winter donors.
- **H300.1 FAIL.** Median |summer - winter offset|: frontier 0.25 h, internal 0.30 h; T -0.05 h, permutation p 0.804
  (500). Largest offset changes: liver -0.69 h, breast +0.69 h, left ventricle +0.59 h; lung -0.43 h, sun-exposed skin
  +0.27 h.
- **Descriptive.** log z2 is lower in the summer arm in 18 of 19 tissues (-0.50 to -1.79; a difference of the arms'
  host scale, not interpretable per tissue); sun-exposed skin is the one tissue with the opposite sign (+1.00).
- **Reading (registered).** No season-dependent frontier effect at this power; the frontier pattern (lower coupling,
  more own variation) stays a post hoc observation.
- Code: `exploratory/A300_frontier_season/scripts/score_a300.py --n-perm 500 --jobs 4`.
<!-- src: paper/exploratory/A300_frontier_season/output/{tests,tissues}.csv -->

### A299 post hoc (2026-11-17; not registered, at the user's question whether the cortex effect is the converse of the medulla hypothesis with clearer markers): the cortex effect is not a medulla effect
- Same 75 samples. Marker consistency (mean pairwise correlation of log fractions) and spread (SD of the log score):
  medulla 0.879 and 1.68 (a fivefold range), cortex 0.704 and 0.50, fat 0.788 and 1.20, capsule / stroma (DCN, LUM,
  COL1A2, COL3A1) 0.731 and 0.68, blood 0.952 and 1.16. Medulla content is well measured and varies widely.
- Score correlations: cortex with medulla -0.13, with stroma -0.44, with fat +0.07; medulla with stroma +0.32.
- Phase change per SD (90 %), single scores: fat +0.35 (-0.38 to +0.86), stroma +0.35 (-0.34 to +0.97), blood -0.05.
  Joint models: cortex -1.29 h (-2.07 to -0.69) with medulla (medulla -0.48, -0.98 to +0.23); cortex -1.23 h (-2.08 to
  -0.51) with medulla, fat and stroma.
- Reading (post hoc). The cortex score is not an inverse medulla score (correlation -0.13) and its effect survives
  adjustment for medulla, fat and capsule content; medulla content, measured with consistent markers over a fivefold
  range, does not advance the phase (if anything the opposite sign). The cortex score varies little (SD 0.50 of the
  log) and its genes (STAR, CYP11A1, CYP11B1, CYP17A1, CYP21A2) are induced by ACTH: it may index steroidogenic activity
  more than cortex content. Under that reading, adrenal samples under stronger ACTH drive sit closer to the body's time
  and express it more strongly, a hormonal route consistent with the fetal cortex result (A264, A266); not tested.
- Code: `exploratory/A299_adrenal_medulla/scripts/posthoc_a299.py`.
<!-- src: paper/exploratory/A299_adrenal_medulla/output/posthoc/{scores,correlations,effects}.csv -->

## A301 (scored 2026-11-17; NOT BLIND for A299 and its post hoc): neither acute ACTH response (H301.1 FAIL) nor cortex identity (H301.2 FAIL) carries the adrenal phase effect; the steroidogenic enzyme score keeps it on its own (-1.10 h per SD)

- Design as registered; one run (seconds). 75 adrenal samples (Hardy 1-4); scores standardised within the adrenal.
- Score consistency (mean pairwise correlation; SD of the log score): ACTH / cAMP response 0.322, 0.82; cortex identity
  0.153, 0.35 (its genes hardly co-vary: a weak content measure); steroidogenic 0.704, 0.50. Correlations between the
  three scores 0.05 to 0.21.
- **H301.1 FAIL.** ACTH score, joint with identity: phase change -0.04 h per SD (90 % -0.57 to +0.63), p 0.483.
- **H301.2 FAIL.** Identity score, joint with ACTH: +0.55 h per SD (-0.04 to +1.03), outside the +-0.5 h bound and of
  the sign opposite to the steroidogenic effect.
- **Three-score model (reported):** steroidogenic -1.10 h per SD (-1.64 to -0.58), gain +0.33 (+0.17 to +0.49); ACTH
  -0.10 h; identity +0.61 h (+0.03 to +1.11). Fast deaths alone: steroidogenic -1.00 h (-1.57 to -0.40), identity
  +0.61 h.
- **Reading (registered).** Neither component reproduces the steroidogenic effect: it is not the acute ACTH / cAMP
  response, and not cortex content as measured by identity genes (which move the phase, if at all, the other way; their
  low consistency limits that reading). The predictor is the level of the steroidogenic enzymes themselves: adrenal
  samples with a higher steroidogenic program sit about 1 h per SD closer to the body's time and express it more
  strongly. Whether that program reflects chronic ACTH drive, zonation or sampling depth within the cortex is not
  decided.
- Code: `exploratory/A299_adrenal_medulla/scripts/score_a301.py`.
<!-- src: paper/exploratory/A299_adrenal_medulla/output/a301/{tests,effects,scores,correlations}.csv -->

### A294b correction note (2026-11-17; interpretation only, post hoc, no new statistic)
The registered reading of H294b.2 named slope 1 as the signature of individual clock periods. Corrected: for a
limit-cycle clock the angle by which an input turns it scales as 1 / amplitude, so the displacement of the state in the
clock plane is independent of amplitude. Three sources of a person's tissue departure then predict different slopes of
log departure variance on log shared variance: person-specific timing of cues (meals, activity, light) or additive
noise, slope 0; a fixed person-specific angle, slope 1; individual intrinsic periods (locking angle growing with
amplitude), slope 2. The observed slopes (0.07 for the clock's own variance, 0.22 for its excess over matched genes)
exclude intrinsic periods as the main source and do not separate person-specific cue history from noise; A302 tests
the first.

## A302 (scored 2026-11-17; NOT BLIND for A291, A292, A294b): the liver's personal departure is not shared preferentially with gut and adipose tissues (H302.1 FAIL)

- Design as registered; one run (1440 s, 4 workers; a check printing shared-donor counts only was run first). 21 pairs
  with >= 30 shared donors (stomach, kidney cortex and coronary artery below); G has 5 tissues.
- **H302.1 FAIL.** T = median r (G) - median r (N) = 0.189 - 0.154 = 0.034; random groups p 0.377; matched sets median T
  0.008, q95 0.093.
- **Descriptive.** Host-error-corrected |correlation| of the liver's departure with each tissue: 0.03 (visceral adipose)
  to 0.26 (thyroid); transverse colon 0.25, pituitary 0.24, tibial artery 0.22, subcutaneous adipose 0.20, skin 0.20;
  the median over the 100 matched sets' pairs 0.12 (set medians 0.05 to 0.20).
- **Reading (registered).** No group structure of departures at this power: the liver's departure is not shared more
  with meal-driven tissues than with others. The personal departure looks tissue-independent with respect to feeding;
  the person-specific cue-history reading of the constant departure (A294b correction note) is not supported here, and
  additive noise remains equally possible. The spread of T among matched sets (SD 0.048) means a group difference
  below about 0.1 could not be detected.
- Code: `exploratory/A302_departure_groups/scripts/score_a302.py --jobs 4`.
<!-- src: paper/exploratory/A302_departure_groups/output/{tests,pairs,sets}.csv -->

## A303 (scored 2026-11-18; descriptive, fixed criterion): the whole adrenal is not offset from the other organs in mouse (pooled -1.0 h) or sheep (-0.6 h); the human +4.4 h is not seen in either species

- Three plumbing failures before any statistic (minimum donors per age; bootstrap ids containing "#"; sheep strata), each
  fixed and committed; during the second, a debug call printed the young-mouse per-tissue phases. One complete run
  (about 25 min).
- **Mouse (Wolff 2023, timed on CT, adrenal against heart, kidney, lung, muscle).** Pooled over ages (age as stratum, 3
  libraries per CT): -1.00 h (90 % -1.87 to -0.39): not the human structure (criterion |offset| >= 2 h). Per age: young
  +0.03 h, aged -3.78 h, old -0.44 h, **without valid intervals**: each age has one library per CT, so resampling within
  CT cannot vary (a design flaw of the script; the "aged" row meets the criterion only through its degenerate interval
  and is not counted).
- **Sheep (untimed atlas, 21 animals with adrenal):** -0.57 h (90 % -1.64 to +0.41); median |offset| of the other
  non-brain tissues 0.45 h. Not the human structure.
- **Reading (registered).** A whole-gland offset near 0 in both species, with the baboon cortex at +0.9 h: the human
  +4.4 h has no counterpart in healthy animals and fits being specific to GTEx (death; the adrenal is the organ most
  engaged by dying) or to humans; the data cannot separate the two.
- Code: `exploratory/A296_genewise_offsets/scripts/score_a303.py`.
<!-- src: paper/exploratory/A296_genewise_offsets/output/a303.csv -->

## A304 (scored 2026-11-18; NOT BLIND for the A138 per-cancer correlations and the tumour-only coherence results): proliferating tumours express less of the body's rhythm (gain loss, H304.1 PASS, clock-specific) and do not drift further from it (departure rise not detected, powered); TP53 and luminal A show no departure rise (FAIL); ER descriptive

- **Design** (registered, staged): TCGA paired tumour-normal, 10 cancers, 517 pairs. Per arm, on the cancer's projector,
  gain G = |C| / V_N (the tumour's slope on its own normal sample: |g| x the normal's reliability) and departure R =
  (V_T - |C|^2 / V_N) / V_N^all; pooled Q_G, Q_R = log(case / reference) with weights n; 2000 within-cancer label
  permutations; 200 matched (decisive) and 200 responsive (reported) pseudo-clock sets. One run (20 s, 4 workers).
- **Power stage (A2, simulation at the real design):** the readouts separate. A pure gain loss (x 0.5) gives Q_G -0.509,
  Q_R -0.100 (it lowers R through the gain-times-host-error term); a pure departure rise (x 2) gives Q_G +0.087, Q_R
  +0.625. Type I error 0.044-0.068. Powered at the registered points (alpha 0.025): H304.1 gain 0.840, departure 1.000;
  H304.3 and H304.4 departure only (0.998, 0.854); H304.2 neither (0.120, 0.698).
- **H304.1 gain loss PASS (primary).** Q_G -0.462 (gain ratio 0.63), permutation p 0.0055, Holm 0.011; null mean 0.005
  (SD 0.178); below the matched q05 (-0.379; 1.5th percentile) and at the 0.5th percentile of the responsive sets.
  Lower G in the high-proliferation arm in 7 of 10 cancers (BLCA, BRCA, KIRC, LIHC, PRAD, THCA, UCEC); higher in COAD,
  KIRP, LUAD.
- **H304.1 departure rise FAIL (primary, powered).** Q_R -0.190, p 0.944: the departure is if anything lower, as a pure
  gain loss predicts (-0.100 in the simulation); a doubling would have been detected with power 1.000.
- **Host check:** the two arms' normal samples have the same clock-plane variance (Q_N -0.089, p 0.40): the gain ratio is
  not a difference of host reliability.
- **Purity:** the high-proliferation arm has more stromal / immune marker expression (lower purity) in 8 of 10 cancers,
  which by A295 would raise gain, not lower it. Sensitivity with cells = cancer x purity half (18 cells, 498 pairs): Q_G
  -0.266, p 0.037; Q_R -0.232, p 0.96 (smaller cells inflate both arms' |C| bias, compressing the ratio).
- **Secondary departure readouts (Holm over 2):** H304.3 TP53 altered vs intact (5 cancers, 299 pairs) Q_R -0.522, p 0.995,
  Holm 1.0: FAIL, the departure is lower in TP53-altered tumours, at the matched q05 (5.5th percentile, not
  clock-specific). H304.4 luminal A vs other (BRCA, 105 pairs) Q_R -0.216, p 0.86: FAIL.
- **Descriptive (not powered):** gain, TP53 -0.049 (p 0.33); luminal A +0.924 (higher gain in luminal A, as A219; inside
  the matched range, 93rd percentile); ER negative -0.451 (p 0.21; matched 8th percentile). ER departure +0.135 (p 0.24).
  Pooled offset of the high- vs low-proliferation tumours +0.60 h.
- **Reading (registered).** "Gain loss detected, departure rise not": in proliferating tumours the body's time is still
  there, expressed at about two thirds of the gain; the tumour clock does not drift further from the patient's time. This
  is the dilution / weakening pattern of lymphocytes in normal tissue, not a desynchronised tumour clock. Bulk cannot say
  whether malignant cells carry a weaker clock or a smaller share of the tumour RNA comes from clock-carrying cells (A212:
  damping and dephasing are not separable in bulk); the purity imbalance argues against more stroma as the cause.
- **A3 (descriptive).** Tumour coupling to the body's time c_T = corr^2 / rel (rel the GTEx coupling of the matching
  tissue; 90 % patient bootstrap): LUAD 0.54 [0.33, 0.83], LIHC 0.38 [0.22, 0.58], KIRC 0.27 [0.11, 0.50], PRAD 0.26
  [0.06, 0.57], BRCA 0.18 [0.07, 0.33], KIRP 0.12 [0.03, 0.39], THCA 0.11 [0.03, 0.25], COAD 0.06 [0.01, 0.21]; BLCA and
  UCEC without a matching tissue (corr 0.52 and 0.08). Gain on body time |g| = G / rel: LUAD 1.21, LIHC 0.79, KIRC 0.64,
  PRAD 0.56, BRCA 0.48, KIRP 0.47, THCA 0.39, COAD 0.28. Per-cancer correlations reproduce A138. The ranking leans on
  GTEx reliabilities from deaths on a multi-tissue host (lung 0.44 is low in GTEx, which lifts LUAD above 1), and the
  normal is adjacent tissue at surgery; intervals overlap except at the ends (LUAD and LIHC against COAD and THCA).
- **Limits.** One normal per patient: G is attenuated by its reliability, equal across arms within a cancer (host check)
  but not between cancers. Proliferation split at the median within cancer; arms of 9 to 56 pairs.
- Code: `exploratory/A304_cancer_gain_departure/scripts/{power_a304.py, score_a304.py} --jobs 4`.
<!-- src: paper/exploratory/A304_cancer_gain_departure/output/power/{calibration,design,power,mde}.csv; paper/exploratory/A304_cancer_gain_departure/output/score/{tests,readouts,per_cancer,companions,sensitivity_purity,purity_imbalance,a3_per_cancer}.csv -->

## A305 power note (2026-11-18; placebo draws only, no observed statistic computed): malignant-stroma pairwise sharing is not powered; A305 closed at the power stage

- 34 single-cell tumour cohorts (A227 list); only 2 have a (malignant, stromal-vascular) pair passing the shared-donor
  rule (min_stratum 3, min_df 15): GSE178341 colorectal (malignant-endothelial, 24 donors) and GSE183904 gastric
  (malignant-endothelial and malignant-fibroblast, 18 donors each): 46 patients, 3 pairs. The other cohorts have 3 to 11
  patients with both compartments.
- 5000 placebo draws: smallest attainable p_cal 0.0036 (matched), 0.0062 (responsive); D = T_MS minus matched median:
  placebo mean 0.0015, SD 0.0616, q95 0.106; **MDE 0.156 > 0.10**: rule (ii) fails, H305.1 not scored.
- Reading (registered): no reading of malignant clock sharing with the tumour stroma; the A149 covariation stays
  unresolved for clock specificity. A powered design needs single-cell tumour cohorts with >= 15 patients carrying both
  malignant and stromal-vascular cells in one stratum.
- Code: `exploratory/A305_malignant_stroma_pairwise/scripts/power_a305.py` (58 s; run by a Sonnet subagent; numbers
  checked against `power.csv` and `design.csv`).
<!-- src: paper/exploratory/A305_malignant_stroma_pairwise/output/power/{power,design,design_pairs}.csv -->

## A306 (scored 2026-11-18; NOT BLIND for d* and offsets of these conditions, A164 / A233 / A235): T2D, hypertension and smoking change neither the gain nor the departure of their organ relative to the rest of the body (5 powered readouts, all FAIL); smoking lung gain lower, not powered (descriptive)

- Design as registered (A233 matching, A164 anchors, 1000 within-set permutations shared across tissues; vent and fast);
  Delta = module minus companion (tibial nerve, both skins) of the per-tissue log z2 ratio (gain) or log psi ratio
  (departure), inverse-null-variance pooling. Power stage 2218 s, scoring 40 s.

| condition | readout | Delta (module / companion) | p | Holm | MDE_80 | verdict |
|---|---|---|---|---|---|---|
| T2D | gain | -0.029 (+0.021 / +0.050) | 0.83 | 1.0 | 0.327 | not detected |
| T2D | departure | -0.139 (-0.050 / +0.090) | 0.24 | 0.98 | 0.288 | not detected |
| hypertension | gain | +0.213 (+0.116 / -0.097) | 0.13 | 0.63 | 0.346 | not detected |
| hypertension | departure | -0.012 (+0.033 / +0.045) | 0.92 | 1.0 | 0.264 | not detected |
| smoking | gain | -0.514 (-0.455 / +0.060) | 0.068 | - | 0.700 | descriptive (not powered) |
| smoking | departure | -0.029 (-0.103 / -0.073) | 0.86 | 1.0 | 0.418 | not detected |

- **Arm-prior sensitivity** (observed only): every Delta within 0.011 of the pooled-prior value (e.g. smoking gain -0.516,
  T2D departure -0.141).
- **Descriptive per tissue:** smoking lung log z2 ratio -0.647 (vent), -0.292 (fast): lower in smokers in both cohorts;
  T2D liver -0.629 (vent, 15 sets), +0.080 (fast, 25 sets). Vent coronary artery +5.19 (47 sets; null SD 0.91, low
  weight) is an unstable estimate.
- **Reading (registered).** With d* known unchanged (A164, A233), no condition scales the organ's clock up or down:
  in T2D (metabolic organs), hypertension (arteries, heart) and smoking (lung, departure only), neither the gain nor the
  departure moves by a factor of about 1.3-1.5 relative to the companions (MDE 0.26-0.42 on the log scale). With the
  earlier offsets (equal within about 1 h) and couplings, chronic disease leaves the organisation of the clock intact at
  the organ level in GTEx. A lower lung gain in smokers (about 0.6 of non-smokers on z2) is the one suggestive signal;
  not powered and not tested.
- **Limits.** GTEx deaths; comorbidity not removed (hypertension with T2D, SMD 0.64-0.66); BMI imbalance in T2D fast
  (0.70); the companion tissues are not disease-free (skin and nerve in T2D).
- Code: `exploratory/A306_disease_gain_departure/scripts/compute_a306.py --stage power | score --jobs 4`.
<!-- src: paper/exploratory/A306_disease_gain_departure/output/{power,power_tissues,balance_summary,score_tests,score_tissues,score_sensitivity}.csv -->

## A307 (scored 2026-11-18; NOT BLIND for A243-A269): the fetal signal-induced clock genes keep the waveform's phase relations to the core genes in adrenal cortex and lung (H307.2 PASS x2); they also follow glucocorticoid-target activity clock-specifically in lung (H307.1 PASS), not beyond matched genes in cortex; liver resolves neither

- Design as registered: rows and adjustment of A264 (cortex, 34 fetuses, 6 cohorts) and A269 (lung 38 / 5, liver 35 / 3);
  L = mean r^2 of PER1, PER2, NR1D1 with the GC score (FKBP5, TSC22D3, ZBTB16, KLF9); W = mean over the 33 (signal-induced,
  core) pairs of r cos(phase difference); 10,000 permutations, 200 matched sets. Power stage (null only) 17 s; scoring 22 s.

| organ | L (GC link) | p (Holm) | matched median / q95 (pct) | W (waveform) | p (Holm) | matched median / q95 (pct) |
|---|---|---|---|---|---|---|
| adrenal cortex | 0.180 | 0.0020 (0.0040) | 0.094 / 0.192 (92.5) | **0.083** | 0.0003 (0.0009) | 0.001 / 0.072 (98.5) |
| lung | **0.204** | 0.0001 (0.0003) | 0.053 / 0.138 (100) | **0.076** | 0.0003 (0.0009) | 0.000 / 0.066 (97.0) |
| liver | 0.061 | 0.117 (0.117) | 0.131 / 0.340 (17.5) | 0.037 | 0.042 (0.042) | 0.000 / 0.062 (81.5) |

- MDE (power note): L 0.096-0.122, W 0.051-0.058.
- **Verdicts.** H307.1: lung PASS; cortex significant but not above the matched q95 (expression-matched genes follow the
  GC score almost as much, median r^2 0.094): not clock-specific; liver fail. H307.2: cortex PASS, lung PASS; liver
  p 0.042 (Holm 0.042) below the matched q95: not clock-specific.
- **Signed correlation with the GC score (descriptive).** PER1 +0.51 / +0.56 / +0.25 (cortex / lung / liver), PER2 +0.49
  / +0.51 / +0.34, NR1D1 +0.21 / +0.20 / -0.04; DBP negative in all three (-0.09, -0.38, -0.39); several core genes are
  also positive (NR1D2 +0.55 / +0.30 / +0.35, CRY2 +0.44 / +0.43 / +0.57, CLOCK +0.44 / +0.20 / +0.57). The GC score is
  a broad axis of the fetal rows, not a PER1 / PER2 / NR1D1-specific one; NR1D1 does not load with the GC targets.
- **Reading (registered).** Cortex: "H307.2 pass, H307.1 fail": PER1, PER2 and NR1D1 covary with the core genes as the
  clock's phase relations predict, beyond expression-matched genes, and their link to GC-target activity is not
  clock-specific: loop-consistent rather than a GC-induced state. Lung: both pass: a GC-driven component on top of
  loop-consistent covariation, not separable. Liver: not resolved.
- **Relation to A266 / A269 (post hoc, interpretation).** A266 and A269 found the core genes alone (without PER1, PER2,
  NR1D1) not significant or not clock-specific; here the cross block between the signal-induced and the core genes has
  the waveform's sign pattern. The two are compatible: the core genes' covariation is weak on its own but aligned with
  PER1 / PER2 / NR1D1 as in an adult clock. The simple reading of A266 ("a hormone-induced response of the
  signal-induced genes alone") is therefore weakened for the cortex and the lung: the signal-induced genes do not move
  alone; whether the shared state is a running loop or an input that moves the whole network in its clock pattern is
  not decided without times (or spliced / unspliced counts).
- **Limits.** Not blind; no adult power (fails read against the MDE only); the cos(phase difference) template uses adult
  GTEx phases; Pearson r over 34-38 fetuses after cohort and age adjustment.
- Code: `exploratory/A307_fetal_gc_pattern/scripts/score_a307.py --stage power | score`.
<!-- src: paper/exploratory/A307_fetal_gc_pattern/output/{power,tests,genes,companions}.csv -->

## A308 (scored 2026-11-18; NOT BLIND for A183 and the archive psoriasis results): the psoriatic lesion delay is beyond the lesion-wide rotation of random expression-matched genes once their estimation noise is removed (H308 PASS, q 0.016); A183's 5.5th percentile reflects noisy companions, and its q05 rule never rejects under the complete null (V308 0.000)

- One run (14 s, 2000 per-person swap draws, the same swaps for every set). Draw 0 reproduces A183 (clock -1.453 h;
  matched q05 -1.560, percentile 0.055; responsive q05 -1.012).
- **V308 (A285 form):** with no lesion effect on any gene, the clock falls below the matched q05 in 0 of 2000 draws
  (responsive 0.024): in this design the percentile rule is conservative, not calibrated. The clock's own offset noise is
  0.33 h; that of the matched sets is mean variance 0.76 h^2 (SD about 0.87 h).
- **H308 PASS.** Matched sets: observed offsets mean -0.04 h, variance 1.08 h^2, of which 0.76 is their own estimation
  noise (noise share 0.71); true lesion-wide rotation of random genes tau 0.56 h. The clock (-1.45 h, noise 0.33 h) lies
  at q = 0.016 of noise-free random sets observed with the clock's noise (90 % over resampled sets 0.0007 to 0.046).
- **Responsive family (reported):** the spread is all noise (noise share 1.23, tau 0); q < 0.0001.
- **Reading (registered).** The matched sets' wide spread in A183 was mostly their own noise (few, weakly covarying
  genes in 81 persons); the lesion does rotate random genes, but by about 0.6 h SD around 0, and the psoriatic delay of
  1.45 h exceeds it. The delay is clock-specific against both companion families. Limits: method-of-moments tau with a
  normal approximation (upper interval bound 0.046, close to the line); not blind; composition not addressed.
- Code: `exploratory/A308_psoriasis_calibration/scripts/score_a308.py`.
<!-- src: paper/exploratory/A308_psoriasis_calibration/output/{tests,sets}.csv -->

## A309 (scored 2026-11-18; NOT BLIND for A236 and the human fetal results): in the late mouse embryo (E15-E18.75, embryos at 00, 06, 12, 18 h) hepatocytes carry a clock-specific covariation between embryos (H309.1 liver PASS) whose signal-induced genes keep the waveform pattern (H309.2 W PASS); ventricular cardiomyocytes, nephron progenitors, lung progenitors and adrenocortical cells do not pass (heart: three quarters of the adult size excluded)

- Design as registered: Qiu 2024 late-window genome-wide pseudobulks (>= 50 cells), batch (embryonic day) residualised,
  10,000 per-gene permutations within batch, 200 matched and 200 responsive sets. Power stage 421 s; two plumbing
  failures before the reported run (CIART below 1 count in a readout's matching table, no statistic computed; then an
  import path of the A307 functions, after the H309.1 table had printed); the final run reproduces H309.1 exactly (64 s).

| readout (cell type) | embryos | C | p | Holm | matched q95 (pct) | responsive q95 (pct) | adult pass 1 / 0.75 / 0.5 | verdict |
|---|---|---|---|---|---|---|---|---|
| liver (hepatocytes) | 16 | 0.642 | 0.0023 | 0.0115 | 0.425 (99.0) | 0.434 (99.5) | 1.000 / 1.000 / 0.945 | **PASS** |
| adrenal (adrenocortical) | 13 | 0.407 | 0.068 | 0.270 | 0.417 (94.0) | 0.377 (96.0) | 1.000 / 0.975 / 0.645 | fail |
| heart (ventricular cardiomyocytes) | 12 | 0.239 | 0.157 | 0.471 | 0.411 (78.0) | 0.351 (85.5) | 1.000 / 0.970 / 0.615 | fail |
| lung (lung progenitors) | 16 | 0.117 | 0.399 | 0.797 | 0.457 (55.5) | 0.526 (56.0) | 1.000 / 1.000 / 0.960 | fail |
| kidney (nephron progenitors) | 16 | -0.516 | 0.946 | 0.946 | 0.416 (7.0) | 0.423 (5.5) | 1.000 / 1.000 / 0.930 | fail |

- **H309.2 (liver only):** L (GC link of PER1, PER2, NR1D1) 0.347, p 0.0016, below the matched q95 0.491 (73.5th
  percentile): not clock-specific; W (waveform cross-agreement) 0.105, p 0.0005, above the matched q95 0.074 (98.0th):
  PASS. Signed correlation with the GC score: PER1 +0.77, NR1D1 +0.67, NR1D2 +0.79, CRY2 +0.80, TEF +0.74, CLOCK -0.52,
  CRY1 -0.45, PER2 +0.03: the GC score lies along the clock's own phase axis in these embryos.
- **All ages (descriptive, quadratic E-day, global permutations, no companions):** liver 54 embryos (E9-E18.75) C 0.526,
  p 0.0006; lung progenitors 36 (E10.75-) 0.437, p 0.013; adrenocortical 15 0.404, p 0.059; heart 12 0.300, p 0.096;
  nephron progenitors 26 -0.012.
- **Reading (registered).** The heart branch: fail, with an adult pass rate 0.970 at 0.75 and 0.615 at 0.5: a ventricular
  clock synchronised across embryos at three quarters of the adult size is excluded at E15-E18.75, a smaller one is not.
  The registered "organ split replicates" branch required adrenal / liver / lung to pass; only the liver does. Read
  together: at the stage where the mouse heart clock is reported to function (explants, Umemura 2017), cardiomyocytes
  of embryos sampled across the day show no clock-shaped covariation, as in the human fetal heart; the kidney is
  negative in both species; the liver is the organ synchronised across individuals in both species. Lung (epithelial
  progenitors here, whole lung in human) and adrenal cortex (13 embryos, p 0.068, C above both companions' q95) are
  not replicated in the late window; lung progenitors covary over all ages (descriptive). In the liver, as in the
  human adrenal cortex (A307), PER1, PER2 and NR1D1 move with the core genes in the clock's phase pattern beyond matched
  genes, and the GC targets lie on the same axis without being clock-specific: a synchronised, clock-patterned state
  of which GC-target activity is part.
- **Limits.** 12-16 embryos per cell type; heart and adrenal powered for three quarters, not half, of the adult size;
  stand-ins from adult pools (Wolff all ages; Manella liver shallower than the fetal hepatocytes, no thinning); cell
  types are not the human organs (progenitor and parenchymal classes vs whole-organ or class pseudobulks); not blind.
- Code: `exploratory/A309_mouse_embryo_organs/scripts/{power_a309.py, score_a309.py}`.
<!-- src: paper/exploratory/A309_mouse_embryo_organs/output/{power_design,power_passrates,l1,gc_pattern,gc_genes,all_ages}.csv -->

## Scout (2026-11-18; metadata only, no data downloaded, no statistic): spliced / unspliced counts and unseen cohorts for the fetal section
Two Sonnet subagents (GEO / ENA / BioStudies / ENCODE / cellxgene metadata); key numbers spot-checked by the lead
(ENA filereports of PRJEB58854 and PRJEB26969).
- **Processed spliced / unspliced counts:** none found in any cohort used (A189-A269); RAW.tar contents and the Cao 2020
  loom layers (12 GB, named "gene_count") not opened.
- **Raw reads open:** Cardoso-Moreira E-MTAB-6814 (bulk polyA, ENA PRJEB26969: prenatal heart about 38, liver 37,
  kidney 30 individuals; about 82 / 73 / 53 GB), Kameneva GSE147821 (10x, about 11 fetal adrenal donors, about 415 GB),
  Wang 2022 GSE167860 (10x, 10 adrenal donors, about 172 GB), Cheng 2022 GSE193672 (about 6 donors, 3-8 wpf, about
  115 GB), Dong 2020 GSE137804 (about 6 fetal samples among 670 GB), Knight-Schrijver GSE216019 (10x, 6-7 hearts, 231 GB),
  Cui 2019 (plate, 20 hearts, 1.5 TB), Hochane GSE114530 (5 kidneys, 59 GB), ENCODE fetal (small sets). R1 of the 10x
  sets is listed in SRA as original files, availability without requester-pays not tested. **Closed:** Cao 2020
  (dbGaP), Sountoulidis 2023 lung (federated EGA), Kildisiute (EGA). Lung cannot reach 15 fetuses with open raw reads.
- **Unseen cohorts (blind replication):** adrenal: del Valle 2023, E-MTAB-12492 / PRJEB58854, bulk RNA-seq of 32 whole
  fetal adrenals (8 per stage at 7, 8, 9, 11.5 wpc, 4 XX / 4 XY), paired-end 43 bp, 46 runs (with 14 non-adrenal
  controls), 49 GB FASTQ, no count matrix deposited. Liver: Popescu 2019 E-MTAB-7407 (14 fetuses, CD45+ / CD45- sorted,
  probable overlap with Suo 2022). Heart, lung: no unseen cohort with >= 8 fetuses found. Kidney: <= 9 (low depth).
<!-- src: subagent reports of 2026-11-18 (session transcript); ENA filereports PRJEB58854 (46 runs, 49.3 GB) and PRJEB26969 (313 runs) -->

## A310 Stage 2a (scored 2026-11-18; NOT BLIND, descriptive fit): the driven population-of-oscillators model is INADEQUATE by the registered criterion (blocks A and B fail; C fails); Stage 2b not run

- Design as registered; block A intervals from 200 donor resamples (`boot_a310.py`, replicate 0 reproduces A294b's
  own variances exactly, 6 min). Fit in one run after a plumbing failure: the first run printed the block A-C tables
  and stopped in the identifiability profiles (overflow at the edge of a parameter range; guarded); the rerun is
  deterministic and printed the same tables. Deviation (same model, not a change of it): the liver's locked state
  is solved exactly (the equation is rotation-covariant, w* = R(|F|) e^{i arg F}; checked against integration)
  instead of integrated.
- **Block A FAIL.** No organ variant meets both rules. The departure E_o (2 psi over the matched floor, 0.37 to 3.26)
  falls inside the tissue's 90 % interval in 33 % (dispersion), 37 % (coupling) and 37 % (drive) of the 27 tissues,
  against 80 % required. The model's slope of log E on log S is inside the A294b interval for dispersion (0.091) and
  coupling (0.084), not for drive (0.000). Residual SD of log E around the best fit 0.42, against a median bootstrap SE
  of 0.19: tissues differ in departure size beyond what their synchrony predicts. Each tissue would need its own
  noise intensity, one parameter per observation.
- **Block B FAIL.** 76.6 % of the 47 atlas relative gains are within 1.645 SE of their class value (78.7 % around the
  A293 random-effects values), against 80 %; myeloid 0.67, epithelial 0.73. Phases: the fitted endothelial lead is
  0.67 h; the 12 class-pair offsets scatter around it with RMS 0.91 h (1.10 h around a shared phase), against 0.6 h.
- **Block C FAIL** (4 of 8 inside their intervals; 6 required). The population model reproduces the removal
  amplitude (predicted 0.715, observed 0.724) and both inversion shifts. It cannot give the two laboratories'
  inversion amplitudes (one value, 0.525, against 0.712 and 0.438) or the 1.6 h removal shift (predicted 0.18 h).
  Linear damped alternative: 3 of 8 (removal amplitude 0.297); weighted SSE 82.8 (population) against 369.3 (linear).
  The best liver fit needs no coupling (K -> 0): the amplitude is held by the cells' fixed amplitude (saturation), not
  by synchronisation.
- **Identifiability.** The best-fitting block A variant is degenerate (drive -> 0, compensated by the noise scale), and
  the liver coupling is unbounded below. No parameter value is quoted.
- **Alternative (b), composition alone (descriptive).** With classes spanning a factor 1.87 in gain (A293 Part 3),
  composition alone needs less than half the clock-gene mRNA from clock-carrying cells only in whole blood (bound
  0.20) and testis (0.31). Skin is next (0.60 and 0.76). The bound therefore does not require an organ synchrony factor
  in any tissue whose dominant cells carry a clock.
- **Global.** 14 parameters for 95 observations (rule met); the verdict is set by the blocks: **INADEQUATE**.
- **Reading (registered).** The model is not worth it as a unifying fit. With one departure intensity and organ
  synchrony as input, it does not reproduce the size of each tissue's personal departure or the cell-class gains
  within their errors, and it adds nothing on organ gain. Which parts survive as qualitative statements, none of them
  a fit:
  - a fixed-size departure turns weak clocks more (A294b correction note, Abraham 2010);
  - the liver keeps its amplitude without its feeding drive because single-cell amplitude saturates (A297);
  - the arithmetic bound leaves composition able to explain organ gain except in blood and testis.
  Stage 2b (orientation of the departure) is not run: the registration conditions it on blocks A and B.
- Code: `exploratory/A310_clock_model/scripts/{model.py, boot_a310.py, fit_a310.py}`.
<!-- src: paper/exploratory/A310_clock_model/output/{boot,blockA_tissues,blockA_variants,blockB_gains,blockB_phases,blockC,profiles,composition_bound,verdict}.csv -->

## A311 (scored 2026-11-18; NOT BLIND for the A104 table and A310; powered): the personal departure of GTEx peripheral tissues is not a phase departure (H311.1 FAIL, share -0.03), with no significant gradient with gain (H311.2 FAIL, rho 0.30, p 0.064): the driven-Hopf model's sharpest prediction is not supported

- Design as registered; one scoring run (2723 s, 4 workers, peak 2.7 GB); clock and 100 matched sets, no errors.
  Primary share = excess tilt / A104 2 psi (deviation note); the script's ML-based share gives the same verdicts
  (median -0.029, 90 % -0.064 to -0.007; rho 0.288, p 0.079).
- **Power (registered, before scoring):** a clock-specific departure 30 % more tangential than radial would give a
  median share near 0.16 (detection 1.00; type I 0.05). A cycle above its bifurcation predicts more than 0.5 of the
  clock-specific departure (model grid).
- **H311.1 FAIL.** Median share over the 27 tissues -0.028 (tissue bootstrap 90 % -0.064 to -0.007); share > 0 in 9 of
  27 tissues; range -0.232 to +0.180. Outside the null at |z| > 2: subcutaneous adipose (-0.191) and lung (+0.180).
  Share of the clock-specific departure (over E_o): median -0.045, above 0.5 in no tissue.
- **H311.2 FAIL.** Spearman rho(share, S) 0.296, one-sided p 0.064.
- **Matched sets.** Their median share is -0.213 (q95 of the set medians 0.016); the clock lies at their 90th
  percentile. Random genes read strongly radial at the same design, so the isotropic simulation does not reproduce
  real non-clock variation (a sample-level scale component of expression is the likely source; not examined). The
  absolute share may therefore be biased downwards for every gene set. Relative to the matched genes the clock is
  not significantly more tangential (90th percentile).
- **Reading (registered).** The upper bound is below 0, but by 0.007 and with a median of -0.03. Of the registered
  readings, "isotropic, no evidence of limit-cycle dynamics" fits better than "radial dominates". The personal
  departure is not the phase-dominated response of a cycle above its bifurcation; a 30 % tangential excess would have
  been detected. On the model grid, a share near 0 places a compartment far below the bifurcation, which contradicts
  the self-sustained liver (A297). The reading is therefore that the departure GTEx measures is not a dynamical
  perturbation of each tissue's cycle. It is variation added to the state (sample-level technical or compositional
  variation, or state at death), consistent with A302 (not shared by meal-driven tissues). The Hopf model stays an
  interpretation of gain and lag; its departure prediction fails in GTEx.
- **Limits.** Deaths: a person's tissue state at death may not be a stationary fluctuation. The calibration is not
  matched to real non-clock variation (matched sets at -0.21). Fast deaths only.
- Code: `exploratory/A311_hopf/scripts/{hopf.py, power_a311.py, score_a311.py}`; `clocksync.orientation` (f9ea9f1).
<!-- src: paper/exploratory/A311_hopf/output/{tests,tissues,tissues_primary,sets,power}.csv -->

## A312 (scored 2026-11-18; NOT BLIND for the A104 table, A310, A311; powered): lower-gain tissues follow the persons' shared clock amplitude slightly more steeply (H312.1 PASS, rho -0.34, p 0.037), but every tissue follows it nearly in proportion (elasticities 0.83 to 1.77, low/high-gain ratio 1.11): the gradient has the Hopf sign, at a small fraction of the predicted size

- Design as registered; power stage (leverage 0.783, power 0.984, gate passed), then one scoring run (171 s).
- **Persons differ strongly in shared amplitude.** The two half-host radii correlate at 0.75 to 0.86 over donors in
  every tissue (power note): in GTEx deaths, how loudly a person's tissues express the shared rhythm is a strong
  person-level property.
- **H312.1 PASS.** Spearman rho(elasticity, S) -0.341, one-sided permutation p 0.037 (27 tissues).
- **Size (descriptive).** Elasticity relative to the other tissues' composite: median 1.19, range 0.83 (aorta) to 1.77
  (adrenal). The 90 % interval lies above 1 in 11 tissues and below 1 in none. Median of the 9 lowest-S tissues 1.19,
  of the 9 highest 1.07, ratio 1.11. The registered alternative was 2.0 against 1.0. On the model grid, mu -1 against
  +1 gives a ratio of 3.4, and a ratio of 1.11 corresponds to nearly the same distance from the bifurcation.
- **Reading (registered: PASS).** The sign is the Hopf model's: tissues with less gain follow the persons' shared
  amplitude more. The size says that most person-level amplitude variation scales every tissue together, nearly in
  proportion. That is a linear or common-factor co-scaling (systemic amplitude, or donor-level technical scale in
  deaths), not receivers at very different distances from their bifurcation. Organ gain differences are therefore at
  most weakly a matter of bifurcation distance, and composition or linear gains remain the main account (A310
  alternative (b) bound). No companion gene set is possible for this statistic (pseudo-clock sets have no host
  amplitude); a donor-level technical scale that acts more on noisy low-gain tissues is not excluded.
- Code: `exploratory/A312_elasticity/scripts/score_a312.py`; `clocksync.elasticity` (eed9786).
<!-- src: paper/exploratory/A312_elasticity/output/{power,power_tissues,tissues,tests}.csv -->

## A314 (scored 2026-11-18; descriptive, post hoc split): ventilated donors with a living brain and brain-dead organ donors have lost alignment to the hour of death alike (R 0.158 and 0.174; 0.180 at equal n)

- Design as registered (A27 procedure per group; one run, minutes).
- **Ventilated post-mortem (living brain), 79 donors:** R 0.158 (permutation p 0.142; null median 0.096, q95 0.198).
- **Brain-dead organ donors, 301 donors:** R 0.174 (p < 0.001; null median 0.046, q95 0.103); at the post-mortem n
  (20 draws of 79): median R 0.180 (range 0.091 to 0.330), median p 0.076, null median 0.093.
- Pooled Hardy 0 for reference (A27): R 0.142, 380 donors. Fast deaths (A27): 0.631, 244 donors.
- **Reading (descriptive).** At equal n the two ventilated groups have the same weak alignment to the hour of death,
  about a quarter of fast deaths. They differ in the organ arrangement (A293 Part 1: 0.690 vs 0.161) but not in the
  loss of external time: losing the external day does not by itself remove the arrangement. No contrast tested.
- Code: `exploratory/A314_vent_walltime/scripts/score_a314.py`.
<!-- src: paper/exploratory/A314_vent_walltime/output/{groups,draws}.csv -->

## A313 (scored 2026-11-18; NOT BLIND for the organ offsets and the liver pairs of A302; powered): adjacent tissues do not share a person's clock departure beyond distant tissues once random genes are taken into account (H313.1 FAIL); composition does not either (H313.2 FAIL); the adjacency gradient exists for random expression-matched genes as well

- Design as registered; power stage (power 0.91 at an excess of 0.10), then one scoring run (6329 s, 16 workers, peak
  2.8 GB); clock and 100 matched sets; 261 of 300 pairs with >= 30 shared donors (the stomach pairs and colon
  transverse-stomach fall below).
- **H313.1 FAIL.** T_adj = median r(adjacent, different) - median r(distant, different) = 0.162 - 0.129 = 0.033;
  region reassignment p 0.098; matched sets median 0.024, q95 0.069; the clock at their 63rd percentile.
- **H313.2 FAIL.** T_fam = 0.131 - 0.129 = 0.002 (p 0.48; matched median 0.037, q95 0.100).
- **Descriptive.** Median r per cell, clock against the median over the matched sets: adjacent-similar 0.265 (3 pairs:
  oesophageal muscularis-junction 0.620, atrial-ventricle 0.265, aorta-coronary 0.191) against 0.210; adjacent-different
  0.162 against 0.137; distant-similar 0.131 against 0.150; distant-different 0.129 against 0.115. Highest single
  adjacent-different pair: colon sigmoid-transverse 0.519 (31 donors).
- **Reading (registered: FAIL, powered).** Departures are shared more by adjacent tissues, but random genes of the same
  expression show the same gradient (adjacent-similar > adjacent-different > distant-different). The regional sharing is
  a property of the samples (shared handling, region-level physiology or composition captured by any gene set), not of
  the clock; the clock-specific part is below the 0.10 the test could detect. With one systemic drive and no regional
  term, the uncoupled entrained-ensemble reading is not contradicted; local coupling or regional cues acting on clocks are
  not detected in GTEx. The 3 adjacent-similar pairs are too few to test; the oesophageal muscularis-junction pair (0.62)
  is effectively two samples of one tissue.
- Code: `exploratory/A313_local_coupling/scripts/score_a313.py --stage score --jobs 16`.
<!-- src: paper/exploratory/A313_local_coupling/output/{pairs,sets,tests,power,power_sets}.csv -->

## A315 (scored 2026-11-18; NOT BLIND for the arrangement results A132 / A293 Part 1 and A314; powered): brain-dead organ donors keep the clock amplitude of ventilated donors with a living brain (H315.1 FAIL, clock-specific ratio 1.05): after loss of the central drive the tissue clocks lose their arrangement, not their amplitude

- Design as registered; power stage (power 0.995 at a 30 % loss), then one scoring run (39 s). 298 organ donors, 79
  ventilated post-mortem donors, 12 tissues (42 to 269 donors per group and tissue).
- **H315.1 FAIL.** Median log ratio of the debiased clock-plane variance, organ donors over post-mortem, Q -0.136 (ratio
  0.87; donor bootstrap 90 % -0.299 to +0.014); random expression-matched genes fall more (median -0.182, q05 -0.303):
  the clock lies at their 69th percentile, clock-specific ratio exp(Q*) 1.05; responsive sets median -0.162 (clock at
  their 62nd percentile). L_t is negative in all 12 tissues and inside the matched spread (tibial artery -0.332, 90 %
  -0.555 to -0.054, matched median -0.174, is the largest).
- **Descriptive (registered).** Fast deaths over ventilated post-mortem: Q +0.174 against matched -0.230 (Q* +0.404; above
  all 200 sets): ventilated deaths carry clearly less clock amplitude than fast deaths, clock-specifically, whether the
  brain lives or not. Within organ donors, longest over shortest ventilation tertile (99 vs 101 donors): Q -0.182 against
  matched +0.053 (Q* -0.236; 1st percentile of the matched sets): clock amplitude falls with time on the ventilator.
- **Reading (registered: FAIL, powered).** Loss of the central drive (brain death) removes the organ arrangement (A132,
  A293 Part 1) but not the tissue clocks' amplitude, relative to donors who are equally ventilated: tissue clocks hold
  their amplitude without the drive and drift as organs, which the registration reads as coupled or self-sustained
  tissue clocks, against the uncoupled entrained-ensemble reading. Amplitude is lost with intensive care itself
  (fast vs ventilated) and, descriptively, with its duration, in both readings a property of the ICU state rather than of
  brain death.
- **Post hoc (not tested).** Arrangement lost with amplitude kept means cells within a tissue stay coherent while tissues
  drift apart: coherence is held inside tissues (local coupling or a tissue-level oscillator) and between tissues only
  by the drive. A313 found no shared departures between adjacent tissues, which concerns coupling between tissues, not
  within them. Time since brain death is not recorded; if it is short against the cells' dephasing time, an uncoupled
  ensemble with small period dispersion would also keep its amplitude, so the result bounds the dephasing rate rather
  than excluding the uncoupled reading outright.
- Code: `exploratory/A315_braindeath_gain/scripts/score_a315.py --stage score`.
<!-- src: paper/exploratory/A315_braindeath_gain/output/{power,tests,tissues,sets,descriptive}.csv -->

## A317 (scored 2026-11-18; NOT BLIND for the arrangement (A132b, A293 Part 1); the group split of C not seen): brain death changes the clock arrangement between organs clock-specifically (H317.1 PASS, D_rho 0.621 above all matched and responsive sets), but organs of brain-dead donors still share a clock state (H317.2 FAIL, excess C 0.269 vs 0.317): the order between organs is rearranged, not dissolved

- Design as registered; power stage (MDE D_rho 0.554, D_C 0.130), then one scoring run (59 s). 298 organ donors drawn
  to 79 (20 draws), 79 ventilated post-mortem donors, 244 fast deaths as reference, 65 pairs.
- **H317.1 PASS.** D_rho = rho(post-mortem) - rho(organ donors) = 0.621 (median over draws); relabelling p 0.0025; above
  all 200 matched (q95 0.092) and all 200 responsive sets (q95 0.152).
- **H317.2 FAIL.** Excess C over matched sets: post-mortem 0.317 (C 0.489), organ donors 0.269 (C 0.430); D_C 0.048,
  p 0.189 (MDE 0.130). Both groups are above all matched and all responsive sets (responsive excess 0.264 and 0.230).
- **Sensitivities (registered, descriptive).** Long-term steroid users excluded (74 / 295 donors): D_rho 0.575, D_C 0.037.
  RIN >= 6: D_rho 0.591, D_C 0.022. Common ischaemic-time range (5th-95th percentile of the post-mortem samples): D_rho
  0.239 on 36 pairs (65 before), D_C 0.005; companion medians 0.030 (matched) and -0.014 (responsive). The ischaemic
  restriction halves the pairs and the difference; whether ischaemic time explains part of the contrast or the loss of
  pairs adds noise is not separable here: ischaemic time is a candidate confound, named as registered.
- **Reading (registered: H317.1 passes, H317.2 fails).** Organs of brain-dead donors keep a shared, clock-specific
  covariance of nearly the size of donors with a living brain, but in a different arrangement: rearranged, not
  drifted. With A315 (tissue amplitude kept) and A314 (alignment to the hour of death lost equally in both groups):
  after brain death the organs keep their clocks and still follow something shared, at different relative phases.
- **Post hoc (not tested).** A shared drive of another kind (donor management: ventilation, infusions, temperature
  control, hormone replacement) would give exactly this: entrained receivers lock to whatever drives them, and their
  relative lags depend on the drive. The living brain would then set the body's organ order not by holding organs
  together, but by being the drive whose timing produces this order.
- Code: `exploratory/A317_braindeath_confirm/scripts/score_a317.py --stage score`.
<!-- src: paper/exploratory/A317_braindeath_confirm/output/{power,tests,sensitivity,sets}.csv -->

## A318 (scored 2026-11-18; NOT BLIND for A317): the brain-death arrangement change survives adjustment for ischaemic time and RNA quality (H318.1 PASS, D_rho 0.584); within organ donors, those with longer ischaemic times keep more of the arrangement (H318.2, +0.410, p 0.040, descriptive): ischaemic time is not a sufficient explanation, but not excluded as a contributor

- Design as registered; power stage (65 pairs kept, MDE 0.582), then one scoring run (46 s). Ischaemic time medians
  258 min (organ donors) and 561 min (post-mortem); RIN 7.6 and 7.2.
- **H318.1 PASS.** With every Hardy 0 state (clock and companions) residualised on ischaemic time and RIN, D_rho 0.584
  (A317 unadjusted 0.621); relabelling p 0.003; above all 200 matched (q95 0.090) and all 200 responsive sets (q95
  0.161). Adjusted D_C 0.088 (A317: 0.048).
- **H318.2 (descriptive).** Organ donors in the longest vs the shortest third of ischaemic time (99 vs 100 donors, 66
  pairs): rho difference against fast deaths +0.410, two-sided relabelling p 0.040: the direction a confound predicts.
- **Reading (registered).** A linear dependence on ischaemic time or RIN does not explain the contrast: the arrangement
  change of brain-dead donors stands, clock-specifically. Within organ donors the arrangement still varies with
  ischaemic time in the confound's direction, so a non-linear or threshold effect of ischaemia (or what it proxies in
  organ donors, e.g. procurement sequence or time since brain death) is not excluded. The brain-death reading is
  supported in GTEx, with ischaemic time stated as an unresolved contributor; an independent replication is needed.
- Code: `exploratory/A317_braindeath_confirm/scripts/score_a318.py --stage score`.
<!-- src: paper/exploratory/A317_braindeath_confirm/output/{a318_power,a318_tests,a318_sets}.csv -->

### A318 post hoc (2026-11-19; registered note before computing; descriptive, NOT BLIND): same-pairs decomposition of the common-ischaemic-range sensitivity: the restriction drops every pair with heart, lung or whole blood; on the remaining 36 pairs the contrast falls from 0.43-0.53 to 0.239, half by donor noise and half by organ donors' rho rising beyond noise (H318.2 selection); R36 sits at the 5th percentile of the noise-matched draws (1 of 20 below)

- Same 36 pairs, A317 estimator; rho with fast deaths (post-mortem / organ donors, median of 20 draws) and D_rho:
  | readout | rho_pm | rho_od | D_rho |
  |---|---|---|---|
  | U65 (A317 H317.1, reproduced) | 0.678 | 0.057 | 0.621 |
  | U36, A317 main draws | 0.424 | -0.109 | 0.533 |
  | U36, sensitivity draws | 0.424 | -0.009 | 0.433 |
  | S36, noise-matched (20 draws: q05 / median / q95) | 0.234 / 0.369 / 0.440 | -0.057 / -0.018 / 0.027 | 0.244 / 0.379 / 0.473 |
  | R36 (A317 sensitivity, reproduced) | 0.423 | 0.184 | 0.239 |
  | U29, the dropped pairs, unrestricted | 0.706 | 0.172 | 0.535 |
- Per-pair estimator check: 0.435 on U36 (sensitivity draws) against 0.433 for the A317 estimator; R36 0.250 with it.
- **Dropped pairs (29):** every pair containing Heart - Atrial Appendage, Lung or Whole Blood (Heart-Lung is not among
  the 65). Their organ-donor samples fall mostly below the post-mortem ischaemic range (201-1271 min): shared donors per
  pair in one organ-donor draw of 79 drop from 20-62 to 3-19 (`dropped_pairs.csv`).
- **Reading as registered.** R36 against S36: 1 of 20 noise-matched draws below 0.239 (2 of 20 below the like-for-like
  per-pair R36, 0.250); interpolated q05 0.244. The registered boundary ("below the 5th percentile") is met for the A317
  estimator and missed for the per-pair one: **borderline, not decided by 20 draws.** The decomposition is clearer than
  the boundary: post-mortem rho does not fall beyond noise (0.423 against S36 median 0.369; 17 of 20 draws lower);
  organ donors' rho rises above every noise-matched draw (0.184 against at most 0.066). The fall of the contrast under
  the restriction is therefore the H318.2 selection (organ-donor samples with post-mortem-range ischaemia are more
  post-mortem-like), not noise alone. U36 is not near R36 (0.43-0.53), and the dropped pairs carry the contrast as
  strongly as the kept ones (U29 0.535): the fall from 0.621 is not located in the dropped pairs.
- **Note for B3.** The three dropped tissues are those the restriction removes, i.e. the organ-donor tissues recovered
  first (heart and lung at procurement, blood): the tissues most likely reached by perfusion or recovered at cross-clamp.
  The organ-donor samples kept in the 9 remaining tissues are the later-recovered, longer-ischaemia ones.
- **Limits.** 20 noise-matched draws (registered); the organ-donor medians move by about 0.1 between two sets of 20
  draws (U36: 0.533 vs 0.433), the same order as the gap to the boundary.

### B4 (2026-11-19, no run): donors after circulatory death: closed, not available
No GTEx variable records donation after circulatory death. Subject data dictionary
(`phs000424.v8.pht002742.v8.GTEx_Subject_Phenotypes.data_dict.xml`): the candidate fields are TRTPTREF (tissue ischaemic
time from actual death, presumed death, or cross-clamp "in the case of surgical donors") and DTHTPTREF. Aggregates in
the analysis universe: all 298 Hardy 0 organ donors carry TRTPTREF and DTHTPTREF "Actual Death" and a cross-clamp time
(TRCCLMP recorded in 298 of 298); post-mortem ventilated donors 77 "Actual Death", 2 "Presumed Death", none with a
cross-clamp time. In the whole universe 300 organ donors, 298 Hardy 0 (one Hardy 2, one Hardy 3). GTEx's Surgical cohort
(15 donors) is not in the analysis universe. The direct control (same collection, no brain death) is not in GTEx.
- Code: `exploratory/A317_braindeath_confirm/scripts/same_pairs_a318.py`.
<!-- src: paper/exploratory/A317_braindeath_confirm/output/{same_pairs,same_pairs_summary,dropped_pairs}.csv -->

### Task 3 precondition (2026-11-19; aggregate run before registering; stopped): the living surgical cohorts cannot decide brain death vs collection mode
Prompt `docs/prompts/cell_route_2026-11-19.md` Task 3. GTEx pair offsets (phase of the pair covariance, h; A317
estimator) of the pairs living surgical cohorts carry; organ donors as the circular mean over the A317 draws of 79;
donor bootstrap (500) of the organ-donor (298) minus post-mortem (79) difference:
| pair | fast | Hardy 4 | post-mortem vent | organ donors | OD - PM (95 % bootstrap) |
|---|---|---|---|---|---|
| subcutaneous - visceral fat | -0.86 | -1.68 | -0.45 | 0.77 | 1.22 (-0.25 to 2.86) |
| muscle - subcutaneous fat | 1.01 | 1.81 | 0.58 | 0.58 | -0.01 (-0.92 to 1.08) |
| subcutaneous fat - lower-leg skin | -1.04 | -1.06 | -0.45 | -0.48 | -0.03 (-1.72 to 1.35) |
| subcutaneous fat - suprapubic skin | -1.87 | -2.39 | -1.05 | 0.63 | 1.68 (-0.24 to 3.44) |
Against the site bound (A251 same-layer site offsets: largest 1.54 h, interval to 2.04 h) no pair's organ-donor minus
post-mortem offset exceeds the bound with its interval (one point estimate, 1.68 h, exceeds 1.54 h; every bootstrap
interval contains 0). **By the stated rule Task 3 cannot decide and is stopped before registration.** The brain-death
contrast is spread over the 65 pairs; these four carry little of it.
- Code: `exploratory/A317_braindeath_confirm/scripts/living_ref_precondition.py`.
<!-- src: paper/exploratory/A317_braindeath_confirm/output/living_ref_precondition.csv -->

## A325 (scored 2026-11-19; NOT BLIND for A132b, A293 Part 1, A314-A318): post-mortem ischaemia does not reshape the reference organ arrangement (H325.2 stable, rho 0.862, p 0.62, power 0.99 at the organ-donor effect size; H325.1 stable, rho 0.778, p 0.68, uninformative at power 0.78): the post-mortem side carries no collection-mode component of the organ-donor size

- Donor ischaemic time (median over samples); tertiles cut within hour-of-death strata (08-20, 20-08 h; x death group
  for H325.2), against 2000 random disjoint splits within the same strata. Night deaths wait longer (fast deaths 844 vs
  1039 min; Spearman with cos(hour) 0.260): stratification was needed.
- **H325.1 (fast deaths, 82 / 82, 63 pairs; ischaemic medians 1214 vs 659 min):** rho_tb 0.778 against null median
  0.730, q05 0.449; p 0.680. Companions: matched p q05 0.134, responsive 0.127 (the clock is not below them). Power at the
  full organ-donor shift 0.78: **stable, uninformative by the registered rule.**
- **H325.2 (all post-mortem deaths, 161 / 161, 66 pairs; 1166 vs 581 min):** rho_tb 0.862 against null median 0.845,
  q05 0.705; p 0.620. Power 0.99 at f = 1, 0.40 at f = 0.5: **stable, informative.**
- Unstratified fast deaths (descriptive): 0.754, p 0.559.
- **Reading as registered.** A 555-585 min contrast of post-mortem ischaemia, larger than the organ-donor minus
  post-mortem 300 min, leaves the arrangement as stable as random splits. Post-mortem ischaemia (slow cooling, hours of
  anoxia) does not reproduce a contrast of the organ-donor size; a component of half that size is not excluded (power
  0.40). This argues against collection mode acting on the post-mortem side. It does not test the organ-donor side
  (perfusion, recovery order, short ischaemia at procurement), which H318.2 and the A318 post hoc note point to.
- **Limits.** Ischaemic time stands for the whole post-mortem collection process; the organ-donor procurement itself
  has no counterpart among post-mortem donors.
- Code: `exploratory/A317_braindeath_confirm/scripts/score_a325.py --stage power | score`.
<!-- src: paper/exploratory/A317_braindeath_confirm/output/{a325_power,a325_tests,a325_sets,a325_hour}.csv -->

### B1 power scout (2026-11-19; before registering; stopped): the baboon cannot referee the GTEx organ arrangement
Prompt Task 4 / REVIEW_2026-11-19 check B1. Shared A109 tissues with the baboon atlas (A36 timed phases): skeletal
muscle, thyroid, visceral fat, suprapubic skin, subcutaneous fat, lung (15 pairs). Fisher-Lee rho of GTEx pair phases
against baboon pair phases, exact tissue-label permutation null (720 orders), both sign conventions read: fast deaths
|rho| 0.310 (p 0.235 in the agreeing orientation), Hardy 4 0.345 (p 0.218), post-mortem ventilated 0.171 (p 0.383),
organ donors 0.060 (p 0.425). Over 25 tissues A36 itself found weak cross-species agreement of organ phases (r_FL 0.175,
p 0.199; `results/figure_1/output/species/a36_readouts.csv`). The fast-death arrangement does not clearly agree with the
baboon's: **B1 is uninformative and stopped** (orientation not chosen; the scout does not support either reading).
- Code: `exploratory/A317_braindeath_confirm/scripts/b1_scout.py`.
<!-- src: paper/exploratory/A317_braindeath_confirm/output/b1_scout.csv -->

## A335 (scored 2026-11-19; NOT BLIND for the 14-gene originals, BLIND for every core-loop, matched, cross-fitted and equivalence statistic except the H293.3 interval, known from A293; one run per part; every 14-gene original reproduced before any new value was read): the night-shift loss and shift, the ICU loudness loss and the order-reading of time of death survive the review's checks; the psoriatic shift survives at the line; living person loudness does not survive the core loop; brain-death loudness is equivalent (+-0.357); the D3 remainder is not estimable

**Part 1, core loop** (11 genes, without PER1, PER2 and NR1D1; companions restricted to the core positions, as A329).

| check | 14 genes (original) | core loop | registered verdict |
|---|---|---|---|
| 1a night shift, muscle (A257): shift s | 2.68 h (90 % 1.81-3.59) | 4.06 h (2.57-5.68) | shift holds |
| 1a amplitude ratio R | 0.234 (0.060-0.414) | 0.239 (0.075-0.398); swap p <= 1/2048; rho_ctl - rho_mis 1.45 > matched q95 0.38 and responsive q95 0.81 | loss holds |
| 1b psoriasis (A308): pooled offset | -1.45 h, q 0.0156 | -1.00 h, q 0.0486 (90 % 0.0003-0.130); 5.5th pct of matched, 1.0th of responsive | shift holds, at the line |
| 1c MuTHER person loudness (A331): leverage | 0.177 (98.5th matched, 93.5th responsive) | 0.107, p 0.0025; matched q95 0.130 (89th pct), responsive 69.5th | falls |
| 1d ICU loudness, fast vs ventilated post mortem (A315): Q* | +0.404, above all matched | +0.457, above all matched and responsive (Q 0.228 > q95 -0.086) | ICU loss holds |
| 1d brain-death loudness, H315.1: Q* | +0.046 (69th pct) | -0.006 (47.5th pct) | FAIL as registered (no loss) |
| 1e brain-death arrows (D2) | lung -1.58, tibial artery +1.09, tibial nerve +0.99 h | lung -1.47 (-2.92 to +0.21), artery +1.20 (0.17-2.03), nerve +1.54 (0.19-2.73), lower-leg skin +1.33 (0.18-2.64) | artery and nerve hold; lung does not; lower-leg skin added |

- 1a per subject (descriptive, 11 subjects): angles -1.39 to +6.57 h (14 genes) and +0.66 to +9.30 h (core; all 11
  positive); amplitude ratios |dm| / |dc| 0.40-0.94 and 0.34-0.92. The night-shift damping is not carried by the
  signal-inducible genes; on the core loop the shift is larger and uniform in sign.
- The H257.3 verdict was read with the registered swap test and R interval plus the stricter companion comparison; both pass.
- Note (2026-11-20, independent review section 9).
  - R = rho_mis / rho_ctl projects the misaligned chord on the expected one, so it contains the rotation. Shift-corrected,
    the core R is 0.239 / cos 60.9 deg = 0.49 (derived), matching the per-subject median |dm| / |dc| of 0.47: the
    misaligned muscle clock keeps about half its amplitude.
  - The swap p is 1/2048. The CSV's 0 is a floating-point tie at the unswapped labelling, which an exact swap test
    counts.

**Part 2, covariate balance** (standardised mean differences, post mortem minus comparator, before / after matching 79
post-mortem cases to up to 3 controls by the A164 rule).
- **Organ donors:**
  - age 0.78 / 0.22, sex 0.02 / 0.10;
  - RIN -0.82 / -0.78 and ischaemic time 1.23 / 1.26 (collection mode; not balanceable).
- **Fast deaths:**
  - age -0.10 / -0.10, sex -0.27 / -0.18;
  - RIN 0.12 / 0.14, ischaemic time -1.09 / -1.06.
- **Matched results:**
  - ICU Q* is +0.389 (14 genes) and +0.459 (core), both above every matched set: **ICU loss robust** to age and sex.
  - Organ donors against post mortem: Q* +0.074 and +0.017 (80th and 54.5th percentiles).

**Part 3, equivalence (90 % intervals).**
- **H315.1 brain-death loudness, Q* within +-0.357:**
  - 14 genes: -0.135 to 0.181 pooled, -0.124 to 0.206 matched;
  - core: -0.147 to 0.137 pooled, -0.157 to 0.161 matched.
  - All four intervals are inside the margin: **equal within a 30 % change of clock-plane variance**, relative to
    expression-matched genes. The raw clock-plane variance is lower in organ donors (core Q -0.186, 90 % -0.307 to
    -0.006), as is the matched sets' (median -0.181).
- **H293.3 class phase order:** the endothelium-minus-fibroblast offset is 0.46 h (0.09-0.90; reproduces A293). It lies
  within +-1 h (equivalent), but not within +-0.6 h. Not blind: the interval was known from A293 before the margin was
  set.
- **H295.1 composition (margin post hoc):** -0.010 per SD (-0.029 to +0.010), within +-0.10.

**Part 4, A270 with a held-out waveform.**
- Halves of 122 fast deaths each, with 12-27 deaths in every 4-h bin.
- The clock R with time of death through a waveform refitted on the other half is 0.562 and 0.714, against 0.546 and 0.716
  with the shipped waveform. It is above all 400 companion sets in each half (maxima 0.211 and 0.291); p 0.001.
- The refitted and shipped acrophases agree (phase concordance 0.987 and 0.966).
- **Holds:** the clock's reading of time of death is not an artefact of fitting the waveform on the same donors.

**Part 5, D3 remainder: not estimable as registered.**
- All 79 post-mortem donors have a ventilation duration.
- Tertiles were cut within day and night strata: top 24 donors (median 372 h), bottom 28 (18.5 h).
- No tissue pair has >= 20 donors shared by fast, top and bottom, so the order rho is undefined. Loudness rests on one
  tissue (Q* -0.259, core -0.595); this is uninformative. No relaxed rule was applied.

**Reading.**
- **Survives:**
  - the night-shift shift (larger on the core loop) and loss (about half the amplitude kept per subject);
  - the ICU loudness loss (core and age / sex matched);
  - the clock's reading of time of death (held-out waveform);
  - brain death leaving loudness unchanged relative to matched genes (now an equivalence, not a failed test);
  - class time-sharing within 1 h (not blind);
  - the psoriatic shift, only at the line (q 0.049).
- **Does not survive:**
  - living person loudness beyond matched genes (A331): carried by PER1, PER2 and NR1D1, or lost precision (see the
    last note);
  - the lung arrow of brain death.
- **Not estimable:** maintenance against inertia in the ventilated post-mortem donors (D3).
- By the rule written before scoring, the Article route stands: the night-shift loss and the psoriatic shift survive.
  The psoriatic shift should be stated with its marginal core q.
- Note (2026-11-20, independent review): the registration's limits line ("the core loop loses the three
  highest-amplitude genes") is wrong. In the shipped `cs.gtex()` waveform PER1, PER2 and NR1D1 rank 14th, 6th and 3rd
  in amplitude. They hold 16 % of sum (amp / sd)^2 but fill the clock plane's weaker direction (second singular value of
  the sd-scaled design 0.98 with 14 genes, 0.64 on the core loop; derived). A core-loop fall of a magnitude statistic
  can therefore be lost precision; the gene-drop null of `prompts/two_readouts_2026-11-20.md` (part 1) separates the
  two.
- Code: `exploratory/A335_review_batch/scripts/{a335_1a_nightshift.py, a335_1b_psoriasis.py, a335_1c_muther.py,
  a335_gtex.py, a335_3_equivalence.py, a335_4_a270_crossfit.py}` (seconds to 9 min each).
<!-- src: paper/exploratory/A335_review_batch/output/{a335_1a_nightshift,a335_1a_per_subject,a335_1b_tests,a335_1c_tests,a315_core,d2_core,balance,matched,equivalence,part3_equivalence,part4_tests,part4_hours,d3_tertiles,repro_14gene}.csv -->

## A336 (scored 2026-11-20; NOT BLIND for the 14-gene values and the A335 / A329 core values; BLIND for every other core value, the gene-drop nulls, the GTEx person leverage on the core loop and the brain-death D_rho on the core loop; one scoring run per part and statistic; every 14-gene value reproduced before any new value was read): the paper's measurements survive without PER1, PER2 and NR1D1; person loudness in GTEx deaths holds on the core loop (it is a clock property in deaths); the brain-death order holds on the core loop; a detection limit now sits beside every "kept"; tier 2: the Figure 2 cell-class contrasts hold on the core loop against expression-matched genes and, in half of the tests, against responsive genes (three PASS -> at the line, cortex not decidable), tumour immune pairing at the line, chronic-disease nulls (A233) hold

**Part 1, two readouts** (`two_readouts.csv`; tier 1: 47 inventory statistics, 55 rows; tier 2: all 31, below; tier 3 not run,
listed in the table). Values: 14 genes -> core; trio contribution = 14 minus core (90 %); drop % = percentile of
the trio drop among the 364 three-gene drops (or 100 + trio) / among the 77 drops of similar singular values (24 of
101 where subsampled).

| figure | statistic | 14 genes -> core | trio contribution (90 %) | drop % all / similar | reading |
|---|---|---|---|---|---|
| 1b | A252 H252.1 rotation rho, exact times | 1.265 -> 1.073 | +0.192 (0.132 to 0.253) | 3.6 / 14.3 | holds |
| 1b | A252 H252.2 seven studies / H252.3 adipose | 2.970 -> 3.902; 2.221 -> 3.861 | -0.93; -1.64 | 95 / 77; 86 / 39 | holds; holds |
| 1b | A256 H256.2 rho / H256.3 cells within 2 h | 2.533 -> 3.159; 17 -> 18 of 20 | -0.63; -1 (-3 to 2) | 90 / 58; 96 / 91 | holds; holds |
| 1c | A270 clock R with hour of death, human / baboon | 0.631 -> 0.627; 0.987 -> 0.982 | +0.003; +0.006 | 49 / 42; 14 / 42 | holds; holds |
| 1d | A257 H257.1 rho / H257.2 shift / H257.3 loss 1 - R | 1.336 -> 1.912; 2.68 -> 4.06 h; 0.767 -> 0.761 | -0.58; -1.38 h (-2.22 to -0.59); +0.005 | 96 / 81; 98 / 92; 51 / 81 | holds x3 |
| 1d | A261 H261.3 programs (core constants 4.06 h, 0.239) | behaviour -> behaviour | - | - | holds |
| 2b | A275 clock sharing SV_X / LY_X (descriptive) | 0.333 -> 0.295; 0.172 -> 0.159 | +0.038 (0.010 to 0.065); +0.013 | 1.1 / 2.6; 16 / 6.5 | holds (rule "both") |
| 4b | A295 H295.2 phase per SD (half-extent, h) | 0.130 -> 0.151 | -0.02 | 8.2 / 10.4 | holds |
| 4b | A295 H295.1 lymphocyte gain per SD | -0.010 (FAIL) -> -0.020 (p 0.040) | +0.010 | 89 / 53 | verdict changes: PASS on the core |
| 4c | A293 H293.1 / H293.2 organ order rho (slow, ventilated) | 0.784 -> 0.831; 0.690 -> 0.700 | -0.05; -0.01 | 97 / 88; 72 / 73 | holds; holds |
| 4e | A291 coupling ~ loudness | 0.717 -> 0.490 | +0.23 (0.09 to 0.46) | 13 / 52 | holds |
| 4e | A294b H294b.1 departure / matched (median log ratio) | 0.969 -> 0.250 (0.14 to 0.36) | +0.72 (0.61 to 0.83) | 84 / 32 | holds |
| 4e | A294b H294b.2 slope | 0.223 -> 0.464 | -0.24 | 64 / 25 | falls: not decidable (precision) |
| 4f | A312 H312.1 elasticity ~ S (rho) | -0.341 -> -0.291 (p 0.063) | -0.05 | 35 / 79 | at the line |
| 4f | A331 MuTHER person leverage | 0.177 -> 0.107 (< matched q95 0.130; p 0.0025) | +0.070 (0.008 to 0.135) | 13 / 43 | at the line (fails the matched condition) |
| 4f | GTEx person leverage (part 2) | 0.765 -> 0.714 | +0.051 (0.020 to 0.087) | 13 / 40 | holds |
| 5b | A315 ICU Q* / H315.1 brain-death Q* | 0.404 -> 0.457; 0.046 -> -0.006 | -0.05; +0.05 | 81 / 78; 15 / 58 | holds; holds (equivalence) |
| 5b | A317 H317.1 / A318 H318.1 D_rho (part 4) | 0.621 -> 0.465; 0.584 -> 0.496 | +0.16 (-0.11 to 0.35) | 16 / 39; 34 / 39 | holds; holds |
| 5c | A297 removal / A322 Q / A324 Q, H324.2 (mouse) | 0.724 -> 0.940; 31.9 -> 27.2; 185 -> 95; FAIL -> FAIL | -0.22; +6.0; +116 | 98 / 88; 29 / 62; 2.5 / 12 | verdicts unchanged |
| 5d, 6b | A304 gain / departure H304.1, H304.3, H304.4 | -0.462 -> -0.668; -0.190 -> 0.089; -0.522 -> 0.017; -0.216 -> -0.762 | +0.21; -0.28; -0.54; +0.55 | 98; 81; 94; 9.6 | holds x4 (departure not detected, powered) |
| 6a | A235 H235.1 metabolic offset, T2D / obesity (equivalence +-1 h) | +0.02 -> +0.14 h (-0.19 to 0.47); +0.20 -> +0.14 h (-0.32 to 0.59) | -0.12; +0.07 | 82 / 29; 8.9 / 33 | holds; holds |
| 6a | A164 H164.1 / H164.2 d* (equivalence +-0.15) | 0.019 -> -0.034; 0.004 -> 0.014 | +0.053; -0.009 | 94; 55 | holds; holds |
| 6a, 5d | A306 gain / departure, T2D, HTN, smoking | see table | - | - | departure HTN holds; five fall: interval crosses the registered MDE, nothing detected (HTN gain "carried by the trio" mechanically, core = 14-gene value) |
| 6b | A138 H138.4 tumour-normal pairing | 0.365 -> 0.252 (< matched q95 0.275, responsive 0.313) | +0.112 (0.073 to 0.159) | 6.9 / 29 | at the line |
| 6b | A138 H138.8 purity D | 0.025 -> 0.109 (p 0.049) | -0.08 | 97 / 88 | falls: not decidable |
| 6b | A271 H271.2 GTEx coupling rho | 0.734 -> 0.501 (> matched q95 0.236, < responsive 0.553) | +0.23 | 15 / 54 | at the line |
| 6c | A254 H254.1 lesion offset / H254.3 lesion r | +0.01 -> +0.31 h (-0.14 to 0.75); -0.417 -> -0.325 | -0.29; -0.09 | 13 / 38; 20 / 65 | holds; at the line (not beyond responsive) |
| 6d | A183 H183.1 psoriasis offset / A308 q | -1.45 -> -1.00 h; q 0.016 -> 0.049 | -0.45 h (-1.14 to -0.06) | 19 / 83 | at the line (matched q05); holds (q < 0.05) |
| 6d | A183 H183.2 psoriasis minus eczema | -1.09 -> -0.72 h | -0.37 | 20 / 74 | falls: not decidable |
| 6d | A115 D115.1 uninvolved skin | -1.81 -> -2.11 h | +0.30 | 85 / 55 | holds |
| 6d | A332 D1 psoriasis r (descriptive) | +0.499 -> +0.129 | +0.37 (0.18 to 0.59) | 3.8 / 17 | descriptive |

- **No statistic reads "carried by PER1, PER2 and NR1D1"** under the registered rule except A306 HTN gain, where the
  core equals the 14-gene value (0.216 against 0.213) and only the interval crosses the margin: a mechanical reading of
  a null. Every other fall is "not decidable: precision" (A294b H294b.2, A138 H138.8, A183 H183.2, four A306 nulls).
- **Trio drops at an extreme of the gene-drop null** (descriptive): low (the trio carried more than three random genes)
  for A252 H252.1 (3.6), A275 SV_X (1.1), A332 D1 (3.8), A324 Q (2.5), A138 H138.4 (6.9); high (dropping the trio
  strengthens the statistic) for A257 shift (98), A297 (98), A304 gain (98), A293 H293.1 (97), A256 (96).
- **Trio alone: never reportable.** With three genes the acrophase gates admit six arrangements (p >= 1/6), so every trio
  gate failed on p_acrophase (A304's pairing test is the exception). The night-shift trio-alone shift (the one new layer
  statistic) is therefore not read (A257 trio gate p 0.994).
- **Night-shift swap test:** compared with a tolerance, core p = 1/2048 (A335 printed 0.0 from a floating-point tie).
- Gene-drop subsampling (100 + trio) by the registered cost rule: A138, A164, A271, A306, A312 (A312 planned at 364 in
  the inventory; the rule applied at 3-4 min per waveform).


**Part 1, tier 2** (scored 2026-11-20, same rules; 14-gene verdict from each entry's stored test table, column
`verdict_14` of `two_readouts.csv`; a core "fall" on a statistic that already fails on 14 genes is no change of verdict).
A283 is the calibrated restatement of A272 H272.1-2 and A274 H274.1-3 (same clock values; the uncalibrated q95 readings
agree and are in the table).

| figure | statistic | 14-gene verdict | 14 genes -> core | trio contribution (90 %) | drop % all / similar | reading |
|---|---|---|---|---|---|---|
| 2a | A283 H283.1 T_SV, A272 atlases | FAIL | 0.499 -> 0.374 (matched p_cal 0.010 -> 0.181) | +0.126 (0.050 to 0.147) | 1.0 / 4.2 | falls (no change of verdict; the matched margin is lost, trio drop low) |
| 2a | A283 H283.2 D, A272 atlases | PASS | 0.277 -> 0.137 (matched p_cal 0.041, responsive 0.280) | +0.139 (0.025 to 0.153) | 4.0 / 16.7 | at the line |
| 2a | A283 H283.3 D, A274 atlases | PASS | 0.152 -> 0.131 (p_cal 0.040, 0.039) | +0.020 (-0.127 to 0.084) | 39 / 33 | holds |
| 2a | A283 H283.4 T_SV, A274 atlases (matched decisive) | PASS | 0.478 -> 0.416 (matched p_cal 0.040) | +0.061 (-0.009 to 0.101) | 13 / 50 | holds |
| 2a | A283 H283.5 T_CX cortex | PASS | 0.497 -> 0.365 (p_cal 0.632, 0.548) | +0.132 (0.091 to 0.152) | 18 / 67 | falls: not decidable (precision) |
| 2a | A283 V283 calibration on the core projector | PASS | max rate 0.063 -> 0.060 | - | - | holds |
| 2a | A280 H280.1 D replication / H280.2 cortex | PASS; FAIL | 0.413 -> 0.308 (p_cal 0.002, 0.010); 0.279 -> 0.225 | +0.105 (0.044 to 0.150); +0.054 | 5.0 / 21; 19 / 63 | holds; falls (no change) |
| 2a | A282 pooled D living / post mortem (descriptive) | both p_cal < 0.05 | 0.225 -> 0.081 (p_cal 0.221, 0.075); 0.261 -> 0.208 (0.006, 0.017) | +0.144 (-0.009 to 0.166); +0.053 | 3.0 / 12.5; 39 / 50 | living: falls, trio drop low (descriptive); post mortem: holds |
| 2a | A282 pooled T_SV living (descriptive) | matched only | 0.450 -> 0.362 (matched p_cal 0.036) | +0.088 (0.014 to 0.132) | 2.0 / 8.3 | at the line (unchanged) |
| 2a | A288 H288.1 T_MAC / H288.2 D_MAC | FAIL; PASS | 0.409 -> 0.331; 0.148 -> 0.098 (matched p_cal 0.011, responsive 0.143) | +0.077 (0.036 to 0.110); +0.050 (0.004 to 0.107) | 4.0 / 17; 5.0 / 21 | falls (no change); at the line |
| 2a | A288 H288.3 alveolar vs monocyte-derived | FAIL | 0.165 -> 0.048 | +0.117 (-0.017 to 0.206) | 5.0 / 21 | falls (no change) |
| 2b | A326 H326.1 lymphocyte bound r* (A281 / A287 / A321) | FAIL (underpowered at r 0.25) | r* 0.7 / 1 / 0.7 -> 0.7 / 1 / 0.7 | - | - | holds (same bound) |
| 2c | A279 Bmal1 knockout (descriptive) | conditions met | max timed coupling 0.200 -> 0.204; SV pairs \|r\| 0.68-0.80 | - | 97 / 88 | at the line (rounding rule) |
| 4a | A273 V273.1 / V273.2 mouse colon offsets | PASS | median error 0.31 -> 0.41 h; slope 1.01 -> 0.97 | - | 15; 22 | holds; holds |
| 4a | A279 V279.1 / H279.1 mouse aorta | PASS; FAIL | 0.30 -> 0.18 h; -0.055 -> -0.071 | - | 89; 43 | holds; falls (no change) |
| 4a | A293 H293.3 class offset, equivalence +-1 h | PASS | 0.46 h (0.09-0.90) -> 0.29 h (-0.19 to 0.80), 4 readable atlases | +0.17 | 14 / 33 | holds |
| 6b | A227 H227.1 myeloid / H227.2 T-NK pairing with normal | PASS; PASS | 0.581 -> 0.473; 0.539 -> 0.502 (pairing p 0.002, 0.0001; below matched q95 0.543, 0.518) | +0.108 (0.001 to 0.137); +0.036 | 6.9 / 29; 29 / 54 | at the line; at the line |
| 6b | A271 H271.1 immune pairing beyond responsive genes | FAIL | myeloid 0.473, T-NK 0.502, CD8-neg T 0.496 (all below responsive q95) | - | 6.9; 29; 34 | falls (no change) |
| 6a | A233 H233.c.1 body-wide d*, five conditions (equivalence +-0.15) | PASS x5 | HTN -0.001 -> +0.024 (90 % -0.000 to 0.048); IHD 0.009 -> 0.005; COPD -0.017 -> -0.026; depression 0.010 -> 0.008; smoking 0.008 -> 0.008; core MDE 0.029-0.059 | HTN -0.025 (-0.041 to -0.002); others within +-0.06 | HTN 7.9 / 29; others 18-44 | holds x5 |
| 6a | A233 H233.c.2-3 reordering Q, module d* and module offset, five conditions | null kept x15 | core Q p 0.13-0.80 (Holm >= 0.65); module p >= 0.08 (Holm >= 0.40); module offset MDE 0.37-3.68 h | - | - | holds x15 (nulls kept) |

- **Verdicts that change on the core loop:**
  - PASS -> at the line: A283 H283.2 (stromal minus lymphocyte, A272 atlases), A288 H288.2 (macrophage minus
    lymphocyte), A227 H227.1-2 (tumour immune cells pair with normal, no longer beyond matched genes). Each keeps the
    matched (A283, A288) or pairing (A227) condition and loses the other.
  - PASS -> falls, not decidable: A283 H283.5 (cortex; cortex was already unresolved, A280 H280.2).
  - Descriptive: A282 living D loses both margins, with the trio drop at the 3rd percentile.
- **Hold:** A283 H283.3-4 and A280 H280.1 (the stromal-lymphocyte contrast and stromal sharing in the A274 / A280
  atlases), A282 post-mortem D, the mouse offsets (A273, A279 V279.1), the human class offset (A293 H293.3), the
  lymphocyte bound (A326).
- **The trio carries a disproportionate share of cross-class sharing in single-cell atlases** (descriptive): its drop
  lies at the 1st-7th percentile of the gene-drop null for A283 H283.1-2, A280 H280.1, A282 living D and T_SV, A288
  H288.1-3 and A227 H227.1, as for A275 SV_X (1.1, tier 1). Among drops of similar singular values (24) the
  percentiles are 4-29: lost precision explains part of it.
- **Process notes:**
  - `p1_a227.py` ran twice to completion (08:37 and 09:04, same seeds); the first run's log was overwritten. Values are
    those of the second run.
  - A279 Bmal1 KO: the "at the line" status was set by the subagent after the run under the registration's rounding
    clause (0.2043 against the stated 0.20); numbers unchanged.
  - A326: the registered H326.1 (power >= 0.8 at r 0.25) fails on both readouts; "holds" means the same bound.
  - A288 and A233 were run by the lead agent after the subagents stopped (usage limit); A233's check mode needed
    `--check-donors` (60 donors gave no estimable cohort); its first scoring launch was stopped during the core nulls
    (no value computed) and relaunched with 12 workers (per-replicate seeds: identical draws).

**Part 2, person loudness in GTEx deaths on the core loop** (`person_gtex.csv`; 244 fast deaths, 27 tissues, 20 halvings).
- 14 genes (not blind; A331 power cv2 0.1239 reproduced): leverage 0.765.
- Core: leverage 0.714, persons-permuted p 0.001 (null q95 0.030), above the core matched q95 0.621 (98.5th percentile)
  and the responsive q95 0.692 (96th). **Holds.**
- Trio contribution +0.051 (0.020 to 0.087); trio drop at the 13th / 40th percentile of the drops.
- The matched sets are high (median 0.119, q95 0.621): a donor-level scale shared by the tissues is large in deaths, and
  the clock is above it.
- Strata (descriptive): core leverage 0.69-0.74 in every tertile of age, ischaemic time and RIN.
- Reading (registered): person loudness is a clock property in deaths, not shown in living tissue (A331 at the line on
  the core loop); Figure 4f keeps GTEx, labelled post-mortem.

**Part 3, detection limits** (`detectability.csv`; power only, no verdict).
- Time kept in intensive care (A293 H293.2): on fast-death halvings, one organ must shift by a median 3.9 h (2.25 h
  tibial nerve to 9.75 h whole blood; lung 8.75 h) to be detected; all organs dispersed, SD 1.5 h.
- Chronic disease time: metabolic offset MDE 0.37 h (T2D), 0.64 h (obesity); module offsets 0.28 h (HTN) to 3.07 h
  (depression).
- Chronic disease coupling: d* MDE 0.083 / 0.130 (A164), 0.028-0.048 (A233); no conversion to hours or fold.
- Chronic disease loudness and departure: fold x1.30 to x1.52 (A306).
- Tumour time (pairing): >= 2.4 h dispersion SD (lower bound).
- Lesion time 0.69 h; class time 0.61 h; composition 0.12 h and x1.03 per SD.
- Brain-death loudness: amplitude x0.90 / x1.11 (registered power 0.995 at Q* -0.357).
- Tumour departure: x1.5.

**Part 4, brain-death order on the core loop** (`braindeath_core.csv`; power stage: core MDE 0.622, A317 0.554).
- 14 genes reproduced (D_rho 0.621, p 0.0025, q95 0.092 / 0.152; A318 0.584).
- Core: D_rho 0.465, relabelling p 0.036, above the core matched q95 0.084 and responsive q95 0.180. **Holds.**
- Adjusted for ischaemic time and RIN: core 0.496, p 0.030, above both q95 (0.091, 0.193). Holds.
- Trio contribution +0.16 (-0.11 to 0.35). Trio drop at the 16th / 39th percentile (adjusted 34th / 39th).
- Trio gate fails (p_acrophase 0.83).
- The core D_rho is below its own MDE (0.622) yet significant.
- Reading (registered): the rearrangement is in the core clock, which answers the leakage objection. Donor-management
  glucocorticoids acting on the core are still not excluded (A136's gate failed). The A318 adjustment removes mean
  effects only, and ischaemic time and RIN cannot be balanced across collection modes (A335 part 2).

**Reading.**
- The split is informative but not destructive. Every Figure 1 statistic holds on the core loop, and the night-shift
  shift is larger there (not blind).
- Person loudness is a clock property in deaths (GTEx core) but not shown beyond matched genes in living women (MuTHER
  core at the line).
- The brain-death order is a core-clock rearrangement.
- Weakened to "at the line" on the core loop: the tumour-normal pairing (A138 H138.4), GTEx strong-clock coupling
  (A271 H271.2), lesion loudness (A254 H254.3), the psoriatic offset against matched sets (A183 H183.1), and the elasticity
  gradient (A312).
- Five chronic-disease loudness / departure nulls (A306) lose their equivalence on the core loop (wider intervals),
  without any detection.
- Tier 2: the cell-class structure of Figure 2 holds on the core loop against expression-matched genes (A283 H283.2-4,
  A280 H280.1, A288 H288.2). Against responsive genes the stromal-lymphocyte contrast holds in three tests (A283
  H283.3, A280 H280.1, A282 post mortem) and is lost in three (A283 H283.2 and A288 H288.2, at the line; A282 living,
  descriptive). Cortex (A283 H283.5) is not decidable. The trio carries a disproportionate
  share of cross-class sharing (drops at the 1st-7th percentile; descriptive). Tumour immune cells still pair with
  their normal tissue but not beyond matched genes (A227 at the line). Class time (A273, A279, A293 H293.3), the
  lymphocyte bound (A326) and the five A233 chronic-disease nulls hold.
- Code: `exploratory/A336_two_readouts/scripts/` (a336_braindeath.py, a336_person_gtex.py, a336_detectability.py,
  p1_*.py, a336_common.py, merge_part1.py).
<!-- src: paper/exploratory/A336_two_readouts/output/{inventory,two_readouts,person_gtex,person_gtex_strata,person_gtex_companions,detectability,detectability_order_sim,braindeath_core,braindeath_core_power,braindeath_drops}.csv, output/part1/*.csv -->

## A333 (power stage 2026-11-19; BLIND, no observed Phi computed): N1 not decidable at this design: with 79 living-brain ventilated donors (42-66 per tissue), the median MDE of log Phi over tissues is 0.726 (a loss to 0.48 of the reference), against the registered log(1 / 0.7) = 0.357; not scored

- **Design as registered.**
  - Fast deaths with a clock table: 244, split once into halves A and B (122 each; seed [333, 1]). V: 79 donors (Hardy 0,
    post-mortem cohort).
  - All 12 body tissues are estimable: fast B 62-114 donors, V 42-66, program genes 58 (whole blood) to 765 (atrial
    appendage), selected in half A with the A91 model and the responsive families excluded.
- **Power (null construction).** The V labels were permuted among the 201 B and V donors 200 times, with each group's
  leave-one-out host refitted every time:

| tissue | SE log Phi | MDE | detectable ratio |
|---|---|---|---|
| tibial artery | 0.207 | 0.516 | 0.60 |
| subcutaneous adipose | 0.224 | 0.557 | 0.57 |
| thyroid | 0.236 | 0.587 | 0.56 |
| tibial nerve | 0.259 | 0.646 | 0.52 |
| suprapubic skin | 0.287 | 0.714 | 0.49 |
| skeletal muscle | 0.290 | 0.723 | 0.49 |
| lung | 0.293 | 0.728 | 0.48 |
| oesophagus mucosa | 0.301 | 0.749 | 0.47 |
| lower-leg skin | 0.358 | 0.891 | 0.41 |
| visceral adipose | 0.384 | 0.955 | 0.39 |
| atrial appendage | 0.449 | 1.119 | 0.33 |
| whole blood | 0.556 | 1.384 | 0.25 |
| **median** | | **0.726** | **0.48** (gate 0.357: FAIL) |

- **Bias under the null.** The permuted log Phi has a positive mean in every tissue, 0.007 (muscle) to 0.318
  (oesophagus mucosa). The smaller group's noisier loadings inflate the median |coefficient| of the program, so any scoring
  would need this bias removed (descriptive).
- **Reading.** Whether adult tissue outputs keep following their clock without behaviour is not decidable in GTEx at
  this design. Only losses to about half would be detectable per tissue. The archive A91 lock (pooled ventilated donors,
  brain-dead included, not blind) stays the only indication, and Figures 1 and 7 keep the N1 panel as an open question.
  A pooled-over-tissues statistic would have more power, but tissues share donors and it is not the registered gate; it
  would need its own registration.
- Code: `exploratory/N1_output_retention/scripts/score_a333.py --stage power --jobs 4` (411 s); the score stage is not
  written (stop rule).
<!-- src: paper/exploratory/N1_output_retention/output/{program,power,power_perm}.csv -->

## A334 / A334b (scored 2026-11-19; BLIND for every per-cell statistic; one calibration run, one power run; neither aorta contrast scored): dephasing and damping of cell clocks cannot be told apart with public mouse 10x data: the lambda / R split fails the skin calibration (V334.1 FAIL: mixing-induced dephasing read 59 % as dephasing) and the dispersion D is not clock-specific in the aorta (A334b: Bmal1 knockout D 2.7 times aligned)

**A334, skin calibration (GSE223109 dermal fibroblasts, 18 samples, 12,224 cells, median 13 clock UMI per cell).** The
estimator is the revised D + rho^2 route; the registered Gamma route is reported only as a column.
- **V334.0 PASS (fibroblasts):**
  - D 9.79, above all 20 permuted-acrophase planes (max 7.58, mean 0.08).
  - lambda 6.56 and R 0.899 (GTEx-waveform units; timed mouse rhythms are larger than the human post-mortem waveform).
  - Myeloid cells fail (D 3.74 against null max 9.79; 4 clock UMI per cell).
- **V334.1 FAIL (decisive).** Mixing cells from samples 8 h apart lowers the coherent amplitude as expected (rho ratio
  0.489, expected 0.500), and D rises from 10.2 to 14.8. But the estimator assigns 59 % of the loss to synchrony (share_R
  0.591, delete-one-ZT 90 % 0.384-0.799; required >= 0.75 with lower bound >= 0.5). lambda falls by dlog -0.29 although
  the per-cell amplitude is unchanged by construction. At 4 h: rho ratio 0.838 (expected 0.866).
- **V334.2 PASS (simulation):** ideal von Mises cells at the skin's parameters give a median share_R of 0.009 under damping
  and 0.987 under dephasing, with no replicate on the wrong side.
- **V334.3 FAIL (reported):** lambda and R vary with time of day beyond their cell-bootstrap errors (Q 37.5 and 23.5, 5 df;
  p < 1e-6 and 0.0003). D runs from 3.4 (ZT10) to 26.6 (ZT18) and R from 0.77 (ZT14) to 0.96 (ZT10). The estimator is
  not rotation-invariant on these data, or the composition changes with time of day.
- **V334.4 PASS:** thinning to 50 % changes log lambda by 0.025 and log R by 0.032.
- **Gate: FAIL.** The aorta lambda / R decomposition is not run ("estimator not calibrated").
- **Reading.** The model-free sign held in the one controlled dephasing on data: D rose while rho fell. The von Mises split
  does not hold, and the time-of-day dependence of D shows that the skin D is not a pure clock quantity either.
- **Run note.** The first calibration run crashed in V334.3 (a negative bootstrap seed). The seed was fixed and the run
  repeated. All earlier checks reproduce identically (same seeds).

**A334b, aorta power (aligned 14 and Bmal1 KO 10 libraries; the misaligned arm was not loaded): clock check FAIL, not
decidable, not scored.**

| class | D aligned (null max) | D Bmal1 KO | KO / aligned | MDE dlog D | expected dlog D, dephasing / damping (rho 0.8) | MDE dlog rho^2 |
|---|---|---|---|---|---|---|
| smooth muscle (50,219 cells, 7.3 UMI) | 3.14 (2.27) | 8.42 | 2.68 | 0.92 | +1.36 / -0.47 | 0.16 |
| fibroblasts (7,582 cells, 5.7 UMI) | 12.8 (15.8) | 22.2 | 1.73 | 0.92 | +0.67 / -0.55 | 0.98 |

- Without a clock, the knockout spreads more than the clocked aorta. In the aorta, D is mostly covariance among clock genes
  that is not driven by the clock, so its sign cannot be read as dephasing or damping.
- D would also be underpowered: an MDE of 0.92 against an expected damping change of 0.47.
- The coherent amplitude itself would be powered (MDE 0.16 against the registered 0.446), but it is not a registered test
  on its own.

**Overall reading.** Whether a clock-loudness loss is carried by cells dephasing or by quieter cells is not decidable
with public mouse 10x data:
- the per-cell split fails a controlled dephasing on data;
- the dispersion D carries non-clock covariance that a clock knockout exposes.

A decisive design needs a clock knockout of the same tissue whose D is near zero, and the depth to resolve it (plate-based
chemistries). The text keeps "quieter" without naming the cell mechanism (PAPER_PLAN section 3).
- Code: `exploratory/A334_dephase_vs_damp/scripts/{score_a334.py --stage calibrate, score_a334b.py --stage power}`
  (1,673 s and 30 s).
<!-- src: paper/exploratory/A334_dephase_vs_damp/output/{calibration,calib_mix,calib_sim,calib_zt,b_power,b_sim,b_placebo}.csv -->

## A329 (scored 2026-11-19; BLIND for the core loop, NOT BLIND for the 14-gene values; one run): the personal departure is mostly carried by PER1, PER2 and NR1D1 (H329.1 falls: excess 2.64 -> 1.28 times matched), but it is not explained by the glucocorticoid / immediate-early signal state of the sample (H329.3 FAIL) and coupling does not rise (H329.2 FAIL): "partly leakage"

> **Erratum 2026-11-23:** the readings of H329.1 and H329.2 ("carried by PER1, PER2 and NR1D1", "coupling falls")
> are projector artefacts; see "A329 erratum" at the end of this file. The statistics stand as computed.

GTEx fast deaths, 27 peripheral tissues, A294b fit. Validation passed (A104 rho 0.734; full median log ratio 0.969).
- **H329.1 (excess) holds.**
  - The core-loop (11-gene) median log ratio to matched sets is 0.250 (90 % 0.143-0.361), against 0.969 for the 14 genes.
    So the departure goes from 2.64 to 1.28 times that of matched sets.
  - The median per-tissue change is -0.765; Wilcoxon p < 1e-7. It falls in all 27 tissues, and the core value is at most
    0.485 (the registered half) in 20 of 27.
  - Highest core values: colon sigmoid 1.475, lung 0.886. Five tissues are at or below matched sets (kidney cortex -0.344,
    breast -0.150).
- **H329.2 (coupling rises) FAIL.** Coupling falls on the core loop instead: median change -0.076, against -0.001 for the
  matched sets at the same positions; it falls in 25 of 27 tissues; one-sided p (rise) 1.0. Median coupling is 0.675 full
  and 0.605 core.
- **H329.3 (departure explained by the signal state) FAIL.**
  - The per-sample departure regressed on the glucocorticoid-target and immediate-early scores gives a median R^2 of 0.048
    (matched 0.056). The excess is -0.009 (90 % -0.013 to 0.022), and no tissue is above the matched q95.
  - On the core loop the excess is 0.001 (-0.007 to 0.030), with 1 tissue above q95 (descriptive).
- **Reading (registered): "partly leakage".**
  - Most of the clock-specific excess of the personal departure sits on the three signal-inducible genes; the 11-gene
    core loop keeps a smaller excess of about 1.3 times matched sets.
  - The departure is not the sample's measured acute signal state: the glucocorticoid and immediate-early scores do not
    explain it beyond matched genes.
  - Coupling falls when the three genes are removed. That fits their large amplitudes carrying part of the shared clock
    state, so the leakage reading's prediction of a cleaner core-loop coupling does not hold.
  - Model part 6 names which part holds: the departure is carried mainly through signal-inducible clock genes, by a
    cause that the two responsive scores do not capture. No core-loop coupling readout is added to clocksync.
- Code: `exploratory/A329_core_departure/scripts/score_a329.py --jobs 4` (7,803 s).
<!-- src: paper/exploratory/A329_core_departure/output/{tests,tissues,sets}.csv -->

## A332 (scored 2026-11-19; DESCRIPTIVE, NOT BLIND): the descriptive companions of the two-axis map: psoriasis gains clock-plane variance but not beyond its companions (D1); in brain death no organ loses loudness and three move by about an hour (D2); within organ donors loudness falls with ventilation time while the order is already lost (D3, organ-donor part only); developmental loudness rises after birth in brain and kidney, with a non-clock exception (D4)

- **D1 (skin; `exploratory/A332_descriptive/output/d1_series.csv`).** Lesion over non-lesion clock-plane variance, r (amplitude
  exp(r / 2)), two-sided swap p; matched / responsive percentiles. GSE154200 psoriasis 81 pairs r +0.499 (1.28; p 0.0004;
  94th / 90.5th) (reproduces A183's +0.50); nummular eczema 36 pairs +0.387 (1.21; p 0.10; 64.5th / 85th); atopic
  dermatitis 37 pairs +0.034 (1.02; p 0.88). A115 psoriasis lesional vs uninvolved: pooled 185 pairs +0.304 (1.16;
  p 0.0004; 56th / 70th); per series -0.142 to +0.410. **No group is outside the central 90 % of both companion families:**
  psoriasis lesions carry more clock-plane variance, as their matched and responsive genes partly do; inflamed skin moves
  time (A183, A308), its loudness change is not clock-specific. Deviations: lesional rows of the A115 series from
  `skin_lesional/processed/` and `processed_heldout/GSE13355` (`processed_a115/` holds no lesional rows; uninvolved rows
  identical); eczema batches by the A183 rules (nummular 4;6, atopic 1;3;5); companions pooled with the clock's weights,
  not tied to the swaps (as A254).
- **D2 (brain death per organ; `exploratory/A317_braindeath_confirm/output/d2_per_organ.csv`; 65 pairs).** Phase change
  organ donors minus ventilated post-mortem (tissue minus its partners, h; donor-bootstrap 90 %, not adjusted over 12
  organs): lung -1.58 (-3.11 to -0.14), whole blood -1.22 (-3.73 to 0.78), suprapubic skin -0.63, muscle -0.27, thyroid
  -0.19, visceral fat -0.21, oesophagus mucosa -0.06, atrial appendage +0.32, subcutaneous fat +0.62, lower-leg skin +0.80,
  tibial nerve +0.99 (0.17 to 1.80), tibial artery +1.09 (0.19 to 2.01). Loudness (A315, amplitude against matched):
  0.92 (tibial artery) to 1.13 (thyroid); every organ within 0.92-1.13. The arrows of Figure 5 are tangential: three organs
  move by about an hour (lung earlier, tibial artery and nerve later), none loses loudness.
- **D3 (ventilation time; organ-donor part from existing outputs).** Organ donors, top vs bottom ventilation tertile
  (100 vs 100; A132): order against fast deaths rho -0.089 vs 0.247 (both p > 0.16; `A132_brain_death/output/a132_groups.csv`);
  loudness Q* -0.236 (below every matched set; A315 `descriptive.csv`). In brain-dead donors the order is already lost in
  both tertiles, so "order held" cannot be read there; loudness is spent with ventilation time. **The ventilated
  post-mortem tertile part (79 donors) was not run** (script `d3_ventilation.py` not written in this session).
- **D4 (development; `exploratory/two_axis_map/output/d4_development.csv`).** log(V / median matched V) with the clock's
  percentile among 200 matched sets, beside C (p). Brain (Kang): >= 12 wpc 0.235 (n 19, 3 covariates; 65th; C -0.234), birth-12 y 0.460
  (75th; C 0.735, p 0.003), >= 13 y 0.683 (89.5th; C 0.662, p 0.007). C reproduces A239 for the postnatal bins only; the prenatal bin differs
  from A239's 13-22 pcw (n 14, C 0.047; -0.500 with RIN) (correction 2026-11-19, independent review). Cardoso-Moreira hindbrain:
  prenatal < 12 wpc 0.316 (68th; C 0.578, p 0.016), >= 12 wpc 0.380 (65th), birth-12 y 1.800 (99.5th; C 0.772, p 0.0005),
  >= 13 y 1.468 (98.5th; C 0.645); forebrain >= 13 y 1.893 (98th; C 0.557). Kidney: < 12 wpc -0.232 (39.5th), >= 12 wpc
  1.865 (99th; **C -0.115**), birth-12 y 1.439 (98th; C 0.719, p 0.0005). Heart prenatal 0.530 and 0.988 (79th, 90.5th;
  C -0.153, 0.058). Liver prenatal 0.173, 0.133 (63rd, 60th; C 0.148, 0.101). Gonads prenatal near 0. Reading
  (descriptive): clock-plane loudness beyond matched genes appears after birth in brain and kidney, together with clock
  shape (C); before birth it is within the matched range, except kidney from 12 wpc, where the variance is large but
  not clock-shaped (C below 0): V alone is not a clock measure, and the curve is read with C. Caveat (A240): the spread
  of sampling times differs between bins. Postnatal Cardoso-Moreira age covariate log2(years + 0.1) (A237 used an ordinal
  stage): those rows do not reproduce an A237 entry; no adult pass rate exists for these bins.
- Code: `exploratory/A332_descriptive/scripts/d1_skin_loudness.py`; `exploratory/A317_braindeath_confirm/scripts/d2_per_organ.py`;
  `exploratory/two_axis_map/d4_development.py`.
<!-- src: paper/exploratory/A332_descriptive/output/d1_series.csv; paper/exploratory/A317_braindeath_confirm/output/d2_per_organ.csv; paper/exploratory/A132_brain_death/output/a132_groups.csv; paper/exploratory/A315_braindeath_gain/output/descriptive.csv; paper/exploratory/two_axis_map/output/d4_development.csv -->

## A331 (scored 2026-11-19; BLIND for the clock statistic; powered): living women differ in a clock loudness shared by fat and skin (H331.1 PASS: leverage 0.177, above every matched set but one in twenty, 98.5th percentile), at the size the GTEx deaths predict; not beyond stress-responsive genes (93.5th percentile)

- Design as registered: MuTHER, 680 women with fat and skin from one punch biopsy; A255 states; `clocksync.elasticity.
  radius_covariance` (correlation of the two tissues' clock-plane radii after the phase harmonics are removed). Power
  stage (power note): GTEx person spread sigma_A 0.352; power 0.947 against the matched q95. One scoring run.
- **H331.1 PASS:** leverage 0.177 (twin-pair bootstrap 90 % 0.119 to 0.236; persons-permuted p 0.0005); matched sets
  median 0.026, q95 0.117 (clock at the 98.5th percentile); responsive sets median 0.121, q95 0.179 (93.5th percentile).
  The simulation at the GTEx spread expected 0.177 (5th percentile 0.116).
- Descriptive: monozygotic 0.195 (267 women), dizygotic 0.152 (413). Shared cv2 0.043 (biased low here: noisy array
  radii inflate the means; the leverage is the registered statistic).
- **Reading (registered):** persons differ in a loudness shared by two tissues in living people, beyond what expression-
  matched genes share through one biopsy and one array batch, and at the size A312 found in deaths: A312's person factor
  is not specific to deaths. Stress-responsive genes share radius almost as much (median 0.121): a shared handling or
  systemic state is part of what the two tissues co-vary in, and the clock is not separated from it here (one biopsy).
  Figure 4f: "persons differ in clock loudness shared across tissues, in deaths (A312) and in living women (A331),
  beyond expression-matched genes".
- **Limits.** One punch biopsy and one array batch for both tissues; women only; arrays; zygosity contrast descriptive.
- Code: `exploratory/P1_person_loudness/scripts/score_a331.py --stage power | score`; `clocksync.elasticity.radius_covariance`
  (bf26544).
<!-- src: paper/exploratory/P1_person_loudness/output/{power,companions,tests}.csv -->

## A330 (power stage 2026-11-19; registered re-test, NOT BLIND for the A36 baboon gains): G1 not decidable at this design: the baboon atlas has 17 passing organs but power 0.60 at a true organ-level Spearman of 0.6 (registered gate 0.8); not scored; the mouse referee (at most 8 organs, power 0.45) is not downloaded

- Design as registered (18 anatomically matched organs; baboon timed gain |g| on the GTEx waveform against GTEx sqrt S).
  The estimator reproduces the A36 `g_abs` of every organ exactly.
- **Precondition met:** 17 of 18 organs read ZT (|g| above the q95 of 2,000 ZT permutations); pituitary fails (|g| 0.81
  against q95 0.85). Residual-bootstrap SE of log |g|: 0.05 to 0.37 (median about 0.09).
- **Power:** 1.000 if the organ ranks were shared exactly (true Spearman 1); **0.604 at a true organ-level Spearman of 0.6**
  (registered criterion >= 0.8); 0.455 at 8 organs, which decides the mouse referee (Zhang 2014 maps to at most 8 GTEx
  tissues): not downloaded.
- **Reading (registered):** not decidable at this design. With 17 organs the test detects only a near-perfect shared organ
  ranking; "organ" stays "open" in Figure 4g. No cross-species correlation was computed. G1b (absolute class gains in the
  GTEx single-nucleus atlas) remains a scout.
- Code: `exploratory/G1_organ_gain_species/scripts/score_a330.py --stage power`.
<!-- src: paper/exploratory/G1_organ_gain_species/output/{organs,power}.csv -->

## A327 (scored 2026-11-19; descriptive, NOT BLIND for A317-A318, A325; underpowered): the brain-death contrast is positive in every perfusion class, leg tissues included (D_rho N 0.244, P 0.265, mixed 0.456), with no perfusion dependence (Delta -0.173, p 0.66); not decided at MDE 1.0

- Perfusion split as an unverified anatomical hypothesis (registration): possibly reached atrial appendage, lung,
  visceral fat, oesophagus mucosa; not reached tibial artery, tibial nerve, lower-leg skin, gastrocnemius, leg
  subcutaneous fat; thyroid and suprapubic skin excluded; blood apart. A317 data, groups, draws, estimator.
  | class | pairs | rho post-mortem | rho organ donors | D_rho |
  |---|---|---|---|---|
  | P (both possibly reached) | 5 | 0.925 | 0.660 | 0.265 |
  | N (both not reached, leg) | 10 | 0.650 | 0.406 | 0.244 |
  | X (mixed) | 20 | 0.858 | 0.402 | 0.456 |
  | P u X | 25 | 0.871 | 0.455 | 0.416 |
- Delta = D_rho(N) - D_rho(P u X) = -0.173, two-sided relabelling p 0.664 (MDE 0.997). D_rho(N) 0.244, p 0.253
  (MDE 0.850).
- **Reading as registered: neither** (no perfusion dependence detected, no significant leg-only contrast): not decided
  at this power. Descriptively the organ donors' arrangement departs from the post-mortem one in each class, including
  the leg tissues that no flush reaches; nothing points to a perfusion-only effect, but 10 leg pairs cannot show it.
- **Limits.** Perfusion and recovery order are aliased (heart and lung are perfused and recovered first); the class
  sizes (5, 10, 20 pairs) leave the test underpowered.
- Code: `exploratory/A317_braindeath_confirm/scripts/score_a327.py --stage power | score`.
<!-- src: paper/exploratory/A317_braindeath_confirm/output/{a327_power,a327_tests}.csv -->

## A326 (scored 2026-11-19; NOT BLIND for the human results and the mouse ratio): the human lymphocyte designs could not have read a mouse-colon-sized T-cell clock (H326.1 FAIL in A281, A287, A321 at r = 0.25: design power 0, 0, 0.05); they read a clock at 0.70 of the fibroblast amplitude (A281 0.996, A321 1.00; A287 only at 1.0): "no readable clock" bounds the lymphocyte clock below about 0.7 of the fibroblasts', not below the mouse ratio

- **Mouse ratio (L1 scout).** Colon chrono-atlas, undamaged arm (8 mice, ZT 1/7/13/19, 14 clock genes, debiased cosinor
  amplitude): T cells 1.12, fibroblasts 4.55 log2-CPM, ratio 0.246 (bootstrap q95 0.404); damaged arms unresolved
  (R-DC 0, q95 0.93; R-PC 0.39, q95 0.89). The 2026 small-intestine immune atlas (CD45+ only, no processed matrix found)
  has no stromal reference: not used.
- **Insertion.** Partner fibroblast-stromal rows, clock-plane states scaled by r (`clocksync.simulate.scale_clock`),
  thinned to the lymphocyte clock-gene count, each design's own statistic; 20 replicates per unit and r.
- **Design power** (A281: P(> 2 of 7 powered atlases pass); A287: P(>= 1 of 3 powered); A321: Stouffer over 3
  informative compartments):
  | r | A281 | A287 | A321 |
  |---|---|---|---|
  | 0.10 | 0.000 | 0.000 | 0.00 |
  | 0.25 | 0.000 | 0.000 | 0.05 |
  | 0.40 | 0.003 | 0.000 | 0.15 |
  | 0.55 | 0.450 | 0.000 | 0.75 |
  | 0.70 | 0.996 | 0.288 | 1.00 |
  | 1.00 | 1.000 | 1.000 | 1.00 |
  Per-unit pass rates at r = 0.25: at most 0.10 in any powered or informative unit (`pass_rates.csv`).
- **Reading as registered.** H326.1 fails in every design at r = 0.25 and at the mouse upper bound 0.40. Bounds x
  (smallest grid r with power >= 0.8): A281 0.70, A321 0.70, A287 > 0.70 (1.0). Wording: "tissue lymphocytes carry no
  clock readable between donors at depths where fibroblast clocks are read; such a clock would have to be weaker than
  about 0.7 of the fibroblast amplitude". A mouse-colon-sized T-cell clock (a quarter of the fibroblasts') would not
  have been seen: the human negative is **not** evidence of a species, residency or tissue difference.
- **Limits.** Scaling the clock also scales its in-plane biological noise (detection rates are lower bounds; the true
  bound may be somewhat lower than 0.70); the inserted clock carries the fibroblasts' non-clock structure and host
  coupling; the mouse ratio is one colon atlas of 8 mice.
- **Code error found and fixed before reporting:** the first scoring run keyed the registered powered-atlas rates by
  atlas name only, so A281 took A287's rates (3 instead of 7 powered atlases); per-unit pass rates were unaffected; the
  reported design powers are from the rerun with the fixed key (`design_power_run1_keybug.csv` kept, not used).
- Validation (r = 1) and its note: registry, A326 validation note.
- Code: `exploratory/L1_lymphocyte_power/scripts/score_a326.py --stage validate | score`.
<!-- src: paper/exploratory/L1_lymphocyte_power/output/{validation,validation_diag_200reps,pass_rates,design_power,replicates}.csv; scout/output/amplitude_ratios.csv -->

## A319 (scored 2026-11-19; NOT BLIND for the published tables): in forebrain Bmal1 KO mice the organs' loss of "within-mouse coherence" is mostly a new, reproducible organ order, not independent drift (H319.2 PASS: under restricted feeding all 6 organs cluster across mice and the new order accounts for all of the within-mouse spread); without any shared drive (DD) no significant shared phase is left beyond the organ order (H319.1 PASS, E 0.069, p 0.070, borderline)

- Data: Izumo et al. 2014 eLife source data (organ means and circular variances across mice; per-mouse circular
  variance over six organs); 11-15 mice per arm; one run, independence null 20000 replicates.
- **Order share** (V_order / mean within-mouse V): WT 0.90 (LD), 0.83 (DD), 0.31 (FR, CT frame); BKO 0.91 (LD), 0.79
  (DD), 1.09 (FR). Mean within-mouse V: WT 0.173 / 0.129 / 0.143; BKO 0.226 / 0.405 / 0.512.
- **H319.2 PASS (BKO FR).** 6/6 organs cluster across mice (Rayleigh p <= 0.005); V_order 0.558 against observed
  0.512: the high within-mouse variance the paper reads as lost coherence is the spread of a consistent new order
  (liver 150 deg, kidney 72, heart 46, lung 9, spleen 73, pituitary 227; shifts from WT FR 50-150 deg; circular
  correlation of organ means with WT 0.02). Excess coherence over independence E 0.076, p 0.035: a mouse's organs
  also share part of their phase beyond the group order (mouse-specific entrainment).
- **H319.1 PASS (BKO DD).** E 0.069, p 0.070: not significant, at the edge; 6/6 organs still cluster across mice
  (liver weakest, R 0.445, p 0.049) and the order explains 0.79 of the within-mouse spread; liver, kidney and heart
  are advanced 47-56 deg relative to WT DD, lung and spleen not (correlation of organ means with WT 0.91).
- **LD:** BKO keeps the WT order (shifts -10 to +19 deg, correlation 0.96), as the paper reports.
- **Amplitude (quoted, not reanalysed):** lower in BKO DD, kept in BKO LD, restored by FR only in liver; heart
  imaging: lower cell amplitude and cell dephasing together.
- **Reading (registered).** The mouse shows the human pattern for order: without the central clock, a substitute
  drive (feeding) sets a new, reproducible organ order rather than leaving organs independent; with light, the order
  is kept. Without any drive for > 30 days, organs keep a reproducible order across mice (possibly the dissection reset
  or a residual cue) and only a borderline shared phase within the mouse. Tension for amplitude: tissue amplitude
  falls after weeks without a shared drive (cell dephasing), which A315's "amplitude kept" in brain-dead donors (hours
  to days, ICU cues) does not contradict but does not generalise to: the human claim is stated for the short term.
- Code: `exploratory/A319_izumo_mouse/scripts/score_a319.py`; card `manuscript/literature/notes/izumo2014differential.md`.
<!-- src: paper/exploratory/A319_izumo_mouse/output/{conditions,organs,tests}.csv -->

## Gap 3 note (2026-11-19; post hoc development on the A163 validation rows, no new simulation): a host-reliability threshold does not rescue the c intervals; c intervals are stated as approximate

- Restricting the df-corrected intervals to hosts with estimated reliability R >= r0 lowers coverage instead of
  raising it (cells with coverage >= 0.90: 57/64 at r0 0, 51/64 at 0.8, 44/60 at 0.85, 28/41 at 0.9; minimum coverage
  0.802 -> 0.000): replicates whose R is overestimated are under-corrected, so selecting on R selects the biased ones.
  An observable domain rule is not available; a different interval (e.g. bias-corrected bootstrap near the boundary)
  would need bootstrap replicates that A163 did not store and a fresh validation (about 19 h at 16 workers).
- **Statement for the paper (gap 3 closed as "state").** Absolute c is reported df-corrected with donor-bootstrap
  intervals labelled approximate: in simulation at the paper's designs (A163), coverage 0.80-0.98 (median 0.94),
  57 of 64 cells >= 0.90; undercoverage is confined to c >= 0.6 (width near the boundary at c 0.9; bias with 30 donors
  and 2 anchors at anchor coupling 0.6). Claims rest on contrasts at a common host (d, d*; calibrated, A154 V154.2), on
  gain and departure, and on ranks; none rests on an absolute c interval.
- Code: `exploratory/A163_reliability_correction/scripts/dev_reliability_rule.py`.
<!-- src: paper/exploratory/A163_reliability_correction/output/validate/{summary,dev_reliability_rule}.csv -->

## A320 (scored 2026-11-19; calibration, NOT BLIND for the clock statistics): a broad stress-responsive family (376 genes) is calibrated at every design (V320 PASS x8) and the living internal-time results are clock-specific against it (H320 PASS x5; A138 PASS); but tumour immune pairing also passes it, so it is less stringent than the glucocorticoid / immediate-early family against a hormonal pulse

- Family: Hallmark NFkB signalling and hypoxia (MSigDB 2024.1) plus the GC targets and immediate-early genes, minus the
  clock genes; sets drawn by the responsive rule. Mean pairwise Jaccard overlap 0.05-0.13 (narrow family 0.65-1.00;
  matched 0.001). Placebo as A285 (2000 draws); draw 0 reproduces every scored statistic and the matched and narrow
  q95. Runs: A251 and A255 under a minute each, A138 41 s, A227 1290 s (peak 31.2 GB).
- **V320 PASS in all 8 tests:** q95-rule rejection 0.0495-0.0635 (band 0.0646); clock above every set under no pairing
  in 0.35-0.85 % of draws (bound 2 %). The narrow family in the same runs: 0.076-0.204 and 1.0-8.4 %.
- **H320 (p_cal against the broad family):**

| test | clock | broad q95 | p_cal broad | p_cal narrow (A285) | p_cal matched |
|---|---|---|---|---|---|
| A251 H251.1 skin sites | 0.628 | 0.379 | 0.0075 | 0.075 | 0.004 |
| A251 H251.2a | 0.813 | 0.518 | 0.0090 | 0.084 | 0.0065 |
| A251 H251.2b | 0.594 | 0.276 | 0.0045 | 0.0635 | 0.0075 |
| A255 MuTHER fat-skin | 0.533 | 0.283 | 0.0065 | 0.0565 | 0.0105 |
| A255 surgical pooled | 0.686 | 0.449 | 0.0060 | 0.017 | 0.0045 |
| A138 H138.4 tumour-normal | 0.365 | 0.246 | 0.0055 | 0.011 | 0.0040 |
| A227 H227.1 tumour myeloid (descriptive) | 0.581 | 0.549 | 0.0040 | 0.439 | 0.0185 |
| A227 H227.2 tumour T/NK (descriptive) | 0.539 | 0.512 | 0.0080 | 0.200 | 0.0350 |

- **Reading (registered) and its qualification (post hoc, from the A227 rows).** By the registered reading, the living
  internal-time results (skin sites, fat-skin, surgical pairs) and the tumour-normal pairing are clock-specific against
  a calibrated stress-responsive family. But the tumour immune pairings, which the narrow family does not support
  (A271), pass the broad family too: diluting the glucocorticoid and immediate-early genes among NFkB and hypoxia genes
  makes the family calibrated and less stringent against a hormonal pulse. The narrow family stays the stringent
  control and is not decidable at one to six pairs. Wording for the paper: "clock-specific against expression-matched
  genes and against a calibrated stress-responsive family; against glucocorticoid and immediate-early genes alone,
  decided only where enough pairs exist (surgical pairs, tumour-normal)". Gap 2 is closed with this wording.
- Code: `exploratory/A285_companion_calibration/scripts/placebo_a285.py --families matched,responsive,responsive_broad
  --out a320`; clocksync c2203c3 (`responsive_broad`).
<!-- src: paper/exploratory/A285_companion_calibration/output/a320/{a251,a255,a138,a227}/{repro,tests}.csv, {a251,a255}/overlap.csv -->

## A322 (scored 2026-11-19; NOT BLIND for the shifts and the liver amplitude): the "cells split between drives" chord does not fit liver, WAT and lung together (H322.1 rejected, Q 31.9, p 0.0005); conflicting feeding costs amplitude in all three tissues, but lung (and slightly liver) keep more than the mix predicts at their shift

- Manella Alb-Cre controls, daytime feeding vs ad libitum, 24 + 24 libraries per tissue, 12 clock genes after QC; one
  projector per tissue over both arms; library bootstrap within ZT (2000). One run.
- **Complex ratios DF / AL:** liver amplitude 0.712 (90 % 0.660-0.768), shift 11.69 h; WAT 0.319 (0.272-0.371),
  10.60 h; lung 0.495 (0.407-0.588), 1.90 h.
- **Fit** (shared Delta, one weight per tissue): Delta 11.28 h (11.09-11.48); b liver 0.844, WAT 0.652, lung 0.371.
  Model points: liver 0.692 at 11.12 h (Mahalanobis 5.8), WAT 0.316 at 10.49 h (0.08), lung 0.274 at 0.98 h (26.0).
  Q 31.86; chi-square (2 df) p < 0.001; parametric bootstrap p 0.0005. **H322.1: rejected.**
- **Reading (registered: rejection with |r| above the chord).** The tissue moved least (lung) keeps about half its
  amplitude where splitting its cells between the two drives would leave about a quarter: self-sustained cells
  re-synchronise at an intermediate phase, so conflict costs less amplitude than the mix. WAT, whose weight is near
  one half, loses most (to a third), as the mix predicts. Federation part 5 is reworded: conflicting drives move a
  tissue by its feeding weight and cost amplitude, most where the weights are balanced, less than a pure mix of
  drives.
- **Limits.** Three tissues, one design, one laboratory; tissues bootstrapped separately though from the same mice; the
  DF arm is not exactly 12 h from AL (Delta fitted at 11.3 h).
- Code: `exploratory/A322_conflicting_drives/scripts/score_a322.py`.
<!-- src: paper/exploratory/A322_conflicting_drives/output/{ratios,fit,tests}.csv -->

## A321 (scored 2026-11-19; BLIND for the atlases): on four fresh atlases tissue lymphocytes show no core-loop clock (H321.1 FAIL with 3 informative compartments, combined p 0.47): "no readable clock" holds with the core loop included; one uninformative exception (skin T cells)

- One scoring run (600 s, 6 workers, peak 3.1 GB). Six lymphocyte compartments, 15-21 shared donors.
- **Primary compartments (core loop C, p; power control; percentile among 200 matched 11-gene sets):**
  - lung NK (Natri controls): C -0.357, p 0.78; power 0.80 (informative); 27.5th percentile.
  - liver CD8 T: C 0.081, p 0.35; power 1.00; 44th.
  - adipose T: C 0.172, p 0.30; power 1.00; 65.5th.
  - skin T: C 0.604, p 0.012; power 0.00 (thinned fibroblasts never pass: not informative); 95th percentile.
- **H321.1 FAIL:** Stouffer over the three informative primary compartments z 0.087, p 0.465. H321.2 not read.
- Descriptive (secondary compartments): lung T C -0.223 (p 0.69, power 0), liver NK -0.394 (p 0.84, power 1.00). The
  14-gene C: negative or near 0 in five of six (skin T 0.469, p 0.010).
- **Reading (registered).** With three informative compartments, the core loop does not reveal a lymphocyte clock in
  fresh tissues: A294's post hoc core-loop pattern does not replicate, and the wording "tissue lymphocytes carry no
  readable clock" holds with the core loop included. Skin T cells cohere on both gene sets (95th percentile of matched
  sets), but the skin's own fibroblasts, thinned to the T cells' counts, never pass, so the test cannot separate a
  lymphocyte clock from other donor-shared covariation there: an exception to check, not a result.
- Code: `exploratory/A321_lymphocyte_core_blind/scripts/score_a321.py --stage score`.
<!-- src: paper/exploratory/A321_lymphocyte_core_blind/output/{design,compartments,tests}.csv -->

## A323 (scored 2026-11-19; NOT BLIND for the GTEx gains and class gains): composition and conserved class gains do not predict the GTEx organ gains (H323.1 FAIL, rho 0.52, p 0.098, 8 tissues); the composition-predicted range (1.29-fold) is far narrower than the observed one (2.89-fold): an organ-level factor scales all cell classes together

- Reference: Eraslan 2022 GTEx single-nucleus atlas (Census 4ed927e9; 209,126 nuclei, 16 donors, 8 tissues),
  clock-gene pseudobulk per donor x tissue x cell type built 2026-11-19 (Census symbol BMAL1 for ARNTL). One run.
- **Code error found and fixed before the reported run:** the first run gave epithelial cells no weight (the family
  is matched by keywords and was missing from the gain map); the registered g (-0.220) is used in the reported run,
  with an assertion that every family has a gain. The first run's numbers are not used.
- **Clock-mRNA shares:** epithelial 0.64-0.74 in breast, lung, esophagus mucosa, prostate, skin; fibroblast-stromal
  0.45-0.49 in heart and esophagus muscularis (mural 0.37 there); parenchyma 0.73 in skeletal muscle, 0.38 in heart;
  lymphocytes <= 0.044 everywhere.
- **Predicted G vs observed sqrt S:** heart LV 1.047 / 2.677; esophagus muscularis 1.062 / 2.995; skeletal muscle
  0.958 / 1.891; breast 0.877 / 1.620; esophagus mucosa 0.860 / 2.768; skin (leg) 0.852 / 1.037; prostate 0.848 /
  2.048; lung 0.824 / 1.760.
- **H323.1 FAIL:** Spearman 0.524, exact one-sided p 0.098. Log R^2 0.35, Theil-Sen slope 1.92. Ranges: predicted
  1.29-fold, observed 2.89-fold. Sensitivities: "other" renormalised 0.381 (p 0.18); nucleus shares 0.381 (p 0.18).
- **Reading (registered FAIL) and the range (descriptive).** Composition does not predict organ gain over these 8
  tissues. Because the conserved class gains (A293: within a factor of about 1.4 against endothelium) and the measured
  compositions can only spread organ gains 1.3-fold, while the organs spread 2.9-fold (esophagus mucosa and skin, both
  epithelial, differ 2.7-fold), an organ-level factor scales every class of an organ together: the A298 alternative
  ("organs drive their cells at different strengths"). Within a tissue, composition still moves gain (A295,
  stromal-vascular +0.041 per SD). Model part 2 is reworded accordingly.
- **Limits.** 8 tissues (low power for the rank test); class gains from other atlases; snRNA nuclei, not cells.
- **Wording note (2026-11-19, REVIEW_2026-11-19 section 7 item 5).** The result reads "composition and conserved class
  gains do not explain the organ gain range". The "organ-level factor" above is an inference by elimination over 8
  tissues; a tissue-specific technical attenuation in GTEx would produce the same pattern. It leaves the results until
  check G1 (organ gains against timed clock amplitudes in other species) separates the two.
- Code: `exploratory/A323_gain_decomposition/scripts/score_a323.py`.
<!-- src: paper/exploratory/A323_gain_decomposition/output/{shares,tissues,tests}.csv -->

## A324 (scored 2026-11-19; BLIND for heart, kidney, VAT): in Xin 2021 the pure mix of drives is rejected again (H324.1 PASS, Q 185, p 0.0005), but not every departure is above it (H324.2 FAIL: visceral fat loses slightly more): conflict costs amplitude; the tissues that barely move keep half or more, adipose (balanced weights) loses most in both datasets

- Xin 2021, DRF vs NRF, 28 + 28 mice per tissue, 14 clock genes; A322's estimation and fit; one run.
- **Complex ratios DRF / NRF:** liver 0.438 (0.413-0.465), 11.29 h; heart 0.589 (0.548-0.638), 1.59 h; kidney 0.648
  (0.604-0.697), 2.49 h; VAT 0.111 (0.091-0.133), 23.52 h (-0.48 h).
- **Fit:** Delta 11.42 h; b liver 0.718, heart 0.233, kidney 0.278, VAT 0.447. Model amplitudes: liver 0.441
  (Mahalanobis 1.2), heart 0.537 (16.8, above), kidney 0.449 (148.5, above), VAT 0.131 (18.9, below).
- **H324.1 PASS** (rejected again: parametric p 0.0005). **H324.2 FAIL:** VAT lies below the chord.
- **Reading (registered: H324.1 passes, H324.2 fails).** Conflict can cost more amplitude than a pure mix (extra
  dephasing in visceral fat); model part 5 is worded "conflicting drives cost amplitude". Descriptive, consistent over
  both datasets (Manella A322, Xin): the tissues that follow the feeding drive least (lung, heart, kidney) shift 1.6 to
  2.5 h and keep half to two thirds of their amplitude; the liver follows almost fully (11.3-11.7 h) and keeps 0.44 to
  0.71; adipose, the tissue with the most balanced weights, loses most (WAT 0.32, VAT 0.11). Xin VAT does not match
  Manella WAT (0.11 vs 0.32; different depots and protocols).
- Code: `exploratory/A322_conflicting_drives/scripts/score_a324.py`.
<!-- src: paper/exploratory/A322_conflicting_drives/output/xin/{ratios,fit,tests}.csv -->

## A328 part a, power stage (2026-11-21; no del Valle clock statistic computed): both del Valle readouts are estimable
- Adult stand-ins at the del Valle design (32 fetal adrenals, stage and sex, permutations within stage;
  `A328_fetal_running_loop/output/power_a/passrates.csv`): pass rate at the full adult clock 1.000 for dv:all and
  dv:core; at 0.75 1.000 for both; at 0.5 0.615 (dv:all) and 0.850 (dv:core). False-positive rate 0.085 (dv:all, CI
  0.050-0.133) and 0.070 (dv:core).
- A FAIL would exclude an adult-size covariation down to 0.75 (14 genes) and 0.5 (core loop). Scoring (one run) waits
  for the user's go.
- Delivery QC (A328 note 1, `output/qc/`): design as registered; Cardoso-Moreira exon counts agree with the authors'
  matrix (median Spearman 0.979; the 7 non-matching libraries are those the authors excluded).

## A328 part a (scored 2026-11-21; BLIND): the fetal adrenal covariation replicates blind in del Valle 2023, and the core loop covaries on its own (H328.1 PASS, H328.2 PASS)
- 32 whole fetal adrenals, 7-11.5 wpc, stage and sex removed, 10,000 permutations within stage
  (`A328_fetal_running_loop/output/score_a/{l1,verdicts}.csv`):

  | test | C | p | matched q95 | responsive q95 | leave-one-gene-out | verdict |
  |---|---|---|---|---|---|---|
  | H328.1 dv:all (14 genes) | 0.592 | 0.0021 | 0.381 | 0.409 | 0.492 (NPAS2) to 0.786 | PASS |
  | H328.2 dv:core (11 genes) | 0.803 | 0.0002 | 0.470 | 0.605 | 0.683 (ARNTL) to 0.841 | PASS |
  | dv:-PER1-NR1D1 (descriptive) | 0.736 | 0.0011 | 0.492 | 0.485 | 0.553 to 0.813 | clock-specific |
  | dv:early, 7-8 wpc (16) | 0.645 | 0.0016 | 0.422 | 0.448 | 0.572 to 0.754 | clock-specific |
  | dv:late, 9-11.5 wpc (16) | 0.450 | 0.025 | 0.409 | 0.363 | 0.231 to 0.592 | clock-specific |
- Both p are below 0.01, so the stand-in false-positive rate (0.085 for dv:all) does not qualify them.
- **Reading (registered).** H328.1: the fetal adrenal covariation replicates blind in a tenth laboratory, before 12
  weeks. H328.2: the core loop covaries without the signal-induced genes, clock-specific against both companions:
  "an entrained loop is supported". This is against the stated expectation (H328.2 FAIL, from A266 on the
  single-cell cortex, where the core loop was not significant) and against claim 3(b) of `PAPER_PLAN_2026-11-19.md`
  ("it rests on the signal-induced genes"). In del Valle the core loop is more coherent than the 14 genes (0.803 vs
  0.592). The cause of the difference with A266 is open (whole bulk adrenal at 32 fetuses vs single-cell cortex
  pseudobulks of fewer fetuses per cohort; A266 power at the adult size).
- Figure 3: "replicated" stands for the adrenal; the "driven state" wording of 3(b) / 3e needs the PI's review
  (Q to `PI_QUESTIONS.md`). Part b (running-loop test) still waits for the 10x deliveries.

## A337 (scored 2026-11-21; BLIND for the unspliced counts): in timed tissues the untimed lead statistic finds the right sign but recovers a lead in only 3 of 7 rhythmic compartments (H337.1 below the registered 2/3); where it does, it equals the timed lead (H337.2 3/3); the responsive condition of H328.3 is met in none (H337.3 0/7)
- Baboon (12 per tissue) and Greenwell liver (AL + NR, 36), 14 genes, subtraction introns
  (`A337_timed_lead/output/baboon_greenwell/{verdicts,summary,compartments}.csv`). All 7 primary compartments are
  rhythmic (spliced state, p <= 0.02).

  | compartment | gate (Holm) | delta_U h (90 %) | delta_T h (90 %) | H337.1 | H337.2 |
  |---|---|---|---|---|---|
  | baboon adrenal cortex | 0.001 | 0.34 (-0.23 to 1.02) | 0.61 (0.02 to 1.26) | fail (interval spans 0) | - |
  | baboon heart | 0.001 | 1.09 (0.42 to 1.76) | 1.27 (0.55 to 1.96) | PASS | PASS |
  | baboon kidney cortex | 0.060 | 3.55 (0.32 to 6.76) | 3.00 (0.81 to 4.94) | fail (gate) | - |
  | baboon kidney medulla | 0.018 | -0.08 (-2.22 to 1.16) | 2.03 (0.12 to 3.63) | fail (interval spans 0) | - |
  | baboon liver | 0.001 | 0.45 (-0.03 to 1.07) | 0.92 (0.21 to 1.95) | fail (interval spans 0) | - |
  | baboon lung | 0.001 | 0.48 (0.08 to 0.98) | 0.95 (0.46 to 1.45) | PASS | PASS |
  | Greenwell liver AL + NR | 0.001 | 2.35 (1.80 to 2.84) | 2.17 (1.63 to 2.67) | PASS | PASS |

- The timed lead is positive in all 7 (median 1.27 h, range 0.61-3.00); the untimed point estimate in 6 of 7. The
  failures are precision at n = 12 (three intervals span 0, one gate at Holm 0.060), not a wrong sign.
- Direct introns: H337.1 4/7 (liver passes, interval 0.09-1.16); core loop: 3/7 on both definitions.
- H337.3: the lead never exceeds the responsive companions' q95 on the 14 genes (0/7; q95 0.67-5.66 h); on the core
  loop it does in lung and Greenwell (2/7).
- Descriptive: adrenal medulla delta_U 0.63 h (0.22-0.99), timed 0.84; Greenwell arrhythmic arms keep a lead (AR 2.87 h,
  DR 2.24 h; timed 3.50 and 2.34): the liver clock runs under arrhythmic feeding. Per-gene timed lag of intron over
  exon (genes with spliced amplitude > 0.15, threshold set after the run): median 1.02 h, positive in 69 % of 96
  gene x compartment pairs.
- **Readings (registered), carried to A328 part b:**
  - H337.1 fails in most rhythmic compartments, mostly with the gate passing and delta_U spanning 0: the statistic lacks
    precision at n = 12; a fetal "no consistent lead" is reported as uninformative, not "driven". It recovers the lead
    at n = 36 (Greenwell) and where the clock is strong (heart, lung).
  - H337.2 PASS where defined: the expected lead for A328 is the median timed lead, 1.27 h (0.61-3.00); it replaces the
    literature midpoint of +2 h in A328's fallback power.
  - H337.3 0/3 among the H337.1 passes: the responsive-set condition of H328.3 misses running clocks; A328 reports
    H328.3 with and without that condition (registered criterion primary; the second reading set by A337 before any
    fetal lead was computed).
- Process: the first verdict table counted the adrenal medulla (descriptive in the registration) among the primary
  compartments; verdicts were re-derived from the stored statistics with the registered set (`score_a337.py --reread`),
  no statistic recomputed.

## A328 part b and A338 (scored 2026-11-22; BLIND for every unspliced count, lead and exon / intron ratio; one run): no consistent pre-mRNA lead in the two decidable compartments (H328.3 FAIL: del Valle, liver), uninformative rather than "driven": the running postnatal liver clock shows the same picture and the real intervals are 1.6-2.9 times wider than the power model's; A338 timescale not resolved in liver
- Design as registered with the power-note decision (matched estimator primary, pergene companion, plain beside;
  scored: del Valle and prenatal liver, Holm over the two; 10,000 permutations, 2,000 bootstrap resamples, 200 + 200
  companions). Code check first (`score_a328b.py --check`, unspliced := thinned spliced: no H328.3 pass, leads -0.72 to
  +0.60 h). Run 287 s (`A328_fetal_running_loop/output/score_b/{compartments,verdicts,companions,genes,a338}.csv`).

  | unit (n) | genes | gate p (Holm), matched | lead h, matched (90 %) | lead h, pergene (90 %) | lead h, plain (90 %) | H328.3 |
  |---|---|---|---|---|---|---|
  | del Valle (32) | 14 | 0.0017 (0.0017) | 0.90 (-0.78 to 2.86) | 0.54 (-1.44 to 1.88) | 0.30 (-1.26 to 1.92) | FAIL |
  | del Valle (32) | core | 0.0001 (0.0002) | 0.04 (-1.10 to 0.98) | -0.11 (-1.26 to 1.32) | -0.03 (-0.86 to 0.83) | FAIL |
  | liver (37) | 14 | 0.0001 (0.0002) | 0.60 (-0.11 to 1.08) | -0.15 (-3.13 to 3.24), gate 0.22 | 0.79 (0.08 to 1.27) | FAIL |
  | liver (37) | core | 0.13 (0.13) | 1.48 (-0.48 to 5.13) | -1.66 (-8.22 to 1.09) | 1.12 (-0.36 to 2.85) | gate fails: uninformative |

  Direct introns agree (del Valle 14 genes 1.19, -0.86 to 3.90; core -0.05, -1.18 to 0.99; liver 14 genes 0.51, -0.49
  to 1.28). With the responsive condition: FAIL in all (responsive q95 0.63-1.37 h).
- *Running-clock reference, same estimator:* postnatal liver (12 libraries, ages 0.04-58 y) matched 0.33 h (-0.30 to
  0.78), gate p 0.0001; pergene 0.47 (0.01 to 1.15); postnatal heart (10) 0.79 (-0.00 to 1.45); postnatal kidney (9)
  gate fails. The statistic does not resolve the lead of a running postnatal clock at n = 10-12 either.
- *Precision against the power stage:* observed matched interval widths del Valle 3.64 h (14 genes) and 2.08 h
  (core), liver 1.18 h, against the GTEx-model medians 1.27, 1.29 and 0.74 h: the model, which shares all biological
  noise between the two count types, overstated precision by 1.6-2.9 times, as the power note anticipated.
- *The gate is not clock-specific:* expression-matched pseudo-clocks pass it in 46-92 % of sets and responsive sets in
  74-100 % (`companions.csv`); spliced and unspliced share any between-fetus variation. Its passes in heart (matched
  0.34 h, 0.005 to 0.84) and kidney (0.12 h, -0.29 to 0.50), the expected negatives, carry no clock meaning.
- *A338 (liver, 14 genes, matched):* phi_f 0.60 h (-0.11 to 1.08); phi_ref = postnatal liver 0.33 h (its gate
  passed); lifetime ratio r 0.81 (0.64 to 1.05; 10 genes); T_eff 10.7 h, 90 % interval unbounded. Registered reading:
  "no daily-scale input of the reference size" (phi_f interval includes 0, compartment decidable). Post hoc caution:
  the same interval also contains the reference lead (0.33 h), and the postnatal reference is itself not resolved, so
  in substance steady and daily inputs both remain compatible: the timescale is not resolved.
- *Fast side read alone (user, 2026-11-22; registered threshold and quantity, read one-sided after scoring):* A338
  fixed "T_eff entirely below 6 h: a fast, pulsatile input" before scoring. With the registered postnatal reference
  (0.33 h, itself unresolved: -0.30 to 0.78) the liver lead's upper bound (1.08 h) gives T_eff >= about 5-6 h at the
  reference's point value (r 0.64-0.81) and no lower bound once the reference's uncertainty is included (tau could be
  near 0). With the literature lifetime registered in A339 (tau = 1 h), T_eff = 2 pi tau / tan(phi): liver >= about
  22 h, del Valle core >= about 24 h, del Valle 14 genes >= about 7 h; allowing for the matched estimator's under-reading
  of the lead (about 2-fold on Greenwell): >= about 10, 11 and 3.6 h. Reading: the pre-specified fast-input reading is
  not reached; **ultradian (hourly) pulses are excluded in liver and the adrenal core, conditional on fetal clock-mRNA
  lifetimes near the adult 1-2 h**; the postnatal calibration is too imprecise to fix tau. A combination of registered
  elements (A338 threshold, A339 lifetime), labelled as such.
- *Descriptive (not decidable):* 10x cortex matched 0.86 h (0.33 to 1.62; core 1.64, 0.81 to 5.41), `cortex_minus1`
  1.36 (0.38 to 2.17); pergene intervals span 0 (14 genes 0.72, -1.68 to 1.81); below the responsive q95 (1.28).
  Per-gene lags (pergene, 14 genes) are erratic at these depths: median del Valle 0.72 h, liver 0.20, postnatal liver
  0.36; positive in 57 %, 50 % and 71 % of genes.
- **Readings.**
  - H328.3 FAILS in both decidable compartments. By A337's carried reading ("a fetal 'no consistent lead' is reported
    as uninformative, not 'driven'"), and because the running postnatal liver clock gives the same unresolved lead with
    the same estimator, the result is **uninformative about loop vs drive**, not evidence for a driven state. The
    registered "with power" branch ("a driven state; a running loop of the reference's size excluded") does not apply:
    the power stage overstated precision (above) and the reference itself is not resolved.
  - Figure 3e stays "the open test" for the running-loop question; claim 3 parts 1, 2 and 4 of the proposed framing
    (Q43) are unaffected; part 3 (daily timescale) is not established by A328 / A338.
  - Not done (registered descriptive): CoPhaser phases on the spliced / unspliced layers, the cycling-share covariate
    sensitivity, the damping rho.

## A339 (scored 2026-11-22; BLIND for rho; one run): fast input not excluded (H339.1, H339.2 FAIL), uninformative: the running postnatal clocks also read rho well below 1, so extra-Poisson intronic variation biases rho down beyond what the power stage captured
- Matched estimator, 20 thinning replicates, exact log-count noise term (power note), 2,000 bootstrap resamples
  (`A328_fetal_running_loop/output/score_a339/{rho,verdicts}.csv`). Code check first (unspliced := thinned spliced:
  rho 0.73-1.85, above 1 where the shared counting noise is over-subtracted, as expected for that check).

  | unit (n) | genes | rho (90 %) | unspliced noise share | H339.1 (> 0.43) Holm p | H339.2 (> 0.69) Holm p |
  |---|---|---|---|---|---|
  | del Valle (32) | 14 | 2.29 (0.87 to 13.54) | 0.85 | 0.36 FAIL | 0.37 FAIL |
  | del Valle (32) | core | 0.63 (0.41 to 1.07) | 0.36 | 0.16 FAIL | - |
  | liver (37) | 14 | 0.45 (0.28 to 0.72) | 0.26 | 0.42 FAIL | 0.94 FAIL |
  | liver (37) | core | 0.45 (0.18 to 0.85) | 0.43 | 0.43 FAIL | - |

- *Running-clock references, same estimator:* postnatal liver 0.74 (0.58 to 0.92), heart 0.89 (0.79 to 1.01), kidney
  0.53 (0.23 to 1.36), against about 0.97 expected for a daily input at tau = 1 h; the Greenwell check of the power stage
  predicted about 0.9 at fetal depths. Descriptive fetal: 10x cortex 0.85 (0.49 to 1.39), heart 0.57, kidney 0.95.
- **Reading (registered):** "a fast component is not excluded"; not "the input is fast" (that reading needed the
  postnatal references near 1, and they are not). Post hoc: real unspliced counts carry variation beyond counting noise
  (intron retention, nascent-RNA bursts), which lowers rho in every unit; fetal liver below postnatal liver cannot be
  separated from noisier fetal introns (unspliced clock-gene median 39). The damping route is uninformative here.

## A340 (scored 2026-11-22; BLIND for the fetal geometry; one run): the fetal clock-plane states are not a single slow displacement (H340.1 PASS: del Valle 14 genes and core, liver), but more elongated than a daily clock sampled in hospital hours; in liver the elongation lies along the glucocorticoid-target direction
- Improperness iota (counting-noise corrected), references noise-matched to each unit (power note), 2,000 bootstrap
  resamples (`A340_fetal_geometry/output/score/{iota,verdicts}.csv`). Code check first (genes permuted within cohort:
  iota inside the no-clock reference ranges).

  | unit (n) | genes | iota (90 %) | line q05 | daily 5-95 % | no clock 5-95 % | p vs line (Holm) | H340.1 | main axis to GC direction |
  |---|---|---|---|---|---|---|---|---|
  | del Valle (32) | 14 | 0.41 (0.26 to 0.71) | 0.97 | 0.13-0.37 | 0.56-0.83 | 0.001 (0.002) | PASS | 32 deg |
  | del Valle (32) | core | 0.55 (0.28 to 0.78) | 0.91 | 0.10-0.42 | 0.72-0.91 | 0.001 | PASS | 88 deg |
  | liver (37) | 14 | 0.68 (0.50 to 0.83) | 0.95 | 0.14-0.53 | 0.42-0.74 | 0.001 (0.002) | PASS | 8 deg |

  Descriptive: liver core 0.83 (0.71 to 0.90; not decidable), axis 2 deg from GC; 10x cortex 0.58 (0.39 to 0.81), inside
  its daily range (0.13-0.60), p vs line 0.008; `cortex_minus1` 0.65.
- **Reading.** H340.1 PASS in all three: the fetal clock-shaped covariation is **not a single slow factor** moving
  fetuses along one line (age error, glucocorticoid tone or termination protocol acting alone). The registered readings
  did not foresee the position found: every point estimate lies above the daily 5-95 % range (del Valle 14 genes just
  above, 0.41 vs 0.37), and far below the line. Post hoc wording: "two-dimensional like a clock sampled across the day,
  with an additional axial component", or a narrower collection window than the 08-18 h GTEx frame; the bootstrap
  intervals overlap the daily range in del Valle (14 genes and core) and touch it in liver. In liver the main axis lies
  8 deg (core 2 deg) from the glucocorticoid-target direction, consistent with a glucocorticoid-like axial component on
  top of the turning state; in the del Valle adrenal it does not (32 and 88 deg). Del Valle 14 genes and core lie below
  the no-clock range, liver inside it (the covariation itself is clock-specific in these units, A328 part a, A243 / A269).
- **For claim 3 (Q43 framing):** part 3 ("daily timescale") now has a lead-free argument against the single-slow-factor
  alternative; the ultradian alternative is not excluded (A339 uninformative, A328 lead upper bounds post hoc only).
- *Output addendum (2026-11-23, no statistic changed):* `score_a340.py` now writes the power stage's no-clock 5-95 %
  range (`noclock_q05`, `noclock_q95`) into `iota.csv` / `verdicts.csv`; the rerun reproduces every earlier column
  exactly (max difference 0). No-clock ranges: del Valle 14 genes 0.56-0.83, core 0.72-0.91, liver 0.42-0.74.
<!-- src: paper/exploratory/A340_fetal_geometry/output/score/{iota,verdicts}.csv -->

## A341 (computed 2026-11-22; descriptive by the power rule): neither glucocorticoid responsiveness nor clock maturity ranks the fetal organs that carry the clock-shaped state
- 7 organ readouts as scored (A244 cortex and medulla, A243 lung, liver, heart, kidney, A239 prenatal brain); measures
  residualised as each readout (`A341_gr_ranking/output/{measures,spearman}.csv`).

  | organ (n) | C | NR3C1 log CPM | GC-score variance (noise-corrected) | core-gene log CPM | ARNTL - NR1D1 |
  |---|---|---|---|---|---|
  | adrenal cortex (34) | 0.608 | 2.13 | 0.129 | 2.04 | 0.17 |
  | liver (35) | 0.583 | 3.86 | 0.113 | 2.57 | 0.92 |
  | lung (38) | 0.530 | 5.05 | 0.057 | 3.19 | 0.61 |
  | heart (72) | 0.032 | 4.53 | 0.174 | 3.08 | 0.37 |
  | kidney (53) | 0.017 | 4.17 | 0.052 | 3.18 | 0.64 |
  | adrenal medulla (30) | -0.048 | 3.75 | 0.175 | 2.29 | 0.24 |
  | brain (20, arrays) | -0.166 | - | 0.007 (3 GC genes, uncorrected) | - | - |

  Spearman with C (one-sided exact p): NR3C1 -0.20 (0.67, 6 organs), GC-score variance 0.25 (0.30, 7), core-gene level
  -0.26 (0.72, 6), ARNTL - NR1D1 -0.03 (0.54, 6).
- **Reading (descriptive):** no measure orders the organs. Heart and adrenal medulla have the largest between-fetus
  glucocorticoid-target variation yet no clock-shaped covariation; the adrenal cortex has the lowest NR3C1 level and the
  strongest covariation. "Organs are selected by glucocorticoid responsiveness" is not supported at the organ level (with
  the power caveat: 7 organs, detectable rho 0.90), nor is a simple clock-maturity ranking. Not done: the mouse
  cell-type replication (descriptive in the registration).

## Scout (2026-11-22; metadata only): no open fetal cohort supports the host-referenced organ-offset test
- Sonnet subagent (web metadata) plus a check of the GSE66302 sample names on disk. Roost 2015 (DeepSAGE, 18 fetuses,
  21 organs): only 5 fetuses have >= 2 of adrenal, liver, lung (2 with all three: CC2, CF2; liver + lung: CI1, CR1,
  DJ1). Cao 2020 (GSE156793, open processed counts, on disk): 28 fetuses, but each organ was taken from only 3-11 of
  them; with >= 1,000 cells, adrenal + liver 4 fetuses, adrenal + lung 6, liver + lung 6, all three 4 (organ
  pseudobulk depths 1-190 M UMI, adequate): the limit is the number of fetuses sharing organs, not depth or access
  (A241: 7 fetuses with an adrenal, adult pass rate 0.14). Suntharalingham
  2025 (45,X and controls, HDBR): no adrenal, lung only pooled. Gerrard 2016: pooled embryos. Han 2020 HCL: same-donor
  overlap not shown. ENCODE / Roadmap fetal donors: not resolved (no shared adrenal-lung-liver donor found). No source
  reports collection time; termination type, where stated, is elective surgical.
- Consequence: the coupling and phase offsets between fetal organs of one individual, clocksync's own readout, cannot be
  measured on open data; it stays a prediction for a fetal cohort in which >= 20 fetuses each give two clock-positive
  organs (Cao's dbGaP raw reads would not add fetuses).

## A343 (computed 2026-11-22): an acute dexamethasone input raises pre-mRNA far ahead of mature mRNA in living human fat (H343.1 PASS) and liver (descriptive, same size)
<!-- src: paper/exploratory/A343_dex_acute/output/score/{verdicts,genes,h343_2,h343_3,two_tissue}.csv; code paper/exploratory/A343_dex_acute/scripts/score_a343.py -->
GSE320593, pre-operative (N) vs intra-operative (S) biopsies; primary units = patients whose two runs share library type.
GR set FKBP5, TSC22D3, ZBTB16, KLF9 against matched controls of the same runs; 200 matched four-gene companions.
- **P343.0 PASS (fat):** DEX patients (38), S - N mRNA change of the GR set (log CPM, vs controls) +0.377, sign-flip
  p < 0.0001, companion q95 0.110. Liver (13, descriptive): +0.990, q95 0.129.
- **H343.1 PASS (fat):** intron-to-exon ratio change D = +1.934 (log), p < 0.0001, companion q95 0.226 (median 0.011).
  Liver (descriptive): +2.393, q95 0.305. Direct-intron sensitivity smaller but above the bar in both: fat +0.963
  (q95 0.184), liver +1.277 (0.224). The pre-mRNA excess (D) exceeds the mRNA rise (E) about five-fold (fat) and
  2.4-fold (liver): the biopsies catch transcription before mature mRNA has caught up.
- Per gene (DEX, primary; D subtraction / direct): FKBP5 2.01 / 1.95 (fat), 2.38 / 2.40 (liver) agree across intron
  definitions; TSC22D3 (1.38 / 0.21 fat, 4.24 / 0.80 liver) and ZBTB16 (4.01 / 1.26, 2.69 / 1.72) are larger by
  subtraction, read as DEX-induced transcripts outside the annotated exons counted as intron by subtraction (GILZ
  isoforms of TSC22D3); the direct definition is the conservative figure. KLF9 small (0.31, 0.40).
- *H343.2 (descriptive):* DEX minus NO_DEX in D: fat +1.709 (permutation p 0.002; 38 vs 7), liver +1.973 (p < 0.001;
  13 vs 5). On all patients the NO_DEX liver D (+1.171, above its q95) carries the library-type artefact the design
  anticipated (6 of 11 NO_DEX surgical runs stranded); primary NO_DEX: fat +0.226, liver +0.420.
- *H343.3 (descriptive):* PER1 S - N change, DEX minus NO_DEX: fat +0.495 (90 % 0.240 to 0.770), liver +0.522 (0.071
  to 0.945). Core-loop displacement difference not resolved (fat primary 1.37, 90 % upper bound 2.26 state units;
  liver 0.76, 1.57): the S - N displacement is about one to two daily amplitudes in both arms (time of day and surgery
  differ between the two biopsies), so the design cannot say whether acute DEX moves the core loop.
- *Two-tissue (descriptive, 21 patients, code 53 excluded):* liver-fat complex correlation of the 14-gene state 0.754
  at N (offset liver minus fat +2.85 h) and 0.491 at S (+3.60 h); no companions, no claim.
- Run note: a first launch crashed before computing anything (`ARNTL` absent: GENCODE v45 names it BMAL1); symbols are
  now mapped through `clocksync.genes.canonical`; the run above is the single scoring run.
- **Reading (as registered):** an acute transcriptional input of known kind appears in pre-mRNA well before mature mRNA
  in living human tissue at bulk depth: the premise of the lead readouts (A328, A337, A342) holds with a known input.
  Acute glucocorticoid raises PER1 (+0.5 log) within the surgical window; whether the core loop moves is not decidable
  here.

## A342 (computed 2026-11-22; unspliced count corrected to direct introns, see erratum): in timed living skin the untimed pre-mRNA lead equals the timed lead (H342.1, H342.2 PASS); the E-P offset equals the timed offset but the pairing is not clock-specific (H342.3 FAIL); no sun-exposure amplitude difference (descriptive)
<!-- src: paper/exploratory/A342_skin_timed/output/score_direct/{verdicts,leads,pairing,genes,sensitivity}.csv (scored run); output/score/ = the bugged first run (subtraction introns); code paper/exploratory/A342_skin_timed/scripts/score_a342.py --direct -->
GSE326530, 20 volunteers x photoprotected (P) / photoexposed (E) x 4 times; one waveform (GTEx not-sun-exposed skin)
for all statistics (power note); volunteer-cluster bootstraps. Unspliced = `counts_intron_direct` (A342 erratum: the
registered subtraction count is 0 for CLOCK and PER2 in this unstranded delivery; user ruled it a bug).
- **P342.0 PASS:** spliced state rhythmic at both sites (p < 0.0001, 14 genes and core).
- **H342.1 PASS (both sites, both gene sets):** untimed lead, matched estimator: P +1.19 h (90 % 0.81 to 1.64), E +1.41
  (1.08 to 1.79); core P +1.19 (0.74 to 1.76), E +1.92 (1.34 to 2.52); gates p < 0.0001. Plain P +1.19, E +1.45;
  per-gene P +0.80 (0.45 to 1.11), E +1.03 (0.65 to 1.45).
- **H342.2 PASS (all four):** timed lead P +1.18 h (0.59 to 1.86), E +1.48 (1.08 to 1.93) (core +1.21, +2.25);
  untimed minus timed P +0.01 (-0.43 to 0.56), E -0.07 (-0.47 to 0.37), core -0.02 (-0.57 to 0.60), -0.33 (-0.77 to
  0.16), all within +-1.5 h. The timed lead (1.2-1.5 h) matches A337's median 1.27 h.
- Sensitivities: QC (58 visits) and site waveforms keep every lead within +-1.25 h of the timed one. On subtraction
  introns (the bugged first run) the timed lead is not resolved from 0 (P +0.23, E +0.57) and the untimed lead is
  inflated (P +2.89, E +1.81): four of 14 genes carry no unspliced signal there.
- **H342.3 FAIL (on clock-specificity only):** untimed E-P offset +1.33 h (0.35 to 2.46) vs timed +1.57 (0.42 to
  2.95), inside the timed interval (core +1.02 vs +1.19): the offsets agree. The pairing gate passes (coherence 0.467,
  p < 0.0001) but sits below both companions' q95 (matched 0.549, responsive 0.644): pseudo-clock genes pair as
  strongly across the two sites of one person at one time. The untimed pairing has no volunteer covariate (as
  registered), so person-, visit- and processing-level covariance shared by the two biopsies of a visit enters every
  gene set. Offsets agree in every variant (QC 1.10 vs 1.13; volunteers 1-11 0.87 vs 0.73; 12-20 1.99 vs 2.79).
- *H342.4 (descriptive):* untimed log relative amplitude E/P -0.253 (-0.429 to -0.081) vs timed +0.087 (-0.319 to
  0.514); difference -0.34 (-0.73 to 0.00). The untimed ratio reports a quieter exposed site that the timed fit does
  not see (likely more non-daily in-plane variance at the protected site).
- *H342.5 (descriptive):* timed log amplitude ratio E/P +0.087 (-0.319 to 0.514; core -0.017), inside the matched
  companions' 5-95 % range (-1.45 to 1.56): no detectable sun-exposure change in clock amplitude. Timed E-P offset
  +1.57 h (0.42 to 2.95; sign as clocksync's offset, case minus reference), not within +-1 h: the exposed forearm and
  the protected buttock differ in phase by about 1.5 h (A251 untimed epidermis forearm-buttock +0.73 h); body site and
  exposure cannot be separated here.
- Run notes: the bugged first run crashed once in the sensitivity block and its rerun reproduced its main tables byte
  for byte (md5); the scored direct-intron run reproduces the pairing table of the first run exactly.
- **Readings:** (1) pre-mRNA leads mRNA by about 1.2-1.5 h in timed living human skin, and the untimed lead statistic
  recovers it at both sites: the first human ground-truth validation of the lead readout. (2) The host-referenced
  offset equals the timed offset in living people (1.3 vs 1.6 h), but in this design the pairing cannot be shown
  clock-specific against companions. (3) Chronic sun exposure leaves skin-clock amplitude unchanged within the
  precision of 20 people.

### 2026-11-22 post hoc per-gene check: why subtraction and direct introns differ (descriptive; no hypothesis)
<!-- src: paper/exploratory/A342_skin_timed/output/intron_definitions/genes.csv; code paper/exploratory/A342_skin_timed/scripts/intron_definitions.py -->
- Skin (GSE326530, all 160 runs unstranded), totals over runs, subtraction vs direct: CLOCK 0 vs 16,143; PER2 0 vs
  10,529; RORA 2,989 vs 176,924; PER3 6,979 vs 29,498. The other ten clock genes agree (direct / subtraction 0.51-0.83:
  direct is the stricter region). All four affected gene bodies overlap another gene (CLOCK: TMEM165; RORA: RORA-AS1;
  PER2, PER3: lncRNAs or UTS2). Mechanism (provenance caveat of the delivery): featureCounts without `-O` drops reads in
  the overlap at gene-body level but keeps them at exon level, so gene body minus exon is clipped to 0. In unstranded
  runs genes on opposite strands overlap too, hence the loss. The direct count uses intron-only regions with every
  other gene's exons removed, and is unaffected.
- Timed lags of intron over exon per gene agree between definitions where both counts are sound (protected site:
  ARNTL 3.46 / 3.46 h, PER1 1.07 / 1.12, CRY1 1.49 / 1.56, CRY2 0.57 / 0.64, TEF -1.73 / -1.46, DBP 0.57 / 0.34) and
  are undefined by subtraction for CLOCK and PER2 (constant zero). The A342 subtraction-intron lead was therefore read
  on a state in which four of 14 genes carry no unspliced signal.
- Dexamethasone set (125 unstranded, 47 stranded runs): CLOCK and PER2 reduced (direct / subtraction 2.5 and 2.3),
  not lost; the GR genes of A343 are not affected in this way (FKBP5 agrees; TSC22D3 and ZBTB16 are larger by
  subtraction, the opposite direction, not explained by this mechanism).
- Every other delivery is stranded (Cardoso-Moreira, del Valle, Greenwell, baboon, oral mucosa, Lundell, Perrin, Tsoi):
  no clock gene shows the loss (direct / subtraction <= 2 for all 14), so A328 part b, A338, A339 and A337 are not
  affected.
- Consequence: direct introns as primary for unstranded libraries; for stranded libraries the two agree for clock
  genes. A342 is scored on direct introns (erratum, user decision).

## A337 Perrin (scored 2026-11-23; BLIND for the unspliced counts): in timed human muscle the untimed lead equals the timed lead (H337.1 PASS, H337.2 PASS on the 14 genes; core loop H337.2 FAIL); responsive condition not met (H337.3 FAIL)
<!-- src: paper/exploratory/A337_timed_lead/output/perrin/{compartments,verdicts,summary,companions,genes}.csv; code paper/exploratory/A337_timed_lead/scripts/score_a337.py --sets perrin -->
Perrin 2018 vastus lateralis, adapter-trimmed re-map, 57 libraries (10 subjects x 6 times, one subject 3), two lane runs
summed, subject as covariate. Stranded: subtraction and direct introns agree (values below: subtraction / direct).
- Rhythmic (p < 0.0001, both gene sets).
- **H337.1 PASS:** untimed lead +1.41 h (90 % 1.17 to 1.66) / +1.42 (1.18 to 1.66); core +1.66 (1.15 to 2.19) /
  +1.63; gate p < 0.0001 (kappa 0.93; core 0.78).
- **H337.2 PASS (14 genes):** timed lead +1.75 h (1.47 to 2.05); untimed minus timed -0.34 (-0.59 to -0.10), TOST
  p < 0.001 within +-1.5 h. **Core loop FAIL:** timed +3.00 (2.31 to 3.80), difference -1.34 (-2.19 to -0.70),
  TOST p 0.33: on the core loop the untimed lead under-reads the timed one by about 1.3 h.
- **H337.3 FAIL:** delta_U below the responsive companions' q95 (3.67 h; core 7.77): rhythmic responsive genes in a
  running clock lead as much, the H328.3 sensitivity problem A337 anticipated.
- Per gene (cosine fits, lag of pre-mRNA over mRNA): ARNTL 1.87 h, CLOCK 3.03, NPAS2 2.52, PER1 1.41, PER2 0.67, PER3
  2.13, CRY2 2.43, NR1D1 0.60, NR1D2 0.74; CRY1, DBP, TEF and RORA small amplitude or long genes (RORA -7.65,
  unspliced mean 22,492: intronic reads dominated by its 700 kb of introns).
- With baboon and Greenwell (A337 entry above: H337.1 3/7, H337.2 3/3 where H337.1 passed) and timed skin (A342: lead
  1.2-1.5 h recovered at both sites), the untimed lead recovers the timed pre-mRNA lead in both living human tissues
  tested (n 57 and 80), and in fewer animal compartments at n 12. Timed leads: muscle 1.75 h, skin 1.18-1.48 h,
  A337 animal median 1.27 h.

## A344 (scored 2026-11-23; NOT BLIND for the full cycle): on hospital-hours samples alone the untimed pre-mRNA lead is recovered (H344.2 PASS); the offset keeps its value but its interval is too wide to decide (H344.1 undecidable); two daytime points are not enough (H344.3); per-visit phases degrade as the arc shrinks
<!-- src: paper/exploratory/A342_skin_timed/output/score_a344/{subsets,verdicts}.csv; code paper/exploratory/A342_skin_timed/scripts/score_a344.py -->
GSE326530, direct introns; subsets of clock times, truth = full-cycle timed values (A342); paired volunteer bootstraps.
The full set reproduces A342 (untimed offset 1.33 vs timed 1.57 h; leads P 1.12, E 1.41).
- **H344.2 PASS (D3 = 06, 12, 18 h; 60 visits):** untimed lead minus full-cycle timed lead P +0.10 h (90 % -0.44 to
  0.63), E -0.09 (-0.57 to 0.53); core -0.13 (-0.75 to 0.54), -0.48 (-0.93 to 0.07).
- **H344.1 undecidable (D3):** untimed offset 1.37 h vs truth 1.57 (difference -0.20, 90 % -1.68 to 0.88; core
  +0.06, -1.28 to 1.01): the point estimate is kept but the interval exceeds +-1 h, largely the uncertainty of the
  timed truth itself (full set: -1.00 to 0.51).
- **H344.3 (D2 = 06, 12 h; 40 visits):** lead E passes (-0.38, -1.06 to 0.33); lead P fails (+0.85, 0.14 to 1.64);
  core fails at both sites (P +1.05, -0.02 to 2.28; E -1.16, -2.07 to -0.24); offset undecidable (+0.26, -2.08 to
  1.98; point 1.83 vs 1.57).
- *Amplitude (descriptive):* the untimed log ratio stays at -0.20 to -0.26 on every subset against the timed +0.087:
  the H342.4 offset neither worsens nor improves on partial arcs.
- *Per-visit phase (descriptive):* concentration R and median absolute error of the visit's state angle against clock
  time: full 0.555 / 2.03 h; D3 0.453 / 2.59 h; D2 0.358 / 3.27 h; N2 0.259 / 3.79 h (all p <= 0.015). Sample times
  are read coarsely and worse as the arc shrinks, while the lead holds on D3: relative readouts survive hospital hours
  better than sample times.
- *N2 night pair (descriptive):* offset -1.27 h from truth, lead P -1.03, E -0.46.
- **Reading:** with three daytime points per person (hospital hours) the untimed lead equals the full-cycle truth at
  both sites; the offset keeps its value (1.37 vs 1.57 h) without a decidable interval; with two points the lead
  drifts by up to about 1 h. Assigning times to samples degrades steadily on partial arcs, which is why the method
  reports relationships rather than sample times.

## A345 (scored 2026-11-23; NOT BLIND for the spliced psoriasis delay, BLIND for introns): the psoriatic lesional delay is carried by pre-mRNA (H345.1 PASS) and is clock-specific at the transcription level (H345.3 PASS, narrowly); equality with the mRNA delay not shown (H345.2 fails as registered: interval too wide), the point estimate puts pre-mRNA at least as late as mRNA
<!-- src: paper/exploratory/A345_tsoi_premrna/output/score/results.csv, output/power/summary.csv; code paper/exploratory/A345_tsoi_premrna/scripts/a345.py -->
Tsoi 2019 (GSE121212), `tsoi2019_skin_trimmed` (delivered 2026-11-23 within the registered window); psoriasis 27
lesional / non-lesional pairs, AD 21 (descriptive). Spliced `counts_exon`, unspliced `counts_intron_direct`; offsets
lesional minus non-lesional (negative = lesional delayed), 90 % patient-bootstrap intervals.
- **P345.0 PASS:** spliced offset -1.38 h (90 % -2.25 to -0.67; delay p 0.001, gate p 0.005); core -1.20 h (p 0.006).
  Reproduces the authors'-count value (-1.45 h).
- **H345.1 PASS:** unspliced offset -2.44 h (-4.35 to -1.05; delay p 0.001, gate p 0.005); core -1.79 h (-2.88 to
  -0.71; p 0.006). Subtraction introns (sensitivity, stranded): -2.26 / -1.50 h, p 0.001 / 0.005.
- **H345.2 fails as registered:** D = unspliced - spliced -1.05 h (-2.64 to +0.11); core -0.59 h (-1.54 to +0.36).
  Neither interval lies inside +-1 h. The realised half-widths (1.37 / 0.95 h) are about three and two times the
  power stand-in's (0.31 / 0.46 h): real pre-mRNA states are noisier than thinned mRNA (the power-note caveat). The
  interval does not exclude equality and excludes the post-transcriptional reading (D > 0) only at its upper end
  (+0.11 / +0.36): pre-mRNA is not less delayed than mRNA.
- **H345.3 PASS (narrow on matched):** coherence-weighted unspliced delay |c| sin(arg c) -0.345 vs matched q05 -0.316
  and responsive q05 -0.229; core -0.299 vs -0.293 and -0.092.
- *AD (descriptive):* spliced -1.15 h (p 0.001), core -1.17 h; unspliced -1.30 h (p 0.007), core -0.81 h (p 0.10,
  gate p 0.015); D -0.15 (-0.77 to +0.47; inside +-1), core +0.36 (-0.52 to +1.35); companions not passed (-0.257 vs
  matched q05 -0.222 but responsive q05 -0.411). On these counts AD lesions also show a spliced delay, unlike A115
  (authors' counts, uninformative) and the A254 reading (no delay in other chronic lesions); descriptive, one series.
- *Not computed:* the registered descriptive lead change within lesional and non-lesional samples (the A342 lead
  statistic) is not in `a345.py`; it can be added as a descriptive follow-up.
- **Reading:** the lesional delay is present in the clock's current transcription, not only in mature mRNA, and beats
  matched and stress-responsive gene sets at the pre-mRNA level, which A183 could not show on mRNA. The registered
  "transcriptional" reading also needs H345.2, which is not decided: the data favour a pre-mRNA delay at least as
  large as the mRNA delay, and rule out a delay confined to mRNA turnover or composition only at the interval's edge.
  Limits: one series; a transcriptional delay can still be compositional at cell-type level (bulk).

## A348 (scored 2026-11-23; NOT BLIND for the 14-gene values, BLIND for every core value; one run per part): tissues of the same living person share one clock state on the core loop, beyond matched and calibrated broad responsive genes, in skin sites, skin layers and the surgical fat-skin pairs (H348 PASS x4); MuTHER fat-skin falls to the line (matched p_cal 0.057)

- Part 1: draw 0 on 14 genes reproduces every stored A320 statistic, q95, p_cal and placebo rate exactly (max |difference|
  0.0). Part 2 (core loop: 11 genes; GSE133979 lacks RORA, 10 genes), 2,000 placebo draws, pairing p from 10,000
  permutations (0.0001 in every test, the floor). Runs under a minute each (`--jobs 4`); one code-bug crash in the repro
  block before any output was written, fixed and rerun.

| test | 14 genes -> core | matched q95 / p_cal | broad responsive q95 / p_cal | narrow responsive q95 / p_cal (floor) | H348 |
|---|---|---|---|---|---|
| A251 H251.1 skin sites | 0.628 -> 0.573 | 0.472 / 0.022 | 0.398 / 0.013 | 0.470 / 0.026 (0.025) | PASS |
| A251 H251.2a | 0.813 -> 0.696 | 0.458 / 0.008 | 0.638 / 0.013 | 0.712 / 0.127 (0.017) | PASS |
| A251 H251.2b | 0.594 -> 0.510 | 0.262 / 0.024 | 0.318 / 0.013 | 0.388 / 0.021 (0.020) | PASS |
| A255 H255.1 MuTHER fat-skin | 0.533 -> 0.311 | 0.328 / 0.057 | 0.278 / 0.019 | 0.348 / 0.197 (0.031) | FAIL (matched) |
| A255 H255.3 surgical pooled | 0.686 -> 0.719 | 0.406 / 0.004 | 0.450 / 0.003 | 0.507 / 0.013 (0.012) | PASS |

- **V348:** matched and broad families calibrated in every test (q95-rule placebo rate 0.037-0.061); the narrow family
  is not (0.080-0.105), as in A285, and is read by p_cal only. Its pi = 1 floor is below 0.05 in all five tests, so on
  the core loop it is decidable: it passes in H251.1, H251.2b and H255.3 and fails in H251.2a and MuTHER.
- **Reading (registered).** Skin sites and layers (A251) and the surgical fat-skin pairs (A255 H255.3, 83 persons) share
  one clock state without PER1, PER2 and NR1D1, beyond expression-matched and calibrated broad responsive genes; the
  surgical pairs also beyond glucocorticoid / immediate-early genes. MuTHER fat-skin (one punch biopsy, one array batch)
  keeps the pairing and the broad family but falls below matched genes: at the line on the core loop. Gap 2 of the
  story audit is closed with this wording.
- Code: `exploratory/A285_companion_calibration/scripts/placebo_a285.py --test a251|a255 --core --families
  matched,responsive,responsive_broad --out a348` (14-gene check: `--out a348_check14`).
<!-- src: paper/exploratory/A285_companion_calibration/output/a348/{a251,a255}/{repro,tests}.csv, output/a348_check14/ -->

## A349 (run 2026-11-23; POWER AND DESCRIPTIVE, no verdict; BLIND for every value below except the D349.4 reproduction): at the ventilated design an organ-order change is detectable only from about 6-7 h in one organ (not 3.9 h) or a dispersion SD of 2.25 h; the ventilated host spread is about half that of fast deaths; host-referenced amplitude falls in most organs, most in organs with a strong fast-death clock; removing the A109 non-temporal axis relaxes the A132b compression from 0.58 to 0.77

GTEx fast deaths (244) against the 79 ventilated post-mortem donors with a living brain, 12 body tissues; each part run
once (part 2: 166 s, peak 1.9 GB). Aggregates only.

- **Part 1, detection limit (200 draws; 66 organ pairs per draw).** Ventilated arm resampled with replacement against
  the full fast-death map; companion = 79 fast deaths drawn without replacement against the remaining 165 (the full map
  cannot contain the test donors).

| gene set | arm | unrotated rho median (5th pct) | one organ: median limit (range) | dispersion SD limit |
|---|---|---|---|---|
| 14 genes | ventilated | 0.577 (0.395) | 6.75 h (2.25-12.0; skin not sun-exposed never) | 2.25 h |
| 14 genes | fast, 79 | 0.778 (0.550) | 5.13 h (2.75-10.5) | 2.00 h |
| core | ventilated | 0.593 (0.421) | 5.63 h (2.5-11.25) | 2.25 h |
| core | fast, 79 | 0.813 (0.627) | 4.13 h (2.75-8.75) | 1.50 h |

  The A336 fast-death limit (3.9 h, SD 1.5 h) is optimistic for intensive care: at 79 donors and their weaker shared
  signal, one organ must move about 6-7 h (14 genes; 5.6 h on the core) to be seen; nerve and subcutaneous fat about
  2.25-2.5 h, blood only at 12 h.
- **D349.1, host spread V** (separate fits, median over targets): 14 genes fast 3.61, ventilated 1.79; core 4.43 against
  2.09 (sigma_f trace 3.61 / 1.66 and 4.39 / 1.99). Median log V ratio fast / ventilated in the joint contrast 0.53
  (pooled) and 0.52 (age / sex matched; core 0.50 / 0.48). The ventilated spread is about half that of fast deaths,
  not the tenth recorded in the A315 registration (cause of the difference not checked).
- **D349.2, host-referenced log z2 fast over ventilated at a common spread** (median over 12 tissues; donor bootstrap
  500): 14 genes 0.142 pooled, 0.111 matched; core 0.201, 0.205. Above 0 in 7 (14 genes) and 8 (core) of 12 tissues;
  90 % interval excluding 0 in 6 / 12 (14 genes) and 3 / 12 (core). Largest losses thyroid +1.06 (0.73-1.44), oesophagus
  mucosa +0.77, heart +0.53; gains in whole blood -1.95 (-2.97 to -0.99), skin not sun-exposed -0.61 (-1.10 to -0.10).
  Median d* 0.075 / 0.059 (14 genes), 0.067 / 0.059 (core); coupling fast 0.686 against ventilated 0.437 (14 genes
  pooled). Beside A315 / A335 Q* +0.404 (pooled), +0.389 / +0.459 (matched, 14 genes / core): the same direction.
- **D349.3, relative gain under the loss:** per-tissue ratio minus the median spans -2.10 (blood) to +0.92 (thyroid);
  the loss rises with the tissue's fast-death z2 (Spearman 0.818 pooled, 0.839 matched; core 0.797, 0.755): the organs
  with the strongest clocks lose most, blood and skin least.
- **D349.4, compression after removing the non-temporal axis.** Reproduction: 0.585 (stored 0.585, 34 pairs). Host
  component along the A109 axis (3.90 h, eigenvalue share 0.66) removed per group: ratio 0.767 (30 pairs), same pole
  0.909, order Spearman 0.219. The removal overshoots (share along the axis -0.10 fast, -0.34 post mortem), so 0.767 is
  an approximate, operator-dependent value (the registration did not define the projection; chosen in the script).
- **Reading (descriptive).** For the paper: intensive care keeps the organ order at a detection limit of about 6-7 h for
  one organ (SD 2.25 h for dispersion), not 3.9 h; the compression toward 0 / 12 h is partly the shared non-temporal
  axis (0.58 -> about 0.77 after its removal). The loudness loss is host-referenced too (log z2 > 0 in most tissues) and
  is not uniform: it scales with clock strength.
- Code: `exploratory/A349_icu_limits/scripts/{a349_common,a349_part1,a349_part2,a349_d4}.py`.
<!-- src: paper/exploratory/A349_icu_limits/output/{part1_limits,part1_sim,part2_spread,part2_summary,part2_tissues,d4_compression}.csv -->

## A347 (scored 2026-11-23; NOT BLIND on 14 genes, BLIND on the core loop; one scoring run, 167 s): psoriatic lesions run late on both readouts and more than atopic lesions (H347.3 PASS); atopic lesions show a small delay on 14 genes (-0.26 h, p 0.034) that is not clock-specific and not detected on the core loop (-0.21 h, p 0.097, MDE 0.37 h): by the registered reading the psoriatic exception stands

17 gated series (authors' counts): psoriasis 9 / 465 pairs, AD 8 / 298. Pooled offset lesional minus non-lesional
(negative = delay), 10,000 swaps; companions calibrated on 2,000 swap draws (A308 q).

| | 14 genes | core loop |
|---|---|---|
| P347 psoriasis offset (p) | -1.507 h (< 1e-4); q matched < 0.001 | -0.958 h (< 1e-4); q matched < 0.001 |
| H347.1 AD offset (p; MDE) | -0.256 h (0.034; 0.35 h): PASS | -0.205 h (0.097; 0.37 h): FAIL |
| H347.2 AD q matched / narrow / broad | 0.114 / 0.101 / 0.103: FAIL | 0.128 / 0.078 / 0.185: FAIL |
| H347.3 D = PSO - AD (90 %) | -1.251 h (-1.613 to -0.920): PASS | -0.753 h (-1.166 to -0.368): PASS |
| H347.4 equivalence +-0.75 h | FAIL | FAIL |

- V347: calibrated (pass rate <= 0.0646) in 11 of 12 family x disease x readout cells; psoriasis core narrow
  responsive 0.070 (read with care; its q < 0.001).
- Per series (core): psoriasis 9 of 9 negative (-0.21 to -2.68 h); AD 7 of 8 negative (-0.06 to -1.70 h; GSE58558
  +0.15 h); the 1/se^2 weights (variance-ratio precision) give the two series nearest 0 (GSE58558, GSE99802) the largest
  AD weights. The Tsoi AD series reads -0.94 h (A345 mapped counts: -1.17 h).
- The 14-gene AD value reproduces the archive A79 pooled value (-0.257 h, p 0.040; not blind).
- **D347 (descriptive).** GSE154200 atopic eczema (3 gated batches, 37 pairs): -0.20 h (14 genes; 90 % -2.04 to +1.68,
  MDE 2.19 h), core +0.43 h: uninformative. Hidradenitis (GSE151243, 20 pairs, gate p 0.005 on both readouts): -1.07 h
  (-1.80 to -0.14, p 0.047), core -0.88 h (-1.58 to -0.05, p 0.053), q 0.31-0.38 against every family: delayed, not
  clock-specific. Tsoi lead (A342 matched estimator) within condition, plain form: psoriasis lesional 0.69 h against
  non-lesional 1.34 h (difference -0.65, 90 % -1.80 to +0.40; core 1.41 against 1.76 h), AD 0.37 against 0.76 h (core
  0.56 against 0.65 h): no lead change resolved. The registered patient-covariate form is degenerate within one
  condition (one sample per patient; it returns the same lead in both conditions) and is not reported. Not fully blind:
  the smoke check on shuffled labels showed the approximate overall lead per disease.
- **Reading (registered; the core loop decides wording).** H347.1 fails on the core with an MDE (0.37 h) below 0.75 h:
  "the psoriatic exception stands; atopic lesions do not run detectably late". On 14 genes a small AD delay is seen,
  not beyond companions; psoriasis runs later than AD on both readouts (by 1.25 h on 14 genes, 0.75 h on the core).
  Proposed wording: "psoriatic lesions run about 1-1.5 h late (core loop 0.96 h), more than atopic lesions, which lean
  late by about 0.2 h without a clock-specific or core-loop delay; hidradenitis leans late (descriptive, not
  clock-specific)". The Tsoi AD series alone (-0.94 to -1.17 h) is the largest AD value of the eight.
- Code: `exploratory/A347_inflamed_skin/scripts/{a347,d347}.py`.
<!-- src: paper/exploratory/A347_inflamed_skin/output/{power/power.csv,score/{families,difference,hypotheses,per_series}.csv,d347/*.csv} -->

## A351 calibration (2026-11-23; stop rule): V351 FAIL, the waveform screen does not separate per-tumour rewiring from damping (H vs D AUC 0.801 in SCAN-B, 0.890 in METABRIC; bar 0.90); "not separable at gene level", no condition scored

- A212 plasmode reproduced with its own generators, seeds and matched strengths (B1: rho 0.832, W 1.806 h, H 2.152 h;
  B2: rho 0.834, H 2.501 h, W at its strongest 2.0 h); 200 replicates per mechanism. Per gene: OLS complex coefficient
  on a fixed projector score per half, a_g = log amplitude ratio, phi_g = phase shift; 50 unit bootstraps for the
  variances; Q_a, Q_phi Cochran heterogeneity; joint AUC = 5-fold cross-validated LDA on (log Q_a, log Q_phi). 25 s.

| base | readout | W vs D joint (Q_phi alone) | H vs D joint (Q_a alone, Q_phi alone) |
|---|---|---|---|
| B1 SCAN-B | 14 genes | 0.999 (0.999) | 0.801 (0.212, 0.566) |
| B1 SCAN-B | core | 0.975 (0.977) | 0.804 (0.203, 0.615) |
| B2 METABRIC | 14 genes | 1.000 (1.000) | 0.890 (0.112, 0.341) |
| B2 METABRIC | core | 0.978 (0.974) | 0.682 (0.314, 0.519) |

- **V351 FAIL** (registered: all four 14-gene AUCs >= 0.90). Shared rewiring (W) is separated from uniform damping
  almost perfectly through Q_phi; per-tumour rewiring (H) is not: random per-tumour phase shifts average out in the
  pooled per-gene coefficient and read as a uniform attenuation (median Q_phi 16.5 under H against 14.2 under D, B1).
- **Reading (registered stop rule).** "Not separable at gene level": a uniform-looking loss cannot distinguish damping
  or dephasing from rewiring that differs between tumours; no condition of A351 is scored. Post hoc note, not a
  result: a screen restricted to "shared rewiring against uniform scaling" would be calibrated (AUC >= 0.975 on both
  bases and readouts); it would need its own registration.
- Code: `exploratory/A351_waveform_screen/scripts/a351_calib.py`.
<!-- src: paper/exploratory/A351_waveform_screen/output/calibration/{auc,replicates}.csv -->

## A346 (scored 2026-11-23; NOT BLIND on 14 genes, BLIND on the core loop; one scoring run): tumour fibroblasts keep the patient's normal fibroblast clock state on the core loop, unchanged (0.686 -> 0.683) and beyond all three companion families; tumour endothelium does not (0.677 -> 0.510, below matched genes); on 14 genes both pass the calibrated rule (narrowly against the narrow responsive family)

A227 cohorts and tables; 10,000 within-cohort permutations; 200 matched, 200 narrow responsive, 200 broad responsive
sets per cohort (core: restricted); 2,000 placebo draws; Holm over the two compartments. Power note: H346.1 was not
decidable at the registered delta (narrow family), so H346.1 is descriptive; pairing, matched and broad are decidable.

| compartment (patients / cohorts) | readout | clock | pairing p (Holm) | matched q95 / p_cal Holm | narrow responsive q95 / p_cal Holm | broad q95 / p_cal Holm |
|---|---|---|---|---|---|---|
| fibroblast (93 / 16) | 14 genes | 0.686 | 0.0002 | 0.539 / 0.010 | 0.638 / 0.0500 (0.04998) | 0.542 / 0.010 |
| fibroblast | core | 0.683 | 0.0002 | 0.558 / 0.009 | 0.635 / 0.012 | 0.540 / 0.012 |
| endothelial (101 / 14) | 14 genes | 0.677 | 0.0002 | 0.512 / 0.010 | 0.659 / 0.0500 (0.04998) | 0.529 / 0.010 |
| endothelial | core | 0.510 | 0.0009 | 0.529 / 0.104 | 0.663 / 0.841 | 0.529 / 0.091 |
| D346 malignant-normal epithelium (190 / 27) | 14 genes | 0.433 | 0.017 | 0.451 / 0.178 (raw) | 0.503 / 0.631 | 0.458 / 0.134 |
| D346 epithelium | core | 0.361 | 0.47 | 0.461 / 0.866 | 0.495 / 0.987 | 0.445 / 0.859 |

- **V346:** matched and broad calibrated everywhere (placebo rate 0.047-0.063) except epithelium core (0.076 / 0.079);
  the narrow family is not (0.071-0.091) and is read by p_cal only.
- **Trio contribution (14 minus core, 90 %, 200 patient bootstraps):** fibroblast +0.003 (-0.054 to +0.072);
  endothelial +0.167 (+0.058 to +0.208); epithelium +0.072 (-0.032 to +0.112).
- **H346.3 (14 genes) PASS** for both compartments: pairing, matched and narrow responsive (Holm 0.04998, at the line)
  calibrated; the A149 note 5 / A271 stromal result now has a calibrated registered form on 14 genes.
- **H346.1 (core; descriptive by the power note):** fibroblasts meet every condition (narrow p_cal Holm 0.012);
  endothelium keeps the pairing (p 0.0009) but falls below matched genes (0.510 against q95 0.529): its 14-gene
  specificity rests on PER1, PER2 and NR1D1 (trio contribution 0.17). **H346.2 (broad, core):** fibroblast PASS (0.012),
  endothelial FAIL (0.091).
- **Reading (registered).** "Tumour fibroblasts keep the patient's normal fibroblast clock state, clock-specifically,
  without PER1, PER2 and NR1D1, in the cohorts available" (beyond matched and broad responsive genes as a registered
  test; beyond glucocorticoid / immediate-early genes as a descriptive, since that condition was not decidable at the
  power stage). The endothelial pairing is not clock-specific on the core loop. The stroma clause can be worded for
  fibroblasts only; no independent cohort (A226). Malignant-normal epithelium stays not clock-specific (D346).
- Code: `exploratory/A346_tumour_stroma_core/scripts/a346.py` (power 2,510 s, peak 21.3 GB; score from the cache).
<!-- src: paper/exploratory/A346_tumour_stroma_core/output/{power/{repro14,power},score/{tests,trio}}.csv -->

## A350 (scored 2026-11-23; BLIND on the core loop; one refit, 4,555 s): relative class gain stays conserved on the core loop in four of five families, but the mural family exceeds the bound (tau 0.097 -> 0.389): H350.1 FAIL

A17 healthy atlas fits refitted on the core loop (21 atlases of A293 Part 3; gate passes in 14), learned anchors
frozen, 500 bootstraps; A293 Part 3 code unchanged.

| family | atlases 14 / core | pooled g, 14 -> core | tau, 14 -> core | Q p core | SE ratio core / 14 |
|---|---|---|---|---|---|
| fibroblast-stromal | 13 / 11 | +0.228 -> +0.192 | 0.181 -> 0.212 | 0.072 | 1.30 |
| epithelial | 11 / 9 | -0.220 -> -0.146 | 0.305 -> 0.283 | 0.030 | 1.11 |
| myeloid | 12 / 10 | -0.399 -> -0.214 | 0.302 -> 0.245 | 0.038 | 1.01 |
| parenchyma | 4 / 4 | -0.126 -> -0.058 | 0.145 -> 0.093 | 0.346 | 1.20 |
| mural | 7 / 7 | -0.079 -> -0.023 | 0.097 -> **0.389** | 0.001 | 2.39 |

- **H350.1 FAIL** (registered: tau <= 0.35 in every tested family): mural cells exceed the bound on the core loop; the
  mural atlas values spread from -0.83 (heart) to +0.93 (skeletal muscle), and the pooled SE more than doubles.
- **Reading (registered).** Relative gain is conserved across organs on the core loop for fibroblasts, epithelium,
  myeloid cells and parenchyma; mural cells' conservation is not shown without the trio (or with the core's lost
  precision: SE ratio 2.39). Wording: "class gains are conserved across organs (14 genes; on the core loop in four of
  five classes, not mural cells)".
- Code: `exploratory/A350_class_gain_core/scripts/a350.py`.
<!-- src: paper/exploratory/A350_class_gain_core/output/{compare,gates_core,a17_healthy_fits_core}.csv, output/part3/{atlases,families}.csv -->

## A352 power note (2026-11-23; no offset computed): the cytokine dose-response of the skin delay is not powered; A352 closed at the power stage

Signature tertiles within series (IL-17 / TNF, Th2; no series lost to platform coverage). The registered primary
(psoriasis, IL-17 / TNF, core loop, 465 patients) has an MDE of 1.19 h for the T3 - T1 difference in delay, against
0.75 h required; atopic dermatitis 0.95-1.58 h; only the pooled Th2 core trend (0.73 h) would be decidable, and it is
not a registered primary. Tertiles split patients three ways, so each stratum's pooled offset carries about the noise of
a third of the cohort (swap-null SD 0.15-0.59 h). Reading: not a negative result; the delay's dependence on IL-17 / TNF
activity cannot be read with tertile strata on these series. The score stage's null was not verified (smoke p 0.005 on
swapped labels). Code: `exploratory/A352_cytokine_delay/scripts/a352.py`.
<!-- src: paper/exploratory/A352_cytokine_delay/output/power/{power,coverage}.csv -->

## A353 (scored 2026-11-23; NOT BLIND for the full-cycle A337 Perrin values and the A344 skin outcome, BLIND for every subset value; one scoring run after one code check): the muscle pre-mRNA lead is recovered from hospital-hours biopsies (H353.1 PASS on 08, 12, 16 h; H353.2 PASS on 08-20 h and on 08, 12 h); reading each biopsy as a time of day degrades as the arc shrinks
<!-- src: paper/exploratory/A337_timed_lead/output/a353/{subsets,verdicts}.csv; code paper/exploratory/A337_timed_lead/scripts/score_a353.py -->
Perrin 2018 muscle, direct introns, 14 genes; truth = full-cycle timed lead (all 57 libraries) +1.75 h (90 % subject
bootstrap 1.49 to 2.02); paired 90 % intervals from 2,000 subject-cluster bootstraps; gate p from 10,000 within-subject
permutations.
- **Full cycle (reference, reproduces A337 Perrin):** untimed +1.42 h (1.19 to 1.66; A337's within-subject interval 1.19
  to 1.65), minus truth -0.32 h (-0.54 to -0.13).
- **H353.1 PASS (08, 12, 16 h; 9 subjects, 27 samples):** untimed +1.74 h (1.21 to 2.38), gate kappa 0.94 (p 0.0001);
  minus truth -0.01 h (-0.50 to 0.63), inside +-1.5 h.
- **H353.2 PASS x2:** 08-20 h (36 samples) +1.40 h, minus truth -0.35 (-0.57 to -0.10); 08, 12 h (18 samples) +1.90 h,
  minus truth +0.16 (-0.55 to 0.89), gate p 0.002.
- *Per-sample phase (descriptive):* R / median absolute error against clock time: full 0.79 / 1.42 h; 08-20 h 0.88 /
  1.30 h; 08, 12, 16 h 0.74 / 2.51 h; 08, 12 h 0.45 / 4.00 h; night 20, 00, 04 h 0.42 / 3.18 h.
- *Night 20, 00, 04 h (descriptive):* +1.05 h, minus truth -0.70 (-1.14 to -0.31).
- *Core loop (descriptive; full-cycle H337.2 core FAIL):* full minus truth -1.34 h (-2.48 to -0.59); 08, 12, 16 h -0.37
  (-1.74 to 0.64); 08-20 h -1.04 (-1.85 to -0.61); 08, 12 h -0.77 (-2.25 to 0.47).
- Code check (`--check`: thinned spliced counts as unspliced, times shuffled within subject): untimed leads -0.36 to
  +0.08 h, timed truth uninformative, every verdict undecidable. One fix after the check, before the scoring run: the
  full-cycle row kept only subjects with all six times (54 libraries); it now keeps all 57, as registered. The
  registration's list of subsets the three-time subject drops from omitted D4, which contains 08 h; the rule (subjects
  with every time of the subset) decides and drops it.
- **Reading (registered).** The muscle lead is recovered from three daytime biopsies per person, a second living tissue
  after skin (A344). Unlike skin, two daytime points (08, 12 h) also pass here, with six full-day times and a
  well-entrained protocol (hourly isocaloric meals, controlled light): the most favourable design. Sample times degrade
  as the arc shrinks (median error 1.42 h full, 2.51 h on 08-16 h, 4.00 h on 08-12 h) while the lead holds: relative
  readouts survive partial sampling better than sample times, as in skin.

## A354 (scored 2026-11-23; NOT BLIND for every clock-gene |r| of A279, BLIND for every companion value; one scoring run after one check run that computed no companion): in the clock-intact aorta every cell-type pair is above both companion families (H354.2 PASS, 6 of 6); in the Bmal1 knockout the stromal pairs fall to the top edge of the expression-matched sets, one above them (H354.1 FAIL: EC-fibroblast 99th percentile)
<!-- src: paper/exploratory/A279_aorta_pairwise_truth/output/a354/{pairs,design,verdicts}.csv; code paper/exploratory/A279_aorta_pairwise_truth/scripts/score_a354.py -->
Auerbach 2026 aorta, genome-wide library x cell type pseudobulks (integer counts), A279 filters; aligned 14 libraries,
4 classes; misaligned 8, 3; knockout 10, 4. 200 matched and 200 responsive sets per arm (28 responsive genes expressed).
Check run (`--check`): the clock |r| reproduces A279 `pairs.csv` in every arm (assert), no companion computed.
- **H354.2 PASS (aligned):** clock |r| 0.83-0.98, above the matched q95 (0.58-0.78) and the responsive q95 (0.71-0.94)
  in all 6 pairs (100th percentile of both families in 5 of 6; EC-macrophage 99.5th matched).
- **H354.1 FAIL (knockout, stromal pairs):** clock |r| against matched q95 (percentile): EC-fibroblast 0.753 against
  0.618 (99.0th; above the responsive q95 0.664 too, 99.5th), EC-SMC 0.629 against 0.636 (94.5th), fibroblast-SMC
  0.661 against 0.666 (93.5th; responsive 33.5th, q95 0.926). Two of three at or below the matched q95, as registered
  all three were required.
- *D354.3 (descriptive):* knockout pairs with macrophages 0.11-0.35, at the matched median (10.5th-58th percentile).
  Misaligned arm (3 stromal pairs): 0.89-0.96, at the 100th percentile of both families. Matched medians are similar in
  every arm (0.28-0.48): the generic co-variation between cell types of a library does not change with the clock.
- **Reading (registered).** The companion contrast separates a running clock from none: with the clock intact every
  pair is far above both families; without BMAL1 the stromal pairs drop to the top edge of the matched range. H354.1
  fails on one pair (EC-fibroblast, 99th percentile of matched and 99.5th of responsive): by the registered reading,
  clock genes of endothelium and fibroblasts covary beyond generic co-variation even without BMAL1 (a non-clock program
  shared by clock genes, or residual clock-gene regulation), so a single pair above the matched q95 is not proof of a
  clock. Caveat for the human cell-type results: they rest on pooled contrasts over many atlases and calibrated
  p-values (A283), not on one pair above q95; the knockout bounds what a non-clock pair can reach at the top of the
  companion range. Limits: one atlas, 10-14 libraries of pooled aortas.

## A355 power note (2026-11-23; gates and swap nulls only, no gain ratio computed): the pooled test and the skin pair are decidable, the fat and oesophagus pairs are not; the pooled weight is about 86 % skin

- Gates (`Dataset.gate`, all compartments, 1000 permutations): GSE119717 p 0.001 (clock share 0.274, null 0.083; 12
  genes), GSE39491 p 0.001 (0.532, null 0.105; 12 genes), Wu 2020 p 0.001 (0.665, null 0.086; 14 genes). All enter.
- Swap-null SD of Delta (1000 swaps) and MDE (2.49 SD) against the GTEx log ratio: GSE119717 omental / subcutaneous 0.250,
  MDE 0.622 vs 0.465: **not decidable**; GSE39491 cardia / squamous 0.196, MDE 0.489 vs 0.293: **not decidable**; Wu
  exposed / protected 0.062, MDE 0.155 vs 0.211: **decidable**; pooled T 0.057, MDE 0.143 vs 0.232: **decidable**.
- Weights 1 / SD^2: Wu 260, GSE39491 26, GSE119717 16: the pooled T is about 86 % the skin pair, whose contrast (sun
  exposure and site) is the weakest same-composition case; the reading states this.
- Choices not fixed by the registration, made before any observed Delta: (i) a clock gene absent from a cohort is
  dropped from its waveform (A255 rule): GSE119717 lacks NR1D1 and NPAS2 (the submitters' detection filter removed the
  probes from the deposited matrix), GSE39491 lacks PER1 and CIART on U133A 2.0; (ii) GSE119717 subjects 50 and 51 carry two
  arrays labelled subcutaneous each (flagged by the scout): both dropped as ambiguous (subcutaneous 26 persons, omental
  29); (iii) GSE39491 has one host compartment and `Dataset.fit` needs two anchors: its ratio is the cross-moment ratio
  log |sum z_cardia conj z_barrett| - log |sum z_squamous conj z_barrett| against the same single host (the gain ratio
  under the one-factor model, free of each target's own noise).
- Code: `exploratory/A355_living_gain_ratios/scripts/a355.py --stage power` (295 s).
<!-- src: paper/exploratory/A355_living_gain_ratios/output/{gates,power}.csv, swap_nulls.npz -->

## A355 (scored 2026-11-23; NOT BLIND for the GTEx gains, A251 and the A254 GSE39491 variance ratios; BLIND for GSE119717 and every living gain; one scoring run, 279 s): the living gain ratios do not follow the GTEx organ ratios as a set (H355.1 FAIL, pooled T 0.003, p 0.50; H355.2 FAIL for skin); fat goes the GTEx way at about two thirds of the size (not decidable), the cardia / squamous pair goes against it: the organ factor stays open

| pair (a / b) | living Delta (fold) | GTEx (fold) | one-sided p | 90 % person bootstrap | status |
|---|---|---|---|---|---|
| omental / subcutaneous fat, GSE119717 (12 genes) | +0.323 (1.38) | +0.465 (1.59) | 0.193 | -0.147 to 0.657 (GTEx inside) | not decidable |
| cardia / squamous oesophagus, GSE39491 (12 genes) | +0.395 (1.49) | -0.293 (0.75) | 0.980 | 0.087 to 0.786 (GTEx outside) | not decidable |
| exposed / protected skin, Wu 2020 (14 genes) | -0.023 (0.98) | -0.211 (0.81) | 0.354 (Holm 0.354) | -0.583 to 0.095 (GTEx inside) | **H355.2 FAIL** |
| pooled T (weights fat 0.053, oesophagus 0.087, skin 0.860) | 0.003 | 0.232 | 0.498 | -0.112 to 0.463 | **H355.1 FAIL** |

- **D355.1, core loop (10-11 genes):** fat +0.367 (-0.366 to 0.847), oesophagus +0.216 (-0.385 to 0.612), skin +0.186
  (-0.351 to 0.327).
- **D355.2, total clock-plane log variance ratio a / b against 200 matched sets:** fat 0.672 (83.5th percentile; matched
  median 0.185), oesophagus 0.027 (68th; -0.212), skin -0.381 (4.5th; -0.058): not clock-specific in fat or oesophagus;
  in skin the exposed sites carry less clock variance than matched genes would.
- **The power rule was optimistic.** For skin the person-bootstrap interval (width 0.68) is about three times wider than
  the swap-null SD (0.062) implied: swapping labels within persons keeps the persons and misses between-person sampling
  variance. The "decidable" skin test therefore carried less information than planned (its interval holds both 1 and
  the GTEx value). Future gain-ratio registrations: power from a person bootstrap or a plasmode, not from within-person
  swaps.
- **Reading (registered: H355.1 fails with T decidable).** The living ratios do not follow the GTEx organ gain ratios as
  a set, so tissue-specific post-mortem attenuation stays a live explanation, as do cohort and anatomical differences;
  "organ" stays open (Figure 4g). Descriptively, the cleanest same-composition pair, omental over subcutaneous fat, goes
  the GTEx way at about two thirds of its size (1.38 vs 1.59), not decidable at 26-29 men, and its total-variance
  difference is not clock-specific; the cardia / squamous pair goes against GTEx, where the anatomy does not match
  (GTEx Stomach is full-thickness stomach body with its muscle layers, GSE39491 mucosal biopsies of the cardia; PER1 and
  CIART absent on its array); the skin contrast is near 1 with a wide interval.
- Code: `exploratory/A355_living_gain_ratios/scripts/a355.py --stage score`.
<!-- src: paper/exploratory/A355_living_gain_ratios/output/{score_pairs,score_pooled,score_descriptive}.csv, boot.npz -->

## A356 power note (2026-11-23; permutation nulls only, no split gain computed): both pooled within-tissue tests are decidable (MDE 0.042 ischaemic time, 0.048 RIN, against the 0.10 bound); across 25 tissues the Spearman test has power 0.80 at rho -0.5

- 25 peripheral tissues (A263 fit, whole blood excluded); median 128 donors per tissue with an ischaemic time, 111 with
  a RIN; median of the halves' medians: ischaemic time 19.0 h (long) vs 12.2 h (short), RIN 7.0 vs 6.0.
- Per-tissue permutation SD of Delta: median 0.097 (ischaemic time), 0.105 (RIN). Pooled (inverse variance): SD 0.017,
  MDE 0.042 (ischaemic time); SD 0.019, MDE 0.048 (RIN): **decidable** (bound 0.10).
- Labels are permuted among different donors of a tissue, so the null carries between-donor sampling variance (unlike
  the within-person swaps of A355).
- Code: `exploratory/A356_gain_attenuation/scripts/a356.py --stage power` (16 s).
<!-- src: paper/exploratory/A356_gain_attenuation/output/{power,power_tissues}.csv (per-tissue aggregates) -->

## A356 (scored 2026-11-23, check G2; NOT BLIND for the GTEx organ gains, BLIND for every split gain and gain-covariate relation; one scoring run, 28 s): post-mortem attenuation does not make the GTEx organ gains: within tissues, gain does not fall with ischaemic time or rise with RIN (H356.1 FAIL, H356.2 FAIL, both decidable); across tissues, ischaemia and RIN could produce at most about 4 % of the organ gain spread

| test | estimate | 90 % (permutation-calibrated) | one-sided p | verdict |
|---|---|---|---|---|
| H356.1 pooled Delta, long minus short ischaemia (6.7 h apart, 25 tissues) | +0.010 (1.01-fold) | -0.018 to +0.038 | 0.756 | FAIL |
| H356.2 pooled Delta, high minus low RIN (0.9 units apart) | -0.006 (0.99-fold) | -0.038 to +0.026 | 0.631 | FAIL |
| H356.3a organ log gain vs tissue median ischaemic time | Spearman -0.366 | - | 0.038, Holm 0.076 | FAIL |
| H356.3b organ log gain vs tissue median RIN | Spearman -0.235 (against the prediction) | - | 0.875, Holm 0.875 | FAIL |

- **Code check:** per-tissue gain of all donors against the A263 `log_amp_rel`: Spearman 0.893.
- **D356.4:** organ log gains spread 0.732 (10th-90th percentile, 25 tissues); tissue median ischaemic time spreads
  0.94 h and median RIN 1.05 units; at the within-tissue slopes (0.0015 per hour, -0.0065 per RIN unit) attenuation
  produces a spread of 0.001 (ischaemia) and 0.007 (RIN); at the 90 % bounds of the slopes at most about 0.003 and 0.030,
  about 4 % of the observed spread. Per-tissue Delta: ischaemia median +0.015 (range -0.305 to +0.145), RIN +0.023
  (-0.497 to +0.347). Core loop pooled: ischaemia +0.011 (SD 0.020, p 0.72), RIN +0.000 (0.021, p 0.50).
- **Reading (registered: all FAIL with the pooled tests decidable).** No sign that post-mortem ischaemia or RNA
  degradation lowers clock gain within tissues or ranks organs; the GTEx organ gains stand with this check. The
  across-tissue rank correlation with median ischaemia (-0.37, at the line before Holm) cannot be attenuation: the
  within-tissue effect is null and organs differ by under an hour in median ischaemia; it more likely reflects
  procurement order confounded with organ. Not tested: agonal state and autolysis (fast deaths only, where agonal time is
  short), and handling routes other than ischaemia and RIN. With A355 (living pairs not reproduced as a set; fat in the
  GTEx direction, not decidable), the organ factor reads: organ gains are not composition (A323) and not post-mortem
  attenuation by ischaemia or RNA integrity (A356); their reproduction in living tissue is not shown (A355, A330).
- Code: `exploratory/A356_gain_attenuation/scripts/a356.py --stage score`.
<!-- src: paper/exploratory/A356_gain_attenuation/output/{score_pooled,score_across,score_descriptive,score_tissues}.csv (per-tissue aggregates) -->

## A357 power note (2026-11-23; matched sets only, no clock-gene residual computed): H357.2 is decidable (power 1.00 at a median R^2 excess of 0.03); H357.1 is not (power 0.08 on the registered SD, 0.66 on a robust SD): not scored

- One power run on the 100 matched sets of A294b (7,442 s, 12 workers, peak 2.7 GB; no errors). A first launch failed on a
  path before loading any data. A check run on 2 sets (timing only) preceded it.
- **H357.2 (decidable).** Median over tissues of the R^2 of the clock-plane departure on the person signal factor:
  matched median 0.008, q95 0.0153; SD over sets 0.0032, tissue-bootstrap SD 0.0027; power 1.00 against an excess of 0.03.
- **H357.1 (not decidable, not scored).** SD of Delta over sets 0.394 and of rho_sig 0.375 (98 sets with a value; 2 sets
  lack a sig-position value): power 0.08 against Delta = 0.10. The SD is inflated by the estimator: the host-error
  correction leaves a near-zero corrected variance for some matched genes, giving 14 of 1,321 per-gene medians outside
  [-1, 1] (9 sets; e.g. -3.30, +2.98) and 79 undefined. On a robust SD (IQR / 1.349: Delta 0.049, rho_sig 0.035) the
  power is still 0.66, below the 0.8 gate. With three signal genes, the per-set statistic is too noisy to resolve a 0.10
  difference in GTEx. By the registered rule H357.1 is not scored: the score stage computes no clock pair
  correlation output (`score_a357.py`, `h1_ok`), and its descriptive companions (per-gene rho of the glucocorticoid and
  immediate-early genes) are dropped with it.
- The H357.2 reading without H357.1 follows the registered line "H357.2 PASS without H357.1: reported as such; not read
  as the signal-layer mechanism" only if H357.1 had failed; with H357.1 unscored, an H357.2 PASS is read as "a
  person-level state of the signal-gene residuals explains part of the departure", without the sharing step.
<!-- src: paper/exploratory/A357_signal_departure/output/{power,power_sets,power_genes}.csv -->

## A357 (scored 2026-11-23; BLIND for every per-gene residual and the person signal factor; one scoring run, 76 s; H357.1 not scored, not powered): a person-level state of the signal-inducible clock genes, read from a person's other tissues, explains a small, clock-specific part of each tissue's departure (H357.2 PASS: median R^2 0.027 against 0.008 for matched genes); that state tracks the person's glucocorticoid-target state (descriptive r 0.77)

GTEx fast deaths, 25 peripheral tissues (Whole Blood and Testis excluded), 14-gene waveform, every tissue an anchor.
- **H357.2 PASS (decidable).**
  - Median over tissues of the R^2 of the clock-plane departure on the person signal factor s(-k): 0.027, above all
    100 matched sets (median 0.008, q95 0.0153).
  - Median excess over matched sets 0.019 (tissue bootstrap 90 % 0.0085 to 0.029).
  - Above the tissue's matched q95 in 6 of 25 tissues: skeletal muscle 0.195, liver 0.088, sun-exposed skin 0.061,
    oesophago-gastric junction 0.063, non-exposed skin 0.052, pituitary 0.044.
    Adrenal 0.096 and transverse colon 0.065 are below their (wider) matched q95.
- **H357.1 not scored** (power note: not resolvable in GTEx).
- **Descriptive (no companion).**
  - The person signal factor correlates with the mean glucocorticoid-target score of the person's other tissues (median
    over tissues r 0.77) and much less with their immediate-early score (0.19).
  - The same factor explains the core-loop (11-gene) departure at least as well (median R^2 0.034).
  - A factor built from the core genes' residuals explains little (0.011 on the 14-gene departure, 0.010 on the
    core-loop departure).
- **Reading (registered: H357.2 PASS with H357.1 unscored; see power note).**
  - A person-level state, carried by the residuals of PER1, PER2 and NR1D1 across a person's tissues, explains part of
    each tissue's departure from the body's clock state beyond matched genes. Descriptively, this state is the person's
    glucocorticoid-target state.
  - It also reaches the core-loop departure, so it is not only these three genes entering the projector.
  - The part is small: about 3 % of the departure (median), up to about 20 % in skeletal muscle. Most of the departure
    remains unexplained.
  - This fits the user's reading at a small scale: one systemic signal, read unequally by tissues, displaces their clock
    states differently (internal desynchrony at the signal layer). In deaths, a daily cortisol history and an agonal
    cortisol surge are not separated.
  - With A329 (the tissue's own glucocorticoid score does not explain its departure): the person-level estimate, averaged
    over many tissues, is the less noisy measure of the signal state; this reconciliation is post hoc.
- Code: `exploratory/A357_signal_departure/scripts/score_a357.py --stage score`.
<!-- src: paper/exploratory/A357_signal_departure/output/{score_tests,score_tissues}.csv (per-tissue aggregates) -->

## A329 erratum (2026-11-23; post hoc diagnosis, on the user's question why removing PER1, PER2 and NR1D1 moves the departure so close to random genes when the glucocorticoid state explains little of it): the fall of the departure excess on the core loop (2.6 -> 1.28 times matched) and of core-loop coupling is a change of projector noise gain, not the trio carrying the departure; the departure is a coherent clock-state shift, about 2.0 times matched noise in every tissue

- **Mechanism.** The projector satisfies L X = I, so a coherent clock-state departure (a shift of the whole waveform in
  phase or amplitude) projects the same with 14 or 11 genes, while gene-level variation is amplified by the projector's
  noise gain. Without NR1D1 (amplitude 0.32) and PER2 (0.22) the 11-gene projector has 3.1 times the noise gain of the
  14-gene one (unweighted trace(L L')). Matched genes, nearly pure gene-level variation, rise by about that factor
  (median own departure 0.709 -> 2.331); the clock genes rise less (1.893 -> 3.145). A ratio (C + g N) / (g N_matched)
  shrinks toward 1 as g grows even when the coherent part C is unchanged, so log ratios on projectors with different
  noise gains cannot be compared. The same added gene-level variation lowers coupling S / (S + D) on the core loop.
- **Split (per tissue, solving D_14 = C + N and D_core = C + g N with g from the tissue's matched sets; median g 3.13).**
  - Coherent clock-state departure C: median 1.42 (0.62 to 3.69), positive in 27 of 27 tissues.
  - Gene-level part of the clock genes N: median 0.47, against 0.71 for matched genes: the clock genes are not noisier
    gene by gene.
  - C / matched gene-level variation: median 2.0.
- **Corrected reading.** The clock-specific personal departure is a coherent shift of each tissue's clock state, present
  in every tissue; it is not carried by PER1, PER2 and NR1D1. This agrees with the A336 gene-drop null (the trio's drop of
  the H294b.1 ratio at the 32nd percentile of three-gene drops of similar weight) and with A357 (a person-level
  glucocorticoid-like state explains about 3 % of it). H329.3 is unaffected (within one projector).
- **Scope.** Within-projector tests against their own matched sets stand (A346, A348, A350, A321 verdicts). Differences
  between the 14-gene and core-loop values of noise-sensitive statistics (departure ratios, coupling, pairing r) are not
  attributable to the trio without a gene-drop null: the A346 endothelial "14-gene specificity rests on PER1, PER2 and
  NR1D1 (trio contribution 0.17)" and the A348 MuTHER fall (0.533 -> 0.311) are read as the core projector's loss of
  precision unless a gene-drop null says otherwise; A336 rows that carry the gene-drop null already make this
  comparison.
- Code: `exploratory/A329_core_departure/scripts/coherent_split.py` (from `tissues.csv`).
<!-- src: paper/exploratory/A329_core_departure/output/{coherent_split,coherent_summary}.csv -->

## A358 (scored 2026-11-24; DESCRIPTIVE; one run): no sign that human organ gain follows how evenly the organ weighs the drives in mouse (rho -0.6 against the -0.7 criterion, 5 organs): "not seen"

- Mouse weights b (mean of A322 / A324 where two exist) against the GTEx log gain log sqrt(S), 5 organs:

  | organ | b | b(1 - b) | amplitude kept under conflict | GTEx log gain |
  |---|---|---|---|---|
  | liver | 0.781 | 0.171 | 0.575 | 0.891 |
  | adipose | 0.549 | 0.248 | 0.215 | 0.785 |
  | lung | 0.371 | 0.233 | 0.495 | 0.566 |
  | heart | 0.233 | 0.179 | 0.589 | 0.971 |
  | kidney | 0.278 | 0.201 | 0.648 | 0.655 |

- Balance against gain: Spearman -0.6, exact one-sided p 0.175: below the registered criterion (<= -0.7), "not seen".
  Amplitude kept under conflict against gain (reported): +0.2 (p 0.39).
- **Reading.** The cross-species panel is not worth showing: the balance spans a narrow range (0.17-0.25) and the ranking
  does not follow it. Drive coherence is not excluded as a source of organ gain; this design cannot test it.
- **Use (user decision 2026-11-24, after scoring).** Mentioned as a suggestive pattern in section 4 (one line, with
  rho, n and the missed criterion), because it fits "conflicting drives cost gain"; never as evidence.
- Code: `exploratory/A358_A360_gain_sources/scripts/a358.py`.
<!-- src: paper/exploratory/A358_A360_gain_sources/output/{a358_organs,a358_tests}.csv -->

## A359 (scored 2026-11-24; BLIND for every reception slope; one run): louder organs do not read the person's systemic glucocorticoid state more strongly (H359.1 FAIL, rho 0.20, p 0.16); the immediate-early reception slope does track gain (rho 0.44, p 0.016; companion, not a registered hypothesis)

- GTEx fast deaths, 25 tissues; GC genes TSC22D3, KLF9, KLF15, ZBTB16, FKBP5, PDK4, TXNIP, SGK1, DDIT4 (NFIL3 and BHLHE40
  excluded as registered); 17 IEG genes; tissue score and person factor residualised on the leave-one-out host.
- The systemic GC state is strongly shared: the tissue's GC score correlates with the person factor of its other
  tissues at r 0.73-0.94 in every tissue except the adrenal (0.25).
- **H359.1 FAIL.** Spearman rho(log gain, GC reception slope) 0.199, one-sided permutation p 0.164.
- **H359.2 (companion).** IEG reception slope against gain: rho 0.435, p 0.016; difference GC - IEG -0.236 (90 % -0.735
  to +0.315).
- **Descriptive.** Adult NR3C1 level against gain: rho -0.07. GC reception against the departure D: -0.15; IEG: +0.03.
- **Reading (registered: FAIL).** Reception of the acute systemic glucocorticoid state does not track organ gain. The
  immediate-early result is a lead, not a finding: louder organs respond more to the person's acute
  stress / immediate-early state (in deaths, partly agonal); it was the specificity companion, chosen before scoring
  but not as a hypothesis. A registered re-test would need an independent design.
- Code: `exploratory/A358_A360_gain_sources/scripts/a359.py`.
<!-- src: paper/exploratory/A358_A360_gain_sources/output/{a359_tissues,a359_tests}.csv (per-tissue aggregates) -->

## A360 (scored 2026-11-24; BLIND for every person-covariate relation; one run): the person's clock volume falls with the person's glucocorticoid-target state (H360.1 PASS, partial R^2 0.061, p 0.0006) and is not a technical scale (H360.2 FAIL, 0.005, p 0.53)

- GTEx fast deaths, 237 donors with >= 5 of 25 tissues and complete covariates; log host radius from all tissues.
- **H360.1 PASS.** The systemic block (person GC and IEG states) explains a partial R^2 of 0.061 beyond the number of
  tissues, time of death, RIN and ischaemic time (permutation p 0.0006). Standardised coefficients: GC state -0.235 (90 %
  -0.352 to -0.109), IEG state +0.047 (-0.054 to +0.150).
- **H360.2 FAIL.** Technical block (median RIN, ischaemic time) partial R^2 0.005, p 0.53: no technical scale detected.
- Descriptive: age, sex and BMI add 0.007 (p 0.65). Full model R^2 0.081. Time of death: cos -0.042, sin +0.028 (both
  intervals include 0). Number of tissues +0.107 (-0.005 to +0.219).
- **Reading (registered: H360.1 PASS, H360.2 FAIL).** The person volume follows the body's systemic state: persons with
  a higher glucocorticoid-target state have quieter shared clocks across their organs. It is a small part (6 %) of the
  person volume, most of which stays unexplained; RIN and ischaemic time do not explain it, which answers A312's
  technical-scale caveat for these two measures. In deaths an acute (agonal) glucocorticoid surge and a chronic state
  are not separated. Not checked: whether the effect holds on the core loop (an acute glucocorticoid push moves PER1,
  PER2, NR1D1 and can change the host radius directly).
- Code: `exploratory/A358_A360_gain_sources/scripts/a360.py`.
<!-- src: paper/exploratory/A358_A360_gain_sources/output/{a360_coefficients,a360_tests}.csv (coefficients only) -->


## A312 calibration note (2026-11-24; post hoc, on the user's question why nearly every tissue's elasticity is above 1): the elasticity level is not anchored at 1, and a gradient of the observed size appears with no true gradient

- Simulation at the GTEx design (`exploratory/A312_elasticity/scripts/calib_level.py`): 27 tissues with the A294b S and
  D, 250 persons with a volume of SD 0.352 in log (A331's GTEx spread), every tissue scaling exactly in proportion
  (true elasticity 1), Poisson counts; six target tissues, 10 splits, 3 draws.
- Estimated elasticity median 1.11 (tissue means over draws 0.97 to 1.36; quiet tissues testis 1.35 and sun-exposed skin
  1.36, loud tissues 0.97 to 1.17); leverage median 0.75 (observed 0.78).
- **Reading (post hoc).** The estimator reads about 1.1 for proportional scaling (the yardstick is a host radius; my
  reading of the cause is the upward bias of a noisy radius, not shown), so the observed median 1.19 does not mean
  "more than in proportion". The H312.1 gradient (rho -0.34; low-gain 1.19, high-gain 1.07) has the size the
  simulation gives with no gradient: it is not shown beyond the estimator. A full null (all tissues, about 20 draws,
  the H312.1 rho under true elasticity 1) would decide whether H312.1 stands; it needs its own registration. What
  stands: persons differ in a volume shared by their organs (leverage), and no organ ignores it.
<!-- src: paper/exploratory/A312_elasticity/output/calib_level.csv -->

## A361 / A363 part 1 power note (2026-11-24; companions and analytic power only, no core-loop clock value computed): A361 is decidable at the A360 size (power 0.89), not at half of it (0.36); A363 part 1 is not powered at the registered size (0.76 < 0.8): not scored

- A360 reproduced exactly first (14 genes, registered model; 237 persons). One power run (1,442 s, peak 3.5 GB).
- Companions (200 matched, 200 responsive pseudo-clock hosts, same persons and covariates): GC coefficient on the core
  positions, median +0.081 matched (q05 -0.041) and +0.087 responsive (q05 -0.035); 14 positions +0.104 (q05 -0.039)
  and +0.089 (q05 -0.055). A generic donor scale does not fall with the GC state: it tends to rise.
- A361: SE 0.074 (A360's bootstrap); threshold -0.144 (the -1.96 SE bound, below the matched q05); power 0.891 at
  -0.235 (A360's size), 0.358 at half. Gate passed.
- A363 part 1: inflammatory score 176 genes (Hallmark TNF-alpha via NF-kB minus immediate-early, GC-target and clock
  genes); its R^2 on the A360 covariates 0.250 (collinear with the GC / IEG states); analytic SE 0.072; threshold -0.141;
  power 0.760 at partial R^2 0.03 (coefficient -0.192), 0.470 at 0.015. **Below the registered 0.8: not scored.**
  Companion inflammatory coefficients: median +0.123 (core matched, q05 -0.056), +0.164 (core responsive, q05 +0.001).
<!-- src: paper/exploratory/A361_A363_volume_mechanism/output/{systemic_power,systemic_companions}.csv -->

## A362 power note (2026-11-24; simulation only, no GTEx elasticity computed): decidable (power 1.00 for the pure and the half conflict pattern)

- 40 null cohorts (every tissue in proportion), 20 pure, 20 half conflict pattern; core waveform; one run (1,081 s).
- Null slope c: median 0.035, SD 0.236, q95 0.266. Pure pattern: median c 1.657, power 1.00; half: 0.865, power 1.00.
- Null bias of the mapped tissues' elasticity (true 1): 1.12 (liver) to 1.25 (kidney): the A312 calibration level
  (about 1.1) again; kidney reads highest under no gradient.
- Limit: the simulator's tissue noise reproduces the GTEx leverage (A312 calibration note: 0.75 against 0.78), so the
  null spread of c should be realistic.
<!-- src: paper/exploratory/A361_A363_volume_mechanism/output/{a362_power,a362_sims}.csv -->

## A363 part 2 power note (2026-11-24; simulation only, no GTEx gene-class elasticity computed): decidable on the simulator (power 1.00 at 0.5 and 0.25), but the null is probably too narrow for GTEx

- 40 null cohorts (uniform scaling), 20 with the suppressed class at a_p^0.5 extra, 20 at a_p^0.25; core classes
  (PER3, DBP, TEF against ARNTL, NPAS2, CLOCK); one run (338 s; a first run failed in the summary step, a column-access
  bug, after identical seeded simulations; fixed and rerun).
- Null T: median -0.001, q05 -0.012, q95 0.018. Alternatives: median T 0.304 (0.5), 0.144 (0.25); power 1.00 for both.
- **Caveat before scoring.** The simulator's per-gene noise (SD 0.1 on log scale, deep Poisson) was calibrated only at
  the level of the tissue state (leverage), not gene by gene; GTEx gene-level residuals are larger and gene-specific, so
  the real null spread of T is likely wider than 0.03. As registered, H363.2 would read T > 0.018 as cytokine-like; this
  threshold should be recalibrated (gene noise set to the GTEx gene residual variance, or a data-based null) in a dated
  addendum before the scoring run.
<!-- src: paper/exploratory/A361_A363_volume_mechanism/output/{a363_power,a363_sims}.csv -->

## A363 part 2 recalibrated power note (2026-11-24, addendum; GTEx noise calibration and simulation only, no gene-class elasticity computed): decidable at a_p^0.5 (power 1.00), not at a_p^0.25 (0.65); the null band is about six times wider than the registered one

- Calibration (core waveform, A360 rows): within-sample noise SD 0.29-0.47 per gene (median over tissues; DBP highest),
  person offset shared by a person's tissues SD 0.18-0.37 (DBP 0.37), against the registered simulator's independent
  SD 0.1. One run (363 s).
- Null T: median -0.006, q05 -0.087, q95 0.092 (registered simulator: -0.012 to 0.018). Alternatives: median T 0.279 at
  0.5 (power 1.00), 0.124 at 0.25 (power 0.65). Gate passed (power >= 0.8 at 0.5).
- Reading for scoring: only a strong cytokine-like distortion (the suppressed genes following the volume about 1.5 times
  as steeply) is reliably detectable; "uniform scaling" excludes that size, not a weaker one.
<!-- src: paper/exploratory/A361_A363_volume_mechanism/output/{a363_power_recal,a363_sims_recal,a363_noise_summary}.csv -->

## A361 (scored 2026-11-24; BLIND for every core-loop relation and companion coefficient; one run, 45 s): the person's clock volume falls with the glucocorticoid-target state in the core loop too, beyond every expression-matched and responsive set (H361.1 PASS)

- A360 reproduced exactly first (14 genes, 237 persons).
- **H361.1 PASS (core loop).** Systemic block partial R^2 0.057, permutation p 0.0014. GC coefficient -0.184 (90 %
  -0.297 to -0.063), below all 200 matched sets (q05 -0.041) and all 200 responsive sets (q05 -0.035); IEG +0.124.
- 14 genes -0.235; 14 minus core -0.052 (90 % -0.116 to +0.016): PER1, PER2 and NR1D1 add a little, not significantly.
- **Reading (registered).** The A360 link does not come from the glucocorticoid-target clock genes, and a generic donor
  scale does not show it (the companions' coefficients are centred on about +0.08). In deaths, acute (agonal) and
  chronic glucocorticoid states are still not separated.
- Code: `exploratory/A361_A363_volume_mechanism/scripts/a361_a363_systemic.py --stage score`.
<!-- src: paper/exploratory/A361_A363_volume_mechanism/output/{systemic_tests,systemic_companions}.csv -->

## A362 (scored 2026-11-24; BLIND for the core-loop tissue elasticities, NOT BLIND and post hoc in origin for the 14 genes; one scoring run, 121 s, after a first run failed in its summary step (the label "null" parsed as missing; nothing written or read); powered, 1.00): organs that lose most amplitude under conflicting drives in mouse do not follow the person's volume more steeply (H362.1 FAIL)

- A312 reproduced exactly first (14 genes, mapped tissues).
- **H362.1 FAIL (core).** Slope of the bias-corrected organ elasticity on the mouse sensitivity k_o: c -0.088 (simulated
  null q95 0.266; p 0.73); intercept 1.19; conflict share -0.06; Spearman -0.3. Corrected organ elasticities: liver
  1.08, adipose 1.19, lung 0.65, heart 0.95, kidney 1.75 (k 0.55, 1.54, 0.70, 0.53, 0.43). Muscle (raw) 1.32.
- 14 genes (reported, not blind): c 0.139 (p 0.29), conflict share 0.11; Spearman -0.4.
- **Reading (registered).** At this power (a pure conflict pattern gives c about 1.7, half of it 0.86) the organ pattern
  of cue conflict is not seen: the person volume is a factor common to the organs. Kidney (41 donors) follows it most
  and lung least, an order the mouse conflict data do not predict.
- Code: `exploratory/A361_A363_volume_mechanism/scripts/a362_organ_conflict.py --stage score`.
<!-- src: paper/exploratory/A361_A363_volume_mechanism/output/{a362_tests,a362_sims,a362_power}.csv -->

## A363 (scored 2026-11-24; BLIND for every inflammatory relation and gene-class elasticity; one run, 56 s; part 1 at power 0.76 (user, addendum), part 2 against the recalibrated null): no inflammatory correlate of the core-loop volume (H363.1 FAIL) and uniform scaling of the clock genes (H363.2): a low volume is not a cytokine-like damping at the detectable size

- **H363.1 FAIL (core).** Inflammatory state partial R^2 0.0046, p 0.30; coefficient -0.074 (90 % -0.243 to +0.083),
  3rd percentile of matched sets. GC coefficient in this model -0.160.
- 14 genes (reported): partial R^2 0.020, p 0.034, coefficient -0.156 (-0.282 to -0.045), below every matched set.
  The core interval (-0.243 to +0.083) contains the 14-gene value: whether PER1, PER2 and NR1D1 carry the correlate
  or the core estimate is only less precise is not decided (no paired difference or gene-drop null for part 1).
- **H363.2: uniform scaling (core).** T = e_supp - e_comp = -0.006 (recalibrated null -0.087 to 0.092); e_supp 0.659
  (PER3, DBP, TEF), e_comp 0.665 (ARNTL, NPAS2, CLOCK). 14 genes (reported) T 0.084. Per-gene elasticities
  (descriptive, core): 0.59-0.81 for nine genes, CRY2 1.16, RORA 1.60.
- **Reading (registered).** A cytokine-like damping that makes the suppressed genes follow the volume 1.5 times as
  steeply is excluded (power 1.00); a weaker distortion is not (0.65 at 1.25 times). The volume scales the core clock
  genes alike, as dephasing, composition or a uniform damping would. The 14-gene inflammatory correlate is not shown on
  the core; that the trio carries it (TNF-alpha and glucocorticoids act on it directly) is not shown either.
- Code: `exploratory/A361_A363_volume_mechanism/scripts/a361_a363_systemic.py`, `a363_gene_classes.py --stage score`.
<!-- src: paper/exploratory/A361_A363_volume_mechanism/output/{systemic_tests,a363_tests,a363_power_recal}.csv -->

## A361 / A363 gene-drop addendum (scored 2026-11-24; BLIND for the 364 drops and A363's paired difference; one run, 81 s): the trio drop is an ordinary three-gene drop for both coefficients: neither the smaller A361 nor the smaller A363 core coefficient is carried by PER1, PER2 and NR1D1 (not decided: precision or the removed genes)

| term | 14 genes | core (trio drop) | 364 drops: median (q10 to q90) | trio percentile, all / similar (77) | predicted by reliability loss |
|---|---|---|---|---|---|
| A361 GC | -0.235 | -0.184 | -0.221 (-0.244 to -0.175) | 85.2 / 51.9 | -0.231 |
| A363 inflammatory | -0.156 | -0.074 | -0.156 (-0.199 to -0.032) | 82.4 / 26.0 | -0.153 |

- A363 paired 14-minus-core difference -0.082 (90 % -0.187 to +0.036; 237 persons).
- **Reading (registered).** Removing any three genes moves both coefficients as much as removing the trio; the trio
  drop sits inside the random drops (below the 90th percentile, mid-range among drops of similar singular values).
- The reliability loss (factor 0.98) predicts far less movement than the drops show: which genes are left, not only how
  precise the projection is, moves a coefficient of this size. The yardstick for a core readout is the drop null, not
  the reliability factor (post hoc observation).
- Consequence for A363 part 1: the core test fails as registered (H363.1 FAIL), and the core value is a typical
  11-gene value, so the inflammatory correlate is neither shown on the clock nor shown to be a trio effect.
<!-- src: paper/exploratory/A361_A363_volume_mechanism/output/{drops_tests,drops_coefs}.csv -->

## A364 power note (2026-11-24; noise, simulation and companions only, no per-volunteer amplitude correlation computed): not decidable at this design (power 0.54 at the GTEx person spread, registered 0.8): not scored

- GSE326530, 20 volunteers x 2 skin sites x 4 times; one run (10 s).
- Core loop: per-sample noise variance 2.57 (P) and 1.38 (E) against a clock radius of 1.48 and 1.17 (noise SD above
  the clock radius at both sites). Expected r at the GTEx spread 0.419; threshold max(matched q95 0.401, simulated null
  q95 0.386); power 0.539 at log SD 0.352, 0.150 at half. 14 genes: expected 0.405, power 0.544.
- Companions: matched median 0.011 (q95 0.401), responsive median 0.120 (q95 0.416) on the core.
- **Reading (registered).** With 20 persons and four noisy samples per site, a GTEx-sized person volume would be seen
  about half the time: the living timed test is not decidable here; the living link stays open (and rests on del Olmo
  2022 and the Saint-Antoine & Anafi preprint).
<!-- src: paper/exploratory/A364_living_volume/output/{power,companions,sim}.csv -->

## A365 power note (2026-11-24; noise, simulation and companions only, no per-subject amplitude correlation computed): not decidable at this design (power 0.51 at the GTEx person spread): not scored

- GSE205155 (A37 download), 154 arrays, 11 subjects, 7 times; one run (15 s).
- Core loop: per-sample noise variance 0.50 (epidermis), 0.31 (dermis) against clock radii 0.78 and 0.70 (signal above
  noise, unlike A364). Expected r at the GTEx spread 0.674, but the matched sets share radius between the two layers
  of one punch (median 0.206, q95 0.668; responsive q95 0.679): threshold 0.668, power 0.512 (0.131 at half). 14 genes:
  expected 0.685, matched q95 0.626, power 0.611.
- **Reading (registered).** With 11 subjects and two layers of one biopsy, a GTEx-sized volume would be separated from
  the shared biopsy scale about half the time: not decidable. The preprint's R^2 0.649 on these data has no such
  companion; our power stage shows that expression-matched genes reach r 0.67 in one punch one time in twenty.
<!-- src: paper/exploratory/A364_living_volume/output/{a365_power,a365_companions,a365_sim}.csv -->

## A366 power note (2026-11-24; geometry, timing-only simulation and companions only, no clock radius statistic computed): MuTHER's clock states lie on an arc (combined elongation 3.4; phase spread about 1.9 h); decidable (power 1.00 at the GTEx spread, 0.92 at half)

- Geometry: elongation fat 3.78, skin 1.72, combined 3.39 (a ring gives about 1): the biopsies sit in a window of the
  day, as expected. Matched by a phase SD of 0.50 rad (1.9 h) in the ring model.
- Timing only (no loudness differences between women), at this geometry: the raw A331 statistic reads median 0.279
  (q95 0.338); the circle-centred statistic median 0.151 (q95 0.209). A331's observed raw 0.177 is below what shared
  timing alone produces in this model (not a test; the ring model's noise and phase distribution are assumptions).
- Companions (circle-centred): matched median 0.029 (q95 0.100), responsive median 0.086 (q95 0.132).
- Threshold max(0.100, 0.209) = 0.209; expected circle-centred leverage with a GTEx-size volume 0.314; power 0.998,
  0.918 at half. Gate passed. One run (10 s).
<!-- src: paper/exploratory/A366_muther_timing/output/{power,geometry,companions,sim}.csv -->

## A366 (scored 2026-11-24; BLIND for the circle-centred leverage and the phase-radius relation; one run, 10 s; powered, 1.00): A331's living fat-skin loudness is not separated from the shared biopsy time (H366.1 FAIL); in the ring model the observed value is what timing alone produces, and well below a GTEx-size or half-size volume

- **H366.1 FAIL.** Circle-centred leverage 0.165 (raw A331 0.177), persons-permuted p 0.0005, above 99.0 % of matched
  and 99.5 % of responsive sets, but below the timing-only q95 0.209. It sits at the 64th percentile of the timing-only
  simulations (median 0.151) and below 99 % of the half-GTEx-size simulations (median 0.262; GTEx size 0.314).
- **Phase and amplitude (the user's question).** A tissue's centred radius rises with the other tissue's distance from
  the cohort's mean phase: Spearman 0.214 (fat radius vs skin position), 0.186 (skin radius vs fat position): the
  timing signature. Circle centres sit 0.19 (fat) and 0.08 (skin) of the mean radius away from the cohort mean.
- **Reading (registered).** A331's leverage is explained by shared biopsy time in this model: A331 is withdrawn as
  person-volume evidence. Its pass over random genes reflects that clock (and cortisol-driven) genes move with the time
  of day, not a shared loudness. The living person volume is not shown in public data (A364, A365 not decidable; A366
  fails); the living link rests on del Olmo 2022 and Saint-Antoine & Anafi 2025 (no matched-gene control).
- **Limits.** The null assumes a sinusoidal ring, normal phase spread and isotropic noise at A331's clock fractions; a
  volume below half the GTEx size is not excluded.
- Code: `exploratory/A366_muther_timing/scripts/a366_muther_timing.py --stage score`.
<!-- src: paper/exploratory/A366_muther_timing/output/{tests,sim,power,companions}.csv -->

## A366 note (2026-11-24; post hoc, synthetic, on the user's remark that a circle centre is hard to find on a narrow arc): the circle fit did not find the ring's centre at MuTHER's noise; the A366 conclusion stands, its label overstates

- At MuTHER's spread (phase SD 1.9 h) the true ring centre lies 1.37 mean radii from the cohort mean. Without noise the
  algebraic fit finds it exactly. With A331's noise it removes only a small part of the distance (remaining error 0.87 of
  it; 0.71 at 3.8 h, 0.66 at 5.7 h). This matches the observed centre shifts (0.19 fat, 0.08 skin).
- The A366 statistic was therefore close to A331's own: radii centred on the cohort mean, residualised on two harmonics
  of both tissues' angles. What A366 decided is whether the observed shared radius exceeds what phase position alone
  produces at the observed geometry. The timing-only and volume simulations went through the same estimator, so
  H366.1 FAIL (powered) stands. "Circle-centred" should read "timing-referenced": no true phases or amplitudes were
  recovered.
- General: absolute phase and amplitude from a hospital-hours arc need a wide spread or low noise; around the clock
  (phase SD 9.5 h) the cohort mean is already the centre (0.04 mean radii away). Relative quantities stay the readout
  for narrow windows.
- Code: `exploratory/A366_muther_timing/scripts/a366_kasa_check.py`.
<!-- src: paper/exploratory/A366_muther_timing/output/kasa_check.csv -->

## A368 power note (2026-11-24; no clock Q computed; 6 s): decidable (power 1.00 at GTEx size, 0.87 at half size); the threshold is set by the matched genes, not by timing

- Timing-only Q at A366's spread (phase SD 0.5 rad, 680 women): median 0.102, q95 0.163; at 0.75 x / 1.25 x the spread
  q95 0.116 / 0.239; uniform phases q95 0.125.
- Companions: matched median 0.085, q95 0.253; responsive median 0.149, q95 0.263. Threshold max(matched q95, largest
  timing-only q95) = 0.253 (matched).
- Volume simulations: GTEx size median 0.458 (power 1.00), half size median 0.292 (power 0.87).
- Code: `exploratory/A368_muther_radial/scripts/a368_muther_radial.py --stage power`.
<!-- src: paper/exploratory/A368_muther_radial/output/{power,sim,companions}.csv -->

## A368 (scored 2026-11-24; BLIND for Q; one run, 3 s; powered, 1.00 / 0.87): at matched arc position, living fat and skin share amplitude beyond biopsy time, but not beyond expression-matched random genes (H368.1 FAIL)

- **H368.1 FAIL.** Q 0.228, persons-permuted p 0.0005, below the threshold 0.253 (matched q95): 93rd percentile of
  matched and 87th of responsive sets.
- **Against timing.** Q is above all 500 timing-only simulations at the fitted spread (q95 0.163). It is below the
  timing-only q95 at 1.25 x the spread (0.239).
- **Against a volume.** Q sits at the 4th percentile of the half-GTEx-size simulations (median 0.292) and below all
  GTEx-size ones (median 0.458).
- **Reading (registered).** No shared clock amplitude at matched phase beyond random genes. A GTEx-size volume is
  excluded, and a half-size one nearly so. Post hoc: what fat and skin share at matched arc position is mostly a
  within-woman scale that random genes share too (array batch, biopsy day, composition). A clock-specific part smaller
  than half the GTEx size is not excluded.
- **Status of the living volume.** Not shown in MuTHER by either statistic (A366, A368). In deaths, the volume sits far
  above its companions (A336 core 0.714 against matched q95 0.621). Its absence at this size in living women is now an
  observation, not only a power limit; A367 asks whether the agonal phase makes it.
- Code: `exploratory/A368_muther_radial/scripts/a368_muther_radial.py --stage score`.
<!-- src: paper/exploratory/A368_muther_radial/output/{tests,sim,power,companions}.csv -->
