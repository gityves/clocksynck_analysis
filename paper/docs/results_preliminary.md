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
