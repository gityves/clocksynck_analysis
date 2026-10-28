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
