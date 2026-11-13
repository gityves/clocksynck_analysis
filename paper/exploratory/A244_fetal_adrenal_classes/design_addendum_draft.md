# A244 design addendum (DRAFT, written 2026-11-13 by the preparing agent; not committed; no fetal clock statistic computed)

Nothing below uses the 14 clock genes in a classifier, a coherence, a pairing, a phase or a companion on fetal data. The
only fetal clock-gene numbers are sums and medians of counts (detection), as in the A242 detection table.

Code: `scripts/prep_a244.py` (stages markers, hist, classes, pseudobulk, design), `scripts/design_a244.py`,
`scripts/feasibility_a244.py`. Outputs: `output/` (tables copied below); the per-cell table and the class pseudobulks are in
`$SCCC_DATA/natexp/fetal_adrenal_open/a244_{cells,class_pseudobulk}.parquet`, `a244_class_samples.csv`.

## 1. Per-cell rule and thresholds (set on the marker counts alone)

Per cell: total UMI count and the counts of the registered marker lists (cortex STAR CYP11A1 CYP17A1 CYP21A2 CYP11B1
NR5A1; medulla CHGA CHGB TH DBH PHOX2B STMN2 ELAVL4 SOX10 PLP1; stroma COL1A1 COL1A2 COL3A1 PDGFRB PECAM1 CDH5; blood
HBG1 HBG2 HBA1 HBB PTPRC); all 26 symbols are present in every source. The marker fraction of a class is its marker
UMIs divided by the total UMIs of the cell. **A class is met when its marker count is >= 3 and its marker fraction is
>= 0.005 (cortex, medulla, stroma) or >= 0.1 (blood).** A cell is assigned to the one of cortex / medulla / stroma that
it meets, provided it meets no other of the four classes; a cell that meets two classes, none, or blood is left out.
One rule, the same constants for every cohort.

Why these constants. (i) The histograms of log10 marker fraction (`output/marker_hist.csv`) of the medulla and stroma
markers show a background mode (ambient RNA, <= 0.002) and a cell-type mode (0.01-0.02) with the valley near 0.005
in the cohorts with a valley (kildisiute, kameneva, wang, sasaki). The cortex fraction has no valley in the pooled
histogram (the modes are 0.003 in dong, 0.02 in kildisiute: platform sensitivity differs). (ii) Two- and three-component
Gaussian mixtures on log10 fraction (cells with count >= 3, pooled) put the class boundary at 0.0026-0.0029 (two
components) and 0.0069 (cortex), 0.0047 (medulla), 0.0042 (stroma) (three components); 0.005 lies in that range for the
three classes, so one constant is used. (iii) The blood fraction is bimodal with erythroid cells at 0.5-0.8 and a
plateau of ambient haemoglobin up to 0.3; 0.1 leaves out every cell above the plateau's lower edge (conservative: more
cells are left out). (iv) The minimum count of 3 keeps a handful of ambient UMIs in a shallow droplet from calling a
class. The published annotations were looked at only after these constants were fixed (the confusion tables below were
the first time the classes were compared with them) and did not change them.

Cells after the rule (all cohorts, 448,931 cells): cortex 71,257, medulla 29,227, stroma
40,736, left out 307,711 (68.5%).
Per cohort:

| cohort | left_out | cortex | medulla | stroma | cells |
|---|---|---|---|---|---|
| dong | 34953 | 18710 | 2941 | 2856 | 59460 |
| han | 22112 | 3371 | 1968 | 2549 | 30000 |
| kameneva | 95059 | 12326 | 4088 | 7246 | 118719 |
| kildisiute | 32097 | 17353 | 5562 | 2960 | 57972 |
| sasaki | 50135 | 8671 | 2514 | 5704 | 67024 |
| wang | 73355 | 10826 | 12154 | 19421 | 115756 |

Whole-adrenal check: the sum of all cells of each fetus (classes plus left out) equals the A242 pseudobulk (asserted in
the script, rtol 1e-3).

### Agreement with the published annotation (Kildisiute: every cell; Dong: 83 % of the cells; kameneva, wang, sasaki, han: no per-cell annotation on disk)

Kildisiute (rows: published label, columns: this rule):

| published | (left out) | cortex | medulla | stroma |
|---|---|---|---|---|
| Bridge | 239 | 0 | 524 | 0 |
| Chromaffin | 18 | 0 | 2035 | 0 |
| Cortex | 10361 | 17189 | 1 | 1 |
| Endothelium | 8001 | 120 | 1 | 32 |
| Erythroblasts | 7910 | 10 | 0 | 0 |
| Leukocytes | 2407 | 26 | 2 | 3 |
| Mesenchyme | 1841 | 3 | 2 | 2897 |
| Other | 711 | 2 | 0 | 1 |
| SCPs | 591 | 3 | 21 | 26 |
| Sympathoblastic | 18 | 0 | 2976 | 0 |

Purity (share of cells of a class carrying the matching published label; cortex = Cortex; medulla = Chromaffin,
Sympathoblastic, SCPs, Bridge; stroma = Mesenchyme, Endothelium): cortex 0.991, medulla 0.999, stroma 0.990. Recall of
published labels: cortex 0.624, medulla 0.861, stroma 0.227 (the endothelial cells mostly carry too few of the six
stroma markers). Dong (rows: published label):

| published | (left out) | cortex | medulla | stroma |
|---|---|---|---|---|
| Capsular cells | 2978 | 894 | 10 | 8 |
| Chromaffin cells | 375 | 70 | 1474 | 2 |
| Endothelium | 6341 | 1878 | 2 | 13 |
| Fibroblast | 2202 | 598 | 17 | 2000 |
| Myeloid cells | 3014 | 424 | 47 | 79 |
| SCPs | 459 | 116 | 44 | 31 |
| Steroidogenic cells | 11151 | 13419 | 0 | 7 |
| Sympathoblasts | 132 | 30 | 932 | 0 |
| T cells | 331 | 97 | 18 | 16 |

Dong purity: cortex 0.766 (the cortex class holds 1,878 endothelial, 894 capsular, 598 fibroblast, 424 myeloid cells with
ambient steroidogenic RNA: Dong's counts are not ambient-corrected), medulla 0.963 (Chromaffin cells, Sympathoblasts,
SCPs), stroma 0.934 (Fibroblast, Endothelium); recall: cortex 0.546, medulla 0.668, stroma 0.154. Cells left out are not
a bias in themselves (the rule is conservative), but the classes are the cells with enough marker UMIs.

## 2. Design (fetuses with >= 200 cells of the class)

39 fetuses have at least one class; cortex 34, medulla 30, stroma 33; both cortex and medulla: 26 (H244.2). Depth =
median pseudobulk depth in millions of UMIs; clock_gene_median_count = median over the 14 genes of the cohort median
count (`output/clock_detection_classes.csv` has the per-gene medians); clock_total_median = median clock-gene sum; a =
median foreign fraction of the other class's markers (section 3).

| compartment | cohort | n_fetuses | cells_total | cells_median | cells_min | depth_median_M | clock_gene_median_count | clock_total_median | a_median |
|---|---|---|---|---|---|---|---|---|---|
| cortex | dong | 3 | 18617 | 7.2e+03 | 323 | 37.2 | 296 | 7.2e+03 | 0.0178 |
| cortex | han | 2 | 3371 | 1.69e+03 | 667 | 1.54 | 19.5 | 271 | 0.0182 |
| cortex | kameneva | 10 | 12324 | 920 | 244 | 7.74 | 87.8 | 2.38e+03 | 0.0036 |
| cortex | kildisiute | 5 | 17353 | 1.53e+03 | 357 | 29.6 | 334 | 7.35e+03 | 0.0018 |
| cortex | sasaki | 5 | 8314 | 872 | 317 | 20.3 | 143 | 4.06e+03 | 0.0057 |
| cortex | wang | 9 | 10694 | 888 | 308 | 12.8 | 94.5 | 2.58e+03 | 0.01 |
| medulla | dong | 4 | 2941 | 717 | 359 | 4.68 | 32 | 514 | 0.0806 |
| medulla | han | 2 | 1968 | 984 | 828 | 0.66 | 4.75 | 77.5 | 0.0393 |
| medulla | kameneva | 7 | 3702 | 302 | 208 | 2.97 | 26.5 | 561 | 0.016 |
| medulla | kildisiute | 3 | 5307 | 2.22e+03 | 383 | 23.1 | 118 | 1.87e+03 | 0.0046 |
| medulla | sasaki | 4 | 2437 | 610 | 253 | 7.87 | 66 | 1.33e+03 | 0.0278 |
| medulla | wang | 10 | 12154 | 832 | 216 | 11 | 102 | 1.6e+03 | 0.0158 |
| stroma | dong | 3 | 2670 | 1.12e+03 | 279 | 3.66 | 28.5 | 655 |  |
| stroma | han | 2 | 2549 | 1.27e+03 | 532 | 1 | 5 | 140 |  |
| stroma | kameneva | 8 | 6826 | 736 | 337 | 6.65 | 51 | 1.16e+03 |  |
| stroma | kildisiute | 4 | 2838 | 478 | 285 | 3.21 | 20.8 | 638 |  |
| stroma | sasaki | 6 | 5596 | 764 | 537 | 7.75 | 73 | 1.78e+03 |  |
| stroma | wang | 10 | 19421 | 2.15e+03 | 222 | 14.8 | 99.8 | 2.21e+03 |  |

Readouts as registered: pooled over cohorts with >= 2 fetuses, cohort indicator and quadratic age trend, permutations
within cohort: `cortex` 34 fetuses in 6 cohorts, `medulla` 30 in 6 (7 and 6 covariates).

## 3. Foreign RNA (the A190 limit), per cohort and class

a = (other-class marker counts / depth in the class pseudobulk) / (the same in the other class's pseudobulk of the fetus,
any size); n_ref_ge_min = fetuses whose reference pseudobulk also has >= 200 cells. share = expected foreign clock-gene
counts / clock-gene counts, median over the 14 genes of the pooled cohort (`foreign_cohort.csv`; per fetus in
`foreign_fetus.csv`). Row `all` pools every fetus of the class.

| cohort | compartment | n_fetuses | n_ref_ge_min | a_median | a_min | a_max | share_median_over_genes | share_max_gene | share_max_gene_name |
|---|---|---|---|---|---|---|---|---|---|
| dong | cortex | 3 | 3 | 0.0178 | 0.0158 | 0.0249 | 0.0171 | 0.0479 | PER3 |
| dong | medulla | 4 | 3 | 0.0806 | 0.0121 | 0.168 | 0.0603 | 0.22 | NR1D1 |
| han | cortex | 2 | 2 | 0.0182 | 0.0128 | 0.0237 | 0.0109 | 0.0437 | PER3 |
| han | medulla | 2 | 2 | 0.0393 | 0.0225 | 0.0562 | 0.0321 | 0.693 | NR1D1 |
| kameneva | cortex | 10 | 7 | 0.0036 | 0.0026 | 0.0196 | 0.0064 | 0.0236 | PER2 |
| kameneva | medulla | 7 | 7 | 0.016 | 0.0078 | 0.0323 | 0.0124 | 0.0726 | DBP |
| kildisiute | cortex | 5 | 3 | 0.0018 | 0.0011 | 0.0053 | 0.0021 | 0.0136 | PER3 |
| kildisiute | medulla | 3 | 3 | 0.0046 | 0.0038 | 0.0063 | 0.004 | 0.0283 | DBP |
| sasaki | cortex | 4 | 2 | 0.0057 | 0.0045 | 0.0283 | 0.0099 | 0.0337 | PER3 |
| sasaki | medulla | 4 | 2 | 0.0278 | 0.0243 | 0.0567 | 0.0243 | 0.102 | DBP |
| wang | cortex | 9 | 9 | 0.01 | 0.0057 | 0.0196 | 0.02 | 0.044 | PER2 |
| wang | medulla | 10 | 9 | 0.0158 | 0.0088 | 0.0466 | 0.01 | 0.0608 | NR1D1 |
| all | cortex | 33 | 26 | 0.0064 | 0.0011 | 0.0283 | 0.009 | 0.0295 | PER3 |
| all | medulla | 30 | 26 | 0.0164 | 0.0038 | 0.168 | 0.0127 | 0.0711 | NR1D1 |

No class is near the confounding limit (median foreign share 0.5): `all` cortex 0.009, medulla 0.013. The largest
single-gene share is 0.693 (NR1D1, han medulla, 2 fetuses); the registered rule uses the median.

## 4. Adult feasibility (`output/feasibility_passrates.csv`, 200 subsamples, 500 permutations, 533 s, 8 workers)

Adult stand-ins: 93 GTEx adrenal glands (death 08-18 h, any Hardy class), thinned to the clock-gene counts of each
class and cohort (share of capped targets: cortex 0.143, medulla 0.000, stroma 0.000; two_classes 0.040). Pass rate =
share of subsamples with p < 0.05.

| test | a1 | a0.75 | a0.5 | shuffled | mispaired | a1 CI | shuffled CI | mispaired CI |
|---|---|---|---|---|---|---|---|---|
| cortex | 1 | 0.995 | 0.52 | 0.055 |  | 0.982-1.0 | 0.028-0.096 |  |
| medulla | 1 | 0.91 | 0.335 | 0.055 |  | 0.982-1.0 | 0.028-0.096 |  |
| stroma | 1 | 0.985 | 0.395 | 0.045 |  | 0.982-1.0 | 0.021-0.084 |  |
| cortex:late | 0.91 | 0.67 | 0.18 | 0.04 |  | 0.861-0.946 | 0.017-0.077 |  |
| cortex:early | 1 | 0.92 | 0.315 | 0.045 |  | 0.982-1.0 | 0.021-0.084 |  |
| medulla:late | 0.57 | 0.29 | 0.135 | 0.06 |  | 0.498-0.64 | 0.031-0.102 |  |
| medulla:early | 0.99 | 0.845 | 0.26 | 0.045 |  | 0.964-0.999 | 0.021-0.084 |  |
| stroma:late | 0.65 | 0.365 | 0.165 | 0.045 |  | 0.58-0.716 | 0.021-0.084 |  |
| stroma:early | 0.995 | 0.925 | 0.33 | 0.05 |  | 0.972-1.0 | 0.024-0.09 |  |
| two_classes | 0.395 | 0.34 | 0.19 |  | 0.04 | 0.327-0.466 |  | 0.017-0.077 |

Design rows (`feasibility_design.csv`):

| test | level | n_fetuses | n_samples | cohorts | n_covariates | pool_donors | share_capped |
|---|---|---|---|---|---|---|---|
| cortex | L1 | 34 | 34 | 6 | 7 | 93 | 0.143 |
| cortex:late | L1 | 13 | 13 | 5 | 5 | 93 | 0.093 |
| cortex:early | L1 | 19 | 19 | 4 | 5 | 93 | 0.06 |
| medulla | L1 | 30 | 30 | 6 | 7 | 93 | 0 |
| medulla:late | L1 | 9 | 9 | 4 | 4 | 93 | 0 |
| medulla:early | L1 | 21 | 21 | 5 | 6 | 93 | 0 |
| stroma | L1 | 33 | 33 | 6 | 7 | 93 | 0 |
| stroma:late | L1 | 9 | 9 | 4 | 4 | 93 | 0 |
| stroma:early | L1 | 21 | 21 | 4 | 5 | 93 | 0 |
| two_classes | L2 | 26 | 52 | 6 | 7 | 80 | 0.04 |

**Rule (>= 0.8 at the full clock): `cortex` 1.000 and `medulla` 1.000 are scored; `stroma` is descriptive.** The
false-positive rate with every clock gene shuffled is 0.055 (0.028-0.096) cortex and 0.055 (0.028-0.096) medulla (0.05
inside the interval). At 0.75 of the adult clock: cortex 0.995, medulla 0.910; at 0.5: 0.520, 0.335.

## 5. H244.2, second adult tissue and calibration

Second tissue: **Muscle - Skeletal**. Rule fixed before any fetal statistic: among the GTEx tissues of the A241 stand-in
list (Lung, Heart - Left Ventricle, Liver, Brain - Cerebellum, Brain - Frontal Cortex, Muscle - Skeletal) the one
with the most donors sharing an adrenal gland in the 08-18 h window (muscle 80 of 93; lung 67, left ventricle 63; whole
blood also 80, excluded: no coherent core clock in blood). Class cortex stands for the adrenal gland, medulla for muscle
(the roles are arbitrary: only the pairing within a donor matters). Gate: `development.shared` on the cohort- and
age-adjusted counts, the cohorts as strata (pairing permuted within cohort, `Dataset.gate` stratum column; no change to
`development.py` needed), p = the larger of p_pairing and p_acrophase, 500 permutations; 26 fetuses in 6 cohorts, 80 adult
donors with both tissues, distinct donors per subsample. Mis-paired: the second tissue taken from another drawn donor.
Rows `two_classes` above: pass rate at the full clock 0.395
(95 % CI 0.327-0.466), 0.75: 0.340, 0.5: 0.190; false-positive rate with mis-paired
donors 0.040 (0.017-0.077; 0.05 inside). **The pass rate is below 0.8: H244.2 is "not estimable" and is not computed.**

## 6. Readouts

Scored (adult pass rate at the full clock >= 0.8): H244.1 `cortex`, `medulla` (Holm over the two; the corrected readouts
`cortex:corrected`, `medulla:corrected` follow the rule of the registration). Not estimable: H244.2. Descriptive as
registered: `stroma`; the windows `late` / `early` of each class (adult pass rate at the full clock: cortex late 0.910,
early 1.000; medulla late 0.570, early 0.990; stroma late 0.650, early 0.995); without each cohort in turn; unsigned R;
leave-one-gene-out; position among the adult subsamples; depth-matched comparison. The medulla `late` window (9
fetuses in 4 cohorts) has an adult pass rate of 0.57 and is to be read with that.

The corrected readouts remove expected foreign counts from the class pseudobulks; the adult calibration above is on the
uncorrected counts (the median foreign share of the clock-gene counts is 0.002-0.060 per cohort and class).
