# A138 (DRAFT for PI review, 2026-10-25; not registered, nothing scored): what can be estimated about the clock in tumours?

Replaces WP5 (A114). Motivation (PI, Q9 e): published untimed orderings of tumours (CoFE, Ananthasubramaniam et al.,
PLoS Biol 2025; CYCLOPS / CYCLOPS 2.0 applications) may reflect the method rather than tumour clocks. Inputs:
`literature.md` (claims checked by the lead against the CoFE and CYCLOPS full texts, 2026-10-25) and `inventory.md`.
Written after A52, A83, A86, A87, A97, A134 and the GTEx method benchmark were scored (their results motivate the
design; none of the readouts below has been computed).

Checked citations:
- **CoFE:** "Successful ordering is signalled by COFE with a clear minimum in the CV output"; sample times "restricted to
  a fraction of the cycle" are "flagged by a CV output with no non-trivial minimum" (synthetic data, S2 Fig).
- **CoFE on TCGA:** "predicted time labels lay throughout the circadian cycle in all ACs"; "the clocks in tumors do not
  have a fixed timing relationship with clocks in healthy tissue"; 2,000-6,000 rhythmic genes; S and G2M scores
  oscillate. Biopsies are assumed taken by day, and no non-circadian (no-clock) control was run on real data.
- **CYCLOPS (Anafi 2017):** ">250 samples are required"; "It requires data from the entire periodic cycle to form an
  ellipse"; seeded ordering of hepatocellular carcinoma "was not able to generate a quality ordering".

Our GTEx benchmark (timed truth, n = 80; CoFE as published genome-wide in 3 draws, CoFE on the clock genes only in 10):
median absolute per-gene acrophase error 3.10 h (full day) / 3.30 h (08-18 h window) genome-wide vs 1.53 / 2.30 h on
the clock genes only; tissue amplitude error 0.258 / 0.395 vs 0.115 / 0.116; failure rate 0.19 / 0.62 vs 0.24 / 0.14.
Genome-wide CoFE recovers core-clock acrophases poorly even when samples cover the whole day, and does better when
restricted to the clock genes (PI, 2026-10-25). The TCGA tumour clock-gene acrophases CoFE reports (e.g. NR1D2 about
8 h late) are read against this error.
<!-- src: paper/results/figure_1/output/method/gtex_benchmark/{genes,summary}.csv (M3 run 2026-10-25) -->

CoFE's tumour core clock (their Fig 2D, consensus over the 11 TCGA cancers; the authors: "arranged ... in the expected
order", "like in a healthy functional clock"). Lead's reading of the published image (approximate, low resolution; exact
values in their S3 Table, not yet read): ARNTL near 17 h, with PER2, PER3, TEF, CRY2, BHLHE40, NR1D2, CLOCK, RORA within
about 3 h of it; DBP near 7 h, NR1D1 near 10 h. In healthy human tissue (GTEx reference waveform, `cs.gtex()`) PER2,
PER3, TEF, DBP and NR1D2 peak 9.3-14.2 h from ARNTL (ARNTL 21.0 h; PER3 9.7, TEF 9.2, DBP 8.5, PER2 11.7, NR1D2 6.8 h).
The PI recalls DBP peaking with BMAL1 in a CoFE figure (Fig 1?); in Fig 2D as read here DBP is about 10 h from ARNTL,
but the E-box / D-box targets PER2, PER3 and TEF sit with ARNTL instead of opposite it. To be confirmed from S3 Table
before any statement.

## Q1. Does the malignant clock carry an estimable, host-coupled signal? (single cell, gate-passing cohorts)
- **H138.1:** pooled anchored coupling of the malignant compartment (copy-number-called malignant cells where available,
  else tumour epithelium) > 0: lower bound of the 95 % patient-bootstrap interval > 0 (DerSimonian-Laird over cohorts,
  as A17). If it fails, the malignant clock state is reported as not estimable at these designs, not as "decoupled".
- **D138.1 (descriptive):** decomposition of the malignant clock-plane variance per cohort, from the anchored fit and
  the per-cell data, into (i) host-coupled (c^2), (ii) patient-specific (psi minus Poisson), (iii) between cells within
  patient (per-cell dispersion D on the same scale, as variance of cell states around the patient mean), (iv) Poisson.
  The same for normal epithelium and stroma of the same patients (specificity companion). Estimand of interest: (ii)
  malignant minus normal epithelium.
- **Power (before scoring):** MDE of the coupling at each cohort's design (`clocksync.power.mde`) and of the pooled
  estimate; the smallest attainable p of each test is stated.

## Q2. When do untimed orderings report a clock that is not there? (simulation, then TCGA placement)
Simulator (new, local first; clocksync once validated): genome-wide log-expression of tumours, n in {50, 100, 300, 500},
p = 2,000 genes (clock module of the 15 CoFE clock genes at GTEx waveforms plus 15 % rhythmic followers), non-clock
structure calibrated to TCGA (purity mixing with stroma, a proliferation programme, subtype mixture, copy-number blocks:
leading non-clock principal components of each TCGA cancer). Arms:
- **null:** no clock (clock genes at their mesor plus the non-clock structure);
- **window:** an intact clock sampled in 08-18 h;
- **full day:** an intact clock sampled over 24 h, with patient phase SD in {0, 2, 4 h};
- **coupled tumour:** a patient's tumour clock = its normal-tissue clock + a common offset + patient SD in {0, 2, 4 h}.
50 replicates per cell; CoFE as published (CV, 5 restarts) and CYCLOPS 2.0 with its default seeds.
- **H138.2 (false validity):** in the null arm at n = 300, CoFE reports an interior CV minimum in > 20 % of replicates
  (one-sided binomial test against 0.20, p < 0.05). Also reported: the same in the window arm, and whether predicted
  phases cover the cycle (Rayleigh R of predicted phases) in each arm.
- **H138.3 (uncoupled inference):** in the coupled-tumour arm (window 08-18 h, patient SD 2 h, n = 300), the
  tumour-minus-normal predicted-time difference is indistinguishable from uniform (Rayleigh p > 0.05) in > 50 % of
  replicates. Pass means CoFE's "no fixed timing relationship" does not follow from its readout.
- **H138.6 (genome vs clock genes):** in the full-day and window arms (n = 300, patient SD 2 h), the median per-gene
  acrophase error of the clock genes is larger with CoFE genome-wide than with CoFE on the clock genes only (paired over
  replicates, one-sided Wilcoxon p < 0.05), replicating the GTEx benchmark on tumour-like structure; and in the null and
  window arms, the share of clock genes whose acrophase deviates >= 8 h from the truth (the size of CoFE's NR1D2 shift)
  under genome-wide CoFE is reported.
- **D138.3 (TCGA, descriptive):** per cancer, agreement between the genome-wide CoFE ordering and the clock-gene-only
  ordering (circular correlation), and the clock-gene acrophase order of each against the GTEx healthy reference
  (`benchmark/tcga.py` `GTEX_MATCH`); which published deviations (NR1D2 +8 h, TEF, BHLHE40, PER2) survive the clock-only
  ordering; the ARNTL-to-target acrophase differences (DBP, TEF, PER2, PER3, NR1D1, NR1D2) against GTEx, and in the
  simulation the same differences when the truth is a healthy clock (does genome-wide CoFE pull targets toward ARNTL?).
- **D138.2 (descriptive):** where each TCGA cancer falls on the simulation map (CV curve depth, phase coverage, gene
  rhythmic fraction), and the same for CYCLOPS 2.0 on the TCGA cancers where it runs; recovery error with known truth
  in every arm (mean absolute error, `evaluate.py`).

## Q3. Does the malignant clock follow the patient's own normal tissue? (paired designs)
- **H138.4:** TCGA paired tumour-normal (10 cancers, 517 patients): circular correlation of the tumour clock state with
  the same patient's normal state exceeds the correlation obtained when normals are permuted across patients within
  cancer (10,000 permutations, pooled over cancers, one-sided p < 0.05). Specificity companion on the same permutations:
  pseudo-clock gene sets matched on expression (A52 rule). Gate: the cancer's normals pass the clock gate.
- **H138.5:** the same in single-cell tumour vs adjacent-normal epithelium (and malignant vs non-malignant tumour
  epithelium), host-referenced, cohorts with >= 6 paired patients; companion: tumour stroma vs adjacent stroma.
- **Power (before scoring):** null SD and MDE of the correlation at each design.

## Order and compute
Q3 and Q1 first (data and code exist; hours on 4 workers); Q2 after the simulator passes a check that, in the full-day arm
without heterogeneity, CoFE and CYCLOPS recover truth as published. CoFE at p = 2,000 and CYCLOPS take minutes to half an
hour per fit, so the Q2 grid (about 3,000 fits) is days of compute within the 8-worker budget: the grid is cut to the
cells named in H138.2-3 first, the rest descriptive.

## Open choices for the PI
1. The thresholds of H138.2 (20 %) and H138.3 (50 %), and whether the null arm should be calibrated per TCGA cancer
   (stronger, slower) or on one pooled structure.
2. Whether CYCLOPS 2.0 enters the primary hypotheses (no standalone paper found; Li et al., PNAS 2024) or stays
   descriptive.
3. Whether to contact the CoFE authors before publication (outward-facing: draft only).
