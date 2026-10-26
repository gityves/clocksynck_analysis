# A140 (DRAFT for PI review, 2026-10-26; not registered, nothing scored): does a CYCLOPS 2.0 donor internal phase across immune compartments report a clock that is not there?

Complements A138 Q2 (bulk tumours, CoFE and the CYCLOPS 2.0 healthy projection). Here: the multi-compartment single-cell
design of the bioRxiv preprint "Circadian misalignment underlies immune escape in breast cancer" (doi
10.64898/2026.05.26.726543, posted 2026-05-29; PDF and supplement in the archive `literature/circadian_breast_cancer{,_supp}.pdf`).
Written after reading its Figure 7 and Figure S8 and after A17 (cancer) was scored; nothing below has been computed.

## What the preprint did (methods, main text lines 275-299 of the methods section)
Human TNBC scRNA-seq (Zhang et al., Cancer Cell 2025; 76 samples). Pseudobulks of CD4+ T, CD8+ T and myeloid cells per
sample, ordered jointly by CYCLOPS 2.0 with cell type as a covariate. The seed is the union of "strongly correlated circadian
genes" of the same three compartments in their timed mouse tumours, mapped to human orthologs (Table S5, not in the PDF;
its size and whether it contains core clock genes are unknown). Per sample, the donor internal phase (DIP) is the circular mean of the three
phases, after iteratively removing cell types more than pi/2 from the mean (up to 5 iterations); a 12 h shift maps mouse to
human time. Samples are binned into six DIP groups; DiscoRhythm on the binned pseudobulks. No healthy reference, projection,
GTEx or timed human sample is mentioned (Li et al. 2024 state their projection explicitly): the reading here is a free fit
on the 228 tumour pseudobulks, the cell-type covariate's offsets learned from the cohort itself. Data: GEO GSE266919
(whether the 76 samples are 76 patients or repeated biopsies of fewer patients is to be checked).
Internal consistency: the DIP is read as clock time (daytime immune activation, link to morning immunotherapy), which needs
immune clocks to follow the host; then DIPs should cluster in surgical hours. If immune clocks in tumours were instead
desynchronised from the host, a spread DIP would be possible but would no longer be clock time.

## Claims, and why each can arise without a clock (lead's reading of the figures)
1. **Cross-compartment concordance** (S8B: circular r 0.81 myeloid-CD8, 0.80 myeloid-CD4, 0.92 CD8-CD4). The three
   pseudobulks of a sample share its non-clock variation (patient immune state, composition, depth, batch) and are ordered
   in one model; the outlier rule pulls the phases together by construction.
2. **Sample phases** (S8C): DIPs spread over about three quarters of the cycle, with a gap; surgical samples are expected
   in a 10-12 h window. A centred, partial-day arc is spread around the circle (A138 note 2 (iii)); archive: CYCLOPS
   ordering skill <= 0.11 on 08-18 h GTEx draws (`figure_1/scripts/method/README.md`, M3).
3. **Gene acrophases** (S8D): the rhythmic genes of all three compartments peak in one narrow sector (about CT6-CT10), not
   spread over the cycle as clock outputs are (GTEx followers). One dominant non-clock axis turned into a circle gives this:
   every gene loading on it peaks at the same end. Thousands of rhythmic genes (Fig 7C, F, I: 5,083 / 4,566 / 3,014) and the
   composition rhythm (Fig 7B: myeloid vs CD8+ T fractions) follow from the same axis.
4. **Clock genes in order** (S8E, S8F): ARNTL, NR1D2, DBP, PER ordered as in primates, ARNTL about 1.6-fold in T cells.
   If they are in the seed, the ordering is fitted to them (circular); archive: a correct clock-gene arrangement does not
   validate the ordering (breast case: GTEx breast ordering skill 0.07 with clock-gene peak times within 1.8 h).

Our own data (A17, cancer, learned anchors, gate-passing cohorts, scored): tumour T/NK coupling to the host is near zero in
the two largest cohorts (esophageal SRP327758, n 50: 0.026 [0, 0.257]; ovarian SRP329496, n 39: 0.005 [0, 0.198]); tumour
myeloid 0.198 [0.014, 0.422] and 0.229 [0.040, 0.450]. These are couplings to the stromal host, not T-myeloid concordance,
and not breast; they motivate, they do not test, the claims above.
<!-- src: paper/results/figure_5/output/cancer/a17_cancer.csv (rule learned, cls T_NK / myeloid, site tumour) -->

## Q1. Simulation of the design
Simulator: `clocksync.simulate`-based, 76 patients x 3 compartments (CD4-like, CD8-like, myeloid-like) at single-cell
pseudobulk depth; p = 2,000 genes: clock genes at GTEx waveforms, rhythmic followers, and a non-clock patient factor
(K = 1-3 shared axes with loadings on 10-30 % of genes, variance calibrated to the leading non-clock principal components of
our gate-passing tumour immune pseudobulks). Sampling 08-18 h (and whole day as a control). Arms:
- **A (clock, coupled):** one patient clock shared by the three compartments (coupling 0.8), plus the non-clock factor;
- **B (no clock):** clock genes at their mesor, non-clock factor only;
- **C (independent clocks):** each compartment its own clock phase (uniform over patients), plus the factor;
- **D (whole day):** arm A sampled over 24 h (positive control).
Seed grid (Table S5 unknown): 50, 150, 500, 1,500 genes, drawn from the genes most correlated with the factor and the clock
in a simulated timed "mouse" reference, with and without the core clock genes. 50 replicates per cell.
Pipeline as published: CYCLOPS 2.0 with the cell-type covariate, DIP with the pi/2 rule, six DIP bins, rhythm detection
(cosinor on bin means, as DiscoRhythm's default; p < 0.05). Same data through `clocksync` (free host on the three
compartments, gate; no anchors exist in this design).
Simulator gate (before scoring): in arm D, CYCLOPS recovers the true phases (circular correlation >= 0.7 at the largest seed);
if not, Q1 stops and the pipeline is fixed first.

- **H140.1 (false concordance):** arm B, seed 500: median cross-compartment circular r of the CYCLOPS phases >= 0.5 in more
  than 50 % of replicates (one-sided exact binomial against 0.5, p < 0.05: pass at >= 32 of 50).
- **H140.2 (false rhythms):** arm B, seed 500: the DIP pipeline calls >= 500 rhythmic genes per compartment in more than 50 %
  of replicates (same test).
- **H140.3 (specificity of clocksync):** arm B and arm C, seed irrelevant: `clocksync` gate fails in >= 90 % of replicates
  (false gate pass <= 10 %), and in arm A at 08-18 h the gate passes in >= 80 %.
- **Descriptive:** DIP coverage (share of the circle, Rayleigh R) per arm; acrophase concentration of the called genes
  (Rayleigh R) per arm; clock-gene order (A138 `GTEX_MATCH`) per arm, with and without clock genes in the seed; all over
  the seed grid.

## Q2. Real data (descriptive unless the PI wants a hypothesis)
On our gate-passing tumour cohorts: the published pipeline (CYCLOPS 2.0 on T/NK and myeloid pseudobulks, cell-type covariate,
DIP), its readouts 1-4, and `clocksync` on the same pseudobulks (free host on T/NK + myeloid with the gate; anchored coupling as
A17). Breast/TNBC data: not in our inventory; Zhang et al. 2025 to be checked for public access first.

## Open choices for the PI
1. Thresholds of H140.1-H140.3, and whether the seed-500 cell is the right primary.
2. Whether Q2 stays descriptive.
3. Whether to ask the authors for Table S5 (outward-facing: draft only).
