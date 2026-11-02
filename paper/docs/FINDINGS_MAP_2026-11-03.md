# Findings map against the Cell criterion (2026-11-03, draft for PI review)

Cell criterion (official): results that provide significant conceptual advances into, or raise provocative questions and
hypotheses regarding, an interesting and important biological question. Sections follow the figure plan in
`provisory_plan_paper.md`. Numbers are copied from `VERDICT_TABLE.md` and `results_preliminary.md` (output tables); archive
rows (A1-A108c) are as flagged there, not rescored. "P" = registered primary; "post hoc" as labelled in the docs.
Status: strong = registered pass (clock-specific where tested); qualified = pass with a stated caveat; weak / open = see note.

## 1. Clock relations and internal time are readable from samples taken once, in hospital hours (what only this method does)

| finding | numbers | source | status |
|---|---|---|---|
| GTEx organ map from 08-18 h deaths only (40 tissues, 111 donors) | offsets r 0.945, amplitudes 0.976, donor order 0.613; CoFE / CYCLOPS offsets −0.28 to −0.42 | A123 H123.1-4 | strong (truth (i) shares donors) |
| Same map vs a donor-disjoint truth | offsets 0.679 (registered bar 0.7: FAIL), amplitudes 0.965; post hoc: half maps agree 0.873, truth (iii) halves 0.544, corrected about 0.84 | A123c, post hoc 2026-11-02 | qualified (shortfall is truth noise; blood explained below) |
| Whole Blood is the outlier because blood has no core-clock time signal | clock genes at the 52nd gene percentile (solid 99.5-99.96th); phase concordance 0.43 vs 0.85-0.94 | post hoc 2026-11-03 | descriptive |
| Heart surgery, morning vs afternoon (living patients) | gate pairing 0.0005; PM − AM displacement p 0.0004, angle 19.73 h (expected 18 h) | A174 V174.1-3 | strong (confirmation) |
| Breast surgery: internal time vs recorded resection time | Spearman 0.56, median error 1.1 h | V25.5 (archive, post hoc note) | weak (needs a registered rerun) |
| Living timed adipose (12 men) | JS 0.888, Fisher-Lee 0.822, p 0.001 | A113 | strong; muscle fails as registered (post hoc FL 0.576) |
| Baboon timed truth | offsets median error 0.455 h, R 0.987 | A36 (archive) | qualified (registered JS statistic degenerate) |
| Regenerating mouse colon, blind | epithelial shift −8.95 h untimed vs −9.03 h timed | A21 R21.4 (archive) | qualified (anchor rule V17.2 failed) |
| Genetic positive control: hepatocyte Bmal1 KO decouples the liver, specifically | d* −0.287 [−0.411, −0.118]; WAT −0.009, lung +0.019 | A111 / A129 H111.3 | strong |
| Feeding shift: untimed offset change tracks the timed change | Fisher-Lee 0.988 (4 tissues, exact p 0.042) | A111 H111.1 Xin | qualified (minimum attainable p) |
| Mouse embryo: early clock structure lets clock time be inferred | **source to identify** (not found in RP; A108 / A121 cover the late window) | ? | to locate |

## 2. Organisation in health: a shared internal time, with the brain (and partly the digestive tract) departing

| finding | numbers | source | status |
|---|---|---|---|
| Brain departs from the periphery, not by artefact | brain r* 0.639 vs random 13-tissue blocks 0.978 (q05 0.936), p 0.001; adjusted for quality 0.627; no batch effect | A142 H142.1, A143 | strong |
| Conserved in sheep | sheep r* 0.625 [0.36, 0.87] | A42 / A103 (archive) | qualified (baboon readout 1 fails) |
| Digestive departure | r* 0.883 vs 0.970, p 0.028 | A142 H142.2 | qualified (not on the bounded r) |
| Coupling is clock-specific across atlases | mean c 0.393 vs pseudo-clock mean 0.242, q95 0.313 (53 compartments) | A150 H150.1 | strong, but the pseudo baseline is large (Δc readout needed) |
| Stronger clocks are more coupled, not shared noise | cross-half rho 0.690 | A144 H144.1 (A104 post hoc direction) | strong |
| Fetal adrenal carries a coherent, clock-specific clock at 13-17 weeks; shared across its classes | R 0.734, Holm 0.0495; within-adrenal gate p_acrophase 0.0018 | A189, A190 | strong / qualified (A190 narrow, depth confound) |
| Mouse embryo: shared clock late in gestation | H108.1 R 0.52, p 0.016; onset not shown (A121) | A108, A121 | weak |
| Ageing: prefrontal amplitude loss above 60; no peripheral age effect | −0.532 log2 (p 0.022); A19 γ −0.029/decade, MDE 0.20 | A3, A94, A19 | qualified (age not equivalent, A180) |

## 3. When daily cues are lost, organs lose wall time but keep their relations

| finding | numbers | source | status |
|---|---|---|---|
| Inter-organ coherence above clock-specific nulls | C 0.463 vs pseudo 0.154 (excess 0.310); vs responsive genes 0.214 | A109 H109.1, -R | strong |
| Ventilated / slow deaths: wall time lost, offsets kept | Hardy 4 wall-time R 0.097 (null 0.101); offsets rho 0.784, residual 0.750 | A109c, H109.2 sec | strong (H109.2 primary fails on ventilated donors) |
| Brain death: preserved order, not through PER1/PER2 | H136.3 0.458, Holm 0.036 | A136 | qualified (GR gate fails) |
| Not temperature; feeding not tested well | drop in C 0.017 (bound 0.05) | A110 H110.3, H110.5 | exclusion |

## 4. Chronic conditions rarely move the network (quantified nulls)

| finding | numbers | source | status |
|---|---|---|---|
| T2D and obesity equivalent to no coupling change (metabolic organs, matched GTEx) | d* +0.019 [−0.021, 0.060]; +0.004 [−0.081, 0.089]; bound ±0.15 | A164 | strong (provocative vs the literature) |
| Disease phase / amplitude equivalence where testable | phase 7 of 8, amplitude 7 of 10 rows equivalent | A180 (A165) | qualified (most rows not testable) |
| Liver fibrosis, other disease atlases | 0 of 3 / 0 of 89 detected; MDEs 0.29-0.43 | A167, A165, A153 | not estimable at d* 0.2 |

## 5. Where clocks break, they break locally

| finding | numbers | source | status |
|---|---|---|---|
| Psoriasis lesion delay, replicated in a fresh series | −1.45 h (Holm 0.0002); matched 5.5th, responsive 0.5th pct; vs eczema −1.09 h | A54, A183 | qualified (not clock-specific by the matched rule; NR1D1-dependent) |
| Tumour clocks follow their own patient's normal tissue | 0.365 vs null 0.122 (q95 0.157), 517 pairs; pseudo q95 0.261 | A138 H138.4 | strong |
| Tumour fibroblasts pair with normal | 0.860 vs 0.504 (98 patients) | A151 | qualified (specificity uninterpretable) |
| Proliferating tumours follow their normal tissue less | TCGA D +0.158, clock-specific; CPTAC D +0.113, not clock-specific | A175, A186 | qualified |
| Coherence loss with proliferation and ER−, by histology | SCAN-B −0.123 / −0.343; METABRIC ER− −0.229; LUAD −0.242, −0.293; gastric −0.187, −0.533; not pan-cancer (A172), not squamous, not HCC | A127, A166, A171, A173, A172 | strong (breast ER−, LUAD, gastric); A178 robustness partial |
| The loss is a swamped clock plane, not disorder; not cell-cycle composition | 7 of 7 Δs below pseudo q05; 0 of 5 composition | A176, A179 | strong / resolved negative |
| Stable: does not recover when proliferation or ER signalling stops | POETIC −0.010 (MDE 0.196); 4 series pooled +0.019 (SD 0.062) | A168, A169, A196 | qualified (moderate rise not excluded) |
| No prognostic value beyond proliferation; published claims reduce to proliferation | A191 HR 0.977, 0.956; A192 score at the 7th pct of proliferation-matched sets; A194 HR 0.966 / 0.957 | A191-A194 (A193 running) | strong negatives |
| Tumour microenvironment programs follow recorded surgery time | S 1.339, Holm 0.0012; adjusted 1.234 | A146 H146.2, H146.5 | strong (internal-time version fails) |
| Fibrotic-lung T cells decouple | HLCA 0.14 → 0.00; Natri 0.45 → 0.00 | A5 / A17 / A98 (archive) | conditional on procurement |

## Not supported or not estimable (to state, not hide)

- Perturbations in blood (night shift A124, forced desynchrony A156): MDE 7.3 h and 13.1 h; blood has no core-clock
  signal (post hoc 2026-11-03).
- Immune cue-following (A85) is not cue-specific (A122).
- Coupling interval coverage fails (A154, A162, A163): report point estimates with permutation nulls, not bootstrap CIs.
- Single-cell tumour atlases do not power within-tumour contrasts (A187, A188, A195).

## Candidate conceptual claims for Cell

1. **Human tissues keep internal time with each other, not with the wall clock, and it can be read in hospital hours.**
   Sections 1 and 3: hospital-hours map, living surgery cohorts, lost cues with relations kept.
2. **The brain (and partly the gut) runs apart from the body**, conserved in sheep (section 2).
3. **Chronic metabolic disease leaves the network intact; breaks are local** (sections 4-5: equivalence vs psoriasis and
   tumours).
4. **Tumours stay entrained to their host, but proliferating adenocarcinomas lose clock coherence stably and without
   prognostic cost** (section 5). Provocative question: what does clock loss do in tumours?
5. **The clock is organised before birth** (fetal adrenal; embryo pending its source).

Main reviewer risk: the pseudo-clock baseline of coupling (A150: 0.242 vs 0.393), so main panels use the Δc readout or
show the matched null.
