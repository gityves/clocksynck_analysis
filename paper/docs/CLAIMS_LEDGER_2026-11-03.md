# Claims ledger (2026-11-03, draft for PI review)

One row per sentence the paper may say. Each row gives the wording that the evidence supports, the evidence (amendment and
number, from `VERDICT_TABLE.md`, `results_preliminary.md`, `FINDINGS_MAP_2026-11-03.md`), the limit the claim rests on, the
wording to avoid, and the objection we expect with its answer. Wording is a draft: the PI decides; outward use (sharing with
the groups whose claims we test) is the PI's call.

## Rules for every sentence

1. **Bounds, not absence.** "No effect larger than X (interval)" or "equivalent within ±X", never "no effect" or "not
   prognostic" alone. A not-estimable test is reported as such, never as a negative.
2. **Scope to what was measured.** Bulk or pseudobulk clock-gene states of organs or cell types; mostly post-mortem
   (GTEx) or surgical material; the cohorts named. Not behavioural rhythms, not the SCN, not cell-to-cell variation within
   a tissue, not transient misalignment within hours of a shift.
3. **Define the term once.** *Internal desynchrony* = a change in the phase offset between compartments of the same
   individual (organs, cell types). Not noise around the host phase (call that loss of clock structure) and not amplitude
   loss (damping).
4. **A positive control next to every null.** The method detects breaks when they exist (row 1.5).
5. **Label provenance.** Registered primary / secondary / post hoc / seen data / archive, as in the verdict table.
6. **Convergence, not refutation,** for published claims: "the association reproduces in direction and is carried by
   proliferation", not "they were wrong".

Status: **A** registered pass, clock-specific where tested; **B** pass with a stated caveat or descriptive; **C** weak,
not estimable or post hoc only (not in the main text).

## 0. Headline (provocative, PI and user's choice)

| | |
|---|---|
| **Sentence** | Across development, ageing, loss of daily cues and chronic disease, human organs keep a shared internal time; where clocks break, they break locally and in a structured way, not as body-wide internal desynchrony. |
| **Rests on** | Sections 2-6 below. Each clause must be traceable to an A or B row. |
| **Avoid** | "Internal desynchrony does not exist"; "circadian disruption is a myth"; "the clock is unaffected by disease". |
| **Main objection** | "You showed robustness where nobody claimed desynchrony; the cases that matter (shift work, jet lag) are absent." **Answer:** state it as a limit. Human misalignment data are not readable (blood, row 7.1); the mouse aorta (row 7.2) is supporting only; the decisive test (multi-organ shift with different baseline phases) is pending (in-house jet-lag data, owner and PI agreement). |
| **Status** | Holds only with the scope sentence of rule 2 and the limit above in the abstract or first paragraph. |

## 1. Internal time is readable from samples taken once, in hospital hours (method)

| # | Sentence | Evidence | Limit | Avoid | Objection → answer | Status |
|---|---|---|---|---|---|---|
| 1.1 | From donors who died between 08 and 18 h, the method recovers the organ phase map that whole-day timing gives. | A123: offsets r 0.945, amplitudes 0.976 (truth shares donors). A123c donor-disjoint truth: 0.679 (registered bar 0.7: FAIL); post hoc split halves 0.873, truth halves 0.544, corrected about 0.84. | The independent-truth bar was missed; the explanation (truth noise) is post hoc. Whole Blood is the outlier (no core-clock signal, post hoc). | "Validated against an independent truth" without the 0.679. | "Circular: truth built from the same donors." → A123c with the miss and the post hoc noise ceiling, both labelled. | B |
| 1.2 | Untimed orderers lose the organ map in working hours; this method keeps it. | A123b: CoFE genome −0.283, CoFE clock −0.417, CYCLOPS −0.296, CHIRAL −0.165 vs 0.945 (Holm 0.0015); with whole-day sampling CHIRAL 0.877, CoFE clock 0.856 recover. Reflection or stretch does not rescue them (post hoc). | The advantage is partial-day coverage, not coupling estimation as such (A145: with full days two-step pipelines do about as well). | "Only our method can infer internal time." | "Unfair settings for competitors." → default settings, whole-day arm shown, A145. | A |
| 1.3 | In living patients, the internal time of breast resections follows the recorded surgery hour, and this ordering is carried by the clock genes. | First analysis (V25.5, 2026-10-15, post hoc): clock-shape gate passed (p_pairing 0.001, p_acrophase 0.002) and single clock genes tracked the hour in their expected directions (normal: ARNTL +0.69, NR1D1 −0.63, NR1D2 −0.67, PER3 −0.71). A200 (registered rerun) adds the matched-set control: Spearman 0.558, above all 200 matched pseudo-clock sets (q95 0.215), Holm 0.006; error 1.42 h without recorded time vs 2.28 h shuffled. Without PER1 0.560. | One cohort of 42, not replicated; leans on NR1D1 (0.233 without PER1 and NR1D1). | "Replicated." | "Surgery hour tracks fasting or anaesthesia, not the clock." → clock-shaped from the first look; the pseudo-clock control and the PER1 removal. | A (one cohort; replication pending) |
| 1.4 | Other living timed material agrees: heart surgery morning vs afternoon; adipose of 12 men; baboon. | A174: displacement p 0.0004, angle 19.73 h (expected 18 h). A113: JS 0.888, Fisher-Lee 0.822, p 0.001 (muscle fails as registered). A36 baboon: offsets error 0.455 h, R 0.987. | Baboon registered statistic degenerate (archive). Muscle fails. | Listing baboon as a registered pass. | — | A (A174, A113), B (baboon) |
| 1.5 | The method detects breaks when they exist (positive controls). | Hepatocyte Bmal1 KO decouples the liver specifically, d* −0.287 [−0.411, −0.118] (A111/A129). Feeding shift: untimed offsets track timed, Fisher-Lee 0.988 (exact p 0.042, the minimum attainable). Brain departure (row 2.2). Psoriasis delay (row 6.1). | Mouse for the genetic and feeding controls. | — | "A null-heavy paper with an insensitive readout." → this row, placed next to every null. | A (KO), B (feeding) |

## 2. Organisation in health: a shared internal time, with the brain apart

| # | Sentence | Evidence | Limit | Avoid | Objection → answer | Status |
|---|---|---|---|---|---|---|
| 2.1 | Organs share a clock state beyond what matched non-clock genes share. | A150: mean c 0.393 vs pseudo-clock 0.242 (q95 0.313), 53 compartments. A144: stronger clocks are more coupled, cross-half rho 0.690. | The pseudo-clock baseline is large (0.242): main panels show Δc or the matched null, never raw c. Coupling intervals under-cover (A154, A162, A163): permutation nulls only. | Raw coupling values as "synchrony". | "Coupling is shared noise or composition." → A150 null, A144 direction. | A (with Δc) |
| 2.2 | The brain runs apart from the body, also in sheep. | A142: brain r* 0.639 vs random 13-tissue blocks 0.978 (q05 0.936), p 0.001; adjusted for quality 0.627; A143 no batch effect. Sheep 0.625 [0.36, 0.87] (archive). | Sheep: baboon readout fails. | "The brain has its own time" (mechanism not shown). | "Post-mortem brain artefact." → quality adjustment, A143. | A (human), B (sheep) |
| 2.3 | Digestive organs also depart from the rest of the body as a group, about a third as much as the brain. | A142 H142.2 (registered): digestive (7 tissues) vs the other 20 peripheral tissues, r* 0.883 [0.819, 0.971] vs random 7-tissue blocks median 0.970 (q05 0.900), Holm 0.028; 1 − r* 0.12 vs 0.36 for the brain. | Holds on the registered reliability-corrected r* only: on the uncorrected r 0.782 vs reference 0.821, p 0.207 (erratum 2026-11-01); without PER1/PER2 p 0.056. Feeding as the cause not tested well (A110 H110.5, MDE 0.20). | "The gut runs on its own clock"; feeding as the explanation. | "Depends on the reliability correction." → state both statistics. | B |
| 2.4 | Ageing: prefrontal clock amplitude falls after 60; we detect no peripheral change. | A3: −0.532 log2 (p 0.022). A19: γ −0.029 per decade, MDE 0.20. A180: age not equivalent (90 % CI [−0.147, 0.089] vs ±0.0375). | Peripheral ageing is "not detected", not "unchanged". | "The peripheral clock does not age." | CHIRAL reports age damping → our archive reanalysis attributes most of it to death type (not re-verified here). | B |

## 3. Development: the clock is organised before birth

| # | Sentence | Evidence | Limit | Avoid | Objection → answer | Status |
|---|---|---|---|---|---|---|
| 3.1 | In the late mouse embryo (E15-E18.75), the clock state alone tells the time of day at which embryos were collected. | A108 H108.1: untimed host phase vs collection clock time R 0.52 (16 embryos, 100 cell types), p 0.016, error 2.02 h; matched pseudo-clock sets carry much less (R 0.191 / 0.334 / 0.194). Earlier windows show no clock-shaped state (A121). | n 16; the model-free clock-shape gate is borderline (p_acrophase 0.08), but the time-of-day inference itself passes and beats matched sets; archive result; onset within the window not shown; organs cannot be traced. | "The clock starts at E15." | "Maternal signals drive expression, not a fetal clock." → open; we show a body-wide time-of-day state readable from clock genes, not its origin. | A/B (time-of-day content A; gate borderline) |
| 3.2 | The human fetal adrenal carries a coherent, clock-specific clock at 13-17 weeks. | A189: R 0.734, Holm 0.0495; within-adrenal gate p_acrophase 0.0018 (A190, narrow, depth confound). | One organ; A190 caveats. | "Human fetal organs are synchronised." | — | A (A189), B (A190) |

## 4. When daily cues are lost, organs lose wall time but keep their relations

| # | Sentence | Evidence | Limit | Avoid | Objection → answer | Status |
|---|---|---|---|---|---|---|
| 4.1 | Without daily cues, organs keep a clock-shaped shared state. | A109 H109.1: C 0.463 vs pseudo-clock 0.154 (Holm 0.025); vs responsive genes 0.214 (q95 0.257). | 0.70 of the shared covariance lies on the non-temporal axis and 0.68 on the glucocorticoid direction: "clock-gene-shaped shared state", not "running clocks". | "Clocks keep running in ventilated donors." | — | A (with the limit) |
| 4.2 | Without daily cues, organs lose wall time but keep the order of their phase relations; the exception is brain death. | Slow deaths (Hardy 4): wall-time R 0.097 (null 0.101); offsets rho 0.784 (H109.2 secondary, p 0.0002). Ventilated donors not declared brain dead keep the order (rho 0.678) while brain-dead organ donors lose it (0.161; difference 0.517, p 0.008; A132b post hoc); the difference is registered in A136 H136.3 (0.458 without PER1/PER2, Holm 0.036). The registered H109.2 on all ventilated donors failed (rho 0.261) because it pooled the two groups. | The ventilated split is post hoc (A132b), its registered confirmation is A136 H136.3; the GR gate of A136 failed, so the mechanism of the brain-death loss is open. Fatal COVID lungs (A109c, descriptive) show no 0 / 12 h push but were too small to show it reliably. | "Ventilation preserves relations" as a registered result; a mechanism for brain death. | "The primary failed." → it pooled brain-dead with other ventilated donors; the split is shown and registered (A136). | A/B |
| 4.3 | Not explained by body temperature. | A110 H110.3: drop in C 0.017 (bound 0.05). | Feeding not tested well. | — | — | B |

## 5. Chronic disease rarely moves the network

| # | Sentence | Evidence | Limit | Avoid | Objection → answer | Status |
|---|---|---|---|---|---|---|
| 5.1 | Type 2 diabetes and obesity leave inter-organ coupling unchanged within ±0.15. | A164 re-read after A202: d* +0.019, 90 % [−0.037, +0.076] (T2D); +0.004 [−0.083, +0.091] (obesity); design-based intervals validated at the A164 design (A202: the A164 T2D bootstrap interval [−0.021, 0.060] was too narrow). | Coupling, metabolic organs, GTEx; phase offsets per organ not separately bounded here. Attenuation of c (A162) not covered by A202; dividing the intervals by 0.8 keeps them inside ±0.15 (descriptive). | "Metabolic disease does not disrupt the clock." | "Contradicts animal and human literature." → bounded estimate at common host spread; scope to inter-organ coupling. | A |
| 5.2 | Where testable, disease atlases show equivalent phase and amplitude. | A180 TOST: phase 7 of 8, amplitude 7 of 10 rows within ±1.5 h / 1.5-fold. | 24 phase and 22 amplitude rows not testable; A165 0 of 89 detected with phase MDE median 2.20 h. | Counting not-estimable rows as negatives. | — | B |

## 6. Where clocks break, they break locally

| # | Sentence | Evidence | Limit | Avoid | Objection → answer | Status |
|---|---|---|---|---|---|---|
| 6.1 | The psoriatic lesion runs late relative to the same person's unaffected skin, replicated in a fresh series. | A54 discovery −1.43 h (0 of 600 pseudo sets reach it, post hoc robustness); A183 fresh: −1.45 h, Holm 0.0002, matched pseudo 5.5th percentile, responsive 0.5th. | Replication sits at the specificity line (5th percentile). Not specific vs eczema (A183 H183.2: matched 20th percentile, NR1D1-carried). | "Psoriasis-specific clock defect." | — | A/B |
| 6.2 | Tumours stay entrained to their own patient's normal tissue. | A138 H138.4: 0.365 vs null 0.122 (q95 0.157), pseudo q95 0.261, responsive q95 0.310; 517 pairs. | — | — | — | A |
| 6.3 | Proliferating and ER-negative adenocarcinomas lose clock coherence: less clock signal in the clock plane (cause not resolved; correction 2026-11-06). | A127, A166, A171, A173: SCAN-B −0.123 / −0.343, METABRIC ER− −0.229, LUAD −0.242, −0.293, gastric −0.187, −0.533. A176: in-plane share falls (7 of 7); disordered relations not tested (H176.2 not informative). Not pan-cancer (A172), not squamous, not HCC. | "Less clock signal" is the supported wording; "swamped" overstated (audit; A212 tests separability). Within ER+ cohorts no proliferation-coherence relation (A196). Tumour-vs-own-normal pairing loss clock-specific in TCGA (A175) but not in CPTAC (A186). Dedifferentiation does not explain it (A197). | "Clock loss drives tumour progression"; "the clock breaks in cancer". | "Composition." → A179, 0 of 5. | A (breast ER−, LUAD, gastric), B (pairing) |
| 6.4 | The loss does not reverse within two weeks of stopping proliferation. | A168 POETIC −0.010 (MDE 0.196); A196 pooled change +0.019 (SD 0.062), about [−0.10, +0.14]; NeoPalAna arrest acted (−0.54 SD) but H196.1 not estimable. | A moderate rise is not excluded; observational treatment series. | "Irreversible tumour state." | — | B |
| 6.5 | Published clock-survival associations in breast cancer reproduce in direction but carry no information beyond proliferation in independent cohorts. | A192 Cadenas genes: score at the 7th percentile of proliferation-matched random sets; 0.88 (p 0.12) with G2M; RORC partial exception (0.80, p 0.019). A194 TimeTeller: HR 0.966 / 0.957 per SD (MDE 1.15). A193 CYCLOPS magnitude: HR 0.920 [0.761, 1.113] top tertile, 0.978 per SD; TCGA reproduction OR 0.81 (p 0.62). A191 our own readout: 0.977, 0.956. A198 no clock x chemotherapy or endocrine interaction beyond proliferation: ratios 1.088 [0.937, 1.264], 1.068 [0.913, 1.249]. | METABRIC and SCAN-B were seen for outcome before some of these tests. A193 used our refit (CYCLOPS seed list, no normal-only eigenspace). Observational treatment. | "Clock genes are not prognostic"; "Li et al. / Cadenas et al. were wrong". | "Your refit is not their pipeline." → list deviations; offer the authors the comparison (PI's call). | A |
| 6.6 | Tumour microenvironment programs follow the recorded surgery hour. | A146 H146.2: S 1.339, Holm 0.0012; adjusted 1.234. | Internal-time version fails; one cohort (GSE233242); fasting confounded with hour. | "Immune escape follows the clock in patients." | — | B (needs a second cohort) |

## 7. Misalignment

| # | Sentence | Evidence | Limit | Avoid | Objection → answer | Status |
|---|---|---|---|---|---|---|
| 7.1 | Human misalignment studies in blood cannot address tissue clocks. | A124 night shift, A156 forced desynchrony: MDE 7.3 h and 13.1 h. Blood time signal follows cell composition (post hoc 2026-11-03: clock genes at the 52nd gene percentile, concordance 0.43 vs 0.85-0.94 in solid tissues). | Composition analysis post hoc. | "Shift work does not desynchronise the body." | — | C (limit statement) |
| 7.2 | Five days after a 6 h advance, aortic cell types sit at the same partial shift. | A201: shift −2.86 h [−3.95, −1.77] (H201.0 pass); SMC −2.85, fibroblast −2.87, EC −2.29 h; D +0.02 h [−0.76, +0.80]. Baseline offset +0.18 h (post hoc). | H201.1 not estimable (1.645 SD 0.78 > 0.75 h); pools of 2 mice, n = 1 per sex x ZT; cell types already in phase at baseline, so preservation of a non-zero offset is not shown; no SCN. | "Jet lag does not cause internal desynchrony." | The authors' central-vs-peripheral misalignment is real and untested here. | B (supporting) |

## 8. Method comparisons with published tools (supplement)

| # | Sentence | Evidence | Limit | Avoid | Status |
|---|---|---|---|---|---|
| 8.1 | In GTEx tissues with recorded time of death, CYCLOPS magnitude is nearly unrelated to a timed measure of clock strength, with or without whole-day sampling. | A199: median Spearman 0.055 (whole day, Holm 0.043) and 0.020 (08-18 h); coverage effect not detected (0.075, p 0.26); truth split-half 0.61 (post hoc). | 58-98 samples per fit; fits unstable across seeds (0.43), unlike tumour cohorts (≥ 0.976). Says nothing direct about magnitude in 500-1,600-sample cohorts. | "CYCLOPS cannot be used in hospital hours"; "CYCLOPS magnitude is meaningless". | B |
| 8.2 | CYCLOPS 2.0 can order samples along a shared non-clock factor. | A140 simulator gate: median circular r 0.060 at full factor; fails to start without it. | Simulator; our port with fewer models in one arm. | — | B |

## Do not claim (not supported or not estimable)

- Any effect of night shift or forced desynchrony on tissue clocks (row 7.1).
- Immune cue-following (A85) as cue-specific (A122 failed).
- Bootstrap confidence intervals for coupling (coverage fails).
- Within-tumour single-cell contrasts (A187, A188, A195 not estimable).
- Dedifferentiation as the cause of tumour coherence loss (A197 failed: SCAN-B +0.092, METABRIC −0.067; MDE 0.11-0.13).
- A clock x treatment effect on breast cancer outcome (A198 failed) or any chronotherapy statement (dosing time not recorded).
- Fibrotic-lung T-cell decoupling beyond "conditional on procurement" (archive).

## To close before submission (affects the ledger)

1. A fresh cohort with recorded resection times (row 1.3 to status A; row 6.6).
2. A multi-organ shift with different baseline phases (headline limit; rows 7.1-7.2): in-house jet-lag data with the
   owner's and PI's agreement, or published multi-organ shift data if found.
3. Main-panel coupling readouts as Δc or with the matched null (row 2.1).
4. State which tests were on fresh cohorts and which on cohorts seen before (rule 5; METABRIC, SCAN-B, GTEx).
5. Embryo (row 3.1): n 16 with a borderline clock-shape gate, although the time-of-day inference passes; a second
   developmental dataset would make "before birth" robust.
6. PI decisions: wording of rows 6.5 and 8.1, and whether to contact the authors before submission.

## Addendum 2026-11-06 (A212-A218; PI review pending; rows above unchanged)

| # | Sentence | Evidence | Limit | Avoid | Status |
|---|---|---|---|---|---|
| 6.3a | The proliferation loss is less in-plane signal without added non-clock variance (swamping excluded). | A212 V212.2 AUC 1.000; observed SCAN-B point outside the swamping cloud. | Damping, dephasing and rewiring not separable in bulk (V212.1 AUC 0.89-0.91) nor in public single cell (A214 scout). | "Swamped"; "rewired". | A |
| 6.7 | Genomic lesions (TP53, copy-number burden, MYC amplification, RB pathway, subclonality) do not lower tumour clock coherence beyond molecular subtype and proliferation. | A217 at matched subtype: TP53 +0.035, FGA −0.015, MYC +0.033 (MDE 0.13, 0.08, 0.17); TCGA BRCA post hoc: losses shrink 59-96 % within PAM50; A215 H215.3 RB −0.022 (MDE 0.068); A213 H213.1 subclonal +0.021 (MDE 0.058); METABRIC (PAM50 strata) no MYC or FGA loss. | Bounds are 0.06-0.17, not 0; 11 TCGA projects have no subtype labels and are unchecked. | "Genomic instability disrupts the clock"; citing A210 H210.2, A213 H213.2, A215 H215.1 as lesion effects. | A (as a bound) |
| 6.9 | In breast cancer, luminal A tumours keep more clock signal than other subtypes at equal proliferation. | A219 H219.2 +0.276 (METABRIC + TCGA BRCA, above all pseudo sets), H219.3 +0.222 (SCAN-B, 97.5th pct); with ER− within PAM50 (A166 H166.2). | Breast only: no subtype effect elsewhere (H219.1, bound 0.025); SCAN-B ordering seen before; ER and PAM50 collinear; causality not shown. | "Estrogen signalling maintains the tumour clock." | A (breast) |
| 6.8 | Physiological proliferation leaves the clock coherent. | A211 regenerating liver T −0.027 (MDE 0.064, powered fail); A126 (underpowered). | One mouse dataset; G2M rise not on the tumour scale. | "Proliferation never affects the clock." | B |

**Do not claim (added):** per-cell rewiring or intact gene relations in malignant cells (A214 scout: not estimable;
WP0 item 2 erratum); checkpoint loss as the link between proliferation and clock loss (A215 H215.4 opposite sign); subtype
as the carrier of the proliferation loss (A216, A218 not scored); a general lineage effect across cancer types (A219 H219.1).

## Addendum 2026-11-12 (cancer synthesis after A224, A227-A229, A245-A250; PI review pending; rows above unchanged)

**Tally of the cancer rows of `VERDICT_TABLE.md` (about 80, by verdict text).** 29 null at the stated MDE; 16 pass; 16
not estimable (power rule or failed precondition); 14 detected but not clock-specific. Of the 16 passes, those still
standing are of two kinds: the tumour sample shares the patient's clock state (H138.4, H151.1, H227.1, H228.4, H200.1;
H229.5 at the line) and clock-gene levels differ in tumours (H207.1, H221.1). Passes that were contrasts between tumour
groups (H127.1, H166.2, H175.1, H210.2, H215.1, H219.2) later failed a robustness check or lost clock specificity (A178,
A185, A186, A217, A224); H213.2 has no follow-up.

| # | Sentence | Evidence | Limit | Avoid | Status |
|---|---|---|---|---|---|
| 6.2a | Within the tumour, myeloid cells and CD8-negative T cells share the clock state of the patient's normal tissue, as fibroblasts do. | A227 H227.1: 0.581 vs null 0.361, pseudo q95 0.548 (169 patients, 24 cohorts); A228 H228.4: 0.532 vs 0.322, pseudo q95 0.501 (214, 26); A151 H151.1 0.860. | CD8 T, NK and the cross-compartment pairings are detected but not clock-specific (A228 H228.1-3, H228.5). Malignant cells: coupling 0.076 (H138.1), pairing not clock-specific (A149 H149.1). | "Malignant cells are entrained"; "malignant cells are uncoupled". | A |
| 6.10 | No tumour program tested follows the tumour's or the patient's clock state: microenvironment programs, endothelial adhesion genes, drug-handling genes, drug targets, cell-cycle and repair genes, MGMT. | A246 H246.1 T −0.121 (MDE 0.85); H246.2 not clock-specific; A229 H229.1 not clock-specific, H229.2 fail; A249 set C S_T 0.117, F 0.122 (null 0.120; MDE 0.018, 0.011); A250 set C 0.097, F 0.097 (null 0.086, p_clock 0.57, 0.58); A247 H247.3-4 (MDE 0.13-0.14). Exceptions: T/NK share (A229 H229.5, borderline); NFIL3 in two cohorts (A249, A250 H250.4: 0.155 vs 0.086, Holm 0.002, implied peak 18.9 and 19.0 h). | A bound at each MDE, on untimed surgical samples. The positive controls of the same designs fail (A249 / A250 precondition: clock-controlled genes of normal tissue not read from the tumour state; A248 blood T cells in GTEx), so these are "not readable", not "not rhythmic". | "Tumour programs are not under clock control"; "there is no chronotherapy window". | B (as a limit statement) |
| 6.11 | Across a patient's two samples, clock-controlled genes of the tumour follow the clock state of the normal tissue in TCGA; this is not clock-specific in CPTAC. | A249 descriptive: set O S_T 0.162 vs 0.120, pseudo q95 0.147, matched q95 0.136 (517 pairs); A250 H250.1: 0.109 vs 0.086, p 0.006, pseudo median 0.107, q95 0.131 (686 pairs). | Descriptive in TCGA, not replicated as clock-specific: not for the main text. | Any wording of an asymmetry between clock output and drug handling (A250 H250.2: Delta 0.012, p 0.099). | C |
| 6.6a | Row 6.6 stands on recorded time only: the hour of surgery is seen in normal breast tissue more than in the tumour, for every gene set tested. | A249 H249.5: set C modulus 2.416 normal vs 1.458 tumour, p 0.030; K 0.065, O 0.082, F 0.047 (42 patients). | One cohort (GSE233242), already used in A146 and A200; no matched-gene companion; not set-specific. | "Drug-handling genes lose their rhythm in tumours." | C |
| 6.9a | Correction to row 6.9: on cis-copy-number-adjusted expression the luminal A excess is 0.101 against 0.276 and is not clock-specific. | A224 H224.1: attenuation 63 %, pseudo 92nd percentile; H224.2 TP53 −0.052 (p 0.097). | See the verdict table for the status of A224. | Citing 6.9 without this row. | B → to be decided by the PI |

**Why the negatives (reading, not a result).** (1) Between patients, time of day is a small axis in tumours: purity,
stroma, proliferation and copy number move clock genes and thousands of others together, so pairings pass against
permutation and fail against pseudo-clock sets (the 14 rows above). Surgery hours compress body time further. (2) The
reader of body time is too noisy for genes downstream of the clock: tumour and normal states pair at 0.365 (H138.4) and
clock-controlled genes share about 0.22-0.24 with the state of their own sample (same-sample readouts of A249 / A250,
confounded upward), which puts their expected cross-sample excess near the MDE of A249 (0.025-0.027); weaker rhythms are
below it by design. This back-of-envelope figure is not a registered number. (3) Where the positive control fails, a
negative is uninformative (rule 1). (4) Copy number and proliferation confound contrasts between tumour groups (A224,
A217, A186). (5) Three published human claims do not hold beyond proliferation (rows 6.5; A194, A220). (6) The premises
are from mouse designs with 24 h sampling.

**Do not claim (added):** a chronotherapy window or its absence (A249, A250: the normal half is unreadable in two
cohorts); a clock in glioma cells tied to the patient's glia (A247: precondition fails); adhesion-gene or T-cell rhythms
in healthy human tissue at the hours implied by A229 (A245: fail; H245.2 uninformative after A248); any statement that
tumour-microenvironment programs follow the internal clock (A246).

**What would change this:** recorded sampling times. The only hour-of-day signals in tumours come from the one timed
cohort (A146 H146.2, A200, A249 H249.5). The 108 REMAGUS biopsy times (Vasilyev et al. 2026; files
`Remagus_Actual_Time.csv` and `GSM_to_number.xlsx`, not public) would give a second one; asking for them is outward-facing
(user's decision).

## Addendum 2026-11-13 (A251-A255: living cohorts, within-person validation, lesions; PI review pending; rows above unchanged)

| # | Sentence | Evidence | Limit | Avoid | Status |
|---|---|---|---|---|---|
| 1.6 | Within a living person the readout advances with the day: between two biopsies of one day the clock state moves along the circle in the predicted direction. | A252 H252.1 (adipose, exact clock times, 15 men x 2 visits): rho 1.265, null SD 0.215, above all 200 matched and 200 responsive sets; H252.2 (seven studies, muscle and adipose, assumed morning times): 2.970; without PER1 and NR1D1 1.058 and 3.566; implied time of day within 1 h of the protocol's at 13 of 14 points. | Every interval but one small arm holds a meal or an intervention (A256 registered to test "nothing done" arms and to replicate); seven studies use protocol clock times; rho mixes amplitude and rotation and exceeds 1. | "The method measures elapsed time"; any hour estimate from rho. | A |
| 2.5 | In living people sampled once in clinic hours, tissues share one clock state: skin layers at three body sites, fat and skin, two fat depots, muscle and fat. | A251 H251.1 (154 women, same layer at different sites): r 0.628 vs null 0.072, matched q95 0.380, responsive q95 0.420; H251.2 0.813 / 0.594. A255 H255.1 (MuTHER, 680 women, fat-skin): 0.533 vs 0.034, matched q95 0.304, responsive q95 0.325; H255.3 (two surgical cohorts, 83 persons): 0.686, responsive q95 0.530. | Women only in the two large cohorts; arrays; fat and skin of MuTHER come from one punch biopsy; ages not public (no age contrast); time spread unknown. | "Organs of living people are synchronised" without the scope (skin, fat, muscle). | A |
| 2.6 | The shared state is body time, not the person's constitutive clock-gene expression: her cultured cell line does not carry it. | A255 H255.2: r(fat, skin) minus mean r(tissue, cell line) 0.479 (0.432 to 0.510); fat-cell line 0.042 (p 0.24), skin-cell line 0.066 (p 0.055); the cell line fails the clock gate (p 0.75). Co-twins pair across persons alike in MZ and DZ pairs (fat 0.522 / 0.520; skin 0.587 / 0.591), above their cell lines by 0.135 and 0.205. | Without PER1 and NR1D1 the skin-cell line pairing is 0.248 (descriptive, no companion): the control is not flat under every gene subset. Co-twin cell lines pair at 0.384 (shared handling). | "Clock phase is not heritable" (not tested: this is state sharing at a visit). | A |
| 1.7 | An untimed offset agrees with a timed one in the same study: dermis runs behind epidermis by 1.6 h (untimed women) and 2.0 h (timed men). | A251 H251.4: -1.58 h (-1.81 to -1.34) vs -1.98 h (-3.15 to -0.95). | One layer pair, one study; the men's control was seen in a code test before the power note (disclosed). | — | A/B |
| 2.7 | Skin sites differ by at most about 1.5 h; the cheek runs latest. | A251 H251.3: offsets -1.54 to +0.73 h; registered equivalence within 1.5 h fails for one of six pairs. | Each compartment is one array batch; exposure, site and composition not separable. | "Sun exposure delays the skin clock." | C (descriptive) |
| 6.1a | The psoriatic delay is not a general property of chronic lesions: lesions of gum, cartilage, artery and uterus stay on the time of the person's own tissue. | A254 H254.1: pooled offset +0.01 h, MDE 0.69 h (184 pairs, all gated); per disease -0.15, +0.04, +0.53 h; hidradenitis (skin) -1.07 h (-1.79 to -0.14), Holm 0.34. | Bound 0.7 h pooled, 0.8 to 2.2 h per disease; arrays; "non-lesional" is tissue of a diseased organ. | "Inflammation delays the clock"; citing hidradenitis as a replication. | A (as a bound) |
| 6.12 | Lesional tissue carries less clock signal than the person's own unaffected tissue. | A254 H254.3: pooled log variance ratio -0.417 (SE 0.088, p 0.0001), below the q05 of matched (-0.404) and responsive (-0.272) sets; periodontitis -0.453 (replicated in a second series, -0.481), plaque -0.545. | Clock-specific only when pooled and at the line; composition not separated; no single disease is clock-specific. | "Lesions lose their clock." | B |

**Do not claim (added):** altered clock relations between brain regions in depression or psychosis, or their
preservation (A253: not scored; in control brains stress-responsive genes pair across regions as strongly as clock
genes, R 0.734 vs responsive q95 0.758; the depression patients of Li 2013 are not public); an age effect in the
living (ages of the skin and twin cohorts are not public); a graft adopting the recipient's time (no usable public
data, second search); a metabolic contrast in living multi-tissue cohorts (no cohort with phenotype and person key).

**Effect on the headline (section 0) and on "Missing for Cell" (REVIEW_2026-11-10 §3.1).** Item 2 (a living cohort
for the robustness arm) is met for the shared-state claim by two large living cohorts and two surgical ones, with a
genetic control; not for the equivalence claims (sex, age, disease), which remain GTEx-only. Fig. 1 gains a
within-person validation (row 1.6). Item 1 (a confirmation of the break arm or a human desynchrony positive) is
unchanged: the brain route is closed (A253).
