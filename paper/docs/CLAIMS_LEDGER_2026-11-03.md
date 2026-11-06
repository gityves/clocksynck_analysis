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
