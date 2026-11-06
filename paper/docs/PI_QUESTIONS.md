# Questions for the PI (open first, answered below)

Each entry: date, question, what is on hold until it is answered, the evidence (with `src` pointers in
`results_preliminary.md`). Work that does not depend on an answer continues. When answered, move the entry to
"Answered" with the answer and its date.

## Open

### Still with the PI (updated 2026-11-06)
- **Q24** (below): the cancer chapter after A214-A226: structure and claims for the paper.
- **Q12** (below): brain time-of-death request. The PI reviews and sends the draft themselves; nothing is sent by the agent.
- **A119:** lung T-cell procurement scout (outward-facing / procurement).
- **Outcome analyses (C8):** scored as A191 on the user's instruction (2026-11-02): both primaries fail (HR 0.98 and 0.96 per SD); draft in
  `paper/exploratory/C8_internal_time_outcomes/DRAFT_C8.md`; results in `results_preliminary.md` (A191). PI to review.
- **Replication of published clock-survival claims (C8b, A192-A194; user's instruction 2026-11-02):**
  - Cadenas genes: carried by proliferation, with RORC as a partial exception.
  - TimeTeller Theta: not prognostic.
  - CYCLOPS magnitude (Li et al. 2024): not prognostic (HR 0.92 top tertile), not reproduced in TCGA with our refit.
  - Clock x treatment (A198): no interaction with chemotherapy or endocrine therapy beyond proliferation (MDE 1.19 / 1.27).
  - Outward-facing use (sharing with these groups, wording in the paper) is the PI's call.
- **Phase 8** (naming, scorecard, figure regeneration, hardening, tutorials, METHOD/CHANGELOG update): PI.
- **Closed 2026-11-01:** Q21, Q22, A136 (registered). **Closed 2026-10-31:** Q8, Q13, Q16, Q17, Q18, Q19, Q20, and the
  A121 and A135 placements.

### Q24 (2026-11-06): the cancer chapter after A214-A226 — structure and claims for the paper
The user wants **one paper** that includes the cancer chapter, held together by "how the body keeps its time and when it
loses it". This needs your review (Phase 8: framing and figures). All numbers: `results_preliminary.md` (sections named).
- **Positive, model-based:**
  - The tumour clock state matches the same patient's normal tissue (H138.4: 0.365 vs pseudo q95 0.261; 517 pairs, 7 of 10
    cancers p < 0.01).
  - The match does not depend on purity (H138.8: low − high +0.025, MDE 0.16), so malignant cells contribute. The shared
    state partly reflects the shared surgery hour, so the claim is that both tissues keep body time, not that a coupling
    strength is measured.
  - High-proliferation tumours pair less (A175 H175.1, +0.158, clock-specific).
  - Tumour fibroblasts and endothelium follow the patient's normal stroma (A149 note 5, post hoc, about 100 patients each).
    A confirmatory test in independent cohorts is registered (A226, running).
  - HCC is less coherent than its own adjacent liver (GSE14520, −0.33).
- **Corrective (the method's specificity and design discipline applied to the literature):**
  - Clock "prognosis" is proliferation (A191-A198).
  - Circadian LASSO signatures are no better than random gene lists; they "validate" exactly where random lists do (A220;
    KIRC, LGG, METABRIC 100 % of random lists).
  - Clock-gene mutations and losses occur at background rates; PER1 / PER3 losses ride on 17p (TP53, 98.9 % co-loss) and
    1p36 (A222).
  - "Clock genes are downregulated" is partly real in paired samples (PER1-3, CRY2, TEF, RORA; L −0.33) but not for ARNTL,
    NPAS2 or NR1D1, which GTEx comparators turn into apparent losses, partly through sampling time (A221).
  - Genomic-lesion effects on coherence are subtype composition (A215-A217).
  - Literature tally: 0 of 79 clinical claims adjusted for proliferation, 1 for subtype, 0 used random-set controls (draft,
    abstract level).
- **Open:**
  - Breast LumA tumours are more coherent at matched proliferation (A219, 3 cohorts, clock-specific).
  - Cis copy-number adjustment removes 63 % of that gap (A224), but a simulation shows additive copy-number noise cannot
    produce it (A225). The copy-number entanglement is real but not resolvable in bulk (biology vs tumour state).
  - TP53-mutant tumours without 17p loss have lower PER1 / PER2 (A224 H224.3, −0.16 log2), not adjusted for tumour state.
- **Not estimable:**
  - Host desynchrony (needs >= 2 normal tissues per cancer patient).
  - Single-cell per-cell questions (A214: 5-10x more depth needed).
  - Level-based damping readout (A223).
- **Proposed shape:**
  1. Tumours and their stroma keep body time.
  2. What the field read as clock disruption or prognostic clock information is proliferation, subtype composition,
     genomic position, design or sampling time; one main figure, the rest in the supplement.
  3. Open: tumour-state differences in coherence.
- **Questions:**
  - Is this structure acceptable?
  - How strongly should the literature correction be worded? It names published signatures, CYCLOPS magnitude and
    TimeTeller.
  - Should the A220-A222 tests appear as main or Extended Data?
- **On hold:** the cancer figure and its text. A226 continues.

### Q12 (2026-10-27): request brain time-of-death fields? (A141)
A141 (registered after your Q6 answer):
- 8 of 12 brain contrasts are noise-driven, so not estimable.
- 4 are anchor-dependent beyond noise: PD_mssm, VaD_mssm, ASCVD_hbcc and the PD atlas.
- Your Q6 answer foresaw a time-of-death request in this case. A draft (not sent) covers PsychAD on Synapse and the PD
  atlas on AMP PD, plus optionally ROSMAP: `paper/exploratory/A141_brain_stability/drafts/tod_request_draft.md`.
- **Question:** send the request(s), edit them, or not request?
- **On hold:** nothing else. The brain panels stay in Extended Data with the result stated.
- Evidence: `results_preliminary.md` (A141).

## Answered
### Q21 (2026-11-01): cancer L2 after A175/A176, and the Q8 wording
**Evidence** (`results_preliminary.md`, A175, A176, C1 replication scout):
- A175 H175.1 (TCGA, secondary): high-proliferation tumours pair less with their own normal tissue (D +0.158, p 0.021,
  clock-specific). ER split not estimable (MDE 0.41).
- A176: the coherence loss is "swamped", not "disordered", in 7 of 7 cohorts. In proliferating tumours the clock-plane
  variance across tumours falls clock-specifically. In ER− tumours it is mostly a generic rise in off-plane variance.
- So the lower pairing in A175 may be less clock signal rather than less host following. The two cannot be separated
  without a spread-matched pairing readout.
- Primary replication: GSE14520 + GSE66229 pooled MDE 0.214 > TCGA D 0.158. Stopped by the stop rule.
**Questions:**
1. Q8 wording: replace "clock loss" by "less clock signal relative to non-clock variation (swamped)" in Figure 5?
2. Register a spread-matched pairing readout (pairing at common tumour clock-plane variance, d* logic) on TCGA as a
   secondary? This would be a second amendment on the C1 question, which the one-amendment rule forbids without your OK.
3. Look for CPTAC (LUAD/STAD/BRCA with normal adjacent tissue) for a powered replication? This means a PDC/GDC download
   (open access), size to be checked.
**On hold:** any L2 wording for cancer beyond "secondary, TCGA". Nothing else is paused.
**Update (scout 2026-11-01, `exploratory/scouts_2026-11-01/SCOUT.md`):**
- A powered C1 replication may exist in GEO: the Nurses' Health Study breast tumour/adjacent pairs (623 unique patients with
  tumour + adjacent normal: 445 on GPL22920, 178 on GPL17585 (GSE115577 re-analyses the GSE93601 arrays); ER− 117 /
  ER+ 503; Affymetrix HTA arrays, all 14 clock genes annotated; `scouts_2026-11-01/prep_nhs/output/nhs_counts.csv`).
- For H1, GSE154200 has 186 persons with lesional + non-lesional skin (psoriasis 88, nummular 52, atopic eczema 46).
- Proposal: download both (open GEO, small), check clock-gene coverage and gates, run the MDE step, and bring the
  primary registrations to you before scoring. Go ahead?

**Update (2026-11-01, after the downloads):** the NHS replication is stopped: 445 usable pairs, tumour–normal pair
gate fails (p_acrophase 0.098) and H175.1 MDE 0.198 > TCGA D 0.158 (`results_preliminary.md`, pre-registration prep).
CPTAC (point 3) is the only route left to a powered C1 replication.

### Q22 (2026-11-01): two primary registrations to approve before scoring (C6 blind types, H1 skin)
**A183 (H1, GSE154200 skin)**, draft `exploratory/scouts_2026-11-01/prep_gse154200/REGISTRATION_DRAFT.md`: blind
replication of the psoriatic lesion delay (A54 held-out −0.93 h; MDE 0.81 h, 87 persons) and its specificity against
nummular and atopic eczema in the same pipeline (MDE 0.69 h). d* is not estimable (two compartments), so the HANDOFF
target "MDE on d* ≤ 0.15" is replaced by the paired offset. **Question:** approve as is?

**A182 (C6)**, draft `exploratory/C6_blind_types/REGISTRATION_DRAFT.md` (headers only read). My defaults in brackets:
1. GSE17536 (177 colon, Moffitt; GSE14333 could not be downloaded): [separate test in the weak/null family].
2. GSE39582 and GSE65858 are ComBat-corrected by the authors: [keep, caveat in the text].
3. Prostate GSE116918 (custom FFPE array, 13 clock genes, no NR1D1; the only prostate cohort ≥ 150): [include with 13].
4. ESCC GSE53625 (no symbol annotation on GEO; relation line in an earlier scout header): [drop].
5. Clock genes: [the pipeline's 14 (with RORA, CIART), as in A171-A173], not the HANDOFF list (BHLHE40/41).
6. Pancreatic, endometrial and cervical have no blind cohort ≥ 150: [accept the narrower family (colon, prostate);
   no near-miss cohorts].
7. HNSC and ESCA lean negative in TCGA (−0.17, −0.15): [keep HNSC as a test; state the TCGA value next to it].
8. Null reading: TCGA values of the weak types (−0.03 to −0.11) are below the expected MDE (about 0.2), so a
   non-significant cohort cannot say "no loss". [Add an equivalence test per cohort: T > −0.19 (the smallest
   replicated cohort effect, GSE62254) with a 90 % interval.]
9. GSE157009 / GSE157010 (two arms of one Alliance study): [separate tests].
**D2 (GSE27272, mother–fetus at term), optional third registration:** 60 mother + cord blood pairs; the cord blood ×
mother gate passes (p 0.022 / 0.013), placenta × mother fails. The only estimable question is one-group: is cord blood
coherent with its own mother's blood beyond shuffled pairs, clock-specifically and beyond the responsive-gene
companion (delivery stress)? [Register it as a secondary, marginal design (n 60); smoking not estimable.] Draft on
your yes.
**On hold:** scoring of A182 and A183. Their power steps are done (A183) or run right after the commit (A182).
- **Answer (PI, 2026-11-01):** Q21.1 yes, Figure 5 uses "less clock signal relative to non-clock variation
  (swamped)". Q21.2 yes, the spread-matched pairing readout is registered as A185 (second C1 amendment, allowed).
  Q21.3 scout CPTAC (metadata, commands for the user). Q22: A183 approved as drafted; A182 with the equivalence bound
  −0.19, registered tests GSE39582, GSE17536, GSE116918, GSE157009, GSE157010, GSE65858 descriptive, GSE53625 dropped,
  all other defaults; D2 registered as a secondary (A184). Also: Q12 the PI reviews and sends the draft; C8 not yet;
  A136 proceeds as drafted.

### Q13 (2026-10-28): roadmap phase 0 decisions (ROADMAP_2026-10-28.md)
I am going ahead with phases 1, 2 and 4 (code, data integrity, simulations) and the phase 6 scouts. None of them
depends on the answers below. The defaults in brackets apply when phase 3 starts, unless you say otherwise.
- **0.1 Reported coupling readout:** excess over the pseudo-clock median (Δc, with pseudo q95), or raw c?
  [default: Δc as main, raw c in the supplement]
- **0.2 Eligibility rule for disease primaries:** a compartment enters only if its control Δc exceeds the pseudo q95 and
  the MDE on d* is <= 0.2, registered as A153 before the contrasts are rerun. [default: adopt]
- **0.4 Anchor self-inclusion (C7): A155 meets the decision rule.** Anchor amplitudes are inflated by 9-20 % with 3-4
  anchors (4-6 % with 6), while offsets (<= 0.06 h), coupling (already leave-one-out) and the A123 map (r unchanged) are
  unaffected. Proposal: leave-one-out becomes clocksync's anchor amplitude and offset readout (opt-in now, default
  after your OK), and Figure 2 amplitude panels are rerun in phase 3. [default: adopt]
- **0.6 Parked items:**
  - A136: register now? [default: I draft the registration for your review, not scored]
  - A147 reference-size sensitivity: [default: run, labelled post hoc, 2 min]
  - A139 grid: [default: skip]
  - A138 Q2 simulator: [default: option (b), fix the clock signal-to-noise with a dated note]
  - Q12 and A119: yours (outward-facing / procurement).
- **On hold:** phase 3 onwards (anchors, backbone, A150 extension, disease reruns).
- **Answer (PI, 2026-10-31, review of STATE_2026-10-31 / HANDOFF_2026-10-31):** defaults accepted: Δc is the main readout (raw c in the supplement); family-level eligibility (as A153 note, Q15 a); leave-one-out anchor amplitude and offset become the default readout; 0.6 items as applied (A136 stays a draft for the PI).

### Q14 (2026-10-28): COARSE ordering choices (bone marrow, muscle)
The data audit found three ordering choices in the shared COARSE rules. They are not bugs of the "T cell" kind:
- "granulocyte monocyte progenitor cell" (21,500 cells) goes to `monocyte`;
- "mesenchymal stem cell" (28,275) and "skeletal muscle satellite stem cell" (13,531) go to `progenitor`, because
  progenitor comes before fibroblast_stromal and muscle.
- **Question:** keep them, or move GMP to progenitor and the two stem labels to fibroblast_stromal / muscle?
- **Default:** keep (the archive A9 choice). Changing them would rebuild organ_bone_marrow and organ_muscle and
  relearn their anchors.
- **On hold:** nothing. If you change them, I rebuild before phase 3.
- **Answer (PI, 2026-10-30):** keep.

### Q11 (2026-10-27): A140 Q1 simulator gate fails; rerun with the published CYCLOPS settings and a cleaner factor calibration?
The registered gate (arm D, whole-day sampling, CYCLOPS circular correlation with truth >= 0.7, seed 1,500) fails: median
0.067, 0 of 50 at >= 0.7. In the same replicates the clocksync gate passes 43 of 50, and the three CYCLOPS compartment phases
agree with each other (median |r| about 0.97) while ignoring the clock: CYCLOPS follows the shared non-clock factor. Per the
registration Q1 stops; H140.1-3 are not scored.
Two choices of the run may cause the failure, so the gate does not yet show that the published pipeline cannot work:
- CYCLOPS was trained with 10 models per fit instead of the template's 80, for cost: 80 models take about 80
  worker-minutes per fit, about 4-5 h on 16 workers for the 50 gate replicates.
- The non-clock factor was calibrated on the GSE266919 Pre samples, because the A17 tables hold clock genes only. These
  samples fail the gate, and their leading eigenvalues include Poisson noise, so the factor may be overstated.
  Post hoc diagnostic (8 replicates): factor scaled to 0.1 gives median truth r 0.30 (1 of 8 >= 0.7); without the
  factor, CYCLOPS crashed in its Julia code.
- **Question:**
  - (a) Accept the stop: report "the published pipeline does not recover a coupled clock at this factor size", as
    descriptive.
  - (b) Or register a revision before any rerun: 80 models, the factor calibrated on gate-passing tumour immune
    pseudobulks with all genes (needs those pseudobulks built), and the factor size as a grid. The gate is then rerun,
    and H140.1-3 only if it passes.
- **On hold:** A140 Q1 scoring and Q2 (ii) (the published pipeline on GSE266919).
- Evidence: `paper/exploratory/A140_dip_cyclops/README.md` (Q1), `output/q1/{simulator_gate,diag_summary}.csv`, `results_preliminary.md` (A140 Q1).

- **Answer (PI, 2026-10-27):** (b) revise and rerun: a dated revision registered before any rerun (80 CYCLOPS models; factor calibrated on genome-wide pseudobulks of gate-passing tumour immune compartments; factor size as a grid); gate, then H140.1-3 only if it passes.
### Q10 (2026-10-26): the A19 age model fails its null calibration in the first full count-level run
The archive validated `Dataset.covariate` (A19 part 2) with a model-based simulation (null rejection 0.025-0.08). The
migrated acceptance (M3 full run) keeps the real count tables and permutes ages within stratum: rejection at alpha 0.05
exceeds the registered 0.07 in 4 of 6 compartments for locked amplitude and for gamma (lung AT1 locked amplitude 0.220;
PsychAD excitatory neurons gamma 0.125, OPC 0.145). Lung AT1 and NK/ILC fits mostly sit at coupling 1 on the age subset.
- **Question:** (a) may I investigate (boundary fits on age subsets; the package's simulated null against the permutation
  null), and (b) if the count-level failure holds, should A19 part 2 (H19.1-2 nulls, the neuron gamma secondary, p <= 0.01
  per compartment against the simulated null) be re-read with a permutation p (ages permuted within stratum), as an
  erratum to the archive?
- **On hold:** any use of the A19 age readouts (neuron gamma) in figures or text; the erratum. M3 of resilience runs as is.
- Evidence: `results/extended_data/scripts/calibration/README.md` (M3 full run), `output/calibration/age_acceptance_*.csv`.
- **Answer (PI, 2026-10-26), part (a):** investigate. Diagnostic `paper/exploratory/A19_acceptance_check/check_null.py` (post hoc, real ages not read). Part (b) open.
- **Diagnostic result (2026-10-26):** the simulated null is too narrow in brain (SD ratio permutation / simulation
  1.10-1.28; 8-13 % outside its 95 %) and off-centre in the lung atlas (AT1 locked amplitude t median 1.84 vs 0.17; 34 %
  outside); the age fit is degenerate (own noise at floor) in 78-84 % of AT1 and NK/ILC permutations
  (`paper/exploratory/A19_acceptance_check/README.md`). Proposal for (b): add a within-stratum permutation null to
  `Dataset.covariate` (`null="permutation"`, default; refit per permutation; fits with own noise at the floor reported as
  not estimable for gamma), rerun the acceptance with it, then re-read A19 part 2 on the real ages as an erratum. Needs
  your go (clocksync change).
- **Answer (PI, 2026-10-26), part (b):** add the within-stratum covariate shuffle. Done: clocksync 613e224
  (`Dataset.covariate(null="permutation")`, default; gamma not estimable when own noise is at its bound). Running: the
  acceptance with it, then A19 part 2 on the real ages (erratum to follow). Note: the M3 port of A19 part 2 has no design
  check, so even under the simulated null the archive's neuron gamma secondary is not reproduced (-0.104, SE 0.077,
  p 0.173 vs -0.124, p 0.013; six compartments the legacy check excluded enter; `resilience/README.md`).
- **Result (2026-10-26):** with the permutation null the acceptance holds within Monte Carlo error (null rejection
  0.015-0.080; one cell of 24 at 0.080). A19 primaries are unchanged (null). Neuron gamma -0.116 (SE 0.054, p 0.031),
  with 2 of the 6 archive compartments at p <= 0.01. Erratum written (`results_preliminary.md`). The hold on the neuron
  gamma readout is lifted; it is a nominal secondary. Closed unless you object.

- **Answer (PI, 2026-10-27):** closed; the neuron gamma result goes in the main text as a nominal, uncorrected secondary (with the erratum).
### Q2 (2026-10-23): figure-table numbers used against their archive verdict (WP0 item 1)
- A104 rho +0.73 quoted as support; H104.1 predicted rho < 0 (archive: not supported; H104.2 not replicated).
- A54 psoriasis "replicated in discovery and held-out": the delay replicates (H54.2) but H54.1 and H54.3 fail.
- Mouse embryo "from E15": H108.1 passes in the late window but the window gate does not (p_acrophase 0.08).
- Quoted values that differ from their table: benchmark offset r 0.934 / 0.900 (quoted 0.91 / 0.88); baboon offsets
  0.455 h (0.5 h); healthy atlases 22 targets / 17 studies ("about 30").
- **Question:** remove, restate or keep each claim (with the verdict stated)?
- **On hold:** edits of the Figure 1, 2 and 5 tables.

- **Answer (PI, 2026-10-27):** go through each claim with the proposed text before editing.
  Item by item (PI, 2026-10-27): A104 restated as not supported (post hoc positive rho labelled); embryo restated with the
  borderline gate; psoriasis 'the delay replicates' with H54.1/H54.3 not met; numbers from the tables (0.934 / 0.900; 0.455 h;
  22 atlases, 17 studies, 14 tissues). Applied to the figure tables of `provisory_plan_paper.md`.
### Q3 (2026-10-23): panels not on the anchored A17+ estimator (WP0 item 4)
13 of 30 figure-table rows flagged (`paper/exploratory/WP0_integrity/output/estimator_consistency.csv`). By design
(the text must name the readout): r* (A103, A42), A104, A85 index, A92/A93 coherence C, A33 tissue host. Needing a
decision: Fig 1 A21 colon (an anchor decouples, V17.2 fails); Fig 2 GTEx organ departures (pre-A17 timed) and
neuron/retina panels (pre-A17 free fit); Fig 4 A20 "damping" (withdrawn in A22) and COVID A105 (2 anchors, no A22
confirmation); Fig 5 fibrotic-lung T cells (A5 free-fit 0.45).
- **Question:** drop, recompute on the anchored estimator, or keep with the readout named, per panel?
- **On hold:** recomputation of those panels.

- **Answer (PI, 2026-10-27):** lead's default per panel: by-design rows kept with the readout named; GTEx departures, neuron/retina and fibrotic-lung T cells recomputed on the anchored estimator; A20 'damping' and A21 colon dropped; COVID A105 to Extended Data with '2 anchors' stated.
- **Done (2026-10-27):** neuron/retina and fibrotic-lung T cells recomputed (results_preliminary.md, "WP0 / PI Q3"). PI (2026-10-27): GTEx departures kept with the timed readout named (no anchored version); neuron panel restated in ED (retina under both rules; cortex anchor-defined), message 'the brain's backbone is neuronal' dropped.
### Q6 (2026-10-24): brain disease contrasts are anchor-dependent (WP0 item 8)
With an all-LOO host instead of the A17 anchors, brain contrasts (Figure 4: PD, AD, ASCVD, VaD, SCZ ...) change sign in
1-7 of 9 compartments (block agreement r* 0.42-0.67); the peripheral headline contrasts (ILD, fibrotic lung, T2D, A44,
AD skin) keep their sign. Per-cohort cancer tumour-minus-normal coupling changes sign in 3 of 9 cohorts (near zero).
- **Question:** keep the brain disease panels with the anchor dependence stated, move them to Extended Data, or drop?
- **PI (2026-10-24):** decide evidence- and literature-driven whether one compartment has a good reason to be the brain reference (covariance, clock strength, literature). -> A135 (registered before any brain contrast is recomputed).
- **A135 result (2026-10-24):** no compartment family qualifies (excitatory neurons first in 2 of 5 control sets; rule 4 of 5). By the registered rule the brain disease panels go to Extended Data with the anchor dependence stated. Confirm, or keep open?
- **Time of death for the brain cohorts (checked 2026-10-24, column names only):** absent from every local donor table
  (PsychAD MSSM/HBCC/RADC, SEA-AD DLPFC/MTG, PD atlas) and from the public SEA-AD donor metadata (PMI only). Controlled
  portals not checked (no access from here): PsychAD on Synapse syn52160016; ROSMAP clinical (`ROSMAP_clinical.csv`, AD
  Knowledge Portal); PD atlas (Prashant et al. 2024, AMP PD Knowledge Platform). PI to decide whether to request it.
- **On hold:** Figure 4 brain panels.
- Evidence: results_preliminary.md, WP0 item 8.

- **Answer (PI, 2026-10-27):** brain donors die at any hour, so the brain cohorts likely cover the full cycle and the
  inferred amplitudes and phases are close to absolute. Compare them with the GTEx and simulation amplitudes to check
  whether the amplitude and phase spread is large enough for the model's estimates to be stable. If not, nothing can be
  said. If yes, describe first which compartments agree, and which should agree or are noise-driven. If they are not
  noise-driven, the PI may request the time-of-death fields to learn what follows it. -> to be registered (descriptive
  amendment) before any brain readout.
### Q9 (2026-10-24): direction after WP0-WP4 (step back)
Proposals (see the chat summary of 2026-10-24; each would be registered or written into the plan once agreed):
- (a) Figure 3 claim: from "a self-sustained peripheral network" to "organs stay coupled clocks when cues are disrupted
  (Hardy 4); perturbation acts in two separable ways: a shared displacement that keeps the relations (ventilated, not
  brain dead) and a loss of the relations (brain dead)". Drop the mechanism row of WP2 from the main text.
- (b) Figure 1: add the genetic positive control (A129: knockout liver decouples, specific) and the living cohorts
  (A113, with the post hoc estimator stated); drop "cue/locked split validated in mouse" (A112 failed).
- (c) Method: make "pull vs order" (signed pole-deviation Spearman and shrink ratio) a registered clocksync readout,
  validated on the A133 simulation grid before use.
- (d) Process: every registration states the attainable minimum p and the null-free ceiling of each statistic
  (H111.1 and A133 set unreachable thresholds).
- (e) Next tests, in order: A136 drug vs brain death within organ donors (glucocorticoid-response score), an external
  brain-dead / multi-organ replication, A134-A135 (running), then WP5-WP8 with WP5 reframed by A134.
- (f) Journal: the disease positive for Cell rests on cancer (being reframed) and the anchor-robust peripheral
  contrasts; Nature Methods is the stronger fit for the validated method.
  Correction (2026-10-24, PI comment): the psoriatic local-clock delay is the strongest pre-specified disease positive
  (archive: A54 H54.2 held-out -0.93 h, p 0.004; A78 core minus systemic cues -1.97 h, p 0.0002; A81b GTEx cue set
  -1.69 h, p 0.0003; other skin diseases -0.75 h, p 0.22). It is a delay (lesional clock later), not an advance. Open:
  trait vs state and treatment specificity (WP6) -> proposed to run WP6 before WP5.
- **Question:** which of (a)-(f) do you agree with?
- **Answer so far (PI, 2026-10-25):** (a) agreed; (b) agreed; (c) asked what "pull vs order" means (explained in chat,
  decision pending); (d) agreed (in the spirit of stating power); (e) cancer is a priority for the literature: CoFE
  (Ananthasubramaniam 2025) may have been misled by its method, and many papers used CYCLOPS 2.0 although it works only
  when samples cover the whole day. Determine whether anything can be estimated in tumours, or whether the variance
  between malignant cells and patients is too large. -> WP5 to be redesigned around this (proposal in chat, to register).
  Follow-up (PI, 2026-10-25): (c) agreed: pull vs order becomes a clocksync readout, validated on the A133 grid first.
  Both WP5 (redesigned) and WP6 will run; cancer first (A138). Literature agents may start now.
  A138 choices (PI, 2026-10-25): thresholds left to the lead (20 % / 50 %); CYCLOPS 2.0 tested both plain and with the
  healthy whole-day anchoring of Li et al. 2024 (assumes tumour = healthy Fourier coefficients); the CoFE authors were
  already told their tumour clock is wrong (no email needed). A138 registered (0162cc2).
- **On hold:** plan rewrite; A136 registration.


- **Answer (PI, 2026-10-27):** (f) journal decided after WP6. Plan rewrite and A136 registration: write both now as drafts for review.
### Q4 (2026-10-23): clocksync change
`clocksync.genes.gene_columns` returns -1 silently for a missing symbol (indexing the last column; the A85 bug).
Proposed: `missing="ignore"|"warn"|"raise"` and an all-absent guard
(`paper/exploratory/WP0_integrity/output/gene_audit/proposed_fix.diff`, with a test). clocksync holds your uncommitted
work, so it was not applied.
- **Question:** apply it (own commit, only these files)?
- **Answer (PI, 2026-10-24):** apply it. Done: clocksync 0a288ac (`genes.py` and a new test only; default behaviour kept, `missing="raise"` used in `paper/common/gtex.py` and `natexp_common`).
- **On hold:** nothing (no current result affected).

### Q7 (2026-10-24): what A110 can still say about the SCN
GTEx organ donors have no brain samples (0 of 301), so the registered SCN test (hypothalamus coupling, Hardy 3-4 vs
brain-dead donors) is impossible; only fast vs Hardy 3-4 remains (descriptive, being scored). The illness-duration test
of memory vs maintained coupling (H110.1) is infeasible: the interval is recorded for 16 of 103 Hardy 4 donors.
- **Question:** accept that GTEx cannot separate SCN loss from donor drugs, or look for another dataset (e.g. a
  brain-dead vs non-brain-dead ICU cohort with multi-organ or blood time series)?
- **Answer (PI, 2026-10-24):** accept that GTEx gives no statistically robust SCN result; the SCN reading stays a hypothesis.
- **On hold:** nothing now; the SCN sentence of Figure 3 stays a hypothesis.

### Q8 (2026-10-24): what the cancer "coherence loss" means (WP0 item 2)
A86/A87's lower share, coherence and alignment in malignant cells are not a depth effect (malignant pseudobulks are
deeper; a 4-fold shallower normal epithelium keeps coherence 0.68 vs malignant 0.23). But per-cell D is equal in
malignant and normal cells (3.89 vs 3.90), which rejects scrambled gene phases within cells (simulated 0.42). A87 also
shows 3-fold higher off-plane (non-clock) variance of the clock genes across tumours, with clock-plane variance kept.
- **Proposal:** restate the main cancer claim as "malignant cells keep a coherent clock in each cell; across tumours,
  clock-gene expression is dominated by non-clock variation that breaks the shared-phase structure", after a test:
  (a) A87 coherence within patients across cells (or clones), (b) coherence after removing the copy-number-explained
  part of each clock gene (WP5 H114 dosage). Register as an amendment before scoring.
- **Question:** agree with testing this before WP5, and with the possible restatement of the Figure 5 claim?
- **Answer (PI, 2026-10-24):** agreed -> A134.
- **A134 result (2026-10-24):** inferred copy number explains a small part (malignant coherence 0.23 -> 0.29, p 0.017 vs normal; about 1/10 of the gap); oncogenic programmes nothing; per-cell test uninformative. Figure 5 wording: keep "coherence lost" with the per-cell caveat and the copy-number share, or restate? (open)
- **On hold:** the Figure 5 cancer narrative; WP5 (A114) framing. The A86/A87 numbers themselves stand.
- Evidence: results_preliminary.md, WP0 item 2 (second and third parts).
- **Answer (PI, 2026-10-31, review of STATE_2026-10-31 / HANDOFF_2026-10-31):** the Figure 5 wording sub-question waits for C3 (in-plane vs off-plane variance: "swamped" vs "disordered").

### Q5 (2026-10-23): A132 outcome and the Figure 3 wording
Within ventilated (Hardy 0) donors, the push of pair phases toward 0 / 12 h is seen in both brain-dead organ donors and
ventilated post-mortem cases (not brain dead), with no separation (H132.1, MDE 0.44-0.54 h). Descriptively, only the
non-brain-dead group keeps the fast-death phase relations (rho 0.69, p 0.0006 vs 0.16, p 0.19 at equal n).
- **Question:** (a) word Figure 3's perturbed group as "ventilated donors" rather than "brain-dead organ donors"?
  (b) register a post hoc test (labelled post hoc) of the offset-correlation difference, organ donors vs ventilated
  post-mortem, as support for your SCN reading, or leave it descriptive?
- **On hold:** Figure 3 narrative wording; any A132 follow-up. A110 continues (H110.4 already uses organ donors).
- Evidence: results_preliminary.md, section A132.
- **Answer (PI, 2026-10-24):** describe both groups (the 0 / 12 h pull in both, phase relations kept only without
  brain death); a post hoc p for the difference is fine, labelled post hoc -> A132b (rho difference 0.517, p 0.008).
  PI reading: the pull is likely drug-induced; the kept offsets are hard to explain.

### Q1 (2026-10-23): Figure 3 framing after A109
H109.1 passes (organs share a clock-shaped state in every death group); H109.2 fails in ventilated donors (phase
relations collapse to 0 / 12 h; about half of the shared covariance is a displacement) but passes in Hardy 4
(secondary). Registered reading: "organs share a donor state", WP2 out of the main text. Proposed alternative: coupled
running clocks without cues, shown in Hardy 4 / Hardy 3, with ventilated donors as a shared displacement overriding
the clock relations. PI comment (2026-10-23): check first that Hardy 4 clocks are decorrelated from wall time, and
whether the 0 / 12 h push is seen in COVID ventilated patients -> registered as A109c (dab55ab), being scored.
- **A109c H109c.1-2 (2026-10-23):** no detectable wall-time alignment in Hardy 4 (R 0.097, p 0.56; MDE 0.145-0.208);
  offsets preserved after removing wall time (rho 0.750, p 0.0002), also without ventilator cases and in inpatients
  only. COVID (D109c.3-4, descriptive): Delorey organs uninformative (too few shared donors); Melms lung shows no
  push in COVID vs controls. The circularity ratio is confounded by host-phase concentration (sampling window), so the
  offset test (H109.2) is the decisive readout; the ventilated collapse to 0 / 12 h rests on it (p_closer 0.001).
- **Question:** after A109c, which framing (registered default, or the Hardy 4 framing as a registered amendment
  labelled post hoc relative to A109)?
- **On hold:** WP2 (A110); edits of the Figure 3 table and narrative.
- Evidence: results_preliminary.md, sections A109, A109b (and A109c when scored).
- **Answer (PI, 2026-10-23):** reframed claim accepted: when daily cues are lost through illness (Hardy 4), organs
  remain coupled clocks with their phase relations; in ventilated (mostly brain-dead) organ donors a shared perturbation
  overrides them. PI note: brain-dead donors have probably lost the SCN, Hardy 4 donors keep it with perturbed
  zeitgebers; caveat: drugs in both groups, not only cues. Follow-up: A132 (SCN loss vs ventilation), A110 rewritten
  for Hardy 4 / Hardy 3.


## Q15 (2026-10-29): A153 has no primaries as registered — which eligibility rule?
**Finding (A159, results_preliminary "Phase 3 ... A159"):**
- No Figure 4 control arm has a genome-wide table on disk.
- The registered per-compartment flag (c above both pseudo q95) is underpowered: a single compartment's pseudo q95 is
  0.6-0.9, so only 3 of 31 healthy compartments pass, all in breast. Heart fibroblast at c 0.72 fails.

So A153 as registered has an empty primary family. Changing the rule after seeing D159.1 is a design change (rule d).
**Options:**
- (a) **Family-level eligibility from healthy atlases** [recommended]. A family is eligible if its coupling, averaged over
  compartments and atlases, exceeds the matched and responsive q95 of the same average under pseudo sets (set j pooled,
  as H150.1 / H159.2). Disease primaries are then the compartments of eligible families with MDE on d* ≤ 0.2.
  - Uses A150 + A159 atlases (5).
  - Needs no disease genome tables.
  - Likely admits stroma and endothelium; immune would be decided by the test.
- (b) Build genome-wide tables for the non-brain disease datasets (about 18 contrasts; Census/h5ad sources), then apply (a)
  per control arm, or a per-arm family version of it. This is 1-2 days of builder work and disk, and also makes H159
  scorable on new atlases. It can follow (a) as a confirmation.
- (c) Drop eligibility: every estimable compartment with MDE ≤ 0.2 is a primary, and Δc is reported descriptively.
  Larger Holm family, so less power.

**Answer (PI, 2026-10-30):** start with (a); (b) is also good and can follow later if it seems important, after the
other tests. Confirmed in session on 2026-10-30: register the A153 note with rule (a) and run A153.

## Q16 (2026-10-29): fresh tumour-normal stromal pairing (roadmap 7.2) — worth registering on what is available?
**Scout** (`exploratory/scouts_2026-10-29/a151_fresh_scout.md`, metadata only). Few fresh paired single-cell cohorts with
stroma exist:
- Kyushu gastric (16 untreated pairs, verified; Seurat .rds of 18-22 GB, parsed in Python only);
- Li 2022 ccRCC (up to 12, unverified);
- OPSCC (≤ 7, unverified);
- oesophageal GSE173950 (conditional).

That is 16 verified patients and about 45 if the rest hold. Colorectal, liver and breast stay not estimable. The bounded
statistic flags small cohorts (n 4-11), so a fresh A151 would likely be an underpowered descriptive.
**Larger alternative (not fresh):** marker-based fibroblast re-annotation of cohorts already in A149/A151 that have 0
fibroblast pairs: GSE150290 23, GSE199654 13, GSE144240 10, GSE178481 9, GSE242299 8. These have been seen at the
malignant level but not for fibroblasts.
**Options:**
- (a) Register a fresh A151 on Kyushu + ccRCC + OPSCC, as the main-text candidate.
- (b) Register the re-annotation of the used cohorts, labelled "not fresh".
- (c) **[default, recommended]** Keep the pairing descriptive in Figure 5 (lead with A97/A127 coherence vs
  proliferation), and spend no more compute on 7.2.
- **Answer (PI, 2026-10-31, review of STATE_2026-10-31 / HANDOFF_2026-10-31):** (c): single-cell stromal pairing stays descriptive. The bulk host-referenced tumour–normal tests (HANDOFF_2026-10-31, C1) take over.

## Q17 (2026-10-29): coupling intervals undercover (A154 V154.3 fails) — how to report c?
**Finding:**
- With every compartment at the same nominal c, ĉ is biased low (c 0.3 → median 0.09-0.26), and 95 % bootstrap
  intervals cover in 22 of 24 cells only at 0.00-0.89.
- A noise-free diagnostic shows the cause is a host built from weakly coupled anchors: with precise anchors ĉ is 0.288
  for 0.3.
- Contrasts at a common host stay calibrated (V154.2 passes). Real anchors are the stromal-vascular compartments, with c
  about 0.6-0.7.

**Options:**
- (a) **[recommended]** A registered follow-up simulation (A162) at realistic anchor coupling (anchors c 0.6-0.8,
  targets 0.1-0.9), measuring coverage and bias. If coverage is adequate there, report c with the caveat "conditional on
  the anchor host".
- (b) Develop an attenuation correction (host reliability) in clocksync and validate it by simulation (several days).
- (c) Report c only as ranks and contrasts (Δc, d*), not absolute values with intervals.

**Default if no answer:** (a). I will register and run it.

**Outcome (2026-10-30):** A162 fails (V162.1 19/32 cells, V162.2 16/32; results_preliminary A162). Coverage holds at
c ≤ 0.3 but not at 0.6-0.9; the estimate is shrunk by a roughly constant factor (0.87-0.92 with anchors at 0.8), as an
errors-in-variables attenuation by the host's reliability. Per the registered reading, absolute c is reported without
intervals (ranks, Δc, d*), and a correction is developed (option b): registered as A163.

**Outcome of A163 (2026-10-31):** the frozen correction fails validation (V163.1 30/32, V163.2 27/32, V163.3 58/64
cells; V163.4 passes). It removes most of the attenuation: coverage at c_target 0.9 goes from 0.02-0.19 to 0.90-0.94
(median), but it misses at c_target ≥ 0.6 with n = 30, K = 4, and with an interval-width shortfall at c_target 0.9.
Per the registered reading, c stays reported without intervals (ranks, Δc, d*). **Decision for you:**
- (1) accept this, which is my default;
- (2) register a restricted claim (c_anchor ≥ 0.8 with n ≥ 100 or K = 8) on fresh seeds;
- (3) develop a better interval (e.g. a bootstrap that re-estimates R per replicate, or BCa) before a new validation.
- **Answer (PI, 2026-10-31, review of STATE_2026-10-31 / HANDOFF_2026-10-31):** option (1): absolute c is reported without intervals (ranks, Δc, d*); no further interval development now.

## Q18 (2026-10-30): A153 under rule (a) still has no primaries — where should the disease positive come from?

- Rule (a) makes fibroblast_stromal the eligible family (endothelial too in the 3-atlas sensitivity).
- Every eligible compartment of the 13 estimable disease targets has an MDE on d* > 0.2. The smallest are
  dcm_healthy heart fibroblast 0.232 and hlca_fibrotic EC_arterial 0.214. So H153.1-3 are not scorable, and the A17
  islets beta effect does not survive (results_preliminary, A153).
- **Nothing depending on this is paused.** A164 (GTEx T2D/obesity, matched) runs as registered, and roadmap 6.6
  (within-person designs) is next in the handoff order.
- **Options for the PI:**
  - (i) Relax the MDE bound to 0.25 in a new registration. That would admit only dcm_healthy heart fibroblast; it is
    post hoc in spirit, so I do not recommend it.
  - (ii) Option (b) of Q15. It does not change the MDEs.
  - (iii) Larger single-cell disease cohorts (an outward action: data requests).
  - (iv) Rely on within-person and GTEx designs for the disease positive.
- My recommendation is (iv) now, with (iii) as a longer-term ask.
- **Answer (PI, 2026-10-31, review of STATE_2026-10-31 / HANDOFF_2026-10-31):** (iv): the disease positive is sought in larger bulk within-person designs from routine clinical sampling (HANDOFF_2026-10-31, Priority 2, H1), with cancer C1 as the host-referenced tumour test. Scouts list only estimable bulk designs: at least 2 compartments (tissues, sites or sorted cell types) per person, since one compartment per person has no host and gives L1 coherence only (PI, 2026-10-31).

## Q19 (2026-10-30): within-person designs (roadmap 6.6) — the on-disk paired cohorts are already analysed; acquire Smillie 2019?

- Scout (metadata only, no clock statistic): `exploratory/scouts_2026-10-30/within_person/SCOUT.md`.
- **On disk, paired, but not blind:**
  - TAURUS (UC + CD): 29 patients with inflamed and non-inflamed biopsies; 23-24 have both eligible families
    (fibroblast_stromal and endothelial, >= 50 cells per state). Archive A45 and A53 already read it within patient:
    stromal coupling 0.80-0.91 in both states, phases unchanged, no amplitude loss. Inflamed vs non-inflamed is also
    confounded with pre/post treatment (pre 67/37, post 37/63 samples).
  - Psoriasis single cell (GSE173706, GSE228421): 12 patients with both families in lesional and non-lesional
    skin; used in A54/A58.
  - Kong 2025 (Crohn vs normal, same donor): 3-4 donors, too few.
- **New, independent:** Smillie 2019 UC atlas (18 UC patients with inflamed and adjacent non-inflamed biopsies,
  12 healthy; 366,650 cells; not on disk, pairing to be verified). Hosted on the Broad Single Cell Portal, which needs
  a login: an access step for the PI.
- **Nothing is paused.** I have not registered a TAURUS or psoriasis re-read: any new readout there would be chosen
  after A45/A53/A54 were seen, so it could only be a secondary.
- **Options:** (i) acquire Smillie 2019 and register a blind within-patient contrast (stroma and endothelium:
  d*, phase, z2, psi), with TAURUS as the not-blind secondary; (ii) register the TAURUS readouts as an explicitly
  non-blind secondary only; (iii) drop 6.6, since gut stroma looks resilient to inflammation (A45).
- My recommendation is (i).
- **Answer (PI, 2026-10-31, review of STATE_2026-10-31 / HANDOFF_2026-10-31):** bulk first: closed for now. H1 scouts paired bulk IBD series; Smillie 2019 is revisited only if no bulk design reaches MDE ≤ 0.15 on d* (the portal login stays a PI step). Scouts list only estimable bulk designs: at least 2 compartments (tissues, sites or sorted cell types) per person, since one compartment per person has no host and gives L1 coherence only (PI, 2026-10-31).

## Q20 (2026-10-31): proliferation coherence loss in gastric malignant cells — worth our own malignant-cell calling?
- **Context:** the proliferation clock-coherence loss is now replicated in breast (A97 discovery, A127 replication; A166
  H166.3 failed, underpowered), lung adenocarcinoma
  (A171) and gastric cancer (A173, two of three cohorts, clock-specific). It is not pan-cancer (A172 fails, power OK).
  In breast single cells, ER− malignant cells are less coherent (A170, p 0.044) but not clock-specifically, so
  the registered cell-intrinsic reading is not supported.
- **Gastric single cell (scout `scouts_2026-10-31/sc_gastric/`):** no public set has malignant-cell labels together
  with raw counts. The candidates are GSE270680 (25 tumour donors), GSE183904 (26) and GSE150290 (23). All three are
  public 10x / csv matrices, about 2.8 GB in total, and none has cell annotations. A test would need our own epithelial
  annotation and CNV-based malignant calling, a methodological choice made before scoring. By analogy with Wu 2021,
  perhaps 10-20 patients per set would have ≥ 50 malignant cells.
- **Nothing is paused.**
- **Options:**
  - (i) register a pooled three-set gastric single-cell test with a pre-specified CNV caller;
  - (ii) stop the gastric line at bulk;
  - (iii) instead, test whether "adenocarcinoma vs squamous" holds where both are large. This would mean the TCGA
    ESCA / HNSC / CESC rows, which are seen, so only as a secondary.
- **My recommendation is (ii) for now.** The bulk replication is the paper-relevant result. Single-cell n would be
  marginal, and the CNV step adds a degree of freedom.
- **Answer (PI, 2026-10-31, review of STATE_2026-10-31 / HANDOFF_2026-10-31):** (ii): gastric stops at bulk. Single-cell cancer work (C4) uses sets with author malignant labels (Wu + Pal breast, a labelled LUAD atlas after a scout).

## Q23 (2026-11-04): ARNTL was zero in the A18-era per-cell files — rebuild all and rerun the per-cell readouts?

- **Finding** (results_preliminary, integrity 2026-11-04): seven per-cell files (hlca_full, islets_t2d, natri2024_ild,
  parkinson_regions, psychad_mssm, psychad_radc, seaad_dlpfc) have ARNTL = 0 (Census symbol BMAL1). Pseudobulk results,
  the A152 files and the cancer per-cell files are not affected.
- **Affected readouts:** A18 part 2 per-cell amplitude / synchrony, A144 part 2 (descriptive, inconclusive), the
  Extended Data per-cell calibration. None is in the claims ledger.
- **Done:** builder fixed; HLCA rebuilt for A203.
- **Nothing else is paused.**
- **Options:** (i) rebuild the six other files (about 3 h of Census streaming, about 9 GB) and rerun the three readouts,
  with the old and new numbers side by side; (ii) drop the per-cell readouts from the paper and keep the erratum.
- **My recommendation is (i)** if any per-cell number stays in the Extended Data, else (ii).
- **Answer (user, 2026-11-04):** (ii) drop the per-cell readouts (A18 part 2 per-cell amplitude and synchrony, A144 part 2,
  the Extended Data per-cell calibration) and keep the erratum. Amended the same day: the malignant per-cell D of WP0 item 2
  stays (cancer files have ARNTL; malignant depth 4.3 clock UMI per cell, where D resolves gene-phase coherence). The
  six other per-cell files are not rebuilt. Basis: the dropped readouts rest on the ARNTL-zero files and were inconclusive; none is in the claims ledger. A203
  bears on donor-level precision, not on within-sample readouts such as D.
