# Questions for the PI (open first, answered below)

Each entry: date, question, what is on hold until it is answered, the evidence (with `src` pointers in
`results_preliminary.md`). Work that does not depend on an answer continues. When answered, move the entry to
"Answered" with the answer and its date.

## Open

### Q2 (2026-10-23): figure-table numbers used against their archive verdict (WP0 item 1)
- A104 rho +0.73 quoted as support; H104.1 predicted rho < 0 (archive: not supported; H104.2 not replicated).
- A54 psoriasis "replicated in discovery and held-out": the delay replicates (H54.2) but H54.1 and H54.3 fail.
- Mouse embryo "from E15": H108.1 passes in the late window but the window gate does not (p_acrophase 0.08).
- Quoted values that differ from their table: benchmark offset r 0.934 / 0.900 (quoted 0.91 / 0.88); baboon offsets
  0.455 h (0.5 h); healthy atlases 22 targets / 17 studies ("about 30").
- **Question:** remove, restate or keep each claim (with the verdict stated)?
- **On hold:** edits of the Figure 1, 2 and 5 tables.

### Q3 (2026-10-23): panels not on the anchored A17+ estimator (WP0 item 4)
13 of 30 figure-table rows flagged (`paper/exploratory/WP0_integrity/output/estimator_consistency.csv`). By design
(the text must name the readout): r* (A103, A42), A104, A85 index, A92/A93 coherence C, A33 tissue host. Needing a
decision: Fig 1 A21 colon (an anchor decouples, V17.2 fails); Fig 2 GTEx organ departures (pre-A17 timed) and
neuron/retina panels (pre-A17 free fit); Fig 4 A20 "damping" (withdrawn in A22) and COVID A105 (2 anchors, no A22
confirmation); Fig 5 fibrotic-lung T cells (A5 free-fit 0.45).
- **Question:** drop, recompute on the anchored estimator, or keep with the readout named, per panel?
- **On hold:** recomputation of those panels.

### Q6 (2026-10-24): brain disease contrasts are anchor-dependent (WP0 item 8)
With an all-LOO host instead of the A17 anchors, brain contrasts (Figure 4: PD, AD, ASCVD, VaD, SCZ ...) change sign in
1-7 of 9 compartments (block agreement r* 0.42-0.67); the peripheral headline contrasts (ILD, fibrotic lung, T2D, A44,
AD skin) keep their sign. Per-cohort cancer tumour-minus-normal coupling changes sign in 3 of 9 cohorts (near zero).
- **Question:** keep the brain disease panels with the anchor dependence stated, move them to Extended Data, or drop?
- **PI (2026-10-24):** decide evidence- and literature-driven whether one compartment has a good reason to be the brain reference (covariance, clock strength, literature). -> A135 (registered before any brain contrast is recomputed); still open.
- **On hold:** Figure 4 brain panels.
- Evidence: results_preliminary.md, WP0 item 8.

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
- **Question:** which of (a)-(f) do you agree with?
- **On hold:** plan rewrite; A136 registration.

## Answered

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
- **On hold:** the Figure 5 cancer narrative; WP5 (A114) framing. The A86/A87 numbers themselves stand.
- Evidence: results_preliminary.md, WP0 item 2 (second and third parts).

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

