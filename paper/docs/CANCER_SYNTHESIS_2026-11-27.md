# Cancer and clock organisation: what we learned, what to re-read, what to test (2026-11-27)

Scope: every registered and post hoc cancer analysis to date (A97-A405; archive A17, A83-A87). Numbers are from
`results_preliminary.md` under the amendment named; claim wording follows the claims-ledger addendum of 2026-11-27.
Status: **A** registered, clock-specific, holds on the core loop or calibrated companions; **B** registered with a stated
limit, or post hoc; **C** descriptive or a lead. Draft for the PI; not reviewed.

## 0. Working picture

A tumour keeps its patient's circadian time. The tumour stroma carries it in every cancer resolved, fibroblasts most
faithfully (in phase with the patient's normal stroma, about 30 % quieter). Cancer cells join only in some cancer
types (lung adenocarcinoma clearly; thyroid and kidney lean), and those are the cancer types whose bulk tumours look
clock-structured. Proliferation dilutes how loudly the tumour expresses the patient's time (gain about two-thirds) without
moving it away from that time. Immune cells share a body-wide state with the stroma that stress-responsive genes share
too; lymphocytes carry no readable clock of their own. Time-of-day signals in tumour immune programs are seen against the
hospital clock, but whether they follow the patient's internal clock cannot be decided with the time references we have.
No clock readout predicts outcome beyond proliferation. The bottleneck of the next round is the time reference, not the
sample size.

## 1. What we learned

### 1.1 The tumour keeps the patient's time (bulk, tumour-normal pairs)

| Finding | Evidence | Limit | Status |
|---|---|---|---|
| Tumour clock state pairs with the same patient's normal tissue | A138 H138.4: 0.365 vs null 0.122, matched q95 0.261 (517 TCGA pairs, 10 cancers); calibrated (A285, A320) | Core loop 0.252 < matched q95 0.275 (A336: at the line); only a dispersion of 2.4 h or more is detectable | A (14 genes), at the line (core) |
| Not carried by the non-malignant fraction | A138 H138.8: low- minus high-purity D +0.025 (MDE 0.16) | Falls to "not decidable" on the core (A336) | B |
| Published "uncoupling" is a sampling-window artefact | A265 (COFE says "uncoupled" in 7 of 10 coupled simulations), A267 (COFE's pairing is generic patient similarity) | - | A |
| "Disruption" scores coexist with pairing | A207 H207.1-2 PASS (CCD, Grabe order) in the tumours that pair | ΔR² negative for any gene set (in-sample phase) | B |
| Clock genes are lower in tumours, not uniformly | A221: PER1 -1.07, PER3 -0.74 log2; NR1D1 +0.38, NPAS2 +0.30; GTEx-as-normal designs exaggerate and invert ARNTL, NPAS2, NR1D1 | Levels, not rhythm | A |
| Per-cancer coupling to body time | A304 A3 (descriptive): LUAD 0.54, LIHC 0.38, KIRC 0.27, PRAD 0.26, BRCA 0.18, KIRP 0.12, THCA 0.11, COAD 0.06 | Leans on GTEx reliabilities | C |

### 1.2 Proliferation lowers how loudly the tumour expresses that time, not the time itself

| Finding | Evidence | Limit | Status |
|---|---|---|---|
| Gain falls, departure does not rise | A304: gain ratio 0.63 (Q_G -0.462, p 0.0055, below matched and responsive); core -0.668; departure Q_R -0.190 (power 1.00 for a doubling) | TCGA only; one normal per patient; bulk cannot say where the loss sits | A |
| Pairing drops with proliferation | A175: D +0.158 (TCGA, clock-specific); A186 CPTAC D 0.113, not clock-specific (90th / 85th percentile) | Replicates in size, not in specificity | B |
| Tumour-only coherence falls with proliferation or ER-negativity in some cancers | SCAN-B -0.123 / ER -0.343 (A127); METABRIC ER -0.229 (A166); LUAD -0.242 / -0.293 (A171); gastric -0.187 / -0.533 (A173) | Not pan-cancer (A172: -0.032, MDE 0.087); not METABRIC proliferation, not liver; damping vs rewiring not separable in bulk (A212) | A (where replicated) |
| Proliferation per se does not swamp a clock | Regenerating liver -0.027 (MDE 0.064, A211); no recovery on arrest (A168 -0.010, MDE 0.196; A196) | Moderate rises not excluded | B |
| Genomic states do not add beyond subtype | TP53, MYC, copy-number losses in TCGA (A210, A213, A215) vanish at matched subtype (A217); luminal A advantage +0.276 -> +0.101 after cis copy-number adjustment, not clock-specific (A224) | - | B |
| Where the gain loss sits is open | Tumour stroma not quieter in proliferating tumours, composition unchanged (A391, descriptive); malignant-cell tests not estimable (A187b, A188, A195, A214, A305) | Power and gates, not negatives | - |

### 1.3 Who carries the time inside the tumour (single cell)

| Compartment | Evidence | Status |
|---|---|---|
| Fibroblasts | A346: 0.686 / 0.683 (core), beyond matched, narrow and broad responsive. A387 full model: c 0.439 / 0.403, peak -0.07 h (-0.77 to 0.63), loudness 0.69 (0.51-0.86). Coupled in every cancer resolved (0.37-0.76, A387 note 1) | A |
| Endothelium | Beyond matched on 14 genes only (A346 core 0.510; A387 core c 0.481 < matched q95 0.515). Peaks before fibroblasts in normal organs, -0.49 h (-0.82 to -0.16), inside the matched band; tumours -0.26 / -0.87 h, descriptive (A389) | B / C |
| Myeloid, T cells | Pair with the patient's normal clock beyond matched (A227: 0.581, 0.539), not beyond responsive (A271). Myeloid coupling to normal stroma 0.229, amplitude 0.31 of fibroblasts, offset +0.07 h, not beyond responsive (A404) | B |
| Lymphocytes | No readable own clock (A281; A404: T/NK and B couple weakly, far below companions) | A (as a limit) |
| Cancer cells | Barely coupled overall: c 0.13-0.15, not beyond matched (A388); c 0.100 / 0.032, loudness 0.32 (A387). Per cancer: thyroid 0.392, kidney 0.330, LUAD 0.290 (z 2.37); breast 0.055, ovary 0.011, colorectal 0.004 (A401). LUAD 0.66 beyond matched (A388, descriptive) | A overall, C per cancer |

### 1.4 Linking cells to bulk

- Cancer types whose cancer cells keep the tumour's time are those whose tumours keep a clock-shaped structure in bulk:
  A401 H401.2 PASS (secondary) rho 0.659, p 0.043 over 8 cancers; leave-one-out p 0.045-0.118 (fragile); standardised
  coupling against bulk pairing rho 0.667 (descriptive). **B**
- The purity-weighted model (stroma plus cancer cells) does not rank the cancers (H401.1 FAIL, rho 0.02): equal stroma
  terms dilute the cancer-cell term.

### 1.5 Immune context and time of day

| Finding | Evidence | Status |
|---|---|---|
| Louder stromal clocks go with colder, less myeloid tumours | A403 slope -0.102 (p 0.025, Holm 0.050: FAIL), beyond both companion families. A403b: shared across cancers (I² 0.13), not proliferation (-0.105 adjusted), core loop not decided (38.6th percentile of drops), not stromal inflammation (-0.019); myeloid share -0.083 with BRCA -0.419 and CESC -0.556 beyond companions | C (lead) |
| Composition does not follow the stroma clock phase beyond companions | A229 T/NK share (0.378 vs 0.302, borderline) not replicated: A402b excess 0.025 (MDE 0.073). Myeloid leans (excess 0.044, p 0.07); its shared peak hour across cancers is inside the matched sets (92.5th percentile, A402b addendum 1). Programs: A246 omnibus -0.121 (MDE 0.85) | B (bound) |
| Breast tumour programs follow the recorded hour of surgery | A146 H146.2: S 1.339, Holm 0.0012 (T-cell activation t 3.48, PD-1 3.22, MHC-I 2.15); T-cell and myeloid abundance do not | B (one cohort) |
| Hour vs internal time not decidable | A146 addendum 1: estimator error 1.86 h on a 10 h window; programs that follow the hour exactly would give S_internal 0.269 (90 % -0.317 to 0.844); observed -0.056 (18th percentile) | - |
| Stromal glucocorticoid state follows the clock phase | A402b D2: r 0.385 (p 1e-4), beyond matched, inside responsive (expected for rhythmic GC targets) | C |
| The preprint's human TNBC phase result is not supported | A140: immune compartments share a non-clock state (gate fails, clock share 0.025 vs null 0.028) | A (as a limit) |

### 1.6 Clinical readouts

- No clock readout informs outcome beyond proliferation: published signatures (A192), TimeTeller (A194), CYCLOPS
  magnitude (A193), our tumour signal (A191), treatment interaction (A198); hazard ratios 0.92-0.98. **A** (bounds)
- Endocrine therapy does not lower coherence (A168, A169). Anti-PD-1 (BIOKEY): no visible difference at about one SD
  (A405). The chronotherapy window cannot be read from TCGA or CPTAC pairs: the positive control fails (A249, A250).

### 1.7 Lessons of method that every next test must carry

1. **The time reference is the bottleneck.** Per-patient internal time (A200: Spearman 0.558; residual 1.86 h) and the
   single-cell host phases are noisy; a null against them bounds what the reader sees, not the rhythm (A146 addendum 1).
   State the expected attenuation with every internal-time test.
2. **Matched genes are necessary, not sufficient.** Immune pairings pass matched and fail responsive (A271, A404); the
   narrow responsive q95 needs placebo calibration (A285); the broad family is calibrated but lets a hormonal pulse
   through (A320).
3. **Core loop and the gene-drop null.** The tumour pairing is at the line on the core (A336). Random gene sets also share
   a phase across cancers (A402b addendum 1): cross-group phase agreement needs common-axis companions.
4. **Covariates.** Subtype and cis copy number confound tumour-group contrasts (A217, A224); check the direction of any
   purity imbalance (A304).

## 2. Earlier results that tools built since could update

Tools named here exist in `clocksync` (dates from its log) or in the analysis scripts (A304, A401, A402b, A146 addendum).

| # | Earlier result | Read then with | Tool now available | What the re-read would settle | Cost | Priority |
|---|---|---|---|---|---|---|
| U1 | A403 / A403b colder-louder (lead) | log host radius, GTEx waveform | phase harmonics removed from the radius (`elasticity.radius_covariance` with sampling-time covariates, A336c; ellipticity check of A400c) | whether the loudness differences are phase differences seen through a mismatched waveform (A402b: the myeloid share also leans with phase) | hours | high |
| U2 | A304 gain / departure (TCGA only) | A304 scripts, 517 pairs | the same readout on CPTAC-3 (713 pairs on disk, A186) | replication of the main-text proliferation statement in a second paired cohort | hours | high |
| U3 | A246, A402b, A146 part 2 (internal-time nulls) | host phase taken as exact | split-half host-phase reliability (`deviation` split genes), noise injection (A146 addendum 1), `simulate` with a mouse-sized infiltration rhythm at working-hours sampling | what these designs could detect on the true-phase scale: "no relation" vs "unreadable" | 1 day | high |
| U4 | A138 H138.4 tumour pairing (at the line on the core) | pooled pairing, 100 + trio drops | full model (`Dataset.fit`, attenuation `df`, d*), all 364 drops, broad responsive family, `deviation.incoherence` (patient dispersion in hours) | firm or settle the abstract's tumour sentence; a dispersion bound in hours | 1 day | high |
| U5 | A212 / A391 damping vs dephasing (not separable) | bulk plasmode; underpowered contrasts | `orientation` (radial vs tangential own variation, calibrated by simulation, A311) | whether cancer-cell clocks vary between patients in amplitude (damped) or in phase (robust oscillators, dephased) | 1 day | medium-high |
| U6 | A187b, A188, A195 cancer-cell clock and cell cycle (not estimable) | 8-17 early cohorts, malignant gate | A401 pooled table (12 cancers, 574 patients) restricted to cancers with coupled cancer cells; CoPhaser cell-cycle phase where median UMI > 1.5k; `contrast(host_prior='arm')` | whether cycling cancer cells are quieter or less coupled where cancer cells keep time | days | medium |
| U7 | A249 / A250 chronotherapy window (positive control failed) | bulk pairs | single cell with the stroma host in cancers with coupled cancer cells; fibroblast clock-controlled genes as positive control | whether malignant clock output is readable where cancer cells keep time | days | medium |
| U8 | A271 / A404 immune coupling (not beyond responsive) | narrow responsive q95 | placebo-calibrated p (A285), broad family (A320) | little expected change; completes the companion set | hours | low |
| U9 | Tumour-only coherence family (A127-A224) | coherence, matched sets only | core loop with 364 drops, responsive families, `decompose.signed_coherence` (sample-level offsets fitted out) | supplement robustness only; the main text uses A304 | 1 day | low |
| U10 | A191-A198 prognosis | tumour-only signal | per-pair gain or contribution on own normal | expected null beyond proliferation | - | not now |

**For `clocksync`** (useful beyond this project): the noise-injection null for a time reference, common-axis companions for
phase agreement across groups, and the A304 gain / departure readout live in scripts; moving them into the package is a
small task (part of the code hardening the PI schedules).

## 3. Hypotheses and analyses to test, ranked

Each is registered (with a power stage) before its data are scored.

**H1. Cancer cells keep time while they keep their lineage identity.** Links A401 to cell identity.
- Within coupled cancers, the per-patient cancer-cell coupling contribution against a differentiation score (LUAD: AT2
  program; THCA: thyroid differentiation score; KIRC: proximal-tubule program).
- Across cancers, cancer-cell coupling against the similarity of malignant cells to their normal epithelium.
- Precondition: split-half reliability of the per-patient contribution. A197 (bulk coherence vs an ES-like index, null)
  is a different readout.
- A positive result reads "lineage-dependent participation"; a null bounds it.

**H2. The patient's clock loudness reaches the tumour.** Connects the cancer chapter to the GTEx person volume (A360,
A361).
- A387 table, tumour and normal stroma per patient: amplitude elasticity of tumour fibroblasts to the patient's shared
  amplitude (`elasticity`), person-level radius covariance of normal and tumour stroma, glucocorticoid-target state of
  the normal stroma against tumour loudness.
- Also decides whether the A403 myeloid link is systemic (present in normal stroma) or local (tumour only).

**H3. Cancer-cell clocks are damped, not dephased** (U5). Orientation tilt of cancer cells against tumour and normal
fibroblasts. A radial tilt fits a gain loss without departure (A304); a tangential tilt points to phase scatter.

**H4. Proliferation lowers the amplitude of cancer cells where they keep time.** LUAD, THCA, KIRC: cancer-cell loudness or
contribution against the G2M index or the CoPhaser cycling fraction (U6). LUAD now has about 80 patients (A188 failed at
about 27). A single-cell counterpart of A304.

**H5. A304 replicates in CPTAC** (U2). Registered replication: gain ratio and departure, high vs low proliferation, both
companion families, core loop. The CPTAC pairing was scored in A186; its gain and departure were never computed.

**H6. The programs that follow the surgery hour follow the stroma clock phase.** T-cell activation, PD-1 and MHC-I (the
A146 leaders) against the tumour-stroma host in the pooled table, per program in r units (per-program MDE 0.053-0.075,
A402 power), after U3 has said what this time reference can detect.

**H7. Tumour stroma does not desynchronise locally.** Equivalence bound on tumour-fibroblast incoherence in hours against
normal fibroblasts (`deviation.incoherence`), with the psoriatic delay (-1.45 h, A183) as the size a local break reaches.
A clean statement for "where clocks break, they break locally".

**H8. The cell-to-bulk ranking replicates in new cancer types.** Extend the A401 ranking (fragile, 8 cancers) with at least
four cancer types that have an A83 value; data cap raised to 700 GB for tumour data (2026-11-26). Register before
downloading.

**H9. Hour vs internal time, decisively** (needs outside input).
- A second timed cohort: the REMAGUS biopsy times (outward request: user's decision).
- An internal inventory of public tumour cohorts with recorded collection times.
- The lab's Bieler 2014 time-lapse data for damping vs dephasing with division (PI question Q44a).

**Not recommended now:** prognosis re-reads (U10); further tumour-only group contrasts (subtype and copy-number confounding);
glioma (A247: precondition fails).

## 4. Suggested order and decisions

1. **Order:**
   - First the cheap checks of current claims: U1, U2 / H5, U3.
   - Then the new biology on data on disk: H2, H3, H1, H4.
   - H6 and H7 follow U3 and U4.
2. **PI:**
   - The cancer chapter wording (claims-ledger addendum 2026-11-27, rows 6.6, 6.10, 6.13-6.18).
   - Whether the stroma clause returns to the abstract. The core-loop re-test the story audit asked for is done (A346,
     A387): fibroblasts hold beyond all companions; endothelium does not.
3. **User:**
   - The REMAGUS request (outward-facing).
   - Downloads for H8.
   - Access to the Bieler data (Q44a).
