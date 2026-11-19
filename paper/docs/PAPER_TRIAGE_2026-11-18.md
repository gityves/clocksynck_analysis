# Paper triage (2026-11-18): which results go where, which do not

Built from the inventory (`triage/inventory_2026-11-18.tsv`, Sonnet extraction, one row per heading) of the 208 scored entries of `results_preliminary.md` (A109 onwards; the archive A1-A108 is
not re-triaged here), STATUS_2026-11-18, the claims ledger, PI_QUESTIONS and the session of 2026-11-18 (A310-A316,
literature cards). IDs refer to `results_preliminary.md`; numbers are in those entries. Draft for the user and the PI;
venue and the split are their decision.

## Two papers, or one

- **B (biology):** "The human body keeps time as a federation": each organ runs its own clock, built from
  cell-autonomous clocks whose loudness is set by cell type and composition; the order between organs exists while a
  living brain drives it. Target: Cell (Nature / Science as alternatives).
- **M (method + cancer):** untimed, host-referenced inference in hospital hours, with its demonstrated edge on real data
  (orderers lose the GTEx organ map in 08-18 h deaths, clocksync keeps it, A123b) and cancer as the application (tumours
  follow their patient's clock, lose loudness not timing as they proliferate; COFE's uncoupling is a sampling artefact).
- One paper = B with cancer as its perturbation figure and M's validation as Extended Data. The tags below allow either.

## Rules used

- **Main:** registered before scoring, powered, clock-specific beyond expression-matched genes (and stress-responsive
  genes where that rule is calibrated: A285 says it is not in bulk), not superseded, answers one question of the paper.
- **Supplement / Extended Data:** calibrations and validations; nulls that bound a claim; descriptive results labelled
  as such; tested limits of the model.
- **Exclude:** failed gates or not estimable (uninformative by rule); superseded or contradicted; not clock-specific
  where specificity is the claim; per-cell readouts dropped by decision (2026-11-04); exploratory dead ends.

## Paper B

Order revised 2026-11-19 (`SYNTHESIS_2026-11-19.md` section 4): cells first, climbing to the body; no organ map
or tissue ranking as a main result (CHIRAL overlap); brain death as the high point, not the opening. Sections are
figures.

| section | main | supplement | notes |
|---|---|---|---|
| 1. Who keeps time inside organs (opening figure) | A272 + A274 + A283 + A280 H280.1 + A282 (stromal-vascular hub, replicated, calibrated); A281 + A287 + A284 + A152 (tissue lymphocytes: no readable clock); A288 (macrophages); A273, A279 (cell classes share the organ phase); A252, A256, A251 (surgical part), A255 (one internal time advancing with elapsed hours in living tissue: last panel) | A275 (hormone and immediate-early states shared, clock not), A294 (heart lymphocytes, breast hormonal; core loop open), A201 (aorta shift together, descriptive), A276-A278 (calibration history), A113/A130 (living cohorts, sign note); method validation: A123b, A145, A209, A285, A258, A200 (or M) | A150 / A263 Part B: anchored atlas couplings not specific against responsive genes: use the pairwise results, not the anchored atlas figure. A285: living fat-skin not specific against responsive genes; the surgical parts stand. Method shown as one schematic panel only |
| 2. What sets loudness: cell class and composition | A293 Part 3 (relative class gain conserved across organs); A293 Part 2 + A295 (organ sets phase; composition moves gain, not phase); A291, A294b (+ correction note), A271 Part B, A144 part 1: gain over a clock-specific departure of near-constant size | A123 / A123c (organ offsets and amplitude ranking: concordance with CHIRAL, Extended Data only); A312 (persons differ in shared amplitude; gradient small), A302, A313 (departures not shared by meal-driven or adjacent tissues), A259 (tension with A294b, report) | coupling = S / (S + D) stated as gain over an added floor, not oscillator stiffness (A311). No organ map or tissue ranking as a main result (overlap with CHIRAL Figure 1) |
| 3. The order between organs needs a living brain | A293 Part 1 (H293.1 slow deaths, H293.2 ventilated living brain), A109 H109.1; A317 H317.1 (brain death rearranges the order, registered, clock-specific); A318 H318.1 (survives ischaemic time and RIN adjustment); A315 (amplitude kept, registered, powered); A319 (mouse central-clock loss: substitute drive, new order; explants, not blind) | A132, A132b, A314 (descriptive / post hoc), A317 H317.2 (shared clock kept) and sensitivities, A318 H318.2 (ischaemic gradient in the confound direction); A231 (age), A208 (sex), A300 | **gap 1** remains as replication outside GTEx. Brain-body split (A142, A143, A263 Part A, A289): two domains, not clock-specific: supplement as a limit, not a result |
| 4. Drives and conflicting cues | A297 (liver self-sustained), A257 (night shift: 2.7 h, amplitude to about a quarter), A262 D1 (forced desynchrony, blood) | A261, A260, A270, A268; model reading (central + behavioural drive, weight b) as a box | A124, A156 superseded by A257 / A262 (MDE too large) |
| 5. Perturbations | chronic disease nulls: A164, A233, A235, A306 (powered); psoriasis: A183 + A308 (delay clock-specific), A254 (other lesions do not carry it) | A202, A180 (equivalence), A115 (uninvolved skin), A234 (age advance: inflated, living cohorts opposite sign) | cancer here only if one paper (see M) |
| 6. Model and its tested limits | one panel: autonomous cellular clocks, weak permissive local coupling, systemic drive keeps organ order | A310 (population fit inadequate), A311 (departure not phase-dominated), A312, A313, A315, A319 | literature cards: Sinturel 2021, Guan 2020, Noguchi 2013, Finger 2021, Liu 2007, Manella 2021, Izumo 2014, Husse 2014, Tahara 2012, Saini 2013 |
| ED. Before birth (Extended Data) | A243, A244, A242 with A264 / A266 / A269 / A307: a driven state carried by signal-induced genes before an adult-size loop; A309 (mouse heart and kidney negatives) | A236, A237-A239, A240, A241, A184 (negatives and bounds) | A189 / A190 superseded by A264 / A266 |

## Paper M (or Extended Data of a single paper)

| section | main | supplement | notes |
|---|---|---|---|
| Method in hospital hours | A123b (orderers fail in 08-18 h, clocksync holds), A123, A174 (heart V3), A206, A209 / A209b (no false clock claims), A155 | A145 (simulation: no significant edge after Holm; report honestly), A139, A147, A140 gates, A199, A204 / A204b (TimeTeller), A205, A157/A158 | **gap 3**: A162 / A163, coverage of c fails at realistic anchor coupling and the df correction did not validate: intervals of c must be stated as approximate or fixed |
| Calibration | A285 (matched rule calibrated, responsive not, in bulk), A283, A258 (unequal host spread: compare at d*) | A154, A203, A180, A273, A279 | the d* rule and the matched-set rule are the method's specificity backbone |
| Tumours follow their patient | A138 H138.4; A265 + A267 (COFE's uncoupling is sampling) | A138 H138.8 (purity), A149 / A151 notes (stroma follows normal stroma; malignant not clock-specific) | the method-exclusive positive result of METHOD_EDGE |
| Loudness, not timing | A304 (gain falls with proliferation, departure does not rise; powered, clock-specific), A221 (clock genes lower in tumours, clock-specific), A127 + A127 post hoc (SCAN-B coherence loss with proliferation, ER-negative within PAM50) | A166, A171, A173 (partial replications), A172 (not general), A176, A212 (swamping excluded; damping vs rewiring not separable in bulk) | |
| Microenvironment | A227 against matched genes + A271 Part A (fibroblasts and endothelium keep the patient's clock; immune pairing not specific against responsive genes) | A228, A229, A246, A245 | |
| Clinical claims tested and rejected | A220 (published circadian risk signatures do not beat random signatures), A222 (no selection on clock genes) | A191, A192 (prognostic, not clock-specific), A193, A194, A198, A249, A250 | useful negative results for the field |

## Exclude (reason)

- **Gate failed / not estimable / not powered:** A132 H132.1, A136, A109b, A109c D109c.3-4, A140 Q2, A148, A153, A160,
  A165 (most rows), A167, A182, A185, A187, A187b, A188, A195, A196, A214, A223, A225, A226, A230, A247, A253, A286,
  A305, A292 (power below 0.8), A316 (feasibility failed), A181.
- **Superseded or corrected:** A126 (by A161), A189 / A190 (by A264 / A266), A210 / A213 / A215 / A217 (subtype post
  hoc), A219 / A224 / A219 qualification (luminal A advantage withdrawn), A166 erratum use, A176 "swamped" wording,
  A142 digestive claim (holds only reliability-corrected), A296-A303 adrenal offset (an effect of death: one line at most).
- **Not clock-specific where that is the claim:** A122, A170, A178 (GSE62254), A186, A207 H207.2, A263 Part A as a brain
  result, A268, A232.
- **Dropped by decision:** A144 part 2 and all per-cell lambda / R readouts (2026-11-04).
- **Dead ends:** A121, A134, A146 parts, A168, A169, A179, A197, A211, A161 as a result (keep as the A126 refutation).

## Gaps before submission (in order)

1. **Brain death: replication.** Confirmatory GTEx test done (A317 H317.1 PASS; A318 adjusted PASS; A319 mouse
   reanalysis supports "rearranged"). Open: an independent human replication (Tabula Sapiens or other organ-donor
   atlases: within-organ cell-type sharing kept, organ order unlike living-brain donors), the ischaemic gradient
   within organ donors (A318 H318.2), time since brain death (not in GTEx). Decides B's high point and venue.
   **Tabula Sapiens scout (2026-11-19): not feasible.** The held Census release (24 donors, 1.14 M cells) has at most
   6 donors per organ with a clock-carrying class (>= 50 cells), 6 donors with >= 4 organs and no living-brain
   reference group, against 79 donors per group in A317; A109b already found the organ arrangement unresolvable here
   (power 0.06-0.13). Donor type (brain death vs circulatory death) is not in the Census metadata (donor summary page
   of the TS portal: cause of death and downtime only). No public multi-organ, same-donor transcriptome of brain-dead
   donors or of experimental brain death in animals was found (searches 2026-11-19; DBD vs DCD datasets are single-organ
   transplant biopsies). An independent human replication needs new data (organ-procurement or autopsy cohorts:
   outward-facing); the mouse central-clock loss (A319) is the available second line of evidence.
2. **Responsive-gene specificity in living tissue.** Closed (A320, 2026-11-19): clock-specific against a calibrated
   broad stress-responsive family (p_cal 0.0045-0.009); the narrow glucocorticoid / immediate-early family stays
   stringent and decidable only with enough pairs (surgical, tumour-normal). Wording in the A320 entry.
3. **Coverage of c (A162 / A163).** Closed as "state" (2026-11-19, `results_preliminary.md` gap 3 note): intervals
   approximate (coverage 0.80-0.98, median 0.94); no claim rests on an absolute c interval.
4. **Gain decomposition, blind.** Predict sample gain from composition x conserved class gains in unseen atlases
   (replaces the closed A298).
5. **Conflicting drives, blind.** One curve of amplitude loss and shift per tissue weight b, on an unseen multi-tissue
   misalignment dataset (Manella liver / lung / fat as the not-blind check first).
6. **Age advance.** Replicate in living subjects or keep in the supplement with its caveats (A130 opposite sign).
7. **Lymphocyte core loop.** Closed (A321, blind, 2026-11-19): no core-loop clock in 3 informative fresh compartments
   (combined p 0.47); "no readable clock" holds with the core loop included; skin T cells an uninformative exception.

## Counts (inventory)

208 entries: PASS 31, partial 66, FAIL 73, descriptive 11, not estimable 14, gate failed 7, not scored 6. Proposed:
about 45 entries carry the main text of B and M together, about 70 go to supplements, about 60 are excluded; FAIL
entries enter mainly as nulls that bound a claim (chronic disease, prognosis, model limits).
