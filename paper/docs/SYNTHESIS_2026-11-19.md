# Synthesis (2026-11-19): what the paper says, what it adds, what is missing, how to present it

Draft for the user and the PI. Builds on `STATUS_2026-11-18.md`, `PAPER_TRIAGE_2026-11-18.md`, the results entries
A310-A319 and the literature cards (`manuscript/literature/notes/`). IDs refer to `results_preliminary.md`. Venue and
the split are the user's and the PI's decision.

## 1. The message in one paragraph

The human body keeps time as a federation. Each organ runs its own clock, built from cell-autonomous clocks in
stromal-vascular, epithelial and parenchymal cells and macrophages (not tissue lymphocytes). A cell class has the same
loudness relative to endothelium in every organ; composition moves a sample's loudness, not its timing, and
composition with these conserved class gains does not explain the range of organ loudness (A323). Organs
express one internal time per person, with gains differing 9-fold (shared variance 85-fold) over a personal,
clock-specific departure of nearly constant size; this departure is why strong clocks look more coupled. The order
between organs does not require alignment to the external day or meals: it is kept in slow deaths and in ventilated
donors with a living brain. In GTEx, brain-dead organ donors carry a different organ order (not dissolved; amplitude
equal to that of equally ventilated donors with a living brain), and brain death cannot be separated there from
collection at procurement. Mouse knockouts of the central clock show the same thing:
a substitute drive (light that bypasses the SCN clock, or feeding) imposes its own organ order, and with no drive
organs drift apart and dampen over days to weeks. Conflicting drives cost amplitude more than phase (night shift:
muscle 2.7 h, excursion to about a third). Chronic disease leaves this organisation intact; psoriasis delays the lesion
clock; tumours keep their patient's time and lose loudness as they proliferate. All of it is read from single,
untimed samples taken in hospital hours, against expression-matched and stress-responsive gene sets.

## 2. What we add to the literature

Legend: **new** = not found in the carded literature or the searches (limited; repeat before submission);
**human first** = known in mouse, first in human tissue; **quantified** = known qualitatively; **qualifies** /
**contradicts** = changes a published reading.

| # | our result (entries) | what was known (cards) | what we add |
|---|---|---|---|
| 1 | One internal time per person across tissues, advancing with elapsed hours; clock-specific (A251-A256) | Donor phases assumed by orderers on full-day GTEx (CHIRAL, Talamanca 2023; CYCLOPS, Ruben 2018); timed multi-site skin (Wu 2018, del Olmo 2022) | **new** as a test: untimed, hospital-hours samples, against matched gene sets |
| 2 | Within organs, time is carried by stromal-vascular, epithelial, parenchymal cells and macrophages; tissue lymphocytes carry no readable clock (A272-A288) | Mouse: all liver populations aligned (Veltri 2025); lymphocyte clocks drive trafficking (Druzd 2017); hepatocytes shape other cells' outputs, not their clocks (Guan 2020) | **human first** at cell-class resolution across organs; **qualifies** the mouse lymphocyte clock ("no readable clock", core loop open, gap 7) |
| 3 | Cell classes share their organ's phase (about 0.6 h); composition moves gain, not phase (A273, A279, A295) | Timed mouse scRNA (aorta, Auerbach 2026; colon atlas); skin cell-type phases (tauFisher, Duan 2024) | **human first**; the composition result is **new** |
| 4 | Each cell class keeps its gain relative to endothelium in every organ (A293 H293.4); composition and conserved class gains do not explain the organ gain range (A323) | none found | **new**: relative loudness is a property of the cell class |
| 5 | Gain over a constant departure: coupling = S / (S + D); gains 9.2-fold apart (shared variance S 85-fold); departure random in direction, not phase-like, not shared by meal-driven tissues (A291, A294b, A302, A311, A312) | Our own earlier observation that stronger clocks look more coupled had no mechanism, and oscillator theory predicts the opposite | **new** explanation; a tested limit: the departure is added variation, not the phase response of a limit cycle (A311 FAIL) |
| 6 | The organ order holds without alignment to the external day (slow deaths) and without meals (ventilated living-brain donors), clock-specifically (A293 Part 1, A109) | Mouse: feeding sets peripheral phases (Damiola 2000; Manella 2021; Greco 2021); human: not tested | **new** in human: the order does not require alignment to the external day, nor meals |
| 7 | In GTEx, brain-dead organ donors carry a different organ order (A317 H317.1 PASS, D_rho 0.621, above all companions; after ischaemic time and RIN adjustment, A318, 0.584; registered re-tests of a post hoc finding, NOT BLIND); aliased with collection at procurement; the shared clock between organs is kept (A317 H317.2 FAIL); amplitude equal to equally ventilated living-brain donors (A315; both ventilated groups below fast deaths, descriptive) | Mouse: SCN lesion loses inter-organ coherence (Yoo 2004; Sinturel 2021); organs persist in vivo at lower amplitude (Tahara 2012); central-clock knockouts keep the order under LD, drift in DD (Husse 2014; Izumo 2014); clinical: brain injury flattens hormone and temperature rhythms | **human first on the order between organs** (single-organ clock-gene differences in brain-dead donor hearts against failing hearts, pseudo-time from collection hour: McTiernan 2021, which reports itself as the first such study); **qualifies** the mouse reading: Izumo's "lost coherence" is mostly a reproducible new order under feeding (A319 H319.2, order share 1.09) |
| 8 | Peripheral clocks are self-sustained and integrate drives with saturation: removing feeding keeps 0.72 of liver amplitude, inverting it moves the liver 11.7 h (A297); the central drive slows feeding-driven shifts | Manella 2021; Saini 2013 (liver re-phases in 1-2 d without the SCN, about 6-7 h after 7 d with it); Kornmann 2007 | **quantified**: additive cues are excluded; one weight b per organ |
| 9 | Night shift: muscle moves 2.7 h after a 12 h behavioural shift, excursion to about a third (A257); forced desynchrony in blood (A262) | Wefers 2018: clock genes "had not adapted"; Archer 2014: loss of rhythmic transcripts | **quantified** (hours and amplitude, clock-specific): conflicting drives cost amplitude more than phase |
| 10 | Fetal organs: a systemically driven clock state before an adult-size loop, carried by signal-induced genes; adrenal cortex, lung, liver, not heart or kidney (A242-A269, A307, A309) | Animal: maternally driven before autonomous; explants rhythmic while in vivo flat (fetal_clocks_animal) | **human first** at organ resolution |
| 11 | Chronic disease (T2D, hypertension, smoking, obesity) leaves couplings, offsets, gain and departure unchanged (A164, A233, A235, A306; powered nulls) | Mouse: diet and obesity rewire peripheral rhythms; human: few multi-organ data | **new** powered nulls in human |
| 12 | Psoriasis delays the lesion clock, clock-specifically against matched genes (A308); five other lesions (periodontitis, osteoarthritis, carotid plaque, fibroid, hidradenitis; one of them skin) do not (A254); against eczema, psoriasis-specificity was not supported (A183 H183.2, not clock-specific) | none found at the clock-state level | **new** |
| 13 | Tumours pair with their patient's clock (A138); COFE's uncoupling is a sampling artefact (A265, A267); proliferating tumours lose gain, not timing (A304); fibroblasts and endothelium keep the patient's clock, immune cells not beyond responsive genes (A227, A271) | COFE (Ananthasubramaniam 2025): tumour-normal uncoupling; TimeTeller; Grabe 2024 (tumour vs own normal, no coupling estimate) | **contradicts** COFE; **new** gain/timing separation in tumours |
| 14 | Age and sex: coupling and offsets equal (A208, A231); CHIRAL's age dampening shrinks with death-type matching (-0.73 to -0.13) | Talamanca 2023: sex-dimorphic, damped with age | **qualifies** Talamanca 2023 (death-type confounding) |
| 15 | Method: untimed, host-referenced coupling, relative amplitude and offset; orderers lose the organ map in 08-18 h deaths, clocksync keeps it (A123, A123b); specificity against matched and responsive sets | CYCLOPS, CHIRAL, COFE, TimeTeller, tauFisher, scRitmo, CoPhaser (synthesis section 2) | **new** estimand and specificity standard; with full-day sampling, two-step pipelines do about as well (A145): the edge is hospital-hours data |

Where the literature corrects us (stated in the paper):
- Mouse tissue amplitude falls after days (Husse: by day 2 in DD) to weeks (Izumo: > 30 days) without any drive,
  through cell dephasing and weaker cell rhythms. A315 is a short-term statement for ventilated, ICU-cued donors.
- Explant phases include a culture reset (Tahara 2012): the mouse support (A319) is from explants and is not blind.
- Weak, permissive coupling between peripheral cells (Finger 2021, TGF-beta; Noguchi 2013, density) exists; the
  federation model keeps it as local and weak, and it is not tested in human (A313 FAIL: the adjacency gradient is not
  clock-specific).

## 3. What is missing

Ranked by how much the headline depends on it.

1. **Brain death: independent replication.** GTEx only (A317, A318). No human cohort found with brain-dead and
   living-brain donors and several organs. Options: (a) Tabula Sapiens or other organ-donor atlases as a within-donor
   check (cell types of an organ still share its phase; organ order differs from GTEx living-brain donors);
   (b) mouse as the second species (A319 done, explants, not blind); (c) per-mouse organ phases from the Takahashi
   lab (outward-facing). Also open: ischaemic time within organ donors runs in the confound's direction (A318 H318.2,
   +0.41, p 0.04); a PubMed check of brain-stem-dead donor endocrinology before claiming priority.
   **Tabula Sapiens scout (2026-11-19): not feasible.** The held Census release (24 donors, 1.14 M cells) has at most
   6 donors per organ with a clock-carrying class (>= 50 cells), 6 donors with >= 4 organs and no living-brain
   reference group, against 79 donors per group in A317; A109b already found the organ arrangement unresolvable here
   (power 0.06-0.13). Donor type (brain death vs circulatory death) is not in the Census metadata (donor summary page
   of the TS portal: cause of death and downtime only). No public multi-organ, same-donor transcriptome of brain-dead
   donors or of experimental brain death in animals was found (searches 2026-11-19; DBD vs DCD datasets are single-organ
   transplant biopsies). An independent human replication needs new data (organ-procurement or autopsy cohorts:
   outward-facing); the mouse central-clock loss (A319) is the available second line of evidence.
   **MODEL prediction 1 (drift timescale after brain death).** The model predicts disorder grows with time since loss
   of the drive; GTEx has ventilation duration, not time since brain death (A132 tertiles in the drift direction, not
   significant). State as a limit, unless a time-since-brain-death variable is found.
2. **Responsive-gene specificity in living tissue: closed (A320).** Clock-specific against a calibrated broad
   stress-responsive family (p_cal 0.0045-0.009); the narrow GC / immediate-early family stays stringent and decidable
   only with enough pairs. Wording in the A320 entry.
3. **Coverage of c: closed as "state".** Intervals approximate (coverage 0.80-0.98, median 0.94); a reliability
   threshold makes coverage worse (gap 3 note); no claim rests on an absolute c interval.
4. **Gain decomposition: done (A323, FAIL).** Composition x conserved class gains spread the organ gains 1.3-fold
   against 2.9-fold observed (8 GTEx tissues, GTEx single-nucleus shares): composition and conserved class gains do
   not explain the organ gain range. Figure 2 and model part 2 reworded. Open: what sets the organ gain (check G1,
   REVIEW_2026-11-19 section 5; biological or technical).
5. **Conflicting drives: done in two mouse datasets (A322, A324).** A pure mix of drives is rejected in both; the
   organs that follow feeding least keep half to two thirds of their amplitude, adipose loses most. Open: a per-organ
   weight predicting a new, partial conflict (no public multi-tissue partial-conflict time series found: jet lag
   designs are liver-only; GSE154797 is fasting / refeeding, i.e. removal).
6. **Age advance** (A234; living cohorts opposite sign): supplement with caveats, or replicate.
7. **Lymphocyte core loop: closed (A321, blind).** No core-loop clock in three informative fresh compartments
   (combined p 0.47); "no readable clock" holds with the core loop included. Skin T cells cohere but are not
   informative (power control 0): a check, not a result. Check L1 (A326): the designs read a lymphocyte clock only
   above about 0.7 of the fibroblast amplitude; a mouse-colon-sized T-cell clock (0.25) would not have been seen, so
   the negative is a bound, not a species difference.
8. **Model document**: done, `MODEL_federation_2026-11-19.md` (with its tested limits A310, A311, A312, A313,
   A315, A319).
9. **Literature before submission**: repeat the inter-organ / inter-tissue coherence search (PubMed, bioRxiv 2026);
   brain-death endocrinology; missing cards (Leproult 2014, Morris 2015 for human misalignment; Janich 2011 for the
   epidermal clock).

Numbering (2026-11-19): gaps 1-7 are numbered as in `PAPER_TRIAGE_2026-11-18.md` and the results entries; the drift
timescale, formerly synthesis gap 2, is MODEL prediction 1.

## 4. How to present it

### One paper (recommended if gap 1 holds) or two

- **One paper, biology first** (Cell): the federation story, with the method's validation in Extended Data and cancer
  as one perturbation figure. The method is the means; the brain-death and cell-class results are what a timed design
  cannot give.
- **Fallback split**: B (biology, figure plan below) and M (method in hospital hours + cancer: A123b, A145, A209,
  A285, A138, A265, A267, A304, A220). The tags in `PAPER_TRIAGE_2026-11-18.md` support either.

### Figure plan (one paper; revised 2026-11-19 after the CHIRAL overlap discussion)

**Principle.** CHIRAL (Talamanca 2023, our lab) already shows that human tissue clocks are synchronised with offsets of
a few hours and that amplitudes rank by tissue. Opening with a GTEx organ map or "one internal time per person" would
repeat its Figure 1 and invite the "reinterpretation only" reading. The paper asks the next question (how that
organisation is built, from cells to organs, and what holds it together) and climbs from the cells to the body. Rules:
- No organ phase map or tissue amplitude ranking as a main result: both go to Extended Data as a concordance check
  with CHIRAL (our offsets vs their tissue phases), credited as its replication.
- Organ-level main panels show only what CHIRAL cannot estimate (it assumes fixed tissue offsets, needs full-day
  sampling, has no cell resolution): conserved class gain, composition effects, the personal departure, and contrasts
  of the order (day, meals, brain death).
- The method is one schematic panel (Fig 1a); its validation goes to Extended Data.
- Figures 1-2 are not GTEx-only (atlases, living surgical cohorts).
- Brain death is the high point (Fig 3), not the opening: it rests on GTEx alone until replicated.

| fig | message | main panels (entries) | Extended Data |
|---|---|---|---|
| 1 | Who keeps time inside human organs | (a) host-referenced readout from untimed samples (schematic); (b) stromal-vascular, epithelial, parenchymal cells and macrophages carry the organ's time, tissue lymphocytes no readable clock (A272, A274, A283, A280, A282, A288); (c) cell classes share the organ phase (A273, A279); (d) one internal time advancing with elapsed hours in living surgical tissue (A251-A256) | per-atlas replicates (A281, A287, A284, A152); hormone and immediate-early states shared, clock not (A275); method validation: orderers fail on 08-18 h deaths, clocksync holds (A123b), simulations (A145), gates (A209), calibration (A285, A258) |
| 2 | Loudness: cell classes keep their relative gain in every organ | conserved class gain across organs (A293 H293.4); composition moves gain within a tissue, not phase (A295); composition and conserved class gains do not explain the organ gain range (A323); gain over a constant personal departure, coupling = S / (S + D) (A291, A294b) | concordance with CHIRAL (organ offsets and amplitude ranking, A123); A302, A311, A312, A313 (limits of the departure) |
| 3 | The order between organs and a living brain | order kept without alignment to the external day and without meals (A293 Part 1, A109); different in brain-dead organ donors, clock-specific, adjusted for ischaemic time and RIN, aliased with collection at procurement (A317, A318; registered re-tests, not blind); amplitude equal to equally ventilated donors (A315); mouse central-clock loss: a substitute drive imposes a new order (A319, with Izumo, Husse, Tahara) | A132, A132b, A314, A317 sensitivities; age and sex equal (A231, A208) |
| 4 | Drives and conflicts | liver self-sustained, saturating integration (A297); night shift: muscle 2.7 h, excursion to about a third (A257); forced desynchrony (A262) | A260, A261, A268, A270; Saini 2013 kinetics |
| 5 | Perturbations | chronic disease nulls (A164, A233, A235, A306); psoriasis delay (A183, A308, A254); tumours keep the patient's time and lose gain (A138, A304; COFE artefact A265, A267) | A202, A180, A115, A227, A271 |
| 6 | Model | the federation: autonomous cellular clocks, class x composition gain, weak local coupling, central drive sets organ order, conflicting drives cost amplitude; tested predictions and failures listed | A310 fit, model grid (A311) |
| ED | Before birth | driven state before an adult-size loop (A243, A244, A266, A269, A307, A309) | negatives and bounds (A236-A241) |

**Introduction sentence.** "CHIRAL established that human tissue clocks are synchronised with offsets of a few hours.
We ask how that organisation is built, from cells to organs, and what holds it together."

### Wording rules (from the results)

- "Population-level estimates of within-person relations, referenced to each donor's own host": never "per person".
- Brain death: "In GTEx, brain-dead organ donors carry a different organ arrangement from every post-mortem group;
  brain death and collection at procurement cannot be separated in this resource"; amplitude "equal to equally
  ventilated donors with a living brain", in the first hours to days; mouse support "explant phases, reanalysis of published summaries".
- Lymphocytes: "no clock readable between donors at depths where fibroblast clocks are read (weaker than about
  0.7 of the fibroblast amplitude)"; never "no clock" or "unlike mouse" (A326).
- The departure: "added variation", not "oscillator stiffness".
- Tumours: "express less of the patient's time", not "lose their clock".
