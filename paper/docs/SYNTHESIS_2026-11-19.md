# Synthesis (2026-11-19): what the paper says, what it adds, what is missing, how to present it

Draft for the user and the PI. Builds on `STATUS_2026-11-18.md`, `PAPER_TRIAGE_2026-11-18.md`, the results entries
A310-A319 and the literature cards (`manuscript/literature/notes/`). IDs refer to `results_preliminary.md`. Venue and
the split are the user's and the PI's decision.

## 1. The message in one paragraph

The human body keeps time as a federation. Each organ runs its own clock, built from cell-autonomous clocks in
stromal-vascular, epithelial and parenchymal cells and macrophages (not tissue lymphocytes). A cell class has the same
loudness relative to endothelium in every organ, and composition changes an organ's loudness, not its timing. Organs
express one internal time per person, with gains differing 85-fold over a personal, clock-specific departure of nearly
constant size; this departure is why strong clocks look more coupled. The order between organs does not come from the
day or from meals: it is kept in slow deaths and in ventilated donors with a living brain, and it is rearranged
(not dissolved, amplitude kept) when the brain dies. Mouse knockouts of the central clock show the same thing:
a substitute drive (light that bypasses the SCN clock, or feeding) imposes its own organ order, and with no drive
organs drift apart and dampen over days to weeks. Conflicting drives cost amplitude more than phase (night shift:
muscle 2.7 h, amplitude to a quarter). Chronic disease leaves this organisation intact; psoriasis delays the lesion
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
| 4 | Each cell class keeps its gain relative to endothelium in every organ (A293 H293.4) | none found | **new**: loudness is a property of the cell class |
| 5 | Gain over a constant departure: coupling = S / (S + D); 85-fold gain range; departure random in direction, not phase-like, not shared by meal-driven tissues (A291, A294b, A302, A311, A312) | Our own earlier observation that stronger clocks look more coupled had no mechanism, and oscillator theory predicts the opposite | **new** explanation; a tested limit: the departure is added variation, not the phase response of a limit cycle (A311 FAIL) |
| 6 | The organ order holds without the day (slow deaths) and without meals (ventilated living-brain donors), clock-specifically (A293 Part 1, A109) | Mouse: feeding sets peripheral phases (Damiola 2000; Manella 2021; Greco 2021); human: not tested | **new** in human: the order does not need the day or meals |
| 7 | Brain death rearranges the organ order (A317 H317.1 PASS, D_rho 0.621, above all companions; survives ischaemic time and RIN adjustment, A318, 0.584); the shared clock between organs is kept (A317 H317.2 FAIL); amplitude kept (A315) | Mouse: SCN lesion loses inter-organ coherence (Yoo 2004; Sinturel 2021); organs persist in vivo at lower amplitude (Tahara 2012); central-clock knockouts keep the order under LD, drift in DD (Husse 2014; Izumo 2014); clinical: brain injury flattens hormone and temperature rhythms | **human first** (no clock-gene or inter-organ study of brain-dead donors found); **qualifies** the mouse reading: Izumo's "lost coherence" is mostly a reproducible new order under feeding (A319 H319.2, order share 1.09) |
| 8 | Peripheral clocks are self-sustained and integrate drives with saturation: removing feeding keeps 0.72 of liver amplitude, inverting it moves the liver 11.7 h (A297); the central drive slows feeding-driven shifts | Manella 2021; Saini 2013 (liver re-phases in 1-2 d without the SCN, about 6-7 h after 7 d with it); Kornmann 2007 | **quantified**: additive cues are excluded; one weight b per organ |
| 9 | Night shift: muscle moves 2.7 h after a 12 h behavioural shift, amplitude to a quarter (A257); forced desynchrony in blood (A262) | Wefers 2018: clock genes "had not adapted"; Archer 2014: loss of rhythmic transcripts | **quantified** (hours and amplitude, clock-specific): conflicting drives cost amplitude more than phase |
| 10 | Fetal organs: a systemically driven clock state before an adult-size loop, carried by signal-induced genes; adrenal cortex, lung, liver, not heart or kidney (A242-A269, A307, A309) | Animal: maternally driven before autonomous; explants rhythmic while in vivo flat (fetal_clocks_animal) | **human first** at organ resolution |
| 11 | Chronic disease (T2D, hypertension, smoking, obesity) leaves couplings, offsets, gain and departure unchanged (A164, A233, A235, A306; powered nulls) | Mouse: diet and obesity rewire peripheral rhythms; human: few multi-organ data | **new** powered nulls in human |
| 12 | Psoriasis delays the lesion clock, clock-specifically; other skin lesions do not (A183, A308, A254) | none found at the clock-state level | **new** |
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
2. **Timescale of the drift after brain death.** The model predicts disorder grows with time since loss of the
   drive; GTEx has ventilation duration, not time since brain death (A132 tertiles in the drift direction, not
   significant). State as a limit, unless a time-since-brain-death variable is found.
3. **Responsive-gene specificity in living tissue** (A285): the living internal-time results stand on matched genes;
   calibrate the responsive rule or state it.
4. **Coverage of c** (A162, A163): intervals undercover at realistic anchor coupling; fix or state as approximate.
5. **Blind gain decomposition** (composition x conserved class gains in unseen atlases): turns findings 3-4 into a
   prediction.
6. **Blind conflicting-drives curve** (weight b per organ): Manella liver, lung and fat as the not-blind check;
   Saini 2013 kinetics as a mouse calibration of b for liver.
7. **Lymphocyte core loop** (A294): blind test before "no clock" loses "readable".
8. **Age advance** (A234; living cohorts opposite sign): supplement with caveats, or replicate.
9. **Model document**: rewrite `MODEL_hopf_2026-11-18.md` as the federation model, with its tested limits
   (A310 inadequate fit, A311 departure not phase-like, A312, A313, A315, A319).
10. **Literature before submission**: repeat the inter-organ / inter-tissue coherence search (PubMed, bioRxiv 2026);
    brain-death endocrinology; missing cards (Leproult 2014, Morris 2015 for human misalignment; Janich 2011 for the
    epidermal clock).

## 4. How to present it

### One paper (recommended if gap 1 holds) or two

- **One paper, biology first** (Cell): the federation story, with the method's validation in Extended Data and cancer
  as one perturbation figure. The method is the means; the brain-death and cell-class results are what a timed design
  cannot give.
- **Fallback split**: B (biology, sections below) and M (method in hospital hours + cancer: A123b, A145, A209,
  A285, A138, A265, A267, A304, A220). The tags in `PAPER_TRIAGE_2026-11-18.md` support either.

### Figure plan (one paper)

| fig | message | main panels (entries) | Extended Data |
|---|---|---|---|
| 1 | Reading the body's clock from untimed hospital samples | host-referenced estimand (c, gain, offset); orderers fail on 08-18 h deaths, clocksync holds (A123b); specificity companions; one internal time per person (A251-A256) | simulations (A145), gates (A209), calibration (A285, A258) |
| 2 | Who keeps time inside organs | stromal-vascular hub vs lymphocytes (A272, A274, A283, A280, A282); macrophages (A288); cell classes share the organ phase (A273, A279) | per-atlas replicates (A281, A287, A284, A152); hormone and immediate-early states shared, clock not (A275) |
| 3 | Loudness is set by cell class and composition | conserved class gain across organs (A293 H293.4); composition moves gain, not phase (A295); 85-fold gain over a constant departure, coupling = S / (S + D) (A291, A294b) | A302, A311, A312, A313 (limits of the departure) |
| 4 | The order between organs needs a living brain | order kept without the day and without meals (A293 Part 1); rearranged in brain death, clock-specific, confound-adjusted (A317, A318); amplitude kept (A315); mouse knockouts: substitute drive, new order (A319, with Izumo, Husse, Tahara) | A132, A132b, A314, A318 sensitivity |
| 5 | Drives and conflicts | liver self-sustained, saturating integration (A297); night shift: muscle 2.7 h, amplitude to a quarter (A257); forced desynchrony (A262) | A260, A261, A268, A270 |
| 6 | Perturbations | chronic disease nulls (A164, A233, A235, A306); psoriasis delay (A183, A308, A254); tumours keep the patient's time and lose gain (A138, A304; COFE artefact A265, A267) | A202, A180, A115, A227, A271 |
| 7 | Model | the federation: autonomous cellular clocks, class x composition gain, weak local coupling, central drive sets organ order, conflicting drives cost amplitude; tested predictions and failures listed | A310 fit, model grid (A311) |
| ED | Before birth | driven state before an adult-size loop (A243, A244, A266, A269, A307, A309) | negatives and bounds (A236-A241) |

### Wording rules (from the results)

- "Population-level estimates of within-person relations, referenced to each donor's own host": never "per person".
- Brain death: "rearranged" and "in GTEx, adjusted for ischaemic time and RNA quality"; amplitude kept "in the first
  hours to days"; mouse support "explant phases, reanalysis of published summaries".
- Lymphocytes: "no readable clock" until gap 7 is closed.
- The departure: "added variation", not "oscillator stiffness".
- Tumours: "express less of the patient's time", not "lose their clock".
