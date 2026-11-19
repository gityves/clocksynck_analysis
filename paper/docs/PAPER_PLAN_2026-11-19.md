# Paper plan (2026-11-19): one time, many volumes

Proposal for the user and the PI. Prompt (user, 2026-11-19), after the discussion of desynchrony as phase and
amplitude: what this means for the paper (figures, analyses done and to perform, novelty), how gain and phase
integrate at organ and cell-type level, where the diseases and the other results go, and (follow-up, same day) how the
development of the clock fits.

Builds on SYNTHESIS_2026-11-19 (section 2, what we add; section 4, figure plan), MODEL_federation_2026-11-19,
MODEL_signal_clock_2026-11-19, REVIEW_2026-11-19, STATUS_2026-11-18 section 2 (development) and PI_QUESTIONS Q25, Q32,
Q40 and Q41. Constraints kept:
- one Cell paper, with the method as the enabler (user, 2026-11-19);
- brain death is not the climax without new data (Q41);
- the claim set stays frozen except for the checks of section 6 (REVIEW 2.4).

Nothing is registered or run. IDs refer to `results_preliminary.md`. Numbers come from the entries; numbers marked
"derived" are conversions of entry values (section 10). If adopted, this replaces the figure plans of SYNTHESIS section 4
and MODEL_signal_clock section 7. The agent prompt for section 6 is `prompts/two_axis_checks_2026-11-19.md`.

## 1. The message

Human tissues keep one time at many volumes. Every clock-carrying cell of an organ keeps the same time, set for the
whole organ by the drives it receives from the brain and from behaviour. How loudly a sample expresses that time is set
separately, by the person, the organ and the cell class. The tissue clock is slow: it carries the person's time
through the acute states of surgery and dying that reset the fast signal layer every cell reads, which is why single
untimed hospital samples can be read. The clock layer is built around birth. Before birth, adrenal cortex, lung and
liver carry the signal layer in the clock's pattern without an adult-size clock. Heart, kidney and brain show no
clock-shaped covariation, and the clock is present in postnatal tissue. Perturbations then fall on two axes:
- **loudness moves often:** body-wide in intensive care, locally in tumours as they proliferate and in chronic lesions,
  and under conflicting drives;
- **time moves rarely:** locally in inflamed skin, partly after a behavioural shift, and between organs only in GTEx
  brain-dead donors (aliased with procurement);
- **the personal departure** from the shared state did not grow where tested (tumours, chronic disease), and chronic
  systemic disease moved neither axis.

What the field calls desynchrony is, in human tissue, mostly a loss of loudness, and the usual synchrony metrics read
loudness.

**Title candidates.**
- "Human tissues keep one circadian time at many volumes"
- "Human tissue clocks lose loudness, not time"
- "A slow clock keeps one time across human tissues"

**The question it raises** (Cell's criterion). If organ time is this robust, and daily programs follow the signals
rather than the clock (A261; mouse), is the clock's loudness, not its phase, the variable that physiology and disease
act on? Development sharpens it: organs run on signals before they have a clock, and lymphocytes all their life.

## 2. Integration: the person's clock plane

### 2.1 One coordinate system for organs and cell classes

Each compartment (an organ's bulk sample or a cell class's pseudobulk) is a point in its person's clock plane:
loudness as the radius, time relative to the person's internal time as the angle. clocksync estimates the point from
untimed samples (relative amplitude |z| and offset arg z against the host). For person p, organ o, cell class k:

    cell class:  z_pok = A_p * G_o * a_k * exp(i phi_o) + delta_pok
    bulk:        z_po  = sum_k s_ok z_pok = A_p * G_o * (sum_k s_ok a_k) * exp(i phi_o) + delta_po

- **Time has one index, the organ (phi_o).**
  - Clock-carrying classes share it: untimed offsets match known time in mouse colon and aorta (A273, A279); human
    classes show no conserved order (A293 H293.3).
  - Composition does not move it (A295, within +-0.13 h per SD).
  - It is set by the organ's mix of drives (A293 Part 1; mouse A297, A319, A322, A324).
- **Loudness has three factors.**
  - Person A_p: persons differ strongly, and every tissue follows the person's amplitude nearly in proportion (A312,
    elasticities 0.83 to 1.77). In living women, fat and skin share it at the size the deaths predict (A331: leverage
    0.177, beyond matched genes, not beyond stress-responsive genes).
  - Cell class a_k: conserved relative to endothelium across organs within a factor of about 1.4 (A293 H293.4);
    lymphocytes below 0.7 of fibroblasts (A272-A288, A321, A326).
  - Organ G_o: gains 9.2-fold apart (A294b). Composition and class gains do not explain them (A323: 1.29-fold
    predicted, 2.89-fold observed). Biological or technical: open; G1 is not decidable with the baboon atlas (A330:
    power 0.60 at a true organ-level Spearman of 0.6).
- **The departure has one size and no direction.**
  - It is clock-specific, about 2.6 times that of expression-matched genes and nearly constant across organs (A294b).
  - It is isotropic (A311) and hardly shared between organs (A302, A313).
  - Coupling c = S / (S + 2 psi) therefore ranks compartments by loudness (A291). A quiet organ's time is read less
    precisely; it is not kept less well.
- **The plane is the slow clock layer.**
  - The fast signal layer (glucocorticoid, feeding and immediate-early states) is read by every cell, lymphocytes and
    fetal cells included (A275, A279, A243-A269).
  - It enters the clock genes through PER1, PER2 and NR1D1; check N2 tests whether this is the departure.
- **Loudness is built in development: a_k depends on the stage (section 2.4).**
  - Before birth, the clock layer is below the sizes the negatives exclude, while the signal layer is present.
  - After birth, the clock is present.

### 2.2 One visual grammar for every figure

Polar plot: angle = time relative to the person's internal time; radius = loudness.
- **Health is points.**
  - Organs spread in radius at offsets of a few hours.
  - The classes of an organ sit at its angle, at different radii; lymphocytes near the origin.
  - Every point carries a departure cloud of the same size, which is why quiet organs look "less coupled".
- **Development is the radius growing from near the origin** (signal layer only) to adult loudness after birth.
- **Perturbations are arrows.** Radial: loudness. Tangential: time. None: nothing moved.
- **One number, if needed.** Delta^2 = (1 - R)^2 + 4 R sin^2(phi / 2) for an amplitude ratio R and a shift phi. In
  clocksync terms, Delta = |z_cond - z_ref| / |z_ref|. The first term is loudness; the second is time, weighted by R.

### 2.3 Where the health results sit

| level | time (angle) | loudness (radius) | departure |
|---|---|---|---|
| cell class | shares the organ's time (A273, A279); no conserved class order (A293 H293.3) | conserved class gains (A293 H293.4); lymphocytes no readable clock, bound 0.7 (A272-A288, A321, A326) | - |
| composition | does not move time (A295) | stromal-vascular content raises gain (A295, descriptive) | - |
| organ | order held without alignment to the day and without meals (A293 Part 1); concordant with CHIRAL (A123) | 9.2-fold (A294b); not composition (A323) | one size, isotropic, not shared (A294b, A311, A313) |
| person | internal time advances with elapsed hours in living tissue (A252, A256) and reads the hour of death (A270) | shared by all tissues (A312) | - |
| development | before birth, the signal-induced genes keep the clock's phase pattern (A307) | below adult size before birth, present after (section 2.4) | - |

### 2.4 Development: the signal layer before the clock

| stage and organ | signal layer | clock layer | entries | status |
|---|---|---|---|---|
| fetal adrenal cortex, lung, liver | clock-shaped covariation between fetuses at half to three quarters of the adult size (adrenal: A242 whole gland, A244 cortex before 12 wpc only; lung, liver: A243); carried by PER1, PER2, NR1D1: not clock-specific on the core loop (A264, A266, A269); these genes keep the clock's phase pattern with the core genes (A307 H307.2, adrenal cortex and lung) and follow glucocorticoid activity in lung (H307.1) | a core loop of adult-like size excluded down to 0.75 of the adult size in adrenal cortex and 0.5 in lung (pooled; lung before 12 wpc is clock-specific on the core loop, descriptive, A269); liver core loop covaries (C 0.686) but not beyond responsive genes: driven or entrained, not separable (A269); a smaller entrained loop not excluded | A242-A244, A264, A266, A269, A307 | not blind; no sampling times |
| fetal adrenal medulla, heart, kidney, brain | none shown | none: a clock of 0.75 to 1 times the adult one excluded (A240); brain 13-22 pcw C 0.047 | A237-A240, A243, A244 | - |
| postnatal brain, heart, kidney, liver | - | present: brain C 0.735 in infancy and childhood, 0.663 from 13 years (A239); postnatal controls C 0.54 to 0.84 (A237) | A237, A239 | - |
| mouse late embryo (E15-E18.75) | signal-induced genes keep the waveform pattern in hepatocytes (A309 H309.2) | hepatocytes clock-specific (A309 H309.1); cardiomyocytes (three quarters of adult size excluded), nephron, lung and adrenocortical progenitors not | A309, A236 | not blind |

**Reading.**
- The systemic drive reaches fetal organs before a core loop of adult size is detectable.
- The clock's loudness is built around birth: in the brain the covariation appears in infancy.
- Liver comes first in both species. Heart and kidney are negative before birth in both.
- Open: is the fetal state a small entrained loop, or an input moving the whole network in the clock's pattern?
  Spliced and unspliced counts decide it (check S1, STATUS_2026-11-18 section 2).
- Precedent (animal): maternally driven before autonomous; explants rhythmic while flat in vivo (fetal_clocks_animal).
  New here: organ and cell-class resolution in human, and the postnatal appearance.

## 3. Perturbations and disease on the two axes

| condition (reference) | time | loudness | departure | entries | status |
|---|---|---|---|---|---|
| intensive care, living brain (fast deaths) | order kept (Fisher-Lee rho 0.690) | lower: clock-specific log variance ratio 0.404, amplitude about 0.82 (derived); within organ donors, lower with ventilation time (descriptive) | - | A293, A315 | registered re-tests, not blind |
| brain death (ventilated, living brain) | order changed (D_rho 0.621; 0.584 adjusted); shared clock state kept; per organ (A332 D2): lung -1.58 h, tibial artery +1.09 h, tibial nerve +0.99 h, others within 1 h | equal: clock-specific ratio 1.05 (amplitude 1.02, derived); every organ 0.92-1.13 of matched (A332 D2) | - | A315, A317, A318, A325, A327, A332 | aliased with procurement; not blind |
| human muscle, day 3 of a 12 h behavioural shift (aligned) | +2.68 h | excursion about 0.31 (derived; projection 0.233; correction note under A257) | - | A257 | not blind; 11 men |
| mouse, daytime feeding (night feeding) | liver 11.3-11.7 h; lung, heart, kidney 1.6-2.5 h; fat 10.60 h (Manella), -0.48 h (Xin) | liver 0.44-0.71; lung, heart, kidney 0.50-0.65; fat 0.32 and 0.11 | - | A322, A324 | A324 blind for heart, kidney, fat |
| mouse liver, feeding rhythm removed | +1.57 h | 0.724 | - | A297 | not blind |
| tumours, high vs low proliferation (own normal tissue) | kept: tumours pair with their patient (A138); COFE's uncoupling is an artefact (A265, A267) | gain 0.63, clock-specific | not raised (power 1.000 at a doubling) | A138, A304 | A304 not blind for A138 |
| tumour cell classes | fibroblasts and endothelium carry the patient's time; immune cells not beyond responsive genes | - | - | A227, A271 | - |
| chronic lesions: periodontitis, osteoarthritis, carotid plaque, fibroid, hidradenitis (own tissue) | kept: pooled +0.01 h, MDE 0.69 h | lower: variance 0.66 pooled (amplitude about 0.81, derived), clock-specific only pooled, at the line | - | A254 | powered for time |
| psoriasis lesion (own skin) | -1.45 h, clock-specific | higher, not clock-specific (r +0.50, amplitude 1.28; 94th / 90.5th percentile of matched / responsive; A115 lesional vs uninvolved pooled +0.30, 56th / 70th) | - | A183, A308, A332 D1 | not blind |
| psoriasis, uninvolved vs healthy skin | later (descriptive) | - | - | A115 | trait component |
| T2D, obesity, hypertension, smoking | offsets and coupling equal (equivalence +-0.15 for T2D and obesity) | equal (powered) | equal (powered) | A164, A233, A235, A306 | - |
| age, sex | equal | CHIRAL's age dampening shrinks with death-type matching | - | A208, A231 | - |

**Reading.**
- Loudness moves often: intensive care, tumours, chronic lesions, conflicting drives.
- Time moves where the configuration of drives changes: brain death (aliased), inflamed skin locally, and partly
  after a behavioural shift.
- The departure did not grow where tested. Chronic systemic disease moved nothing.

**Model rules (MODEL_signal_clock section 3).**
- Time follows which drives set the organ. Loudness follows how strongly and coherently they agree, and how many
  cells follow them.
- Time where you resolve, loudness where you average. Disagreement among the units a measurement averages (cells, cell
  classes in bulk, the drives converging on an organ, people in a group fit) appears as lost loudness.
- So a bulk loss of loudness can be desynchrony below the resolution. This is classic at the cell level
  (welsh2004bioluminescence, nagoshi2004circadian; Ukai 2007, to card) and only single cells can decide it. Public
  data cannot, in human (A160, A305).

**Delta for the map** (derived; loudness share of Delta^2 in brackets).
- Human: muscle 0.79 (77 %); tumours 0.37 (100 %); intensive care 0.18 (about 100 %); chronic lesions 0.19 (about
  100 %).
- Mouse: liver 1.71 (3 %), lung 0.61 (68 %), visceral fat 0.89 (about 100 %).
- A phase-only index scores intensive care, tumours and visceral fat near zero. A loudness-only index scores brain
  death zero.

## 4. What is new

1. **One time, many volumes.**
   - In human tissue, time is set per organ and shared by its clock-carrying classes, while loudness factorises into
     person, organ and class.
   - Composition moves loudness, not time.
   - Conserved class gains and the composition result were not found in any species (SYNTHESIS section 2, rows 3-4).
   - Person loudness shared across tissues is new in human (A312); no study of a donor's amplitude across tissues was
     found.
   - Searches are limited; repeat before submission (SYNTHESIS gap 9).
2. **Coupling measures loudness.** The departure has one size, so synchrony metrics rank loudness (A291, A294b, A311).
   This is a correction the field can apply to its tumour, ageing and disease "desynchrony" claims.
3. **Human perturbations move loudness, rarely time.**
   - The first sorting of human perturbations on both axes, with the same companions:
     - intensive care against brain death, in one cohort;
     - tumours and lesions lose loudness and keep time;
     - inflamed skin moves time;
     - chronic systemic disease gives powered nulls.
   - This answers the paper's disease question: in human tissue, mostly quieter, rarely desynchronised.
4. **A slow clock behind a fast signal layer, in human tissue.**
   - The clock's inertia is shown against the signal layer in the same samples: living interventions, deaths and a
     behavioural shift (A256, A270, A257, A261).
   - It is why untimed hospital tissue can be read; the companions are the operational separation of the two layers.
   - The concept is the mouse consensus, credited as such. The human demonstration and its use are new.
5. **The human clock layer is built around birth.**
   - Organ by organ: the signal layer comes before the clock in adrenal cortex, lung and liver.
   - Heart, kidney and brain show no clock-shaped covariation before birth; the brain clock appears in infancy.
   - Liver comes first, as in the mouse embryo.
   - The animal literature has the maternally driven fetal state; the human organ and cell-class resolution and the
     postnatal appearance are new.
   - Caveats: not blind, no sampling times, signal or small entrained loop undecided until S1.
6. **Corrections of published claims.**
   - COFE's tumour uncoupling (A265, A267).
   - CHIRAL's age dampening (death type).
   - Prognostic clock signatures reduce to proliferation (A191-A194, A220).
   - Human T2D and obesity "disruption" from cultured cells, against powered in vivo nulls (A164, A306).

**Human firsts (not headlines).**
- Cell-class participation, with the lymphocyte bound.
- The order held without meals while the brain lives.
- The loudness cost of a behavioural shift.

**Credited to mouse and other animals.**
- Two or three layers: kornmann2007system, koronowski2019defining, greenwell2019rhythmic.
- Drives set peripheral phases: damiola2000restricted, saini2013realtime.
- Bulk damping as cell dephasing: welsh2004bioluminescence, nagoshi2004circadian.
- Maternally driven fetal rhythms: fetal_clocks_animal.

## 5. Figures (Cell Article, seven)

| fig | message | main panels (entries) | to the Supplement |
|---|---|---|---|
| 1 | A slow clock, read from untimed hospital tissue | (a) the clock plane and the readout; companions as the separation of clock and signal layers (schematic); (b) living tissue, two biopsies hours apart: the clock rotates with elapsed hours, responsive sets do not (A252, A256); (c) deaths: the clock reads the hour of death, responsive and matched sets do not; baboon as the timed reference (A270); (d) night shift: clock 2.68 h at about 0.31, program 11.65 h at 0.854 (A257, A261) | method validation (A123b, A145, A209, A258, A285, A320); blood forced desynchrony (A262) |
| 2 | Structural cells keep the clock; every cell reads the signals | (a) participation across atlases and living donors (A272, A274, A280, A282, A283, A288); (b) lymphocytes share glucocorticoid and immediate-early states, not the clock (A275), with the bound and the core loop (A326, A321); (c) Bmal1 knockout: a shared signal state without a clock (A279) | per-atlas replicates (A281, A287, A284, A152); cortex unresolved (A280, A283) |
| 3 | Development: the signal layer comes before the clock | (a) before birth, adrenal cortex, lung and liver carry a clock-shaped covariation at half to three quarters of the adult size; medulla, heart, kidney and brain none (A237-A239, A242-A244); (b) it rests on the signal-induced genes: not clock-specific on the core loop (A264, A266, A269); they keep the clock's phase pattern and follow glucocorticoid activity in lung (A307); (c) after birth the clock is present: brain 13-22 pcw C 0.047 against 0.735 in infancy and childhood (A239), postnatal controls (A237), the loudness curve (D4); (d) mouse late embryo: hepatocytes first, heart and kidney not (A309); (e) the time-free test: spliced and unspliced phase relations, and the blind adrenal replication (S1, registered as A328; counts not yet delivered), or stated as the open test; loudness after birth (A332 D4) | negatives and bounds (A236-A241); Cao design not estimable (A181); mother-fetus not resolved (A184); developmental loudness table (T5) |
| 4 | One time, many volumes | (a) classes share the organ's time (A273, A279; A293 H293.3); (b) composition moves loudness, not time (A295); (c) the organ order held without alignment to the day and without meals (A293 Part 1); (d) conserved class gains (A293 H293.4); (e) organ gains over a departure of one size: coupling reads loudness; information in bits (A291, A294b, A311, A313; N2); (f) the person's loudness, shared by every tissue (A312), in living women too (A331); (g) the organ range is not composition (A323); biological or technical open (A330 not decidable) | CHIRAL concordance (A123); order sensitivities (A109, A132, A314); A302; A310 fit; A312 per tissue |
| 5 | Perturbations move loudness, rarely time | (a) the map: every condition of section 3 as an arrow; (b) intensive care lowers loudness and keeps the order, brain death changes the order and keeps loudness, in GTEx, aliased with procurement (A293, A315, A317, A318); (c) drives removed and in conflict, mouse: feeding removed (A297), central clock lost (A319), daytime feeding by organ weight (A322, A324); (d) the departure does not grow (A304, A306) | brain-death checks (A132b, A318 H318.2, A325, A327; Q41 Tasks 2-3, B1, B4); ventilation time (A132, A315; D3) |
| 6 | Disease: chronic disease spares the clock, local disease quiets it, inflamed skin shifts it | (a) chronic systemic disease: equivalence and powered nulls (A164, A233, A235, A306); (b) tumours keep the patient's time, their stroma carries it, gain falls with proliferation, the departure stays flat (A138, A227, A271, A304), and the COFE correction (A265, A267); (c) chronic lesions keep time and lose loudness (A254); (d) psoriasis: the lesion later (A308), eczema alike (A183), uninvolved skin late (A115), its loudness (D1) | per-cancer ranking (A3); prognostic nulls (A191-A194, A198, A220); smoking lung gain (A306) |
| 7 | Model | signal and clock layers; the clock layer built after the signal layer; time from drives, loudness from cells; time where resolved, loudness where averaged; tested parts and failures (A310, A311, A312, A313); predictions (N1, tumour single cells, per-organ drive weights, S1) | equations (MODEL_signal_clock section 3) |

**Changes from SYNTHESIS section 4.**
- Development becomes a main figure (Figure 3). The former figures on the order between organs and on loudness merge
  into Figure 4, so the paper stays at seven.
- Brain death moves from its own figure to one panel of Figure 5, stated with its alias (Q41).
- "Drives and conflicts" splits: inertia goes to Figure 1; drive removal and conflict go to Figure 5.
- "Perturbations" splits: the map goes to Figure 5, disease to Figure 6.
- The opening keeps CHIRAL as the starting point: "CHIRAL established that human tissue clocks are synchronised with
  offsets of a few hours. We ask how that time is kept, from cells to organs and from before birth, and what moves it."

**Other results.**
- Age advance (A234): Supplement, with its caveats.
- Brain-body two domains (A263, A289): Supplement; not clock-specific.
- Adrenal offset (A296-A303): out, read as an effect of death.
- Outcome analyses (A191): Supplement.
- Night-shift PBMC (A262 D2, uninformative): Supplement.
- Mouse embryo timing (A108, A236): Figure 3 Supplement, wording per Q30.

## 6. Analyses

### 6.1 Done

The frozen claim set: every entry cited in sections 2, 3 and 5 is scored.

### 6.2 To perform

**Tier 0: assembly (no new claim; needed for the figures).**
- T1. The two-axis table: every condition of section 3, with intervals, on one scale:
  - amplitude = square root of the variance ratio;
  - A304 gain as reported (a slope);
  - A257 excursion = R / cos s;
  - Delta and its split.
  Each row names its source CSV.
- T2. Polar-plot data: organ gains and offsets (A291 and A104 tables), class radii (A293 Part 3), person spread
  (A312), departure size (A294b).
- T3. The text fixes of section 8.
- T4. The registry table: every amendment with its registration and scoring dates, blind status and verdict
  (REVIEW 2.4).
- T5. The developmental loudness table (section 2.4 in numbers):
  - per organ, cell class and stage: the estimate or the bound relative to the adult size;
  - the layer (signal or clock);
  - the core-loop result;
  - the source CSV.
  From the A237-A244, A264-A269, A307 and A309 outputs.

**Tier 1: decisive checks (register, power first; each settles one wording of the claim).**

| check | question | decides | cost |
|---|---|---|---|
| G1 (REVIEW section 5) | is the organ factor biological? GTEx organ gains against timed clock amplitudes of the same organs in baboon and mouse | Figure 4g: "organ" as a loudness factor, or "open" | hours |
| P1 | is person loudness biology or a scale of death? MuTHER fat and skin (680 women, `living_multitissue/E-TABM-1140`): correlation of the two tissues' clock-plane radii across persons, against matched and responsive sets (one punch biopsy, so the companions decide) | Figure 4f in living people; the A312 caveat | hours |
| N2 (MODEL_signal_clock section 6) | is the departure the signal layer read through PER1, PER2 and NR1D1? A294b on the core loop | model part 6; Figure 4e | one 25 min run |
| N1 (MODEL_signal_clock section 6) | without behaviour, do outputs follow the clock? Intensive-care donors, output program against the clock. Partly answered in the archive: A91 (genes following internal time in fast deaths keep a higher ventilated lock than wall-time genes, H91 +1.23) and A27 D5; read them first, register only what they leave open (the per-tissue retention ratio) | the one human test of the two-layer claim; Figures 1 and 7 | hours; power first (79 donors) |
| S1 (STATUS_2026-11-18 section 2; scout of 2026-11-18) | is the fetal state a running loop? The time-free test on spliced and unspliced counts of the open fetal cohorts, plus the blind replication of the adrenal covariation in del Valle 2023 (unseen) | Figure 3e: "a small entrained loop" or "a driven state"; the adrenal replication | registration now (blind); scoring when the server counts (`prompts/server_mapping_fetal_2026-11-18.md`) arrive |

**Tier 2: descriptive companions (complete the map; not blind; labelled descriptive).**
- D1. Psoriasis and eczema loudness with companions: the A254 H254.3 statistic on GSE154200. Places inflamed skin on
  both axes.
- D2. Brain death per organ: phase change and loudness ratio per organ, from the A315 and A317 outputs. Shows where
  the arrows are; no test.
- D3. Ventilation time: the order by tertile (A132) beside loudness by tertile (A315), as "loudness spent, order held".
- D4. The developmental loudness curve. Per organ and developmental bin (prenatal windows; infancy and childhood;
  adolescence and adult), the clock-plane variance against matched sets, beside C with its adult pass rate. Data:
  Cardoso-Moreira 2019 (A237: seven organs, prenatal and postnatal series) and Kang 2011 (A239: brain). The
  sampling-time spread caveat as in A240. Shows when loudness appears.

**Order.**
1. T1-T5 now.
2. Register S1 at once, so it stays blind whenever the counts arrive.
3. G1, P1 and N2 next (hours; data on disk).
4. N1 if its power scout passes.
5. D1-D4 alongside.

### 6.2b Status of the section 6 checks (2026-11-19, prompt `prompts/two_axis_checks_2026-11-19.md`)

| check | amendment | status | what it changes |
|---|---|---|---|
| T1-T5 assembly | - | done: `exploratory/two_axis_map/` (two-axis table, polar points, development table), `docs/triage/registry_2026-11-19.tsv` (219 entries), text fixes | 2.4: adrenal band from A242; liver core loop not excluded (driven or entrained); lung before 12 wpc core-loop positive (descriptive) |
| S1 fetal running loop | A328 | registered blind; waits for the server counts; timed positive control missing (user's decision) | Figure 3e stays "the open test" until scored |
| G1 organ factor | A330 | not decidable at this design (baboon power 0.60 at true rho 0.6; mouse 0.45 at 8 organs, not downloaded) | Figure 4g: "organ" stays open |
| P1 person loudness, living | A331 | PASS (blind, powered): leverage 0.177 beyond matched genes (98.5th), not beyond responsive (93.5th) | Figure 4f holds in living women, with the responsive caveat |
| N2 departure on the core loop | A329 | partly leakage: excess 2.64 -> 1.28 times matched without PER1, PER2, NR1D1 (H329.1); not the GC / IEG signal state (H329.3 FAIL); coupling falls (H329.2 FAIL) | model part 6: departure carried mainly through signal-inducible clock genes, cause open |
| N1 output retention | A333 | registered (staged); power stage not run; archive A91 answers partly (pooled ventilated lock of clock-following genes above 1 in 23 of 26 tissues, liver highest) | Figures 1 and 7: pending |
| D1 skin loudness | A332 | descriptive: psoriasis variance up, not beyond companions | 3, Figure 6d: loudness change not clock-specific |
| D2 brain death per organ | A332 | descriptive: no organ loses loudness; lung, tibial artery and nerve move about 1 h | Figure 5b arrows |
| D3 ventilation time | A332 | organ-donor part only (loudness falls, order already lost); post-mortem tertiles not run | Figure 5 supplement |
| D4 developmental curve | A332 | descriptive: loudness beyond matched after birth in brain and kidney | Figure 3c; V read with C |

### 6.3 Not to do (closed, not powered, or outside the freeze)

- Psoriasis within keratinocytes (A160: MDE above 1.7 h).
- Tumour cell clocks, weaker or fewer (A305: not powered).
- Brain-death replication in public data (Tabula Sapiens scout, SYNTHESIS gap 1).
- Paediatric multi-compartment designs (the D3 scout found none; A181 reading).
- Fetal cohorts behind dbGaP or EGA (Cao 2020, Kildisiute, Sountoulidis): access requests are outward-facing.
- New disease atlases.

### 6.4 Outward-facing (PI; drafts only on request)

- Procurement biopsies after brain death and after circulatory death (Q41 option 1).
- The in-house jet-lag data, as a blind test of per-organ drive weights (Q41 option 3).
- Human single-cell tissue under conflicting drives (night shift), to test whether lost loudness is cell dephasing.
- Access to closed fetal cohorts (Q32), if S1 needs more fetuses than the open ones give.

## 7. Risks at review in this frame

| risk | answer in the paper |
|---|---|
| "time is robust" read as alignment to the day | the untimed design reads relative time between compartments; alignment to the external day cancels into the host, and the paper says so |
| loudness has non-clock causes (composition, attenuation, death) | every loudness claim against matched and responsive sets and referenced to the person; matched genes also lose variance (A254, A315), and the clock is read beyond them |
| brain death | one panel, with its alias and the checks (Q41) |
| GTEx deaths carry Figures 4 and 5 | living cohorts carry Figures 1-2, the lesions and P1; registry table for blind status |
| fetal results: not blind, no sampling times, signal or small loop undecided | state the bounds (A240) and call it a driven state; S1 is the blind, time-free test |
| the model is structural | A310 failed as a joint fit, and the paper says so; each part has its own test |
| human dynamics are not timed | time constants from mouse; one human point (A257) |

## 8. Text fixes now

- **A257 excursion.** R = 0.233 is the misaligned change projected on the aligned chord. The excursion ratio is
  R / cos s, about 0.31 (derived). "Amplitude to a quarter" should read "excursion to about a third (a 2.7 h rotation
  alone would leave 0.77)".
  - To fix: the A257 headline in `results_preliminary.md`, SYNTHESIS lines 21, 41 and 136, and MODEL_federation
    line 63.
  - MODEL_signal_clock was corrected on 2026-11-19.
- **Wording.**
  - "Loudness" (or "gain") for the clock's amplitude relative to the person.
  - "Time" for relative timing between compartments.
  - "Desynchrony" only for phase relations.
  - "Quieter", not "disrupted", for a loss of loudness.
  - "Tumours express less of the patient's time" stays (SYNTHESIS wording rules).
  - Before birth: "a driven state carried by signal-induced genes", not "a fetal clock". Also "an adult-size core loop
    excluded down to x of the adult size; a smaller entrained loop is not excluded".

## 9. Decisions for the user and the PI

1. "One time, many volumes" as the message, and the clock plane as the grammar of the figures.
2. The seven-figure plan (section 5), with development as Figure 3 and time and loudness merged into Figure 4.
3. Which checks to register: S1 now (blind); G1, P1 and N2 next; N1 after its power scout; D1-D4 as descriptive.
4. Title.
5. Q25 recut: the two kinds of desynchrony become time (phase, rare) and loudness (common).

## 10. Conversions

- Variance to amplitude: square root of exp(log ratio) (A315, A254, A183).
- A304 gain: a slope (|C| / V_N), on the amplitude scale.
- A257: excursion ratio = R / cos s = 0.233 / cos 40.2 deg = 0.31.
- Delta^2 = (1 - R)^2 + 4 R sin^2(phi / 2), with phi in radians (hours x 2 pi / 24).
- Intensive care: exp(-0.404 / 2) = 0.82. Brain death: exp(0.046 / 2) = 1.02. Lesions: exp(-0.417 / 2) = 0.81.
