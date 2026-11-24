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

**Wording fixes (2026-11-20).** The independent review's wording fixes (`reviews/plan_review_2026-11-19.md`, sections 2
and 9) and the A335 results are applied in sections 1-5, 6.2b, 7, 8 and 10. Structural recommendations wait in section
9; the follow-up checks are `prompts/two_readouts_2026-11-20.md`. Steps that need data or access not available now:
`DATA_ROADMAP_2026-11-20.md`.

**Revision (2026-11-23).** After the story audit (`STORY_AUDIT_2026-11-23.md`) and the discussion with the user, the
wording fixes are applied in place (sections 1-4 and 10), and section 11 records the framing changes: the mechanism,
the mouse-predicted desynchrony found absent, the positive-control figure, the record of blind tests, and the
mechanism additions. The abstract is revised in `ABSTRACT_DRAFT_2026-11-22.md`.

## 1. The message

At the resolution of tissues and cell classes, human tissues keep one time at many volumes. The clock-carrying classes
of an organ share its time within about an hour; in mouse, that time is set for the whole organ by the drives it
receives from the brain and from behaviour. How loudly a sample expresses that time is set separately, by the organ and
the cell class; persons differ in a body-wide expression scale that the clock follows more reliably than random
genes in GTEx deaths, instantaneous deaths included (A336, A367), not shown as a clock property in living tissue
(A366, A368; revised 2026-11-24). Evidence base (2026-11-24): the organ level (organ gains, order, coupling, the person
scale) is shown in post-mortem donors only, defended against post-mortem attenuation (A356) and by its structure (a
body-wide death process cannot create differences between the organs of one person); the cell level is also shown in
living tissue (participation ranking in five living atlases, A282); in living people, a behavioural shift costs the
muscle clock loudness (A257). Living two-biopsy and two-organ results (A251, A252, A255, A256) are positive controls of
the method, not findings about organisation. The tissue clock is slower than the signal layer: it carries the person's time through the acute states of
surgery and dying that reset the fast signal layer every cell reads, which is why single untimed hospital samples can
be read. Before birth (revised 2026-11-23), the adrenal already carries a clock-like organisation: its core loop
covaries between fetuses on its own, replicated blind (A328 part a), and its states turn in two dimensions (A340).
Liver and lung carry the signal layer in the clock's pattern, not shown to be a clock, and not explained by
glucocorticoid activity alone (A269, A307, A341). Heart, kidney and adrenal medulla show no clock-shaped covariation;
brain shows none at 13-22 pcw but covaries before 12 wpc in one cohort (descriptive). The
clock is present in postnatal tissue. Perturbations then fall on two axes:
- **loudness moves often:** body-wide in intensive care, locally in tumours as they proliferate and in chronic lesions,
  and under conflicting drives;
- **time moves rarely:** partly after a behavioural shift (the core loop by about 4 h), perhaps locally in inflamed skin
  (at the line on the core loop), and between organs only in GTEx brain-dead donors (aliased with procurement and donor
  management);
- **the personal departure** from the shared state did not grow where tested (tumours, chronic disease), and chronic
  systemic disease moved neither axis.

What the field calls desynchrony is, at tissue resolution, mostly a loss of loudness, and the usual synchrony metrics
read mostly loudness. Below that resolution a loss of loudness can itself be desynchrony among cells, which public data
cannot decide (A334).

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
  - Clock-carrying classes share it: untimed offsets match known time in mouse colon and aorta (A273, A279; both hold
    on the core loop, A336). Human classes share it within 1 h, with a small endothelium-minus-fibroblast offset of the
    mouse sign (0.46 h, 90 % 0.09-0.90; A293, A335 equivalence, not blind; core loop 0.29 h, -0.19 to 0.80, four
    readable atlases, A336); the six-atlas sign test failed (H293.3).
  - Composition does not move it (A295, within +-0.13 h per SD).
  - In mouse it is set by the organ's mix of drives (A297, A319, A322, A324); in human only the order is shown (A293
    Part 1).
- **Loudness: organ and class factors, approximately separable; the person factor is shown in deaths.**
  - Person A_p (revised 2026-11-24, `docs/PERSON_VOLUME_2026-11-24.md` section 12): persons differ in a body-wide
    expression scale; in GTEx deaths the clock follows it more reliably than random genes (core half-host leverage
    0.714 against matched q95 0.621 and responsive q95 0.692; A336 part 2), but its size is what random genes share
    (cv2 0.120 against a matched median of 0.071; A336 size note). It is present at full size in violent,
    instantaneous deaths (core 0.726, 34 donors; A367), scales every core clock gene alike (A363) and every tissue
    nearly in proportion (A312), and does not follow the mouse cue-conflict organ pattern (A362). In living people a
    clock-specific excess is not shown: the A331 fat-skin leverage is the shared biopsy time (A366), and at matched
    phase fat and skin share no more than random genes (A368; a GTEx-size clock volume excluded); STAGE and AAGMEx
    are not decidable (A369, A370).
  - Cell class a_k: conserved relative to endothelium across organs within a factor of about 1.4 (A293 H293.4; on the core loop in four
    of five classes, mural cells not, A350);
    lymphocytes do not share the tissue's clock state (A272-A288); an own clock is not excluded below about 0.7 of the
    fibroblasts', a level other resident classes such as myeloid (0.53) also sit below (A281, A321, A326; the bound is
    unchanged on the core loop, A336) (revised 2026-11-23).
  - Organ G_o: gains 9.2-fold apart, 3.1-fold without whole blood and testis (27 peripheral tissues; A104 table,
    derived). Composition and class gains do not explain them (A323: 1.29-fold predicted, 2.89-fold observed). Not
    post-mortem attenuation: within tissues gain does not move with ischaemic time (6.7 h apart: 1.01-fold) or RIN
    (0.99-fold), decidable; ischaemia and RIN could produce at most about 4 % of the organ spread (A356, check G2).
    Reproduction in living tissue not shown: living pairs do not follow the GTEx ratios as a set, omental over
    subcutaneous fat in the GTEx direction at two thirds of its size, not decidable (A355); the baboon ranking test is
    not decidable (A330). In a timed primate organ amplitudes span 2.8-fold (10th-90th percentile, A36 table,
    descriptive) (revised 2026-11-23).
- **The departure varies little and has no dominant direction.**
  - It is clock-specific, about 2.6 times that of expression-matched genes (A294b). It is a coherent shift of the
    tissue's whole clock state, about 2.0 times matched gene-level variation and present in all 27 tissues; it is not
    carried by PER1, PER2 and NR1D1. The core-loop fall to 1.28 times matched (A329) is the 11-gene projector's 3.1-fold
    noise gain (A329 erratum, post hoc; A336 gene-drop null, 32nd percentile) (revised 2026-11-23).
  - A small part of it is a person-level glucocorticoid-like state read unequally by tissues: the signal-gene
    residuals of a person's other tissues explain a median R^2 of 0.027 of a tissue's departure, against 0.008 for
    matched genes (A357, registered; skeletal muscle 0.195); this state tracks the person's glucocorticoid-target
    state (r 0.77, descriptive). Daily cortisol history and an agonal surge are not separated (revised 2026-11-23).
  - Across the 27 peripheral tissues 2 psi spans 4.6-fold, against 85-fold for S (A104 table, derived).
  - It is not phase-dominated (A311) and hardly shared between organs (A302, A313).
  - Coupling c = S / (S + 2 psi) therefore ranks compartments mostly by loudness (A291). Gain is the transfer of the
    body's time to an organ; the departure is added at the organ, after the transfer (slope 0.07, not 1), so coupling
    rises with gain and louder organs stay closer to the body's time in hours (implied 3.0 h against 1.1 h from
    A294b, matching A259): gain and coupling are one organ property seen two ways. What sets the gain is open (not
    composition, A323; not post-mortem attenuation, A356) (revised 2026-11-24).
- **The plane is the slow clock layer.**
  - The fast signal layer (glucocorticoid, feeding and immediate-early states) is read by every cell, lymphocytes and
    fetal cells included (A275, A279, A243-A269).
  - It enters the clock genes through PER1, PER2 and NR1D1, but neither the departure's excess nor the living person
    leverage is carried by them beyond what three random genes of similar weight carry (A336 gene-drop null; revised
    2026-11-23).
  - In single-cell atlases they also carry a disproportionate share of what classes share: dropping them puts the
    cross-class statistics at the 1st-7th percentile of all three-gene drops (A275, A280, A282, A283, A288, A227;
    A336, descriptive; 4th-29th among drops of similar precision). The structural-lymphocyte contrast survives on the
    core loop against expression-matched genes, and against responsive genes in three of six tests.
- **Loudness is built in development: a_k depends on the stage (section 2.4).**
  - Before birth, the clock layer is below the sizes the negatives exclude, while the signal layer is present.
  - After birth, the clock is present.

- **Where coupling to the body's time comes from (revised 2026-11-24; `COUPLING_SOURCES_2026-11-24.md`).**
  - Coupling c_o = S_o / (S_o + D_o) has two sources: the transfer G_o of the body's time to the organ, and the
    departure D_o added after it. D_o is nearly constant, so coupling ranks organs by G_o; the question is what sets
    G_o (and phi_o).
  - Set by the brain: the arrangement, not the amplitude. Order kept without alignment to the day in slow deaths (A293 Part 1; in intensive care only "no large disruption": one organ about 6-7 h, dispersion SD 2.25 h, A349; revised 2026-11-24), changed
    after brain death (A317 / A318), while amplitude equals that of equally ventilated living-brain donors (A315).
  - Set by the body's state: intensive care costs amplitude (A315 / A335, falls with ventilation time); persons differ
    in a body-wide scale that every organ follows, shown in deaths (A312, A336, A367); a person-level glucocorticoid-like state is read unequally by organs and
    explains a small part of the departure (A357).
  - Drives, quantified in mouse: removing a drive keeps amplitude (A297, 0.724); conflicting drives cost it, most where
    an organ's weights are balanced, less than a pure mix (A322, A324); a behavioural shift moves the human muscle clock
    and costs amplitude (A257, A335).
  - Not set by cells or post-mortem attenuation (A323, A356). What sets G_o is open: drive coherence across species
    is a suggestive pattern only (A358, descriptive: organs that weigh one drive in mouse are the loud ones in human,
    adipose the exception; rho -0.6, 5 organs, p 0.175, below the registered -0.7: mention, not a claim); reception of the systemic glucocorticoid state does not track gain (A359
    FAIL; the immediate-early reception slope does, rho 0.44, a lead only); the in-house jet lag remains the
    out-of-sample test of organ drive weights (data needed).
  - In deaths with a terminal phase the person scale falls with the glucocorticoid-target state (A360 PASS, partial
    R^2 0.06; core loop -0.184 beyond all 400 companions, A361); not in instantaneous deaths, where the coefficient is
    +0.22 and the GC-target state lower (A367, descriptive, 34 donors): probably an agonal relation. Not a technical
    scale of RIN or ischaemic time (A360 H360.2 FAIL).

### 2.2 One visual grammar for every figure

Polar plot: angle = time relative to the person's internal time; radius = loudness.
- **Health is points.**
  - Organs spread in radius at offsets of a few hours.
  - The classes of an organ sit at its angle, at different radii; lymphocytes near the origin.
  - Every point carries a departure cloud of the same size, which is why quiet organs look "less coupled".
- **Development is the radius growing from near the origin** (signal layer only) to adult loudness after birth.
- **Perturbations are arrows.** Radial: loudness. Tangential: time. None: nothing moved.
  - A change shared by every organ (the person's overall loudness, as in intensive care) is drawn apart from the
    relative pattern: organ gains relative to their mean, beside organ phases relative to the person.
  - Each condition carries its detection limit as a circle around the origin: an arrow inside it reads "not
    detectable", not "no change" (`prompts/two_readouts_2026-11-20.md`, part 3).
- **One number, if needed.** Delta^2 = (1 - R)^2 + 4 R sin^2(phi / 2) for an amplitude ratio R and a shift phi. In
  clocksync terms, Delta = |z_cond - z_ref| / |z_ref|. The first term is loudness; the second is time, weighted by R.

### 2.3 Where the health results sit

| level | time (angle) | loudness (radius) | departure |
|---|---|---|---|
| cell class | shares the organ's time (A273, A279); human classes within 1 h, endothelium minus fibroblast 0.46 h (A293, A335; not blind) | conserved class gains (A293 H293.4); lymphocytes do not share the tissue's clock state (A272-A288); own clock not excluded below about 0.7 of fibroblasts, where myeloid (0.53) also sits (A281, A321, A326; revised 2026-11-23) | - |
| composition | does not move time (A295) | stromal-vascular content raises gain (A295, descriptive); lymphocyte content dilutes relative gain slightly (-0.020 per SD on the core loop, p 0.040, A336; within the post hoc ±0.10 margin, A335) | - |
| organ | order held without alignment to the day in slow deaths; in intensive care no large disruption at the ventilated detection limit (one organ about 6-7 h, A349) with offsets compressed toward 0 / 12 h (A132b) (A293 Part 1); concordant with CHIRAL (A123) | 9.2-fold, 3.1-fold without blood and testis (A104 table, derived); not composition (A323) | varies little, not phase-dominated, not shared (A294b, A311, A313); mostly PER1, PER2, NR1D1 (A329) |
| person | internal time advances with elapsed hours in living tissue (A252, A256) and reads the hour of death (A270; held-out waveform, A335) | shared by all tissues in deaths (A312); in living tissue not on the core loop (A331, A335) | - |
| development | before birth, liver and lung: the signal layer in the clock's pattern (A269, A307) | adrenal before birth: a clock-like organisation, core loop clock-specific (A328 part a, A340); other organs below adult size before birth, present after (section 2.4) | - |

### 2.4 Development: a clock-like adrenal, the signal layer elsewhere (revised 2026-11-23)

| stage and organ | signal layer | clock layer | entries | status |
|---|---|---|---|---|
| fetal adrenal | clock-shaped covariation of the 14 genes at about two thirds of the adult size (A244 cortex cells, C 0.61; A242 whole gland), replicated blind before 12 wpc (del Valle, C 0.59, A328 part a H328.1) | core loop clock-specific on its own in del Valle (C 0.80, A328 part a H328.2); two-dimensional states, not one slow factor (A340 H340.1); fast pulsatile input excluded (post hoc, A328 part b, A338, A339 at tau = 1 h); single-cell cortex core loop 0.42, three quarters of the adult size excluded, half not (A266): the datasets differ, cause open | A242, A244, A264, A266, A328, A340 | del Valle blind; no sampling times |
| fetal lung, liver | clock-shaped covariation of the 14 genes (lung 0.53, liver 0.58; A243 pooled; before 12 wpc lung passes, C 0.679, liver does not, C 0.169); PER1, PER2, NR1D1 keep the clock's phase pattern with the core genes (A307 H307.2, lung) and follow glucocorticoid activity in lung (H307.1); the covariation survives adjustment for glucocorticoid-target activity (A269 H269.2); glucocorticoid responsiveness does not rank the organs (A341) | lung: core loop small (0.47; adult size excluded down to 0.5, A269); liver: core loop covaries (0.69) but not beyond responsive genes, states in the no-clock elongation range with the main axis along the glucocorticoid-target direction (A340, descriptive): driven or entrained, not separable | A243, A269, A307, A340, A341 | not blind; no sampling times |
| fetal adrenal medulla, heart, kidney, brain | none shown, except brain before 12 wpc | none: a clock of 0.75 to 1 times the adult one excluded (A240; heart also at half, A243); brain 13-22 pcw C 0.047 (A239); hindbrain and forebrain before 12 wpc covary, C 0.578 and 0.465 (p 0.016, 0.027; Cardoso-Moreira, D4, descriptive) | A237-A240, A243, A244, A332 D4 | the negatives exclude 0.75 to 1 times the adult size; the positives sit at 0.5-0.75 |
| postnatal brain, heart, kidney, liver | - | present: brain C 0.735 in infancy and childhood, 0.663 from 13 years (A239); postnatal controls C 0.54 to 0.84 (A237); Cardoso-Moreira birth-12 y: hindbrain 0.772, forebrain 0.044 (D4) | A237, A239, A332 D4 | - |
| mouse late embryo (E15-E18.75) | signal-induced genes keep the waveform pattern in hepatocytes (A309 H309.2) | hepatocytes clock-specific (A309 H309.1); cardiomyocytes (three quarters of adult size excluded), nephron, lung and adrenocortical progenitors not | A309, A236 | not blind |

**Reading.**
- Before birth the adrenal carries a clock-like organisation (core loop, blind); in liver and lung a systemic input
  moves the clock genes in the clock's pattern, glucocorticoids a candidate, not shown.
- The clock's loudness is built around birth. In brain the covariation appears in infancy in Kang 2011 (A239); in
  Cardoso-Moreira forebrain it is present before 12 wpc and absent at birth-12 y (D4): the cohorts disagree.
- In mouse, hepatocytes come first (A309); in human, lung is shown before 12 wpc and liver is not (A243). Heart and
  kidney are negative before birth in both species.
- Open: does the adrenal clock run on its own or follow a daily maternal / ACTH signal? The spliced / unspliced lead
  did not decide it (A328 part b, A338, A339: uninformative, the running postnatal liver is not resolved either).
- Precedent (animal): maternally driven before autonomous; explants rhythmic while flat in vivo (fetal_clocks_animal).
  New here: organ and cell-class resolution in human, and the postnatal appearance.

## 3. Perturbations and disease on the two axes

| condition (reference) | time | loudness | departure | entries | status |
|---|---|---|---|---|---|
| intensive care, living brain (fast deaths) | no large disruption of the order (Fisher-Lee rho 0.690; registered re-test, A293, of the descriptive A132 readout after H109.2 failed); offsets compressed toward 0 / 12 h (0.58 of fast deaths, A132b), read as a shared non-temporal state on running clocks (A109); detection limit at the ventilated design: one organ about 6-7 h (5.6 h core), dispersion SD 2.25 h (A349): "no large disruption" | lower: clock-specific log variance ratio 0.404, loudness about 0.82 (derived); holds on the core loop (+0.457) and matched on age and sex (+0.459), with raw clock variance higher in fast deaths while matched genes go the other way (A335); within organ donors, lower with ventilation time (descriptive) | - | A293, A315, A335 | registered re-tests, not blind; core loop and matching blind (A335) |
| brain death (ventilated, living brain) | order changed (D_rho 0.621; 0.584 adjusted); shared clock state kept; per organ (A332 D2): lung -1.58 h, tibial artery +1.09 h, tibial nerve +0.99 h (90 % intervals exclude 0); whole blood -1.22 h (-3.73 to +0.78); the other eight within 0.8 h. On the core loop artery (+1.20 h) and nerve (+1.54 h) hold, lung does not, lower-leg skin (+1.33 h) appears (A335). 3 of 12 intervals exclude 0 on either gene set, against 1.2 expected by chance: descriptive | equal relative to expression-matched genes: clock-specific ratio 1.05 (loudness 1.02, derived), equivalent within +-0.357 on 14 genes and the core loop, pooled and matched (A335); raw clock variance lower in organ donors, as for matched genes (core Q -0.186, matched median -0.181); every organ 0.92-1.13 of matched (A332 D2) | - | A315, A317, A318, A325, A327, A332, A335 | aliased with procurement and donor management; not blind |
| human muscle, day 3 of a 12 h behavioural shift (aligned) | core loop +4.06 h (90 % 2.57-5.68; all 11 subjects the same sign); 14 genes +2.68 h | core loop: excursion about 0.49 (derived; projection 0.239), per subject 0.34-0.92 (median 0.47, all 11 below 1); 14 genes: about 0.31 (derived; projection 0.233; correction note under A257) | - | A257, A335 | 14 genes not blind; core loop blind (A335); 11 men |
| mouse, daytime feeding (night feeding) | liver 11.3-11.7 h; lung, heart, kidney 1.6-2.5 h; fat 10.60 h (Manella), -0.48 h (Xin) | liver 0.44-0.71; lung, heart, kidney 0.50-0.65; fat 0.32 and 0.11 | - | A322, A324 | A324 blind for heart, kidney, fat |
| mouse liver, feeding rhythm removed | +1.57 h | 0.724 | - | A297 | not blind |
| tumours, high vs low proliferation (own normal tissue) | kept: tumours pair with their patient, at the line on the core loop (0.252 against matched q95 0.275, A336) (A138; pairing is blind to a constant tumour-normal offset); COFE's uncoupling is an artefact (A265, A267) | gain 0.63, clock-specific; within purity halves Q_G -0.266 (gain about 0.77, derived; A304 sensitivity) | not raised (power 1.000 at a doubling) | A138, A304 | A304 not blind for A138 |
| tumour cell classes | fibroblasts keep the patient's normal fibroblast clock state on the core loop, beyond matched and broad responsive genes (A346; 0.686 → 0.683); endothelium pairs but is not beyond matched genes on the core loop (A346); malignant epithelium not clock-specific (A346 D346); no independent cohort (A226); immune cells pair with normal tissue, not beyond responsive genes, and on the core loop not beyond matched genes (A227 at the line, A336) | - | - | A227, A271, A336 | - |
| chronic lesions: periodontitis, osteoarthritis, carotid plaque, fibroid, hidradenitis (own tissue) | kept: pooled +0.01 h, MDE 0.69 h | lower: variance 0.66 pooled (loudness about 0.81, derived), clock-specific only pooled, at the line | - | A254 | powered for time |
| psoriasis lesion (own skin) | 14 genes -1.45 h, q 0.016 (A308 rule, adopted after A183's rule failed); core loop -1.00 h, q 0.049, at the line (5.5th percentile of matched sets, as on 14 genes; A335): suggestive | higher, not clock-specific (r +0.50, loudness 1.28; 94th / 90.5th percentile of matched / responsive; A115 lesional vs uninvolved pooled +0.30, 56th / 70th) | - | A183, A308, A332 D1, A335 | not blind; core loop blind (A335) |
| psoriasis, uninvolved vs healthy skin | later (descriptive) | - | - | A115 | trait component |
| T2D, obesity, hypertension, smoking | offsets and coupling equal (equivalence +-0.15 for T2D and obesity; A233 and A235 hold on the core loop, A336) | equal (powered) | equal (powered) | A164, A233, A235, A306 | - |
| age, sex | equal | equal: clock 1.038 (0.995-1.083), outputs 1.030 (0.965-1.078), old / young after death-type matching (A373, equivalence) | - | A208, A231, A373 | - |

**Reading.**
- Loudness moves often: intensive care, tumours, chronic lesions, conflicting drives.
- Time moves where the configuration of drives changes: partly after a behavioural shift (the core loop by about
  4 h), between organs in brain death (aliased), and perhaps locally in inflamed skin (at the line).
- The departure did not grow where tested. Chronic systemic disease moved nothing.

**Model rules (MODEL_signal_clock section 3).**
- Time follows which drives set the organ. Loudness follows how strongly and coherently they agree, and how many
  cells follow them.
- Time where you resolve, loudness where you average. Disagreement among the units a measurement averages (cells, cell
  classes in bulk, the drives converging on an organ, people in a group fit) appears as lost loudness.
- So a bulk loss of loudness can be desynchrony below the resolution. This is classic at the cell level
  (welsh2004bioluminescence, nagoshi2004circadian; Ukai 2007, to card) and only single cells can decide it. Public
  data cannot, in human (A160, A305), nor in timed mouse 10x tissue (A334, A334b).

**Delta for the map** (derived; loudness share of Delta^2 in brackets).
- Human: muscle 0.79 (77 %) on 14 genes and 0.87 (34 %) on the core loop, where the night shift is mostly a time
  move (A335); tumours 0.37 (100 %); intensive care 0.18 (about 100 %); chronic lesions 0.19 (about 100 %).
- Mouse: liver 1.71 (3 %), lung 0.61 (68 %), visceral fat 0.89 (about 100 %).
- A phase-only index scores intensive care, tumours and visceral fat near zero. A loudness-only index scores brain
  death zero.

## 4. What is new

1. **One time, many volumes.**
   - In human tissue, at tissue and cell-class resolution, organs keep their order and their clock-carrying classes
     share the organ's time within about an hour; loudness is approximately separable into organ and class factors.
   - Composition does not move time; its effect on loudness is descriptive.
   - Conserved class gains and the composition result were not found in any species (SYNTHESIS section 2, rows 3-4).
   - Person loudness shared across tissues is shown in GTEx deaths, on the core loop and beyond expression-matched and
     stress-responsive genes (A312, A336); in living women it is at the line on the core loop (A331, A335, A336). No
     study of a donor's amplitude across tissues was found.
   - Searches are limited; repeat before submission (SYNTHESIS gap 9).
2. **Coupling measures loudness.** The departure varies 4.6-fold where loudness varies 85-fold (derived), so synchrony
   metrics rank mostly loudness (A291, A294b, A311). This is an empirical warning the field can apply to its tumour,
   ageing and disease "desynchrony" claims.
   - Coupling to the body is the transfer of the body's time to the organ (gain) over a departure of fixed size added
     after it: louder organs stay closer to the body's time in hours (3.0 h to 1.1 h implied, matching A259).
2b. **Where coupling comes from, in human (revised 2026-11-24; `COUPLING_SOURCES_2026-11-24.md`).**
   - The brain sets the organ arrangement but not the amplitude (A293, A315, A317); the body's state sets amplitude
     (intensive care, A315 / A335; the person volume, A312).
   - Removing a drive and conflicting drives act differently on amplitude, quantified (A297, A322, A324): the
     literature states it verbally (manella2021liverclock).
   - A person-level glucocorticoid-like state, read unequally by organs, is a small source of person-specific internal
     desynchrony (A357): the human counterpart of balsalobre2000resetting.
   - The organ gain is an organ property beyond cells (A323, A356); the literature ranks amplitudes
     (talamanca2023sexdimorphic) without decomposing them.
3. **Human perturbations move loudness, rarely time.**
   - The first sorting of human perturbations on both axes, with the same companions:
     - intensive care against brain death, in one cohort;
     - tumours and lesions lose loudness and keep time;
     - inflamed skin may move time (at the line on the core loop);
     - chronic systemic disease gives powered nulls.
   - Every "kept" is stated with its detection limit (A336 part 3): in intensive care the order excludes organ phases
     dispersed by SD 1.5 h and a single-organ shift of about 4 h (median over organs; 2.3 h nerve to 9.8 h blood);
     lesions 0.69 h; class time 0.61 h; brain-death loudness x0.90 to x1.11; tumour time only a dispersion >= 2.4 h.
     "Kept" means "not detected at that limit".
   - Brain death's rearrangement of the order holds on the core loop (D_rho 0.465; 0.496 adjusted; A336 part 4):
     it is in the core clock, not leakage through PER1, PER2 and NR1D1; donor-management glucocorticoids acting on the
     core are not excluded (A136).
   - Five chronic-disease loudness and departure nulls lose their equivalence on the core loop (wider intervals),
     with nothing detected (A306, A336).
   - This answers the paper's disease question: in human tissue, at tissue resolution, mostly quieter, rarely
     desynchronised.
4. **A slow clock behind a fast signal layer, in human tissue.**
   - The clock's inertia is shown against the signal layer in the same samples: living interventions, deaths and a
     behavioural shift (A256, A270, A257, A261).
   - On the core loop the behavioural shift moves the clock 4.06 h while PER1, PER2 and NR1D1 move less (A335,
     descriptive); the output program moves 11.65 h (A261).
   - It is why untimed hospital tissue can be read; the companions are the operational separation of the two layers.
   - The concept is the mouse consensus, credited as such. The human demonstration and its use are new.
5. **Before birth, a clock-like adrenal; the signal layer elsewhere (revised 2026-11-23).**
   - The adrenal core loop covaries between fetuses on its own, replicated blind before 12 wpc (A328 part a), with
     two-dimensional states (A340). Lung and liver carry the signal layer in the clock's pattern, not shown to be a
     clock; not explained by glucocorticoid activity alone (A269, A307, A341). Liver passes only when stages are pooled
     (not before 12 wpc; A243).
   - Heart and kidney show no clock-shaped covariation before birth. Brain shows none at 13-22 pcw but covaries before
     12 wpc in one cohort (D4, descriptive); when it appears after birth differs between cohorts (A239, D4).
   - In the mouse embryo hepatocytes come first (A309); in human, lung, not liver, is shown before 12 wpc (A243).
   - The animal literature has the maternally driven fetal state; the human organ and cell-class resolution and the
     postnatal appearance are new.
   - Caveats: no sampling times; the adrenal running on its own or driven daily is undecided (S1 lead not decidable);
     the two adrenal datasets disagree on the core loop (A266 vs A328); the positives sit at 0.5-0.75 of the adult
     size while the negatives exclude 0.75 to 1 (A240), and only heart is powered at half (A243).
6. **Corrections of published claims.**
   - COFE's tumour uncoupling (A265, A267).
   - CHIRAL's age dampening (revised 2026-11-24, A373, registered, blind; PI to see before the text): between 21-49
     and 60-70 years, once death type is balanced, neither the peripheral clock (absolute amplitude ratio old / young
     1.038, 90 % 0.995-1.083) nor the rhythmic output genes on CHIRAL's own DIP (1.030, 0.965-1.078; Whole Blood
     excluded) lose amplitude; equivalence within +-10 % passes for both. CHIRAL's analogue (per-gene cosinor on the
     DIP) gives 0.995 matched and 0.978 unmatched. The published dampening is most likely a count effect of death type
     and noise (its genome-wide rhythmic-gene count shrinks under matching, WP0 item 5, unregistered). Timed mouse data
     show a modest fall at 27 months that is not clock-specific (A372, descriptive).
   - Published clock prognostic signatures add nothing beyond proliferation or random gene sets (A192-A194, A220).
   - Human T2D and obesity "disruption" from cultured cells, against powered in vivo nulls (A164, A306).

**Human firsts (not headlines).**
- Cell-class participation, with the lymphocyte bound.
- No large disruption of the order in ventilated donors while the brain lives (detection limit one organ about 6-7 h,
  A349; offsets compressed toward 0 / 12 h, A132b; feeding unrecorded; maintenance by the brain is not separated from
  inertia, D3 not estimable, A335).
- The loudness cost of a behavioural shift.

**Credited to mouse and other animals.**
- Two or three layers: kornmann2007system, koronowski2019defining, greenwell2019rhythmic.
- Drives set peripheral phases: damiola2000restricted, stokkan2001entrainment, balsalobre2000resetting,
  brown2002rhythms, buhr2010temperature, saini2013realtime; organ-specific weights: manella2021liverclock,
  izumo2014differential; the SCN holds organs together, not cells within an organ: yoo2004periodluciferase,
  sinturel2021circadian, husse2014lightdark; jet lag: yamazaki2000resetting, kiessling2010adrenal.
- Bulk damping as cell dephasing: welsh2004bioluminescence, nagoshi2004circadian.
- Maternally driven fetal rhythms: fetal_clocks_animal.

## 5. Figures (Cell Article, seven)

| fig | message | main panels (entries) | to the Supplement |
|---|---|---|---|
| 1 | A slow clock, read from untimed hospital tissue | (a) the clock plane and the two readouts (14 genes and the core loop without PER1, PER2 and NR1D1, with the gene-drop null; A336); companions as the separation of clock and signal layers (schematic); (b) living tissue, two biopsies hours apart: the clock rotates with elapsed hours, responsive sets do not, on 14 genes and on the core loop (A252, A256, A336); (c) deaths: the clock reads the hour of death, responsive and matched sets do not, also with a waveform fitted on other donors; baboon as the timed reference (A270, A335); (d) night shift: the core loop moves 4.06 h and keeps about half its amplitude (14 genes 2.68 h at about 0.31), PER1, PER2 and NR1D1 move less (not blind), the program 11.65 h at 0.854 and still 'behaviour' with the core constants (A257, A261, A335, A336); the trio alone is not readable (three-gene gates cannot pass) | method validation (A123b, A145, A209, A258, A285, A320); blood forced desynchrony (A262) |
| 2 | Structural cells keep the clock; every cell reads the signals | (a) participation across atlases and living donors (A272, A274, A280, A282, A283, A288); (b) lymphocytes share glucocorticoid and immediate-early states, not the clock (A275), with the bound and the core loop (A326, A321); (c) Bmal1 knockout: a shared signal state without a clock (A279) | per-atlas replicates (A281, A287, A284, A152); cortex unresolved (A280, A283); core loop (A336): the structural-lymphocyte contrast holds in the A274 / A280 atlases and post mortem, and only against matched genes in the A272 atlases and for macrophages (at the line), cortex falls (not decidable), living pooled D falls (descriptive); Bmal1 knockout at the line |
| 3 | Development: before birth the adrenal carries a clock-like organisation; liver and lung carry the clock genes' signal layer (revised 2026-11-23, user) | (a) before birth, the 14 clock genes covary between individuals in the clock's pattern in adrenal cortex, liver and lung at about two thirds of the adult size (A243, A244), replicated blind in the whole adrenal before 12 wpc (A328 part a H328.1, C 0.59); adrenal medulla, heart, kidney none at adult power; brain none (A237-A240, A244); (b) adrenal: the core loop covaries on its own, clock-specific, in the blind cohort (A328 part a H328.2, C 0.80); the states are two-dimensional, not one slow factor (A340 H340.1); fast pulsatile input excluded (post hoc: lead upper bounds at tau = 1 h, A328 part b, A338, A339); in the single-cell cortex cohorts the core loop is 0.42 (three quarters of the adult size excluded, half not; A266): the two adrenal datasets differ, cause open; (c) liver and lung: the covariation sits in the signal-responsive genes. Lung core loop small (0.47, below the half-size median 0.70, power 0.84; A269); liver core loop covaries (0.69) but not beyond stress-responsive genes (q95 0.71), and its states lie in the no-clock elongation range with the main axis along the glucocorticoid-target direction (A340, descriptive). Consistent with a systemic input; glucocorticoids are a candidate, not shown: adjusting for glucocorticoid-target activity keeps the covariation (A269 H269.2), the PER1 / PER2 / NR1D1 link to glucocorticoid targets is clock-specific in lung only (A307), glucocorticoid responsiveness does not rank the organs (A341); (d) mouse late embryo: hepatocytes, heart and kidney not (A309); (e) after birth the clock is present (A237, A239, D4); open: does the adrenal clock run on its own or follow a daily maternal / ACTH signal (the lead is not decidable: A328 part b, A338, A339) | negatives and bounds (A236-A241); phase-pattern pairs (A307); organ ranking (A341); core-loop power tables (A266, A269, A328); Cao design not estimable (A181); mother-fetus not resolved (A184); developmental loudness table (T5) |
| 4 | One time, many volumes | (a) classes share the organ's time, in human within 1 h with a small endothelial offset of the mouse sign (A273, A279; A293, A335; offsets descriptive, H293.3 FAIL; atlases 3 living, 1 mixed, 2 post-mortem); (b) composition does not move time, and its loudness effect is descriptive (A295, A335); (c) the organ order held without alignment to the day in slow deaths; in intensive care no large disruption at its detection limit (A293 Part 1, A349; post-mortem); (d) conserved class gains (A293 H293.4; 13 atlases: 6 living, 1 mixed, 3 post-mortem, 3 not classified); (e) organ gains over a departure that varies little: coupling reads mostly loudness; information in bits (A291, A294b, A311, A313; N2, A329; post-mortem); (f) the person scale (revised 2026-11-24): a body-wide expression scale the clock follows more reliably than random genes, in all fast deaths and in violent, instantaneous deaths (A336, A367), not larger than random genes' (size note); post-mortem only (A331 withdrawn, A366; no clock excess at matched phase in living women, A368); (g) the organ range is not composition (A323) and not post-mortem attenuation by ischaemia or RIN (A356); reproduction in living tissue not shown (A355, A330 not decidable) (revised 2026-11-23) | CHIRAL concordance (A123); order sensitivities (A109, A132, A314); A302; A310 fit; A312 per tissue |
| 5 | Perturbations move loudness, rarely time | (a) the map: every condition of section 3 as an arrow, each with its detection limit drawn as a circle around the origin (an arrow inside its circle reads 'not detectable', not 'no change'; A336 part 3), overall loudness drawn apart from the relative pattern; (b) intensive care lowers loudness and keeps the order (robust on the core loop and matched); brain death changes the order, also on the core loop and after ischaemic-time / RIN adjustment (A336; glucocorticoids on the core not excluded), and keeps loudness relative to matched genes; in GTEx, aliased with procurement and donor management; per-organ arrows descriptive (A293, A315, A317, A318, A335); (c) drives removed and in conflict, mouse: feeding removed (A297), central clock lost (A319), daytime feeding by organ weight (A322, A324); (d) the departure does not grow (A304, A306) | brain-death checks (A132b, A318 H318.2, A325, A327; Q41 Tasks 2-3, B1, B4); ventilation time (A132, A315; D3) |
| 6 | Disease: chronic disease spares the clock, local disease quiets it, inflamed skin may delay it | (a) chronic systemic disease: equivalence and powered nulls (A164, A233, A235, A306); (b) tumours keep the patient's time, their stroma carries it (immune cells at the line on the core loop, A227, A336), gain falls with proliferation, the departure stays flat (A138, A227, A271, A304), and the COFE correction (A265, A267); (c) chronic lesions keep time and lose loudness (A254); (d) psoriasis: the lesion later, at the line on the core loop (A308, A335; rule revised after A183), eczema alike (A183), uninvolved skin late (A115), its loudness (D1) | per-cancer ranking (A3); prognostic nulls (A191-A194, A198, A220); smoking lung gain (A306) |
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
| S1 fetal running loop | A328, A337, A338 | part a scored 2026-11-21: del Valle replicates blind, core loop PASS (Q43); timed control A337 scored (untimed statistic underpowered at n = 12); part b + A338 scored 2026-11-22: registered estimator depth-biased (power note), matched estimator: no consistent lead in del Valle or liver (H328.3 FAIL), uninformative (the running postnatal liver clock is not resolved either); A338 timescale not resolved; 10x cortex not decidable | Figure 3e stays "the open test"; 3(b) wording paused (Q43) |
| G1 organ factor | A330 | not decidable at this design (baboon power 0.60 at true rho 0.6; mouse 0.45 at 8 organs, not downloaded) | Figure 4g: "organ" stays open |
| G1b living pairs | A355 | scored 2026-11-23: living gain ratios do not follow GTEx as a set (H355.1 FAIL); omental / subcutaneous fat 1.38 vs GTEx 1.59, not decidable; cardia / squamous opposite (anatomy mismatch); swap-null power optimistic | Figure 4g: living reproduction not shown |
| G2 post-mortem attenuation | A356 | scored 2026-11-23: no effect of ischaemic time or RIN on gain within tissues (H356.1-2 FAIL, decidable); at most about 4 % of the organ spread | Figure 4g: "not post-mortem attenuation" |
| N2b departure as a systemic signal read unequally | A357 | scored 2026-11-23: person-level signal-gene state explains a small, clock-specific part of each tissue's departure (H357.2 PASS, median R^2 0.027 vs 0.008); tracks the GC-target state (r 0.77, descriptive); sharing across tissues not resolvable (H357.1 not powered) | model part 6: departure mostly a coherent clock-state shift of unknown cause; a small GC-like component read unequally |
| P1 person loudness, living | A331, A335 1c; A364-A370 | A331 withdrawn: its leverage is the shared biopsy time (A366); at matched phase no clock excess beyond random genes (A368); A364, A365, A369, A370 not decidable | Figure 4f: deaths only (all fast deaths and instantaneous deaths, leverage and size against companions); no living panel |
| N2 departure on the core loop | A329 (+ erratum 2026-11-23) | the core-loop fall (2.64 -> 1.28 times matched) and coupling fall are the 11-gene projector's noise gain (3.1-fold), not leakage: the departure is a coherent clock-state shift, 2.0 times matched gene-level variation, in 27 of 27 tissues; not the sample's GC / IEG state (H329.3 FAIL) | model part 6: departure carried mainly through signal-inducible clock genes, cause open |
| N1 output retention | A333 | not decidable at this design: median MDE log Phi 0.726 (only a loss to 0.48 detectable; gate 0.357); archive A91 (pooled, not blind) the only indication | Figures 1 and 7: N1 stays an open question |
| D1 skin loudness | A332 | descriptive: psoriasis variance up, not beyond companions | 3, Figure 6d: loudness change not clock-specific |
| D2 brain death per organ | A332, A335 1e | descriptive: no organ loses loudness; on the core loop tibial artery (+1.20 h) and nerve (+1.54 h) hold, lung does not, lower-leg skin appears (+1.33 h) | Figure 5b arrows: artery and nerve, descriptive (3 of 12 intervals exclude 0 on either gene set; 1.2 expected by chance) |
| D3 ventilation time | A332, A335 5 | organ-donor part only; the post-mortem remainder is not estimable (no tissue pair with 20 shared donors) | Figure 5 supplement: organ-donor part only |
| D4 developmental curve | A332 | descriptive: loudness beyond matched after birth in brain and kidney | Figure 3c; V read with C |
| cell mechanism of lost loudness (user, 2026-11-19) | A334, A334b | not decidable with public mouse 10x data: the lambda / R split fails the skin calibration (8 h mixing read 59 % as dephasing); aorta D not clock-specific (Bmal1 KO D 2.7 times aligned) | section 3: "a bulk loss of loudness can be desynchrony below the resolution" stays; decisive design needs a same-tissue clock knockout with D near 0 and plate-based depth |
| review checks (independent review 2026-11-19) | A335 | night-shift loss and shift hold on the core loop (s 4.06 h, R 0.239); ICU loudness loss holds on the core and age / sex matched (Q* +0.457, +0.459); psoriatic shift holds at the line (core q 0.049); brain-death loudness equivalent within +-0.357, relative to matched genes; classes share time within 1 h (not blind: interval known from A293); time-of-death reading holds with a held-out waveform | Q41: Article route by the pre-set rule; psoriasis stated with its core q; 'equal' for brain-death loudness now an equivalence; wording applied 2026-11-20 |
| two readouts, person loudness, detection limits, brain-death core (prompt `prompts/two_readouts_2026-11-20.md`) | A336 | tier 1 (55 rows): Figure 1 holds on the core loop; no fall attributable to PER1, PER2 and NR1D1 (A306 HTN gain mechanical); at the line on the core: tumour-normal pairing, GTEx strong-clock coupling, lesion loudness, psoriatic offset vs matched, elasticity gradient, MuTHER person leverage; GTEx person leverage holds on the core (0.714 > matched q95 0.621); brain-death D_rho holds on the core (0.465; adjusted 0.496); detection limits in `detectability.csv`; tier 2 (31 statistics): Figure 2 cell-class contrasts hold against matched genes, against responsive genes in three of six tests (A283 H283.2, A288 H288.2 at the line; cortex H283.5 not decidable), tumour immune pairing at the line (A227), class time, lymphocyte bound and A233 nulls hold | claim 1 keeps the person factor (deaths); claim 3 and Figure 5a: every 'kept' with its detection limit; Figure 5b: rearrangement in the core clock; Figure 1a: two readouts |
| review follow-up (2026-11-20) | A336 | registered 2026-11-20 (cc76c34), scoring pending: two readouts for every main-figure statistic, with a gene-drop null; GTEx person loudness on the core loop; detection limits; brain-death order on the core loop | `prompts/two_readouts_2026-11-20.md`; decides claim 1's person factor and Figures 1d, 4f, 5a, 5b |

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
| brain death | one panel, with its aliases (procurement, donor-management glucocorticoids) and the checks (Q41, A335) |
| GTEx deaths carry Figures 4 and 5 | living cohorts carry Figures 1-2, the lesions and P1; registry table for blind status |
| fetal results: no sampling times; adrenal running or driven undecided | state the bounds (A240); adrenal: "a clock-like organisation" (core loop, blind, A328 part a; geometry A340); liver and lung: "the signal layer in the clock's pattern"; the lead is not decidable (A328 part b, A338, A339) |
| the model is structural | A310 failed as a joint fit, and the paper says so; each part has its own test |
| human dynamics are not timed | time constants from mouse; one human point (A257: the core loop 4.06 h on day 3, A335) |
| PER1, PER2 and NR1D1 carry signals into the clock genes | every main-figure statistic on the core loop and as the three genes' contribution, with a gene-drop null (`prompts/two_readouts_2026-11-20.md`) |

## 8. Text fixes now

- **A257 excursion.** R = 0.233 is the misaligned change projected on the aligned chord. The excursion ratio is
  R / cos s, about 0.31 (derived). "Amplitude to a quarter" should read "excursion to about a third (a 2.7 h rotation
  alone would leave 0.77)".
  - To fix: the A257 headline in `results_preliminary.md`, SYNTHESIS lines 21, 41 and 136, and MODEL_federation
    line 63.
  - MODEL_signal_clock was corrected on 2026-11-19.
  - On the core loop (A335): projection 0.239 at a 4.06 h shift; excursion 0.239 / cos 60.9 deg = 0.49 (derived),
    matching the per-subject median 0.47: "keeps about half its amplitude".
- **Wording.**
  - "Loudness" (or "gain") for the clock's amplitude relative to the person.
  - "Time" for relative timing between compartments.
  - "Desynchrony" only for phase relations.
  - "Quieter", not "disrupted", for a loss of loudness.
  - "Tumours express less of the patient's time" stays (SYNTHESIS wording rules).
  - Before birth (revised 2026-11-23): adrenal, "a clock-like organisation (its core loop covaries on its own)", not
    "a running fetal clock"; liver and lung, "the signal layer in the clock's pattern", not "glucocorticoid-driven".
- **Review wording (2026-11-20), applied in sections 1-5.**
  - Psoriasis: "suggestive", with its core-loop q (0.049) and the rule history (A308 adopted after A183's rule failed).
  - Brain-death loudness: "equal relative to expression-matched genes".
  - Cell classes: "within 1 h, with a small endothelial offset"; that equivalence is not blind.
  - "Desynchrony is mostly loudness" only with "at tissue resolution", and not in the title or abstract.
  - No sentence claims that the brain maintains the order (D3 not estimable).
  - "Liver first" and "brain none" are dropped for human (A243, D4).

## 9. Decisions for the user and the PI

1. "One time, many volumes" as the message, and the clock plane as the grammar of the figures.
2. The seven-figure plan (section 5), with development as Figure 3 and time and loudness merged into Figure 4.
3. Which checks to register: S1 now (blind); G1, P1 and N2 next; N1 after its power scout; D1-D4 as descriptive.
4. Title.
5. Q25 recut: the two kinds of desynchrony become time (phase, rare) and loudness (common).
6. Figure 3 (development): main text only if S1 is scored blind with a timed positive control; otherwise the
   Supplement (independent review, change 3).
7. The person factor of claim 1: decided by the GTEx core-loop test and the gene-drop check
   (`prompts/two_readouts_2026-11-20.md`, parts 1-2).

## 10. Conversions

- Variance to loudness: square root of exp(log ratio) (A315, A254, A183). A loudness ratio, not a daily-amplitude ratio
  (METHOD section 9); do not set it on one scale with the timed amplitude ratios (A257, mouse).
- A304 gain: a slope (|C| / V_N), on the amplitude scale.
- A257: excursion ratio = R / cos s = 0.233 / cos 40.2 deg = 0.31; core loop (A335): 0.239 / cos 60.9 deg = 0.49.
- Night shift on the core loop: Delta = 0.87, loudness share 34 % (R 0.49, phi 4.06 h).
- Tumour purity sensitivity (A304): exp(-0.266) = 0.77.
- Organ gains: square root of the S range over the 27 peripheral tissues, S = 2 psi c / (1 - c) from the A104 table:
  9.24, and 3.10 without whole blood and testis; 2 psi spans 0.90-4.14 (4.6-fold).
- Delta^2 = (1 - R)^2 + 4 R sin^2(phi / 2), with phi in radians (hours x 2 pi / 24).
- Intensive care: exp(-0.404 / 2) = 0.82. Brain death: exp(0.046 / 2) = 1.02. Lesions: exp(-0.417 / 2) = 0.81.

## 11. Framing after the story audit (2026-11-23)

### 11.1 The mechanism, as tested in human

The paper has a systems mechanism, tested in independent examples:

| part | tested in | outcome |
|---|---|---|
| time is set by the drives an organ follows | mouse: A322 / A324 (daytime feeding), A319 (forebrain Bmal1 knockout), A129 (hepatocyte knockout); human: A293 (order without daily alignment), brain death (aliased) | holds qualitatively; the quantitative chord was rejected (A322 H322.1) |
| the clock's inertia carries time through acute states | A252 / A256 (through surgery), A270 (hour of death), A257 / A261 (core loop +4 h, program 11.65 h), A343 (glucocorticoid reaches PER1 transcription within hours) | holds |
| a fast signal layer read by every cell, a clock layer in structural cells | A275, A279, A281-A287, A321 | holds |
| loudness follows drive strength, conflict and composition | A315 / A335 (intensive care), A322 and A257 (conflicting drives), A304 (proliferation), A295 (composition, small) | holds as a pattern; the cell-level cause is undetermined (A334) |
| quantitative oscillator model | A310, A311, A312, A323 | failed or much smaller than predicted (Supplement) |

Credit the mouse consensus for the parts. The new content is the human demonstration of each part, and two human
rules: conserved class gains and the personal departure.

### 11.2 The conceptual claim: mouse-predicted desynchrony is absent in human disease

Each mouse expectation, its human test and the result:
- **High-fat feeding advances the liver clock by hours.** Excluded in T2D and obesity, equivalent within ±1 h (A235;
  MDE 0.37 h).
- **Metabolic disease "disrupts" peripheral clocks** (often from cultured cells). Coupling equivalent within ±0.15;
  no gain or departure change (A164, A306).
- **Ageing desynchronises organs.** Coupling equivalent; no organ shift larger than 0.20 h (A231).
- **Chronic disease and lesions shift local clocks.** No reordering (A233); lesions within 0.69 h (A254).

The mechanism explains the difference: in mouse, desynchrony needs conflicting drives (rest-phase feeding), which
human chronic disease rarely imposes. Where human drives conflict, time moves (night shift, A257), and locally in
inflamed skin. Wording: "the same mechanism predicts robust human time under chronic disease and ageing; the
desynchrony extrapolated from mouse models and cultured cells is absent at our detection limits". Do not write
"human clocks resist desynchrony".

Together with "the usual synchrony metrics read mostly loudness" (A291, A294b), this is the claim about the field:
true desynchrony in human tissue is rare and local, and most reported "disruption" is loss of loudness.

### 11.3 Figures

- **New supplementary figure: "the readout sees time move when it moves".**
  - Positive controls: A111 (Xin feeding, against timed truth), A129, A319, A273 / A279, A251 H251.4, A342 (offset
    against timed truth), A257.
  - Beside them, every null with its detection limit (A336 part 3; intensive care at the ventilated design, A349).
- **Figure 1** also shows A348: living tissues share one time on the core loop (4 of 5 tests; MuTHER fat-skin at the
  line).
- **Figure 1** adds A251 and A255 (living organs share one time; a person's cell line does not) beside A252 / A256, as
  positive controls of the method (revised 2026-11-24: expected from the definition of a circadian clock, not findings
  about organisation).
- **Figure 6d**: psoriasis is the exception.
  - Psoriasis delay: A183 blind, A308, A345 pre-mRNA, A347 pooled -1.51 / -0.96 h.
  - Psoriasis runs later than atopic lesions (A347: -1.25 / -0.75 h). The atopic delay is small, not clock-specific
    and not on the core loop; it is not shown to be zero (equivalence failed).
  - Hidradenitis is descriptive.
  - Non-skin lesions are kept (A254).
- **Intensive care** moves out of the "order kept" claims into the loudness claims.
  - The order statement is "no large disruption": one organ about 6-7 h, dispersion SD 2.25 h (A349).
  - Intensive care serves as the reference group for brain death.
  - The loudness loss is host-referenced too (A349 log z2 > 0) and largest in the organs with the strongest clocks
    (A349, Spearman about 0.82; descriptive).
- **Tumours**: "stay paired with the patient's time, carried at least by their fibroblasts" (A346). Do not write
  "only fibroblasts": endothelium, cancer cells and immune cells are not resolved on the core loop.
- **Supplementary registry table** with blind status (70 blind, 42 partly blind of 219; registry 2026-11-19) and a
  methods paragraph: no new human or animal sampling, every result from public or consented data, pre-registered.

### 11.4 Mechanism additions (to register; prompts on disk)

1. A waveform screen of uniform scaling against gene-specific change, across tumours, lesions, intensive care and
   psoriasis. Calibrate it on the A212 plasmode first. It separates "intact but quieter" from "rewired", not
   dephasing from per-cell damping.
2. A cytokine dose-response of the skin delay: IL-17 / TNF against Th2 strata; pre-mRNA in Tsoi.
3. In-house single-cell time-lapse of clock and cell cycle: does faster division dephase or damp the clock? This
   needs the PI (Q44).

Prompts: `prompts/audit_followups_2026-11-23.md` (gaps) and `prompts/mechanisms_2026-11-23.md` (items 1-2).

**Note (2026-11-23, after the mechanism runs; proposed, not applied).**
- Item 1 (A351): calibration failed (V351: per-tumour rewiring against damping AUC 0.80 / 0.89, bar 0.90): "intact but
  quieter" cannot be told from rewiring that differs between tumours at gene level in bulk; no condition was scored.
  Shared rewiring is separable (AUC >= 0.975): a narrower screen ("shared rewiring against uniform scaling") would be
  calibrated and needs its own registration (user decision).
- Item 2 (A352): not powered (primary MDE 1.19 h against 0.75 h); closed at the power stage. The mechanism of the
  psoriatic delay stays open; A345 (pre-mRNA delay, clock-specific narrowly) remains the only mechanistic step.
- Wording: keep "loudness, not time" as a pattern; do not state "intact but quieter" as shown.

## 12. Evidence base by level (2026-11-24, after A361-A373 and the user's remark that the living anchors are positive controls)

| level | finding | living evidence | post-mortem evidence |
|---|---|---|---|
| method | the clock turns with elapsed hours; organs of a person share one state; a cell line does not | A251, A252, A255, A256 (positive controls) | A270, A335 (hour of death) |
| cell | structural cells carry the clock, lymphocytes do not (ranking) | A282: five living atlases (HLCA living donors, breast, colon, Kong ileum, Kong colon), beyond both companion families | A272, A274, A282 post-mortem arm |
| cell | classes share the organ's time within about 1 h; endothelium slightly early | A293 part 2: kidney, breast, Kong ileum (living), HLCA (mixed) | optic nerve, anterior eye; the order test failed (H293.3) |
| cell | class gains conserved across organs | A293 part 3: breast, muscle, colon, kidney, Kong ileum, Kong colon (living), HLCA (mixed) | heart, optic nerve, anterior eye; endometrium, islets, Natri lung, skin AD control not classified in A282 |
| organ | gains about ninefold, order, coupling reads loudness | not shown (A371 not decidable: the organs in living cohorts share one phase; A355, A330 not decidable) | A104, A291, A293 part 1, A323, A356 |
| person | body-wide scale, clock follows it most reliably | not shown (A366, A368 fail; A364, A365, A369, A370 not decidable) | A336, A367 |
| perturbation | behavioural shift costs loudness and moves the core loop | A257 (timed, living) | - |
| age | clock and outputs keep amplitude to 70 | - | A373 (matched; equivalence) |

Living organ-level evidence needs STARNET (dbGaP application, user and PI). Atlas donor type from A282 `arms.csv`;
the four unclassified atlases are to be checked before their donor type is written in the paper.
