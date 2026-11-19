# Signal and clock: a two-layer reading of human circadian time (proposal, 2026-11-19)

Proposal for the user and the PI. Prompt (user, 2026-11-19): after PI_QUESTIONS, SYNTHESIS_2026-11-19,
REVIEW_2026-11-19 and MODEL_federation_2026-11-19, propose a mathematical framework or another way to assemble the
findings into an updated view of the human circadian system. The user also observed that the clock acts like a watch
with inertia, while the tissue's own output genes shift without inertia, as if the tissue did not read its clock.

The frame keeps the six parts of the federation and adds three things: a fast signal layer, the clock's time constants
and one readout term. Nothing is registered and no analysis was run. Numbers were re-read in the output tables listed
in section 10, and IDs refer to `results_preliminary.md`.

## 1. The view in one paragraph

Human tissues keep two kinds of time.
- **Signal layer (fast).** Glucocorticoid, feeding, activity, temperature and immediate-early responses reach every
  cell, lymphocytes and fetal cells included, and follow the drives within hours.
- **Clock layer (slow).** Stromal-vascular, epithelial and parenchymal cells and macrophages integrate the same drives
  over days.

The clock's phase is set by the mix of drives an organ receives, and the living brain sets the body's organ order
through the drives it makes. The clock's loudness is set by the cells. Because it is slow, the clock carries the body's
time through acute states that reset the signal layer: surgery, interventions, the hours of dying. That is why single
tissue samples taken in hospital hours can be read for internal time. When drives conflict or fade, the clock pays in
amplitude before it moves in time. Most of a tissue's daily program follows the signal layer: in human muscle after a
12 h behavioural shift, the program moved with behaviour while the clock moved 2.7 h. So the peripheral clock acts less
like a watch the tissue reads than like a flywheel. While the signals are present, it does not time the daily program;
when they fail, it holds the body's time and the order between organs.

## 2. The observation: a clock with inertia that its own tissue does not follow

| design | clock layer | signal layer or daily program | entries |
|---|---|---|---|
| human muscle, 3 days of a 12 h behavioural shift | moves 2.68 h (90 % 1.81-3.59), excursion about 0.31 (0.233 projected on the aligned chord) | genome-wide program moves 11.65 h (10.54-12.76), excursion 0.854 ("shifted or evoked") | A257, A261 |
| living people, two biopsies hours apart, 16 studies | rotates with elapsed hours (rho 2.533); interventions do not add (H256.4) | responsive sets do not rotate (median -0.380, q95 0.779) | A252, A256 |
| GTEx fast deaths (hours of dying) | reads the hour of death (R 0.631) | responsive sets do not (median R 0.073; matched 0.054) | A270 |
| baboon, living, timed (reference) | reads ZT (R 0.987) | responsive and matched sets read ZT too (0.537, 0.622) | A270 |
| mouse liver, feeding rhythm removed | keeps 0.724 of amplitude, +1.57 h | > 70 % of cycling outputs lose rhythm (literature) | A297; greenwell2019rhythmic |
| tissue lymphocytes | no readable clock (bound: below about 0.7 of fibroblasts) | share glucocorticoid-target (0.390) and immediate-early (0.392) states like other cells; clock genes at the matched level (0.172 vs 0.167) | A275, A281, A321, A326 |
| fetal adrenal cortex, lung, liver | no core loop shown | covariation carried by the signal-induced PER1, PER2, NR1D1 | A243-A269, A307 |
| mouse aorta, Bmal1 KO | no timed clock (coupling 0 to 0.2) | stromal-vascular classes still share clock-gene variation (abs r 0.63-0.75) | A279 |

**Already known.** The mouse literature describes three layers: a small clock-autonomous core, a large layer driven
by feeding, activity, temperature and hormones, and a gating layer in which the local clock integrates systemic signals
(kornmann2007system, koronowski2019defining, greenwell2019rhythmic, weger2021systematic, greco2021integration;
`literature_positioning.md` section 11). Wefers 2018 reported the muscle dissociation on the arrays used in A257 and
A261. "Outputs follow signals" is not new as a statement.

**New here.**
1. In human tissue, at body scale, the two layers separate by cell class (lymphocytes carry the signal layer only) and
   by development (fetal organs carry the signal layer first).
2. The clock's inertia is shown with one readout in three human settings (living interventions, deaths, a behavioural
   shift), each time against the signal layer.
3. Inertia is the property that makes untimed hospital tissue readable. The specificity standard (clock beyond
   responsive genes) is the operational separation of the two layers, so the method's companions become biology.
4. Across perturbations, the clock pays in amplitude before it moves in time (section 4).

**The surprise, restated.** People read a watch; nobody reads a flywheel. A flywheel turns with its engine, resists
sudden changes, keeps turning when the engine stutters, and slows when two engines push against each other. The data
show each of these:
- the clock lags a behavioural shift (A257);
- it keeps its amplitude when one drive is removed (A297, A315);
- it loses amplitude under conflicting drives (A257, A322, A324).

The flywheel picture also predicts something not yet shown in human: when the signals are absent, the outputs that
persist should follow the clock. Check N1 tests this.

## 3. Formal version

First (24 h) harmonic, complex amplitudes, envelopes varying over days. Notation as in MODEL_federation where shared.

    drives          U_d = A_d exp(i Phi_d)      d: central (while the brain lives), meals, activity/sleep,
                                                other (ICU care, light by non-SCN routes)
    signal layer    S_ok = sum_d beta_kd U_d + eta_o
                    every class k, lymphocytes included; follows U within hours; eta_o: acute state at sampling
    clock layer     Z_o = r_o exp(i phi_o), carried by clock-carrying classes only (class gain a_k > 0)
      phase         kappa_o dphi_o/dt = |V_o| sin(arg V_o - phi_o),       V_o = sum_d w_od U_d
      amplitude     tau_r dr_o/dt = g_o h_o(rho_o) - r_o,                rho_o = |V_o| / sum_d w_od A_d  (drive coherence)
      gain          g_o = A_p G_o sum_k s_ok a_k                         person x organ factor x composition x class gain
    outputs         Y_og = alpha_g Z_o + beta_g S_o      (clock gating folds into alpha and beta at the first harmonic)
    readouts        clock genes x_o = Z_o + lambda S_o + eps;  responsive genes S_o + eps';  matched genes eps''

- **Two time scales.** The signal layer follows within hours. The clock's phase relaxes with tau_phi = kappa_o / |V_o|
  (days; slower when the drives are weak or in conflict), and its amplitude relaxes with tau_r.
- **Shape of h_o.** h_o(1) = 1.
  - h_o(rho) >= rho where cells hold their amplitude and resynchronise: lung, heart and kidney sit above the chord.
  - h_o(rho) <= rho where cells split between drives: adipose sits at or below it.
  - With no drive at all, r decays at the dephasing rate (mouse: days to weeks) and phi drifts.
- **Leakage term.** lambda is not 0 because PER1, PER2 and NR1D1 are induced by signals (A264, A266). PER1 is the most
  shared clock gene and a direct glucocorticoid and immediate-early target (A275).

**What clocksync measures, in the model** (normal life, steady state):

| clocksync | model |
|---|---|
| host | phase of the person's shared Z (internal time) |
| offset arg z | arg V_o relative to the host: the organ order is the order of the weighted drive mixes |
| gain | g_o h_o(rho_o), that is g_o in normal life |
| departure 2 psi | Var(Z_o around its mean) + lambda^2 Var(eta_o) + technical variation |
| coupling c | S / (S + 2 psi) |
| "clock-specific" | Z beyond S (responsive family) and beyond eps (matched family) |
| readable from one untimed sample | Z at sampling is Z over the last tau: an acute state moves S, not Z |

**Relation to the federation.** Parts 1-6 stay unchanged. Added:
- the signal layer S, a measured state shared by all cells;
- the time constants (tau_phi and tau_r, days, against hours for S);
- outputs as the mixture alpha Z + beta S;
- the leakage term lambda S in the clock readout.

Not added: a population fit. A310 failed as a joint fit, so h_o, kappa_o, alpha and beta are not fitted jointly; each
part is tested on its own, as in the federation.

**Box (optional, for the method chapter): coupling as information.** For a complex Gaussian channel, the information
a tissue's clock state carries about the person's internal time is I = log2(1 + S / 2 psi) = -log2(1 - c) bits per
sample. Over the GTEx periphery this is 0.10 to 2.77 bits (c 0.067 to 0.853; A291). Departures are only weakly
correlated between organs (median r 0.13 to 0.16; A313), so signal-to-noise ratios roughly add across organs. That gives
the precision of a person's internal time read from k organs. It also turns "coupling mostly measures loudness"
(REVIEW 2.3) into a positive statement: loud organs carry more of the body's time.

## 4. The results as perturbations of the model

| perturbation | model | clock layer observed | signal layer / outputs observed | entries |
|---|---|---|---|---|
| acute state at sampling (hours) | eta moves S; Z unchanged | internal time read in surgery and deaths; rotation with elapsed hours | does not read time in deaths and does not rotate in living tissue | A252, A256, A270 |
| meals and activity removed, brain alive (ICU); mouse arrhythmic feeding | arg V and rho kept: order kept; amplitude lower, slowly | order kept (rho 0.690 against fast deaths); clock-plane variance 0.84 of fast deaths (0.67 against matched genes); lower with longer ventilation (organ donors, descriptive); mouse liver 0.724, +1.57 h | human not measured (N1); mouse: most liver outputs lose rhythm | A293, A315, A297 |
| a drive displaced (night shift; daytime feeding) | arg V moves by the weight on the displaced drive; rho < 1: amplitude falls, phase moves slowly | muscle 2.68 h at about 0.31; mouse liver 11.3-11.7 h at 0.44-0.71; lung, heart, kidney 1.6-2.5 h at 0.50-0.65; adipose 0.11-0.32 | muscle program 11.65 h at 0.854; blood program 3.96 h at 0.355 | A257, A261, A262, A322, A324 |
| central drive replaced (brain death; SCN-clock loss) | new arg V per organ from the remaining drives: new order; amplitude held while a drive remains | GTEx: different organ order (D_rho 0.621), shared clock state kept, amplitude equal to ventilated living-brain donors (clock-specific ratio 1.05); aliased with procurement. Mouse: new reproducible order under feeding, drift in DD | not measured | A317, A318, A325, A327, A315, A319 |
| receiver changed (proliferation, cell class, composition) | g changes, phi does not | proliferating tumours: gain 0.63, no drift; class gains conserved; composition moves gain, not phase | - | A304, A293 Part 3, A295 |
| cells without a clock layer (lymphocytes, fetal, Bmal1 KO) | a_k near 0, beta > 0 | no readable clock (bound 0.7); no core loop before birth; no timed clock in the KO | shared signal states | A275, A281, A321, A326, A264-A269, A307, A279 |
| nothing changed (chronic disease, age, sex) | no parameter moves | couplings, offsets, gain and departure equal (powered) | - | A164, A208, A231, A233, A235, A306 |
| local change (psoriasis lesion) | phi of the lesion moves | -1.45 h, clock-specific against both families | - | A183, A308 |

### The phase-amplitude map (one panel for the whole paper)

x axis: shift in hours, or the organ order where a shift is not defined. y axis: amplitude relative to the reference
(amplitude = square root of a variance ratio where the entry reports variance).

| layer | condition (reference) | shift | amplitude ratio | entry |
|---|---|---|---|---|
| clock | human muscle, day 3 of a 12 h behavioural shift (aligned) | +2.68 h | about 0.31 (projection 0.233) | A257 |
| program | same subjects, genome-wide rhythmic program | 11.65 h | 0.854 | A261 |
| program | human blood, forced desynchrony (day 1) | 3.96 h | 0.355 | A262 |
| clock | mouse liver, feeding rhythm removed (ad libitum) | +1.57 h | 0.724 | A297 |
| clock | mouse liver, daytime feeding (Manella; Xin) | 11.69; 11.29 h | 0.712; 0.438 | A322, A324 |
| clock | mouse lung (Manella), heart and kidney (Xin), daytime feeding | 1.90; 1.59; 2.49 h | 0.495; 0.589; 0.648 | A322, A324 |
| clock | mouse WAT (Manella), VAT (Xin), daytime feeding | 10.60; -0.48 h | 0.319; 0.111 | A322, A324 |
| clock | human ventilated, living brain (fast deaths) | order kept (rho 0.690) | 0.92 (0.82 against matched genes) | A293, A315 |
| clock | human brain-dead organ donors (ventilated, living brain) | order changed (D_rho 0.621) | 0.93 (1.02 against matched genes) | A317, A315 |
| clock | proliferating tumours (less proliferating) | no drift | 0.63 | A304 |
| clock | psoriasis lesion (own skin) | -1.45 h | - | A308 |
| clock | women vs men, old vs young, T2D, obesity | 0.04, -0.05, 0.02, 0.20 h | equal | roadmap A1, A164, A208 |

How the model reads the map:
- **Signal layer and daily programs** sit at the shift of their drive, with amplitude near 1 (muscle) or lose amplitude
  when they mix drives (blood).
- **Clocks** sit lower, at the shift set by the organ's weight on the displaced drive.
- **Time moves without much amplitude loss** only when the configuration of the drives changes: brain death, or
  inverted feeding for the high-weight liver.
- **Receiver changes** move amplitude only.
- **Chronic disease, age and sex** sit at the origin.

## 5. The updated view

1. **Two layers.** Every cell reads the body's signals within hours. Only structural cells (stromal-vascular,
   epithelial, parenchymal) and macrophages keep a clock. Tissue lymphocytes and fetal organs carry the signal layer
   without a readable adult-size clock.
2. **The clock is a flywheel.** It carries the body's time through acute states and lags behavioural shifts. The daily
   program follows the signals (muscle, mouse liver).
3. **Time from the drives, loudness from the cells.** An organ's phase is the weighted mix of its drives. Its loudness
   is the product of class gain, composition, an organ factor (open: G1) and the person's shared amplitude.
4. **The brain sets the order through the drives it makes.** With a living brain, the order holds without meals and
   does not require alignment to the external day. In GTEx brain-dead donors it is different (aliased with
   procurement); in mouse, a substitute drive imposes a new order.
5. **Perturbations are paid in amplitude first.** Conflict, intensive care and proliferation lower the clock's
   amplitude. Time moves when the configuration of drives changes, or locally (psoriasis). Chronic disease, age and sex
   change neither.

| textbook (mouse-derived) | human tissue, untimed (this paper) |
|---|---|
| the SCN sets peripheral clocks, which time their tissue's physiology | the brain and behaviour make systemic signals that every cell reads within hours; structural cells integrate them into a slow clock that holds time and order |
| every cell type keeps a clock | stromal-vascular, epithelial, parenchymal cells and macrophages do; tissue lymphocytes read signals with no clock readable at fibroblast depth |
| desynchrony means clocks shifting apart | perturbations cost amplitude first; time moves when the drives change |
| organ phases follow feeding | the order holds without meals while the brain lives |
| disease disrupts clocks | chronic disease leaves the organisation intact; tumours express less of the patient's time as they proliferate |

**The question the paper raises** (Cell's criterion): why do human tissues keep a clock that their own daily program
does not follow? The frame's answer is a hypothesis: the clock is the body's memory of time. It is read when the
signals are absent, noisy or in conflict (night, fasting, illness, intensive care), and it keeps the order between
organs.

## 6. Checks (not registered; power first; ranked by what each decides per hour of work)

**N1. Intensive care as the human version of arrhythmic feeding (direct test of the flywheel).**
- *Setting.* Ventilated donors with a living brain have no daily rhythm of meals or activity, while the central drive
  and the clock remain (A293 H293.2, A315).
- *Statistic per tissue:* Phi_o = (L_out ventilated / L_out fast) / (L_clock ventilated / L_clock fast). L is the
  loading on the leave-one-out host of:
  - (a) an output program: non-clock, non-responsive genes that follow the host in fast deaths (A91's
    "clock-following" class), selected on a held-out half of the fast-death donors;
  - (b) the tissue's own clock state.

  Host reliability cancels in the ratio.
- *Readings.*
  - Phi near 1: the program follows the clock, or the central drive, without behaviour.
  - Phi well below 1: it follows behaviour, as most mouse liver outputs do (Greenwell 2019).
  - A gradient with meal-driven tissues lowest. The archive's A33 found the clock itself losing most coupling in liver,
    pancreas and muscle.
- *Companions.* Matched non-rhythmic sets (floor) and the responsive family.
- *Extension.* In brain-dead donors, the host is rearranged, so use the loading of the program on the tissue's own
  clock state. It is kept if the tissue reads its clock and lost if the program reads central signals.
- *Limits.*
  - The central drive persists in (a), so Phi counts centrally driven outputs with the clock.
  - Continuous feeding is not recorded.
  - 79 ventilated post-mortem donors: power first.
  - Donor steroids confound the extension.
- *Cost.* A91 machinery and the A315 groups: hours.

**N2. Is the personal departure the signal layer read through PER1, PER2 and NR1D1?**
- *Setting.* A294b rows (GTEx fast deaths, 27 peripheral tissues), refitted on the core loop: the 11-gene waveform
  without PER1, PER2 and NR1D1, as in A321. Each gene set is compared against 100 matched sets of its own size.
- *Statistics.*
  - The H294b.1 excess (log of clock own variance over the matched median) on the core loop, against the full set
    (full set: 0.969, about 2.6-fold).
  - Coupling on the core loop against the full set, paired over tissues.
  - Within tissue, the projection of each sample's departure on its responsive state (A275 glucocorticoid and
    immediate-early classes), against matched sets.
- *Readings.*
  - If the departure is leakage, the core-loop excess falls toward 0, coupling rises, and the projection is positive.
    Model part 6 then reads "the acute signal state seen through signal-inducible clock genes", and clocksync users get
    a core-loop readout for coupling.
  - If the excess stays, part 6 is a person-specific clock departure, with its cause open.
- *Blind status.* The core-loop departure of these tissues has not been computed; check before registering (A296 used
  the core loop for the adrenal offset only).
- *Cost.* One A294 run, about 25 min on 4 workers.

**N3. Night shift: does each signal family follow its own drive?**
- *Setting.* GSE106800, as in A257 and A261, with glucocorticoid-target and immediate-early families (A275 lists)
  placed against the clock and the program.
- *Prediction.* Immediate-early and meal-responsive families move with behaviour at full excursion. Glucocorticoid
  targets follow cortisol (central drive, partly sleep-linked) and sit in between.
- *Power first.* A261 needed about 600 reference-rhythmic genes per program to clear its gate, so families of 11 to 17
  genes may not be estimable.

**N4. Within intensive care, phase held and amplitude spent (descriptive).**
- *Setting.* Ventilated donors with a living brain, by ventilation tertile.
- *Readouts.* Clock-plane variance (A315 statistic) and the organ order against fast deaths (A293 statistic).
- *Prediction.* Amplitude falls with duration while the order holds. Within organ donors amplitude already falls with
  ventilation (A315 descriptive, clock-specific log ratio -0.236). This would give the order of magnitude of tau_r in
  human.
- *Power.* 79 donors: descriptive only.

**Recommendation.** Run N1 and N2. Both use machinery and data on disk. N2 decides model part 6; N1 is the one human
test of the frame's new claim. Run N3 and N4 only if their power scouts are cheap. The rest of the claim set stays
frozen (REVIEW 2.4).

## 7. How the paper reads in this frame

**Title candidates.**
- "Human tissues keep a slow clock behind a fast signal layer"
- "The human peripheral clock is a flywheel"
- "Two layers of time in the human body, read from hospital tissue"

**Opening.** "CHIRAL established that human tissue clocks are synchronised with offsets of a few hours. We ask which
cells keep that time, what sets it, and what tissues do with it."

| fig | message | main panels (entries) |
|---|---|---|
| 1 | Two layers in every organ: all cells read the signals, structural cells keep the clock | readout schematic with both layers; stromal-vascular, epithelial, parenchymal cells and macrophages vs lymphocytes (A272, A274, A283, A281, A321, bound A326); gene classes shared by cell class (A275); Bmal1 KO: shared signal state without a clock (A279); classes share the organ phase, composition moves gain (A273, A279, A295) |
| 2 | The clock is slow: it keeps the body's time through acute states and lags shifts; the daily program follows the signals | elapsed hours in living tissue beyond responsive genes (A252, A256); hour of death read by the clock, not by the signal layer, baboon as the timed reference (A270); night shift, clock 2.68 h at about 0.31 against program 11.65 h at 0.854 (A257, A261), blood (A262); mouse removal and conflict (A297, A322, A324); N1 if run |
| 3 | Loudness comes from the cells | class gains (A293 Part 3); composition and class gains do not explain the organ range (A323); gain over a departure of near-constant size (A291, A294b, A311; N2 if run); person amplitude (A312) |
| 4 | The order between organs is set by the brain's drives | order without meals and without alignment to the external day (A293 Part 1); GTEx brain-dead donors: different order, shared state and amplitude kept, aliased with procurement (A317, A318, A325, A327, A315); mouse (A319) |
| 5 | Perturbations are paid in amplitude | the phase-amplitude map (section 4): chronic disease, age, sex; psoriasis; tumours (A138, A304); ICU; brain death; night shift; mouse conflicts |
| 6 | Model | signal and clock, tested parts and failures (A310, A311, A312, A313) |
| ED | Before birth: signal before clock | A243, A244, A264, A266, A269, A307, A309 |

**Q41.** The frame gives the Article route a high point that does not rest on one aliased cohort. Figure 2 rests on:
- living interventions (16 studies);
- GTEx deaths, with a timed baboon reference;
- the night-shift design;
- two mouse datasets.

Brain death becomes one consequence (Figure 4), stated with its alias. The frame also suits the Resource route:
inertia is why hospital tissue is readable, and the specificity standard is the separation of the two layers.

## 8. Limits of the frame

1. **"Outputs follow signals" in human rests on one design**: A261, on the Wefers 2018 arrays, with a masking limit
   ("shifted or evoked"). The rest is mouse work; blood follows neither (A262). Until N1, the claim is stated for muscle
   and mouse liver only.
2. **The responsive family is not a clean signal readout.** It mixes systemic and handling states (A275). Its q95 rule
   is not calibrated in bulk (A285), while the broad family of A320 is.
3. **Inertia is not timed in human.** tau_phi and tau_r are mouse numbers (saini2013realtime, husse2014lightdark,
   izumo2014differential) plus one human point: 2.68 h after three days.
4. **The equations are a frame, not a fit.** A310 failed as a joint population fit, and the same discipline applies.
5. **Two layers is the mouse consensus**, which has three layers counting gating. What is new is the human, cell-class
   and body scale, the link to untimed readability, and the amplitude-first rule across perturbations.
6. **Person amplitude.** Tissues follow the person's shared amplitude nearly in proportion (A312: elasticities above
   1 in 11 tissues), which a saturating h does not produce. Either that variation is a donor-level technical scale in
   deaths, or human drives are not saturating. The frame does not decide.
7. **Brain death.** Donor management includes high-dose methylprednisolone (A132 registration), and glucocorticoids
   reset peripheral clocks (balsalobre2000resetting). In the frame this is one more substitute drive. It does not
   separate loss of the central drive from donor management or procurement.

## 9. Decisions for the user and the PI

1. Adopt signal and clock as the paper's organising model (it extends the federation; MODEL_federation stays as parts
   1-6), or keep the federation as is.
2. Register N1 and N2 (power first), or keep the freeze and state the flywheel prediction as untested.
3. The figure plan of section 7, with Figure 2 (inertia and the two layers) as the high point in place of brain death.
4. Title and the question the paper raises.

If adopted, SYNTHESIS section 1 and MODEL_federation section 1 get a pointer here; nothing else changes until N1 and
N2 are scored.

## 10. Sources (re-read 2026-11-19)

`exploratory/` paths:
- `A257_nightshift_muscle/output/score_tests.csv`
- `A261_nightshift_programs/output/tests.csv`
- `A262_blood_misalignment/output/tests.csv`
- `A256_elapsed_replication/output/score_tests.csv`
- `A270_responsive_time/output/{tests_human,tests_baboon}.csv`
- `A275_shared_gene_classes/output/classes.csv`
- `A279_aorta_pairwise_truth/output/pairs.csv`
- `A293_organisation_levels/output/part1/groups.csv`
- `A294_floor/output/tests.csv`
- `A297_feeding_model/output/arms.csv`
- `A304_cancer_gain_departure/output/score/tests.csv`
- `A311_hopf/output/tests.csv`
- `A312_elasticity/output/tests.csv`
- `A315_braindeath_gain/output/{tests,descriptive}.csv`
- `A317_braindeath_confirm/output/tests.csv`
- `A322_conflicting_drives/output/{ratios,xin/ratios}.csv`
- `A323_gain_decomposition/output/tests.csv`

From results entries: A308 (-1.45 h), roadmap A1 offsets (`roadmap_A1_offset_bound/output/offsets_hours.csv`), A326
bound, and A291 coupling range. Archive A33: `scCircadianCancer/docs/RESULTS_preliminary.md`.

Conversions:
- variance to amplitude: square root of exp(log ratio);
- A257 excursion: R (0.233) is the misaligned change projected on the aligned chord; the excursion ratio is
  R / cos s = 0.233 / cos 40.2 deg = 0.31 (correction 2026-11-19; earlier versions of this file plotted 0.233);
- bits: -log2(1 - c).
