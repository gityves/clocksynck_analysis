# Results section: slide by slide

Order (agreed 2026-11-18): R1 method (description and validation, short); R2 organisation in health (organs,
composition, cells, age and sex); what holds the organisation (loss of daily cues, central drive); development (only if
re-scored); internal desynchrony in disease; synthesis. Limits inside each chapter plus one final slide. Every number below was re-read from the table named in its slide's **Source** line. GTEx figures show
per-tissue aggregates only, never per-donor times.

Status labels: **registered** (primary or secondary, scored as written), **confirmation** (archive result known,
re-run on current code), **seen data** (registered after the data were looked at), **post hoc**.

---

# R1. Model validation

**Chapter message:** where the true time is known, the untimed method recovers it. Organ phase map from
working-hours deaths, patient order in living surgery, organ offsets in a full-day primate atlas, and the two
perturbations that should move it (a liver without its clock, a feeding shift) are detected.

---

## R1.1 GTEx: the organ map from deaths between 08 and 18 h

**Message:** from donors who all died in working hours, the method recovers the organ phase map and the organ
amplitudes that whole-day timing gives.

**Design:** GTEx v10 fast deaths (Hardy 1-2). Test set: the 111 donors who died 08-18 h, 40 tissues with ≥ 20 such
donors. Truth (i): cosinor on the recorded time of death of all fast deaths. Truth (iii): the same on the 133 donors who
died outside 08-18 h only (disjoint from the test set).

**Bullets**
1. Offsets against truth (i): circular r 0.945 [0.907, 0.972], median error 0.35 h. Amplitudes: Pearson (log2) 0.976.
2. Against the disjoint truth (iii): amplitudes 0.965 (pass). Offsets 0.679 [0.467, 0.839], below the registered bar
   of 0.7 (**fail as registered**).
3. Why (post hoc): the estimate reproduces across disjoint donor halves (0.873), the night-only truth does not (its
   halves agree at 0.544). The shortfall comes from the truth's noise, not from the estimate.

**Figure:** `figures/gtex_organ_map/gtex_organ_map.png` (`scripts/gtex_organ_map.py`), the CHIRAL-paper organ map
(angle: tissue offset, radius: relative amplitude, upper half of the clock face, -6 to +6 h) for A. the timed truth (i),
B. clocksync and C. CHIRAL, both untimed on the 08-18 h donors; in B and C each tissue is joined to its truth (grey).
To add if wanted: the split-half agreement (estimate halves vs truth (iii) halves).

**Speaker notes**
- Why the best direction and rotation: an orderer's cycle has an arbitrary direction and origin; scored as run, a
  reversed draw negates every tissue offset (CYCLOPS was reversed in all full-day draws), which made the orderers look
  worse with the full day than in hospital hours. Aligned for every method alike, the full day helps everyone, and the
  advantage of clocksync is specific to hospital hours (post hoc scoring, 2026-11-17).
- The orderers' working-hours failure is not a reflection or stretch artefact: allowing the better direction and a
  least-squares stretch, CHIRAL stays at |r| 0.165 (post hoc check).
- CoFE genome-wide and CYCLOPS also fail with the whole day (0.139, 0.115): they order donors along non-clock
  variation.
- The distinctive property on real data is partial-day coverage, not coupling estimation as such: on idealised
  whole-day simulations two-step orderer pipelines recover coupling nearly as well (A145).
- DIP vs internal time (post hoc diagnostic, `scripts/gtex_dip_diag.py`, full day): the two disagree at DIP 12-18 h,
  where CHIRAL places only 35 donors (62 died in that interval). These bins hold the weakest shared clock states (0.75-0.80
  of the median donor) and their recorded times are scattered (R 0.32-0.44); at DIP 15-18 h the DIP is > 6 h off the
  recorded time for 32 % of donors (clocksync 21 %). Weak clocks are harder for both (Spearman of strength with error
  −0.27 clocksync, −0.40 DIP); strongest third 1.41 vs 1.44 h, weakest third 2.87 vs 4.18 h.
- Why not circular r: over the full day the DIP's |circular r| with time of death is higher (0.677 vs 0.600), but the
  difference is not resolved (paired bootstrap 95 % CI [−0.46, +0.40]): with deaths spread over the day the
  Jammalamadaka-SenGupta mean direction is ill-defined. The median error difference is (DIP − clocksync +0.53 h
  [0.16, 0.86]).
- Status: H123.4 registered, pass.

**Source:** `exploratory/A123_working_hours/output/a123b/{tests,h123_4,posthoc_reflect}.csv`.

---

## R1.3 Living patients: internal time follows the surgery hour

**Message:** in living humans, each patient's untimed internal time follows the recorded hour, and this is carried by
the clock genes.

**Bullets**
1. Breast surgery (42 patients, tumour + adjacent normal, resection 08-18 h): Spearman of internal time with the
   recorded hour 0.558; 200 expression-matched non-clock gene sets: median −0.053, none reaches it (p 0.005).
   Median error with no recorded time used 1.42 h (shuffled hours 2.28 h).
2. Heart surgery (16 patients, 8 operated 06-10 h, 8 at 14-18 h; four cell types): the afternoon host sits away from
   the morning host (p 0.0004), in the direction forward time predicts (19.7 h against 18 h, within the registered
   ±3 h).
3. Subcutaneous adipose (12 men, 48 biopsies over the day): circular r with clock time 0.822 (p 0.001), median error
   1.11 h.

**Figure:** `figures/living_patients/living_patients.png` (`scripts/living_breast_compute.py`, `scripts/living_patients.py`).
A. breast: internal time (13:00 rule, no recorded time used) vs resection hour, one dot per patient, same range on both
axes. B. heart: each patient's internal time by surgery window (06-10 h vs 14-18 h), circular means. C. adipose:
single-biopsy internal time vs zeitgeber time (common offset removed), circular mean and SD per time point. D. per-sample
absolute error after removing one common offset: breast (08-18 h) 1.12 h, adipose (full day) 1.11 h.

**Speaker notes**
- Breast: the ordering does not rest on glucocorticoid-induced PER1 (0.560 without it); it leans on NR1D1 (0.233
  without PER1 and NR1D1). This cohort was seen before the registration (V25.5): it is a robustness test, not a
  replication. A fresh cohort with recorded resection times is the open step.
- Breast, why most patients look early (post hoc, `output/living/breast_calibration.csv`): the 13:00 rule assumes the
  cohort's mean hour; the recorded mean is 13.8 h, so every estimate starts 0.86 h early (64 % underestimated). The
  internal time is also compressed (slope 0.59 on the recorded hour): late patients come out about 2 h early. The order is
  right (Spearman 0.558); with the centre fixed, the per-patient error (1.12 h) equals the full-day adipose biopsies (1.11 h).
- Compartments: breast = two bulk compartments per patient (tumour, adjacent normal), internal time = phase of their shared
  host. Heart = four snRNA-seq pseudobulks (cardiomyocyte, fibroblast, endothelial, neuron; anchors fibroblast,
  cardiomyocyte, neuron). Adipose = one bulk biopsy, no host model: the biopsy's phase on the GTEx clock plane.
- Heart: the truth is a 4-h window per patient, so the error is on the group displacement, not per patient. Holds in
  each surgical arm (p 0.032, 0.033), with fixed anchors and with a left-ventricle waveform. Archive result known:
  confirmation on current code.
- Adipose: muscle (GSE108539) fails on the registered statistic; it passes only with an estimator chosen after the
  failure (Fisher-Lee 0.576), so the muscle pass is post hoc and not shown.
- Status: A200 H200.1-2 registered, pass (seen data); A174 V174.1-3 pass (confirmation); A113 H113.1 adipose
  registered, pass.

**Source:** `exploratory/A200_breast_surgery_rerun/output/{tests,secondary}.csv`;
`exploratory/A174_v3_heart/output/v3_scores.csv`; `exploratory/A113_living_time/output/a113_summary.csv`.

---

## R1.4 A primate sampled around the clock: 63 tissue offsets

**Message:** in a baboon atlas sampled every 2 h over a day, the untimed tissue offsets match the timed ones within
half an hour.

**Bullets**
1. 63 tissues, untimed offset against the timed offset: circular r (Fisher-Lee) 0.912, median error 0.45 h
   (peripheral tissues 0.38 h, neural and ocular 0.50 h).
2. The untimed host phase orders the 12 sampling times: R 0.987, median error 0.39 h.

**Figure:** `figures/baboon_offsets/baboon_offsets.png` (`scripts/baboon_compute.py`, `scripts/baboon_offsets.py`; wide).
A. organ map as the GTEx one: untimed offset (angle) and relative amplitude (radius) per tissue, coloured by organ system,
grey rings: the timed values. B. 63 tissues, untimed vs timed offset after one common shift; dashed: identity. C. the 12
animals' internal time (host phase) vs zeitgeber time after one common shift. Tissue codes: Mure 2018
(`figures/baboon_offsets/tissue_mapping.jpeg`).

**Speaker notes**
- Archive amendment (A36), re-run on current code with leave-one-out anchors (median error 0.455 → 0.452 h). The
  registered statistic (Jammalamadaka-SenGupta on the host phase) is degenerate for an even 12-point grid: the numbers
  here are descriptive.
- One animal per time point: the 12 host phases are 12 animals.
- In 8-h windows of the same atlas (12 windows of 5 time points): order |Spearman| ≥ 0.9 in every window, direction
  right in 11 of 12; offsets degrade to median r 0.55, median error 1.2 h (A40).

**Source:** `results/figure_1/output/species/a36_readouts.csv`; `results/figure_1/output/windows/baboon_windows_summary.csv`.

---

## R1.5 Positive controls: the readout moves when the clock is broken or shifted

**Message:** removing the liver's clock decouples the liver, and only the liver; shifting the liver with daytime
feeding moves its untimed offset as the timed one moves.

**Bullets**
1. Hepatocyte Bmal1 knockout vs control (24 vs 24 mice, ad libitum): liver coupling 0.967 → 0.676, d* −0.287
   (p 0.001). Same permutations: white fat −0.009, lung +0.019; the liver drop exceeds both (p 0.001).
2. Night to day restricted feeding (Xin, 4 tissues): untimed offset changes track the timed ones, Fisher-Lee 0.988
   (exact p 0.042, the minimum attainable with 4 tissues). Liver moves +8.8 h untimed, +9.2 h timed.

**Figure:** `figures/positive_controls/positive_controls.png` (`scripts/positive_controls.py`; half slide). A. change in
coupling at the control's host spread (d*) per tissue, knockout minus control, 95 % interval: liver −0.287 [−0.411, −0.118],
white fat −0.009 [−0.089, 0.055], lung +0.019 [−0.045, 0.087]. B. inverted feeding (food in the day instead of the night): each
tissue's shift relative to the other three, untimed vs timed; dashed: identity.

**Speaker notes**
- The knockout liver keeps a clock-shaped shared state (its gate still passes) and sits −1.31 h from the fat + lung
  host (timed −0.98 h).
- Not shown but stated: daytime feeding in the Manella cohort (3 tissues) cannot reach p < 0.05 by design (minimum
  exact p 0.167; r 0.986); and the registered prediction that daytime feeding lowers liver coupling failed (white fat
  and lung drop as much). The readout detected the shift, not a liver-specific decoupling.
- Mouse controls: the whole day is sampled, so these test the readout, not the hospital-hours design.
- Status: A111/A129 H111.3 registered, pass; H111.1 Xin registered, pass; Manella fail by design; H111.2 fail.

**Source:** `exploratory/A111_feeding_controls/output/{a111_h111_3,a111_h111_1,a111_h111_2,a111_clockloss_liver_offset}.csv`.

---

## R1 what is not claimed

- The disjoint-truth bar for GTEx offsets was missed (0.679 < 0.7); the noise explanation is post hoc.
- Living-patient validation rests on cohorts seen before registration (breast) or known from the archive (heart);
  no fresh timed resection cohort yet.
- Gate passes alone are not claims: the gate makes 13 to 26 false claims in 150 no-clock simulations, the gate plus
  the matched pseudo-clock companion none (A209, A209b; method slide M6).

---

# R2. Organisation in health

**Chapter message:** a person's organs share one clock state, each at its own gain over a departure of constant size,
so the stronger an organ's clock, the more it reads as coupled. Composition does not set an organ's phase; inside tissues
each cell class keeps its gain relative to the others, and lymphocytes carry no readable clock. Age and sex leave the
organisation unchanged; with age the whole body clock runs earlier.

---

## R2.1 Organs share one clock state: each expresses the body's time with its own gain over a departure of constant size

**Message:** every organ carries the body's clock state, at its own gain, plus a departure of its own that is about
the same size in every organ and specific to the clock genes. Coupling is the ratio of the two, so strong clocks read
as the coupled ones: not because they are held more tightly, but because the same departure weighs less against a
larger shared rhythm. The brain follows the same rule.

**Bullets**
1. Across 27 GTEx peripheral tissues, the clock variation shared with the body (S) spans 85-fold (0.12, whole blood,
   to 10.3), the tissue's own departure (D) only 4.6-fold (0.90 to 4.14), and D barely grows with S (slope of log D on
   log S 0.07). Coupling = S / (S + D): stronger clocks are more coupled (Spearman rho +0.734; registered prediction
   from oscillator theory, rho < 0, not supported).
2. Not measurement: a perfectly locked tissue would read a coupling of 0.97 at every tissue's real design (A291).
   The departure is clock-specific: 2.6 times the own variation of expression-matched genes (median over tissues; 90 %
   2.4 to 2.9), above it in all 27 tissues (A294b, registered).
3. The brain follows the rule. In the baboon, where each animal's time is known, coupling rises with each tissue's
   rhythm amplitude against known time (Spearman 0.82 over 57 tissues; periphery 0.77, brain 0.74), and weakly
   rhythmic brain regions fall on the peripheral line (A289: median residual −0.001). Human brain couples at 0.406
   against 0.679 in the periphery.
4. Source of the departure: open. Not individual clock periods (it would grow with amplitude), not shared by
   meal-driven organs (A302, H302.1 fail); noise or other personal cues remain. Calibration against stress-responsive
   genes: A290, running.

**Figure:** `figures/organs_health/organs_health.png` (`scripts/organs_health.py`; wide). Human panels: GTEx fast
deaths, 27 peripheral tissues, one dot per tissue coloured by organ system. S = z2 V_h: between-donor clock-plane
variance of a tissue carried by the body's clock state (gain squared times the host spread); D = 2 psi: the tissue's
clock-plane variance at a fixed body state, net of counting noise, i.e. its own departure from the body's clock state
in amplitude and phase (both in clock-plane log-expression units squared, common to all tissues). A. S against D, log
axes of the same range; grey: Theil-Sen line of log D on log S. B. coupling
against S; grey: S / (S + D) at the median D (1.9), the definition of coupling drawn at a constant departure. C. own
departure of the clock genes against that of expression-matched genes (median of 100 sets); dashed, 1:1; dotted, the
median ratio 2.6. D. baboon: rhythm amplitude against known time (A36 timed gain, estimated without the host; log
scale) against coupling, 35 peripheral tissues and 22 brain regions (pink); grey: the peripheral line.

**Speaker notes**
- Status: A104 registered (H104.1, rho < 0, not supported; the positive relation is a post hoc reading made registered
  by A144 and A271 Part B). A291 registered (H291.1 pass: a perfectly locked tissue reads 0.966 to 0.977 at every
  peripheral tissue's real design; the relation holds after that ceiling, rho 0.717). The S / D decomposition is the
  A291 post hoc table; A294b registered (H294b.1 pass; H294b.2 slope 0.22 of the excess over matched genes, 90 % 0.04
  to 0.33). A302 registered (H302.1 fail). A289 H289.3 registered (fail: brain on the line, the outcome shown).
- Definitions (`clocksync` `_core/host.py`, `on_host`; A104 fit, leave-one-out peripheral host). Each donor x tissue
  pseudobulk is projected onto the reference clock plane: a complex state beta_i with known Poisson covariance V_i.
  The host factor model over the other tissues gives each donor's body clock state, posterior mean f_s and covariance
  S_u. Per tissue k: beta_i = a_stratum + z_k f_s + e_i, e_i ~ N(0, psi_k I + V_i), with |z|^2 debiased by tr cov(z)
  (z2) and psi corrected for its degrees of freedom. V_h = tr(cov f + mean S_u), the spread of the body's state over
  the tissue's donors. S = z2 V_h, D = 2 psi_k, c = S / (S + D).
- What S is: gain squared times host spread, invariant to the host's scale convention (z2 and V_h separately are not).
  V_h is 3.2 to 3.8 in 26 of 27 tissues; skeletal muscle 1.11 (with z2 3.23), so compare S, not z2. Over the 27
  tissues rho(S, z2) 0.90.
- What D is: the scatter of the tissue's clock-plane state around the point the body's state predicts, isotropic
  (psi I), so a departure in amplitude and phase from the body's clock state; "departure from the body's time" is a
  shorthand. It contains (i) the tissue's person-specific departure; (ii) non-count technical variation (dissection,
  composition, RNA quality), not separated by A291; (iii) a small host-error term: regressing on posterior means puts
  about |z|^2 tr(S_u) into the residual; A291 sizes it, a perfectly locked tissue reads c 0.966 to 0.977 (median
  0.972), so about 3 % of S. Counting noise is excluded by construction (V_i modelled). The clock-specific excess over
  matched genes (A294b) argues against generic technical noise, not against technical variation specific to clock
  genes (none known). In hours: A259, SD 1.0 to 3.2 h, a different estimator.
- Panel B is not a test: c = S / (S + D) is how coupling is defined. What the panel shows is that one constant D
  places every tissue near one curve; the finding is panel A (D nearly constant while S spans 85-fold).
- Reading of the slope (A294b correction note, post hoc): for a limit-cycle clock, the angle an input turns it scales as
  1 / amplitude, so a departure independent of amplitude (slope 0) means person-specific cue timing or additive noise;
  a fixed angle would give slope 1, individual intrinsic periods slope 2. The observed 0.07 (0.22 for the excess over
  matched genes) excludes intrinsic periods as the main source. A302 tested the cue-history reading for feeding (the
  liver's departure is not shared more with gut and adipose; differences below about 0.1 undetectable).
- Clock-specificity: the departure exceeds the matched median in all 27 tissues and the matched q95 in 10. The matched
  sets share little between tissues (median S 0.016 to 0.193), so their own variance is read on a small shared state;
  non-count technical variation specific to clock genes (none known) cannot be excluded.
- The departure in hours (A259, registered): SD across donors 1.0 h (thyroid) to 3.2 h (lung), 11 body tissues; an
  upper bound, not above the companion gene sets on that statistic (H259.2). Different estimator from A294b.
- Gene-set companions of the rho (former panel C, A271 Part B, H271.2 pass on the q95 rule): 200 expression-matched
  and 200 stress-responsive gene sets give rho medians −0.160 and 0.109 (q95 0.258 and 0.453). The calibrated form
  (A290) runs the weekend of 2026-11-21; update when its `tests.csv` exists.
- Why panel D uses the timed amplitude (post hoc choice of display): the untimed "own amplitude" counts every
  clock-gene variation between animals, not only the rhythm, and shares a term with coupling. Baboon tissues at
  coupling 0 with sizeable untimed amplitude (lateral hypothalamus, thalamus, testis, pineal; smooth muscle 0.06) have
  weak rhythms against known time (0.29 to 0.76, smooth muscle 1.67; periphery median 2.46): their clock-gene
  variation does not follow the body's clock. Whether no clock or a clock of their own with a phase that differs between
  animals is not tested (underpowered with 12 animals). The baboon (12 animals) was not simulated as in A291.
- The brain's lower coupling is not specific to the clock genes (A263): every gene set drops in the brain, and the
  clock is the program the brain shares most with the body (brain coupling 0.406 against stress-responsive and
  matched set medians 0.064 and 0.029; post hoc). Brain vs periphery as blocks (A142, registered): r* 0.639 [0.530,
  0.728] against random 13-tissue peripheral blocks (median 0.978). Wording: "the brain is less coupled to the body's
  shared state", never "the brain's clock is the exception". In the brain block, own variation, not the shared rhythm,
  varies between regions (S 2-fold, D 4-fold; A291 post hoc; units not comparable with the periphery).
- Not the blood-brain barrier (A289, registered): pineal and pituitary, outside it, couple less than brain regions in
  baboon (0.12 vs 0.42) and sheep (0.08 vs 0.20), the pineal not at all (post hoc reading: its clock is set by the
  SCN's sympathetic input).
- Compared with CHIRAL (Talamanca 2023): CHIRAL gives each tissue a mean amplitude and phase relative to one donor
  phase, assuming every tissue follows it fully. Coupling estimates how fully; D is the spread between people around a
  tissue's mean offset. With donors over the whole day, a two-step CHIRAL pipeline could approach coupling (A145);
  hospital hours and the gene-set comparison are what it lacks.

**Source:** `exploratory/A294_floor/output/{tissues,tests}.csv`; `exploratory/A291_attenuation/output/tests.csv`;
`results/figure_2/output/backbone/a104_gtex_tissues.csv`; `results/figure_2/output/organs/r2_couplings.csv`;
`results/figure_1/output/species/a36_tissues.csv`; `exploratory/A302_departure_groups/output/tests.csv`;
`exploratory/A271_responsive_companions/output/{a104_sets,a104_tests}.csv`; `exploratory/A289_brain_barrier/output/tests.csv`;
`exploratory/A263_companions/output/sets_human.csv`; `exploratory/A142_organ_blocks/output/tests.csv`;
`exploratory/A259_tissue_deviation/output/score/{targets,tests}.csv`.

---

## R2.2 Back to the organs: composition does not set an organ's phase and explains only part of its gain

**Message:** an organ's phase is not a mix of its cells: changing a sample's composition leaves its phase in place,
except in the adrenal. Composition does change gain, but by a few percent per SD, while organs made mostly of the same
cell class differ 1.2- to 2.6-fold: on top of its cell classes, each organ scales how strongly it expresses the body's
time.

**Bullets**
1. Phase (A295, registered, H295.2 pass): one SD more stromal-vascular cells moves a sample's phase by −0.05 h
   (90 % −0.12 to +0.02), one SD more lymphocytes by +0.05 h (−0.03 to +0.13), against organ offsets from −1.2 h
   (testis) to +4.9 h (adrenal) over 26 peripheral tissues.
2. The exception is the adrenal: one SD more cortex moves a sample −1.15 h (90 % −1.87 to −0.56) towards the body's
   time (A299, descriptive arm); part of the adrenal offset may come from its non-cortex part. The medulla markers do
   not carry it (H299.1 fail), nor ACTH response or cortex identity (A301).
3. Gain: one SD more stromal-vascular cells raises gain 1.04-fold (90 % 1.02 to 1.06; descriptive), lymphocytes do not
   detectably lower it (0.99; H295.1 fail); adrenal cortex 1.30-fold. Organ gains span 0.26 (blood) to 2.28.
4. Tissues dominated by the same cell class differ in gain: squamous epithelium 2.6-fold (oesophagus mucosa over
   sun-exposed skin) and 2.1-fold (over unexposed skin), adipocytes 1.6-fold (visceral over subcutaneous adipose and
   breast), smooth muscle 1.2-fold (oesophagus muscularis over sigmoid colon), while a class's gain relative to
   endothelium varies between organs only by a factor 1.10 to 1.36 (R2.3). Reading: an organ-wide drive scales its
   classes together. Not tested directly (A298 closed).

**Figure:** `figures/organs_composition/organs_composition.png` (`scripts/organs_composition.py`; wide). A. phase
relative to the body [h]: top, each peripheral tissue's offset (one dot per tissue, coloured by organ system); below,
phase change per SD of a sample's stromal-vascular or lymphocyte marker score (dots, one per tissue; black, pooled
value with 90 % interval) and of the adrenal cortex score (purple, 90 % interval). B. gain [fold], log scale: top, each
tissue's gain on the body's time; below, gain change per SD of the same scores (fold = 1 + change; pooled with 90 %
interval). C. tissues made mostly of one cell class (squamous epithelium, adipocytes, smooth muscle; post hoc choice):
each tissue's gain on the body's time, one labelled dot per tissue coloured by organ system; right, highest over
lowest gain of the row.

**Speaker notes**
- Status: A295 registered (NOT BLIND for the organ map, A291, A293): H295.2 pass (phase), H295.1 fail (lymphocyte gain);
  the stromal-vascular gain effect is descriptive. A299 registered (H299.1 fail, medulla); the cortex arm is its
  descriptive specificity arm. The same-class pairs are chosen post hoc (A291 post hoc table), the display is post hoc.
- What the within-tissue slopes do not show: composition varies far more between organs than within one (blood
  against liver), so a per-SD slope cannot be extrapolated to say how much of the between-organ gain composition
  explains. The argument for "only part" is the same-class pairs read with R2.2: classes keep their relative gains in
  every organ, yet organs made of the same class differ 1.2- to 2.6-fold. The pairs differ in their minor classes too, so
  this is a reading, not a partition.
- The direct partition (A298: predicted gain from atlas cell shares against GTEx gain) was closed before scoring: the
  marker-based deconvolution explained little of the GTEx profiles (weighted R^2 0.008 to 0.403) and several atlases
  lack families. It needs bulk-calibrated references or several organs of the same donors.
- Within-tissue composition effects per tissue (descriptive): stromal-vascular gain 0.65 (adrenal) to 1.42 (testis)
  fold per SD; the largest phase changes are in whole blood (−0.76 h stromal-vascular, +2.35 h lymphocyte), where the
  marker scores are least meaningful.
- Adrenal: A296 found the +4.4 h offset a whole-clock property (with and without PER1, PER2, NR1D1); in the baboon the
  offset sits in the medulla (+5.9 h), not the cortex (+0.9 h); mouse and sheep whole adrenals show no offset (A303).
  The cortex effect (towards the body's time) and the baboon medulla agree in direction; the human medulla markers do
  not carry it. Leave the adrenal as an open organ-level observation.
- Organ rows: 26 peripheral tissues of the organ map (pituitary in neither group, as R2.1 panel A); one value per
  tissue, no per-donor values.

**Source:** `presentation/output/organs_composition/rows.csv` (from `presentation/output/organs/couplings.csv`,
`exploratory/A295_composition/output/{tissues,tests}.csv`, `exploratory/A299_adrenal_medulla/output/tests.csv`);
`exploratory/A293_organisation_levels/output/part3/families.csv`; `exploratory/A298_organ_drive/output/design.csv`;
`results_preliminary.md` A296, A301, A303.

---

## R2.3 Inside tissues: the organ sets the time; each cell class sets how strongly it expresses it

**Message:** within a tissue, every resident cell class runs at the organ's time, and each class expresses that time
with its own gain, the same in every organ: fibroblasts loudest, endothelium and mural cells next, epithelium and
macrophage-containing myeloid cells softer. Lymphocytes are the zero-gain end of this scale: no readable clock of
their own. The organ also scales the gain of all its classes together (R2.2).

**Bullets**
1. Gain relative to endothelium of the same atlas (registered, A293 H293.4 pass): fibroblast 1.26 [1.07, 1.48],
   mural 0.92, parenchyma 0.88, epithelial 0.80, myeloid 0.67 [0.53, 0.84] (pooled over 4 to 13 body atlases). A
   class's gain varies between organs by a factor of 1.10 to 1.36, less than the same class differs between GTEx
   organs (1.3 to 2.7): the organ differences come from composition and an organ-wide drive, not from classes
   responding differently in different organs.
2. The organ sets the phase: offsets between stromal-vascular classes stay within about 1.5 h; endothelium leads
   fibroblasts by +0.46 h over six atlases (90 % 0.09 to 0.90), the mouse sign, but not in a consistent order
   (H293.3 fail): a lead, not a result. A sample's composition moves its phase by at most 0.13 h per SD (A295, H295.2
   pass); it moves only how loudly it expresses the body's time.
3. Lymphocytes share no more of the tissue's clock state than random genes (stromal-vascular minus lymphocyte sharing
   beyond both companion families: five atlases 0.277, replication 0.152, IBD cohort 0.413, living donors 0.225;
   p_cal 0.004 to 0.020).
4. No clock of their own: lymphocyte clock genes are coherent between donors in 2 of 9 atlases, fibroblasts of the
   same donors thinned to the lymphocytes' depth in 7 of 7 powered atlases (A281); where a clock is readable at their
   depth, no clock-gene structure beyond random genes (A287): "no readable clock on the full clock-gene set".

**Figure:** `figures/celltypes_health/celltypes_health.png` (`scripts/celltypes_health.py`; wide). A. gain of each
cell class relative to the endothelium of the same atlas (ratio of shared-rhythm amplitudes, log scale): coloured
dots, one per atlas; black, random-effects pooled value with 95 % interval; lymphocytes: no gain estimated, the arrow
marks the zero-gain end. B. phase offsets between stromal-vascular classes (positive = first class earlier): dots,
human atlases (untimed pairwise); black, circular mean over atlases with atlas bootstrap 90 % interval; open
diamonds, mouse colon and aorta against known time. C. own clock: signed clock-gene coherence between donors of the
lymphocytes (cyan) and of stromal-vascular cells of the same donors thinned to the lymphocytes' clock-gene counts
(brown, median of 20); grey bar, 5-95 % range of the lymphocytes' 200 expression-matched gene sets (A294), nine
atlases.

**Speaker notes**
- Status: A293 Parts 2 and 3 registered (NOT BLIND: the per-family gains and the mouse offsets were seen in earlier
  work); H293.4 pass (tau <= 0.35 in all five families), H293.3 fail (4 of 6 atlases with the mouse sign, sign test
  p 0.344); A295 registered (H295.2 pass, H295.1 fail).
- The organ factor of gain is inferred, not tested: relative gains cancel any organ-wide drive, and the direct test
  (A298, composition against organ gain in GTEx) was closed before scoring because the atlas shares were implausible.
  Say "composition and an organ-wide drive", not which one.
- Heterogeneity: epithelial (tau 0.305) and myeloid (0.302) are close to the bound and Q shows spread beyond the
  bootstrap SEs (p 0.001, < 0.001); "myeloid" mixes macrophages with monocytes and dendritic cells (only macrophages
  carry a clock, A288), which may explain part of it. Atlases have one or two compartments per family; the SE ignores
  the covariance with endothelium.
- Lymphocytes are not given a number on the gain axis: with no readable clock, their "gain" would be noise. The arrow
  says where they sit on the scale, the evidence is panel C. Do not suggest a function for this.
- Phase, mouse calibration (timed): untimed pairwise offsets between cell types match known time in the colon (A273:
  median error 0.31 h, slope 1.005) and the aorta (A279: 0.30 h, slope 1.04). In mouse the order endothelial >
  mural > fibroblast holds in both tissues; the human pooled endothelial lead is about a third of the mouse size.
  The breast endothelial-mural offset (+3.2 h) is one compartment pair.
- Composition (A295): stromal-vascular content raises gain (+0.041 per SD, positive in 21 of 27 tissues;
  descriptive), lymphocyte content does not detectably lower it (H295.1 fail: within-tissue variation too small).
- Panel C, read honestly (A294, registered): the grey bar is the 5-95 % range of random expression-matched genes in
  the lymphocytes. Lymphocytes sit inside it in seven of nine atlases. Breast CD4 T (0.38) is just above it, and is a
  hormonal pulse: without PER1, PER2 and NR1D1 it is not significant (p 0.081) and it is not above stress-responsive
  genes (H294.1, H294.2 fail). Heart lymphocytes (0.50, 36 donors) keep a coherence without those genes (p 0.014) and
  above stress-responsive genes (both pass), but not above random genes: a weak signal in one atlas.
- Caveat, post hoc (A294): on the core loop without the hormone-induced genes, lymphocyte coherence rises in most
  atlases (6 of 9 at p < 0.10, e.g. IBD ileum −0.43 to 0.40). The hormone-induced genes may mask a weak lymphocyte
  core loop; "no readable clock" holds on the full clock-gene set, and a blind core-loop test on fresh atlases is
  needed before the slide states it without qualification.
- Sharing per cell family (post hoc, descriptive; median |r| over compartments of ten body atlases; the panel of the
  previous version): fibroblast 0.63, parenchyma 0.60, mural 0.58, epithelial 0.58, endothelial 0.58, glia /
  melanocyte / progenitor 0.55, myeloid 0.48, lymphocyte 0.22 (random genes 0.26 to 0.37 for resident classes, 0.21
  for lymphocytes). Every resident class shares the clock state; immune cells are the exception. Macrophages join
  the resident cells: macrophage minus lymphocyte 0.148 (p_cal 0.003 / 0.006, 8 atlases; A288).
- Readout: `clocksync.pairwise`, anchor-free: for two cell classes of the same donors, the complex correlation over
  donors of their clock-plane states (|r|). No host, no anchor choice, so the A17 anchor question does not arise.
- Status: A272 (H272.2), A274 (H274.1 replication, blind for its atlases), A280 (H280.1, blind), all restated on the
  calibrated test of A283 (5000 placebo draws; rejection 0.031 to 0.063 under no pairing); A282 (living arm,
  descriptive question fixed in advance); A288 (H288.2 pass). The earlier q95 rule was anti-conservative (A276); no
  verdict depended on it.
- What is clock-specific is the contrast. The level of sharing is not shown above stress-responsive genes everywhere:
  stromal-vascular sharing passes both families in the A274 atlases (0.478; p_cal 0.009 / 0.013) but not against
  stress-responsive genes in the A272 atlases (0.499 vs median 0.434, p 0.299); macrophage sharing likewise (0.409 vs
  0.367, p 0.133). A shared body-wide stress state carries part of the level; the gap to lymphocytes is clock.
- Living donors: in HLCA, the one atlas with both arms, stromal-vascular sharing is the same living and post-mortem
  (0.523 vs 0.516): the backbone is not an agonal artefact.
- Not residency (A284, registered, none passes): resident T cells do not share more than circulating ones (Delta
  0.010), plasma cells share less than B cells in every atlas, alveolar over monocyte-derived macrophages misses
  (A288 H288.3, p 0.051).
- Cortex: the cell classes share a clock state pooled over three cohorts (A283 H283.5 pass) but not in the largest,
  PsychAD MSSM (A280 H280.2 fail): not claimed.

**Source:** `exploratory/A293_organisation_levels/output/part3/{atlases,families}.csv`;
`exploratory/A293_organisation_levels/output/part2/{atlases,tests}.csv`; `exploratory/A295_composition/output/tests.csv`;
`results_preliminary.md` A293 (GTEx scale, A291 post hoc table); `exploratory/A272_backbone_pairwise/output/pairs.csv`; `exploratory/A274_celltype_map/output/pairs.csv`;
`exploratory/A276_pairwise_calibration/output/a283/{per_atlas,tests}.csv`;
`exploratory/A280_celltype_replication2/output/calibrated/tests.csv`; `exploratory/A282_pooled_living/output/{tests,per_atlas}.csv`;
`exploratory/A288_macrophage_pairwise/output/{per_atlas,tests}.csv`; `exploratory/A284_residency/output/tests.csv`;
`exploratory/A281_lymphocyte_clock/output/{pairs,tests}.csv`; `exploratory/A287_lymphocyte_covariance/output/tests.csv`;
`exploratory/A294_lymphocyte_hormonal/output/{pairs,tests}.csv`.

---

## R2.4 Age and sex: the core clock's organisation is the same; ageing shifts the body's clock as a whole

**Message:** as CHIRAL found, the core clock is the same in old and young, women and men: here, coupling and organ
offsets are equivalent to no difference. The reported changes in output rhythms (fewer rhythmic genes in older donors,
more in women) largely follow the type of death and shrink or vanish when it is matched. What changes with age is the
timing of the whole body clock: it runs earlier, the brain with the body.

**Bullets**
1. Output rhythms (CHIRAL's statistic, re-read with its TMM normalisation; log ratio of rhythmic-gene counts): old over
   young −0.67 on CHIRAL's donor subsample, −0.58 on all donors, −0.20 matched on death type (95 % −0.98 to 0.51),
   0.00 with RNA quality also matched; female over male +0.61, +0.51, +0.23 (−0.08 to 0.60), +0.17. Young GTEx donors
   are mostly violent fast deaths.
2. Core-clock organisation (registered): coupling old − young +0.003 (90 % −0.027 to 0.033), female − male +0.020
   (−0.003 to 0.043), both equivalent within ±0.15; organ offsets −0.05 h and +0.04 h, equivalent within ±1 h (A231,
   A208), in fast and ventilated donors.
3. Timing (registered, timed readout): internal time runs ahead of the recorded time of death by +0.64 h per decade in
   the body (90 % 0.17 to 1.11) and +0.99 h in the brain (0.49 to 1.49), death type and sex held fixed; brain and body
   do not drift apart (−0.01 h per decade, A234).

**Figure:** `figures/ageing_sex/ageing_sex.png` (`scripts/ageing_sex.py`; wide). A. CHIRAL's statistic: total log ratio
of rhythmic-gene counts (TMM), old / young (filled circles) and female / male (open diamonds), on CHIRAL's donor
subsample, all donors at equal n, donors matched on death type (Hardy class), and also on RNA quality (RIN tertile);
tissue bootstrap 95 % interval. B. coupling difference at a common host spread (d*), old − young and female − male,
90 % interval; dashed: equivalence bounds ±0.15. C. organ offset difference [h], same; dashed: ±1 h. D. internal time
minus recorded time of death, change per decade of age, body and brain hosts, 90 % interval.

**Speaker notes**
- What Talamanca et al. 2023 claim (re-read, `literature/talamanca_science_2023.pdf`): "The relative phases and
  amplitudes of clock genes were conserved in males and females"; "the overall amplitudes, phases, and relative
  relationships of the clock genes were conserved with age"; the number of rhythmic genes "was higher in females by
  about twofold", and rhythms were "overall strongly damped in the older donors" (< 50 vs > 60 y). Supplement
  (`scCircadianCancer/literature/chiral_supp.pdf`, methods): "we regressed out ischemic time, sex, age and type of
  death as categorical variables with a linear regression model" on log CPM, per gene and tissue; for the sex and age
  contrasts, "when the number of samples in the two groups differed, we subsampled the group with more samples to match
  the sample size of the smaller group". Death type is therefore removed from mean expression but not balanced between
  the groups compared. In the main text it appears only for time-of-death concordance (better in fast deaths, fig.
  S2E); they list "complex covariate structures, or sampling bias" among limits. So panel B agrees with them; panel A
  concerns their output-gene claims.
- Why matching and not covariate removal: regressing death type out of mean expression does not remove its effect on
  rhythm amplitude (fast deaths carry stronger rhythms); matching equalises the death-type mix of the two groups. Our
  re-run keeps CHIRAL's covariate regression (ischaemic time, age, sex, death type) and adds the matching, so the drop
  in panel A is what matching adds on top of their correction. Reading of the mechanism, not a separate test.
- Status of panel A: confirmation of archive A101b / A102b (WP0 item 5, 2026-10-23), CHIRAL's statistic, not
  clocksync. Is the drop itself significant (matched minus pooled, same tissues, tissue bootstrap)? Sex: −0.27 (95 %
  −0.39 to −0.14, p < 0.001); with RNA quality −0.34 (p < 0.001). Age: +0.39 (−0.10 to 0.84, p 0.092), +0.59 with RNA
  quality (0.02 to 1.25, p 0.040). Age matching keeps 29 of 37 tissues. With CPM, the age drop is +0.57 (p 0.024).
  Wording: "with death type matched, the age dampening is no longer seen and the female excess halves and is not
  significant"; not "CHIRAL's results are artefacts". Their within-tissue findings (e.g. coronary artery, liver
  programmes) are not re-tested here.
- Panels B, C: A231 (registered, H231.1 pass; H231.2, H231.3 fail: no shift detected, organ-specific shifts with SD
  above 0.45 h would have been detected); A208 (registered, H208.1 pass, H208.2 fail; offset equivalent within about
  ±0.2 h; pseudo-clock |offset| q95 1.29 h, so no clock-specific offset). MDEs: d* 0.037 (age), 0.032 (sex). The old
  arm carries more T2D (not adjusted; T2D leaves coupling unchanged, A164 / A202). A shared amplitude loss would leave
  coupling and offsets unchanged by construction: B and C do not test dampening.
- Panel D: A234 (registered, H234.1 and H234.3 pass; H234.2 post hoc). A timed readout (recorded time of death), fast
  deaths only. Likely inflated (winner's curse): power about 0.29 at the archive prior 0.35 h per decade, which the
  interval includes; young fast donors are few (19 under 40 y). Age is not associated with time of death (post hoc,
  R 0.05). Direction as chronotype epidemiology (Roenneberg 2004; to verify before citing). Sex and internal time
  versus the clock was not tested here (A99 archive, timed, +0.35 h female − male, SE 0.31, post hoc).
- Host spread with age (A231 D231.6, descriptive): no consistent loss of the shared clock's spread in fast donors
  (median log ratio −0.03), confounded by time of death.

**Source:** `exploratory/WP0_integrity/output/chiral_tmm/{count_comparison_tmm_summary,sex_comparison_matched_summary}.csv`;
`exploratory/A231_age_coupling/output/{score_tests,score_d231_6_spread}.csv`;
`exploratory/A208_sex_coupling/output/{tests,pooled}.csv`; `exploratory/A234_internal_time_age/output/score_tests.csv`;
`literature/talamanca_science_2023.pdf`; `scCircadianCancer/literature/chiral_supp.pdf`.

---

# R3. What holds the organisation

**Chapter message:** the organ arrangement does not need the daily cues: it survives slow deaths and ventilation
without meals or activity, a peripheral clock keeps its amplitude when its main cue is removed, and a shifted behaviour
moves the outputs and barely the clock. It is lost when the brain is lost. The central drive, not the external day, holds the organs together (a hypothesis for the brain part:
GTEx is the only data).

---

## R3.1 The organ arrangement survives the loss of external time, not the loss of the brain

**Message:** in slow deaths and in ventilated donors the body's clock no longer reads the hour of death, yet the organs
keep the arrangement of fast deaths, on the clock genes only. Brain-dead organ donors, ventilated in the same way, lose
it.

**Bullets**
1. External time is lost: the body's clock reads the recorded hour of death in fast deaths (R 0.63, 244 donors) but
   barely in slow deaths (0.16, p 0.079, 106) and ventilated donors (0.14, 380; post hoc readout).
2. The arrangement is kept: slow deaths (Hardy 4, 103 donors) share the fast-death organ arrangement at Fisher-Lee rho
   0.784, above all 200 expression-matched (median 0.038) and all 200 stress-responsive gene sets (median 0.234)
   (registered, H293.1); ventilated donors with a living brain (79), without meals or activity, 0.690 (H293.2).
3. Brain-dead organ donors (79, same ventilation and intensive care): 0.161, not different from tissue-label
   permutation (p 0.18), at the level of random genes (descriptive).

**Figure:** `figures/holds_organisation/holds_organisation.png` (`scripts/holds_organisation.py`; wide). Rows in the
same order in both panels: fast deaths (reference), slow, intermediate, ventilated with a living brain, brain-dead
organ donors. A. correlation of each group's organ-pair phase differences with those of fast deaths (Fisher-Lee, 12
body tissues): bars for the clock genes (blue) and the medians of 200 expression-matched (light grey) and 200
stress-responsive gene sets (dark grey), whisker to their 95th percentile, same donors and pairs; fast deaths are the
reference (1 by construction, outlined). B. alignment of the body's clock (host phase from the clock genes of 15
reference tissues, fitted without time) to the recorded hour of death, R = |mean exp(i (host phase − ω t_death))|, per
death type, with its permutation p and donor count; the ventilated readout pools living-brain and brain-dead donors
(one light bar across both rows).

**Speaker notes**- Status: A293 Part 1 registered, NOT BLIND (the clock rho of every group was known from A109 / A132; the companion
  sets were new); H293.1 and H293.2 pass. Hardy 3 and brain-dead donors descriptive by decision. Brain death vs
  ventilated as a difference: 0.517 (p 0.008, A132b), post hoc.
- What the groups remove (A109c, A132): slow deaths are not aligned to the hour of death (R 0.097, null median 0.101;
  entrainment as strong as in ventilated donors not excluded); their arrangement survives removing the wall-time
  component (rho 0.750). Ventilated donors and brain-dead donors share intensive care, ventilation and the absence of
  meals and activity; they differ in brain death. Both are pushed toward 0 / 12 h relative to fast deaths (A132).
- Not the drugs (A136): donor-management glucocorticoids do not explain the loss (exposure not robustly different,
  no dose-response at MDE 0.52), and the loss stays without PER1 and PER2 (0.726 vs 0.268, Holm p 0.036). The gate
  failure is not robust when long-term steroid users are removed: a caveat.
- Why "hypothesis" for the brain: GTEx organ donors have no brain samples, so the SCN cannot be read; the registered
  SCN test was impossible (PI, Q7, 2026-10-24: the SCN reading stays a hypothesis). Consistent with SCN lesions
  destroying inter-tissue phase synchrony in mice (Yoo 2004); contrary report: Tahara 2012. Wording: "organs keep a
  shared internal time unless the central drive is lost".
- Panel B (A27, post hoc readout, reused as context in A46, A92, A93): one host for all four groups. The registered
  test for slow deaths (A109c H109c.1, 12 body tissues, 103 donors) gives R 0.097 against a null median 0.101 (p 0.56):
  no detectable alignment, but alignment as strong as in ventilated donors is not excluded (80 %-power R 0.208).
  Ventilated donors keep a small alignment (R 0.14, p 0.0005) and are pushed toward 0 / 12 h (A132); the brain-dead
  and living-brain subgroups are not separated in this readout. Intermediate deaths: R 0.25 (p 0.039, 51 donors).

**Source:** `exploratory/A293_organisation_levels/output/part1/{groups,sets,reproduction}.csv`;
`results/figure_3/output/cues_lost/walltime_hosts.csv`; `exploratory/A109_coupling_nulls/output/a109c_walltime.csv`;
`results_preliminary.md` A109c, A132, A136; `docs/PI_QUESTIONS.md` Q7.

---

## R3.2 A peripheral clock holds without its cue; its outputs follow behaviour

**Message:** removing or shifting one cue in one organ leaves its core clock in place: without its feeding rhythm the
mouse liver clock keeps three quarters of its amplitude, and after a 12 h shift of behaviour the human muscle clock
stays near central time, damped, while its rhythmic outputs follow behaviour.

**Bullets**
1. Mouse liver (A297, registered): removing the feeding rhythm keeps 0.72 of the clock's amplitude (90 % 0.66 to 0.79)
   and moves it 1.6 h; if the clock were the sum of its cues it would fall to 0.15 to 0.29.
2. Human muscle, three days of 12 h-shifted sleep, light, activity and meals (11 men): the core clock moves 2.7 h
   (90 % 1.8 to 3.6) and its daily change falls to 0.23 (A257); the 1,587 rhythmic output genes move 11.6 h (10.5 to
   12.8) at 0.85 of their change (A261): the clock and its outputs come apart.

**Figure:** `figures/holds_organisation/cue_removal.png` (`scripts/holds_organisation.py`; wide). A. mouse liver:
core-clock amplitude relative to ad libitum feeding when the feeding rhythm is removed (blue, 90 % bootstrap interval,
phase shift beside it), and the amplitude expected if the clock state were the sum of its cues (grey open; calibrated
on an inverted feeding rhythm in two datasets, Manella upper, Xin lower). B. human muscle: shift of the 07:00-to-19:00
change after a 12 h shift of behaviour against the aligned condition (dashed: 0, central time; 12 h, behaviour), 90 %
interval, for the core clock (14 clock genes, blue) and the rhythmic outputs (black); right, the size of the change
relative to aligned.

**Speaker notes**
- Mouse (A297, registered, NOT BLIND for the A111 shifts; H297.1 fail, above the interval): different laboratories, strains and protocols; liver only; amplitudes read in each dataset's own
  projector. The additive prediction (1 + rho) / 2 assumes the clock's state is the sum of its cues. The
  night-restricted control is 1.04 (0.93 to 1.14). Reading for humans: removing one cue (as in intensive care, A46)
  leaves a peripheral clock running near its phase; conflicting cues move it.
- The inversions behind panel A's prediction (registered calibration): day feeding against ad libitum (Manella) moves the
  liver clock 11.7 h at 0.71 of its amplitude; day- against night-restricted feeding (Xin) 11.3 h at 0.44.
- Panel B (A257 registered, NOT BLIND, H257.2 partial shift, H257.3 damping pass; A261 registered, H261.3
  "behaviour"): within-person, Wefers 2018 (GSE106800), biopsies at 07:00 and 19:00, day 3 of each protocol. The core
  clock's 12 h ambiguity: a 12 h shift has no direction of its own, so "towards behaviour" cannot be read from the
  sign; 10 of 11 men shift positively. Masking limit (A261): with two biopsy times, a reversed 19:00 - 07:00 difference
  of the outputs cannot be told from an acute response to fasting and waking (both biopsies are at the behavioural
  morning and evening in both conditions); a full cycle under a constant routine would separate them. Say "outputs
  follow behaviour (shifted or evoked)". Wefers 2018 gave the same reading (clock genes kept their pattern, fatty acid /
  PPAR genes changed); A257 / A261 add hours, intervals and companions. Light moved with behaviour, so "central time"
  is the expected position of a clock that did not follow the 12 h shift, not a measured SCN phase.

**Source:** `exploratory/A297_feeding_model/output/{tests,arms}.csv`; `exploratory/A257_nightshift_muscle/output/score_tests.csv`;
`exploratory/A261_nightshift_programs/output/programs.csv`; `results_preliminary.md` A257, A261, A297.
