# Results section: slide by slide

Order (agreed 2026-11-11): R1 model validation, R2 organisation in health (with ageing), R3 development, R4 loss of
daily cues, R5 chronic disease, R6 where clocks break: locally (psoriasis, cancer), R7 misalignment and limits,
R8 synthesis. Every number below was re-read from the table named in its slide's **Source** line. GTEx figures show
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
- Truth (i) shares donors with the test set (the window donors are part of the cosinor); truth (iii) does not, which is
  why it was registered as the stricter test.
- The two truths agree with each other at only 0.821 [0.657, 0.925]: the ceiling of any comparison with truth (iii)
  is low.
- Donor order within the window (tangent coordinate vs recorded time): Spearman 0.613, direction correct.
- Status: H123.1-3 registered, pass; H123.5 registered, **fail** (offsets); H123.6 registered, pass; split halves post
  hoc (PI decision 2026-10-28: keep the panel with the verdict as scored).

**Source:** `exploratory/A123_working_hours/output/tests.csv`; `a123c/tests.csv`; `a123c/diag/estimate_split_half.csv`,
`a123c/diag/split_half.csv`.

---

## R1.2 Orderers need the whole day; this method does not

**Message:** untimed orderers lose the organ map when samples come from working hours only; with the whole day, the
ones built on clock genes recover it. The method keeps the map in both designs.

**Bullets**
1. Working hours (111 donors), offsets vs truth (i), circular r: clocksync 0.945; CHIRAL −0.165, CYCLOPS 2.0 −0.296,
   CoFE genome-wide −0.283, CoFE clock genes −0.417. Each difference with clocksync Holm p 0.0015 (registered).
2. Whole day (all fast deaths): CHIRAL 0.877, CoFE clock genes 0.856, clocksync 0.966. The failure is the design, not
   the orderers.
3. In working hours the method also beats the supervised TimeTeller (0.754) and a timed oracle using the window
   donors' recorded times (0.738).

**Figures**
- `figures/gtex_sweep/gtex_sweep.png` (`scripts/gtex_sweep.py`): median offset and amplitude error per method against the
  number of donors (10-80), whole day vs 08-18 h, from the Figure 1 GTEx benchmark (8 tissues, 10 draws per n; truth:
  cosinor on the recorded times of all fast deaths). In hospital hours clocksync's offset error stays 0.35-0.56 h from
  10 to 80 donors; every orderer, TimeTeller and the timed oracle start near 1.6-2.2 h at 10 donors. Amplitudes: all
  methods comparable, except the timed oracle in hospital hours (a 10 h window cosinor).
- `figures/gtex_dip/gtex_dip.png` (`scripts/gtex_dip_compute.py`, `scripts/gtex_dip.py`): CHIRAL's published donor internal
  phase (DIP) against clocksync's internal time, binned (3 h, bins with < 5 donors hidden), clocksync fitted on the full day
  (241 donors; median |difference| 1.7 h) or on hospital hours only (110 donors; 1.8 h). Panel C: median absolute error
  against the recorded time of death, donor bootstrap 95 % CI. Full day: clocksync 2.03 h [1.71, 2.31], DIP 2.58 h
  [2.19, 2.91]. Hospital hours: clocksync 1.71 h [1.19, 2.22], DIP 2.93 h [2.28, 3.43]. Within hospital hours the DIP
  barely orders the donors; clocksync's hospital-hours internal time does.

**Caveats for the figures:** the sweep and the DIP comparison are post hoc and descriptive (presentation figures). clocksync's
sweep runs were refitted on current code (c4dde99, leave-one-out anchors; 2026-11-14); the other methods' runs are
unchanged from 2026-10-25. CYCLOPS and CoFE genome-wide were run only at 80 donors (3 draws) for cost.

**Speaker notes**
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
A. breast: internal time (13:00 rule, no recorded time used) vs resection hour, one dot per patient, equal axes. B. heart:
each patient's internal time by surgery window (06-10 h vs 14-18 h), circular means. C. per-sample absolute error after
removing one common offset: breast (08-18 h) 1.12 h, adipose (full day) 1.11 h.

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

**Figure (to build):** untimed vs timed offset per tissue on a clock face, peripheral and neural tissues in two colours.

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

**Figure (to build):** A. coupling per tissue, control vs knockout, paired dots (liver falls, fat and lung flat).
B. untimed vs timed offset change per tissue under day feeding, the diagonal.

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
