# Results section outline (2026-11-23, built slide by slide with the user)

Prompt: `docs/prompts/presentation_results_2026-11-23.md`. Phase convention on slides (user, 2026-11-23): peak time,
+ = later, as CHIRAL Fig. 1B (clocksync offsets, > 0 advanced, are shown negated); offsets between two compartments
state their direction in words. Framing: "loudness, not time" (`docs/ABSTRACT_DRAFT_2026-11-22.md`,
`docs/PAPER_PLAN_2026-11-19.md` section 11). Slides kept from `results.md` are referenced, not repeated. Tags: blind / not
blind / descriptive / post hoc. Every number below was read from the output table named in its **Source** line.

---

# 1. Validation: the readout sees time, and sees it move

| slide | status | figure |
|---|---|---|
| 1.1 GTEx organ map from working-hours deaths | kept (`results.md` R1.1) | `figures/gtex_organ_map/` |
| 1.2 Living patients: breast, heart, adipose | kept as is (`results.md` R1.3; user, 2026-11-23); speaker note below | `figures/living_patients/` |
| 1.3 Baboon atlas, 63 tissue offsets | kept (`results.md` R1.4; user, 2026-11-23) | `figures/baboon_offsets/` |
| 1.4 Timed living skin and muscle, and hospital hours | **new** (below) | `figures/living_timed/` |
| 1.5a Positive controls, cell types: offsets and the knockout | **new** (below) | `figures/celltype_controls/` |
| 1.5b Positive controls, organs: the readout moves when the clock is broken or shifted | **revised** (below) | `figures/positive_controls/` |
| 1.6 Negative control: tissues of one living person share one time; her cell line does not | **new** (below) | `figures/living_tissues/` |


**1.2 speaker note (within-person rotation, not shown):** between two biopsies of one day the state moves along the
clock circle in the protocol's direction (13 of 14 cells within 1 h, A252; 17 of 20 within 2 h, A256 blind replication),
beyond matched and stress-responsive genes. Source: `exploratory/A252_elapsed_time/output/score_datasets.csv`,
`exploratory/A256_elapsed_replication/output/score_cells.csv`.

## 1.4 In timed living skin and muscle, the untimed readouts match the timed ones, and survive hospital hours

**Bullets**
1. Offset between compartments of one person: skin dermis minus epidermis, clocksync -1.58 h (90 % -1.81 to -1.34; 148
   untimed women) against timed -1.98 h (-3.15 to -0.95; 20 men of the same study) (A251 H251.4, pass). Exposed minus
   protected skin site: 1.33 h against 1.57 h timed (A342).
2. Pre-mRNA lead over mRNA: skin protected 1.19 h against 1.18 h timed, exposed 1.41 against 1.48 (A342 H342.1-2, pass);
   muscle 1.42 against 1.75 (A337 Perrin H337.1-2, pass on 14 genes).
3. Hospital hours (same skin, 06, 12, 18 h only): the lead stays within the bound at both sites (untimed minus truth
   protected +0.10 h, -0.44 to 0.63; exposed -0.09, -0.57 to 0.53; H344.2 pass); the offset keeps its value (1.37 against
   1.57 h) but its interval is too wide to decide (difference -0.20, -1.68 to 0.88; H344.1 undecidable). With two
   daytime points (06, 12 h) the protected lead drifts +0.85 h (0.14 to 1.64; H344.3 fail).
4. Hospital hours in muscle (08, 12, 16 h only): the lead is recovered, untimed minus truth -0.01 h (-0.50 to 0.63;
   A353 H353.1 pass); also on 08-20 h (-0.35, -0.57 to -0.10) and on two points, 08, 12 h (+0.16, -0.55 to 0.89; H353.2
   pass).

**Figure:** `figures/living_timed/living_timed.png` (`scripts/living_timed.py`; wide). A. offset between compartments,
clocksync (blue circles) against timed (black diamonds). B. pre-mRNA lead, same encoding. Filled: 14 clock genes; open:
core loop. C. skin (A344) and muscle (A353) cut to subsets of clock times (full day, three day times, two day times):
untimed minus the full-cycle timed truth, skin lead at each site, muscle lead and the skin offset; grey bands, registered bounds (lead +-1.5 h, offset +-1 h). 90 % intervals throughout.

**Tags:** A251 not blind (registered, one cohort); A342 registered, scored on direct introns after the erratum; A337
Perrin blind for the unspliced counts; A344 and A353 not blind for the full cycle, blind for every subset value.

**Caveats (speaker notes)**
- The exposed-protected pairing is not clock-specific against companions (H342.3 fails on that condition only;
  coherence 0.467 against matched q95 0.549): the offsets agree, the pairing could be shared visit-level variation.
- Dermis-epidermis compares different people: untimed women against timed men.
- Core loop: muscle lead under-read by about 1.3 h (1.63 against 2.97 h; H337.2 core FAIL). Skin hospital hours on the
  core: pass on 06-18 h (leads -0.13, -0.48 h); on 06-12 h the exposed lead also fails (-1.16 h).
- Reading each visit as a time of day degrades as the arc shrinks (median error 2.03 h full day, 2.59 h 06-18 h, 3.27 h
  06-12 h; night pair 3.79 h; descriptive): the method reports relationships, not sample times.
- The H344.1 interval is wide largely because the timed truth itself is (full set -1.00 to 0.51).
- Muscle passes on two daytime points where skin does not: six full-day times and a well-entrained protocol (hourly
  isocaloric meals, controlled light), the most favourable design. Its per-biopsy time error also grows as the arc
  shrinks (1.42 h full day, 2.51 h on 08-16 h, 4.00 h on 08-12 h; A353, descriptive). Muscle has one compartment: no
  offset there.

**Source:** `exploratory/A251_skin_sites/output/{score_offsets,score_descriptive}.csv`;
`exploratory/A342_skin_timed/output/score_direct/{pairing,leads}.csv`;
`exploratory/A337_timed_lead/output/perrin/verdicts.csv`; `exploratory/A342_skin_timed/output/score_a344/subsets.csv`;
`exploratory/A337_timed_lead/output/a353/subsets.csv`; `docs/results_preliminary.md` A251, A337 Perrin, A342, A344,
A353.

## 1.5a Between cell types of one sample, untimed offsets equal timed ones; without a clock they stop being clock-specific

**Bullets**
1. Mouse colon (A273): offsets between cell classes of one sample, clocksync against timed, median error 0.31 h over 43
   pairs (undamaged 28, regenerating after injury 15), on the identity line also where the classes move apart.
2. Mouse aorta (A279, Auerbach 2026): median error 0.30 h over 6 pairs of EC, SMC, fibroblast, macrophage.
3. Clock intact, every aorta pair correlates above both companion families: clock |r| 0.83-0.98 against matched q95
   0.58-0.78 (A354 H354.2 pass, 6 of 6).
4. Bmal1 knockout: macrophage pairs fall to the companion median (|r| 0.11-0.35); stromal pairs fall to the top edge of
   the matched range (EC-SMC 0.629 vs q95 0.636; fib-SMC 0.661 vs 0.666), EC-fibroblast stays above (0.753 vs 0.618;
   H354.1 fail).

**Figure:** `figures/celltype_controls/celltype_controls.png` (`scripts/celltype_controls.py`; half). A. offset between
cell types, clocksync against timed: colon (green filled), regenerating colon (green open), aorta (red); dashed:
identity. B. aorta, |r| between the clock states of two cell types, clock intact (filled, upper row) and Bmal1 knockout
(open, lower row), beside the 5-95 % range and median of 200 expression-matched (light grey) and 200 stress-responsive
(dark grey) gene sets.

**Tags:** A273 registered calibration (not blind for the per-class timed phases); A279 registered (not blind for
A201's per-cell-type shifts); A354 registered, not blind for the clock |r|, blind for every companion value.

**Caveats (speaker notes)**
- One pair above the matched q95 is not proof of a running clock: EC and fibroblasts keep clock-gene covariation beyond
  random genes without BMAL1 (shared non-clock regulation of the same genes). The human cell-type results rest on pooled
  contrasts over many atlases with calibrated p-values (A283), not on single pairs.
- Stress-responsive genes covary strongly between stromal classes (q95 up to 0.94): the clock is separated from them
  only with the clock intact.
- Mouse, timed designs: calibration of the offset readout, not a test in human tissue.

**Source:** `exploratory/A273_pairwise_offset_truth/output/pairs.csv`; `exploratory/A279_aorta_pairwise_truth/output/pairs.csv`;
`exploratory/A279_aorta_pairwise_truth/output/a354/{pairs,verdicts}.csv`; `docs/results_preliminary.md` A273, A279, A354.

## 1.5b The readout moves when a clock is broken or shifted (mouse organs)

**Bullets**
1. Liver Bmal1 knockout (Manella 2021, ad libitum, 24 vs 24 mice): liver coupling falls by 0.287 at the control's host
   spread (95 % -0.411 to -0.118, p 0.001); white fat (-0.009) and lung (+0.019) do not move (A129 / H111.3 pass;
   companions on the same permutations).
2. Daytime feeding moves the liver clock against the other tissues; the untimed offset changes follow the timed ones:
   Xin 2021 (night to day feeding, 4 tissues), Fisher-Lee 0.988, exact p 0.042 (H111.1 pass); Manella 2021 (ad
   libitum to day feeding, 3 tissues), 0.986, exact p 0.167, the minimum with 3 tissues (H111.1 fail by design).

**Figure:** `figures/positive_controls/positive_controls.png` (`scripts/positive_controls.py`; half). A. change in
coupling, knockout minus control, per tissue (95 % interval from the permutations). B. per tissue, offset change after
the feeding switch, clocksync against timed, centred on the dataset's circular mean; circles Xin, squares Manella;
dashed: identity. Colours: organ systems.

**Tags:** registered (A111, A129), not blind (scored on timed data); Manella H111.1 fails by design (registration flaw
noticed at scoring).

**Caveats (speaker notes)**
- 3-4 tissues per dataset: the exact p cannot be small.
- A whole-body shift is not identifiable untimed: only relative shifts are compared.
- Manella lung reads about 1.6 h short (untimed -5.76 h, timed -7.37 h).
- Night shift (A257), forebrain knockout (A319), conflicting drives (A322, A324) and feeding removal (A297) go to
  section 4 (drives).

**Source:** `exploratory/A111_feeding_controls/output/{a111_h111_1,a111_h111_3}.csv`; `docs/results_preliminary.md` WP3
(A111, A129).

## 1.6 Tissues of one living person, sampled once without a recorded time, share one clock state; her own cell line does not

**Bullets**
1. Two tissues of a living person share one clock state: skin sites and layers r 0.59-0.81 (A251, 146-150 women per pair), fat and
   muscle or visceral fat at surgery 0.69 (A255 H255.3, 83 persons), fat and skin 0.53 (A255 H255.1, 680 women); each above
   200 expression-matched and 200 stress-responsive gene sets.
2. The woman's own cell line, which carries her genome but not her body time, shares none of it: r 0.04 with fat, 0.07 with
   skin (H255.2). The shared state is body time, not constitutive clock-gene expression.
3. On the core loop (11 genes, A348 blind) it holds in skin sites (0.57), skin layers (0.70, 0.51) and the surgical pairs
   (0.72); fat-skin from one punch biopsy falls to the line (0.31 against matched genes q95 0.33, p_cal 0.057).

**Figure:** `figures/living_tissues/living_tissues.png` (`scripts/living_tissues_compute.py`, `scripts/living_tissues.py`;
wide). A, B: MuTHER, pointwise mutual information [bits] of the clock-state angle in fat against skin (A) and against
the own cell line (B); von Mises kernel density on the 24 h torus (SD 1.5 h), grey where fewer than 0.5 persons are
expected per 1 h cell. Mean on the diagonal 0.87 bits (fat-skin) against 0.03 (fat-cell line). C: r per tissue pair, 14
genes (ink) and core loop (blue), with matched and stress-responsive q95 lines; open = not clock-specific.

**Tags:** A251, A255 registered, not blind for the earlier GTEx pairing; A348 blind for every core value; A and B
descriptive (kernel width and grey cut-off chosen for display; the compute script reproduces the scored r exactly).

**Placement:** after the positive controls (1.5a/b) as the negative control: same genome, no body time, no pairing.
Moved from section 3 (user, 2026-11-23): it validates the premise (one time per person) rather than adding biology,
and it closes the post-mortem objection before the GTEx results of section 3.

**Caveats (speaker notes)**
- These are living people: the GTEx pairing between organs is not a property of death.
- The angle is the direction of a person's deviation from the cohort mean, not her time of day. With biopsies at clinic
  hours the states lie on a short arc (almost a chord); a woman seen earlier than average points one way, one seen later
  the opposite way. Hence the two red blobs about 12 h apart on the diagonal (fat angles: 189 women at 9-12 h, 185 at
  21-24 h, 51-84 in the other 3 h bins). It does not mean persons are shifted by many hours.
- r = |sum za conj(zb)| / sqrt(sum |za|^2 sum |zb|^2) on the clock-plane projection (no model fit): 1 when each person's
  state in one tissue is her state in the other turned by one common angle.
- MuTHER fat and skin come from one punch biopsy and one array batch; the responsive companion bounds what shared
  handling gives. Women only, arrays.
- Companion families differ: narrow responsive sets on 14 genes (A251, A255), calibrated broad sets on the core loop
  (A348; the narrow family is not calibrated there).
- Co-twins (H255.4) not shown.

**Source:** `exploratory/A251_skin_sites/output/score_tests.csv`; `exploratory/A255_living_tissues/output/score_tests.csv`;
`exploratory/A285_companion_calibration/output/a348/{a251,a255}/tests.csv`;
`presentation/output/living_tissues/{muther_phases,check}.csv`; `docs/results_preliminary.md` A251, A252, A255, A256, A348.

---

# 2. Before birth: gates and coherence only

| slide | status | figure |
|---|---|---|
| 2.1 Where the fetal clock-shaped state is | **revised** from `results.md` R4 panel A (below) | `figures/development/` |
| 2.2 Is it a clock? Only the adrenal holds | **new** (below) | `figures/fetal_clock/` |
| 2.3 Liver and lung: not current glucocorticoid activity under another name | **new** (below) | `figures/fetal_signal/` |

## 2.1 Before birth, some organs share a clock-shaped state between individuals; others, at adult power, do not

**Bullets**
1. Human fetus: a clock-shaped covariation of the 14 clock genes between individuals in the adrenal cortex (C 0.61, 34
   fetuses), liver (0.58, 35) and lung (0.53, 38), above expression-matched and stress-responsive gene sets, at about
   two thirds of the adult value at the same design (0.85, 0.86, 0.84) (A243, A244). Whether it is a clock is tested
   on the next slide (only the adrenal holds).
2. The adrenal replicates blind in a new laboratory before 12 weeks: whole adrenal, del Valle 2023, C 0.59, p 0.002, 32
   fetuses at 7-11.5 wpc (A328 part a, H328.1 pass).
3. None in the heart (0.03, 72 fetuses), kidney (0.02, 53) or adrenal medulla (-0.05, 30), where an adult-size clock
   gives 0.80-0.90 and passes in 100 % of subsamples; none in the brain (-0.17, 20; A239).
4. Late mouse embryo (E15-E18.75, A309): hepatocytes replicate (0.64); ventricular cardiomyocytes (0.24) and nephron
   progenitors do not covary.

**Figure:** `figures/development/development.png` (`scripts/development.py`; half). Organs along x; bars: measured C,
human fetus blue, mouse embryo orange (right slot empty when not measured), filled when significant and above both
companion families, open when not; open diamond: median C of adult tissue subsampled to the same design (brain:
postnatal brain at 20, A240); grey line: matched genes, 95th percentile. The R4 pair panels (B-E, PER1 / PER2 / NR1D1
with the core genes, A307) are off this slide.

**Tags:** A243, A244, A239 registered, blind; A328 part a registered, blind; A309 registered, not blind for the earlier
fetal results. Aggregates only.

**Caveats (speaker notes)**
- No sampling times: C reads a shared state between individuals, not a rhythm.
- "Two thirds of adult": C is not linear in amplitude, a rough reading.
- Brain: the negative excludes only a full-size clock (adult pass rate 0.95 at full size, 0.365 at three quarters; A240).
  Hindbrain at 4-8 wpc: C 0.41, unadjusted p 0.040, Holm 0.278, not clock-specific (A237 H237.3: suggestive, not
  established).
- Mouse adrenocortical cells (0.41, 13 embryos) and lung progenitors (0.12) do not replicate in this window; mouse
  kidney -0.52 is one gene (NPAS2; 0.07 without it).
- Mouse heart clocks run in explants at E15-E18 (Umemura 2017) yet show no shared state in vivo: cells able to
  oscillate, not synchronised.

**Source:** `exploratory/A243_fetal_organs_pooled/output/{l1,verdicts,position_in_adult}.csv`;
`exploratory/A244_fetal_adrenal_classes/output/{l1,verdicts_a244,position}.csv`;
`exploratory/A328_fetal_running_loop/output/{score_a/l1,score_a/verdicts,power_a/passrates}.csv`;
`exploratory/A239_fetal_brain/output/{l1,verdicts}.csv`; `exploratory/A240_fetal_amplitude/output/curve_array.csv`;
`exploratory/A309_mouse_embryo_organs/output/{l1,power_passrates}.csv`; `docs/results_preliminary.md` A237-A240, A243,
A244, A309, A328.

## 2.2 Only the fetal adrenal carries a clock-like organisation

**Bullets**
1. In the blind adrenal cohort the core loop (11 genes, without PER1, PER2, NR1D1) covaries between fetuses on its own:
   C 0.80 against 0.59 on the 14 genes, above both companion families (A328 part a H328.2, pass). In liver (0.69) it is
   significant but not above stress-responsive genes (q95 0.71); in lung (0.47) and the single-cell adrenal cortex
   (0.42) it is not significant (A266, A269).
2. The fetal states are not one slow factor moving fetuses along a line (A340 H340.1, pass in the adrenal on 14 genes
   and the core loop, and in liver), but more elongated than a daily clock sampled in hospital hours.
3. No fast pulsatile drive: the pre-mRNA lead stays below what an input every 6 h or faster gives (3.09 h at a 1 h mRNA
   lifetime), in the adrenal core loop (upper bound 0.98 h; input slower than 24 h) and liver (1.08 h; 22 h) (post hoc).

**Figure:** `figures/fetal_clock/fetal_clock.png` (`scripts/fetal_clock.py`; wide). A. C, 14 genes (black) and core loop
alone (blue), filled when clock-specific; lines across each bar: matched (light grey) and stress-responsive (dark grey)
95th percentiles. B. elongation of the fetal states (0 circle, 1 line), 90 % interval, against a daily clock sampled in
hospital hours (blue band), no clock (grey band) and the lower bound of one slow factor (dark line). C. pre-mRNA lead
over mRNA [h], 90 % interval, with the running postnatal liver (grey); dashed: expected lead for a daily input; red
band: input period <= 6 h.

**Tags:** A328 part a, A340 registered, blind; A266, A269 registered, not blind for the earlier fetal results; panel C
post hoc (A338 threshold, A339 lifetime; the minimum-period companion of A339 not blind).

**Caveats (speaker notes)**
- Not lost power: an adult clock at half size gives core-loop C 0.69 in liver and 0.70 in lung and passes in 91 % and 84 %
  of subsamples (A269 power stage); a three-quarter-size loop is excluded in the single-cell cortex (0.995), a half-size
  one not (0.63) (A266). Liver is a specificity problem, lung a small core loop.
- The two adrenal datasets disagree on the core loop (whole adrenal 0.80, cortex cells 0.42); cause open (whole tissue vs
  cortex cells, one cohort of 32 vs several smaller).
- Liver states lie inside the no-clock elongation range (0.42-0.74) with their main axis 8 deg from the
  glucocorticoid-target direction (adrenal 32 and 88 deg; descriptive).
- The matched estimator reads leads about 2-fold low (A337): minimum periods may be overestimated up to 2-fold (adrenal
  core and liver still about 11 and 10 h; adrenal 14 genes not excluded). Conditional on fetal clock-mRNA lifetimes near
  the adult 1-2 h. The damping test (A339) was uninformative.
- The lead does not separate a running loop from a daily drive: the running postnatal liver is unresolved too (A328
  part b, A338).

**Source:** `exploratory/A328_fetal_running_loop/output/{score_a/l1,score_b/verdicts}.csv`;
`exploratory/A244_fetal_adrenal_classes/output/{l1,verdicts_a244}.csv`; `exploratory/A266_fetal_adrenal_drive/output/{l1,feasibility_passrates}.csv`;
`exploratory/A243_fetal_organs_pooled/output/l1.csv`; `exploratory/A269_fetal_organs_drive/output/{l1,feasibility_passrates}.csv`;
`exploratory/A340_fetal_geometry/output/score/verdicts.csv`; `presentation/output/fetal_clock/lead_bounds.csv`;
`docs/results_preliminary.md` A266, A269, A328, A338, A339, A340.

## 2.3 The fetal liver and lung covariation is not current glucocorticoid activity under another name

**Bullets**
1. Regressing out current glucocorticoid-target activity (four direct receptor targets, FKBP5, TSC22D3, ZBTB16, KLF9;
   linear, own slope per clock gene) leaves the covariation intact: liver 0.58 to 0.59, lung 0.53 to 0.57, still
   clock-specific (A269 H269.2, pass). A single shared glucocorticoid axis is excluded; a rhythmic glucocorticoid input
   acting with delays is not.
2. Glucocorticoid activity does not rank the fetal organs: heart and adrenal medulla vary most between fetuses and show
   no clock-shaped covariation; lung varies little and covaries (Spearman 0.25, one-sided p 0.30; A341, descriptive).

**Figure:** `figures/fetal_signal/fetal_signal.png` (`scripts/fetal_signal.py`; half). A. 14-gene C before (filled) and
after the glucocorticoid score is regressed out (hatched); grey line: matched genes, 95th percentile. B. seven fetal
organs, 14-gene C against the between-fetus variance of the glucocorticoid-target score; filled: clock-specific; brain
open (three genes, not noise-corrected). Colours: organ systems.

**Tags:** A243 registered, blind; A269 registered, not blind for the earlier fetal results; A341 descriptive by its power
rule.

**Caveats (speaker notes)**
- Four genes are a noisy measure of glucocorticoid activity: the adjustment may remove less than the true effect.
- PER1, PER2 and NR1D1 keep the clock's phase pattern with the core genes in lung, and follow glucocorticoid activity
  there only (A307).
- A341 has seven organs: only a very strong ranking (rho about 0.90) is detectable.
- Other systemic inputs (cAMP, feeding, temperature) are not addressed.

**Source:** `exploratory/A243_fetal_organs_pooled/output/l1.csv`; `exploratory/A269_fetal_organs_drive/output/l1.csv`;
`exploratory/A341_gr_ranking/output/{measures,spearman}.csv`; `docs/results_preliminary.md` A269, A307, A341.

## Section 2 conclusion

Before birth, a clock-like organisation is established in one organ, the adrenal: its core loop covaries between
fetuses on its own, replicated blind before 12 weeks, and its states turn in two dimensions rather than along a line.
In liver and lung the 14 clock genes covary in the clock's pattern, but the signal sits in, or cannot be separated
from, signal-responsive genes, consistent with a systemic input; glucocorticoids are a candidate, not shown. Heart,
kidney, adrenal medulla and brain show none at adult power. Whether the adrenal clock runs on its own or follows a daily
maternal or ACTH signal remains open. No sampling times: these are shared states, not rhythms.

---

# 3. Adult health: one time per person, many volumes

| slide | status | figure |
|---|---|---|
| 3.1 The reference organisation: GTEx fast deaths, whole day | **new** (below) | `figures/reference_map/` |
| 3.1b Whole blood is the exception; testis and brain, as quiet, keep time | **new** (below) | `figures/blood_exception/` |
| 3.2 Inside tissues: cell classes run at the organ's time with their own loudness; lymphocytes do not share it | **revised** from `results.md` R2.3 (below) | `figures/celltypes_health/` |
| 3.3 Organs differ in loudness: not their cells, not post-mortem attenuation | **new** (below; composition figure to backup) | `figures/organ_gain/` |
| 3.4a Louder organs are more coupled (human, baboon) | **revised** from `results.md` R2.1 (below) | `figures/organs_health/` |
| 3.4b Because each organ's departure has a fixed size, not a fixed angle | **new** (below) | `figures/organs_departure/` |
| 3.5 Persons differ in a body-wide expression scale; in deaths the clock follows it most reliably | **revised 2026-11-24** (below) | `figures/person_loudness/` |

CHIRAL concordance (A123) is slide 1.1. The ventilated and brain-dead rows of `results.md` R3.1 go to sections 4 and 7
(`figures/holds_organisation/`, unchanged).

## 3.1 The reference organisation: each organ keeps its place on the body's time, and organs differ mostly in loudness

**Bullets**
1. Most peripheral organs peak up to about 1.5 h before the person's body time, with gains of about 1.2 to 2.3; a
   few peak clearly earlier: adrenal −4.86 h (95 % −5.57 to −4.23), lung −2.88 h, unexposed skin −2.19 h.
2. Brain regions read on the body's host peak later (+0.97 to +2.58 h) and are quieter (gain 0.40 to 0.57).
3. Gains span 0.26 (whole blood) to 2.28 (gastro-oesophageal junction); testis 0.38. Loudness, not timing, is where
   organs differ most (sections 3.4 onwards).

**Figure:** `figures/reference_map/reference_map.png` (`scripts/reference_map.py`; half). Clock face of slide 1.1: angle =
peak time of the tissue's clock genes relative to the person's body time (0 h; later peaks clockwise, as CHIRAL
Fig. 1B: the clocksync offset, > 0 advanced, negated; host from the peripheral tissues, peripheral tissues read
leave-one-out), radius
= relative amplitude (gain); crosses: 95 % donor-bootstrap intervals (200 draws). 40 tissues, 41 to 220 fast-death donors
each; organ-system colours; direct labels only on tissues discussed later.

**Tags:** A263 R2 fit (registered readout), per-tissue intervals from `organs_compute.py`; display descriptive (labelled
tissues chosen for the deck).

**Caveats (speaker notes)**
- 0 h here is the person's body time; on slide 1.1 it is the mean over tissues.
- Sign convention (user, 2026-11-23): organ maps show peak times (+ = clock genes peak later), as CHIRAL; numbers on
  these slides are peak times. `couplings.csv` and the clocksync API keep the method's convention (offset > 0 =
  advanced = peaks earlier): the adrenal is +4.86 h there.
- Whole blood has no readable core clock (post hoc entry 2026-11-03): its phase interval spans 1.09 to 3.72 h.
- The brain is read on a host built from peripheral tissues: its offset is relative to the periphery, and its low gain
  partly reflects that reference (R2.1, A289).
- Fast deaths are the reference because their body clock reads the hour of death (R 0.63, slide 4.1).

**Source:** `presentation/output/organs/couplings.csv` (A263 R2 fit with donor bootstrap, `scripts/organs_compute.py`).

## 3.1b Whole blood is the exception: its core clock genes keep no better time than random genes, while other quiet tissues keep it

**Bullets**
1. In GTEx fast deaths, the 14 core clock genes of whole blood (212 donors) follow the hour of death no better than a
   typical gene: median cosinor R^2 0.0096, 52nd percentile of all expressed genes; in five solid tissues (179-220
   donors) they sit at the 99.5th to 99.96th percentile.
2. Their phases agree with the body waveform at 0.43 (solid tissues 0.85 to 0.94). Removing cell composition changes
   neither (0.0088; 0.48).
3. Not a matter of loudness: the other quiet tissues keep time. Testis (gain 0.38; 135 donors): 92nd percentile,
   phase agreement 0.92; the 13 brain regions (gains 0.40 to 0.57): median 99.85th percentile (90.8 to 99.96), agreement
   0.84 (0.79 to 0.90).
4. Consistent with lymphocytes not sharing the organ's time (slide 3.2): an expected exception, not a failure of the
   readout.

**Figure:** `figures/blood_exception/blood_exception.png` (`scripts/blood_exception.py`; wide, bar plots). A. median
cosinor R^2 of the core clock genes against the hour of death, as a percentile of all expressed genes of the tissue
(dashed: a typical gene). B. phase agreement of the clock genes with the GTEx waveform, as measured (filled bars) and
after regressing out composition scores (open bars; blood-cell markers: neutrophil, T, B, NK, monocyte, eosinophil, erythroid markers). Brain: bar at the median over the 13 regions, one dot per region. Organ-system colours.

**Tags:** post hoc, descriptive (user observation stated before running, 2026-11-03; testis and brain added post hoc
2026-11-25 as low-gain comparisons).

**Caveats (speaker notes)**
- Shown in post-mortem bulk blood only.
- Timing signatures built from living blood are not a readout of a blood core clock (user's reading, 2026-11-25; gene
  lists checked, not tested). Their selectors, trained on timed samples, never pick the loudest core-clock outputs
  (ARNTL / BMAL1, DBP). They pick glucocorticoid targets (TimeSignature, whole blood: DDIT4 and PER1 selected in every
  run; BodyTime, monocytes: FKBP4, KLF9, PER1), the fast signal layer every cell reads (A275, A279, A343), and
  cell-type markers (TimeSignature: IL1B, CD38, GZMB, CD1C, CLEC10A). NR1D1 / NR1D2 (both) and PER2, PER3, CRY1
  (BodyTime) also appear. This is the same picture as GTEx blood. Braun 2018 PNAS Table 1 (41 genes); Wittenbrink 2018
  JCI (2-gene predictor NR1D2 + PER2; 12- and 13-gene predictors).
- Blood carries little time-of-death signal at all: the top 100 time genes (chosen in one donor half) reach held-out
  R^2 0.021 vs 0.014 for random genes, and composition removes 6 % of it. That its time signal is diurnal composition
  is therefore only weakly supported here: there is little signal to explain.
- Living forced desynchrony (A262, backup): the blood program damps to about a third and is displaced about 4 h; the
  night-shift PBMC design does not pass the gate.
- Possible test (not registered): whether blood's weak time-of-death signal sits in glucocorticoid targets rather than
  in the core loop.

**Source:** `exploratory/blood_composition_posthoc/output/{blood_composition,blood_composition_low_gain}.csv`; `docs/results_preliminary.md`, post hoc
entry "GTEx Whole Blood has no core-clock time signal" (2026-11-03); `presentation/output/organs/couplings.csv` (gain).

## 3.2 Inside tissues: every resident class runs at the organ's time with its own loudness; lymphocytes do not share it

**Bullets**
1. Each cell class has its own gain relative to endothelium of the same atlas, the same across organs: fibroblast 1.26
   (95 % 1.07 to 1.48), mural 0.92, parenchyma 0.88, epithelial 0.80, myeloid 0.67 (0.53 to 0.84) (A293 H293.4, pass).
   On the core loop (A350, blind) conserved in four of five classes; mural cells not (tau 0.10 -> 0.39, H350.1 fail).
2. Classes run at the organ's time: stromal-vascular classes peak within about 1 h of each other; endothelium peaks
   0.46 h before fibroblasts (90 % 0.09 to 0.90), the mouse direction, but the order is not shown (H293.3 fail).
3. Lymphocytes do not share the tissue's clock state (pairwise, beyond matched and stress-responsive genes, A272-A288;
   coupling median 0.13 in the atlas fits). Whether they keep a clock of their own is open: the own-clock test would
   only see one about as loud as the fibroblasts'.

**Figure:** `figures/celltypes_health/celltypes_health.png` (`scripts/celltypes_health.py`; wide, two panels). A. gain
relative to endothelium per class: atlases (dots), pooled 14 genes (black, 95 %), pooled core loop (open blue; mural
"not conserved"); lymphocytes at the label "not coupled to the tissue's clock", no position on the gain scale. B. lead of the first class between
stromal-vascular classes [h]: human atlases, circular mean with 90 % atlas bootstrap, mouse against known time
(diamonds). Lymphocyte panels removed (user, 2026-11-23): too much weight on a negative.

**Tags:** A293 parts 2-3 registered, NOT BLIND (per-family gains and mouse offsets seen earlier); A350 blind on the core
loop; A281 registered; A321 blind for the atlases; A326 registered, not blind for the human results.

**Caveats (speaker notes)**
- Why this slide before composition: classes differ in loudness and lymphocytes have none, so composition could move an
  organ's gain; classes sit within about 1 h of each other, so it can barely move its phase (3.3 tests both).
- "Lead of the first class" is in words: positive = the first class peaks earlier (unaffected by the sign convention).
- Heart lymphocytes (0.50) sit just below matched genes (0.58); breast CD4 T (0.38) above (0.30) is a hormonal pulse:
  without PER1, PER2, NR1D1 not significant (A294).
- Epithelial and myeloid tau are close to the bound (0.30); "myeloid" mixes macrophages (which carry a clock, A288)
  with monocytes and dendritic cells.
- Sharing per family, the calibrated stromal-vascular minus lymphocyte contrast (A272-A288), living donors (HLCA) and
  residency (A284): as in `results.md` R2.3 speaker notes.
- Why myeloid gain is read at 0.53 of the fibroblasts' while a lymphocyte clock below 0.7 is not: different tests. The
  gain is host-referenced (a class's states regressed on the tissue's shared state, estimated from all compartments);
  the own-clock test (A281, A321) asks the class's clock genes to cohere among themselves with no reference, and A326
  calibrated it on fibroblasts. Classes relative to fibroblasts: endothelial 0.79, mural 0.73, parenchyma 0.70,
  epithelial 0.63, myeloid 0.53; a lymphocyte clock of myeloid size is not excluded by the own-clock test.
- In the host-referenced A17 atlas fits (descriptive, post hoc; `results/figure_2/output/atlas_fits/a17_healthy_fits.csv`,
  medians over compartments): coupling fibroblast 0.65, endothelial 0.53, myeloid 0.27, lymphocyte 0.13 (66
  compartments, 24 % fitted at the coupling boundary); width of the log-amplitude interval 0.92, 0.87, 1.58, 2.09. A
  lymphocyte gain on panel A's scale would need its own registration (A293 part 3 method, power note first).
- Lymphocytes (detail for questions, not on the slide):
  - No readable clock: their coherence between donors is above matched genes in 1 of 9 atlases (breast,
    a hormonal pulse, A294), fibroblasts of the same donors thinned to their depth in 7 of 9; blind on the core loop in
    fresh atlases, none (A321, 3 informative compartments, combined p 0.47).
  - What that excludes (A326): a lymphocyte clock at 0.7 of the fibroblasts' would have been read (design power 0.996
    and 1.00); a mouse-colon-sized one (0.25 of the fibroblasts', q95 0.40) would not (0.00 and 0.05). A bound, not a
    species or residency difference.
- Skin T cells (A321) cohere but their thinned fibroblasts never pass: uninformative, an exception to check.

**Source:** `exploratory/A293_organisation_levels/output/part3/families.csv`, `part2/tests.csv`;
`exploratory/A350_class_gain_core/output/compare.csv`; `exploratory/A281_lymphocyte_clock/output/pairs.csv`;
`exploratory/A294_lymphocyte_hormonal/output/pairs.csv`; `exploratory/L1_lymphocyte_power/output/design_power.csv`,
`scout/output/ratios_vs_fibroblast.csv`; `docs/results_preliminary.md` A293, A321, A326, A350.

## 3.3 Organs differ in loudness, and the difference is neither their cells nor post-mortem attenuation

**Bullets**
1. Organs made mostly of one cell class differ in gain: squamous epithelium 2.6-fold (oesophagus mucosa over sun-exposed
   skin), adipocytes 1.6-fold (visceral over subcutaneous and breast), smooth muscle 1.2-fold (oesophagus muscularis over
   sigmoid colon), while a class's relative gain varies between organs within about 1.4 (3.2).
2. Composition and the conserved class gains predict a 1.29-fold range over eight tissues; the organs span 2.89-fold
   (A323, registered; the rank test fails, rho 0.52, p 0.098, 8 tissues).
3. Not post-mortem attenuation: within each tissue, donors with 6.7 h longer ischaemia have the same gain (1.01-fold,
   90 % 0.98 to 1.04), donors with higher RIN too (0.99-fold, 0.96 to 1.03); ischaemia and RIN could produce at most
   about 4 % of the organ spread (A356, registered, blind, H356.1-2 FAIL as the artefact needs them to pass).
4. On top of its cells, each organ scales how loudly it expresses the body's time. Reproduction of the human organ
   ratios in living tissue is not shown (A355).

**Figure:** `figures/organ_gain/organ_gain.png` (`scripts/organ_gain.py`; wide). A. tissues made mostly of one class:
gain per tissue, highest over lowest. B. eight tissues: gain predicted from composition and class gains against
observed, each relative to its geometric mean (A323). C. gain ratio within tissue, long over short ischaemia and high
over low RIN: tissues (dots, organ-system colours) and pooled with 90 % (A356). The former composition figure
(`figures/organs_composition/`) moves to backup.

**Tags:** A323 registered, not blind for the GTEx and class gains; A356 registered, blind for every split gain; panel A
post hoc display (same-class rows chosen post hoc).

**Caveats (speaker notes)**
- Within tissues, composition does not move time (A295 H295.2: +0.05 h per SD stromal-vascular, -0.05 h lymphocytes,
  against organ peak times from -4.86 h, adrenal, to +1.20 h, testis) and moves gain only a few percent (1.04-fold per
  SD stromal-vascular, descriptive; lymphocytes 0.99, H295.1 fail); adrenal cortex +1.15 h towards the body's time,
  gain 1.30 (A299, descriptive). Backup: `figures/organs_composition/`.
- A323 uses eight tissues of the GTEx single-nucleus atlas (16 donors), class gains from other atlases, nuclei not
  cells; the "organ factor" is an inference by elimination.
- Living pairs (A355): omental over subcutaneous fat 1.38 against GTEx 1.59 (not decidable, 26-29 men); cardia over
  squamous oesophagus opposite to GTEx (mucosal biopsies against full-thickness stomach); skin near 1; the pooled test
  fails. In the timed baboon, organ amplitudes span 2.8-fold (10th-90th percentile, descriptive).
- Not tested: agonal state and autolysis (fast deaths only).
- Loudness is per-cell amplitude times synchrony among the cells; bulk and pseudobulk cannot separate the two, and
  public 10x data cannot either (A334: the split fails a controlled dephasing; the dispersion is not clock-specific in
  the aorta). Either way it is a property of the organ.
- Panel C excludes whole blood (ischaemic time undefined for a draw). Across tissues, organ gain vs median ischaemia
  rho -0.37 (Holm 0.076): organs differ by under 1 h in median ischaemia, so this cannot be attenuation.

**Source:** `exploratory/A323_gain_decomposition/output/{tissues,tests}.csv`;
`exploratory/A356_gain_attenuation/output/{score_tissues,score_pooled,score_across,score_descriptive}.csv`;
`exploratory/A355_living_gain_ratios/output/{score_pairs,score_pooled}.csv`; `presentation/output/organs/couplings.csv`;
`presentation/output/organs_composition/rows.csv`; `docs/results_preliminary.md` A295, A299, A323, A355, A356.

## 3.4a Louder organs are more coupled to the body's clock state, in human and baboon

**Bullets**
1. Across 27 GTEx peripheral tissues, coupling rises with the amplitude a tissue shares with the body (Spearman rho
   +0.734). Organs differ 9.2-fold in that amplitude (S spans 85-fold: 0.12, whole blood, to 10.3). Not measurement: a
   perfectly locked tissue would read 0.97 at every tissue's real design (A291).
2. The same holds in the baboon, where each animal's time is known, brain included: coupling rises with each tissue's
   rhythm amplitude (Spearman 0.82 over 57 tissues; periphery 0.77, brain 0.74), and weakly rhythmic brain regions fall
   near the peripheral curve (A289: median residual -0.001 on its line). Human brain couples at 0.406 against 0.679 in the periphery.

**Figure:** `figures/organs_health/organs_health.png` (`scripts/organs_health.py`; half). A. human: coupling against the
shared amplitude a = sqrt(S), one dot per tissue; grey: a^2 / (a^2 + D) = S / (S + D) at the median D (1.89), the definition of coupling at a
constant departure, not a fit (explained in 3.4b). B. baboon: amplitude against known time (A36 timed gain, estimated
without the host) against coupling, 35 peripheral tissues and 22 brain regions; grey: the same curve, a^2 / (a^2 + D),
D fitted on the peripheral tissues (2.30; one parameter; amplitude units differ from A). Right axes: the phase SD a coupling implies,
(24 / 2 pi) sqrt((1 - c) / (2 c)) hours (isotropic departure, small angle, as 3.4b panel C; approximate beyond about
3 h): coupling 0.8 is about 1.3 h, 0.5 about 2.7 h.

**Tags:** A104 registered (H104.1, rho < 0, not supported; the positive relation is post hoc, made registered by A144
and A271 Part B). A291 registered (H291.1 pass). A289 H289.3 registered (fail: brain on the line, the outcome shown);
panel B display post hoc (the fitted curve replaces A289's peripheral line; on it the brain's median residual is
+0.08, periphery +0.04, descriptive; A289's registered residual, -0.001, is on its linear line).

## 3.4b Why: each organ's departure from the body's clock state has a fixed size, not a fixed angle, so a louder organ stays closer to the body's time

**Bullets**
1. A tissue's own departure from the body's clock state (D) barely grows with its shared variation S (slope of log D on
   log S 0.07; D spans 4.6-fold, 0.90 to 4.14, while S spans 85-fold). A phase jitter of fixed angle would give slope 1,
   individual clock periods slope 2.
2. So the same displacement turns a loud clock by a small angle and a quiet one by a large angle: the phase SD it
   implies falls with amplitude, from 3.0 h (sun-exposed skin) to 1.1 h (oesophagus muscularis), and matches the
   departure measured in hours by a different estimator (A259; Spearman 0.945 over 11 tissues, 1.0 to 3.2 h).
3. The departure is as much amplitude as phase (fraction along phase 0.49; matched genes 0.38; A311): not a phase
   jitter. It is clock-specific (2.6 times matched genes, A294b) and a coherent shift of the whole clock state (A329
   erratum).
4. Coupling, the share of an organ's person-to-person clock variation that follows the body's clock state, rises with
   the organ's gain, the transfer of the body's time to the organ. The organ's own departure is added after that
   transfer and is about the same size in every organ, so louder organs receive more of the body's time and
   stay closer to it in hours (1.1 h against 3.0 h): gain and coupling rank organs alike. What sets an organ's gain is
   open: not its cell composition (A323), not post-mortem attenuation (A356).

**Figure:** `figures/organs_departure/organs_departure.png` (`scripts/organs_departure.py`, from
`scripts/organs_departure_compute.py`; wide). A. schematic to scale with real values: sun-exposed skin (quiet) and
thyroid (loud) on circles of radius sqrt(S), departure disc of 1 SD (sqrt(D / 2)), wedge the implied phase SD.
B. D against S per tissue; grey lines through the geometric-mean S and median D: fixed size (0), fixed angle (1), own
clock periods (2). C. phase SD in hours against amplitude: model estimate sqrt(D / (2 S)) (dots) and the gene-split estimate (A259)
with 90 % interval (diamonds); grey: fixed-size prediction and fixed-angle line. Whole blood and testis omitted. D. fraction of
the departure along the phase direction, clock genes and matched genes (A311). Reference lines are named in legends.

**Tags:** A294b registered (H294b.1 pass; slope descriptive); A259 registered; A311 registered (H311.1 FAIL: not
phase-dominated, the outcome shown); A329 erratum post hoc; panel A and the implied phase SD (C) post hoc displays.

**Caveats (speaker notes)**
- S = z2 V_h (gain squared times host spread; invariant to the host's scale); D = 2 psi (scatter around the state the
  body predicts, in amplitude and phase, net of counting noise). Both in clock-plane log-expression units squared.
  Definitions as in `results.md` R2.1 speaker notes.
- 3.4a: the curve in panel A is not a test: it is the definition of c. The finding is 3.4b panel B (D nearly
  constant while S spans 85-fold).
- 3.4b panel C: small-angle approximation (sqrt(D / 2) is the tangential SD, sqrt S the RMS radius of the shared state).
  Dots, "model estimate": from the host factor model (A294b), Poisson noise modelled, half of D counted as phase.
  Diamonds, "gene-split estimate" (A259): the 14 genes and the host tissues split in halves, noise cancels in the
  cross-moments, the whole deviation counted as phase (registered as an upper bound); 11 body tissues, each against the
  other ten. Same GTEx donors, so the agreement (Spearman 0.945) is between estimators, not datasets. In A259 the
  deviation was not above expression-matched genes in any tissue (matched genes alone reach 30-74 % of the clock's
  value): the hour values include sample-level variation, not only clock. Whole blood (10 h implied) and
  testis (5.2 h) are omitted: no shown clock, and the small-angle approximation fails.
- 3.4b panel D: matched genes lean towards amplitude (a sample-level scale component, A311); relative to them the clock
  departure is more tangential, but not significantly (90th percentile of matched sets).
- Specificity (former panel): own departure of the clock genes 2.6 times that of expression-matched genes (90 % 2.4 to
  2.9, all 27 tissues; A294b, registered).
- Why gain is a transfer (post hoc reading): a departure that arose before the gain (e.g. the body's signal reaching
  each organ differently in each person) would be amplified with it (D proportional to S, slope 1, equal coupling in
  every organ). Slope 0.07 places the departure after the transfer, whatever the clocks' dynamical regime. "Added after
  the transfer" covers both a real fluctuation of the organ's own and variation in the sample's state at death; GTEx
  cannot separate them. Do not read A311 as "tissues are driven below their self-sustained regime": A311's own reading
  is that the departure is variation added to the state, not a dynamical perturbation of the cycle (a share near 0 would
  contradict the self-sustained liver, A297; matched genes read radial, so the absolute share is biased).
- Slope reading (post hoc): 0.22 for the excess over matched genes (A294b H294b.2). A kicked limit cycle would also
  give a departure of fixed size, but mostly along the phase (A311 fails that), so the departure is variation added to
  the state: person-specific cue state (the A357 glucocorticoid part is small), sample-level variation or the state at
  death. Feeding cue history not supported (A302, H302.1 fail).
- Why not the core loop (A329 erratum): the projector maps a coherent clock-state shift the same way with 14 or 11
  genes, but the 11-gene projector amplifies gene-level variation 3.1-fold. So the core-loop ratio (1.28) and coupling
  fall without the coherent part changing. Split per tissue: coherent departure median 1.42 (0.62 to 3.69, positive in
  27 of 27), clock genes' gene-level part 0.47 against 0.71 for matched genes. Agrees with the A336 gene-drop null
  (the trio's drop at the 32nd percentile of similar three-gene drops).
- The glucocorticoid-like person state read unequally by organs (A357) is a drive: shown in section 4, not here.
- The departure in hours (A259): SD across donors 1.0 h (thyroid) to 3.2 h (lung); an upper bound, not above the
  companion sets on that statistic.
- The stress-responsive calibration of "stronger clocks are more coupled" (A290) is registered but not scored (only
  placebo draws on disk).
- Brain: lower coupling is not clock-specific (A263; every gene set drops); not the blood-brain barrier (A289: pineal
  and pituitary couple less than brain regions in baboon and sheep). Wording: "the brain is less coupled to the body's
  shared state", never "the brain's clock is the exception".

**Source:** `exploratory/A294_floor/output/{tissues,tests}.csv`; `exploratory/A329_core_departure/output/coherent_split.csv`;
`presentation/output/organs_departure/{tissues,summary}.csv` (from A294b, A311 `sets.csv`, A259 `score/targets.csv`);
`exploratory/A291_attenuation/output/tissues.csv`; `results/figure_2/output/organs/r2_couplings.csv`;
`results/figure_1/output/species/a36_tissues.csv`; `presentation/results.md` R2.1; `docs/results_preliminary.md` A289,
A259, A291, A294b, A302, A311, A329 erratum.

## 3.5 Persons differ in a body-wide expression scale; in deaths the clock follows it most reliably

(Revised 2026-11-24 after A361-A371; wording from `docs/PERSON_VOLUME_2026-11-24.md` section 12.)

**Bullets**
1. Persons differ in a body-wide scale of expression, and in GTEx deaths the clock follows it across a person's tissues
   more reliably than random genes: core-loop leverage 0.714, above 98.5 % of expression-matched gene sets (q95 0.621)
   and 96 % of stress-responsive ones (q95 0.692) (A336 part 2, registered). Its size is what random genes share too
   (cv2 0.120; matched median 0.071, responsive 0.262; post hoc size note).
2. It does not need a slow death: in violent, instantaneous deaths it is as large (core 0.726, 34 donors; 99th
   percentile of matched sets; A367, registered, gate relaxed by the user).
3. It is a plain scale: every core clock gene follows it alike (A363) and every organ nearly in proportion (A312); its
   organ pattern is not the one conflicting cues produce in mouse (A362).
4. In deaths with a terminal phase it falls with the glucocorticoid-target state (core loop -0.18, 90 % -0.30 to
   -0.06, beyond all 400 companion sets; A361); in instantaneous deaths it does not (+0.22; A367, descriptive): an
   agonal relation, most likely.

**Figure:** `figures/person_loudness/person_loudness.png` (`scripts/person_loudness.py`; wide). A. how reliably the
loudness is shared across a person's tissues (core loop), all fast deaths and violent deaths: clock genes (blue)
against 200 matched (light grey) and 200 stress-responsive (dark grey) gene sets. B. how large the shared loudness is
(cv2), same sets, all fast deaths. E (third panel). the clock's pairwise leverage under controls for batch, collection
centre, cell-type markers and the main expression axes, against the matched sets' 95th percentile under the same
control (A375). D (fourth panel). standardised effects on the person's scale (log host radius, 14 genes, A360) with 90 % intervals; black:
interval excludes 0.

**Tags:** A336 part 2 registered (core loop blind); A367 registered (blind; H367.1 at power 0.775 by the user's
decision); A361, A362, A363 registered; A366, A368 registered (living); A312, A360 registered; cv2 comparison post hoc.

**Caveats (speaker notes)**
- Leverage: correlation of the clock-plane radii of two halves of a person's tissues after the phase harmonics are
  removed (`clocksync.elasticity.radius_covariance`); cv2: their shared squared coefficient of variation.
- Deaths only. Processes acting within minutes of death or after it are not excluded. Not explained by RIN, ischaemic
  time, time of death, age, sex or BMI (A360), nor by batch, collection centre or body-wide composition (blood, immune,
  endothelial, fibroblast and smooth-muscle markers; top genome PCs): with all of these controlled the clock keeps 81 %
  of its pairwise leverage, at the 98th percentile of matched sets (A375).
- The GC contrast between violent and other fast deaths rests on 34 donors (power 0.24) and Hardy 1 differs in age,
  sex and cause of death.
- Living (removed from the bullets 2026-11-24): a clock-specific excess is not shown in living people. At matched
  phase, fat and skin of living women (MuTHER) share no more than random genes and a GTEx-size clock volume is excluded
  (A368); the earlier living leverage was the shared biopsy time (A366, A331 withdrawn).
- Living: STAGE (66 bypass patients) and AAGMEx (251 adults) are not decidable (A369, A370: arc width, noisy muscle
  clock). STARNET (about 600 patients, up to seven tissues) would decide it.
- Timed living skin (post hoc, descriptive): the published person-level amplitude shared by epidermis and dermis
  (Saint-Antoine & Anafi 2025, preprint; del Olmo data) reproduces (r 0.78, R^2 0.61) but random non-rhythmic genes
  reach r 0.70 (q95 0.72): mostly a shared scale of one biopsy. Between two separate skin sites (GSE326530) no
  amplitude is shared (their statistic r 0.13; our core clock r -0.04, 41st percentile of matched sets). Cite the
  preprint with this caveat.
- Elasticity (panel removed 2026-11-24, backup slide; every organ follows the scale nearly in proportion, A312): level not anchored at 1 (simulated proportional scaling reads median 1.11; grey band); the
  registered gradient has the size the same simulation gives with no true gradient (A312 calibration note).
- What the scale physically is (uniform damping of cell clocks, cells dephasing, body-wide composition) is not
  separable in bulk.

**Source:** `exploratory/A336_two_readouts/output/{person_gtex,person_gtex_companions}.csv`;
`exploratory/A367_hardy_agonal/output/{tests,companions,gc}.csv`; `exploratory/A312_elasticity/output/{tissues,calib_level}.csv`;
`exploratory/A358_A360_gain_sources/output/a360_coefficients.csv`; `exploratory/A368_muther_radial/output/tests.csv`;
`exploratory/A364_living_volume/output/{posthoc_clock,posthoc_global}.csv`;
`docs/results_preliminary.md` A312, A336, A360-A363, A366-A370, the size note and the post hoc skin entry.

---

# 4. Drives: where coupling to the body's time comes from (to build; revised 2026-11-24)

Section 3 describes the organisation and explains coupling without drives (transfer over a fixed departure, 3.4a-b).
Section 4 asks what sets the transfer and the organ's place in time. Source: `docs/COUPLING_SOURCES_2026-11-24.md`.

| slide | content | entries |
|---|---|---|
| 4.1 When the body no longer reads the hour, the organ arrangement holds; only brain death moves it (former 3.last) | **new** (below) | A27, A293, A314, A317, A374 |
| 4.1b Intensive care quiets the clocks; brain death keeps their amplitude | **new** (below) | A315, A335, A349 |
| 4.2 Removing a drive keeps amplitude, conflicting drives cost it | **new** (below): mouse liver without a feeding rhythm; daytime feeding across organs, weights and amplitude kept; one-line suggestive mention: organs that weigh one drive in mouse (liver, heart) are the loudest in human, the balanced lung the quietest, adipose the exception (A358, descriptive, rho -0.6, 5 organs, below its criterion) | A297, A322, A324, A358 |
| 4.3 A behavioural shift moves the human clock | night shift, muscle core loop +4 h, amplitude cost | A257, A335 |
| 4.4 A systemic signal read unequally by organs | glucocorticoid-like person state explains a small part of each organ's departure | A357 |
| 4.5 What sets an organ's gain | drive coherence (descriptive), reception, person volume; jet lag (in-house, to come) | registered 2026-11-24 |

## 4.1 When the body's clock no longer reads the hour, the organ arrangement holds; only brain death moves it

(User, 2026-11-25: slide 3.last moved here, as the opening of the drives section: dying and intensive care are not
physiological; brain death shown here rather than at the end of the life-course order. Amplitude: slide 4.1b.)

**Bullets**
1. The body's clock reads the hour of death in fast deaths (R 0.63, 244 donors), barely in slow (0.16, 106),
   intermediate (0.25, 51) or ventilated deaths (0.16, 79; brain death 0.17, 301) (A27 post hoc readout; A314,
   descriptive).
2. Slow deaths keep the fast-death organ arrangement: Fisher-Lee rho 0.78 (103 donors, 12 body tissues), above the 95th
   percentile of matched (0.13) and stress-responsive gene sets (0.44) (A293 H293.1, pass); organs moved by at most
   0.9 h (core 0.6 h), against an arrangement of about 1 h (median organ-pair offset 1.02 h) (A374).
3. Ventilated donors with a living brain keep it too: rho 0.690 (79 donors), above the matched (0.108) and responsive
   (0.318) 95th percentiles (A293 H293.2, pass), within sampling noise (expected 0.785, 5-95 % 0.519-0.879); organs
   moved by at most 1.8 h (core 1.5 h) (A374).
4. Brain death moves it: rho 0.161 at 79 donors, below the responsive 95th percentile; organs moved by at least 1.6 h
   (core 1.2 h) (A317, A374; descriptive). The arrangement is internal and held by the living brain, not by the day.

**Figure:** `figures/intensive_care/intensive_care.png` (`scripts/intensive_care.py`; wide; the panels of the former
slide 3.last with two rows added, user 2026-11-25). Rows: fast deaths (reference), slow, intermediate, ventilated with a
living brain, brain death. A. organ arrangement shared with fast deaths (clock genes; matched and responsive 95th
percentiles; diamond: expected from sampling noise at the group's donor count, brain death at 79). B. body clock
aligned to the hour of death, R. C. 90 % compatibility interval of the organ scatter (arrow: upper bound beyond the
simulated grid); dashed: the size of the arrangement. 14 genes; core loop in the notes.

**Tags:** A293 part 1 registered, NOT BLIND for the clock rho (known from A109 / A132; companions new): H293.1 slow and
H293.2 ventilated pass; intermediate deaths and brain death descriptive. Panel B: A27 post hoc readout, A314
descriptive. A374 registered (intervals blind, correlations not). A317 registered (brain death, H317.1 pass).

**A374 note (speaker notes):** a detection limit ("a scatter of SD 2.25 h would be seen", A349) exceeds the
arrangement itself; the compatibility interval says how far the organs can have moved. Upper bounds absorb any
systematic distortion as scatter (conservative).

**Caveats (speaker notes)**
- Slow deaths are not aligned to the hour (A109c H109c.1: R 0.097 against a null median 0.101, p 0.56), but alignment
  as strong as in ventilated donors is not excluded (80 %-power R 0.208).
- The arrangement survives removing the wall-time component (rho 0.750, A109c).
- Donor counts differ between A (12 body tissues, A293) and B (all donors of the group, A27): 103 against 106.
- Intermediate deaths look less similar than slow deaths (rho 0.624 against 0.784; core 0.721 against 0.831), but they
  have half the donors (50, 48 pairs, against 103, 66 pairs) and noise alone lowers rho at small n (79 fast deaths
  against the rest: median 0.778, A349). Both compatibility intervals include no movement and overlap (0-1.8 h against
  0-0.9 h; A374); the groups were not compared (Hardy 3 descriptive in A293). Not a difference. Expected from noise
  alone (fast deaths at the same n, post hoc): 0.796 (5-95 % 0.637-0.901) at 103 donors, 0.751 (0.491-0.871) at 50;
  both observed values lie inside (panel A, diamonds).

- Detection limit in intensive care: one organ must move about 6-7 h (6.75 h, 14 genes; core 5.6 h), or all organs
  with an SD of 2.25 h, to be seen (A349). "Kept" means no large disruption, not unchanged.
- Offsets in ventilated donors are compressed toward 0 / 12 h by a shared non-temporal axis (A109, A132b); removing it
  relaxes the compression from 0.58 to about 0.77 (A349 D349.4, approximate).
- Brain death: descriptive, inseparable from procurement in GTEx (A318 H318.2, A325, A327). Organ donors differ in RIN
  and ischaemic time (not balanceable).

**Source:** `exploratory/A314_vent_walltime/output/groups.csv`; `exploratory/A374_order_interval/output/noise_reference.csv`; `exploratory/A293_organisation_levels/output/part1/groups.csv`;
`results/figure_3/output/cues_lost/walltime_hosts.csv`; `exploratory/A109_coupling_nulls/output/a109c_walltime.csv`;
`exploratory/A374_order_interval/output/intervals.csv`; `docs/results_preliminary.md` A109c, A293, A374. A314, A317, A349.

## 4.1b Intensive care quiets the tissue clocks, most where they were loudest; brain death keeps their amplitude

**Bullets**
1. Ventilated donors' clocks are quieter: clock-specific amplitude 0.82 of fast deaths (Q* +0.404, above all 200
   matched sets), core 0.80; unchanged by age / sex matching (Q* +0.389, core +0.459) (A315 descriptive, A335).
2. The loss rises with the organ's loudness on the reference map (Spearman 0.77 over 11 tissues, gain of slide 3.1,
   post hoc; A349, descriptive: 0.83 with blood on the scale-free S = z2 V, erratum to its z2 readout): thyroid, oesophagus mucosa and heart lose most.
3. Coupling to the body's clock falls too, compared at a common host spread: median d* 0.08 over the 11 tissues, clear
   in thyroid (0.28, 90 % 0.18-0.38) and visceral fat (0.20, 0.06-0.32), and larger in louder organs (Spearman 0.79,
   post hoc; A349, descriptive).
4. Brain death keeps the amplitude of ventilated donors (1.02, 90 % 0.94-1.10, inside the registered equivalence margin
   0.84-1.20; A315 H315.1 FAIL, A335 part 3).
5. Amplitude falls with time on the ventilator among organ donors (longest vs shortest tertile Q* -0.236, 1st
   percentile of matched sets; A315, descriptive).

**Figure:** `figures/icu_amplitude/icu_amplitude.png` (`scripts/icu_amplitude.py`; wide). A. clock-specific clock
amplitude, one row per comparison: ventilated vs fast deaths (0.82; matched and responsive 95th percentiles 0.93 and
0.90), brain death vs ventilated (1.02, 90 % 0.94-1.10; shaded: equivalence margin 0.84-1.20), longest vs shortest
ventilation tertile among organ donors (0.89; companion percentiles not stored for this contrast). exp(Q* / 2) from the
clock-plane variance ratio. B. per tissue, host-referenced loudness lost in intensive care (log z2, fast over
ventilated, at a common spread, 90 % interval) against the tissue's gain on the reference map (slide 3.1); C. per
tissue, coupling lost (d*, fast minus ventilated at a common spread, 90 % interval), same x. Whole blood omitted (no
readable clock, slide 3.1b); x = reference gain rather than z2 (user, 2026-11-25; post hoc). 14 genes.

**Tags:** A315 descriptive contrasts (registered) and H315.1 FAIL (powered); A335 registered (core loop blind); A349
descriptive; amplitude display exp(Q* / 2) post hoc.

**Caveats (speaker notes)**
- The ventilated host spread is about half that of fast deaths (A349); per-tissue losses are compared at a common
  spread.
- Why not z2 on the x axis (as registered in A349): z2 is read in each tissue's host scale, and an anchor's
  leave-one-out host can differ. Muscle's host spread is 0.18 against about 3.5 for the other tissues, which inflates its
  z2 to 19.5 and made it look like the loudest organ losing little. On the gauge-free reference gain (1.34) it sits with
  lung and loses as much (0.24 vs 0.23). On the core loop muscle loses 0.41 and its coupling falls more clearly (d* 0.21,
  0.03-0.40).
- Raw couplings are not compared between groups (e.g. muscle 0.59 to 0.38): the host spread itself falls under
  ventilation, and coupling falls with spread alone; d* compares them at the pooled spread.
- Whole blood, omitted from panels B and C, "gains" (-1.95; A349): it has no readable clock (slide 3.1b). The registered Q* of
  panel A is a median over the 12 tissues including blood (a median, so one tissue moves it little).
- Raw clock-plane variance is lower in organ donors, as is the matched sets' (core Q -0.186): the brain-death
  equivalence is clock-specific, relative to random genes; also equivalent after age matching (A335).
- Deaths only; intensive care withdraws light-dark, meal and activity cycles together.

**Source:** `exploratory/A335_review_batch/output/{a315_core,matched,equivalence}.csv`;
`exploratory/A315_braindeath_gain/output/descriptive.csv`; `exploratory/A349_icu_limits/output/part2_tissues.csv`;
`presentation/output/organs/couplings.csv` (gain); `docs/results_preliminary.md` A315, A335, A349.

## 4.2 Removing a drive keeps a clock's amplitude; conflicting drives cost it

**Bullets**
1. Mouse liver without a feeding rhythm (Greenwell, 1/8 of the food every 3 h, 18 + 18 libraries) keeps 0.72 of its
   amplitude (90 % 0.66-0.79) and moves +1.6 h; if cues added, removal would leave 0.15-0.29 (calibrated on inverted
   feeding). H297.1 FAIL, above the interval: inputs set the liver clock's phase, the oscillator holds its amplitude
   (A297).
2. Inverted feeding moves each organ by its weight on the feeding drive: the liver almost fully (11.3-11.7 h), lung,
   heart and kidney 1.6-2.5 h (A322, A324).
3. Conflict costs amplitude, but not as a pure mix of the two drives: lung, heart and kidney keep 0.50-0.65, above the
   mix (H322.1 rejected, Q 31.9, p 0.0005; H324.1 PASS, Q 185); visceral fat loses more (0.11; H324.2 FAIL). Adipose, the
   organ with the most balanced weights, loses most in both datasets (white fat 0.32, visceral fat 0.11).
4. Suggestive only: the organs that weigh one drive in mouse (liver, heart) are the loudest in human, the balanced lung
   the quietest, adipose the exception (Spearman -0.6, 5 organs, below the registered -0.7; A358, "not seen").

**Figure:** `figures/holds_organisation/cue_removal.png` (`scripts/holds_organisation.py`, `cue_removal`; wide). Clock
faces: angle = phase shift against the reference arm (clockwise = later), radius = amplitude ratio; 90 % library-
bootstrap intervals. A. mouse liver: night-restricted (control), feeding inverted (Manella, Xin), feeding rhythm removed
(Greenwell), and the removal expected if cues add (open diamonds). B. mouse organs under inverted feeding: Manella
day-fed vs ad libitum (circles: liver, white fat, lung) and Xin day- vs night-restricted (squares: liver, heart, kidney,
visceral fat); dashed: the pure mix of the two drives, a straight chord from the reference to the fully shifted state
(Delta 11.28 and 11.42 h). Organ-system colours.

**Tags:** A297 H297.1 FAIL (registered, NOT BLIND for the shifts); A322 H322.1 rejected (NOT BLIND); A324 H324.1 PASS,
H324.2 FAIL (BLIND for heart, kidney, visceral fat); A358 descriptive, below its criterion.

**Caveats (speaker notes)**
- On the chord, an organ's expected state is the point at its fitted weight b (liver 0.84 / 0.72, white fat 0.65, lung
  0.37, heart 0.23, kidney 0.28, visceral fat 0.45); lung sits 26 Mahalanobis units above its point, kidney 148, heart 17;
  visceral fat 19 below.
- Different laboratories, strains and protocols; amplitudes are read in each dataset's own projector; Manella's
  reference is ad libitum, Xin's night-restricted feeding. Xin visceral fat and Manella white fat are different depots.
- Removal is tested in the mouse liver only; the cross-organ panel is inversion only.
- A358 is one line, never evidence: the balance b(1 - b) spans only 0.17-0.25 over the 5 organs.
- Human muscle under 12 h-shifted behaviour (the former panel B) is slide 4.3.

**Source:** `exploratory/A297_feeding_model/output/{tests,arms}.csv`;
`exploratory/A322_conflicting_drives/output/{ratios,fit,tests}.csv` and `xin/{ratios,fit,tests}.csv`;
`exploratory/A358_A360_gain_sources/output/{a358_organs,a358_tests}.csv`; `docs/results_preliminary.md` A297, A322,
A324, A358.
