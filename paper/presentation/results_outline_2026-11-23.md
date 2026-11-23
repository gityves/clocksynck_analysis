# Results section outline (2026-11-23, built slide by slide with the user)

Prompt: `docs/prompts/presentation_results_2026-11-23.md`. Framing: "loudness, not time" (`docs/ABSTRACT_DRAFT_2026-11-22.md`,
`docs/PAPER_PLAN_2026-11-19.md` section 11). Slides kept from `results.md` are referenced, not repeated. Tags: blind / not
blind / descriptive / post hoc. Every number below was read from the output table named in its **Source** line.

---

# 1. Validation: the readout sees time, and sees it move

| slide | status | figure |
|---|---|---|
| 1.1 GTEx organ map from working-hours deaths | kept (`results.md` R1.1) | `figures/gtex_organ_map/` |
| 1.2 Living patients: breast, heart, adipose | kept as is (`results.md` R1.3; user, 2026-11-23) | `figures/living_patients/` |
| 1.3 Baboon atlas, 63 tissue offsets | kept (`results.md` R1.4; user, 2026-11-23) | `figures/baboon_offsets/` |
| 1.4 Timed living skin and muscle, and hospital hours | **new** (below) | `figures/living_timed/` |
| 1.5a Positive controls, cell types: offsets and the knockout | **new** (below) | `figures/celltype_controls/` |
| 1.5b Positive controls, organs: the readout moves when the clock is broken or shifted | **revised** (below) | `figures/positive_controls/` |

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
