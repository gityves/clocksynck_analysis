# Registered experiments (append only)

Amendments A1-A108c are registered in the archive (`$SCCC_ARCHIVE/docs/PLAN_internal_desync.md`); numbering continues
here from A109. Drafts of A109-A128 are in `provisory_plan_paper.md`; an amendment counts as registered only once it is
copied here, completed and committed, before its data are scored.

Each entry: `## A1xx (date, written before scoring | post hoc): title`, then data and data rule, estimand, hypothesis ids
with the pass rule (one- or two-sided, alpha, multiplicity), null and its size, power or MDE, reading of each outcome,
code folder (`paper/exploratory/A1xx_<slug>/`). Changes after registration are dated notes below the entry, labelled
post hoc when written after any result was seen.

## A109 (2026-10-23, written before scoring): is the coupling of organs without cues a clock?

Figure 3 rests on inter-organ coherence C staying at 0.41-0.46 in ventilated and slow deaths (archive A92, A93). C has
no clock-specific null, and A25 found a shared non-temporal post-mortem state inside the clock plane. A109 tests whether
C is clock-shaped (H109.1) and whether it reflects clocks that keep their phase relations rather than a shared push
(H109.2). Code: `paper/exploratory/A109_coupling_nulls/`.

**Data and data rule.** GTEx v10, the 12 A92 tissues; groups fast (Hardy 1-2), Hardy 3, Hardy 4, ventilated (Hardy 0).
Clock states as A92 (GTEx waveform projection, Poisson variance, centred within tissue and group). Equal n per group by
subsampling as A93. Group sizes, equal n and MDE: completed after the scout (note below), before scoring.

**Fixed confounds.** (1) Procurement: ventilated donors are beating-heart organ donors, fast deaths post-mortem
collections; every fast-vs-ventilated contrast is reported beside fast-vs-Hardy 4 (post-mortem, no daily cues).
(2) Amplitude: every C is reported with the per-tissue clock-plane variance (z2 V + 2 psi, A104) of the same group.

**Nulls (fixed).** (1) Pseudo-clock genes: 200 expression-matched 14-gene sets (`clocksync.paired.pseudo_waveforms`,
matched on mean log CPM per tissue, clock and D6 clock-output genes excluded) carrying the GTEx waveform loadings; C
recomputed for each. (2) Permuted-acrophase planes (secondary, conservative): 200 waveforms with gene acrophases
permuted (`Waveform.permuted`). (3) Shuffled pairing: donor labels shuffled within tissue (sanity, expected C about 0).

**H109.1 (primary, confirmatory family of five: H109.1, H109.2, H114.1, H115.1, H116.1; Holm, family-wise alpha 0.05).**
In ventilated donors, C_clock exceeds the pseudo-clock null. Statistic: C_clock minus the pseudo-clock null median;
one-sided empirical p (add-one, 200 sets); raw and Holm-adjusted p reported. Secondary: the same against the acrophase
null.

**H109.2 (primary, confirmatory family; Holm as above).** Offset preservation: for each tissue pair, the phase of the
noise-corrected complex covariance arg E[z_a conj(z_b)] in ventilated donors vs fast deaths. Statistic: Fisher-Lee
circular correlation over the 66 pairs; null: tissue labels permuted jointly in the ventilated group (10,000). Pass:
correlation > 0 with one-sided p (Holm), AND the median absolute pair phase of ventilated donors, folded to [0, 12 h],
is not closer to 0 or 12 h than in fast deaths (one-sided permutation p > 0.05 for "closer").

**Secondary.** (a) H109.1 in Hardy 3, Hardy 4 and fast deaths (fast = positive reference). (b) Excess coherence
(C_clock minus pseudo-clock median) per group with donor bootstrap; change fast -> ventilated. (c) Non-temporal
component: in fast deaths, regress each tissue's state on (cos, sin) of time of death; the residual cross-tissue
covariance defines the non-temporal shared direction; share of ventilated C along it vs orthogonal. (d) Processing
check: C on tissue pairs collected by different procedures, if GTEx annotations allow. (e) Procurement check: H109.1 and
H109.2 repeated in Hardy 4. (f) Circularity: pooled |E[z_a z_b]| / |E[z_a conj(z_b)]| over pairs, compared with 200
simulations each at the ventilated design (same donors, tissues, depths): rotating clocks with uniform phases at the
observed coupling; pure displacement of the observed size along the A27 glucocorticoid direction (17.7 h); mixtures at
25/50/75 % displacement; report the best-matching mixture fraction. Descriptive: major-axis share of host covariance;
Rayleigh and axial R of host phases; angle major axis vs glucocorticoid direction (mod 12 h).

**Reading.** H109.1 and H109.2 pass: organs share a clock-shaped state without cues and keep their phase relations;
Figure 3 stands. H109.1 passes, H109.2 fails: organs share a donor state, not a running clock; Figure 3 reframed, WP2
dropped from the main text. H109.1 fails: coupling-without-cues claim removed. Pass in ventilated with a failure in
Hardy 4: reported as procurement-dependent. GTEx time of death: aggregates only.

## A109b (2026-10-23, written before scoring): replication in Tabula Sapiens (secondary)

Tabula Sapiens donors are brain-dead ventilated organ donors with many organs each. Per donor x organ pseudobulks of
stromal-vascular cells (fibroblasts, endothelium, mural cells; A17 anchor families), organs with >= 5 donors.
Statistics: C against the pseudo-clock null (as H109.1) and offset preservation against the GTEx fast-death pair
relations for the organs both share (as H109.2). One-sided p < 0.05 each, secondary (not in the family). Feasibility
rule: >= 6 organs shared by >= 5 donors, otherwise descriptive. Gene matching on `feature_name`. Builder: a new
`common/datasets/tabula_sapiens_organs` (or `cells_census`). Organ and donor counts: completed after the scout, before
scoring. Code: `paper/exploratory/A109_coupling_nulls/`.

**A109b note (2026-10-23, written after the scout, before scoring).** Data: the held Census Tabula Sapiens group
pseudobulk (`$SCCC_DATA/atlas/tabula_sapiens/`, PROVENANCE.md there); builder
`paper/common/datasets/tabula_sapiens_organs/build.py` -> `$SCCC_DATA/pseudobulk/tabula_sapiens_organs.parquet`
(10x only, majority assay per donor, stromal-vascular families pooled per donor x organ, >= 50 cells; 14/14 clock genes
on `feature_name`). Attrition: 24 donors -> 21 with a stromal-vascular pseudobulk -> 12 donors, 49 pseudobulks in the 9
organs with >= 5 donors (adipose, bladder, exocrine gland, heart, thymus ["endocrine gland" in Census], musculature,
respiratory system, small intestine, spleen). The feasibility rule (>= 6 organs with >= 5 donors) is met, but 5 donors
carry 7-9 organs each and the others 1-2; 7 of 36 organ pairs share >= 5 donors. Variants (separate 10x strata, 30
cells, Smart-seq included) recover no donor; the registered rule is kept. Organs shared with the GTEx A92 tissues for
the offset test: adipose (subcutaneous; visceral as sensitivity), heart (atrial appendage), musculature (skeletal
muscle), respiratory system (lung): 6 pairs. With 5 multi-organ donors and 6 pairs the replication has low power; it is
reported as secondary with its MDE (added before scoring, from the A109 null at this design).

**A109 note (2026-10-23, written after the scout, before scoring).** Scout: `paper/exploratory/A109_coupling_nulls/output/scout_{counts,summary,h109_2_power}.csv`
(no C on the real clock genes was computed in the scout; the archive A92/A93 C values, already published, were reproduced
to validate the estimator).
- *Donors* (12 A92 tissues, >= 2 tissues, time of death present): fast 244, Hardy 3 50, Hardy 4 103, ventilated 377.
  Equal n = 103 (Hardy 4; A93 rule; Hardy 3 keeps its 50 donors, not equal-n). Tissue pairs with >= 20 shared donors at
  equal n: 66/66 (fast, Hardy 4, ventilated), 48 (Hardy 3). Universe = the 12 tissues (A93 used 107 donors including
  liver, adrenal and brain); C reproduces `figure_3/output/gtex_cues/a93_groups.csv` (Hardy 3 0.4494, Hardy 4 0.4127
  identical; fast 0.608 vs 0.600, ventilated 0.452 vs 0.462, the difference from universe and draws).
- *Estimator*: C as A92 (mean over pairs of the Poisson-corrected |complex correlation|, centred over shared donors,
  >= 20 shared), averaged over 20 equal-n draws; the same draws for data and every null set.
- *MDE, H109.1* (pseudo-clock null at n = 103; median, sd, excess needed for p < 0.05): ventilated 0.154, 0.010, 0.018
  (0.028 for p < 0.01); fast 0.140, 0.009, 0.018; Hardy 4 0.158, 0.013, 0.021; Hardy 3 0.203, 0.018, 0.032. The
  shuffled-pairing sanity null is expected at the estimator's finite-n floor of mean |correlation| (about 0.15), not 0.
- *Power, H109.2*: 66 pairs, label-permutation test: 80 % power at Fisher-Lee rho about 0.29 (p < 0.05).
- *Design choices*: H109.2 on full groups (primary), equal-n rows as sensitivity; group-centred states; no noise
  correction of the cross-covariance (noise independent across tissues). "Closer" statistic: median distance of the
  pair phases to the nearest of 0 / 12 h, p from 1,000 permutations of donor group labels (ventilated vs fast; Hardy 4
  vs fast). Sizes: 200 pseudo-clock, 200 acrophase, 200 shuffled-pairing sets; 10,000 joint tissue-label permutations;
  200 simulations per circularity scenario. Pseudo-clock genes exclude the 14 clock genes, `CLOCK_LIKE`
  (`clocksync/_core/controls.py`) and the 18 D6 clock-output genes. Bootstrap: donors resampled, then one equal-n
  draw. GC direction 17.68 h (`diag6c_gc_summary.csv`). Processing check (d): only whole blood uses a different
  collection kit (SMCAT Yellow vs Aqua), so it compares the 11 blood pairs with the 55 others. Holm adjustment over the
  family is applied in `results_preliminary.md` when all five are scored.
- *A109b MDE*: computed with the same null at the Tabula Sapiens design before A109b is scored.

**A109b MDE note (2026-10-23, written after A109 was scored, before A109b is scored).** Pseudo-clock null of C at the
Tabula Sapiens design (12 donors, 9 organs, 200 sets; pair C clipped at 1): pairs with >= 3 shared donors (primary, 36
pairs) null median 0.675, sd 0.060, MDE 0.104 (p < 0.05); >= 5 shared donors (sensitivity, 7 pairs) 0.531, 0.102,
0.187. The A92 minimum of 20 shared donors cannot apply here. Offset test: 4 GTEx-shared organs (adipose ->
subcutaneous, heart -> atrial appendage, musculature -> skeletal muscle, respiratory -> lung), 6 pairs, 24 exact label
permutations (minimum p 1/24 = 0.042); power at p <= 0.042 is 0.06-0.13 for couplings 0.3-0.8 (0.41 at 0.95). The offset
test is therefore reported but read as uninformative whatever its p; the C test is the A109b readout. Since A109 is
known, this replication is read for ventilated donors in a second cohort (H109.1-type C only). Script
`paper/exploratory/A109_coupling_nulls/scripts/compute_tabula_sapiens.py`.
<!-- src: paper/exploratory/A109_coupling_nulls/output/a109b_c_null.csv; a109b_power.csv -->

## A109c (2026-10-23, written after A109 was scored, before scoring; motivated by A109): is the Hardy 4 pass inherited from wall time, and is the 0 / 12 h push found in fatal COVID?

A109 found the organ phase relations of fast deaths preserved in Hardy 4 deaths (secondary) but collapsed towards
0 / 12 h in ventilated donors, with about half of the shared covariance behaving as a displacement along one axis
(PI comment 2026-10-23: as expected for a stopped clock, e.g. Bmal1 or Per knockout, but the push acts identically in
all organs). Hardy 4 donors are mostly hospital inpatients (archive A93: 85 of 107; 12 hospice, 4 ventilated) and may
keep meals (oral, enteral or intravenous) and ward light; their wall-time alignment R was 0.10 (fast 0.53, ventilated
0.15), never tested against a null. All tests below are secondary (own alpha 0.05, outside the confirmatory family).
Code: `paper/exploratory/A109_coupling_nulls/`.

**H109c.1 (Hardy 4 wall time).** R (A92/A93 wall-time alignment of the host phase with time of death, equal-n 103)
against 1,000 time-of-death permutations within Hardy 4. One-sided p; MDE reported. Reading: p >= 0.05 with an MDE
below 0.15 means no detectable entrainment; p < 0.05 means residual entrainment, and H109c.2 decides whether it
explains the offsets.

**H109c.2 (offsets without wall time).** H109.2 repeated in Hardy 4 after removing the wall-time component: each
tissue's state regressed on (cos, sin) of time of death within Hardy 4, residuals used. Same statistic, null (10,000
joint tissue-label permutations) and pass rule as H109.2 (rho > 0 one-sided p < 0.05 AND p_closer > 0.05), fast-death
relations as reference (unchanged). Sensitivities: without the 4 ventilated Hardy 4 donors; hospital inpatients only
(if the death-place variable allows; aggregates only). Reading: pass -> the Hardy 4 preservation is not inherited from
linear wall-time entrainment; fail -> it may be.

**D109c.3 (descriptive, fatal COVID, Delorey 2021, `covid_delorey`).** Four organs (heart, liver, kidney, lung), 6
pairs: pooled circularity ratio |E[z_a z_b]| / |E[z_a conj z_b]|, median folded pair phase distance to 0 / 12 h, and
Fisher-Lee correlation with the GTEx fast-death relations of the matching tissues; each compared with 200 simulations
at the Delorey design (rotating vs pure displacement vs 50 % mixture, as A109) and with the GTEx fast and ventilated
values on the same tissues. Descriptive, because the archive organ gate is mixed (H105.2: p_pairing 0.284,
p_acrophase 0.032) and the coupling is 0 in 3 of 4 organs; a failed gate is uninformative.

**D109c.4 (descriptive, fatal COVID lung, Melms 2021, `covid_melms`).** The same circularity ratio and folded distance
across the lung's cell classes (gate passes, H105.1), COVID donors vs the 7 surgical controls, with a donor
permutation p for the difference in ratio (1,000) reported but not tested. It asks whether a shared push across cell
types appears within a lung of fatal COVID (mostly ventilated).

**A109 note (2026-10-23, post hoc, after A109 and A109c were scored): framing decided by the PI.** The registered
reading for "H109.1 passes, H109.2 fails" was "organs share a donor state". The PI chose instead: organs keep a shared
clock-shaped state without daily cues (H109.1, confirmatory); when cues are lost through illness (Hardy 4, registered
secondary check (e); H109c.2) the organs keep their phase relations; in ventilated, mostly brain-dead organ donors a
shared perturbation overrides those relations (H109.2 fails, reported as a failure of a confirmatory hypothesis). This
framing is post hoc; the Holm table keeps H109.2 as failed.

## A132 (2026-10-23, written after A109/A109c were scored, before scoring): brain death or ventilation? The collapse of organ phase relations in ventilated donors

PI hypothesis: brain-dead donors have probably lost the SCN, whereas Hardy 4 donors keep it with perturbed zeitgebers.
Brain death co-occurs with donor management (methylprednisolone, thyroid hormone, vasopressin, vasopressors), so loss of
the SCN and these drugs cannot be separated in GTEx; this amendment separates brain death (with its management) from
ventilation and ICU care. All hypotheses secondary (own alpha 0.05). Data: GTEx v10, A109 states and tissues; protected
phenotypes, aggregates only. Code: `paper/exploratory/A132_brain_death/`.

**Groups.** Within Hardy 0: organ donors (`COHORT` organ donor / OPO, brain dead) vs ventilated post-mortem cases
(`COHORT` post-mortem, not declared brain dead). Group sizes and pairs with >= 20 shared donors: after the scout.
Feasibility rule: >= 40 donors and >= 40 of 66 pairs in the smaller group, otherwise descriptive.

**H132.1 (collapse in brain death).** The A109 readouts per group, fast-death relations as reference: Fisher-Lee rho of
the pair offsets (10,000 joint tissue-label permutations) and the folded distance to 0 / 12 h (p_closer, 1,000
donor-label permutations vs fast). Statistic: difference in folded distance, organ donors minus ventilated post-mortem
cases, with a donor-label permutation p (1,000, one-sided, prediction < 0: organ donors closer to 0 / 12 h), at equal n.
Reading: organ donors collapse and ventilated post-mortem cases do not -> the collapse goes with brain death and its
management (SCN loss or donor drugs); both collapse -> it goes with ventilation / ICU care; neither differs -> no
separation at this power (MDE reported).

**H132.2 (dose).** Within Hardy 0 organ donors, ventilation duration (protected phenotype) top vs bottom tertile:
difference in folded distance, donor-label permutation p (1,000, one-sided, prediction: longer -> closer to 0 / 12 h).
MDE after the scout.

**D132.3 (descriptive).** Long-term steroids (protected phenotype) yes vs no within Hardy 0 and within Hardy 4: folded
distance and the glucocorticoid-direction share (A109 non-temporal output), if >= 15 donors per arm.

**A132 note (2026-10-23, written after the scout, before scoring).** Scout:
`paper/exploratory/A132_brain_death/output/scout_{variables,cohort_hardy,groups,tertiles,steroids,mde}.csv` (full data;
null spreads only, no observed contrast). An end-to-end `--subsample` check also ran; its observed values were deleted
and not read by the lead. Variables: `COHORT`, `DTHVNT`, `DTHVNTD` (hours), `MHSTRDLT`. Analysis universe (>= 2 tissues,
time of death): Hardy 0 organ donors 298, Hardy 0 post-mortem 79 (all `DTHVNT` = 1), fast reference 244.
- **Feasible.** Smaller group (ventilated post-mortem) 79 donors, 66 / 66 pairs with >= 20 shared donors (organ donors at
  equal n = 79: 64 pairs). Fixed choices: pairs with < 20 shared donors in either group are masked for that comparison;
  the equal-n contrast uses one fixed draw of 79 organ donors (seed in the script).
- **MDE (80 %, one-sided 0.05, 2.49 x null SD).** H132.1: 0.46 h (null SD 0.186 h). The A109 ventilated-minus-fast
  difference in median folded distance was 0.42 h (0.60 vs 1.02 h), so a full collapse in organ donors with none in
  post-mortem cases is detectable with power about 0.7; a partial split is not. A null H132.1 is uninformative, not
  evidence that the two groups collapse alike.
- **H132.2.** Tertile cuts 61.0 h and 101.0 h (median 79.2 h, no missing durations); bottom 101, top 100 donors, 66 pairs
  each; MDE 0.46 h.
- **D132.3 not reportable:** long-term steroids yes in 8 Hardy 0 and 5 Hardy 4 donors (< 15 per arm). Dropped.

## A110 (2026-10-23, written after A109/A109c were scored and the PI chose the framing (Q1), before scoring): what keeps organs coupled when daily cues are disrupted by illness?

Replaces the WP2 draft in `provisory_plan_paper.md`, which assumed ventilated donors kept their phase relations
(H109.2 failed there). Groups follow the PI's reading: Hardy 4 (slow death after long illness; SCN presumably intact,
zeitgebers perturbed, drugs present) is the main group, Hardy 3 a secondary group, fast deaths (Hardy 1-2) the reference
and Hardy 0 organ donors (brain dead, SCN presumably lost) the contrast for H110.4. Drugs and illness cannot be
separated from cue loss in GTEx; every reading below is about this state, not about cues alone. Data, tissues, states and
C as A109 (`paper/exploratory/A109_coupling_nulls/scripts/_helpers.py`); protected phenotypes aggregates only. Code:
`paper/exploratory/A110_mechanism/`. H110.1 is the primary of A110 (secondary for the paper); the others are secondary,
with Holm within H110.2.

**H110.1 (does coupling decay with time ill?).** Hardy 4 donors split by the interval from onset of the immediate cause
of death to death (`DTHCODD`, unit `DTHCODDU`; the scout checks completeness; if < 60 % of Hardy 4 donors have it, the
first underlying cause interval `DTHFUCODD` is used instead). Statistic: C in the top minus the bottom tertile (median
split if either tertile has < 40 pairs with >= 20 shared donors), equal n, one-sided donor-label permutation p (1,000;
memory predicts top < bottom). MDE after the scout; a memory-model C drop at the top-tertile median interval for a period
spread of 0.3 h/day is simulated before scoring (as A92's power note). Reading: decay -> coupling decays without regular
cues (memory); no decay with the MDE below the simulated drop -> coupling is maintained; otherwise uninformative. The
interval is not the time without cues (hospitalisation dates are not recorded): stated as a limitation.

**H110.2 (signatures of a coupled network; Holm over a-b).** In Hardy 4, fast deaths as reference.
- (a) Strong clocks carry the network: Spearman(own amplitude A, coupling c) over peripheral tissues (A104 statistic),
  one-sided rho > 0, tissue-label permutation p.
- (b) Stroma: retained coupling per tissue (Hardy 4 / fast, A33 machinery, at common spread d*) against the Eraslan
  stromal-vascular fraction; one-sided Spearman > 0, permutation p.
- (c, descriptive) Hubs: the WP2-draft split-half hub removal; Hardy 4 is too small for a confirmatory split.

**H110.3 (temperature, exclusion).** As the WP2 draft (CIRBP, RBM3; HSPA1A, HSPA1B, HSPH1, DNAJB1; residualised on
ischaemic time and its square per tissue), in Hardy 4: excluded if the drop in C after regressing the scores out is
< 0.05 with upper bound < 0.05 and not above 200 expression-matched random 6-gene scores.

**H110.4 (SCN).** Coupling of Brain - Hypothalamus minus the median other brain region, each to the peripheral host
(anchored estimator, common d*), in Hardy 3-4 pooled vs Hardy 0 organ donors. PI hypothesis predicts the difference is
larger in Hardy 3-4; one-sided donor-label permutation p (1,000). Feasibility rule: >= 20 donors with hypothalamus and
>= 3 other brain regions in each group, otherwise descriptive (fast vs Hardy 3-4 as a fallback, labelled descriptive).
Caveats in the text: GTEx hypothalamus dissections need not contain the SCN, and brain death alters brain tissue
itself, so a pass supports but does not show an SCN-driven body.

**H110.5 (feeding, exclusion).** Digestive (colon transverse and sigmoid, stomach, oesophagus mucosa and muscularis,
liver, pancreas) vs other A109 tissues: within-digestive minus digestive-to-other coherence, Hardy 4 minus fast,
group-label permutation; excluded if not > 0 at p < 0.05.

**Reading.** The Figure 3 panel reports each explanation with its verdict. The abstract sentence may say that organs of
people dying after long illness remain coupled clocks with their phase relations, and (if H110.1 shows no decay and
H110.2a or b passes) that the coupling is maintained rather than remembered. No sentence claims coupling without the
brain: Hardy 4 donors have it.

## A132b (2026-10-24, POST HOC: chosen after seeing A132's per-group readouts; requested by the PI): do ventilated post-mortem cases keep the fast-death phase relations better than brain-dead organ donors?

A132 showed both Hardy 0 groups pushed toward 0 / 12 h, with offset preservation descriptively only in the ventilated
post-mortem cases (rho 0.690 vs 0.161). Statistic: Fisher-Lee rho vs fast deaths, ventilated post-mortem minus organ
donors, at equal n (the same 79 organ donors as A132, seed 132), over the pairs with >= 20 shared donors in fast deaths
and in both groups (fixed mask); null: donor labels permuted between the two groups (1,000), two-sided p. Reported as
post hoc, with both per-group readouts (offsets and folded distance) described side by side. Code:
`paper/exploratory/A132_brain_death/scripts/compute_posthoc_rho.py`; output `a132b_rho_difference.csv`.
