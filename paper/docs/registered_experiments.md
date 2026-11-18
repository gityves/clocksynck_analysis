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

## A133 (2026-10-24, written after A109/A132b, before any GTEx scoring; method): separating running clocks from a shared displacement

Motivation (post hoc relative to A109/A132): ventilated pair phases are pulled toward 0 / 12 h. Clock time spread over
the cycle contributes to the Hermitian cross-covariance E[z_a conj z_b] only (its pseudo-covariance E[z_a z_b] vanishes
for uniformly covered phases); a shared donor-level displacement s_d u_a (s real, varying between donors; a constant
shift cancels on centring) contributes s-variance u_a conj(u_b) to the first and u_a u_b to the second. Estimator: fit
the axial part from the pseudo-covariance matrix (rank-1 Takagi factor u), subtract it from the Hermitian
cross-covariance, read coupling and offsets on the remainder ("circular part").
- **V133.1 (simulation, gate for the rest).** Running clocks with known offsets plus displacements (mixture fraction 0 to
  0.8), at GTEx ventilated sizes (n 79 and 298, 12 tissues, real noise): the circular part recovers the offsets
  (Fisher-Lee rho with truth >= 0.8 at fraction 0.5, n 79) where the raw estimator does not; with no displacement it
  agrees with the raw estimator (rho >= 0.95). Confounder: host phases concentrated (Rayleigh R up to 0.3) with no
  displacement must not be "corrected" into different offsets (rho >= 0.9). Fails -> A133 stops, the method stays.
- **H133.1 (secondary; only if V133.1 passes).** On the circular part, Fisher-Lee rho vs fast deaths in ventilated
  post-mortem cases and in organ donors (A132 groups, equal n, same permutations as A132). Reading: organ donors' offsets
  recovered (one-sided p < 0.05) -> a displacement masked running clocks with kept relations; not recovered with the
  post-mortem cases recovered -> relations lost in brain death (A132b). Folded distance on the circular part reported.
- Code in `paper/exploratory/A133_axial_decomposition/`; a clocksync function is proposed to the PI only after V133.1.

**A110 note (2026-10-24, written after the scout, before scoring).** Scout:
`paper/exploratory/A110_mechanism/output/scout_*.csv` (counts and null spreads; subsample smoke outputs deleted unread).
- **H110.1 infeasible, not scored.** The interval is recorded for 16 of 103 Hardy 4 donors (`DTHCODD` 0.155,
  `DTHFUCODD` 0.146; Hardy 3 0.30, fast 0.32, Hardy 0 0.85); a median split gives 8 vs 8 donors and no pair with >= 20
  shared donors. With complete data the median split would have had MDE 0.143. A110 has no primary; memory vs
  maintained coupling is not testable in GTEx Hardy 4. (Hardy 0 has the interval but carries the A132 displacement.)
- **H110.2 feasible.** 32 peripheral tissues; 23 with >= 20 Hardy 4 and fast donors (floor 20, as A33); Eraslan
  stromal-vascular fraction for 8 tissues, 7 at >= 20 donors, so H110.2b is a Spearman over 7 tissues (low power).
  Retained coupling = Hardy 4 / fast at d* from `Dataset.contrast`. H110.2c descriptive.
- **H110.3 feasible.** All six genes present; ischaemic time complete for >= 98.8 % of Hardy 4 samples per tissue
  (median 820 min). Fixed design: cold and heat scores regressed out jointly; donor bootstrap with the fit fixed on the
  full sample; random sets matched within 0.25 log CPM per gene.
- **H110.4: the registered contrast is impossible:** GTEx organ donors have no brain samples (0 of 301; 10 Hardy 0
  donors with hypothalamus in all). The SCN-loss reading cannot be tested on GTEx brain. Registered fallback, descriptive:
  fast (145 donors with hypothalamus, >= 3 other regions, >= 2 peripheral tissues) vs Hardy 3-4 (78); statistic at d*:
  hypothalamus minus the median other region, cervical spinal cord excluded.
- **H110.5 feasible, weak.** Digestive set as registered (7 tissues; beyond the 12 A109 tissues); pairs with >= 20 shared
  donors in Hardy 4 and fast at equal n (103): 4 within digestive, 40 digestive-to-other; null SD 0.088, MDE 0.218.

## A111, A112, A129 (WP3; 2026-10-24, written before any data of these studies were downloaded or scored): feeding and genetic controls

Text as drafted in `provisory_plan_paper.md` WP3, copied verbatim below. Secondary for the paper. Scout notes (accessions, designs, group sizes, MDEs) follow before scoring; A129 is registered as a hypothesis only after its scout.

A49 (Deota) could not serve as a positive control because night-time feeding produced almost no desynchrony. Daytime and inverted feeding do. Arrhythmic feeding gives a timed truth for the cue/locked distinction.

*Base code:* mouse orthologue convention (A48) `results/figure_4/scripts/resilience/compute_mouse_ageing.py`; Deota builder `common/datasets/mouse_ageing_deota` (A49); timed truth `clocksync.timed` (`Dataset.timed`, `timed_loadings`, `cosinor`). New builders go to `common/datasets/<name>/build.py` with a `PROVENANCE.md` under `$SCCC_DATA/<name>/`.

**Datasets (scout first; take accessions from the papers).**

| dataset | design | tissues | role |
| --- | --- | --- | --- |
| [Manella et al. 2021, Nat Metab](https://pmc.ncbi.nlm.nih.gov/articles/PMC7611072) | ad libitum vs daytime feeding, liver-clock mutants, around the clock | liver, WAT, lung, quadriceps, kidney, heart | primary positive control |
| [Xin et al. 2021, iScience](https://pmc.ncbi.nlm.nih.gov/articles/PMC8050734) | day- vs night-restricted feeding, 1 week, 24 h | several peripheral tissues (TPM only) | replication; pseudo-counts as V25.5 |
| [Greenwell et al. 2019, GSE118967](https://omicsdi.org/dataset/geo/GSE118967) | ad libitum, night-only, arrhythmic (1/8 of food every 3 h) | liver | lock validation |

**A111, positive control (Manella primary, Xin replication).**

- Units: one mouse per time point and condition; tissues as compartments; host anchored on all tissues except liver (all except the tissue read, leave-one-out). Human GTEx waveform with mouse orthologues (A48 convention); ZT folded to 24 h.
- Timed truth first: per tissue and condition, cosinor phase of the clock-plane state on ZT; truth change = phase(daytime feeding) minus phase(ad libitum) per tissue, relative to the tissue mean.
- H111.1 (offsets): the untimed change in tissue offsets relative to the host correlates with the timed change, circular correlation (Fisher-Lee) > 0, exact permutation p < 0.05 over tissues.
- H111.2 (coupling): liver coupling to the non-liver host falls under daytime feeding (untimed d < 0), in the direction of the timed liver-vs-rest R^2 change. Permutation of condition labels within ZT.
- Informativeness rule (as A49): the timed truth must show a liver phase shift >= 3 h relative to the other tissues; otherwise report uninformative.
- H111.3 (secondary, a positive control for a break): in the liver-clock mutant arm under ad libitum feeding, liver coupling to the non-liver host is lower than in wild type (untimed d < 0 at a common host spread, d\*), or the liver fails the cross gate while the other tissues pass; condition labels permuted within ZT, one-sided p < 0.05. Check the mutant's design in the paper first (which clock gene, which cells) and write it into the amendment.
- Descriptive: does the rest of the body stay together without a liver clock (C among non-liver tissues, mutant vs wild type)?

**A112, lock validation (Greenwell).**

- Single tissue, so the host is the liver's own clock-plane state per sample. For each gene, lock = |coefficient on the host| in arrhythmic feeding divided by the same in ad libitum (A81 definition, both timed and untimed versions).
- Gene sets fixed now: mouse orthologues of the GTEx A81b cue set and locked set (40 each), core clock genes as reference.
- H112.1: cue-set orthologues lose more lock under arrhythmic feeding than locked-set orthologues (median log lock ratio, cue minus locked < 0), gene-label permutation p < 0.05.
- Truth check: Greenwell's own rhythmicity calls (rhythmic in ad libitum, lost in arrhythmic) should be enriched among cue orthologues.
- Reading: a pass shows the GTEx-derived cue/locked split corresponds to feeding-driven vs clock-driven rhythms in a controlled experiment, across species.

**A129, a tissue without its clock (scout, then register).** No work package so far shows that the method reports a break when a tissue's own clock is removed.

- Scout public multi-tissue transcriptomes of tissue-specific clock knockout or reconstitution mice with >= 2 tissues per mouse: e.g. Koronowski et al. 2019 (liver-only Bmal1 reconstitution), Welz et al. 2019 (epidermis-only reconstitution), and tissue-specific Bmal1 knockouts. Verify accessions and designs; use only open matrices.
- If found (>= 6 mice per genotype with the tissues of interest): H129.1, the tissue lacking its clock has lower coupling to the host of the other tissues than in controls (d\* < 0), or fails the cross gate, while a companion tissue with an intact clock does not change (specificity on the same permutations). Timed truth where ZT is recorded. Register after the scout, before scoring.
- If none: state it as a limitation; A111 H111.3 is the only genetic control.

**Addition (2026-10-24, after A109/A132/A133): clock-loss readouts on the knockout tissues (A129 and the A111 liver-clock mutant arm, descriptive).** Where a tissue lacks its clock, report the A109 offset readouts against the wild-type relations (Fisher-Lee rho, folded distance to 0 / 12 h, Spearman of the signed pole deviations, as A132b): the expected signature of lost clock relations (order lost) for comparison with GTEx brain-dead organ donors.

## A113, A123, A130, A131 (WP4; 2026-10-24, written before any data of the new cohorts were downloaded or scored; A123: the 8-tissue benchmark window results have been seen): internal time in living humans

Text as drafted in `provisory_plan_paper.md` WP4, copied verbatim below. A113 is primary of WP4 (secondary for the paper); A123 primary for the Figure 1 organ map; A130 secondary; A131 a scout (exploratory). Scout notes follow before scoring.

V25.5 (breast surgery, 42 patients) is the only living-hospital validation of per-patient internal time. A113 adds living timed cohorts, including one where the feeding rhythm is flattened.

*Base code:* A28 E1 projector `results/extended_data/scripts/gate/compute_extended_clock.py`; A35 within-subject arm `results/figure_2/scripts/gene_rhythms/compute_blood_clock.py`; A25 window estimator `results/figure_1/scripts/windows/{compute_internal_time,_internal_time}.py`; V25.5 `.../windows/compute_breast_timed.py`.

**Datasets (scout first).**

- GSE108539 (timed skeletal muscle) and GSE104674 (timed adipose): listed in the plan as timed human bulk sets, never used. Confirm subjects, sampling times, repeated biopsies per subject, platform and count type.
- Geneva human muscle with hourly isocaloric feeding (Perrin et al. 2018; Loizides-Mangold et al. 2017): serial biopsies every 4 h. Use only if a transcriptome matrix is public.
- Smith et al. 2024 (JCEM, hourly liquid meals, muscle every 4 h): use only if genome-wide or at least the 14 clock genes are public.

**Estimators.**

- *Single sample, single tissue:* clock-plane phase from the 14-gene GTEx projector (A28 E1); error after a common offset.
- *Within subject:* states centred on the subject mean (A35 arm c); ordering of each subject's biopsies.
- *Cohort window:* A25 axis and curvature on subsets mimicking hospital windows (08-18 h, 09-17 h) built from the timed samples, flip resolved by the hospital prior (V25.5 rule).

**Primary (H113.1), per dataset.** Circular correlation of estimated phase with recorded clock time > 0.3, permutation p < 0.05 (subject-level permutation when biopsies are repeated). Pass in >= 2 datasets.

**Secondary.**

- Median per-sample error after the common offset, against V25.5 (1.1 h) and GTEx single-tissue (2.6 h).
- Window subsets: Spearman of the tangent coordinate with clock time; direction correct under the hospital prior.
- Hourly-feeding cohort: does the muscle clock order samples as well as under normal meals? Descriptive, the human counterpart of ICU feeding.
- Common offset (internal frame minus wall clock) per dataset, next to V25.5 (-2.2 h).

**Reading.** Passing in two or more living cohorts turns internal-time estimation from biopsies into a validated readout and supports the translational paragraph of Figure 3.

**A123, the GTEx organ map from working-hours deaths only (written before scoring; the 8-tissue window results of the benchmark have been seen).** Reproduce the CHIRAL organ figure using only donors who died in working hours, with times hidden.

- *Base code:* `results/figure_1/scripts/method/compute_gtex_benchmark.py` (benchmark truth, CHIRAL runner, windows), `common/benchmark/chiral.py`, `results/figure_2/scripts/organs/compute_organs_timed.py`.
- Donors: GTEx v10 fast deaths (Hardy 1-2) who died 08-18 h; every tissue with >= 20 such donors.
- Estimates: untimed reference-free tissue phase offsets and relative amplitudes; per-donor internal time (A25 tangent coordinate, flip by the hospital prior).
- Truth (i): per-tissue cosinor on recorded time of all fast deaths, reduced by complex SVD (the benchmark truth). Truth (ii): CHIRAL's published per-tissue phases and amplitudes (naef-lab/CHIRAL, commit 1f78a9a).
- H123.1: offsets vs truth (i), circular correlation > 0.7, tissue-permutation p < 0.05. H123.2: relative amplitudes vs truth (i), Pearson > 0.7. H123.3: donor ordering, Spearman of the tangent coordinate with recorded time > 0.4 within the window, with the direction correct.
- Secondary: the same against truth (ii); CHIRAL (pychiral) run on the same 08-18 h donors as the comparator that assigns phases.
- Deliverable: a CHIRAL-style organ figure built from working-hours donors only, for Figure 1.

**A130, internal time against known chronotype epidemiology.** Chronotype advances with age after early adulthood and is earlier in women until about 50 (Roenneberg et al. 2004; verify before citing). In GTEx the age advance was already seen (A94: about 0.2-0.5 h per decade, peripheral tissues) and the sex phase difference was estimated (A99: +0.35 h, SE 0.31): both are reported as consistent with that epidemiology, labelled post hoc.

- H130.1 (secondary, pre-registered replication outside GTEx): per-subject internal minus wall time against age, pooled over the living cohorts of A113 and V25.5 that record age (cohort as stratum, A25 internal time with the hospital prior where sampling is in working hours, timed phase where times span the day). Slope < 0 (earlier with age), one-sided permutation of age within cohort, p < 0.05. Feasibility: >= 30 subjects with age and an age range >= 25 years, otherwise descriptive.
- Descriptive: women minus men, same estimand, by age above or below 50.

**A131, sun time within time zones (scout, exploratory).** Human chronotype follows sun time within a time zone (Roenneberg et al. 2007; verify). GTEx samples carry a collection-centre code (`SMCENTER`). Scout whether public GTEx documentation maps the codes to locations, and whether the donor's place of death is near the centre (OPO catchment). If yes: centre longitude offset from its time-zone meridian against the donor's host phase minus time of death (fast deaths; A93 H93.3 machinery), centre as the unit, descriptive with interval. Centres differ in procedures, so any trend is a candidate, not a finding. Aggregates per centre only.

**WP4 note (2026-10-24, written after the scout, before scoring).** Scout: `paper/exploratory/WP4_scout/output/scout_wp4.csv`,
`gtex_window_counts.csv` (aggregates; no hypothesis statistic computed).
- **A113 datasets: two independent new cohorts.** Perrin 2018 (GSE108539), Loizides-Mangold 2017 (lipidomics on the same
  biopsies) and Smith 2024 (same cohort summary; matrix not public, not requested) are one Geneva cohort.
  - *Muscle, GSE108539:* 10 subjects (9 M / 1 F), 6 biopsies 4-hourly over 24 h (one subject 3), hourly isocaloric
    meals; public matrix log2 edgeR RPKM. Fixed: exon layer, array-type input (as GSE39445), 13 clock genes (ARNTL
    absent from the public matrix; the projector uses the 13 present genes, weights of the 14-gene projector without
    ARNTL); subject-level permutation; the 3-biopsy subject kept.
  - *Adipose, GSE104674:* 12 men (6 lean, 6 T2D), 4 biopsies within 20 h (ZT 15:30, 0:15, 5:45, 11:15), raw counts, 14
    genes. Fixed: scored against ZT (the common offset absorbs the unknown wake-up clock time); lean and T2D pooled,
    permutation within subject group; the internal-minus-wall offset is not available for adipose.
  - H113.1 "pass in >= 2 datasets" therefore requires both cohorts to pass; V25.5 is the prior validation, not counted.
- **A123 feasible.** 244 fast deaths with time of death, 111 died 08-18 h; 40 of 48 tissues have >= 20 window donors
  (the 8 benchmark tissues 79-98; whole blood 99, included by the registered >= 20 rule). Truth (ii): CHIRAL's per-tissue
  tables are not local (only on a Drive link); re-derived by harmonic regression of the clock genes on CHIRAL's
  per-sample CT (`sample_metadata.txt`, commit 1f78a9a), as `compute_chiral_age.py`; no download.
- **A130 feasible.** Pool V25.5 (42 women, 33-93 y, per-subject age from PNAS SI Table S1) and GSE104674 (12 men,
  40-71 y): 54 subjects, range 60 y; cohort as stratum. Muscle has no per-subject age and is excluded. The sex
  descriptive is dropped: sex is fully confounded with cohort. V25.5 ages are added in its builder from Table S1.
- **A131 dropped (infeasible).** Public GTEx dictionaries define `SMCENTER` only as a code (B1, C1, D1); no locations;
  and fast deaths with time of death come from effectively two centres (B1 128, C1 115).

**WP3 note (2026-10-24, written after the scout, before scoring).** Scout: `paper/exploratory/WP3_scout/output/scout_wp3.csv`
(designs and counts; no statistic, informativeness not computed). Fixed choices:
- **A111 primary, Manella 2021 (GSE159135):** liver, WAT and lung only (the other tissues were qPCR); 4 arms (Alb-Cre
  controls and hepatocyte Bmal1 KO "BLKO", ad libitum / daytime feeding) x 12 ZT x 2 mice. Labels `<ZT><A|B>` are
  assumed to be the same mouse across tissues (not stated by the authors). Host for liver: WAT + lung (minimum host 2
  tissues). Informativeness: liver shift >= 3 h relative to WAT and lung. Depth = sum over the authors' retained genes
  (pre-filtered 3' counts).
- **A111 replication, Xin 2021:** TPM pseudo-counts (x 30 M, no length correction); liver, heart, kidney, VAT; ZT24 folded
  to 0 and kept as separate mice; the kidney-only mouse `N24_5` enters kidney only. Mouse pairing across series assumed.
- **H111.3 = A129 primary (merged).** Manella BLKO vs Alb-Cre within ad libitum (24 vs 24 mice): liver coupling to the
  WAT + lung host, d* < 0 one-sided, condition labels permuted within ZT (1,000). Specificity companion on the same
  permutations: WAT coupling to lung and lung to WAT (hosts without the knockout liver), BLKO minus Alb-Cre, reported
  with intervals; specific if the liver drop exceeds both companions'. Daytime-feeding arms descriptive. Arntl mRNA is
  not a knockout readout in these 3' data (still detected).
- **Expression QC:** a gene failing the median-count >= 5 rule is dropped in that tissue and dataset (CIART in Manella
  liver; CIART and NPAS2 in Greenwell liver); tissues are kept.
- **A112, Greenwell 2019 (GSE118967):** primary contrast arrhythmic (AR) vs ad libitum (LB); the partly adapted arm (DR)
  and night-restricted (NR) descriptive. The rhythmicity-call truth check uses the paper's supplementary tables if public,
  otherwise it is dropped.
- **Clock-loss readouts (addition, descriptive):** Manella has 3 pairs only, so instead of the A132b pair statistics:
  liver offset to the WAT + lung host, BLKO vs Alb-Cre, untimed and timed. Further descriptive sets: VIP-Cre Bmal1 KO
  (GSE296764; SCN-specific knockout, SCN + liver + heart, CT0 / CT12 only, 6-8 mice per arm): liver-heart coupling and
  offset, knockout vs wild type, as a mouse analogue of the brain-dead question; whole-body Bmal1 KO (Jouffe 2022,
  GSE150102, single ZT12): coupling only. Other candidates (Koronowski 2019, Welz 2019, liver-only sets) are
  single-tissue and not used.

## A134 (2026-10-24, written after WP0 item 2 and the PI's answer to Q8, before scoring): is the malignant "coherence loss" non-clock variation across tumours?

WP0 item 2: malignant pseudobulks are deeper than normal epithelium, yet A87 coherence is lower (0.23 vs 0.78) while
per-cell D is equal (3.89 vs 3.90) and off-plane clock-gene variance across tumours is about 3-fold higher. Candidate:
each malignant cell keeps a coherent clock; across tumours, clock-gene expression carries large non-clock variation
that hides the shared phase structure. Same 79 patients (A87 single-cell set). Secondary; the Figure 5 statement
depends on it (PI, Q8). Code: `paper/exploratory/A134_cancer_coherence/`.
- **H134.1 (primary of A134): copy number.** Per patient and clock gene, the copy-number state of the gene's locus in
  malignant cells (the inference used by `compute_cancer_split.py --cnv`, or author annotations; fixed in the scout
  note). Regress each clock gene's between-patient malignant pseudobulk (log, centred within study x stratum) on its
  own copy-number state, recompute A87 coherence on the residuals. Statistic: coherence gain (adjusted minus raw) in
  malignant cells minus the gain in normal epithelium given the same regression (specificity companion; normal
  epithelium's own copy number is expected near neutral), patient bootstrap (1,000), one-sided > 0.
- **H134.2 (secondary): non-clock programmes.** Same as H134.1 with patient-level malignant scores of proliferation
  (A24 cycling fraction), MYC targets and hypoxia (MSigDB hallmark sets, fixed versions in the scout note) instead of
  copy number. One-sided > 0.
- **D134.3 (descriptive): tumour-specific gene-phase relations.** Per-cell D with the common GTEx gene phases vs with
  patient-specific gene phases fitted on half of the patient's cells and evaluated on the other half (cross-fit),
  malignant vs normal.
- Reading: H134.1 or H134.2 passes -> restate Figure 5 as "malignant cells keep a coherent clock; tumour-to-tumour
  variation from copy number / oncogenic programmes breaks the shared phase structure"; both fail -> the A87 statement
  stands with the per-cell D caveat, the reconciliation stays open. Feasibility (scout): copy-number state for >= 40
  patients and >= 10 of 14 clock-gene loci.

## A135 (2026-10-24, written after WP0 item 8 showed the brain contrasts are anchor-dependent; PI request, Q6): which brain compartment, if any, should carry the host?

The disease contrasts under A17 anchors and all-LOO hosts have been seen (WP0 item 8). The choice of reference is
therefore made only from control donors and literature, with the rule fixed here, before any contrast is recomputed.
Datasets: the brain sets of WP0 item 8 (snRNA-seq, compartments as in A17). Code: `paper/exploratory/A135_brain_reference/`.
- **Criteria per compartment, control donors only:** (1) clock strength: Poisson-corrected clock-plane share against
  its permuted-acrophase null; (2) gate: share of other compartments with which the model-free cross gate passes;
  (3) centrality: mean pairwise coherence (Poisson-corrected |complex correlation|, as A109 C) with the other
  compartments; (4) stability: each of (1)-(3) per dataset.
- **Rule.** Rank compartments within each dataset by the sum of ranks of (1)-(3). A compartment family (e.g. excitatory
  neurons, astrocytes) is the reference if it ranks first in >= 2/3 of the datasets and passes the gate with >= 80 % of
  the others in each; ties -> the family with the larger median (3). Otherwise no single reference is justified and the
  brain panels are reported with their anchor dependence (Extended Data).
- **Literature (recorded, not used by the rule).** Neurons carry the brain host in the archive (02, 04, A18 3.5);
  astrocyte clocks can drive SCN timekeeping (Brancaccio et al. 2017, 2019; verify before citing).
- **Then:** recompute the brain disease contrasts with the chosen reference (single-family host, leave-one-out for its
  own members); report per contrast the sign agreement with the A17-anchored and all-LOO values. Descriptive; the
  disease verdicts of A17 are not re-tested.

**A134 note (2026-10-24, written after the scout, before scoring).** Scout: `paper/exploratory/A134_cancer_coherence/output/scout_*.csv`
(raw A87 coherence reproduced: malignant 0.228, normal 0.782; no observed gain computed).
- **Copy number:** no locus-level state is stored or annotated, so it is derived: per patient and clock-gene locus, the
  mean over malignant cells of the inferCNV-like profile of `cancer_common.cnv` (log1p CP10k minus the other patients'
  normal epithelium of the cohort, clipped at +-3, 101-gene running mean in genome order, own gene excluded, nearest
  expressed gene if the locus is unexpressed, <= 392 kb). Continuous score, stratum-centred, primary (no gain / loss
  cut). >= 50 malignant cells: all 79 patients, 14 of 14 loci. Expression-derived and relative: stated as inferred copy
  number. Normal epithelium: the same score of its own cells (companion, near neutral).
- **Estimator:** per gene, OLS on its own-locus score; A87 coherence recomputed on the residuals with the Poisson
  diagonal kept. Gain = adjusted minus raw. Statistic: malignant gain minus normal gain, patient bootstrap (1,000) for
  the interval; p from 1,000 permutations of whole-patient copy-number rows across patients; one-sided > 0. Null SD
  0.042, MDE 0.104 (the gap to explain is 0.55).
- **H134.2:** MSigDB Hallmark v2024.1 (MYC_TARGETS_V1, HYPOXIA without RORA; identical to v2026.1 for these sets);
  score = set minus expression-matched controls (`cells_cancer.scores`) on malignant cells from whole-transcriptome
  counts; cycling = Seurat S or G2M phase call (A24 definition) on malignant cells. The three scores enter one joint
  regression per gene; same statistic and nulls as H134.1. Not independent of H134.1 (shared neighbourhood genes).
- **D134.3:** the 64 patients with >= 100 cells in both groups (halves >= 50).

## A138 (2026-10-25, replaces WP5 / A114; written after A52, A83, A86, A87, A97, A134 and the GTEx method benchmark were scored, before any A138 readout; PI, Q9 e): what can be estimated about the clock in tumours?

Motivation: published untimed orderings of tumours (CoFE, Ananthasubramaniam et al., PLoS Biol 2025, 11 TCGA
adenocarcinomas; CYCLOPS / CYCLOPS 2.0 applications) may reflect the method rather than tumour clocks. CoFE's validity
signal is an interior minimum of its cross-validation curve; it reports that tumour clocks "do not have a fixed timing
relationship with clocks in healthy tissue", and its consensus tumour core clock places PER2, PER3, TEF and NR1D2 near
ARNTL (Fig 2D, lead's reading of the image), where healthy tissue has them 9-14 h apart (GTEx reference). CYCLOPS
requires ">250 samples" and "data from the entire periodic cycle" (Anafi 2017). Our GTEx benchmark (n = 80): CoFE
genome-wide per-gene acrophase error 3.10 h / 3.30 h (full day / 08-18 h) vs 1.53 / 2.30 h on the clock genes only;
CoFE flags its ordering as valid in all 11 TCGA cancers and in every GTEx run, including the 08-18 h windows. Sources and
quotes: `paper/exploratory/A138_cancer_estimability/{literature,inventory,A138_draft}.md`. A failed test is reported as
"not estimable at this design", never as "decoupled". Code: `paper/exploratory/A138_cancer_estimability/`.

**Q1. Does the malignant clock carry an estimable host-coupled signal? (single cell, gate-passing cohorts)**
- **H138.1:** pooled anchored coupling of the malignant compartment (copy-number-called malignant cells where available,
  else tumour epithelium; learned A17 anchors) > 0: lower bound of the 95 % patient-bootstrap interval (1,000) > 0,
  DerSimonian-Laird over cohorts as A17.
- **D138.1 (descriptive):** per cohort, the malignant clock-plane variance split into host-coupled (c^2 of the anchored
  fit), patient-specific (psi minus the Poisson expectation), between cells within patient (per-cell dispersion on the
  same variance scale) and Poisson; the same for normal epithelium and stroma of the same patients (specificity
  companion). Estimand of interest: patient-specific share, malignant minus normal epithelium.

**Q2. When do untimed orderings report a clock that is not there? (simulation, then TCGA placement)**
Simulator: genome-wide log-expression, p = 2,000 genes: the 15 CoFE clock genes at GTEx waveforms plus 15 % rhythmic
followers (GTEx gene phases), and non-clock structure from TCGA BRCA (782 tumours; the cancer of Li et al. 2024): its
leading non-clock principal components (purity, proliferation, subtype, copy-number blocks) with their variances. Arms
(50 replicates per cell; n in {100, 300, 500}):
- null: no clock (clock genes at mesor) plus the non-clock structure;
- window: intact clock, samples 08-18 h;
- full: intact clock over 24 h, patient phase SD in {0, 2, 4 h};
- coupled: each patient's tumour clock = own normal-tissue clock + common offset + patient SD in {0, 2, 4 h}, sampled
  08-18 h, tumour and normal ordered as in CoFE;
- tumour coefficients: tumour clock-gene amplitudes and phases differ from healthy (amplitude x 0.5, PER/TEF/NR1D2
  phases shifted 6 h), with whole-day healthy samples added as in Li et al. 2024 (CYCLOPS 2.0 with dataset and tumour
  covariates).
Methods: CoFE as published (genome-wide, CV, 5 restarts), CoFE on the clock genes only, CYCLOPS 2.0 with default seeds
(without healthy samples, and with them in the tumour-coefficients arm). Check before any scoring (simulator gate): in the
full arm without heterogeneity, CoFE and CYCLOPS recover the truth to within their published GTEx-like error (CoFE clock
genes median acrophase error <= 2 h at n = 300); if not, Q2 stops and the simulator is fixed first.
- **H138.2 (false validity):** null arm, n = 300: CoFE reports an interior CV minimum in > 20 % of replicates
  (one-sided exact binomial test against 0.20, p < 0.05: pass at >= 16 of 50; smallest attainable p 1e-35; power
  0.85 at a true rate of 0.38). Also reported: window arm; CYCLOPS 2.0's own fit statistic (fixed from its output
  in a dated note before scoring) in the same arms.
- **H138.3 (uncoupled inference):** coupled arm, patient SD 2 h, n = 300: the tumour-minus-normal predicted-time
  difference passes as uniform (Rayleigh p > 0.05) in > 50 % of replicates (one-sided exact binomial against 0.50: pass at >= 32 of 50; smallest
  attainable p 9e-16).
  Pass: CoFE's "no fixed timing relationship" does not follow from its readout.
- **H138.6 (genome vs clock genes):** full and window arms, n = 300, patient SD 2 h: median per-gene acrophase error of
  the clock genes larger with CoFE genome-wide than with CoFE on the clock genes only (paired over replicates, one-sided
  Wilcoxon p < 0.05).
- **H138.7 (healthy anchoring):** tumour-coefficients arm, n = 300: CYCLOPS 2.0 with healthy whole-day samples recovers
  the tumour clock-gene acrophases with median error > 3 h (one-sided, replicate bootstrap lower bound > 3 h), i.e. the
  healthy-coefficient assumption biases tumour readouts toward the healthy clock.
- **D138.2 / D138.3 (descriptive, TCGA):** where each of the 11 TCGA cancers falls on the simulation map (CV depth,
  phase coverage, rhythmic fraction); per cancer, agreement of the genome-wide and clock-only CoFE orderings; ARNTL-to-
  target acrophase differences (DBP, TEF, PER2, PER3, NR1D1, NR1D2) against GTEx, and the same differences in the
  simulation when the truth is a healthy clock (does genome-wide CoFE pull targets toward ARNTL?); share of clock genes
  shifted >= 8 h under the null and window arms.

**Q3. Does the malignant clock follow the patient's own normal tissue? (paired designs)**
- **H138.4:** TCGA paired tumour-normal (10 cancers, 517 patients): circular correlation of the tumour clock state with
  the same patient's normal state exceeds that with normals permuted across patients within cancer (10,000
  permutations, pooled over cancers, one-sided p < 0.05; smallest attainable p 1e-4). Gate: the cancer's normals pass
  the clock gate. Specificity companion on the same permutations: expression-matched pseudo-clock gene sets (A52 rule).
- **H138.5:** single cell, malignant (or tumour epithelium) vs adjacent-normal epithelium of the same patient,
  host-referenced states, cohorts with >= 6 paired patients; same permutation scheme; companion: tumour stroma vs
  adjacent stroma.
- **Power (before scoring, dated note):** null SD and MDE of the H138.1, H138.4 and H138.5 statistics at their designs.

Order: Q3 and Q1 (existing data and code), then Q2 after the simulator gate. Holm over H138.1-H138.7 within A138.
Reading: H138.2 / H138.3 / H138.6 / H138.7 pass -> the published tumour orderings and clock structures are not evidence
of tumour clocks; H138.1 / H138.4 / H138.5 then state what is estimable (or that nothing is, at these designs).

**A138 note (2026-10-25, before any A138 readout): the Li et al. 2024 design in the tumour-coefficients arm.** Checked in
the PNAS text (PMC10873596): healthy samples were not added to fill hours; the luminal A tumours (TCGA 193, Manchester 18)
were "projected ... onto the eigengene space computed from the noncancerous samples" (GTEx breast 167, TCGA normal 106,
26 timed Manchester normals), with tumour status as a covariate; tumour phases were not validated against times
(only the 26 timed normals, CorrFisher 0.69). The arm therefore simulates: healthy reference ordered with CYCLOPS 2.0
(GTEx-like coverage), tumours projected onto its eigengene space with a tumour covariate, tumour clock coefficients
differing from healthy as specified. H138.7 is unchanged (tumour acrophase error > 3 h). Other CYCLOPS tumour orderings
found: Anafi 2017 (HCC, projected on normal-liver eigenvectors after seeded ordering failed), Schwartz 2023 (PDAC, no
projection, coverage presumed 06-18 h), both from the same group (`literature.md` section 6).

**A138 note 2 (2026-10-25, before any A138 readout; PI discussion of the CYCLOPS projection).** (i) `clocksync` enters
every Q2 arm (fairness: the same simulated data, its anchored fit with the stroma-like compartments as anchors where the
arm has compartments, else reported as not applicable). (ii) The tumour-coefficients arm is split into a global change
(tumour clock amplitude x 0.5 and phase +3 h for all genes: clocksync's free z absorbs it, CYCLOPS projection does not
by construction) and gene-specific shifts (PER2, PER3, TEF, NR1D2 +6 h: violate clocksync's conserved waveform too;
clocksync's gate share and, if A139 passes, its gene refit are reported). (iii) The healthy reference is sampled with
GTEx-like coverage, tumours 08-18 h; CYCLOPS 2.0 is run with the tumour covariate's offset and scale free (as Li et al.)
and fixed to the healthy values. Rationale: with free tumour means, the mean of window-sampled tumours is not the
centre of the healthy circle, so tumour phases and magnitudes depend on tumour coverage; with fixed means, tumour and
healthy mesors are assumed equal. H138.7 is read for both variants.

**A138 note 3 (2026-10-27, new lead; before any A138 readout or power computation): operational statistics of H138.1,
H138.4, H138.5 and the power step.** Known at registration (stated, not new): the A52 per-cancer tumour-normal
correlation `corr` (archive "up to 0.52"); H138.4 adds the pairing null, not a new readout.
- **H138.4 statistic.** Per cancer, states of tumour and normal on one projector (`clocksync.paired.cohort_states`,
  counts with depth), centred within sample type; per-cancer statistic = Poisson-corrected |complex correlation|
  (`clocksync.paired.cross`, `corr`). Pooled = mean over the gated cancers (normals pass `reference_gate`, A52), weights
  n_pairs. Null: normals' rows permuted across patients within cancer (pairing broken, both marginals kept), 10,000
  permutations, same weights; p = (1 + #null >= obs) / 10,001. Specificity companion on the same permutations: 200
  expression-matched pseudo-clock gene sets (A52 rule, `pseudo_waveforms`); reported as the pseudo-set pooled statistic's
  95th percentile next to the clock value (descriptive, not a pass condition).
- **H138.5 statistic.** Units: patients with malignant (copy-number-called, else tumour epithelium) and adjacent-normal
  epithelium and a stromal host (`cancer_sc.paired_epithelium`, main tables, `cnv=False` except for the malignant call).
  Host = the patient's stromal compartments (`Dataset.host()` on the stromal anchors of the cohort). "Host-referenced
  state" = the residual of the compartment's state after its complex least-squares regression on the host over patients
  (removes the shared host component). Statistic: |complex correlation| of the malignant and normal-epithelium residuals
  (Poisson-corrected as `cross`), per cohort with >= 6 paired patients, pooled with weights n. Null: normal residuals
  permuted across patients within cohort (10,000). Companion on the same permutations: tumour vs adjacent fibroblasts,
  host = the remaining stromal compartments; descriptive.
- **H138.1** unchanged (anchored coupling of the malignant compartment, learned A17 anchors, 1,000 patient bootstraps,
  DerSimonian-Laird).
- **Power step (written to `paper/exploratory/A138_cancer_estimability/output/power/` and recorded before scoring).**
  H138.4/H138.5: the permutation null distribution on the real data (no observed statistic computed); reported: null
  mean, SD, 95th percentile, smallest attainable p (1/10,001), MDE = null 95th percentile + 0.84 x null SD (normal
  approximation; one-sided 0.05, power 0.8). H138.1: `clocksync.simulate` at the observed design of each Q1 cohort
  (patients, malignant and host depth), true coupling c in {0, 0.1, 0.2, 0.3, 0.5}, 100 replicates, 200 bootstraps;
  MDE = smallest c with pooled lower bound > 0 in >= 80 % of replicates; false-positive rate at c = 0 reported.
  If an MDE exceeds the largest grid value, the hypothesis is reported as "underpowered at this design" alongside its
  verdict.

**A138 note 4 (2026-10-27, after the H138.4/H138.5 null distributions and before any observed statistic or H138.1
simulation).** (i) H138.1 power simulation reduced for cost (registered design: about 31 h at 4 workers): 30 replicates
per c, 100 bootstraps per fit (about 5 h); MDE rule unchanged (>= 24 of 30). (ii) H138.5 data rule: a cohort whose
Poisson-corrected residual variance is <= 0 in either compartment, or whose null mean exceeds 1 (the correction has broken
down), is not estimable and leaves the pool; applied identically to the observed data. (iii) The fibroblast companion
uses its own permutations (its patients differ from the primary's), so it is not "on the same permutations"; it is
descriptive. (iv) The H138.4 pseudo-clock companion is computed only at scoring (it is an observed quantity).

**A138 note 5 (2026-10-27, after the power step and before any observed A138 statistic; PI decisions).**
(i) **H138.1 decision rule calibrated.** The power simulation shows that the registered rule (pooled lower bound > 0)
passes in 4 of 30 null replicates (false-positive rate 0.133; coupling is bounded at 0, so the pooled estimate is
biased upward: mean 0.030 at c = 0). New rule: H138.1 passes if the observed pooled coupling exceeds 0.083, the 95th
percentile of the simulated null pooled estimate (`output/power/h138_1_sim.csv`, c = 0). One-sided p = share of the
30 null replicates whose pooled estimate is >= the observed value, (1 + k) / 31 (smallest attainable 0.032). The
registered lower-bound rule is reported next to it. MDE about 0.2 under either rule.
(ii) **H138.5 is descriptive** (MDE 0.677 against a null of 0.443 ± 0.092): reported with its detection limit, out of
the A138 Holm family. Holm runs over the remaining confirmatory A138 hypotheses (H138.1-H138.4, H138.6, H138.7).

## A139 (2026-10-25, method; written before any simulation of it; PI request): can clock-gene acrophases be refitted under hospital hours?

Motivation: `clocksync` fixes the 14-gene GTEx waveform (relative amplitudes and acrophases) because with all loadings
free and partial-day sampling the loadings are identified only up to a 2 x 2 linear map (anisotropic, unknown host
covariance). Fixing the other genes pins that map, so one gene's loading is identified by its regression on the host
built without it. Never simulated under hospital hours (archive: A28 refitted with times of death; A31 window check:
follower genes R 0.91-0.96 with 3 consecutive ZTs, 0.15-0.33 with 2). Code: `paper/exploratory/A139_waveform_refit/`.
Simulator: `clocksync.simulate` (GTEx waveform, human depth), 4 anchor compartments (coupling 0.8-1.0) and one target
compartment, `gene_shift_h` for shifted genes; windows 24 h, 08-18 h, 10-16 h; n in {40, 80, 200}; 100 replicates per
cell. Scenarios: S0 no shift; S1 one target gene +6 h (PER2); S3 three target genes +6 h (PER2, TEF, NR1D2); SW the same
three shifted in every compartment (waveform itself wrong).
Estimators: (E-a) target-compartment gene refit: each clock gene of the target regressed on the anchored host
(`gene_rhythms` machinery), deviation = arg(c_g) - arg(z_target) - psi_g; (E-b) iterative leave-one-gene-out waveform
refit for SW: host fitted on 13 genes over all compartments, the left-out gene's loading regressed on it, all genes in
turn, repeated until the loadings change < 0.1 h (max 10 rounds). A gene is flagged when the 95 % donor-bootstrap
interval (200) of its deviation excludes 0.
- **V139.1 (refit, primary; 08-18 h, n = 80):** E-a in S1 and S3, E-b in SW: median absolute error of the refitted
  deviation <= 1.5 h for shifted genes and <= 1.0 h for unshifted genes. Pass -> the refit becomes a clocksync option
  (`refit_genes`), with the flag.
- **V139.2 (flag, fallback):** S0: per-gene false flag rate <= 0.05 (each gene, binomial upper bound); S1/S3/SW: power to
  flag a 6 h shift >= 0.8. V139.1 fails and V139.2 passes -> clocksync flags misfitting genes and refits the compartment
  without them (genes dropped, waveform otherwise fixed). Both fail -> the conserved waveform stays an untestable
  assumption, stated as such.
- **D139.3 (descriptive):** all cells of the grid (window x n), the gate's clock share under S3 / SW, and the bias of
  the compartment's z when shifted genes are not dropped.
Order: run before A138's Q2 (its clocksync arm uses the outcome). Compute: 14-gene fits, minutes per cell, run when the
M3 lanes free workers.

**A139 note (2026-10-25, before any A139 run).** Estimator E-a is clocksync's existing E4, `Dataset.gene_loadings`
(per gene, free loading on the host of the other compartments; `dev_h` against the conserved waveform; weighted
likelihood-ratio test of conserved vs free with a residual-permutation null). It was pre-registered in the archive
cancer plan (2026-10-03) and used in the cancer pipeline (`summary_*_e4_genes`, descriptive), the cycling follow-ups
(C3 half-life gradient; per-gene hi - lo timing, e.g. PER3 -2.1 h, no gene passes FDR) and A24-related checks. Its
only validation was a null calibration (0 of 40 false positives, "low power", archive ANALYSIS_PLAN); its accuracy and
power for known gene shifts under hospital hours were never simulated. A139 therefore also scores E4's likelihood-ratio
test as a flag (false-flag rate in S0, power in S1 / S3) next to the bootstrap flag; until A139 is scored, every E4
readout is unvalidated and is not quoted as a gene-specific finding.

## A140 (2026-10-26; written after reading the preprint's Fig 7 and Fig S8 and after A17 (cancer) was scored, before any A140 run; PI, 2026-10-26): does a CYCLOPS 2.0 donor internal phase across immune compartments report a clock that is not there?

Complements A138 Q2 (bulk tumours). Target: the multi-compartment design of the bioRxiv preprint "Circadian misalignment
underlies immune escape in breast cancer" (doi 10.64898/2026.05.26.726543; archive `literature/circadian_breast_cancer{,_supp}.pdf`):
TNBC scRNA-seq (Zhang et al., Cancer Cell 2025, GEO GSE266919, 76 samples), pseudobulks of CD4+ T, CD8+ T and myeloid
cells ordered jointly by CYCLOPS 2.0 with cell type as a covariate (free fit on the cohort: no healthy reference or
projection is described), seed = union of mouse immune "strongly correlated circadian genes" (Table S5, not available),
donor internal phase (DIP) = circular mean of the three phases after iterative removal of phases > pi/2 from the mean,
12 h mouse-to-human shift, six DIP bins, DiscoRhythm. Claims read as possible artefacts (S8B cross-type circular r
0.80-0.92; S8C DIPs over about three quarters of the cycle for surgical samples; S8D acrophases of thousands of rhythmic
genes in one sector; S8E/F clock genes in primate order) and the reasons are in `paper/exploratory/A140_dip_cyclops/A140_draft.md`.
Code: `paper/exploratory/A140_dip_cyclops/`. A failed clocksync gate is reported as not estimable, never as decoupled.

**Q1. Simulation of the design.** 76 patients x 3 compartments (CD4-like, CD8-like, myeloid-like) at single-cell
pseudobulk depth (`clocksync.simulate`); p = 2,000 genes: clock genes at GTEx waveforms, rhythmic followers (GTEx gene
phases), and a non-clock patient factor (K in {1, 3} shared axes loading on 10-30 % of genes, variance calibrated to the
leading non-clock principal components of our gate-passing tumour immune pseudobulks). Sampling 08-18 h unless stated.
Arms: A coupled clock (one patient clock shared by the three compartments, coupling 0.8) plus the factor; B no clock
(clock genes at mesor) plus the factor; C independent clocks (each compartment its own phase, uniform over patients)
plus the factor; D arm A sampled over 24 h (positive control). Seed grid: 50, 150, 500, 1,500 genes, the genes most
correlated with the factor and the clock in a simulated timed reference, with and without the core clock genes; 50
replicates per cell. Pipeline as published (CYCLOPS 2.0 with the cell-type covariate, DIP with the pi/2 rule, six bins,
cosinor on bin means, p < 0.05). Same data through `clocksync` (free host on the three compartments, gate).
Simulator gate (before scoring): in arm D at seed 1,500, CYCLOPS phases reach circular correlation >= 0.7 with the truth;
if not, Q1 stops and the pipeline is fixed first.
- **H140.1 (false concordance):** arm B, seed 500, K = 1: median cross-compartment circular r of the CYCLOPS phases
  >= 0.5 in more than 50 % of replicates (one-sided exact binomial against 0.5, p < 0.05: pass at >= 32 of 50).
- **H140.2 (false rhythms):** arm B, seed 500, K = 1: >= 500 rhythmic genes called per compartment in more than 50 % of
  replicates (same test).
- **H140.3 (clocksync specificity and sensitivity):** gate passes in <= 10 % of replicates in arm B and in arm C
  (one-sided exact binomial against 0.10 in each, p < 0.05: pass at <= 1 of 50), and in >= 80 % in arm A (exact binomial
  against 0.80 not rejected downward at p < 0.05: pass at >= 35 of 50).
- **Descriptive:** DIP coverage (share of the circle, Rayleigh R), acrophase concentration of the called genes (Rayleigh
  R) and clock-gene order (A138 `GTEX_MATCH`) per arm, over the whole seed grid and K = 3.
Holm over H140.1-H140.3. Reading: H140.1 and H140.2 pass -> the preprint's concordance and rhythm counts are expected
without a clock at its design; H140.3 states whether clocksync separates the arms.

**Q2. The preprint's cohort (GSE266919).** Data: processed matrices from GEO (PROVENANCE.md). Readouts, fixed in a dated
note after the data inventory (samples per patient, compartments, cell counts, available metadata) and before any clock
readout (PI, 2026-10-26): (i) gate and candidate anchor compartments (A17 rules); (ii) the published pipeline and
`clocksync` on the same pseudobulks; (iii) directional replication of the preprint's mouse findings on TME composition
(myeloid and CD8+ T fractions) and antigen presentation (myeloid) vs T-cell recognition along the host phase within the
daytime window (human daytime ~ mouse CT16-CT22 after the 12 h shift); with daytime sampling only, the expectation is a
small change; the readout is its direction.

**A140 Q2 readouts (dated note 2026-10-26, after the data inventory `paper/exploratory/A140_dip_cyclops/inventory_GSE266919.md`
and before any clock readout on GSE266919).** Mouse directions read from the preprint figures (no tables available):
from CT16 to CT22, immune fractions (S1H) CD4+ T about 48 -> 52 %, CD8+ T about 13 -> 17 %, Macro/Mono about 11 -> 6 %;
Macro/Mono antigen-presentation score (5F) high -> about neutral; T-cell interferon-response and effector scores (5C) rise.
After the 12 h shift, CT16 -> CT22 is human CT4 -> CT10 (late morning to late afternoon). The mouse fractions change
monotonically over CT4-CT22 (CD4+ T 26 -> 52 %), so the mouse "rhythm" may itself be a trend.
- **Data and units.** Compartments = the five deposited objects (CD4_T, CD8_T, myeloid, B, NK); no stromal, vascular or
  epithelial compartment exists. Unit = biopsy (sample). Pseudobulk = raw UMI sums per sample x compartment, kept at
  >= 20 cells. Primary set: the 44 Pre-treatment samples (one per patient), strata = biopsy site (breast vs other).
  Sensitivity: all 78 samples, strata = Group x site, patient as permutation block where the call supports it
  (descriptive).
- **(i) Gate and anchors.** `Dataset.gate(n_perm=1000)`, pass at p < 0.05. Anchors by the A17 learned rule
  (`host_identity/anchor_rule.py`, frozen learning tables, target study zhang2025, tissue breast). From the learning
  table alone (no GSE266919 data read): myeloid ranks first (pooled 0.202, lower 0.136) and lymphocyte qualifies (0.091,
  lower 0.027), so all five compartments are candidate anchors, subject to the target's own pair confirmation. The
  fixed rule gives no anchor. The host is therefore an immune host, not comparable with A17's stromal-host couplings.
  Reported: per-compartment leave-one-out coupling, `stability` (identified at >= 0.8) and `free_agreement`. Gate
  failed or host unstable: (iii) is not estimable.
  <!-- src: paper/results/extended_data/output/host_identity/a17_anchor_families.csv (target organ_breast, all tissues) -->
- **(ii) Published pipeline.** Run after the Q1 simulator gate passes, on the CD4_T, CD8_T and myeloid pseudobulks of the
  primary set. Table S5 is unavailable, so two seeds are used: the CYCLOPS 2.0 default gene selection (500 genes), and
  the same plus the 14 GTEx clock genes. Readouts 1-4 of the draft (cross-type circular r, DIP coverage, acrophase
  concentration, clock-gene order), plus the DIP and the clocksync host phase compared (circular r). Paired Pre - Post
  host phase within the 32 patients with both is descriptive (treatment confounds DIP bins).
- **(iii) Directional replication.** Host phase per sample from `Dataset.host()`. Its orientation (larger = later) is
  checked on a simulated timed cohort before scoring. The phase is unwrapped around the circular mean (within +/- 12 h).
  Per readout: OLS slope on host phase (per h) with one intercept per stratum; 95 % CI by 1,000 sample bootstraps with
  host refit; p from 1,000 within-stratum permutations of the host phase (the specificity companion: an arbitrary phase
  with the same spread). Predicted signs (mouse CT16 -> CT22):
  - **P140.1:** CD8_T fraction of all deposited immune cells, logit, slope > 0.
  - **P140.2:** myeloid fraction, logit, slope < 0.
  - **P140.3:** myeloid antigen-presentation score, slope < 0. The score is the mean z-scored log CPM of HLA-DRA,
    HLA-DRB1, HLA-DPA1, HLA-DPB1, HLA-DQA1, HLA-DQB1, CD74, CIITA, HLA-A, HLA-B, HLA-C, B2M, TAP1, TAP2, TAPBP, PSMB8 and
    PSMB9 (class II and class I sub-scores descriptive).
  Holm over P140.1-P140.3, one-sided in the predicted direction. "Same direction" = point estimate of the predicted
  sign; "replicated" = Holm-adjusted p < 0.05. Descriptive: CD4_T fraction (> 0), T-cell interferon-response and
  recognition genes (CD3D, CD3E, LCK, TRBC1), and the same readouts on the DIP of (ii) for comparison.

**A140 Q1 revision (2026-10-27, after the simulator gate failed and before any rerun; PI, Q11 (b)).** The registered gate
failed (CYCLOPS truth r median 0.067, 0 of 50; clocksync gate 43 of 50; `results_preliminary.md`). Two run choices may
cause it, so the simulator is revised before anything is rerun. Everything else in Q1 is unchanged: arms, H140.1-3,
thresholds and Holm.
- **CYCLOPS:** the template's 80 models per fit (the failed run used 10).
- **Factor calibration:** genome-wide pseudobulks (all genes, raw UMI sums, >= 50 cells) of the tumour T_NK and tumour
  myeloid compartments of the 21 gate-passing A17 cancer cohorts (`cancer_sc`, main tables, `cnv=False`). Per
  compartment and cohort, log CPM is centred within cohort. The non-clock factor is the leading principal components
  after removing the projection on the 14 clock genes' host plane. Their eigenvalues are corrected for Poisson noise
  (expected sampling variance of log CPM subtracted per gene before the PCA). The median over cohorts sets the factor
  sd at scale 1.
- **Factor size as a grid:** scales 0, 0.25, 0.5 and 1 of the calibrated sd. Scale 1 is the primary.
  - Scale 0 runs only if the CYCLOPS port completes on it; a crash is reported as "not run".
  - The simulator gate is evaluated per scale (arm D, seed 1,500, 50 replicates, median circular r >= 0.7).
  - H140.1-3 are scored at scale 1 only if the gate passes there.
  - If it fails at scale 1 but passes at a smaller scale, H140.1-3 are not scored. Reported instead: the largest scale
    at which CYCLOPS recovers the clock, against the calibrated tumour-immune factor size (descriptive).
- **Cost check before the full run:** 2 replicates at 80 models. If one gate cell (50 replicates) would exceed 6 h on
  16 workers, replicates drop to 30 per cell for the gate only. H140.1-3 keep 50 (thresholds as registered).

**A140 Q1 revision note (2026-10-27 13:20, POST HOC: written after the scale-1 gate was seen; PI decision).** Scale 1
fails (50 of 50 replicates: median truth r 0.060, 0 of 50 >= 0.7), so H140.1-3 are not scored, as registered. The
remaining queue (scales 0.5, 0.25, 0 at 50 replicates, about 16 h) is stopped. Scale 0 (no factor) runs next with 20
replicates: it separates a factor-driven failure from a failure of the CYCLOPS pipeline at this pseudobulk design. Scale
0.25 runs (20 replicates) only if scale 0 passes. The 13 finished scale-0.5 replicates are reported as descriptive.
Gate rule per scale unchanged (median truth r >= 0.7).

## A141 (2026-10-27, written after WP0 item 8 and A135, before any readout below; PI, Q6): are the brain disease contrasts stable enough to read, given that brain donors die at any hour?

Brain donors (PsychAD MSSM / HBCC, SEA-AD, PD atlas, RADC) die at any hour. A host fitted on them should then cover the
full cycle, so compartment amplitudes and phases are close to absolute. The question is whether the anchor dependence of
the brain contrasts (WP0 item 8: 1-7 of 9 compartments change sign between the A17 anchors and the all-LOO host; r*
0.42-0.67) is what noise alone produces at these amplitudes and donor numbers. Descriptive; no disease hypothesis. Data:
the brain contrast sets of WP0 item 8. Code: `paper/exploratory/A141_brain_stability/`. A failed gate is uninformative.
- **Step 1, coverage.** Per cohort (control donors, then all donors): Rayleigh R of the host phases and the host spread.
  Reference: the same statistics in GTEx fast-death brain tissues (timed, any hour) and in simulations with uniform
  24 h sampling. "Full cycle" = R within the 95 % range of the uniform-sampling simulation at the cohort's n.
- **Step 2, absolute amplitude.** Per compartment: host-referenced log amplitude plus the host amplitude, compared with
  the timed GTEx brain-region amplitudes (`organs/tissues.csv`, fast cohort, brain tissues). Same scale = within the
  range of the GTEx brain regions. Reported, not tested.
- **Step 3, stability by simulation (decision).** Per cohort, 200 simulations (`clocksync.simulate`) at the fitted
  amplitudes, phases and own noise of every compartment, the cohort's donor numbers per group, uniform 24 h sampling, and
  no case effect. Each is refitted with the A17 anchors and with the all-LOO host, at the contrast's design.
  - Statistic: the number of compartments whose case-minus-control d changes sign between the two hosts.
  - Noise-driven: the observed count is <= the 95th percentile of the simulated counts. Then nothing can be said about
    that contrast's compartments (reported as not estimable at this design). Otherwise: not noise-driven.
- **Step 4, descriptive map (only for contrasts not noise-driven).** Which compartments agree in sign and size across
  the two hosts, read against expectations fixed now:
  - expected to agree: the neuron classes among themselves, and astrocytes with oligodendrocytes;
  - expected noise-driven: compartments with the lowest coupling or < 20 donors per group.
  If some contrast is not noise-driven, the lead drafts a request for time-of-death fields to the PI (outward-facing:
  PI decides).

## A123b (2026-10-27, written before any run of the added methods on the A123 design; PI request): every benchmark method on the full working-hours organ map

Motivation (PI, 2026-10-27): Figure 1 keeps one organ-map panel, A123 (40 tissues, 111 working-hours donors), and drops
the 8-tissue benchmark row; the method comparison moves onto the A123 design. **Known at registration:** A123 (clocksync
offsets vs truth (i) 0.945; CHIRAL -0.165 on the same donors) and the 8-tissue benchmark at n = 80 (`gtex_benchmark/
summary.csv`: in the window design, offset r clocksync 0.900, CHIRAL 0.350, CoFE clock genes 0.289, CoFE genome-wide
0.099, CYCLOPS -0.495, TimeTeller 0.644, timed oracle 0.165; failure rates up to 0.62). The prediction below is
therefore informed by the benchmark; what is new is the 40-tissue design. Code: `paper/exploratory/A123_working_hours/`.

- **Designs.** (W) the A123 design: GTEx v10 fast deaths (Hardy 1-2) who died 08-18 h, the 40 tissues with >= 20 such
  donors, all window samples (no draws). (F, descriptive) the same 40 tissues with all fast deaths (whole day).
- **Methods, run as in `results/figure_1/scripts/method/compute_gtex_benchmark.py` (`run_order`, default settings):**
  - CoFE genome-wide (`cofe_filter`, 5 restarts) and CoFE on the clock genes;
  - CYCLOPS 2.0 with its human seed genes and template hyperparameters;
  - CHIRAL (as A123);
  - TimeTeller trained on recorded times, 5-fold out-of-fold over all fast deaths (a supervised reference, not untimed);
  - the timed oracle (per-tissue cosinor on the recorded times of the design's own donors, rank-1 complex SVD);
  - clocksync as A123. In (W): waveform cross-fitted on fast deaths outside the window. In (F) no donor is left out, so
    the GTEx waveform is used (A123 sensitivity: unchanged results); stated as a mild advantage.
- **From orderings to the organ map (as A123 did for CHIRAL).**
  - Each method orders each tissue's samples separately.
  - Donor times are a consensus: every tissue is aligned (rotation, reflection) to the reference tissue on shared donors,
    then averaged on the circle. The reference is Muscle - Skeletal, or the tissue with the most donors whose fit
    succeeded if Muscle fails.
  - Per tissue: clock-gene cosinor on the consensus times, rank-1 complex SVD, giving offsets and log2 relative
    amplitudes.
  - A failed fit (error, or all phases missing) drops that tissue for that method; the failure rate is reported.
- **Readouts against truth (i) and (ii), the A123 `compare` (10,000 tissue permutations, 200 tissue bootstraps):**
  - offset circular correlation and median absolute error;
  - log2 amplitude Pearson;
  - for the orderings, donor order: |circular correlation| of the consensus times with recorded time, best rotation and
    reflection (favourable to them), and median absolute error.
- **H123.4 (secondary, design W, truth (i)):** clocksync's offset circular correlation exceeds that of each untimed
  ordering method other than CHIRAL (CHIRAL's value is already known): CoFE genome-wide, CoFE clock genes, CYCLOPS 2.0.
  - Statistic: difference in circular correlation, paired tissue bootstrap (2,000 resamples of the tissues both methods
    estimate).
  - One-sided p is the share of resamples with difference <= 0 (smallest attainable p 1/2001). Holm over the three.
  - A method with estimates in < 20 tissues is "not estimable" for this comparison.
  - Pass per method at Holm-adjusted p < 0.05.
- **Descriptive:** the same readouts in design F; TimeTeller and the oracle as references; per-method runtime and failure
  rate.
- **Compute.** CoFE genome-wide and CYCLOPS take about 30 min per tissue fit at n = 80 (4 threads). (W) is about 80 such
  fits, roughly 11 h on 12 workers, and runs first. (F) has up to 244 donors per tissue, about 2-3x longer, and runs
  after (W). Detached `run.sh` job, at most 16 workers and the RAM budget.
- **Figure consequence:** Figure 1 organ-map panel = A123 with all methods (design W main, design F in ED or as a second
  column).

## A142 (2026-10-27, written before any comparator below was computed; PI request): is the brain-periphery split specific to the brain, and do digestive organs carry a departure of their own?

Motivation (PI, 2026-10-27): r* < 1 for the brain says little if any group of organs gives r* < 1 against the rest.
The PI also expects the digestive organs to depart from the body's internal time, because the time of the last meal
resets their clocks and is partly decoupled from internal time. **Known at registration:**
- brain vs periphery r* 0.64 [0.56, 0.73] (A103 base);
- the archive factor model of the departures (01 GTEx organs): outside the brain, 2 factors fit hardly better than 1
  (held-out 0.349 vs 0.337), and all peripheral tissues, including liver and gut, load on one factor. That model was
  not a test of a digestive block.

None of the r* values below have been computed. Code: `paper/exploratory/A142_organ_blocks/`.

- **Data and statistic, as A103 base** (`results/figure_2/scripts/organs/compute_brain_artefact.py`):
  - GTEx v10 fast deaths (Hardy 1-2), tissues with >= 40 donors: 13 brain regions, 27 peripheral tissues, 244 donors.
  - Per-sample linear phase departures from the timed fit of the 14-gene plane, standardised per tissue.
  - A block score is the donor's mean departure over the block's tissues (>= 2 tissues).
  - r = correlation over donors of block A vs block B scores. r* = r / sqrt(rel_A rel_B), with split-half
    (Spearman-Brown) reliabilities over 50 random tissue splits.
  - The statistic is the same function generalised from "Brain*" to any tissue set.
- **Blocks:**
  - brain = the 13 brain regions, against the 27 peripheral tissues (A103 base);
  - digestive = Stomach, Colon - Transverse, Colon - Sigmoid, Esophagus - Gastroesophageal Junction, Esophagus -
    Muscularis, Esophagus - Mucosa, Liver (7; pancreas and small intestine have < 40 donors), against the other 20
    peripheral tissues.
- **Reference distribution:** 1,000 random blocks of k peripheral tissues, each against the remaining 27 - k peripheral
  tissues, with k = 13 for the brain comparison and k = 7 for the digestive one (seed 142). This is what an arbitrary
  group of peripheral organs gives.
- **H142.1 (brain specificity):** brain r* is below the k = 13 reference.
  - One-sided p = (1 + #{reference r* <= brain r*}) / 1,001 (smallest attainable 0.001).
- **H142.2 (digestive departure, the PI's hypothesis):** digestive r* is below the k = 7 reference, with the same
  p definition.
  - Also reported: 200-replicate donor-bootstrap interval of the digestive r*.
- **Decision rule:** Holm over H142.1-2, alpha 0.05. Power is stated as the reference's 5th percentile (the largest r*
  that could pass) and its SD.
- **Reading rules:**
  - H142.1 pass: the brain separates more than an arbitrary organ group, so "partly separate" is brain-specific
    (at tissue level, in GTEx).
  - H142.1 fail: the brain's r* is what any organ group gives, and the Figure 2 brain claim is dropped or restated.
  - H142.2 pass: digestive organs carry a shared departure of their own (consistent with meal timing; not shown to be
    meal timing). H142.2 fail: no digestive departure beyond what any organ group shows at this design.
- **Sensitivities (descriptive), on both blocks and on their references:**
  - the plane without PER1/PER2;
  - departures adjusted for ischaemic time and Hardy class (the A103 `adjust`).
- **Descriptive blocks**, each against the rest of the periphery with its own size-matched reference:
  - food-exposed epithelium and liver (Stomach, Colon - Transverse, Esophagus - Mucosa, Liver);
  - cardiovascular (Artery - Aorta, Artery - Coronary, Artery - Tibial, Heart - Atrial Appendage, Heart - Left
    Ventricle);
  - endocrine and adipose (Adipose - Subcutaneous, Adipose - Visceral, Adrenal Gland, Thyroid, Pituitary).
- **Figure consequence:** the Figure 2 brain-periphery row cites H142.1. A digestive departure, if found, becomes a
  Figure 2 row (main if it passes, else not shown).

## A143 (2026-10-27, written before any readout below; PI request): is the brain-periphery split a processing artefact, and does the hypothalamus follow the body?

Motivation (PI, 2026-10-27): A103 excluded collection delay, death type and acute PER1/PER2 induction. It did not
exclude a per-donor technical factor shared by the brain samples, which are collected and processed together, apart
from the other tissues (for example RNA quality, autolysis, nucleic-acid isolation batch). Baboon does not show the split
(r* 0.97, 12 animals; raw brain-periphery correlation 0.46); sheep does (0.63).

**Known at registration:**
- A103 (r* 0.64 base; 0.64 adjusted for ischaemic time and Hardy);
- A110 H110.4 (slow vs fast deaths: the hypothalamus keeps its coupling to the peripheral host while the other brain
  regions lose it, difference 0.208, p 0.002, descriptive);
- the batch structure, counts only: per donor a median 8 brain samples in 4 isolation batches (SMNABTCH) and 20
  peripheral samples in 7;
- the join coverage: 85 % of fast-death samples have v8 sample attributes.

No departure has been read against batch or quality. Code: `paper/exploratory/A143_brain_artefact/`. Same cohort,
departures and r* as A103 base and A142 (14-gene plane, fast deaths, tissues with >= 40 donors). Sample attributes from
`paper/data/Gtex_phenotypes` (controlled access: aggregates only). Samples without attributes are dropped from H143.1-2;
coverage is reported; < 50 % coverage of brain samples makes H143.1-2 not estimable.

- **H143.1 (quality adjustment).**
  - Per tissue, departures regressed on RIN (SMRIN), autolysis score (SMATSSCR, as a factor), ischaemic time, Hardy 2 vs
    1 and collection centre (SMCENTER) where it varies. Residuals plus intercept, as A103 `adjust`.
  - r* recomputed with 200 donor bootstraps.
  - Pass (split not explained by sample quality): adjusted r* < 1 with bootstrap upper bound < 1.
  - Also reported: the change from base r* on the same samples.
- **H143.2 (isolation batch carries a departure).** For each pair of brain regions:
  - the correlation over donors of the two standardised departures, separately in donors whose two samples share an
    isolation batch and in donors where they do not;
  - statistic: same-batch minus different-batch correlation, averaged over region pairs with >= 10 donors in each class,
    weighted by the harmonic mean of the two counts;
  - null: within each region pair, same / different labels permuted across donors (2,000 permutations);
  - one-sided p, smallest attainable 1/2001.
  - A pass is evidence for a batch-borne departure, so an artefact. The same is computed for peripheral pairs as a
    reference (descriptive).
  - Power: null SD and 95th percentile.
- **H143.3 (the hypothalamus follows the body more than the other brain regions; fast deaths).**
  - Per brain region k: alignment ratio = [r(k, P) / sqrt(rel_P)] / [r(k, B_-k) / sqrt(rel_B-k)]. Here P is the
    peripheral score, B_-k the score of the other brain regions, and rel their split-half reliabilities. The region's own
    noise cancels in the ratio.
  - Statistic: hypothalamus ratio minus the median ratio of the other 12 regions, with 1,000 donor bootstraps.
  - Pass: one-sided bootstrap p (share of replicates <= 0) < 0.05.
  - Reading: a technical factor shared by all brain samples predicts equal ratios; a pass means at least one region
    departs differently, as biology would allow.
  - Descriptive: every region's ratio, and the same in Hardy 3-4 deaths.
- **Decision rule:** Holm over H143.1-3.
- **Reading rules:**
  - H143.1 and H143.3 pass, H143.2 fails: the split is not explained by quality or batch; Figure 2 keeps "partly
    separate", with the hypothalamus panel if H143.3 passes.
  - H143.2 passes, or H143.1 fails: the brain claim is restated as "not separable from brain processing" and moves to
    Extended Data.

**A143 note 1 (2026-10-27, after a 150-donor check run, before the full run; H143.1 gave no estimate in the check).**
Data rule:
- The autolysis score (SMATSSCR) is present for 15 % of fast-death brain samples (84 % peripheral), against 80 % / 89 %
  for RIN, centre and isolation batch. Brain tissues are mostly not scored.
- H143.1 therefore adjusts for RIN, ischaemic time, Hardy 2 vs 1 and centre, the same model in every tissue.
- The registered model with autolysis is reported as a descriptive sensitivity in the tissues where >= 70 % of samples
  have a score. The brain is then represented only by regions that meet that threshold, if any.

Other parts are unchanged. A code fix that does not change any statistic: verdicts of H143.2-3 are set after Holm.

## A144 (2026-10-27, written before any readout below; PI request): why are strong clocks the coupled ones?

Motivation (PI, 2026-10-27): oscillator theory predicts that low-amplitude clocks entrain more easily, but A104 found
the opposite (GTEx periphery rho(A, c) +0.734 over 27 tissues; brain -0.115). Three explanations:
1. **Algebraic and noise sharing:** c = s / (s + u) and A = sqrt(s + u), with s = z2 V the locked and u = 2 psi the
   unlocked variance, estimated from the same donors.
2. **Tissue amplitude is a population amplitude:** synchrony of cells within the sample, which body-wide signals
   maintain, raises both amplitude and coupling.
3. **Cue exposure:** the tissue's input, rather than its oscillators' rigidity.

**Known at registration:** A104 (above); A18 part 2 per-cell readouts (`figure_2/output/backbone/a18_cell_split.csv`,
not read against coupling); the eligible-row counts of that table. Code: `paper/exploratory/A144_strong_clocks/`.

**Part 1, GTEx (explanation 1).** A104 design: fast deaths, the 27 peripheral tissues with >= 40 donors, all tissues
anchors, leave-one-out host, `Dataset.fit`.
- 100 random splits of the donors into halves H1 and H2 (seed 144). In each half, per tissue: c, A, s and u.
- Cross-half rho: mean of Spearman(A_H1, c_H2) and Spearman(A_H2, c_H1). Same-half rho: mean of Spearman(A_H1, c_H1)
  and Spearman(A_H2, c_H2). Each is averaged over the splits.
- **H144.1:** cross-half rho > 0.
  - Null: tissue labels of A permuted against c within each split, the same permutation in all splits (10,000).
  - One-sided p, smallest attainable 1e-4.
  - Pass: the amplitude-coupling association is not noise shared within the same donors.
- **D144.1 (descriptive):**
  - same-half minus cross-half rho (the noise-shared part), with its spread over splits;
  - the SD over tissues of log s and of log u, full sample: which part varies more;
  - Spearman of log s and of log u with c.

  Reading: if log s varies more than log u, strong clocks are strong by their shared part (explanations 2-3). The
  algebraic share is stated, not removed.

**Part 2, single cell (explanations 2 vs theory), descriptive (D144.2).**
- **Why only descriptive:**
  - The A18 per-cell amplitude lambda and within-sample synchrony R satisfy lambda R = rho, with rho^2 = s / (1 - R_W^2),
    built from the same fit as c. A regression of c on log lambda and log R therefore carries an equal mechanical
    coefficient on both through log rho.
  - The split between them is set by the within-sample dispersion D, which is computed without the host.
  - The informative estimand is the coefficient difference beta_R - beta_lambda in logit(c) ~ log lambda + log R +
    target. Population view: > 0 (synchrony, not per-cell amplitude, carries coupling). Theory, at fixed synchrony:
    beta_lambda < 0, so also > 0. The two differ in the sign of beta_lambda - beta_rho, where rho is the
    mechanical-only reference. Both coefficients are reported.
  - Eligible rows are few, mostly brain: about 40 healthy or control rows.
- **Data:** `a18_cell_split.csv` rows with `eligible` True in control groups (and `lung_hlca` "all").
  - MSSM controls are counted once (AD_mssm control; age_mssm dropped).
  - c is the coupling of the same group's anchored fit (A18's host refitted per group, frozen learned anchors),
    recomputed with the A18 loaders.
- **Statistic:** OLS with target fixed effects; 1,000 bootstraps of compartments within target. Reported: beta_lambda,
  beta_R, beta_R - beta_lambda with 95 % intervals. No pass / fail.

**Reading:**
- Part 1 pass with log s dominating: the strong-clock association is real and lies in the shared part. That is
  consistent with population synchrony or cue input, not with single-oscillator rigidity, which these data do not test
  directly.
- Part 1 fail: the A104 association is noise shared within donors, and the Figure 2 row is dropped.

## A115 (2026-10-27, written before any readout on healthy skin; PI decision): trait component of the psoriatic delay, descriptive

**Status: descriptive (PI, 2026-10-27).** H115.1 leaves the confirmatory family: the design cannot detect a plausible
trait effect (power below). Reported with its detection limit, no pass / fail wording beyond the interval.

**Known at registration** (blind preparation, `paper/exploratory/A115_A116_psoriasis/output/`; no clock statistic was
computed on healthy samples):
- Series with healthy controls from other persons (raw files on disk; the processed matrices hold lesional and
  uninvolved rows only): psoriasis GSE13355 (58 pairs, 64 healthy), GSE121212 (27, 38), GSE201827 (65, 17),
  GSE268958 (35, 19); atopic dermatitis GSE121212 (21, 38). GSE66511 (12, 12) failed the A54 gate and is excluded.
  `a115_inventory.csv`.
- Estimand check (15 gated A54 series, 594 pairs): the paired A78/A81b estimator gives lesional minus uninvolved
  -1.703 h (archive -1.689); the within-group estimator used here gives -1.045 h (14 of 15 series negative; coherence
  0.22 lesional, 0.41 uninvolved). The A78 systemic set is unusable within group (group offsets -6.5 and -3.2 h).
  `estimand_check_pooled.csv`.
- Null from random splits of uninvolved samples, pooled over the 4 psoriasis series: SD 0.882 h at full group sizes
  (with replacement), 1.163 h for the proportional split; MDE (one-sided 0.05, 80 %) 2.19-2.89 h. Power at a true
  -1.0 h: 0.22-0.31. Smallest attainable p 1e-4 (10,000 permutations). Atopic dermatitis: SD 4.5-6.5 h.
  `a115_null_pooled.csv`, `power_table.csv`.

**Data rule.** Per series, the healthy columns are added from the raw file on disk (`a115_inventory.csv`,
`healthy_source`), processed exactly as that series' builder processes its lesional and uninvolved rows (same platform,
gene mapping, scale and depth). Acceptance check before reading healthy rows: the rebuilt uninvolved rows reproduce
the processed ones (per-gene correlation > 0.99, or exact equality), else that series is dropped and reported.
Uninvolved = patients with a lesional-uninvolved pair, one sample per patient (as A54); healthy = one sample per person.

**Estimand.** Per series, one projector for the uninvolved and healthy rows together (`_lib.states_for`, GTEx core
waveform; GTEx cue set A81b restricted to the genes present). Group offset = arg E[c_core c_cue*] over the group's
samples (`S.ccorr`). Pooled offset per group = arg of the A54-weighted sum of the series' complex group statistics;
D115.1 = uninvolved minus healthy, in hours (wrapped).

**Reported (D115.1, primary descriptive):** D115.1 with a 95 % bootstrap interval (samples resampled within group and
series, 2,000), and the one-sided group-label permutation p within series (10,000), stated with MDE 2.19 h.
Reading rule: a negative estimate whose interval excludes 0 is "a trait component of at least ...", anything else is
"not detected; effects below about 2.2 h not excluded". Per-series values are shown.
**Descriptive:** the same for atopic dermatitis (GSE121212, one series); the keratinocyte differentiation score (A61
sets) uninvolved minus healthy, per series and pooled (standardised mean difference); coherence per group.
No systemic-set sensitivity (unusable within group, above).

## A116 (2026-10-27; PI decision): not registered as a new test

The only GEO series meeting the placebo rule (>= 10 placebo patients with lesional and uninvolved biopsies at both
visits) are GSE117468 (22 placebo, 83 treated) and GSE201827 (22, 33): the data on which A56 already reported the
treated-minus-placebo contrast post hoc (-0.50 h). A test registered now on these data adds no evidence. Arm-label null
SD 0.864 h, MDE 2.15 h; PASI75 within treated (GSE117468, 72 vs 11) MDE 5.64 h (`a116_null.csv`). Near misses
(GSE136757 9, GSE252029 8, GSE283313 8 placebo patients) have uninvolved biopsies at baseline only.
**Decision:** H116.1 and H116.2 are withdrawn; H116.1 leaves the confirmatory family. The A56 interaction stays post
hoc and is reported with the 2.15 h detection limit; the response analysis (A57) stays descriptive.

## A145 (2026-10-27, written before any readout; PI request): what is specific to clocksync against two-step orderers?

**Question (PI, 2026-10-27).** CHIRAL fits a shared donor phase assuming every tissue follows it; CYCLOPS and CoFE
order samples. Any orderer can be used in two steps: (1) order donors on the anchor compartments, (2) regress each
target compartment on the inferred phase. That two-step pipeline yields offsets, relative amplitudes and a coupling
proxy. Which clocksync readouts does it reproduce, and where does clocksync differ? No clocksync advantage is assumed;
the method description will state only what this amendment supports.

**Known at registration:** A123 (working hours: clocksync organ offsets r 0.945, CHIRAL -0.165); A140 Q1 (CYCLOPS
orders by a shared non-clock factor; cannot start without it at 1,500 seed genes); A123b running (phase recovery
only, no coupling). No two-step coupling has been computed anywhere in the project.

**Simulator.** `clocksync.simulate`, GTEx 14-gene waveform, n = 100 donors, depth and tau at the defaults, 4 anchor
compartments (coupled, psi chosen for c = 0.9) and 5 targets with true coupling c in {0, 0.25, 0.5, 0.75, 1.0} (psi
set from c; c = 0 is `coupled=False`), true offsets spread over 0-12 h, relative amplitudes in {0.5, 1}. 30 replicates
per arm (seed 145).
- **Arm F (full day):** window 0-24 h, shared_sd 0.
- **Arm W (working hours):** window 8-18 h, jitter 1 h, shared_sd 0.
- **Arm S (shared non-clock factor):** full day, shared_sd in {0.5, 1.0} (the A140 sizes on the log scale, matched by
  the share of clock-plane variance it takes; the value is fixed by a check on simulated data before any method runs).
- **Arm N (no clock):** amplitude 0, shared_sd 1.0.

**Methods.**
- **clocksync:** `Dataset.fit`, the 4 anchors, leave-one-out host; gate reported; coupling c, offset, relative
  amplitude per target.
- **Two-step, for CHIRAL, CYCLOPS 2.0 (clock-gene seed) and CoFE (clock genes):** step 1, the orderer on the donor x
  (gene, anchor) matrix of log2 normalised counts (concatenated anchors, as a multi-tissue CHIRAL input); step 2, per
  target, the same clock-plane state as clocksync (same projector, depth and Poisson variance) regressed on (cos, sin)
  of the inferred phase over donors: coupling proxy = explained variance / (total minus mean Poisson variance), clipped
  to [0, 1]; offset = arg of the complex coefficient; relative amplitude = its modulus over the anchors' mean.
- **Timed oracle:** step 2 on the true times (upper reference).
A failed orderer fit is data (counted as failure, not dropped).

**Readouts per method and arm:** Spearman(estimated c, true c) over targets, pooled over replicates; mean absolute
coupling error; circular correlation of offsets with truth; in arm N, the share of targets with c > 0.2 (false
coupling) and, for clocksync, the gate pass rate.

**Hypotheses (secondary family A145, Holm over the three; each against the best two-step method of that arm):**
- **H145.1 (working hours):** clocksync mean absolute coupling error < that of the best two-step method. Paired over
  replicates, one-sided sign test, p < 0.05 (smallest attainable 9e-10 at 30 replicates).
- **H145.2 (shared non-clock factor, shared_sd 1.0):** the same comparison as H145.1 in arm S.
- **H145.3 (no clock):** clocksync's false-claim rate, counting a target as claimed only if the gate passes and
  c > 0.2, is lower than the best two-step method's false-coupling rate. One-sided exact McNemar-type test over
  targets x replicates, p < 0.05.
**Descriptive D145.1 (arm F):** all readouts. Expected: two-step methods reproduce coupling, offsets and amplitudes on
full-day data. If they do, the paper states that coupling is recoverable by a two-step orderer pipeline when sampling
covers the day, and clocksync's specificity is limited to what H145.1-3 support.
**Descriptive D145.2 (real known answer):** A129 hepatocyte Bmal1 knockout vs control (timed mouse, full day): liver
coupling to the host of the other tissues, clocksync vs the three two-step pipelines vs the timed oracle; does each
see the drop (0.967 -> 0.676 in clocksync)?

**Reading.** Each passed H145 is a stated property of the method (partial coverage; robustness to shared non-clock
structure; an explicit test of clock presence). A failed one is removed from the method's claims. D145.1 sets how
"coupling" is presented: as a readout clocksync makes explicit, not one only it can compute, unless the data say
otherwise. Code: `paper/exploratory/A145_two_step/`.

## A146 (2026-10-27, written before any readout; PI request): does the tumour microenvironment follow the time-of-day programs of the mouse TNBC atlas in human breast cancer?

**Source of the predictions.** Preprint "Circadian misalignment underlies immune escape in breast cancer" (bioRxiv
10.64898/2026.05.26.726543; archive `literature/circadian_breast_cancer.pdf`), 4T1 mouse tumours at CT4/10/16/22. Its
human arm (CYCLOPS on immune pseudobulks of TNBC) is questioned in A140 and not used. Mouse acrophases are taken from
its main text and mapped to human wall time by activity onset: mouse CT12 = human 07:00, so human h = mouse CT - 5 h
(the preprint uses a 12 h shift in CT; this is the same mapping expressed in wall time).

**Family:** secondary, own Holm over the A146 hypotheses (PI, 2026-10-27). Part 1 below; later parts (TNBC single
cell, pan-cancer) are registered separately before their data are read.

**Known at registration.** GSE233242 (Li et al. 2024): 42 patients (primary join; builder `breast_timed_gse233242`),
tumour and adjacent-normal bulk TPM, recorded resection times 08-18 h; V25.5 per-patient internal time `est_day`
(`figure_1/output/windows/v25_5_patients.csv`; Spearman with recorded time 0.56, MAE 1.12 h). Archive D8: a
glucocorticoid score was read against these times (not one of the programs below). No program below has been
computed on these data.

### Part 1: GSE233242 (bulk, recorded and internal time)

**Programs** (tumour libraries; score = mean z of log2(TPM + 1) over the genes, z across patients; human orthologs of
the genes named in the preprint), with mouse acrophase -> predicted human acrophase:
| program | genes | mouse CT | human h |
|---|---|---|---|
| P1 proliferation (cancer epithelium) | MKI67 CCNB1 CENPA TOP2A AURKB PCNA MCM2 BUB1B BRCA2 MRE11 | 1.5 | 20.5 |
| P2 MHC-I antigen presentation (cancer epithelium) | B2M PSMB8 TAPBP TAP1 HLA-A HLA-B HLA-C | 22 | 17 |
| P3 interferon response | STAT1 IRF1 IFNGR1 IFNGR2 IFNAR2 STAT3 | 22 | 17 |
| P4 T-cell activation (residual on T abundance) | TCF7 LEF1 RUNX3 BCL6 | 18 | 13 |
| P5 myeloid antigen presentation (residual on myeloid abundance) | CIITA HLA-DRA HLA-DRB1 CD74 | 18 | 13 |
| P6 phagocytosis (residual on myeloid abundance) | DNM2 CDC42 SYK FCGR1A FCGR2A | 18 | 13 |
| P7 T-cell abundance | CD3D CD3E CD2 TRAC | 16 | 11 |
| P8 myeloid abundance | CD68 CD14 CSF1R LYZ | 4 | 23 |
| P9 PD-L1 | CD274 | 22 | 17 |
| P10 PD-1 (residual on T abundance) | PDCD1 | 4 | 23 |
Residual = residual of the program score on the abundance score (OLS over patients), re-z-scored. Genes absent from
the matrix are dropped and listed; a program with fewer than half of its genes is not scored.

**Statistic.** Per program k, OLS over patients of the score on cos(2 pi (t - h_k) / 24); t_k = t statistic of the
slope (positive = the predicted timing). Omnibus S = mean of t_k over the scored programs. Null: patient times
permuted (10,000; the same permutation for every program), p = (1 + #{S_null >= S}) / 10,001, one-sided.
- **H146.1 (primary of part 1):** S > 0 with t = internal time (`est_day`).
- **H146.2:** S > 0 with t = recorded resection time.
- **Holm** over H146.1-2 within part 1 (and over the whole A146 family once later parts are scored).
- **D146.1 (descriptive):** per-program t_k and permutation p for both times; S(internal) - S(recorded) with a patient
  bootstrap (2,000); the same omnibus in adjacent-normal libraries (programs P2-P10; P1 descriptive); S as a function
  of the mapping offset (human h = CT - 5 h +- 0, 2, 4, 6 h), shown, not tested.
**Power (stated before scoring):** a single program at n = 42 needs |r| about 0.40 (one-sided 0.05, 80 %). The omnibus
null SD is measured from the permutations and its MDE reported with the result. Small-sample, bulk data: a fail is
"not detected at this design".
**Caveats fixed in advance:** bulk mixes composition and expression (residual scores address immune programs only
partly); the window covers about 10 h, so a program whose predicted peak lies outside it is read mainly through its
slope; the internal time comes from the same libraries' clock genes (the programs exclude them).
**Reading.** H146.1 pass: human breast tumours show the predicted time-of-day organisation against internal time; a
larger S for internal than recorded time supports a clock rather than a scheduling origin. Fail: not detected; the
single-cell part (TNBC) is then registered with its own power before running.

**A146 note 1 (2026-10-27, after part 1 and before any part 2 readout).** The planned part 2 (TNBC, Zhang et al. 2025,
GSE266919) is not estimable: the deposit holds only immune objects (44 patients, 78 samples; no epithelial, fibroblast,
endothelial or mural cells), so no stromal host or internal time can be built, and the immune-only gate failed (A140
Q2). No recorded times. Blind inventory: `paper/exploratory/A146_tme_time/output/part2_prep_*.csv`. Part 2 is
re-planned on the on-disk single-cell cancer cohorts with a stromal host (pan-cancer, breast where present) and will be
registered after its own blind inventory and power.

## A121 (2026-10-27, written before any readout below; PI request to bring the embryo into Figure 2): the embryo clock against the known timing of fetal clock emergence

Transcribed from the WP7 plan (`provisory_plan_paper.md`, A121), with the published reference fixed here.
**Known at registration:** A108 (Qiu 2024 mouse embryo atlas: shared clock only late in gestation, E15-E18.75, H108.1
R 0.52, p 0.016, n 16; gate borderline p 0.08), A108b, A108c (per-trajectory amplitude, phase and p per gestational
half, `figure_2/output/embryo/a108c_timed_halves.csv`; cardiomyocytes not rhythmic in either half, p 0.22 and 0.28).
The pooled host statistic R per half (H121.1) has not been computed. Code base: `figure_2/scripts/embryo/`.
**Published reference (checked 2026-10-27):** fetoplacental PER2::LUC rhythms in utero emerge during pregnancy and
stabilise around E15.5 with peaks at dusk to early night (bioRxiv 10.1101/2023.06.30.547166; later version
10.64898/2025.12.17.695051). Reference interval: ZT12-ZT14 (lights-on ZT0). No cell-autonomous oscillation in E10
hearts, clear at E18 (Umemura et al. 2017; verify before citing).
- **H121.1 (onset within the window):** alignment with clock time is stronger in E17.0-E18.75 than in E15.0-E16.75.
  Per half: the H108.1 statistic R (8 embryos each, batch intercepts, host refitted per half). Statistic R(late) minus
  R(early); null: clock times permuted within batch in both halves (10,000); one-sided p < 0.05. Secondary: the same
  per lineage (hepatocytes, mesenchyme), descriptive.
- **H121.2 (phase against the known fetal timing):** the common offset of the embryo host (host phase minus clock
  time) is converted to mouse ZT with a human-to-mouse frame offset calibrated on the timed mouse datasets that use the
  same orthologue projection (A21 colon, A48 mouse ages, A72 mouse skin): inverse-variance weighted circular mean of
  the per-dataset offsets, SE by bootstrap over datasets and samples. The light schedule is not recorded: predicted ZT
  of the embryo clock peak reported for lights-on at 06:00, 07:00 and 08:00. Pass if, under lights-on 07:00, the 95 %
  interval of the predicted peak ZT (offset SE and embryo host SE combined) overlaps ZT12-ZT14. Uninformative (reported
  descriptively) if the calibrated offset's SE exceeds 1.5 h or the datasets disagree by more than 3 h.
- **Heart, stated in advance:** "consistent with" the heart literature, not "replicates".
- **Family:** secondary; Holm over H121.1-2.
- **Reading:** H121.1 passes and H121.2 passes (or is uninformative): Figure 2 embryo panel, labelled exploratory,
  framed as consistent with and extending known late emergence. H121.1 or H121.2 fails: the embryo goes to Extended
  Data.

### A146 part 2 (2026-10-27, written after the blind preparation and before any program-time readout): pan-cancer single cell, programs against internal time within compartments

**Known at registration** (blind preparation, `paper/exploratory/A146_tme_time/output/part2pc_prep_*.csv`; no program
score was related to any time): 43 on-disk cohorts (`cancer_sc`), 576 tumour patients; stromal host in 23 cohorts, gate
passes in 17; A25 internal times (hospital prior) for 270 patients; circular SD of times median 3.0 h in gate-passing
cohorts. Pool = gate-passing cohorts with >= 6 timed patients: 16 cohorts, 226 timed patients (malignant programs 166
patients in 13 cohorts; immune programs 212 in 14). Omnibus null (2,000 within-cohort time permutations): SD 0.398,
q95 0.658, MDE 0.99 (S units). Part 1 result known (tumour programs follow recorded time, not internal time).
**Design (as the preparation script `part2pc_prep.py`, frozen):** same 10 programs and predicted human acrophases as
part 1, scored within compartment (P1, P2, P3, P9 in malignant cells, CNV call where available, else tumour
epithelium; P4, P10 in T_NK; P5, P6 in myeloid; P7, P8 as log fractions of T_NK and myeloid cells among all tumour
cells, >= 100 cells); gene score = mean within-cohort z of log2(CPM + 1); residual programs as in part 1; >= 50 cells
per compartment, >= 6 patients per cohort x program. Statistic per program: slope t of the stacked OLS of the
within-cohort z-scored score on the within-cohort-centred cos(2 pi (t - h_k)/24), t = per-patient internal time;
S = mean of the 10 t values.
- **H146.3:** S > 0. Null: times permuted within cohort (10,000, seed 146), one-sided p. Holm over H146.1-3 (A146
  family, secondary).
- **D146.2 (descriptive):** per program and per cohort t values; omnibus with P9 also in myeloid; the clock-keeping
  subset (3 cohorts, 31 timed patients: PRAD, LIHC, LUAD) and the non-keeping subset (2 COAD cohorts, 13 patients),
  not powered (|r| 0.49-0.66 per program); the breast cohorts (ERP130387, SRP273482; 17 timed patients); the mapping
  offset sweep (+-2, 4, 6 h).
**Caveats fixed in advance:** no recorded times, so the internal-versus-recorded question of part 1 cannot be asked;
P7 and P8 depend on each cohort's cell capture (CD45-sorted samples are already excluded by the builder); cancer labels
are cohort-level (SRP329496 mixes ovarian and breast).
**Reading:** a pass shows the predicted organisation against internal time across cancers (which part 1 did not find
in bulk breast); a fail with S below its MDE says the internal-time association is absent or below about 1 S unit
across cancers, consistent with part 1.

**A145 note 1 (2026-10-27, after the simulator check and a 2-replicate code check, before the full run; no H145
verdict read).** (i) Shared factor fixed by the registered check (`output/simulator_check.csv`): the A140 factor takes
0.22 (scale 1) and 0.066 (scale 0.5) of the clock-plane variance; the matching simulator values are shared_sd 0.207
(arm S_large, used for H145.2 and arm N) and 0.103 (arm S_small, descriptive); the literal 0.5 / 1.0 were labels of the
A140 sizes. (ii) Technical choices (as `scripts/lib.py`): psi = z2 V (1 - c) / (2c); offsets {0, 3, 6, 9, 12} h and
amplitudes assigned to targets by a per-replicate permutation shared across arms; orderer phase reflected to the clock
plane's rotation sense, offsets relative to the anchors' circular mean; orderers see the 56 clock-gene x anchor
features; CYCLOPS with the template hyperparameters (80 models; not reduced, the published pipeline is what is
tested); a failed two-step fit is excluded from the pairing with a sensitivity counting it as a clocksync win.
(iii) D145.2 (real data) was computed in the code check; it is reported from the full run unchanged.

## A147 (2026-10-28, method; written before any simulation of it; PI decision Q9 c): validation of the pull-versus-order readout

**Readouts (as A132b, made general).** Against a reference group (fast deaths), for every compartment pair with >= 20
shared donors in both groups: the signed deviation of the pair phase arg E[z_a conj z_b] from its nearest pole (0 or
12 h). *Order* = Spearman of the group's against the reference's signed deviations over pairs (null: tissue labels
permuted in the group, 1,000). *Shrink ratio* = median of group / reference deviation over same-pole pairs with
|reference deviation| >= 1 h (1 = no pull, 0 = full pull). Known at registration: A132b values (ventilated post-mortem
order 0.649, ratio 0.58; organ donors 0.257 / 0.40 at n 79, 0.095 / -0.01 at n 298; Hardy 4 0.788 / 0.77) and the A133
grid (raw Fisher-Lee rho 0.70 at displacement fraction 0.44, n 79).
**Simulator:** the A133 generative model (`paper/exploratory/A133_axial_decomposition/scripts/compute_v133_sims.py`:
fast-death offsets, ventilated amplitudes and noise, 12 tissues), 200 replicates per cell, n in {79, 298}:
- pull arms: common-direction displacement fraction f in {0, 0.2, 0.44, 0.8}, offsets intact;
- loss arms: a share q in {0.25, 0.5, 1} of tissues get offsets redrawn uniformly (relations lost), f in {0, 0.44}.
- **V147.1 (order separates loss from pull):** AUC of *order* between the pull arm f = 0.44 and the loss arm q = 0.5,
  f = 0 is >= 0.9 at n 79; and at f = 0.8 with offsets intact, order stays above the q = 0.5 median in >= 80 % of
  replicates.
- **V147.2 (shrink ratio tracks the pull):** median ratio at f = 0 within [0.85, 1.15], and Spearman of the cell
  medians with f over the pull arms <= -0.9, at both n.
- **Pass both:** a clocksync readout (`pull_order`, with a unit test) is added, and the Figure 3 "two modes of
  perturbation" panel uses it with its validation; the A132b group values are reported as computed by the validated
  readout (data unchanged, still post hoc). **Fail:** the readout stays descriptive in Extended Data.
**D147 (descriptive):** all cells: order and ratio distributions, and the order null false-positive rate in the loss arm
q = 1. Code: `paper/exploratory/A147_pull_order/`.

## A122 (2026-10-28, written before any readout below; transcribed from the WP7 plan): controls for the corrected A85 immune cue-following result

**Known at registration:** A85 corrected: immune cue-following index 0.182 vs structural 0.095, difference +0.087,
one-sided p 0.002 (`paper/results/figure_5/output/cancer/`). Immune classes have a small host loading, so the index
|b| / (|a| + |b|) can rise mechanically. A146 (2026-10-27) suggests tumour immune programs follow the hospital clock,
which motivates the control but does not enter it. Gene matching on `feature_name`.
- **H122.1 (locked-set control, primary):** rerun A85 with the A81b locked set in place of the cue set. Statistic:
  (immune minus structural, cue index) minus (immune minus structural, locked index), pooled over atlases as A85.
  Pass: > 0, one-sided p < 0.05 (atlas bootstrap, 2,000).
- **H122.2 (pseudo-cue sets):** 200 expression-matched random 40-gene sets (excluding clock, cue and locked genes).
  Pass: the observed immune-minus-structural difference exceeds their 95th percentile (smallest attainable 1/201).
- **Descriptive:** sign of the difference in each atlas; within atlas, the class index against log pseudobulk depth
  and the difference after adjusting for depth; variance of the cue state's own expression by class.
- **Family:** secondary; Holm over H122.1-2.
- **Reading:** both pass: Figure 2 panel, "immune cells take relatively more of their timing from body-wide signals".
  Otherwise Extended Data, described as a consequence of low host coupling.

**A138 note 6 (2026-10-28; PI request; written after H138.1 and H138.4 were scored and before any purity readout):
H138.8, does the tumour-normal clock match come from the non-malignant fraction?** Motivation: H138.4 passed (pooled
0.365) while malignant host coupling was not shown (H138.1, 0.076) and malignant clock variance is mostly
patient-specific (D138.1). If the bulk match comes from stroma and immune cells, it should weaken with tumour purity.
Secondary, outside the A138 confirmatory family (registered after seeing H138.1/H138.4; the purity association itself
has not been computed).
- **Purity:** PanCanAtlas ABSOLUTE purity of the tumour aliquot (public GDC PanCanAtlas supplement, downloaded and
  parsed only, PROVENANCE.md); fallback where missing: none (pairs without a call are dropped and counted). Sensitivity:
  ESTIMATE-derived purity computed from the same counts (Yoshihara 2013 gene sets and formula), descriptive.
- **Statistic:** within each H138.4 cancer, pairs split at the cancer's median purity; the H138.4 statistic
  (`clocksync.paired.cross` on one projector per cancer, as note 3) computed separately in the low- and high-purity
  halves; D = pooled(low) - pooled(high), pooled over cancers with weights n_pairs, cancers with >= 10 pairs per half.
- **H138.8:** D > 0. Null: purity labels permuted within cancer (10,000), one-sided p = (1 + #{D_null >= D}) / 10,001.
- **Descriptive:** per-cancer D; pair-level association (Spearman of purity with the pair's contribution
  Re(z_T conj z_N e^{-i theta_c}) / (sd_T sd_N), theta_c the cancer's pooled phase); the ESTIMATE sensitivity; the
  null SD and MDE of D, reported with the result.
- **Reading:** pass: the bulk tumour-normal match is carried mainly by non-malignant cells, and Figure 5 says so.
  Fail with D near 0 and the MDE small: the match does not depend on purity, so malignant cells contribute to it.
  Fail with a large MDE: not estimable at this design.

## A148 (2026-10-28, written after A115 was scored and before any decomposition readout; PI request): how much of the psoriatic clock delay does the keratinocyte state mix explain, step by step?

**Known at registration:** A115 (uninvolved minus healthy core-cue offset -1.81 h; uninvolved more differentiated,
Hedges g +0.47); the paired lesional delay (A78: -1.70 h core-cue, core -1.43 h); archive A61/A69 (single-cell state
mix predicts -0.20 to -0.36 h with healthy phases, up to -0.5 to -1.0 h with coarse states); A63/A64 (healthy
keratinocyte clock advances along differentiation: cycling ~0, basal +0.6 h, most differentiated +1.4-1.6 h). No
deconvolution of the A115 bulk data has been done.
**Data:** the four A115 psoriasis series with healthy controls (`skin_lesional/processed_a115/`, lesional rows from the
original processed matrices), one sample per person per group.
**Estimator.**
1. Per sample, keratinocyte state shares by non-negative least squares on CPM: reference profiles = the A63/A65 bin
   pseudobulks (kc_cycling, kc_d1-d4; `pseudobulk_bins_both.csv`, healthy donors) plus non-keratinocyte references
   (fibroblast, endothelial, immune, melanocyte) from the same single-cell sources; genes = the 500 most bin-variable
   genes excluding clock, cue and locked genes. Shares renormalised within keratinocytes.
2. Composition-only predicted keratinocyte clock state of a group, as A69: sum_k g_k z_k, g_k = mean share of bin k x
   the bin's clock-gene gain (A60 projector), z_k = the A63 healthy per-bin complex phase and amplitude.
3. Predicted step offset between groups A and B: arg(sum g_k^A z_k / sum g_k^B z_k) in hours; per series, pooled with
   the A54 weights (arg of the weighted sum of unit complex ratios).
4. Interval: 2,000 resamples of samples within group and series, each paired with a draw from the A63 donor bootstrap
   (500 draws, recycled).
**Acceptance check (before reading any predicted offset):** the deconvolved lesional-minus-uninvolved shares move as in
single cell (A61): cycling share up and most-differentiated share down, pooled. If not, A148 stops (deconvolution not
valid) and is reported as such.
- **H148.1:** predicted lesional-minus-uninvolved offset < 0 (composition predicts a delay); one-sided, interval upper
  bound < 0.
- **H148.2:** predicted uninvolved-minus-healthy offset > 0 (composition predicts an advance, opposite to the observed
  -1.81 h); one-sided, interval lower bound > 0.
- **D148 (descriptive):** share explained per step, predicted / observed (observed: A115 within-group offsets on the
  same samples, lesional-uninvolved and uninvolved-healthy), residual = observed minus predicted; per-series values;
  shares per bin per group.
- **Family:** secondary; Holm over H148.1-2.
- **Reading:** H148.1 and H148.2 pass: the state mix explains part of the lesional step but cannot explain the
  uninvolved shift (wrong sign), which then needs another cause (site or batch of the healthy controls, or a trait
  mechanism outside keratinocyte differentiation). H148.2 fails with a negative prediction: the uninvolved shift is
  partly a differentiation effect, contrary to the transit model's sign.
Code: `paper/exploratory/A148_psoriasis_decomposition/`.

**A148 note 1 (2026-10-28, PI decision):** A148 is closed after its acceptance check failed; the post hoc gene-set
variant is not registered. The psoriasis differentiation argument rests on the single-cell results (A61, A63, A69)
and the A115 differentiation descriptive.

## A149 (2026-10-28, written after A138 H138.1/H138.4/H138.5/H138.8 were scored and before any A149 readout; PI request): does the malignant clock follow the patient's own normal epithelium in single cells, and does the link run through the tumour stroma? Paired tumour / adjacent-normal cohorts pooled, per cancer type.

**Known at registration:** H138.4 (TCGA bulk, tumour vs own normal 0.365, null 0.122; per cancer BLCA 0.52 to UCEC
0.08); H138.8 (no purity dependence, D +0.025, MDE 0.16); H138.5 (malignant vs adjacent-normal epithelium after
removing the stromal host, 52 patients, 6 cohorts with >= 6 pairs: 0.583, null 0.443 +- 0.092, p 0.069, MDE 0.677);
H138.1 (malignant coupling to the in-sample stromal host 0.076 [0.014, 0.139], MDE about 0.2); A87 (alignment of
malignant vs normal epithelium to the stromal host 0.27 vs 0.64, 79 patients); erratum of 2026-10-28 (0.365^2 = 0.13
shared variance, back-of-envelope). Counted today, no statistic computed: the 43 `cancer_sc` cohorts hold 137
patients with tumour (malignant or tumour) epithelium and adjacent-normal epithelium in 25 cohorts, 85 of them with a
tumour stromal compartment. The unreferenced single-cell pairing (H149.1) has not been computed on any cohort.
**Data:** the `cancer_sc` cohorts plus new public cohorts with tumour and adjacent-normal samples of the same patients
(search list `paper/exploratory/A149_paired_tumour_normal/candidates.md`; the cohorts added are fixed in a dated note
before any readout), processed counts only, built with the unchanged `cancer_sc` rules (patient = donor, untreated,
not CD45-sorted, median UMI >= 1,500, majority chemistry, marker compartments, >= 50 cells per pseudobulk). Copy-number
calls for new cohorts with `cancer_common.studies.calls`; malignant compartment as A138 (`q3_lib.study_table`: called
malignant cells where >= 5 patients are called, else tumour epithelium). Eligible cohort: >= 4 paired patients (the
A138 limit of 6 is lowered for pooling; per-cohort centring and permutation keep the null exact).
**States:** per cohort, tumour and adjacent-normal epithelium on one projector (`clocksync.paired.projector`, A138
`sc_pair`), centred within cohort x chemistry; Poisson-corrected |complex correlation| (`corr_stat`); pooled =
mean over cohorts weighted by paired patients; null = rows of the second compartment permuted across patients within
cohort (10,000; the same permutation indices for the pooled and per-cancer statistics); one-sided p, floor 1e-4.
- **H149.1 (single-cell, purity-free H138.4):** malignant vs own adjacent-normal epithelium, unreferenced states, >
  null. Specificity companion: 200 expression-matched pseudo-clock gene sets (A52 rule); the clock statistic is read
  as clock-specific only if it exceeds their 95th percentile.
- **H149.2 (beyond the tumour stroma):** the same after removing the patient's tumour stromal host from both states
  (H138.5 statistic, host = `q3_lib.stromal_anchors`), patients with a host, > null.
- **H149.3 (through the tumour stroma):** malignant vs the patient's tumour stromal host state, > null (host rows
  permuted within cohort). Companion: adjacent-normal epithelium vs the same host (expected positive; reported).
- **D149.1 (descriptive):** retained share = (S2 - null mean 2) / (S1 - null mean 1) on the patients with a host, with
  a 2,000-draw patient bootstrap (within cohort) interval.
- **D149.2 (descriptive, PI: cancers need not behave alike):** per cancer type with >= 12 paired patients over its
  cohorts: S1, S2, S3, null mean / SD, p, excess over null with bootstrap interval; I^2 of the excess across cancer
  types. Pre-specified contrast: cancers whose tumours keep a clock-shaped signal in TCGA (A83: kidney, liver, lung,
  thyroid, prostate, breast) minus those at null (colorectal, bladder, endometrial, ovarian), excess difference with
  bootstrap interval; other types are reported but not in the contrast.
- **Power (before scoring, dated note):** null mean / SD and MDE of S1-S3 on the final data (permutations only, no
  observed statistic), as the A138 power step.
- **Family:** secondary (outside the confirmatory family and the A138 family); Holm over H149.1-H149.3.
- **Reading:** H149.1 pass: the malignant clock follows the patient's own normal tissue in single cells, without
  stromal admixture (time of surgery and patient-specific clock traits are not separated: both are "the body").
  Then D149.1 near 1 with H149.2 pass: the link does not run through the tumour stroma; D149.1 near 0 with H149.3
  pass: it runs through the stroma; intermediate: both. Any failed test is "not detected at this design" with its MDE,
  never "uncoupled". D149.2 states where cancers differ; no cancer-type claim is confirmatory.
Code: `paper/exploratory/A149_paired_tumour_normal/`.

**A149 note 1 (2026-10-28, before any A149 data were downloaded or read): cohorts added.** From `candidates.md`
(search of 2026-10-28), the cohorts with a verified count of >= 4 paired patients, droplet or well-based UMI data
(cells or nuclei), patient and tissue type in the metadata, no stated neoadjuvant treatment and no fraction sorting
that removes epithelium: GSE164789, GSE229253, GSE333596, GSE153935, GSE310916, GSE200972, GSE117570 (lung);
GSE206785, GSE150290, GSE167297 (gastric); GSE199654, GSE160269 (oesophageal squamous; GSE199654 is epithelium only
and enters H149.1 only); GSE178481, GSE242299 (kidney); GSE282701, GSE242889, GSE290298 (liver); GSE293189 (bladder;
the 2 paragangliomas are dropped, no epithelium of origin); GSE144240 (cutaneous SCC); GSE193337, GSE176031 (prostate).
Data rules: a cohort whose supplied matrix fails the integer check of the legacy GEO builder (raw UMI) is dropped;
multiple tumour samples or lesions of a patient are pooled as `cancer_sc` does; GSE164789 keeps all lesion stages
(AIS / MIA / invasive), malignant cells taken from copy-number calls where called. Not added (unverified pairs or
format, Seurat-only objects, treated subsets mixed): Joanito 2022, Moorman 2025, OPSCC (CELLxGENE), GSE173950,
GSE201347, GSE185344, GSE202642, GSE210347 (re-compiled, overlaps). Copy-number calls of the new cohorts are written
to `$SCCC_DATA/cnv_a149/` so that the cohort set of the frozen A83-A138 analyses (`CNV_STUDIES`) is unchanged.
Cancer classes for D149.2: lung, liver, kidney, prostate, thyroid, breast = TCGA clock-keeping; colorectal, bladder =
at null; gastric, oesophageal, skin, pancreatic, head and neck, ampullary = not in the contrast.

**A149 note 2 (2026-10-28, before any A149 statistic or null was computed): operational choices.**
(i) Host = stromal compartments of the tumour sample only (`tumour:` fibroblast / endothelial / mural): learned A17
anchors among them where >= 2, else all of them (the `q3_lib.stromal_anchors` fallback would also admit adjacent-normal
stroma, which is not the "tumour stromal host" of H149.2/H149.3). (ii) Estimability per cohort and statistic as A138
note 4 (ii) (corrected variances > 0, null mean <= 1); H149.2/H149.3 need >= 5 patients with a host (`Dataset.host`).
(iii) Existing cohorts use the stored `cancer_sc` tables (as A138); new cohorts are built from counts with
`cancer_sc.build`. (iv) Pseudo-clock companion of H149.1: genome-wide sums of the two epithelial compartments per
patient (same cells and labels, recomputed from the counts), 200 expression-matched sets per cohort
(`clocksync.paired.pseudo_waveforms`, matched on the cohort's normal epithelium), pooled per set with the H149.1
weights; 95th percentile of the pooled sets. (v) Bootstrap (D149.1, D149.2): patients resampled within cohort,
statistic recomputed, excess = statistic - the unresampled null mean; I^2 from the per-cancer excess and its bootstrap
SE. (vi) Seed stream 149.

**A146 note 2 (2026-10-28; PI request; written after part 1 was scored and before any of the quantities below was
computed): origin of the recorded-time pattern in GSE233242 tumours.** Systemic cue tied to a shared wake time vs
procedural or scheduling confounds. Data and programs as part 1 (42 patients, tumour TPM, recorded resection time
08-18 h; GEO records only age and subtype). Secondary family "A146 checks", Holm over H146.4-H146.5; the rest is
descriptive. Permutations: the part 1 stream (seed 146, 10,000).
- **H146.4 (systemic glucocorticoid cue, positive control of the wall-time reading):** glucocorticoid-response score
  (mean z of log2 TPM+1: TSC22D3, FKBP5, DDIT4, ZBTB16, KLF9; PER1 left out, it is a clock gene of the internal-time
  estimate) falls with recorded time (cortisol peaks after waking and declines over the day): t of the slope < 0,
  one-sided permutation p. Normal tissue: descriptive.
- **H146.5 (scheduling):** the part 1 omnibus S against recorded time with every program residualised on age and
  subtype (dummies) stays > 0 (same permutations). Descriptive: Spearman of age, Kruskal-Wallis of subtype, and
  Spearman of the T-cell and myeloid abundance scores (P7, P8) with resection hour.
- **D146.2 (handling):** immediate-early score (FOS, FOSB, JUN, EGR1, ATF3, IER2) and hypoxia score (VEGFA, SLC2A1,
  PGK1, CA9, ADM, BNIP3) against resection hour (t, permutation p); S with programs residualised on both.
- **D146.3 (trend vs predicted phases):** observed S against S with the 10 predicted hours shuffled among programs
  (10,000 random assignments): a generic monotone trend scores under any assignment, the predicted structure only
  under the registered one. Plus S restricted to the immune-cell programs (P4-P8, P10) and to the others (P1-P3, P9).
- **Reading:** H146.4 pass and H146.5 pass: a wall-time systemic cue is present in the same samples and scheduling by
  age or subtype does not explain the pattern, consistent with a shared-wake-time systemic cue (fasting is not
  separated: it is also tied to the hour). H146.4 fail: the positive control for a wall-time systemic signal is
  absent (power stated), weakening that reading. H146.5 fail: scheduling explains part. D146.2 large: handling
  artefact candidate. D146.3 shuffled assignments reaching the observed S: a generic trend, not the predicted phases.
Code: `paper/exploratory/A146_tme_time/scripts/checks_part1.py`.

**A149 note 3 (2026-10-28, after the cohort builds and before any A149 statistic or null): build outcomes.**
Dropped by the note 1 integer rule: GSE206785 and GSE290298 (log-normalised matrices only). Dropped as a data rule:
GSE333596 (S1/S3 tumour and normal libraries are not separable; for the remaining A/B libraries GEO does not state
whether A and B are different patients or two regions of one patient, so the patient unit is not identifiable).
Kept with their builder notes: GSE199654 has 8 pairs (not 13: the authors' epithelial matrix holds 59 adjacent-normal
cells); GSE293189 8 pairs (paragangliomas "other"); GSE178481 tumour thrombus "other" and one treated patient
excluded by the treatment rule; GSE160269 CD45- libraries only enter (the builder drops CD45+-sorted samples).
Several cohorts (GSE242299, GSE176031, GSE153935, GSE293189) sit largely below the unchanged median UMI >= 1,500
sample rule and may contribute few or no patients. Build scripts:
`paper/exploratory/A149_paired_tumour_normal/scripts/build_geo_a149_part{1,2}.py`.

**A149 note 4 (2026-10-28, power step; nulls only, no observed statistic computed).** Units (`output/units.csv`):
S1 190 patients in 27 cohorts (130 from existing cohorts, 60 new: GSE164789 10, GSE229253 9, GSE150290 10, GSE144240
8, GSE310916 5, GSE282701 5, GSE242889 5, GSE200972 4, GSE167297 4); S2 61 patients / 9 cohorts; S3 51 / 7 estimable.
Lost to the unchanged builder rules: kidney (no epithelial compartment called in GSE178481 / GSE242299), GSE160269,
GSE199654, GSE153935, GSE293189, GSE193337, GSE176031 (depth or adjacent-normal epithelium below the cell and patient
rules); GSE117570 fails in the state computation (ZeroDivisionError; left out). Pooled nulls (mean +- SD, q95, MDE):
S1 0.446 +- 0.042, 0.517, 0.552; S2 0.539 +- 0.092, 0.691, 0.768; S3 0.410 +- 0.078, 0.542, 0.608. Per cancer type
(>= 12 patients): S1 colorectal 82, lung 42, gastric 21, liver 14 (MDE 0.57-1.02); S2 and S3 colorectal and lung
only. D149.2 classes: keeping 71 (S1), at null 82 (colorectal only). H149.1 is powered for an excess of about 0.11
over its null; H149.2 / H149.3 only for about 0.2-0.23, and the per-cancer estimates are descriptive at best.

## A150 (2026-10-28, written after A149 was scored and before any quantity below was computed; PI request): is the healthy single-cell coupling of Figure 2 clock-specific, and is it inflated by ambient RNA?

**Motivation.** A149: within-patient pairing between compartments was matched or exceeded by expression-matched
pseudo-clock gene sets, so pairing in single-cell data can be carried by non-clock structure. The Figure 2 atlases
have faced the gate (permuted acrophases of the clock genes, shuffled pairings), not expression-matched non-clock genes;
ambient RNA would add clock-shaped shared covariance that the gate cannot reject. **Known at registration:** the A17
healthy fits (`results/figure_2/output/atlas_fits/a17_healthy_fits.csv`); A109 organ-level pseudo-clock excess (GTEx
0.463 vs 0.154). No pseudo-clock or ambient readout of any atlas has been computed.
**Data:** the gate-passing A17 healthy atlases with a genome-wide twin (same rows, `_atlas/genome.py` `CLOCK_TWIN`):
lung_hlca, organ_breast, organ_muscle, organ_colon, and heart if `heart_global_genome` reproduces the A17 `heart` rows
(checked before any readout; else left out). organ_adipose fails the gate. The other 19 gate-passing atlases have no
genome-wide table and are not covered (stated as a limitation). Anchors and rows: the A17 learned rule, unchanged.
**Part 1, pseudo-clock control.**
- Per atlas, 200 expression-matched pseudo-clock waveforms (`clocksync.paired.pseudo_waveforms`, matched on mean log
  fraction over all rows of the atlas, clock-like genes excluded); for each, the full A17 anchored fit (host and targets
  on the pseudo genes, `Dataset(..., waveform=w).fit(gate=False)`, point estimates).
- **H150.1:** the mean anchored coupling over all compartments of all covered atlases with the clock waveform exceeds the
  same mean under the pseudo sets (set j pooled across atlases as set j): p = (1 + #{pseudo >= clock}) / 201, floor
  0.005; pass at p < 0.05.
- **D150.1:** per compartment and per cell family (the Figure 2 family classes, `backbone/compute_family_summary.py`):
  clock coupling, pseudo median and q95, excess = clock minus pseudo median. A family whose clock coupling is below its
  pseudo q95 is "not clock-specific at this design".
**Part 2, ambient RNA (marker-based over-correction, SoupX-like).**
- Soup profile of a donor s_d = its summed counts over all compartments, as fractions. For compartment k, exclusion genes
  N_k = the 30 non-clock genes with the largest ratio of atlas-level soup fraction to k's atlas-level fraction (soup
  fraction >= 1e-5). Contamination rho_dk = sum_{g in N_k} x_dkg / (depth_dk sum_{g in N_k} s_dg), capped at 0.5: an upper
  bound, since true expression of N_k in k also counts.
- Corrected clock counts x' = max(x - rho depth s, 0), depth' = depth (1 - rho); the A17 learned-rule fit is redone.
- **D150.2 (descriptive, reading pre-specified):** per compartment and family: median rho, coupling before and after,
  change; Spearman of compartment couplings before vs after. Reading: median |change| < 0.05 in every family and the
  family order unchanged -> ambient RNA does not drive the backbone; otherwise the affected families are named.
- **Family and reading:** secondary; H150.1 alone is tested. H150.1 fail or a family below its pseudo q95 is brought to
  the PI before any Figure 2 wording changes. Code: `paper/exploratory/A150_atlas_specificity/`.

## A151 (2026-10-28, written before any fibroblast pairing readout; PI request; run after A150): do tumour fibroblasts follow the patient's own normal-tissue fibroblasts, clock-specifically?

**Known at registration:** tumour fibroblasts keep a clock-shaped signal (A86 clock share about 0.46-0.50; D138.1 stroma
host-coupled share 0.53); H138.5 fibroblast companion (1 cohort, 6 patients: 0.240, p 0.87, uninformative); A149
(epithelial pairing within patients matched or exceeded by pseudo-clock sets; H149.1 not clock-specific). No tumour vs
adjacent-normal fibroblast pairing has been computed on the A149 cohort set.
**Data and states:** the A149 cohort set and tables (notes 1-3), compartments `tumour:fibroblast` and
`adjacent_normal:fibroblast`, patients with both; cohorts with >= 4 such patients. States, centring, Poisson-corrected
|complex correlation|, n-weighted pooling, within-cohort permutation null (10,000) and estimability as A149 S1.
- **H151.1:** pooled tumour vs own adjacent-normal fibroblast pairing > null (one-sided, floor 1e-4).
- **H151.2 (clock specificity):** the pooled clock statistic exceeds the pooled statistic of 200 expression-matched
  pseudo-clock sets (A149 note 2 (iv), genome-wide sums of the two fibroblast compartments, matched on the cohort's
  adjacent-normal fibroblasts): p = (1 + #{pseudo >= clock}) / 201, floor 0.005.
- **Companion (descriptive):** the same for endothelial cells (`tumour:endothelial` vs `adjacent_normal:endothelial`).
- **Power (before scoring, dated note):** units and null of H151.1 from permutations only.
- **D151.1 (descriptive):** per cancer type with >= 12 patients: statistic, null mean, p. No bootstrap intervals (the
  A149 patient bootstrap is invalid for |correlation|).
- **Family:** secondary; Holm over H151.1-H151.2. **Reading:** both pass -> tumour stroma keeps a clock-specific,
  patient-specific relation to the patient's normal-tissue stroma (the stromal counterpart of H138.4). H151.1 pass and
  H151.2 fail -> pairing carried by non-clock structure, as A149. H151.1 fail -> not detected at this design (MDE).
Code: `paper/exploratory/A151_fibroblast_pairing/`.

**A151 note 1 (2026-10-28, power step; nulls only).** Fibroblast pairs: 98 patients in 17 cohorts (colorectal 42 / 7,
lung 24 / 3, others < 12 per type); null 0.504 +- 0.058, q95 0.601, MDE 0.650 (detectable excess about 0.15).
Endothelial companion: 103 patients / 15 cohorts, null 0.521 +- 0.068, MDE 0.695. Scoring follows A150.

**A149 note 5 / A151 note 2 (2026-10-28, POST HOC: written after H149.1 and H151.2 were scored; PI decision).** The
Poisson-corrected |correlation| exceeds 1 for many pseudo-clock sets at stromal pseudobulk depths (A151), which makes
the registered specificity comparison uninterpretable. Post hoc bounded comparison: the uncorrected |complex
correlation| of the centred states, |sum za conj(zb)| / sqrt(sum |za|^2 sum |zb|^2), for the clock waveform and for the
same 200 pseudo-clock sets per cohort (same seeds, patients and genome tables as the registered companions), pooled
with the same weights; p = (1 + #{pseudo >= clock}) / 201. Applied to A149 S1 (malignant vs adjacent-normal epithelium)
and A151 (fibroblasts; endothelial descriptive). Per-cohort clock and pseudo values are saved, with the share of
corrected values above 1. Reported as post hoc; the registered verdicts stay as scored.

## A152 (2026-10-28, written after A150 and before any residency split of heart or breast; PI request): is lymphocyte coupling carried by tissue-resident rather than circulating cells?

**Known at registration:** A150 per compartment: lymphocyte coupling heart 0.556, breast CD4 / CD8 T 0.442 / 0.456,
plasma 0.810, lung T 0.159, B 0.100, colon mostly 0. Archive A6 (pre-clocksync estimator, A4 class rules; healthy HLCA,
198 donors): resident-like 0.61 [0.35, 1.00] vs circulating-like 1.00 (boundary), difference -0.39 [-0.63, 0.00]; fibrotic
lung the reverse in direction. The A6 values are not comparable with the A17 lung T value (0.159); the reason is not
established. Blood atlases fail the gate.
**Residency state per cell (A6 rule, unchanged):** x = log1p(1e4 count / total UMI); score = mean x of CD69, ITGAE,
CXCR6 minus mean x of SELL, CCR7, KLF2, S1PR1; tertiles within block (study x assay x suspension type), ties at random;
resident-like and circulating-like tertiles kept, intermediate dropped.
**Part 1, lung (HLCA healthy, existing `hlca_tstates` builder):** mode "state" (T class replaced by
`lung:T_resident_like`, `lung:T_circulating_like`, >= 20 cells), A17 learned anchors of lung_hlca, current estimator.
D152.0 (reconciliation, descriptive): mode "T" under the A17 setup vs A17 lung:T; then the A6 split under the A17 setup.
**Part 2, heart (heart_global) and breast (breast_kumar):** per-cell marker counts streamed from CELLxGENE Census
(release 2025-11-08, as the pseudobulks) for the T / NK-containing classes (heart: `lymphocyte`, `natural killer cell`;
breast: every T-cell, NK T and NK class); cells scored as above; the lymphocyte classes of the A17 table are replaced by
`<tissue>:lymph_resident_like` and `<tissue>:lymph_circulating_like` pseudobulks (>= 20 cells); other classes, donors,
strata and A17 learned anchors unchanged. A tissue enters if >= 15 donors have both states.
- **H152.1:** pooled over the part 2 tissues that enter, resident-like minus circulating-like coupling > 0. Statistic:
  mean over tissues of c_res - c_circ (A17 fit, point estimates). Null: within each donor the two state labels swapped
  at random (1,000 permutations, refit), one-sided p; MDE from the null. Pass at p < 0.05.
- **D152.1 (descriptive):** per tissue (lung, heart, breast): c_res, c_circ, difference with the permutation p, relative
  amplitudes, pseudobulk depth per state (low depth can push either estimate to the boundary), resident-like share of
  lymphocytes per donor; the CD8-only split where >= 15 donors.
- **Specificity:** genome-wide state pseudobulks do not exist; only the gate is applied (stated limitation).
- **Reading (PI):** pass -> residency carries lymphocyte coupling in heart / breast. A positive but non-significant
  difference stays descriptive. A negative difference (circulating > resident) in heart would point to intravascular
  cells following the body's clock. Secondary family. Code: `paper/exploratory/A152_lymphocyte_residency/`.

## A17 note: anchors relearned after the class-rule fix (2026-10-28, PI; dated data-rule fix, before any refit under the new tables is read)

**Why:** commits 281bb26 and 588d37d fixed the cell-class rules: case-insensitive "T cell" / "B cell" put mast, goblet,
tuft, myofibroblast and club cells into T_other / B / heart lymphocyte. New classes were added at the same time: heart
mast / adipocyte / mesothelial, and melanocyte. Several classes also changed content, for example colon:T_other
16,444 -> 3,388 cells and ileum:T_other 24,091 -> 6,244. The frozen learned anchors of `common/anchors_a17.csv` were
computed on the contaminated classes. A first patch added the new compartments as non-anchors and kept every frozen
flag. That contradicts the learned rule, which selects whole families (an independent review found ia:mast unflagged
while myeloid is an anchor family there). The patch is discarded, unread: no fit made under it was looked at.
**Rule (unchanged code):** `extended_data/scripts/host_identity/compute_coherence.py` then `compute_anchors.py` are
rerun on the fixed tables. The learning set, family scores, leave-one-study-out selection and the "< 3 -> no anchor"
rule are as registered in A17 step 2. The relearned table replaces `common/anchors_a17.csv`, and the pre-fix table is
kept as `common/anchors_a17_pre_classfix.csv`. The fixed stromal-vascular rule is unchanged apart from compartments that
appear or disappear. Every result anchored on the A17 table is then rerun and reported old vs new, verdict changes
first: the A17 healthy fits, the family summary, A150, the Figure 4 contrasts, A152 heart, cue_units (A85 / A122), and
Tabula Sapiens rebuilt. The old values stay in the docs as superseded.

## A123c (2026-10-28, written before any truth is computed from out-of-window donors alone; PI request after the independent review): the working-hours organ map against a truth that shares no donor with the test

**Why:** the review of 2026-10-28 (C9) noted that truth (i) is a cosinor on all fast deaths, so it includes the 111
window donors whose untimed estimates are tested. After the addendum, the method's real-data advantage rests on A123,
so the overlap is removed.
**Known at registration:** A123 H123.1-3 (offsets r 0.945 vs truth (i), amplitudes 0.976, donor Spearman 0.613); A123b
(CoFE genome -0.283, CoFE clock -0.417, CYCLOPS -0.296, CHIRAL -0.165 vs truth (i) in the window design). No truth from
out-of-window donors has been computed.
- **Truth (iii).** GTEx v10 fast deaths (Hardy 1-2) who died outside 08-18 h. Per-tissue cosinor on recorded time,
  reduced by rank-1 complex SVD (A123 truth (i) code, same reference tissue and centring), for the A123 tissues with
  >= 20 such donors. Tissues below that are dropped and counted. The truth and test donor sets are disjoint.
- **Estimates:** the stored A123 clocksync estimates, unchanged. Their waveform was cross-fitted on the out-of-window
  donors, which overlap with the truth (iii) donors (stated).
- **Sensitivity (split):** the out-of-window donors are split at random into two halves (seed 123, stratified by tissue
  count). Half 1 gives the waveform, clocksync is refitted on the window donors, and the estimates are scored against a
  truth (iii) built from half 2 only.
- **H123.5:** offsets vs truth (iii), circular correlation > 0.7 and tissue-permutation p < 0.05 (A123 `compare`,
  10,000 permutations, 200 tissue bootstraps).
- **H123.6:** log2 relative amplitudes vs truth (iii), Pearson > 0.7.
- Each passes on its own, as H123.1-2 did.
- **Descriptive:**
  - truth (i) vs truth (iii) agreement;
  - the A123b orderers (CHIRAL, CoFE genome and clock, CYCLOPS) and the timed oracle vs truth (iii), from their
    stored consensus estimates;
  - the split sensitivity;
  - a prior-only baseline for donor timing: median absolute error of the A123 tangent-coordinate hour on the 111 window
    donors vs always predicting the window centre (13:00), after the same common offset.
- **Reading:** H123.5-6 pass -> the A123 organ map does not depend on the test donors entering the truth. A fail would
  remove or reframe the Figure 1 panel (PI decision).
Code: `paper/exploratory/A123_working_hours/scripts/compute_truth_iii.py`.

## A140 Q1 note 2 (2026-10-28, a design change after a result, approved by the PI after the independent review): the registered core-clock seed

**Why:** the review addendum notes that every Q1 CYCLOPS run so far used the seed without the 14 core clock genes.
CYCLOPS is normally seeded with clock-informed genes, so the failure is open to a straw-man charge. The registered seed
grid includes the core-clock variant ("with and without the core clock genes"), but the simulator gate was only run
without it.
**Known at registration:** the revised gate without core genes:
- at scale 1, CYCLOPS vs truth median 0.060, 0 of 50 >= 0.7;
- at scale 0, every CYCLOPS fit fails to start;
- clocksync gate 24 of 50 at scale 1 and 20 of 20 at scale 0;
- on the GTEx window design (A123b), CYCLOPS with its human seed: offsets -0.296.
- **Gate with core genes:** arm D, K = 1, factor scale 1, seed 1,500 = the 14 core clock genes plus the top 1,486 others
  (`q1_lib.seed_set(core=True)`), 80 models, calibrated factor (revision settings), 50 replicates. The worlds are the
  same as the gate without core genes (base seed 1500), so the comparison is paired. Pass: median CYCLOPS-vs-truth
  circular correlation >= 0.7 (as registered).
- **Descriptive:** scale 0 with the core seed, 20 replicates. Does CYCLOPS start and recover the clock when there is no
  factor?
- **If the scale-1 gate passes:** H140.1-3 are run and scored as registered (arm B, K = 1, seed 500 with core genes; the
  clocksync gate in arms A, B, C; exact binomial thresholds; Holm over the three).
- **If it fails:** H140.1-3 are not scored. Reading: at this design CYCLOPS 2.0 does not recover a coupled clock even
  when seeded with the core clock genes.
Code: `A140_CORE=1` in `q1_run.py` / `q1_score.py`; outputs `output/q1rev/{gate,primary}_core/`,
`simulator_gate_core.csv`, `primary_scores_core.csv`.

## A154 (2026-10-28, method; written before any simulation of it; roadmap 4.1, review S1): gate calibration under shared non-clock structure

**Why:** every "not estimable" in the paper rests on the gate. A145 arm N gave a gate pass in 4 of 30 runs without a
clock, while the review's isolated rank-1 test was calibrated (0.040, 0.051). A140 showed gate sensitivity falling from
20 of 20 to 24 of 50 at full factor size.
**Simulator:** `clocksync.simulate`, current code and defaults of the paper (gate n_perm 1,000).
- **Designs** (8): donors n in {30, 100} x compartments K in {4, 8} x depth in {single-cell pseudobulk (median A17 clock
  UMI per compartment), bulk}.
- **Arms per design:**
  - N0: no clock, no factor.
  - N1: no clock, plus a shared non-clock donor factor. Rank in {1, 3}; loading on 10 % or 30 % of genes, with clock
    genes eligible; size scale in {0.25, 0.5, 1, 2} x the A140 Q1 calibration (gate-passing tumour immune eigenvalues).
  - P: one coupled clock with c in {0.3, 0.6, 0.9}, without a factor and with the factor at scale 1, rank 1, 30 %.
- **Replicates:** 200 per cell, seed 154.
- **Readouts:**
  - (i) gate pass rate;
  - (ii) under a global null (two arms with the same c = 0.6), uniformity of the contrast p (KS, and the share
    < 0.05);
  - (iii) coverage of the 95 % bootstrap interval of c in arm P without a factor;
  - (iv) coverage of the boundary lower bound when the fit is at the boundary.
- **V154.1:** gate pass rate in N0 <= 0.075 in every design (one-sided exact binomial against 0.05 not rejected at
  p < 0.01).
- **V154.2:** the share of contrast p < 0.05 under the global null is <= 0.075 in every design.
- **V154.3:** coverage of the 95 % interval of c >= 0.90 in every design and every c.
- **Descriptive:** the N1 map (pass rate vs rank, loading and size), and P sensitivity with vs without the factor.
- **Reading:** V154.1 fails -> the gate is anti-conservative and every gate-based claim needs a calibrated threshold.
  Large N1 pass rates at realistic sizes -> a passed gate does not by itself show a clock; claims lean on the
  pseudo-clock companions. Loss of P sensitivity with the factor -> "not estimable" is worded as "clock not separable
  from shared non-clock structure", not "no clock". Code: `paper/exploratory/A154_gate_calibration/`.

## A155 (2026-10-28, method; written before any simulation of it; roadmap 1.5, review C7): bias from reading anchors on a host that contains them

**Why:** `Dataset` reads an anchor's amplitude |z| and offset arg z on the full host, which includes the anchor itself
(`dataset.py` ~323). The anchor's own noise then enters its regressor. In A123 every tissue is an anchor.
**Simulator:** `clocksync.simulate`, known per-compartment amplitude, offset and coupling.
- **Designs:**
  - anchors m in {3, 4, 6, 10} plus 4 non-anchor targets, n in {30, 100}, depth single-cell and bulk;
  - an A123-like all-anchor design: 40 compartments, n 111, GTEx bulk depth, amplitudes and offsets from A123 truth (i).
- **Replicates:** 200, seed 155.
- **Readouts:** for anchors, bias and RMSE of log |z|, arg z (h) and c under (a) the current full host and (b) a
  leave-one-out host; non-anchors as reference. For the A123-like design, the circular correlation of offsets and the
  Pearson of log2 amplitudes against truth, under (a) and (b).
- **Decision rule (method):** if under (a) the anchor bias exceeds 5 % in |z| or 0.25 h in offset, or the A123-like
  offset r changes by > 0.02 between (a) and (b), at A17-typical designs (m 3-6, n 30-100), leave-one-out becomes the
  anchor readout in clocksync. Every anchored result is then rerun in roadmap phase 3, and A123 in phase 8. The PI
  confirms (Q13, 0.4). Otherwise the current readout stays, with the measured bias stated in METHOD.
Code: `paper/exploratory/A155_anchor_self_inclusion/`.

## A124 (2026-10-28, written before any download of GSE106800; roadmap 6.1, from WP8): does the readout detect misalignment imposed on human muscle?

**Why:** the review asks for a pre-specified human perturbation positive. A simulated night shift is a within-subject,
imposed misalignment.
**Data:** Wefers et al. 2018, GSE106800. 14 healthy men, crossover: a control protocol, and a protocol with the
behavioural (sleep, activity, meal) cycle shifted by 12 h. Vastus lateralis biopsies at 2 clock times per condition
(Methods 08:00 and 20:00; Results 07:00 and 19:00, an unexplained discrepancy). Affymetrix HuGene 2.1 ST, public series
matrix (47 arrays after the authors' QC). Only the processed matrix is used.
- **Data rule:** fixed in a dated note after the inventory and before any clock readout: probe-to-gene mapping, the
  arrays per subject, and subjects with all 4 biopsies.
- **Estimator (the A115 / A78 estimator, A53 array convention):**
  - log2 intensity, equal gene weights;
  - core state = GTEx 14-gene waveform; cue state and locked state = the A81b sets (`figure_3/output/gtex_cues/a81b_sets.csv`)
    restricted to the genes present;
  - one projector per gene set over all arrays, states centred within subject;
  - per condition, the core-cue offset = arg E[c_core c_cue*] over that condition's arrays (`S.ccorr`); the core-locked
    offset likewise.
- **H124.1 (primary):** the core-cue offset changes between control and misalignment.
  - Statistic: |Δ| (h, wrapped) = |offset(misaligned) − offset(control)|.
  - Null: condition labels swapped within subject (all 2^n assignments if n <= 16, else 10,000).
  - Two-sided p < 0.05.
- **H124.2 (specificity, same permutations):** |Δ core-cue| − |Δ core-locked| > 0, one-sided p < 0.05.
- Holm over H124.1-2.
- **MDE:** from the null distribution of |Δ|, reported.
- **Descriptive:**
  - Δ core-cue and Δ core-locked with bootstrap intervals (subjects resampled, 2,000);
  - the core-state difference between the two biopsy times per condition (does the muscle core clock follow the shifted
    behaviour?);
  - the same for the cue state;
  - coherence per condition.
- **Reading:**
  - H124.1 passes -> the readout detects imposed misalignment in human tissue (the human perturbation positive), and
    the WP8 cue analyses (A126-A127) can become claims.
  - H124.1 fails -> "not detected at n = 14, MDE ...". A126-A127 stay descriptive.
  - H124.2 fails with H124.1 passing -> the change is not specific to cue-driven genes.
Code: `paper/exploratory/A124_night_shift/`.

**A124 data rule (dated note 2026-10-28, after the blind inventory `exploratory/A124_night_shift/inventory.md` and before any clock readout):**
- **Arrays and subjects:** the public series matrix has 47 arrays from 12 subjects, not 14.
  - 11 subjects have all 4 biopsies and form the primary set. The exact test has 2^11 = 2,048 assignments, so the
    smallest p is 4.9e-4.
  - RAS9 lacks the misaligned evening array and is descriptive only.
- **Biopsy times:** as recorded, 07:00 and 19:00.
- **Probe mapping:** main-category probe sets whose gene_assignment lists the symbol; probe sets listing the symbol
  alone, where they exist; mean of log2 per gene. NR1D1 is kept through its two THRA;NR1D1 read-through probe sets.
- **Gene names:** mapped with the project's canonical alias table (`clocksync.genes.canonical`), as everywhere else.
  Expected coverage: core 14 / 14, cue 40 / 40, locked 39 / 40 (CLUHP3 absent).

## A156 (2026-10-28, written after A124 was scored and before any download of GSE48113; roadmap 6.2): forced desynchrony in human blood, a timed within-subject misalignment

**Why:** A124 (muscle, 2 biopsy times per condition) had MDE 7.3 h and did not detect a change. GSE48113 (Archer et al.
2014, forced desynchrony on 28 h days) samples blood 7 times on day 1 (sleep aligned with the central clock) and on day 4
(sleep displaced by about 12 h) in 22 subjects. It is the largest timed within-person misalignment series found.
**Known at registration:**
- Archer 2014 reports far fewer rhythmic transcripts relative to the melatonin phase on day 4, and a reduced amplitude
  of rhythmic genes.
- A124: the muscle core clock keeps its direction under a behavioural shift and loses amplitude.
- No GSE48113 data have been read.
- **Data rule:** a dated note after a blind inventory (subjects with both days, samples per day, sample clock times or
  hours since wake, probe mapping), before any clock readout. Only the processed series matrix is used.
- **Estimator (A124 / A115, A53 array convention):** core, cue (A81b) and locked (A81b) states. One projector per set
  over all arrays. States centred within subject x day. Per subject and day, core-cue offset = arg Σ c_core c_cue*
  over that subject-day's samples; core-locked likewise.
- **H156.1 (primary):** the core-cue offset changes from day 1 to day 4.
  - Statistic: |circular mean of Δ_i| in hours, unweighted, where Δ_i = offset(D4) − offset(D1) per subject.
  - Null: day labels swapped within subject (10,000 random assignments, or exact if n <= 16).
  - Two-sided p < 0.05.
- **H156.2 (specificity, same permutations):** |mean Δ core-cue| − |mean Δ core-locked| > 0, one-sided.
- Holm over H156.1-2. The MDE is reported from the null, defined as in A124.
- **Descriptive (timed):**
  - per subject and day, the core-state phase against sample time (a cosinor on the 7 samples), the change from D1 to
    D4, and the same for the cue state. Does the blood core clock follow the displaced sleep, or stay with the central
    clock?
  - the untimed Δ core-cue against the timed core-minus-cue phase change (circular correlation over subjects): a check
    of the untimed readout against timing;
  - coherence per day.
- **Reading:**
  - H156.1 passes -> the untimed core-cue readout detects imposed misalignment in human tissue, a perturbation positive
    (blood).
  - H156.1 fails -> not detected, with the MDE.
  - The timed descriptive states whether the readout tracks the known displacement either way.
Code: `paper/exploratory/A156_forced_desynchrony/`.

**A156 data rule (dated note 2026-10-28, after the blind inventory `exploratory/A156_forced_desynchrony/inventory.md` and before any clock readout):**
- **Data:** GSE48113 has 287 arrays from 22 subjects; 21 have both days.
- **Primary set:** subjects with >= 5 arrays on each day (19 subjects); 10,000 random within-subject day swaps.
  Sensitivities: the 15 subjects with 7 + 7 arrays (exact test), and all 21.
- **Day labels:** samples are assigned by the recorded day field. BB0121's day-1 point 7 is labelled "out of phase" but
  stays day 1, as recorded.
- **Probe mapping:** non-control probes whose symbol maps to the gene; mean of the log2 values present. A gene with all
  probes missing in an array is imputed with its mean within that subject-day; the frequency is reported.
- **Gene names:** the inventory showed that `clocksync.genes.canonical` misses the old platform names known from A124.
  The primary therefore uses canonical plus the 8 A124 aliases (core 14, cue 39, locked 38); canonical alone is a
  sensitivity. This is fixed before the readout, unlike in A124.
- **Cosinor (descriptive):** recorded clock time, 24 h period. BB0264 day-1 point 7 (no time) is dropped from the
  cosinor only.

## A109 / A138 note (2026-10-28, roadmap 5.1-5.2, written after the independent review and before any responsive-gene readout): second specificity companion for H109.1 and H138.4

**Why:** the pseudo-clock pool excludes glucocorticoid and immediate-early genes (review C2). A shared acute or GC pulse
across organs (ventilated donors), or across the two samples of one surgery (TCGA), is therefore not in the null.
`clocksync` now provides responsive-gene pseudo sets (expression-matched draws from GC targets and immediate-early
genes, `pseudo_waveforms(kind="responsive")`) and `drop_genes`. These readouts are secondary. They do not change the
registered verdicts, and they qualify the readings.
- **H109.1-R:**
  - In ventilated donors, C_clock exceeds the median of 200 responsive-gene sets carrying the GTEx waveform loadings.
    One-sided add-one p < 0.05.
  - Sensitivity: H109.1 with PER1 and PER2 dropped, against the registered matched pseudo-clock null.
- **H138.4-R:** tumour-normal pairing statistic vs 200 responsive-gene sets, with the same rule. Plus:
  - (i) PER1/PER2-free pairing vs the matched pseudo null;
  - (ii) an immediate-early positive-control set (FOS, FOSB, JUN, EGR1, ATF3, IER2, DUSP1, NR4A1), whose pairing
    shows how strongly shared surgery stress pairs the two samples;
  - (iii) the pairing null permuted within cancer x tissue source site (and within plate where annotated) instead of
    within cancer.
- Holm over H109.1-R and H138.4-R (secondary family).
- **Responsive pool:** small (about 27 genes present per dataset), so responsive sets overlap heavily. The overlap (mean
  Jaccard) is reported.
- **Reading:**
  - C_clock (pairing) above the responsive sets -> the shared state is not a GC or immediate-early response.
  - Not above -> "clock-gene-shaped shared state, not separable from a shared acute response".
  - A drop in the H138.4 pairing under TSS/plate stratification -> part of the pairing is batch.
Code: the A109 and A138 script folders (new `*_responsive.py`).

## A157 (2026-10-28, method; written before any simulation of it; roadmap 4.5): is the Jammalamadaka-SenGupta correlation usable on an evenly spread 4-h grid?

**Why:** A113 muscle gave JS -0.854 with the registered estimator and Fisher-Lee 0.576 (chosen afterwards). The muscle
pass is post hoc unless the JS failure is a property of the design.
**Simulator:** the A113 muscle design (GSE108539: 10 subjects, their actual sampling grid of clock times, 57 samples).
True internal times = clock time + a subject offset (SD 1 h) + noise. Estimated times are generated by the A113
estimator on `clocksync.simulate` counts at the observed muscle depth and coupling, so that the estimate error matches
the observed median error (1.30 h, calibrated before any correlation is read). 500 replicates, seed 157. The adipose
design (GSE104674) is the comparison.
- **Readouts:** distribution of JS and Fisher-Lee between estimated and true times; share of replicates with JS < 0
  while Fisher-Lee > 0.5.
- **V157.1:** in the muscle design, JS < 0 in >= 20 % of replicates while Fisher-Lee stays > 0.5 in >= 80 %; in the
  adipose design, JS < 0 in < 5 %.
- **Reading:** V157.1 passes -> JS is ill-defined on this design; the Fisher-Lee muscle value is reported as the design-
  appropriate estimator, with the post hoc label and this simulation. V157.1 fails -> muscle goes to Extended Data as a
  registered fail.

## A158 (2026-10-28, method; written before any simulation of it; roadmap 4.6): why is the clocksync coupling drop in hepatocyte Bmal1 KO larger than the timed oracle's?

**Why:** A129: clocksync d* -0.293 vs timed oracle -0.186 (D145.2).
**Simulator:** the Jouffe 2022 design (genotypes, tissues and timepoints as in `jouffe2022_bmal1ko`), with the liver's
true coupling lowered from its Alb-Cre value to values {0.5, 0.6, 0.7, 0.8} x control, and other tissues unchanged.
200 replicates per value, seed 158.
- **Readouts:** bias of clocksync d* and of the timed oracle's coupling change against the true change; the effect of
  the liver being an anchor or not, and of the KO liver's lower amplitude on the truncation (z2 floor).
- **Descriptive:** the true change at which each estimator reproduces the observed values.
- **Reading:** states which estimator is biased and why. It is a sentence in Methods, not a test.

## A159 (2026-10-28, written after A150 and A155 were scored and before any pseudo-gene fit of the atlases below; roadmap 3.3, Q13 defaults 0.1-0.2): clock-specific coupling for every healthy atlas, the reported Δc readout, and the eligibility table for disease primaries

**Why:** A150 found that raw coupling is not clock-pure (pseudo-clock mean 0.249 vs clock 0.404, 5 atlases), and that
lymphocyte, myeloid and mural coupling is not clock-specific. Under the Q13 defaults, the paper reports Δc = c − pseudo
median, and disease primaries are restricted to compartments with clock-specific coupling in controls (0.2).
**Data:** every gate-passing healthy A17 atlas, and the control arm of every Figure 4 contrast, for which a genome-wide
table at the A17 rows can be built from data on disk (`atlas_wg`, Census group pseudobulks, `_atlas/genome.py`).
- **Row check:** before any readout, the clock-gene columns of the genome-wide table must reproduce the A17 table rows
  exactly (donors, compartments, counts). Atlases that fail or cannot be built are listed and not covered.
- **Anchors:** the relearned A17 table (2026-10-28). Code: current clocksync, `anchor_readout="loo"` (Q13 0.4 default).
- **Fits:** per atlas, 200 expression-matched pseudo-clock waveforms and 200 responsive-gene waveforms
  (`pseudo_waveforms(kind="matched"|"responsive")`, matched on mean log fraction over the atlas rows). Each gets the full
  anchored fit (`fit(gate=False)`, point estimates), as A150 part 1.
- **H159.1 (replication of H150.1 on new atlases):** over the covered atlases not in A150, the mean anchored coupling
  under the clock waveform exceeds the same mean under the matched sets (set j pooled across atlases), with
  p = (1 + #{pseudo >= clock}) / 201.
- **H159.2:** over the same atlases, the pooled stromal-vascular coupling (families fibroblast_stromal and endothelial,
  compartment mean) exceeds its matched-pseudo q95, by the same rule.
- Holm over H159.1-2.
- **D159.1 (the deliverable table):** per atlas x compartment:
  - c_clock, the matched median and q95, the responsive median and q95;
  - Δc = c_clock − matched median;
  - flag `clock_specific` = c_clock > matched q95 AND c_clock > responsive q95.
  For contrast control arms, the same per compartment. Per family: median Δc, and the share of clock-specific
  compartments.
- **Reading:**
  - H159.1-2 pass -> the A150 result generalises, and Figure 2 reports Δc with families ranked on it.
  - D159.1 is the eligibility table that A153 (disease contrasts) uses: a compartment enters as a primary only if
    clock-specific in the control arm, and the MDE on d* is <= 0.2.
Code: `paper/exploratory/A159_specific_coupling/`.

## A153 (2026-10-29, written before any contrast below is rerun and before the lead read D159.1 (computed by the Phase 3 runner on 2026-10-28 21:36; the eligibility rule itself was registered in A159); roadmap 6.4-6.5, Q13 defaults 0.1-0.2): corrected disease contrasts, eligibility-restricted primaries and equivalence tests

**Why:** the review of 2026-10-28 found three problems with the A17 disease contrasts:
- several were read through a host that failed its gate in one arm (A111 H111.2);
- truncated estimates hid boundary rows (S3);
- the primaries sat on classes whose coupling is not clock-specific (A150: lymphocyte, myeloid, mural).

As a result, "no change" was never a tested claim. A153 reruns every Figure 4 contrast with the corrected code and
defines one family of primaries before it is seen. The registered A17 results stay reported as registered.
**Data:** the 36 targets of the A17 registry (`datasets/_atlas/a17_targets.csv`) at their original labels, blocks and
rules, rebuilt under the 2026-10-28 class rules (rule-stamped tables).
**Code:** current clocksync; the relearned A17 anchors (learned rule); `anchor_readout="loo"`; seeds from the target name.
Permutation and bootstrap counts are those of `compute_contrasts_a17.py`.

**Per target:**
1. **Gates.** Run the pooled model-free gate, then `contrast(require_host_gate_each_arm=True)`.
   - If either fails, every row of the target is **not estimable**.
2. **Estimates per compartment.**
   - Untruncated, debiased c per arm, d, d* (common spread), p, p*, se_d_star, the boundary type per arm, and the
     stromal-vascular companion. The companion uses the same permutations as before; neurons stand in for
     stromal-vascular in brain.
   - The MDE on d* is computed at 80 % power with `power.mde` (n_rep 100, n_perm 200) at the real design. Damping s
     scales the host-locked variance by s², so the MDE on d* is c*0 − s*² c*0 / (s*² c*0 + 1 − c*0), where:
     - c*0 is the compartment's pooled coupling at the common spread;
     - s* is the damping at which power reaches 0.8, interpolated linearly on the scale grid.
     - If power does not reach 0.8 at s = 0, the MDE is set to "> max".
3. **Pseudo-clock companion (descriptive).** Where a genome-wide table exists (as in A159), run the same contrast for 200
   expression-matched pseudo-clock waveforms (point estimates):
   - Δc per arm;
   - the quantile of the clock d* among the pseudo d*.

**Eligibility (the A153 primary family; decided now, applied mechanically).** A compartment of a target is a primary if
all of the following hold:
- it is `clock_specific` in that target's control arm in D159.1 (A159);
- its MDE on d* is ≤ 0.2;
- the target is estimable.

Targets without a D159.1 control-arm row (no genome-wide table) have no A153 primaries. They are reported with the
reason. Brain targets go to the Extended Data with the A141 stability result stated; they are not in the family.

**Tests:**
- **H153.1 (effect).** Run Holm over all A153 primaries of all targets together, on the two-sided p*. An **effect** is
  a Holm-adjusted p* < 0.05.
  - The specificity qualifier, reported next to every effect, is the companion p (target minus stromal-vascular) and the
    pseudo-d* quantile.
  - An effect whose companion p ≥ 0.05 is labelled "not target-specific".
- **H153.2 (equivalence, TOST).** A primary without an effect is **equivalent** (an estimable null) if its 90 % interval
  d* ± 1.645·se_d_star lies within ±0.15. Otherwise it is **inconclusive**.
- **H153.3 (pooled family equivalence).** Run this for each of the families fibroblast_stromal and endothelial,
  separately:
  - Per target, take the inverse-variance mean of its A153-primary compartments in the family, with SE equal to their
    mean se_d_star (no credit for correlation).
  - Pool over targets with `meta(method="hk")`.
  - The family is equivalent if the 90 % HK interval lies within ±0.15, and shows an effect if the 95 % interval
    excludes 0.
  - This needs at least 3 targets; otherwise it is reported as not estimable.

**Readout:**
- one row per target × compartment, with class (not estimable / not eligible / effect / equivalent / inconclusive) and
  the reason;
- a summary of each class per disease;
- a comparison with the A17 registered verdicts, with changed verdicts first.

**Reading:**
- An effect in H153.1 is the pre-specified disease positive the paper needs.
- Equivalence (H153.2-3) turns "coupling is preserved in disease" into a tested claim, bounded at |d*| < 0.15.
- Inconclusive and not-eligible rows are reported as such, never as "unchanged". A failed gate is uninformative.

**Code:** `paper/exploratory/A153_corrected_contrasts/`. Promoted to `figure_4/scripts/resilience/` if Figure 4 uses it.

**A153 note (2026-10-30, PI Q15 answer (a); a design change made after D159.1 was seen, approved by the PI; written
before the family statistic below is computed and before any A153 contrast is run).** The per-compartment eligibility
above is replaced by a family-level rule. Everything else in A153 is unchanged.
- **Source.** The A159 fits of the 5 healthy atlases with genome-wide tables (lung_hlca, heart, organ_breast,
  organ_muscle, organ_colon), as on disk: `exploratory/A159_specific_coupling/output/d159_compartments.csv` (c_clock),
  `pseudo_matched.csv` and `pseudo_responsive.csv` (200 sets each). No refit.
- **Family.** `family()` of `results/extended_data/scripts/host_identity/anchor_rule.py`, applied to the healthy and the
  disease compartment names alike.
- **Statistic.** For family F, C_F = the mean of c_clock over every compartment of F in the 5 atlases (equal weight per
  compartment). For pseudo set j, C_F,j = the mean of set j's coupling over the same compartments (set j pooled across
  atlases, as H150.1 / H159.2). Compartments with a failed pseudo fit are dropped from that set's mean.
- **Eligible family.** C_F > the q95 of C_F,j under the matched sets AND > the q95 under the responsive sets, with
  at least 3 compartments in at least 2 atlases. p = (1 + #{C_F,j ≥ C_F}) / 201 is reported per kind; no Holm (an
  eligibility screen, not a test).
- **A153 primary.** A compartment of an estimable target whose family is eligible, with MDE on d* ≤ 0.2. Brain targets
  stay out of the family (as above). The D159.1 per-compartment flag is reported next to each row, descriptively.
- **Sensitivity (descriptive, not used for eligibility).** The same rule on the 3 atlases that passed the A159 exact row
  check (heart, organ_breast, organ_colon); lung_hlca (467 rows with count differences) and organ_muscle (1 row) did not.
- **Option (b)** (genome-wide tables for the disease control arms) may follow later as a confirmation (PI); it is not
  part of A153.

## A160 (2026-10-29, written before any lesional statistic below is computed; roadmap 6.9): is the psoriatic core-cue delay present within keratinocytes (composition-controlled)?

**Why:** the psoriatic delay is read in bulk skin:
- paired lesional minus uninvolved core-cue offset of −1.70 h (A78);
- core offset of −1.43 h, robust to gene subsets and 600 pseudo sets (psoriasis robustness, 2026-10-28).

Lesions change cell-type and keratinocyte-state composition. Archive A61/A69 attribute −0.2 to −1.0 h to the state
mix. A delay read within keratinocyte-only pseudobulks, and within keratinocyte differentiation states, removes these
two routes.
**Data:** single-cell psoriasis with paired lesional (PP) and uninvolved (PN) skin of the same patients: GSE173706 (on
disk; 11 PN-PP pairs, 8 healthy). GSE228421 (5 patients) is added only if it passes the same checks; no other series.
**Pseudobulk:** cells labelled keratinocyte by the existing annotation (`_skin_sc`, `psoriasis_sc_kc_states`), summed per
library, all genes, CPM, then the A78 processing.
- **Composition-controlled variant:** the same per library × A63 differentiation bin (`psoriasis_sc_bins`).

**Estimand:** as A78 (paired, A81b GTEx cue set restricted to genes present, GTEx core waveform; `_cues.py`):
- the core-cue offset arg E[c_core c_cue*] in PP minus PN, on the 11 pairs, in hours (wrapped);
- the bin variant uses the inverse-variance (bootstrap) weighted circular mean over eligible bins.

**Feasibility and power, computed blind before any PP statistic.** Only PN and healthy libraries are used:
- **Data rule:** a library enters with at least 200 keratinocytes, and a library × bin with at least 50. A bin is
  eligible with at least 8 complete pairs. All core genes and at least 60 % of the cue set must be detected.
- **Null and MDE:** split PN-vs-PN plus healthy by random assignment into pseudo-pairs, with the same pair count and
  estimator. Null SD and MDE (one-sided 0.05, 80 % power) per variant.
- **Gate:** a variant is confirmatory only if its MDE is ≤ 1.7 h (the bulk paired delay). Otherwise it is descriptive,
  reported with its MDE, and has no pass/fail wording.

**Tests (confirmatory variants only; Holm over those that qualify):**
- **H160.1:** the keratinocyte-only PP − PN core-cue offset is < 0, by a one-sided exact within-patient label swap
  (2^11 = 2,048), p < 0.05.
- **H160.2:** the same for the differentiation-bin variant.

**Descriptive:**
- the PN − healthy offset (trait, as A115);
- the whole-skin pseudobulk of the same libraries, giving the bulk-like reference on identical samples;
- the core offset alone;
- the NR1D1 leave-out.

**Reading:**
- H160.1 passes → the delay is within keratinocytes, not cell-type composition.
- H160.2 passes → it is not differentiation-state mix either.
- A whole-skin delay with a null keratinocyte offset at adequate MDE → composition is the route.

Code: `paper/exploratory/A160_psoriasis_keratinocyte/`.

## A127 (2026-10-29, written before SCAN-B is downloaded; WP8, roadmap 7.5; H127.3 single-cell part not registered here): breast tumour clock coherence vs proliferation and ER status at scale (SCAN-B)

**Why:** A97 found lower clock-gene coherence in high-proliferation TCGA tumours (H97.2) and lower coherence in
ER-negative breast tumours (exploratory). SCAN-B (GSE96058, about 3,000 primary breast tumours, RNA-seq, with PAM50, ER
and follow-up) replicates both at about 3x the TCGA breast size. A124 failed, so cue results stay descriptive (WP8 rule);
H127.1-2 are not cue results.
**Data:** GSE96058 processed expression matrix and sample characteristics, public. One sample per patient: technical
replicates are dropped, keeping the first listed. Only primary tumours with PAM50 and ER recorded.
- **Data rule:** if the processed matrix is not on a log-FPKM / TPM-like scale that can be converted to the A97 scale (log2
  CPM-like within cohort), use the published scale with z-scores within cohort, and state it. All 14 clock genes and
  ≥ 90 % of the G2M set must be present.

**Definitions** (frozen from A97, `figure_5/scripts/cancer/compute_proliferation.py`):
- proliferation score = mean z of Hallmark G2M genes, clock genes removed;
- purity proxy = mean z of COL1A1, COL1A2, COL3A1, DCN, LUM, PDGFRB, PTPRC, CD3E, CD14, LYZ, CD74 (high = low purity);
- coherence = `decompose.components` gene coherence (A87) of the 14 clock genes over a group of tumours.

Tertiles are by rank within stratum.

**Tests:**
- **H127.1.** Strata: PAM50 subtype × purity-proxy tertile, keeping strata with ≥ 30 tumours.
  - T = mean over strata of coherence(high proliferation tertile) − coherence(low).
  - Null: tumour order permuted within stratum, 2,000 times. One-sided p (loss) < 0.05.
- **H127.2.** Coherence(ER−) − coherence(ER+) within purity-proxy tertile, mean over tertiles; ER labels permuted within
  tertile, 2,000 times. One-sided p < 0.05.
  - Sensitivity, not tested: the same within proliferation tertile, to ask whether ER adds anything beyond proliferation.
- Holm over H127.1-2.

**Specificity companion (same permutations):** 200 expression-matched pseudo-clock sets (14 genes, matched on mean and
variance, clock and G2M genes excluded; `clocksync` controls). T_pseudo per set.
- A pass is labelled "clock-specific" only if T is also below the 5th percentile of T_pseudo.
- Otherwise the result is reported as "coherence loss not specific to the clock".

**Exploratory (pre-specified):**
- survival as A97: Cox with age, tumour size, node status and PAM50, using the per-tumour off-plane score; overall
  survival;
- coherence per PAM50 subtype;
- the purity-proxy tertile companion as A97.

**Reading:**
- H127.1 passes, clock-specific → A97 replicates at scale, within subtype and purity; Figure 5 panel.
- H127.1 passes, not specific → proliferating tumours lose coherence in any co-regulated set; the A97 wording changes to
  "general transcriptional decoherence".
- A failure at this size (MDE reported) bounds the A97 effect.

**Code:** `paper/exploratory/A127_scanb/`. Data: `$SCCC_DATA/scanb_gse96058/` with PROVENANCE.md.

## A126 (2026-10-29, written after a metadata-only scout (`exploratory/scouts_2026-10-29/a126_scout.md`, no clock statistic) and before any clock readout; WP8, roadmap 7.5): breast epithelial clock coherence across the menstrual cycle

**Why:** tumour clock coherence falls with proliferation (A97; replicated at scale in A127). Each luteal phase drives
breast epithelial proliferation. If physiological proliferation lowers epithelial clock coherence as it does in tumours,
the clock gives a mechanistic reading of lifetime menstrual cycles as a breast cancer risk factor. Progesterone also
regulates clock genes directly, so a proliferation control within tissue comes first (H126.0). A124 failed, so all cue
readouts are descriptive (WP8 rule).

**Definitions** (fixed in the WP8 plan):
- **luteal score:** mean z of PGR, TNFSF11, WNT4, AREG in HR+ luminal cells (per-donor log CPM pseudobulk; z within
  study);
- **proliferation:** the A97 G2M set; cycling cells by the A24 top-quartile rule;
- **coherence:** A87 `decompose.components` gene coherence of the 14 clock genes over a group of donors;
- **pseudo-clock companion:** 200 expression-matched 14-gene sets (mean and variance; clock, G2M and luteal genes
  excluded), the same permutations, as A127. A pass is "clock-specific" only if T is below the pseudo 5th percentile.

**H126.0 (proliferation control; gates the reading of H126.1).**
- Data: endometrium GSE266265 (60 donors, cycle phase and serum progesterone; per-cell export on disk).
- Per donor, an epithelial pseudobulk of cycling cells and one of non-cycling cells.
- Donor eligible with ≥ 30 cycling epithelial cells; the test needs ≥ 15 eligible donors, else infeasible.
- T0 = coherence(cycling) − coherence(non-cycling) over eligible donors. Null: within-donor label swaps, exact or
  10,000. One-sided p < 0.05 (loss).
- Descriptive: endometrial coherence, proliferative minus secretory phase; the same breast test on Kumar + Reed per-cell
  extractions where ≥ 15 donors reach 30 cycling cells.

**H126.1 (primary).**
- Data: premenopausal donors pooled, study as stratum:
  - Kumar HBCA, rebuilt by the authors' `broad_cell_type` so that HR+ luminal (lumhr) is separate;
  - Reed HBCA, donors aged ≤ 45 (no menopause field);
  - Murrow GSE198732 premenopausal donors (public raw counts, parsed in Python only).
- Donor eligible with ≥ 50 HR+ luminal cells.
- Tertiles of the luteal score within study.
- T1 = the n-weighted mean over studies of coherence(top tertile) − coherence(bottom), on HR+ luminal pseudobulks.
- Null: luteal labels permuted within study, 2,000 times. One-sided p < 0.05 (loss).
- Feasibility: ≥ 10 donors per pooled tertile.
- **Score validity, checked and reported before any clock reading:**
  - spread of the score;
  - Spearman correlation with HR+ cell number (if |ρ| > 0.3, tertiles are cut on the residual after regressing on log cell
    number);
  - direction with age (expected to fall).

**H126.3 (obesity).** Over the eligible premenopausal donors with BMI (Kumar, Reed, Murrow), coherence of HR+ luminal
pseudobulks in BMI ≥ 30 vs < 30. Labels permuted within study, two-sided p < 0.05. Infeasible if there are fewer than
10 obese donors.

**Multiplicity:** Holm over H126.1 and H126.3. H126.0 is a gate, not part of the family.
- If H126.0 fails or is infeasible, H126.1 is reported but read as "luteal state, cause open (hormone or proliferation)".

**Descriptive:**
- luminal progenitor and basal coherence by luteal tertile;
- HR+ phase relative to the stromal host (A23 y4);
- core-cue offset vs luteal score (H126.2 of the plan, descriptive here);
- parity and BRCA carriers (Reed);
- negative control: Kumar postmenopausal donors split by the same score (no luteal biology expected);
- Murrow without hormonal-contraception users;
- GTEx breast, women < 45, aggregates only.

**Reading:**
- H126.0 passes and H126.1 passes, clock-specific → physiological proliferation lowers epithelial clock coherence as in
  tumours. This is a Figure 5 panel, linked in the Discussion to menstrual cycles as a risk factor.
- H126.1 passes without H126.0 → luteal state matters, mechanism open.
- A failure is reported with its MDE.

**Code:** `paper/exploratory/A126_menstrual/`.

## A161 (2026-10-29, written after a metadata-only scout (`exploratory/scouts_2026-10-29/endometrium_scout.md`) and before any expression of the series below is read): is endometrial clock coherence lower in the secretory phase? (confirmation of the A126 descriptive)

**Why:** in GSE266265, endometrial epithelial clock coherence across donors was 0.653 in the proliferative phase and
0.126 in the secretory phase (p 0.0005; A126 descriptive, not registered). Cycling cells showed no loss, which points to
a hormone (progesterone) effect on the clock rather than proliferation. An independent series is needed before any
claim.
**Data (primary):** GSE51981 (GPL570, whole endometrium).
- Phase from the tissue label: proliferative 64 vs secretory 82 (early + mid + late); samples with unknown phase are
  dropped.
- **Data rule:** use the processed series-matrix values; if there are none, the GEO SOFT values. No CEL processing in R.
  Gene = the probe with the highest mean over all samples. All 14 clock genes must be mapped.
- Patient overlap with GSE6364 / GSE4888 is checked from sample metadata before any readout. An overlapping
  replication series is dropped.

**Estimand:** A87 coherence (`decompose.components`) of the 14 clock genes over the samples of a phase group, on log2
values z-scored per gene over all samples (no Poisson term; microarray).

**H161.1:** coherence(secretory) − coherence(proliferative) < 0.
- Phase labels permuted within endometriosis status (case / control), 10,000 times. One-sided p < 0.05.
- **Specificity:** 200 matched pseudo-clock 14-gene sets (mean and variance; clock-like genes excluded), the A127
  variance-matched draw, on observed labels. "Clock-specific" only if T is below the pseudo 5th percentile.
- MDE = 2.49 × null SD.

**Sensitivities (not tested):**
- (i) **Composition:** every gene (clock and pseudo) residualised on an epithelial score (mean z of EPCAM, KRT8, KRT18,
  CDH1) and a stromal score (mean z of VIM, PDGFRB, COL1A1, DCN), fitted over all samples, then the same T.
- (ii) Controls without endometriosis or pelvic pathology only.
- (iii) Early vs mid secretory, each against proliferative.
- (iv) Replication series: GSE6364 (11 vs 26) and GSE141549 (24 vs 43, platform stratum, "medication" samples
  excluded), the same T and p within each.

**Reading:**
- H161.1 passes, clock-specific, and survives (i) → the secretory-phase loss of clock coherence replicates in an
  independent cohort. This is a candidate physiological-perturbation panel: a hormone state that decouples a tissue
  clock.
- Passes but not clock-specific → a general secretory transcriptional reorganisation.
- Fails → the A126 observation is not confirmed (MDE reported).

**Code:** `paper/exploratory/A161_endometrium_phase/`.

## A162 (2026-10-29, method; written after A154 V154.3 failed and before any simulation below; Q17 default (a)): coverage and bias of c at realistic anchor coupling

**Why:** A154 put every compartment at the same nominal c, so the host was built from weakly coupled anchors. ĉ was then
biased low and intervals undercovered (V154.3). Real anchors (stromal-vascular) have c about 0.6-0.8. The question is
whether c and its 95 % bootstrap interval are adequate at the anchor coupling the paper actually uses.
**Simulator:** as A154 arm P: the same 8 designs, depths, seeds scheme (`SeedSequence([162, design, cell, rep])`), no
factor, 200 replicates per cell.
- Anchors: K/2 compartments with c_anchor in {0.6, 0.8}.
- Targets: K/2 with c_target in {0.1, 0.3, 0.6, 0.9}.
- Truth c = the nominal value, defined as in A154.

**Readouts:**
- target coverage of the 95 % interval;
- median bias;
- the boundary share;
- anchor coverage, read on the leave-one-out host (the Q13 0.4 default).

**V162.1:** target coverage ≥ 0.90 in every design and cell with c_anchor = 0.8.
**V162.2:** the same with c_anchor = 0.6.
**Descriptive:** bias vs c_anchor, with the A154 cells (c_anchor = c_target) as reference.

**Reading:**
- V162.1-2 pass → c intervals are reported as they are, with the caveat that they are conditional on a well-coupled
  anchor host. The A154 failure is stated as the weak-anchor limit.
- V162.1 passes but V162.2 fails → intervals are reported only where the anchors' estimated c is ≥ 0.8, else as
  contrasts.
- V162.1 fails → absolute c is reported without intervals (ranks, Δc, d*), and a correction is developed (Q17 b).

**Code:** `paper/exploratory/A154_gate_calibration/` (`--kinds P162`), output `output/a162/`.

## A163 (2026-10-30, method; written after A162 was scored and before any validation simulation below; Q17 option (b), required by the registered A162 reading): a host-reliability correction of c, validated on fresh seeds

**Why:** A162 failed (V162.1 19/32 cells). ĉ is shrunk by a roughly constant factor (0.87-0.92 with anchors at
0.8, 0.80-0.88 at 0.6) that does not shrink with n, as errors-in-variables attenuation by the reliability of the
estimated host. Per the A162 reading, absolute c is reported without intervals until a correction is validated.
**Development (not confirmatory).** Candidate corrections are developed on the A162 seeds
(`SeedSequence([162, ...])`) and on the A154 cells only. Candidates include:
- (i) removing the host posterior uncertainty from the residual variance, psi_c = psi − |z|² tr(S̄_u) / 2;
- (ii) dividing c by the host reliability R = 1 − tr(S̄_u) / V;
- (iii) a reliability read from the anchors' leave-one-out couplings (Spearman-Brown on the anchor composite);
- (iv) design-calibrated rescaling by simulation at the real design (`host.simulator`), with the host re-simulated.

The interval is the donor-bootstrap percentile interval of the corrected c (host refitted, the correction applied
within each replicate). One candidate is frozen: committed to clocksync as an opt-in readout (`coupling_corrected`,
`coupling_corrected_lo/hi`), with the commit hash recorded in a dated note here before the validation runs. The default
readout is unchanged.
**Validation (confirmatory; fresh seeds `SeedSequence([163, design, cell, rep])`, 200 replicates per cell):** the A162
simulator and its 8 designs, with anchors at c_anchor ∈ {0.6, 0.8} and targets at c_target ∈ {0.1, 0.3, 0.6, 0.9}
(32 cells per c_anchor), no factor.
- **V163.1:** coverage of the corrected 95 % interval ≥ 0.90 in every design × c_target cell at c_anchor 0.8.
- **V163.2:** the same at c_anchor 0.6.
- **V163.3:** |median bias| of the corrected c ≤ 0.05 in every cell at both c_anchor values.
- **V163.4 (no inflation where the default is fine):** at c_target 0.1, coverage of the corrected interval ≥ 0.90 in every
  design (both c_anchor values).

**Descriptive:**
- the A154 equal-c cells (weak anchors; 100 replicates), as the weak-anchor limit;
- the anchors' own corrected c on the leave-one-out host;
- the share of corrected values > 1.

**Reading:**
- V163.1-4 pass → c is reported corrected, with intervals, conditional on a gate-passing anchored host with c_anchor
  ≥ 0.6 (Q17 b closed).
- V163.1 and V163.3 pass at 0.8 only → corrected intervals only where the anchors' corrected c is ≥ 0.8; elsewhere
  ranks, Δc and d*.
- Otherwise → c stays reported without intervals (A162 reading), and the correction is not used.
- Contrasts (d, d*) are not changed by A163: at a common host they stay calibrated (V154.2).

**Code:** `paper/exploratory/A163_reliability_correction/`; the clocksync change has a test.

**A163 note (2026-10-30, freeze; written before any validation seed (163) is simulated).** The frozen candidate is
**(i) with a degrees-of-freedom-adjusted host posterior** (`attenuation="df"`):
- R = 1 − tr S̄ / (tr cov(f̂) + tr S̄), where S̄ is recomputed with the anchors' ψ inflated by m / (m − 1) (m = number
  of host compartments).
- c_corrected = q / (qR + 1), with q = c / (1 − c): the exact errors-in-variables correction of ψ, not c / R.
- The interval is the donor-bootstrap percentile interval, with the host refitted and R recomputed in every replicate.
- Code: clocksync commit **9db4ed8** (`clocksync._core.host.reliability`, `corrected`;
  `Dataset.fit(attenuation="df")`). The default readout is unchanged.

Development (A162 seeds, 20 replicates per cell, 64 target cells; `output/dev/summary.csv`, `cells.csv`):

| readout | cells with coverage ≥ 0.90 | minimum coverage | max \|median bias\| |
|---|---|---|---|
| raw | 32 | 0.000 | 0.271 |
| "model" (fitted S) | 49 | 0.425 | 0.165 |
| "df" | 55 | 0.712 | 0.139 |

- The basic (reverse-percentile) interval of "df" covered worse: 25 of 64 cells.
- "df" is frozen.
- Expected weak spots, from development only: c_target 0.9 in n100_K8_bulk (0.71 at c_anchor 0.8) and n30_K4_bulk
  (0.78 at 0.6), and small-n bias at n30_K4 (−0.10 to −0.14 at c_target 0.6-0.9).
- The validation runs as registered: `compute_a163.py --mode validate --kind df --reps 200 --eq-reps 100`.

## A164 (2026-10-30, written before any matched statistic below is computed; roadmap 6.7; A47 and A96 part 1 were scored in the archive and are known: T2D module d* −0.043 [−0.091, 0.004], BMI H47.1 +0.025, p 0.79): do T2D and obesity leave the metabolic organs coupled? Matched design, power and equivalence

**Why:** A47 compared T2D cases with every non-diabetic donor of a cohort, and BMI as a slope. Cases differ from
controls in age, Hardy class and sex, and "no loss" was not a tested claim. A164 matches on the death circumstances and
tests equivalence, as A153 does for Figure 4.
**Data:** GTEx v8 clock-gene pseudobulks of the A47 Datasets: the `fast` (Hardy 1-2) and `vent` (Hardy 0) cohorts of
`organs/tissues.csv`, always kept separate, with the A47 anchors (all GTEx tissues of the cohort; a target anchor on its
leave-one-out host). Protected phenotypes (`common/gtex.py`): T2D, BMI, exact age, sex, Hardy class, ischaemic time.
Controlled access: aggregates only.
- **Module:** liver, subcutaneous and visceral adipose, pancreas, skeletal muscle (A47 MET).
- **Companion:** tibial nerve and the two skins (A47 COMP).
- A tissue enters a cohort's module only with ≥ 10 matched cases (fast pancreas, 5 cases, is out).

**Matching (within cohort, blind to clock data):**
- **T2D:** each case is matched to up to 3 non-diabetic controls, without replacement.
  - Exact on sex and Hardy class.
  - Nearest on age, with a caliper of 5 years, then on ischaemic time.
  - Unmatched cases are dropped and counted.
- **Obesity:** BMI ≥ 30 vs 18.5 ≤ BMI < 25, T2D donors excluded. Each obese donor is matched 1:1 to a normal-weight
  donor in the same way.
- **Balance** is reported before any clock readout: standardised mean differences of age, ischaemic time and BMI.

**Estimand:** per tissue, d* = c*_case − c*_control at the common host spread (`Dataset.contrast`, untruncated), then:
- the module mean d*;
- the module minus companion mean d*;
- pooled over the two cohorts by inverse variance.

The variance comes from a donor bootstrap that resamples matched sets (500 replicates), with the host refitted.
Permutations swap labels within matched sets (2,000).
**Tests:**
- **H164.1 (T2D, equivalence):** the 90 % interval of the pooled module d* lies within ±0.15 (the A153 bound).
- **H164.2 (obesity, equivalence):** the same for BMI ≥ 30 vs < 25.
- **H164.3 (effect):** the two-sided permutation p of the pooled module d*, Holm over T2D and obesity.
  - A Holm p < 0.05 is an effect.
  - It is "organ-specific" only if the module-minus-companion p, on the same permutations, is also < 0.05.

**Power, computed before the tests from the matched design only:**
- MDE = 2.49 × SD of the permutation null of the pooled module d* (80 % power, two-sided 0.05).
- The expected 90 % half-width = 1.645 × the bootstrap SE.

If the expected half-width exceeds 0.15, the equivalence test is stated "not powered" and reported without a verdict
word.
**Sensitivities (descriptive):**
- T2D additionally matched on BMI (caliper 3);
- without long-term steroid users;
- per cohort;
- per tissue;
- the A47 unmatched estimates side by side.

**Reading:**
- H164.1-2 pass → "T2D and obesity do not loosen the metabolic organs from the body clock" becomes a tested
  equivalence at |d*| < 0.15, on matched donors.
- An effect in H164.3 that is organ-specific → a disease positive.
- Not powered → the claim stays at its MDE.

**Code:** `paper/exploratory/A164_gtex_metabolic_matched/` (model: `results/figure_3/scripts/gtex_cues/compute_metabolic.py`).

## A165 (2026-10-30, written after A153 was scored and before any phase, amplitude or noise estimate of its contrasts was looked at; only their bootstrap SEs were read, to judge power; PI request 2026-10-30): do diseases shift, damp or blur the clocks of normal compartments? Phase, amplitude and noise readouts of the A153 contrasts

**Why:** A17 and A153 tested only the coupling c. c mixes the amplitude of the host-locked part (z2) with the
compartment's own noise around the host (ψ), and it is blind to phase. A153 found no compartment able to detect
|Δc| ≥ 0.2. Three disease mechanisms are separable in the same fits:
- **Internal desynchrony:** a compartment no longer on the same internal time as the body's other clocks. This
  extends the usual peripheral-clock meaning to cell types. Readout: the phase offset of the compartment against the
  host, arg z, case minus control (`offset_diff_h`).
- **Damping:** a smaller host-locked amplitude. Readout: log(z2_case / z2_control).
- **Loss of clock structure:** more clock-gene variance not locked to the host. Readout: log(ψ_case / ψ_control).
  Expected in malignant cells; in normal cells of diseased tissue it is the less expected and more interesting case.

The bootstrap SEs in the A153 output (fibroblast_stromal and endothelial rows, read before this note) give median
approximate MDEs of about 2.2 h for phase, about 2-fold for amplitude and about 3-fold for ψ. Phase is the
best-powered readout.
**Data and code:** the 13 estimable A153 targets, at the same tables, anchors, seeds, blocks, gates and
`anchor_readout="loo"`. The contrast is rerun with 1,000 permutations and 200 bootstrap draws (A153 counts). The
permutation p of the log z2 and log ψ ratios is computed on the same permutation draws as d* and the offset.
clocksync's `Dataset.contrast` gains the columns `p_log_z2` and `p_log_psi` for this; existing columns are unchanged.
Brain targets are out (as in A153).
**Primaries:** every estimable compartment of the fibroblast_stromal and endothelial families. These are the two
families eligible in either scope of the A153 family screen and clock-specific in A150 D150.1. Unlike A153, no MDE
bound is applied: no disease-effect size can be anticipated (PI), so the MDE is reported next to every row instead.
Per readout, a row is excluded when:
- **phase and z2:** either arm is at the boundary with z2 truncated, or z2 ≤ 0 in either arm;
- **ψ:** either arm is at the ψ floor.

Other families are reported descriptively with the same readouts, without tests.
**Tests:**
- **H165.1 (per compartment, effect).** For each primary row and each of the three readouts, take the two-sided
  permutation p. Holm runs over all primary rows × readouts together. An effect is a Holm p < 0.05.
  - **Specificity qualifier** (reported, not a gate), on the same permutations: the readout of the compartment minus
    the mean of the stromal-vascular companion set (leave-self-out, as in A153).
  - A phase effect also requires the shift not to be shared by the companion set, because offsets are already
    relative to the host.
- **H165.2 (pooled, per family, directional).** These are the two mechanisms with a predicted sign:
  - damping: pooled log z2 ratio < 0;
  - loss of clock structure: pooled log ψ ratio > 0.

  Per target, take the inverse-variance mean of its primary rows in the family, with SE equal to their mean bootstrap
  SE. Pool over targets with `meta(method="hk")`, one-sided at 0.05, Holm over the 4 tests (2 families × 2 readouts).
  This needs ≥ 3 targets per family, otherwise it is not estimable. Targets sharing a control arm are not independent;
  this is stated and not adjusted, as in A153.
- **Phase is not pooled.** Its sign is not predicted, and shifts of opposite sign across diseases would cancel.
**Power:** approximate MDE = 2.80 × the bootstrap SE per row and readout (80 %, two-sided 0.05), reported in hours,
fold-change in amplitude (exp(MDE/2)) and fold-change in ψ. A non-significant row is "not detected (MDE …)", never
"unchanged". No equivalence test is run: no margin is pre-specified for these readouts.
**Descriptive:**
- the gate clock share per arm (shared non-clock structure), case minus control, per target;
- the host spread ratio log(V_case / V_control);
- all readouts for the non-primary families.
**Reading:**
- An H165.1 phase effect that the companion does not share is internal desynchrony of that compartment in disease.
- An H165.2 ψ effect is a general loss of clock structure in normal stromal-vascular cells of diseased tissue.
- Nothing detected: the disease clocks of these compartments are not shown to shift by more than about 1-2 h, damp by
  more than about 2-fold, or blur by more than about 3-fold (per-row MDEs).
**Code:** `paper/exploratory/A165_disease_readouts/`; clocksync `Dataset.contrast` (p_log_z2, p_log_psi).

## A19 note (2026-10-30, roadmap 3.5; written before the rerun below): age effects rerun on the relearned anchors, with MDEs; H19.2 under the family rule

**Why:** the A19 part 2 re-read (2026-10-26, permutation null) predates the relearned anchors (2026-10-28) and reports
no power. H19.2 contrasts lymphocyte + myeloid with stromal-vascular compartments. Under the A153 family rule
(note of 2026-10-30), lymphocyte and myeloid are not clock-specific, so a locked-amplitude difference built on them has
no demonstrated clock content.
**Rerun:** `figure_4/scripts/resilience/compute_age.py`, unchanged:
- permutation null;
- the current rule-stamped tables and relearned learned anchors;
- `--jobs 4`.

The 2026-10-26 outputs are kept for comparison as `*_2026-10-26.csv`.
**Additions (descriptive, from the pooled table):**
- an approximate MDE = 2.80 × pooled SE (80 %, two-sided 0.05) for H19.1, H19.2 and every family secondary, in the
  estimand's units per decade;
- **H19.2 is labelled uninformative** whatever its p: its immune arm is not eligible under the family rule. It is
  reported, never read as "immune clocks do not age differently".
- The eligible-family readouts, fibroblast_stromal (all5) and endothelial (rowcheck3 sensitivity), are reported
  per estimand with their MDEs. These are secondaries, not tests.

**Reading:** the H19.1 verdict stands as registered and gains an MDE. "No age effect" is stated only down to that MDE.

## A166 (2026-10-30, written after a metadata-only scout (`exploratory/scouts_2026-10-30/metabric/SCOUT.md`: counts, gene coverage, no clock statistic) and before any clock readout of METABRIC or of TCGA BRCA within PAM50; PI request 2026-10-30): independent replication of the breast-tumour ER and proliferation coherence losses (METABRIC)

**Why:** in SCAN-B, ER-negative tumours had less coherent clocks (H127.2, T −0.343). After scoring, this held within
PAM50 subtype × proliferation tertile (A127 post hoc, T −0.268, p 0.0045). The post hoc result needs an independent,
pre-specified replication before it can be a paper result. TCGA BRCA is not independent for the ER effect: A97 first
saw it there, exploratorily. METABRIC (about 2,000 tumours, a different platform and country) is the confirmatory cohort.
**Data:**
- **METABRIC** (cBioPortal `brca_metabric`, ODbL 1.0; `$SCCC_DATA/metabric/`, PROVENANCE.md).
  - `data_mrna_illumina_microarray.txt` (Illumina HT-12 v3, log2 intensity), 1,980 tumours with expression.
  - ER = `ER_IHC` (Positive / Negative; missing out). PAM50 = `CLAUDIN_SUBTYPE` restricted to the A127 five (LumA,
    LumB, Her2, Basal, Normal); claudin-low and NC are out of the tests.
  - Stratifier: `COHORT` (5 levels).
- **TCGA BRCA** (secondary, not independent): recount3 primary tumours, one per patient (first aliquot), log CPM as
  the repo convention. ER from the GDC XML (`xml_breast_carcinoma_estrogen_receptor_status`), PAM50 from PanCanAtlas
  (`$SCCC_DATA/tcga/pancan_subtypes/`).

**Definitions (frozen from A127, `exploratory/A127_scanb/scripts/a127_lib.py`):**
- Y = ln 2 × log2 expression of the 14 clock genes.
- Proliferation = mean z of Hallmark G2M (clock genes removed); purity proxy = the A97 11-gene mean z.
- Both scores are z-scored within cohort (METABRIC: within COHORT); tertiles are by rank within the stratum named below.
- Coherence = A87 `decompose.components` gene coherence over a group of tumours, groups centred.

**Tests (METABRIC, 2,000 permutations each, one-sided for a loss, Holm over the three):**
- **H166.1 (ER, = H127.2):** coherence(ER−) − coherence(ER+) within COHORT × purity-proxy tertile, mean over strata
  with ≥ 10 tumours per arm; ER labels permuted within stratum.
- **H166.2 (ER within PAM50 × proliferation, = the A127 post hoc):** cells PAM50 × proliferation tertile (tertiles
  within PAM50), with ≥ 10 tumours per ER arm; T = mean over cells; ER labels permuted within cell.
- **H166.3 (proliferation, = H127.1):** strata PAM50 × purity-proxy tertile with ≥ 30 tumours; T = mean of
  coherence(high proliferation tertile) − coherence(low); tumour order permuted within stratum.

**Specificity companion (same statistic, observed labels):** 200 pseudo-clock 14-gene sets matched on mean and variance
in METABRIC (clock and G2M genes excluded; `wf_for` of A127). A pass is "clock-specific" only if T is below the 5th
percentile of the pseudo T. Otherwise it is "coherence loss not specific to the clock".
**Power, before the tests:** MDE = 2.49 × SD of each permutation null (one-sided 0.05, 80 %), reported next to the SCAN-B
effect sizes (H127.2 −0.343, post hoc −0.268, H127.1 −0.123).
**Sensitivities (descriptive):**
- ER by `ER_STATUS`;
- claudin-low included as its own stratum;
- without COHORT stratification;
- per COHORT;
- per-gene mean intensity against the array floor (4.56).
**TCGA BRCA (secondary, descriptive; not confirmatory for H166.1):** H166.2 and H166.3 with the same rules. Cells with
< 10 tumours per ER arm are dropped.
**Not registered here:** survival or endocrine-therapy response. These need PI approval, and endocrine response is
confounded by indication.
**Reading:**
- H166.1-2 pass, clock-specific → the ER effect, and its independence from subtype and proliferation, replicate in an
  independent cohort and platform. It is a pre-specified cancer positive (Figure 5).
- Fail at a powered design → the SCAN-B within-PAM50 result stays post hoc, and the effect is bounded at the MDE.
- Not specific → the wording changes to "general transcriptional decoherence".
**Code:** `paper/exploratory/A166_metabric/`.

## A167 (2026-10-30, written before any power or contrast statistic below is computed; NOT BLIND: the same four cohorts were scored in archive A38, A43, A44 and A96; PI request 2026-10-30): advanced liver fibrosis and the stellate and endothelial clocks, without living-donor controls

**Why:** the liver was not in A153/A165. With about 110 patients, the pooled open MASLD cohorts are the largest
single-cell disease set on disk, and stellate cells (fibroblast_stromal) and endothelial cells are the two families
eligible for clock readouts (A153 family screen, A150 D150.1). Stellate activation drives MASH fibrosis.
**Known before writing (archive):**
- A38 (GSE202379, SAF): hepatocyte coupling 0.68 → 0.54, p 0.195; endothelial 0.89 → 0.65, p 0.075.
- A43/A44: coupling falls in every class, anchors included, and the host spread narrows; the pooled effect is carried
  by the 20 GSE289173 living-donor controls and is near null without them. D2 (F3-4 vs F0-2, with the living donors)
  hepatocyte minus cholangiocyte −0.29 (p 0.007).
- The stellate and endothelial rows of the 2026-10-26 M3 output (`a44_contrast.csv`) were not read for this note.

**Data:** `data.load("masld_pooled")` (GSE202379, GSE289173, GSE244832, GSE256398; unit = patient). The 20 GSE289173
living-donor controls are excluded (procurement). Metadata read for this note (no clock statistic):

| | GSE202379 | GSE244832 | GSE256398 | GSE289173 (no living donors) | total |
|---|---|---|---|---|---|
| D2 case (F3-4) | 21 | 6 | 7 | 10 | 44 |
| D2 control (F0-2) | 26 | 9 | 9 | 18 | 62 |

Patients with the compartment (all studies, before exclusion): stellate 112, endothelial 124, hepatocyte 126.
**Design:** D2, F3-4 (case) vs F0-2 (control). Labels are permuted within study. Anchors as A43/A44
(endothelial, stellate, macrophage) with `anchor_readout="loo"`, so each anchor is read against the host of the
others. Pooled gate and host gate in each arm required (`require_host_gate_each_arm=True`); failing either makes the
contrast not estimable (uninformative, never "decoupled"). 1,000 permutations, 200 bootstrap draws; seed crc32("A167").
No A22 own-d* anchor drop (anchors are targets here). The A163 correction is not applied (not yet validated).
**Power step (run first, recorded before the contrast is run):** MDE on d* of stellate and endothelial by
`power.mde` with the A153 conversion (100 replicates, 200 permutations, within study). No MDE bound is applied
(as A165); the MDE is reported next to each verdict.
**Primaries:** liver:stellate and liver:endothelial.
**Tests:**
- **H167.1 (coupling):** d* (case minus control at a common host spread), two-sided permutation p, Holm over the 2
  primaries. Effect if Holm p < 0.05.
- **H167.2 (mechanism):** phase offset (`offset_diff_h`, internal desynchrony), log z2 ratio (damping) and log ψ
  ratio (loss of clock structure) of the 2 primaries, two-sided permutation p on the same draws, Holm over the 6.
  Same per-readout exclusions as A165.
- **Specificity qualifier (reported, not a gate), same permutations:** each readout of the primary minus the mean of
  all other estimable liver compartments (leave-self-out, "liver-wide companion"). A liver-wide loss (A44) would leave
  this near 0. An effect whose companion p ≥ 0.05 is labelled "liver-wide, not compartment-specific".
- A phase effect also requires the companion phase p < 0.05 (as A165).

**Secondaries (descriptive, no test):**
- D1 (MASH vs no MASLD / steatosis only), living donors excluded, same readouts;
- leave-one-study-out of H167.1-2;
- hepatocyte, cholangiocyte, macrophage, lymphocyte rows (not eligible families);
- host spread ratio log(V_case / V_control) and gate clock share per arm.

**Power statement:** approximate MDE = 2.80 × bootstrap SE for H167.2 readouts (units as A165). A non-significant
row is "not detected (MDE …)", never "unchanged".
**Reading:**
- H167.1 or H167.2 effect, compartment-specific: advanced fibrosis changes the stellate (or endothelial) clock beyond
  the rest of the liver. Reported as a non-blind re-read of seen cohorts, to be replicated.
- Effect but liver-wide: consistent with A44; disease cannot be separated from procurement and host-spread change.
- Nothing detected: stated down to the MDEs.

**Code:** `paper/exploratory/A167_liver_fibrosis/`.

**A167 amendment 1 (2026-10-30, PI; written while the power step for stellate / endothelial was running, before any
power or contrast output was read): the hepatocyte is a third primary.**
- **Why:** the A153 family screen marks "other" not eligible on 2 compartments from 2 atlases (thin data, not a
  negative); D150.1 put "other" above the pseudo-clock q95; hepatocytes are the liver clock (A129: hepatocyte Bmal1 KO
  decouples the liver, d* −0.287).
- **Even less blind** than the stellate / endothelial rows: A38-A44 targeted hepatocytes (D1 without living donors
  −0.14, p 0.066; D2 with living donors hepatocyte minus cholangiocyte −0.29, p 0.007).
- **Changes:** primaries liver:stellate, liver:endothelial, liver:hepatocyte. H167.1 Holm over 3; H167.2 Holm over 9.
  The hepatocyte gets the same power step (MDE on d*) and the same liver-wide companion. Anchors unchanged (the
  hepatocyte is not an anchor).
- **Caveat stated in the reading:** hepatocytes dominate the nuclei (median about 2,000 per patient vs 166 stellate),
  so ambient hepatocyte RNA may enter the other pseudobulks (D150.2 flags the epithelial family); a hepatocyte effect
  shared by the liver-wide companion is read as liver-wide.

**A167 power step (2026-10-30, recorded before the contrast is run):** D2 without living donors, 62 control / 44 case,
`power.mde` 100 replicates × 200 permutations within study, A153 conversion.

| primary | c0 (pooled) | s* | MDE on d* |
|---|---|---|---|
| liver:stellate | 0.753 | 0.528 | 0.293 |
| liver:endothelial | 0.715 | 0.400 | 0.429 |
| liver:hepatocyte | 0.561 | 0.459 | 0.349 |

All three exceed 0.2; as registered, no bound is applied and each verdict carries its MDE. A null is "not detected
(MDE …)". <!-- src: paper/exploratory/A167_liver_fibrosis/output/power.csv -->

## A168 (2026-10-30, written after a metadata-only scout (`exploratory/scouts_2026-10-30/endocrine/SCOUT.md`: sample tables, probe existence; no clock statistic or clock-gene value) and before any clock readout of these series; PI approval 2026-10-30): does blocking oestrogen signalling reduce breast-tumour clock coherence? Neoadjuvant aromatase-inhibitor pairs

**Why:** ER-negative tumours have less coherent clocks than ER-positive ones (H127.2; replicated in METABRIC, H166.1),
also within PAM50 × proliferation (A127 post hoc; H166.2, Holm 0.047). Both are between-tumour comparisons. If ER
signalling sustains the tumour clock, two weeks of oestrogen deprivation should reduce coherence in the same tumours.
This turns the association into an intervention.
**Data** (`$SCCC_DATA/scout_endocrine/`, PROVENANCE.md; GEO series matrices as deposited; all tumours ER+):
- **POETIC, GSE105777** (Illumina HT-12 v4, RSN-normalised): 157 patients with a baseline core biopsy (B) and a
  surgical sample after 2 weeks of an aromatase inhibitor (S); 56 untreated patients with a diagnostic core (D) and a
  surgical sample (S). Pairing from the sample titles.
- **GSE20181** (Edinburgh letrozole, Affymetrix U133A): 58 patients with pre and 10-14-day samples (A, B); day-90 (C)
  samples are descriptive only (submitted 3.5 years later: batch). Contains GSE5462 (never counted twice).
- **GSE87411** (Z1031B, Agilent two-colour log ratios to a universal reference; baseline vs 2-4 weeks), descriptive.

**Definitions (frozen from A127 / A166):**
- Expression on the log2 scale (log2-transformed if the deposited values are linear, i.e. maximum > 100); probes
  averaged per gene symbol; Y = ln 2 × log2 expression of the 14 A127 clock genes. GSE20181 has no CIART probe and uses
  the other 13 (its pseudo sets have 13 genes).
- Coherence = A87 `decompose.components` gene coherence over a group of samples, genes centred within the group
  (`a127_lib.coh`). Groups: each arm × time point.
- Proliferation = mean z of Hallmark G2M (clock genes removed); purity proxy = the A97 11-gene mean z; z-scored within
  series.

**Tests (one-sided for a loss, 2,000 permutations each, Holm over H168.1-3):**
- **H168.1 (primary, POETIC difference-in-differences):** T = [coh(AI, S) − coh(AI, B)] − [coh(untreated, S) −
  coh(untreated, D)] < 0. Null: patients reassigned between the two arms at random (sizes 157 / 56 kept), each
  patient keeping both samples. This removes the core-biopsy-vs-surgery sampling effect, which is shared by both arms.
- **H168.2 (POETIC, treated change):** coh(AI, S) − coh(AI, B) < 0. Null: the two time labels swapped at random within
  each treated patient. A pass without H168.1 is read as a sampling effect, not a treatment effect.
- **H168.3 (replication, GSE20181):** coh(day 10-14) − coh(pre) < 0 over the 58 pairs; time labels swapped within
  patient. No untreated arm: supports H168.1 in direction only.

**Specificity companion (same statistic, observed labels):** 200 pseudo-clock sets per series matched on mean and SD
(`a166_lib.pseudo_sets`; pool excludes clock-like, G2M, purity-proxy genes and, new here, the Hallmark oestrogen
response early and late genes, which the treatment moves directly). A pass is "clock-specific" only if T is below the
5th percentile of the pseudo T; otherwise "coherence loss not specific to the clock".
**Power, before the tests:** MDE = 2.49 × SD of each permutation null (one-sided 0.05, 80 %), reported next to the
METABRIC ER effect beyond its null mean (about −0.165) and the SCAN-B one (about −0.315).
**Sensitivities (descriptive):**
- without PER1 and without the immediate-early-responsive PER1, NR1D1 pair (sampling-sensitive genes);
- Y residualised on proliferation and purity proxy within arm × time point (proliferation falls on AI, and in A127
  higher proliferation means lower coherence, so the proliferation route predicts a coherence *rise*);
- POETIC within thirds of the change in proliferation score (a response proxy: no Ki67 in GEO);
- GSE20181 within responders and non-responders (title labels), and pre vs day 90;
- GSE87411 baseline vs 2-4 weeks (all pairs, with and without the "re-analyzed" baselines).
- Change in proliferation, purity proxy and the oestrogen-response score per arm (to show the treatment acted).

**Reading:**
- H168.1 pass, clock-specific → oestrogen deprivation reduces tumour clock coherence beyond the sampling effect: ER
  signalling sustains the clock in ER+ tumours, the interventional counterpart of the ER− loss. With H168.3, replicated
  across platform and drug.
- H168.1 fail at a powered design → the ER− loss is not reproduced by short-term ER blockade: it is a stable property
  of the tumour type (or needs longer than 2 weeks), bounded at the MDE.
- H168.2 pass, H168.1 fail → a sampling effect.
- Not specific → "general transcriptional decoherence under treatment".

**Not registered here:** outcome or response prediction (PI).
**Code:** `paper/exploratory/A168_endocrine/`.

**A168 power step (2026-10-31, recorded before any observed statistic is computed):** 2,000 permutations per null;
MDE = 2.49 × null SD (`run_a168.py --stage power`).

| test | patients | null mean | null SD | MDE | vs METABRIC (0.165) / SCAN-B (0.315) |
|---|---|---|---|---|---|
| H168.1 POETIC difference-in-differences | 157 AI / 56 untreated | −0.044 | 0.173 | 0.431 | underpowered for both |
| H168.2 POETIC treated change | 157 | 0.000 | 0.079 | 0.196 | powered for SCAN-B, not METABRIC |
| H168.3 GSE20181 | 58 | −0.005 | 0.203 | 0.504 | underpowered for both |

The 56 untreated pairs limit H168.1. Data used: POETIC all 157 + 56 pairs; GSE20181 58 pairs, 13 clock genes (no
CIART; linear MAS5 values log2-transformed); GSE87411 81 BL→2-4 wk pairs (97 with the "re-analyzed" baselines).
Implementation choices: `paper/exploratory/A168_endocrine/README.md`. <!-- src: paper/exploratory/A168_endocrine/output/{power,sample_counts}.csv -->

## A169 (2026-10-31, written after A168 was scored and before any clock readout of GSE111563, GSE59515 or GSE55374; GSE20181 day 90 is already seen in A168: −0.143, p 0.32; PI request 2026-10-31): does tumour clock coherence fall progressively with the duration of oestrogen deprivation?

**Why:** A168 found no coherence loss after 2 weeks of aromatase inhibition. The ER− loss may need longer deprivation
to develop (tumours adapt, clones are selected). Edinburgh letrozole series sample the same tumours up to and beyond
3 months.
**Data** (`$SCCC_DATA/scout_endocrine/`, GEO series matrices; all ER+, neoadjuvant letrozole; each series was
submitted in one batch):
- **GSE111563** (HT-12 v4, probe-filtered matrix; 12 clock genes, no NR1D1 or CIART): time windows tp1 (< 0 d,
  pre), tp2 (0-120 d), tp4 (> 120 d). 37 patients with tp1 and tp4; 24 with all three.
- **GSE59515** (HT-12 v4, 14 genes): pre, 2 weeks, 3 months; 25 patients with all three.
- **GSE55374** (HT-12 v4): 13 patients pre and 3 months (10 with 2 weeks); descriptive.
- Patient overlap between the three Edinburgh series is unresolved (SCOUT.md): they are never pooled, and H169.1-2 are
  not independent replications of each other.

**Definitions:** as A168 (`A168_endocrine/scripts/a168_lib.py`: log2, probes averaged per symbol, coherence
`a127_lib.coh` per time-point group, scores z-scored within series).
**Tests (one-sided for a loss, 2,000 permutations, time labels permuted within patient; Holm over the two):**
- **H169.1 (GSE111563):** coh(tp4, > 120 d) − coh(tp1) over the 37 patients with both.
- **H169.2 (GSE59515):** coh(3 months) − coh(pre) over the 25 patients.
- **Progression (descriptive, both series):** coherence at each time point on the patients with all three, and the
  trend slope over the time-point rank (null: the three labels permuted within patient), one-sided p reported.

**Specificity companion:** 200 pseudo-clock sets per series, as A168 (oestrogen-response genes excluded from the pool).
"Clock-specific" only if T is below the 5th pseudo percentile.
**Power, before the tests:** MDE = 2.49 × null SD, next to the METABRIC (0.165) and SCAN-B (0.315) ER effects.
**Sensitivities (descriptive):** without PER1 (and NR1D1 where present); Y residualised on proliferation and purity
proxy within time point; responders vs non-responders (GSE59515 labels); the change in proliferation and
oestrogen-response scores per time point (to show the treatment acts throughout).
**Reading:**
- Pass, clock-specific → coherence falls with prolonged oestrogen deprivation: the ER− state can be approached by
  long deprivation, though not by 2 weeks (A168).
- Fail at these sizes → no progressive loss down to the MDE; with A168, the ER− loss is a stable tumour-type
  property. Small series: a fail is weak evidence; MDEs are stated.

**Code:** `paper/exploratory/A169_endocrine_duration/` (reuses the A168 library).

**A169 power step (2026-10-31, recorded before any observed statistic):** 2,000 permutations; MDE = 2.49 × null SD.
H169.1 GSE111563 (37 patients, 12 genes): null SD 0.263, MDE 0.653 (underpowered for METABRIC and SCAN-B).
H169.2 GSE59515 (25 patients, 14 genes): null SD 0.124, MDE 0.307 (about the SCAN-B size).
<!-- src: paper/exploratory/A169_endocrine_duration/output/power.csv -->

## A170 (2026-10-31, written after a metadata-only scout (`exploratory/scouts_2026-10-31/sc_breast_er/SCOUT.md`, parts 1-2: cell counts per patient and class, ER labels; no expression value read except the two archive rows noted in part 1, which are not from Wu or Pal) and before any expression of Wu 2021 or Pal 2021 is read; PI request 2026-10-31): is the ER− loss of breast-tumour clock coherence carried by the malignant cells themselves?

**Why:** ER− breast tumours have lower bulk clock-gene coherence (A127 SCAN-B, A166 METABRIC, clock-specific), and
oestrogen deprivation does not reproduce it (A168, A169). Bulk cannot say whether the loss sits in the malignant cells
or in the microenvironment. Single-cell malignant pseudobulks can.
**Data** (`$SCCC_DATA/sc_breast_er/`):
- **Wu 2021** (GSE176078; raw counts, author cell types `celltype_major`): malignant = "Cancer Epithelial". ER from the
  paper's Supplementary Table 1 (`wu_reconciled.csv`, `final_ER`).
- **Pal 2021** (GSE161529; 10x matrices of the figshare deposit 17058077; malignant = the authors' `Tum` clusters of the
  ERTotal, TNBC and HER2 objects, cell names `<sample>_<barcode>`). ER+ = the ER objects; ER− = TNBC (BRCA1 carriers
  included). HER2+ tumours (ER not reported) are not used.
- **Units:** one tumour per patient with ≥ 50 malignant cells. Excluded: Pal ER-0319 (ER-low), patient 0114 (one tumour
  in each arm: both dropped). Expected from the scout: ER+ 22 (Wu 11, Pal 11), ER− 16 (Wu 9, Pal 7).

**Definitions:** malignant pseudobulk = summed raw counts of the patient's malignant cells; y = log((n + 0.5) / depth)
of the 14 A127 clock genes present (`decompose.logfrac`, depth = total counts of the pseudobulk), centred per gene
within cohort. Coherence = `a127_lib.coh` (A97 / A127 statistic) over the tumours of an arm.
**Test:**
- **H170.1 (primary, one-sided for a loss):** T = coh(ER−) − coh(ER+) on the pooled Wu + Pal malignant pseudobulks;
  2,000 permutations of ER labels within cohort. Pass: p < 0.05.
- **Specificity companion:** 200 pseudo-clock sets (A168 construction: mean / SD matched in the pooled malignant
  pseudobulks; pool without clock-like, G2M, purity-proxy, oestrogen-response or clock genes; genes with ≥ 1 count per
  1e5 in every unit). "Clock-specific" only if T is below the 5th pseudo percentile.
**Power, before the test:** null SD of T from the permutations, MDE = 2.49 × SD, next to the bulk ER effects (METABRIC
−0.165 beyond the null, SCAN-B −0.315).
**Descriptive (no test):** T per cohort; whole-tumour pseudobulk (all cells) T in the same patients, to show whether the
bulk effect is visible at this size; Wu CAF and endothelial coherence by ER (patients with ≥ 50 cells of the class);
without CID3963 (borderline ER-low); without the 4 Pal BRCA1 TNBC.
**Reading:**
- Pass, clock-specific → the ER− coherence loss is carried by the malignant cells (cell-intrinsic).
- Fail with the whole-tumour T also near 0 → the cohort is too small to see the bulk effect; uninformative.
- Fail with the whole-tumour T negative → the loss is not in the malignant cells: the microenvironment or composition
  carries it (descriptive only).
**Code:** `paper/exploratory/A170_sc_breast_er/`.

**A170 power step (2026-10-31, recorded before any observed statistic):** 38 malignant pseudobulks (ER+ 22: Wu 11, Pal 11;
ER− 16: Wu 9, Pal 7), 14 clock genes, 2,000 within-cohort permutations: null mean +0.022, SD 0.155, MDE 0.385 (above the
SCAN-B effect 0.315: a fail cannot exclude a bulk-sized loss).
<!-- src: paper/exploratory/A170_sc_breast_er/output/power.csv -->

## A171 (2026-10-31, written after a metadata-only scout (`exploratory/scouts_2026-10-31/luad_lihc/SCOUT.md`: sample tables, platform gene coverage, whole-matrix scale; no clock statistic or clock-gene value) and before any clock readout of the series below; TCGA LUAD / LIHC are already seen (A97); PI request 2026-10-31): does tumour clock coherence fall with proliferation in lung adenocarcinoma and hepatocellular carcinoma? Independent replication of the A97 per-cancer losses

**Why:** in TCGA, clock-gene coherence falls from the low to the high proliferation tertile most strongly in LUAD
(−0.515, 516 tumours) and LIHC (−0.359, 371), more than in BRCA (−0.163), where the loss replicated (A127, A166). These
per-cancer values were descriptive, unstratified and are already seen; they need independent, pre-specified cohorts.
**Data** (`$SCCC_DATA/scout_luad_lihc/`, GEO series matrices and platform annotations; primary tumours only):
- **H171.1 LUAD GSE72094** (Moffitt; GPL15048 custom Affymetrix, log2 as deposited; 442 tumours; 14 clock genes).
- **H171.2 LUAD GSE68465** (Director's Challenge; GPL96, MAS5 linear → log2; 443 adenocarcinomas; 13 clock genes, no
  CIART; 4 sites).
- **H171.3 LIHC GSE14520** (GPL3921 tumours only, 225; log2 as deposited; 13 clock genes, no CIART).
- **H171.4 LIHC GSE36376** (GPL10558, log2 as deposited; 240 tumours; 14 clock genes).
- **Descriptive:** GSE31210 (LUAD, 226, stage I-II; MAS5 → log2), GSE41271 (adenocarcinomas only), GSE63898 (LIHC,
  228); TCGA LUAD and LIHC with the same statistic (not independent).
**Definitions** (A166 / A168 conventions, `A168_endocrine/scripts/prep_a168.py`): linear matrices (max > 100) → log2
with the smallest positive value as floor; probes averaged per symbol (a multi-symbol probe counts for each);
Y = ln 2 × log2 expression of the clock genes present (the waveform restricted to them); proliferation = Hallmark G2M
mean z, purity proxy = the A97 11-gene mean z, both z-scored within series (GSE68465: within site, and Y centred per
gene within site before any coherence). Coherence = `a127_lib.coh`.
**Tests (one-sided for a loss, 2,000 permutations each, Holm over the four):** the H166.3 statistic. Strata = purity-proxy
tertiles of the series (≥ 30 tumours); within each stratum, proliferation tertiles by rank; T = mean over strata of
coh(high) − coh(low); null: proliferation scores permuted within stratum.
**Specificity companion:** 200 pseudo-clock sets per series (A168 construction; pool without clock-like, G2M,
purity-proxy or clock genes). "Clock-specific" only if T is below the 5th pseudo percentile.
**Power, before the tests:** MDE = 2.49 × null SD per test, next to the TCGA values with the same statistic (computed in
the power stage; seen data).
**Descriptive (no Holm, no claim):**
- the four statistics in GSE31210, GSE41271, GSE63898 and TCGA LUAD / LIHC;
- LUAD GSE72094: coherence by KRAS, EGFR, STK11, TP53 mutation and smoking (ever / never): coh(group) − coh(rest)
  within purity tertile, and within proliferation tertile, with permutation p;
- LIHC: tumour vs paired non-tumour coherence in GSE14520 (GPL3921, labels swapped within patient) and tumour vs
  adjacent tissue in GSE36376 (unpaired permutation);
- each test without PER1 and NR1D1.
**Reading:**
- A cancer's tests pass, clock-specific → the proliferation coherence loss replicates in that cancer: with breast, a
  pan-cancer property of proliferating tumours (Figure 5).
- Pass, not specific → proliferating tumours lose coherence in co-regulated gene sets in general.
- Fail at a powered design → the TCGA value of that cancer is not confirmed; bounded at the MDE.
**Code:** `paper/exploratory/A171_luad_lihc/`.

**A171 power step (2026-10-31, recorded before any statistic of the four test series):** 2,000 permutations; three
purity strata in each series. MDE (2.49 × null SD): H171.1 GSE72094 (442 tumours, 14 genes) 0.167; H171.2 GSE68465 (443,
13) 0.304; H171.3 GSE14520 (225, 13) 0.318; H171.4 GSE36376 (240, 14) 0.363. TCGA with the same statistic (seen data):
LUAD T −0.361 (516, p 0.0015), LIHC −0.512 (371, p 0.0005). All four tests are powered for the TCGA effect of their
cancer, H171.4 narrowly for LUAD-sized effects only.
<!-- src: paper/exploratory/A171_luad_lihc/output/{power,tcga}.csv -->

## A172 (2026-10-31, written after A171 was scored and before any of the TCGA projects below is downloaded; the 11 A97 projects are seen and excluded; PI request 2026-10-31 "work on the different cancer types"): is the proliferation clock-coherence loss a general property of tumours? The unseen TCGA projects

**Why:** proliferating tumours lose clock coherence in breast (A127, A166) and lung adenocarcinoma (A171), not in HCC
(A171). In the 11 TCGA projects of A97 (seen; `A171_luad_lihc/output/tcga_all.csv`, post hoc) the loss is clear in
LUAD, LIHC and KIRP and absent in BLCA, KIRC and THCA. The remaining TCGA projects have never been read and form an
independent pan-cancer test with the same pipeline.
**Data:** recount3 TCGA (G029; same mirror, scaling and loader as `paper/common/datasets/tcga_paired/build.py`) for every
TCGA project not among the 11 (ACC, CESC, CHOL, DLBC, ESCA, GBM, HNSC, KICH, LAML, LGG, LUSC, MESO, PAAD, PCPG, READ,
SARC, SKCM, STAD, TGCT, THYM, UCS, UVM). **Eligible:** ≥ 150 patients with a primary tumour (sample type "Primary Tumor";
LAML "Primary Blood Derived Cancer - Peripheral Blood"; one per patient, lowest barcode, any histology), counted from the
metadata before any expression is read. Expected from GDC: CESC, ESCA, GBM, HNSC, LAML, LGG, LUSC, PAAD, PCPG, SARC,
STAD, TGCT.
**Definitions (= A171 TCGA arm, `A171_luad_lihc/scripts/run_a171.py` `tcga()`):** Y = log((n + 0.5) / depth) of the 14
clock genes; proliferation = Hallmark G2M mean z, purity proxy = A97 11 genes mean z of log2 CPM, z-scored within project.
Statistic per project = H166.3 / H171 (within purity-proxy tertiles of ≥ 30 tumours, coh(high proliferation tertile) −
coh(low)).
**Tests (one-sided for a loss, 2,000 permutations; proliferation permuted within stratum):**
- **H172.1 (primary):** T = unweighted mean of the per-project statistics over the eligible projects; null = the mean of
  the per-project null draws of the same index.
- **H172.2:** per eligible project, Holm over the eligible projects.
**Specificity companion:** per project 200 pseudo-clock sets (A168 construction on log2 CPM; pool without clock-like, G2M,
purity-proxy or clock genes); pooled pseudo T = mean over projects of the k-th set's T. "Clock-specific" only if T is
below the 5th pseudo percentile (H172.1; per project for passing H172.2 rows).
**Power, before the tests:** MDE = 2.49 × null SD for H172.1 and each project.
**Descriptive:** H172.1 over carcinomas only (CESC, ESCA, HNSC, LUSC, PAAD, STAD); without LAML (no solid stroma: the
purity proxy is not interpretable); each project without PER1 and NR1D1.
**Reading:**
- H172.1 passes, clock-specific → the proliferation coherence loss generalises to unseen tumour types; the per-project
  rows say where.
- Fail at a powered design → the loss is cancer-type specific (breast, LUAD), not a general tumour property.
**Code:** `paper/exploratory/A172_tcga_unseen/`.

**A172 power step (2026-10-31, before the tests; `A172_tcga_unseen/output/{eligibility,power}.csv`):** 12 eligible
projects: CESC 304, ESCA 184, GBM 155, HNSC 502, LGG 514, LUSC 501, PAAD 178, PCPG 179, READ 166, SARC 259, STAD 415, TGCT
150. Not eligible: LAML (126 primary blood samples), and ACC, CHOL, DLBC, KICH, MESO, SKCM (103 primaries; most SKCM
samples are metastatic), THYM, UCS, UVM. READ was not among the expected projects; LAML drops out, so "without LAML" =
H172.1. MDE (80 %): H172.1 0.087, carcinomas only 0.122; per project 0.175 (LGG) to 0.370 (GBM): CESC 0.259, ESCA 0.278,
HNSC 0.287, LUSC 0.317, PAAD 0.282, PCPG 0.340, READ 0.318, SARC 0.330, STAD 0.357, TGCT 0.286. For scale, the seen TCGA
effects are −0.51 to +0.15 (`A171_luad_lihc/output/tcga_all.csv`); the breast / LUAD effects are −0.17 to −0.36.

## A173 (2026-10-31, written after A172 was scored and the gastric / lung-histology scout (metadata only), before any expression value of the series below is read; GSE41271 adenocarcinomas were read in A171 (descriptive, −0.101), its squamous tumours were not; PI request 2026-10-31 "work on the different cancer types"): does the proliferation clock-coherence loss replicate in gastric cancer, and does it follow histology in lung?

**Why:** A172 (unseen TCGA) failed pan-cancer, with STAD the strongest project (−0.305, 2nd pseudo percentile, Holm
0.18) and the carcinomas leaning negative; LUSC was positive (+0.21) while LUAD is the strongest seen loss (−0.36) and
replicated in A171. Two follow-ups: (A) is the gastric loss real in independent cohorts; (B) is the lung loss specific to
adenocarcinoma.
**Data (scout `scouts_2026-10-31/gastric_lung_histology/`; parsed with the A171 parser, probes averaged per symbol):**
- Gastric: GSE62254 (ACRG, 300 tumours, GPL570; deposited values are log10-like (0.57-4.41) and are multiplied by log2 10),
  GSE84437 (Yonsei, 483, GPL6947, linear → log2; NaN cells as in A171 `finish`), GSE15459 (Singapore, 200 arrays, GPL570,
  linear → log2; all 200 arrays, the authors' 8 exclusions are not identifiable in GEO).
- Lung: adenocarcinoma and squamous tumours of GSE37745 (106 / 66, GPL570), GSE41271 (183 / 80, GPL6884), GSE30219 (85 /
  61, GPL570; ADC / SQC labels; other histologies and normals dropped).
**Scores:** G2M and purity proxy (A171 definitions) z-scored over the tumours of the analysis (gastric: series; lung: the
ADC + SQC tumours of the series).
**Tests (one-sided, 2,000 permutations):**
- **H173.1-3 (gastric):** the H171 statistic (purity-proxy tertile strata ≥ 30, coh(high proliferation tertile) − coh(low);
  proliferation permuted within stratum) < 0 in GSE62254, GSE84437, GSE15459; Holm over the three.
- **H173.4 (lung histology):** D = mean over the three lung series of (T_ADC − T_SQC) < 0, where T_h is the H171 statistic
  within histology h on purity-proxy halves (median split within histology; tertile strata would leave squamous strata
  below 30). Null: proliferation permuted within histology × half; D from the same permutation. Single test (α 0.05).
**Specificity companion:** 200 pseudo-clock sets per series (A168 construction; pool without clock-like, G2M, purity or
clock genes); for H173.4 the same sets give pseudo D. "Clock-specific" only if the observed value is below the 5th pseudo
percentile; evaluated for passing tests.
**Power, before the tests:** MDE = 2.49 × null SD per test.
**Descriptive:** T_ADC and T_SQC per lung series and their means; gastric GSE26253 (FFPE DASL, 432, PER3 missing); each
test without PER1 and NR1D1; gastric pooled mean of H173.1-3.
**Reading:**
- H173.1-3 pass (≥ 2 of 3) → the gastric loss replicates: a third cancer type after breast and lung adenocarcinoma.
- H173.4 passes → in lung the loss is a property of adenocarcinoma, not of the organ.
- Failures at a powered design bound the claim to breast and lung adenocarcinoma.
**Code:** `paper/exploratory/A173_gastric_lung_histology/`.

**A173 power step (2026-10-31, before the tests; `A173_gastric_lung_histology/output/power.csv`):** gastric strata
GSE62254 100/100/100, GSE84437 161 ×3, GSE15459 67/66/67; lung halves (ADC; SQC) GSE37745 53/53; 33/33, GSE41271
92/91; 40/40, GSE30219 43/42; 31/30. All series have 14 clock genes (GSE26253 13, no PER3). MDE (80 %): H173.1 0.224,
H173.2 0.229, H173.3 0.332 (TCGA STAD −0.305); H173.4 0.355 (TCGA LUAD − LUSC = −0.361 − 0.212 = −0.57); per-series D
0.60-0.63, per-histology T 0.41-0.49 (descriptive only).

## A174 (2026-11-01, written after the pseudobulk build (`A174_v3_heart/output/cells_per_class.csv`: cell counts only, no clock-gene value read) and before any gate, host or AM/PM readout on current code; the archive result is known (V3), so this is a **confirmation on current code**, not a new test; HANDOFF_2026-10-31 T1): does the V3 heart-surgery validation hold on the current estimator?

**Why:** V3 (archive, pre-A17 estimator: `scCircadianCancer/legacy/scripts/v3_heart.py`, `outputs/v3_heart/checks.csv`) is the
only true-time human validation that has not been ported: gate p 0.0005 (minimum at 2,000 permutations), PM − AM host
displacement p 0.0004 (each surgical arm 0.032 / 0.033), displacement angle 19.3 h against 18 h predicted for forward time.
**Data:** GSE270258 (Jia 2024), left-ventricle snRNA-seq of 16 cardiac-surgery patients, 8 operated 06-10 h (AM) and 8
14-18 h (PM), 4 + 4 in each arm (NCR, ACR); window only, no per-patient time. No author labels: cell classes from the
archive marker rule, frozen in `build_v3.py` (no clock gene among the markers), floor 20 cells.
**Design (fixed now):**
- Compartments: non-immune classes with ≥ 12 patients: cardiomyocyte (16), fibroblast (16), endothelial (16), neuron (15).
- Anchors: the learned A17 heart anchors present here (`common/anchors_a17.csv`, heart_global: cardiomyocyte, fibroblast,
  neuron; mural is learned but present in 4 patients only); clocksync c266e20 defaults (`anchor_readout="loo"`).
- Host state: `Dataset.host()` on these anchors; reference anchor fibroblast (as in the archive).
**Predictions (restated from the archive):**
- **V174.1** the gate passes (pairing and acrophase p < 0.05, 2,000 permutations) on the 16 patients.
- **V174.2** the PM − AM mean host-state displacement is longer than under window-label permutation (5,000, archive
  statistic), p < 0.05.
- **V174.3** its angle is within ±3 h of 18 h (forward time; 6 h if reversed).
**Sensitivity (descriptive):** each surgical arm separately (8 patients, as in the archive); archive floor of 50 cells; the
fixed stromal-vascular anchors available here (fibroblast, endothelial; 2 anchors, so only non-anchor targets are read);
myeloid added (12 patients).
**Reading:** V174.1-3 pass → V3 enters the true-time panel (T3) on current code. A failure is reported as such against the
archive values; no further amendment (stop rule).
**Code:** `paper/exploratory/A174_v3_heart/`.

## A175 (2026-11-01, written after a counts-only scout (pairs per cancer and per ER status; no pairing statistic by stratum computed) and before any stratified pairing readout; TCGA pairs and per-cancer coherence are seen (A97, H138.4, H138.8), so every test here is **secondary**; HANDOFF_2026-10-31 C1): does the tumour follow its own normal tissue less when it is proliferating or ER−? Readout level **L2** (tumour to the patient's normal tissue as host)

**Why:** the coherence loss in proliferating and ER− tumours (L1; A127, A166, A171, A173) does not use host referencing.
H138.4 showed that tumours follow their own patient's normal tissue (0.365 vs pairing null 0.122). C1 asks whether that
following is weaker in the tumour states that lose coherence.
**Readings stated in advance (both publishable):** lower pairing in high-proliferation (or ER−) tumours → the L1 loss
comes with a loss of host following. Pairing not lower, with an MDE at or below the TCGA H138.4 margin over its null
(0.24), → coherence is lost inside a tumour that still tracks the person.
**Data:** the H138.4 set (517 tumour-normal pairs, 10 gated cancers, `A138_cancer_estimability`, states and A52 gate
unchanged; the script first reproduces H138.4 = 0.3647).
**Strata:**
- Proliferation: G2M score of the tumour sample (A97 definition, `figure_5/scripts/cancer/_shared.py`: mean z of
  HALLMARK_G2M_CHECKPOINT without clock genes, on `tcga_full_counts`), tertiles by rank within cancer; high vs low
  tertile, cancers with ≥ 10 pairs per arm (expected to drop BLCA and UCEC).
- ER (BRCA only): `xml_breast_carcinoma_estrogen_receptor_status` (Negative vs Positive; scout: 21 vs 78 pairs).
**Tests (one-sided, label permutation within cancer, 10,000, p = (1 + #{D_null ≥ D}) / 10,001):**
- **H175.1:** D_prolif = pooled pairing(low tertile) − pooled pairing(high tertile) > 0 (pairs-weighted pooling, as
  H138.8).
- **H175.2:** D_ER = pairing(ER+) − pairing(ER−) > 0 in BRCA. Run only if its MDE ≤ 0.3; otherwise "not estimable at
  this n" (stop rule).
- Holm over the tests that run.
**Power, before the tests:** MDE = q95 + 0.84 × SD of the permutation null of D (labels only; the observed D is not
evaluated in this step).
**Specificity companions (same permutations):** D on 200 expression-matched pseudo-clock 14-gene sets per cancer (A52
rule) and on the 200 responsive pseudo sets (H138.4-R); a loss is "clock-specific" only if D exceeds the 95th pseudo
percentile of both.
**Descriptive:** D per cancer; proliferation residualised on ABSOLUTE purity within cancer (475 pairs); BRCA PAM50
basal + HER2 vs luminal; mid tertile.
**Reading:** see above; a failed H175.1 with MDE > 0.24 is "not resolved at this n", not "unchanged".
**Replication (primary, registered separately after this result):** GSE14520 HCC (213 pairs), GSE66229 gastric (98 paired
patients, scout header only), any CPTAC LUAD/STAD; MDE step first, stop if MDE exceeds the TCGA D.
**Code:** `paper/exploratory/A175_c1_pairing_strata/`.

## A176 (2026-11-01, written before any clock-plane or off-plane variance is read in these cohorts (the A127/A166/A171/A173 outputs hold coherence only; scout `A127` off-plane descriptive is a different per-tumour measure); the coherence contrasts are seen, so the coherence column here is a reproduction; HANDOFF_2026-10-31 C3): is the coherence loss in proliferating and ER− tumours a swamped clock or a disordered one? Readout level **L1**

**Why:** coherence R is a signal-to-noise quantity. Lower R can mean (a) disordered gene relations, or (b) intact
relations with less clock-plane variance relative to non-clock variance in the clock genes ("swamped"). The variance
decomposition of `clocksync.decompose.components` separates them: `log_plane` (clock-plane variance) and `log_off`
(off-plane variance), from the same within-stratum-centred covariance as R.
**Data and contrasts:** exactly the registered coherence contrasts, same strata, centring and arms: SCAN-B ER− vs ER+
(H127.2) and high vs low proliferation (H127.1); METABRIC ER− vs ER+ (H166.1); GSE72094, GSE68465 (H171), GSE62254,
GSE15459 (H173) high vs low proliferation tertile. Seven contrasts.
**Statistics (low-coherence arm minus reference arm, combined over strata as the registered coherence statistic):**
Δlog_plane, Δlog_off, and the in-plane log-odds Δs = Δlog_plane − Δlog_off; ΔR reproduced.
**Tests (one-sided, the registered within-stratum label permutations of each contrast, 2,000):**
- **H176.1 (swamped):** Δs < 0 in each of the seven contrasts; Holm over seven. Supported if ≥ 5 of 7 pass.
- **H176.2 (disordered, per contrast):** "disordered" is called where Δs is not < 0 (H176.1 fails) **and** the MDE of
  Δs is ≤ |Δs| expected under swamping (the log-odds change that alone reproduces the observed ΔR, from the pseudo-clock
  calibration below); otherwise "unresolved".
**Calibration (before reading Δs):** on the 200 pseudo-clock sets of each cohort, the regression of ΔR on Δs gives the Δs
that a pure swamping would need to produce the observed ΔR; and the permutation null of Δs gives its MDE (2.49 × null
SD).
**Specificity companion:** Δs on the 200 pseudo-clock sets of each cohort (same permutations); swamping is reported as
"clock-gene specific" when the observed Δs is below the 5th pseudo percentile, otherwise as a general feature of these
tumours' expression (e.g. more variance in every gene set).
**Descriptive:** Δlog_plane and Δlog_off separately; the same decomposition without PER1 and NR1D1; a proliferation-
responsive pseudo pool (genes in the top decile of |correlation| with G2M, matched on mean and SD) as a second companion,
which is the C4 pool.
**Reading:** H176.1 supported → "clock signal swamped by non-clock variance" replaces "clock loss" wording (Q8);
disordered calls in ≥ 5 of 7 → "clock relations disordered"; mixed → stated per cohort.
**Code:** `paper/exploratory/A176_c3_plane_variance/`.

## A178 (2026-11-01, written before any leave-one-gene-out or enlarged-pseudo-pool statistic is computed; the full-gene statistics are seen; HANDOFF_2026-10-31 C7): are the breast, lung and gastric coherence losses robust to single genes and to the pseudo pool size? Readout level **L1**

**Why:** two passing cohorts sit at the pseudo 5 % line (GSE72094 5.5th, GSE62254 5.0th percentile, 200 sets), and only
PER1 and NR1D1 have been dropped so far.
**V178.1 (leave-one-gene-out):** SCAN-B H127.1 and H127.2, METABRIC H166.1 and H166.2, each registered statistic recomputed
with each of the 14 clock genes dropped (waveform subset as in the existing without-PER1/NR1D1 sensitivity, same strata and
2,000 permutations). Robust if, per contrast, all 14 statistics are negative and ≥ 12 of 14 have one-sided p < 0.05.
**V178.2 (pseudo pool size):** GSE72094 (H171.1) and GSE62254 (H173.1) against 1,000 new pseudo-clock sets (A168
construction, new seed; the 200 earlier sets not reused). The observed statistic's percentile is reported with its exact
binomial 95 % interval. Clock-specific if the percentile is < 5; "at the line" if the interval contains 5.
**Descriptive:** the gene whose removal moves each statistic most; the same 1,000-set percentile for the other five passing
contrasts (H127.1, H127.2, H166.1, H171.2, H173.3).
**Reading:** V178.1 and V178.2 pass → the losses are not carried by one gene or by the pseudo pool size. A V178.1 failure
names the carrying gene.
**Code:** `paper/exploratory/A178_c7_robustness/`.

## A179 (2026-11-01, written before any cell-cycle-conditioned coherence is computed; the unconditioned contrasts are seen; A176 is scored (swamped: clock-plane variance falls in proliferating tumours); HANDOFF_2026-10-31 C4, bulk part): is the proliferation coherence loss carried by clock-gene variation that tracks the cell cycle? Readout level **L1**

**Why:** the mammalian clock can phase-lock to the cell cycle, and PER/CRY/NR1D1 respond to it. If the clock genes'
variation across proliferating tumours follows cell-cycle composition, conditioning on it should restore coherence.
**Data and contrasts:** the five registered proliferation contrasts (H127.1, H171.1, H171.2, H173.1, H173.3), with their
strata, arms and centring.
**Conditioning:** within each arm × stratum, every clock gene is regressed on the tumours' G1/S (HALLMARK_E2F_TARGETS
without clock genes) and G2M (A97) scores, and coherence is recomputed on the residuals. Shrinkage S = ΔR_conditioned −
ΔR (> 0 when conditioning reduces the loss).
**Test H179.1 (one-sided):** S > 0, against 200 random score pairs (scores built like G1/S and G2M from expression-matched
random gene sets of the same size, A168 pool), per contrast; p = share of random S ≥ observed. Holm over five. Supported
if ≥ 3 of 5 pass **and** the median shrinkage is ≥ 50 % of |ΔR|.
**Descriptive:** fraction of each clock gene's variance explained by the two scores, high vs low arm; Δlog_plane after
conditioning.
**Readings:** H179.1 supported → the loss is carried by cell-cycle-associated clock-gene variation (compatible with
locking or with cell-cycle regulation of clock genes; bulk data cannot tell them apart; the single-cell part of C4
would). Not supported, with S's MDE below 50 % of |ΔR| → not explained by linear cell-cycle composition.
**Code:** `paper/exploratory/A179_c4_cellcycle/`.

## A180 (2026-11-01, written after the A19 and A165 estimates and SEs were scored and seen, so both tests are **secondary**; the bounds are fixed by existing conventions, not by the estimates; HANDOFF_2026-10-31 "Turning nulls into tested claims"): are the A19 age and A165 disease nulls equivalent to no effect?

**Why:** the paper must never say "unchanged" without an equivalence bound or an MDE.
**Bounds:**
- Age (A19 H19.1, γ per decade, all compartments): the A164 bound of ±0.15 on coupling, spread over the 40 years that
  separate GTEx's youngest and oldest decades: ±0.0375 per decade.
- Disease (A165): phase ±1.5 h; amplitude 1.5-fold (log z2 ratio within ±2 ln 1.5, since z2 is amplitude squared). These
  are the "best-powered rows detect about 1-1.5 h or 1.3-1.5-fold" scale of the A165 entry. TOST is applied only to rows
  whose MDE is at or below the bound (the others cannot be equivalent at this n).
**Test:** TOST at α 0.05 each side, i.e. the 90 % interval (estimate ± 1.645 SE) inside the bounds.
**Reading:** equivalence passes → "no change larger than the bound" for that row. Fails → the null stays "not detected,
MDE x". No multiplicity correction for equivalence (each row is its own claim). Rows are listed with their verdict.
**Code:** `paper/exploratory/A180_equivalence/`.

## A181 (2026-11-01, written after the archive A107b fetal gates on Cao 2020 were seen (no primary passes after Holm; p_acrophase 0.002-0.038, p_pairing 0.028-0.746) and before any downsampled-adult gate or responsive-gene readout; HANDOFF_2026-10-31 D1): is the mid-gestation failure informative? Depth control and the responsive-gene shared state

**Why:** a failed gate is uninformative unless the design could have passed. A107b used about 1,000 counts per cell
and pseudobulk depths of 1e6-1.4e8 (`legacy/.../build_cao.md`), with 8-23 fetuses.
**V181.1 (depth control, the precondition):** GTEx bulk, deaths at 08-18 h (fetal material is collected in hospital
hours), the A107b H107b.1 organs (brain cortex for cerebrum, heart left ventricle, kidney cortex, lung). Draw 200
subsamples of 23 donors with all four organs. Thin each sample's clock-gene counts binomially to the Cao median
pseudobulk depth of that organ (whole-organ depths, `build_cao.md`), scaled by the sample's clock-gene share of the
library. Run the same gate (GTEx waveform, 1,000 permutations, pass = p_pairing and p_acrophase < 0.05). Repeat with 11
donors and the lung within-organ analogue: no adult within-organ bulk exists, so H107b.3-4 get no control and stay
uninformative. **Sensitivity = pass rate.** The A107b failure counts as "no shared clock state at a depth and n where
adults show one" only if the pass rate is ≥ 0.8. Below 0.8 it stays uninformative (stop rule: no further amendment).
GTEx is used here as a benchmark for fetal material, not as a GTEx-only analysis. Only aggregate pass rates are
reported.
**H181.2 (responsive-gene shared state; one-sided):** in the A107b whole-organ and stromal designs, the same gate on the
immediate-early and glucocorticoid-responsive genes of the Cao panel. Gene set: FOS, FOSB, JUN, JUNB, EGR1, EGR3, IER2,
ATF3, NR4A1, DUSP1, FKBP5, TSC22D3, ZBTB16, DDIT4, KLF9. The pairing statistic is the cross-organ covariance of the
responsive-gene state against shuffled fetal pairing (`Dataset.gate` p_pairing on these genes with a flat waveform:
pairing only, no acrophase test, since these genes have no waveform). Pass: p_pairing < 0.05 in either design (Holm
over two). Companion: 200 expression-matched pseudo sets from the 68-gene panel's housekeeping/other genes are too few,
so the companion is the same pairing statistic on the 14 clock genes (A107b p_pairing 0.746 / 0.269, seen).
**Reading:** V181.1 ≥ 0.8 and H181.2 passes → "fetal organs share a cue-driven (responsive-gene) state but no clock
state", the mirror image of A109. V181.1 ≥ 0.8 and H181.2 fails → no shared state of either kind at this depth.
V181.1 < 0.8 → the fetal clock result stays uninformative, and H181.2 is reported descriptively.
**Code:** `paper/exploratory/A181_fetal_depth/`.

**A181 design note (2026-11-01, before any GTEx thinning or gate and before H181.2; only the A107b reproduction has run):**
- "23 donors with all four organs" cannot be drawn. GTEx 08-18 h deaths with all four organs: at most 12 donors (BA9).
- The fetal design does not have all four organs either: its 23 fetuses carry 8-11 organs each of Cerebrum, Heart,
  Kidney and Lung.
- V181.1 therefore matches the fetal **organ-presence pattern**. Each subsample draws 23 donors from the 08-18 h deaths
  (any Hardy class) that have ≥ 2 of the 4 organs. The 23 H107b.1 fetal organ sets are assigned to them at random,
  redrawing a donor that lacks an organ of its assigned set. Organs outside the set are masked. Pairs and per-organ n
  then equal the fetal design.
- The brain tissue is the one with more donors in the window (Frontal Cortex BA9).
- The 11-donor arm is dropped: it was the lung within-organ analogue, and adult bulk has no within-organ compartments.
- H181.2 uses an identity projector on the 15 responsive genes. A flat waveform leaves no clock plane, so
  `Dataset.gate` cannot run; the pairing null is computed in the script.
- Criteria and readings unchanged.

## A136 (registered 2026-11-01 from the 2026-10-27 draft `exploratory/A136_glucocorticoid/A136_draft.md`, PI approval 2026-11-01 "proceed with A136"; nothing computed before this commit): is the loss of organ phase relations in brain-dead organ donors a drug effect? Readout level **L3**

Registered exactly as the draft (exposure score, H136.1-H136.3, specificity companions, Holm over three, readings).
The draft's open choices take their defaults: the five-gene GR set (FKBP5, TSC22D3, ZBTB16, KLF15, DDIT4; no PER1, MT2A
or CDKN1C); H136.1 is the gate for H136.2; no external replication in this amendment. The H136.2 power step (null SD
and MDE from tertile-label permutations) is run and written before the tests. Controlled access: aggregates only.
**Code:** `paper/exploratory/A136_glucocorticoid/`.

## A182 (2026-11-01, written after a metadata-only scout (`exploratory/scouts_2026-11-01/C6_blind_cohorts.md`) and a header-only summary (`exploratory/C6_blind_types/output/header_summary.csv`: header lines, the ID_REF column and platform annotation; no expression value) and before any expression value of the series below is read; the TCGA values of the same types are seen; design choices approved by the PI 2026-11-01 (Q22); HANDOFF_2026-10-31 C6): is the proliferation clock-coherence loss specific to the types where it replicated (breast, lung adenocarcinoma, gastric)? Readout level **L1**

**Why:** no new type enters the claim unless scored blind. C6 scores types whose TCGA value was weak or null, and lung
squamous carcinoma, in cohorts never read.
**Data** (`$SCCC_DATA/c6_*`, GEO series matrices, A171 parser; tumour filters from `header_summary.csv`):
- Family 1 (weak/null TCGA value): **H182.1** GSE39582 colon (566 primary; 19 non-tumoral dropped; authors' RMA +
  ComBat), **H182.2** GSE17536 colon (177), **H182.3** GSE116918 prostate (248; 13 clock genes, no NR1D1; purity proxy
  10 of 11, no CD3E).
- Family 2 (squamous): **H182.4** GSE157009 and **H182.5** GSE157010 lung squamous (249, 235; two arms of one study,
  separate tests).
- Descriptive only: GSE65858 head and neck (253 primary; TCGA HNSC leans negative), two-sided, with its pseudo
  percentile. GSE53625 (oesophageal) is dropped (no symbol annotation). GSE14333 could not be downloaded.
**Definitions (= A171 / A173):** pipeline 14 clock genes (`a168_lib.GENES`, those present); linear matrices (max > 100)
→ log2; probes averaged per symbol; G2M (Hallmark without clock genes) and purity proxy (A97 11 genes) mean z over the
cohort's tumours; coherence `a127_lib.coh`.
**Statistic:** H171: purity-proxy tertile strata (≥ 30 tumours), T = mean over strata of coh(high G2M tertile) −
coh(low); proliferation permuted within stratum, 2,000 permutations, one-sided for a loss; Holm within each family.
**Seen TCGA values (H171 statistic):** COAD −0.108 (p 0.146), PRAD −0.101 (0.159), LUSC +0.212 (0.95), HNSC −0.170
(0.079) (`A171_luad_lihc/output/tcga_all.csv`, `A172_tcga_unseen/output/per_project.csv`).
**Equivalence (per cohort, PI 2026-11-01):** a cohort shows "no loss of the replicated size" if T − 1.645 × null SD >
−0.19 (the smallest replicated cohort effect, GSE62254 −0.187), where null SD is the SD of the permutation null.
**Specificity companion:** 200 pseudo-clock sets per cohort (A168 construction), same permutations; "clock-specific"
only below the 5th pseudo percentile.
**Power (after this commit, before the tests):** MDE = 2.49 × null SD per cohort; MDE > 0.30 → "not estimable".
**Predictions:** no Family 1 or Family 2 cohort loses coherence clock-specifically; each is equivalent at −0.19.
**Reading:**
- All cohorts fail and are equivalent → the loss is specific to breast, lung adenocarcinoma and gastric cancer.
- A Family 1 cohort passes clock-specifically → that type's TCGA value was a false negative; the type joins the claim.
- A Family 2 cohort passes clock-specifically → the loss is not adenocarcinoma-specific in lung; the histology wording
  is dropped (cf. A173 H173.4).
- Neither pass nor equivalence → "not resolved at this n".
**Descriptive:** each test without PER1 and NR1D1; T per stratum; GSE39582 by dataset (discovery / validation).
**Code:** `paper/exploratory/C6_blind_types/` (renamed in the code header as A182).

## A183 (2026-11-01, written after the GSE154200 prep, gate and label-permutation nulls (`scouts_2026-11-01/prep_gse154200/output/{gates.csv,power.csv}`) and before any lesional-minus-non-lesional or between-diagnosis value is computed; the archive psoriasis result is seen (A54-A59); PI approval 2026-11-01 (Q22); HANDOFF_2026-10-31 H1): does the psoriatic lesion delay replicate blind in a single-centre cohort, and is it specific to psoriasis? Readout level **L3**

**Why:** A54 replicated a lesional clock delay in psoriasis on held-out GEO series (0.93 h, p 0.004; discovery 1.50 h);
A59 found no psoriasis-sized delay in other lesional skin diseases. GSE154200 samples psoriasis, nummular eczema and
atopic eczema the same way at one centre, so diagnosis is not confounded with series.
**Data:** `$SCCC_DATA/skin_gse154200/GSE154200_countMat.txt.gz`; persons with both biopsies in one batch: psoriasis 87,
nummular 52, atopic 46. Gate passes in each diagnosis and pooled (p ≤ 0.002). 14 of 14 clock genes.
**Estimator:** `clocksync.paired` (counts, GTEx waveform, one projector, cohort = batch), offset = arg(z_lesional /
z_nonlesional) in hours, negative = lesion delayed. d* is not estimable on two compartments.
**Tests (one-sided, Holm over two):**
- **H183.1:** psoriasis lesional offset < 0 (per-person label-swap null, 10,000). MDE 0.81 h.
- **H183.2:** psoriasis offset − mean offset of the two eczemas < 0 (diagnosis labels permuted within batch, 10,000).
  MDE 0.69 h.
**Specificity companions (same permutations):** 200 expression-matched pseudo-clock sets (`pseudo_means` from the
genome-wide counts) and the responsive-gene companion (`pseudo_responsive=True`); "clock-specific" only beyond the q05 of
both.
**Descriptive:** eczema offsets with bootstrap intervals; coherence and r per diagnosis; leave-one-gene-out; without PER1
and NR1D1; offset by batch.
**Reading:** H183.1 passes clock-specifically → the delay replicates blind in a single-centre cohort. H183.2 passes → it
is a psoriasis property, not one of eczematous inflammation. A fail bounds the delay below about 0.8 h here and does not
retract A54.
**Code:** `paper/exploratory/A183_skin_gse154200/`.

## A184 (2026-11-01, written after the GSE27272 prep, gates and nulls (`scouts_2026-11-01/prep_d2_gse27272/output/{gates.csv,power.csv}`; the cord blood × mother gate is seen, p_pairing 0.022, p_acrophase 0.013) and before any pairwise coherence or offset is computed; PI approval 2026-11-01 as a **secondary** (marginal design); HANDOFF_2026-10-31 D2): does the term fetus follow its mother's internal time? Readout level **L2** (pairwise)

**Data:** GSE27272 (69 vaginal deliveries, GPL6883, 13 clock genes, no PER3); 60 deliveries with maternal and cord blood.
**Estimator:** `clocksync.paired` on (maternal blood, cord blood) per delivery, log2 intensities; coherence = length of the
mean of the per-delivery unit offset vectors (the pairwise complex coherence of the power step).
**Test H184.1 (one-sided):** coherence > its null from shuffling cord-blood states across deliveries (5,000). Null q95
0.224 (`power.csv`).
**Companions (same shuffles):** 200 expression-matched pseudo-clock sets and the responsive-gene sets
(`pseudo_responsive=True`; delivery stress and glucocorticoid response). "Clock-specific" only above the q95 of both.
**Descriptive:** mean offset (cord minus mother) with bootstrap interval; the same without PER1 and NR1D1. Placenta ×
mother failed its gate and is not read. Smoking is not estimable (arm gates fail).
**Reading:** H184.1 passes clock-specifically → the term fetus shares its mother's clock state at delivery (blood only;
shared delivery-time structure not fully excluded). A pass beaten by the responsive companion → shared delivery state,
not a shared clock. A fail is "not resolved at n 60".
**Code:** `paper/exploratory/A184_mother_fetus/`.

## A185 (2026-11-01, written after A175 and A176 were scored (A175 H175.1 D +0.158; A176 swamped 7 of 7) and before any tumour-amplitude-matched pairing is computed; second amendment on C1, allowed by the PI 2026-11-01 (Q21.2); **secondary**): is the lower tumour–normal pairing in proliferating tumours explained by less tumour clock signal? Readout level **L2**

**Why:** A176 showed that proliferating tumours have less clock-plane variance. Pairing (`q3_lib.corr_stat`) falls when the
tumour clock state is small even if the tumour still follows its normal tissue, so A175's D can be "less signal" rather
than "less following". Comparing pairing at a common tumour amplitude separates the two (the d* logic).
**Data:** the A175 set and arms (517 H138.4 pairs, cancers with ≥ 10 pairs per arm, G2M tertiles within cancer).
**Matching:** tumour amplitude a = |z_tumour| in the H138.4 state space (the states entering `corr_stat`). Within each
cancer, pairs of the low and high arms are pooled and cut into amplitude tertiles; D_match = pairs-weighted mean over
cancer × amplitude bins of pairing(low) − pairing(high), bins with ≥ 5 pairs per arm.
**Test H185.1 (one-sided, D_match > 0):** proliferation labels permuted within cancer × amplitude bin, 10,000.
**Power (before the test):** MDE = q95 + 0.84 SD of the null; MDE > 0.3 → "not estimable"; the reading uses the A175 D
0.158 as reference.
**Companions:** D_match on the 200 pseudo-clock and 200 responsive sets of A175 (same permutations).
**Descriptive:** the share of A175's D removed by matching, (D − D_match) / D; matching on normal-tissue amplitude
instead; amplitude distributions per arm (aggregates).
**Reading:** D_match > 0 (passes) → proliferating tumours follow their normal tissue less beyond having less clock signal.
D_match fails with MDE ≤ 0.158 → A175's pairing drop is explained by less clock signal (swamped), and the C1 wording is
"less clock signal, host following not shown to change". Otherwise "not resolved".
**Code:** `paper/exploratory/A185_c1_spread_matched/`.

## A186 (2026-11-01, written after the CPTAC metadata scout (`exploratory/scouts_2026-11-01/CPTAC_scout.md`: pair counts and file sizes only) and before any CPTAC expression file is read (download running); approved by the PI 2026-11-01 ("download + register anyway"); the A175 replication clause): does the proliferation pairing drop (A175 H175.1) replicate in independent CPTAC-3 patients? Readout level **L2**, **primary**

**Data:** CPTAC-3 open STAR gene counts (`$SCCC_DATA/cptac_<type>/`, manifests and sample sheets in
`scouts_2026-11-01/cptac_manifests/`): cases with a Primary Tumor and a Solid Tissue Normal file, LUAD 215, CCRCC 151,
LSCC 96, SKCM 66, HNSCC 61, PDAC 52, UCEC 45, STAD 27 (scout type mapping; cases without a clean type are dropped).
Several aliquots of one sample type in one case: unstranded counts summed. No TCGA patient is included.
**Pipeline (= H138.4 / A175):** per type, clock states of tumour and normal on the A138 projector and the per-type gate of
H138.4 (`A138_cancer_estimability`); only gated types enter. G2M = A97 score (Hallmark G2M without clock genes) of the
tumour sample, tertiles by rank within type; types with ≥ 10 pairs per arm.
**Test H186.1 (one-sided):** D = pooled pairing(low tertile) − pooled pairing(high tertile) > 0, pairs-weighted pooling
over types (`q3_lib.corr_stat`, as A175); proliferation labels permuted within type, 10,000; p = (1 + #{D_null ≥ D}) /
10,001.
**Power and stop rule (before the test):** MDE = q95 + 0.84 SD of the label-permutation null (labels only). If MDE >
0.158 (the TCGA D) the test is not run and the result is "not estimable in CPTAC".
**Specificity companions (same permutations):** 200 expression-matched pseudo-clock 14-gene sets per type and 200
responsive sets (A175 construction); "clock-specific" only above the 95th percentile of both.
**Descriptive:** D per type; the A185 tumour-amplitude-matched D_match; the reproduction of the gates per type.
**Reading:** H186.1 passes clock-specifically → the drop replicates in independent patients and C1 becomes a replicated
L2 result (proliferation only; ER cannot be tested, no breast pairs). A fail at MDE ≤ 0.158 → not replicated. Gate
failures in a type are uninformative and only remove that type.
**Code:** `paper/exploratory/A186_cptac_c1/`.

## A187 (2026-11-01, written after the C4 scout (`scouts_2026-11-01/c4_single_cell/SCOUT_C4.md`: cell counts, depth, clock-gene detection) and the CoPhaser cell-cycle refits (`C4_single_cell/pilot`, `C4_single_cell/fits`: phases and QC only), before any clock-gene value is read across cell-cycle phase; the archive C2 result on overlapping cells is known: CoPhaser G2/M − non-cycling offset −0.02 h [−0.55, 0.51], 9 studies, `figure_5/output/cancer/decision.csv`): is the clock of malignant cells locked to their cell cycle? Readout level **L3** (within-patient phase offset between cell-cycle states, `clocksync.paired`), **primary** (PI 2026-11-01: primary despite the archive C2 result on overlapping cells; bins NC / S / G2M and stop rule MDE > 1.0 h approved)

**Why:** coupled-oscillator work predicts that the mammalian clock phase-locks to the cell cycle (Bieler 2014, Feillet
2014, Droin 2019). Locking would make the proliferation coherence loss (A97, A175, A186) a mechanism. The bulk test
(A179) found no linear cell-cycle composition effect; the archive single-cell test (C2) used older CoPhaser fits, two of
them degenerate (GSE183904, SRP327758: 100 % cycling, S→G2M order failed). Both are valid in the current fits.
**Data:** CNV-malignant cells of every study whose current CoPhaser fit is valid (validation plot and S→G2M order,
`fits/output/qc.csv`; fit chosen without clock genes, which are excluded from rhythmic, context and score sets), plus Wu
and Pal breast epithelial tumour cells. Studies whose kept seed is unstable (circular correlation with another usable
seed < 0.8 over cycling cells) are excluded.
**Bins (per patient):** NC = cycling probability < 0.5; among cycling cells, each cell goes to the nearer (circularly) of
the study's S-score peak and G2M-score peak along theta: S bin, G2/M bin. Clock-gene counts are summed per patient × bin;
depth = total UMI of the bin. A patient enters a contrast with ≥ 50 cells in each of its two bins.
**Estimator:** `clocksync.paired` (counts, GTEx waveform, unit = patient, cohort = study, reference = NC), offset =
arg(z_bin / z_NC) in hours; reference gate on the NC pseudobulks per study (`gated_only=True`).
**Power and stop rule (before the test):** MDE = 2.49 × SD of the per-patient label-swap null of the pooled offset. If MDE >
1.0 h the test is not run ("not estimable"). (Locking 1:1 to a cycle of about a day separates G2/M and G1 clocks by several
hours.)
**Test H187.1 (two-sided):** pooled G2/M − NC offset ≠ 0; per-patient label-swap null, 10,000.
**Specificity (same swaps):** 200 expression-matched pseudo-clock sets (`pseudo_means` from the malignant pseudobulks) and
the responsive companion. "Clock-specific" only if |offset| exceeds the pseudo 95th percentile of |offset|.
**Artefact checks (pre-specified, descriptive):** (a) without the short-lived clock mRNAs NR1D1 and PER1 (`drop_genes`,
as in A183; chosen a priori, not from the archive per-gene table), i.e. mitotic transcript loss; (b) the same contrast in cycling T/NK and myeloid cells where
fits allow (no locking expected to show there if it is a malignant-cell property; an equal offset means artefact);
(c) library-size-matched subsampling of the G2/M bin to the NC depth.
**Secondary H187.2 (runs only if MDE ≤ 1.0 h):** |pooled offset| in high minus low within-study tertile of patient cycling
fraction > 0 (tertile labels permuted within study, 10,000).
**Descriptive:** S − NC offset; r (log variance ratio) per contrast; per-study offsets; leave-one-study-out; breast by ER.
**Reading:**
- H187.1 passes clock-specifically and survives (a)-(c): the malignant clock is phase-locked to the cell cycle. That
  makes the proliferation coherence loss a candidate mechanism, and H187.2 tests the dose.
- A fail at MDE ≤ 1.0 h bounds any G2/M–G1 clock offset below about 1 h, so phase locking does not explain the
  coherence loss. This replicates archive C2 on corrected fits.
**Code:** `paper/exploratory/C4_single_cell/` (A187).

**A187 inclusion, fixed before scoring (2026-11-01; after the CoPhaser fits and seed checks, no clock value read):** "stable"
= the kept seed's theta has circular correlation ≥ 0.8 with at least one other usable (finite, non-vetoed) seed over the
kept seed's cycling cells; a study with a single usable seed gets seeds 5-9 first (done for GSE131907, GSE176078, GSE161529,
ERP149257). Included (valid plot and stable): GSE178341, GSE132465, GSE131907, GSE176078 (Wu), GSE161529 (Pal), SRP327758,
GSE181919, SRP416180. Excluded: unstable GSE183904 (max 0.33), GSE294300 (0.77), GSE315534 (0.74), GSE166555, GSE200997,
GSE335070; no second usable seed ERP149257; invalid fits (all seeds vetoed) GSE159115, GSE212966, GSE270680
(`C4_single_cell/fits/output/{qc.csv,seed_stability.csv}`, `FITS.md`). Wu and Pal malignant cells are the authors' tumour
epithelial labels (no CNV calls).

## A187b (2026-11-01, POST HOC: written after the A187 gate p-values were seen, `A187/output/gates.csv`, and before any offset was read; requested by the PI): A187 with the gate taken over all patients with ≥ 50 NC cells

**Change (only this):** the per-study reference gate is computed on the NC pseudobulks of every patient with ≥ 50 NC cells
(not only the patients in the contrast), same `n_gate` and seed as A187. Known at registration: this gate passes in GSE131907
(p 0.025, 14 patients) and SRP327758 (p 0.035, 22 patients) and in no other study. Everything else as A187: same 8 studies,
bins, ≥ 50 cells per bin per contrast, `paired` with min_pairs = 10 (as in the A187 scoring; min_pairs = 4 reported as a
sensitivity), power step first with stop if MDE > 1.0 h, then H187.1, companions, artefact checks (a)-(c), H187.2 and the
descriptives. Status: **post hoc**; any pass is exploratory and is not counted among the pre-specified results.
**Code:** `paper/exploratory/C4_single_cell/A187/` (A187b outputs with suffix `_b`).

## A188 (2026-11-01, written before any C5 contrast, power value or per-arm gate is computed; the A138 Q1 single-cell set and H138.1 are seen: malignant host coupling 0.076 [0.014, 0.139], MDE about 0.2; the bulk proliferation pairing drop is seen (A175, A186); requested by the PI 2026-11-01 with the rule "test if powered, else descriptive only"; HANDOFF_2026-10-31 C5): do malignant cells follow their own tumour stroma less in high-proliferation (and ER−) tumours? Readout level **L2** (anchored coupling at a common host spread, d*), **secondary**

**Why:** A175/A186 show that proliferating tumours follow the patient's normal tissue less (in CPTAC not clock-specific).
C5 asks the within-tumour version: whether malignant cells couple less to the stromal host of their own tumour.
**Data:** the A138 Q1 single-cell set: cohorts and patients of H138.1 (`cancer_sc`, copy-number-called malignant cells,
else tumour epithelium). Host = the cohort's learned A17 stromal anchors, as in H138.1; ≥ 2 anchors are needed for d*.
**Groups:** proliferation = the A97 G2M score (Hallmark G2M without clock genes) of the patient's malignant pseudobulk,
tertiles by rank within cohort; high vs low (middle dropped). ER (H188.2) = breast cohorts with ER labels (Wu, Pal; labels as
in A170 and `scouts_2026-10-31/sc_breast_er`), ER− vs ER+.
**Estimator:** per cohort `Dataset.contrast(groups, levels=(low, high), targets=malignant, companion=<myeloid, else the mean
of the non-anchor non-malignant compartments>, require_host_gate_each_arm=True, mde=True, anchor_readout="loo")`. A cohort
enters with ≥ 5 patients per arm and both arm gates passing. Pooled d* = inverse-variance (se_d_star) mean over the
entering cohorts; the pooled null is the same weighting over the per-cohort permutation draws (labels permuted within
cohort).
**Power and branch (binding, decided now):** pooled MDE on d* = 2.49 × pooled se of d*, with pooled se = (Σ 1/se_d_star²)^−½.
The per-cohort `power.mde` values are reported alongside.
- **MDE ≤ 0.3: test.** H188.1 (one-sided): pooled d*(malignant, high − low) < 0, p from the pooled permutation null
  (≥ 1,000 per cohort). Specificity: d*(malignant) − d*(companion) on the same permutations; "malignant-specific" only if
  that difference is also < 0 at p < 0.05.
- **MDE > 0.3: descriptive only.** The same quantities are computed and reported with intervals. There is no p-value
  verdict, and the result is "not estimable at this n".
- H188.2 (ER− vs ER+) follows the same branch rule on its own MDE.
**Descriptive (either branch):** per-cohort d*, d, c per arm, z2 and psi ratios, offsets; gates per arm; leave-one-cohort-out.
**Reading:** a powered pass means proliferating malignant cells follow their own stroma less, the within-tumour
counterpart of A175. A powered fail bounds the drop below the MDE. Descriptive-only output is a feasibility statement,
not "no change".
**Code:** `paper/exploratory/A188_c5_within_tumour/`.

## A189 (2026-11-02, written after the D1b feasibility calibrations on adult material (`exploratory/D1b_fetal_feasibility/FEASIBILITY.md`, commits fef03e4, 3226f03) and before any fetal coherence value, kidney within-organ gate or genome-wide fetal pseudobulk is computed; seen: the archive A107b gates on the same fetuses (whole organs, stromal cells per organ, within lung p_pairing 0.028 / p_acrophase 0.003, within heart 0.093 / 0.023), the A181 depth control, and fetal clock-gene counts next to adult expectations (`markers_expected.csv`)): which mid-gestation fetal organs and cell types carry clock-shaped gene organisation? Readout levels **L1** (coherence within one compartment) and **L2** (gate across cell classes within one organ); **primary** for the blind tests, **secondary** for the A107b re-reads

**Why:** human tissue clocks cannot be timed in utero, and imaging is not possible in humans. Snapshot atlases are the
only human route to where and when a peripheral clock is present. A181 showed the cross-organ gate cannot pass at the Cao
design (few fetuses share organs; adult pass rate 0.20). Readouts within one compartment or one organ avoid that limit.
Inclusion was fixed on adult calibrations at the fetal design (sensitivity = pass rate over 200 adult subsamples cut to
the fetal n and depth), threshold ≥ 0.5 (user, 2026-11-02: keep underpowered compartments; a fail is uninformative anyway).
**Data:** Cao 2020 (GSE156793), 28 fetuses, 72-129 days. Organ pseudobulks with ≥ 1,000 cells, class pseudobulks with ≥ 200
cells (`fetal_cao`, classes as `build.CLASSMAP`). Genome-wide pseudobulks of the same rows are rebuilt from the GEO loom
(`fetal_cao` builder with all genes; the 21 GB decompressed loom is temporary and deleted after the build).
**H189.1 (L1, blind, primary; one-sided, Holm over 9):** coherence R (`decompose.components`, GTEx waveform) of one
compartment across fetuses exceeds the per-gene donor permutation null (10,000 permutations; breaks gene covariance,
keeps each gene's distribution). Compartments (calibration pass rate, thinned, `passrates.csv`): Liver 0.805, Lung 0.715,
Adrenal 0.705, Heart 0.51 (whole organs, GTEx 08-18 h calibrators); Heart:stromal 0.935, Lung:stromal 0.815,
Heart:endothelial 0.73, Lung:epithelial 0.62, Kidney:parenchyma 0.55 (classes, adult nucleus-atlas calibrators).
**H189.2 (L2, blind, primary; one test):** within the fetal kidney (parenchyma, stromal; 9 fetuses; calibration 0.76),
`Dataset.gate` (10,000 permutations) passes: p_pairing < 0.05 and p_acrophase < 0.05.
**Specificity (both):** 200 expression-matched pseudo-clock sets (`paired.pseudo_waveforms`, kind "matched", from the
compartment's genome-wide pseudobulk; for H189.2 the mean log expression over both classes), each carrying the GTEx
amplitudes and phases. H189.1: "clock-specific" only if R exceeds the pseudo 95th percentile of R in that compartment.
H189.2: "clock-specific" only if the clock p_pairing is below the pseudo 5th percentile of p_pairing (gate rerun with each
pseudo waveform, 1,000 permutations). The responsive companion (kind "responsive") is reported for both, descriptive.
**H189.3 (secondary, seen, not counted):** the A107b within-heart (cardiomyocyte, stromal, endothelial; 9 fetuses;
calibration 0.94) and within-lung (epithelial, stromal, endothelial; 11 fetuses; no adult calibration) gates rescored on
current code with the same specificity companion.
**Descriptive:** leave-one-gene-out R; R per compartment with fetal ages listed (no age test: 5-11 fetuses); marker table
(seen) with adult expectations.
**Reading:**
- A clock-specific pass places clock-shaped gene organisation in that fetal compartment at 13-18 post-conception weeks
  (L1), or clock-shaped covariance shared across cell types of one fetal organ (L2: a coupled organ clock).
- A fail is uninformative: a weaker or less coherent clock, or too few fetuses (sensitivity 0.51-0.94 even for adult clocks).
- Compartments are not compared with each other and organs are not ranked by onset: their sensitivities differ.
**Stop rule:** no further amendment of D1/A189; whatever the outcome, it is written up.
**Code:** `paper/exploratory/A189_fetal_within/`.

## A190 (2026-11-02, written after A189 was scored (commit 99529d5; seen: fetal adrenal whole-organ coherence passes clock-specifically, R 0.734, Holm 0.0495) and after a proxy power calibration on adult material (`exploratory/D1b_fetal_feasibility`, `--adrenal`); before any fetal adrenal class-level coherence or within-adrenal gate is computed; not scored in the archive (A107b gated heart and lung only)): is the fetal adrenal clock shared across its cell types? **L2 gate within one organ**, **primary** (one test), plus a descriptive L1 per class

**Why:** user request (2026-11-02) after the A189 adrenal positive. A189 shows clock-shaped organisation in the whole organ,
which is dominated by the cortex. A gate across cell classes asks whether the cortex, the stroma and the vessels of the
same fetuses share one clock state (a coupled organ clock), not only whether the cortex carries one.
**Selection:** the organ was chosen after its A189 L1 pass on the same 7 fetuses. The two readouts are not independent
(the parenchyma enters both), so a pass is a follow-up of A189, not a replication of it.
**Data:** Cao 2020 class pseudobulks (≥ 200 cells): Adrenal:parenchyma (adrenocortical cells, 94%, and chromaffin cells,
`build.CLASSMAP`), Adrenal:stromal, Adrenal:endothelial, all present in the same 7 fetuses (89-122 days). Neural
(sympathoblasts, Schwann; 6 fetuses) and myeloid (5) are left out.
**Power (proxy, no adult adrenal nucleus atlas):** adult heart (cardiomyocyte, fibroblast, endothelial) and kidney
(proximal tubule, interstitial, endothelial) atlases cut to 7 donors and thinned to the adrenal class depths: pass rate
0.90 [0.85, 0.94] and 0.875 [0.82, 0.92] over 200 subsamples (`output/passrates_adrenal.csv`).
**H190.1 (L2, blind, primary):** `Dataset.gate` (no anchors, GTEx waveform, 10,000 permutations) on the three classes
passes: p_pairing < 0.05 and p_acrophase < 0.05. Specificity as H189.2: "clock-specific" only if the clock p_pairing is
below the 5th percentile of p_pairing over 200 expression-matched pseudo-clock sets (gate rerun with each, 1,000
permutations; mean log expression over the three classes); responsive sets reported, descriptive.
**Descriptive (not tested):** L1 coherence R per adrenal class against the per-gene donor permutation null (10,000) with
the matched pseudo 95th percentile; leave-one-gene-out range of the gate's p_pairing is not computed.
**Reading:**
- A clock-specific pass: the fetal adrenal cortex, stroma and vessels share one clock state at 13-17 weeks (a coupled
  organ clock), on the same fetuses as A189.
- A fail is uninformative (proxy sensitivity 0.875-0.90, adult sampling times unknown); it does not show the classes
  run separate clocks.
- The class-level L1 values are descriptive: they do not localise the A189 clock to one class (no multiplicity control,
  stromal and endothelial pseudobulks are 30-70 times shallower than the parenchyma).
**Stop rule:** one amendment; no further adrenal sub-question on these fetuses without a second fetal dataset.
**Code:** `paper/exploratory/A190_fetal_adrenal/`.

## A191 C8 (2026-11-02, started on the user's instruction; draft `exploratory/C8_internal_time_outcomes/DRAFT_C8.md`, family chosen by the user's delegation; PI review pending): do a tumour's clock signal and its phase offset from its own normal tissue relate to recurrence? **L1 and L3, primary (Holm over two), exploratory in spirit**

**Seen before writing:** outcome event counts and covariate distributions (metadata only); A166 group-level METABRIC
coherence (ER, proliferation) and A176 group-level in-plane log-odds (swamped clock); A127 SCAN-B off-plane score vs
survival (HR 0.94 per SD, p 0.35); A177 GSE14520 tumour-normal pair gate (passes: p_pairing 0.0005, p_acrophase 0.010)
and its G2M tertiles; A171 GSE14520 tumour vs normal coherence. No per-tumour clock readout has been related to any
outcome in either cohort.
**H191.1 (L1, METABRIC, ER+ by IHC with hormone therapy, relapse-free survival):** per-tumour in-plane log-odds
s_i = log(y_i' P y_i) - log(y_i' (I - P) y_i), y_i = the tumour's 14 clock genes (log2 x ln 2, gene-centred within
METABRIC cohort), P = projector on the GTEx waveform plane (`clocksync.decompose.sample_log_odds`); the per-sample
counterpart of the A176 readout. Cox model (statsmodels PHReg, strata = METABRIC cohort), s standardised, adjusted for age,
tumour size, grade (categorical), log(1 + positive nodes), PAM50 (claudin-low kept as its own level) and the G2M score
(A166 `prolif`); complete cases. Statistic: Wald z of s; two-sided.
**H191.2 (L3, GSE14520 HCC, recurrence):** per-patient folded phase distance of the tumour from its own non-tumour
liver, relative to the cohort offset: states from `paired.cohort_states` (arrays, CIART absent so 13 genes), centred
across patients, d_i = |arg(z_T,i conj(z_N,i) conj(c))| in hours (0-12), c = sum z_T conj(z_N). Cox model, d standardised,
adjusted for age, sex, BCLC (0/A vs B/C), cirrhosis, AFP (high), size (large), multinodular; complete cases. Wald z, two-sided.
**Power step (before any outcome is read, written to `output/power.csv`):** MDE of the HR per SD at 80 % power, two-sided
alpha 0.05, from the number of events and the variance of the readout not explained by the covariates (Schoenfeld with
1 - R^2 of the readout on the covariates). MDE above HR 1.5 per SD: not estimable, not run.
**Specificity (both):** the same Cox z with the readout computed from 200 expression-matched pseudo-clock sets (A166/A127
construction: each clock gene replaced by one of its 20 nearest non-clock, non-G2M, non-purity genes in mean and SD,
waveform amplitudes and phases carried over). "Clock-specific" only if |z_clock| exceeds the 95th percentile of |z_pseudo|.
**Pass:** Holm-adjusted p < 0.05 over H191.1 and H191.2 and clock-specific.
**Secondary (not counted):** H191.1 for overall survival and in all METABRIC tumours (ER as an added covariate); H191.2 for
overall survival; the unadjusted HRs; GSE14520 internal time at surgery (non-tumour liver phase relative to the cohort,
cos and sin terms, likelihood-ratio test, 2 df), read as surgery hour plus individual phase, no causal reading.
**Reading:**
- A clock-specific association independent of proliferation and stage would make the tumour clock a prognostic variable
  (L1) or tumour-host desynchrony a prognostic variable (L3).
- The prediction for H191.1 is no independent association (A176: the loss is swamping by proliferation). A fail is
  reported with its MDE, never as "no effect" without it.
- A pass that is not clock-specific means a matched gene set carries the same prognostic information.
**Stop rule:** one amendment. NHS, SCAN-B, POETIC and the A146 replication are not part of A191.
**Code:** `paper/exploratory/C8_internal_time_outcomes/`.

## A192 (2026-11-02, on the user's instruction "replication/explanation" of published clock-survival claims; literature in `literature_positioning.md` section 10): is the published prognostic value of clock genes in breast cancer (Cadenas et al., Cell Cycle 2014) carried by proliferation? **Fresh cohorts, primary (Holm over two)**

**Seen before writing:** event counts only (Mainz GSE11121 200 tumours, 46 distant metastases; Rotterdam GSE2034 286, 107;
Transbig GSE7390 198, 62; all node-negative, no systemic therapy, U133A); Cadenas' published per-gene HRs and their
correlations with a proliferation metagene (CRY2 -0.612, PER2 -0.440, PER3 -0.409, TIMELESS +0.674). No expression value of
these cohorts has been looked at. Data: `paper/common/datasets/breast_prognosis_geo/build.py` (log2, one probe per gene;
CIART not on U133A, so 13 of the 14 clock genes).
**Definitions.** Endpoint: distant-metastasis-free survival (months). Within-cohort z-scores of log2 expression. G2M:
mean within-cohort z of the hallmark G2M genes present (A127 definition). ESR1 and ERBB2 expression stand in for ER and HER2
(recorded ER is missing in Mainz). Cox (statsmodels PHReg, Efron ties), strata = cohort, readouts standardised.
**H192.1 (reproduction, a pipeline check, not counted):** univariate HR of each gene Cadenas reported (protective: CLOCK,
PER1, PER2, PER3, CRY2, NPAS2, RORC; adverse: TIMELESS). Reproduced if same direction with BH q < 0.05 over the 8.
**H192.2 (primary, mechanism):** clock score S = mean z of the 7 protective genes.
(a) Specificity against proliferation: the unadjusted Cox z of S against 1,000 proliferation-matched pseudo-scores, each
gene replaced by one of its 20 nearest pool genes (non-clock, non-G2M, non-purity) in standardised (mean, SD, correlation
with G2M). Clock-specific if z_S lies outside the 2.5-97.5 percentiles of the pseudo z.
(b) Adjusted: Wald test of S adjusted for G2M, ESR1 and ERBB2.
Reading: "carried by proliferation" (prediction) if (a) is not clock-specific and (b) p >= 0.05; "clock-specific
prognosis" if (a) is clock-specific and (b) p < 0.05 in the published direction; otherwise "mixed", reported as such.
**H192.3 (primary, fresh test of A191):** per-tumour in-plane log-odds (`decompose.sample_log_odds`, 13 genes, GTEx waveform
restricted to them, gene-centred within cohort), adjusted for G2M, ESR1, ERBB2; specificity against 200 expression-matched
pseudo-clock sets as in A191.
**Power step first** (Schoenfeld with 1 - R^2 on the covariates, written before outcomes are read); MDE above HR 1.5 per SD:
not estimable, not run. **Holm** over the p-values of H192.2(b) and H192.3.
**Secondary (not counted):** Cadenas' multivariate model (Mainz and Transbig: grade, size, ESR1, ERBB2) for PER3, RORC and
TIMELESS with and without G2M; MKI67 alone and AURKA alone as proliferation terms; each gene's correlation with G2M against
Cadenas' values; the share of log HR of S absorbed by G2M (1,000 bootstrap resamples); per-cohort HRs; ER+ only;
random-effects pooling of H192.3 with A191 H191.1 and an equivalence test (TOST, margin HR 1.15 per SD).
**Stop rule:** one amendment; no other proliferation definition is tried for the primary.
**Code:** `paper/exploratory/C8b_clock_prognosis_replication/`.

## A193 (2026-11-02, same instruction): does high CYCLOPS magnitude predict worse survival in luminal A breast cancer (Li, Hammarlund, Clarke et al., PNAS 2024)? **Reproduction in TCGA, replication in METABRIC and SCAN-B, primary (Holm over two)**

**Seen before writing:** event counts (TCGA PanCanAtlas luminal A 499 patients, 57 deaths, 31 within 5 years; METABRIC PAM50
LumA 700, 378 deaths, 231 relapses; SCAN-B LumA 1,642, 115 deaths); the authors' public code (GitHub ranafi/tumor-circadian-
clock-strength-..., mortality script and training script, read only). METABRIC and SCAN-B outcomes have been related to other
clock readouts (A127 SCAN-B off-plane, A191 METABRIC ER+ in-plane, both null); CYCLOPS magnitude has never been computed in
any of these cohorts.
**Readout.** CYCLOPS 2.0 (`common.benchmark.cyclops`, commit 56dbacb, template hyperparameters, 80 models), seeds = the
CYCLOPS human seed list (`cyclops_seed_genes.txt`; the authors' Zhang/Meng list is not public), magnitude = `Magnitude` of
`fit_output.csv`; three fits per cohort (sample-order seeds 0-2), readout = mean of the within-fit z-scored magnitudes.
TCGA: luminal A tumours fitted together with TCGA adjacent normal breast, tumour status as the discrete covariate (the
paper's design; the authors' normal-only eigenspace and GTEx/Manchester training samples are not reproduced). METABRIC and
SCAN-B: luminal A tumours only (no normal tissue), batch covariate = METABRIC cohort or SCAN-B platform.
**Reproduction (TCGA, not counted):** the authors' analysis: "died within 5 years" (death at <= 5 y = 1, every other patient
0, as in their script), top tertile vs the rest, logistic regression with age and metastasis at diagnosis, and with MKI67
added; then a Cox model for OS (top tertile, and per SD), adjusted for age, stage and G2M.
**H193.1 (primary):** METABRIC + SCAN-B luminal A, overall survival, top magnitude tertile (within cohort) vs the rest, Cox
with strata = cohort, adjusted for age and G2M (their age + MKI67 model; metastasis at diagnosis is not recorded).
**H193.2 (primary):** the same, magnitude per SD, adjusted for age, tumour size, grade, log(1 + positive nodes) (METABRIC) or
node status (SCAN-B), and G2M.
Two-sided; the published direction is HR > 1. **Replicated** if Holm-adjusted p < 0.05 with HR > 1. Power step first (MDE
of the tertile HR and of the HR per SD); MDE above 1.5 (per SD) or 2.0 (tertile): not estimable.
**Secondary (not counted):** METABRIC relapse-free survival; disease-specific survival where recorded; magnitude vs G2M,
the purity score and the A191 in-plane log-odds (Spearman); the same with MKI67 instead of G2M; seed concordance of the three
fits; their 5-year logistic endpoint in METABRIC and SCAN-B.
**Stop rule:** one amendment. If CYCLOPS fails to fit a cohort (fewer than 2 eigengenes, or no fit), that cohort is dropped
and reported.
**Code:** `paper/exploratory/C8b_clock_prognosis_replication/`.

## A194 (2026-11-02, same instruction): is TimeTeller's clock-dysfunction score Theta an independent prognostic factor in breast cancer (Vlachou et al., bioRxiv 2019), and does it survive a proliferation term? **Primary (Holm over two)**

**Seen before writing:** event counts (SCAN-B all tumours 3,073, 289 deaths; TCGA BRCA PanCanAtlas, OS events to be counted
in the power step); the TimeTeller preprint (REMAGUS: favourable good clock function in grade 1-2, HR 0.24, no proliferation
term; REMAGUS survival data are not public, GSE26639 has none). Theta has never been computed in these cohorts.
**Readout.** TimeTeller (`common.benchmark.timeteller`, the method-comparison install), genes = the published human panel
present (ARNTL, NPAS2, PER1, PER2, PER3, NR1D1, NR1D2, CIART, TEF, DBP), trained on GTEx "Breast - Mammary Tissue" fast
deaths (Hardy 1-2) with recorded time of death (3 h bins, as in the benchmark), log2-normalised counts; applied to tumour
log2 expression (SCAN-B log2 FPKM, TCGA recount3 log2 CPM). Theta per tumour (continuous). Deviation: the authors trained on
oral mucosa arrays; their GCF/PCF thresholds are training-specific, so Theta is used as a continuous score and in tertiles.
GTEx: only aggregates are released.
**H194.1 (primary):** SCAN-B + TCGA BRCA (strata = cohort), overall survival, Theta per SD, adjusted for the clinical model
of the preprint (age, size, nodes, grade where recorded, ER, HER2). Published direction: higher Theta (worse clock), worse
outcome.
**H194.2 (primary):** the same plus G2M.
Two-sided, Holm over two; power step first (MDE above 1.5 per SD: not estimable).
**Secondary (not counted):** grade 1-2 subgroup (the preprint's claim); Theta tertiles; METABRIC (arrays; cross-platform from
an RNA-seq training set) with OS and RFS; Theta vs G2M and vs the A191 in-plane log-odds (Spearman); MKI67 instead of G2M.
**Stop rule:** one amendment; no retraining on another reference.
**Code:** `paper/exploratory/C8b_clock_prognosis_replication/`.

**A194 addendum (2026-11-02, before any outcome is read):** the registered readout (tumour log2 values as they are)
saturates: Theta = 1 for 97 % (SCAN-B), 94 % (TCGA) and 100 % (METABRIC) of tumours, because constant per-gene offsets
between the GTEx training scale and each tumour cohort (platform, depth, gene length) put tumours outside the training
distribution. Primary readout is therefore the "aligned" variant: each panel gene shifted so that the cohort mean equals
the GTEx training mean (shift only; no outcome involved). The raw variant is reported as uninformative. Training check:
128 GTEx donors, 8 bins, 5-fold out-of-fold median absolute time error 4.38 h (weaker than muscle or thyroid in the
benchmark). Covariates per cohort enter as cohort-specific columns of the stratified model (TCGA: age, T stage, N+,
recorded ER, ERBB2 expression; SCAN-B: age, log size, N+, grade, ER, HER2; METABRIC: age, size, grade, nodes, ER, ERBB2
expression). SCAN-B is the A127 analysis set (PAM50 and ER recorded).

## A195 (2026-11-02, on the user's question "is clock desynchrony a byproduct of division, or is division eased by the absence of a coherent clock?"; written before any clock-plane variance is computed on single-cell cell-cycle bins): is the proliferation clock-signal loss confined to cycling cells? **L1 (clock-plane variance), primary (Holm over two)**

**Seen before writing:** bulk proliferation contrasts lower the clock-plane variance (A176: Δlog_plane SCAN-B −0.228,
GSE72094 −0.545, GSE62254 −0.687, GSE68465 −0.232); linear cell-cycle composition does not explain it in bulk (A179); in
malignant single cells the G2/M − NC phase offset is null in the one gated study (A187b, −0.07 h, MDE 0.48 h); the A187
per-study NC reference gates (pass only GSE131907 and SRP327758 over all NC patients); within-tumour L2 coupling by
proliferation is not estimable (A188). Never computed: log_plane of NC or cycling pseudobulks, or of any arm defined by
cycling fraction.
**Why:** a cross-sectional association cannot order clock loss and proliferation. One mechanism can be separated: if the
loss is a byproduct of division in the dividing cells (mitotic transcriptional shutdown, cycle-driven clock-gene
transcription), it sits in cycling cells, and the non-cycling cells of highly proliferative tumours keep a full clock
signal. If it persists in their non-cycling cells, it is a property of the tumour (cell-intrinsic clock loss that may permit
division, or a shared cause), not of the act of dividing. This does not order the two causally; it rules one direction's
simplest mechanism in or out.
**Data:** the A187 pseudobulks (`C4_single_cell/A187/output/pseudobulk_clock_*.csv.gz`, `genome_*.csv.gz`; 8 studies,
malignant cells, CoPhaser bins NC / S / G2M); cycling = S + G2M summed; ALL = NC + S + G2M summed. Cycling fraction per
patient = share of the patient's malignant cells with cycling probability ≥ 0.5 (A187 `cycling_fraction`).
**Statistic:** `clocksync.decompose.components` on `logfrac` (counts, Poisson term removed), GTEx waveform, log fractions
centred within study (and within arm), pooled over studies (k = number of studies); log_plane difference between arms.
**H195.1 (cell state, within patient):** Δ1 = log_plane(cycling) − log_plane(NC) over patients with ≥ 50 cells in both;
null: NC/cycling labels swapped within patient (10,000). Byproduct predicts Δ1 < 0.
**H195.2 (tumour state, decisive):** Δ2 = log_plane(NC; high cycling-fraction tertile) − log_plane(NC; low tertile),
tertiles by rank within study (middle dropped), patients with ≥ 50 NC cells, studies with ≥ 3 patients per arm; null:
tertile labels permuted within study (10,000). The same contrast on ALL pseudobulks (ΔALL) is the reference: it should
reproduce the bulk direction.
Two-sided p, Holm over the two. **Power and stop rule (before the tests):** MDE = 2.49 × null SD; a test with MDE > 0.40
(about the median bulk |Δlog_plane|) is not estimable and is not read.
**Signal check:** NC pseudobulks must carry a clock-shaped signal: pooled gene coherence R (`components`) above the 95th
percentile of 200 expression-matched pseudo-clock sets (`paired.pseudo_waveforms`, means from the NC pseudobulks). If not,
both tests are uninformative.
**Specificity:** Δ1, Δ2, ΔALL of the 200 matched pseudo sets on the observed labels; percentile reported.
**Artefact checks (descriptive):** (c1) NC binomially thinned to the cycling depth per patient (H195.1); (c2) every NC
pseudobulk thinned to the 10th-percentile NC depth of its study (H195.2); (d) per-study Δ; S and G2M separately for Δ1.
**Reading:**
- ΔALL < 0 and Δ2 bounded (fails at MDE ≤ 0.40, |Δ2| < |ΔALL| / 2) with Δ1 < 0: the loss sits in dividing cells, a
  byproduct of division (at this resolution).
- Δ2 < 0 (pass, same sign as ΔALL): non-cycling cells of proliferative tumours also lose the signal; division in the cell
  does not explain it. Direction stays open (clock loss permissive, or a shared cause).
- ΔALL not negative: the single-cell data do not reproduce the bulk loss; uninformative.
**Limitation (registered):** NC = no S/G2M signature, i.e. G0 and G1 cells; G1 cells of cycling lineages are inside NC.
**Stop rule:** one amendment.
**Code:** `paper/exploratory/C4_single_cell/A195/`.

## A196 (2026-11-02, same question as A195; A195 not estimable; A168 and A169 scored and seen (2-week aromatase inhibition: proliferation −0.9 SD, clock coherence unchanged, POETIC T −0.010, GSE20181 +0.191, GSE87411 0.529 → 0.598); written before any clock value of GSE93204 is read): does arresting proliferation restore tumour clock coherence? **L1 within-patient, primary H196.1 (one test); H196.2 post hoc**

**Why:** if the proliferation coherence loss is a byproduct of ongoing division, stopping division should restore
coherence within days; a stable tumour property (or clock loss that permits division) predicts no change. NeoPalAna
(GSE93204) adds palbociclib to anastrozole between C1D1 and C1D15, so C1D1 → C1D15 is a CDK4/6 arrest on a constant
endocrine background, the cleanest cell-cycle intervention in open data (scout: `scratchpad`, verified GEO counts).
**Data:** GSE93204 (Agilent GPL6480 two-colour, log2 ratio to a universal reference, as GSE87411), parsed with the A168
prep functions (probes averaged per symbol, genes missing in > 20 % dropped, A97 scores over the series).
**H196.1 (primary):** T = coh(C1D15) − coh(C1D1) over patients with both samples, registered coherence (A166 `coh`, clock
genes present); null: time labels swapped within patient (2,000); two-sided. Byproduct predicts T > 0.
**Treatment check (first):** paired proliferation change C1D1 → C1D15 negative (Wilcoxon p < 0.05); otherwise
uninformative. **Power:** MDE = 2.49 × null SD; MDE > 0.40 = not estimable, T is then descriptive only.
**Specificity:** 200 matched pseudo-clock sets (A168 `pseudo_sets`), percentile of T.
**Descriptive:** BL → C1D1 (anastrozole), BL → C1D15, C1D15 → surgery; residualised on proliferation and purity.
**H196.2 (POST HOC: A168 estimates seen): reversal calibration.** For POETIC AI B → S (157), GSE20181 A → B, GSE87411
base → w2 and GSE93204 C1D1 → C1D15: predicted rise P = b × (−Δprolif), b = [coh(low) − coh(high proliferation tertile)]
/ [mean prolif(high) − mean prolif(low)] in the series' pre-treatment samples, Δprolif = mean paired change. Observed T
(A168 code) minus P, z with var(T) from the swap null and var(P) from the tertile-label null of the numerator (×
Δprolif / denominator); pooled by inverse variance. Uninformative if the pooled pre-treatment numerator is not above its
null (one-sided p ≥ 0.05). Reading: pooled T − P < 0 (p < 0.05) = proliferation arrest does not restore coherence at the
size the cross-sectional loss implies (not an acute byproduct of division); T ≈ P = compatible with the byproduct.
**Caveats (registered):** aromatase inhibition also blocks ER signalling; POETIC on-treatment samples are surgical
(the untreated arm falls, A168); two-colour ratios; no biopsy clock time.
**Stop rule:** one amendment.
**Code:** `paper/exploratory/A196_proliferation_arrest/`.

## A197 (2026-11-03, on the user's agreement to test the "clock is a property of differentiated tissue" message; written before any ES-like score is computed in any cohort): does tumour clock coherence fall with dedifferentiation beyond proliferation and subtype? **L1, primary H197.1-2 (Holm over two); H197.3-4 secondary**

**Seen before writing:** coherence falls with proliferation (A127 H127.1, A171, A173) and with ER− (A127 H127.2, A166, also
within PAM50 × proliferation, H166.2); not cell-cycle composition (A179); no recovery on arrest (A168, A196). No
dedifferentiation index has been computed in these cohorts.
**Why:** if the clock is a property of differentiated tissue, coherence should fall with dedifferentiation at equal
proliferation and subtype. ER− and histology results are compatible but confounded with proliferation and subtype.
**Index (fixed now):** ES-like score (Ben-Porath et al., Nat Genet 2008) = mean within-cohort gene z of BENPORATH_ES_1
(MSigDB C2 v2024.1, `$SCCC_DATA/reference/c2.cgp.v2024.1.Hs.symbols.gmt`) after removing Hallmark G2M,
BENPORATH_PROLIFERATION, BENPORATH_CYCLING_GENES, the purity genes, clock-like and clock genes (gene count reported).
**Statistic:** within strata, coherence (A127 `coh`) of the high minus the low ES-like tertile (tertiles by rank within
stratum); T = mean over strata; null: ES-like tertile labels permuted within stratum (2,000). One-sided, T < 0.
Breast strata = PAM50 × proliferation tertile (tertile within PAM50), ≥ 30 tumours. LUAD strata = proliferation tertile
× purity tertile (within cohort), ≥ 30 tumours.
**H197.1:** SCAN-B (A127 analysis set). **H197.2:** METABRIC (A166 set). Holm over the two.
**H197.3-4 (secondary, Holm over two):** LUAD GSE72094, GSE68465 (A171 tables).
**Power and stop rule:** MDE = 2.49 × null SD; a test with MDE > 0.15 is not estimable and not read.
**Specificity:** 200 expression-matched pseudo-clock sets (A127 construction; pool also without the ES genes), observed
labels; clock-specific if T < pseudo q05.
**Descriptive:** Spearman of the ES-like score with proliferation and ER; T without proliferation strata; PRC2 targets
score (BENPORATH_PRC2_TARGETS, expected opposite sign) as an alternative index; leave-one-gene-out.
**Reading:**
- H197.1-2 pass, clock-specifically: coherence follows dedifferentiation beyond proliferation and subtype; supports the
  "clock is a property of differentiated tissue" message.
- Fail at MDE ≤ 0.15: the ES-like index adds nothing beyond proliferation and subtype; the message rests on ER− and
  histology only and is stated as a hypothesis.
**Caveat (registered):** the ES-like signature is a transcriptional proxy for dedifferentiation, not lineage
differentiation itself.
**Stop rule:** one amendment.
**Code:** `paper/exploratory/A197_dedifferentiation/`.

## A198 (2026-11-03, on the user's instruction): does the tumour clock signal modify the benefit of chemotherapy or endocrine therapy in breast cancer, beyond proliferation? **Clock x treatment interaction, primary (Holm over two), observational**

**Seen before writing:** treatment counts (METABRIC: chemotherapy 412 yes / 1,568 no, hormone therapy 1,216 / 764, 529
missing; SCAN-B: chemotherapy 1,153 / 1,903, endocrine 2,523 / 532); main-effect outcome analyses of clock readouts in
both cohorts (A127 SCAN-B off-plane HR 0.94; A191 METABRIC ER+ hormone-treated in-plane null; A193 CYCLOPS magnitude and
in-plane log-odds correlations; A192, A194). No clock x treatment term has been fitted in any cohort.
**Readout:** per-tumour in-plane log-odds s (`clocksync.decompose.sample_log_odds`, GTEx waveform, 14 clock genes, log2 x
ln 2), genes centred within METABRIC cohort or SCAN-B platform; s standardised within cohort (A191 readout).
**Model:** Cox, overall survival, METABRIC + SCAN-B pooled, strata = cohort. Terms: s, treatment T, s x T, G2M, G2M x T,
age, size, grade (1 / 3 vs 2), nodes (METABRIC log(1 + positive), SCAN-B status), PAM50 (claudin-low its own level), ER and
HER2 where not fixed by the population, all cohort-specific as in A193 (`outcomes.block_design`); complete cases.
Statistic: Wald z of s x T, two-sided. The G2M x T term is the proliferation companion: proliferation is known to predict
chemotherapy benefit, so the clock term is read only beyond it.
**H198.1:** T = chemotherapy, all tumours.
**H198.2:** T = endocrine therapy, ER+ tumours.
**Power step first (`output/a198_power.csv`):** MDE of the interaction (ratio of HRs per SD of s) at 80 % power, two-sided
alpha 0.05: Schoenfeld with the variance of s x T not explained by all other terms. MDE above 1.5: not estimable, not run.
**Specificity:** the same interaction z from 200 expression-matched pseudo-clock sets (A166 construction). Clock-specific
only if |z| exceeds the 95th percentile of |z_pseudo|.
**Pass:** Holm p < 0.05 and clock-specific. Even then the reading is hypothesis-generating: treatment is not randomised
(confounding by indication), and a pass would need a randomised-trial cohort before any clinical statement.
**Secondary (not counted):** without the G2M x T term; METABRIC RFS; per cohort; within ER+ for H198.1; MKI67 instead of
G2M.
**Reading:** a fail with its MDE closes the last prognosis question that proliferation does not already answer, as far as
these cohorts can. Wording in the paper: "no evidence of an interaction larger than X", never "the clock does not affect
treatment response".
**Stop rule:** one amendment.
**Code:** `paper/exploratory/C8b_clock_prognosis_replication/scripts/score_a198.py`.

## A199 (2026-11-03, on the user's instruction): does CYCLOPS magnitude measure clock strength in samples taken in working hours? **GTEx fast deaths, primary (Holm over two)**

**Question.** Li et al. (PNAS 2024) read CYCLOPS 2.0 magnitude as tumour clock strength in surgical samples, which are taken
in daytime. A123b showed that untimed orderers, CYCLOPS included, lose the organ map on 08-18 h donors, but CYCLOPS also did
poorly with the whole day (offsets 0.115) and its donor order did not degrade (0.47 vs 0.50); magnitude was never tested.
**Seen before writing:** A123b orderings (phases only) and window / whole-day donor counts per tissue; no CYCLOPS magnitude
has been computed in GTEx.
**Data.** GTEx v10 fast deaths (Hardy 1-2) with recorded time of death (`_shared.load`), as A123 (user: fast deaths only;
CHIRAL's donor internal phase not used). Tissues: Breast - Mammary Tissue plus the 11 other tissues with most window donors,
Whole Blood excluded (no core-clock signal, post hoc 2026-11-03). Designs per tissue: **W** = donors who died 08-18 h;
**Fs** = the same number of fast deaths drawn at random from the whole day (seed 20261103), so W and Fs differ only in time
coverage.
**Fits.** CYCLOPS 2.0 as in A123b (`run_order` input: 2^log2-normalised - 1, CYCLOPS human seed list, template
hyperparameters, 80 models, no covariate), three sample-order seeds per (tissue, design); magnitude = mean over seeds of the
within-fit z-scored `Magnitude`.
**Timed clock strength (truth).** Per tissue, cosinor of the 14 GTEx clock genes (log2-normalised) on recorded time over
all fast deaths, 5-fold cross-fitted (the donor's own fold excluded); a donor's strength = projection coefficient of its
centred clock-gene vector on the expected profile at its recorded time, sum_g (y_g - m_g) p_g / sum_g p_g^2.
**H199.1 (primary):** magnitude tracks timed clock strength less in working hours. Per tissue rho = Spearman(magnitude,
timed strength) in W and in Fs; statistic = median over tissues of rho_Fs - rho_W; one-sided Wilcoxon signed-rank over
tissues.
**H199.2 (primary):** with whole-day coverage, magnitude tracks timed clock strength: median rho_Fs > 0, one-sided
Wilcoxon over tissues. (If H199.2 fails, magnitude does not measure clock strength in GTEx whatever the coverage, and
H199.1 is read in that light.)
**Power step first:** MDE = 2.49 x sqrt(mean over tissues of (1/(n_W - 3) + 1/(n_Fs - 3))) x 1.05 / sqrt(n tissues)
(Fisher-z approximation, Wilcoxon efficiency); above 0.20: not estimable.
**Pass:** Holm p < 0.05 for each.
**Secondary (not counted):** R^2 of magnitude on (cos, sin) of the hour of death, W vs Fs; Spearman of magnitude with
sequencing depth and with the first principal component of non-clock genes; seed concordance of magnitude per design;
CYCLOPS phase vs recorded time (|circular r| after rotation and reflection), descriptive only (range restriction makes it
not comparable between designs).
**Reading:** H199.1 and H199.2 pass: magnitude measures clock strength only with whole-day sampling, so its use on daytime
surgical samples is unsupported. H199.1 fails with H199.2 passing: coverage is not the problem for magnitude. H199.2 fails:
magnitude is not a clock-strength readout in GTEx at all. Wording in the paper stays at "in GTEx fast deaths"; tumours are
not tested here.
**Controlled access:** recorded time of death and per-donor values stay in memory or in the git-ignored output cache; only
per-tissue aggregates are written.
**Stop rule:** one amendment. A tissue whose fit fails in any seed or design is dropped and reported.
**Code:** `paper/exploratory/A199_cyclops_working_hours/`.

## A200 (2026-11-03, on the user's instruction): registered rerun of V25.5, is per-patient internal time in breast surgery clock-specific? **Seen data, primary (Holm over two)**

**Seen before writing:** the V25.5 result on these data (GSE233242, 42 patients, tumour + adjacent normal, resections
08-18 h): Spearman of the tangent coordinate with the recorded hour 0.558 [0.297, 0.772], median error 1.12 h with the
error-minimising offset, 1.42 h with the estimates' circular mean set to 13:00, against 2.10 h for "always 13:00" and
1.67 h for the best constant (`A25_prior_baseline`); PER1 removal unchanged (archive); A146 programs on the tumours.
Nothing in this amendment is fresh data: it fixes the analysis plan and adds the specificity control that V25.5 lacked.
**Estimator (unchanged from V25.5, `compute_breast_timed.py`):** two-compartment host (normal, tumour), A25 internal-time
estimate, 12 h flip resolved by the hospital prior, 30 M reads, pooled human waveform, primary join (42 patients).
**H200.1 (clock specificity):** the Spearman of the tangent coordinate with the recorded hour exceeds the 95th percentile of
the same statistic from 200 expression-matched pseudo-clock sets run through the identical pipeline (each clock gene
replaced by one of its 20 nearest non-clock-like, non-G2M genes in mean and SD of log2(TPM + 1) over the 84 libraries;
counts built as for the clock genes; the pseudo genes take the clock genes' waveform amplitudes and phases). One-sided p =
(1 + #pseudo >= clock) / (1 + #pseudo). Fits that fail are excluded and counted; above 20 % failures, H200.1 is
uninformative.
**H200.2 (no use of truth):** median absolute error with the estimates' circular mean set to 13:00 (the hospital prior,
no recorded time used) is below the permutation null (2,000 shuffles of the recorded hours, same rule); one-sided.
**Pass:** Holm p < 0.05 for each.
**Secondary (not counted):** without PER1 (glucocorticoid-induced at surgery); without PER1, PER2 and TSC22D3-responsive
genes is not applicable (14-gene set: drop PER1 only, plus PER1 and NR1D1); per compartment alone (normal, tumour); breast
waveform; error of the pseudo sets under the 13:00 rule (distribution); the pseudo-set percentile of the error.
**Reading:** H200.1 pass = the patient ordering by recorded surgery hour is carried by the clock, not by any matched gene
set that tracks the hour (fasting, anaesthesia, handling). H200.1 fail = the per-patient internal-time claim is withdrawn
from the main text (ordering may be generic time-of-surgery biology). Being on seen data, this is a robustness test, not a
replication; a fresh cohort with recorded resection times remains the decisive step.
**Stop rule:** one amendment.
**Code:** `paper/exploratory/A200_breast_surgery_rerun/`.

## A201 (2026-11-03, on the user's instruction): after a 6 h light advance, do aortic cell types keep their clock phase relations while losing alignment to the light cycle? **Mouse aorta scRNA-seq (Auerbach et al., bioRxiv 2026), primary**

**Data and what was seen:** Zenodo 22847203 (`$SCCC_DATA/aorta_jetlag_auerbach2026`, PROVENANCE.md). Metadata only:
32 libraries, each a pool of 2 aortas (paper Methods); aligned controls 14 (sex x ZT0/6/12/18, 1-3 per cell), misaligned
8 (one per sex x ZT, 5 days after a 6 h advance, ZT on the new schedule), ubiquitous Bmal1 KO (CAGG-CreERT2, tamoxifen)
10; cells per library and type: SMC 672-6,321, fibroblast 36-1,151, EC 34-332, macrophage 8-198. The paper reports reduced
amplitudes, limited adaptation of peak times and, within SMCs, Dbp and Nr1d1 about 1.2-1.3 h closer to Arntl. No
expression value has been read.
**Readout.** Per library x cell type (SMC, fibroblast; EC and macrophage secondary, used only where >= 30 cells in every
library of the comparison), pseudobulk counts of the 14 clock genes (mouse orthologs of the GTEx set), log(CPM + 1) with CPM
on the library's total counts of that cell type. Per condition and cell type, per gene: y = m_g + sex_g + a_g cos(w ZT) +
b_g sin(w ZT); complex amplitude z_g = a_g - i b_g. Phase shift of cell type c, misaligned vs aligned:
Delta_c = arg sum_g z_g^mis conj(z_g^ali) (hours, amplitude-weighted).
**H201.0 (gate, first):** the misalignment acted on the peripheral clock: mean of Delta_SMC and Delta_fib differs from 0
(90 % interval excludes 0). If it fails, H201.1 is uninformative (nothing to desynchronise).
**H201.1 (primary):** relations kept: D = Delta_SMC - Delta_fib lies within +/- 1.5 h (90 % interval inside the bounds,
TOST at alpha 0.05).
**Uncertainty:** the larger of (a) a cell bootstrap (cells resampled within library, 1,000) and (b) a library jackknife
(leave one library out, within condition); intervals as estimate +/- 1.645 SD.
**Power step first:** SD of D; estimable if 1.645 SD <= 0.75 h, else H201.1 is not estimable.
**Secondary (not counted):** amplitude ratio misaligned / aligned per cell type; Delta for EC and macrophage; Delta per sex;
clocksync untimed analysis (libraries as donors, cell types as compartments): gate per condition and arg z of fibroblast
vs SMC, aligned vs misaligned; Bmal1 KO as negative control: the gate should fail and the clock-gene amplitude collapse
(ratio KO / control of sum |z_g|).
**Reading:** H201.0 and H201.1 pass: five days after a shift, aortic cell types moved together, consistent with a shared
internal time that is misaligned with the light cycle rather than desynchronised internally. H201.1 fails with H201.0
passing: cell types re-entrain at different rates (internal desynchrony). One tissue, one time after the shift, pools of
two mice, n = 1 per sex x ZT: an existence test, not a body-wide statement.
**Stop rule:** one amendment.
**Code:** `paper/exploratory/A201_aorta_jetlag/`.

## A202 (2026-11-03, on the user's instruction after the claims-ledger audit; written before any simulation below; A164 scored and seen): is the A164 equivalence interval valid at its own design? **Method, confirmatory (V202.1-2)**

**Why:** H164.1-2 (T2D, obesity) are the paper's equivalence claims: the 90 % interval D ± 1.645 × bootstrap SE of the
pooled module d* lies inside ±0.15. Coupling intervals undercover elsewhere (A154 V154.3, A162, A163), and the coverage of
d* intervals was never checked. **Seen before writing:** `A164.../output/a164_power.csv` gives a pooled T2D bootstrap SE
of 0.0247 against a permutation-null SD of 0.0334 (ratio 0.74; vent 0.77, fast 0.72); obesity 0.0518 vs 0.0523. The
T2D interval may therefore be too narrow. A162 also reports an attenuation of c by about 0.8-0.9 at realistic anchor
coupling, which a thinning design cannot see (below).
**Design (real data, A164 matched sets, both cohorts, module tissues only):** truth injected by binomial thinning of the
real counts (`clocksync.power.damped`, as `power.mde`), with exchangeable labels:
1. Per replicate, labels are reassigned within every matched set exactly as the A164 permutations (one uniform key per
   donor, shared across tissues; the member with the largest key is the case).
2. Per module tissue k, the host-locked clock part of k is damped to s_k times itself in the donors labelled **case**
   (δ < 0) or **control** (δ > 0), with s_k chosen so that c0_k − c′(s_k) = |δ|, c0_k the pooled coupling of k on the
   labelled donors and c′(s) = s² c0 / (s² c0 + 1 − c0) (`power.mde`). Every module tissue gets the same δ. A tissue with
   c0_k < |δ| + 0.02 is damped to s = 0 and its realised δ_k recorded; the module truth is the mean of the realised δ_k.
3. d* per tissue as A164 (`Dataset.contrast`, one call per tissue on its leave-one-out host, untruncated); module mean;
   pooled over cohorts with the **A164 weights** (fixed, `a164_power.csv`).
- Each tissue is read on a host that excludes it, so damping all module tissues at once leaves each tissue's host
  unchanged (as `power.mde`). Under exchangeable labels both arms have the same host spread in expectation, so d* ≈ the
  injected change of c.
- Seeds `SeedSequence([202, design, cell, rep])`. Designs: T2D (1:3) and obesity (1:1), the primary A164 designs.

**Arm P (point estimates, no bootstrap):** cells δ ∈ {0, −0.075, −0.15, +0.15}, 1,000 replicates per cell and design.
Readouts: median bias of the pooled module D̂ minus the module truth; the SD of D̂ over replicates (SD_sim), against the
A164 bootstrap SE and null SD; the coverage of D̂ ± 1.645 × SE_A164 (the registered half-width, fixed at its observed
value); the TOST rate: share of replicates declared equivalent (interval inside ±0.15) at δ = ±0.15.
**Arm B (the registered procedure in full):** cells δ ∈ {0, −0.15}, 100 replicates per cell and design, each with the
A164 matched-set bootstrap (200 replicates instead of 500, host refitted; a cost deviation stated here) and its own
interval D̂ ± 1.645 × SE_boot.
**Validation criteria:**
- **V202.1 (coverage of the registered procedure, arm B):** coverage of the 90 % interval ≥ 0.85 in every cell × design
  (100 replicates: binomial SD 0.03 at 0.90).
- **V202.2 (equivalence size, arm P):** the TOST rate at δ = −0.15 and at δ = +0.15 is ≤ 0.075 in each design (nominal
  0.05; the A154 tolerance).

**Descriptive:** bias per cell; SD_sim / SE_A164 and SD_sim / null SD; coverage with SD_sim and with the null SD in place
of SE_A164; per-cohort results; arm B's distribution of SE_boot. **Attenuation (descriptive, not simulated):** the A164
intervals divided by the A162 attenuation range (0.80 and 0.90), read against ±0.15. This is arithmetic on seen numbers,
reported as a sensitivity, not a test.
**Reading:**
- V202.1-2 pass → H164.1-2 stand as tested equivalences, with "interval validated at the design (A202)"; the attenuation
  sensitivity is stated.
- V202.1 fails or V202.2 fails → the A164 equivalence is re-read (labelled as a re-reading after A202, not a new test)
  with the design-based half-width 1.645 × max(SD_sim, null SD) and, if the bias exceeds 0.02 in magnitude, with D̂ − bias.
  If that interval still lies inside ±0.15, H164.1-2 hold under the corrected interval; otherwise the claim drops to
  "no effect detected (MDE 0.083 / 0.130)".
- Thinning on real data defines the truth on the estimator's own scale, so attenuation of the true c (A162) is not
  covered by A202; it stays a stated limit.
**Compute:** arm P about 25 min on 16 workers; arm B about 4 h on 16 workers (A164 run: about 2.5 CPU-s per cohort evaluation of 8 tissues). Run
through `paper/common/run.sh`, arm P first; arm B only after arm P is logged.
**Stop rule:** one amendment; no new bound, design or test is added after either arm is seen.
**Code:** `paper/exploratory/A202_a164_coverage/` (reuses `A164.../scripts/compute_a164.py` helpers and `clocksync.power.damped`).

## A203 (2026-11-04, on the user's question whether single cells could be used better than pseudobulks; written before any subsampled fit is computed; seen: the A17 `lung_hlca` table (198 healthy donors, 18 compartments, median cells per pseudobulk 113-948 by compartment) and the per-cell file's columns): is the precision of internal-time readouts limited by cells or by donors? **Method, confirmatory (V203.1-2); enters the Methods only if both pass (user, 2026-11-04)**

**Why:** clock time varies between donors only, so donors should set the precision of coupling and offsets, and cells
only each donor's measurement error. If this holds at a real design it gives a design rule (how many cells per
pseudobulk are enough) and states why the method works on pseudobulks.
**Data:** HLCA per-cell clock-gene counts (`$SCCC_DATA/cells/hlca_full/cells.parquet`, Census 2025-11-08), regrouped
into the A17 `lung_hlca` healthy rows (same donors, compartments from `ann_level_3` / `ann_level_2` as
`common/datasets/_atlas/lung.py`, the donor's majority assay, strata study | assay, depth = summed `raw_sum`), A17
learned anchors (`common/anchors_a17.csv`). **Acceptance:** the regrouped full table reproduces the A17 `lung_hlca`
clock counts of every row (else the cause is reported and the scoring waits for a fix).
**Readout compartments:** those with ≥ 40 donors in the full table (rule fixed here; the others stay in the fit). Fit:
`Dataset.fit(gate=False, n_boot=0)`, LOO anchor readout (default): per compartment offset arg z (h) and coupling c.
**Error metric (split halves):** donors split at random into two halves, balanced within stratum; both halves fitted on
the same cell subsample rule; per compartment |offset₁ − offset₂| (circular, h) and |c₁ − c₂|. err = median over
readout compartments and 50 splits. Descriptive: circular correlation of offsets and Pearson of c between halves.
- **Cells arm:** all 198 donors; cells per pseudobulk capped at m ∈ {10, 25, 50, 100, 200, all} by random subsampling
  of cells without replacement (a new subsample per split and half).
- **Donors arm:** all cells; donors first subsampled to N ∈ {50, 100, 150, 198} (within strata), then split.
- Seeds `SeedSequence([203, arm, level, split])`.
**Validation criteria:**
- **V203.1 (cells saturate):** err at m = 100 ≤ 1.10 × err at m = all, for the offset and for c.
- **V203.2 (donors limit):** err at N = 100 ≥ 1.25 × err at N = 198, for the offset and for c (√2 ≈ 1.41 is expected
  if donors set the precision).
**Descriptive:** m₉₅ (cells at which err is within 5 % of err(all), interpolated); per compartment family; binomial
thinning of the full pseudobulks to the expected UMI of m cells (m = 25, 100) against true cell subsampling, which
shows whether cell-to-cell heterogeneity adds error beyond counting noise.
**Reading:**
- V203.1 and V203.2 pass → Methods: "beyond about m cells per pseudobulk, the precision of coupling and offsets is set
  by the number of donors" with the curves (Extended Data).
- Either fails → reported in `results_preliminary.md` only, not in the Methods (user, 2026-11-04).
- One atlas (lung); the rule is not claimed for other tissues or depths.
**Stop rule:** one amendment; no new metric or threshold after scoring.
**Code:** `paper/exploratory/A203_cells_vs_donors/`.

## A204 (2026-11-05, from the literature follow-up list, `paper/docs/prompts/literature_followup_analyses.md` item 1; written before any simulated readout; seen: one COFE and one TimeTeller run per design for runtime only, no output read): how do orderers produce "uncoupled tumour clocks" under working hours? **Simulation, secondary family A204 (Holm over four)**

**Why.** Two published claims of tumour-clock uncoupling rest on readouts that working-hours sampling may produce by
itself. COFE (Ananthasubramaniam & Venkataramanan, PLoS Biol 2025, TCGA): "predicted time labels lay throughout the
circadian cycle in all ACs" and "the difference between the predicted time labels for healthy and tumor tissue was
distributed throughout the circadian cycle". TimeTeller (Vasilyev et al., bioRxiv 2026, REMAGUS breast, n = 108):
"sampling time was not correlated with either TT time (slope -0.027, p = 0.80)". A145 measured coupling error of
two-step orderers (W arm, CoFE 0.206 vs clocksync 0.060), not the orderers' own outputs.
**Relation to A138 Q2 (H138.3).** H138.3 registered the same Rayleigh readout on a genome-scale simulator; its
simulator gate failed in a two-replicate check and Q2 is parked (ROADMAP 4.8). A204 uses the A145 clock-gene simulator,
on which CoFE recovers time with whole-day sampling (A145 D145.1). It does not score H138.3; that stays parked.
Limit stated now: COFE in TCGA used genome-wide features; A204 tests COFE on the 14 clock genes.

**Simulator** (`clocksync.simulate`, GTEx 14-gene waveform, depth, tau, shift and mesor_sd at their defaults; the
per-compartment gene mesor offsets play the tissue difference between tumour and normal). Per patient:
- normal N: z = 1, psi for c = 0.9, presence 0.2 (about 60 matched normals, as TCGA per cancer);
- tumour T: z = 0.7 exp(i 3 h), presence 1; truth levels c_T in {0.9, 0.5} (coupled) and "uncoupled" (`coupled=False`,
  psi as for 0.5).
psi is set from c at whole-day host variance (`A145 lib.psi_of`, V of 0-24 h) in every window: the biology is fixed and
the window changes, so population coupling is lower in W (stated, as CLAUDE.md "coupling scales with host spread").
n = 300 patients. Windows: F (0-24 h, jitter 0) and W (08-18 h, jitter 1 h). 30 replicates per window x truth (seed
`[204, window, truth, rep]`).
**COFE as in the paper:** `common.benchmark.cofe.run_cofe` (standardised features, 3 x 5-fold CV, 5 restarts) trained on
the 300 tumours (log2 normalised clock genes); the matched normals are test samples projected with COFE's own test rule
(centred and scaled by the tumour means and SDs, projected on the tumour loadings, angle = time label), as COFE's
train/test API does. Direction and origin of the labels are COFE's (arbitrary); the readouts below do not depend on them.
**clocksync on the same data:** the A138 H138.4 tumour-normal statistic on the matched patients:
`paired.cohort_states` (one projector, counts and depth) and `paired.cross` corr (Poisson-corrected modulus of the
complex correlation); null: normals permuted across patients (1,000), p = (1 + #null >= obs) / 1,001; offset arg.
**TimeTeller (REMAGUS design):** `common.benchmark.timeteller.run_timeteller`, trained on 144 timed whole-day healthy
samples (compartment z = 1, psi for c 0.9, its own mesor offsets: a tissue mismatch as healthy training vs tumour), tested
on 108 tumours (z = 0.7 exp(i 3 h)) with c_T in {0 (uncoupled), 0.25, 0.5, 0.75, 0.9}, windows F, W and W_narrow
(10-16 h, jitter 1 h; REMAGUS times "cluster around midday", range not reported: the window is an assumption). Readout as
the paper: linear regression of TT time (h) on sampling time (h), slope p (two-sided), plus the circular-linear
correlation. 30 replicates per cell.

**Simulator gate (checked first, from arm F; if it fails A204 stops and is reported as such):** COFE, F, c_T 0.9:
median |circular correlation| of tumour labels with true times >= 0.7; TimeTeller, F, c_T 0.9: median circular error of
TT time minus (true time + 3 h) <= 3 h.

**Hypotheses (secondary family A204, Holm over H204.1-H204.4; 30 replicates; "≥ 20 of 30" is the exact one-sided
binomial test against 0.5, p 0.049):**
- **H204.1 (labels span the cycle):** W, c_T 0.9: COFE's tumour labels leave no empty arc of 4 h or more (largest gap
  between consecutive sorted labels < 4 h) in >= 20 of 30 replicates. The true times' largest gap (about 12 h) is reported.
- **H204.2 (false uncoupling):** W, c_T 0.9: the COFE tumour-minus-normal label differences pass as uniform (Rayleigh
  p > 0.05) in >= 20 of 30 replicates (the H138.3 readout).
- **H204.3 (clocksync on the same data):** W, c_T 0.9: the pairing test gives p < 0.05 in >= 24 of 30 replicates, and
  the specificity companion on the same permutation scheme, W uncoupled, gives p < 0.05 in <= 3 of 30 (both required).
  Test statistic for Holm: the one-sided exact binomial p of the detection count against 0.5.
- **H204.4 (TimeTeller at the REMAGUS design):** W, c_T 0.9, n = 108: the slope test detects (p < 0.05) in >= 24 of 30
  replicates (power >= 0.8). Holm p: one-sided exact binomial of the count against 0.5.
**Descriptive (D204):** per cell, the resultant length R of COFE labels and of the true times; R and Rayleigh p of
the differences; COFE's interior CV minimum share; |circular correlation| of labels with true times; clocksync corr,
offset (truth 3 h) and detection rate; the discrimination AUC between c_T 0.9 and uncoupled in W of COFE's R_diff and of
clocksync's corr, on the same replicates; TimeTeller power per c_T and window and its MDE (smallest c_T with power >=
0.8); the share of TT times inside 07-20 h.

**Reading (fixed now):**
- H204.1 and H204.2 pass → under working hours, clock-gene COFE returns full-cycle labels and uniform tumour-normal
  differences for a coupled tumour: the two published observations are what COFE outputs whether or not the tumour is
  coupled, so they are not evidence of uncoupling. If H204.3 also passes, clocksync detects the coupling in the same
  patients. Genome-wide COFE is not tested.
- H204.2 fails (differences concentrated) → in this simulator COFE would have seen coupling under working hours; the
  TCGA observation is then not explained by the sampling window, and the paper does not say it is.
- H204.4 passes → a coupled tumour clock correlates with sampling time at REMAGUS size if TimeTeller times are as
  accurate as in this simulator; r ≈ 0 then is not explained by the narrow window alone (supports the authors' reading
  under the stated window). H204.4 fails → a narrow window with TimeTeller's error gives r ≈ 0 for a coupled tumour; the
  MDE is reported. W_narrow is a sensitivity, not a test.
**Stop rule:** one amendment; no truth level, window or threshold added after any readout is seen.
**Compute:** about 180 COFE and 450 TimeTeller runs of about 30 s each, 8 workers, about 1 h; through `paper/common/run.sh`.
**Code:** `paper/exploratory/A204_orderer_uncoupling/`.

## A205 (2026-11-05, literature follow-up item 4a; written before any simulation of it): does the gate pass, and where does the offset land, when the genes of one compartment shift by different amounts? **Simulation, method (V205.1-2), secondary**

**Why.** After feeding inversion, liver Bmal1-luc shifts about 6-7 h while Rev-erbα and PER2 shift about 12 h (Saini 2013,
Kiessling 2010 cards). clocksync assumes one shift per compartment (conserved waveform). METHOD.md tests scrambled
phases only; A139 tested refit and flags (both failed) but never computed the gate or the offset bias for this pattern
(its D139.3 was not run).
**Simulator:** A145 design (`clocksync.simulate`, GTEx 14 genes, defaults), 4 anchors (z 1, c 0.9) and one target
(z 1, c 0.9, offset 0) whose genes carry `gene_shift_h`:
- P0: none (reference); P12: all 14 genes +12 h (coherent shift); P_partial: ARNTL, CLOCK, NPAS2, CRY1 (the
  Bmal1 loop, about 6 h) +6 h, the other 10 genes +12 h; P_partial3: the same with +3 h / +12 h (stronger incoherence).
Windows F (0-24 h) and W (08-18 h, jitter 1 h); n = 80 donors; 100 replicates per cell (seed `[205, ...]`).
**Fit:** `Dataset.fit(n_boot=0, force=True)`, 4 anchors, leave-one-out host; gate (`passes`, clock share) and the target's
offset (relative to the anchors' circular mean), coupling and log relative amplitude.
**Validation criteria (W, the primary window):**
- **V205.1 (the gate tests a shared clock, not the waveform):** the gate passes in >= 80 % of P_partial replicates
  (anchors are intact, so a clock is present). Fail → the gate also rejects incoherent compartments; METHOD.md states it.
- **V205.2 (offset stays within the gene shifts):** the target offset lies in [6, 12] h in >= 90 % of P_partial
  replicates (the readout is an average of the gene shifts, not an artefact outside them).
**Descriptive:** per cell, gate pass rate, median offset (and its position between 6 and 12 h as a weight), coupling and
|z| against P0 and P12 (does incoherence read as lower amplitude or coupling?), P_partial3 the same.
**Reading:** V205.1-2 pass → METHOD.md: a compartment whose genes shift unequally passes the gate and returns an offset
between the gene shifts, with lower |z| (if seen); a measured offset is then not attributable to a single loop. Stop
rule: one amendment, no new scenario after scoring.
**Code:** `paper/exploratory/A205_partial_shift_gate/`.

## A204b (2026-11-06, post hoc to A204's failed gate, on the user's request; written before any rerun): TimeTeller at the REMAGUS design with the training tissue's mean levels matched
**Why (post hoc diagnosis, 2 replicates per case, scratch only, not a readout of any H204):** TimeTeller recovers time on
the A204 simulator (median error 0.94-1.21 h, whole day) when training and test samples share their per-gene mean
levels, and fails when they carry independent per-gene offsets (`mesor_sd` 0.5, the A204 setting). The A204 gate
therefore failed on the tissue mismatch, not on TimeTeller itself; the user reports TimeTeller works on GTEx 08-18 h
after training on the same tissue. The tumour offset convention (internal time = t + 3 h) was confirmed (1.2-1.3 h vs
5.7-6.4 h with the opposite sign).
**Design:** A204's TimeTeller arm unchanged except `mesor_sd = 0` for the training and the tumour compartments (tumour
z 0.7 exp(i 3 h) kept: amplitude and phase differ, mean levels do not). Same seeds, windows, c_T levels, 30 replicates.
**Gate:** as A204 (F, c_T 0.9, median error <= 3 h). **H204b.4:** as H204.4 (W, c_T 0.9, n 108: slope p < 0.05 in >= 24
of 30), single test. Readings as registered for H204.4; MDE per window reported. H204.1-H204.3 stay unscored (their
registration made the whole of A204 conditional on the joint gate); the COFE results remain descriptive.
**Code:** `compute_a204.py --only tt --mesor-sd 0 --out output/a204b`.

## A206 (2026-11-06, literature follow-up item 4b; written before any expression readout; seen: sample titles and CTs only): is a feeding-rhythmic but clockless liver clock-shaped? **Mouse, Cry1/Cry2 double knockout (Vollmers et al. 2009, GSE13093), secondary**
**Known at registration:** the cross-compartment gate did not fail in Bmal1 knockouts (A201 aorta: clock share 0.77 -> 0.36,
p 0.003; A129 hepatocyte KO liver passes). GSE13093 has one tissue, so the cross-compartment gate cannot run; the
single-set check of the project is `paired.reference_gate` (A52, A83: clock-plane share of the gene x gene covariance
of gene-centred log expression against 200 permuted-acrophase planes). GTEx human waveform, as for every mouse set.
**Data:** GSE13093 (`$SCCC_DATA/vollmers2009_GSE13093/`, parsed by `A206_cry_dko_gate/scripts/parse_gse13093.py`): liver,
log2 MAS5, 14 genes. Sets: Cry DKO restricted feeding (RF, 12 samples, CT24-46 every 2 h), Cry DKO ad libitum (AL, 12,
same CTs), WT RF hourly CT24-47 split into even CTs (12, the KO design) and odd CTs (12).
**Power companion (must pass for H206.1 to be read):** both WT RF halves pass reference_gate (p < 0.05).
- **H206.1:** Cry DKO RF fails reference_gate (p >= 0.05). Pass → a feeding-driven rhythm without a clock does not
  look clock-shaped to the single-set gate. Fail → feeding alone produces a clock-shaped covariance in liver; the gate
  cannot separate a fed rhythm from a clock in a single tissue (stated in METHOD.md).
**Descriptive:** clock share and p of the four sets; Cry DKO AL (expected no rhythm); per-gene cosinor amplitude on
the known CTs (which genes carry the fed rhythm), as context only.
**Stop rule:** one amendment, 200 permutations (the A52 setting) and seed 206; no other set or gene subset added after
the readout.
**Code:** `paper/exploratory/A206_cry_dko_gate/`.

## A207 (2026-11-06, literature follow-up item 6; written before any CCD or Grabe readout on TCGA; seen: A138 H138.4 per-cancer pairing results, scored 2026-10-27): do co-expression and order-based "tumour clock disruption" readouts coexist with the tumour-normal clock pairing? **TCGA pairs, secondary (Holm over two)**
**Why.** Shilts 2018 (ΔCCD) and Grabe 2024 (tumours ordered by the paired normal's inferred time; R² ARNTL 0.76 normal vs
0.11 tumour in LUAD) read loss of co-expression or order as clock disruption. A138 H138.4 found that tumour clock states
follow the same patient's normal (pooled over 10 cancers, p 1e-4). Nothing in the project has run CCD or the Grabe
procedure on these pairs (`common/benchmark/ccd.py` only used in the GTEx benchmark).
**Data and code reused:** the A138 table (`figure_5/scripts/cancer/_shared.py` load_tcga / pair_frames; 10 cancers,
517 pairs; log expression of the 12 CCD genes and the 14 GTEx genes), `common.benchmark.ccd.ccd` (Shilts reference,
12 genes), H138.4 per-cancer statistics (`A138_cancer_estimability/output/score/per_unit.csv`, not recomputed).
- **ΔCCD** per cancer = CCD(tumours) − CCD(normals of the same patients). Null: tumour/normal labels swapped within
  patient (1,000 swaps, the H138.4 scheme); pooled = mean over cancers weighted by pairs, the same swaps jointly.
- **Grabe procedure** (on the 14 GTEx clock genes; Grabe's 13-gene list is only partly reported): PCA of the normals'
  centred log expression, phase = atan2(PC2, PC1); per gene, cosinor R² on that phase in the normals and in the paired
  tumours; per cancer, ΔR² = median over genes (tumour − normal). Same swaps; pooled as ΔCCD.
- **H207.1:** pooled ΔCCD > 0 (one-sided swap p < 0.05).
- **H207.2:** pooled ΔR² < 0 (one-sided swap p < 0.05). **Specificity companion (on the same swaps):** the same ΔR² for
  200 expression-matched pseudo-clock gene sets (A52 rule, `paired.pseudo_waveforms`, genome counts from
  `tcga_paired/build.py`); H207.2 counts only if the clock ΔR² is below the 5th percentile of the pseudo-set ΔR².
**Descriptive:** per cancer ΔCCD, ΔR², the A138 corr and its p, side by side; Spearman across the 10 cancers (n = 10,
descriptive only).
**Reading:** H207.1 and/or H207.2 pass together with H138.4 → co-expression or order disruption and tumour-normal clock
pairing coexist in the same tumours: these readouts do not measure coupling to the patient's clock, and "disrupted"
must not be read as "uncoupled". Both fail → no disruption at this design; reported as such.
**Stop rule:** one amendment; no gene list change after the readout. **Code:** `paper/exploratory/A207_ccd_grabe_tcga/`.

## A208 (2026-11-06, literature follow-up item 7; written before any sex-stratified fit; seen: donor counts by sex only, A99 archive result): do women and men differ in organ coupling or offsets? **GTEx, A164 procedure with sex as the label, secondary**
**Why.** Chronotherapy shows sex x dosing-time interactions (Lévi 2010). A99 (archive) found no sex difference in
radius (p 0.36) and a timed phase +0.35 h (SE 0.31, female minus male, post hoc). Coupling and offsets by sex at a common
host spread were never estimated.
**Data:** `common.gtex.load` (public sex), the two A164 cohorts kept separate: fast (Hardy 1-2, with time of death; 187 men,
57 women) and ventilated (Hardy 0; 242 / 139). Tissues: >= 20 donors of each sex in the cohort, sex-specific organs out.
**Procedure:** `compute_a164.py` unchanged except the label (female = case, male = control): matching within cohort on
Hardy class and age bracket (1:1, women as the smaller group; `_covariates` helpers), `Dataset.contrast` d* per tissue,
host refitted per matched-set bootstrap (200), label permutations within matched sets (1,000), pooled module d* =
inverse-variance mean over tissues, then over cohorts; interval half-width from max(bootstrap SE, null SD) (A202,
`Contrast.equivalence`). Offset: the same for `offset_diff_h` (circular differences, pooled as a weighted mean).
**Power first (dated note before any d* is read):** MDE_80 = 2.49 x null SD of the pooled d* and of the pooled offset.
- **H208.1 (coupling):** if MDE <= 0.15, the pooled d* (female − male) lies within ±0.15 (90 % interval, as A164).
  If MDE > 0.15: not powered, d* reported descriptively only.
- **H208.2 (offset):** pooled offset difference two-sided permutation p < 0.05; secondary reading against A99's
  +0.35 h (same sign or not). Equivalence within ±1 h reported if the interval allows.
**Specificity companion:** on the same permutations, the pooled d* of the 3 lowest-coupling tissues (expected noise) and
expression-matched pseudo-clock gene sets for the offset (A52 rule, 50 sets) — the clock statistic must exceed their 95th
percentile for a difference to be called clock-specific.
**Reading:** H208.1 pass → organ coupling does not differ by sex (within ±0.15); a sex x schedule interaction would then
not come from coupling. H208.2 pass and companion → a clock-specific sex offset (direction and size stated, no
per-donor phases released).
**Stop rule:** one amendment, no tissue subset chosen after the readout. **Code:** `paper/exploratory/A208_sex_coupling/`.

**A208 power note (2026-11-06, before any observed d* or offset; from `A208_sex_coupling/output/power.csv`, 1,000
permutations).** 55 tissue x cohort units (vent 24, fast 31), 2,981 matched pairs. MDE_80: pooled d* 0.032 (vent 0.041,
fast 0.050) <= 0.15, so H208.1 is powered; pooled offset 0.17 h (vent 0.23, fast 0.24). Implementation choices (agent,
docstring of `compute_a208.py`): pooling weights from the permutation null SD (independent of observed values); tissues
with < 10 matched pairs dropped (none were); sex-specific organs excluded: ovary, uterus, vagina, cervix, fallopian
tube, prostate, testis and breast (breast is sex-differentiated in GTEx, excluded with them); age matched on the public
bracket midpoint.

## A209 (2026-11-06, user request after the novelty assessment; written before any A209 run; seen: one arm-N replicate (rep 0) fitted for runtime: gate failed, and one pseudo-set fit gave target couplings about 0.85-0.91): what is the false-claim rate of the gate plus the matched-gene companion, the rule the real-data claims use? **Simulation, method, secondary (Holm over three)**
**Why.** A145 H145.3 failed: in a no-clock world with a shared non-clock factor, a claim counted as gate pass and c > 0.2
occurred at 0.127 (CoFE 0.06). Real-data claims also require the clock statistic to beat expression-matched pseudo-clock
sets (and responsive sets) on the same data. The false-claim rate of that combined rule was never simulated.
**Simulator (= A145 `lib.cohort`, unchanged design):** 100 donors, 4 anchors (c 0.9), 5 targets (c 0, 0.25, 0.5, 0.75,
1.0; A145 `design(rep)` offsets and amplitudes), GTEx 14-gene waveform extended by 200 non-clock genes (amplitude 0;
mesors resampled from the clock-gene mesors plus N(0, 0.3), so pseudo sets are expression-matched by construction).
30 replicates per arm, seed 209.
- **Arm N:** A145 arm N: amplitude 0, whole day, shared program iid per gene with sd = A145 `shared_sd("large")`.
- **Arm N_resp:** arm N plus a responsive factor: per donor f ~ N(0, 1), common to its compartments, added with loading
  2 x shared_sd("large") to 5 clock genes (PER1, PER2, NR1D1, DBP, TEF) and to a responsive pool of 30 of the 200
  extra genes; the 9 other clock genes and 170 extra genes get none. (`clocksync.simulate` gains an optional per-gene
  shared loading vector; default unchanged.)
- **Arm W:** A145 arm W (clock present, 08-18 h, jitter 1 h, no shared factor), for power.
**Companion (per replicate, same data):** 100 matched pseudo sets (14 genes drawn from the non-responsive extra genes)
and, in N_resp, 100 responsive sets (5 genes from the responsive pool in the responsive positions plus 9 non-responsive
extras); each set renamed to the 14 clock genes (waveform carried, A52 rule) and fitted with the same anchors, gate off.
**Claim rule (as in real data):** a target is claimed if the clock-gene gate passes, c > 0.2, and its c exceeds the 95th
percentile of its own c over the matched sets (N, W), and over both matched and responsive sets (N_resp).
- **H209.1 (arm N):** false-claim rate over 150 targets x replicates <= 0.05 and below A145's 0.127 (one-sided exact
  binomial p < 0.05 against 0.127).
- **H209.2 (arm N_resp):** the same with both companions.
- **H209.3 (arm W, power):** among targets with true c >= 0.5 (90), the claim rate with the companion is >= 0.8 x the
  gate-and-c rate without it.
**Descriptive:** gate-only and matched-only claim rates in every arm; per-c claim rates; the companion percentile
distribution.
**Reading:** H209.1-2 pass -> the rule used for real-data claims is calibrated against shared non-clock factors,
including clock-gene-biased ones; H209.3 says what it costs. A failure is reported as the rule's false-claim rate.
**Stop rule:** one run; loadings, gene lists and the claim rule fixed here. **Code:** `paper/exploratory/A209_gate_companion/`.

## A210 (2026-11-06, user request; written before any copy-number file is read; seen: tumour coherence vs proliferation in METABRIC (A166) and TCGA (A97, A171, A172), MYC never examined): does MYC amplification lower tumour clock coherence at matched proliferation? **METABRIC primary, TCGA secondary (Holm over two)**
**Why.** Proliferating adenocarcinomas lose clock coherence (A127, A166, A171), but physiological proliferation does not
reproduce it (A126) and no cell-cycle locking is seen (A187b). MYC represses the clock directly (Altman et al. 2015,
card missing). Copy number is a natural experiment that separates an oncogenic programme from proliferation as such.
**Data:**
- **METABRIC:** `data_cna.txt` from the existing cBioPortal tarball (`$SCCC_DATA/metabric/`, GISTIC thresholded: 2 =
  amplification, 1 = gain, 0 = neutral); expression, PAM50, scores and pseudo sets exactly as A166.
- **TCGA:** UCSC Xena `TCGA.PANCAN.sampleMap/Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes` (public), joined on
  the 15-character primary-tumour barcode; expression as A172 (recount3 loader) for the 11 A97 and 12 A172 eligible
  projects.
**Definitions:** amplified = MYC GISTIC 2; neutral = 0; gains (1) and losses out of the tests. Proliferation (G2M mean
z) and coherence as A166 / A172 (A87 `decompose.components`, groups centred).
- **H210.1 (METABRIC, primary):** cells PAM50 x proliferation tertile (tertiles within PAM50) with >= 10 amplified and
  >= 10 neutral tumours; T = mean over cells of coherence(amplified) − coherence(neutral); MYC labels permuted within
  cell (2,000); one-sided for a loss.
- **H210.2 (TCGA, secondary):** cells project x proliferation tertile (within project), same rule and statistic.
**Specificity companion (same permutations):** 200 pseudo-clock sets (A166 construction in METABRIC, A172 in TCGA). A
pass is "clock-specific" only if T is below their 5th percentile.
**Power, before the tests (dated note):** MDE = 2.49 x null SD. A hypothesis with MDE > 0.25 is not scored (reported as
not estimable). Reference sizes: METABRIC ER− −0.229, proliferation −0.123 (A127).
**Descriptive:** amplified minus neutral proliferation within cells (residual imbalance); gain vs neutral; MYC expression
top vs bottom tertile within PAM50 x proliferation cells; per-cell T.
**Reading:** pass, clock-specific -> an oncogenic programme lowers coherence beyond proliferation; fail at a powered
design -> MYC amplification adds nothing detectable beyond proliferation (bounded by the MDE).
**Stop rule:** one amendment; no cell rule or threshold change after the power note. **Code:** `paper/exploratory/A210_myc_coherence/`.

## A211 (2026-11-06, user request; written before the expression files are downloaded; seen: GEO sample titles and design text only): does synchronous physiological proliferation swamp the liver clock? **Mouse partial hepatectomy, timed (GSE180012), secondary**
**Why.** Tumour clock coherence falls with proliferation (A127, A166, A171). If dividing alone swamps the clock, the
regenerating liver after a two-thirds partial hepatectomy (PHx), whose hepatocytes enter S phase together, should show
it against sham-operated liver at the same times. A126 (endometrium) did not, at MDE 0.28.
**Data:** GSE180012: liver every 4 h from 4 to 32 h after PHx or sham (3 replicates per point, per titles), plus PHx in
the morning or evening read at 4 and 24 h; processed files `GSE180012_processed_data_file_{1,2}.xls.gz` (parse only).
Main series only for the tests; the morning/evening series descriptive.
**Manipulation gate (both required, else uninformative):** (a) proliferation (Hallmark G2M mean z, clock genes removed)
higher in PHx than sham over 20-32 h (one-sided label-swap p < 0.05 within time point); (b) the sham clock is visible:
`paired.reference_gate` on the sham samples p < 0.05 (200 permuted-acrophase planes).
- **H211.1:** coherence(PHx) − coherence(sham) over the matched time points (A87 gene coherence of the 14 clock genes,
  each arm centred within time point is NOT applied: samples centred per arm, so the time-of-day rhythm is the signal);
  null: arm labels swapped within time point (2,000); one-sided for a loss.
**Specificity companion (same swaps):** 200 expression-matched pseudo-clock sets (A168 construction). "Clock-specific"
only if T is below their 5th percentile.
**Power (dated note before the test):** MDE = 2.49 x null SD; if > 0.3, H211.1 is reported as not estimable.
**Descriptive:** the same split into 4-16 h (injury, pre-replicative) and 20-32 h (replicative); clock-gene cosinor
amplitudes and phases per arm over the 28 h series; the morning/evening series.
**Reading:** H211.1 passes, clock-specific, with the loss in 20-32 h -> synchronous proliferation can swamp a clock in
normal tissue (the tumour loss need not be malignant); a powered fail with the gate passed -> dividing alone does not
swamp the liver clock, supporting a tumour-state explanation (with A126). PHx is also an injury: an early loss is read as
injury, not proliferation.
**Stop rule:** one amendment. **Code:** `paper/exploratory/A211_liver_regeneration/`.

**Not registered (2026-11-06): SCN-lesion prediction of the brain-death reading.** No public dataset with several
peripheral organs from the same SCN-lesioned (or SCN-clock-deleted) animals was found (E-utilities searches
"suprachiasmatic lesion", "SCN lesion", "SCNx"; Sonnet scout). Single-tissue SCNx series exist (GSE243858 choroid plexus,
12 control and 12 SCNx mice every 4 h in DD; GSE300640 brown fat, ZT4/16) but cannot test within-animal inter-organ order.

**A211 gate and power note (2026-11-06, before any observed coherence; `A211_liver_regeneration/output/{gate,power}.csv`).**
Data: file 1 holds log2(RPKM + 1) for 24 sham and 24 PHx samples (3 per time point, 4-32 h; no raw counts, so Y is ln 2 x
log2 expression as A168 arrays); surgery 08:00-09:00 (GEO extract protocol). Gate (a) passes: G2M PHx − sham over
20-32 h +1.05 z, p 0.0005 (2,000 swaps; −1.11 at 4 h rising to +1.39 at 32 h); gate (b) passes: sham clock share 0.673
vs 0.077, p 0.005. MDE (2,000 swaps within time point): registered 4-32 h 0.064; descriptive 4-16 h 0.137, 20-32 h
0.040. H211.1 is estimable.

**A210 power note (2026-11-06, before any observed coherence difference; `A210_myc_coherence/output/{eligibility,counts,power}.csv`).**
METABRIC: 14 of 15 PAM50 x proliferation cells eligible, 1,294 tumours (450 amplified); TCGA: 30 project x tertile cells
in 13 projects, 2,728 tumours (844 amplified). MDE (2,000 permutations): H210.1 0.130, H210.2 0.086; both scored. The
null means are −0.020 and −0.044 (smaller amplified groups read lower coherence), so p values come from the permutation
null, not from T < 0. Choices left open by the registration (README): tertiles over all CNA-called tumours of the
stratum; Holm over the scored hypotheses; METABRIC pseudo sets = the 200 stored by A166, TCGA pseudo sets per A172
(seed crc32("A210")).

## A209b (2026-11-06, post hoc to A209 on a simulator asymmetry found after its scoring; written before any A209b run): the A209 false-claim test with non-clock genes exchangeable with the clock genes
**Why.** In A209 arms N and N_resp the clock amplitude is 0, but the clock genes still receive each compartment's own
clock-plane noise (psi) through their waveform amplitudes; the 200 extra genes had amplitude 0 and received none, so
pseudo sets looked more coupled than the clock genes (median clock percentile 0) and the companion was strict by
construction. In N_resp the gate never passed (0 of 30), so H209.2 did not test the companion.
**Change (only this):** the extra genes get amplitudes resampled from the 14 clock-gene amplitudes and acrophases
uniform on the cycle (seed 209), so in the no-clock arms they carry the same compartment noise and shared programs as the
clock genes; amplitude stays 0 (no rhythm). Arms: N_b (= N), N_resp_b (= N_resp, loading 2 x shared_sd("large")) and
N_resp1_b (loading 1 x shared_sd("large"), added because the 2x factor made the gate fail every time). Arm W is not
rerun (A209 H209.3 stands: pseudo genes there have no rhythm either way). Same 30 replicates, pseudo sets, claim rule.
- **H209b.1 (N_b):** false-claim rate (gate, c > 0.2, above matched q95) <= 0.05 and below 0.127 (one-sided exact
  binomial p < 0.05).
- **H209b.2 (N_resp_b and N_resp1_b pooled, 300 targets):** the same with both companions.
**Descriptive:** gate-only, matched-only (without gate) and full rates per arm; median clock percentile among pseudo sets
(expected about 0.5 under exchangeability: the check that the fix works).
**Reading:** pass -> the real-data rule is calibrated against shared non-clock factors when pseudo genes are
exchangeable; fail -> the rule's false-claim rate is reported as such. **Code:** `compute_a209.py --exchangeable`.

## A212 (2026-11-06, user request; written before any A212 code is run; seen: the A176 components of all seven contrasts, the A210 and A211 results): can bulk coherence data tell a damped or dephased clock from a rewired one or a swamped one? **Method, plasmode simulation**
**Why.** A176 H176.1 shows less in-plane clock signal in proliferating and ER− tumours (7 of 7), but H176.2, the test
for disordered relations, was degenerate, so "swamped, not disordered" was not shown. The reference waveform is shared
by all tumours. Clonal heterogeneity can lower coherence in two ways the pseudo-clock companion cannot catch, because
both act on the clock genes only: clones out of phase with each other (bulk damping) and clones whose gene relations
differ from the reference (rewiring). Within-tumour dephasing and lower per-cell amplitude are the same thing in bulk
to first order (amplitude factor |E e^{i delta}|), so they are one mechanism here (D).
**Bases (real tumours, no mechanism added):** B1 = the reference arm of A176 H127.1 (SCAN-B, lowest proliferation
tertile of each PAM50 x purity stratum); B2 = the reference arm of H166.1 (METABRIC ER+, A166 strata). Y as in A176,
centred within stratum. Per tumour, clock state beta = Y X (X'X)^-1 (GTEx waveform X), clock part C = beta X', rest
E = Y − C. Each replicate splits every stratum's tumours at random into halves; the case half gets the mechanism, and
the A176 statistics (dR, dlog_plane, dlog_off, ds; `decompose.components`, mean over strata) are case minus reference.
**Mechanisms (case half only):**
- **D, damping / dephasing:** Y' = rho C + E, rho in {1, 0.9, ..., 0.1}.
- **S, swamping:** Y' = Y + U L', three non-clock factors, U ~ N(0, 1) per tumour, loadings L ~ N(0, s^2) per gene
  (redrawn per replicate), s on a grid of 10 values from 0 to 2 x the per-gene SD of E.
- **W, homogeneous rewiring:** Y' = beta X_w' + E, X_w the waveform with each gene's acrophase shifted by N(0, sigma^2)
  (one draw per replicate, shared by all case tumours), sigma in {0, 1, ..., 6} h.
- **H, heterogeneous rewiring (the user's clonal scenario):** as W with an independent shift vector per tumour.
**Strength matching:** for each base and mechanism, 50 replicates per grid value; the strength whose mean dR equals the
observed dR of the base's contrast (B1 −0.123, B2 −0.229; A176 `tests.csv`) is found by linear interpolation. If no
grid value reaches it, the mechanism is recorded as "cannot reach the observed loss" and is run at its strongest value.
Then 200 replicates per mechanism at the matched strength (seeds SeedSequence([crc32("A212"), base, mechanism, rep])).
- **V212.1 (primary, B1):** rewiring is separable from damping at matched dR: cross-validated (5-fold) linear
  discriminant AUC on (dlog_plane, dlog_off) >= 0.9 for both W vs D and H vs D.
- **V212.2 (B1):** swamping is separable from damping: AUC S vs D >= 0.9.
- **V212.3 (B2):** V212.1 and V212.2 repeated on B2.
**Predictions written now:** D lowers dlog_plane with dlog_off near 0; S raises dlog_off; W and H lower dlog_plane and
raise dlog_off, and lower dR at a smaller plane loss than D.
**Descriptive (compared with seen data, labelled as such):** the observed (dlog_plane, dlog_off) of H127.1 against the
B1 clouds and of H166.1 against B2 (Mahalanobis distance to each cloud and its percentile within the cloud); mixtures are
not fitted.
**Reading:** V212.1 passes and the observed point sits in the D cloud only -> the coherence loss is less clock amplitude
in bulk (lower per-cell amplitude or intra-tumour dephasing), not rewired relations; observed in the W/H cloud -> rewired
relations; V212.1 fails -> bulk data cannot tell rewiring from damping, and the paper says "less clock signal", cause
open. No mechanism that fits is read as proven: the clouds bound what bulk statistics can distinguish.
**Stop rule:** one amendment; grids and AUC bar fixed here. **Code:** `paper/exploratory/A212_coherence_mechanisms/`.

## A213 (2026-11-06, user request; written before the ABSOLUTE table is downloaded and before any genome-altered fraction is computed; seen: the A210 results and the header and first two rows of the ABSOLUTE file): is tumour coherence lower in heterogeneous or copy-number-burdened tumours at matched proliferation? **TCGA primary, METABRIC secondary**
**Why.** Heterogeneity and genome burden go with proliferation, MYC amplification and outcome. If they lower coherence
at matched proliferation, they are a route by which a shared-waveform model and survival could both depend on tumour
heterogeneity, and they could carry the A210 TCGA result (MYC expression showed no loss).
**Data:**
- **Subclonal genome fraction** (ABSOLUTE, PanCanAtlas `TCGA_mastercalls.abs_tables_JSedit.fixed.txt`, GDC
  api.gdc.cancer.gov/data/4f277128-f793-4354-a13d-30cc7fe9f6b5; `call status` = called only), joined on the
  15-character barcode. Purity and ploidy from the same file for descriptives.
- **Fraction of genome altered (FGA):** fraction of genes with a GISTIC call != 0 in the copy-number matrices already used
  by A210 (Xena TCGA; METABRIC `data_cna.txt`), over genes called in all samples of the matrix.
- Expression, projects, proliferation and coherence exactly as A210.
**Cells:** stratum x proliferation tertile as A210 (TCGA project; METABRIC PAM50), all tumours with the measure. Within
each cell, the measure's top vs bottom tertile; a cell needs >= 10 tumours per arm. T = mean over cells of
coherence(top) − coherence(bottom); labels permuted within cell (2,000); one-sided for a loss.
- **H213.1 (primary):** TCGA, subclonal genome fraction.
- **H213.2:** TCGA, FGA. **H213.3:** METABRIC, FGA. Holm over the three.
**Specificity companion (same permutations):** 200 pseudo-clock sets (A210 construction; seeds crc32("A213")). A pass is
"clock-specific" only if T is below their 5th percentile.
**Power (dated note before the tests):** MDE = 2.49 x null SD; a hypothesis with MDE > 0.25 is not scored.
- **S213.4 (sensitivity of A210 H210.2):** within each A210 TCGA cell, amplified tumours matched 1:1 to neutral ones on
  FGA (nearest neighbour without replacement, caliper 0.05; unmatched dropped); T_MYC on the matched set with labels
  swapped within pairs (2,000). Read: FGA carries the A210 result if (T_matched − its null mean) is less than half of the
  A210 excess (T − null mean = −0.085) and p >= 0.05.
**Descriptive:** Spearman of subclonal fraction with FGA and with proliferation; residual proliferation and purity
imbalance top − bottom within cells; per-project T; the A176 components (dlog_plane, dlog_off, ds) of each H213
contrast, placed against the A212 clouds.
**Reading:** pass, clock-specific -> heterogeneous or burdened tumours carry less clock signal beyond proliferation, a
shared cause candidate for coherence and outcome; powered fail -> no such loss within the MDE, and heterogeneity does
not explain the coherence results at matched proliferation. Neither result says that heterogeneity changes outcome
through the clock: no survival test is registered here (A191-A198 found no clock information beyond proliferation).
**Stop rule:** one amendment; no cell rule or tertile change after the power note. **Code:**
`paper/exploratory/A213_heterogeneity_coherence/`.

**A213 power note (2026-11-06, before any observed coherence difference; `A213_heterogeneity_coherence/output/{eligibility,counts,power}.csv`).**
TCGA: 69 project x proliferation-tertile cells in all 23 projects for both measures (all cells eligible); H213.1
(subclonal genome fraction, 8,216 primaries with an ABSOLUTE "called" value) 5,499 tumours in the two arms, MDE 0.058
(null SD 0.0234); H213.2 (FGA, 9,005 primaries with a call) 6,024 tumours, MDE 0.058 (null SD 0.0235). METABRIC H213.3
(FGA, 1,756 tumours): 15 PAM50 x tertile cells, 1,176 tumours, MDE 0.136 (null SD 0.0545). Null means are about 0 (−0.0001,
−0.0004, 0.0003) because the arms are equal-sized tertiles. All three are scored (MDE <= 0.25). Choices left open by the
registration (README): FGA over genes complete in all samples (TCGA 24,776; METABRIC 22,079); tertiles by rank with ties in
table order (18.7% of called tumours have subclonal fraction 0); proliferation tertile over tumours with the measure; one
pseudo-set family per TCGA project (all primaries) shared by H213.1 and H213.2, METABRIC the 200 stored A166 sets; S213.4
cells need >= 10 matched pairs, pairs taken in increasing FGA distance, swap null within pairs.

## A215 (2026-11-06, user agreement; written before any mutation file is downloaded or read; seen: A171 GSE72094 LUAD descriptive "TP53 within purity −0.132, p 0.062", the A210, A211, A213 and A197 results, the A214 scout): does checkpoint-deregulated proliferation lower tumour clock coherence? **L1, TCGA primary, METABRIC secondary**
**Why.** Physiological proliferation keeps a coherent clock (A211 regenerating liver, powered; A126 endometrium), yet
proliferating tumours carry less in-plane clock signal (A176), and the clock signal adds no survival information beyond
proliferation (A191-A194). In physiological proliferation the RB and p53 checkpoints are intact, and the clock gates the
cycle through them (CDKN1A, WEE1). If the tumour loss belongs to proliferation that has escaped checkpoint control rather
than to division itself, then (i) tumours with p53 or RB-pathway loss carry less clock signal at matched proliferation,
and (ii) proliferation lowers coherence mainly where checkpoints are lost. The LUAD descriptive above is known, so the TCGA
LUAD project is not blind; the pan-cancer and METABRIC tests are.
**Data.** TCGA MC3 public MAF (`mc3.v0.2.8.PUBLIC.maf.gz`, GDC PanCanAtlas, open access) and the A210 GISTIC thresholded
matrix; METABRIC `data_mutations.txt`, `case_lists/cases_sequenced.txt` and `data_cna.txt` from the stored cBioPortal
tarball. Expression, projects, proliferation score and coherence exactly as A210 / A213 (TCGA 23 projects; METABRIC PAM50).
**Definitions (sequenced tumours only: TCGA barcode present in the MAF, METABRIC in `cases_sequenced`):**
- TP53 altered: >= 1 non-silent TP53 mutation (Missense, Nonsense, Frame_Shift_Del/Ins, In_Frame_Del/Ins, Splice_Site,
  Translation_Start_Site, Nonstop) or GISTIC −2; intact: none of these.
- RB pathway altered: a non-silent RB1 or CDKN2A mutation, or GISTIC −2 at RB1 or CDKN2A; intact: none.
**Cells:** stratum (TCGA project; METABRIC PAM50) x proliferation tertile (within stratum, over tumours with the status);
a cell needs >= 10 tumours per arm. T = mean over cells of coherence(altered) − coherence(intact); status labels permuted
within cell (2,000); one-sided for a loss. Arms are unequal, so the null mean is not 0: the excess T − null mean is
reported and p is taken against the permutation null.
- **H215.1 (primary):** TCGA, TP53. **H215.2:** METABRIC, TP53. **H215.3:** TCGA, RB pathway. Holm over the three.
- **H215.4 (secondary, interaction; TCGA TP53):** within project, I = [coherence(top) − coherence(bottom proliferation
  tertile)] among altered minus the same among intact (tertiles within project x status, >= 10 tumours per tertile arm in
  both statuses); mean over projects; status permuted within project x proliferation tertile (2,000); one-sided I < 0.
**Specificity companion (same permutations):** 200 pseudo-clock sets (A213 construction; seeds crc32("A215")); a pass is
clock-specific only if T (or I) is below their 5th percentile.
**Sensitivity S215.5:** H215.1 on 1:1 proliferation-matched pairs within cells (nearest neighbour on the continuous
proliferation score, caliper 0.1 SD of the project's score, without replacement); swap null within pairs.
**Power (dated note before the tests):** MDE = 2.49 x null SD; a hypothesis with MDE > 0.25 is not scored.
**Descriptive:** residual proliferation and purity imbalance altered − intact within cells; per-project T; truncating vs
missense TP53; the A176 components (dlog_plane, dlog_off, ds) placed against the A212 clouds.
**Reading:**
- H215.1 (+ H215.2) pass clock-specifically and H215.4 passes: clock loss goes with checkpoint-deregulated proliferation,
  not with division; this reconciles A211 / A126 with A176.
- H215.1 passes, H215.4 fails: TP53 loss is a state marker with its own clock loss, not a modifier of the proliferation
  effect.
- Powered fails of H215.1-3: checkpoint loss is not the shared cause within the MDE.
No survival test is registered. **Stop rule:** one amendment; no definition, cell or tertile change after the power note.
**Code:** `paper/exploratory/A215_checkpoint_coherence/`.

## A216 (2026-11-06, user agreement; written with A215, before any mutation file is read; seen: the proliferation losses of A97 / A171 / A172 / A176, A210, A213): does proliferation still lower tumour coherence at matched tumour state? **L1, TCGA primary, METABRIC secondary**
**Why.** If tumour proliferation is a marker of an aggressive state rather than the cause of clock loss (A211, A126), its
association with coherence should shrink when that state is matched.
**State strata:** TCGA project x FGA tertile (A213) x TP53 status (A215); METABRIC PAM50 x FGA tertile x TP53 status.
Tertiles within project / PAM50 over sequenced tumours with FGA.
**Statistic.** P1 = mean over state strata of coherence(top) − coherence(bottom proliferation tertile), tertiles within
stratum, >= 10 tumours per arm. P0 = the same with strata = project (METABRIC: PAM50) only, on the same tumours (those of
the eligible P1 strata). Proliferation labels of the top and bottom arms permuted within stratum (2,000); one-sided.
- **H216.1 (primary):** TCGA, P1 < 0. **H216.2:** METABRIC, P1 < 0. Holm over the two.
- **Attenuation A = P1 / P0** with a tumour bootstrap within strata (1,000), percentile interval.
**Specificity companion:** 200 pseudo-clock sets (seeds crc32("A216")); clock-specific if P1 is below their 5th
percentile.
**Power (dated note before the tests):** MDE of P1 = 2.49 x null SD; not scored if > 0.25 or if P0 itself is not below its
own pseudo 5th percentile on these tumours (then there is no proliferation loss to explain).
**Reading:**
- H216.1 passes, clock-specific, A >= 0.5: the proliferation loss persists at matched FGA, TP53 and type; these state
  variables do not carry it (proliferation-linked, or an unmeasured state).
- Powered fail with the upper end of A's interval < 0.5: the proliferation loss is carried by these state variables;
  proliferation is a marker.
- Otherwise inconclusive.
Residual proliferation imbalance within strata is reported. No survival test. **Stop rule:** as A215. **Code:**
`paper/exploratory/A215_checkpoint_coherence/` (shared with A215).

**A215/A216 power note (2026-11-06, before any observed coherence difference; `A215_checkpoint_coherence/output/{eligibility,counts,power}.csv`).**
Status: TCGA MC3 (10,295 sequenced tumours; 8,446 primaries in the 23 projects with a TP53 status, 41.4% altered; 8,423 with an RB-pathway status, 23.1% altered);
METABRIC (1,699 of 1,756 expression tumours; TP53 32.9% altered, RB pathway 5.8%). Null means are not 0 because the arms are unequal and proliferation-defined.
H215.1 (TCGA TP53): 18 of 23 projects, 46 cells, 6,449 tumours, null SD 0.0263, MDE 0.066. H215.2 (METABRIC TP53): 12 cells, 1,519 tumours, null SD 0.0561, MDE 0.139.
H215.3 (TCGA RB pathway): 19 projects, 49 cells, 6,586 tumours, null SD 0.0274, MDE 0.068. H215.4 (interaction): 15 projects, 6,272 tumours, null SD 0.0622, MDE 0.155.
S215.5 (matched, 44 cells, 3,652 tumours in pairs): null SD 0.0300, MDE 0.075. H216.1 (TCGA, 23 projects, 94 strata, 5,336 tumours in the arms, 7,955 in P0): null SD 0.0221, MDE 0.055.
H216.2 (METABRIC, 18 strata, 1,032 tumours in the arms, 1,539 in P0): null SD 0.0536, MDE 0.133. All seven have MDE <= 0.25 and are scored on power; for H216.1-2 the registered
precondition (observed P0 below its pseudo 5th percentile) cannot be evaluated before the tests and is checked there. Choices left open (README): MC3 FILTER kept PASS / wga / native_wga_mix;
sequenced = barcode in the MAF (METABRIC: `cases_sequenced`); intact requires a copy-number call; H215.4 status permuted within project x pooled proliferation tertile with
status-specific tertiles recomputed; S215.5 caliper 0.1 SD of the project's score over tumours with TP53 status; P0 on all tertiles of the eligible P1 strata; "residual
proliferation imbalance" in A216 read as top - bottom FGA, TP53 fraction and purity within strata.

## A217 (2026-11-06, user agreement; written before PanCanAtlas subtype labels are downloaded for any project other than BRCA; seen: A210, A213, A215 results and the post hoc TCGA BRCA PAM50 rerun (losses shrink 59-96 % at matched subtype)): do the TCGA genomic-state losses survive matching on molecular subtype? **L1, TCGA, BRCA excluded from the tests (not blind)**
**Why.** The pan-cancer TCGA passes (A210 H210.2 MYC, A213 H213.2 FGA, A215 H215.1 TP53) used project x proliferation
cells. In breast the altered arms are mostly Basal / LumB, and at matched PAM50 the losses mostly vanish. If the same holds
in the other projects, these passes are subtype composition, not lesion effects.
**Data.** PanCanAtlas molecular subtypes: patient attribute `SUBTYPE` of the cBioPortal studies
`<project>_tcga_pan_can_atlas_2018` (COAD and READ: `coadread`), API or study tarball, provenance recorded. Projects of the
A210 set (23) with a subtype for >= 50 % of their primaries enter; tumours without a subtype are dropped. Lesions,
expression, proliferation, coherence exactly as A210 / A213 / A215.
**Cells:** project x subtype x proliferation tertile (tertile within project x subtype); arms as registered originally
(MYC: GISTIC 2 vs 0; FGA: top vs bottom tertile within cell; TP53: altered vs intact); >= 10 per arm. T = mean over cells of
coherence(arm 1) − coherence(arm 0); labels permuted within cell (2,000); one-sided; excess over the null mean reported.
- **H217.1:** TP53. **H217.2:** FGA. **H217.3:** MYC amplification. Holm over the three. BRCA excluded.
**Specificity companion:** 200 pseudo-clock sets per project (seeds crc32("A217")); clock-specific only below their 5th
percentile.
**Descriptive:** the same tumours with project x proliferation cells (the original design) and the ratio of the two
excesses; per-project and per-subtype values; subtype composition of the arms; BRCA reported separately.
**Power (dated note before the tests):** MDE = 2.49 x null SD; > 0.25 not scored.
**Reading:** pass, clock-specific -> the lesion goes with less clock signal beyond subtype; powered fail -> the original
TCGA pass is read as subtype composition and is withdrawn as a lesion effect. **Stop rule:** one amendment.
**Code:** `paper/exploratory/A217_subtype_matched/`.

## A218 (2026-11-06, user agreement; written with A217, before any subtype label outside BRCA is read; seen: A127 H127.1 (SCAN-B, proliferation loss within PAM50 x purity, pass), A166 H166.3 (METABRIC, fail), A97 / A172 per-project proliferation losses, A216 precondition fail): is the tumour proliferation clock loss carried by molecular subtype? **L1, TCGA projects with subtypes, BRCA excluded**
**Why.** Within breast the proliferation loss was tested within PAM50 (SCAN-B pass, METABRIC fail). Elsewhere it was tested
within project only, where proliferation differs by subtype. If subtype (lineage or differentiation state) carries it,
proliferation is a marker of subtype in these projects.
**Design.** Projects and subtypes as A217. Purity proxy as A127 (mean z of COL1A1, COL1A2, COL3A1, DCN, LUM, PDGFRB, PTPRC,
CD3E, CD14, LYZ, CD74), tertile within project. P_proj = mean over project x purity-tertile cells of coherence(top) −
coherence(bottom proliferation tertile, within cell); P_sub = the same over project x subtype x purity-tertile cells;
>= 10 per arm; both on the tumours of the eligible P_sub cells. Proliferation labels of the two arms permuted within cell
(2,000).
- **H218.1:** P_sub < 0 (one-sided), clock-specific (below the pseudo 5th percentile, 200 sets, seeds crc32("A218")).
- **Precondition:** P_proj below its own pseudo 5th percentile on these tumours; otherwise H218.1 is not scored (nothing
  to explain) and both values are descriptive.
- **Attenuation:** A = P_sub / P_proj, delete-one-project jackknife interval (no resampling with replacement).
**Power (dated note before the tests):** MDE of P_sub = 2.49 x null SD; > 0.25 not scored.
**Reading:** H218.1 passes and A >= 0.5 -> the proliferation loss holds within subtype (as SCAN-B); powered fail with the
upper jackknife bound of A < 0.5 -> it is carried by subtype; otherwise inconclusive. Per-project values descriptive.
**Stop rule:** one amendment. **Code:** `paper/exploratory/A217_subtype_matched/` (shared).

**A217/A218 power note (2026-11-06, before any observed coherence difference; `A217_subtype_matched/output/{projects,eligibility,counts,power}.csv`).**
Subtype labels (cBioPortal PanCanAtlas `SUBTYPE`, API 2026-11-06): 22 of the 23 A210 projects have >= 50 % of their primaries labelled (OV 41.9 % does not), but in 10 of them the label equals the
project name (BLCA, KIRC, KIRP, LIHC, LUAD, LUSC, PAAD, PCPG, PRAD, THCA; OV too), which carries no subtype. Projects entering (>= 50 % and >= 2 levels, BRCA excluded): CESC, COAD, ESCA, GBM, HNSC, LGG, READ, SARC,
STAD, TGCT, UCEC (coverage 71 % to 98 %; BRCA 89.8 %, reported separately). H217.1 (TP53): 7 projects, 17 cells, 1,189 tumours, null mean 0.0162, null SD 0.0514, MDE 0.128. H217.2 (FGA): 9 projects, 39 cells, 1,579 tumours,
null mean -0.0002, null SD 0.0329, MDE 0.082. H217.3 (MYC): 4 projects, 9 cells, 387 tumours, null mean -0.0583, null SD 0.0669, MDE 0.166. H218.1 (P_sub): 11 projects, 45 cells, 1,755 tumours in the arms,
null mean -0.0002, null SD 0.0292, MDE 0.073. All four have MDE <= 0.25 and are scored on power; for H218.1 the registered precondition (observed P_proj below its pseudo 5th percentile) is checked in the
tests. Choices left open (README): entry requires >= 2 subtype levels; purity tertile within project over all primaries; P_proj on the tumours of the eligible P_sub cells with proliferation tertile within the cell;
jackknife interval A +/- 1.96 x jackknife SE; original design on projects with an eligible matched cell.

## A219 (2026-11-06, user agreement; written before any per-subtype coherence at matched proliferation is computed; seen: SCAN-B coherence per PAM50 without proliferation matching (A127 exploratory: LumA 0.754, Normal 0.663, Basal 0.538, Her2 0.382, LumB 0.267), the A217/A218 results (lesion and proliferation contrasts within subtype, not subtype contrasts), subtype counts per TCGA project): does tumour clock coherence differ between molecular subtypes at matched proliferation? **L1; TCGA subtyped projects primary (blind), breast directional secondary**
**Why.** A217 found no lesion effect beyond subtype, and the TCGA BRCA arms differed mainly in subtype. If lineage or
differentiation state (subtype) sets tumour clock coherence, coherence should differ between subtypes of one cancer type at
equal proliferation. This makes "subtype, not lesion" a direct test rather than an inference from fails.
**Data.** TCGA: A217 subtypes and projects (CESC, COAD, ESCA, GBM, HNSC, LGG, READ, SARC, STAD, TGCT, UCEC), expression,
proliferation and coherence as A210 / A217. METABRIC PAM50 (`CLAUDIN_SUBTYPE`, A166 five) and TCGA BRCA PAM50 (A217). SCAN-B
PAM50 (A127).
**Cells:** project (METABRIC, SCAN-B: cohort) x proliferation tertile (within project); in a cell a subtype enters with
>= 15 tumours; a cell needs >= 2 subtypes.
- **H219.1 (primary, TCGA non-BRCA):** S = mean over cells of the tumour-weighted SD of the subtype coherences around the
  cell's weighted mean; subtype labels permuted within cell (2,000); one-sided greater (the null mean is > 0 from sampling
  noise; p against the permutation null, excess reported).
- **H219.2 (secondary, breast, directional):** D = mean over cells of coherence(LumA) − coherence(LumB + Her2 + Basal pooled;
  Normal-like left out); one-sided greater; METABRIC and TCGA BRCA pooled (cells from both); labels permuted within cell.
- **H219.3 (secondary, SCAN-B, not blind for the unmatched ordering):** D as H219.2.
Holm over H219.1-3.
**Specificity companion (same permutations):** 200 pseudo-clock sets (A213 / A217 construction; METABRIC and SCAN-B the
stored A166 / A127 sets; seeds crc32("A219")). Subtypes differ in many genes, so a pass is "clock-specific" only if S (D) is
above the 95th percentile of the pseudo sets.
**Power (dated note before the tests):** MDE = 2.49 x null SD (H219.1: on S − null mean); > 0.25 not scored.
**Descriptive:** subtype coherences per cell and project; the same without proliferation matching; residual proliferation
imbalance between subtypes within cells; ER within LumB / Her2 where n allows.
**Reading:** H219.1 passes clock-specifically -> subtype sets clock coherence beyond proliferation in cancer types other than
breast: lineage or differentiation state is a determinant of tumour clock signal; H219.2-3 pass -> luminal A breast
tumours keep more clock signal than other subtypes at equal proliferation. Passes that are not clock-specific -> subtypes
differ in co-expression generally, not in the clock. Powered fails -> subtype does not set coherence at matched
proliferation, and the A217 fails leave proliferation as the only tumour-state correlate. **Stop rule:** one amendment.
**Code:** `paper/exploratory/A219_subtype_coherence/`.

**A219 power note (2026-11-06, before any observed coherence difference; `A219_subtype_coherence/output/{eligibility,counts,power}.csv`).**
Cells project (cohort) x proliferation tertile; a subtype enters with >= 15 tumours, a cell needs >= 2 entering subtypes. H219.1 (TCGA non-BRCA, S): 9 of the 11 A217 projects have an eligible cell (GBM and READ have none), 23 of 33 cells,
2,388 tumours, null mean 0.0764, null SD 0.0100, MDE 0.025 (on S - null mean). H219.2 (METABRIC + TCGA BRCA, D): 2 cohorts, 6 cells (3 + 3), 2,526 tumours, null mean -0.0011, null SD 0.0623, MDE 0.155. H219.3 (SCAN-B, D): 1 cohort, 3 cells,
2,865 tumours, null mean 0.0152, null SD 0.0488, MDE 0.121. All three have MDE <= 0.25 and are scored. Choices left open (README): proliferation tertile within project over all tumours with a subtype label (before the 15-tumour rule);
for D each of LumA, LumB, Her2, Basal enters on its own count (LumA and >= 1 other required, the entering others pooled); weighted SD in population form, unweighted mean over cells; METABRIC pseudo sets are the stored A166 sets and
SCAN-B the stored A127 mean- and variance-matched sets.

## A220 (2026-11-06, user agreement; written before any survival file is downloaded or any signature is fitted; seen: the literature tally (79 claims; LASSO "circadian risk scores" the most common type), Venet 2011, A191-A194 / A198 (single clock scores carry no outcome information beyond proliferation in METABRIC, GSE14520, SCAN-B), A210-A219 TCGA coherence results): are published-style circadian risk signatures better than random gene signatures? **Outcome; TCGA primary, breast secondary**
**Why.** Most clock-outcome claims are LASSO-Cox signatures built from a "circadian rhythm-related" gene list in one TCGA
cancer type and validated by a held-out log-rank test. Random gene signatures are prognostic in many cancers, mostly through
proliferation (Venet 2011). If the clock list does no better than size- and expression-matched random lists, the genre
reports a property of any gene list, not of the clock.
**Data.** TCGA recount3 primary tumours (one per patient, as A210); overall survival from TCGA-CDR (Liu 2018, GDC
PanCanAtlas supplement); stage and age from TCGA-CDR; proliferation metagene of A210; PanCanAtlas subtype where A217 has
one. METABRIC (cBioPortal, OS), SCAN-B GSE96058 (OS). Eligible project: >= 60 OS events among tumours with expression.
**Gene lists.** Primary: MSigDB 2024.1 GOBP_CIRCADIAN_RHYTHM, restricted to genes expressed in the project (median
log2 CPM above the project's 25th percentile of gene medians). Secondary (descriptive): the 14 waveform clock genes.
**Pipeline (the typical published one).** 10 random 70/30 splits per project stratified by event. In training:
univariate Cox p < 0.05 filter (if < 5 genes pass, the 5 lowest p); L1-penalised Cox, penalty by 5-fold CV;
risk score on the held-out 30 %. Outcome per list: held-out Harrell C, mean over splits.
**Random lists.** 200 per project, same size, matched on expression decile (project median log2 CPM), drawn from
expressed genes outside the circadian list; proliferation genes not excluded (that is the point); seeds crc32("A220").
- **H220.1 (primary, unadjusted, TCGA):** P = mean over projects of the clock list's percentile among its random lists;
  null = the same statistic for each random list ranked among the other 199; one-sided greater.
- **H220.2 (primary, adjusted, TCGA):** the same on incremental C: held-out C of Cox(age + stage + proliferation [+
  subtype]) + risk score minus C of the covariate model (both fitted on the training split).
- **H220.3 (secondary, breast):** H220.2 in METABRIC and SCAN-B (covariates age, PAM50, proliferation; stage / size
  where recorded), mean of the two percentiles.
Holm over H220.1-3.
**Power (dated note before scoring):** MDE = 2.49 x null SD of P; a test with MDE > 0.20 (percentile units) not scored.
**Descriptive:** per project unadjusted C, percentile and "validated" (held-out median-split log-rank p < 0.05) for the
clock list and the share of random lists that validate; the 14-gene list; selected genes and their correlation with
proliferation.
**Reading:** H220.1 fails powered and random lists validate as often -> circadian risk signatures are the Venet effect
(any gene list works). H220.1 passes but H220.2 fails -> the clock list beats random lists only through proliferation /
stage. H220.2 (or H220.3) passes -> clock genes carry outcome information beyond random lists and covariates in those
cancers: a real lead, to be replicated. **Stop rule:** one amendment. **Code:** `paper/exploratory/A220_signature_benchmark/`.

## A221 (2026-11-06, user agreement; written before any tumour-vs-normal level is computed and before GTEx recount3 is downloaded; seen: H138.4 (tumour-normal clock pairing, coupling only), A173 (PER1 and NR1D1 sampling-sensitive in liver), A146 / A204b (surgery-time effects), GTEx time-of-death results of the atlas chapter (clock phase follows recorded time of death)): are "clock genes downregulated in tumours" real or a design artefact? **Levels; TCGA paired and TCGA vs GTEx**
**Why.** Many papers report lower PER / CRY / NR1D / DBP / TEF in tumours, often against GTEx (GEPIA-type designs).
Clock-gene levels depend on sampling time: tumours are taken at daytime surgery, GTEx donors die at any hour and some
are ventilated for days; PER1 is also an immediate-early gene. A tumour-vs-adjacent-normal pair shares surgery time and
ischaemia, so it removes these terms.
**Data.** TCGA recount3 (on disk): primary tumours and solid-tissue normals. GTEx recount3 (same pipeline; downloaded for
the matched tissues only): BRCA-breast, LUAD/LUSC-lung, KIRC/KIRP/KICH-kidney cortex, LIHC-liver, THCA-thyroid,
PRAD-prostate, COAD-transverse colon, STAD-stomach, ESCA-oesophagus mucosa, UCEC-uterus. GTEx time of death (v8,
controlled): held in memory, aggregates only. Eligible project: >= 10 tumour-normal pairs (same patient) and >= 50 GTEx
samples of the matched tissue.
**Measures.** log2(CPM + 1), library = protein-coding counts. Per gene: LFC_paired = mean over pairs of tumour - normal;
LFC_gtex = mean over all tumours - mean over all GTEx samples of the tissue (all Hardy classes, as GEPIA).
Clock set: the 14 waveform genes. Random sets: 200 per project, 14 genes matched on expression decile (mean over
tumours and normals), drawn outside the clock set; the same sets for both designs; seeds crc32("A221").
- **H221.1 (primary, paired):** L = mean over projects of [mean clock LFC_paired - mean matched-set LFC_paired]; percentile
  against the 200 sets; one-sided lower (clock-specific downregulation at matched surgery time).
- **H221.2 (primary, design):** Delta = mean over projects of mean clock (LFC_gtex - LFC_paired), same null; two-sided (the
  GTEx design distorts clock genes specifically).
- **H221.3 (secondary, mechanism):** across clock genes x projects, Spearman between the design difference (LFC_gtex -
  LFC_paired) and the GTEx daytime offset (mean in donors with time of death 08-18 h minus mean over all donors of the
  tissue); genes permuted within project (2,000); one-sided greater.
Holm over H221.1-3.
**Power (dated note before scoring):** MDE = 2.49 x null SD on L and Delta (log2 units); not scored if MDE > 0.5.
**Descriptive:** per gene and project both LFCs; share of clock genes "significantly down" (Welch p < 0.05) in each
design vs random genes; immediate-early set (FOS, FOSB, JUN, JUNB, EGR1, EGR2, EGR3, IER2, ATF3, NR4A1); GTEx by Hardy
class; clock level vs ABSOLUTE purity within tumours (composition).
**Reading:** H221.1 passes -> clock genes are lower in tumours than matched genes at the same surgery time (real
level loss; composition still to separate via the purity descriptive). H221.1 fails powered -> the downregulation is
not clock-specific. H221.2 passes -> claims from GTEx designs depend on the comparator; H221.3 passes -> sampling time
explains part of it. **Stop rule:** one amendment. **Code:** `paper/exploratory/A221_tumour_normal_levels/`.

## A222 (2026-11-06, user agreement; written before any clock-gene mutation or copy-number frequency is computed; seen: A215 / A217 used MC3 and GISTIC for TP53, RB1, MYC and FGA only; textbook positions PER1 17p13.1 (next to TP53) and PER3 1p36): are clock genes selected for alteration in tumours, or passengers? **Genomic; TCGA pan-cancer**
**Why.** Papers report frequent clock-gene mutation or loss. Mutation counts scale with coding length, expression and
local mutation rate; copy-number losses ride on arm-level events and on neighbouring tumour suppressors (PER1 next to
TP53). Selection shows as an excess of non-silent over silent mutations, and as losses more focal at the gene than at
its neighbours.
**Data.** MC3 v0.2.8 public MAF (FILTER PASS, as A215); samples with > 1,000 coding mutations excluded (hypermutators);
GISTIC2 thresholded by gene (Xena, as A215); coding length from GENCODE v29 CDS (downloaded); expression = pan-TCGA
median log2 CPM (recount3). Clock set: the 14 waveform genes.
- **H222.1 (primary, mutation):** R = log[(non-silent + 0.5) / (silent + 0.5)] summed over the clock set; null 1,000
  random 14-gene sets matched on CDS-length decile and expression decile; one-sided greater.
- **H222.2 (primary, deletion focality):** per gene and project, deep-deletion (-2) frequency minus the median over its
  25 flanking genes on each side (same chromosome, gene order of the GENCODE coordinates); mean over clock genes and
  projects; null 1,000 random 14-gene sets; one-sided greater.
- **H222.3 (secondary, amplification focality):** H222.2 on +2; two-sided.
Holm over H222.1-3.
**Power (dated note before scoring):** MDE = 2.49 x null SD of each statistic; reported in the null's units and as a
percentile shift.
**Descriptive:** per clock gene: mutation rate per Mb CDS vs matched genes, truncating share, shallow and deep loss
frequency vs neighbours per project; PER1 loss conditional on TP53 loss (17p); PER3 vs the 1p36 block.
**Reading:** both primaries fail powered -> clock genes are mutated and lost as expected for their length, expression and
position: passengers, and "frequent alteration" claims report background rates. A pass names a clock gene set under
selection (which genes drive it is descriptive). **Stop rule:** one amendment. **Code:**
`paper/exploratory/A222_clock_alterations/`.

**A222 power note (2026-11-06, before any clock-gene statistic; `A222_clock_alterations/output/power.csv`, null draws `null_draws.csv`).**
MC3: 10,178 samples, 448 hypermutators excluded, 9,730 kept; 19,849 protein-coding genes with expression; GISTIC 26 projects
(>= 100 samples), 9,681 samples; all 14 clock genes have a full 25-gene flank. H222.1 (R): null mean 1.136, SD 0.145, MDE
0.360 (a 1.43-fold non-silent / silent excess). H222.2 (deep-deletion focal excess): null mean 0.00027, SD 0.00043, MDE
0.00107 (0.11 percentage points above the flank). H222.3 (amplification): null mean 0.00034, SD 0.00039, MDE 0.00097. All
scored. Choices left open: README.

**A221 power note (2026-11-06, before any clock-gene level is computed; `A221_tumour_normal_levels/output/{eligibility,power,null_draws,random_sets}.csv`).**
13 projects eligible (BRCA 112 pairs, COAD 41, ESCA 13, KICH 25, KIRC 72, KIRP 32, LIHC 50, LUAD 58, LUSC 51, PRAD 52,
STAD 34, THCA 59, UCEC 22); GTEx recount3 samples per tissue 94 (kidney cortex, shared by KICH / KIRC / KIRP) to 706
(thyroid); lung shared by LUAD / LUSC, so projects are not independent in H221.2-3. Random 14-gene sets (200): mean paired
LFC +0.091, mean design difference −0.023. H221.1 (L): null SD 0.065, MDE 0.163 log2. H221.2 (Delta): null SD 0.049, MDE
0.121 log2. H221.3: 182 gene x project points. All scored. Choices left open (README): CPM on recount3 coverage over
protein-coding genes; one primary per patient (lowest barcode); random-set candidates mean CPM >= 1; daytime window
[8, 18) h; Hardy class from the public recount3 GTEx annotation.

## A223 (2026-11-06, user agreement; written after the A221 TCGA per-gene mean paired LFCs were seen (PER1-3, CRY2, TEF, RORA down; ARNTL, NPAS2, NR1D1, DBP, CIART flat or up), before any CPTAC clock-gene level is computed; seen: A186 CPTAC pairs (coupling only, no levels), A27 (surgical glucocorticoids induce PER1)): is the tumour-vs-paired-normal clock-gene level pattern a damped clock read at daytime surgery? **Levels; CPTAC primary (blind), TCGA secondary (not blind)**
**Why.** If a tumour's clock is damped (amplitude a < 1 relative to its normal tissue) and both are sampled at the same
daytime hour t, the expected log level difference of clock gene g is -(1 - a) lam amp_g cos(w t - psi_g): genes near their
peak at surgery fall, genes near their trough rise, and the mean over genes spread around the cycle stays near 0. This
reads damping from levels alone, independently of coherence. Competing reading: stress or glucocorticoid response of the
adjacent tissue (PER1, TEF), which predicts losses of response genes, not trough-gene gains.
**Data.** CPTAC-3 open pairs (A186 cache, 8 types, 713 pairs; types with >= 20 pairs); TCGA A221 pairs (13 projects).
log2(CPM + 1), protein-coding library (CPTAC: all ENSG genes as cached, documented). Per pair and gene: tumour - normal.
**Predictor.** x_g = -amp_g C_g / ln 2 with the GTEx waveform (`clocksync.gtex()`, amp_g and psi_g fitted on time of
death), C_g = mean of cos(w t - psi_g) over t uniform in [8, 18) h (surgery window, as A221 / A204b). Under the model the
slope of LFC_g on x_g is (1 - a) lam.
**Statistic.** Per type, OLS slope beta of the mean paired LFC_g (14 genes) on x_g with an intercept (a uniform shift is not
damping); B = mean of beta over types. Null: phase rotation, psi_g + delta for all genes, delta on a 1-h grid uniform over
the 23 non-zero offsets plus 977 continuous draws (1,000 total), recomputing x and B; one-sided greater. Seeds
crc32("A223").
- **H223.1 (primary, CPTAC):** B > rotation null.
- **H223.2 (secondary, CPTAC without PER1):** the same on 13 genes (PER1 is a glucocorticoid / surgery gene, A27).
- **H223.3 (secondary, TCGA, not blind):** H223.1 on the 13 A221 projects.
Holm over H223.1-3.
**Power (dated note before scoring):** under the null, B's rotation SD; MDE = 2.49 x SD; not scored if MDE > 0.5 (B = 1 is
complete damping at lam = 1).
**Descriptive:** per type beta and R^2; the best-fitting surgery-window centre over a 0-23 h grid (should be daytime if
the model holds); immediate-early gene LFC per type next to beta; implied (1 - a) with the A-chapter lam of the matched
GTEx tissue where available.
**Reading:** H223.1 passes -> the paired level pattern is phase-structured as a damped tumour clock read at daytime
surgery, a level-based readout of tumour damping replicated in an independent cohort; with H223.2 it does not rest on PER1.
Fails powered -> the paired loss is not organised by circadian phase (stress / composition remain). A pass whose best
window centre is at night argues against the surgery-time model. **Stop rule:** one amendment. **Code:**
`paper/exploratory/A223_damping_levels/`.

**A223 power note (2026-11-06; `A223_damping_levels/output/power.csv`): not scored.** Rotation null (1,000 offsets): H223.1
CPTAC (8 types, 14 genes) SD 1.284, MDE 3.20; H223.2 (13 genes) SD 1.261, MDE 3.14; H223.3 TCGA (13 projects) SD 1.127,
MDE 2.81. All above the 0.5 limit (complete damping at lam = 1 gives B = 1). No B at the observed phases was computed.

**A220 power note (2026-11-06, before any clock-list C-index; `A220_signature_benchmark/output/{eligibility,power,null_percentiles}.csv`).**
18 TCGA projects with >= 60 OS events (BLCA, BRCA, CESC, COAD, ESCA, GBM, HNSC, KIRC, LGG, LIHC, LUAD, LUSC, MESO, OV, PAAD,
SARC, STAD, UCEC; 71-262 events); METABRIC 1,756 (1,044 events), SCAN-B 3,073 (289). GOBP_CIRCADIAN_RHYTHM 207 genes, 137-154
expressed per cohort. 200 random lists per cohort, 10 splits. H220.1: null SD 0.071, MDE 0.177 (percentile units); H220.2:
SD 0.071, MDE 0.177; both scored. H220.3 (2 cohorts): SD 0.215, MDE 0.535 > 0.20: not scored (per-cohort percentiles
descriptive). Random lists alone: the share validating (held-out median-split log-rank p < 0.05) is 0-100 % by cohort (100 %
in KIRC, LGG, METABRIC; 99.5 % SCAN-B). Choices left open (README): in-sample training risk score in the adjusted Cox model
(same for clock and random lists); covariates observed in < 50 % dropped; lambda.min on CV partial likelihood; no gene
selected -> C 0.5.

## A224 (2026-11-06, user agreement; written before any clock-gene expression is related to its copy number; seen: A222 (PER1 losses co-occur with TP53 loss in 98.9 %, PER3 with 1p36), A215 H215.1 (TCGA TP53 loss −0.085, clock-specific), A217 (lost at matched subtype), A219 H219.2 (LumA − others +0.276, METABRIC + TCGA BRCA)): does cis copy number of the clock genes carry the TP53 and subtype coherence contrasts, and is the "p53-clock" link positional? **L1 + levels; TCGA, METABRIC**
**Why.** PER1 sits next to TP53 (17p13.1), PER3 in 1p36, NR1D1 in the ERBB2 amplicon (17q12). Copy-number variation of a
clock gene between tumours adds non-clock variance to its level and lowers between-tumour coherence most where the genome is
most altered (TP53-altered, non-luminal-A tumours). A TP53-mutant vs wild-type level contrast also inherits PER1 dosage
through 17p loss.
**Time-dependence (user caution).** Clock-gene levels depend on sampling time, so dosage is estimated with phase
covariates: per tumour, the in-plane phase from the other 13 clock genes (leave-the-gene-out projection on the GTEx
waveform plane, cohort-centred levels), and log2 level of gene g ~ project + beta_g CN_g + cos + sin of that phase. Dosage
enters as an additive log offset; a CN_g x (cos, sin) term is fitted descriptively (dosage on amplitude rather than mean).
CN_g = GISTIC thresholded state (-2..2, numeric; TCGA Xena by gene, METABRIC data_cna.txt).
**Adjusted expression.** y_sg - beta_g CN_sg, with beta_g pooled within cohort (TCGA: over the A219 / A215 projects;
METABRIC alone). Pseudo-clock genes are adjusted the same way with their own CN and their own beta (same phase covariates);
genes without a CN call keep their level.
- **H224.1 (primary, METABRIC + TCGA BRCA):** A219 H219.2 recomputed on adjusted expression (same cells, permutations, 200
  pseudo sets adjusted alike). Outcome D_adj; pass = permutation p < 0.05 and above the pseudo q95. Attenuation (D - D_adj)
  / D with a leave-one-cell-out jackknife, reported.
- **H224.2 (secondary, TCGA):** A215 H215.1 TP53 contrast (project x proliferation cells) on adjusted expression, same
  design and companion; outcome T_adj and attenuation.
- **H224.3 (secondary, TCGA levels):** among tumours with neither PER1 / TP53 loss (GISTIC >= 0 at both) nor PER2 change,
  TP53 non-silent mutant vs wild-type difference in phase-adjusted log2 PER1 and PER2 (mean of the two, project fixed
  effects; labels permuted within project, 2,000); two-sided (p53-biology prediction: a mutation effect without 17p loss).
  Descriptive: the naive mutant vs wild-type contrast over all tumours, and the 2 x 2 (mutation x 17p loss) means.
Holm over H224.1-3.
**Power (dated note before scoring):** MDE = 2.49 x null SD; H224.1-2 not scored if MDE > 0.25 (coherence); H224.3 not
scored if MDE > 0.3 log2. Also before scoring: the fitted beta_g (dosage sanity; expected > 0) and the share of each clock
gene's between-tumour variance carried by CN.
**Reading:** H224.1 passes with small attenuation -> the luminal A coherence advantage is not cis copy-number noise; large
attenuation or a fail -> part of it is genome instability on the clock loci, and the A219 reading weakens accordingly.
H224.2 near 0 -> the TCGA TP53 loss was positional or compositional (with A217). H224.3 fails powered -> no TP53-mutation
effect on PER1 / PER2 beyond 17p dosage: "p53 regulates the clock" claims from tumour level contrasts are positional.
**Stop rule:** one amendment. **Code:** `paper/exploratory/A224_cis_copy_number/`.

**A224 power note (2026-11-06, before any adjusted contrast; `A224_cis_copy_number/output/{power,dosage,counts,eligibility,eligibility_h224_3}.csv`).**
Dosage (phase-adjusted, log2 per GISTIC step, all p < 0.001): TCGA pooled (23 projects) beta 0.236 (CRY1) to 0.462 (PER3);
PER1 0.388, PER2 0.362, NR1D1 0.274; CN share of between-tumour variance 0.016 (NPAS2) to 0.153 (CLOCK), PER1 0.048.
METABRIC (array) beta 0.017 (NR1D1) to 0.411 (PER2), PER1 0.033; CIART has no CN row (not adjusted). CN x phase joint p
< 0.01: TCGA ARNTL, NPAS2, NR1D2, CIART; METABRIC PER3, DBP (descriptive). Power: H224.1 (6 cells, 2,526 tumours) null SD
0.053, MDE 0.133; H224.2 (46 cells, 6,449) SD 0.027, MDE 0.066; H224.3 (14 projects, 2,362 tumours) SD 0.032, MDE 0.079
log2. All scored. Choices (README): seeds of A219 / A215 reused, crc32("A224") for H224.3; TCGA beta pooled over the A215
projects; pseudo genes use the 14-gene phase; H224.1 one-sided greater, H224.2 one-sided less, H224.3 two-sided.

## A225 (2026-11-06, user agreement; written after A224 and its post hoc breakdown were seen (LumA − others +0.276 → +0.101 after cis-CN adjustment; ARNTL and PER3 carry most of it; CN variance 2-4 x higher outside LumA); before any simulation, masked estimate or trans contrast is computed): copy number at the clock loci — readout artefact or clock biology? **Validation (simulation) gating L1 tests; METABRIC + TCGA BRCA, TCGA pan-project**
**Why.** A224 subtracts an estimated dosage offset; this removes a readout artefact (a gene's level shifted at every hour
inflates its between-tumour spread and lowers coherence with an intact oscillator) but could also remove real effects of
copy number on the oscillator. Three pieces separate them: (A) a semi-synthetic test of the estimators on the real copy-number
profiles, which gates (B) a model-free masked estimate, and (C) a trans test of copy-number effects on the rest of the
clock.
**A. Simulation (gate).** Cells, tumours, subtypes and real GISTIC profiles of A224 H224.1 (METABRIC, TCGA BRCA). Per tumour:
internal phase uniform, log level of clock gene g = m_g + lam amp_g cos(phi - psi_g) + beta_g CN_g + (CN x phase terms of
A224 `dosage.csv`) + Gaussian noise with the gene's A224 residual SD; lam calibrated so that the raw estimator reproduces the
observed LumA coherence of each cohort (LumA only, the CN-quiet arm). Pseudo genes: same construction without the clock
term, with their own CN and betas. 200 replicates per scenario, seeds crc32("A225").
- Scenario N (null): identical clocks in all subtypes. Scenario D: non-LumA lam x 0.7 (30 % damping).
- Estimators: raw (A219), adjusted (A224, betas re-estimated in each replicate), masked (B).
- **V225.1 (positive control, descriptive):** mean raw D under N, as a share of the observed +0.276.
- **V225.2 (gate):** an estimator passes if under N |mean D| <= 0.03 and the A219 permutation test rejects in <= 7 % of
  replicates, and under D its mean recovers the planted difference within ±25 % (true D under D from the noise-free
  coherences of the simulation). Only passing estimators are used in B; if none passes, B is not scored.
**B. Masked estimate.** A gene's value is set missing in tumours where its GISTIC state is not 0; coherence from
pairwise-complete gene-gene covariances (tumour weights as A219); pseudo genes masked by their own CN.
- **H225.1 (primary):** A219 H219.2 (LumA − others, METABRIC + TCGA BRCA, same cells and permutations) with the
  gate-passing estimator (masked if it passes, else adjusted); pass = permutation p < 0.05 and above the q95 of the 200 pseudo
  sets treated alike.
**C. Trans test (TCGA, 23 A215 projects).** For clock gene g, arms: GISTIC(g) <= -1 vs = 0; cells = project x FGA tertile x
proliferation tertile (both within project); arm >= 10 per cell. Outcome: coherence of the other 13 clock genes (raw; each
of them masked where its own CN is not 0, so they carry no dosage artefact). D_g = mean over cells of coh(loss) − coh(neutral);
labels permuted within cell (2,000).
- **H225.2 (primary, g = ARNTL):** one-sided less; companion: 200 expression-matched 13-gene pseudo sets with the same
  ARNTL-defined arms (masked alike): pass = p < 0.05 and below the pseudo q05 (a loss specific to the clock, not to the
  genomic state of ARNTL-lost tumours).
- **H225.3 (secondary):** the same for each of the other 13 clock genes with >= 5 eligible cells; Holm within H225.3.
Holm over H225.1-2.
**Power (dated note before B and C are scored, after A):** MDE = 2.49 x null SD; H225.1 not scored if MDE > 0.25,
H225.2 not scored if MDE > 0.15.
**Reading:** V225.1 near 1 -> copy-number noise alone can produce the luminal A gap. H225.1 fails with a passing estimator
-> the luminal A advantage is a readout artefact (A219 withdrawn as a biological claim); passes -> a clock-specific
advantage survives without dosage modelling. H225.2 passes -> ARNTL copy loss disorganises the rest of the tumour clock (a
biological copy-number effect, and the correct estimator is the one validated in A); fails powered -> ARNTL loss acts on its
own level only. **Stop rule:** one amendment. **Code:** `paper/exploratory/A225_cn_artefact_or_biology/`.

**A225 part A result and power note (2026-11-06, before B / C are scored; `A225_cn_artefact_or_biology/output/{sim_calibration,sim_summary,sim_replicates,power,eligibility}.csv`).**
Calibration (LumA raw coherence): METABRIC lam 0.318 (0.646 vs target 0.647), TCGA BRCA lam 0.542 (0.486 vs 0.489). True D
under scenario D defined as the mean D of 200 copy-number-free simulations at the real cell sizes (the population coherence of
a rank-2 clock is about 1 at any amplitude, so "noise-free" was read this way): pooled 0.158.
- **V225.1:** mean raw D under N −0.003 (SD 0.086) = −0.012 of the observed +0.276: additive cis dosage noise at the A224
  betas, CN x phase terms and the real joint CN profiles does not produce a luminal A gap.
- **V225.2 gate:** raw passes (N: −0.003, rejection 3 %; D: recovery 0.97); adjusted passes (−0.007, 3 %; 0.97); masked
  fails (rejection 13 % under N; estimable only in METABRIC: TCGA BRCA arms lack pairwise-complete tumours). In the
  simulation adjusted ≈ raw.
- **H225.1** uses the adjusted estimator (masked failed): identical to A224 H224.1 (already scored, +0.101, not
  clock-specific); null SD 0.053, MDE 0.133.
- **H225.2 (ARNTL) not scored:** 59 cells with both arms >= 10, none estimable after masking the other 13 genes (19-58 % of
  ARNTL-loss tumours CN-neutral at each other gene). **H225.3:** PER1 (9 cells, 662 tumours, MDE 0.182), TEF (7, 592, 0.193),
  NR1D2 (5, 303, 0.259) eligible (>= 5 cells); the others not.

## A226 (2026-11-06, user agreement; confirmatory; written before any new cohort is identified, downloaded or read; seen: A151 H151.1 (pairing 0.860, 98 patients / 17 cohorts), A149 note 5 / A151 note 2 POST HOC (bounded statistic: fibroblasts 0.699 vs pseudo q95 0.542, endothelium 0.701 vs 0.527, p 0.005)): does tumour stroma follow the same patient's normal-tissue stroma clock-specifically, in cohorts not used before? **L3-like pairing; single-cell, independent cohorts**
**Why.** The stromal pairing is one of the two positive cancer claims, but its clock specificity was established post hoc
(bounded statistic chosen after the registered one broke down). This registers that statistic and tests it on new data.
**Cohorts.** Single-cell or single-nucleus studies with, from the same patients, tumour and adjacent-normal (or
paired normal) tissue and annotated or annotatable fibroblasts / endothelial cells, **not among the A149 / A151 cohorts**
(`A151_fibroblast_pairing/output/units.csv`, `A149_paired_tumour_normal/output/units.csv`). Sources: on-disk studies not
yet used (scRNAseq_database, CELLxGENE Census incl. Moorman 2025 and LuCA donors outside A149), public GEO / CELLxGENE
studies found by a dated search (processed matrices only; <= 50 GB new data). Cell typing as A149 (marker rules of its
note 1, or the authors' labels when present, documented per study). The cohort list is fixed in the power note before
any statistic is computed.
**Statistic (A149 note 5, now registered).** Per patient, clock state of each compartment from the pseudobulk on the GTEx
waveform plane, centred within cohort; per cohort the uncorrected |complex correlation| |sum za conj(zb)| / sqrt(sum |za|^2
sum |zb|^2) between tumour and normal compartment states; pooled with patient-count weights. Eligibility per compartment
and cohort as A151 (>= 4 patients with both; pseudobulk depth rules of A149).
- **H226.1 (primary, fibroblasts):** pooled statistic > within-cohort permutation null (patients' normal states shuffled,
  10,000), one-sided.
- **H226.2 (primary, fibroblast specificity):** pooled clock statistic > the pooled statistic of 200 expression-matched
  pseudo-clock sets (A149 note 2 (iv) construction, same patients and weights); p = (1 + #{pseudo >= clock}) / 201.
- **H226.3 / H226.4 (secondary):** the same for endothelial cells.
Holm over H226.1-2 and over H226.3-4 separately.
**Power (dated note before scoring):** patients and cohorts per compartment; permutation null SD and MDE (2.49 x SD).
Not scored if < 30 patients in a compartment or if MDE > the A151-note-2 clock−null difference halved (i.e. the test
must detect half the post hoc effect).
**Reading:** H226.1-2 pass -> the stromal host pairing is confirmed clock-specific in independent data (main-text claim).
H226.1 pass, H226.2 fail -> pairing confirmed, clock specificity not replicated (the claim stays "patient-specific stromal
pairing", specificity post hoc only). Not scored / fails -> the stromal claim stays post hoc and is reported as such.
**Stop rule:** one amendment. **Code:** `paper/exploratory/A226_stroma_pairing_confirm/`.

**A226 power note (2026-11-06; `A226_stroma_pairing_confirm/output/{eligibility,power}.csv`, scout `scout/SCOUT.md`): not scored.**
Dated scout (user database, on-disk inventories, CELLxGENE Census 2025-11-08, GEO E-utilities 6,226 single-cell series): every
study with same-patient tumour and normal samples and enough pairs is already an A149 cohort. New data (9.9 GB, `$SCCC_DATA/a226`)
give, under the A149 rules, fibroblasts 11 patients / 2 cohorts (GSE308103 lung 5, Qian 2020 CRC 6) and endothelium 13 / 2
(GSE308103 5, GSE183219 8); without depth and presence rules 48 / 6 and 47 / 5. Null SD (simulated independent states at the
real cohort sizes, calibrated on the A151 pseudo-set spread) 0.154 / 0.162, MDE 0.383 / 0.403 (sensitivity set 0.184 / 0.212)
against the threshold (half the A151 clock − pseudo-median difference) 0.120 / 0.134. Both compartments below 30 patients:
not scored. Design note: the threshold was too strict for any feasible cohort set (the same method gives MDE 0.127 at
A151's own 98 patients); the claim nevertheless stays post hoc, as registered.

## A227 (2026-11-06, user agreement; written before any tumour-vs-normal pairing of immune or mural compartments is computed; seen: A149 / A151 / A226 for epithelium, fibroblasts and endothelium (pairing and post hoc clock specificity), A152 (resident vs circulating lymphocyte coupling in normal tissues), the blood-composition result (blood carries diurnal composition, no coherent core clock), the compartment label counts of the A149 pseudobulk files): which tumour-microenvironment compartments keep the same patient's normal-tissue clock state? **L3-like pairing; single-cell, all paired patients**
**Why.** H138.4 shows bulk tumours share the patient's normal-tissue clock state; stroma does so post hoc (A149 note 5). The
immune and perivascular compartments have never been tested. Resident compartments should share the host clock more than
recently recruited ones.
**Data.** The A149 cohort set (27 cohorts) plus the A226 cohorts, same cell labels (`cancer_sc` labels / A149 marker rules)
and the same eligibility rules (>= 4 patients with both tumour and normal pseudobulks of the compartment, >= 50 cells,
sample median UMI >= 1,500, compartment present in >= 60 % of patients or >= 6). Compartments: myeloid, T_NK, B_plasma,
mural (pericytes / smooth muscle), mast. Fibroblasts and endothelium recomputed on the same set as descriptive anchors only
(not confirmatory: mostly the same patients).
**Statistic (A226).** Uncorrected |complex correlation| of tumour vs normal compartment states (GTEx waveform plane, centred
within cohort), pooled with patient-count weights; within-cohort permutation null (10,000); 200 expression-matched
pseudo-clock sets (A149 note 2 (iv) construction).
- **H227.1 / H227.2 (primary: myeloid / T_NK):** pairing > permutation null (one-sided) **and** clock statistic > pseudo q95
  (p = (1 + #{pseudo >= clock}) / 201); a pass needs both.
- **H227.3-5 (secondary: B_plasma, mural, mast):** the same.
Holm over the pairing p of the scored primaries, and separately over the scored secondaries; clock specificity read per
compartment.
**Power (dated note before scoring):** patients and cohorts per compartment; MDE = 2.49 x permutation-null SD (null from
shuffled patient labels only). Not scored if < 40 patients or MDE > 0.25.
**Descriptive:** per compartment clock statistic minus pseudo median, ordered; the same in the cancer types with >= 12
patients; fibroblast / endothelial anchors.
**Reading:** a clock-specific pass -> that compartment keeps the patient's clock state in the tumour. Immune passes carry
the caveat that tumour and normal immune cells share the same circulation (diurnal trafficking and composition), so they
mark host time reaching the tumour, not necessarily a local clock. Pairing without clock specificity -> patient-specific
non-clock structure (as A149 epithelium). **Stop rule:** one amendment. **Code:** `paper/exploratory/A227_tme_compartment_pairing/`.

**A227 power note (2026-11-06, before any observed pairing; `A227_tme_compartment_pairing/output/{eligibility,eligibility_by_cancer,power}.csv`).**
34 cohorts (27 A149 with estimable S1 + 7 A226). Permutation null from shuffled pairings only (identity never evaluated),
10,000 per cohort. Myeloid 169 patients / 24 cohorts, null 0.361 SD 0.036, MDE 0.089; T_NK 227 / 27, 0.320, 0.032, 0.079;
B_plasma 172 / 22, 0.339, 0.035, 0.087; mast 50 / 8, 0.390, 0.066, 0.165; mural 37 / 7 (MDE 0.189): **not scored** (< 40
patients). Anchors: fibroblast 93 / 16 (MDE 0.119), endothelial 101 / 14 (0.114). Scored: H227.1-3, H227.5. Patients are not
de-duplicated across A149 and A226 cohorts (Qian 2020 donor overlap with A149 unconfirmed; README).

## A228 (2026-11-07, user agreement; written after A227 was scored (T_NK 0.539 vs pseudo q95 0.510; myeloid 0.581 vs 0.548) and before any cross-compartment or T/NK-subset pairing is computed or any subset label is called; seen: A150 (healthy lymphocyte coupling 0.094 vs pseudo 0.087), A152, A159, A140 (TNBC immune gate fails), fibrotic-lung T-cell loss of host loading, the blood-composition result): does the tumour immune clock state carry body time, and is the T/NK pairing a composition effect? **L3-like pairing; single-cell, A227 cohorts**
**Why.** A227 pairing is positive if a patient's immune clock state is shared between tumour and normal tissue. Three sources
give that: (i) the cells follow body time (surgery hour, host phase); (ii) a patient-specific immune clock state not tied to
body time (recirculating clones, a T-cell trait); (iii) diurnal or patient-specific subset composition (CD4 / CD8 / NK) whose
clock-gene profiles differ. (i) predicts that the tumour immune state also pairs with the same patient's normal-tissue stroma,
which follows body time (A150, A151 post hoc); (ii) and (iii) do not. (iii) predicts that pairing vanishes within subsets.
Tumour immune cells and normal stroma come from different samples, so ambient RNA cannot create cross-compartment pairing.
**Data.** The 34 A227 cohorts, A227 tables, labels and eligibility (>= 4 patients with both pseudobulks, >= 50 cells,
sample median UMI >= 1,500). Statistic, permutation null (10,000, normal-side rows shuffled within cohort), 200
expression-matched pseudo-clock sets and pooling exactly as A227.
**Part 1, cross-compartment (does the immune state carry body time?).** Reference = the patient's adjacent-normal
fibroblast pseudobulk (the strongest host compartment in healthy tissue, A150).
- **H228.1 (primary):** tumour T_NK vs adjacent-normal fibroblast: pairing > permutation null **and** clock statistic >
  pseudo q95.
- **H228.2 (primary):** tumour myeloid vs adjacent-normal fibroblast: the same.
- **Descriptive:** the same against adjacent-normal endothelium (replication); tumour endothelium vs adjacent-normal
  fibroblast (positive control: two compartments that both follow body time); adjacent-normal T_NK / myeloid vs
  adjacent-normal fibroblast (same sample, ambient not excluded; reported only).
**Part 2, within T/NK subsets (is the T/NK pairing composition?).** Per-cell marker rule inside the A227 T_NK compartment,
raw counts, fixed now: CD3 = CD3D + CD3E + CD3G, CD8 = CD8A + CD8B; **CD8 T**: CD3 >= 1 and CD8 >= 1 and CD4 = 0; **CD8-negative
T** (mostly CD4, plus gamma-delta / MAIT): CD3 >= 1 and CD8 = 0; **NK**: CD3 = 0 and KLRF1 + NCAM1 + FCGR3A + GNLY + NKG7 >= 2;
other cells dropped. Subset pseudobulks per sample (>= 50 cells), tumour vs adjacent-normal of the same subset.
- **H228.3:** CD8 T tumour vs normal; **H228.4:** CD8-negative T; **H228.5:** NK. Each: pairing > null **and** clock > pseudo
  q95 (pseudo sets rebuilt on the subset pseudobulks).
- **Descriptive:** per patient and sample, subset fractions (CD8 / CD8-negative / NK / dropped) in tumour and normal.
**Holm** over the scored hypotheses of Part 1 and, separately, of Part 2; clock specificity read per hypothesis.
**Power (dated note before scoring):** patients and cohorts per test; permutation-null SD from shuffled pairings only; MDE
= 2.49 x SD. Not scored if < 40 patients or MDE > 0.25.
**Reading.**
- H228.1 pass: tumour T/NK cells carry body time (source (i)); the circulation is then a route of entrainment, not an
  alternative to a clock. H228.1 not passing while A227 H227.2 passed: not evidence for (ii); reported as "not detected
  at MDE x" (the cross statistic is expected to be smaller than the within-type one even under (i), because T/NK coupling
  to the host is low).
- H228.3 and H228.4 both pass: the T/NK pairing is not subset composition (composition inside a subset, e.g. naive vs
  effector, remains possible). Neither passing with adequate power: composition not excluded; A227 H227.2 is then read
  as "T/NK compartment state", not "T-cell clock".
- Same rules for myeloid in Part 1. No myeloid subsets (stop rule).
**Stop rule:** one amendment; no further marker rules or references after the power note. **Code:**
`paper/exploratory/A228_immune_body_time/`.

**A228 power note (2026-11-07, before any observed pairing; `A228_immune_body_time/output/{eligibility,power,subset_status,subset_fraction_summary}.csv`, commit a4a78d7).**
34 A227 cohorts; permutation null from shuffled pairings only (10,000 per cohort). Part 1: H228.1 tumour T_NK vs normal
fibroblast 118 patients / 17 cohorts, null 0.365 SD 0.043, MDE 0.108; H228.2 tumour myeloid vs normal fibroblast 116 / 17,
0.370, 0.043, 0.107; descriptive vs normal endothelium T_NK 129 / 18 (MDE 0.104), myeloid 124 / 17 (0.102); positive control
tumour endothelium vs normal fibroblast 76 / 13 (0.131); same-sample normal T_NK 135 / 18 (0.098), normal myeloid 101 / 16
(0.117). Part 2: H228.3 CD8 T 191 / 23, null 0.321 SD 0.035, MDE 0.088; H228.4 CD8-negative T 214 / 26, 0.322, 0.034, 0.084;
H228.5 NK 67 / 9, 0.341, 0.062, 0.155. All five hypotheses scored. GSE308103 lacks CD8B (CD8 = CD8A there). Median T_NK
fractions per sample, tumour / normal: CD8 T 0.275 / 0.309, CD8-negative T 0.594 / 0.504, NK 0.044 / 0.051, dropped 0.029 /
0.029. The A149 compartment-presence rule is not applied to subset pseudobulks (>= 50 cells and median UMI rules only).

## A229 (2026-11-07, user agreement; written after A228 was scored and after reading summaries of Wang et al., Cell 2024 (mouse: endothelial clocks time ICAM-1 and gate CD8+ T cell infiltration), before any adhesion gene or T-cell share is related to any clock state; motivated in part by the A228 descriptive endothelial result (tumour T/NK and myeloid vs normal endothelium above pseudo q95), so the choice of endothelium is not blind): do tumour endothelial adhesion genes and the tumour CD8 T share follow the tumour endothelial clock state in human tumours? **L1-like association; single-cell, A227 cohorts**
**Why.** Wang et al. show in mice that endothelial clocks time adhesion molecules and thereby CD8+ T cell entry into
tumours. Without sampling times, the human counterpart is: across patients, adhesion-gene expression of tumour endothelium
and the CD8 share of tumour T cells vary with the tumour endothelial clock state, beyond what matched non-clock gene sets
give. Spread of clock states comes from surgery hour and patient phase (mostly daytime surgery), which limits power.
**Data.** The 34 A227 cohorts and tables; tumour endothelial pseudobulk per patient (>= 50 cells, A227 eligibility);
T/NK subset counts from A228 (`subset_fractions.csv`, tumour samples, summed per patient, patients with >= 50 T/NK cells).
**State.** z = complex clock state of the tumour endothelial pseudobulk (`clocksync.paired.projector` / `states`, GTEx
waveform, counts=True, fitted on the tumour endothelial rows of the cohort), centred within cohort x chemistry.
**Statistic.** Readout y (real), centred within cohort x chemistry: r = |sum y conj(z)| / sqrt(sum y^2 sum |z|^2) (the A227
statistic with a real second argument), pooled over cohorts with patient-count weights. Null: y permuted across patients
within cohort (10,000). Clock specificity: z recomputed from 200 expression-matched pseudo-clock sets on the same tumour
endothelial genome pseudobulks (A227 construction), same y; read as clock-specific only if r > pseudo q95.
- **H229.1 (primary):** y = log2(CPM + 1) of ICAM1 in the tumour endothelial pseudobulk.
- **H229.2 (primary):** y = logit of the CD8 T share of tumour T cells, CD8T / (CD8T + CD8neg_T) (A228 rule; +0.5 to
  both counts).
- **H229.3 / H229.4 (secondary):** VCAM1, SELE (as H229.1).
- **H229.5 (secondary):** y = log of the T/NK share of all cells of the tumour sample (sum of compartment cell counts).
Each: pass = Holm p < 0.05 (primaries and secondaries separately) **and** r > pseudo q95.
**Descriptive.** Percentile of ICAM1 among 200 expression-matched endothelial genes used as y (are adhesion genes more
tied to the endothelial clock state than typical endothelial genes?); the same tests in adjacent-normal endothelium; per
cancer type with >= 12 patients; arg(sum y conj(z)) mapped to the GTEx clock-time axis (implied peak, no prediction).
**Power (dated note before scoring):** patients and cohorts per test; permutation-null SD (shuffled y only, observed never
computed); MDE = 2.49 x SD. Not scored if < 40 patients or MDE > 0.25.
**Reading.** H229.1 pass: human tumour endothelial ICAM1 varies with the endothelial clock state (Wang's gate, in
humans, untimed). H229.2 pass: the CD8 share of tumour T cells follows the endothelial clock state. Both: the mouse
mechanism has a human correlate. Not passing: "not detected at MDE x", not a refutation (daytime surgery compresses the
spread). Association, not causation: a shared systemic cue could move both. **Stop rule:** one amendment. **Code:**
`paper/exploratory/A229_endothelial_gate/`.

**A229 power note (2026-11-07, before any observed association; `A229_endothelial_gate/output/{eligibility,power}.csv`, commit 87ebf4f).**
Null from y shuffled within cohort only (10,000). Tumour side: H229.1 ICAM1 190 patients / 21 cohorts, null 0.306 SD 0.035, MDE
0.086; H229.2 CD8 share 172 / 20, 0.315, 0.036, 0.089; H229.3 VCAM1 190 / 21 (MDE 0.083); H229.4 SELE 190 / 21 (0.087); H229.5
T/NK share 164 / 18 (0.093). All scored. Adjacent-normal (descriptive): 144-155 patients, MDE 0.094-0.097. Exclusions: 10
cohorts without a tumour endothelial pseudobulk, 3 with < 4 endothelial patients; GSE183219 out of H229.2 (1 patient with >= 50
T/NK cells); 3 cohorts out of H229.5 (no T_NK pseudobulk). Implied peak = arg(sum y z) in hours (`clocksync.stats.HOURS`).

## A230 (2026-11-07, user agreement; written before GSE266616 is downloaded or any of its counts are read, and before any epithelial offset is computed in TAURUS; seen: the GEO sample metadata of GSE266616 (14 patients with inflamed and adjacent-to-inflamed or non-inflamed samples at the same site, scout `exploratory/scouts_2026-11-07/ibd_scout.md`), TAURUS stromal results A45 / A53 (stromal coupling 0.80-0.91 in both states, phases unchanged), the mouse regenerating-colon result (Carmona-Alcocer et al., Nat Commun 2026; A21 archive: epithelium −8.95 h untimed vs −9.03 h timed)): is the epithelial clock offset from the patient's stromal host in inflamed human intestine? **Within-person contrast; single cell; GSE266616 blind primary, TAURUS non-blind secondary**
**Why.** Q25: local desynchrony after injury needs no cue conflict in mice (regenerating colon). Its human counterpart is
an epithelial offset against the same tissue's stroma that differs between inflamed and non-inflamed tissue of the same
patient.
**Data.**
- **GSE266616** (Crohn's disease, terminal ileum and ascending colon; epithelial and lamina-propria fractions in one
  library; 10x 3' v3.1; processed counts from GEO; `$SCCC_DATA/gse266616/` with PROVENANCE.md). Arms: inflamed vs the
  same patient's adjacent-to-inflamed or non-inflamed sample at the same site ("reference"); one pair per patient and
  site; if a patient has both sites, both pairs enter with patient as the permutation block.
- **TAURUS** (on disk, `taurus` dataset; inflamed vs non-inflamed, UC + CD): non-blind secondary.
**Labels and pseudobulks (fixed now).** No author annotation is used. Cells classed by the `cancer_sc` marker rules
(epithelial, fibroblast, endothelial, mural, myeloid, T_NK, B_plasma, mast); doublet-like cells (two class scores within
the rule's margin) dropped as in `cancer_sc`. Pseudobulk per sample x class, >= 50 cells; sample kept if median UMI per
cell >= 1,500. Unit (resection) = sample; donor table carries patient, site and arm.
**Model.** `clocksync.Dataset` with fixed stromal-vascular anchors (fibroblast, endothelial; mural if present in >= 60 %
of samples), strata = site; `Dataset.contrast(groups = arm, levels = ("reference", "inflamed"), targets =
["epithelial"], companion = "fibroblast" (read on its leave-one-out host), blocks = patient, n_perm = all within-patient
label swaps if <= 10,000 else 10,000 random, require_host_gate_each_arm = True)`.
- **H230.1 (primary):** epithelial offset change |offset_diff_h| (inflamed − reference) > permutation null (two-sided,
  within-patient swaps), **and** larger than the companion's offset change on the same permutations (`Contrast.module`,
  epithelial − fibroblast offset difference, two-sided p < 0.05).
- **H230.2 (secondary):** epithelial coupling change d* (inflamed − reference) ≠ 0, with the companion difference on the
  same permutations.
- **H230.3 (secondary, non-blind):** H230.1 in TAURUS (same rules, blocks = patient).
Holm over H230.2-H230.3.
**Gate.** Host-only gate in each arm (`require_host_gate_each_arm`); a failed gate makes the contrast not estimable.
**Power (dated note before scoring; no observed contrast computed).** Eligible patients (both arms with epithelium and
>= 2 anchors), median UMI, cells per class; offset MDE = 2.49 x SD of offset_diff_h under within-patient label swaps
(permutation draws only); coupling MDE via `power.mde`. **Not scored** if < 6 eligible patients, or offset MDE > 6 h
(the mouse effect is about 9-12 h), or the gate fails in either arm.
**Reading.** H230.1 pass: inflamed human intestinal epithelium runs at a different phase from its own stroma than
non-inflamed epithelium of the same patient: a human local desynchrony without cue conflict. Fail with MDE <= 6 h: a
mouse-sized shift is not present; smaller shifts not excluded. "Adjacent-to-inflamed" is not strictly non-inflamed,
which biases toward no difference. Not scored: the local clause stays supported by mouse colon and tumours only.
**Stop rule:** one amendment; no other label rule, anchor set or arm definition after the power note. **Code:**
`paper/exploratory/A230_ibd_local_offset/`.

**A230 power note (2026-11-07, before any observed contrast; `A230_ibd_local_offset/output/{power,gate,cells,eligible_samples}_*.csv`, commit 6710885).**
GSE266616: 17 inflamed / adjacent pairs (17 patients; the GEO `patient` field is malformed for 3, patient taken from the
title), all references adjacent-to-inflamed (none non-inflamed); all 85 samples pass median UMI >= 1,500; 8 patients
eligible (epithelium and >= 2 anchors in both samples); mural in 9 % of samples, so anchors = fibroblast + endothelial.
Host-only gate fails in both arms (p_pairing 0.052 reference, 0.170 inflamed; one anchor pair, 8 donors): **H230.1 and
H230.2 not scored** (registered rule). Offset MDE 3.17 h. TAURUS: 25 eligible patients (82 non-inflamed, 60 inflamed
samples), mural 45 % so anchors = fibroblast + endothelial; gate passes in both arms (p_pairing 0.001); offset MDE 1.71 h;
d* MDE 0.19. **Design flaw found at the power step:** with two anchors the fibroblast companion has no leave-one-out
host (needs >= 3 anchors), so the companion clause of H230.1 / H230.3 cannot be evaluated; no observed contrast has
been computed. `Dataset.contrast` draws random within-patient swaps (no exact enumeration). Decision on the companion
pending (user).

**A230 design change (2026-11-07, user decision; after the power note, before any observed contrast).** The fibroblast
companion of H230.1 / H230.3 is not estimable with two anchors. It is replaced by the project's standard specificity
check: the epithelial |offset_diff_h| (inflamed − reference) must exceed the 95th percentile of the same statistic over
200 expression-matched pseudo-clock sets (`clocksync.paired.pseudo_waveforms`, kind "matched", log means from the
genome-wide non-inflamed epithelial pseudobulks, seed 230), each fitted with the same anchors, strata and samples
(gate off for pseudo sets, as A209). H230.3 pass = permutation p < 0.05 (Holm with H230.2 is moot: H230.2 belongs to the
unscored GSE266616 arm, so H230.3 is tested alone) **and** clock > pseudo q95; p_clock = (1 + #{pseudo >= clock}) / 201.
The d* readout gets the same pseudo check (descriptive). GSE266616 stays not scored. This is a design change recorded
before scoring, not a post hoc reading; TAURUS remains non-blind (stroma seen in A45 / A53).

## A231 (2026-11-07, user agreement; literature follow-up after the robustness audit; written before any age-stratified fit; seen: donor counts by age arm and cohort only (phenotype table, aggregates), A19 / A180 / A130 results, A208 design and results): do organ coupling or organ phase relations change with age? **GTEx, A208 procedure with age as the label, secondary (Holm over H231.1-3)**
**Why.** Ageing is in the headline, but our support is a non-detection on own noise (A19 MDE 0.20 per decade; A180
equivalence failed, 90 % [-0.147, 0.089] vs ±0.0375 per decade); coupling and offsets were never contrasted by age.
Talamanca et al. (Science 2023) report rhythms that dampen with age. A130 internal minus wall time drifts -0.339 h per
decade (p 0.054).
**Data and labels:** `common.gtex` with the protected exact age (in memory only). Cohorts as A164 / A208: fast (Hardy
1-2) and vent (Hardy 0), kept separate. Old = 60-70 y (case), young = 21-49 y (control); 50-59 y not used. Counts in the
phenotype table (before tissue filters): fast 126 old / 48 young, vent 105 / 243.
**Procedure:** `compute_a208.py` with the label replaced: 1:1 matching within cohort, exact on sex and Hardy class,
nearest ischaemic time (the smaller arm processed first, as A208); tissues with >= 20 donors per arm in the cohort and
>= 10 matched pairs (all non-cultured tissues; no organ excluded); host = A164 anchors; 1,000 within-set label
permutations shared across tissues and statistics; 200 matched-set bootstraps with host refit; pooling by inverse null
variance within then across cohorts; 90 % interval = ± 1.645 x max(bootstrap SE, null SD) (A202).
**Power first (dated note before any observed statistic):** MDE_80 = 2.49 x null SD for the pooled d* and the pooled
offset; for Q (H231.3), the SD s of independent per-tissue offset shifts detected with 80 % power at the null q95.
- **H231.1 (coupling, equivalence):** pooled d* (old − young) within ±0.15 if MDE <= 0.15; else descriptive only.
- **H231.2 (offset, effect):** pooled offset difference, two-sided permutation p; clock-specific only if |pooled offset|
  exceeds the q95 of 50 expression-matched pseudo-clock sets (A208 companion). Equivalence within ±1 h reported.
- **H231.3 (reordering, effect):** Q = Σ_t w_t (offset_diff_t)² over tissues (w_t = 1 / null var), permutation p
  (one-sided, upper); clock-specific only if Q exceeds the q95 of the same pseudo-clock sets. This is the body-wide
  desynchrony readout (organs shifting by different amounts), which a pooled mean offset can miss.
**Companion for d*:** the 3 lowest-coupling tissues per cohort (A208).
**Descriptive, answering the "GTEx is timed" objection (no pass / fail):**
- **D231.4 timed comparator (fast cohort):** on the same matched donors, per tissue and arm, a cosinor fit of each clock
  gene on time of death with the GTEx waveform (tissue phase = the waveform phase shift maximising the likelihood), the
  tissue phase minus the arm's circular mean over tissues, age difference pooled as H231.2. Reported: the timed pooled
  offset and its bootstrap SE next to clocksync's, their agreement per tissue (circular correlation), and SE ratio.
- **D231.5 working hours:** H231.1-3 statistics and D231.4 recomputed on fast donors who died 08-18 h only: estimate,
  SE and agreement with the full-day values. Shows what each approach keeps when the window is narrow.
- **D231.6 host spread:** host spread (the common-spread denominator of d*) old vs young per cohort. A body-wide
  amplitude loss (Talamanca) does not change coupling or offsets; this readout is the closest, confounded in fast by
  the time-of-death distribution and in vent by the ICU state; descriptive only.
**Reading:** H231.1 pass: inter-organ coupling does not change with age within ±0.15 (the headline "ageing" clause
becomes an equivalence for coupling). H231.2 / H231.3 pass with companion: organs shift relative to the body clock
with age (direction, size); fail with MDE stated: not detected. Results per cohort are reported; the vent cohort is
the one no timed method can use (wall-time entrainment R 0.15 vs 0.53 in fast deaths, A109c).
Controlled access: aggregates only (per tissue, per cohort, pooled); no per-donor ages or phases written.
**Stop rule:** one amendment; no age cut, tissue subset or cohort chosen after the power note. **Code:**
`paper/exploratory/A231_age_coupling/`.

**A231 power note (2026-11-08, before any observed statistic; `A231_age_coupling/output/{power,balance_summary,eligible_summary}.csv`).**
Eligible tissue x cohort units: vent 25 (1,212 matched pairs, median 47 per tissue), fast 29 (649, median 23); pooled 54
units, 1,861 pairs. MDE_80: pooled d* 0.037 (vent 0.051, fast 0.052), far below ±0.15, so H231.1 is powered; pooled
offset 0.20 h (vent 0.28, fast 0.27); Q: per-tissue shift SD 0.45 h (vent 0.51, fast 0.64; null Q mean 54.0, q95 73.5).
Balance after matching (median |SMD| over tissues): ischaemic time 0.05 / 0.09 (fast / vent), BMI 0.13 / 0.15; T2D is
not a matching variable and stays imbalanced (0.46 / 0.79, max 1.03), older donors carrying more T2D. T2D leaves
coupling unchanged within ±0.15 (A164 / A202), so this is stated as a caveat, not adjusted. Implementation choices in
the docstrings of `compute_a231.py` and `matched_label.py` (score stage checked end to end on placebo labels only).
The power stage ran twice concurrently by mistake; a single rerun from the cached null draws gave byte-identical
`power.csv`. Scoring proceeds as registered.

## A232 (2026-11-07, user agreement; literature follow-up item 3; written before any proliferation score is computed or related to an offset; seen: the GTEx organ offsets (Figure 1 organ map, A92 / A123), never against proliferation): do organs with more proliferation run at a different phase? **GTEx, secondary (Holm over H232.1-3)**
**Why.** Droin, Paquet & Naef (Nat Phys 2019) explain shifted circadian phases of tissues by their proliferation rates
(cell cycle coupled to the clock); in mouse, tissue phases correlate with Mki67 and Myc. In humans, no within-donor test
exists. A97 / A126 / A127 / A166 tested coherence vs proliferation, never the offset.
**Data:** GTEx v10 (`common.gtex.load`), cohorts fast and vent separate, all non-cultured tissues with >= 50 donors in the
cohort; clock fit as the organ map (A164 anchors, target anchors on their leave-one-out host), offset_t = arg z_t (h),
centred on the circular mean over tissues.
**Proliferation:** per sample, the mean z (across all samples of the cohort, so that tissues are
comparable) of log CPM of the A97 G2M set (Hallmark G2M, clock genes removed); tissue score = median over its donors.
MKI67 and MYC log CPM medians reported alongside (secondary readouts, not tested).
- **H232.1 (primary of A232, fast cohort):** Spearman rho between offset_t and tissue G2M score across tissues ≠ 0;
  p by 10,000 permutations of the scores over tissues (two-sided). Weighted check (descriptive): weighted slope with
  w = 1 / bootstrap var of offset_t.
- **H232.2 (vent cohort):** the same; replication reading = same sign and p < 0.05.
- **H232.3 (within tissue):** in each tissue of the fast and vent cohorts, donors in the top vs bottom tertile of that
  tissue's G2M score (within cohort), `Dataset.contrast` offset_diff_h (top − bottom), 1,000 label permutations within
  tissue, pooled by inverse null variance (A208 pooling); two-sided p. Reading: same sign as H232.1 or not.
**Specificity (each test, same data):**
- **clock:** 200 expression-matched pseudo-clock sets (`clocksync.paired.pseudo_waveforms`, kind "matched", seed 232),
  each through the whole pipeline (offsets per tissue, same scores); the clock |rho| (H232.1-2) or |pooled offset|
  (H232.3) must exceed the pseudo q95 (a composition artefact would move pseudo-clock offsets too).
- **proliferation:** 200 random gene sets of the G2M set's size, matched on mean expression (deciles), clock and G2M
  genes excluded, scored the same way; the G2M |rho| must exceed their q95 (H232.1-2).
A test is called clock- and proliferation-specific only when p (Holm) < 0.05 and both companions pass.
**Power (dated note before any rho):** number of tissues per cohort, offset bootstrap SEs, MDE of rho at n tissues
(80 %, two-sided 0.05) and of the H232.3 pooled offset (2.49 x null SD).
**Reading:** pass: human organs that proliferate more run at a shifted phase relative to the body clock (direction and
slope in h per unit score), as predicted by cell-cycle-clock coupling; this also bears on the cancer chapter (the
proliferation-linked coherence loss). Fail with MDE stated: not detected at this size. Tissue composition differs
between tissues and is not removed; the pseudo-clock companion guards against expression-driven artefacts only.
**Stop rule:** one amendment; G2M set, cohorts and tissue filter fixed here. **Code:**
`paper/exploratory/A232_offset_proliferation/`.

**A232 power note (2026-11-07, before any rho or observed top-vs-bottom offset; `A232_offset_proliferation/output/power.csv`,
`tissue_offsets.csv`).** Tissues with >= 50 donors: fast 38 (median 150.5 donors), vent 30 (221). Median offset
bootstrap SE 0.195 h (fast) and 0.197 h (vent); between-tissue offset SD 1.61 h and 0.61 h, so offset reliability 0.98
and 0.70. MDE of Spearman rho (80 %, two-sided 0.05; Bonett-Wright approximation, checked by simulation): fast 0.46,
vent 0.52 (0.62 after attenuation by offset noise): only strong associations are detectable, tissues being the unit.
H232.3 pooled within-tissue offset MDE 0.13 h (fast 0.18, vent 0.21; null SD 0.053 h pooled). G2M set: MSigDB Hallmark
G2M v2024.1 (200 genes after removing clock genes, 199 in GTEx); A97 used v2026.1, not on disk (version differs).
Implementation choices in the docstring of `compute_a232.py`. All tests powered as registered; scoring proceeds.

## A233 (2026-11-07, user agreement; robustness audit row 4; written before any condition-stratified fit; seen: case / non-case counts per condition and cohort (phenotype table, aggregates), A164 / A208 designs and results): do chronic conditions other than T2D and obesity change organ coupling or phase relations? **GTEx, A164 procedure, secondary (Holm per hypothesis over conditions)**
**Why.** "Chronic disease" in the headline rests on T2D and obesity (A164 / A202) and per-row atlas TOSTs (A180). The
literature reports clock disruption in hypertension (non-dipping blood pressure, vascular clocks), ischaemic heart
disease, COPD and cigarette smoke (lung clock), and major depression (brain, Li et al. 2013). Citations to be pinned
in `literature_positioning.md` before writing; they do not change the design.
**Conditions (protected phenotypes; case = 1, control = 0, unknown excluded)** with counts in the phenotype table
(vent / fast cases): hypertension MHHTN (188 / 152), ischaemic heart disease MHHRTDIS (56 / 41), COPD MHCOPD (69 / 48),
major depression MHDPRSSN (63 / 18), smoking MHSMKSTS Yes vs No (259 / 169). Sepsis (9 / 1) and others below the A164
eligibility are not used.
**Modules (local readout; fixed here):** hypertension: heart (LV, atrial appendage), arteries (aorta, coronary,
tibial), kidney cortex; ischaemic heart disease: heart (LV, atrial appendage), coronary artery, aorta; COPD and
smoking: lung; depression: all brain regions. Companion: tibial nerve and the two skins (A164 COMP), any of them in the
module (none) excluded.
**Matching (A164 T2D rule, blind to clock data):** each case to up to 3 controls without replacement, exact on sex and
Hardy class, nearest age (caliper 5 y), then ischaemic time; balance (SMD of age, ischaemic time, BMI, T2D share)
reported before any clock readout. A tissue enters with >= 10 matched cases.
**Estimation:** A164 / A208 (`Dataset.contrast` d* and offset_diff_h per tissue, permutations within matched sets
(1,000, shared), bootstrap (200, host refit), inverse-null-variance pooling within then across cohorts, A202 interval).
**Power first (dated note per condition before any observed statistic):** MDE_80 of the body-wide pooled d*, module d*,
module offset; Q as A231.
Per condition c:
- **H233.c.1 (body-wide coupling, equivalence):** all-tissue pooled d* within ±0.15 (if MDE <= 0.15; else descriptive).
- **H233.c.2 (body-wide reordering):** Q over all tissues (A231 H231.3 definition), one-sided permutation p, Holm over
  conditions; clock-specific only above the q95 of 50 pseudo-clock sets.
- **H233.c.3 (local):** module d* and module offset, two-sided permutation p, Holm over conditions (per statistic);
  "organ-specific" only if module − companion is also p < 0.05 (d*) and the offset beats the pseudo-clock q95 (offset).
**Reporting per cohort:** vent and fast separately and pooled; the vent cohort is where untimed reading is the only option.
**Reading:** H233.c.1 pass: body-wide coupling unchanged in condition c within ±0.15; H233.c.2 / c.3 pass with
companions: body-wide reordering or a local organ shift (named). Fail with MDE: not detected. Comorbidity is not
removed (hypertension with T2D, smoking with COPD); conditions are not independent and are read separately.
Controlled access: aggregates only.
**Stop rule:** one amendment; conditions, modules, matching and cohorts fixed here; no condition added after the power
note. **Code:** `paper/exploratory/A233_chronic_conditions/` (shares the matched-contrast code with A231).

**A233 power note (2026-11-08, before any observed statistic; `A233_chronic_conditions/output/{power,balance_summary,eligible_summary}.csv`).**
Eligible tissues vent / fast and matched cases: hypertension 30 / 39 (1,987 / 1,675), ischaemic heart disease 26 / 35
(793 / 864), COPD 29 / 35 (806 / 758), depression 29 / 9 (848 / 98), smoking 30 / 40 (2,036 / 1,547). Pooled MDE_80
(body-wide d*, offset h, Q shift SD h; module d*, offset h; module tissues): hypertension 0.028, 0.14, 0.31; 0.054,
0.28 (10); ischaemic heart disease 0.035, 0.18, 0.44; 0.087, 0.51 (7); COPD 0.035, 0.18, 0.44; 0.222, 1.45 (lung, 2);
depression 0.048, 0.28, 0.58; 0.431, 3.07 (1 brain tissue unit); smoking 0.030, 0.15, 0.34; 0.166, 1.21 (lung, 2).
H233.c.1 is powered for all five conditions. The lung (COPD, smoking) and brain (depression) local readouts are weak
(module d* MDE > 0.15, offset MDE > 1 h); they are scored as registered and read with their MDE. Balance after matching
(median |SMD|): age <= 0.16 except hypertension in vent (0.54, max 0.60; the 5 y caliper does not remove it); BMI 0.49
(COPD fast) and 0.35 (hypertension fast, max 0.97); ischaemic time 0.40 (COPD vent); T2D 0.66 / 0.72 for hypertension
(vent / fast). Comorbidity is not removed, as registered; A231 reads the age contrast directly. Implementation choices
in the docstring of `compute_a233.py`. The power stage ran twice concurrently by mistake (later copy stopped); a single
rerun from the cached null draws gave byte-identical `power.csv`. Scoring proceeds as registered.

## A234 (2026-11-08, user agreement; follow-up of A231; written before any internal-time-by-age readout in GTEx; seen: A231 per-tissue offsets (fast brain regions old − young mean -1.15 h, 9 of 10 negative; other tissues -0.06 h), fast-donor counts by Hardy class x exact-age band and by brain region (aggregates), A130, archive `legacy/desync/01_gtex_organs` (timed phase advance 0.35 h per decade, bracket age, "young fast deaths are mostly violent, Hardy 1")): does internal time drift from external time with age, and does the brain drift differently from the body? **GTEx fast deaths, secondary (Holm over H234.1 and H234.3; H234.2 post hoc)**
**Why.** Chronotype advances with age (epidemiology). The ageing evidence so far: A130 in living humans (-0.339 h per
decade, p 0.054) and an archive GTEx cosinor estimate that is entangled with death type (fast deaths under 40 y: 14
Hardy 1, 7 Hardy 2; 40-70 y: 20 Hardy 1, 203 Hardy 2) and used bracket age. A231 cannot see a shared drift (offsets are
relative to the donor's body clock), and its brain offsets suggest the brain falls behind the body with age.
**Sign convention (fixed here).** Internal time = the clocksync host phase of the donor (h); Δ = internal − recorded
time of death (circular, h). The host phase increases with time of day (waveform convention; A113 H113.1: single-sample
phase vs ZT Fisher-Lee +0.82), so Δ > 0 means the internal clock is ahead of the wall clock, and **a positive slope on
age is a phase advance with age** (earlier chronotype). Checked in the data before any age readout (D234.0). Note: A130
registered "slope < 0 (earlier with age)" with Δ computed the same way (internal − wall), so its label looks inverted
(its -0.339 h per decade would be a delay); to be confirmed in the A130 code and corrected in `results_preliminary.md`.
**Data.** `common.gtex.load`, fast deaths (Hardy 1-2) with time of death; protected exact age (`AGE`) and sex in memory
only. Counts: 244 fast donors (22-70 y), 241 with >= 5 body tissues, 197 with >= 3 brain regions.
**Hosts (two separate fits on the fast cohort, untimed).**
- Body host: `cs.Dataset` on the A81 / A92 body tissues (`_shared.BODY`, 12 peripheral tissues; no brain, no pituitary)
  as anchors; donors with >= 5 of them.
- Brain host: brain regions with >= 100 fast donors, spinal cord excluded (follows the body, A143): amygdala, anterior
  cingulate, caudate, cerebellar hemisphere, cerebellum, cortex, frontal cortex BA9, hippocampus, hypothalamus, nucleus
  accumbens, putamen (11), all as anchors; donors with >= 3 of them.
- Host phase per donor from `Fit.host`; each fit's constant frame offset cancels in the slope.
**Estimator (all three tests).** Circular-linear regression Δ_i = μ_s(i) + β·age_i (von Mises likelihood, i.e. maximise
Σ cos(Δ_i − μ_s − β age_i)), stratum s = Hardy class x sex; β in h per decade. Null: age permuted within strata
(10,000; two-sided p). 90 % interval ± 1.645 x max(donor-bootstrap SE (200, hosts refitted), null SD). MDE_80 = 2.80 x
null SD (two-sided 0.05).
- **D234.0 (direction check, before any age readout):** Fisher-Lee correlation of each host phase with time of death
  must be positive; if negative, that host's phases are negated before Δ. Wall-time alignment R per host reported.
- **H234.1 (body, timed reference):** β of Δ_body = body host − time of death on age ≠ 0. Read against the archive prior
  (+0.35 h per decade if the archive sign holds). Sensitivity: Hardy 2 only (176-210 donors, few under 40); bracket age
  instead of exact age (archive comparability).
- **H234.2 (brain vs body, untimed; POST HOC: the pattern was seen in A231):** β of Δ_bb = brain host − body host on age,
  donors with both hosts, each donor counted once (answers the shared-donor problem of the A231 per-region readout).
  Companions: 50 expression-matched pseudo-clock sets (`clocksync.paired.pseudo_waveforms`, kind "matched", seed 234)
  through both hosts, |β| must exceed their q95; brain host without the two cerebellar regions (anchor dependence, WP0
  item 8). A231 predicts β < 0 (brain behind body with age; -1.15 h over the old - young midpoint gap of about 3
  decades). Reported as post hoc whatever the result.
- **H234.3 (brain, timed reference):** β of Δ_brain = brain host − time of death on age ≠ 0.
**Power (dated note before any observed β):** donors per test, null SD and MDE_80 of β, R of each host with time of death.
**Reading.** H234.1 pass: body internal time drifts with age relative to clock time (advance if β > 0), with death type
and sex held fixed; fail with MDE: not detected at this size. H234.3 with H234.1 and H234.2 decide whether brain and body
drift together, apart (brain keeps clock time while the body advances, or the reverse), or neither. Time of death has
recording errors and fast deaths are not living subjects; constant offsets cancel, age-related reporting errors do not.
Controlled access: aggregates only (β, intervals, nulls); no per-donor phases or ages written.
**Stop rule:** one amendment; hosts, strata, age variable and donor filters fixed here. **Code:**
`paper/exploratory/A234_internal_time_age/`.
**A234 power note (2026-11-08, before any observed β; `A234_internal_time_age/output/{power,direction,counts}.csv`).**
D234.0: both hosts run with time of death in the registered direction (Fisher-Lee body 0.380, brain 0.372; none
negated); wall-time alignment R 0.62 (body, 233 donors) and 0.62 (brain, 196). Donors: H234.1 233 (Hardy 1 / 2: 31 / 202;
19 under 40 y), H234.2 192, H234.3 196. Null SD of β (10,000 age permutations within Hardy x sex) and MDE_80: H234.1
0.249 / 0.70 h per decade, H234.2 0.234 / 0.66, H234.3 0.302 / 0.85. **Underpowered against the expected sizes:** power
about 0.29 at the archive prior (0.35 h per decade) for H234.1, about 0.37 at the A231-implied brain − body drift
(-0.38 h per decade) for H234.2; only drifts of about 0.7 h per decade or more are detectable. Scored as registered and
read with the MDE. The pseudo-clock companion is weak: pseudo hosts carry no time-of-day signal, so their β scatter up to
the estimator bound (±3.05 h per decade; placebo check). Code checked end to end on shuffled ages only (`--placebo`).

## A235 (2026-11-09, user agreement; follow-up of the mouse high-fat-diet literature; written before any metabolic-organ offset is computed for T2D or obesity; seen: A164 coupling results (T2D module d* +0.019, obesity +0.004; `a164_tissues.csv` holds no offsets), A233 offsets for five other conditions, A164 matching and balance tables; the archive A47 outputs were not re-read for offsets): are metabolic organs phase-shifted relative to the body clock in T2D and obesity? **GTEx, A164 matching, A233 estimation; fast deaths primary, vent negative control (Holm over T2D and obesity)**
**Why.** In mice a high-fat diet advances peripheral clocks, the liver most clearly, against an almost unchanged SCN
(Pendergast et al. 2013; also Kohsaka et al. 2007, Eckel-Mahan et al. 2013; from memory, to be pinned in
`literature_positioning.md` before writing; they do not change the design). The shift is carried mostly by food eaten in
the rest phase (Hatori et al. 2012). A164 tested coupling only, so the mouse prediction (an organ offset, not a loss of
coupling) was never tested in humans.
**Expected direction:** module tissues advanced relative to the host in cases (offset_diff_h > 0; `Dataset.contrast`
convention: case minus control, > 0 advanced). The expected human size is unknown and probably far below the mouse
shift: obese people do not systematically eat in their biological night, and terminal illness removes the feeding cue.
**Labels and matching (A164, blind to clock data):** T2D cases vs non-diabetic controls, up to 3 controls per case;
obesity BMI >= 30 vs 18.5 <= BMI < 25, T2D donors excluded, 1:1. Exact on sex and Hardy class, age caliper 5 y,
nearest age then ischaemic time. Balance (SMD of age, ischaemic time, BMI, T2D share, sex) reported before any clock
readout. A tissue enters with >= 10 matched cases in the cohort.
**Module:** liver, subcutaneous and visceral adipose, pancreas, skeletal muscle (A164 MET). **Companion:** tibial nerve
and the two skins (A164 COMP). Host: A164 anchors; a module tissue that is an anchor is read on its leave-one-out host.
**Estimation (A233):** per tissue `Dataset.contrast` offset_diff_h and d*; 1,000 permutations within matched sets,
shared across tissues; 200 matched-set bootstraps with host refit; inverse-null-variance pooling over module tissues;
90 % interval ± 1.645 x max(bootstrap SE, null SD).
**Power first (dated note before any observed statistic):** MDE_80 = 2.49 x null SD of the fast module offset, vent and
pooled, per condition.
- **H235.1 (primary, advance):** fast-cohort module offset > 0, one-sided permutation p, Holm over T2D and obesity.
  Clock-specific only if it exceeds the q95 of |module offset| over 50 expression-matched pseudo-clock sets (fast,
  A233 companion). "Metabolic-specific" only if the fast module − companion offset, on the same permutations, is also
  > 0 at one-sided p < 0.05. The two-sided p and the equivalence within ±1 h are reported.
- **H235.2 (negative control, descriptive):** vent-cohort module offset and its 90 % interval (no feeding rhythm in
  ventilated care; prediction: no shift). A shift of the same size in vent argues against the feeding route.
**Descriptive (no pass / fail):** pooled (vent + fast) module offset; per-tissue offsets (liver named); module log z2
ratio (relative amplitude, the damping readout, case vs control) per cohort.
**Sensitivity (registered, not a test):** H235.1 recomputed with a host that excludes the module tissues (anchors =
A164 anchors minus skeletal muscle and the two adiposes), so a shift shared by all metabolic organs is not partly
absorbed by the host.
**Reading:** H235.1 pass with companions: metabolic organs run advanced relative to the body clock in condition c,
consistent with the mouse result (size stated against the mouse shift). Fail with MDE: a shift larger than the MDE is
not present in fast deaths; a mouse-sized shift is excluded only if the 90 % interval excludes it. The feeding route
is not tested beyond H235.2. Controlled access: aggregates only.
**Stop rule:** one amendment; conditions, module, cohorts and sign fixed here; no tissue subset chosen after the power
note. **Code:** `paper/exploratory/A235_metabolic_offset/` (A231 / A233 matched-contrast code).

**A235 power note (2026-11-09, before any observed statistic; `A235_metabolic_offset/output/{power,power_mc,matching,balance_summary,eligible_summary}.csv`).**
Module tissues per cohort: fast liver, both adiposes and skeletal muscle (pancreas out: 5 T2D and 3 obese cases);
vent all five. Matched cases in the module (fast / vent): T2D 156 / 201, obesity 111 / 240. MDE_80 of the module
offset (A164 host; fast / vent / pooled): T2D 0.37 / 0.53 / 0.30 h, obesity 0.64 / 0.82 / 0.50 h; fast module −
companion offset 0.71 h (T2D) and 1.03 h (obesity). Sensitivity host without the module anchors: fast 0.42 h (T2D) and
0.70 h (obesity). H235.1 is powered for shifts of about 0.4 h (T2D) and 0.6 h (obesity), well below the multi-hour
liver advance of high-fat-fed mice; the expected human size is unknown. Module d* MDE <= 0.16 everywhere (descriptive
here, A164 tested coupling). Balance after matching (median |SMD|): age <= 0.25, ischaemic time <= 0.23; BMI
imbalanced in the T2D design (fast 0.68, vent 0.23; cases heavier, as in A164; T2D and obesity are not separable in
this design and are read together). Code checked end to end on placebo labels only (`--placebo`, subsample). Scoring
proceeds as registered.

## A236 (2026-11-10, user agreement; follow-up of A108 / A108b / A108c; written before any gene-level time-of-day coefficient, any signed-area or polarisation statistic, any responsive pseudo-set and any matched pseudo-set beyond the three of A108 are computed, and before the genome-wide pseudobulks are built; seen: every A108 / A108b / A108c / A121 result (H108.1 R 0.516, p 0.016; gate p_pairing 0.001, p_acrophase 0.083; three matched pseudo sets R 0.19 / 0.33 / 0.19; D1 glucocorticoid score p 0.26) and a post hoc look on 2026-11-10 at the host states per clock time (host-state covariance eigenvalue ratio 0.058; median error 0.4 h at 0 h against 2.0-5.5 h at 6, 12 and 18 h; axial concentration 0.67): is the time-of-day state of the late mouse embryo a running clock? **Qiu 2024, E15.0-E18.75, 16 embryos; secondary (development section); H236.2 is post hoc in motivation**
**Why.** A108 showed that the untimed host phase of late embryos follows the time of collection (R 0.516). The user's
reading: maternal cues act as a zeitgeber on a clock that already runs, as light does in adults. The post hoc look
says the data do not show this yet: the host states lie almost on one line and the alignment is carried by the 0 h
embryos, which is what a single cue-driven axis would give (R ~ 0.5 at four evenly spaced times). Three properties
separate "a time-of-day state" from "a running clock", and none was tested: the gene content of the time-locked
signal, its second dimension (rotation), and its specificity against cue-responsive genes. Entrained and driven
oscillators cannot be separated cross-sectionally: no outcome below shows self-sustained rhythms.
**Data and design (fixed by A108):** late window, 16 embryos in 4 batches (one E-day and its sequencing runs each),
each batch with one embryo at 0, 6, 12 and 18 h; cell types by the A108 rule (>= 50 cells, present in >= 75 % of the
embryos; 100, of which 77 in all 16); 14 clock genes; GTEx global waveform, as in A108.
**Timed model (all three parts).** Per cell type, a response is regressed on batch intercepts, a linear within-batch
drift (t - 9 h) and cos wt, sin wt (w = 2 pi / 24 h). The drift term is there because age and clock time move
together within a batch (A108 limits): a developmental trend plus a one-dimensional oscillation would otherwise
imitate a rotation. With four evenly spaced times the drift is identified by the alternating contrast (0 - 6 + 12 -
18 h), so a 12-h harmonic is assumed absent. The model without the drift term is reported as a sensitivity.
- **H236.1 (primary, gene content):** the 24-h harmonics of the 14 clock genes lie in the clock plane. Response:
  log fractions centred over genes (`centred_log_fraction`); per gene, common cos / sin coefficients (a_g, b_g) over
  cell types with cell-type-specific intercepts and drifts, weights `gene_weights` of the cell type. Statistic:
  S1 = (a'Pa + b'Pb) / (a'a + b'b), P the projector on the clock plane (`_core.study.plane`). Null: 10,000 waveforms
  with acrophases permuted over genes (`Waveform.permuted`), one-sided. S1 does not depend on rotation: an axis in
  the clock plane passes.
  **H236.1b (secondary, beyond the direct glucocorticoid targets):** the same without PER1 and PER2 (12 genes,
  re-centred), p < 0.05; with it the H108.1 statistic without the two genes (`drop_genes`), R and permutation p.
- **H236.2 (primary, rotation; post hoc in motivation):** the time-locked component turns in the adult direction.
  Response: clock-plane state of each cell type as a complex number; fitted c_k (cos) and s_k (sin). Polarisation
  rho = 2 sum_k w_k Im(s_k conj(c_k)) / sum_k w_k (|c_k|^2 + |s_k|^2), w_k = 1 / residual mean square of cell type k
  (full fit, fixed): 1 for a circular rotation in the adult direction, 0 for an axis, -1 for the reverse; the
  numerator has expectation zero under an axis at any angle and under no signal, whatever the noise covariance; it
  needs no phase alignment between cell types. Test: leave-one-batch-out jackknife t of rho (4 batches, 3 df),
  one-sided rho > 0. Time permutations are not used: they are not valid for the composite null "axis".
- **H236.3 (primary, specificity):** H108.1's R exceeds the q95 of both 200 expression-matched and 200 responsive
  pseudo-clock sets (`clocksync.paired.pseudo_waveforms`; glucocorticoid and immediate-early genes), each scored as
  H108.1 on genome-wide pseudobulks. "Clock-specific" only above both. The genome-wide cell type x embryo
  pseudobulks of the late window are rebuilt from the GEO count matrices (six runs, streamed and parsed only, no raw
  matrix stored; clock-gene columns checked for equality with the stored panel; genes with mean count >= 1 enter the
  matching). If the rebuild fails, the 67-gene panel gives a descriptive answer only.
Holm over the three primary tests.
**Power first (dated note before any observed statistic):** gene-level plasmode on the cell types present in all 16
embryos: fitted intercepts and drifts, residuals of the timed model reassigned between embryos as whole blocks
(random permutation and sign, rescaled for the fitted degrees of freedom), and a planted clock-plane signal of the
noise-corrected fitted 24-h power per cell type. Scenarios: no signal and axis (size of H236.2), rotation with rho 1
and 0.5 (power of H236.2), axis (power of H236.1), harmonics of the same size on a permuted-acrophase design (size of
H236.1). If the size of H236.2 exceeds 0.07, its critical value is replaced by the simulated 95 % quantile, stated in
the note. The observed rho, its numerator and S1 are not printed before the note is committed.
**Descriptive (no pass / fail):** per-gene amplitude and acrophase against the reference (BMAL1 arm against the
E-box / D-box arm); leave-one-gene-out range of S1 and R; rho with a noise-corrected denominator, per trajectory
(hepatocytes and cardiomyocytes named, where a clock is expected first) and leave-one-batch-out; share of pseudo sets
with a timed p < 0.05. **Sensitivity:** all three without the drift term; H236.2 on the 77 complete cell types.
**Reading:** all three pass: the late embryo carries a time-of-day-locked state with the gene content of the clock,
turning in the adult direction, not shared by cue-responsive genes: a running clock at population level; entrained
or driven stays open. H236.1 and H236.3 pass, H236.2 fails: a clock-shaped, clock-specific time-of-day axis; rotation
not shown, read against the power note (a failed test is uninformative below the stated power). H236.1 fails: the
wording stays "time-of-day state". H236.1 passes and H236.1b fails: the axis is carried by PER1 and PER2, compatible
with a direct glucocorticoid drive. H236.3 fails: the state is shared with cue-responsive genes.
**Limits:** four batches, one embryo per time and batch; batch, age and sequencing run are confounded; no light
schedule recorded. **Stop rule:** one amendment; statistics, gene sets, weights and the drift model fixed here; no
lineage chosen after the power note. **Code:** `paper/exploratory/A236_embryo_running_clock/`.

**A236 power note (2026-11-10, before any observed rho, numerator or S1; `A236_embryo_running_clock/output/{power,power_info}.csv`, 2,000 data sets per scenario).**
Plasmode on the 77 cell types present in all 16 embryos: fitted batch intercepts and drifts kept, gene-level residual
blocks (embryo x cell type x gene) of the timed model rescaled by sqrt(16 / 9) and reassigned between embryos (random
permutation and sign). **Deviation 1 (planted signal):** the noise-corrected fitted 24-h power is positive in only 35
of the 77 cell types, so planting it per cell type along the fitted orientation would select on noise. The signal is
planted instead with one orientation common to the cell types and a 24-h power of lambda x the expected noise power
of each cell type's harmonics, lambda = 0.25, 0.5, 1, 2 and the observed ratio 0.72 (untruncated, summed over cell
types; the only observed quantity computed).
**Deviation 2 (null of H236.1, forced by the plasmode):** the registered acrophase-permutation null is not valid here.
With no time-locked signal it rejects in 0.41 of the data sets, and in 0.07-0.20 with harmonics on a permuted-acrophase
design: the untimed variation between embryos is itself clock-shaped (the A108 gate), and with 16 embryos it leaks
into the fitted harmonics by chance, so S1 is high without any time-locked signal (null mean 0.37). A permutation of clock times would have the same defect in reverse (no power). The null of S1 is
therefore the plasmode without time-locked signal (10,000 data sets, residual blocks reassigned between embryos), on
the 77 complete cell types where it is defined: critical value 0.651, size 0.054, rejection 0.029-0.042 under
non-clock harmonics, power for a clock-plane axis 0.65 / 0.92 / 0.98 / 1.00 / 1.00 at lambda 0.25 / 0.5 / 0.72 / 1 /
2. H236.1 now reads: the time-locked harmonics are more clock-shaped than the untimed variation of the same embryos
leaves by chance. The same null (rebuilt for each gene set and model) serves H236.1b and the no-drift sensitivity;
the acrophase p is reported beside it, and S1 on all 100 cell types is descriptive. With cell-type-specific
amplitudes the gene weights hold S1 of a pure clock-plane signal slightly below 1 (0.996 for the planted axis).
**H236.2, critical value replaced as registered:** with t_3 the size is 0.047 with no signal but 0.12-0.14 under an
axis and 0.07-0.11 under non-clock harmonics (four batches; the jackknife underestimates the variance of a ratio).
The critical value is the largest simulated 95 % quantile over the null scenarios, t = 4.76 (axis, lambda 0.5),
instead of 2.35; the p value is the rank of the observed t in that scenario's 2,000 draws. Size at this value: 0.007
(no signal), 0.035-0.050 (axis), 0.019-0.033 (non-clock harmonics). Power: rho 1, 0.16 / 0.25 / 0.33 / 0.40 / 0.58;
rho 0.5, 0.07 / 0.12 / 0.14 / 0.16 / 0.19 at the five lambda. **H236.2 is underpowered at the observed signal size
(0.33 for a perfect rotation): a failure is uninformative and will not be read as "an axis".** Expected rho is
attenuated by noise in the denominator (mean 0.37 for a true rho of 1 at lambda 0.72), so the noise-corrected rho is
the estimate to read. **Deviation 3:** because the critical value is defined on the 77 complete cell types, they
carry the primary H236.2 and the 100 cell types become the sensitivity (the registration had them the other way
round); the sensitivities are read against the same critical value (approximate).
**Clarification (H236.3):** R has no drift term, so its "no drift" sensitivity is void; its counterpart is R after
residualisation on embryonic day (the A108 sensitivity), for the clock set and every pseudo set.
Scoring proceeds with these three deviations, all fixed before any observed statistic.

## A237 (2026-11-10, user agreement: "you can go ahead, you can download data if needed ... if you have ideas of datasets / ways to find clocks in published human embryos"; registered together with A238 and A239, before any clock statistic is computed on Cardoso-Moreira 2019, Garcia-Alonso 2022 or Kang 2011, prenatal or postnatal; seen: the sample sheets and design tables of the scout (`scouts_2026-11-10/FETAL_SCOUT.md`), the detection (median expression) of the clock genes per organ and period (`prep_fetal_bulk/output/clock_detection_cardoso.csv`), the linking of libraries to individuals (expression of non-clock genes only), the adult feasibility below, every A189 / A181 / A236 result, a post hoc re-read of A189 with the signed coherence (`A189_fetal_within/output/h189_1_signed_posthoc.csv`) and the shift simulation below): do human prenatal organs carry a clock, alone and shared across the organs of one fetus? **Cardoso-Moreira 2019, bulk RNA-seq, 4-20 wpc, 7 organs, 207 prenatal libraries; secondary (development section)**
**Why.** The human developmental result of the paper is one organ in one cohort (A189, fetal adrenal, 7 fetuses) and
the mouse counterpart was downgraded (A236). The scout found three open designs in other organs with enough
individuals; this one is the only one with several organs per fetus (the paper's question, asked in development:
A181 could not estimate it in Cao 2020) and with postnatal samples of the same organs on the same platform.
**Common design of A237-A239 (`paper/common/development.py`; fixed before scoring).**
1. *Adjustment (A236 lesson).* Log fractions, log((count + 0.5) / depth), are residualised gene by gene within each
   compartment on the covariates of the amendment (age trend, sex, ...); every pseudo-clock gene set receives the
   same adjustment. Age-related change of the clock genes is thereby not read as a clock, and an age-dependent clock
   amplitude is not modelled.
2. *L1, one compartment across individuals:* signed coherence C (`clocksync.decompose.signed_coherence`, added for
   this purpose): the coherence of the clock genes with the reference waveform after a per-sample offset is fitted
   out, sign kept. Null: 10,000 per-gene permutations of the residuals, one-sided. Reason: the unsigned R used so
   far (A189) reads a shift common to the clock genes of a sample (depth, composition, RNA quality) as coherence. In
   simulation without any clock (`feasibility_fetal/output/shift_artefact.csv`, 200 data sets per scenario) a common
   shift of sd 0.2 passes the permutation test of R in 0.625 (n 10) and 1.000 (n 40) of the data sets, against 0.050
   and 0.035 for C; with a clock both detect it (C 1.000 at every shift; R 0.690 at n 10 with the shift). R is
   reported beside C, not tested.
3. *L2, the compartments of one individual:* `Dataset.gate` on the adjusted table, free host, compartments with >= 5
   individuals, strata as stated (cross-covariances centred and pairings permuted within strata), 10,000
   permutations. The p of the readout is the larger of p_pairing and p_acrophase.
4. *Companions (both required):* 200 expression-matched and 200 responsive pseudo-clock sets
   (`clocksync.paired.pseudo_waveforms`; genes with mean count >= 1 in the samples of the test), adjusted and scored
   as the clock set: C for L1; for L2 the gate's `pairing_z` (shared clock-plane covariance standardised by its
   pairing null; 1,000 permutations, same seed). Clock-specific: the clock statistic exceeds the q95 of both.
5. *Verdict (`development.verdicts`, families written in each scoring script):* Holm within the L1 family of a
   hypothesis; an L2 readout is a hypothesis of its own. PASS: Holm p < 0.05 and clock-specific. p < 0.05 without
   clock-specificity is reported as "not clock-specific" (a state shared with other genes, not shown to be the
   clock). A gate that is not estimable has no verdict.
6. *Negatives:* a readout is negative when its unadjusted p is >= 0.05; it is informative only where the adult pass
   rate of the readout is >= 0.8 (D1b rule) and its positive control, where one is defined, passes. Even then it
   reads "no adult-like clock organisation detected" (limits). Unadjusted p < 0.05 that does not survive Holm:
   "suggestive, not established", neither a pass nor a negative.
7. *Feasibility (D1b):* adult material cut to the fetal design. Each fetal individual is replaced by one adult donor
   with the same compartments, who receives its age, sex and stratum, so that the adjustment costs the same degrees
   of freedom; pass = the registered permutation test at alpha 0.05 (500 permutations), without the Holm step and
   without the companions, whose cost is not calibrated.
**Data and design (fixed).** E-MTAB-6814 CPM matrix; counts recovered as CPM / smallest non-zero CPM of the library;
one row per gene symbol; 14 clock genes, GTEx waveform. Organs: forebrain, hindbrain, heart, kidney, liver, ovary,
testis. The matrix does not name its libraries: each prenatal column was assigned to a library of Supplementary
Table 2, hence to an individual, from non-clock genes (`prep_fetal_bulk/scripts/link_cardoso_individuals.py`; 152 of
207 columns linked; leave-one-individual-out: 0 wrong of 1,183 decisions, 0.047 undecided). One linked liver library
has the other sex than its individual (the source disagrees with itself) and is treated as unlinked. L1 uses every
column (organ, age and sex are exact); L2 uses the linked ones. Covariates, per organ: polynomial of age in weeks
(degree 3 where the range exceeds 8 weeks, else 2) and sex. L2 strata: 4-8 | 9-12 | 13-20 wpc. Seed 237.
- **H237.1 (primary, L2):** a clock state shared across the organs of one fetus, `organs:9-20wpc`: 26 individuals
  with >= 2 linked organs (10 with 2, 6 with 3, 7 with 4, 3 with 5), 81 libraries, 7 organs (forebrain 14, heart 17,
  hindbrain 12, kidney 8, liver 14, ovary 5, testis 11 individuals); strata 9-12 wpc 15, 13-20 wpc 11. The range is
  the one where the adult pass rate reaches 0.8 (below); it excludes the first trimester's first half, where no clock
  is expected.
- **H237.2 (primary, L1, Holm over 6):** signed coherence across individuals at 4-20 wpc in heart (34), liver (37),
  forebrain (32), hindbrain (33), kidney (27) and ovary (18).
- **H237.3 (secondary, L1, Holm over 7; onset):** the same before and after the first trimester where the adult pass
  rate is >= 0.8: 4-8 wpc in heart (14), liver (16), hindbrain (16), kidney (9); 13-20 wpc in heart (9), liver (10),
  kidney (10).
- **Positive control (scored with the hypotheses, read first):** postnatal libraries of the same organs (neonate to
  senior, one individual each), L1 with a quadratic of the stage rank and sex: forebrain (21 individuals), hindbrain (25), heart (10),
  kidney (9) and liver (12), each a family of its own (adult pass rates 0.880-1.000, below); testis (13; 0.465) is
  descriptive. Same criteria (p < 0.05
  and clock-specific). The postnatal individuals died at unrecorded times.
**Power first (done; GTEx donors who died at 08-18 h, terminations being done in hospital hours; 200 subsamples;
`scouts_2026-11-10/feasibility_fetal/output/passrates_cm_day_sex.csv`).** Each
fetal individual is stood for by one GTEx donor with its organs (forebrain: frontal cortex BA9; hindbrain: cerebellum;
heart: left ventricle; kidney: cortex). Adjusted pass rates at alpha 0.05: H237.1 `organs:9-20wpc` 0.825 (95 % CI
0.765-0.875: the interval reaches below 0.8). H237.2: heart, liver, hindbrain and kidney 1.000, forebrain 0.965,
ovary 0.940 (at the Holm threshold of the worst rank, 0.05 / 6: 0.995 for the first four, forebrain 0.680, ovary
0.850). H237.3: 4-8 wpc heart 0.990, liver 0.985, hindbrain 0.985, kidney 0.875; 13-20 wpc heart 0.870, liver 0.960,
kidney 0.870 (worst rank, 0.05 / 7: 0.850, 0.860, 0.850, 0.670; 0.580, 0.725, 0.660). Controls: forebrain 0.880,
hindbrain 1.000, heart 0.955, kidney 0.895, liver 0.975. Below 0.8, hence descriptive: testis (prenatal 0.215-0.760,
postnatal 0.465), forebrain at 4-8 and 13-20 wpc (0.640 each), ovary at 4-8 (0.680) and 9-20 wpc (0.775), L2 at 4-20
(0.630), 13-20 (0.165) and 4-8 wpc (0.140). With the unsigned R (first run, no sex covariate,
`passrates_cm_day.csv`) the 4-20 wpc rates were 0.055 in forebrain, 0.495 in kidney and 0.730 in heart, against
0.970, 1.000 and 1.000 for C at the same designs (`passrates_cm_day_L1.csv`).
**Descriptive (no pass / fail):** testis at every range; every L1 window not named above (9-20 wpc; forebrain and
ovary windows); L2 at 4-20, 13-20 and 4-8 wpc (adult pass rates 0.630, 0.165 and 0.140: uninformative if negative); `organs_all` (unlinked
columns at their guess: an upper bound on what linking could add); the RNA quality number as a further covariate
(`:rqn`); unsigned R; leave-one-gene-out range of C; share of pseudo sets with p < 0.05.
**Reading (A237-A239 together).** Primary hypotheses: H237.1, H237.2, H238.1, H238.2, H239.1, H239.2 (14 primary
readouts). No correction is applied across hypotheses, so one isolated pass is read as A189 was: suggestive, one
cohort. The developmental claim of the paper is strengthened only by (i) H237.1, a clock state shared across the
organs of a fetus, or (ii) the same organ system passing in two datasets (brain: H237.2 forebrain or hindbrain with
H239.1 or H239.2; gonad: H237.2 ovary with H238.1 or H238.2), or (iii) an L1 pass that survives Holm with both
companions in an organ whose postnatal control passes. H237.3: late passes and early fails at an informative design:
compatible with an onset during the second trimester, not a test of it (the difference is not tested). All primary
readouts fail where controls pass: no adult-like clock organisation is detectable in human fetal organs at these
designs; the adrenal result of A189 stays a single observation and the paper says so. Any of these outcomes goes to
the PI before the paper plan changes.
**Limits.** (1) Sampling time is not recorded and terminations are done in hospital hours: a positive needs clock
phases that differ between individuals, as they do in adults sampled in such a window (the calibration); a fetal
clock entrained tightly to maternal time, or of low amplitude, would be missed and a negative cannot tell these from
no clock. (2) The companions' cost in sensitivity is not calibrated. (3) Bulk tissue: composition changes with age;
the adjustment removes its smooth part only. (4) L2 rests on the linking; a wrong link can only weaken it. (5)
Feasibility uses adult organs at adult amplitude and GTEx read depth. (6) Organs of one fetus share handling and RNA
quality: a gate passing its permutations but not the companions is "not clock-specific".
**Stop rule:** one scoring run per amendment; organs, ranges, covariates, strata, statistics and families fixed here
and in `FAMILIES` of each script; any further cut is a new amendment. **Code:**
`paper/exploratory/A237_fetal_organs/scripts/score_a237.py`, `paper/common/development.py`.

## A238 (2026-11-10, same agreement and same state of knowledge as A237; written before any clock statistic is computed on Garcia-Alonso 2022; seen: the pseudobulk sheet (donor, class, cells, depth, age, sex, assay), the detection of the clock genes per class (`prep_fetal_gonad/output/clock_detection_*.csv`) and the adult feasibility below): do the cell classes of the human fetal gonad carry a clock, alone and shared within one gonad? **Garcia-Alonso 2022, scRNA 10x 5', 6-21 wpf, 50 donors; secondary (development section)**
**Why.** The largest open human fetal single-cell cohort by donors: three to seven times the n of the adrenal positive
of A189, with stromal and endothelial classes (where the host is carried in adults) and several classes per donor,
which Cao 2020 did not have (A181). Its age range reaches the second trimester, where A189 found the adrenal clock.
**Data and design (fixed; common design of A237).** CELLxGENE collection 661a402a, female and male gonad datasets;
donor x class pseudobulks with >= 200 cells of the four classes present in both sexes: mesenchymal-stromal,
supporting (granulosa / Sertoli lineage), surface epithelial, endothelial (`prep_fetal_gonad/scripts/prep_gonad.py`;
germ cells, immune and other classes are not used: too few donors or one sex only). 14 clock genes, GTEx waveform.
Covariates, per class: quadratic age trend of each sex, sex, assay (5' v1 / v2); 6 columns in the pooled tests. The
per-sex trend was chosen after the first feasibility run (one common trend), on the argument that testis and ovary
diverge from 6-8 weeks, before any fetal statistic; the feasibility was rerun with it. L2 strata: sex x (before | from
13 wpf). Seed 238.
- **H238.1 (primary, L1, Holm over 4):** signed coherence across donors, pooled sexes, in mesenchymal-stromal (47
  donors), supporting (46), epithelial (33; 28 female, 5 male) and endothelial (24).
- **H238.2 (primary, L2):** shared clock state across the four classes of one gonad, `classes4:all`: 47 donors with
  >= 2 of the four classes (11 with 2, 19 with 3, 17 with 4; 147 pseudobulks; per class 45 / 46 / 33 / 23 donors;
  strata female early 16, female late 13, male early 7, male late 11).
**Power first (done, `scouts_2026-11-10/feasibility_fetal/output/passrates_gonad_sextrend.csv`, 200 subsamples).** No
adult gonad atlas is in the project: adult breast nuclei (10x 3' v3 stratum) stand in, fibroblasts for the stroma,
epithelium for the supporting cells, endothelium for endothelium and mural cells for the surface epithelium (which
has no adult counterpart; this fixes the number of compartments and their depth, not the biology). Each fetal donor
is replaced by an adult donor with its classes, thinned to the fetal depths, with its age, sex and assay labels. Pass
rates at alpha 0.05, adjusted as registered: L1 stromal 1.000, supporting 0.995, epithelial 0.985, endothelial 0.930
(at the Holm threshold of the worst rank, 0.05 / 4: 1.000, 0.895, 0.815, 0.840); L2 `classes4:all` 1.000. With heart
nuclei (pool large enough for two L1 jobs only): 1.000 and 1.000. The companions are not part of the calibration.
**Descriptive (no pass / fail):** every test per sex and from 13 wpf (the first run, with one age trend, gave adult
pass rates of 0.650-0.745 for endothelium per sex or from 13 wpf, 0.740 for epithelium from 13 wpf and 0.660 for
`classes4` in males: uninformative if negative; the others 0.870-1.000, `passrates_gonad.csv`); `classes3` (without
the surface epithelium) and the stromal-supporting pair; unsigned R; leave-one-gene-out range of C. The 2024 gonad
and reproductive-tract dataset (11 gonad donors) is not used.
**Reading.** As A237. A pass of H238.1 in a class is a clock-shaped, clock-specific order of the clock genes across
gonads; a pass of H238.2 is a state shared by the cell classes of a gonad beyond what matched and cue-responsive
genes share. If the L2 gate passes its permutations but not the companions, the shared state is not shown to be the
clock (ambient RNA, dissociation and handling are shared by the classes of one library): reported as "not
clock-specific". A negative at these pass rates reads "no adult-like clock organisation detected in the fetal gonad",
with the limits below.
**Limits.** Sampling times are not recorded; a negative is informative only if fetal clock phases differ between
donors as adult ones do (limits of A237). The adult stand-ins are other tissues and nuclei, not cells; the classes of
one donor come from the same libraries. Two assay versions, confounded with nothing by design but unequal (38 / 12
donors). **Stop rule:** one scoring run; classes, covariates, strata and families fixed here (`FAMILIES` in the
script). **Code:** `paper/exploratory/A238_fetal_gonad/scripts/score_a238.py`.

## A239 (2026-11-10, same agreement and same state of knowledge as A237; written before any clock statistic is computed on Kang 2011, prenatal or postnatal; seen: the sample sheet (brain, region, stage, age, sex, PMI, pH, RIN), the presence and median intensity of the clock genes (`prep_fetal_bulk/output/clock_detection_kang.csv`) and the adult feasibility below): do the regions of the human brain share a clock state before birth as they do after? **Kang 2011, exon arrays, 16 regions per brain, 20 prenatal and 32 postnatal brains; secondary (development section)**
**Why.** The only open dataset with many compartments per individual on both sides of birth, on one platform: the
postnatal brains are a positive control for the design, the arrays and the statistic, which neither A189 nor A237 /
A238 have for L2. The gate passes across the brain regions of adult donors cut to this design (feasibility below):
the prenatal test asks whether such a shared state exists before birth. It also reaches the third trimester and the
first postnatal year.
**Data and design (fixed; common design of A237).** GSE25219 gene-level table (log2 intensities), which has
no PER3 and no NR1D1: 12 clock genes, GTEx waveform without the two. The 16 regions sampled at every stage from
stage 3 on (A1C, AMY, CBC, DFC, HIP, IPC, ITC, M1C, MD, MFC, OFC, S1C, STC, STR, V1C, VFC); samples of one brain and
region averaged (hemispheres); brains with >= 6 of the 16 regions; 2^log2 intensity stands for a count at a constant
depth. Covariates, per region: polynomial of log age (degree 2 for the prenatal and the postnatal group, 1 for the
sub-groups) and sex. L2 strata: stages 1-3 | 4-6 | 7 | 8-9 | 10-11 | 12-15. L1 is computed on the mean over the
regions of a brain (one value per brain). Seed 239.
- **Positive control (scored with the hypotheses, read first):** postnatal brains (32; stages 8-15, 4 months to 82
  years; strata 6 / 7 / 19 brains; 481 brain x region means, 26-32 brains per region): L2 `regions:postnatal` and L1
  `brain_mean:postnatal`, same criteria as the hypotheses (p < 0.05 and clock-specific).
- **H239.1 (primary, L2):** shared clock state across the regions of one prenatal brain, `regions:prenatal`: 20
  brains (stages 2-7, 8-37 pcw; 11 female, 9 male; strata 4 / 14 / 2 brains), 285 brain x region means, 14-20 brains
  per region.
- **H239.2 (primary, L1):** signed coherence of the brain means across the 20 prenatal brains, `brain_mean:prenatal`.
**Power first (done, `scouts_2026-11-10/feasibility_fetal/output/passrates_kang_all.csv`, 100 subsamples).** GTEx
brain donors with >= 6 of 13 regions (290 donors, deaths at any hour, the 12 genes), cut to the number of brains of
each group with the ages and strata of the Kang brains (sex was not among the covariates of the calibration: one
parameter more per region in the registered analysis). Pass rates at alpha 0.05: `regions:prenatal` 0.96,
`brain_mean:prenatal` 1.00, `regions:postnatal` 0.98, `brain_mean:postnatal` 1.00; descriptive groups 0.89-0.98.
This calibration is RNA-seq: array noise, RNA quality and the companions are not in it, which is why the postnatal
control decides whether a prenatal negative can be read.
**Descriptive (no pass / fail):** trajectory of both readouts over 13-22 pcw (14 brains), infancy and childhood (13)
and from 13 years (19); both readouts with RIN added to the covariates (prenatal and postnatal); unsigned R;
leave-one-gene-out range of C.
**Reading.** Control passes and H239.1 passes: the regions of the human brain share a clock-specific state before
birth. Control passes, H239.1 and H239.2 fail: no adult-like shared clock state in the prenatal brain at a design
where the postnatal one is detected (informative, with the limits below); the descriptive trajectory then says where
it appears, without a test. Control fails: the arrays or the design do not carry the readout; the prenatal results
are uninformative whatever they are. A gate passing its permutations but not the companions is "not clock-specific"
(post-mortem interval, pH and RNA quality are shared by the regions of a brain).
**Limits.** Post-mortem brains with unrecorded time of death at every age; prenatal and postnatal brains differ in
cause of death and tissue handling; 4 of the 20 prenatal brains are at 8-12 pcw and 2 at 35-37 pcw, so the prenatal
group is essentially 13-22 pcw; 12 genes instead of 14; cell composition changes with age within a region (age
adjustment only). **Stop rule:** one scoring run; regions, groups, covariates, strata and families fixed here
(`FAMILIES` in the script). **Code:** `paper/exploratory/A239_fetal_brain/scripts/score_a239.py`.

## A240 (2026-11-11, user: "it would be very unlikely for human foetus that the collection time are very synchronized ... Register the two others"; a calibration on adult material, post hoc with respect to A237-A239, whose results are all seen; written before any adult pass rate at a reduced amplitude or in a narrower window is computed, except the code check of the script (3 subsamples per readout at scales 1, 0.5 and 0 and in both windows, `output/*_check.csv`: no pass at scale 0 in any readout)): how small, or how synchronised, would a fetal clock have to be to stay undetected in A237-A239? **Amplitude and collection-window bound; adult material at the fetal designs; secondary (development section)**
**Why.** The negatives of A237-A239 read "no adult-like clock organisation detected" because the feasibility used
adult clocks at adult amplitude, with adult donors who died at 08-18 h. Two alternatives were left open (limits of
A237): a fetal clock of lower amplitude, and fetal clock phases closer to each other than adult ones in such a window.
This amendment turns each informative negative into a number: the largest clock, relative to the adult one, that the
design would have missed.
**Data and design (fixed).** No fetal expression value is read. The registered feasibility runs are repeated on the
same subsamples and seeds (`scouts_2026-11-10/feasibility_fetal`: `passrates_cm_day_sex.csv`,
`passrates_gonad_sextrend.csv`, `passrates_kang_all.csv`), with the adult clock scaled: in every tissue or cell class
of the adult pool, the clock-plane state of each sample, centred over the pool, is multiplied by a
(`clocksync.simulate.scale_clock`, added for this purpose, clocksync 721115e; counts by binomial thinning, arrays in
log2 intensity). Grid a = 1, 0.75, 0.5, 0.35, 0.25, 0.15, 0. Pass = the registered permutation test at alpha 0.05 on
the adjusted data (500 permutations), without Holm and without companions, as in the feasibility. 200 subsamples per
readout and scale (Kang RNA-seq arm: the first 100 are the registered ones). Seed 240 for the thinning.
- *Readouts (21):* H237.1 `organs:9-20wpc`; the six H237.2 and the seven H237.3 L1 readouts; the four H238.1 classes
  and H238.2 `classes4:all` (adult breast nuclei, thinned to the fetal depths); H239.2 `brain_mean:prenatal` (GTEx
  brains, 12 genes). H239.1 is left out (its postnatal control was not clock-specific).
- *Array arm:* the 32 postnatal brain means of Kang 2011 (the passed L1 control of A239), 20 per subsample (the
  prenatal n), registered adjustment, clock scaled in log2 intensity: the bound on the platform of the prenatal
  negative.
- *Window arm (the user's question):* the 13 Cardoso-Moreira L1 readouts at a = 1 with the adult pool restricted to
  GTEx deaths at 08-14 h (177 donors; 12 readouts with a pool >= 1.5 x the fetal n: all but kidney 4-20 wpc) and at
  08-12 h (113 donors; 10 readouts: also without liver and hindbrain 4-20 wpc), against the registered 08-18 h (294
  donors).
**Readouts of the amendment (no pass / fail; a calibration).**
- **a80** per readout: the smallest a of the grid with a pass rate >= 0.8 at it and at every larger a. Statement
  licensed for a negative readout whose control passed: "an adult-organised clock of at least a80 times the adult
  between-individual clock variation would have been detected in 8 of 10 such cohorts". Reported with the pass-rate
  curve and its exact 95 % CI; a50 beside it.
- **Window:** pass rate per readout at 08-14 h and 08-12 h against 08-18 h.
- Descriptive: `below_obs`, the share of adult subsamples whose statistic is at or below the fetal C already scored
  (A237-A239 `l1.csv`), per scale, and a_excl, the smallest a with `below_obs` <= 0.05 at it and at every larger a;
  `kept`, the share of the clock-gene counts left after thinning.
**Checks (stop the reading of a readout if they fail).** At a = 1 every subsample reproduces the registered pass /
fail (the script stops otherwise). At a = 0 no clock-plane variation is left beyond counting noise, so the pass rate
must be at or below the nominal 0.05 (the plane is emptier than chance); if it exceeds 0.05 the scaling has failed for
that readout and its bound is not reported.
**Reading.** a is a bound on the product (amplitude) x (spread of phases between individuals), relative to adults
dying at 08-18 h (breast donors: unrecorded times; GTEx brains: any hour): a clock at full amplitude with phases half
as spread, and a clock of half the amplitude with adult spread, are the same point of the curve. The window arm
separates the two partly: if the pass rates hold at 08-12 h, a collection restricted to a morning theatre list does
not explain the negatives. The suggestive readouts (liver 4-20 wpc, hindbrain 4-8 wpc, gonad endothelium) get curves
but no bound claim: they are not negatives. No readout of this amendment makes a fetal negative a positive or the
reverse; the outcome changes one sentence of the limits ("a fetal clock below a80 of the adult one is not excluded").
**Limits.** (1) Conservative by construction, by an unknown margin: the scaling removes the non-clock variation of
the clock plane together with the clock (a real clock of relative size a would sit on the full in-plane noise, which
the coherence counts in its favour), and thinning removes counts (`kept`). The pass rate at a is a lower bound of the
detection rate of a clock of relative size a, and a80 an upper bound of the smallest detectable clock; `below_obs`
is biased the same way. (2) Adult tissues, stand-ins and depths; fetal
waveform assumed to be the GTEx one. (3) The companions are not calibrated. (4) GTEx time of death is not a
collection time, and a 4-6 h window of deaths still holds chronotype spread; the true fetal collection protocols
were not checked. (5) Array subsamples overlap (20 of 32 brains): their CI is not valid. (6) a_excl uses the fetal
statistics already seen and is descriptive.
**Stop rule:** one run; grid, readouts, windows and criteria fixed here. **Code:**
`paper/exploratory/A240_fetal_amplitude/scripts/amplitude_bound.py`.

## A241 (2026-11-11, user: "Register the two others"; written before any clock statistic of an organ x adrenal pair is computed in Cao 2020; seen: A189 (the adrenal is the one L1 positive: this is why it is chosen as host, post hoc), A190, A181 and the archive A107b free-host gates on the same fetuses (not estimable or failed), A237-A239, the sample sheet, the adult calibration below, the adrenal marker fractions of the other organs (median 0.001-0.016 of the adrenal's, maximum 0.077) and a dry run of the script on one expression-matched pseudo gene set): do the organs of a human fetus share the clock state of its adrenal? **Cao 2020, sci-RNA-seq3, the 7 fetuses with an adrenal (89-122 days); secondary (development section)**
**Why.** The fetal adrenal is the only compartment with a clock-specific order of the clock genes (A189). If the
adrenal clock runs and other organs follow it (or the same maternal cue), the clock state of an organ should covary
with that of the adrenal of the same fetus even where the organ alone shows nothing at L1 (3-6 fetuses). A181 asked
this with a free host and could not estimate it; here the host is fixed to the compartment known to carry a clock.
**Data and design (fixed).** Genome-wide whole-organ pseudobulks (`natexp/cao2020/genome/pseudobulk_organ.h5ad`), >=
1,000 cells, fetuses with an adrenal pseudobulk (7). 14 clock genes, GTEx waveform. Host-restricted gate,
`Dataset.gate(targets=organs, hosts=["Adrenal"], min_shared=3)`: only organ x adrenal cross-covariances; pairing null
= the fetuses of each organ permuted; 10,000 permutations; p = the larger of p_pairing and p_acrophase. No age
adjustment (3-6 fetuses per organ; an age-driven shared programme is caught by the companions). Seed 241.
- **H241.1 (L2, the only hypothesis):** `organs6`: lung (6 fetuses with an adrenal), heart (4), liver (4), cerebellum
  (3), cerebrum (3), muscle (3): 23 organ pseudobulks, 6 organ x adrenal pairs. The organs are those with an adult
  stand-in in the calibration; `min_shared` 3 is the variant with the highest adult pass rate at the fetal depths
  (below).
- *Companions (both required):* 200 expression-matched and 200 responsive pseudo-clock sets, same gate (1,000
  permutations, same seed), compared on `pairing_z`; clock-specific: above the q95 of both.
- *Adrenal RNA (required; A190 lesson):* per fetus and organ, the fraction a of adrenal cortex marker counts (STAR,
  CYP11A1, CYP17A1, CYP21A2, CYP11B1, NR5A1) relative to the adrenal of the same fetus (an upper bound); the expected
  adrenal-derived clock-gene counts are subtracted and the gate is rerun (`organs6:corrected`).
- **PASS:** p < 0.05, clock-specific, and p < 0.05 on the corrected counts. p < 0.05 without clock-specificity: "not
  clock-specific". p < 0.05 lost after correction: "not shown (adrenal RNA)".
**Power first (done, `output/passrates.csv`; 200 subsamples, 500 permutations).** Each of the 7 fetuses is stood for
by one GTEx donor (death at 08-18 h) with an adrenal and the stand-in organs of that fetus (frontal cortex BA9 for
cerebrum, left ventricle, skeletal muscle), thinned to the fetal depths. Adult pass rates: `min_shared` 3: 0.140 (95
% CI 0.095-0.196); 4: 0.125; 5 (lung only): 0.060; unthinned 0.130, 0.135, 0.065. p_acrophase < 0.05 in 0.98 of the
subsamples and p_pairing < 0.05 in 0.14 (`min_shared` 3): the shape is seen, the pairing of 7 individuals is not.
**The design is not powered: an adult clock shared with the adrenal would be missed in 6 of 7 such cohorts.** It is
registered and scored once all the same, at the user's request, with the reading fixed in advance.
**Descriptive (no pass / fail):** `organs6+placenta` (placenta 5 fetuses; marker correction not valid there: the
placenta expresses CYP11A1), with companions; one organ at a time (`only:<organ>`); the lung on log fractions
residualised on a linear age trend within organ (`only:Lung:age`); the marker fractions. Phase offsets of the organs
against the adrenal are estimated only in a later amendment, and only if H241.1 passes.
**Reading.** Fail (expected in 86 % of adult-like cohorts): uninformative, never "the organs do not follow the
adrenal"; the question stays open and needs more fetuses with an adrenal. PASS: suggestive evidence that organs of a
mid-gestation fetus share a clock-specific state with its adrenal; one cohort, 7 fetuses, a host chosen after A189:
to the PI before any use, and a replication cohort becomes necessary. Not clock-specific or lost after correction: a
shared non-clock state or adrenal RNA, reported as such.
**Limits.** 7 fetuses; organs of one fetus share dissociation, batch and RNA quality (companions); the host was
chosen after seeing A189; the same 7 adrenals enter every pair, so the pairs are not independent (the permutation
null accounts for it); adult stand-ins at 08-18 h; no age adjustment in the primary.
**Stop rule:** one scoring run; organs, `min_shared`, companions and correction fixed here (`TESTS`, `FAMILIES` in
the script). **Code:** `paper/exploratory/A241_fetal_adrenal_host/scripts/{feasibility_a241.py,score_a241.py}`.

## A242 (2026-11-11, user: "Yes go ahead", to "register the open-data adrenal replication as A242, starting with that power calculation"; written before any clock statistic is computed on these fetal data; seen: the sample sheets and the design tables (`output/design_cohorts.csv`, `design.csv`), the median count of each clock gene per cohort (`output/clock_detection.csv`), the adult feasibility below, a dry run of the companion gene sets (gene lists only, no statistic), every result of A189, A190 and A237-A241; the surrogate components of the Schroder libraries are computed by the script and were not inspected): does the clock of the human fetal adrenal replicate in open data from other laboratories? **Nine open cohorts, 59 fetal adrenals; L1; secondary (development section)**
**Why.** The fetal adrenal of Cao 2020 is the one human prenatal positive (A189: 7 fetuses, 89-122 days; signed
coherence C 0.926, p 0.0001, post hoc). Three other datasets are negative in other organs (A237-A239) and exclude a
clock of 0.75-1 times the adult one (A240). One organ in one cohort cannot carry a section (Q32); the scout found
open fetal adrenal data from eight other laboratories and ENCODE, none of which shares a fetus with Cao 2020.
**Data and design (fixed; `prep_adrenal.py`, `score_a242.py`).** One whole-adrenal pseudobulk per fetus and cohort:
all cells of the fetus summed (the organ-level pseudobulk of A189), bulk libraries as released. Samples of one source
with the same stated age are one fetus unless the source names them as different individuals. 14 clock genes, GTEx
waveform, signed coherence C (`common.development`, design points 1, 2, 4 and 5 of A237). Seed 242.
- *Cohorts (fetuses; age as stated by the source, in weeks):* `encode_pa` (ENCODE polyA RNA-seq, 4; 12.1-15.4, post
  conception), `kameneva` (GSE147821, 10x, 11; 6-14, post conception), `dong` (GSE137804, 10x, 4; 8.4-14.6), `wang`
  (GSE167860, 10x, 10; 7-14, gestational), `han` (GSE134355, Microwell-seq, 2; 12 and 14), `sasaki` (GSE193672 and
  GSE291218, 10x, 7; 5-19, post fertilisation; three samples deposited in both series are counted once), `roost`
  (GSE66302, DeepSAGE, 3; 9.4-22), `kildisiute` (Kildisiute 2021 processed object, 10x, 5; 8-21), `schroder`
  (GSE196700, targeted transcriptome, 13 bulk adrenals of the first and second trimester, no age per sample; three
  sequencing batches of 4, 7 and 2).
- *Pooling.* Log fractions are residualised gene by gene on a cohort indicator and a common polynomial of age; only
  cohorts with >= 2 fetuses in a test enter it (a single fetus has no residual). The 10,000 per-gene permutations of
  the null are done within cohort (`development.coherence(strata=)`, added for this amendment): cohorts differ in
  platform, depth and noise, so their fetuses are not exchangeable. A deviation from A237-A239 (one cohort each),
  fixed before scoring. Companions: 200 expression-matched and 200 responsive pseudo-clock sets among the genes with
  a mean count >= 1 and present in every cohort of the test, same adjustment; clock-specific: C above the q95 of both.
- **H242.1 (primary; the replication; Holm over its 2 readouts):**
  `late`: fetuses from 12 weeks, the age range of A189 (12.7-17.4 weeks): 19 fetuses in 7 cohorts (encode_pa 4,
  kameneva 3, wang 3, sasaki 3, dong 2, han 2, roost 2); cohort indicator and a linear age trend (7 covariates).
  `schroder`: the 13 bulk adrenals; batch indicator and the first principal component of the log fractions of the
  non-clock genes (mean count >= 50, batch regressed out; `development.surrogates`) in place of the unrecorded age;
  permutations within batch.
- **H242.2 (secondary; onset):** `early`: fetuses before 12 weeks: 25 fetuses in 5 cohorts (kameneva 8, wang 7,
  sasaki 4, kildisiute 4, dong 2); cohort indicator and a quadratic age trend (6 covariates).
- **PASS** (`development.verdicts`): Holm p < 0.05 within the hypothesis and clock-specific. p < 0.05 without
  clock-specificity: "not clock-specific". Unadjusted p < 0.05 lost at Holm: "suggestive, not established".
**Power first (done; `feasibility_a242.py`, 538 s; `output/feasibility_passrates.csv`; 200 subsamples, 500
permutations within cohort).** Each fetus is stood for by one distinct GTEx adrenal gland (93 donors who died at
08-18 h) that receives its cohort, age and depth; the adult clock-gene counts are thinned gene by gene to what the
cohort captures (cohort median count of the gene x relative depth of the fetus; where the adult library holds fewer
reads the targets of the gene are lowered over the cohort, which the cohort indicator absorbs: 35 % of the fetus x
gene targets of `late`). Adult pass rates at alpha 0.05: `late` 0.995 (95 % CI 0.972-1.000), `schroder` 0.990
(0.964-0.999), `early` 1.000 (0.982-1.000); at the Holm threshold of the worse rank (0.025): `late` 0.985, `schroder`
0.965. With the adult clock scaled (`simulate.scale_clock`) to 0.75: 0.930, 0.895, 1.000; to 0.5: 0.375, 0.410,
0.560. False-positive rate (every clock gene permuted across the drawn donors on its own): 0.040 (0.017-0.077),
0.075 (0.043-0.121), 0.065 (0.035-0.109): the nominal 0.05 lies in each interval. The other designs: `all` 1.000
(0.75: 1.000; 0.5: 0.670), `late:deep` 0.985, `kameneva` 0.995, `wang` 0.985, `schroder:pc0` 1.000, `schroder:pc2`
0.970. `late` is kept with its two shallow cohorts (han, roost: median clock-gene count < 100): leaving them out
does not raise the adult pass rate.
**Descriptive (no pass / fail):** `all` (every age, 46 fetuses in 8 cohorts, quadratic age trend); `late:deep`;
`late` without each cohort in turn (`late:-<cohort>`); `kameneva` (11) and `wang` (10) alone; `schroder:pc0` and
`schroder:pc2`; unsigned R; leave-one-gene-out range of C; the position of the fetal C in the adult subsamples of
the feasibility (`feasibility_subsamples.csv`).
**Reading.** `late` PASS: the fetal adrenal clock is replicated in fetuses from other laboratories and platforms, at
the age of A189; the adrenal becomes a two-dataset positive, the exception to the prenatal negatives is real, and
the development section can state it (to the PI: Q32). `schroder` PASS alone: suggestive support on a third platform
with unknown ages, not the replication. `late` negative (unadjusted p >= 0.05; informative: adult pass rate 0.995,
0.930 at 0.75 of the adult clock): A189 is not replicated; an adrenal clock of at least three quarters of the adult
one, with phases that differ between the fetuses of a cohort, is excluded at 12-22 weeks in these cohorts, and A189
stays a single-cohort observation that the paper must present as unreplicated (or drop: PI). Not clock-specific: a
state shared with expression-matched or glucocorticoid-responsive genes (the adrenal makes the glucocorticoid), not
shown to be the clock. H242.2: `early` negative with `late` positive is compatible with an onset near the end of
the first trimester, not a test of it; `early` positive: the adrenal clock is present before 12 weeks.
**Limits.** (1) Sampling times are not recorded. The cohort indicator removes what the fetuses of a cohort share: a
laboratory that collects at one hour of the day contributes nothing, and only phase differences within a cohort are
read. (2) Age conventions differ between sources (post conception, post fertilisation, gestational, not stated):
up to two weeks of ambiguity at the 12-week boundary and in the common age trend. (3) Fetus identity rests on the
sources' sample sheets (same-age samples merged). (4) Sex is not adjusted (not given by every source; 11 residual
degrees of freedom in `late`). (5) Whole-adrenal pseudobulks: cortex, medulla and blood proportions vary between
fetuses; the age trend removes their smooth part and the companions guard the rest. (6) Eight of the nine cohorts
count 3' tags or amplicons, one (ENCODE) full-length reads; Kildisiute's counts are recovered from ambient-corrected
log-normalised values. (7) Feasibility uses adult adrenals at adult amplitude; the surrogate component of the
adult libraries may carry time of death, which makes the `schroder` pass rate conservative. (8) The ENCODE fetuses
come from the same tissue programme as Cao 2020 (different fetuses and assay): `late:-encode_pa` shows what rests on
them. (9) Excluded before scoring: ENCODE long-read libraries (1.5-9 reads per clock gene, two clock genes absent),
single-fetus cohorts in a window, GSE64283 (no gene quantification), sources that release reads or alignments only.
**Stop rule:** one scoring run; cohorts, windows, covariates, permutation scheme, companions and families fixed here
and in `FAMILIES`; any further cut is a new amendment. **Code:**
`paper/exploratory/A242_fetal_adrenal_replication/scripts/{prep_adrenal.py,feasibility_a242.py,score_a242.py}`,
`paper/common/development.py`.

## A243 (2026-11-11, user: "okay, register 1 and 2 and write a prompt for a new agent to launch them"; staged: the hypotheses, statistic, eligibility rules, families and readings are fixed here, the cohort list and the adult power go in a design addendum committed before any clock statistic is computed on the new cohorts; seen: every result of A189 (Cao 2020 organs, signed C: heart 0.140, lung 0.116, liver 0.114, p >= 0.29), A237 (Cardoso-Moreira, 4-20 wpc: heart C -0.094, liver 0.387, kidney -0.350), A240 and A242; no other fetal data of these organs): is a clock of the size read in the fetal adrenal present in the heart, kidney, lung or liver before birth when open cohorts are pooled? **Four organs, pooled open cohorts; L1; secondary (development section)**
**Why.** The fetal adrenal carries a clock-shaped covariation between fetuses in ten cohorts (A189, A242), at
between half and three quarters of the adult between-individual clock variation. The negatives of the other organs
(A189, A237-A239) exclude a clock of 0.75 to 1 times the adult one only (A240): a clock of the adrenal's size would
have been missed, so "the adrenal starts first" is not shown (Q32). Pooling cohorts with a cohort indicator and
within-cohort permutations (A242) is the one design that reaches half the adult size (`all` of A242: adult pass rate
0.670 at 0.5 with 46 fetuses). This amendment applies it to four organs with open fetal data from several
laboratories.
**Fixed now.**
- *Organs:* heart, kidney, lung, liver. Other organs met by the scout are reported, not scored.
- *Cohort eligibility* (applied to metadata only): (E1) human embryonic or fetal tissue of the organ, not cultured,
  not an organoid, no stated anomaly of the organ; (E2) open gene-level counts (UMI, tag or read counts; expected
  counts of a quantifier accepted), bulk or single-cell (all cells of a fetus summed); no FASTQ mapping, no
  controlled access, no source that releases normalised values only unless the counts are recoverable exactly as for
  Kildisiute in A242; (E3) an age stated per sample; (E4) the 14 clock genes present and a median count per clock
  gene >= 20 (median over genes of the cohort medians); (E5) >= 2 fetuses in the readout; (E6) a fetus present in two
  sources is counted once, in the deeper one. Cao 2020 (organ pseudobulks with >= 1,000 cells, as A241) and
  Cardoso-Moreira 2019 (prenatal libraries with a resolved individual, as A237) are cohorts of the pool; both were
  scored alone before (negative).
- *Unit, statistic, pooling:* as A242. One whole-organ pseudobulk or bulk library per fetus and cohort (samples of
  one source with the same stated age are one fetus unless the source names different individuals; several parts of
  one organ are summed); 14 clock genes, GTEx waveform, signed coherence C; log fractions residualised gene by gene
  on a cohort indicator and a common age polynomial; 10,000 per-gene permutations within cohort; 200
  expression-matched and 200 responsive pseudo-clock sets among the genes with a mean count >= 1 present in every
  cohort of the readout; clock-specific: C above the q95 of both. Seed 243. `paper/common/development.py`.
- *Readouts per organ:* `<organ>:all` (every fetus; quadratic age), `<organ>:late` (from 12 weeks; linear age),
  `<organ>:early` (before 12 weeks; quadratic age): the windows and covariates of A242.
- *Power first* (as `feasibility_a242.py`): each fetus stood for by one GTEx donor of the adult organ (heart left
  ventricle, kidney cortex, lung, liver; deaths at 08-18 h; distinct donors within a cohort) thinned gene by gene to
  the clock-gene counts of its cohort; adult clock at 1, 0.75 and 0.5 (`simulate.scale_clock`) and without a clock
  (clock genes permuted across the drawn donors); 200 subsamples, 500 permutations within cohort, alpha 0.05. A
  readout is scored only if its adult pass rate at the full clock is >= 0.8 (D1b rule); otherwise it is "not
  estimable" and its fetal statistic is not computed.
- **H243.1 (primary):** `<organ>:all` for each organ scored; Holm over the organs scored.
- **H243.2 (secondary; onset):** `<organ>:late` and `<organ>:early` for each organ scored; Holm over the readouts
  scored.
- **PASS** (`development.verdicts`): Holm p < 0.05 within the hypothesis and clock-specific. p < 0.05 without
  clock-specificity: "not clock-specific". Unadjusted p < 0.05 lost at Holm: "suggestive, not established".
- *Descriptive (no pass / fail):* `<organ>:all:new` (without Cao 2020 and Cardoso-Moreira 2019); `<organ>:all`
  without each cohort in turn; each cohort with >= 10 fetuses alone; unsigned R; leave-one-gene-out range; the
  position of the fetal C among the adult subsamples at 1, 0.75 and 0.5 (as `position_a242.py`), beside that of the
  adrenal (`all` of A242).
**Design addendum (before scoring; `### A243 design addendum`, committed with the prep and feasibility scripts).**
The cohorts found and those excluded, each with the rule that excludes it; per cohort and window: accession, assay,
age convention, fetuses, depth, median clock-gene counts; the design table (fetuses, cohorts, covariates per
readout); the adult pass rates at 1, 0.75, 0.5 and without a clock; the readouts scored; a dry run of the companion
gene lists. Nothing else may be chosen there.
**Reading.** `<organ>:all` PASS: a clock-shaped covariation between fetuses is present in that organ before birth,
in pooled cohorts, as in the adrenal; the developmental result becomes "present in several organs, below the adult
size", and the earlier negatives are read as bounds on size (to the PI, Q32). Negative (unadjusted p >= 0.05) with an
adult pass rate >= 0.8 at 0.5: a clock of half the organ's adult size is excluded; the organ is behind the adrenal
(whose signal lies between 0.5 and 0.75 of its own adult size), and "the adrenal is ahead of <organ>" can be stated
for the ages covered. Negative with >= 0.8 at 0.75 only: three quarters excluded, the comparison with the adrenal
stays open. Negative below that: bound as in A240, nothing new. Not clock-specific: a state shared with matched or
responsive genes, not shown to be the clock. H243.2 orders nothing by itself: a window is read with its own adult
pass rate.
**Limits.** (1) to (4) of A242 (unrecorded times and the cohort indicator; age conventions; fetus identity; no sex
adjustment). (5) Cao 2020 and Cardoso-Moreira 2019 were scored before on their own: `all:new` shows what the new
cohorts give. (6) Sizes are relative to each organ's own adult clock variation in GTEx; the adult stand-in is one
part of the organ (left ventricle, kidney cortex). (7) The fetal liver is a haematopoietic organ at these ages:
its whole-organ pseudobulk is largely blood cells, which carry no coherent core clock in the adult; a liver result
is of that mixture. (8) Whole-organ pseudobulks: cell composition varies between fetuses; the age trend and the
companions are the control. (9) New data: at most 60 GB, within the 500 GB cap.
**Stop rule:** one addendum, one scoring run; any further cut is a new amendment. **Code:**
`paper/exploratory/A243_fetal_organs_pooled/scripts/` (prep, feasibility, score, position), on
`paper/common/development.py` and the scripts of A242.

### A243 design addendum (2026-11-16; before any clock statistic on the new cohorts; seen: the metadata, the clock-gene count medians of `prep_a243.py`, the adult feasibility runs and a dry run of the companion gene lists; no fetal clock statistic)
**Cohorts** (`output/design_cohorts.csv`; E4 = median over the 14 clock genes of the cohort medians of the count;
search and per-source verdicts in `scout/sources.md`):

| organ | cohort | accession | assay | age convention | fetuses (< 12 / >= 12 wk) | weeks | median depth | E4 |
|---|---|---|---|---|---|---|---|---|
| heart | cao2020 | GSE156793 | sci-RNA-seq3 | post conception (ga_pcw) | 9 (0 / 9) | 12.86-17.43 | 5.1 M | 151.5 |
| heart | cardoso_moreira | E-MTAB-6814 | bulk RNA-seq | post conception | 27 (18 / 9) | 4-19 | 23.8 M | 666.5 |
| heart | cui | GSE106118 | plate scRNA-seq (UMI) | not stated | 18 (8 / 10) | 5-25 | 18.1 M | 236 |
| heart | encode_cshl | ENCSR000AEZ | total RNA-seq | not stated (ENCODE) | 2 (0 / 2) | 19-28 | 23.2 M | 655.75 |
| heart | encode_uw | ENCSR047LLJ+ENCSR863BUL | polyA+ RNA-seq | not stated (ENCODE) | 2 (0 / 2) | 13-17.14 | 52.8 M | 1835.85 |
| heart | han_hcl | GSE134355 | Microwell-seq | not stated | 2 (1 / 1) | 11-12 | 5.1 M | 38.5 |
| heart | knight_schrijver | GSE216019 | 10x 3' v3 | not stated | 6 (6 / 0) | 8.29-11.43 | 107.3 M | 570.5 |
| heart | roost | GSE66302 | DeepSAGE tags | gestational | 6 (2 / 4) | 9-22 | 10.5 M | 119.25 |
| kidney | cao2020 | GSE156793 | sci-RNA-seq3 | post conception (ga_pcw) | 9 (1 / 8) | 10.29-17.86 | 9.1 M | 300 |
| kidney | cardoso_moreira | E-MTAB-6814 | bulk RNA-seq | post conception | 21 (12 / 9) | 4-20 | 21.7 M | 641.5 |
| kidney | encode_uw | ENCSR212AMA+ENCSR495UXA+ENCSR896QPD+ENCSR907KDH | polyA+ RNA-seq | not stated (ENCODE) | 4 (0 / 4) | 12.14-15.43 | 45.0 M | 1337.05 |
| kidney | han_hcl | GSE134355 | Microwell-seq | not stated | 4 (2 / 2) | 11-13 | 6.0 M | 34.5 |
| kidney | hochane | GSE114530 | 10x 3' v2 | not stated | 5 (2 / 3) | 9-18 | 20.9 M | 326.5 |
| kidney | lindstrom_combes | GSE102596+GSE124472+GSE139280 | 10x 3' | not stated | 4 (0 / 4) | 14-17 | 39.3 M | 270 |
| kidney | roost | GSE66302 | DeepSAGE tags | gestational | 6 (2 / 4) | 9.3-22 | 6.8 M | 53.5 |
| lung | cao2020 | GSE156793 | sci-RNA-seq3 | post conception (ga_pcw) | 11 (0 / 11) | 12.71-17.86 | 10.8 M | 275 |
| lung | encode_cshl | ENCSR000AFC | total RNA-seq | not stated (ENCODE) | 2 (0 / 2) | 20-24 | 34.1 M | 895.75 |
| lung | han_hcl | GSE134355 | Microwell-seq | not stated | 2 (1 / 1) | 11-12 | 6.0 M | 37 |
| lung | roost | GSE66302 | DeepSAGE tags | gestational | 6 (2 / 4) | 9.1-22 | 7.2 M | 71.5 |
| lung | sountoulidis | GSE215895 | 10x 3' | post conception | 17 (11 / 6) | 5-14 | 173.5 M | 3664.5 |
| liver | cao2020 | GSE156793 | sci-RNA-seq3 | post conception (ga_pcw) | 7 (0 / 7) | 12.71-17.14 | 6.8 M | 93 |
| liver | cardoso_moreira | E-MTAB-6814 | bulk RNA-seq | post conception | 26 (19 / 7) | 4-19 | 20.3 M | 230 |
| liver | encode_cshl | ENCSR000AFB | total RNA-seq | not stated (ENCODE) | 2 (0 / 2) | 20-22 | 41.0 M | 822.25 |
| liver | roost | GSE66302 | DeepSAGE tags | gestational | 6 (2 / 4) | 8.3-22 | 6.9 M | 19 (fails) |

**Excluded, with the rule** (all in `scout/sources.md`): Roost liver (E4: 19.0 < 20); He 2022 lung (E1: the Methods
state CD326/EPCAM MACS enrichment, the enriched libraries are not identifiable); Barnes 2023 lung (E1: enrichment not
verifiable); Renesme 2025 lung (E2: the UMI matrix has no gene names or ids); Quach 2024 lung (E2: an 8 GB Seurat
object only; reading it would execute downloaded code); Popescu 2019 kidney and liver, Stewart 2019 kidney (E1: sorted
fractions; Stewart and Popescu share fetuses, E6); Han HCL liver (E3/E5: one or two samples of 26-27 weeks, ambiguous
donor); sorted, part-of-organ, spatial, TPM/RPKM-only, long-read ENCODE and controlled-access sources (E1, E2, as
tabled). Differences from the scout's counts: Cardoso-Moreira counts only prenatal individuals resolved by the A237 link
tiers (heart 27, kidney 21, liver 26; the scout's 34, 27, 37 counted unresolved libraries); Cui heart 18 fetuses;
Lindstrom / Combes kidney 4. Fetuses of UK HDBR sources (Cardoso-Moreira, Knight-Schrijver) cannot be cross-matched on
metadata; ENCODE UW and Cao fetuses of the same stated age and sex are kept apart as in A242 (limit (3) of A242).

**Design and adult power** (`output/feasibility_design.csv`, `output/feasibility_passrates.csv`;
`feasibility_a243.py`, 200 subsamples, 500 permutations; fetuses of a window counted in cohorts with >= 2 fetuses in
that window, so `heart:early` and `heart:late` drop the one Han fetus each):

| readout | fetuses | cohorts | covariates | adult pass rate at 1 / 0.75 / 0.5 | without a clock |
|---|---|---|---|---|---|
| `heart:all` | 72 | 8 | 9 | 1.000 / 1.000 / 0.800 | 0.055 |
| `heart:late` | 36 | 6 | 6 | 1.000 / 0.995 / 0.605 | 0.065 |
| `heart:early` | 34 | 4 | 5 | 1.000 / 0.985 / 0.520 | 0.050 |
| `heart:all:new` | 36 | 6 | 7 | 1.000 / 0.985 / 0.435 | 0.040 |
| `kidney:all` | 53 | 7 | 8 | 1.000 / 0.995 / 0.250 | 0.045 |
| `kidney:late` | 34 | 7 | 7 | 0.995 / 0.945 / 0.155 | 0.015 |
| `kidney:early` | 18 | 4 | 5 | 0.985 / 0.855 / 0.160 | 0.045 |
| `kidney:all:new` | 23 | 5 | 6 | 0.985 / 0.785 / 0.155 | 0.055 |
| `lung:all` | 38 | 5 | 6 | 1.000 / 0.995 / 0.390 | 0.070 |
| `lung:late` | 23 | 4 | 4 | 1.000 / 0.945 / 0.300 | 0.045 |
| `lung:early` | 13 | 2 | 3 | 0.950 / 0.670 / 0.175 | 0.045 |
| `lung:all:new` | 27 | 4 | 5 | 0.995 / 0.940 / 0.285 | 0.080 |
| `liver:all` | 35 | 3 | 4 | 1.000 / 0.995 / 0.405 | 0.065 |
| `liver:late` | 16 | 3 | 3 | 1.000 / 0.840 / 0.225 | 0.065 |
| `liver:early` | 19 | 1 | 2 | 1.000 / 0.935 / 0.260 | 0.050 |
| `heart:cardoso_moreira` | 27 | 1 | 2 | 1.000 / 1.000 / 0.620 | 0.075 |
| `heart:cui` | 18 | 1 | 2 | 1.000 / 0.960 / 0.545 | 0.060 |
| `kidney:cardoso_moreira` | 21 | 1 | 2 | 1.000 / 0.960 / 0.160 | 0.055 |
| `liver:cardoso_moreira` | 26 | 1 | 2 | 1.000 / 0.995 / 0.285 | 0.045 |
| `lung:cao2020` | 11 | 1 | 2 | 0.915 / 0.705 / 0.185 | 0.040 |
| `lung:sountoulidis` | 17 | 1 | 2 | 0.995 / 0.895 / 0.285 | 0.045 |

**Readouts scored.** All registered readouts have an adult pass rate >= 0.8 at the full clock (lowest `lung:early`
0.950), so H243.1 has four organs (Holm over 4) and H243.2 eight windows (Holm over 8). `liver:all:new` does not exist
(without Cao and Cardoso-Moreira the liver keeps ENCODE CSHL's 2 fetuses; the readout needs >= 3). The descriptive
leave-one-cohort-out readouts (`<organ>:all:-<cohort>`) were not simulated; they are scored with the power of
`<organ>:all` as stated in `score_a243.tests` and read as descriptive only. At half the adult size only `heart:all`
reaches 0.8 (0.800): by the registered reading, a negative heart excludes a clock of half the adult size; a negative in
kidney, lung or liver excludes three quarters (pass rates at 0.75 >= 0.995 for every `all`), and the comparison with the
adrenal stays open there.

**Companion dry run** (`companions_dryrun_a243.py`, `output/companions_dryrun.csv`; no clock statistic): the scored pool
is that of `development.log_mean` as in A242, i.e. genes measured (not missing) in every cohort of the readout with a
mean count >= 1 over its fetuses. Pools of the primary readouts: heart 13,242 genes, kidney 13,402, lung 13,542, liver
17,749; of these, 230, 122, 192 and 2,086 have a zero total in one cohort (measured, not detected; mostly in the
shallow Cao and Han libraries). The first matched and responsive set of each readout is listed in the table (the
responsive sets draw from the immediate-early and glucocorticoid-responsive genes, e.g. EGR2, KLF2, FKBP5, NR4A2 for
`heart:all`). Nothing else is chosen here.

## A244 (2026-11-11, user: "okay, register 1 and 2 and write a prompt for a new agent to launch them"; staged: hypotheses, classes, marker lists, families and readings are fixed here, the per-cell classifier thresholds, the design table and the adult power go in a design addendum committed before any clock statistic is computed on a cell class; seen: every result of A242 on the whole adrenals of the same fetuses, A189, A190 and its post hoc checks; no per-cell data of these cohorts beyond the summed counts): is the clock-shaped covariation of the fetal adrenal carried by the cortex, by the medulla, or by both? **Six single-cell cohorts of A242, 39 fetuses; L1 per cell class; secondary (development section)**
**Why.** A242 reads the adrenal as one sample per fetus. A whole-organ signal can come from the steroidogenic
cortex (the adult adrenal clock), from the neural-crest-derived medulla, from both, or from the proportions of the
classes varying between fetuses (limit 5 of A242). Within Cao 2020 the classes could not be separated: the shallow
classes were about one third cortex RNA (A190, post hoc). The six single-cell cohorts of A242 allow the split with
a cohort indicator, and droplet data may carry less foreign RNA.
**Fixed now.**
- *Cohorts:* the single-cell cohorts of A242: `kameneva` (11 fetuses), `wang` (10), `sasaki` (7), `kildisiute` (5),
  `dong` (4), `han` (2); the three bulk cohorts cannot be split. Fetuses, ages and cells as in `prep_adrenal.py`.
- *Classes, by marker counts per cell; clock genes and pseudo-clock sets are never read by the classifier:*
  `cortex` (steroidogenic; the A190 list: STAR, CYP11A1, CYP17A1, CYP21A2, CYP11B1, NR5A1); `medulla` (neural
  crest-derived: chromaffin cells, sympathoblasts, Schwann cell precursors: CHGA, CHGB, TH, DBH, PHOX2B, STMN2,
  ELAVL4, SOX10, PLP1); `stroma` (descriptive: COL1A1, COL1A2, COL3A1, PDGFRB, PECAM1, CDH5); blood cells (HBG1,
  HBG2, HBA1, HBB, PTPRC) are left out. One per-cell rule on the marker fractions for every cohort; a cell that
  meets two classes or none is left out.
- *Unit:* one pseudobulk per fetus and class with >= 200 cells of the class (A190); a cohort enters a readout with
  >= 2 such fetuses.
- *Statistic, pooling:* as the `all` readout of A242: signed coherence C across fetuses; cohort indicator and a
  common quadratic age trend; 10,000 per-gene permutations within cohort; both companions (200 sets each) among the
  genes with a mean count >= 1 present in every cohort; clock-specific: C above the q95 of both. Seed 244.
- *RNA of the other class (the A190 limit):* per fetus, a = fraction of the other class's marker counts in the class
  pseudobulk, relative to that fraction in the other class's pseudobulk of the same fetus (an upper bound of the
  foreign share); the expected foreign clock-gene counts are subtracted (floored at 0, rounded; as
  `score_a241.contamination`) and the readout is rerun (`<class>:corrected`). Reported per cohort: a and the
  expected foreign share of the clock-gene counts. A class whose median foreign share of the clock-gene counts
  exceeds 0.5 is "confounded" and has no verdict.
- *Power first:* GTEx adrenal glands thinned to the clock-gene counts of each class and cohort
  (`feasibility_a242.py`), adult clock at 1, 0.75 and 0.5 and without a clock; 200 subsamples. A class is scored
  only if its adult pass rate at the full clock is >= 0.8; otherwise "not estimable", not computed.
- **H244.1 (primary; Holm over the classes scored):** `cortex`, `medulla`. **PASS:** Holm p < 0.05, clock-specific,
  and p < 0.05 in `<class>:corrected`; with the last condition failing: "not shown (RNA of the other class)". p <
  0.05 without clock-specificity: "not clock-specific"; unadjusted p < 0.05 lost at Holm: "suggestive, not
  established".
- **H244.2 (secondary; conditional):** cortex and medulla of one fetus share a clock state: gate over the two
  classes on the corrected, cohort- and age-adjusted counts (`development.shared`; p = the larger of p_pairing and
  p_acrophase; companions on `pairing_z`), fetuses with both pseudobulks. Scored only if, on adult material cut to
  the design (two GTEx tissues of one donor, adrenal gland and one other fixed in the addendum, thinned to the class
  counts), the pass rate is >= 0.8 and the false-positive rate with mis-paired donors has 0.05 inside its 95 % CI;
  otherwise "not estimable", not computed.
- *Descriptive (no pass / fail):* `stroma`; each class in the windows of A242 (`late`, `early`); each class without
  each cohort in turn; unsigned R; leave-one-gene-out range; the position of each fetal C among its adult
  subsamples; a depth-matched comparison (the deeper class of each cohort thinned gene by gene to the clock-gene
  counts of the shallower one; median C over 100 thinnings; difference cortex minus medulla with a 95 % interval from
  1,000 resamplings of the fetuses within cohort).
**Design addendum (before scoring; `### A244 design addendum`, committed with the prep and feasibility scripts).**
The per-cell rule and its thresholds, set on the marker counts alone; agreement with the published cell annotation
where a source releases one; per cohort and class: fetuses with >= 200 cells, cells, depth, median clock-gene
counts; a per class and cohort; the design table; the adult pass rates; the readouts scored; for H244.2 the second
adult tissue and its calibration. Nothing else may be chosen there.
**Reading.** A localisation within the fetuses of A242, not a further replication. Both classes PASS: the
covariation is present in the steroidogenic and in the neural-crest lineage; it is not cortex-only and not an effect
of class proportions. Cortex PASS, medulla negative (unadjusted p >= 0.05) at an adult pass rate >= 0.8: the cortex
carries it; a medullary clock of the size given by the power table is excluded. Medulla PASS lost after correction,
or confounded: not separable from cortex RNA, as in A190. Cortex negative at an adult pass rate >= 0.8: the
whole-adrenal result of A242 is not found within the steroidogenic cells and may rest on the proportions of the cell
classes; A242 is then qualified in the docs and goes back to the PI (Q32). H244.2 PASS: the two lineages of one
fetus share a clock state (a coupled organ clock); negative at >= 0.8: no shared state of adult strength; not
estimable: no statement.
**Limits.** (1) The same fetuses as A242 (39 of the 46 of its `all` readout). (2) Marker-defined classes; the
medulla forms during the ages covered, so its cells are few in the youngest fetuses and its composition
(sympathoblasts, chromaffin cells, precursors) changes with age. (3) The foreign-RNA correction is a marker-based
upper bound (A190: approximate); Kildisiute's counts are already ambient-corrected. (4) The adult stand-in is the
whole adrenal gland for both classes. (5) The classes differ in cells and depth: their C are compared only in the
depth-matched descriptive. (6) Limits 1 to 4 of A242.
**Stop rule:** one addendum, one scoring run; any further cut is a new amendment. **Code:**
`paper/exploratory/A244_fetal_adrenal_classes/scripts/` (prep, feasibility, score), on
`paper/common/development.py`, `score_a241.py` (contamination) and the scripts of A242.

## A245 (2026-11-12, user agreement 2026-11-12 ("yes go ahead") to the draft committed as 66fccf5, text unchanged; nothing computed; staged: hypotheses, readouts, statistic, families and readings are fixed here, the tissue list and the power go in a dated note before any observed statistic; seen: every result of A229 (tumour endothelium: VCAM1 clock-specific, ICAM1 not, implied peaks ICAM1 18.5 h, VCAM1 20.5 h, SELE 19.0 h, T/NK share 4.8 h; lung-weighted), A146 parts 1-2 and checks, A227, A228; the card `wang2024circadian`; no adhesion gene or leukocyte marker has been read against time of death or host phase in GTEx by this author, and no genome-wide GTEx rhythm table of the repository was opened for these genes): do endothelial adhesion genes and the T-cell content of healthy human tissues follow the time of day predicted by A229? **GTEx fast deaths, recorded and internal time; other; secondary (cancer section, companion of A229)**
**Why.** In mice the endothelial clock times ICAM-1 and thereby leukocyte entry into tumours (Wang et al., Cell 2024).
The human counterpart so far is A229: two secondaries pass, both primaries fail, in tumours resected in working
hours, without times. A229's implied peaks were descriptive. GTEx gives recorded times over the whole day, an
untimed host per donor and tissues from other people: the implied peaks become predictions.
**Data.** GTEx v10 counts, fast deaths (Hardy 1-2) with a time of death (`common.gtex.load`). Target tissues: every
non-cultured tissue except the brain regions, whole blood and spleen, with >= 40 such donors (list fixed in the power
note). Internal time: the untimed body host of A234 (12 body tissues, >= 5 per donor, GTEx waveform, phase oriented
by D234.0); when the target tissue is one of the 12, the host is refitted without it.
**Readouts** (per donor x tissue, log2(CPM + 1), mean over genes for a score):
- adhesion genes VCAM1, ICAM1, SELE, each as the residual of its regression on an endothelial content score (PECAM1,
  CDH5, VWF) within tissue;
- T-cell content: CD3D, CD3E, CD2, TRAC (P7 of A146); leukocyte content: PTPRC.
**Statistic.** Per tissue, the OLS slope t of y on cos(2 pi (time - h) / 24), h the predicted peak, with Hardy x sex
indicators, age and ischaemic time as covariates; S = mean of t over tissues. Null: times permuted across donors
within Hardy x sex, the same permutation for all tissues (10,000; seed 245), one-sided. Phase-free companion
statistic R = |mean over tissues of the complex cosinor t| on the same permutations (a common phase, not the
predicted one). Each statistic is computed against recorded time and against internal time.
**Specificity (same permutations).** Per readout, S for 200 expression-matched genes (gene readouts) or 200
expression-matched 4-gene sets (T-cell score), same tissues, same h: the readout is "beyond typical genes" only if S
exceeds their q95. Internal time only: the host refitted on 50 expression-matched pseudo-clock sets (A234 rule),
clock-specific only above their q95.
- **H245.1 (primary):** VCAM1, h = 20.5 h, recorded time, S > 0.
- **H245.2 (primary):** T-cell content, h = 4.8 h, recorded time, S > 0.
- **H245.3 (secondary):** ICAM1 (18.5 h), SELE (19.0 h), recorded time.
- **H245.4 (secondary):** H245.1 and H245.2 against internal time; D = S(internal) - S(recorded) with a 90 % interval
  from 2,000 donor bootstraps (the A146 question in healthy tissue: hospital clock or body clock).
Pass: Holm p < 0.05 within the family (primaries; secondaries) and S above the matched q95.
**Descriptive (no pass / fail).** R for every readout; per-tissue amplitude and peak (lung named, as A229 is
lung-weighted); PTPRC; ventilated donors (Hardy 0) against internal time; T-cell content after residualising on HBB
(blood in the sample); offset sweep of h (+-2, 4, 6 h).
**Power first (dated note, no observed statistic):** tissues and donors; null SD of S and R per readout; MDE_80 =
2.49 x null SD; the amplitude (log2 units) that MDE corresponds to in the median tissue. A hypothesis is scored only
with >= 10 tissues.
**Reading.** H245.1 pass: VCAM1 of healthy human tissues varies over the day at the hour the tumour endothelial clock
state implied; the gate has a human, timed correlate outside tumours. H245.2 pass: tissue T-cell content peaks in
the predicted late night. Both with H245.4 favouring internal time: body time, not the hospital clock, carries it
(the opposite of A146 part 1 in tumours). R passes while S fails: a daily pattern at another hour; A229's implied
peaks do not transfer. Fail below the MDE: not detected at that size in bulk tissue, not a refutation of the mouse
mechanism.
**Limits.** (1) Bulk tissue: endothelial expression is not separated from other cells beyond the content score; the
T-cell score includes intravascular blood, whose lymphocyte count is itself circadian, so H245.2 cannot tell entry
into tissue from blood content (HBB sensitivity only). (2) Post-mortem donors; time of death is associated with the
cause of death (strata Hardy x sex only). (3) The predictions come from tumour endothelium, lung-weighted, and are not
blind to A229. (4) Controlled access: aggregates only.
**Stop rule:** one power note, one scoring run. **Code:** `paper/exploratory/A245_trafficking_gtex/scripts/`
(`common.gtex`, the A234 host).

**A245 power note (2026-11-12, before any observed statistic; `A245_trafficking_gtex/output/{power,tissues,direction}.csv`).**
233 fast-death donors with a time of death and a body host; 26 target tissues with 41 to 216 donors (11 read on a
leave-one-out host). Hosts read time of death (Fisher-Lee 0.380 full, 0.360-0.401 leave-one-out; none negated). Null SD of
S (10,000 permutations within Hardy x sex) and MDE_80, recorded / internal time: VCAM1 0.405 / 0.405, MDE 1.01 / 1.01
(amplitude 0.105 / 0.115 log2 units in the median tissue); T-cell content 0.312 / 0.315, MDE 0.78 / 0.78 (0.054 / 0.055);
ICAM1 0.490 / 0.486, MDE 1.22 / 1.21 (0.150 / 0.158); SELE 0.454 / 0.448, MDE 1.13 / 1.11 (0.224 / 0.236). Null q95 of R:
0.76-1.19. All hypotheses are scored (26 tissues >= 10). Code checked on shuffled times only (`--placebo`, subsample):
the fast t agrees with statsmodels OLS, the identity path equals the permutation path, and unstratified shuffles are
rejected at 5.0 % / 5.3 % against the stratified null. Choices not fixed by the registration are listed in the script
docstring (donor set, strata merging, matched companions, bootstrap without host refit, Holm over the four
secondaries). Scoring proceeds as registered.

## A246 (2026-11-12, user agreement 2026-11-12 ("yes go ahead") to the draft committed as 66fccf5, text unchanged; nothing computed; staged: programs, state, statistic, families and readings are fixed here, eligibility and power go in a dated note before any observed association; seen: A146 parts 1-2 and checks (programs against recorded time pass in breast; against internal time not detected, largest t: T-cell activation 1.83, PD-L1 1.31), every result of A227, A228 and A229 on the same cohorts, including H229.5 (T/NK share, pass) and the per-cancer pattern (lung positive, colorectal at the pseudo level); no program below has been related to any clock state): do tumour-microenvironment programs follow the clock state of the tumour, at any phase? **Single cell, A227 cohorts; L1-like association (A229 statistic); other; secondary (cancer section)**
**Why.** A146 part 2 asked whether ten programs follow internal time at phases transferred from a mouse atlas (mouse
CT - 5 h) and on a noisy per-patient time. Its failure is conditional on both. The A229 statistic needs neither: it
asks whether a readout varies with a clock state across patients at any phase, against expression-matched
pseudo-clock sets. A229 applied it to three adhesion genes and two cell shares only.
**Data.** The 34 A227 cohorts and tables; pseudobulks per patient x compartment with >= 50 cells (A227 eligibility).
**Programs** (genes, scores and residual definitions of A146; scored in the compartment of A146 part 2): P1
proliferation, P2 MHC-I, P3 interferon response, P9 PD-L1 in malignant cells; P4 T-cell activation and P10 PD-1
(residuals on T abundance) in T/NK; P5 myeloid antigen presentation and P6 phagocytosis (residuals on myeloid
abundance) in myeloid cells; P8 myeloid share of all tumour cells. P7 (T-cell share) is H229.5 and is not rescored.
A program enters with a median pseudobulk count >= 10 summed over its genes; no program gene is a clock gene.
**States.** (a) *Body time in the tumour:* z of the tumour endothelial pseudobulk, exactly as A229. (b) *Own clock:*
z of the compartment in which the program is scored (same construction, fitted on that compartment's rows; P8 on
myeloid). Both centred within cohort x chemistry.
**Statistic.** r of A229 per program and state, pooled over cohorts with patient-count weights; null: patients of y
permuted within cohort, the same permutation for all programs (10,000; seed 246). Omnibus per state: T = mean over
programs of (r - null mean) / null SD. Clock specificity: z recomputed from 200 expression-matched pseudo-clock sets
(A227 construction); a program or the omnibus is clock-specific only above the pseudo q95.
- **H246.1 (primary):** omnibus T > 0 with state (a).
- **H246.2 (primary):** omnibus T > 0 with state (b). Holm over the two.
- **H246.3 (secondary):** each program with each state; Holm over the programs within a state.
Pass: Holm p < 0.05 and above the pseudo q95.
**Descriptive (no pass / fail).** Percentile of each program among 200 expression- and size-matched gene sets of its
compartment used as y; implied peak per program (`implied_peak_h`, as A229) beside the hour A146 transferred from
the mouse, and their circular correlation over programs; per cancer type with >= 12 patients; the same tests in
adjacent-normal samples; breast cohorts alone (the cancer of A146 part 1).
**Power first (dated note, no observed association):** patients and cohorts per program and state; null SD; MDE =
2.49 x SD; a program is not scored with < 40 patients or MDE > 0.25 (A229 rule), the omnibus not with < 5 programs.
**Reading.** H246.1 pass: programs of the tumour microenvironment vary with the body time read in the tumour's own
endothelium, in untimed human resections; with the implied peaks agreeing with the mouse order, the mouse
organisation holds in humans and A146 part 2 failed on the transferred phases or the time estimate. H246.2 pass
without H246.1: programs follow the clock state of their own cells, not the body time of the tumour. Both fail below
the MDE: with A146 part 2, no clock-state organisation of these programs is detected at any phase; the
recorded-time result of A146 part 1 then stands alone as a hospital-clock effect awaiting a second timed cohort.
**Limits.** (1) Daytime surgery compresses the spread of states. (2) Association: a systemic cue can move a clock
state and a program together. (3) Not blind to A227-A229 on the same patients; lung is expected to carry any signal.
(4) State (b) in malignant cells rests on a weak clock signal (A138 H138.1); a negative there is read with its MDE.
(5) Ambient RNA shared within a sample is controlled by the pseudo-clock sets only.
**Stop rule:** one power note, one scoring run. **Code:** `paper/exploratory/A246_tme_clock_state/scripts/` (on
`A229_endothelial_gate/scripts/compute_a229.py` and the A146 program definitions).

**A246 power note (2026-11-12, before any observed association; `A246_tme_clock_state/output/{eligibility,power}.csv`; script written by a Sonnet agent, numbers re-read in the table).**
Null from y shuffled within cohort only (10,000; one ordering of a cohort's patients per iteration, shared by all
programs). Tumour side, state (a) endothelial: P1-P3 134 patients / 18 cohorts (MDE 0.104-0.106), P9 73 / 9 (0.143), P4 and
P10 164 / 18 (0.092, 0.093), P5, P6 and P8 169 / 19 (0.089-0.092). State (b) own compartment: P1-P3 298 / 31 (0.071-0.073),
P9 162 / 16 (0.097), P4 and P10 318 / 29 (0.070), P5, P6 and P8 314 / 29 (0.070-0.071). All nine programs are scored on
both states. Omnibus T: null SD 0.340 (a) and 0.349 (b), MDE 0.85 and 0.87 (T units; the SD is close to 1 / sqrt(9), so
the programs are nearly independent under the null). Adjacent-normal side (descriptive): 75-263 patients per program,
all scored. P9 (CD274) loses 14 cohorts to the count rule; 10 to 12 cohorts have no tumour endothelial pseudobulk.
Choices not fixed by the registration (script docstring): abundance denominators are the cells of the A227 tables
(as H229.5), patients with < 100 such cells dropped for P4-P6, P8, P10; malignant compartment = copy-number calls where
the cohort has them, else tumour epithelium; a program needs half of its genes; z-scores over the patients with the
compartment. Code tested only on shuffled y (`--placebo`, 3 cohorts). Scoring proceeds as registered.

## A247 (2026-11-12, user agreement 2026-11-12 ("yes go ahead with A247") to the draft committed as 729a317, text unchanged; nothing computed; staged: classes, statistics, hypotheses, families and readings are fixed here, eligibility, the healthy-brain calibration and the power go in a dated note before any clock statistic is computed on GBmap; seen: the scout `exploratory/scouts_2026-11-12/glioma_scout.md` (cells, patients, summed clock-gene UMI and MGMT UMI per class; no clock statistic, no gene-by-gene clock counts), every carcinoma result (A138 H138.1: malignant coupling c 0.076; A149, A151, A227-A229), the brain results A135 / A141 / WP0 Q3 (no principled brain reference; glia against neurons is anchor-dependent), the cards `dong2019targeting` and `gonzalezaponte2025circadian`): do glioma cells share a clock state with the patient's own non-malignant brain cells, and does MGMT follow that state? **Extended GBmap (CELLxGENE), droplet data; L3-like pairing within the tumour sample and L1-like association; secondary (cancer section)**
**Why.** Glioblastoma stem cells keep robust clocks in culture and depend on BMAL1 / CLOCK (Dong et al. 2019), unlike
carcinoma cells, whose clock state is barely tied to the patient's other cells (A138 H138.1). MGMT promoter
methylation calls vary with the hour of surgery between 10:00 and 18:00, and MGMT oscillates in synchronised glioma
lines (Gonzalez-Aponte et al. 2025); no tissue clock readout exists in patients. No glioma cohort has been analysed.
**Data.** Extended GBmap, droplet libraries only (Smart-seq2 out), published annotation. Classes: `malignant`
(annotation_level_1 Neoplastic), `myeloid` (TAM-BDM, TAM-MG, Mono, DC; neutrophils and mast cells out),
`oligodendrocyte` (mature oligodendrocytes only: OPC are left out, since OPC-like malignant cells mislabelled as OPC
would create pairing). One pseudobulk per patient x class with >= 100 cells; strata = study x assay (cells /
nuclei); a study enters a test with >= 4 patients holding both pseudobulks. Duplicated cells (`is_primary_data`
False) once.
**States and statistics.** Complex clock state z per pseudobulk (`clocksync.paired`, GTEx waveform, counts, projector
fitted per class and study, centred within stratum), as A227. *Pairing* of two classes of one patient: the A227
statistic, pooled over studies with patient-count weights; null: rows of the second class permuted across patients
within stratum (10,000; seed 247). *Association* of a real readout y with a state: the A229 statistic r, same null.
**Clock specificity.** 200 expression-matched pseudo-clock sets (A227 construction) and 200 responsive pseudo-clock
sets (the stress- and glucocorticoid-responsive construction of A109 H109.1-R / A242; most patients receive
dexamethasone before surgery): clock-specific only above the q95 of both.
**RNA of the other class (same sample).** As A244: per patient, the share of class b's marker counts found in the
class a pseudobulk relative to class b's own (markers, never clock genes: myeloid PTPRC, CSF1R, C1QA, C1QB, TYROBP;
oligodendrocyte MBP, PLP1, MOG, MOBP; malignant SOX2, PTPRZ1, EGFR) bounds the foreign share; the expected foreign
clock-gene counts are subtracted and the test rerun (`:corrected`). A pass must hold in `:corrected`.
- **P247.0 (precondition, scored first):** myeloid vs oligodendrocyte pairing is above the null (p < 0.05) and
  clock-specific: two non-malignant classes of a tumour-bearing brain share a clock state. If it fails, H247.1-2 are
  reported but read as uninformative.
- **H247.1 (primary):** malignant vs myeloid pairing. **H247.2 (primary):** malignant vs oligodendrocyte pairing.
  Holm over the two; pass = Holm p < 0.05, clock-specific, and p < 0.05 in `:corrected`.
- **H247.3 (primary, MGMT):** y = log2(CPM + 1) of MGMT in the malignant pseudobulk (>= 10 MGMT UMI), centred within
  stratum, against the myeloid state of the same patient (a state not read from the cells that express y). Pass: p <
  0.05 and clock-specific.
- **H247.4 (secondary):** the same against the malignant cells' own state; and the methylation call against the
  myeloid state: T = |mean z (methylated) - mean z (unmethylated)| pooled over strata holding both calls, labels
  permuted within stratum, both companions (the untimed counterpart of the hour-of-surgery effect). Holm over the two.
**Descriptive (no pass / fail).** MGMT percentile among 200 expression-matched malignant genes used as y; implied
peak of MGMT (`implied_peak_h`); H247.3 within methylation status; cells and nuclei separately; primary vs recurrent;
T/NK vs myeloid pairing; clocksync coupling c of the malignant class on a host anchored on myeloid and
oligodendrocyte cells, beside the carcinoma value of H138.1 (only if the gate and its companions pass); per study
with >= 12 patients.
**Power first (dated note, no observed statistic on GBmap).** Patients and studies per test; null SD from shuffled
pairings only; MDE = 2.49 x SD; a test is not scored with < 40 patients or MDE > 0.25 (A227 rule). Calibration on
healthy brain: microglia vs oligodendrocyte pairing in an on-disk control atlas (PsychAD or SEA-AD controls),
pseudobulks thinned to the cells and clock-gene counts of the GBmap design; the pass rate (200 subsamples) says
whether a negative P247.0 would be informative; without such an atlas on disk the note says so and P247.0 stands on
its MDE. The coding of the MGMT labels is fixed there from the metadata alone.
**Reading.** H247.1 or H247.2 pass (with P247.0): glioma cells stay on the clock state of the patient's brain cells,
as their retained clock in culture predicts and unlike carcinoma cells; the loss of host coupling is then a property
of carcinomas, not of malignancy. Fail below the MDE with P247.0 passing: glioma cells are not tied to the patient's
glia at that size, as in carcinomas, although their clock runs in culture. P247.0 failing: no statement; the
calibration says whether tumour-bearing brain differs from healthy brain at this depth. H247.3 pass: MGMT expression
varies with the body time read in the tumour, the first in-patient support for a clock-timed MGMT. H247.4
(methylation) pass: the call depends on the clock state at resection, the claim of Gonzalez-Aponte et al. without
recorded hours. Fails are read with their MDE: daytime surgery compresses the spread of states.
**Limits.** (1) Both classes come from one sample: shared ambient RNA and handling; the pseudo-clock sets and the
`:corrected` rerun are the controls, the marker bound is approximate (A190). (2) Dexamethasone, temozolomide and
radiotherapy are not recorded; stage is unknown for most patients. (3) No host anchor learned for brain
(neurons and endothelium are too rare here): pairing between classes, not coupling to a host, is the readout.
(4) Published malignant calls (copy-number based) and labels; sorted or enriched studies change which patients
qualify, not the statistic. (5) Methylation status is known for about half of the patients, with mixed codes.
(6) Association, not causation; no surgery times.
**Stop rule:** one power note, one scoring run; Wang 2022 (GSE174554, tumour / normal calls only) is not used here.
**Code:** `paper/exploratory/A247_glioma/scripts/` (on `compute_a227.py`, `compute_a229.py`, `score_a241.contamination`).

## A248 (2026-11-12, user agreement 2026-11-12 ("yes go ahead with A248") to the draft committed as e03899a, text unchanged; nothing computed; user: "draft the glucocorticoid T-cell amendment", after remarking that glucocorticoid control of T-cell localisation is solid in the literature and that the A245 T-cell result is unexpected; staged: readouts, predicted hours, statistics, families and readings are fixed here, the power goes in a dated note before any observed statistic; seen: every result of A245 on the 26 tissues, including the per-tissue peaks of T-cell content (pooled 16.4 h, lung 12.6 h, 19 of 26 tissues between 9 and 21 h; phase-free p 0.077) and of the adhesion genes; NOT seen: any readout of GTEx whole blood, any glucocorticoid-response score in GTEx, any naive or effector marker score): does T-cell content of human blood follow the published glucocorticoid-driven rhythm in GTEx, do tissues carry a glucocorticoid-response rhythm, and is tissue T-cell content tied to it? **GTEx fast deaths, recorded time; blood positive control and glucocorticoid score blind, tissue T-cell hour post hoc to A245; other; secondary**
**Why.** A245 tested T-cell content at an hour taken from A229 (4.8 h) and found the opposite sign. That hour did not
come from the redistribution literature. In healthy men blood T-cell counts peak at night and fall in the day, 3 h
behind cortisol, by redistribution of naive cells (to bone marrow, through CXCR4): CD3+ acrophase 02:10 h (SEM 20
min), amplitude 388 on a mesor of 1,902 cells per microlitre, nadir about 14:00 h; naive CD4+ 02:41 h, naive CD8+
02:09 h (night-to-day change -41 % and -46 %); effector CD8+ 15:34 h (SEM 116 min); cortisol maximum about 08:00 h
(cosinor acrophase 10:40 h), correlation with T-cell counts maximal at a 3 h delay (r -0.65) (Dimitrov et al., Blood
2009, Table 1, Table 2 and Results, read in PMC2686184; 14 men, 24 h, sampled every 1.5 to 3 h). Blocking cortisol
synthesis or the glucocorticoid receptor blunts the morning fall of naive T cells (Besedovsky et al., FASEB J 2014;
abstract only); in mice the rhythm is lost without the T-cell glucocorticoid receptor (Shimba et al., Immunity 2018;
abstract only). Cards for the three papers are written before the result entry. A tissue sampled in bulk holds
blood and resident cells: if its T-cell signal is blood, it peaks at night; if cells that leave blood by day are
found in tissue, it peaks by day (what A245 hints at, post hoc).
**Data.** As A245 (GTEx v10, fast deaths with a time of death, the 233 donors of A245, same covariates, strata and
permutation scheme), plus whole blood (excluded from A245) as its own compartment. Code: the functions of
`compute_a245.py`.
**Readouts** (log2(CPM + 1), mean over genes): T-cell content (CD3D, CD3E, CD2, TRAC, as A245); glucocorticoid
response GC (TSC22D3, FKBP5, DDIT4, ZBTB16, KLF9: the A146 set; PER1 left out, a clock gene); naive (CCR7, SELL,
LEF1, TCF7) and cytotoxic effector (GZMB, PRF1, CX3CR1, FGFBP2) scores, descriptive.
**Statistics.** S, R and the permutation null of A245 (10,000 permutations of time within Hardy x sex; seed 248);
blood is one compartment, so S is its t. Matched companion: 200 expression-matched gene sets of the same size
(A245 rule); beyond typical genes only above their q95. Association of two readouts in the same sample: partial
correlation of T-cell content with GC given the covariates, mean of Fisher z over tissues; null: the GC vector
permuted across donors within Hardy x sex (shared across tissues); companion: GC replaced by 200 matched 5-gene
sets, on the same samples.
- **H248.1 (primary; blind; positive control of the readout):** blood T-cell content peaks at 02:10 h (S > 0,
  recorded time).
- **H248.2 (primary; blind):** the GC score of the 26 tissues peaks at 10:40 h, the cortisol acrophase (S > 0,
  recorded time); a transcript lag of up to 3 h costs at most 30 % of the projection. Holm over H248.1-2.
- **H248.3 (secondary; blind):** in blood, T-cell content is lower where GC is higher (partial correlation < 0).
- **H248.4 (secondary; POST HOC to A245, labelled as such wherever reported):** tissue T-cell content peaks at
  14:00 h, the blood nadir (S > 0 over the 26 tissues, recorded time). Holm over H248.3-4.
Pass: Holm p < 0.05 and above the matched q95.
**Descriptive (no pass / fail).** GC in blood and per tissue (peak, amplitude); partial correlation of T-cell content
with GC over the 26 tissues, with and without cos / sin of time among the covariates (does a donor's glucocorticoid
state relate to tissue T-cell content beyond the hour?); naive and effector scores in blood against 02:09-02:41 h
and 15:34 h; every readout against internal time (the A245 host) and D = S(internal) - S(recorded); ventilated
donors; blood T-cell content after residualising on a neutrophil score (CSF3R, FCGR3B, CXCR2).
**Power first (dated note, no observed statistic):** donors with blood; null SD and MDE (2.49 x SD) of every test;
for H248.1 the log2 amplitude at the MDE beside the published one (a peak 20 % above the mesor, about 0.27 log2
units, if the share of T-cell RNA follows the count).
**Reading.** H248.1 pass: the GTEx readout recovers the known human rhythm; A245's failure is then about the A229
hour, not about the ability to see T-cell rhythms in these donors. H248.1 fail at an MDE below the published
amplitude: the readout or these donors (post-mortem, terminal stress, relative RNA shares) do not carry the known
rhythm, and no T-cell statement from GTEx, A245 included, should be read as biology. H248.2 pass: a
glucocorticoid-response rhythm at the cortisol hour is visible across human tissues without sampling times being
controlled; with H248.3, the glucocorticoid link of the blood rhythm is seen in the same samples. H248.4 pass
(post hoc): tissue T-cell content is in antiphase to blood, as cells leaving blood by day would give; it needs an
independent cohort before it is a claim. H248.4 fail with H248.1 passing: bulk tissue content does not mirror the
blood rhythm at that size.
**Limits.** (1) RNA shares, not cell counts: a T-cell score in blood falls when neutrophils rise (neutrophil
sensitivity only). (2) Death is not a resting baseline: cortisol rises with terminal stress, and the hour of death
is tied to its cause (strata Hardy x sex only). (3) Bone marrow, the proposed daytime destination of naive T cells,
is not in GTEx; lymphoid organs are not targets. (4) The published hours come from 14 young men under a fixed sleep
schedule; donors differ in age, sex and habits. (5) H248.4 and the tissue descriptives are not blind. (6) Controlled
access: aggregates only.
**Stop rule:** one power note, one scoring run. **Code:** `paper/exploratory/A248_glucocorticoid_tcells/scripts/`.

**A248 power note (2026-11-12, before any observed statistic; `A248_glucocorticoid_tcells/output/{power,compartments}.csv`).**
208 of the 233 donors have whole blood; 26 tissues as A245. Null SD and MDE_80 (10,000 permutations): H248.1 blood T-cell
content 1.007 (t units), MDE 2.51, i.e. an amplitude of 0.274 log2 units, the size of the published rhythm (about 0.27 if
RNA share follows the count): power is about 0.8 only for the full published amplitude, so a failure of the positive
control is informative only for an undiminished rhythm. H248.2 GC over 26 tissues: null SD 0.641 (the score is strongly
shared between the tissues of a donor), MDE 1.60, amplitude 0.172 log2 in the median tissue. H248.4 tissue T-cell content:
null SD 0.316, MDE 0.79 (0.059 log2). H248.3: null SD of Fisher z 0.071, MDE partial r 0.175. Code checked on shuffled
times and shuffled GC only (`--placebo`, subsample; unstratified GC shuffles fall outside the null's central 95 % in 6.3 %
of 400 draws, null SD 0.084 against 0.085). Choices not fixed by the registration are in the script docstring. Scoring
proceeds as registered.

**A247 power note (2026-11-12, before any observed statistic on GBmap; `A247_glioma/output/{power,eligibility,calibration,foreign_rna,prep_summary,patients}.csv`; scripts by a Sonnet agent, tables re-read). NOT SCORED YET: the calibration below goes to the user first.**
Pseudobulks with >= 100 cells: malignant 129 patients, myeloid 162, mature oligodendrocyte 65, T/NK 76. Null from
shuffled rows only (10,000). P247.0 myeloid vs oligodendrocyte: 53 patients / 8 strata, null 0.373 (SD 0.068), MDE 0.168.
H247.1 malignant vs myeloid: 105 / 11, MDE 0.118. H247.2 malignant vs oligodendrocyte: 42 / 5, MDE 0.186. H247.3 MGMT vs
myeloid state: 85 / 10, MDE 0.140. H247.4 MGMT vs own state: 93 / 10, MDE 0.127; methylation call vs myeloid state: 75
patients / 7 strata (32 methylated, 43 unmethylated), MDE 0.69 on the z scale (0.32 of the state's RMS): the 0.25 rule
was written for correlations and does not apply; by the letter it is not scored. MGMT labels: M, Meth = methylated; UM,
Un, NM = unmethylated; NT and unknown = no call (36 / 60 / 86 of 182 patients; no conflict). Foreign-RNA bound (median a):
0.014-0.027 between malignant, myeloid and oligodendrocyte pseudobulks, except malignant RNA in oligodendrocytes 0.110
(maximum 0.61); myeloid RNA in T/NK 0.258.
**Calibration on healthy brain (microglia vs oligodendrocytes of control donors, thinned to the P247.0 design; 200
subsamples, 500 permutations):** pass rate 0.225 [0.169, 0.289] (PsychAD HBCC, 119 donors) and 0.170 [0.121, 0.229]
(PsychAD MSSM, 222); with donors mis-paired 0.145 [0.099, 0.202] and 0.120 [0.078, 0.173], above the nominal 0.05; at the
atlases' own depth the pass rate is 0.935 (HBCC) and 0.315 (MSSM). SEA-AD has too few control donors (44, 45); no
pseudo-clock companion (the atlases on disk hold a 45-gene panel). Agent's diagnostic (not in the tables): with one
common thinning depth the mis-paired rate is 0.05 and the pass rate 0.05. Two consequences: (i) at GBmap depth a
healthy brain passes P247.0 little more often than mis-paired donors do, so a negative P247.0 is uninformative; (ii)
when the two classes of a patient share their depth, the within-stratum permutation p of the pairing statistic is
anti-conservative (state moduli covary although phases do not), so a positive p alone is not trustworthy either; the
pseudo-clock companions, which share the depth structure, are then the only valid control.
**Choices not fixed by the registration** (docstrings): cells flagged `is_primary_data` False kept (each cell occurs
once in the extended file); projector per class and stratum; pseudo-clock sets drawn once per test; stratum = study x
suspension; foreign-RNA share clipped to [0, 1]. Code run only on shuffled pairings (`--placebo`).

**A247 note 1 (2026-11-12, user: "go ahead with option 1"; before any observed statistic is computed by the tests stage): the calibration of the power note was flawed; corrected, the registered null is calibrated and nothing in the tests changes.**
(i) *Flaw.* `compute_a247.thin` thinned every clock gene of a healthy donor to the count of the same gene in a GBmap
patient. Where the atlas is deeper than GBmap this replaces the donor's clock-gene pattern by the GBmap patient's. The
rates of the power note (`calibration.csv`: pass 0.225 / 0.170, mis-paired 0.145 / 0.120) therefore describe GBmap
patterns with binomial noise, not healthy brain, and its two consequences are withdrawn. (ii) *Disclosure.* Through
that flaw, those rates and a first run of `calibrate_null.py` with the same thinning (mis-paired rates of 0.10-0.19 and
pass rates of 0.18-0.40 for the P247.0 and H247.1 designs) are noisy functions of the real GBmap pairings: they hint
that the observed statistics of P247.0 and H247.1 are not far above their permutation nulls. No choice below was made
from them. (iii) *Corrected calibration* (`calibrate_null.py`, one thinning probability per sample so that the donor
keeps its own pattern; 200 subsamples, 500 permutations; `output/calibration_null_<design>_<assay>.csv`). Registered
statistic and permutation null, rate of p < 0.05 with donors mis-paired, PsychAD HBCC / MSSM: P247.0 design 0.070 / 0.050;
H247.1 design 0.065 / 0.070; H247.2 design 0.075 / 0.060 (each from 200 subsamples; 95 % binomial intervals of about
+-0.035 include 0.05). True pairs at GBmap depth: P247.0 0.610 / 0.240; H247.1 design (microglia at the myeloid depth,
oligodendrocytes at the malignant depth) 0.965 / 0.435; H247.2 design 0.645 / 0.210; at the atlases' own depth 0.950 /
0.345 (P247.0 design), so the two atlases differ in the strength of glial pairing itself. Variants tried on the
corrected calibration (unit-modulus states, states residualised on their clock-gene count, permutation within depth
blocks) are no better calibrated and less powerful; residualising on both counts and a sign-flip null are invalid.
(iv) *Decision.* The registered statistic, null, companions and pass rules stand unchanged. Readings added: a negative
P247.0 is weakly informative (a healthy-brain pass rate of 0.24 to 0.61 at this depth); H247.1 has a power of 0.44 to
0.97 for malignant cells tied to myeloid cells as healthy oligodendrocytes are to microglia, H247.2 of 0.21 to 0.65.
H247.4 (methylation call) stays unscored by the letter of the power rule. No pseudo-clock calibration is possible (45-gene
atlases).

## A249 (2026-11-12, user agreement 2026-11-12 ("yes go ahead with A249") to the draft committed as 4a1644c, text unchanged; nothing computed; staged: gene sets, states, statistics, hypotheses, families and readings are fixed here, eligibility and power go in a dated note before any observed association; seen: A138 H138.4 (tumour and normal clock states of a patient pair, 0.365, clock-specific), A151, A207, A221 (clock-gene levels, tumour vs normal), A229, A246 (programs do not follow the tumour's clock state; in adjacent-normal tissue they do, at the line), A247 (same-sample pairing is not clock-specific), archive A100 (drug targets in healthy GTEx tissue); NOT seen: any non-clock gene of a tumour or normal sample against the clock state of the patient's other sample): are drug-relevant genes on body time in normal tissue and off it in the tumour of the same patient? **TCGA tumour-normal pairs primary, single-cell cohorts and the timed breast cohort secondary; L1-like association across samples (A229 statistic); secondary (cancer section)**
**Why.** Timed chemotherapy rests on a difference between host and tumour: toxicity pathways of normal tissue follow
the day while the tumour's do not, so that a dosing hour can spare the host (Levi et al. 2010, which names the genes
of set C below; Levi et al. 2024). The premise comes from rodents and cell lines, trials are mixed (Sancar and Van
Gelder 2021), and it has not been tested in patients: tissue is taken once, in working hours, without a usable time.
A patient's two samples are taken at the same operation. Body time can therefore be read from one sample and a
gene from the other, which removes the same-sample covariation that made A247 unreadable.
**Data.** `common/datasets/tcga_paired` (recount3 counts; one primary tumour and one solid-tissue normal per patient),
the cancers and pairs of H138.4. Secondary: the A151 single-cell cohorts (patients with fibroblast pseudobulks in
both the tumour and the adjacent-normal sample); GSE233242 (42 breast patients, tumour and normal, recorded
resection hour).
**States.** z_N, z_T: complex clock state of the normal and of the tumour sample (`clocksync.paired`, GTEx waveform,
counts, one projector per cancer as H138.4), centred within cancer.
**Gene sets (fixed now, from the literature; no data used).**
- *C, chronopharmacology* (genes named in Levi et al. 2010 for 5-FU, irinotecan, platinum handling and apoptosis):
  DPYD, TYMS, TOP1, UGT1A1, CES1, CES2, ABCB1, ABCC2, CYP3A4, ERCC1, XPA, BCL2, BAX, TP53, WEE1.
- *K, cell-cycle gating and repair* (Matsuo et al. 2003; Fu et al. 2002; Grechez-Cassiau et al. 2008;
  Gonzalez-Aponte et al. 2025): WEE1, CCNB1, CDK1, MYC, CCND1, GADD45A, MDM2, CDKN1A, MGMT.
- *F, approved drug targets:* the Human Protein Atlas list of A100 (`$SCCC_DATA/drug_targets`).
- *O, clock-controlled genes outside the 14-gene waveform (positive control):* BHLHE40, BHLHE41, NFIL3, HLF, RORC,
  ARNTL2.
A gene enters a cancer with a median count >= 10 in both sample types; a set enters with >= half of its genes
(>= 100 genes for F) in >= 5 cancers.
**Statistic.** Per gene and cancer, r of A229 between y (log2(CPM + 1), centred within cancer) and a state; pooled
over cancers with patient-count weights; a set's S is the mean of r over its genes. Two cross-sample readouts:
S_N = S(y of the normal sample, z_T) (is the normal tissue's gene on the body time read in the tumour?) and
S_T = S(y of the tumour sample, z_N). Null: patients of y permuted within cancer, one permutation for all genes
(10,000; seed 249). Clock specificity: the state recomputed from 200 expression-matched pseudo-clock sets (A227
construction), above their q95. Typical genes: 200 expression-matched gene sets of the same size as y, above their
q95 ("beyond typical genes"). Window: D = S_N - S_T with a 90 % interval from 2,000 patient bootstraps within
cancer (states not refitted).
- **P249.0 (precondition):** set O, S_N: p < 0.05 and clock-specific (body time read in the tumour sample orders
  clock-controlled genes of the normal sample). If it fails, everything below is uninformative.
- **H249.1 (primary):** set C is on body time in normal tissue: S_N, p < 0.05, clock-specific and beyond typical genes.
- **H249.2 (primary; the window):** H249.1 passes, S_T of set C is not detected (p >= 0.05) at an MDE below the
  excess of S_N over its null mean, and the 90 % interval of D excludes 0. z_T is the noisier reading of body time
  (H138.4), which works against S_N, so the asymmetry is conservative for a window.
- **H249.3 (secondary):** the same two tests for sets K and F (Holm over the two sets, per test).
- **H249.4 (secondary; cell-type resolved):** single cell, body time from the fibroblasts of the other sample: set C
  in normal epithelium against the tumour-sample fibroblast state, and in malignant cells against the normal-sample
  fibroblast state; same statistics.
- **H249.5 (secondary; timed):** GSE233242: mean over the genes of set C of the two-term cosinor t modulus against the
  recorded resection hour, normal minus tumour, hours permuted across patients (the window against a clock on the wall).
**Descriptive (no pass / fail).** Per gene S_N, S_T and implied peak; per cancer with >= 30 pairs; the same-sample
readouts (y and z of one sample), labelled as confounded; set O in the tumour (S_T): are clock-controlled genes
themselves off body time in tumours; S_N and S_T against tumour purity tertiles (ABSOLUTE, on disk).
**Power first (dated note, no observed association):** pairs and cancers per set; null SD and MDE (2.49 x SD) of S_N,
S_T; not scored with < 100 pairs (TCGA), < 40 patients (single cell) or an MDE above 0.10; for H249.5 the null SD only.
**Reading.** H249.1 and H249.2 pass: in patients, genes that handle anticancer drugs are on body time in normal
tissue and not in the tumour: the first direct support of the chronotherapy premise in human tissue. H249.1 passes
and S_T passes too: these genes stay on body time in the tumour; there is no window at the transcript level (in
line with the mixed trials), although timing could still matter for toxicity alone. H249.1 fails with P249.0
passing: these genes are not detectably on body time in human normal tissue in working hours; the rodent premise is
not seen at the MDE. The comparison of sets C, K and F says whether any of this is particular to drug handling.
**Limits.** (1) Bulk samples: a gene of the tumour sample is read in a mixture of malignant and other cells (purity
descriptive; H249.4 resolves cell types on fewer patients). (2) Surgery in working hours compresses the spread of
body time; both sides are compressed alike. (3) "Normal" is tumour-adjacent tissue. (4) Transcripts, not enzyme
activity or drug levels. (5) A patient-level state other than time that moves clock genes in one sample and a gene
in the other would pass; the pseudo-clock sets are the control. (6) Set C is a literature list drawn up for
colorectal chemotherapy; several genes are liver- or gut-specific and will drop out elsewhere.
**Stop rule:** one power note, one scoring run. **Code:** `paper/exploratory/A249_chrono_window/scripts/` (on
`common/datasets/tcga_paired`, the H138.4 states and `compute_a229.py`).

**A249 power note (2026-11-12, before any observed association; `A249_chrono_window/output/power.csv`).**
Null from patients of y permuted within block only (10,000). *TCGA:* the 10 gated cancers of H138.4, 517 pairs (19 to 112
per cancer; states reproduce `q3_lib.tcga_states`). Genes scored / eligible in >= 5 cancers: C 15 / 14, K 9 / 9, O 6 / 6, F
761 / 623 of 854. Null mean of S about 0.120; null SD and MDE of S_N / S_T: C 0.0092 / 0.0070, MDE 0.023 / 0.018; K MDE 0.024 /
0.021; O 0.027 / 0.025; F 0.016 / 0.011. All sets scored. *Single cell (H249.4):* 13 cohorts hold patients with both
fibroblast and both epithelial pseudobulks, 61 patients in blocks of 4 to 7; null mean of S about 0.50 (blocks of four
patients), MDE 0.058-0.074 (C, K, O) and 0.025 (F): scored by the letter (>= 40 patients, MDE <= 0.10), but with blocks
this small the statistic is dominated by its small-sample bias and only the companions can carry a claim. The tumour
side is all tumour epithelium (the genome-wide tables have no copy-number split), not called malignant cells: a
deviation from "malignant cells" in the registration, forced by the data. *Timed (H249.5):* 42 patients; genes with TPM
>= 1 in both tissues: C 12, K 9, O 6, F 475; null SD of normal minus tumour modulus 0.460 (C), 0.527 (K), 0.538 (O), 0.307
(F). Code run only on shuffled patients (`--placebo --smoke`). Choices not fixed by the registration are in the script
docstring. Scoring proceeds as registered.

## A250 (2026-11-12, user agreement 2026-11-12 ("Yes go ahead") to the draft committed as b1b5bc0, text unchanged; user: "draft the confirmatory test on CPTAC pairs"; nothing computed; staged: gene sets, states, statistics, hypotheses, families and readings are fixed here, eligibility and power go in a dated note before any observed association; seen: every result of A249 in TCGA, including the descriptive ones that motivate this amendment (set O in the tumour against the normal sample's state: S_T 0.162, p 0.0001, pseudo q95 0.147, matched q95 0.136; set C S_T 0.117, set F 0.122; per gene NFIL3 0.206, RORC 0.195, ARNTL2 0.160; per cancer), A186 (CPTAC pair counts, gates and tumour-normal pairing of the clock states, by proliferation tertile), the A223 power table; NOT seen: any non-clock gene of a CPTAC sample against the clock state of the patient's other sample): are a tumour's clock-controlled genes on the patient's body time while its drug-handling genes are not? **CPTAC-3 tumour-normal pairs, patients independent of TCGA; confirmatory for the descriptive result of A249; L1-like association across samples (A229 statistic); secondary (cancer section)**
**Why.** In A249 the registered window test could not be read (the tumour sample was too poor a reader of body time),
but the other direction gave a descriptive result: read against the body time of the patient's normal sample, the
tumour's clock-controlled genes outside the 14-gene panel follow it, clock-specifically, and the tumour's
drug-handling genes and drug targets do not. Both halves were seen before being stated, so neither is a finding.
CPTAC-3 holds more pairs than TCGA from other patients, centres and pipelines; the hypotheses, sets and thresholds
below are fixed from TCGA and tested once.
**Data.** The CPTAC-3 pairs of A186 (`A186_cptac_c1/output/cache`, raw GDC STAR counts, aliquots summed): the 7
types of its pool (LUAD, CCRCC, LSCC, SKCM, HNSCC, PDAC, UCEC; 686 pairs, all gated). STAD (27 pairs) is descriptive.
**States.** z_N, z_T as `a186_lib.states` (GTEx waveform, one projector per type), centred within type.
**Gene sets.** C, K, F and O of A249, unchanged. A gene enters a type with a median count >= 10 in both sample
types; a set enters with >= half of its genes (>= 100 for F) in >= 4 types.
**Statistic.** As A249: per gene and type, r of A229 between y (log2(CPM + 1), centred within type) and a state,
pooled with patient-count weights; S is the mean of r over a set. S_T = S(y of the tumour sample, z_N), S_N = S(y of
the normal sample, z_T). Null: patients of y permuted within type, one permutation for all genes and sets (10,000;
seed 250). Clock specificity: above the q95 of 200 expression-matched pseudo-clock sets (A227 construction). Beyond
typical genes: above the q95 of 200 expression-matched gene sets of the same size. Contrast of two sets:
Delta = S_T(O) - S_T(X), one-sided p from the same permutations (which carry the set-size bias of each S), and a 90 %
interval from 2,000 patient bootstraps within type (states not refitted).
- **H250.1 (primary; replication):** set O, S_T: p < 0.05, clock-specific and beyond typical genes.
- **H250.2 (primary; the asymmetry):** H250.1 passes; set C, S_T is not detected (p >= 0.05) at an MDE below 0.042
  (the excess of set O over its null mean in TCGA); and Delta(O, C) has p < 0.05 with a 90 % interval excluding 0.
- **H250.3 (secondary):** the same asymmetry for X = F and X = K (Holm over the two).
- **H250.4 (secondary; genes):** NFIL3, RORC and ARNTL2 in the tumour sample against z_N, each p < 0.05 (Holm over
  the three; gene-level null from the same permutations).
- **H250.5 (secondary; the window of A249, second attempt in independent pairs):** P249.0 as registered (set O,
  S_N: p < 0.05 and clock-specific); only if it passes, H249.1 and H249.2 as registered in A249. If it fails again,
  the normal half stays unreadable from tumour-sample states and this is stated as such.
**Descriptive (no pass / fail).** Per gene S_T, S_N and implied peak, beside the TCGA values; per type; STAD; the
same-sample readouts, labelled as confounded; S_T of sets O and C by tertile of tumour G2M score (A186 labels) and,
if a purity table is on disk at the power note, by purity tertile; TCGA and CPTAC combined (weighted by pairs).
**Power first (dated note, no observed association):** pairs and types per set; genes eligible; null SD and MDE
(2.49 x SD) of S_T, S_N and of each Delta. Not scored if the MDE of S_T for set O exceeds 0.042 or that of
Delta(O, C) exceeds 0.045 (the TCGA difference): the replication could then not see an effect of the TCGA size.
**Reading.** H250.1 and H250.2 pass: in two independent sets of patients, clock output in the tumour stays on the
patient's body time while the transcripts that handle anticancer drugs do not follow it. This is half of the
chronotherapy premise (the tumour side offers no timed target at the transcript level); the normal side stays open
unless H250.5 reads it. H250.1 passes and set C is detected too: drug handling in tumours is on body time in CPTAC,
the TCGA null was a power or cohort matter. H250.1 fails at an adequate MDE: the TCGA result does not replicate and
is dropped from the text. H250.1 passes without clock specificity: a patient-level state shared by the two samples,
not shown to be the clock (as H186.1).
**Limits.** (1) Bulk: a gene of the tumour sample is read in malignant and other cells; clock output carried by the
stroma of the tumour sample would pass H250.1 (H249.4 could not resolve this: 61 patients). (2) Tumour and normal of
a patient share the operation, and often the processing batch: the pseudo-clock sets are the control. (3) "Not on
body time" for set C is a bound at the MDE, never proof of absence; several set C genes are liver- or gut-specific.
(4) Set O has six genes. (5) Surgery in working hours compresses body time. (6) Transcripts, not activity.
**Stop rule:** one power note, one scoring run. **Code:** `paper/exploratory/A250_cptac_window/scripts/` (on
`a186_lib` and the functions of `compute_a249.py`).

**A250 power note (2026-11-12, before any observed association; `A250_cptac_window/output/power.csv`).**
Null from patients of y permuted within type only (10,000). The 7 CPTAC types of the A186 pool, 686 pairs (45 to 215
per type; states reproduce `a186_lib.states`). Genes scored / eligible in >= 4 types: C 15 / 14, K 9 / 9, O 6 / 6, F 784 /
687 of 854; all sets enter. Null mean of S about 0.086. Null SD and MDE of S_T: O 0.0084, MDE 0.021 (limit 0.042); C
0.0059, MDE 0.015 (below 0.042, as H250.2 requires); K MDE 0.018; F MDE 0.009. MDE of S_N: O 0.024, C 0.019. Contrasts:
Delta(O, C) null SD 0.0092, MDE 0.023 (limit 0.045); Delta(O, F) MDE 0.019; Delta(O, K) MDE 0.024. Both power rules are
met: an effect of half the TCGA size is detectable. No purity table of CPTAC is on disk, so the purity strata of the
descriptive part are omitted. Code run only on shuffled patients (`--placebo --smoke`). Choices not fixed by the
registration are in the script docstring. Scoring proceeds as registered.

## A251 (2026-11-12, registered under the user's instruction of 2026-11-12 to "work autonomously ... register the tests, download data" on the directions of REVIEW_2026-11-10 §4 (D4) and Q29 (b); staged: compartments, statistics, hypotheses, families and readings are fixed here, gates and power go in a dated note before any observed pairing or offset; seen: the scout `scouts_2026-11-12/brain_skin_scout.md` (sample sheet, batches, probe sets of the clock genes and their median intensity percentile); the 2026-11-08 metadata check (no age in GEO); NOT seen: any clock statistic, any expression value compared between compartments or persons; Wu 2020's own results are known from the card `wu2020population` (a population phase predictor; no within-person pairing is reported)): do skin compartments of distant body sites of a living person share one clock state? **Wu 2020 (del Olmo), GSE139301 + GSE139300, 154 women, 3 sites x 2 layers, untimed clinic biopsies, HG-U219 arrays; L2 pairwise; secondary by registration (robustness arm, living cohort)**
**Why.** The paper's statements on a shared internal time rest on post-mortem GTEx (A109, A150). Review item 2
("missing for Cell"): one living cohort. Here each woman gives epidermis and dermis from forearm, buttock and cheek at
one visit, without a recorded time: the design the method is built for. Each compartment was hybridised in its own
batch, so a pairing across compartments cannot come from a shared array batch.
**Data.** Women: 154 (137 with all six compartments). Timed arm: 20 men, forearm, both layers, at 0, 6, 12 and 18 h.
Probe set of a gene, per layer series: among its probe sets the one with the highest mean over the women's samples,
restricted to those ending in `_at` without `_s_` / `_x_` when any exists. Age is not public: no age contrast.
**States.** Per sample, the clock-plane state of the 14 genes (GTEx waveform, arrays: unweighted projector,
`clocksync.paired.states`), centred within compartment (women) or within compartment and man (timed arm).
**Statistic.** Pair of compartments a, b: r = |sum z_a conj(z_b)| / sqrt(sum |z_a|^2 sum |z_b|^2) over the women with
both (A227 `raw_stat`); a family of pairs is pooled with weights = women. Null: women permuted independently within
each compartment (10,000; seed 251). Companions on the same statistic: 200 expression-matched pseudo-clock sets and 200
responsive sets (`pseudo_waveforms`, means from the layer's genome-wide probe sets collapsed by the same rule);
"clock-specific" only above the q95 of both. Offset of a pair: arg of the cross-moment (h), 90 % interval from 2,000
bootstraps of women.
- **P251.0 (positive control, timed men):** in each layer, the state follows the collection time: modulus of the
  complex cosinor of the state on clock time (states centred within man), against times permuted within man (10,000);
  p < 0.05 in both layers. If it fails in both, the platform does not read the clock and everything below is
  uninformative.
- **H251.1 (primary):** the same layer at different body sites shares a clock state: the six pairs (three site pairs
  x two layers) pooled, p < 0.05 and clock-specific.
- **H251.2 (secondary):** epidermis and dermis share it: (a) at the same site (three pairs), (b) at different sites
  (six pairs); each p < 0.05 and clock-specific (Holm over the two).
- **H251.3 (secondary; one time across the body surface):** every same-layer site-pair offset has its 90 % interval
  within +/- 1.5 h (equivalence; passes only if all six do).
- **H251.4 (secondary; offset against timed truth):** the dermis-minus-epidermis offset of the women's forearm lies
  within the 90 % bootstrap interval of the same offset in the timed men (difference of the two cosinor acrophases,
  men resampled), or the two differ by less than 1 h.
**Descriptive.** Each of the 15 pairs (r, companions, offset); clock-plane variance per compartment; without PER1
and NR1D1; the reference gate of each compartment (`reference_gate`); cheek (sun-exposed) against the two other sites.
**Power first (dated note, no observed pairing):** women per pair; gate per compartment; null SD and MDE (2.49 x SD)
of each pooled r; bootstrap SE of offsets is part of scoring. Not scored if fewer than four compartments pass the gate
or the MDE of H251.1 exceeds 0.10.
**Reading.** H251.1 passes clock-specifically: in living people, skin at sites a metre apart is on one internal time,
the first living-cohort support for the shared-time statement, from clinic-hours biopsies. Passes but not above the
companions: a person-level state shared across sites (age is a candidate: 20-74 y, not available), not shown to be the
clock. Fails at the MDE: site clocks are not detectably shared here, which would restrict the GTEx statement to
post-mortem tissue. H251.3: offsets within 1.5 h support "one time"; H251.4 ties the untimed offset to a timed one.
**Limits.** (1) One cohort, women only, arrays. (2) The spread of collection times is unknown; a narrow window lowers
r but does not create it. (3) Age and other traits vary between women and could move clock genes: the companions are
the control, and age cannot be adjusted. (4) Epidermis and dermis come from one biopsy (laser microdissection): a
cross-layer pairing at the same site shares the biopsy; cross-site pairings do not. (5) Layers were normalised
separately; statistics use centred states only.
**Stop rule:** one power note, one scoring run. **Code:** `paper/exploratory/A251_skin_sites/scripts/`.

**A251 power note (2026-11-12, before any observed pairing or offset of the women; `A251_skin_sites/output/{power,gates}.csv`).**
All six compartments pass the reference gate (clock share 0.59-0.76 against a null mean of 0.09-0.12; p 0.005, the floor
of 200 permutations). Women per pair 146 to 150. Null (women permuted within compartment, 10,000): H251.1 null mean 0.072,
SD 0.017, MDE 0.041 (limit 0.10); H251.2a SD 0.023, MDE 0.058; H251.2b SD 0.017, MDE 0.041. Scored. **Disclosure:** the code
test (`--placebo --smoke`, women shuffled within compartment) does not shuffle the timed men, so P251.0 was seen there
before this note: modulus 1.62 (epidermis) and 0.99 (dermis) against null q95 0.46 and 0.31, p at the floor of 200
permutations; men's dermis-minus-epidermis offset -2.0 h (-3.2 to -1.0). The positive control therefore passes on seen
data; nothing of the women was seen (their statistics in that run are on shuffled women). Choices not fixed by the
registration are in the script docstring. Scoring proceeds as registered.

## A252 (2026-11-12, registered under the user's instruction of 2026-11-12 ("1 looks great ... work autonomously ... register the tests, download data"); staged: datasets, assumed clock times, statistic, hypotheses, families and readings are fixed here, eligibility and power go in a dated note before any observed change of a clock-gene state; seen: the scout `scouts_2026-11-12/repeat_biopsy_scout.md` (sample sheets, time points, presence and median percentile of the clock genes); NOT seen: any expression value compared between time points; none of the datasets below was used before in this project or its archive (GSE106800, used in A124, is excluded)): does the clock state of a living person's tissue advance by the elapsed time between two biopsies of the same day? **Ten repeat-biopsy studies (adipose, muscle; atrium and liver graft descriptive); L1 state change within person; secondary by registration (method validation in the living, Fig. 1)**
**Why.** Every human validation of the readout so far is across persons (GTEx deaths, one breast cohort, the timed men
of A251). Clinical studies that take a biopsy, intervene and take another hours later give a within-person test: the
same tissue of the same person, so gene-specific levels cancel, and the state must move along the clock's circle by
15 degrees per hour. Deviations are of interest (exercise, meals, surgery), but the first question is whether the
rotation is seen at all.
**Data and clock times (fixed now).** Elapsed hours as in the scout sheets. Clock time of the first biopsy: stated for
GSE118280 (adipose, 08:40, then 12:20 and 19:00; 15 men, two diet visits); assumed 09:30 for GSE231509 (muscle,
arrival 09:00 fasted; +1 h after a meal) and 09:00 for the other morning fasted protocols: GSE252357 (muscle; rest
arm and resistance-exercise arm; +1.75, 4.25, 9.25 h; the 24 h point is not used), GSE142401 (adipose, +4 h after a
meal), GSE151066 (muscle, +0.67, 3.67 h, cycling), GSE273041 (muscle, +3, 5 h, 3 h run), GSE303983 (adipose, +2.25 h:
exercise session plus 1.5 h), GSE120862 (muscle, +2, 5 h, exercised and resting leg). Descriptive only (no clock
time, no elapsed time, or not the same body): GSE4386 (atrium, start and end of bypass surgery), GSE75432 (muscle,
before and >= 2 h after knee surgery), GSE151648 (liver graft before and after transplantation).
**State.** Per sample z = L y, with y the log expression of the clock genes present (arrays: log2 x ln 2; counts:
ln(count + 1)) and L the unweighted projector of the GTEx waveform (`clocksync` `projector`; L 1 = 0, so depth and
array scale cancel). Unit = person x visit (or x leg).
**Statistic.** For a unit and a later time point k, delta = z_k - z_0 and the predicted chord
c_k = exp(i w (t0 + D_k)) - exp(i w t0) (w = 2 pi / 24 h, t0 clock time of the first biopsy, D_k elapsed hours).
rho = sum Re(delta conj(c)) / sum |c|^2 over units and time points of a dataset, and with pooled sums over datasets:
the amplitude of the rotation in units of the GTEx waveform (1 = the GTEx amplitude, 0 = no rotation; negative =
moving backwards). Null: time-point labels permuted within unit (10,000; all assignments when fewer). Companions: rho
from 200 expression-matched pseudo-clock sets and 200 responsive sets carrying the waveform's phases (per dataset, from
its genome-wide means); "clock-specific" only above the q95 of both.
- **H252.1 (primary):** GSE118280 (exact clock times): rho > 0, p < 0.05, clock-specific.
- **H252.2 (secondary):** the seven other datasets pooled (assumed clock times): rho > 0, p < 0.05, clock-specific.
- **H252.3 (secondary):** by tissue among those seven: muscle (5 datasets) and adipose (2), Holm over the two.
- **H252.4 (secondary; without the genes a meal or stress moves directly):** H252.1 and H252.2 with PER1 and NR1D1
  removed from the waveform: rho > 0, p < 0.05 (Holm over the two).
**Descriptive (no pass / fail).** rho per dataset and time point; the implied clock time of the interval's midpoint
(direction of the mean delta) against the protocol's; rho against elapsed time; t0 moved by +/- 2 h; rest against
exercise (GSE252357 arms, GSE120862 legs): difference of the mean delta at matched times with a permutation p (arms
are 5 vs 8 persons and 7 persons: no claim); the three descriptive datasets: modulus and direction of the mean delta
against the companions.
**Power first (dated note, no observed change):** units and time points per dataset; clock genes present; null SD and
MDE (2.49 x SD) of rho for H252.1-3. Not scored if the MDE exceeds 0.7 (a rotation with the GTEx amplitude must be
detectable with margin).
**Reading.** H252.1 passes clock-specifically: within a living person the readout advances with the day, the first
within-person validation, and rho is the tissue's amplitude relative to GTEx. H252.2 passes too: it holds across
tissues and laboratories with protocol times only. A pass that the responsive sets match: a state change shared with
stress or feeding genes, not shown to be the clock's rotation. A fail at the MDE: the rotation is not seen within
persons over 1 to 10 h, which would question state changes read as time. H252.4 says whether PER1 and NR1D1 carry it.
**Limits.** (1) Every interval contains an intervention or a meal; only the rest arm of GSE252357 (5 persons) has
none, so a pass is "elapsed time plus what happened in it". (2) Assumed clock times (H252.2-3) are protocol guesses;
an error of 2 h turns the prediction by 30 degrees (cos 0.87). (3) A second biopsy of the same muscle or fat depot
carries a wound response from the first; the responsive companion is the control. (4) rho mixes amplitude and
rotation: a full-amplitude clock that stalls and a damped clock that turns give the same value. (5) Arrays and
normalised counts are used as deposited.
**Stop rule:** one power note, one scoring run. **Code:** `paper/exploratory/A252_elapsed_time/scripts/`.

**A252 power note (2026-11-12, before any observed change of a clock-gene state; `A252_elapsed_time/output/{power,design}.csv`).**
Units (person x visit or leg) and later time points: GSE118280 30 units (15 men x 2 visits), +3.67 and +10.33 h; GSE231509 21,
+1 h; GSE252357 13, +1.75 / 4.25 / 9.25 h; GSE142401 38, +4 h; GSE151066 19, +0.67 / 3.67 h; GSE273041 20, +3 / 5 h; GSE303983 45,
+2.25 h; GSE120862 28 (7 persons x 2 legs x 2 training states), +2 / 5 h. All 14 clock genes present (13 in GSE4386). Null of
rho (time labels permuted within unit, 10,000): H252.1 SD 0.215, MDE 0.53; H252.2 SD 0.237, MDE 0.59; H252.3 adipose SD 0.276,
MDE 0.69; H252.3 muscle SD 0.304, MDE 0.76: above 0.7, **not scored** (its value is printed by the script and is reported
as descriptive only; the Holm family of H252.3 reduces to adipose). Deviation: GSE4386 is deposited on the linear scale
and is log2-transformed (descriptive dataset). Code run only with time labels shuffled within unit (`--placebo --smoke`).
Choices not fixed by the registration are in the script docstring. Scoring proceeds as registered.

## A253 (2026-11-12, registered under the user's instructions of 2026-11-12 ("Draft the amendment for depression"; then "you can register and run the depression amendment if you believe it has a good chance of passing, otherwise don't bother"); the lead's judgement: the cohorts found hold about 67 depressed and 65 control donors with three to six regions each, enough to see a clear loss of region relations if there is one, and both outcomes are usable (REVIEW_2026-11-10 D2 b); staged: cohorts, statistics, hypotheses, families and readings are fixed here, the precondition and power go in a dated note before any case-control value; seen: the scout `scouts_2026-11-12/psych_transplant_scout.md` (sample sheets, donors per group and region, presence and median percentile of the clock genes); the card `li2013circadian` (weaker rhythms in depression by time of death; the patients of that study are not public); WP0 item 8 and A141 (brain disease contrasts depend on the host anchor; the statistics below use no host); NOT seen: any clock statistic of these cohorts): are the clock relations between brain regions of one person changed in major depression or psychosis? **GSE80655, GSE102556, GSE53987 (depression); GSE80655, GSE202537, GSE53987 (schizophrenia and bipolar disorder); post-mortem, 3 to 6 regions per donor; L2 / L3 pairwise, anchor-free; secondary by registration (candidate human desynchrony test, Q29 a)**
**Why.** Depression is the human disease for which internal desynchrony is most often proposed, and the evidence is
weaker rhythms against time of death (Li 2013), which cannot tell a shifted patient from regions out of step with
each other. Region pairs of one donor can: both regions were taken at the same death, so their relation needs no
time. By the ledger's definition (rule 3) internal desynchrony is a change of the phase offset between compartments;
a lower pairing is loss of shared clock structure.
**Data.** GSE80655 (counts; anterior cingulate, DLPFC, nucleus accumbens; 24 donors each of control, depression,
bipolar disorder, schizophrenia); GSE102556 (log FPKM; six regions; 22 controls, 26 depression; the 18 BA25 libraries
without a GEO row take the diagnosis of their donor); GSE53987 (arrays; BA46, hippocampus, striatum; donors linked by
identical age, sex, PMI and pH); GSE202537 (counts; caudate, putamen, nucleus accumbens; 36 controls, 28
schizophrenia, 8 bipolar psychosis; time of death recorded). "Psychosis" = schizophrenia and bipolar disorder pooled.
**States.** Per sample the clock-plane state of the genes present (GTEx waveform, unweighted projector as A252),
centred within cohort x region x group.
**Statistics (per cohort, group and region pair; donors with both regions).** Pairing r (A227 `raw_stat`) and offset
(arg of the cross-moment, h). R_g: mean of r over the cohort's region pairs (weights = donors).
D = R_control - R_case and A = mean over region pairs of |offset_case - offset_control| (wrapped; weights = case
donors), each pooled over the family's cohorts with weights = case donors. Null: diagnosis labels permuted among the
donors of a cohort (cases of the family and controls only; 10,000; seed 253). Companions: D and A from 200
expression-matched pseudo-clock sets and 200 responsive sets; "clock-specific" only above the q95 of both. Interval
of D: 90 %, 2,000 bootstraps of donors within cohort x group.
- **P253.0 (precondition):** the controls' regions share a clock state: R_control pooled over the four cohorts, p <
  0.05 against donors permuted within region, and clock-specific. If it fails, everything below is uninformative.
- **H253.1 (primary; depression, offsets):** A > its null, p < 0.05, clock-specific.
- **H253.2 (primary; depression, pairing):** D > 0, p < 0.05, clock-specific. Holm over H253.1-2.
- **H253.3-4 (secondary; psychosis):** the same two tests, Holm over the two.
- **Equivalence reading of D** (if H253.2 or H253.4 fails): "kept" if the 90 % interval of D lies below 0.15.
**Descriptive.** Per cohort and per diagnosis (schizophrenia and bipolar apart); r and offset per region pair and
group; clock-plane variance per group (a lower spread lowers r without any desynchrony: flagged if the log ratio
exceeds 0.3 in absolute value); GSE202537: modulus and acrophase of the complex cosinor of the state on time of death
per group and region; D by suicide (GSE102556) if both arms hold >= 10 cases.
**Power first (dated note, no case-control value):** donors per group, cohort and pair; P253.0; null SD and MDE (2.49
x SD) of D and of A per family. Not scored if P253.0 fails or the MDE of D for depression exceeds 0.20.
**Reading.** H253.1 passes clock-specifically: regions of depressed (or psychotic) donors keep other phase relations
than controls': internal desynchrony in a human disease, the positive the paper lacks, to be replicated. H253.2
alone: the regions share less clock structure, without a shift between them (noise or damping; the variance
descriptive says which is more likely). Both fail with D "kept": region relations are preserved in depression within
the bound, so the weaker rhythms against time of death are a shift of the whole brain against the day or of the
patients' time of death, which fits the paper's thesis. Fails without the bound: not resolved.
**Limits.** (1) Diagnosis is confounded with suicide, medication, brain pH and agonal state; pH and PMI are
descriptive covariates only. (2) Death at any hour gives the spread; groups may differ in their hours of death. (3)
Three to six regions, none of them the SCN. (4) GSE53987 donors are linked by inference; GSE102556 is FPKM. (5) Pooled
"psychosis" mixes two diagnoses. (6) Lower RNA quality in cases lowers r: the companions are the control.
**Stop rule:** one power note, one scoring run. **Code:** `paper/exploratory/A253_psych_regions/scripts/`.

**A253 power note (2026-11-12; `A253_psych_regions/output/{design,precondition,power}.csv`): NOT SCORED, precondition fails.**
Donors: depression 68 cases and 65 controls in three cohorts; psychosis 122 and 79. Power would have sufficed (label
permutation, 10,000): D null SD 0.054, MDE 0.134 (depression; limit 0.20) and 0.168 (psychosis); A null SD 0.14 and
0.17 h (MDE 0.34 and 0.41 h). **P253.0 fails:** the controls' regions share a state far above the donor-permutation
null (R_control 0.734, null mean 0.182, SD 0.026, p 0.0001) and above every expression-matched set (median 0.387, q95
0.549), but not above the stress-responsive sets (median 0.648, q95 0.758; p_responsive 0.139). By the registered
rule the case-control tests are not scored: a pairing that responsive genes match cannot be read as a clock relation.
**Disclosures.** (1) The code test (`--placebo --smoke`) shuffled the labels of the tests but its descriptive block
used the true groups, so the pairing R per true group and cohort (not the offsets) was printed and seen before this
note, with GSE53987 still on the wrong (linear) scale; no decision above depends on it, and the precondition uses
controls only. (2) GSE53987 is deposited on the linear scale and is now log2-transformed (deviation). The scoring
stage was not run.

## A254 (2026-11-12, registered under the user's instruction of 2026-11-12 ("3 seems great ... work autonomously ... register the tests, download data"); staged: datasets, estimator, hypotheses, families and readings are fixed here, gates and power go in a dated note before any lesional-minus-non-lesional value; seen: the scout `scouts_2026-11-12/lesion_multitissue_scout.md` (sample sheets, pair counts, presence and median percentile of the clock genes); every psoriasis and eczema result (A54-A59, A115, A183: lesion delayed by about 1.4 h, at the specificity line, not specific against eczema), A230 (inflamed gut: not detected); NOT seen: any expression value of these datasets compared between lesion and non-lesion): is a clock delay of the lesion against the same person's unaffected tissue a general property of chronic lesions outside the skin? **Periodontitis, osteoarthritis, carotid plaque, uterine fibroid, hidradenitis; paired samples of one procedure; L3 (offset) and clock-plane variance; secondary by registration (claim "where clocks break, they break locally")**
**Why.** The local-break claim rests on psoriasis, where the lesion runs about 1.4 h behind the person's own
unaffected skin. If chronic lesions of other tissues show a shift against their own tissue, the claim is general; if
none does at an adequate bound, the psoriatic delay is a property of that skin disease.
**Data (pairs = persons with a lesional and a non-lesional sample of the same tissue from one procedure).**
GSE16134 (periodontitis-affected vs healthy gingival sites of one mouth, 69 persons; several affected sites of a
person are averaged on the log scale); GSE57218 (osteoarthritic vs preserved cartilage of one joint, 33); GSE43292
(carotid plaque vs adjacent intact artery, 32); GSE192352 (leiomyoma vs own myometrium, 30; counts; a benign tumour,
kept as its own row); GSE151243 (hidradenitis lesion vs perilesional skin, 20; counts). Descriptive only: GSE10334
(periodontitis, 64 persons; overlap with GSE16134 undocumented), GSE39491 (Barrett mucosa vs the same patient's
squamous oesophagus, 40: a different epithelium). Excluded: GSE26051 (different tendons), GSE119717 (fat depots, two
clock genes missing).
**Estimator.** `clocksync.paired` as A183: both samples of a dataset on one projector (GTEx waveform; counts or
arrays), states centred within sample type; offset = arg of the complex correlation of the paired states (h;
negative = lesion delayed); r = log ratio of the clock-plane variance, lesion over non-lesion. Null: the two labels
swapped per person (10,000; seed 254). Companions: the same from 200 expression-matched and 200 responsive
pseudo-clock sets per dataset.
- **Gate:** a dataset enters with >= 20 pairs and a reference gate of its non-lesional samples at p < 0.05.
- **H254.1 (primary; the psoriasis direction, pooled):** the pooled offset of the gated datasets (weighted resultant,
  weights = pairs) is negative: one-sided swap p < 0.05, and below the q05 of the pooled offsets of both companions.
- **H254.2 (secondary; each disease):** offset different from 0, two-sided swap p, Holm over the gated datasets;
  "clock-specific" if |offset| exceeds the q95 of |offset| of both companions.
- **H254.3 (secondary; clock signal):** r different from 0 per dataset (two-sided swap p, Holm) and pooled
  (inverse-variance); "clock-specific" if outside the central 90 % of both companions.
**Descriptive.** Coherence of the pair, bootstrap interval of each offset, without PER1 and NR1D1, GSE10334 and
GSE39491 (both control tissues) with the same readouts.
**Power first (dated note, no lesional-minus-non-lesional value):** pairs, gate, coherence flag (`ccorr`), swap-null SD
and MDE (2.49 x SD) of each offset and of the pooled one. H254.1 is not scored if the pooled MDE exceeds 1.0 h or
fewer than three datasets pass the gate; a dataset's own test is not scored with an MDE above 2.5 h.
**Reading.** H254.1 passes clock-specifically: lesions of several tissues run late against the person's own tissue,
and the local break is a general feature of chronic lesions. Fails with a pooled MDE below 1 h: no psoriasis-sized
delay outside skin; the claim stays a skin statement. H254.2 can single out one disease in either direction. H254.3:
lesions that lose or gain clock-plane variance (as tumours lose it).
**Limits.** (1) Offsets are read from the covariance over persons, so they need a spread of procedure times; the
gate tests that the non-lesional tissue shows it. (2) Lesion and normal tissue differ in cell composition; a
composition axis that is clock-shaped would mimic an offset; the companions are the control. (3) Arrays, small n per
disease. (4) "Non-lesional" tissue of a diseased organ is not healthy tissue. (5) No disease here is psoriasis-like
in mechanism; a null does not bear on psoriasis itself.
**Stop rule:** one power note, one scoring run. **Code:** `paper/exploratory/A254_paired_lesions/scripts/`.

## A255 (2026-11-12, registered under the same instruction ("4 seems great"); staged: cohorts, statistics, hypotheses and readings are fixed here, gates and power go in a dated note before any observed pairing; seen: the scout `scouts_2026-11-12/lesion_multitissue_scout.md` (sample sheets; no metabolic phenotype is public in any multi-tissue living cohort, and the bariatric cohort of Greenawalt has no person key, so the obesity / diabetes contrast of A164 / A235 cannot be repeated in the living); A251 (skin compartments of living women share a clock state); NOT seen: any expression value of these cohorts): do two different tissues of a living person share a clock state, and is it more than what the person's genome alone would give? **MuTHER twins (E-TABM-1140: subcutaneous fat, skin, lymphoblastoid cell line; 856 women) primary; three surgical two-tissue cohorts secondary; L2 pairwise; secondary by registration (robustness arm, living cohorts)**
**Why.** A251 shows a shared state within one organ (skin sites and layers). The paper's statement is about organs.
MuTHER has fat and skin of 681 living women, and a built-in control no other cohort offers: a lymphoblastoid cell
line of the same woman, grown in culture, carries her genome but not her body time. Clock genes have cis-regulatory
variants that act across tissues, so a pairing of two tissues could be genetic; the cell line measures that part.
**Data.** MuTHER processed matrices (log2, one column per woman and tissue). Surgical cohorts: GSE73439 (subcutaneous
and omental fat, 51 women, caesarean or gynaecological surgery), GSE133979 (rectus muscle and subcutaneous fat, 33
patients, abdominal surgery), GSE24883 (subcutaneous and visceral fat, 32 women; two-colour ratios; pairing by title
index, unverified).
**States and statistic.** As A251: clock-plane state per sample (GTEx waveform, unweighted projector), centred within
tissue (and cohort); pairing r of two tissues over the persons with both; null: persons permuted within tissue
(10,000; seed 255); companions: 200 expression-matched and 200 responsive pseudo-clock sets; "clock-specific" only
above the q95 of both.
- **Gate:** reference gate of fat and of skin (p < 0.05).
- **H255.1 (primary):** fat and skin of a woman share a clock state: r > null, p < 0.05, clock-specific.
- **H255.2 (primary; beyond the genome):** r(fat, skin) exceeds the mean of r(fat, cell line) and r(skin, cell
  line): difference > 0 with a 90 % bootstrap interval (2,000, women) excluding 0. Holm over H255.1-2 for the p of
  H255.1 only (H255.2 is an interval).
- **H255.3 (secondary; surgical cohorts):** the three cohorts pooled (weights = persons): r > null, p < 0.05,
  clock-specific.
- **H255.4 (secondary; co-twins):** in fat and in skin, the states of co-twins pair (r over twin pairs, the order
  within a pair randomised; null: co-twins permuted) more than in the cell line (difference with a bootstrap interval
  over twin pairs), Holm over the two tissues. Twins usually attend the clinic together, so this reads shared visit
  time plus shared genome; the cell line reads the genome.
**Descriptive.** r(tissue, cell line) with companions; co-twin pairing by zygosity (monozygotic, dizygotic) per
tissue; fat-skin offset with a bootstrap interval; without PER1 and NR1D1; each surgical cohort.
**Power first (dated note, no observed pairing):** persons per pair, gates, null SD and MDE of each r. Not scored if
a gate fails (that tissue) or the MDE of H255.1 exceeds 0.10.
**Reading.** H255.1 and H255.2 pass: two organs of a living person share a clock state that her cell line does not
carry: body time, not genotype. H255.1 passes and H255.2 fails: the shared state is at least partly the person's
constitutive clock-gene levels, and the GTEx statement inherits that caveat. H255.4: co-twins sharing a visit share
body time across persons, a between-person positive control.
**Limits.** (1) Fat and skin come from one punch biopsy; time to freezing is shared (the responsive companion is the
control). (2) Women only, arrays. (3) The cell line differs from tissue in more than time (culture, EBV); it bounds
the genetic part only for genes expressed in it. (4) Whether co-twins were sampled together is not recorded per pair.
**Stop rule:** one power note, one scoring run. **Code:** `paper/exploratory/A255_living_tissues/scripts/`.

**A255 power note (2026-11-12, before any observed pairing; `A255_living_tissues/output/{gates,power}.csv`).**
MuTHER: 680 women with fat and skin, 794 with fat and cell line, 679 with skin and cell line; 14 clock genes. Gates:
fat passes (clock share 0.348 against 0.102, p 0.010), skin passes (0.205 against 0.066, p 0.015); the cell line shows no
clock-shaped covariance (0.042 against 0.065, p 0.75), as expected of cultured cells. Null (women permuted, 10,000): fat-skin
null mean 0.034, SD 0.018, MDE 0.044 (limit 0.10); fat-cell line MDE 0.042; skin-cell line MDE 0.045. Surgical cohorts:
GSE73439 (51 women, both gates p 0.005, MDE 0.16) and GSE133979 (32 patients, 13 clock genes, both gates p 0.005, MDE
0.21) enter; GSE24883 fails its gates (p 0.78 and 0.13; 8 clock genes left after missing probes; its person link was
also unverified) and leaves H255.3, which pools 83 persons (null SD 0.052, MDE 0.129). Co-twin pairs are counted at
scoring. Code run only with persons shuffled within tissue (`--placebo --smoke`). Choices not fixed by the registration are
in the script docstring. Scoring proceeds as registered.

**A254 power note (2026-11-12, before any lesional-minus-non-lesional value; `A254_paired_lesions/output/power.csv`).**
All five datasets pass the gate (p 0.005 to 0.015; clock share 0.32 to 0.83 in non-lesional tissue) with unflagged
coherence: periodontitis 69 pairs (13 clock genes), osteoarthritis 33, carotid plaque 32, uterine fibroid 30,
hidradenitis 20. Label-swap null (10,000) of the offset: SD 0.33 h (MDE 0.83 h), 0.66 (1.65), 0.86 (2.15), 1.03 (2.56), 0.66
(1.64). The fibroid's own test exceeds 2.5 h and is **not scored** (it stays in the pooled test). Pooled offset (184
pairs): null SD 0.275 h, MDE 0.69 h (limit 1.0): H254.1 is scored; a psoriasis-sized delay (1.4 h) is well inside reach.
Person links were checked with sex genes at preparation (concordance 0.94 to 1.00 where informative). Code run only
with labels swapped at random per person (`--placebo --smoke`). Choices not fixed by the registration are in the script
docstring. Scoring proceeds as registered.

## A256 (2026-11-12, registered under the user's instruction of 2026-11-12 (autonomous work on direction 1); follow-up of A252, whose results are all seen (rho 1.27 with exact times, 2.97 in seven studies, direction within 1 h at 13 of 14 points; rest arm of GSE252357 rho 2.89 on 5 persons); staged: statistic, hypotheses and readings are fixed here, the dataset list with arms, elapsed hours and assumed clock times goes in a design addendum committed before any expression value of these datasets is read; NOT seen: any expression value of the candidate datasets (rows 16 to 25 of `scouts_2026-11-12/repeat_biopsy_scout.md`: GSE40589, GSE9105, GSE7146, GSE59088, GSE220899, GSE202295, GSE164081, GSE163356, GSE43219, GSE126296; none downloaded yet)): does the within-person rotation of the clock state replicate in fresh studies, and is it seen when nothing is done between the biopsies? **Fresh repeat-biopsy studies; L1 state change within person; secondary by registration (replication of A252)**
**Why.** A252 passed, but (i) seven of its eight studies were scored with assumed clock times after one registration,
and (ii) every interval but one small arm (5 persons) contained a meal, exercise or surgery, so "elapsed time" and
"what happened in it" are not separated. Infusion studies have saline arms in which fasted persons lie still between
two biopsies: elapsed time alone.
**Statistic.** rho of A252, unchanged (same projector, chords, permutation of time labels within unit, matched and
responsive companions, seed 256). Clock time of the first biopsy: as stated by the study, otherwise 09:00.
- **H256.1 (primary; nothing done):** units of arms without exercise, meal, insulin, glucose or drug between the
  biopsies (saline infusion or rest), pooled: rho > 0, p < 0.05, clock-specific.
- **H256.2 (primary; replication):** all fresh units pooled: rho > 0, p < 0.05, clock-specific. Holm over H256.1-2.
- **H256.3 (secondary; direction):** over all dataset x time-point cells, the implied clock time of the interval's
  midpoint lies within 2 h of the protocol's in more cells than under time labels permuted within unit (count as
  statistic, p < 0.05).
- **H256.4 (secondary; does an intervention add to it):** rho of the intervention arms minus rho of the
  nothing-done arms in studies that have both, arm labels permuted over persons within study; two-sided.
**Descriptive.** rho per dataset, arm and time point; without PER1 and NR1D1; rho against elapsed hours.
**Design addendum first (dated, no expression value read):** datasets retained (>= 8 persons with a first and a later
biopsy on the same day, per-sample matrix public, clock genes present), arms and their class (nothing done /
intervention), elapsed hours, clock times; then the power note (null SD and MDE of rho; H256.1 is not scored with
fewer than 10 units or an MDE above 1.0, H256.2 with an MDE above 0.7).
**Reading.** H256.1 passes: the state advances within a person when nothing is done: elapsed time itself. H256.2
passes: A252 replicates on studies chosen before their data were opened. H256.1 fails at an adequate MDE while H256.2
passes: the change needs the intervention (meal or activity as the driver, not time alone), which would change the
reading of A252. H256.4 says whether interventions add to or subtract from the rotation.
**Limits.** As A252 (protocol clock times, amplitude and rotation mixed, repeated biopsy); saline arms are fasted,
so "nothing done" still includes continued fasting.
**Stop rule:** one design addendum, one power note, one scoring run. **Code:** `paper/exploratory/A256_elapsed_replication/scripts/`.

**A256 design addendum (2026-11-13, before any expression value of these datasets is read; table and protocol quotes in `A256_elapsed_replication/design_addendum_draft.md`, prepared by a Sonnet agent with `prep_a256.py`, no clock statistic).**
Retained (7): GSE40589 (subcutaneous fat, 8 persons, six crossover arms, 0 and +4 h, morning fasted; two-colour log
ratios, used as deposited: the base of the logarithm scales rho but no p), GSE9105 (muscle, insulin clamp, +1.5 and +5 h;
12 persons; PER1 and CIART not on the array), GSE59088 (muscle; the same-day pair is 13:00 and 15:30, +2.5 h, clock time
stated; rest-only control and two arms that exercised before the first biopsy; 6 persons each), GSE220899 (muscle,
resistance exercise, +4 and +7 h, two bouts per person), GSE202295 (muscle, first biopsy right after cycling, +3 h of
rest; 35 persons), GSE163356 (muscle, interval exercise, +4 h; two sessions), GSE126296 (muscle, sprints, +3.67 h; 14
persons). Rejected: GSE7146 (6 persons per protocol and platform), GSE164081 (no per-sample matrix), GSE43219 (rest
biopsy on another day). Clock time of the first biopsy 09:00 by the registered default except GSE59088 (13:00).
**Arm classes (decided here).** *Nothing done* (H256.1): saline infusion alone (GSE40589 INF_NACL, 7 units) and rest
alone (GSE59088 control, 6 units): 13 units. *Recovery*: exercise ended before the first biopsy and only fasting rest
lies between the biopsies (GSE59088 endurance and resistance, GSE202295): by the letter of the registration these would
count as nothing done; they are kept out of H256.1 because the earlier exercise may still act, and H256.1 with them
included is reported as a sensitivity (descriptive). *Intervention*: everything else. H256.2 pools all units. H256.4
uses GSE40589 (five clamp or GIP arms against saline alone; the same persons, so arm labels are permuted within
person) and GSE59088 (recovery arms against control; labels permuted over persons), pooled. Person links: verified by
sex genes in GSE9105 and GSE126296; by title or characteristics in the male-only studies; in GSE40589 only by the
"replicate" number of the sample description (unverified).

**A256 power note (2026-11-13, before any observed state change; `A256_elapsed_replication/output/{design,power}.csv`).**
166 units in 17 arms of seven studies. H256.2 (all fresh units): null SD of rho 0.239, MDE 0.60 (limit 0.7): scored.
H256.1 (nothing done: 13 units, saline alone 7 and rest alone 6): null SD 1.03, MDE 2.57, above the limit of 1.0: **not
scored**; its value is printed by the script and reported as descriptive only, and the Holm family reduces to H256.2.
The question "does the state advance when nothing is done" therefore stays open in public data; H256.4 (interventions
against nothing done) is reported with its permutation p but rests on the same 13 units. Code run only with time
labels shuffled within unit (`--placebo --smoke`). Choices not fixed by the registration are in the script docstring.

## A257 (2026-11-13, user: "Go for 1, 2 and 3" (item 2 of REVIEW_2026-11-10 §6.4); **NOT BLIND**: the same data were scored in A124, whose descriptives are seen: the core state moves between 07:00 and 19:00 in the same direction in both conditions (cosine 0.91) with about half the amplitude under misalignment (3.20 -> 1.67); the within-person statistic of A252 / A256 is seen to work in 16 other studies): after a 12 h shift of behaviour, does the human muscle clock stay on central time or follow behaviour? **GSE106800 (Wefers 2018), 12 men, vastus lateralis at 07:00 and 19:00 after an aligned and after a misaligned protocol; L1 within person; secondary by registration; formalises a seen descriptive with intervals and companions**
**Why.** The only public human misalignment design in a solid tissue. A124 asked about a between-gene-set offset
(MDE 7.3 h) and could not resolve it. The within-person change of the clock state between the two biopsy times of a
condition has a predicted direction: if the muscle clock follows behaviour shifted by s hours, the change from 07:00
to 19:00 turns by s hours (s = 12 h: it reverses).
**Data and states.** Series matrix and probe map of A124 (main probe sets, mean log2 per gene, canonical names with
the eight platform aliases of the A124 sensitivity); 11 subjects with all four arrays (RAS9 enters H257.1 only).
State z = L y, unweighted GTEx projector over all arrays, y = log2 x ln 2. Unit = subject x condition; delta = z(19:00) -
z(07:00).
**Statistics.**
- **H257.1 (positive control):** aligned units: rho of A252 with t0 = 07:00 and 12 h elapsed; exact null over the
  swaps of the two time labels within unit (2^12); matched and responsive companions (200 each, from the arrays'
  genome-wide means); pass: p < 0.05 and above the q95 of both. If it fails, H257.2-3 are uninformative.
- **H257.2 (primary; where the muscle clock sits):** shift s = arg(sum over subjects of delta_mis conj(delta_ctl)) /
  w, in hours in (-12, 12]; 90 % interval from 2,000 bootstraps of subjects. Readings fixed now: interval inside
  [-3, 3] h: "stays on central time" (follows behaviour by less than a quarter of the shift); interval excluding 0 and
  12: "partial shift of s h"; interval containing 12 (or -12): "follows behaviour". Otherwise "not resolved".
- **H257.3 (secondary; damping):** R = rho_mis / rho_ctl (both against the aligned chord); exact within-subject swap of
  the condition labels for rho_ctl - rho_mis > 0 (one-sided, 2^11); 90 % bootstrap interval of R.
**Descriptive.** s without PER1 and NR1D1; the companions' rho_ctl - rho_mis; per subject the angle of delta_mis
against delta_ctl.
**Limits.** (1) Not blind (above). (2) Two clock times per condition, 12 h apart: a rotation of the state and a
change of its amplitude cannot both be read from one chord per condition; s is the angle only. (3) 11 subjects. (4)
Light, sleep and meals all shift together in the protocol; "behaviour" means all of them. (5) Biopsies on the third
day of the misaligned protocol: a slow adjustment would read as a small s.
**Stop rule:** one scoring run. **Code:** `paper/exploratory/A257_nightshift_muscle/scripts/`.

## A258 (2026-11-13, user: "Go for 1, 2 and 3" (item 1 of REVIEW_2026-11-10 §6.4: D1, PI question Q27); a calibration, not a biological test; seen: the review's post hoc feasibility table (12 seeds: d* bias -0.11 to -0.19 at log V ratio -0.8 to -1.0 with anchor noise 0.05-0.15, <= 0.07 at -0.3 to -0.45, none at equal spread) and the observed log V ratios of the equivalence families (A231 medians -0.03 / +0.23, A233 -0.32 to +0.24, largest single tissue 0.76); nothing of the correction below has been computed): how biased are group contrasts at unequal host spread, and does an arm-specific host prior remove the bias? **Simulation (`clocksync.simulate`), then a sensitivity column on the equivalence families; method (Q27)**
**Why.** `Dataset.contrast` regresses each arm on host posterior means shrunk by one pooled prior; in the arm with the
narrower host spread the loading, z2 and d* are pulled toward 0 (review M1). The equivalence claims (sex, age,
obesity, T2D, five chronic conditions) are read at a bound of +/- 0.15 on d*, so the size of the bias at their designs
decides whether they hold as stated.
**Candidate correction (fixed now; an option, the default is unchanged).** For the donors of an arm: with the pooled
prior Sigma and each donor's posterior mean f and covariance S_post, the shrinkage K = I - S_post Sigma^-1 gives the
unshrunk anchor estimate r = K^-1 f and its measurement covariance S = K^-1 Sigma - Sigma; the arm's prior Sigma_g =
cov(r, centred within stratum) - mean(S), projected onto positive semi-definite matrices; the arm's posterior means and
covariances follow from Sigma_g (empirical Bayes per arm). The arm's loading, psi and V are read on them; V0 (the
common spread of d*) stays the pooled one. `clocksync`: `on_host(..., prior="arm")`, `Dataset.contrast(...,
host_prior="arm")`, with a unit test (noisy anchors, unequal windows).
**Simulation design (fixed now).** `clocksync.simulate` with the compartments of the review (three anchors, one target
at coupling as there); two arms; reference window 6-18 h, case window 6-18 / 8-16 / 10-14 / 11-13 h; anchor noise psi
0.002 / 0.05 / 0.15; donors per arm 40 / 80 / 200; 100 seeds per cell (seed 258). Truth: equal coupling and loading in
both arms (d* = 0, log z2 ratio = 0). A power arm: the target's coupling lowered by 0.15 in the case arm at 80 donors.
Outputs per cell and estimator (pooled, arm): mean and SEM of d*, d, log z2 ratio and observed log V ratio; RMSE;
size of the label-permutation test (200 permutations, 30 seeds per cell at 80 donors) and power in the power arm.
- **V258.1 (bias):** pooled estimator: mean d* bias per cell, reported against the observed log V ratio.
- **V258.2 (correction accepted if all hold):** with the arm prior, |mean d* bias| <= 0.02 in every cell with anchor
  psi <= 0.15 and log V ratio >= -1.0 (the range of the real designs); RMSE not above the pooled estimator's by more
  than 10 % at equal spread; power in the power arm not below the pooled estimator's by more than 0.05.
- **V258.3 (real designs):** a bias table indexed by the observed log V ratio and anchor reliability; for every row
  of the equivalence families of A164 / A202, A208, A231, A233 and A235, the expected bias at its observed log V ratio
  (interpolated from V258.1), flagged when it exceeds 0.05.
- **Part 2, only if V258.2 passes (a separate dated note before it is run):** the equivalence families recomputed
  with `host_prior="arm"` as a sensitivity column; registered numbers and verdicts unchanged; a row whose 90 %
  interval leaves +/- 0.15 under the correction is reported as "not robust to the spread correction".
**Reading.** V258.1 sets how far the pooled numbers can be trusted at each design; V258.2 decides whether the paper
reports a corrected column; V258.3 lists the rows to flag. A failed V258.2 leaves the flag table as the answer.
**Code:** `paper/exploratory/A258_spread_bias/scripts/`; `clocksync` changes committed in its repository with a test.

### A244 design addendum (2026-11-13, before any clock statistic on a cell class; full tables in `A244_fetal_adrenal_classes/design_addendum_draft.md` and `output/`, prepared by a Sonnet agent with `prep_a244.py`, `design_a244.py`, `feasibility_a244.py`, checked by the lead)
- **Per-cell rule (set on marker counts alone; clock genes never read):** a class is met when its marker count is >= 3
  and its marker fraction (marker UMIs / cell UMIs) is >= 0.005 (cortex, medulla, stroma) or >= 0.1 (blood); a cell is
  assigned to the single class among cortex, medulla and stroma it meets if it meets no other of the four. Constants
  from the valleys of the marker-fraction histograms and mixture boundaries (0.004-0.007). 448,931 cells: cortex
  71,257, medulla 29,227, stroma 40,736, left out 307,711.
- **Agreement with published annotations:** Kildisiute purity cortex 0.991, medulla 0.999, stroma 0.990 (recall 0.624,
  0.861, 0.227); Dong purity 0.766 (cortex holds endothelial and capsular cells with ambient steroidogenic RNA), 0.963,
  0.934. The other four cohorts have no annotation on disk.
- **Design (fetuses with >= 200 cells of the class):** cortex 34 fetuses in 6 cohorts (kameneva 10, wang 9, sasaki 5,
  kildisiute 5, dong 3, han 2; median clock-gene count per cohort 19.5 to 334); medulla 30 in 6 (wang 10, kameneva 7,
  dong 4, sasaki 4, kildisiute 3, han 2; 4.75 to 118); both classes 26.
- **Foreign RNA:** median a 0.0064 (cortex) and 0.0164 (medulla); expected foreign share of the clock-gene counts 0.009
  and 0.013 (largest cohort: dong medulla 0.060): no class is confounded. The reference of a (the other class's
  pseudobulk of the fetus) is taken at any size, as the registration states; for 26 of 33 cortex and 26 of 30 medulla
  fetuses it has >= 200 cells; the `:corrected` readouts restricted to those are reported as descriptive.
- **Adult power (93 GTEx adrenals, 200 subsamples):** pass rate at the full clock cortex 1.000, medulla 1.000 (0.75:
  0.995, 0.910; 0.5: 0.520, 0.335; shuffled 0.055 each); 14 % of the cortex targets capped. **H244.1 scores `cortex` and
  `medulla`** (all ages, cohorts with >= 2 fetuses, cohort indicator and quadratic age, permutations within cohort);
  windows `late` (cortex 0.910, medulla 0.570) and `early` are descriptive as registered.
- **H244.2: not estimable, not computed.** Second adult tissue fixed by rule before any fetal statistic: the A241
  stand-in tissue sharing most 08-18 h donors with the adrenal, excluding blood: skeletal muscle (80 donors). Adult
  pass rate 0.395 (0.327-0.466) at the full clock, mis-paired false-positive rate 0.040 (0.017-0.077): below 0.8.

## A259 (2026-11-15, roadmap item A1 (c) under the lead prompt `prompts/target_line_what_internal_time_controls_2026-11-15.md` (user, 2026-11-15: target plan 1, gate G3); staged: estimator, data, hypotheses and readings are fixed here, the power note (simulation only, `A259_tissue_deviation/scripts/power_a259.py`, running at registration) goes in a dated addendum before any observed statistic; seen: the offset equivalences of A208, A231, A233, A235 (roadmap A1 (a), `roadmap_A1_offset_bound/output/offsets_hours.csv`), A133 (no circular-part estimator), A255 (fat-skin pairing r 0.533), the psi of past fits (noise around the host, not split into true deviation and measurement noise); the estimator was checked on simulated data only; NOT seen: any split-half or deviation statistic of real data): how far does a tissue's clock depart from the clock state its person's other tissues share, in hours, net of measurement noise? **GTEx fast deaths (11 body tissues as targets), MuTHER fat and skin (living); L2 / L3, person-specific deviation; primary for gate G3**
**Why.** Plan 2's headline is a bound in hours: one internal time per person. The group offset equivalences bound
differences between groups of people, not the spread of a tissue around its own person's body time. The same quantity
decides whether the tissue's own clock can be read apart from the body's in humans (roadmap T2): a deviation that is
not larger than measurement noise leaves nothing for a "tissue clock" regressor to read.
**Estimator (fixed now; `clocksync.deviation.incoherence`, committed with a test).** Per target compartment k:
1 + kappa_k = |m(beta_k^A, beta_k^B)| |m(g1, g2)| / sqrt(|m(beta_k^A, g1)| |m(beta_k^B, g2)| |m(beta_k^A, g2)| |m(beta_k^B, g1)|),
m(x, y) = mean x conj(y) over donors (centred), beta^A / beta^B the clock-plane states of disjoint halves of the 14
clock genes (one projector per half over all rows: GTEx counts, log fractions with pooled median-count weights;
arrays, log2 x ln 2, unweighted), g1 / g2 composites of two random halves of the host compartments (each rotated onto
the frame of the best-sampled host). Measurement noise is independent between gene halves and between samples, so
kappa_k = var(person-specific deviation) / (clock variance of k). sigma_k = (12 / pi) sqrt(kappa_k) is its SD in hours
if the deviation were pure phase; amplitude deviations count as phase, so sigma_k is an upper bound. With one host
compartment (MuTHER) the pair statistic (1 + kappa_a)(1 + kappa_b) - 1 is used, sigma_pair = the SD of the relative
phase of the two tissues. 20 gene splits (balanced, each half's centred design with condition number <= 10), 20 host
splits, 1,000 donor bootstraps (seed 259), 90 % percentile intervals.
**Data.** (i) GTEx v10 fast deaths (Hardy 1-2), clock-gene table of `gtex_cues/_shared.clock_table` (depth = library
size), targets: the 11 body tissues of `_shared.BODY` without whole blood (no coherent clock); hosts of a target: the
other ten; a target is eligible with >= 80 donors having it and >= 1 tissue in each host half. (ii) MuTHER (E-TABM-1140),
subcutaneous fat and skin of the same women (A255 matrices, `A254_paired_lesions/output/prep/`), pair mode; the cell
line is not a clock compartment (A255) and is not used. (iii) Descriptive, same pair in GTEx: subcutaneous adipose with
each skin site.
**Companions (same splits).** 100 matched and 100 stress-responsive pseudo-clock waveforms (`paired.pseudo_waveforms`,
from the mean log fraction over all rows): for each, the within-sample split-half moment |m(beta^A, beta^B)| of the
target divided by the target's clock variance (|m(beta_k^A, beta_k^B)| / (1 + kappa_k)) gives kappa_tech, the share
of non-clock sample-level variation that projects onto the plane; and |m(beta_k, g)| of the pseudo set over the
clock's, the share of the cross-tissue moment carried by non-clock variation shared by a donor's samples.
**Hypotheses.**
- **H259.1 (primary, the bound; G3):** the median over eligible GTEx targets of the upper 90 % bound of sigma_k is
  below 1.0 h. Also reported: the number of targets whose upper bound is below 1.0 h.
- **H259.2 (secondary, existence; per target, Holm over targets):** kappa_k > 0 (one-sided bootstrap p < 0.05) and
  kappa_k above the q95 of kappa_tech of both companion families. A pass means a person-specific deviation of the
  tissue's clock from the body time, beyond technical variation.
- **H259.3 (secondary, living):** MuTHER fat-skin upper 90 % bound of sigma_pair below 1.41 h (sqrt 2 x 1 h: both
  tissues within 1 h).
**Readings (fixed now).**
- G3 (with A1 (b)): plan 2's bound is a headline if H259.1 and H259.3 pass and the A1 (b) simulation leaves the four
  scored offset equivalences standing; it is written as a PI question, with the evidence, before anything that depends
  on it.
- H259.1 fails with H259.2 passing in several tissues: tissue clocks depart from the body time by a measurable,
  person-specific amount; that is a result for plan 1 (the L term of T2 has something to read) and against the "one
  time per person" wording.
- H259.2 fails in a tissue at a detection limit delta (power note): a tissue-clock regressor cannot separate the
  tissue's clock from the body's there for deviations below delta; T2's L term is not registered for that tissue.
- A companion family reaching kappa_k (kappa_tech q95 >= kappa_k): the deviation is not shown to be clock; sigma_k
  stays an upper bound only. A cross-tissue companion share above 0.2 is flagged: kappa may be understated.
**Limits.** (1) sigma counts amplitude deviation as phase (conservative for the bound, not for H259.2's "phase"
wording: H259.2 says "deviation"). (2) The hours conversion assumes host phases spread over the day (GTEx deaths
cover it; MuTHER visit times are unknown: a narrow window inflates kappa, conservative for the bound). (3) Non-clock
variation shared by a donor's tissues lowers kappa (companion share reported). (4) GTEx is post mortem: agonal and
ischaemic variation within a sample counts as deviation (conservative). (5) GTEx time of death and phenotypes are
controlled access: aggregates only.
**Stop rule:** one scoring run after the power note is dated. **Code:** `paper/exploratory/A259_tissue_deviation/scripts/`.

## A260 (2026-11-15, roadmap item T0c under the lead prompt of 2026-11-15; a calibration of the readers phase T will use, for gate G2 (ii); seen: the scout `exploratory/scouts_2026-11-15/T0c_calibration/` (timed cosinor class sizes and timed feeding shifts: local 584, systemic 40, knockout-only 26 genes; control liver clock -11.56 h, WAT -9.24 h, lung -2.31 h under daytime feeding), A111 / A129 (liver coupling 0.967 -> 0.676 without its clock); NOT seen: any untimed per-gene loading of these data): does the untimed procedure put mouse liver genes in the right class, with and without a local clock and with feeding moved? **Manella 2021 (GSE159135), Alb-Cre control and hepatocyte Bmal1 knockout x ad libitum (AL) / daytime feeding (DF), liver, WAT, lung, 24 mice per arm; per-gene untimed loadings against timed truth; calibration (G2 (ii))**
**Why.** No reader of the cue line was ever checked where the answer is known (inventory 2026-11-15, point 5). Here
it is known from the recorded times: genes that keep their rhythm in a liver without its clock are driven from outside
the hepatocyte (body or feeding), genes that lose it need the local clock, and daytime feeding moves feeding-driven
genes with the food.
**Truth (timed, recorded ZT; fixed as in the scout).** Per arm and tissue, cosinor of log((count + 0.5) / depth) on ZT
(genes with mean count >= 5), BH q. Liver, non-clock genes (the 14 GTEx clock genes excluded): **local** = q < 0.05 in
control AL, p > 0.2 in knockout AL and amplitude ratio knockout / control < 0.5; **systemic** = q < 0.05 in both and
ratio >= 0.5. Timed shift of a gene, DF minus AL, relative to the lung: its phase difference minus the median
difference of the lung's rhythmic clock genes of that genotype.
**Untimed procedure (times hidden).** `Dataset.gene_rhythms` on the clock-gene table of `common.data.load("manella2021")`
with the genome table of the builder's `read_tissue` (liver), one arm at a time, host anchors as below, `n_null` 50 (host
states permuted): complex loading c, amplitude |c|, acrophase in the anchors' frame, Wald p.
- **V260.1 (body reader):** host from WAT and lung. In the knockout AL liver, -log10 p of the untimed loading separates
  systemic from local genes: AUC >= 0.80 with the lower 90 % bound (gene bootstrap, 2,000) above 0.65.
- **V260.2 (local reader by contrast):** log(|c| knockout AL / |c| control AL) on the same anchors, against the timed log
  amplitude ratio, over genes with timed q < 0.05 in either arm: Spearman >= 0.60; and its AUC for local vs systemic
  >= 0.80 (lower 90 % bound above 0.65).
- **V260.3 (shift reader; the untimed analogue of T1):** host from the lung alone (feeding moves WAT). Per genotype,
  untimed acrophase DF minus AL per gene against the timed shift relative to the lung, over genes with timed q < 0.05 in
  both arms: circular correlation >= 0.70 in the control; the knockout is reported (few genes are rhythmic in both of
  its arms).
- **D260.4 (descriptive; tissue-deviation reader):** A259's kappa of the liver against WAT and lung (target mode, one
  host per half) in control AL and knockout AL; and the share of local genes whose untimed loading on the liver's
  deviation (liver clock state minus its prediction from WAT and lung, both from `clocksync.deviation.states`) has p <
  0.01 in control AL. No pass rule: whether a tissue term can be read when the local clock is coupled is the question
  A259 answers in humans.
**Companions.** Untimed calls are compared with the host-permuted null of `gene_rhythms` (empirical FDR); AUCs are
against the systemic / local truth, so a reader that calls every rhythmic gene scores 0.5.
**Readings.** V260.1 and V260.2 pass: the untimed body reader is calibrated: in a liver without its clock, genes that
still follow the body's state are the ones that keep a rhythm, and the knockout / control contrast of untimed loadings
recovers clock dependence (G2 (ii) met for the body / local split). V260.3 pass: the shift reader recovers feeding-driven
moves (the basis of T1). A failure is reported as a limit of the procedure, and phase T designs that rest on the
failed reader are not scored (G2 (ii) not met for that reader).
**Limits.** Same mice for truth and procedure (unavoidable: the truth is the hidden times). The mouse label is assumed
to be the same animal in the three tissue files (A111 assumption). MARS-seq 3' counts of about 11.6 k genes; depth is
the sum of retained genes. A systemic class of 40 genes.
**Stop rule:** one scoring run. **Code:** `paper/exploratory/A260_mouse_calibration/scripts/`.

## A261 (2026-11-15, roadmap item T1 under the lead prompt of 2026-11-15; **NOT BLIND**: A124 and A257 scored these arrays for the 14 clock genes and the A81b cue / locked sets (A257: muscle clock shift 2.68 h, excursion ratio 0.233); seen: the scout `exploratory/scouts_2026-11-15/T1_nightshift/` (reference cosinor of GSE108539, residual SDs of the GSE106800 differences around their condition means, simulated separability); NOT seen: any program-level or non-clock gene-level contrast between the two conditions): which daily programs of human muscle follow central time, the muscle clock, or behaviour after a 12 h shift of behaviour? **GSE106800 (Wefers 2018), 11 men with four arrays; programs of reference-rhythmic genes; L1 within person; primary for plan 1 (T1)**
**Why.** The one human design that moves behaviour (12 h), the muscle clock (2.68 h, damped to 0.233: A257) and central
time (about 0 h) by known, different amounts. A124 failed with gene-set states (MDE 7.3 h); here each program is read
against gene phases from a timed reference, and the excursion ratio separates central time from the damped tissue
clock.
**Data.** Arrays and probe map of A257 (main probe sets, mean log2 per canonical gene with the eight aliases); 11
subjects with 07:00 and 19:00 arrays in both conditions; per subject and condition d = y(19:00) - y(07:00) (log2).
Reference: GSE108539 exon layer (10 subjects x 6 biopsies), per-gene cosinor with subject intercepts, log2 RPKM > 0 on
average; amplitude A_g, acrophase phi_g (clock time); rhythmic: p < 0.01. The 14 clock genes are excluded from every
program.
**Programs (fixed now).** (i) Hallmark sets (MSigDB v2024.1, `$SCCC_DATA/reference/h.all.v2024.1.Hs.symbols.gmt`) with
>= 20 reference-rhythmic genes on the array; (ii) "all rhythmic": every reference-rhythmic gene on the array; (iii) six
phase programs: reference-rhythmic genes binned by phi_g in 4 h bins (0-4, ..., 20-24 h), descriptive.
**Statistic.** Per subject and condition, z = regression (no intercept) of the program's d_g on (C_g, S_g) = 2 A_g
(cos, sin)(w (19 - phi_g)), as a complex number. Shift s = arg(mean z_mis / mean z_al) / w (h), ratio r = |mean z_mis /
mean z_al|; 90 % jackknife-t intervals over subjects. Hypotheses (s, r): central (0 h, 1), tissue clock (2.68 h,
0.233), behaviour (12 h, 1); each is tested by the Hotelling T2 of the per-subject residual z_mis - r e^{iws} z_al
(F(2, n - 2)), **excluded at p < 0.01** (scout: truth excluded in at most 0.092 of simulated programs; alternatives
excluded with power >= 0.915 at 20 genes).
**Rhythm gate and companion.** A program is read only if |mean z_al| exceeds the q95 of 200 expression-matched random
gene sets of the same size drawn from reference non-rhythmic genes (p > 0.5, matched on mean array intensity in deciles),
each carrying its own reference A_g, phi_g (the specificity companion: a set without a reference rhythm should not
look rhythmic in the aligned condition). The same 200 sets give the null distribution of the class calls.
**Hypotheses.**
- **H261.1 (positive control):** the "all rhythmic" program passes the rhythm gate and its aligned loading has an
  argument within 3 h of 0 (the reference phases transfer to the arrays).
- **H261.2 (primary, the partition):** per gated program, the class: the single hypothesis not excluded (central,
  tissue clock or behaviour), "mixed" (all three excluded) or "unresolved" (two or three not excluded). Reported: the
  counts and shares of gated Hallmark programs per class, with s and r and their intervals per program.
- **H261.3 (primary, the genome-wide program):** the class of "all rhythmic".
**Readings (fixed now).** If H261.1 fails, H261.2-3 are uninformative. "Central" for a program: its daily change kept
both direction and size under a 12 h shift of behaviour: a program of the body's central time (light-entrained,
systemic). "Tissue clock": it moved and damped like the muscle clock. "Behaviour": it reversed with sleep, activity
and meals. "Mixed": none of the three alone (a mixture of times, or another one). The partition is "mostly X" when X
holds more than half of the gated programs; any share is reported with the number of programs, which overlap in genes
(not independent units).
**Limits.** (1) Not blind (clock genes only seen). (2) Two clock times per condition, third day of the protocol. (3) 11
men. (4) Light, sleep and meals move together: "behaviour" is all of them. (5) The reference cohort ate hourly while
awake; meal-driven rhythms are under-represented among reference-rhythmic genes. (6) Programs share genes.
**Stop rule:** one scoring run. **Code:** `paper/exploratory/A261_nightshift_programs/scripts/`.

### A259 power note (2026-11-15, simulation only, before any observed statistic; `A259_tissue_deviation/output/power/summary.csv`, 40 data sets per cell)
- **Target mode** (ten hosts): unbiased without a shared program (mean kappa bias -0.003 to +0.004 at amplitude 1, up to
  +0.041 at amplitude 0.5 and 100 donors). Detection (lower 90 % bound > 0) at amplitude 1: 200 donors 0.575 at 0.75 h,
  0.850 at 1.0 h; 400 donors 0.750 at 0.75 h, 0.950 at 1.0 h; at amplitude 0.5, 200 donors: 0.300 at 1.0 h, 0.725 at
  1.5 h. False detection at sigma 0: 0.000-0.125. Upper bound below the truth (nominal 0.05): 0.000-0.175.
- **The bound depends on clock strength.** At sigma 0 the median upper 90 % bound is 0.70 h (amplitude 1, 200 donors)
  and 1.72 h (amplitude 0.5, 200 donors): H259.1 can pass only where tissue clocks are about as strong as the GTEx
  waveform; a failure in weak-clock tissues is a power limit, not evidence of deviation (read with H259.2).
- **Shared non-clock program across a donor's tissues** (shared_sd 0.1): kappa is biased **up** by 0.125-0.155 (sigma
  1.28 h at a true 0), and the lower bound exceeds 0 in 0.75 of data sets at sigma 0. This corrects limit (3) of the
  registration (written as "lowers kappa"): such variation inflates kappa, which keeps H259.1 conservative but makes
  H259.2 (existence) depend on its companions. Reading added before scoring: an H259.2 pass with a cross-tissue
  companion share (median of the matched sets) above 0.2 is reported as "not separable from shared non-clock
  variation", not as a person-specific clock deviation.
- **Pair mode** (MuTHER design, 680 persons): kappa unbiased (-0.022 to +0.029). At sigma 0 the truth already includes
  the host tissue's own noise (0.57 h at amplitude 1, 1.14 h at 0.5); median upper bound 0.89 h (amplitude 1) and
  1.78 h (0.5) against the registered 1.41 h.

## A262 (2026-11-15, roadmap G2 (iii) under the follow-up of 2026-11-15; seen: the scout `exploratory/scouts_2026-11-15/T5_blood_misalignment/` (designs, a GSE39445 reference cosinor, within subject x condition SDs, simulated power); A156 read GSE48113 for the 14 clock genes and the A81b cue / locked sets only (timed: core state -1.36 h, cue state -6.85 h from day 1 to day 4); NOT seen: any program-level or non-clock gene-level contrast between conditions in GSE48113 or GSE122725): in human blood, does the genome-wide rhythmic program follow central time or the displaced sleep schedule? **GSE48113 (Archer 2014, forced desynchrony) and the dim-light arm of GSE122725 (= GSE107537, Kervezee 2018, simulated night shift); programs of reference-rhythmic genes; L1 within person; the second misalignment design for plan 1**
**Why.** A261 places the rhythmic program of muscle on behaviour (shifted or evoked). The only misaligned human designs
with public genome-wide data beyond GSE106800 are in blood (scout T5: no muscle or adipose design with shifted meals or
sleep and public matrices). The literature suggests blood goes the other way (Kervezee 2018: 676 of 929 baseline-rhythmic
probe sets stay on the habitual schedule; Archer 2014: rhythmic transcripts fall from 6.4 % to 1.0 %). Blood has no
coherent core clock for this project (expected exception), so no host-referenced readout is made: the question is only
which time the rhythmic program follows, with the A261 machinery.
**Data.**
- D1 GSE48113 (series matrix, GPL15331 symbols, mean log2 per symbol): subjects with >= 3 samples on both day 1
  (aligned) and day 4 (sleep and meals displaced by about 12 h; melatonin at about the same clock time on both days,
  Archer 2014); time = recorded clock time of the sample.
- D2 GSE122725, experimental group "Control" (dim light; the 8 subjects of GSE107537): day-oriented (aligned) vs
  night-oriented (sleep delayed by 10 h for 3 days) schedule, subjects with >= 3 samples in both; time = recorded
  relative clock time (mod 24 through the cosine). Symbols from the first RefSeq record of the GPL24324 annotation
  (main probe sets), mean log2 per symbol.
- Reference (central time without sleep or meal rhythms): GSE39445 sleep-extension arm under extended wakefulness;
  per gene (mean log2 per GPL15331 symbol, the lowest-intensity fifth dropped), cosinor on the melatonin-referenced
  circadian phase with subject intercepts and a linear term in hours awake; amplitude A_g, phase phi_g; reference-
  rhythmic: p < 0.01. Phases are mapped to clock time by adding the circular mean of (clock time - circadian phase) over
  the reference samples. The 14 clock genes are excluded from every program.
**Programs (fixed now).** (i) "all rhythmic": every reference-rhythmic gene present in the dataset (primary); (ii) six
phase programs by 4 h bins of the clock-time phase (descriptive). Hallmark programs are not scored (scout: programs of
20 genes are resolved in only 0.01 to 0.37 at a quarter of the reference amplitude).
**Statistic.** Per subject and condition, values centred per gene within the unit; z = a + i b from the regression (no
intercept) of the program's centred values on A_g (cos, sin)(w (t - phi_g)) over genes and samples, so a shift s
(later = positive) and excursion ratio r give z = r e^{i w s}. Shift s = arg(mean z_mis / mean z_al) / w, ratio
r = |mean z_mis / mean z_al|, 90 % jackknife-t intervals over subjects (as A261). Hypotheses, each tested by the
Hotelling T2 of z_mis - r e^{i w s} z_al over subjects (F(2, n - 2)), **excluded at p < 0.01**:
- D1: central (0 h, 1), behaviour (12 h, 1), lost (r = 0, i.e. mean z_mis = 0);
- D2: central (1.5 h, 1: the dim-light melatonin delay of about 1 to 2 h, Kervezee 2019), behaviour (10 h, 1), lost;
  sensitivity: central at 0 h and at 2 h.
Class: the single hypothesis not excluded; "mixed" if all are excluded; "unresolved" if two or three are not.
**Rhythm gate and companion.** A program is read only if |mean z_al| exceeds the q95 of 200 random sets of the same
size drawn from reference non-rhythmic genes (p > 0.5), one per program gene from the same decile of mean intensity in
the dataset, without replacement while the decile has unused genes, then with replacement (the non-rhythmic pool is
smaller than the genome-wide program), each carrying its own reference A_g, phi_g; seed 262. The same sets give the
null distribution of the class calls (specificity companion: no class should be called in sets without a reference
rhythm).
**Hypotheses.**
- **H262.1 (positive control, per dataset):** the "all rhythmic" program passes the gate and its aligned argument
  arg(mean z_al) / w is within 3 h of 0 (the reference phases transfer).
- **H262.2 (primary):** the class of "all rhythmic" in D1.
- **H262.3 (primary):** the class of "all rhythmic" in D2.
**Readings (fixed now).** If H262.1 fails in a dataset, its class is uninformative. "Central": the blood program kept
its timing and size while sleep (and in D1 meals) moved. "Behaviour": it moved with the sleep schedule. "Lost": its
rhythm is abolished under misalignment (the Archer reading). "Mixed": none alone (partial shift and damping).
"Unresolved": not separable. If both D1 and D2 are "central" (or one central and one lost) and A261 is "behaviour",
the result is stated as "the rhythmic programs of muscle and blood follow different times under the same kind of
misalignment" with the limits of each design; the sentence "which time a tissue's outputs follow is tissue-specific"
is written only then, as a hypothesis from two tissues. If blood is also "behaviour", no tissue difference is claimed.
**Limits.** (1) GSE48113 not blind for the clock and cue genes (A156). (2) Blood rhythms include diurnal cell
composition (memory: blood composition clock); a "central" program may be composition driven by central hormones, not
a cellular clock. (3) D1: sleep, meals and posture move together and the sampling day contains the displaced sleep
(masking possible, but a full cycle is sampled); D2: sampling under constant posture with hourly snacks (no meal
response), sleep inside the window, central time known only at group level. (4) Possible subject overlap between the
reference (GSE39445) and D1 (same laboratory; different identifiers). (5) D2 has 8 subjects.
**Stop rule:** one scoring run. **Code:** `paper/exploratory/A262_blood_misalignment/scripts/`.

## A264 (2026-11-15, roadmap item A2 under the lead prompt `prompts/biology_paper_checks_2026-11-15.md`, check 1a; **NOT BLIND**: every fetal adrenal statistic of A189, A242 and A244 has been seen (A244 cortex C 0.608, n 34, Holm p 0.0076, clock-specific; medulla -0.048); seen before writing: the adult power below and two aggregates of the adult adrenal pool; NOT seen: any fetal statistic under the adjustments below): is the clock-shaped covariation of the fetal adrenal cortex a static axis between fetuses? **A244 cortex readout (6 cohorts, 34 fetuses) under three further adjustments; L1; secondary (development section)**
**Why.** A244 localises the fetal adrenal covariation to the steroidogenic cortex, but a static axis between fetuses
(gestational age or cohort beyond the A244 adjustment, steroidogenic activity) would pass the same test, and two of the
14 clock genes (PER1, NR1D1) respond directly to glucocorticoid / cAMP signalling. Each adjustment below removes one such
axis; none can show that the covariation is a running clock.
**Data and statistic (fixed; `design_a264.py`, `score_a264.py`).** Rows, pseudobulks, statistic, permutations and
companions as the `cortex` readout of A244: signed coherence C over the cortex pseudobulks (>= 200 cells) of the cohorts
with >= 2 fetuses (dong 3, han 2, kameneva 10, kildisiute 5, sasaki 5, wang 9), log fractions residualised gene by gene on
the covariates; 10,000 per-gene permutations within cohort; 200 expression-matched and 200 responsive companion sets,
adjusted in the same way and drawn for the waveform of the readout; clock-specific: C above the q95 of both. Seed 264.
- `cortex:age` (i): cohort indicator, common **cubic** age polynomial and a linear age slope of its own for every cohort
  with >= 5 fetuses (kameneva, kildisiute, sasaki, wang): 12 covariates (A244: 7).
- `cortex:steroid` (ii): the A244 covariates plus a steroidogenic activity score: per fetus, the equal-weight mean over
  20 genes of the gene's log fraction standardised over the cortex pseudobulks. Genes, fixed from the literature before
  any statistic of this amendment: those raised >= 4-fold by 48 h of ACTH in **both** adult and fetal primary human
  adrenocortical cells (Xing, Parker, Edwards, Rainey 2010, J Mol Endocrinol 45:59, Table 3, PMC3725639): MRAP, CYP11B1,
  HSD3B2, CYP21A2, MVD, SCARB1, MT3, IDI1, MC2R, AMDHD1, ACSS2, TMEM97, GSTA4, STAR, CYP11A1, EBP, FDX1, CYP17A1, EPB41L1,
  INHA (none is a clock gene). In the adult adrenal pool (93 GTEx glands, deaths at 08-18 h) the 20 genes cohere (mean
  pairwise r 0.485) and the score shares almost no variance with the hour of death (cosinor R^2 0.004;
  `output/adult_scores.csv`).
- `cortex:-PER1-NR1D1` (iii): the A244 covariates; the waveform without PER1 and NR1D1 (12 genes), companions of 12 genes.
- Descriptive (no verdict): `cortex` (A244 reproduced: it must return C 0.608), `cortex:all3` ((i) to (iii) together),
  `cortex:acute` (A244 covariates plus the score of NR4A1, NR4A2, NR4A3, the immediate-early cAMP targets named in the
  prompt; not taken from Xing 2010, whose design was 48 h); leave-one-gene-out range; unsigned R.
**Power first (done; `feasibility_a264.py`, 278 s; `output/feasibility_passrates.csv`; 200 subsamples, 500
permutations).** As A244's feasibility: each fetus stood for by a distinct GTEx adrenal gland thinned to the cortex
clock-gene counts of its cohort, which also brings its own adult activity scores. Pass rate at alpha 0.05, adult clock at
1 / 0.75 / 0.5, and without a clock: `cortex:age` 1.000 / 0.990 / 0.365, 0.055; `cortex:steroid` 1.000 / 1.000 / 0.460,
0.080; `cortex:-PER1-NR1D1` 1.000 / 0.995 / 0.545, 0.025; `cortex` (A244 design) 1.000 / 0.995 / 0.440, 0.065;
`cortex:all3` 1.000 / 0.910 / 0.455, 0.075. Every registered readout is estimable (>= 0.8 at the full clock). The adult
power does not test separability from steroidogenic activity: in the adult the score is independent of the clock (R^2
0.004), so the adjustment costs no adult power by construction.
**Hypotheses (each read on its own: a robustness check of one standing result, not a discovery family; no Holm).**
- **H264.1:** `cortex:age` p < 0.05 and clock-specific.
- **H264.2:** `cortex:steroid` p < 0.05 and clock-specific.
- **H264.3:** `cortex:-PER1-NR1D1` p < 0.05 and clock-specific.
Verdict per readout: PASS; "not clock-specific" (p < 0.05, not above both companions); fail (p >= 0.05).
**Readings (fixed now).**
- All three PASS: the cortex covariation is not carried by a smoother age trend or cohort-specific age slopes, by the
  ACTH-driven steroidogenic program, or by the two signal-induced clock genes. The fetal result stays a candidate for the
  main text, stated as "clock-shaped covariation that survives these adjustments"; a static axis of another kind is not
  excluded, and a running clock is not shown.
- H264.2 fails or is not clock-specific: **the covariation is explained by steroidogenic activity** (a static axis, or a
  clock that cannot be told from it in these data): the fetal result moves to the supplement as a hypothesis, and the
  development section loses "the adrenal starts first" as a main-text claim (to the PI).
- H264.1 fails or is not clock-specific: the covariation rests on the age and cohort structure of the pool; same move.
- H264.3 fails or is not clock-specific: the covariation rests on PER1 and NR1D1; it is not shown to be a clock-network
  state; same move.
- Any one of the three moves the fetal exception to "fallen" in the tally of the prompt (stop rule 1); all three PASS
  keeps it "standing, not shown to be running".
**Limits.** (1) Not blind (above). (2) Five of the 20 score genes are also cortex classifier markers of A244 (STAR,
CYP11A1, CYP17A1, CYP21A2, CYP11B1): within the class, their fraction mixes activity and purity; the score is used as one
covariate, not interpreted. (3) The score is from 48 h ACTH in culture, not from acute or in vivo signalling. (4) If the
fetal clock drives steroidogenesis, adjusting for the score removes clock signal too: H264.2 failing cannot tell
"explained by" from "drives"; the reading above treats both as "not shown". (5) Limits 1 to 6 of A244.
**Stop rule:** one scoring run (`score_a264.py --jobs 6`); any further cut is a new amendment. **Code:**
`paper/exploratory/A264_fetal_adrenal_covariates/scripts/{design_a264.py,feasibility_a264.py,score_a264.py}`.

## A265 (2026-11-15, roadmap item A3 under the lead prompt `prompts/biology_paper_checks_2026-11-15.md`, check 2; seen: A138 H138.4 on the 517 TCGA pairs of the same ten cancers (pooled clocksync 0.365, null mean 0.122, matched q95 0.260, responsive q95 0.310), A204 (COFE on 14 clock genes in a simulator, descriptive), the scout `output/scout_design.csv` (design counts) and one COFE run on gene-permuted BRCA tumours for runtime (134 s, 8 cores); NOT seen: any COFE label of real TCGA data): COFE against clocksync on the TCGA tumour-normal pairs: which reading is right, tumour clocks uncoupled from the host (COFE) or paired with the patient's normal tissue (A138)? **TCGA, ten cancers, 464 pairs; simulation on the real design; other; primary for the cancer chapter (head-to-head)**
**COFE's claim, verified at the source** (Ananthasubramaniam and Venkataramanan 2025, PLoS Biol 23:e3003196; bioRxiv
2024.03.13.584582; S4 Text and S2 Table). COFE (sparse cyclic PCA, sparsity by unsupervised cross-validation) was run on
each of 11 TCGA adenocarcinomas separately, on tumours of one histological type (S2 Table), one tumour per patient,
recount3 counts, protein-coding genes with mean count >= 64, log2 cpm (DESeq2), genes with mean log2 cpm >= 7, defaults
(K 5, 3 repeats, 5 restarts). Labels of tumours "lay throughout the circadian cycle in all ACs"; patient-matched healthy
samples (10 of the 11 cancers, "a small fraction of the original patients") were compared with their tumours: "the
difference between the predicted time labels for healthy and tumor tissue was distributed throughout the circadian cycle
in most ACs", with no statistic reported; control: held-out tumour samples of the same patients cluster near zero.
Conclusion: "the clocks in tumors do not have a fixed timing relationship with clocks in healthy tissue" (abstract:
rhythms "uncoupled with clocks in healthy matched tissue"). How the healthy samples were labelled is not stated; this
amendment uses COFE's rule for samples not used in training (projection on the trained loadings, `run_cofe(test=)`), the
reading implied by the held-out control (assumption, stated as such).
**Data (fixed; `lib_a265.py`).** The ten cancers with normals (BLCA, BRCA, COAD, KIRC, KIRP, LIHC, LUAD, PRAD, THCA, UCEC;
COFE's eleventh, OV, has none). Tumours: primary tumours of the S2 Table histology, one per patient (lowest barcode).
Normals: one solid-tissue normal per patient whose tumour is in that set. 464 pairs (BLCA 19, BRCA 91, COAD 39, KIRC 72,
KIRP 32, LIHC 49, LUAD 43, PRAD 51, THCA 49, UCEC 19), 288 to 782 tumours per cancer, 830 to 1,321 COFE genes
(`output/scout_design.csv`). COFE's transform is log2(cpm + 1) with the protein-coding library size (S4 Text names DESeq2
without the call: deviation). COFE as `common.benchmark.cofe.run_cofe` (COFE's fits, folds and CV error; its driver
replaced for memory and a shared-object bug, A138), seed 265 + cancer index.
**Readouts.**
- (i) *COFE's readout reproduced (descriptive with a criterion):* per cancer, Rayleigh test of the tumour-minus-normal
  label differences; "uniform" if p >= 0.05; largest empty arc of the tumour labels; CV minimum interior. Reproduced if
  the differences are uniform in >= 6 of 10 cancers ("most") and tumour labels have no empty arc > 12 h.
- (ii) **H265.1** *pairing of COFE's own labels:* R_c = |mean e^{i (phi_T - phi_N)}| per cancer, pooled with weights n
  pairs; one-sided p against 10,000 permutations of the normals within cancer. clocksync on the same 464 pairs: the A138
  statistic (`q3_lib.corr_stat`; normals gated by `reference_gate`, pooled over gated cancers), same permutations
  (descriptive: A138 is the registered clocksync result on 517 pairs).
- (iii) **H265.2 (decisive)** *simulation on the real design (`score_a265.py --part sim`).* Per cancer the real tumour
  and normal counts, with the normals' patient labels permuted (real pairing broken); every patient gets a time t, normal
  phase t, tumour phase t + N(0, 1 h) (coupled) or uniform on 24 h (uncoupled); counts multiplied by
  exp(A_g cos(w (phase - phi_g))) for a rhythm program fixed from GTEx (A27 timed cosinor of log CPM on the hour of death:
  every gene with BH q < 0.05 in >= 1 of 13 tissues, plus the 14 clock genes; 365 genes; amplitude the median over its
  significant tissues, acrophase their circular mean). Arms: `coupled_surgical` (t uniform 08-18 h), `uncoupled_surgical`,
  `coupled_wholeday` (t uniform 0-24 h); 10 replicates each, 1,000 permutations per replicate; both methods exactly as
  in (i) and (ii). Rates per arm: clocksync pass (pooled p < 0.05), COFE pairing pass (H265.1 statistic p < 0.05), COFE
  published reading "uncoupled" (uniform differences in >= 6 cancers).
**Readings (fixed now).**
- *Sampling window explains the disagreement:* in `coupled_surgical`, COFE's published reading says "uncoupled" in
  >= 5 of 10 replicates while clocksync passes in >= 8 of 10, and clocksync passes in <= 2 of 10 `uncoupled_surgical`
  replicates. Then COFE's published negative is not evidence against coupling under surgical hours, and A138 stands.
- *COFE right, A138 carried by something else:* clocksync passes in >= 3 of 10 `uncoupled_surgical` replicates (its
  pairing statistic picks up something other than the simulated coupling in these data), or H265.1 fails on real data
  while COFE's pairing statistic passes in >= 8 of 10 `coupled_surgical` replicates (COFE would have seen real coupling).
  Then the A138 companions are checked again (matched and responsive q95 on these 464 pairs, new amendment) and the
  tumour-normal pairing claim is suspended in the roadmap.
- *Real data:* H265.1 PASS (p < 0.05): COFE's own labels pair tumour with normal beyond permuted pairs, contradicting the
  published reading with COFE's output. H265.1 fail is read with the `coupled_surgical` COFE pairing rate: below 0.8,
  uninformative.
- Anything else: reported as observed, with the rates; no reading.
**Limits.** (1) Normal labelling by projection is an assumption (above). (2) The simulated program comes from GTEx
cosinors on the hour of death (attenuated, winner's curse on the amplitudes); tumour amplitude equals normal amplitude.
(3) The background keeps each sample's real clock, unrelated to the simulated time after the normals are re-paired; this
noise is the same in every arm. (4) 10 replicates per arm: rates have wide intervals. (5) COFE's transform deviates from
DESeq2 as stated. (6) A138's 517 pairs include other histologies; here 464.
**Stop rule:** one run per part. **Code:** `paper/exploratory/A265_cofe_head_to_head/scripts/{lib_a265,scout_a265,score_a265}.py`.

## A266 (2026-11-15, user: "Yes run these new amendments for the fetal adrenal gland", to the follow-ups proposed after A264; **NOT BLIND**: every fetal adrenal statistic of A189, A242, A244 and A264 has been seen (A264: without PER1 and NR1D1, C 0.489, p 0.029, above matched sets, below the responsive q95 0.634); seen before writing: the adult power below; NOT seen: any fetal statistic on the CORE waveform, any fetal relation of the responsive genes to the clock state, any fetal receptor expression): is the clock-shaped covariation of the fetal adrenal cortex an entrained clock loop or a hormone-driven response? **A244 cortex readout (34 fetuses, 6 cohorts); L1 and a gene-relation readout; secondary (development section)**
**Why.** After A264 the cortex covariation is not explained by age, cohort or steroidogenic activity, but without PER1 and
NR1D1 it is not above the glucocorticoid / immediate-early companion. A rhythmic hormonal input (fetal ACTH via MC2R /
cAMP; cortisol via GR; maternal melatonin via MT1, reported in non-human primate fetal adrenal) would move the
signal-induced clock genes and the responsive genes together. Three questions separate "driven response" from "entrained
loop" as far as these data allow.
**Fixed now (`design_a266.py`, `score_a266.py`).** Rows, covariates (cohort indicator, quadratic age), permutations
(10,000, per gene within cohort) and companions (200 matched, 200 responsive, same waveform) as A264. Seed 266.
- **H266.1** `cortex:core`: signed coherence C on the CORE waveform = the 14 GTEx clock genes without those with a
  documented direct acute response to glucocorticoid or cAMP (PER1: GRE and CRE; PER2: GRE; NR1D1: glucocorticoid
  repression): ARNTL, CLOCK, NPAS2, PER3, CRY1, CRY2, NR1D2, DBP, RORA, TEF, CIART. PASS: p < 0.05 and above the q95 of
  both companions. Descriptive: `cortex:core:steroid` (plus the A264 steroid score), `medulla:core`.
- **H266.2** `cortex:relation`: per sample, the centred clock state z from the CORE waveform (`clocksync.paired`
  projector and states); per responsive gene (`RESPONSIVE_GENES` with mean count >= 20), its centred log fraction
  regressed on (Re z, Im z): loading beta_g. Agreement with the adult adrenal (GTEx, all 229 glands, untimed, same
  computation): A = Re sum beta_g conj(beta_adult,g) / norms; p against 10,000 permutations of the fetuses' responsive
  rows (genes together) within cohort. Scored only if the adult power is >= 0.8.
- **H266.3** (descriptive) receptors: median log2 cpm of MTNR1A, MTNR1B, MC2R, MRAP, NR3C1, HSD11B2 in the fetal cortex,
  medulla and stroma pseudobulks and the adult adrenal (MTNR1A/B are absent from the GTEx matrix).
**Power first (done).** H266.1 (`feasibility_a266.py`, the A264 feasibility on these tests, 200 subsamples):
`cortex:core` adult pass rate 1.000 / 0.995 / 0.630 at 1 / 0.75 / 0.5, false-positive rate 0.055;
`cortex:core:steroid` 1.000 / 0.975 / 0.650, 0.085; `medulla:core` 0.990 / 0.910 / 0.415, 0.050. H266.2
(`score_a266.py --stage power`: 200 draws of 34 adult glands against the remaining 195, 200 permutations each, 28
responsive genes): pass rate **0.500** (median A 0.344), false-positive rate 0.075. **H266.2 is not estimable** at 34
samples by the D1b rule and its fetal statistic is not computed (a weaker design is not substituted).
**Readings (fixed now).**
- H266.1 PASS: the covariation is carried by the core loop without the signal-induced genes, above both companions:
  an entrained clock loop is supported over a purely acute hormonal response (not a running clock shown; still not
  timed). The fetal result returns to "candidate for the main text" with this qualification (to the PI, Q32).
- H266.1 p < 0.05, not clock-specific: the core loop covaries but not beyond the hormone-response genes: "driven or
  entrained, not separable"; the fetal result stays in the supplement as A264 placed it.
- H266.1 fail (p >= 0.05; adult power 1.000 at the full clock, 0.995 at 0.75): the covariation rests on PER1, PER2 and
  NR1D1: a hormone-driven response is the simpler reading; a core loop of three quarters of the adult size is excluded.
- H266.3 orders nothing; it reports whether the receptors of the candidate inputs are expressed in the fetal cortex.
**Limits.** (1) Not blind. (2) The split of the clock genes into signal-induced and core is from the literature (from
memory at writing; PER2 and NR1D1 regulation by glucocorticoid is documented mainly in rodent liver and cell lines). (3)
The responsive companion with 11 genes may have a high q95 (as in A264 with 12); the verdict is read as registered. (4)
No sampling times: entrainment cannot be shown, only covariation. (5) Limits of A244 and A264.
**Stop rule:** one scoring run (`score_a266.py --stage score`). **Code:**
`paper/exploratory/A266_fetal_adrenal_drive/scripts/{design_a266,feasibility_a266,score_a266}.py`.

## A263 (2026-11-15, roadmap item A4 under the lead prompt `prompts/biology_paper_checks_2026-11-15.md`, check 3; number reserved when the inventory started, registered after A264-A266; **NOT BLIND**: reruns of scored data (A42, A142, A150); seen: the inventory `exploratory/A263_companions/inventory.md` (no statistic), A142 (GTEx brain block r* 0.639 [0.530, 0.728] vs reference median 0.978, timed), A42 (sheep brain vs periphery coupling medians 0.205 vs 0.416, baboon 0.420 vs 0.712; r* sheep 0.625, baboon 0.970, not separate), A150 (pooled atlas coupling 0.393 vs matched q95 0.313); NOT seen: any pseudo-set statistic of the readouts below, any human brain-periphery coupling ratio): are the brain exception and the atlas couplings clock-specific against matched and stress-responsive gene sets? **GTEx, SheepGTEx, baboon; five atlases of A150; L2; primary for the organ-map and brain claims (companions)**
**Why.** A253 showed that pairing between brain regions is matched by stress-responsive sets, and A259 that matched sets
carry 0.30 to 0.74 of a tissue's apparent deviation. The brain exception and the atlas couplings were scored before the
responsive companion existed. Inventory (`inventory.md`, no statistic): no main-figure result has a gene-set companion
on the same permutations; the brain readouts and the organ map have none; the atlas couplings have a matched one only.
One rule for all: **a result is clock-specific if it exceeds the q95 (here, for a drop, lies below the q05) of both
families**. Original verdicts are not changed; a companion column is added; a result that loses specificity gets an
erratum in RP and moves in the roadmap status tables.
**Part A, the brain exception (directional; `scripts/score_a263.py`, `build_genomewide.py`).** Per species one table
(GTEx v10 fast deaths, tissues with >= 40 donors as A142: brain = the "Brain -" regions, Pituitary excluded, periphery
the rest; sheep and baboon as A42, genome-wide rows identical to the A42 builders, checked on the clock genes). For the
clock genes and for each of 200 expression-matched and 200 responsive pseudo-clock sets (`pseudo_waveforms`, GTEx
waveform amplitudes and phases; matching variable = mean log fraction over all rows of the species, genes with mean
count >= 1; seed 263):
- **R2 (primary, three species):** coupling of every tissue to the peripheral-anchor host (`Dataset.fit(n_boot=0,
  gate=False)`, anchors = periphery, peripheral tissues leave-one-out: the A42 readout 2, now also in human);
  L = log(median brain coupling / median peripheral coupling). **H263.1 (human), H263.2 (sheep), H263.3 (baboon):** the
  brain's drop is clock-specific if L_clock < q05 of L over the matched sets and < q05 over the responsive sets (sets
  with an undefined L are counted and left out).
- **R1 (secondary, human and sheep; baboon has 12 animals):** block agreement r* of brain vs periphery hosts
  (`Dataset.agreement`, n_boot = 0), clock-specific if below the q05 of both families (sets whose r* is undefined are
  counted). In human R1 is the untimed counterpart of A142's timed r*.
**Readings (fixed now).**
- H263.1 PASS and at least one of H263.2-3 PASS: "brain couples less to the body's clock state" stands as a clock
  statement, in human and the passing species.
- H263.1 PASS only: stands in human; the cross-species sentence is restated per species.
- H263.1 fail (L_clock not below the q05 of one family): brain's lower coupling is shared by non-clock gene sets of that
  family (for the responsive family: brain responds less to the body's shared stress / hormonal state); **the brain
  exception falls as a clock statement** (tally of the lead prompt), and A142 gets an erratum: "brain is less coupled
  to the body's shared state, not shown to be clock-specific".
- In any species, L_clock not negative: no drop in that readout; reported.
**Part B, atlas couplings (`exploratory/A150_atlas_specificity/scripts/compute_a150.py --kind responsive`, the A150
script with the responsive family added; outputs in `A150_atlas_specificity/output/responsive/`).** Same five atlases,
rows, anchors and fits as A150; 200 responsive sets. **H263.4:** pooled mean coupling of the 53 compartments above the
q95 of the responsive pooled sets (one-sided p = (1 + #{set >= clock}) / (1 + n)), as H150.1 for the matched family
(H150.1: PASS, 0.393 vs q95 0.313). Descriptive: per compartment, the share above both q95.
**Readings.** H263.4 PASS: the atlas couplings are clock-specific against both families. Fail: the atlas couplings are
matched by responsive sets; the organ-map coupling claims of Figure 2 are restated as "shared state", not clock (erratum).
**Not in this amendment (inventory, stated reasons):** timed validations (A36 baboon host vs ZT, A123 organ-map recovery,
A201 aorta shift) test recovery of recorded time, not specificity, and genome-wide diurnal programs would pass too;
offset equivalences (A208, A231, A235) are null results; A132b stays descriptive (user decision); A91 belongs to roadmap
T2; A227 (responsive companion missing) follows in a later amendment with the A138 template.
**Limits.** (1) Not blind. (2) R2 in human is new (A142 used timed departures); R1 is untimed. (3) The responsive pool
holds 28 to 29 genes: its sets overlap heavily. (4) Sheep symbols are NCBI symbols; baboon genes are one-to-one
orthologs with GENCODE v29 names (genes without a one-to-one ortholog cannot be drawn). (5) Pseudo sets carry the GTEx
clock waveform's amplitudes and phases, as in A150.
**Stop rule:** one run per species and one Part B run. **Code:** `paper/exploratory/A263_companions/scripts/`,
`paper/exploratory/A150_atlas_specificity/scripts/compute_a150.py`.

## A267 (2026-11-15, user: "You can" (register the COFE follow-up), with the remark that in GTEx COFE failed on all genes but worked reasonably on the core clock genes with whole-day sampling (RP A123b: organ offsets r 0.856 with the clock genes vs 0.139 genome-wide, design F; both fail in working hours); follow-up of A265; seen: A265 part (i)-(ii) on real data (genome-wide COFE pooled R_c 0.302, null mean 0.245, p 0.0013; differences uniform in 3 of 10 cancers; clocksync 0.361, p 0.0001), one COFE run on a random 14-gene set of KIRP for runtime (output not read); NOT seen: any clock-gene COFE label of TCGA, any companion of COFE, the A265 simulation summary): is the pairing of COFE's tumour and normal labels a clock signal or a patient-shared non-clock effect? **TCGA, ten cancers, 464 pairs; other; primary for the cancer chapter (head-to-head)**
**Why.** In A265 COFE's own labels pair tumour with normal beyond permuted pairs, but COFE orders the strongest circular
structure of about 1,200 genes and has no specificity companion: the pairing may come from shared patient effects
(normal contamination of the tumour, handling, batch, age, sex) rather than the clock. COFE on the clock genes is the
variant that works in GTEx (with whole-day sampling); it can be given the same gene-set companions as clocksync.
**Fixed now (`score_a267.py`; data, transform and pairing statistic of A265; seed 267).**
- **H267.1 (primary):** COFE as published (standardised features, 3 x 5-fold CV, 5 restarts) on the 14 GTEx clock genes
  of each cancer's tumours, log2(cpm + 1), no expression filter; normals projected (COFE's rule for unseen samples).
  Pooled R_c with 10,000 within-cancer permutations of the normals. Companions: COFE on 100 expression-matched and 100
  responsive 14-gene sets per cancer (`pseudo_waveforms`, matching on the mean log fraction over the cancer's tumours,
  genes with mean count >= 1; set j pooled across cancers). PASS: p < 0.05 and above the q95 of both families.
- **H267.2:** generic-axis companion of the genome-wide pairing of A265: on the A265 COFE matrices, the top 10 principal
  components of the standardised tumours (normals projected); for each of the 45 pairs of components, the angles of
  tumours and normals and their pooled R_c. PASS: the A265 genome-wide pooled R_c (0.302) above the q95 of the 45.
- Descriptive: genome-wide COFE without the 365 genes of the A265 simulation program (pooled R_c, 10,000 permutations);
  per cancer, agreement of clock-gene COFE tumour phases with the clocksync tumour state phase (|mean e^{i(phi - psi)}|,
  the better of the two directions).
**Readings (fixed now).**
- H267.1 PASS: COFE on the clock genes also pairs tumour with normal clock-specifically: both methods agree that tumour
  clocks keep a fixed relation to the patient's normal tissue; the published "uncoupled" reading is contradicted by
  COFE itself with its working variant.
- H267.1 p < 0.05, not clock-specific: COFE's clock-gene pairing is matched by non-clock sets (shared patient effects
  or a generic axis); for the clock, the evidence is clocksync's A138 (clock-specific against both families).
- H267.1 fail: clock-gene COFE does not detect pairing in surgical-hour data (consistent with A123b in working hours);
  uninformative about coupling.
- H267.2 fail: the genome-wide COFE pairing of A265 is no larger than what generic expression axes give: it is patient
  similarity, not evidence of clock coupling (A265 H265.1 is then read as "COFE's genome-wide labels carry patient-shared
  structure", in RP). H267.2 PASS: COFE's genome-wide axis pairs more than generic axes; not shown to be clock.
**Limits.** (1) Projection of the normals is an assumption (A265). (2) COFE's sparsity CV on 14 features has a short
grid (s from 1 to about 3.7). (3) 100 sets per family (runtime: about 36 s per 14-gene COFE run). (4) The responsive
pool holds about 29 genes; its sets overlap. (5) Principal-component axes are a generic companion, not a model of a
specific confound.
**Stop rule:** one run. **Code:** `paper/exploratory/A265_cofe_head_to_head/scripts/score_a267.py`.

## A268 (2026-11-16, roadmap item T4 under the lead prompt `prompts/biology_paper_checks_2026-11-15.md`, check 4; reader rule chosen by the user ("Relaxed rule") after the scout `exploratory/scouts_2026-11-15/T4_systemic_signal/`; seen: the scout (mouse candidate table; GTEx power with inserted relations), A27 D6b (GC score R^2 on the host 0.014 in fast deaths, 0.234 in ventilated donors), A136, A248, the MuTHER power below; NOT seen: any observed reader-host relation in GTEx or MuTHER): does a glucocorticoid or temperature reader keep a fixed phase relation to the host clock state across tissues and people, clock-specifically, confirmed in living people, and is the relation looser in brain? **GTEx fast deaths (15 host tissues), MuTHER fat and skin; other; primary for T4 (the systemic signal)**
**Readers (fixed; deviation from the prompt's rule, chosen by the user).** The prompt's rule (rhythmic without a local
clock in mouse liver, Manella 2021 hepatocyte Bmal1 knockout) leaves the glucocorticoid and cold readers empty (scout).
Relaxed rule: literature candidates rhythmic (q < 0.05) in control mouse liver at the expected phase: glucocorticoid
FKBP5 + TSC22D3 (ZT 15.0, 16.2; after the corticosterone peak), cold CIRBP (ZT 6.4; cool rest phase), heat HSP90AA1 (ZT
20.2; warm active phase). Reader value per sample: mean over its genes of the within-tissue z of log expression. The 14
clock genes are never in a reader.
**Host (same construction in both cohorts).** Clock-plane state of the GTEx waveform (unweighted projector on
gene-centred log expression), centred within tissue; the host of a (person, tissue) = the mean state of the person's
other tissues (leave-one-tissue-out). Variant: waveform without PER1 (a glucocorticoid target); every result is reported
with and without it.
**Statistic.** Per tissue k, the complex loading beta_k = sum y conj(z) / sum |z|^2 of the centred reader on the centred
host; consistency K = |sum_k beta_k| / sum_k |beta_k| over the 12 body tissues (one phase relation in every tissue: K
near 1); one-sided p against 10,000 permutations of the reader across persons within tissue. Companions on the same
host: 200 readers of expression-matched genes (the reader's size, `pseudo_clock_genes` bins on the GTEx mean log CPM,
CLOCK_LIKE excluded) and 200 of responsive genes (`RESPONSIVE_GENES` minus the reader's genes).
- **H268.1-3 (glucocorticoid, cold, heat; Holm over the three, per host variant):** PASS = Holm p < 0.05, K above the
  q95 of both companions, AND living confirmation: in MuTHER (fat on the skin host, skin on the fat host), the pooled
  loading |beta_fat + beta_skin| has p < 0.05 against 10,000 permutations of persons within tissue, and its argument
  lies within +-3 h of the GTEx body argument arg(sum_body beta).
- Verdict classes: "PASS (confirmed)"; "GTEx only (not confirmed)" (does not count as a positive); "not specific"
  (Holm p < 0.05, not above both companions); fail.
- Descriptive: brain looser: projection of the cerebellum and cortex loadings on the unit vector of sum_body beta,
  relative to the median body projection, 200 donor bootstraps.
**Power first (done).** GTEx (scout, clocksync host; an inserted relation of fixed phase): pass rate 0.99-1.00 at
R^2 0.014 (the A27 GC level), 1.00 at 0.03; false-positive rate 0.02-0.07. MuTHER confirmation
(`score_a268.py --power`, `output/power_muther.csv`, 200 draws, 1,000 permutations): glucocorticoid 0.96 at R^2 0.014,
1.00 at 0.03, false-positive rate 0.04; cold 0.975, 1.00, 0.065; heat 0.965, 1.00, 0.06.
**Readings (fixed now).**
- A reader PASS (confirmed): that signal keeps a fixed phase relation to the body's clock state in dead and in living
  people, beyond matched and stress-responsive genes: the named candidate for the systemic signal (correlational: a
  marker of the signal, not shown to set the clocks). With PASS in the host without PER1 as well, the relation does not
  rest on the GC-responsive clock gene.
- GTEx only: a relation in post-mortem tissue that living tissue does not confirm: not counted (agonal or post-mortem).
- Not specific: the relation is shared by matched or stress-responsive genes: a general state, not this signal.
- All three fail at the stated power: a clean negative: these public data cannot name the signal by these readers.
- Brain (descriptive): a ratio below 1 with its interval below 1: the relation is looser in brain.
**Limits.** (1) Relaxed reader rule (user decision); the mouse check is in control liver, where the local clock may
contribute. (2) Agonal stress drives glucocorticoid and heat-shock genes and can move PER1 in GTEx: the responsive
companion and the host without PER1 address this, the living confirmation decides. (3) MuTHER has two tissues and no
sampling times; a confirmation is of the phase relation, not of time. (4) Readers of one or two genes.
(5) Extension to GTEx long-term steroid users (protected phenotype, aggregates only) only after a confirmed positive.
**Stop rule:** one scoring run. **Code:** `paper/exploratory/A268_systemic_signal/scripts/score_a268.py`.

## A269 (2026-11-16, user: "yes go ahead", to the follow-up proposed after A243 ("how does lung liver compare to adrenal gland, which might be driven only by hormonal response that don't oscillate with a period of 24h"); **NOT BLIND**: every statistic of A243 has been seen (H243.1 PASS: lung:all C 0.530, Holm p 0.026; liver:all C 0.583, Holm p 0.025; leave-one-gene-out minimum lung 0.439 without NPAS2, liver 0.499 without CRY1; liver without Cardoso-Moreira 0.249, p 0.19), and A264 / A266 on the adrenal (core loop C 0.417, p 0.057); seen before writing: the GC-gene medians per cohort (detection only) and the adult power below; NOT seen: any fetal statistic on the CORE waveform or with the GC score): is the clock-shaped covariation of the fetal lung and liver carried by the core loop without its signal-induced genes, and does it survive adjustment for glucocorticoid-target activity? **A243 `lung:all` (38 fetuses, 5 cohorts) and `liver:all` (35, 3); L1; secondary (development section)**
**Why.** Signed coherence C measures a clock-shaped covariation between fetuses, not a 24 h period: a hormone level
that differs between fetuses and moves the signal-induced clock genes would give it as well as a running clock. In the
adrenal cortex the covariation did not survive the removal of PER1, PER2 and NR1D1 (A266) and a hormone-driven
response became the simpler reading. Glucocorticoids act on the fetal lung (maturation) and liver; the same two checks
decide whether lung and liver carry a loop the adrenal lacks, or the same driven response.
**Fixed now (`design_a269.py`, `feasibility_a269.py`, `score_a269.py`).** Rows, covariates (cohort indicator, common
quadratic age), permutations (10,000, per gene within cohort) and companions (200 matched, 200 responsive, adjusted in
the same way, on the same waveform) as A243. Seed 269.
- CORE: the A266 waveform (the 14 GTEx clock genes without PER1, PER2, NR1D1; 11 genes). GC score: per fetus, the mean
  over FKBP5, TSC22D3, ZBTB16, KLF9 (canonical direct glucocorticoid-receptor targets, not clock genes; fixed before any
  statistic) of the gene's log fraction standardised over the fetuses of the organ (`design_a264.activity`); all four
  genes are measured in every lung and liver cohort.
- **H269.1** `lung:core`, `liver:core`: C on CORE. Holm over the readouts scored.
- **H269.2** `lung:gc`, `liver:gc`: C on the 14 genes, the GC score added as a covariate. Holm over the readouts scored.
- PASS (`development.verdicts`): Holm p < 0.05 and above the q95 of both companions; "not clock-specific" and
  "suggestive, not established" as A243.
- Descriptive: `<organ>:core:gc` (both), `lung:early:core` (13 fetuses before 12 weeks, quadratic age).
**Power first (done; `feasibility_a269.py`, 200 subsamples, 500 permutations).** The A243 stand-ins (GTEx Lung, Liver;
deaths 08-18 h; thinned to each cohort) also receive the adult donor's own GC score (unthinned), so the adjustment
removes from the adult clock what the adult score shares with it:

| readout | fetuses | genes | covariates | adult pass rate at 1 / 0.75 / 0.5 | without a clock |
|---|---|---|---|---|---|
| `lung:core` | 38 | 11 | 6 | 1.000 / 1.000 / 0.840 | 0.065 |
| `lung:gc` | 38 | 14 | 7 | 1.000 / 0.960 / 0.340 | 0.030 |
| `lung:core:gc` | 38 | 11 | 7 | 1.000 / 1.000 / 0.865 | 0.105 |
| `liver:core` | 35 | 11 | 4 | 1.000 / 1.000 / 0.910 | 0.070 |
| `liver:gc` | 35 | 14 | 5 | 1.000 / 1.000 / 0.445 | 0.045 |
| `liver:core:gc` | 35 | 11 | 5 | 1.000 / 1.000 / 0.855 | 0.040 |
| `lung:early:core` | 13 | 11 | 3 | 0.975 / 0.845 / 0.505 | 0.045 |

Adult GC score (`adult_scores.csv`): lung (208 donors) R^2 on the hour of death 0.0041, mean pairwise r of the four genes 0.599; liver (106 donors) R^2 on the hour of death 0.0146, mean pairwise r of the four genes 0.583.
A readout is scored only if its adult pass rate at the full clock is >= 0.8 (D1b rule).
**Readings (fixed now).**
- H269.1 PASS in an organ: the covariation is carried by the core loop without the signal-induced genes, above both
  companions: the organ carries what the adrenal cortex did not (A266); "an entrained loop is supported over a purely
  acute response" for that organ (not a running clock: no times). The developmental line becomes "lung / liver ahead of
  the adrenal on this criterion" (to the PI, Q32).
- H269.1 fail (p >= 0.05): with the adult pass rates above (>= 0.84 at 0.5), a core loop of half the adult size is
  excluded, i.e. a core loop of the size of the 14-gene signal of A243 (between 0.5 and 0.75); the organ reads as the
  adrenal: a response through PER1, PER2, NR1D1 is the simpler reading. If both organs fail, the
  developmental result narrows to "the clock genes covary before birth in the adrenal, lung and liver; a core loop is
  not shown in any".
- H269.1 p < 0.05, not clock-specific: "driven or entrained, not separable", as A266.
- H269.2 PASS: the covariation is not explained by glucocorticoid-target activity. Fail: it is shared with that
  activity. Neither direction orders cause: in adults GC targets are also clock outputs, so a fail does not show that a
  hormone drives the clock genes; it is read together with H269.1.
**Limits.** (1) Not blind. (2) CORE and GC are literature lists (from memory at writing; PER2 / NR1D1 regulation by
glucocorticoid documented mainly in rodent liver and cell lines). (3) The GC score is four genes: a coarse readout of
glucocorticoid action, not of hormone level; other inputs (cAMP, temperature, maternal melatonin) are not adjusted.
(4) Fetal liver is largely haematopoietic at these ages (A243 limit 7); liver rests on pooling (no cohort alone is
significant). (5) The responsive companion is redrawn on CORE (11 genes) and its q95 may be high, as in A264 and A266.
(6) Limits of A243.
**Stop rule:** one scoring run. **Code:** `paper/exploratory/A269_fetal_organs_drive/scripts/`.

## A270 (2026-11-16, user: "yes go ahead", to the test proposed after A263 Part B ("most cells should share a clock with the internal time"); **NOT BLIND** on the clock: the clock-gene host vs recorded time is known in both cohorts (baboon A36 R 0.987; GTEx fast deaths, `walltime_hosts.csv`, R 0.631); seen: A263 Part A (human, baboon) and Part B (H263.4 FAIL, clock 0.393 at the 87.5th percentile of the responsive pooled sets); seen before the full run: the baboon `--check` run (3 responsive and 3 matched sets, 200 permutations; for the code only); NOT seen: any other host of a pseudo gene set read against recorded time): do the responsive companion sets carry internal time? **Baboon (GSE98965, 12 animals, 63 tissues) and GTEx fast deaths (Hardy 1-2, 15 reference tissues); control; primary for the reading of every responsive-companion verdict (A150/A263 Part B, A263 Part A, A269)**
**Why.** The responsive family (`controls.RESPONSIVE_GENES`: glucocorticoid targets and immediate-early genes) was used as
a non-clock control. Glucocorticoid targets follow the daily cortisol rhythm in vivo, so a responsive set may be a
systemic output of internal time rather than a non-clock state. If so, failing against it says the clock-gene readout is
not separable from the systemic daily program, not that the shared state is non-clock.
**Design.** Per species, 200 responsive and 200 expression-matched sets (`clocksync.paired.pseudo_waveforms`, GTEx
waveform amplitudes and phases, matching variable = mean log fraction over the rows used, genes with mean count >= 1;
seeds [270, 1] and [270, 2]) plus the clock. Per set: the untimed host (`Dataset`, every tissue an anchor, as A36 and
`compute_walltime_hosts.py`) read against recorded time with `Dataset.timed` (R, its permutation p over 1000 shuffles of
the times, Fisher-Lee r). Baboon: ZT; human: time of death (fast deaths, the 15 tissues of `V2_TISSUES`, depth = library
size; aggregates only).
**H270.1 (baboon), H270.2 (GTEx):** the share of responsive sets with p_R < 0.05 is >= 0.5 ("the responsive family carries
internal time"). Descriptive: the same share for matched sets, the median R of each family, the percentile of the clock R
in each family, and the median of max(R, R of the reflected host) (a set whose genes peak in another order can read time
mirrored).
**Readings.** H270.1 and H270.2 PASS: the responsive companion is an output of internal time; the verdicts against it
(A150/A263 Part B, A263 Part A, A269) read "not separable from the systemic daily (glucocorticoid) program", not "shared
non-clock state"; specificity of a clock statement then rests on the matched family. Both FAIL: the responsive sets do not
read time; the failures stand as non-clock warnings. Split: per species. If matched sets also pass in a species, internal
time is carried genome-wide there (Mure 2018: most baboon genes are rhythmic in some tissue) and is reported as such.
**Limits.** (1) Not blind on the clock. (2) The pseudo sets carry the clock waveform's phases; a responsive set whose genes
peak in another order is read with a distorted host (hence the reflected descriptive). (3) 12 animals: permutation p has a
floor near 1/1000 and R is noisy per set. (4) GTEx time of death is recorded wall time; the clock R is 0.631.
**Stop rule:** one run per species. **Code:** `paper/exploratory/A270_responsive_time/scripts/score_a270.py`.

## A271 (2026-11-16, user: "Yes go ahead", to the responsive-companion check proposed after A270 (A227 T/NK positive, A104); **NOT BLIND**: every matched-family statistic is known (A227: H227.1 myeloid 0.581 vs matched q95 0.548; H227.2 T/NK 0.539 vs 0.510; A228 H228.4 CD8-negative T 0.532 vs 0.501; A104 periphery rho +0.734, brain -0.115), and A263 Part B / A270 (responsive sets match the clock on atlas coupling and do not read GTEx time of death); NOT seen: any responsive-set statistic of A227, A228 or A104): are the T/NK and myeloid tumour-normal pairings and the amplitude-coupling relation specific against the responsive companion? **Single-cell TME (34 cohorts) and GTEx fast deaths; control; primary for the immune lines of the cancer chapter and the A104 panel of Figure 2**
**Why.** A227/A228 are clock-specific against matched genes only (listed as missing in A263, "Not in this amendment").
Tumour and normal samples of a patient share one operation (anaesthesia, often dexamethasone), a shared glucocorticoid
state that the responsive family captures. A104: coupling c = z2 V / (z2 V + 2 psi) and own amplitude A = sqrt(z2 V + 2 psi)
share z2 V, so rho(A, c) > 0 may hold for any gene set. Figure 3 (H109.1-R) and TCGA pairing (H138.4-R) already passed
the responsive companion (2026-10-28) and are not retested.
**Part A (`compute_a227.py --stage tests --kind responsive`, `compute_a228.py --stage tests --kind responsive`; the
registered scripts with the responsive family passed to `pseudo_one`, otherwise unchanged: same cohorts, units, statistic,
pooling, seeds; outputs in `output/responsive/`).** **H271.1 (Holm over 3):** for H227.1 (myeloid), H227.2 (T/NK) and
H228.4 (CD8-negative T), the pooled pairing is above the q95 of the 200 pooled responsive sets, with p_resp = (1 + #{set >=
clock}) / (1 + n) Holm-adjusted < 0.05. The other A227 / A228 rows are descriptive.
**Part B (`exploratory/A271_responsive_companions/scripts/score_a271_a104.py`).** A104 rows, blocks and fit; 200 matched
and 200 responsive sets (seeds [271, 1], [271, 2]); rho per block, point estimates. **H271.2:** periphery rho of the clock
above the q95 of both families. Brain block descriptive.
**Readings.** H271.1 per compartment: PASS, the pairing is clock-specific against both families (T/NK: the first
lymphocyte line clock-specific against both); fail, "the patient's tumour and normal immune cells share a state" without
the clock. H271.2 PASS: "stronger clocks are more coupled" is a clock property; fail: generic to any waveform-shaped gene set
(the shared z2 V term), and the A104 panel leaves Figure 2 or is reworded as a property of the estimator.
**Limits.** (1) Not blind on the matched family. (2) Responsive sets overlap heavily (Jaccard 0.79 GTEx, 0.91 TCGA in the
A109 note); per cohort, single-cell pools can be smaller. (3) p floor 1/201. (4) A270: in GTEx the responsive shared state
is not time of death; what it is in the surgical samples is not tested here.
**Stop rule:** one run per script. **Code:** the two scripts above and
`paper/exploratory/A271_responsive_companions/scripts/score_a271_a104.py`.

## A272 (2026-11-16, user: "Okay, start this" (development of the cellular backbone: coverage, circularity, who follows whom); **NOT BLIND** for the anchored readout: the A150 couplings by family against both companion families (fibroblast above both; endothelial, epithelial, myeloid not above the responsive sets; lymphocytes at 0.069) and A271 Part A (tumour-normal pairing: fibroblast and endothelial above the responsive sets, immune cells not) are known; NOT seen: any pairwise statistic between compartments of a healthy atlas, on any gene set; seen before registration: the design stage of the script (pairs, donors, degrees of freedom; no state computed)): is the stromal-vascular backbone seen without a host? **Five healthy single-cell atlases; anchor-free; control of the backbone claim; primary H272.1, H272.2**

**Why.** In the anchored fits the host is defined on stromal-vascular compartments (learned anchors), so their high
coupling is partly by construction, even with the leave-one-out readout (the other stromal-vascular anchors stay in the
host). The backbone claim ("a person's shared clock state is carried by the stromal and vascular cells of a tissue,
lymphocytes do not carry it") needs a readout with no host and no anchors, against both companion families.
**Data.** The A150 atlases (gate-passing A17 healthy atlases with a genome-wide twin): lung_hlca, heart, organ_breast,
organ_muscle, organ_colon; rows, strata and families as in A150 (`anchor_rule.family`). Design
(`exploratory/A272_backbone_pairwise/output/design.csv`): 249 pairs of compartments; pairs (stromal-vascular, X) 60 in 5
atlases, pairs (lymphocyte, X) 60 in 4 atlases (muscle has no lymphocyte compartment), X = any family other than
stromal-vascular and lymphocyte (epithelial, myeloid, parenchyma, neuron, other).
**Statistic.** One projector per atlas and gene set (`clocksync.paired.projector`, `states`). For a pair (a, b), over the
donors having both, in strata with >= 3 shared donors, pooled degrees of freedom >= 15: states centred within stratum,
r_ab = sum a conj(b) / sqrt(E_a E_b), with E the summed squared states minus the Poisson variance (the pooled form of
`paired.ccorr`); a pair with a non-positive corrected variance is NaN and left out of medians. |r_ab|: shared state;
arg r_ab: offset (h, a minus b). T_SV = mean over atlases of the median |r| over pairs (stromal-vascular, X); T_LY the same
for (lymphocyte, X); D = T_SV - T_LY over the 4 atlases with lymphocytes. Companions: 200 expression-matched and 200
stress-responsive pseudo-clock sets per atlas (`paired.pseudo_waveforms`, matched on the mean log fraction over the rows
of the atlas, genes with mean count >= 1; seeds [272, atlas, 1 | 2]); set j is pooled across atlases as set j.
**H272.1 (the backbone without a host):** T_SV of the clock genes is above the q95 of both companion families.
**H272.2 (the ranking is a clock property):** D of the clock genes is above the q95 of D in both companion families: the
same partner cell types share more clock state with stromal-vascular cells than with lymphocytes, beyond what
non-clock genes show. No multiplicity correction between the two (H272.2 is read only with H272.1's outcome stated).
**Descriptive.** D272.1 family x family matrix of median |r| (clock, companion medians, share of pairs above both q95) and
per-compartment centrality (median |r| with compartments of other families) with its excess over each companion family;
Spearman of that excess with the A150 anchored excess. D272.2 offsets: for pairs above the q95 of both families, arg r
in hours with a donor bootstrap within stratum (500, 90 %), and the closure of triangles (offset ab + bc - ac). D272.3
sensitivity of T_SV and D: correlation without the Poisson correction; clock counts after the A150 ambient
over-correction.
**Readings.** H272.1 and H272.2 PASS: the backbone holds without a host and against both companions; the section is
written on the pairwise readout, the anchored couplings become the secondary display. H272.1 PASS, H272.2 fail: stromal-
vascular cells share a clock state with the other cells, but lymphocytes are not shown to share less of it than non-clock
genes predict (depth or responsiveness may explain their low coupling): "lymphocytes do not carry it" is dropped.
H272.1 fail against the responsive family only: the shared state of stromal-vascular cells is not separable from the
daily or acute hormone response in these atlases (A270: such sets read time in the baboon, hardly in GTEx): reported as
such, backbone kept as "shared state". H272.1 fail against both: the backbone was an effect of the anchoring; the
section is withdrawn.
**Not in this amendment.** No lead or lag between cell types is interpreted: D272.2 gives the precision only. A
driver-follower reading (lag-amplitude relation, mediation) needs first a calibration of pairwise offsets against timed
truth (mouse colon with known ZT, aorta) and is a separate amendment; offsets below about 1 h are not interpreted (the
project's untimed-vs-timed errors are 0.4 to 0.5 h). Extension to the atlases without a genome-wide table needs new
builds from the Census and is decided after this result.
**Limits.** (1) Compartments of one sample share ambient RNA and dissociation effects, which inflate |r| and pull offsets
to 0; the companions carry the same inflation for |r|, D272.3 bounds it for the clock. (2) The Poisson correction is
unstable at low counts (lymphocytes): NaN pairs are counted, medians used, the uncorrected value reported. (3) Five
atlases, 4 with lymphocytes; lung_hlca contributes the most pairs (per-atlas medians, then the mean, limit its weight).
(4) |r| is positively biased at small n; the companions share the donors, so the comparison is unaffected.
**Code.** `exploratory/A272_backbone_pairwise/scripts/score_a272.py` (`--stage design|score`, `--placebo`: donors permuted
within compartment and stratum, a code test run before scoring).

### A272 design addendum (2026-11-16, before scoring; from the placebo run only)
The placebo run (donors permuted within compartment and stratum, 40 sets per family; pairing destroyed, so no real
pairwise statistic was produced) showed that the Poisson-corrected variance is non-positive in 144 of the 249 clock
pairs and in most pseudo sets (median number of valid matched sets per pair 2 to 29 of 40): the 1 / count variance of
`paired.states` over-corrects at single-cell pseudobulk counts, and the registered statistic is not computable.
**Change, fixed before scoring:** the primary statistic of H272.1 and H272.2 is the complex correlation without the
Poisson correction (bounded by 1, defined for every pair); the Poisson-corrected value becomes the descriptive column of
D272.3. Consequence: |r| is attenuated by counting noise in shallow compartments (lymphocytes most). The matched
companion sets are matched on expression, so they carry the same attenuation: the comparisons of H272.1 and H272.2
(clock against companions on the same compartments) are unaffected, absolute values of |r| are lower bounds and are not
compared between cell families without the companions. Everything else is unchanged.

## A273 (2026-11-16, user: "Okay, work on this autonomously" (a cell-type level understanding of the clock); a calibration of the A272 pairwise readout, no biological hypothesis; **NOT BLIND** for the per-class timed phases of this atlas (A21, A31: per-class gene phases against ZT R 0.98; epithelial shift after injury); NOT seen: any pairwise offset between cell classes of this atlas, timed or untimed; seen before registration: the design stage (mice and compartments per arm)): are untimed pairwise offsets between cell types of one tissue sample unbiased? **Mouse colon chrono-atlas (known ZT); calibration; V273.1, V273.2**

**Why.** A272 reports offsets between cell types as the argument of a pairwise correlation. Compartments of one sample
share ambient RNA and dissociation effects, which would pull such offsets towards 0; before any lead or lag between
cell types is read in the human atlases, the readout is checked where time is known.
**Data.** `colon_chrono` (clock genes, CellBender counts), rows with >= 30 cells, compartments with >= 7 of the 8 mice of
an arm. Arms C-DC (undamaged) and R-DC (regenerating) primary, 12 compartments and 66 pairs each; R-PC (10
compartments, 45 pairs) as sensitivity. ZT 1 / 7 / 13 / 19, two mice per ZT and arm.
**Readouts, per arm.** Timed: `Dataset.timed_loadings("zt")`, the phase of each compartment against ZT and its coupling
to the known time; the timed offset of a pair is the difference of the two phases. Untimed: the argument of the A272
pairwise statistic over the mice (no time). Error = untimed minus timed, with a leave-one-mouse-out jackknife SD.
Primary set: pairs of C-DC and R-DC in which both compartments have a timed coupling >= 0.5 ("timed-rhythmic").
**V273.1:** median |error| over the primary set <= 1.0 h. **V273.2:** slope of the untimed on the timed offsets (through
the origin) between 0.7 and 1.3.
**Readings.** Both pass: pairwise offsets between cell types of one sample are unbiased at the 1 h scale; the human
offsets of A272 are limited by their own intervals only. V273.2 below 0.7: offsets are shrunk towards 0 (shared RNA):
human pairwise offsets are lower bounds, and a lead or lag is read on the anchored fit only. V273.1 fails with V273.2
passing: unbiased but imprecise at 8 animals; no conclusion on bias beyond the slope.
**Descriptive.** The timed offsets themselves (median and largest |offset| between cell classes of the colon): whether
cell types of one tissue differ in clock phase when time is known; per arm; all pairs including weakly rhythmic
compartments; the relation of |error| to |r|.
**Limits.** Eight mice per arm with full 24 h coverage: this tests bias, not the precision reached in a narrow
sampling window (human hospital hours), which the donor bootstrap of A272 and the timed human comparisons (A251, A123)
cover. The human waveform is used on mouse data, as in A21 and A31. Mouse, one tissue.
**Code.** `exploratory/A273_pairwise_offset_truth/scripts/score_a273.py`.

## A274 (2026-11-16, user: "Okay, work on this autonomously" (cell-type level understanding of the clock); staged: hypotheses, statistic, eligibility and readings are fixed here, before the genome-wide tables of the new atlases exist; the design table (atlases, pairs) goes in an addendum committed before scoring; seen: A272 in full (5 atlases), A273, A150, A271; NOT seen: any genome-wide table or pairwise statistic of the atlases below): does the cell-type ranking of A272 replicate in atlases it did not use, and which cell classes share the clock state in the human cortex? **Healthy single-cell atlases not in A272; anchor-free pairwise readout; replication (H274.1, H274.2) and cortex (H274.3)**

**Why.** A272 rests on five atlases, one of them weak (colon), and was not blind for the anchored readout. The atlases
with a clock-gene table only are being rebuilt with all genes (Census 2025-11-08, same grouping as their panel builds:
`exploratory/celltype_map_builds/run.sh`, `celltype_map_builds_gut/run.sh`), which makes companions possible there.
**Atlases.** Body: organ_optic_nerve, organ_eye_anterior, kidney_sc, islets, kong_ileum, kong_colon,
lung_natri_control. Cortex: psychad_hbcc, seaad_dlpfc, seaad_mtg. Rows, strata, classes and healthy donors: those of the
clock-gene dataset of each atlas (A17 healthy target). **Eligibility, fixed:** an atlas enters only if its genome-wide
table reproduces the clock-gene table on the twin check (clock counts, cells, depth, strata identical on the common
rows); pairs as in A272 (strata with >= 3 shared donors, pooled degrees of freedom >= 15). An atlas that fails the twin
check or is not built is left out and named; no atlas is added.
**Statistic.** As A272 after its design addendum (complex correlation without the Poisson correction, one projector
per atlas and gene set; 200 expression-matched and 200 stress-responsive sets, seeds [274, atlas, 1 | 2]; set j pooled
across atlases as set j). Families from `anchor_rule.family`; SV = fibroblast_stromal, endothelial, mural; LY =
lymphocyte; X = any other family.
**H274.1 (replication of H272.2, primary):** D = T_SV - T_LY over the eligible body atlases that have at least 3 pairs
(SV, X) and 3 pairs (LY, X), above the q95 of D in both companion families.
**H274.2 (replication of the level):** T_SV over the eligible body atlases with at least 3 pairs (SV, X), above the q95
of the matched family; the responsive family is reported with the same rule (A272: p 0.149).
**H274.3 (cortex):** T_CX = mean over the eligible cortex atlases of the median |r| over pairs of compartments of
different families (neuron, glia, endothelial, myeloid), above the q95 of both companion families: the cell classes of
the cortex share a clock-gene state beyond non-clock genes. Descriptive: the family x family matrix with companions;
neuron-neuron against neuron-glia against glia-glia pairs; offsets with donor-bootstrap intervals for pairs above both
companions (A273: unbiased at the 1 h scale), and their agreement between the three cortex atlases.
**Readings.** H274.1 PASS: the ranking (stromal-vascular cells share the clock state with the tissue, lymphocytes do
not) is replicated in independent atlases and is written as a result. Fail with the same sign and the clock D above the
companion medians: not replicated at this power, reported with A272 as the only support. Fail with D at or below the
companion medians: the A272 ranking is not general; the section states it for the A272 atlases only, as a lead.
H274.2 as H272.1. H274.3 PASS: a shared clock state links the cortical cell classes (the brain is internally coherent,
whatever its coupling to the body, A263); fail against the responsive family only: shared state not separable from the
agonal or hormonal response (post-mortem tissue); fail against both: no cell-class clock relation is read in cortex.
**Not in this amendment.** No driver-follower reading. Gene-level outputs per cell type (A50 maps) are not rescored.
**Limits.** Post-mortem cortex (agonal state shared by all classes of a donor: the responsive family is the relevant
control); nuclei, not cells, in the cortex atlases; eye and optic nerve come from post-mortem donors too; the body
atlases with lymphocytes may be only the two gut atlases and one or two others (few pairs: the addendum gives them).
**Code.** `exploratory/A272_backbone_pairwise/scripts/score_a272.py` generalised to a list of atlases and a family
grouping (`exploratory/A274_celltype_map/scripts/score_a274.py`), placebo before scoring.

## A275 (2026-11-16, autonomous follow-up of A272; **post hoc and descriptive**: A272 is known in full (the responsive sets reach T_SV 0.434 against the clock's 0.499; lymphocytes with the tissue: clock 0.220, responsive 0.349); NOT seen: any per-gene sharing): which genes are shared between the cell types of a donor: clock genes, glucocorticoid targets, or immediate-early genes? **The five A272 atlases; per-gene decomposition of the responsive companion; no test**

**Why.** The responsive companion pools two kinds of genes (`controls.GC_TARGETS`: glucocorticoid targets, several of
them clock outputs such as NFIL3, BHLHE40, KLF9, TSC22D3; `controls.IMMEDIATE_EARLY`: FOS, JUN, EGR1 and the like, the
core of the tissue-dissociation response). Failing against it means different things for the two: a daily hormonal
program shared by the cells of a donor is internal time by another route; an immediate-early response shared by the
cells of one sample is a handling state. No separate dissociation companion is added: the immediate-early pool is that
signature.
**Readout.** Per gene and pair of compartments (the A272 pairs): correlation over donors of log((count + 0.5) / depth),
centred within stratum, pooled over strata. Sharing of a gene in a pair class (SV_X, LY_X, ...) = median over pairs;
of a gene class (clock, glucocorticoid target, immediate-early, expression-matched genes of 200 matched sets) = median
over genes per atlas, then mean over atlases.
**Reported, no verdict.** The sharing of the four gene classes in (SV, X) and (LY, X) pairs, and the per-gene values
of the clock, glucocorticoid and immediate-early genes. **Reading rule fixed in advance:** if immediate-early genes
are shared at least as much as the clock genes while glucocorticoid targets are not, the responsive level of A272 is
a handling state and the matched family is the relevant control for the clock; if glucocorticoid targets carry it, the
responsive level is a systemic hormonal state, not separable from internal time in these data (A270); if both, both.
Per-gene sharing is not corrected for counting noise: genes are compared with expression-matched genes only.
**Code.** `exploratory/A275_shared_gene_classes/scripts/score_a275.py`.

### A274 design addendum (2026-11-16, before any pairwise statistic of these atlases; from `--stage design`, which computes eligibility and pair counts only)
- **Eligible (genome table reproduces every clock count of the clock-gene table, all rows):** body: `organ_optic_nerve`
  (83 donors), `organ_eye_anterior` (67), `kidney_sc` (26), `kong_ileum` (30), `kong_colon` (21); cortex: `psychad_hbcc`
  (145), `seaad_dlpfc` (44), `seaad_mtg` (47).
- **Not eligible, by the registered rule:** `islets` (3 of 118 rows differ from the clock-gene table, largest
  difference 17 counts; cause not traced); `lung_natri_control` (no exact twin: its T class is restricted by a
  per-cell cycling label that the genome-wide group table does not carry). Neither is approximated.
- **Pairs per class (body):** SV_SV 23 (optic nerve 13, eye 8, kidney 1, ileum 1), SV_LY 20 (ileum 14, colon 3, optic
  nerve 3), LY_LY 31 (ileum 21, colon 10), SV_X 92, LY_X 50, X_X 66. H274.1 (D = T_SV - T_LY) therefore rests on the
  gut for its lymphocyte side and on the eye tissues for its stromal-vascular side: stated here as a limit, the
  statistic is unchanged. **Cortex:** CX_diff 71 (27 + 22 + 22), CX_same 21.
- **Dependence:** `seaad_dlpfc` and `seaad_mtg` are two regions of the same cohort and share donors; they count as
  one cohort in the reading of H274.3 (replication = PsychAD and SEA-AD agree), and per-atlas values are reported.
- Source: `exploratory/A274_celltype_map/output/{twin,design}.csv`. Next: placebo, then one scoring run.

## A276 (2026-11-16, user: "you can work autonomously on paper/docs/prompts/celltype_clock_program_2026-11-16.md" (block A1); a calibration of the A272 pairwise test, no biological hypothesis; seen: A272 and A274 in full, including their single placebo runs (A274: all six placebo clock values above the companion medians, one of six at its q95, p 0.050); NOT seen: any placebo draw beyond those two): is the pairwise test (clock against the q95 of both companion families) calibrated when the pairing of donors is destroyed? **The thirteen atlases of A272 and A274; calibration; V276.1, V276.2**

**Why.** H272.2, H274.1-3 compare a clock statistic with the q95 of 200 companion sets. Under no pairing the clock
statistic should fall above that q95 in about 5 % of placebo draws. Heavier tails or a different shape of the clock-plane
states than the companion states would make the clock's null |r| larger and the test anti-conservative; one placebo run
per amendment cannot show this.
**Design.** Rows, atlases, pairs, families, statistics and companion sets exactly as scored (A272: lung_hlca, heart,
organ_breast, organ_muscle, organ_colon, seeds [272, atlas, kind]; A274: the eight eligible atlases of its design
addendum, seeds [274, atlas, kind]; 200 matched and 200 responsive sets each). States are computed once per row and gene
set (`paired.projector`, `states`: row-wise, so unchanged by relabelling donors). Placebo draw d (seed [276, d],
d = 1..500): donor labels permuted within compartment x stratum in every atlas (`score_a272.placebo`), pairs rebuilt
(`pair_index`), the ten statistic-family combinations recomputed: H272.1 T_SV, H272.2 D, H274.1 D, H274.2 T_SV, H274.3
T_CX, each against the matched and the responsive family. Draw 0 (no permutation) must reproduce the scored
`tests.csv` of A272 and A274 (code check).
**V276.1 (each test):** the rejection rate (clock above the companion q95) over 500 draws is <= 0.075 for each of the
ten combinations (binomial 99 % upper limit for a true rate of 0.05 is 0.073).
**V276.2 (overall):** the mean rejection rate over the ten combinations is <= 0.065.
**Descriptive.** Mean percentile of the clock among its companions per combination (0.5 expected); the same rate per
atlas and pair class (median |r| over the class, clock against the class q95), which A2 uses to read single atlases;
the clock's null |r| against the companions' null |r| per pair.
**Readings.** Both pass: the test is calibrated, H272.2 and H274.1-3 stand as scored, and the readout moves to
`clocksync`. A rate above the band: the scored verdicts of that statistic are suspended (not reversed) and the cause is
traced before any further test of this line is scored: the per-pair null |r| of the clock against the companions, the
anisotropy and tails of the clock-plane states (fraction of the plane variance on its first axis, kurtosis of |z|), and
the number of genes carrying the projection. A corrected test is a new amendment.
**Limits.** The placebo removes all pairing, so it checks the null only; power is not tested here. Correlated draws
(same rows, same companion sets) make the ten rates dependent; the band is per combination.
**Code.** `exploratory/A276_pairwise_calibration/scripts/calibrate_a276.py`.

## A277 (2026-11-16, autonomous, block A1 of the cell-type program after A276 failed; **NOT BLIND**: every clock and companion statistic of A272 and A274 on the real pairing is known, and A276 in full (null shift of the clock 0.12 to 0.55 SD); NOT seen: any statistic in the corrected form below, any draw with the seeds below): the pairwise test with each gene set read against its own no-pairing null. **The thirteen atlases of A272 and A274; calibration V277.1-2, then the restated tests H277.1-5**

**Why.** A276: under no pairing the clock's |r| statistics sit above those of its companions because their states
differ in shape (isotropy, tails, genes carrying the projection), not because of any pairing. Reading every gene set
against its own placebo null removes that shift and keeps the statistic.
**Corrected statistic.** For each gene set j (clock, 200 matched, 200 responsive per atlas, exactly as scored) and each
of the five statistics (T_SV and D of A272; D, T_SV, T_CX of A274): excess_j = T_j - mu_j, with mu_j the mean of T_j over
placebo draws (donors permuted within compartment x stratum; states unchanged). The test: excess of the clock above the
q95 of the companions' excess, in each family; p = (1 + #{companion excess >= clock excess}) / 201.
**Calibration, before any restated test is read.** Null set N1: 500 draws, seeds [277, 1, d]; calibration set N2: 500
draws, seeds [277, 2, d]. mu_j from N1; in each N2 draw the clock's excess against the companions' excess.
**V277.1:** rejection rate <= 0.075 for each of the ten statistic x family combinations. **V277.2:** mean <= 0.065.
If either fails, the restated tests are not read and the program stops for a report.
**Restated tests (read only if V277.1 and V277.2 pass), mu_j from N1 and N2 pooled (1000 draws):** H277.1 = H272.1
(T_SV, A272 atlases), H277.2 = H272.2 (D), H277.3 = H274.1 (D), H277.4 = H274.2 (T_SV), H277.5 = H274.3 (T_CX), each
against both families with the registered rule of its original (H272.1, H272.2, H274.1, H274.3 above the q95 of both;
H274.2 against the matched family, the responsive one reported). **Descriptive:** the same per atlas and pair class
(A2 reads it), the companion q95 of the excess, and the clock's excess in |r| units.
**Readings.** A restated test that passes stands as before, now on a calibrated test. One that passed and now fails is
reversed, and the text drops it. H272.1 against the responsive family failed before; if it passes now, that is reported
as a change caused by the correction (not blind).
**Limits.** Not blind: the margins of A272/A274 and the size of the null shift were known when this was written; the
correction was chosen because it keeps the statistic and is fixed by A276's diagnosis, not by its effect on the
verdicts. The placebo null removes the pairing between compartments only; a shared non-clock state (ambient RNA,
handling) is controlled by the companions, as before.
**Code.** `exploratory/A276_pairwise_calibration/scripts/corrected_a277.py` (outputs in `output/a277/`).

## A278 (2026-11-16, autonomous, after A277's calibration failed; **NOT BLIND**: the clock and companion statistics of A272 and A274 on the real pairing are known (they set the observed percentiles below), A276 and A277 calibrations in full; NOT seen: any calibrated p-value of the form below): the pairwise test calibrated by construction. **The thirteen atlases of A272 and A274; calibration V278 now; restated tests H278.1-5 read only after the user has seen the A276/A277 report**

**Why.** A276 (location shift) and A277 (correlated responsive companions, wider clock null) show that a q95 rule
against companions is not calibrated for this statistic. The clock's percentile among its companions, read against the
distribution of that same percentile under no pairing, absorbs shift, spread and correlation together.
**Test.** For each statistic T (T_SV, D of A272; D, T_SV, T_CX of A274) and family: pi = share of the 200 companion values
of T below the clock's T (raw statistics, as scored). Null: pi in each of the 1000 placebo draws of A277 (N1 and N2,
seeds [277, 1 | 2, d]; donors permuted within compartment x stratum). p_cal = (1 + #{null pi >= observed pi}) / 1001. A
family is passed when p_cal < 0.05.
**V278 (calibration, null only, run now):** with N1 as the null, the rejection rate over the N2 draws (p_cal < 0.05) is
<= 0.075 for each of the ten combinations. Exchangeable draws make this a code check.
**Restated tests (read after the user has seen the report):** H278.1-5 = H272.1, H272.2, H274.1, H274.2 (matched family
decides, responsive reported), H274.3, each with its original rule, the family criterion replaced by p_cal < 0.05.
Descriptive: the same per atlas and pair class, and the null rate of pi = 1 (clock above all 200 companions), which
bounds the smallest attainable p_cal.
**Readings.** As A277: a test that passed and still passes stands on a calibrated test; one that passed and fails is
reversed; H272.1 against the responsive family is reported as changed by the correction if it changes (not blind).
**Code.** `exploratory/A276_pairwise_calibration/scripts/calibrated_a278.py` (reads `output/a277/null_draws.npz`).

## A279 (2026-11-16, autonomous, block A4 of the cell-type program; calibration of the pairwise offsets in a second mouse tissue plus one ordering question; does not use the companion q95 test held by Q36; **NOT BLIND** for the per-cell-type timed phases of this dataset (A201: aligned arm clock phases by cell type, SMC, fibroblast and endothelium within about 0.6 h of each other; Bmal1 KO amplitude 0.09-0.24 of aligned) and for A273 (colon); NOT seen: any pairwise statistic or timed coupling per cell type of this dataset): are untimed pairwise offsets unbiased in the aorta, and do macrophages couple less to known time than the stromal-vascular cells? **Auerbach 2026 aorta (pools of 2 aortas, ZT 0/6/12/18, both sexes); V279.1, H279.1**

**Data.** `$SCCC_DATA/aorta_jetlag_auerbach2026` (A201 loader): library x cell type pseudobulks of the 14 clock genes,
cell types SMC, fibroblast, EC, macrophage (the A201 MAIN and EXTRA classes; "Other" and T cells left out), >= 30 cells,
compartments in >= 10 libraries of an arm. Arms: aligned (14 libraries, primary), misaligned (8) and Bmal1 KO (10)
descriptive. One stratum per arm (as A273); sex as stratum is a sensitivity.
**Readouts, as A273.** Timed: `Dataset.timed_loadings("zt")` per compartment (phase, coupling to known time); timed
offset of a pair = difference of phases. Untimed: arg of the A272 pairwise statistic over libraries. Error = untimed
minus timed (wrapped), leave-one-library-out jackknife SD.
**V279.1:** median |error| <= 1.0 h over the aligned pairs whose two compartments have timed coupling >= 0.5.
**H279.1:** in the aligned arm, the timed coupling of macrophages is below that of each stromal-vascular class (SMC,
fibroblast, EC), with the jackknife: macrophage minus the smallest stromal-vascular coupling < 0 and its jackknife 90 %
interval excludes 0. Not tested if macrophages are in fewer than 10 libraries (named, not scored).
**Descriptive.** Slope of untimed on timed offsets (only if the timed offsets span >= 1 h; A201 suggests they do not);
|r| of every pair; the KO arm (no clock: low untimed |r| and timed coupling expected) and the misaligned arm.
**Readings.** V279.1 PASS: pairwise offsets are unbiased in a second tissue (with A273, two mouse tissues). FAIL: the
calibration is colon-specific and human offsets get no support beyond it. H279.1 PASS: the myeloid cells of the aorta
follow the tissue clock less than its stromal-vascular cells, as in the human atlases; fail: no order shown here.
**Limits.** Pools of 2 aortas (pool-to-pool variation dominated A201); 1-3 pools per sex x ZT; SMC-dominated tissue;
timed offsets between these classes are small, so V279.1 tests bias near zero, not the slope.
**Code.** `exploratory/A279_aorta_pairwise_truth/scripts/score_a279.py`.

## A280 (2026-11-17, autonomous, blocks A3 and C3 of the cell-type program; staged: hypotheses, atlases and eligibility fixed here, design addendum before scoring; **scoring waits for Q36** (the calibrated form of the companion test, chosen by the user before any statistic of these atlases is computed); seen: A272, A274, A276-A279; NOT seen: any pairwise statistic of the atlases below): second replication of the lymphocyte ranking, and the third and fourth cortex cohorts. **IBDverse ileum and rectum (healthy), PsychAD MSSM (healthy) and RADC (normal); H280.1-3**

**Atlases and eligibility.** Body: `ibd_rectum`, `ibd_ileum` (A17 healthy fits; genome twins exact, data note
2026-11-17). Cortex: `psychad_mssm` (A17 healthy fit, 257 donors) once its genome table passes the twin rule;
`psychad_radc` has no A17 healthy fit (its 26 "normal" donors are reported descriptively only). Rows, strata, families,
pairs (strata with >= 3 shared donors, pooled df >= 15) and companions (200 matched, 200 responsive; seeds [280, atlas,
kind]) as A274. An atlas that fails the twin rule is named and left out.
**H280.1 (second replication of H272.2 / H274.1):** D = T_SV - T_LY over the IBD atlases with >= 3 pairs (SV, X) and
>= 3 pairs (LY, X), above both companion families by the test adopted under Q36.
**H280.2 (cortex, replication of H274.3):** T_CX over the eligible healthy cortex atlases (MSSM), above both families
by the same test.
**H280.3 (offset-sign tie-break):** the inhibitory_CGE minus vascular offset (arg r, h) in the healthy MSSM donors with
its donor-bootstrap 90 % interval (500, within stratum): negative and excluding 0 agrees with PsychAD HBCC (-1.53 h),
positive and excluding 0 with SEA-AD DLPFC (+1.33 h); an interval including 0 decides nothing. Read only if the pair's
|r| is above both companion q95 (A274 rule for offsets).
**Descriptive.** Per atlas and pair class with companions; the IBD family matrix (lymphocytes split T, B, plasma as
labelled: resident-plasma against B is a preview of B2, not a test); RADC normal donors as MSSM.
**Readings.** H280.1 PASS: the ranking holds in a third, independent body source with many more lymphocyte donors
(61 rectum, 53 ileum). Fail with D above the companion medians: not replicated at this power. H280.2 as H274.3. H280.3:
MSSM shares processing with HBCC (same consortium), so agreement with HBCC is weaker evidence than agreement with SEA-AD;
stated with the result.
**Limits.** IBDverse healthy donors are biopsies of a Crohn's study's controls (living tissue). MSSM and HBCC share
protocols; cortex classes carry no lymphocytes.
**Code.** `exploratory/A280_celltype_replication2/scripts/score_a280.py` (the A274 script with these atlases).
**Design addendum (2026-11-17, before scoring; `--stage design`, output `twin.csv`, `design.csv`).** All three atlases
pass the twin rule exactly (0 rows with a count difference): `ibd_rectum` 1,659 rows, 274 donors, 8 compartments;
`ibd_ileum` 1,569 rows, 242 donors, 8 compartments; `psychad_mssm` 1,930 rows, 257 donors, 9 compartments (the genome
table, 33,994 genes, also reproduces the full 7,523-row clock table). Pairs: rectum SV_X 4 (median n 58, min df 47),
LY_X 12 (205, 115); ileum SV_X 4 (52.5, 51), LY_X 12 (191.5, 161): H280.1 is estimable in both. MSSM CX_diff 27
(median n 213, min df 77), CX_same 9; inhibitory_CGE and vascular are both present (H280.3 estimable). RADC (no A17
healthy fit; genome twin exact, 1,173 rows) stays descriptive. Nothing scored; `--stage score` stays blocked until Q36.
**Addendum after Q36 (2026-11-17, before scoring; the user chose option (b), A283).** Readout of H280.1 and H280.2:
the A283 form (clock's pi among its 200 companions, p_cal against the clock's pi in 5000 placebo draws of these
atlases, seeds [280, 9, d], d = 0..4999; the first 1000 already computed, null only); passed when p_cal < 0.05 in both
families. H280.3 gate: the inhibitory_CGE-vascular pair's |r| pi, with p_cal < 0.05 in both families against the same
5000 draws of that pair (replaces "above both companion q95"); the offset interval as registered (500 donor bootstrap
draws within stratum, seed [280, 7]). Read only if A283's V283 passes; if it fails, A280 waits with the line. The A274
descriptive outputs (`--stage score`: family matrix, per-atlas rows, pairs) are run for the tables; their q95 columns
are not used for any verdict.

## A281 (2026-11-17, autonomous, block B1 of the cell-type program; uses no pairwise |r| test (the within-compartment coherence of `common/development.py`, calibrated by per-gene permutation, A237 onwards); **NOT BLIND** for: lymphocyte couplings to the anchored host (A150, 0.069), lymphocyte pairwise |r| with other classes (A272, A274: at the matched median) and per-gene sharing (A275); NOT seen: any within-compartment coherence of these atlases' compartments, any thinned table): do lymphocytes carry a clock of their own that they do not share with the tissue, or no readable clock? **Healthy body atlases of A272, A274 and A280 (IBD); B1; H281.1, H281.2**

**Why.** Lymphocytes share no more clock-gene state with the other cells of a tissue than non-clock genes do. Either
their clock state varies between donors but independently of the tissue (own time: a clock that keeps the time of
where the cell came from), or their clock genes carry no readable clock variation (weak clock, or too shallow).
**Readout.** Per compartment (>= 15 donors): signed clock-gene coherence C across donors
(`decompose.signed_coherence` as in `development.coherence`: log((count + 0.5) / depth), centred within stratum;
per-sample offset fitted out, sign kept), one-sided permutation p against 1000 per-gene permutations within strata. Not used here: any host or partner.
**Depth control (power).** In each atlas with a lymphocyte and a fibroblast-stromal compartment: the fibroblast-stromal
compartment (largest number of donors) is thinned binomially, row by row, to the median depth of the lymphocyte
compartment's clock genes (sum of the 14 clock-gene counts per row; rows of donors present in both), 20 replicates, and
C with its p computed on each (200 permutations per replicate); thinned pass rate = share of replicates with p < 0.05.
The lymphocyte compartment is read on the same donors. Atlases: lung_hlca, heart, organ_breast, organ_colon,
organ_optic_nerve, kong_ileum, kong_colon, ibd_ileum, ibd_rectum (healthy rows; organ_muscle, eye and kidney have no
lymphocyte compartment with >= 15 donors and are named, not scored).
**H281.1 (own time):** across these atlases, the share of lymphocyte compartments (the largest per atlas, on the shared
donors) with C p < 0.05 is >= 0.5. **H281.2 (weak clock, read only with H281.1 failed):** in the atlases where the
thinned fibroblast-stromal pass rate is >= 0.8 (powered), the lymphocyte compartment fails (p >= 0.05) in at least
two thirds.
**Readings.** H281.1 PASS: lymphocytes carry clock-coherent variation between donors that the tissue does not share:
own time (with A272/A274). H281.1 fail and H281.2 PASS: at depths where a stromal clock is read, the lymphocyte clock is
not: a weak or absent clock signal in tissue lymphocytes. Both fail: not decided (power or mixed). Descriptive: C and p
of every compartment by family (all classes), C of 200 expression-matched pseudo sets per lymphocyte compartment
(percentile of the clock), T, B and plasma separately where labelled (IBDverse).
**Limits.** Donors of an atlas are sampled in a limited window of the day (surgery, biopsies; post-mortem for optic
nerve), which limits between-donor clock variation in every class alike; the thinned fibroblasts carry the same window.
Signed coherence reads coherent variation in the waveform pattern; a clock with altered gene phase relations in
lymphocytes would read low.
**Code.** `exploratory/A281_lymphocyte_clock/scripts/score_a281.py`.

## A282 (2026-11-17, autonomous, block A2 of the cell-type program; staged: question, split and table fixed here; **scoring waits for Q36** (the companion readout of each cell uses the form the user adopts); descriptive, no new test; **NOT BLIND**: every clock and companion statistic of A272 and A274 on all donors is known (A272 T_SV 0.499, lymphocytes 0.220; A276-A278 calibrations); NOT seen: any statistic restricted to living or post-mortem donors; seen before registration: the per-donor source classes below (metadata only)): does the cell-type ranking hold in living tissue alone, where no agonal state is shared by the cells of a donor? **The thirteen atlases of A272 and A274, pooled; living against post-mortem donors; descriptive**
**Source classes.** `exploratory/A282_pooled_living/meta/living_postmortem.csv` (per donor, from each atlas's metadata;
1,298 donors over the 15 atlases of A272/A274 builds; collected by a subagent 2026-11-17, counts checked against the
file). Living: surgical resection, biopsy, living donor. Post-mortem: deceased organ donor, autopsy, brain bank.
Calls that rest on recall of the source paper and not on its text (heart: Litvinukova 2020 and Tucker 2020; muscle:
Tabula Sapiens and He 2020) are set to unknown in the primary split; unknown donors are left out of both arms.
Expected arms (healthy donors only, read at scoring): living = organ_breast, organ_colon, kidney_sc, kong_ileum,
kong_colon, organ_muscle (Perez, De Micheli), lung_hlca (biopsy, resection); post-mortem = heart (Koenig, Kuppe),
organ_optic_nerve, organ_eye_anterior, lung_hlca (donor lung, autopsy), the three cortex atlases.
**Table.** Per arm, pooled over atlases with the A272 statistic unchanged (pairs of compartments with >= 15 shared
donors of that arm, strata as scored): family by family median |r| matrix with the clock value, both companion
medians and the companion readout adopted under Q36; per-atlas rows (T_SV, T_LY, D, T_CX where defined) per arm; an
atlas contributes to an arm only if that arm keeps >= 3 (SV, X) and >= 3 (LY, X) pairs.
**Question fixed in advance.** In the living arm: is D = T_SV - T_LY > 0 with T_SV above both companion medians? Yes:
the stromal-vascular-over-lymphocyte ranking is not a product of a shared agonal state. No: stated with the arm's
power (pairs, median df). The post-mortem arm, and the living minus post-mortem difference of D, are descriptive
(atlases differ between arms: organ and donor type are confounded with the split). Anything chosen after seeing the
table is labelled post hoc.
**Limits.** Living tissue is surgical (ischaemia time, anaesthesia, disease-adjacent margins) and biopsies; the arms
share no atlas except HLCA. The cortex is post-mortem only and has no lymphocytes.
**Code.** to be written as `exploratory/A282_pooled_living/scripts/score_a282.py` (A276 loaders, donor filter per arm).

## A283 (2026-11-17, user's answer to Q36: "(b) 5000 draws, then read"; a calibration and the restatement of the cell-type tests; **NOT BLIND**: the clock and companion statistics of A272 and A274 on the real pairing, A276-A278 in full (A278: rates 0.018 to 0.082 with 500 null and 500 evaluation draws; null rate of pi = 1 0.005 to 0.028); NOT seen: any draw with the seeds below, any calibrated p-value on the real pairing): the A278 test with a null of 5000 placebo draws. **The thirteen atlases of A272 and A274; V283, then H283.1-5**
**Test.** As A278: for each statistic T and family, pi = share of the 200 companion values of T below the clock's T;
p_cal = (1 + #{null pi >= observed pi}) / (1 + draws). Null: 5000 new placebo draws (seeds [283, d], d = 0..4999; donors
permuted within compartment x stratum, states unchanged; A277 `draw`). The A277 draws (N1, N2) are not reused.
**V283 (calibration, null only, read before any real pairing):** split the draws into halves (d < 2500, d >= 2500);
with one half as the null, the rejection rate (p_cal < 0.05) over the other half is <= 0.0685 for each of the ten
combinations (five statistics x two families), in both directions. The band is 0.05 + 3 SD of a difference of two
independent proportions at 0.05 over 2500 draws each (A278's band, 0.075 over 500 draws, ignored the error of the null
half). If V283 fails, nothing is read: the cause is reported and the line waits for the user.
**Restated tests (read only if V283 passes; null = all 5000 draws):** H283.1-5 = H272.1 (T_SV), H272.2 (D), H274.1 (D),
H274.2 (T_SV; matched family decides, responsive reported), H274.3 (T_CX), each passed when p_cal < 0.05 in the family or
families of its original rule. Smallest attainable p_cal is the null rate of pi = 1 (reported per combination).
**Descriptive.** Per atlas and pair class: the clock's pi and p_cal against the same 5000 draws of that class median.
**Readings.** As A278: a test that passed and still passes stands on a calibrated test; one that passed and now fails
is reversed and named as such in RP and VERDICT_TABLE. The calibrated form becomes the readout of every later test of
the line (A280, A282, B2, C1-C3).
**Code.** `exploratory/A276_pairwise_calibration/scripts/calibrated_a283.py` (A277 `load_all`, `draw`).
**Design addendum (2026-11-17, after Q36, before scoring; `score_a282.py --design`, output `arms.csv`; no statistic
computed).** Readout: the A283 form (clock's pi among its 200 companions; p_cal against 5000 placebo draws of the arm,
seeds [282, arm, d]); scoring only after V283 passes. Pairs re-derived on each arm's rows with the A272 pair rule
(strata with >= 3 shared donors, pooled df >= 15; stricter than the ">= 15 shared donors" written above). Every row has
a source class (0 missing). Living: HLCA 55 donors (SV_X 3, LY_X 3), breast 126 (9, 20), muscle 27 (4, 0), colon 36
(6, 15), kidney 26 (6, 0), Kong ileum 30 (8, 28), Kong colon 21 (3, 19): D over five atlases (HLCA, breast, colon, Kong
ileum, Kong colon), T_SV over seven. Post-mortem: HLCA 124 (22, 17), heart 27 (Koenig, Kuppe; 12, 4), optic nerve 83
(34, 3), anterior eye 67 (41, 0), cortex HBCC 145, DLPFC 44, MTG 47 (CX_diff 27, 22, 22): D over three atlases (HLCA,
heart, optic nerve), T_SV over four, T_CX over three. HLCA is the only atlas in both arms.

## A285 (2026-11-17, user: "check whether the other companion tests share the bias ... I would register it as a calibration before anything is written up"; a calibration and calibrated restatement, no new biology; **NOT BLIND**: every clock and companion statistic of the tests below is known; seen before registration: a read-only code inventory (`exploratory/A285_companion_calibration/scout/inventory.md`) and the smoke-test placebos of A251/A255 (5 sets, code checks); NOT seen: any placebo draw with the seeds below, any calibrated p): are the companion rules of the bulk and tumour headline results calibrated under no pairing, and do their verdicts stand on the calibrated percentile test of A283? **A251 skin, A255 fat-skin and surgical, A138 tumour-normal, A227 tumour myeloid and T/NK; A263 Part A not run (below)**
**Tests covered.** H251.1, H251.2a, H251.2b; H255.1, H255.3; H138.4 (matched q95 descriptive, responsive H138.4-R);
H227.1, H227.2 (matched, the registered family; responsive from A271 reported). Each exactly as scored: statistic,
pooling weights, 200 + 200 companion sets with the scored seeds, states computed once.
**Placebo.** The pairing unit is permuted within the scored strata (A251 women within compartment; A255 persons within
tissue and cohort; A138 normal samples across patients within cancer; A227 normal-side rows across patients within
cohort x chemistry stratum), the clock and every companion set re-paired identically in a draw; 2000 draws per test,
seeds [285, test, d]. States unchanged (no refit).
**V285 (per test x family):** the placebo rejection rate of the scored rule ("clock above the companion q95") <= 0.0646
(0.05 + 3 SD at 2000 draws). A test that exceeds it carries the A276 bias.
**Restatement H285 (read for every test, whatever V285):** pi = the clock's percentile among its 200 companions; p_cal =
(1 + #{placebo pi >= pi}) / 2001; the clock-specific part of each verdict stands when p_cal < 0.05 in each family its
rule used (A251, A255: both; A138: responsive (H138.4-R) and matched as reported; A227: matched, responsive reported).
The pairing p-values against each test's own permutation null are unaffected (a permutation test of the clock itself)
and are not restated. Reported: rate, pi, p_cal, the placebo rate of pi = 1 (the floor of p_cal).
**Readings.** V285 pass: the scored rule was calibrated for that design (large bulk cohorts, as expected by the user).
V285 fail with H285 pass: the verdict stands on the calibrated test. H285 fail: the clock-specific claim is reversed
and named as such in RP and VERDICT_TABLE before any write-up.
**A263 Part A not run.** Its coupling L is a ratio of medians from a full `Dataset.fit` per gene set: a placebo needs a
refit of 401 sets per draw (baboon 49 min, human 4.4 h, sheep 13.8 h per pass), and permuting donors drives couplings to
0 so L is undefined. Its verdicts are failures by wide margins (clock L above the companion medians, at the 30th to
88th percentile); a calibration bias of the A276 size cannot turn them into passes. Stated as a limit.
**Code.** `exploratory/A285_companion_calibration/scripts/placebo_a285.py` (one module per result, importing the scored
scripts' functions).

## A286 (2026-11-17, user: "Lung fibrosis: the one powered disease contrast, with a power calculation first"; staged: power stage now (permuted labels only, no observed contrast computed), statistic and readings fixed here, the primary chosen by the rule below from the power table, scoring after a power note; seen: anchored fibrosis contrasts on `hlca_fibrotic` (A4/A5/A12 lung T cells; A153: every MDE > 0.2), A272/A274/A282/A283 (healthy HLCA pairwise statistics); NOT seen: any pairwise statistic of a fibrotic donor): does fibrosis change how much clock state the lung's cell classes share? **HLCA, 200 healthy vs 71 fibrotic donors (`hlca_fibrotic_genome`); anchor-free pairwise readout (`clocksync.pairwise`)**
**Data.** `hlca_fibrotic_genome` (rules of `lung_hlca_genome`; its healthy rows identical, 1,799 rows, 0 differences);
case = `lung.FIBROTIC` (IPF 52, SSc-ILD 8, HP, NSIP, end-stage fibrosis, myositis-ILD, ILD), control = Healthy and
Healthy (tumour adjacent). Pairs re-derived on each arm's rows (A272 rule: strata study|assay with >= 3 shared donors,
pooled df >= 15); classes and families as A272 (SV = fibroblast_stromal, endothelial, mural; LY = lymphocyte).
**Candidate statistics (per arm, clock |r| per pair; contrast = case minus control).**
  S1  T_SV: median |r| over SV_X pairs.
  S2  D = T_SV - T_LY.
  S3  F: median |r| of fibroblast-X pairs minus median |r| of endothelial-X pairs (X = non-SV classes).
**Design problem stated in advance.** Cases are 65 of 71 lung explants (taken at surgery, daytime), controls mostly
donor lungs (unknown times): a narrower spread of sampling times lowers every clock |r| of an arm without biology. S1
carries that confound fully; S2 and S3 are within-arm contrasts in which a common loss of shared clock variance cancels
to first order. The A282 within-HLCA comparison (living vs post-mortem SV_X 0.523 vs 0.516) bounds the donor-type effect
on SV sharing but not the time-spread effect.
**Power stage (now).** Null of each contrast: arm labels permuted across donors within stratum (study|assay; strata
without both arms keep their labels), 2000 draws, seeds [286, d]; per statistic: null SD and MDE = 2.49 SD, pairs and
donors per arm. Also with controls restricted to living tissue (surgical resection, biopsy, explant) as a sensitivity
design. No observed contrast is computed.
**Primary rule (fixed now).** Primary = S3 if its MDE <= 0.10, else S2 if its MDE <= 0.10, else no primary (descriptive
only, stated as not powered). S1 is reported descriptively whatever its power (confounded). Specificity: each scored
contrast against the same contrast of the 200 matched and 200 responsive companion sets on the same label
permutations (A283 form: the clock's percentile among its companions, p_cal against its placebo distribution).
**Readings (for the scoring addendum).** A PASS of S3 (fibroblast sharing changes relative to endothelium, beyond both
families): fibrosis remodels the fibroblast's clock coupling to the tissue. FAIL with MDE stated: no change detected at
that size. A difference in S1 alone is not read as biology (sampling-time spread).
**Code.** `exploratory/A286_lung_fibrosis_pairwise/scripts/score_a286.py --stage power`.
**Power-stage amendment (2026-11-17, permuted labels only; no observed contrast computed).** (i) HLCA: only 4 study|assay
strata hold both arms (97 donors); elsewhere cases and controls come from different studies, so arm = study and the
registered "strata without both arms keep their labels" leaves that confound in the null (null mean of S2 0.08). The
only unconfounded HLCA design is within study (donors of the 4 strata with both arms; 48 controls, 49 cases); its MDE
are S1 0.238, S2 0.274, S3 0.387 (2000 draws, `output/power.csv`): above 0.10, so HLCA has no primary. (ii) A second
fibrosis cohort is added to the power stage before any observed value: Natri 2024 ILD (`lung_natri_genome`, exact
twin of `lung_natri`: 1,027 rows, 0 differences; 40 controls, 61 cases, one study, strata = processing site), same
statistics, nulls and rule (within-stratum label permutation; strata with both arms). The primary rule is applied to
whichever cohort reaches it; if both do, both are scored and the HLCA one is the replication.

## A284 (2026-11-17, user: "Test residency (B2) first ... Fix in advance which labels count as resident and which as circulating ... then score once"; block B2 of the cell-type program; **NOT BLIND** for: the pairwise |r| of the coarse classes (A272/A274/A280/A283: lymphocytes at the companion median; plasma and B are both "lymphocyte"), A281 (lymphocyte clocks not coherent), and the label scout (`exploratory/A284_residency/scout/residency_labels.csv`: labels, donors, cell counts only); NOT seen: any statistic of a resident or circulating class separately): does sharing of the tissue's clock state increase with residency? **Healthy atlases; resident against circulating T, plasma against B, macrophages against monocytes; H284.1-3**
**Labels (fixed now; exact strings).**
- *T (H284.1).* Resident: IBDverse `author_cell_type` in {CD8 TRM TRGC2, Tissue resident GD T}; breast (Kumar)
  Census `cell_type` "CD8-positive, alpha-beta memory T cell" (= author CD8-Trm, same 19,213 cells); Kong ileum
  "intraepithelial lymphocyte" (= IELs ID3 ENTPD1, 5,004 cells). Circulating: IBDverse {CD4 CM/naive, CD4 CM, CD8
  naive}; breast "effector memory CD8-positive, alpha-beta T cell", "effector memory CD4-positive, alpha-beta T cell",
  "CD4-positive, alpha-beta T cell" (= CD8-Tem, CD4-Tem, CD4-naive); Kong ileum "naive T cell" (= T cells Naive CD4).
  Atlases: IBDverse ileum, IBDverse rectum (healthy), breast, Kong ileum. No TRM label in HLCA, heart, optic nerve,
  muscle, Pelka colon, Kong colon (scout): not used.
- *B lineage (H284.2).* Resident = plasma, circulating = B, as the built compartments: breast, Pelka colon, Kong ileum,
  Kong colon, IBDverse ileum and rectum.
- *Myeloid (H284.3).* Resident = macrophage, circulating = monocyte, as built: HLCA (lung:macrophage, lung:monocyte),
  breast (macrophage_DC, monocyte). Pelka colon, Kong and IBDverse lack a usable monocyte class (scout): not used.
**New tables (built before any statistic; all genes, same donors and strata rules as the source genome tables).**
IBDverse T subclasses from the local raw file (`author_cell_type`, healthy donors); breast and Kong ileum T classes from
the genome-wide group pseudobulks (`atlas_wg`, Census `cell_type`). Each new table also carries the source's non-immune
classes unchanged (checked against the existing genome tables: identical rows).
**Statistic.** Per atlas and contrast, on donors having both the resident and the circulating class (>= 15 donors):
T_res = median clock |r| (A272 statistic, strata as the atlas) of the resident class with each tissue class X
(non-immune families: fibroblast_stromal, endothelial, mural, epithelial, parenchyma, other non-immune), T_circ the same
with the circulating class on the same X and donors; Delta = T_res - T_circ. **Depth control:** per donor, the deeper
of the two immune rows is thinned binomially (all genes) to the other's clock-gene count; 10 replicates (seeds [284,
atlas, rep]), Delta averaged. Pooled Delta = mean over the contrast's atlases.
**Tests.** Null: the resident and circulating rows of each donor swapped at random (2000 draws, seeds [284, 9, d]; exact
under "no difference between the two classes"). H284.k passes when pooled Delta > 0 with one-sided p < 0.05, Holm over
the three. Specificity (reported, A283 form): Delta of the 200 matched and 200 responsive companion sets on the same
thinned tables; the clock's pi among them, p_cal against the swap distribution of pi.
**Expectation (mouse colon, A21/A31):** macrophage 0.514, plasma 0.495, B 0.167, T 0.000 (couplings to known time).
**What each outcome means for the health chapter (stated now).** All three pass: "cells that live in the tissue share
its clock state; circulating immune cells do not" becomes the explanation of the lymphocyte gap (A272-A283), and the
lymphocyte class is a mixture of the two. T fails with plasma or macrophage passing: residency matters for the B and
myeloid lineages but resident T cells still do not share; the T-cell gap is a property of T cells (with A281: no
readable clock of their own). None passes: the gap is not about residency; the chapter keeps "lymphocytes are apart"
without a mechanism. A Delta that passes the swap test but not the companions is a shared non-clock state of resident
cells, not a clock statement.
**Limits.** Labels are the authors' (marker-defined); circulating T in tissue may be partly recirculating memory;
breast Tem may be tissue effectors (stated by the scout). Kong ileum IELs come from the epithelial fraction.
**Code.** `exploratory/A284_residency/scripts/build_a284.py`, `score_a284.py`.
**Design addendum (2026-11-17, before scoring; `score_a284.py --design`, outputs `design.csv`, `table_checks.csv`; no
statistic computed).** Tables: IBDverse T subclasses built from the raw file (`$SCCC_DATA/atlas_wg/ibdverse_tsub`,
healthy cells; 2,433 ileum and 1,677 rectum source rows unchanged); breast and Kong ileum from the genome builder with
the registered labels first. Breast: only T rows change (CD8_T and part of CD4_T move to the new classes). **Kong
ileum: the source class rule had put "intraepithelial lymphocyte" into `ileum:epithelial`** (IELs are a median 9.6 %,
max 23.2 %, of that class's cells per donor); A284 uses the epithelial rows without them (the only non-immune rows that
differ). Donors with both classes: T: IBDverse ileum 100, rectum 57, breast 94, Kong ileum 22; plasma-B: breast 34,
Pelka colon 26, Kong ileum 25, Kong colon 16, IBDverse ileum 149, rectum 173; macrophage-monocyte: HLCA 107, breast
65. All twelve are eligible (>= 15). Rows: the A17 healthy targets (as A272/A274/A280), new T rows on the same donors
and strata.

## A287 (2026-11-17, user: "The unsigned lymphocyte covariance test, to separate 'no clock' from 'a different clock'"; follow-up of A281 (B1); **NOT BLIND** for: A281 in full (lymphocyte signed coherence C fails in 7 of 9 atlases, IBDverse T C -0.43 / -0.45 at the 12th / 6th percentile of matched sets; thinned fibroblasts pass in 7 atlases), A275 (per-gene sharing); NOT seen: any waveform-free covariance statistic of any compartment): do tissue lymphocytes have structured clock-gene covariance between donors that is not the tissue's waveform (a different clock), or none (no clock)? **The nine atlases of A281 with a lymphocyte and a fibroblast-stromal compartment on >= 15 shared donors; H287.1, H287.2**
**Why a new statistic.** A281's signed C reads the clock genes against the GTEx waveform's phase relations: a clock whose
genes keep other relations reads low. The unsigned R of `development.unsigned` is no remedy: it reads a shift common to
the clock genes of a sample as coherence (A237 simulation: false positives 0.625 to 1.000). The statistic here uses no
waveform and removes the common shift.
**Statistic U (waveform-free).** Per compartment, on the donors shared by the atlas's lymphocyte and fibroblast-stromal
compartments (the A281 pairs): y = log((count + 0.5) / depth) of the 14 clock genes, centred per gene within stratum,
then each row centred over the genes (the per-sample shift removed); U = share of the variance of the gene correlation
matrix of y in its two leading eigenvalues ((lambda1 + lambda2) / trace). Any clock with one phase per gene puts its
covariance in a plane (rank 2), whatever the phase relations; no clock gives the noise level.
**References.** (i) 200 expression-matched gene sets (the A281 sets, `A272.pseudo_sets`, seeds as A281), each treated
identically; pi = share of matched sets with U below the clock's. (ii) 1000 per-gene permutations of y within strata:
p_perm. (iii) 200 responsive sets (`clocksync.pairwise.companions`, seeds [287, atlas, 1]): reported, not decisive (they
covary by design). A compartment **carries clock-gene structure** when pi >= 0.95 and p_perm < 0.05 (no pairing is
involved, so the A276 bias of the clock plane does not arise; under exchangeability with matched genes pi is uniform).
**Power.** The fibroblast-stromal rows thinned binomially (all genes of the sets, as A281 but every gene) to the
lymphocyte compartment's median clock-gene count, 20 replicates (seeds [287, atlas, 100 + rep]), 200 permutations each:
an atlas is **powered** when the thinned fibroblasts carry structure in >= 80 % of replicates.
**Tests.**
- **H287.1 (a different clock):** in the powered atlases, the lymphocytes carry clock-gene structure in >= half.
- **H287.2 (no clock; read only if H287.1 fails):** in the powered atlases, the lymphocytes do not in >= 2/3.
- Neither: undecided. Fewer than 3 powered atlases: not estimable.
**Descriptive.** U, pi, p_perm of every compartment of the nine pairs, unthinned; the leading-plane phase of each clock
gene (angle of its loading on the two leading eigenvectors) and its circular correlation with the waveform acrophases,
for compartments that carry structure; responsive pi.
**What each outcome means (stated now).** H287.1 passes: lymphocytes keep a clock-like covariance of their own whose
gene relations are not the tissue's waveform (A281: signed C low), a candidate "different clock"; the text then says
the lymphocyte clock is unread by the tissue waveform, not absent, and the descriptive phases say how it differs.
H287.2 passes: no clock-gene covariance beyond expression-matched genes where depth-matched fibroblasts show it; "no
readable clock" holds for any phase pattern. Neither: the lymphocyte clock stays undecided.
**Limits.** A non-clock programme acting on several clock genes (PER1 is a glucocorticoid and immediate-early target,
A275) also gives structure: a pass is clock-gene covariance, read with the responsive pi and the phases. The between-
donor clock spread depends on the atlas's sampling window (shared by both compartments).
**Code.** `exploratory/A287_lymphocyte_covariance/scripts/score_a287.py` (A281 tables and sets; `--design` first).

## A288 (2026-11-17, user: "A larger myeloid design ... since the macrophage lead is the only residency-like signal left", then "okay go ahead" to the plan: macrophages against both companion families in every atlas with a macrophage class, IBDverse relabelled; resident against monocyte-derived macrophages as a secondary; macrophage against monocyte dropped (monocytes too rare in tissue, scout); **NOT BLIND** for: A272/A274 family matrices (myeloid with fibroblasts 0.559, matched 0.290, responsive 0.360; myeloid-X 0.409 in A272), A281 (macrophages coherent in most atlases, post hoc), A282 post hoc (fibroblast-myeloid lower in living tissue), A284 H284.3 (macrophage over monocyte, Delta 0.092, FAIL after Holm), the label scout of A284; NOT seen: any statistic of a macrophage class on the calibrated test, any IBDverse macrophage row, any macrophage subclass): do macrophages keep the tissue's shared clock state as stromal-vascular cells do? **Healthy body atlases with a macrophage class; H288.1-2 (primary), H288.3 (secondary)**
**Classes.** MAC = the macrophage class of each atlas as built: lung_hlca `lung:macrophage`; breast, Pelka colon,
muscle, Kong ileum, Kong colon `*:macrophage_DC`; optic nerve (head, nerve) and anterior eye (ciliary, trabecular)
`*_macrophage_DC` (two classes, both MAC); IBDverse ileum and rectum: new rows `<tissue>:macrophage` = the healthy
cells of the author labels {Macrophage CD163hi, Macrophage CD163lo, Macrophage CCL3/4 CXCL9/10, Macrophage MMP9,
Macrophage proliferating} (= Census "tissue-resident macrophage"), from the raw file as A284. Heart (`heart:myeloid`,
mixed) is not used. `macrophage_DC` classes include dendritic cells (limit). Kong ileum: epithelial rows without the
IELs (A284 data note).
**Statistic (A283 form, `clocksync.pairwise`).** X = non-immune classes (A284). Per atlas: T_MAC = median clock |r| over
(MAC, X) pairs; T_LY = median over (lymphocyte, X) pairs; D_MAC = T_MAC - T_LY. Pooled = mean over the atlases where the
statistic has >= 1 pair (an atlas without a pair is dropped from that mean, fixed now; A284 lesson). Companions: 200
matched + 200 responsive sets per atlas (seeds [288, atlas, family]); pi of the pooled clock value among the pooled
companion values; null: 5000 placebo draws (donors permuted within compartment x stratum, seeds [288, d, atlas]);
p_cal as A283.
- **H288.1:** pooled T_MAC is clock-specific: p_cal < 0.05 against both families (the matched level is already seen
  at the family level, so the responsive family decides in practice).
- **H288.2:** pooled D_MAC > 0 beyond both families (p_cal < 0.05 for each): macrophages share more clock state with the
  tissue than lymphocytes do (the H272.2 contrast with macrophages in place of stromal-vascular cells).
- Holm over H288.1-2 (max of the two family p_cal per test).
**Secondary H288.3, resident against monocyte-derived macrophages (A284 design: donors with both, deeper row thinned to
the other's clock-gene count, 10 replicates, Delta = T_res - T_md over X, swap null 2000 draws, seeds [288, 9, d]).**
Decided on HLCA only, where origin is labelled: resident = Census `cell_type` "alveolar macrophage" (HLCA "Alveolar
macrophages", "Alveolar Mph CCL3+", "Alveolar Mph MT-positive", "Alveolar Mph proliferating"), monocyte-derived =
"elicited macrophage" (HLCA "Monocyte-derived Mph"); "lung macrophage" (interstitial perivascular) not used. Passes when
Delta > 0 with one-sided p < 0.05. Descriptive, IBDverse ileum and rectum: CD163hi against CD163lo (a maturity proxy
from markers, not origin). Eligibility: >= 15 donors with both and >= 1 pair per class.
**What each outcome means (stated now).** H288.1 and H288.2 pass: macrophages join the stromal-vascular cells as
keepers of the tissue's shared clock state, and the cell-type section reads "tissue cells and macrophages, not
lymphocytes" (with mouse colon 0.514 and aorta 0.810 to known time). H288.1 passes, H288.2 fails: macrophages share
clock state but not measurably more than lymphocytes. H288.1 fails: macrophages are not shown to keep the shared state
beyond non-clock genes (as stromal-vascular sharing in A272, H283.1); the section names stromal-vascular cells only.
H288.3 passes: long-lived resident macrophages share more than recruited ones (a residency effect within one lineage).
**Limits.** Labels are the authors'; most MAC classes include DCs; the matched-family level is not blind.
**Code.** `exploratory/A288_macrophage_pairwise/scripts/build_a288.py`, `score_a288.py` (`--design` first, a design
addendum with atlases and pair counts before scoring).
**Design addendum (2026-11-17, before scoring; `score_a288.py --design`, outputs `design.csv`, `design_h3.csv`; no
statistic computed).** IBDverse macrophage table built from the raw file (`$SCCC_DATA/atlas_wg/ibdverse_mac`, healthy
cells; donor x label groups: ileum CD163hi 228, CD163lo 231, proliferating 169, MMP9 53, CCL 39; rectum 183, 193, 97,
1, 25), rows >= 20 cells; HLCA split through a class-refinement hook of the genome builder (`genome.build(...,
refine=)`; other classes unchanged). Pairs (MAC_X / LY_X): lung 9 / 20, breast 4 / 28, Pelka colon 4 / 20, muscle 3 / 0,
optic nerve 19 / 5, anterior eye 23 / 0, Kong ileum 4 / 28, Kong colon 3 / 13, IBDverse ileum 4 / 12, rectum 4 / 12. All
ten atlases enter T_MAC; eight enter D_MAC (muscle and anterior eye have no lymphocyte pair). H288.3: HLCA 92 donors
with both (alveolar 124, monocyte-derived 105): eligible if both classes have a pair; IBDverse ileum 60 (descriptive);
rectum 5 (not eligible).

## A289 (2026-11-17, user: "draft the baboon amendment", after the R2.1 discussion of why the brain is less coupled to the body's clock state; **NOT BLIND** for: the baboon and sheep brain and periphery medians of the peripheral-anchor coupling (A42 readout 2: baboon 0.420 vs 0.712, sheep 0.205 vs 0.416), A263 for both species (every gene set drops in the brain), the baboon couplings of AMY 0.60, ARC 0.84 and CER 0.67 (seen in a table head, 2026-11-17), the sheep couplings of Brainstem 0.31, Cerebellum 0.13 and Cerebral cortex 0.30 (same table head), GTEx: pituitary 0.605 between brain 0.406 and periphery 0.679, and GTEx brain regions all above the peripheral amplitude-coupling line (13 of 13, median residual +0.19; post hoc, presentation R2.1); NOT seen: any baboon or sheep coupling of the pineal or pituitary, any other baboon or sheep brain coupling, any baboon or sheep amplitude other than those of the six tissues above (correction before scoring: the same table head showed their own amplitudes, baboon AMY 1.04, ARC 2.11, CER 1.77; sheep Brainstem 0.76, Cerebellum 0.80, Cerebral cortex 0.89), any peripheral amplitude of either species): is the brain's lower coupling explained by the blood-brain barrier, and do brain clocks couple more than their amplitude predicts? **Baboon (Mure 2018, 12 animals) and SheepGTEx, the A42 tables and fits; H289.1-3**
**Data.** `results/figure_2/output/organs/r2_couplings.csv` (A42 readout 2, current code: coupling of every tissue to
the peripheral-anchor host, peripheral tissues read leave-one-out, `own_amplitude` per tissue). No refit. Groups as
A42: brain = the CNS regions (baboon 22, sheep 8), periphery = the rest (35 each); the excluded tissues are pineal
and pituitary (both species) and the ocular tissues (baboon COR, ONH, RET, RPE).
**Barrier (H289.1, H289.2).** Outside the barrier: pineal and pituitary (OUT, 2 per species; neural or neuroendocrine
tissue with fenestrated capillaries). Behind it: the brain regions as A42 (IN; baboon ARC stays in IN, as defined,
although it borders the median eminence and was seen). Statistic: Delta = median coupling(OUT) - median coupling(IN);
p = exact one-sided permutation over the assignment of the OUT labels among OUT + IN tissues (baboon 276, sheep 45
assignments).
- **H289.1 (baboon, primary):** Delta > 0, p < 0.05.
- **H289.2 (sheep, replication):** Delta > 0, p < 0.05 (smallest attainable p 0.022).
**Amplitude (H289.3).** Per species, over the peripheral tissues: least squares coupling = a + b log(own amplitude).
Residual of each brain region from that line; statistic: median residual of the brain regions; p = one-sided sign
test (residual > 0). **H289.3:** in baboon, median residual > 0 and p < 0.05. Sheep: the same, reported (8 regions,
smallest p 0.004). Extrapolation check, stated with the verdict: the number of peripheral tissues within the brain
regions' amplitude range; fewer than 3 labels the verdict "extrapolated".
**Descriptive.** All couplings and amplitudes of OUT, IN, periphery and ocular; Delta with the ocular tissues as a
second IN group (behind the blood-retina barrier, outside the brain); GTEx: pituitary against brain (seen) for
reference; for H289.1 the same Delta on the A263 gene sets is not computed (the barrier reading predicts more sharing
for any gene set, so the companions are no falsifier here).
**What each outcome means (stated now).** H289.1 and H289.2 pass: tissues outside the barrier couple more than
brain regions behind it, in two species: the barrier is a candidate cause of the brain's lower coupling, read with
the limit below. H289.1 fails: the barrier reading is not supported in the timed atlas; the R2.1 notes keep both
readings open. H289.3 passes: brain clocks couple more than their amplitude predicts in a second species: the brain
keeps more of the body's clock state than its weak clocks imply (the R2.1 post hoc note becomes a registered result).
**Limits.** Pineal and pituitary are glands: function and cell composition differ from the brain regions (the pineal
is driven by the SCN through a sympathetic route), so H289.1 cannot separate the barrier from gland type. 12 baboon
animals give noisy per-tissue couplings; sheep strata are mixed. ARC and three regions per species are seen.
**Code.** `exploratory/A289_brain_barrier/scripts/score_a289.py` (reads the A42 table; one scoring run).

## A290 (2026-11-17, user: "run the calibration check first", before H271.2 enters the merged R2.1 slide; a calibration and calibrated restatement in the A285 form, no new biology; **NOT BLIND**: A104 and A271 Part B in full (periphery clock rho 0.734; responsive median 0.109, q95 0.453, p 0.005, i.e. the clock above all 200 responsive sets; matched q95 0.258, p 0.005), A285 (the responsive q95 rule is anti-conservative in bulk, 0.076 to 0.204, where it rests on few pairs; the matched rule is calibrated in every test); NOT seen: any placebo draw of the A104 statistic): is "stronger clocks are more coupled" clock-specific against the responsive family on a calibrated test? **GTEx v10 fast deaths, the 27 peripheral tissues of A104; V290, H290**
**Statistic (A271 Part B, unchanged).** Per gene set (the clock and the 200 responsive sets of A271 Part B, rebuilt
from the gene lists of `a104_sets.csv` with the clock's amplitudes and phases, as `pseudo_waveforms`): Spearman
rho(own amplitude A, coupling c) over the peripheral tissues, every tissue an anchor, leave-one-out host
(`compute_oscillator_gtex.per_tissue`). pi = share of responsive sets with rho below the clock's (sets with undefined
rho count as below).
**Placebo (A285 form).** Donor labels permuted within each tissue (the pairing unit within its stratum), one
permutation per draw applied to every set; seeds [290, d], d = 1..400; draw 0 = the real pairing, must reproduce the
clock rho 0.734 and the responsive median and q95 of A271 Part B (to 3 decimals).
- **V290 (calibration of the q95 rule, descriptive):** share of placebo draws in which the clock rho exceeds the
  responsive q95 of the same draw; calibrated if <= 0.0685 (A283 band).
- **H290 (restated H271.2, responsive family):** p_cal = (1 + #{d: pi_d >= pi_0}) / (1 + 400) < 0.05.
The matched family is not re-run (calibrated in every A285 test; its A271 verdict stands as the A285 reading).
**What each outcome means (stated now).** H290 passes: "stronger clocks are more coupled" is clock-specific against both
families on a calibrated test; the merged R2.1 slide states it. H290 fails: clock-specific against matched genes only;
the slide says so ("not decidable against stress-responsive genes" if the floor of p_cal is above 0.05).
**Limits.** 400 draws (p_cal floor 0.0025; rejection-rate SE about 0.011 at 0.05), fewer than A285's 2000, for cost
(201 fits of about 4 s per draw).
**Code.** `exploratory/A290_a104_calibration/scripts/placebo_a290.py` (`--check`: draw 0 with 5 sets).
**Deviation (2026-11-18, before any placebo draw was read; draw 0 only exists, it reproduces A271 Part B).** Placebo fits
take up to about 5 min (not 4 s): 400 draws need about 11 days on 8 workers. The run goes over one weekend (user) and is
stopped at a time fixed by the calendar, not by the draws; N = the draws completed (expected about 100). H290 and V290 are
scored on those N draws with the same formulas (p_cal floor 1 / (N + 1); with N = 100, 0.0099; V290 SE about 0.022 at
0.05, so "calibrated" is read as <= 0.0685 with that uncertainty stated). If N < 60, H290 is reported as not decidable.

## A291 (2026-11-17, user: "okay, register them", after the discussion of why weaker clocks share less of the body's time (A104, A271 Part B, A289) and whether oscillator theory is contradicted; a calibration of the amplitude-coupling relation at the real GTEx design; **NOT BLIND**: A104 per-tissue table (periphery rho 0.734 over 27 tissues; brain -0.115 over 13), A271 Part B and its post hoc count check, A289; NOT seen: any simulated coupling at these designs): is "stronger clocks are more coupled" what measurement alone produces? **GTEx fast deaths, the A104 rows; coupling of a perfectly locked tissue at each tissue's real design**

**Why.** Coupling is c = S / (S + 2 psi) and own amplitude A = sqrt(S + 2 psi) (S the shared clock-plane variance, psi
the tissue's own variance per dimension; checked on the A104 table). A tissue with a weak clock could read as less
coupled because its fitted coupling is attenuated by counting noise and by the reliability of its host, even if it
were perfectly locked (method limit, roadmap 1.4; A154, A162). A271's post hoc check excluded counting noise as the
cause of the amplitude variable, not attenuation of c. Three scenarios explain the line: (S1) perfect locking, the
relation from measurement only; (S2) locking with own variance that does not scale with the rhythm (weak clocks
noisier relative to their rhythm, biology or non-count technical noise); (S3) transient tracking: low-amplitude clocks
follow recent external cues and depart from the inertial consensus of the strong clocks (oscillator theory away from
steady state). This amendment tests S1; A292 separates S2 from S3.
**Rows and fit.** As A104 (`results/figure_2/scripts/backbone/compute_oscillator_gtex.py`: GTEx fast deaths, tissues
with >= 40 donors, periphery and brain blocks, every tissue an anchor, leave-one-out host, `fit(n_boot=0, gate=False)`).
Reproduction check before any simulation: periphery rho 0.734 and brain -0.115 to 3 decimals.
**Simulation.** Per tissue k, `clocksync.simulate.like(ds, k)` (the tissue's real donors, host means and posterior
covariances, measurement covariances, strata, intercepts): 100 replicates at true coupling 1 (psi = 0, perfect
locking) and 100 at true coupling = the tissue's estimate (design bias, descriptive). Seeds [291, block, tissue index,
1 | 2]. c_ceil,k = median fitted coupling at truth 1.
**H291.1 (primary, periphery):** the relation survives the measurement ceiling: Spearman rho(A_k, c_k / c_ceil,k) > 0,
one-sided permutation p < 0.05 (A shuffled over tissues, 10000 permutations, seed 291). Reported with it: rho(A_k,
c_ceil,k) (how much of the line measurement alone produces) and the median c_ceil.
**Readings, fixed now.** H291.1 PASS: counting noise and host reliability do not produce the relation; S2 or S3 remain
(A292). FAIL with rho(A, c_ceil) at or above 0.5: the amplitude-coupling relation is what the measurement produces;
the organ-level rule ("weaker clocks share less") is dropped, and the brain reading of A289 (brain on the periphery's
line) becomes "brain is measured like a weak peripheral tissue", not a statement about its clocks. FAIL otherwise:
not decided at 27 tissues; the rule stays out of the results.
**Descriptive.** Brain block the same (13 regions; A104 found no relation there); per-tissue c, c_ceil, design bias.
**Limits.** The simulator's noise is counting plus host uncertainty; non-count technical noise (dissection,
composition, RNA quality) sits in psi and cannot be told from biological own variation here: S1 tests measurement
in the narrow sense only. Aggregates only (GTEx controlled access).
**Code.** `exploratory/A291_attenuation/scripts/score_a291.py` (`--check`: reproduction and one tissue).

## A292 (2026-11-17, user: "okay, register them", with A291; **NOT BLIND**: A104 / A271 per-tissue amplitudes and couplings, A259 (tissue departures from body time bounded at 1 to 3 h, mostly carried by matched non-clock genes), A91 (rhythmic genes follow internal time more often than the hour of death), A270 (in GTEx the clock host reads time of death, R 0.631; few companion sets do), A263 and A289 (brain); NOT seen: any regression of a tissue's clock state on the hour of death beyond its host): do weak clocks depart from the body's time toward the external time of day? **GTEx fast deaths, the A104 periphery rows; per tissue, the part of its clock state that follows the hour of death beyond its leave-one-out host**

**Why.** S3 of A291 (transient tracking) predicts that the departures of low-amplitude tissues from the consensus of
the strong clocks are structured: pulled toward recent external time. S2 (own variance that does not scale with the
rhythm) predicts departures unrelated to anything external. The hour of death is the one external time recorded for
every GTEx donor (controlled access, used in aggregate only).
**Statistic.** Per periphery tissue k: clock-plane states x_dk of its samples (`clocksync.paired.states`, GTEx
waveform), its leave-one-out host h_d (the A104 host without k, `Dataset.host(anchors=others)`), and the external
time u_d = exp(i 2 pi tod_d / 24) rotated into the host's frame by delta = arg sum_d h_d conj(u_d) (full host, all
donors). Complex least squares x_dk = mu_k + a_k h_d + b_k u_d (one intercept per stratum). Tracking weight
w_k = Re(b_k / a_k) / sd(h): the external-time component of the tissue's state per unit of its host component, in
the same direction. Leakage note (stated now): errors in h attenuate a_k and move a share of the host component
into b_k; that share is proportional to a_k and therefore the same fraction for every tissue on a common host, so it
shifts all w_k equally and does not create an amplitude dependence.
**H292.1 (primary):** weak clocks track external time more: Spearman rho(A_k, w_k) < 0 over the 27 tissues,
one-sided p < 0.05 against the null of the hour of death permuted over donors (2000 permutations, seeds [292, d];
every w_k recomputed).
**V292 (calibration and power, before the observed statistic is read).** State-level simulation at the real design:
truth h = the full host, tissue states a_k h + N(0, psi_k) + the sample's measurement covariance, the leave-one-out
host perturbed by its posterior covariance, the real hours of death. (i) S2 (no tracking, b = 0): rejection rate of
H292.1 over 500 simulations must be <= 0.0685 (A283 band); if not, H292.1 is not read. (ii) S3: b_k = w0 a_k sd(h)
(A_min / A_k), w0 in {0.1, 0.25}: power. If power at w0 = 0.25 is below 0.8, a failure of H292.1 is uninformative.
**Specificity (reported; decisive only for the wording).** The same rho for 100 expression-matched and 100
stress-responsive pseudo-clock sets (`pseudo_waveforms`, matching as A271 Part B; seeds [292, 1 | 2]); the clock's
percentile pi in each family. A pass of H292.1 with the clock below the 95th percentile of the responsive family is
"weak tissues track a stress state timed with death", not a clock statement.
**Readings, fixed now.** H292.1 PASS (V292 calibrated): low-amplitude clocks depart from the body's time toward
external time; with A291 PASS, the amplitude-coupling relation reads as transient tracking (S3): strong peripheral
clocks hold an inertial internal time, weak clocks follow the present. FAIL at adequate power: no evidence of
transient tracking; with A291 PASS, S2 is the remaining reading (not shown, by elimination). FAIL at low power:
uninformative.
**Descriptive (no verdict).** Mean w over tissues (contains the leakage, not interpretable alone); brain block the same
(rho and the median w of brain regions against the periphery: does the brain follow the hour of death more?);
ventilated donors (Hardy 0) the same.
**Limits.** The hour of death is the end of an agonal period, not the person's recent schedule; fast deaths keep that
period short. Amplitudes A_k are estimated from the same data. Aggregates only; no per-donor value leaves the script.
**Code.** `exploratory/A292_external_tracking/scripts/score_a292.py` (`--stage calibrate | score`).

## A293 (2026-11-17, user: "Yes go ahead", after the synthesis of the clock's organisation into three levels: alignment to the day, organ arrangement, cell participation): which level of the organisation is set by what? **Part 1 GTEx death types with companions (level B); Part 2 conserved stromal-vascular phase order across organs (level C); Part 3 conserved relative gain of cell classes across organs (level C)**

**Seen before registration (NOT BLIND).** Part 1: A109, A109c, A132, A132b in full (clock-gene offset preservation:
Hardy 4 rho 0.784, ventilated post-mortem 0.690 / 0.678, organ donors 0.161; Hardy 4 alignment to the hour of death at
its null), A46. Part 2: the mouse timed phases of A273 (colon) and A279 (aorta), which motivated the hypothesis
(endothelium earliest, then smooth muscle, then fibroblasts in both tissues and both aorta arms); of the human pairs, one
endothelial-fibroblast offset (optic nerve head, A274 list of pairs above both companions) was seen; the counts of
pairs per class were read for this design, no other offset. Part 3: the A17 atlas fits (`a17_healthy_fits.csv`) have
been used throughout Figure 2; the ratios to endothelium defined below have not been computed.

### Part 1 (level B): is the preserved organ arrangement of slow deaths and of ventilated living-brain donors clock-specific?
**Data and statistic.** As A109 / A132 (`exploratory/A109_coupling_nulls/scripts/_helpers.py`; GTEx v10, the 12 A92
body tissues, donors with >= 2 tissues and a time of death). Groups: fast deaths (Hardy 1-2, reference, all donors);
slow deaths (Hardy 4, all); Hardy 3 (all); Hardy 0 ventilated post-mortem cases and Hardy 0 brain-dead organ donors
(COHORT), each at the equal-n draw of A132 (`default_rng(132 + seed)`, seed 0). Per group: Fisher-Lee correlation rho
of the pair phases (pairs with >= 20 shared donors in the group and in fast deaths) with those of fast deaths, for the
clock genes and for the 200 expression-matched pseudo-clock sets per tissue of A109 (`Tissues(200, seed 0)`) and the
200 stress-responsive sets of the A109 note (`compute_h109_responsive.add_sets`), on the same donors and pairs.
Reproduction check first: the clock rho of A109 / A132 (Hardy 4 0.784; ventilated post-mortem 0.690; organ donors
0.161) to 3 decimals, else stop.
**H293.1 (slow deaths):** clock rho above the q95 of the matched family (decisive: the A285 calibration found this
family calibrated, the responsive one anti-conservative), with the tissue-label permutation p < 0.05 (10,000, as
A109). The responsive family is reported (percentile pi of the clock).
**H293.2 (ventilated post-mortem, living brain):** the same.
**Descriptive.** Hardy 3; organ donors (brain death stays descriptive by decision: replication of A132b with companions,
no mechanism reading); the group's companion median rho (does non-clock covariation also carry the fast-death pattern?).
**Readings, fixed now.** H293.1 PASS: in slow deaths, whose internal time is not aligned to the hour of death (A109c),
the organs keep the fast-death arrangement on clock genes beyond non-clock genes: the organ arrangement (level B) does
not depend on alignment to the day (level A). FAIL with companions at or above the clock: the preserved relations are
a shared non-clock tissue structure; the claim is withdrawn for slow deaths. H293.2 likewise for the living-brain
ventilated donors: PASS supports "the arrangement holds without behavioural and feeding cues".

### Part 2 (level C): do stromal-vascular cell classes keep the same phase order in every organ?
**Data.** Human: the pair tables of A272 and A274 (pairs.csv; pairs of compartments of one atlas, untimed offsets
unbiased against known time in mouse colon and aorta, A273, A279). Mouse: the untimed offsets of the same class pairs in
A273 (colon, arm C-DC) and A279 (aorta, aligned arm), from their pairs.csv.
**Statistic.** Class pair endothelial-fibroblast_stromal (primary), endothelial-mural and mural-fibroblast_stromal
(secondary). Each pair's offset oriented as the first class named minus the second. Readable pairs: |r| above the
pair's matched q95 (as stored). Per atlas: circular mean of the readable pairs' offsets, weights |r|^2. Reference sign:
the sign of the mouse untimed offset of the same class pair (colon and aorta must agree; if they disagree the class
pair is not tested).
**H293.3 (primary):** the per-atlas endothelial-minus-fibroblast offset has the mouse sign in more atlases than chance:
one-sided sign test p < 0.05 over the atlases with >= 1 readable pair. Reported: the circular mean over atlases (atlas
bootstrap 90 % interval, 2000 draws, seed 293) against the mouse values.
**Readings.** PASS: a small phase order of stromal-vascular classes, with the mouse sign, is shared across human organs:
a cell-class property of timing exists beside the organ's phase (with endothelium at the interface with blood-borne
signals, the order would fit signals reaching endothelium first; interpretation, not tested). FAIL: no conserved
cell-class phase order is shown; within organs the classes share the organ's phase (A272-A274) and timing is an organ
property only, at this resolution.
**Limits.** Few atlases (6 with a readable endothelial-fibroblast pair by the count read for this design); offsets of
1 h are near the resolution of single pairs.

### Part 3 (level C): does each cell class carry the same gain relative to endothelium in every organ?
**Data.** A17 healthy atlas fits (`results/figure_2/output/atlas_fits/a17_healthy_fits.csv`, rule "learned", status
ok), body atlases only (blood, bone marrow, lymph nodes, spleen, brain, retina and Tabula Sapiens excluded); families
by `clocksync.anchors.family`.
**Statistic.** Per atlas and family f with a clock (fibroblast_stromal, mural, epithelial, parenchyma, myeloid):
g_f = mean log_amp_rel of the family's compartments minus mean log_amp_rel of the atlas's endothelial compartments
(the atlas's host scale and sampling spread cancel); SE from the 95 % bootstrap intervals ((hi - lo) / 3.92, combined
in quadrature; covariance ignored, which overstates the SE). Per family over atlases: random-effects pooled g, tau
(between-atlas SD) and Cochran's Q (`clocksync.stats.meta`).
**H293.4:** the gain relative to endothelium is a cell-class property: for the families with >= 4 atlases, tau <= 0.35
(an amplitude factor of 1.4) in every family. Reported with it: Q p per family, and the GTEx same-class spread (log
amplitude ratios of the tissues dominated by one class, A291 post hoc table) for scale.
**Readings.** PASS: relative gain is conserved across organs; organ differences in gain then come from composition and
an organ-wide drive (both cell classes of an organ scaled together), not from cell classes behaving differently in
different organs. FAIL (tau > 0.35 in a family): that class's responsiveness depends on the organ.
**Limits.** Relative gains, not absolute: an organ-wide drive is invisible to this test by construction and needs the
GTEx decomposition (not in this amendment). Lymphocytes are excluded (no readable clock, A281, A287).

**Code.** `exploratory/A293_organisation_levels/scripts/part1_death_types.py`, `part2_phase_order.py`,
`part3_relative_gain.py`. Aggregates only for GTEx (controlled access).

## A294 (2026-11-17, user: "register the test", after asking whether the breast and heart lymphocyte coherence (A281's 2 of 9 passes) is a clock or a glucocorticoid / stress response; **NOT BLIND** for: A281 in full (breast CD4 T C 0.384, p 0.013, 97th percentile of matched sets; heart lymphocytes C 0.502, p 0.012, 85th percentile; the other seven), A287 in full (breast CD4 T: U at the 59th percentile of matched and below every responsive set; heart: 19th / 60th), A275 (PER1 is the most shared clock gene and a glucocorticoid / immediate-early target), A266 (fetal adrenal coherence vanished without PER1, PER2, NR1D1); NOT seen: any coherence on the 11-gene core loop or against responsive sets for any lymphocyte compartment, any distribution of the matched sets beyond A281's median and percentile): is the lymphocyte coherence of breast and heart carried by the hormone-induced clock genes? **The nine A281 lymphocyte / fibroblast-stromal pairs; breast CD4 T and heart lymphocytes decide; H294.1-2**
**Statistic (A281, unchanged).** Signed clock-gene coherence C across donors (`decompose.signed_coherence`, per-sample
offset fitted out) of the lymphocyte compartment on the donors it shares with the fibroblast-stromal compartment;
log((count + 0.5) / depth) centred within stratum; per-gene permutation within strata, 1000 permutations.
**Core loop.** The 14-gene waveform without PER1, PER2 and NR1D1 (the glucocorticoid-induced genes, as A266), via
`Waveform.subset` (11 genes, acrophases and amplitudes of the GTEx waveform). Seeds [294, k] (k = atlas index of A281).
**Power control (read first).** The same donors' fibroblast-stromal rows thinned to the lymphocytes' clock-gene count
(A281 rule, 20 replicates, 200 permutations each) read on the core loop: the core-loop test of a compartment is
informative only where the thinned fibroblasts pass (p < 0.05) in >= 80 % of replicates.
- **H294.1 (core loop, breast CD4 T and heart lymphocytes, each):** in an informative compartment, the core-loop C has
  p < 0.05.
- **H294.2 (responsive family, each):** the 14-gene C is above the 95th percentile of 200 stress-responsive pseudo-clock
  sets (`A272.pseudo_sets(A, k, 200, "responsive")`, the A281 lymphocyte rows).
**What each outcome means (stated now), per compartment.** H294.1 and H294.2 fail: the coherence is carried by the
hormone-induced genes and matched by a stress-responsive pulse: "hormone-driven, not a lymphocyte clock". Both pass: a
clock-pattern coherence survives without those genes and beyond a hormonal pulse: "a weak lymphocyte clock signal in
that tissue" (an exception to A281/A287, stated as such). Mixed, or the power control fails: undecided.
**Descriptive (same run, no test).** For all nine pairs: the 5th, 50th and 95th percentiles of the 200 matched sets'
C of the lymphocytes (A281's sets and seeds; must reproduce A281's median and percentile), for the slide's random range;
the core-loop C and p and the responsive percentile of every lymphocyte compartment.
**Limits.** Removing three genes lowers the signal of a real clock too (the power control bounds this); two
compartments, so no pooled statement; A281's per-compartment p of these two is not corrected (Holm over nine: 0.11).
**Code.** `exploratory/A294_lymphocyte_hormonal/scripts/score_a294.py` (one scoring run).

## A294b (renumbered 2026-11-18 from a second "A294": the number was already taken by the lymphocyte entry above, registered first; 2026-11-17, user: "Okay you can launch the remaining tests"; the floor test proposed after A291's post hoc decomposition; **NOT BLIND** for the A104 table (S 85-fold, 2 psi 5-fold over 27 tissues) and A259; NOT seen: own variance of any non-clock set at these rows): is the own-variation floor of peripheral tissues specific to the clock genes? **GTEx fast deaths, A104 rows; clock own variance against 100 expression-matched pseudo-clock sets fitted the same way**

**Statistic.** A104 rows and fit (every tissue an anchor, leave-one-out host, `fit(n_boot=0, gate=False)`), periphery
block. For the clock and each of 100 matched sets (`pseudo_waveforms(kind="matched")`, GTEx waveform amplitudes and
phases, matching as A271 Part B, seed [294, 1]): per tissue, own variance 2 psi and shared variance S = z2 V.
Reproduction check: clock periphery rho 0.734.
**H294b.1 (primary):** the clock genes carry own variation beyond the generic floor: median over the 27 tissues of
log(2 psi_clock / median over sets of 2 psi_set) > 0, with a tissue-bootstrap 90 % interval (2000, seed 294) above 0.
**H294b.2 (conditional on H294b.1 PASS):** the clock's excess own variance E = 2 psi_clock - median 2 psi_set scales with
the shared rhythm as person-specific phase or amplitude variation would: Theil-Sen slope of log E on log S_clock
(tissues with E > 0) with a bootstrap 90 % interval including 1 and excluding 0. Interval including 0 and excluding 1:
additive own variation.
**Readings.** H294b.1 FAIL with the interval including 0 or below: the floor is generic (technical or sample-level,
shared with non-clock genes); peripheral tissues are locked to the body's time up to the measurement floor and the
amplitude-coupling relation is the gain over that floor. PASS: the clock carries own variation beyond the floor; H294b.2
says whether it has the signature of individual clock periods (slope 1) or of additive noise (slope 0).
**Limits.** Matched sets have their own shared state (non-clock between-person variation); their psi is the floor
for genes of that expression read through the same projector. Aggregates only.
**Code.** `exploratory/A294_floor/scripts/score_a294.py`.

## A295 (2026-11-17, same instruction; **NOT BLIND** for the GTEx organ map and A291 / A293; NOT seen: any per-sample composition score): within a tissue, does a sample's composition change the gain of its clock state but not its phase? **GTEx fast deaths, A104 periphery rows; marker scores fixed here**

**Scores (fixed now; log((sum of marker counts + 0.5) / depth), standardised within tissue).** Lymphocyte: CD3D,
CD3E, CD2, CD247, MS4A1, CD79A, CD19. Stromal-vascular: DCN, LUM, COL1A2, PECAM1, VWF, CDH5.
**Model.** Per tissue (A292 design: states x, leave-one-out host h, strata intercepts): complex least squares
x = mu + c s + a h + b (s h), s one score at a time. Gain change per SD = Re(b / a); phase change per SD (h) =
Im(b / a) x 24 / (2 pi). SE by donor bootstrap (200, seed [295, tissue]). Pooled over tissues by random effects
(`clocksync.stats.meta`).
**H295.1 (gain, positive control):** pooled Re(b / a) for the lymphocyte score < 0, one-sided p < 0.05 (lymphocytes
carry no readable clock, A281, A287: more lymphocyte RNA dilutes the rhythm).
**H295.2 (phase):** the pooled phase change per SD of each score lies within +-0.5 h (90 % interval inside the bound),
for both scores.
**Readings.** H295.1 and H295.2 PASS: composition changes how loudly a sample expresses the body's time, not when:
phase is not a composition property. H295.2 PASS with H295.1 FAIL: composition varies too little within tissues to
move the gain detectably; the phase bound holds but is a weak statement. H295.2 FAIL: composition moves the phase;
organ offsets may partly be composition.
**Code.** `exploratory/A295_composition/scripts/score_a295.py`. Aggregates only.

## A296 (2026-11-17, same instruction; **NOT BLIND** for the timed organ offsets (adrenal +5.1 h, lung +2.5 h, skin +1.8 / +2.4 h, A123 organ map) and A264 / A266 (PER1, PER2, NR1D1 signal-induced); NOT seen: any offset on a reduced waveform): are the large organ offsets whole-clock properties or carried by signal-induced clock genes? **GTEx fast deaths, timed; baboon replication**

**Statistic.** GTEx fast deaths, A104 periphery rows, timed against the hour of death (`Dataset.timed_loadings("tod")`):
each tissue's offset relative to the circular mean of all tissues, with the full 14-gene waveform and with PER1, PER2
and NR1D1 dropped (`Waveform.drop`). Difference per tissue (reduced minus full), donor bootstrap 90 % interval (200,
seed 296).
**H296.1 (adrenal):** whole-clock property: |reduced - full| <= 1 h for the adrenal gland and its reduced offset >= +3 h.
Secondary, the same rule: lung, skin (sun-exposed, not sun-exposed). Descriptive: per-gene timing on the host
(`Dataset.gene_loadings`) for these tissues; baboon (GSE98965, timed, all tissues an anchor): offsets of ADC, ADM, LUN,
SKI with full and reduced waveforms.
**Readings.** PASS: the adrenal sits about 4 to 5 h from the other organs on its core clock, not only on hormone-
induced genes; an organ-level property to explain (direct neural input is the hypothesis, not tested). FAIL: the
adrenal offset is carried by PER1 / PER2 / NR1D1 induction (waveform distortion), not by the organ's clock phase.
**Code.** `exploratory/A296_genewise_offsets/scripts/score_a296.py`. Aggregates only.

## A297 (2026-11-17, same instruction; **NOT BLIND** for the Manella and Xin timed shifts (A111: liver +11.7 h, WAT +10.6 h, lung +1.9 h under daytime feeding; Xin liver +10.1 h relative to the rest); the WP3 script also wrote a Greenwell per-arm host table (`a112_host_state_by_arm.csv`) that is not in RP and was not read for this design): does removing the feeding rhythm do what the additive cue model, calibrated on inverted feeding, predicts? **Mouse liver: Manella (inversion, calibration), Xin (inversion, second calibration), Greenwell (removal, test)**

**Model.** The liver's timed clock state per arm is a complex amplitude A (state = mu + A e^{iWt}). Additive cues:
A_AL = d + f (direct cues d, feeding cue f, in phase with normal nocturnal feeding); inversion (food in the light phase)
A_inv = d - f. With rho = A_inv / A_AL: f / d = (1 - rho) / (1 + rho), and removal of the feeding rhythm predicts
A_removed / A_AL = (1 + rho) / 2 (closed form: amplitude ratio |1 + rho| / 2, phase shift arg(1 + rho)).
**Estimation.** One projector per dataset and tissue over all its arms (clock genes passing the A111 QC); per arm, OLS of
the state on cos, sin of ZT; A = complex amplitude. rho from Manella Alb-Cre controls (AL vs DF; liver); mouse
bootstrap within arm (2000, seed 297) gives the prediction interval of (1 + rho) / 2. Observed: Greenwell liver,
arrhythmic (adapted) vs ad libitum, A_AR / A_AL with its bootstrap interval; night-restricted vs ad libitum as the
internal control (predicted ratio 1, shift 0).
**H297.1:** the observed removal ratio's amplitude lies inside the 90 % prediction interval of the additive model. Reported
with it: the observed phase shift against the predicted one; rho from Xin (NRF vs DRF) as a second calibration.
**Readings (fixed now).** PASS: an additive driven model explains both inversion and removal; the liver clock's
amplitude is set by its inputs. FAIL with the observed amplitude ratio above the interval (near 1) and a small phase
shift: the liver clock keeps its amplitude when a cue is removed (self-sustained): inputs set its phase, the
oscillator sets its amplitude. FAIL below the interval: removal costs more amplitude than additive cues predict.
**Limits.** Different laboratories, strains and protocols; inversion in Manella is daytime feeding against ad libitum
(not exactly 12 h); Greenwell has the liver only. Not a human result: it decides how removal (human intensive care,
A46) and inversion relate.
**Code.** `exploratory/A297_feeding_model/scripts/score_a297.py`.

## A298 (2026-11-17, same instruction; staged: hypothesis now, reference mapping and organs in a design addendum before scoring; **NOT BLIND** for GTEx gains (A291 post hoc) and A293 Part 3 (relative gains: fibroblast +0.228, mural -0.079, parenchyma -0.126, epithelial -0.220, myeloid -0.399 log units against endothelium)): does composition explain the gain differences between GTEx organs, or is an organ-wide drive needed? **GTEx periphery organs with a matching healthy atlas**

**Statistic.** Per GTEx tissue t with a matching atlas: shares of clock-gene mRNA per cell family, p_c(t), from a
non-negative least-squares fit of the tissue's mean expression (marker genes per family, fixed in the addendum) onto
the atlas's family profiles; predicted relative gain G_t = (sum_c p_c(t) exp(g_c))^2 with g_c from A293 Part 3
(endothelium 0; lymphocyte and other families without a readable clock at exp(g) = 0). Organ drive D_t = S_t / G_t (S
from the A104 table).
**H298.1:** composition explains part of the gain differences: Spearman rho(G_t, S_t) > 0, one-sided permutation p < 0.05.
**H298.2:** an organ-wide drive is needed: the spread of log D_t over tissues exceeds that of log S_t explained by G_t,
read as: the variance of log D_t > 0.5 x the variance of log S_t.
**Readings.** H298.1 PASS, H298.2 FAIL: composition and cell-class gains explain tissue gain; no organ drive needed.
H298.2 PASS: organs drive their cells at different strengths (with A293 Part 3: cell classes keep their relative
gains, the organ scales them together). Both FAIL: not decided.
**Limits.** About 10 organs; marker-based NNLS shares; atlas and GTEx donors differ.
**Code.** `exploratory/A298_organ_drive/scripts/` (design addendum first).

### A298 design addendum (2026-11-17, before any statistic; composition only, `scripts/build_shares.py`): not scored
The composition step does not meet the quality the registered statistic needs. Marker-based NNLS of the GTEx tissue
profiles on the atlas family profiles explains little of them (weighted R^2 0.008 to 0.403 over the 11 tissues;
unweighted R^2 negative in 9), and several atlases lack families the tissue contains (subcutaneous adipose: no
endothelial, fibroblast or mural family; skeletal muscle: no parenchyma; transverse colon: fibroblast and myeloid
dropped for < 2,000 cells), giving implausible shares (skin 42 % myeloid mRNA). A298 is closed before scoring; S and
the gain predictions were not computed. A test of the organ-wide drive needs a deconvolution with complete,
bulk-calibrated references (or single-cell data of several organs of the same donors), not available here.
<!-- src: paper/exploratory/A298_organ_drive/output/{shares,design}.csv -->

## A299 (2026-11-17, user: "Yes start these suggestions"; **NOT BLIND** for the GTEx adrenal offset (+4.42 h, A296), its own variation and coupling (4.14, 0.46; A294b / A104), the baboon cortex +0.91 h and medulla +5.90 h (A296), and A295 (composition does not move phase in other tissues); NOT seen: any medulla or cortex score in GTEx): does the medulla content of a GTEx adrenal sample move its clock phase towards the medulla's? **The one tissue where composition should move phase: two compartments at different phases**

**Why.** If the adrenal cortex follows the body's time and the medulla (a sympathetic ganglion on the direct neural
route) sits about 5 h earlier, as in the baboon, a bulk sample's phase is a mixture: samples with more medulla RNA should
be more advanced. A295 found no phase effect of composition in tissues whose cell types share one phase; here the
prediction is positive.
**Scores (fixed now; log((sum of marker counts + 0.5) / depth), standardised within the adrenal).** Medulla: TH, DBH,
PNMT, CHGA, CHGB. Cortex (specificity): STAR, CYP11A1, CYP11B1, CYP17A1, CYP21A2.
**Model.** A295 model on the adrenal only: x = mu + c s + a h + b (s h) with the adrenal's leave-one-out host (every
peripheral tissue an anchor, A104 rows); phase change per SD = Im(b / a) x 24 / (2 pi) h (positive = advanced, the
convention of the offsets); gain change Re(b / a); donor bootstrap (2000, seed 299).
**Donors.** Primary: Hardy 1 to 4 (fast and slow deaths, whose organ arrangement is preserved, A293); ventilated
donors excluded (organ arrangement pushed towards 0 / 12 h, A109). Descriptive: fast deaths only (A104 rows).
**H299.1:** medulla score, phase change per SD > 0 (towards the medulla's advance), one-sided bootstrap p < 0.05.
Reported with it: the cortex score's phase change (predicted opposite sign or null) and both gain changes.
**Readings.** PASS: composition moves phase where compartments differ in phase; the GTEx adrenal offset is at least
partly a medulla contribution, as in the baboon (with A295: phase is an organ-compartment property, and a bulk tissue
mixing two compartments reads in between). FAIL: the human adrenal phase is not explained by its medulla content at
this power; the adrenal offset stays an organ-level observation.
**Limits.** About 100 adrenal samples; markers fixed from adrenal biology, not checked for cross-expression; the
medulla's own phase in humans is not measured.
**Code.** `exploratory/A299_adrenal_medulla/scripts/score_a299.py`. Aggregates only.

## A300 (2026-11-17, same instruction; **NOT BLIND** for the post hoc look at "frontier" tissues (lower coupling 0.51 against 0.69, higher own variation 2.80 against 1.87, more hour-of-death weight in A292, offsets +0.24 against -0.22 h; list chosen after seeing the table); NOT seen: any season contrast): do the tissues facing the external environment change their clock relation to the body with the season? **GTEx fast deaths, A104 periphery rows; summer against winter deaths**

**Groups (fixed now, by anatomy).** Frontier, tissues whose epithelium faces air, light or the gut lumen: Lung, Skin -
Sun Exposed, Skin - Not Sun Exposed, Esophagus - Mucosa, Colon - Transverse, Stomach. Internal: the other periphery
tissues except Whole Blood (circulating cells) and Testis (germ cells). Tissues with >= 20 donors in each season arm.
**Statistic.** Death season (`DTHSEASON`, protected): Summer and Winter donors (Spring and Fall not used). Per arm, the
A104 periphery fit (every tissue an anchor, leave-one-out host); per tissue the offset (h, relative to the circular
mean) and log z2 (relative amplitude, comparable between arms). Delta_k = offset_summer - offset_winter (wrapped).
T = median |Delta| over frontier tissues minus median |Delta| over internal tissues.
**H300.1:** T > 0, one-sided p < 0.05 against season labels permuted over the summer and winter donors (500
permutations, seeds [300, d], both arms refitted).
**Descriptive.** The same for log z2 changes; skin and lung separately; the direction of the skin and lung offsets
(summer minus winter).
**Readings.** PASS: frontier tissues move their clock relative to the body with the season more than internal tissues:
an environmental input beyond the body's time reaches them (with the post hoc pattern: lower coupling, more own
variation). FAIL: no season-dependent frontier effect at this power; the frontier pattern stays a post hoc observation.
**Limits.** About 60 donors per season arm; season also changes day length, temperature, infections and activity; GTEx
death season is not the person's exposure history; offsets per arm carry the fit's noise (the permutation null includes
it).
**Code.** `exploratory/A300_frontier_season/scripts/score_a300.py`. Aggregates only.

## A301 (2026-11-17, user: "Okay go ahead", after the A299 post hoc: the steroidogenic cortex score moves adrenal phase -1.15 h per SD and survives adjustment for medulla, fat and capsule content; **NOT BLIND** for A299 and its post hoc in full; NOT seen: any ACTH-response or cortex-identity score): is the adrenal phase effect carried by ACTH activity or by cortex content? **GTEx adrenal, the A299 rows and model**

**Scores (fixed now; log((sum of counts + 0.5) / depth), standardised within the adrenal).**
- ACTH / cAMP response, without steroidogenic enzymes: NR4A1, NR4A2, NR4A3, SIK1, CREM, MC2R, MRAP.
- Cortex identity, genes not known as acute ACTH targets: NR5A1, NR0B1, KCNK3, VSNL1, CCN3.
- Steroidogenic (the A299 cortex score, for reference): STAR, CYP11A1, CYP11B1, CYP17A1, CYP21A2.
**Model.** The A299 model with the ACTH and identity scores fitted jointly (x = mu + sum_j c_j s_j + a h + sum_j b_j
s_j h); phase change per SD Im(b_j / a) x 24 / (2 pi) h, gain change Re(b_j / a); donor bootstrap 2000 (seed 301);
Hardy 1-4 primary, fast deaths descriptive.
**H301.1:** ACTH activity moves the phase towards the body's time: ACTH-score phase change < 0, one-sided bootstrap
p < 0.05, in the joint model.
**H301.2:** cortex content does not: identity-score phase change within +-0.5 h (90 % interval inside the bound), in
the joint model.
**Readings.** H301.1 PASS, H301.2 PASS: the adrenal effect is hormonal activity, not composition: samples under stronger
ACTH drive sit closer to the body's time (a hormonal input to the adrenal clock, with the fetal cortex result as
context). H301.1 FAIL, identity score moving phase: the effect is cortex content (composition of compartments at
different phases, the non-cortex part being advanced). Both moving: not separable. Neither: the steroidogenic effect of
A299 is not reproduced by either component. Reported: all three scores' correlations; the three-score joint model.
**Limits.** 75 samples; NR4A and CREM are also general stress / immediate-early genes (agonal state); marker sets
fixed from adrenal biology, consistency checked descriptively, not used to select genes.
**Code.** `exploratory/A299_adrenal_medulla/scripts/score_a301.py`. Aggregates only.

## A302 (2026-11-17, user: "Yes go ahead", after the reading that person-specific cue timing (e.g. meal times) would give a departure of constant size in the clock plane (A294b slope near 0); **NOT BLIND** for A291, A292, A294b (per-tissue S, own variation, coupling) and the A142 organ blocks; NOT seen: any correlation of departures between tissues; the trait-or-state test proposed with it is not registered: every repeat-biopsy cohort samples one tissue per person, so a tissue's departure cannot be separated from the person's own offset): do the departures of meal-driven tissues share a direction within a person beyond the body's time? **GTEx fast deaths, A104 periphery rows; the liver's departure against gut and adipose tissues versus other tissues**

**Groups (fixed now).** Anchor: Liver (the tissue that follows feeding in both mouse inversion designs: Manella +11.7 h,
Xin +11.3 h, A111). Meal-driven group G: Stomach, Colon - Transverse, Colon - Sigmoid, Esophagus - Mucosa, Adipose -
Subcutaneous, Adipose - Visceral (Omentum) (digestive and adipose tissues; the mouse data on disk are mixed for adipose
and heart / kidney, so the group is fixed from digestive physiology, not from a mouse threshold). Other group N: the
remaining periphery tissues except Liver, Whole Blood and Testis.
**Departures.** For the liver and a tissue j: a host built from all periphery tissues except both (`clocksync._core.host`,
the A104 rows), each tissue's state regressed on it with strata intercepts (complex least squares), residuals e. Over the
donors with both tissues: cross-moment M = mean e_liver conj(e_j), minus the part induced by the common host error,
L = a_liver conj(a_j) x mean trace of the host's posterior covariance; same correction for each tissue's own variance;
r_j = |M - L| / sqrt((M_ll - L_ll)(M_jj - L_jj)) (NaN when a corrected variance is <= 0). Pairs with >= 30 shared donors.
**Statistic.** T = median r_j over j in G minus median r_j over j in N.
**H302.1:** T > 0 with (i) p < 0.05 against random groups of the same size drawn from G and N together (10,000 draws,
seed 302), and (ii) T above the q95 of the same statistic for 100 expression-matched pseudo-clock sets (A271 Part B
matching; seed [302, 1]) on the same donors (specificity: shared handling of gut samples would also move matched genes).
**Readings.** PASS: meal-driven tissues depart from the body's time together within a person: part of the personal
departure is a shared, cue-driven displacement (person-specific feeding history), not independent noise. FAIL (i):
no group structure of departures at this power: the departure looks tissue-independent. PASS (i) but FAIL (ii): the
shared departure is not clock-specific (shared sampling of gut tissues). Reported: r_j per tissue, the corrected and
uncorrected values, the same for the matched median.
**Limits.** About 50 to 180 shared donors per pair; the leakage correction assumes isotropic host error; fasting
state before death is unknown; adipose sensitivity to feeding is mixed in the mouse data.
**Code.** `exploratory/A302_departure_groups/scripts/score_a302.py`. Aggregates only.

## A303 (2026-11-17, user: "for the adrenal can we see the same structure we found in mouse or sheep?"; descriptive with a fixed criterion; **NOT BLIND** for the human (+4.4 h) and baboon (cortex +0.9 h, medulla +5.9 h) adrenal offsets; NOT seen: any adrenal offset in mouse or sheep): is the adrenal offset from the other organs in mouse and sheep?
**Data.** Mouse: Wolff 2023 (GSE201207; adrenal, heart, kidney, lung, muscle, hypothalamus; 12 circadian times in
constant darkness x 3 ages; pools of 2 mice), timed (`timed_loadings("ct")`), per age and pooled with age as stratum.
Sheep: the sheep atlas (untimed; 21 animals with adrenal), `Dataset.fit` with every non-brain tissue an anchor.
**Readout.** The adrenal's offset (h) relative to the circular mean of the non-brain tissues (mouse: heart, kidney,
lung, muscle; sheep: all non-brain tissues), with a bootstrap 90 % interval (mouse: libraries within CT, 1000; sheep:
donors, 500).
**Criterion (fixed now).** "The same structure as in humans" = the adrenal's offset >= 2 h from the mean with its
interval excluding 0 (sign reported against the human +4.4 h, positive = advanced). A whole-gland offset near 0 in
both species, with the human result, would fit the adrenal offset being specific to GTEx (death) or to primates.
**Code.** `exploratory/A296_genewise_offsets/scripts/score_a303.py`.

## A304 (2026-11-18, lead prompt `prompts/organisation_followups_2026-11-18.md`, task A1-A3; staged: this entry fixes design, statistics, readings and the power rule; a design addendum with the power results and the contrasts to be scored follows before any observed contrast is computed; **NOT BLIND** for the A138 per-cancer tumour-normal correlations (all pairs), the A127-A224 tumour-only coherence results and the H138.8 purity split; **seen before writing:** design counts only (pairs per cancer; per arm: TP53 altered / intact, ER- / ER+ and LumA / other in BRCA; no clock statistic of any arm); NOT seen: any tumour-normal statistic computed within an arm): is the loss of tumour clock coherence with tumour state a loss of gain or a rise of departure? **TCGA paired tumour-normal (`tcga_paired`, 10 cancers, 517 pairs)**

**Framework.** Normal sample N = h + e_N (h the patient's clock state at surgery, variance S; own departure 2 psi_N);
tumour T = g h + e_T (gain g, departure 2 psi_T). With one normal per patient the host is N itself, of reliability
rel = S / (S + 2 psi_N). Per arm, on the states of the cancer's one projector (`q3_lib.tcga_states`, as A138), centred
within arm and sample type, Poisson-corrected: V_N, V_T (clock-plane variances), C = mean z_T conj(z_N) (n - 1).
- **Gain** G = |C| / V_N = |g| rel (the tumour's regression slope on its own normal; psi_T does not enter it).
- **Departure** R = (V_T - |C|^2 / V_N) / V_N^all = (g^2 S (1 - rel) + 2 psi_T) / V_N^all, with V_N^all the cancer's
  normal variance over all its pairs (fixed across arms; makes R unitless and comparable between cancers).
- **Offset** arg(C_case / C_ref) (h), descriptive.
Within a cancer rel cancels from the gain ratio if the two arms' normals are alike; the pure gain loss lowers G and
also R (through g^2 S (1 - rel)); a pure departure rise raises R and leaves G. The host check below tests "normals alike".

**Contrasts (case vs reference; cells = cancer; arms with >= 8 pairs each, else the cancer leaves that contrast).**
- H304.1 (primary): proliferation high vs low: median split, within cancer, of the A127 score (mean z of Hallmark G2M
  without clock genes, z within the project over all its primaries; `figure_5/.../_shared.load_tumours`, column
  `prolif`).
- H304.2: ER negative vs positive (BRCA; clinical XML field, `load_tumours` column `er`).
- H304.3: TP53 altered vs intact (A215 MC3 status, `A215_checkpoint_coherence/output/status_tcga.csv`).
- H304.4: luminal A vs other PAM50 subtypes (BRCA; PanCanAtlas subtypes, `A217_subtype_matched/output/subtypes_tcga.csv`).
**Pooled statistics** (weights n_c = pairs of the cancer in the contrast): Q_G = log(sum n_c G_c,case / sum n_c
G_c,ref) and Q_R = log(sum n_c R_c,case / sum n_c R_c,ref). Null: labels permuted within cancer (2000 draws, seed
[304, contrast index]). One-sided p: gain loss P(null <= Q_G); departure rise P(null >= Q_R) (add-one).
**Specificity (same labels):** 200 expression-matched pseudo-clock sets per cancer (`paired.pseudo_waveforms`,
genome-wide counts of the pairs as A138 responsive note; seed [304, 7, cancer index]), pooled per set index; a readout
is clock-specific if Q_G is below the matched q05 (gain) or Q_R above the matched q95 (departure) (A285: calibrated in
bulk). 200 responsive sets reported, not decisive (A285: anti-conservative).
**Host check (same permutations):** log ratio of V_N case / ref pooled as Q_G; reported with its p. If |log ratio|
has p < 0.05 the arms' normals differ and the gain readout is read with that caveat (rel may differ).
**Tests.** Per contrast two readouts: gain loss (Q_G < 0) and departure rise (Q_R > 0), each "detected" if
permutation p < 0.05 after Holm over the readouts of its family (primary: H304.1's two; secondary: the readouts of the
powered contrasts among H304.2-4) and clock-specific.
**Readings.** Gain loss detected, departure rise not: malignant cells dilute or weaken the body's rhythm in the tumour
(as lymphocytes in normal tissue), they do not drift from it. Departure rise detected, gain loss not: the tumour clock
drifts from the body's time. Both: a mixture. Neither, powered: no change beyond the MDE in either; underpowered:
uninformative. Not clock-specific: the change is shared by matched genes (tumour-state shared variance, not the clock).
**Power stage (A2, run and reported before the addendum).** Simulation at the real design of each contrast (real arm
sizes and per-sample Poisson variances; cancer moments S, g, psi_T calibrated from all pairs of the cancer, the A138
data already seen, with rel set to the GTEx coupling of the matching tissue from the A104 table: BRCA breast, COAD
colon transverse, KIRC / KIRP kidney cortex, LIHC liver, LUAD lung, PRAD prostate, THCA thyroid; BLCA and UCEC the median
of these). `clocksync.simulate.like` needs an anchored host of >= 2 compartments and cannot represent this design; the
simulation is written for the paired model (complex Gaussian states, Poisson noise of each sample). Scenarios for the
case arm: null; pure gain loss g x k_G; pure departure rise 2 psi_T x k_D; 500 replicates, 500 permutations each.
**Power rule (fixed now):** a contrast is scored as a test only if power >= 0.8 at k_G = 0.5 (the tumour expresses half
the body's rhythm, e.g. malignant cells without a clock forming half of its clock RNA) for the gain readout, or at
k_D = 2 (departure doubled) for the departure readout; readouts below the rule are reported descriptively. Reported:
power curves, MDE (k at 80 %), the type I error of each readout, and the cross-readouts (Q_R under a gain loss, Q_G
under a departure rise).
**Sensitivity (reported, H304.1):** cells = cancer x purity half (A127 purity proxy, `load_tumours` column `purity`).
**A3 (descriptive, after A1).** Per cancer, the tumour's coupling to the body's time c_T = corr^2 / rel (corr the A138
Poisson-corrected tumour-normal correlation, rel as above; BLCA and UCEC without a matching tissue are left out) and the
tumour gain on body time |g| = G / rel, with patient-bootstrap 90 % intervals (1000, seed 304) on corr and G; ranked.
Caveats stated with it: rel is from GTEx deaths on a multi-tissue host, the normal is adjacent tissue at surgery.
**Code.** `exploratory/A304_cancer_gain_departure/scripts/{power_a304.py, score_a304.py}`.

### A304 design addendum (2026-11-18, after the power stage; no observed contrast computed)
`power_a304.py --jobs 4` (38 s; `output/power/{calibration,design,power,mde}.csv`). The power rule was applied at
alpha 0.025 (Holm's first step over the two readouts of a family). Type I error at 0.05 over 500 null replicates: 0.044
to 0.068 for every readout and contrast (Monte Carlo SE about 0.01). Cross-readouts: a pure gain loss raises Q_R in
<= 0.03 of replicates, a pure departure rise lowers Q_G in <= 0.07: the two readouts separate. Signatures at the
rule points (H304.1): gain x 0.5: mean Q_G -0.509, Q_R -0.100; departure x 2: Q_G +0.087, Q_R +0.625.
| contrast | cells (pairs, case / ref) | power gain at k_G 0.5 | power departure at k_D 2 | scored |
|---|---|---|---|---|
| H304.1 proliferation | 10 cancers (517, 257 / 260) | 0.840 | 1.000 | gain and departure (primary, Holm over 2) |
| H304.2 ER (BRCA) | 1 (99, 21 / 78) | 0.120 | 0.698 | descriptive only |
| H304.3 TP53 | 5 (BRCA, COAD, LIHC, LUAD, PRAD; 299, 100 / 199) | 0.524 | 0.998 | departure (secondary) |
| H304.4 LumA (BRCA) | 1 (105, 58 / 47) | 0.200 | 0.854 | departure (secondary) |
Secondary family: H304.3 and H304.4 departure readouts, Holm over 2. Every other readout is reported descriptively.

## A305 (2026-11-18, lead prompt task A4; staged: power stage only now (placebo draws, no observed statistic), scored only if the power rule below passes, after a dated power note; **NOT BLIND** for A149 H149.3 (malignant cells covary with the tumour stromal host, specificity not tested), A138 H138.1, A227 / A271 (tumour fibroblasts and endothelium keep the patient's normal state, clock-specifically) and the healthy pairwise results A272-A283; NOT seen: any pairwise |r| of malignant cells with stromal-vascular compartments of the same tumour sample): do malignant cells share the clock state of their own tumour's stromal-vascular cells beyond the companion families? **Single-cell tumour cohorts (A149 / A227 tables with genome-wide counts), anchor-free pairwise readout (`clocksync.pairwise`, A283 form)**
**Pairs.** Within each cohort (stratum = cohort x chemistry, as A149), compartment pairs (malignant, X) with X in the
tumour sample's stromal-vascular compartments (fibroblast, endothelial, mural); the malignant compartment as A138
(`q3_lib.study_table`). Shared-donor rule of `pairwise.pair_index` (min_stratum 3, min_df 15). Companion pairs
reported: (tumour SV, tumour SV) (positive control, expected high) and (malignant, T/NK) where present.
**Statistic.** T_MS = median over cohorts of the cohort median clock |r| (uncorrected) of (malignant, SV) pairs; read as
the clock's percentile pi among 200 matched and 200 responsive companion sets, p_cal against 5000 placebo draws
(donor labels permuted within compartment x stratum), as A283.
**H305.1:** T_MS above both companion families (p_cal < 0.05 against each; the matched family decides if the
responsive one cannot reach 0.05, A285).
**Power stage (now; placebo only).** Per cohort: patients, pairs, df; the placebo rate of pi = 1 per family (smallest
attainable p_cal); the placebo SD of T_MS minus its matched-companion median and MDE = q95 + 0.84 SD - placebo mean
of that difference (one-sided 0.05, power 0.8).
**Power rule (fixed now).** H305.1 is scored only if (i) the smallest attainable p_cal is < 0.05 for the matched family
and (ii) the MDE of T_MS minus the matched median is <= 0.10 (A286 bound). Otherwise A305 closes at the power stage
("not powered"; no reading of malignant clock sharing).
**Readings (for scoring).** PASS: malignant cells share their tumour's stromal clock state clock-specifically (the
tumour's own clock environment reaches them). FAIL with MDE: not detected at that size; sharing matched by companions:
the A149 covariation is a shared non-clock state.
**Code.** `exploratory/A305_malignant_stroma_pairwise/scripts/power_a305.py`.

## A306 (2026-11-18, lead prompt task B; staged: power stage (permutation nulls only) then a dated power note, scoring only of powered readouts; **NOT BLIND** for: A164 / A202 / A233 / A235 (module d* equivalent to 0 within about +-0.05 for T2D, hypertension and smoking; offsets equal within about 1 h), A235's descriptive module log z2 ratio (T2D -0.003 fast, -0.047 vent); per-tissue z2 of each arm was written to the A233 / A235 output tables (not read for this registration); **seen before writing:** long-term steroid counts (MHSTRDLT, aggregates: 6 fast, 8 vent cases); NOT seen: any psi of a condition arm, any module-minus-companion z2 or psi contrast): do chronic conditions change the gain or the departure of the organ they affect? **GTEx, A233 matched procedure; readouts log z2 ratio and log psi ratio**
**Conditions and modules (fixed now).** T2D (A164 label; module = A164 MET: liver, both adiposes, pancreas, skeletal
muscle), hypertension (MHHTN; module = arteries (aorta, coronary, tibial) and heart (left ventricle, atrial
appendage)), smoking (MHSMKSTS; module = lung). Companion: A164 COMP (tibial nerve, both skins). Long-term steroids
are not used (6 / 8 cases, below the A233 minimum of 10 matched cases per tissue).
**Matching and estimation.** A233 exactly (`A231_age_coupling/scripts/matched_label.py`: up to 3 controls, exact sex and
Hardy class, age caliper 5 y, then ischaemic time; A164 anchors; `Dataset.contrast` per tissue, pooled host prior;
1000 permutations swapping the case label within matched sets, shared across tissues; vent and fast cohorts). Only
module and companion tissues are fitted.
**Readouts.** Per tissue the contrast table's log z2 ratio (gain: case / control loading on the common host) and log psi
ratio (departure). Per condition and readout: Delta = (module pooled) - (companion pooled), inverse-null-variance
pooling over tissues within cohort, then over cohorts (A233 module-minus-companion). Subtracting the companion
removes what the arms' host does to every tissue alike (host gauge and spread, A258), so Delta is the organ's change
relative to the rest of the body. Two-sided permutation p.
**Power rule (fixed now).** MDE_80 = 2.49 x null SD of Delta (A233). A readout is scored if MDE_80 <= log 2 = 0.693
(a doubling or halving of the organ's z2 or psi relative to the companions); otherwise descriptive. Holm over all scored
readouts (at most 6).
**Tests.** H306.c.g (gain) and H306.c.d (departure), c in {t2d, htn, smoking}: Delta != 0, Holm p < 0.05.
**Readings.** d* is known to be unchanged (A164, A233), so z2 and psi can only move together (coupling = z2 V / (z2 V + 2
psi) fixed): both readouts detected in the same direction = the organ's clock expresses more or less of everything
(gain and departure scaled together, a change of amplitude, not of synchrony). Gain alone or departure alone detected
would contradict the d* equivalence and is reported as such (power of the two readouts differs). Neither, powered: no
change of the organ's gain or departure beyond the MDE. Not powered: uninformative.
**Reported.** Per tissue and cohort log z2 and log psi ratios; module and companion pooled values separately; arm-prior
sensitivity (`host_prior="arm"`, observed only); matching balance.
**Code.** `exploratory/A306_disease_gain_departure/scripts/compute_a306.py --stage power | score`.

### A306 power note (2026-11-18, permutation nulls only; no observed statistic computed)
`compute_a306.py --stage power --jobs 4` (2218 s, peak 2.9 GB; `output/{power,power_tissues,matching,balance_summary}.csv`).
Eligible tissues: T2D vent 8 (module 5 incl. pancreas) / fast 7 (pancreas out); hypertension 8 / 8; smoking 4 / 4.
MDE_80 of Delta (module minus companion; bound log 2 = 0.693): T2D gain 0.327, departure 0.288; hypertension gain 0.346,
departure 0.264; smoking gain 0.700 (not powered: descriptive), departure 0.418. Scored: five readouts, Holm over 5.
Balance after matching (median |SMD|): BMI 0.70 (T2D fast), 0.32 (T2D vent), 0.38 (hypertension fast); T2D share 0.66 / 0.64
in hypertension (fast / vent); age 0.53 (hypertension vent); others <= 0.28. Comorbidity not removed, as in A233.

## A307 (2026-11-18, user: "Run the glucocorticoid-pattern check for the fetal organs" (lead prompt task C, optional check); staged: null-only power stage, then scoring; **NOT BLIND**: A243, A244, A264, A266, A269 seen in full (in particular: the 14-gene covariation survives adjustment for the GC score in lung and liver, H269.2 PASS; without PER1, PER2, NR1D1 it fails in cortex and lung and is not clock-specific in liver); **seen before writing:** detection and mean counts of FKBP5, TSC22D3, ZBTB16, KLF9 in the three organs' rows (all detected in >= 94 % of rows); NOT seen: any correlation of a single gene with the GC score or between clock genes in the fetal rows): is the fetal clock-gene covariation the glucocorticoid-induction pattern or the waveform's phase relations? **Fetal adrenal cortex (A244 `cortex`, 34 fetuses, 6 cohorts), lung and liver (A243 `lung:all` 38 / 5, `liver:all` 35 / 3)**
**Rows, adjustment.** As A264 / A269: log fractions residualised within each organ on cohort indicators and common
quadratic age (`common.development.adjust`); genes: the 14 GTEx clock genes and the A269 GC targets (FKBP5, TSC22D3,
ZBTB16, KLF9). SI = PER1, PER2, NR1D1 (the signal-induced clock genes of A266); CORE = the other 11.
GC score: mean of the four GC targets' adjusted log fractions standardised over the organ's rows.
**Statistics (per organ).**
- GC link L = mean over the three SI genes of r^2 (Pearson, over fetuses, adjusted) with the GC score. Sign-free: a
  hormone pulse that represses NR1D1 and induces PER1 / PER2 counts as the same pattern. Signed r reported.
- Waveform cross-agreement W = mean over the 33 (SI, CORE) pairs of r_gh cos(phi_g - phi_h) (phi: GTEx waveform
  acrophases): positive when the SI genes covary with the core genes as the clock's phase relations predict (e.g.
  against ARNTL, with DBP); near 0 when the SI genes move on their own.
**Nulls and companions.** L: GC score permuted over fetuses within cohort (10,000; seed [307, organ, 0]); companion: 200
expression-matched non-clock 3-gene sets in place of SI (`paired.pseudo_waveforms` on the SI sub-waveform; sets
containing a GC target redrawn; seed [307, organ, 1]). W: every gene permuted on its own within cohort (10,000; seed
[307, organ, 2], as `development.coherence`); companion: 200 expression-matched sets on the full waveform (A243 matched
family, seed [307, organ, 3]) with the same SI / CORE positions. The responsive family is not used (its pool is the GC
targets themselves).
**Tests (Holm over the three organs per statistic).** H307.1 (GC pattern): L > null (p < 0.05) and above the matched
q95. H307.2 (waveform pattern): W > null (p < 0.05) and above the matched q95.
**Readings (per organ).** H307.1 pass, H307.2 fail: the covariation is the glucocorticoid-induction pattern (the
signal-induced clock genes move with GC targets, not with the core loop): a driven response. H307.2 pass, H307.1 fail:
the signal-induced genes keep the clock's phase relations to the core genes, independent of GC-target activity: a loop
(entrained or running) rather than a GC-driven state. Both pass: a GC-driven component on top of loop-consistent
covariation (not separable). Neither: not resolved at the stated MDE.
**Power stage (null only, before scoring).** Per organ and statistic: null mean, SD, q95 and MDE = q95 + 0.84 SD; matched
companion medians and q95 are part of scoring (they are on the true rows) and are not computed at this stage. No adult
power is computed (the user asked for the cheap check): a fail is read only against the MDE, not as an exclusion of a
loop of a given adult size.
**Code.** `exploratory/A307_fetal_gc_pattern/scripts/score_a307.py --stage power | score`.

### A307 power note (2026-11-18, nulls only; no observed statistic computed)
`score_a307.py --stage power` (17 s; `output/power.csv`, 10,000 draws). MDE (q95 + 0.84 SD of the null): L (mean r^2 of
SI with the GC score) cortex 0.122, lung 0.096, liver 0.104 (null means 0.028-0.035); W cortex 0.058, lung 0.051, liver
0.053 (null means 0.000). Both statistics are scored in all three organs as registered.

## A308 (2026-11-18, user: "yes register it and run the calibration", after re-reading the psoriasis caveats against the newer tests; a calibration and a variance decomposition of an existing statistic, no new data; **NOT BLIND**: A183 in full (psoriasis lesional offset -1.45 h, p 0.0001; matched q05 -1.56, clock at the 5.5th percentile, 11 of 200 sets at or below it; responsive q05 -1.01, 0.5th percentile; matched offsets SD 1.04 h; gene sensitivity, without PER1 and NR1D1 -1.08 h), A54-A82 (archive), A254, A285 (the matched q95 rule is calibrated in bulk, the responsive one is not where few pairs); NOT seen: any swap-placebo offset of a companion set in GSE154200): is the psoriatic lesion delay beyond the lesion-wide rotation of random expression-matched genes, once the companions' own estimation noise is removed? **GSE154200, the 81 psoriasis persons and batches of A183 H183.1**
**Statistic (A183 H183.1, unchanged).** Pooled offset over the used batches of `clocksync.paired` (per-batch
`ccorr` of lesional vs non-lesional states on one projector, weights = the clock fit's inverse variances, fixed), for
the clock and the 200 matched and 200 responsive sets of A183 (same seeds, so the same gene sets). States centred within
sample type per batch as in `paired.swap_offset_null`.
**Placebo.** Per-person label swaps (lesional <-> non-lesional), 2000 draws, seed [308, d]; each draw applies the same
swaps to the clock and to every companion set. Draw 0 = no swap, must reproduce A183 (clock -1.453 h, matched q05
-1.560, percentile 0.055; responsive q05 -1.012) to 3 decimals.
- **V308 (A285 form, descriptive):** share of draws in which the clock offset is below the matched q05 of the same
  draw (calibrated if <= 0.0685; conservative if far below 0.05). Responsive family reported. This checks the rule under
  the complete null (no lesion effect on any gene); it does not test specificity, since the swaps remove the lesion's
  rotation of the companion sets too.
- **H308 (specificity, deconvolved).** Per matched set j: observed offset o_j (draw 0) and its noise SD s_j (SD of its
  wrapped placebo offsets); for the clock, o_c and s_c. True lesion-wide rotation of random genes: mean m = mean_j o_j,
  variance tau^2 = max(0, var_j(o_j) - mean_j s_j^2). Statistic: q = Phi((o_c - m) / sqrt(tau^2 + s_c^2)), the
  percentile of the clock among noise-free random sets observed with the clock's own noise. **H308 passes if q < 0.05.**
  Interval: 90 % over 2000 bootstrap resamples of the 200 sets (reported, not decisive). Also reported: noise share
  mean_j s_j^2 / var_j(o_j); the same for the responsive family.
**What each outcome means (stated now).** H308 passes: the matched sets' spread is mostly their own estimation noise;
the psoriatic delay exceeds what the lesion does to random genes of the same expression: clock-specific, and A183's
5.5th percentile reflects noisy companions. H308 fails with tau^2 well above 0: the lesion rotates random genes'
projections by about as much as the clock; the delay is a lesion-wide transcriptional rotation that the clock genes
share, not shown to be a clock property (the slide says so). V308 well below 0.05: the q05 rule is conservative in this
design under the complete null (stated; it does not by itself change the A183 verdict).
**Limits.** tau^2 is a method-of-moments estimate with a normal approximation; circular wrapping is ignored beyond
|offset| < 12 h (the observed sets lie within about +-3 h). Composition is not addressed.
**Code.** `exploratory/A308_psoriasis_calibration/scripts/score_a308.py` (imports A183's loaders).

## A309 (2026-11-18, user: "yes go ahead, register and run it", after the question why the fetal heart is negative when mouse heart Per1 / Per2 are reported to oscillate; staged: adult power first, then scoring of powered readouts; **NOT BLIND**: A236 (whole-embryo late window: time-locked harmonics not in the clock plane; H108.1 R 0.516 not clock-specific; descriptive rotation estimates hepatocytes 0.67 (se 0.91), cardiomyocytes 0.52 (se 0.31)), A108 / A108c lineage gates (mostly "too few"), the human fetal results A243-A269 and A307; **seen before writing:** embryo counts, cell counts and median clock-gene counts per cell type (counts only); NOT seen: any between-embryo clock statistic of a single mouse cell type): does the late mouse embryo show the human fetal organ pattern (adrenal cortex, lung, liver covary; heart, kidney do not), at the stage where the mouse heart clock is reported to function? **Qiu 2024 late window (E15.0-E18.75: 4 batches = embryonic days, one embryo at 00, 06, 12, 18 h each; `natexp/qiu2024/genome_late`)**
**Readouts (fixed now; pseudobulks with >= 50 cells).** `adrenal` = Adrenocortical cells (13 embryos), `liver` =
Hepatocytes (16), `lung` = Lung progenitor cells (16), `kidney` = Nephron progenitors (16), `heart` = Ventricular
cardiomyocytes (12). Atrial cardiomyocytes (6 embryos) and alveolar cells (<= 8 in the window) are not used.
**Statistic.** As A243 (`common.development.l1_test`): signed coherence C of the 14 clock genes across embryos after
residualising every gene on batch indicators (the embryonic day; within a batch the four embryos span 18 h of age),
10,000 per-gene permutations within batch; 200 matched and 200 responsive companion sets from the genome-wide late
table (gene symbols upper-cased, clock genes by the A236 panel; duplicated symbols: most expressed kept), adjusted in the
same way. Sampling times are not used (the readout is the human one; times enter only through the design).
**Power first (adult stand-ins, before any embryo statistic).** Every embryo of a readout stood for by one timed adult
mouse sample of the organ (distinct samples): Wolff 2023 (`mouse_ageing_wolff`, all ages; adrenal, lung, kidney,
heart) and Manella 2021 liver (`manella2021`, AlbCre, ad libitum; 24 samples); stand-ins assigned to the readout's batch
structure at random; clock-gene counts thinned gene by gene to the cell type's median count (p_g = min(1, fetal median_g /
adult median_g)); adult clock at 1, 0.75, 0.5 (`simulate.scale_clock`) and shuffled (false-positive rate); 200 subsamples,
500 permutations within batch. A readout is scored only if its adult pass rate at the full clock is >= 0.8 (D1b rule).
**H309.1 (per readout, Holm over the scored readouts).** PASS: Holm p < 0.05 and above both companions' q95; "not
clock-specific": p < 0.05 below a companion q95; fail otherwise (read against the pass rate at 0.75 and 0.5).
**Readings.** `heart` PASS: the mouse heart carries the covariation at E15-E18.75, so the human heart negative (4-25
weeks) is a matter of stage (or species), not an organ property. `heart` fail with a pass rate >= 0.8 at 0.5 while
adrenal / liver / lung pass: the organ split replicates in a second species at a stage where heart cells can oscillate
(explants): in vivo the heart is not synchronised between individuals the way adrenal, liver and lung are. `heart` not
powered: the question stays open. Kidney is read in the same way as the heart.
**H309.2 (secondary; only for readouts passing H309.1; Holm over them per statistic).** The A307 statistics (L: mean r^2
of PER1, PER2, NR1D1 with the GC score of FKBP5, TSC22D3, ZBTB16, KLF9; W: (signal-induced, core) waveform agreement),
nulls and matched companions as A307 (within batch). Readings as A307.
**Descriptive.** All-age readouts of the same cell types from the 14-gene table (`mouse_qiu`; quadratic E-day, global
permutations, no companions): C and p only.
**Code.** `exploratory/A309_mouse_embryo_organs/scripts/{power_a309.py, score_a309.py}`.

### A309 power note (2026-11-18, adult stand-ins only; no embryo statistic computed)
`power_a309.py --n-sub 200 --jobs 4` (421 s; `output/power_{design,subsamples,passrates}.csv`). Embryos (batches):
adrenal 13 (4), liver 16 (4), lung 16 (4), kidney 16 (4), heart 12 (4; one batch with a single embryo). Adult pass rate
at the full clock / 0.75 / 0.5 (shuffled): adrenal 1.000 / 0.975 / 0.645 (0.030); liver 1.000 / 1.000 / 0.945 (0.040);
lung 1.000 / 1.000 / 0.960 (0.050); kidney 1.000 / 1.000 / 0.930 (0.065); heart 1.000 / 0.970 / 0.615 (0.050). All five
readouts are scored (Holm over 5). Heart and adrenal: a fail excludes a clock of three quarters of the adult size, not
one of half. Choice made in the script, not in the registration: stand-ins carry the adult sample's own library size
(gene-wise thinning is a constant log offset per gene).

## A310 (2026-11-18, user: "yes go ahead with Stage 1" of `docs/prompts/clock_model_2026-11-18.md`; a model synthesis, staged: fit and adequacy (Stage 2a, **not blind, descriptive**), then one prediction on a GTEx statistic not yet computed (Stage 2b), then the single-cell power stage (Stage 3, separate amendment); **NOT BLIND** for every input of the fit: A104 table, A291 and its post hoc, A294b (slope of log excess own variance on log S 0.223, 90 % 0.044 to 0.332), A302, A293 Parts 1-3, A295, A297 (all four contrasts), A304, A183 / A308, A164-A306, A132; literature `docs/literature/clock_population_models_2026-11-18.md`; **seen before writing (model only, no data):** the model's slope of log departure variance on log R^2 and its radial / tangential split at arbitrary parameters (organs differing in dispersion: slope about 0.3 and tangential excess growing with R; in coupling: about 0.1, same sign; in drive strength: negative slope); NOT seen: the radial / tangential split of any GTEx tissue's own variation): can one driven population-of-oscillators model reproduce the organisation results with far fewer parameters than observations, and does it predict a new property of the personal departure?

**Model (code `exploratory/A310_clock_model/scripts/model.py`).** Cells of class c in tissue o are phase oscillators
of fixed single-cell amplitude a_c, Lorentzian intrinsic frequencies (centre omega0, half-width gamma_o), mean-field
coupling K_o and a systemic drive F (the body's time theta_p). In the drive's frame the tissue order parameter
w_o = R_o e^{i phi_o} follows the Ott-Antonsen equation dw/dt = (i nu - gamma) w + (H - conj(H) w^2) / 2, H = K w + F.
The bulk clock state of class c in tissue o of person p is

    x_ocp = a_c w_o* e^{i theta_p} + a_c delta_op + measurement,

with w_o* the stable locked state (gain a_c R_o, phase phi_o) and delta_op the person's departure: the stationary
fluctuation of w around w_o* under additive white noise of one intensity q for all tissues (J C + C J^T + q I = 0).
A sample mixes classes: G_o = R_o sum_c s_oc a_c, s_oc the share of clock-gene mRNA from class c (lymphocytes a = 0).
Consequences: coupling = S / (S + 2 psi) with S = G_o^2 V; the departure's angle shrinks with gain; its size and
radial / tangential split follow from (gamma_o, K_o, F), not from a free term. Liver feeding (A297) adds a second
drive F_f e^{i alpha_f} to the liver equation (removed: F_f = 0; inverted: alpha_f + pi).
**Organ variants (one choice, counted as a parameter).** Tissues differ in one of: dispersion gamma_o (primary), coupling
K_o, drive F_o; the other two are shared. In each, R_o is the tissue's input, set from S_o.
**Alternatives it must beat.** (a) Linear damped driven oscillator dw/dt = (i nu - lam) w + F (A297's additive-cue
reading): its departure is isotropic, its slope of log departure variance on log S is 0.5 if tissues differ in damping
and 0 if in drive. (b) Composition alone (R_o = 1): organ gain only through s_oc. (c) Free per-organ gains without
structure.

**Observations, parameters (fit, Stage 2a).**
| block | observations | model parameters |
|---|---|---|
| A. GTEx periphery (A104 rows, 27 tissues) | 27 excess own variances E_o = 2 psi_clock - median 2 psi_matched (A294b sets), each with a donor-bootstrap 90 % interval (200 resamples, matched floor fixed); the slope of log E_o on log S_o | q, the amplitude scale a (R of the strongest tissue), the ratio of the two shared constants, the variant: 4. Inputs, not fitted: R_o from S_o (27), phi_o (27) |
| B. Cell classes (A293 Parts 2-3) | 47 atlas x family log gains relative to endothelium with their SEs; 12 atlas x class-pair phase offsets | 5 class gains log(a_c / a_endo); 1 endothelial lead (others share the organ phase): 6 |
| C. Liver drive (A297) | 8: amplitude ratio and shift of Manella DF / AL, Xin DRF / NRF, Greenwell arrhythmic / AL and night-restricted / AL, each with its 90 % interval | F_f / F_s, alpha_f, gamma / F_s, K / F_s (nu = 0): 4. Night-restricted = ad libitum in the model |
Total 14 parameters for 95 observations (inputs excluded). Organ gain and organ phase are inputs: on them the model is
alternative (c) and claims no compression; stated, not tested.

**Adequacy criterion (registered).** Global: parameters < observations / 3. Block A: the prediction lies in the
tissue's 90 % interval of E_o in >= 80 % of tissues, and the model's slope of log E on log S over the observed S range
lies in 0.044 to 0.332; the variant is chosen as the one meeting both with the smallest weighted squared error on
log E_o (if none, block A fails). Block B: >= 80 % of the 47 relative gains within pooled +- 1.645 x their own SE; RMS
of the 12 class-pair offsets around the model (endothelium leads by one fitted lag, the others 0) <= 0.6 h. Block C:
>= 6 of the 8 predictions inside their 90 % intervals. Fits by weighted least squares (log scale for gains, variances
and amplitude ratios; hours for shifts); identifiability by profile: a parameter whose profile does not rise by 2.71
(chi-square, 1 df, 90 %) on both sides within a factor 100 of its estimate is reported not identifiable and never
quoted as a biological value.
**Comparison with the alternatives (descriptive).** (a) Fitted to blocks A and C with its own parameters (3 for C;
its slope fixed at 0 or 0.5 for A): number of block-C hits and whether its slope lies in the A294b interval. (b) For
each tissue, the largest share of clock-carrying mRNA compatible with composition alone, s_o <= (G_o / G_max) x
(a_max / a_min) with a_max / a_min from A293 Part 3 (the strongest tissue fully clock-carrying): tissues whose bound
falls below 0.5 while their dominant cell class carries a clock (muscle, hepatocytes, adipocytes, squamous epithelium,
A293) would need more than R_o = 1. (c) Compression counted per block.
**Verdicts.** ADEQUATE: global rule and blocks A, B, C pass. PARTIAL: A and B pass, C fails: an organisation model
without the liver's drive dynamics. INADEQUATE: A or B fails: the model is not worth it; stop (reported as such).

**Stage 2b, prediction (only if A and B pass; scored once).** The model predicts how the personal departure is
oriented: excess own variation across the tissue's state (tangential, a phase departure) exceeds that along it
(radial, a gain departure), the more so the higher R; alternative (a) predicts no orientation; a person x organ gain
departure (as the person-level amplitude of Saint-Antoine and Anafi 2025) predicts a radial excess.
Statistic: A104 rows and fit (fast deaths, periphery, every tissue an anchor, leave-one-out host). Per tissue o and
donor p, residual r = beta_op - z_o u_p (u_p the host posterior mean), radial and tangential parts on the direction
of z_o u_p; T_o = var(tangential) - var(radial), minus the mean T_o of 100 `clocksync.simulate.like` replicates at the
tissue's estimate (isotropic own variation, real design); share sh_o = T_o / E_o. Companions: the same on 100
matched sets (their own fits), reported per tissue and pooled.
- **H310.1:** median over tissues of sh_o > 0, tissue-bootstrap 90 % lower bound > 0, and above the q95 of the same
  median over the matched sets.
- **H310.2:** Spearman rho(sh_o, S_o) > 0, one-sided permutation p < 0.05 (2000).
- Reported with them: the chosen variant's predicted median share (computed and committed before scoring) and
  whether it lies in the observed interval; the radial excess alone.
- **Power first (committed before scoring):** detection rate of H310.1 in 200 simulations at the real design with a
  tangential excess of the variant's predicted size added to the `like` states; below 0.8 the test is reported not
  powered and not scored.
**Readings.** H310.1 PASS: the departure is mostly a phase departure, as a synchronising population predicts and a
linear damped oscillator does not: the model makes a confirmed prediction beyond re-description (with H310.2 PASS, one
that scales with synchrony as predicted). FAIL, powered, share near 0: departures are isotropic: the population
dynamics are not supported in GTEx and the model stays a description. Radial excess (upper bound < 0): person x organ
gain variation dominates the departure, outside the model.
**Disease and brain death (descriptive, parameter changes, no fit).** Tumour gain loss with no departure rise (A304):
lower a or R at fixed q predicts a departure that does not rise (A304's pure-gain simulation, Q_R -0.10). Psoriatic delay
with cue genes on time (A183, A308): a change of phi_o, from detuning or coupling; whether it comes with a gain
change is read from the A183 states. Chronic conditions (A164-A306): no parameter change. Brain death: F -> 0; a tissue
with K_o > 2 gamma_o keeps R and free-runs (organ arrangement drifts, gain kept); with K_o < 2 gamma_o it dephases
(gain lost). Correction to the task prompt, from `A132_brain_death/output/a132_groups.csv`: in A132 H132.2 the
longest-ventilation tertile has the lower organ-arrangement rho (-0.089, p 0.64) and the shortest the higher (0.247,
p 0.16), neither significant: the direction a drift since loss of drive predicts, not against it; ventilation time is
not time since brain death.
**Limits.** Phase oscillators fix the single-cell amplitude; coupling-dependent cell amplitude (Webb 2012) is outside
the model. The departure is one stationary noise of one intensity; static person-specific cue timing is a different
process. GTEx deaths; the fit is a synthesis, not evidence; liver drive from three mouse laboratories.
**Code.** `exploratory/A310_clock_model/scripts/{model.py, fit_a310.py, predict_a310.py, power_a310.py, score_a310.py}`.

## A311 (2026-11-18, user: "This fit well the lab typical view of the circadian cycle, let's see if we can justify this model", on the driven-Hopf model statement `docs/MODEL_hopf_2026-11-18.md`; staged: power first, then one scoring run; **NOT BLIND** for the A104 table (S, 2 psi, couplings per tissue), A291, A294b (slope 0.223), A302, A310 (fit and verdict); **seen before writing (model only):** the forced Hopf normal form's orientation share (0.015 at mu -4, 0.500 at mu 0, 0.904 at mu +4 in units of F^{2/3}), and that the new `clocksync.orientation` estimator reads no excess under isotropic own variation and recovers an injected tilt on simulated cohorts; NOT seen: the radial / tangential split of any real compartment): is the personal departure of peripheral tissues mostly a phase departure, as for driven cycles above their Hopf bifurcation, and more so in tissues with more gain? **GTEx fast deaths, the 27 peripheral tissues of A104**

**Statistic (`clocksync.orientation`, clocksync f9ea9f1).** A104 rows and fit (every peripheral tissue an anchor,
leave-one-out host). Per tissue: residual of each donor's state on the host posterior mean (loading and intercepts
refitted), split on the direction of the predicted state; tilt T = mean(tangential^2) - mean(radial^2), minus the mean
of 200 simulations at the tissue's real design with isotropic own variation of the estimated size (seed [311, tissue
index]); **share_o = excess tilt / 2 psi_o** (tangential minus radial share of the tissue's own variation). Companion:
100 expression-matched pseudo-clock sets (A294b construction, seed [294, 1], so the same sets), each on its own fit,
50 simulations per tissue.
- **H311.1:** median over the 27 tissues of share_o > 0: tissue-bootstrap (2000) 90 % lower bound > 0, and the
  median above the q95 of the same median over the 100 matched sets.
- **H311.2:** Spearman rho(share_o, S_o) > 0, one-sided permutation p < 0.05 (2000).
- **Reported (descriptive):** per tissue share_o and its null SD; the share of the clock-specific departure,
  share_o x 2 psi_o / E_o (E_o the excess over the matched floor, A294b), read on the model grid as a position
  relative to the bifurcation (> 0.5: above, for small detuning); radial and tangential excess separately; the
  matched sets' per-tissue shares.
**Power (run and committed before scoring).** Per tissue, `orientation.simulator` at the real design with a true tilt
of 0.3 x E_o (a clock-specific departure 30 % more tangential than radial, mu about -0.7 F^{2/3} on the model grid,
i.e. slightly below the bifurcation) and with 0; 200 replicates of the 27-tissue median and its bootstrap bound. Power of
the bootstrap part of H311.1 >= 0.8 at 0.3 x E_o is required to score; otherwise the test is reported not powered and
not scored. Type I rate under tilt 0 reported (must be <= 0.10). The matched-set part is not simulated.
**Readings (stated now).**
- **H311.1 and H311.2 PASS:** departures are mostly in phase, more so where the gain is high, as driven cycles above
  their bifurcation predict and as neither a damped responder (isotropic) nor phase noise in the input (tangential in
  every tissue, departure growing with S at slope 1, excluded by A294b's 0.223) does. The Hopf model makes a confirmed
  prediction; the share of the clock-specific departure places tissues relative to the bifurcation (descriptive).
- **H311.1 PASS, H311.2 FAIL:** phase departures without the predicted gradient; consistent with the model only if
  all tissues sit at a similar distance above the bifurcation; phase noise in the input not excluded.
- **H311.1 FAIL (powered), share near 0:** isotropic departures: no evidence of limit-cycle dynamics in the personal
  departure; the Hopf statement stays an interpretation, not a supported model of human tissue clocks.
- **Median share < 0 (upper bound < 0):** radial departures dominate: person x organ variation of gain, outside the model.
**Limits.** The direction is the predicted state's, from an estimated loading and host; the simulation calibrates its
error at the real design but assumes Gaussian own variation. Non-count technical variation is part of 2 psi and
dilutes the share (isotropic if technical). GTEx deaths. The model grid assumes additive isotropic perturbations of
one compartment; a mixture of compartments at different mu reads as their weighted mixture.
**Code.** `exploratory/A311_hopf/scripts/{hopf.py, power_a311.py, score_a311.py}`.

### A311 power note (2026-11-18, before scoring; `power_a311.py`, 42 s; only simulated statistics computed)
Bootstrap part of H311.1 at the real design: rejection 0.05 under isotropic own variation (type I, <= 0.10 required);
detection 1.00 with a true tilt of 0.3 x E_o (200 replicates each). The estimated median share averages 0.164 against
0.191 put in (attenuated by about 14 %; the scored shares are read as slightly conservative). Powered: H311.1 and
H311.2 are scored.
<!-- src: paper/exploratory/A311_hopf/output/{power,power_tissues}.csv -->

### A311 deviation note (2026-11-18, during the scoring run, before any orientation statistic was printed)
`clocksync.orientation` reports own = 2 x the ML psi; A104 / A294b's 2 psi is small-sample corrected (`on_host`), so
the script's check printed max |diff| 0.159 against A294b. The registered denominator is the A104 2 psi: the primary
share is recomputed as excess tilt / A294b `own` (the script's ML-based share is reported alongside). H311.1 and
H311.2 are scored on the primary share; the null simulations use the ML psi as the isotropic own size (as in the power
stage).

## A312 (2026-11-18, user: "yes register prediction 1 and run it", prediction 2 of `docs/MODEL_hopf_2026-11-18.md`; staged: leverage and power first, then one scoring run; **NOT BLIND** for the A104 table, A291, A294b, A310, A311 (the personal departure is not phase-dominated); **seen before writing (model and simulated cohorts only):** the forced Hopf normal form's elasticity of gain to drive (0.971 at mu -4, 0.333 at mu 0, 0.050 at mu +4, F^{2/3} units; `A311_hopf/output/hopf_grid.csv`), and that `clocksync.elasticity` reads 1.06 / 0.57 / 0.07 for simulated elasticities 1 / 0.5 / 0, with leverage near 0 when persons do not differ in amplitude; NOT seen: any elasticity, leverage or host radius of GTEx): do tissues with less gain follow the persons' shared clock amplitude more steeply than tissues with more gain, as receivers nearer or below their Hopf bifurcation would? **GTEx fast deaths, the 27 peripheral tissues of A104**

**Why this tests the model.** Near or below its bifurcation a receiver's amplitude is sensitive to any person-level
systemic change (drive strength, or a shared shift of mu); far above it the amplitude saturates. If organ gains
differed only by composition or by linear gains on one drive (alternatives (b), (c) of A310, and a linear receiver
model), every tissue would follow the persons' shared amplitude with the same relative elasticity, whatever its gain.
The Hopf model predicts a gradient: relative elasticity falling with gain.
**Statistic (`clocksync.elasticity`, clocksync eed9786).** A104 rows (every peripheral tissue an anchor). Per target
tissue k: 20 random halvings of the other 26 tissues (seed [312, k]); per halving and role swap, host states of the two
halves (`Dataset.host`), the target's amplitude a_k along its predicted direction on host B, and the radii r_A, r_B,
all residualised on the first two harmonics of both hosts' phases; elasticity_k = mean cov(a_k, r_A) / mean(a_k) over
mean cov(r_B, r_A) / mean(r_B) (relative to the shared amplitude of the other tissues); leverage_k = median
corr(r_A, r_B). Donor bootstrap (200 draws of the common donors, hosts fixed) for each tissue's SE.
- **H312.1:** Spearman rho(elasticity_k, S_k) < 0 over the 27 tissues, one-sided tissue-label permutation p < 0.05
  (2000).
- **Reported (descriptive):** per-tissue elasticity with bootstrap 90 % interval and leverage; median elasticity of
  the 9 lowest- and 9 highest-S tissues and their ratio; on the model grid, the ratio places the two groups relative
  to the bifurcation (e.g. ratio 3.4 for mu -1 against +1).
**Stage 1, leverage and power (run and committed before scoring; writes leverages and SEs only, no elasticity).**
Gate: median leverage over tissues >= 0.15 (persons differ in shared amplitude beyond estimation noise) and power of
H312.1 >= 0.8 against a relative elasticity falling log-linearly in S from 2.0 (lowest S) to 1.0 (highest), drawn with
each tissue's bootstrap SE (2000 simulated sets of 27). If the gate fails: reported not testable, not scored.
**Readings (stated now).** PASS: low-gain tissues follow the persons' shared amplitude more steeply than high-gain
tissues, as the Hopf model predicts and composition or linear gains do not: gain differences between organs are, at
least in part, distance to the bifurcation. FAIL (powered): a flat profile; organ gain differences are compatible with
composition or linear gains; the bifurcation stays an interpretation (fetal, lymphocytes, liver saturation) without
support from adult organ gains. A positive rho: high-gain tissues follow more steeply, against the model.
**Limits.** The persons' shared amplitude in deaths mixes drive strength, shared mu shifts and donor-level technical
scale (RNA quality, agonal state); a technical scale acting proportionally on every tissue does not create a gradient,
one acting more on noisy low-gain tissues could; the bootstrap holds the hosts fixed. Elasticities are relative to the
other tissues' composite, which is dominated by high-gain tissues.
**Code.** `exploratory/A312_elasticity/scripts/{power_a312.py, score_a312.py}`.

### A312 power note (2026-11-18, before scoring; `score_a312.py --stage power`; leverages and SEs only)
Leverage (corr of the two half-host radii over donors): median 0.783 over the 27 tissues (0.747 to 0.857): persons
differ strongly in the amplitude with which their tissues jointly express the shared rhythm. Bootstrap SE of the
elasticity: median 0.164 (0.088 thyroid to 0.670 whole blood). Power of H312.1 against the registered alternative
0.984. Gate passed: H312.1 is scored.
<!-- src: paper/exploratory/A312_elasticity/output/{power,power_tissues}.csv -->

## A313 (2026-11-18, user: "yes register A313 and run the power stage", after the question whether cell populations of one organ are more synchronised than populations of another organ at the same physical distance; staged: power first (matched sets only), then scoring on the user's go; **NOT BLIND** for the GTEx organ offsets (A296; descriptive look of 2026-11-18: adjacent same-type pairs 0.00-0.15 h apart, oesophageal mucosa vs muscularis 1.10 h, median of all pairs 0.85 h), A302 (the liver's 21 pair correlations, 0.03 to 0.26), A293 Part 2, A310-A312; NOT seen: the departure correlation of any pair not involving the liver): do anatomically adjacent tissues share a person's departure from the body's time more than distant tissues, at equal composition? **GTEx fast deaths, A104 periphery rows, Whole Blood and Testis excluded (as A302): 25 tissues**

**Why.** In the entrained-ensemble reading (cell-autonomous clocks entrained by one systemic drive, no coupling), a
person's departures in two tissues share only what the systemic drive and sample handling share. Coupling between
neighbouring cells, or a cue local to a body region (innervation, perfusion, temperature), makes adjacent tissues
depart together. Similar composition can also make tissues depart together (shared technical or cell-type
behaviour). The 2 x 2 design separates adjacency from composition. Coupling and local cues are not separated.
**Regions (adjacent = same region; fixed now).** leg: Skin - Sun Exposed (Lower leg), Adipose - Subcutaneous, Muscle -
Skeletal, Artery - Tibial, Nerve - Tibial. heart: Heart - Atrial Appendage, Heart - Left Ventricle, Artery - Coronary,
Artery - Aorta. upper GI: Esophagus - Mucosa, Esophagus - Muscularis, Esophagus - Gastroesophageal Junction, Stomach.
abdomen: Colon - Transverse, Colon - Sigmoid, Adipose - Visceral (Omentum). Every other tissue is its own region.
**Families (similar = same family; fixed now).** squamous: both skins, Esophagus - Mucosa. adipose: both adiposes,
Breast - Mammary Tissue. artery: the three arteries. heart: the two heart tissues. GI smooth muscle: Esophagus -
Muscularis, Esophagus - Gastroesophageal Junction, Colon - Sigmoid. GI mucosa: Stomach, Colon - Transverse. Every
other tissue is its own family.
**Cells.** AS adjacent-similar (3 pairs: atrial-ventricle, aorta-coronary, muscularis-GEJ); AD adjacent-different
(22: 10 leg, 4 heart, 5 upper GI, 3 abdomen); DS distant-similar (11); DD distant-different (the remaining pairs).
Pairs with >= 30 shared donors.
**Pair statistic (A302 estimator).** Host from the periphery tissues except both members; each state regressed on it
with strata intercepts; corrected |correlation| r of the residuals (host-error leakage removed).
- **H313.1 (local factor):** T_adj = median r(AD) - median r(DD) > 0 with (i) p < 0.05 against 10,000 random
  reassignments of the 25 tissues to regions of the same sizes (families kept), and (ii) T_adj above the q95 of the
  same statistic over 100 expression-matched pseudo-clock sets (A302 construction, seed [313, 1]) on the same donors.
- **H313.2 (composition, secondary):** T_fam = median r(DS) - median r(DD) > 0, same two criteria with families
  reassigned (regions kept).
- **Reported:** r per pair and cell; AS; the matched medians per cell; the phase differences of the A296 offsets per
  cell (descriptive, seen).
**Power stage (before scoring; matched sets only, no clock pair correlation computed).** On 20 matched sets: SD of
T_adj over sets and SD of its region-reassignment null within sets; power of H313.1 against an excess of 0.10 in the
corrected correlation for AD pairs (the scale of A302: liver pairs 0.03 to 0.26, matched median 0.12), power =
Phi(0.10 / SD - 1.645) with SD the larger of the two. Gate: power >= 0.8; otherwise reported as not resolvable in
GTEx and not scored.
**Readings (stated now).** H313.1 PASS: adjacent tissues of different composition share a person's departure beyond
distant tissues and beyond matched genes: a local factor (coupling or a regional cue) acts on peripheral clocks;
the uncoupled entrained-ensemble model needs a regional term. FAIL (powered): no regional sharing above 0.10; departures are
regionally independent, as the uncoupled ensemble with one systemic drive predicts. H313.2 PASS with H313.1 FAIL:
sharing follows composition, not location.
**Limits.** Deaths; the departure includes sample-level technical variation (A291, A311), which can be shared by
tissues processed together, e.g. tissues of one region collected in sequence (no GTEx field orders sample collection;
not controlled). The correction assumes isotropic host error. Region and family sizes are small.
**Code.** `exploratory/A313_local_coupling/scripts/score_a313.py --stage power | score`.

## A314 (2026-11-18, user: "we need to split the ventilated apart in B because they don't behave at all the same way in A", slide R3.1; descriptive, **post hoc split**: chosen after A293 Part 1 / A132 showed the organ arrangement kept in ventilated post-mortem donors (rho 0.690) and lost in brain-dead organ donors (0.161); seen: the pooled Hardy 0 alignment (A27 `walltime_hosts.csv`: R 0.142, p 0.0005, 380 donors) and A132's folded distance to 0 / 12 h (organ donors 0.75 h, post-mortem 0.77 h, fast 1.01 h); NOT seen: the alignment to the hour of death of either subgroup alone): does the body's clock read the hour of death in ventilated donors with a living brain and in brain-dead organ donors? **GTEx Hardy 0 split by COHORT (Postmortem; Organ Donor (OPO)); descriptive, no hypothesis**

**Why.** R3.1 panel B reports the pooled Hardy 0 alignment, which mixes the two groups that differ in panel A. A
descriptive readout per group lets the slide show both on the same rows; it tests nothing.
**Procedure (as A27 `compute_walltime_hosts.py`, unchanged).** Per group: `gtex.load` on the 15 reference tissues
(`_shared.V2_TISSUES`), Hardy 0, clock genes, library-size depth; host fitted without time with every tissue an anchor
(`cs.Dataset`); `Dataset.timed(tod, n_perm=10,000)`: R = |mean exp(i (host phase − ω t_death))| and its permutation p.
Groups: ventilated post-mortem (all donors) and organ donors (all donors); because R under no alignment shrinks with
n, also organ donors at the post-mortem n (20 random draws, seed 314, median R and the range).
**Reported.** R, p, donors per group; the null median of R per group (from 1,000 of the permutations) so the two can
be compared at their own n. No contrast between groups is tested.
**Reading (stated now).** Descriptive only: the slide shows each group's R beside the fast-death value. Whether the
living-brain donors keep more alignment than brain-dead donors is not a registered question; any difference is
labelled post hoc.
**Limits.** Same as A27 (one host per group; recorded time of death); organ donors are about four times as many as
post-mortem ventilated donors.
**Code.** `exploratory/A314_vent_walltime/scripts/score_a314.py` (aggregates only).

### A313 power note (2026-11-18, before scoring; `score_a313.py --stage power`, 20 matched sets only, 4518 s)
SD of T_adj over the 20 sets 0.033, median SD of the region-reassignment null within a set 0.021; power of H313.1
against an excess of 0.10 for adjacent-different pairs 0.91: gate passed. Usable pairs per cell (>= 30 shared donors):
AS 3, AD 19, DS 10, DD 229. Matched sets already show T_adj above 0 on average (mean 0.017, 12 of 20 sets positive), and
spread more than the region-reassignment null: the reassignment criterion alone would be anti-conservative, and the
matched-set criterion (ii) carries the specificity, as registered. Scoring needs the clock and 100 sets (about 6 h at 4
workers).
<!-- src: paper/exploratory/A313_local_coupling/output/{power,power_sets}.csv -->

## A315 (2026-11-18, user: "check if brain-death gains were looked at, then register it"; check done: no amplitude or gain comparison of brain-dead organ donors against ventilated post-mortem donors exists (A132, A132b, A293 Part 1 compared the organ arrangement; A314 the alignment to the hour of death; PI Q28 (i) proposed this table and is on hold); staged: power first (null spreads only), then one scoring run; **NOT BLIND** for A132 / A132b / A293 Part 1 (arrangement 0.161 vs 0.690), A314, H109.1 (shared-covariance excess lower in pooled ventilated donors than in fast deaths) and A231 D231.6 (host spread of pooled ventilated donors about a tenth of fast deaths, separate fits); NOT seen: any within-tissue clock-plane variance of either ventilated group): when the central drive is lost (brain death), do tissue clocks lose amplitude (dephasing of uncoupled cellular clocks) or keep it while their arrangement is lost (coupled or limit-cycle tissue clocks drifting)? **GTEx Hardy 0: brain-dead organ donors (COHORT "Organ Donor (OPO)") vs ventilated post-mortem donors (COHORT "Postmortem"), the 12 A109 body tissues**

**Prediction (entrained-ensemble reading).** Cell-autonomous clocks kept in phase only by the systemic drive dephase
when it stops (order parameter decays at the rate of their period dispersion), so each tissue's bulk clock-plane
amplitude falls; tissues whose clocks hold their amplitude by coupling (a tissue limit cycle) keep it and drift in
phase. Both lose the arrangement (seen). Brain death is the loss of the drive in this reading; ventilation, ICU care
and the absence of meals are common to both groups.
**Data.** A109 `Tissues(200, seed 0)` (GTEx, 12 body tissues, first sample per donor, donors with >= 2 tissues and a
recorded time of death; one projector per tissue fitted over Hardy 0-4), clock states and Poisson variances; the 200
matched sets of that construction; responsive sets (`compute_h109_responsive.add_sets`, 200). Groups as A132 (all
donors of each group, no equal-n draw). Tissues with >= 20 donors in each group.
**Statistic.** Per tissue t and group g: debiased clock-plane variance W_tg = mean |z - mean_g z|^2 - mean Poisson
variance (trace). L_t = log(W_t,OD / W_t,PM) (tissues with both W > 0). Q = median over tissues of L_t. The same for each
matched and responsive set: Q_k.
- **H315.1:** Q < 0 with (i) Q below the q05 of the 200 matched Q_k, and (ii) the 90 % upper bound of Q < 0 (donor
  bootstrap within groups, 1000, the matched and clock rows resampled together).
- **Reported:** L_t per tissue with intervals; Q* = Q - median matched Q_k (clock-specific log ratio) and exp(Q*);
  responsive Q_k distribution; fast deaths against ventilated post-mortem (descriptive, the ICU effect); within organ
  donors, top vs bottom tertile of ventilation duration (A132 cuts 61 / 101 h; descriptive, a weak proxy of time since
  loss of drive, protected phenotype, aggregates only).
**Power stage (before scoring; no group difference of the clock computed).** SD of Q over the 200 matched sets and SD
of the clock's Q under 1000 random relabellings of Hardy 0 donors into groups of the real sizes; power against a
clock-specific variance loss of 30 % (Q* = log 0.7) = Phi(0.357 / SD - 1.645), SD the larger of the two. Gate: power >=
0.8, else reported not resolvable and not scored.
**Readings (stated now).** PASS: tissue clocks of brain-dead donors carry less amplitude than those of ventilated
donors with a living brain, beyond matched genes: loss of the central drive dephases peripheral clocks (uncoupled
entrained ensemble). FAIL with Q* near 0 (powered): amplitude kept while the arrangement is lost: tissue clocks keep
their own rhythm without the drive and drift (coupled or limit-cycle tissue clocks), against the uncoupled reading.
Q > 0 (clock-specific): brain-dead donors carry more clock amplitude, outside both readings.
**Limits.** Procurement differs (organ donors: surgical procurement with perfusion; post-mortem: autopsy), handled only
as far as matched genes capture it; donor management (steroids, hormones) differs; time since brain death is not
recorded; a variance loss can also come from organ-specific additive noise falling (not expected); one projector per
tissue fitted on all death types.
**Code.** `exploratory/A315_braindeath_gain/scripts/score_a315.py --stage power | score`.

### A315 power note (2026-11-18, before scoring; `score_a315.py --stage power`, 1 min; no clock contrast computed)
298 organ donors, 79 post-mortem donors, 12 tissues with >= 20 donors per group. SD of Q over the 200 matched sets
0.082; SD of the clock's Q under 1000 relabellings of Hardy 0 donors 0.084; power against a clock-specific variance loss
of 30 % 0.995. Gate passed: H315.1 is scored.
<!-- src: paper/exploratory/A315_braindeath_gain/output/power.csv -->

## A316 (2026-11-18, user: "yes register it with the feasibility stage first ... is CoPhaser really better than this method, the projection should probably identical": the method's own per-cell estimator is used, not CoPhaser; staged: feasibility on the aligned arm only, then one scoring run; **NOT BLIND** for A201 (pseudobulk shifts SMC -2.85, fibroblast -2.87, EC -2.29 h; cosinor amplitude ratio misaligned / aligned SMC 0.875, fibroblast 0.888) and A279; NOT seen: any within-library cell-to-cell clock statistic of this dataset): after a 6 h light advance, do the cells of a tissue shift as a block, drift apart (dephasing of autonomous cells), or lose amplitude each (damping)? **Auerbach 2026 mouse aorta (A201 data; aligned 14, misaligned 8, Bmal1 KO 10 libraries, pools of 2 aortas)**

**Why.** A315 found that tissue clocks keep their amplitude after loss of the central drive while organs drift apart,
and A201 that aortic cell types shift together. Whether a tissue is a block of coherent cells or an ensemble of
autonomous cells is read in single cells: during re-entrainment, autonomous cells with different periods move at
different speeds and spread; a coherent tissue moves as a block; damped cells lose amplitude each. A component that
paces the others would lead and its followers would spread more.
**Estimator (`clocksync.cells.dispersion`, no per-cell phases).** Per library x cell type with >= 100 cells: D =
within-library cell-to-cell in-plane second moment of the clock state, from cross-gene moments of split gene halves
(free of counting noise; log scale; 30 splits; the 14 clock genes, mouse symbols, GTEx waveform design). Per
condition x cell type: D_bar = mean D over libraries; coherent power A^2 = mean over libraries of |b_l - mean b|^2,
b_l the library's pseudobulk clock state (same projector on the library's log mean expression; ZT-balanced within
condition, sex centred), minus its Poisson variance. Coherent fraction F = A^2 / (A^2 + D_bar); per-cell power
P = A^2 + D_bar. Contrasts misaligned / aligned on the log scale, per cell type and pooled (mean of SMC and fibroblast);
library jackknife within condition for the SE.
- **H316.1 (dephasing):** pooled log(F_mis / F_al) < 0, 90 % upper bound < 0.
- **H316.2 (block):** pooled log(F_mis / F_al) within +-0.26 (90 % interval; 0.26 = the drop if the whole 0.88
  amplitude ratio were dephasing at constant P).
- **Reported:** log P ratio (damping if P falls with F kept); per cell type F and P ratios, EC and macrophages
  descriptive (a pacemaker type would keep F while its followers lose it); Bmal1 KO F and D against the aligned arm
  (negative control: coherent power collapses); clock UMI per cell.
**Feasibility stage (aligned arm only; written before scoring).** (a) median clock UMI per cell >= 2 in SMC and
fibroblasts; (b) aligned D above its permuted-acrophase null (`dispersion_null`, 20 permutations) in both; (c) placebo:
the aligned libraries split into two ZT- and sex-balanced halves (all balanced splits, or 50 random ones), the pooled
log F ratio and its jackknife SE computed per split; SE (median) <= 0.105, so a ratio of exp(-0.26) is detected at about
80 %; placebo mean within +-0.1 of 0. Any criterion failing: power note, not scored.
**Readings (stated now).** H316.1 PASS with P kept: autonomous cells drift apart during re-entrainment (uncoupled
ensemble within the tissue). H316.2 PASS: the tissue shifts as a coherent block (coupling or a tissue oscillator), and
with P down, the 0.88 is per-cell damping. Neither: not resolved. A cell type keeping F while others lose it: candidate
pacemaker (descriptive).
**Limits.** One tissue, one time after the shift (day 5), pools of 2 aortas; D includes non-clock covariance of the
clock genes (null-corrected only in the feasibility check); coherent power mixes rhythm amplitude with pool-to-pool
variation in each arm.
**Code.** `exploratory/A316_aorta_cells/scripts/score_a316.py --stage feasibility | score`.

### A316 feasibility note (2026-11-18, aligned arm only; `score_a316.py --stage feasibility`): gate FAILED, not scored
(a) Clock UMI per cell: median 6 (SMC), 4 (fibroblasts): passes. (b) Aligned D above its permuted-acrophase null in
both (SMC 0.047 vs -0.521; fibroblasts 5.418 vs -0.136): passes, but the null means are negative and the two cell
types differ by two orders of magnitude, so D is not stable enough here to be read as a dispersion. (c) Placebo (50
ZT- and sex-covering halvings of the 14 aligned libraries): mean log F ratio 0.033 (SD 0.209), median jackknife SE
0.288, above the 0.105 required: a drop of exp(-0.26) in the coherent fraction cannot be resolved; the coherent power,
estimated from a handful of pooled libraries per arm, dominates the error (as the pool-to-pool variation in A201).
Registered consequence: H316.1 and H316.2 are not scored; the misaligned arm was not read.
<!-- src: paper/exploratory/A316_aorta_cells/output/{feasibility,feasibility_gate,placebo}.csv -->

## A317 (2026-11-18, user: "yes go ahead with gap 1" of `docs/PAPER_TRIAGE_2026-11-18.md`; confirmatory for the biology paper's headline; staged: power first (label-permutation nulls only), then one scoring run; **NOT BLIND** for H317.1: A132b (post-mortem 0.678 vs organ donors 0.161 against fast deaths, difference 0.517, p 0.008, post hoc), A293 Part 1 (organ donors 0.161: 99th percentile of matched and 75th of responsive sets), A133 post hoc, A136 (glucocorticoid gate failed), A314 (equal loss of alignment to the hour of death), A315 (tissue amplitude kept); H109.1 pooled over all ventilated donors (C_clock 0.463 above the pseudo-clock null); **NOT seen**: C or its excess split by organ donors vs ventilated post-mortem donors, and any of the confound sensitivities below): does brain death remove the clock relation between organs, clock-specifically and beyond procurement confounds? **GTEx Hardy 0, organ donors (298) vs ventilated post-mortem donors (79); fast deaths (244) as the arrangement reference; the 12 A109 body tissues**

**Data.** A109 `Tissues(200, seed 0)`: clock states, 200 expression-matched pseudo-clock sets per tissue, 200
stress-responsive sets (`compute_h109_responsive.add_sets`). Groups as A132 (COHORT). Equal n: organ donors drawn to 79,
20 draws (seed 317); statistics are medians over draws. Pairs with >= 20 shared donors in every group compared (A132 rule).
- **H317.1 (arrangement, not blind).** D_rho = rho(post-mortem, fast) - rho(organ donors, fast), Fisher-Lee correlation of
  pair phases with fast deaths (A132 / A293 estimator). Passes if (i) D_rho > 0 with donor-label permutation p < 0.05
  (2000 relabellings of Hardy 0 donors into groups of 79 and 79, same draws), and (ii) D_rho above the q95 of the same
  difference for the matched sets and for the responsive sets (same donors and draws).
- **H317.2 (shared clock between organs, not seen by group).** Excess C = C_clock - median C over the 200 matched sets
  (A109 H109.1 estimator: mean |Poisson-corrected complex correlation| over pairs). D_C = excess(post-mortem) -
  excess(organ donors). Passes if D_C > 0 with donor-label permutation p < 0.05 (same relabellings). Reported: each
  group's excess with its matched-set percentile, and the responsive-set excess.
- **Confound sensitivities (pre-specified, descriptive; D_rho and D_C with their matched and responsive medians):**
  (a) donors with long-term steroid use excluded (MHSTRDLT = 1); (b) samples restricted to the common range of
  ischaemic time of the two groups (5th to 95th percentile of the smaller group's range); (c) samples with RIN >= 6.
**Prediction (federation reading).** Without the central drive organs drift apart while each keeps its amplitude
(A315): H317.1 and H317.2 pass, robust to (a)-(c).
**Readings (stated now).** Both pass: brain death removes the clock relation between organs, beyond matched and
responsive genes; with A315, organs lose their order, not their clocks. H317.1 passes, H317.2 fails: organs keep a
shared clock covariance in a different arrangement (rearranged, not drifted). H317.1 fails on (ii): the arrangement
loss is not specific to clock genes (a generic organ-donor state). A sensitivity reversing the sign: the confound named
is a candidate explanation and is reported as such.
**Power stage (committed before scoring).** SD of D_rho and D_C under the 2000 relabellings (observed values not
computed); MDE = 2.49 SD. Reported with the A132b difference for reference.
**Limits.** Not blind for the arrangement; procurement (surgical vs autopsy), donor management and time since brain death
differ or are unknown; one cohort. Replication outside GTEx: Tabula Sapiens cannot resolve the arrangement (A109b, 5
multi-organ donors, power 0.06-0.13); a mouse counterpart (organs after loss of the central clock) is a separate step.
**Code.** `exploratory/A317_braindeath_confirm/scripts/score_a317.py --stage power | score`.

### A317 power note (2026-11-18, before scoring; `score_a317.py --stage power`; relabelling nulls only)
65 pairs. SD of D_rho under 2000 relabellings 0.223, MDE 0.554: the A132b difference (0.517, one equal-n draw) lies just
below it, so H317.1 has roughly three-quarter power at the effect already seen. SD of D_C 0.052, MDE 0.130. No stop
rule was registered for this stage; both hypotheses are scored. Choice fixed in the script before scoring: the
sensitivity rows carry companion medians over 20 sets (every 10th) and 5 organ-donor draws, to bound runtime.
<!-- src: paper/exploratory/A317_braindeath_confirm/output/power.csv -->

## A318 (2026-11-18, user: "yes do both in parallel" (the ischaemic-time follow-up of A317 and the mouse data search); staged: power first, then one scoring run; **NOT BLIND** for A317 (D_rho 0.621; sensitivities: steroids 0.575, RIN >= 6 0.591, common ischaemic range 0.239 on 36 of 65 pairs) and everything A317 lists; NOT seen: any clock state adjusted for ischaemic time, any within-group ischaemic gradient of the arrangement): is the brain-death change of the organ arrangement explained by ischaemic time (or RNA quality) rather than by brain death? **A317 data, groups, draws and pairs**

**Why.** Cutting samples to a common ischaemic range halved the pairs and the effect (A317), which cannot separate a
confound from noise. Adjusting each sample's state keeps every pair.
- **H318.1 (primary).** Per tissue, the clock state of every Hardy 0 sample (both groups pooled) is regressed (complex
  least squares, with intercept) on standardised ischaemic time and RIN; residual states replace the states (samples
  without both attributes dropped); the same adjustment for every matched and responsive set. D_rho as A317 H317.1 on
  the adjusted states, same draws, pairs (>= 20 shared donors after the drop) and criteria: D_rho > 0, relabelling p <
  0.05 (2000), above the matched and the responsive q95. Pooled fitting removes any part of a true group difference
  that is collinear with ischaemic time: the test is conservative.
- **H318.2 (descriptive).** Within organ donors, donor-level median ischaemic time in tertiles: rho(top tertile, fast) -
  rho(bottom tertile, fast), unadjusted states, relabelling null within organ donors (2000). If ischaemic time drove the
  contrast, organ donors with longer ischaemic times (closer to the post-mortem range) would keep more of the arrangement
  (positive difference).
- **Reported:** ischaemic time and RIN medians per group (aggregates); the adjusted D_C (as A317 H317.2).
**Readings (stated now).** H318.1 PASS: the arrangement change survives the adjustment: not explained by ischaemic
time or RNA quality. FAIL with D_rho near 0: ischaemic time (or RNA quality) is a sufficient explanation of the
contrast at this resolution, and the brain-death reading loses its GTEx support. FAIL with D_rho clearly > 0 but below
the companions: not decided. H318.2 is read only beside H318.1.
**Power stage.** SD of the adjusted D_rho under 2000 relabellings and the MDE (2.49 SD); pairs remaining. No stop rule:
both are scored, the MDE reported.
**Code.** `exploratory/A317_braindeath_confirm/scripts/score_a318.py --stage power | score`.

### A318 power note (2026-11-18, before scoring; `score_a318.py --stage power`, adjusted relabelling null only)
All 65 pairs kept (298 organ donors, 79 post-mortem). SD of the adjusted D_rho 0.234, MDE 0.582. Confound aggregates
(medians over samples): ischaemic time 258 min (organ donors) vs 561 min (post-mortem); RIN 7.6 vs 7.2.
<!-- src: paper/exploratory/A317_braindeath_confirm/output/a318_power.csv -->
