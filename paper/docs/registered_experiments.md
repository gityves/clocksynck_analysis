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
