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
