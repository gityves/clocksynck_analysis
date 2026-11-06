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
