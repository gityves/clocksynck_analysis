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
