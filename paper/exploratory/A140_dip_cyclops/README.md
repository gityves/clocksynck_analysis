# A140 (CYCLOPS 2.0 donor internal phase vs clocksync): code and results

## Q2 (i), (iii)

Registered in `paper/docs/registered_experiments.md` (A140 Q2 readouts, dated note 2026-10-26). Part (ii) (CYCLOPS
pipeline) is not run here (waits for Q1).

**Run.** `build_q2.py` (`--jobs 4`, log `paper/results/logs/build_q2.log`, peak 1 GB) streams the five GSE266919 objects
with `rds_parse.py` and sums raw UMIs per sample x compartment (`output/q2/pb_*.csv.gz`, `cells_per_sample.csv`,
`samples.csv`; git-ignored). `q2_readouts.py` (`--jobs 4`, log `paper/results/logs/q2_readouts.log`) runs (a) the orientation
check, (b) part (i) on the primary set (44 Pre-treatment biopsies; 43 enter, one has no compartment with >= 20 cells;
strata breast vs other) and on the all-78 set (77 enter), (c) part (iii) only if gate passes and stability >= 0.8.

**Orientation (a).** Simulated timed cohort, 5 compartments, 60 donors, 24 h window, 5 seeds: the resultant length of
(phase - t) is 0.984-0.989 and of (phase + t) 0.025-0.188; a larger `Dataset.host()` phase_h is later (no sign flip).
<!-- src: paper/exploratory/A140_dip_cyclops/output/q2/orientation.csv -->

**Part (i).**

| set | samples | gate p (shared signal) | gate p (clock shape) | clock share (null) | gate | anchors | stability | free agreement |
|---|---|---|---|---|---|---|---|---|
| primary | 43 | 0.000999 | 0.4585 | 0.0246 (0.0283) | fails | B, CD4_T, CD8_T, NK, myeloid (learned, status ok) | 1.000 (median similarity 0.980) | 1.000 |
| all 78 | 77 | 0.000999 | 0.5405 | 0.0249 (0.0311) | fails | same | 1.000 (0.989) | 1.000 |

<!-- src: paper/exploratory/A140_dip_cyclops/output/q2/part_i_summary.csv -->

The pairing p is at the permutation floor (1,000 permutations) but the clock-shape p is not: the shared covariance is
not clock-shaped (clock share below its null mean). The gate fails in both sets, so `Dataset.fit` couplings are withheld
(uninformative; not "decoupled"). The fixed stromal-vascular rule gives no anchor (no such compartment deposited).
Learned families (tissue scope, one breast study in the learning set): myeloid pooled 0.199 (lower 0.113), lymphocyte
0.123 (lower 0.040); both anchor families.
<!-- src: output/q2/primary_families.csv -->
These differ from the note's values (0.202 / 0.136 and 0.091 / 0.027), which came from `a17_anchor_families.csv` (target
organ_breast, all-tissue scope: the breast study itself was the target there); the anchor outcome is the same.
Stability and free agreement are computed on the anchored host regardless of the gate; the free host equals the anchored
host here because all five compartments are anchors, so the agreement of 1.0 is by construction. Primary pair
coherence excess is 0.46 (CD4_T-CD8_T) down to negative (B with CD4_T, NK, myeloid): `output/q2/primary_pairs.csv`.

**Part (iii).** Not estimable (gate failed, both sets): P140.1-3 rows read "not estimable: gate failed".
<!-- src: output/q2/part_iii_readouts.csv -->
The code path (host refit bootstrap, within-stratum permutations, Holm) was exercised on a check run only (40 bootstraps,
100 permutations on the primary set, forced past the gate; not a result and not kept).

### Choices not in the registration
- Compartments with >= 20 cells enter the clock table; the clock genes absent from an object would be set to 0 (none were).
- Depth = sum of all genes of the object (gene sets differ between objects).
- Stability uses the package default threshold 0.7 (share of bootstrap hosts above it); "identified" = share >= 0.8.
- Learning tables rebuilt as `compute_anchors.learn` (gate-passing healthy fits, tissue not "multi"); target pairs from
  `Dataset.coherence(n_perm=200, n_boot=200)` of each set (the all-78 set uses its own pairs).
- Gate seed 0, 1,000 permutations; fit would use n_boot 200, seed 0 (not run).
- Orientation: 60 donors, uniform 24 h window, no jitter, psi 0.05, z = 0.8-1.2 real amplitudes, anchors = all five; sign
  decided by the larger resultant length.
- Part (iii) specifics not exercised (fixed in code): fraction logit = log((n+0.5)/(N-n+0.5)) with N all deposited immune
  cells of the sample; scores z-scored once over the readout's samples (not per bootstrap replicate); bootstrap hosts
  unwrapped about their own circular mean; percentile CIs; interferon-response gene set (ISG15, IFIT1, IFIT2, IFIT3, IFI6,
  IFI44L, MX1, OAS1, STAT1, IRF7) chosen by me, T-cell scores in CD4_T and CD8_T separately, predicted sign > 0; the
  all-78 permutations run within Group x site strata without patient blocks.
