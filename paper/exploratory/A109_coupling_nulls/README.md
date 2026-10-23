# A109: coupling nulls (GTEx v10)

Registered in `paper/docs/registered_experiments.md` (A109; A109b Tabula Sapiens is a separate builder and script).
`scripts/compute_coupling_nulls.py` scores H109.1 (C_clock against the pseudo-clock null, ventilated donors) and H109.2
(Fisher-Lee correlation of pair phases, ventilated vs fast) with the secondaries (acrophase null, Hardy 3 / 4 / fast,
excess bootstrap, non-temporal direction, kit-differing pairs, circularity simulations). `scripts/_helpers.py`: data
access (`Tissues`) and statistics not yet in `clocksync`.

```
run.sh --scout           group sizes, pair counts, null spread (MDE), power of the Fisher-Lee test (no real C of ventilated donors)
run.sh --subsample 500   end-to-end check (20 sets, 200 permutations), outputs in output/subsample/
run.sh                   full scoring (200 sets, 10,000 label permutations, 200 bootstraps)
```
Outputs (`output/`, aggregates only): `h109_1_null.csv`, `clock_plane_variance.csv`, `h109_2_offsets.csv`,
`pair_relations.csv`, `excess_bootstrap.csv`, `nontemporal_share.csv`, `circularity_sims.csv`, `procedure_check.csv`,
`h109_2_power.csv`, `scout_*.csv`.

Choices not fixed by the registered text: see the report of the scouting step (donor universe, n, folded statistic,
equal-n bootstrap, circularity simulation model, GC direction constant). Controlled access: time of death and
phenotypes stay in memory; kit annotations enter as per-tissue majority only.

Proposed `clocksync` additions (with tests): `paired.ccorr_matrix` (donor x tissue complex covariance and
pseudo-covariance with pairwise-complete donors), `stats.circularity_ratio`, a batched `paired.coherence` with
pair masks, `controls.pseudo_clock_states` (states of pseudo-clock sets), `simulate.shared_state` (rotating / displaced
mixture of the observed design), and a permutation Fisher-Lee test over tissue labels.
