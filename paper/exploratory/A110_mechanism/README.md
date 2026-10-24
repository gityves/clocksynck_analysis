# A110: what keeps organs coupled when daily cues are disrupted by illness?

Registered in `paper/docs/registered_experiments.md` (A110). `scripts/compute_mechanism.py` holds the scout and the
scoring of H110.1 (decay of coupling with the interval ill, Hardy 4, `DTHCODD`), H110.2 (amplitude and stroma
signatures of the coupled network, Hardy 4 vs fast), H110.3 (temperature, exclusion), H110.4 (hypothalamus,
Hardy 3-4 vs Hardy 0 organ donors) and H110.5 (digestive tissues, exclusion). Reuses `A109_coupling_nulls/scripts`
(`_helpers`: states and C) and `A132_brain_death` conventions (equal n, fixed pair mask, MDE = 2.49 x null SD);
the clock-plane fits use the `clocksync` public API (`Dataset.fit`, `.contrast`, `.subset`).

```
run.sh --scout           counts, completeness, pair structure, null spreads, memory-model drop (no observed statistic)
run.sh --subsample 500   end-to-end check, outputs in output/subsample/ (delete without reading before scoring)
run.sh --only h1 h4      score only the listed hypotheses (default: all)
run.sh                   full scoring (1,000 donor-label permutations)
```
Outputs (`output/`, aggregates only): `scout_*.csv`, `a110_*.csv`. Controlled access: phenotypes, time of death and
the interval ill stay in memory.
