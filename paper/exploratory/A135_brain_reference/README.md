# A135 part 1: control-only reference choice for the brain host

Registered in `paper/docs/registered_experiments.md` (A135). Control donors only; no disease contrast is computed.
`scripts/compute_brain_reference.py` scores per dataset and compartment (1) clock strength (Poisson-corrected
clock-plane share against 1,000 permuted-acrophase planes, as A86), (2) the share of other compartments with which the
pairwise model-free gate passes (`Dataset.gate` on the two-compartment subset), (3) mean pairwise C with the others
(`Dataset.coherence`), then ranks by the rank sum and applies the registered rule. The docstring lists the choices the
registration leaves open (dataset deduplication to 5 independent sets on the largest control set, minimum 15 control
donors per compartment, 15 shared donors per pair, ties on the rank sum broken by centrality).

```
run.sh --subsample 40   check (fewer permutations, 40 control donors per dataset), outputs in output/subsample/
run.sh                  full run (1,000 acrophase planes, 200 gate and 100 coherence permutations)
```
Outputs (`output/`, git-ignored): `a135_criteria.csv`, `a135_ranks.csv`, `a135_decision.csv` (registered rule, then a
post hoc neurons-merged sensitivity), `a135_family_map.csv`, `a135_datasets.csv`, `a135_gate_pairs.csv`.
