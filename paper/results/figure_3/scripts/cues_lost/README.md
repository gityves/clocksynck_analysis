# cues_lost: organs keep their own time when daily cues are lost

Ventilated (Hardy 0), slow (Hardy 3-4) and fast (Hardy 1-2) GTEx deaths. Narrative: `docs/WORKPLAN.md` (Figure 4). The
amendments A46, A47, A81/A81b, A89-A93, A100-A102b (gradient, ventilation dose, cue/locked genes, drug targets, covariates) are in
`gtex_cues` and `gene_rhythms`; they feed this figure. New work packages land here:
WP1 A109 (clock-specific nulls and rotating clocks vs shared displacement), WP2 A110 (mechanism), WP3 A111/A112 (feeding, mouse),
WP4 A113 (living humans). Outputs: `paper/results/figure_3/output/cues_lost/`. GTEx time of death is controlled access: aggregates only.

`_shared.py` re-exports `organs/_shared.py`.

| script | amendment, hypothesis | question and estimand | legacy script replaced |
|---|---|---|---|
| `compute_ventilator_tissues.py` | A33 | coupling to the host of the other 14 tissues in ventilated vs fast deaths, food-entrained vs vascular / stromal tissues, D = mean loss (food) - mean loss (vascular), donor bootstrap and tissue-label permutation; with and without the glucocorticoid score | `ventilator_tissues/scripts/a33_ventilator_tissues.py` |
| `compute_ventilator_gc.py` | A33 ("with GC" readout, target-specific adjustment) | the A33 loss and D with the target's state adjusted on its own GC score and the mean GC score of the donor's host tissues (`Dataset.adjust(['gc', 'gcH'])`, one dataset per target), as the legacy covariate | `a33_ventilator_tissues.py` (`coupling_gc`) |
| `compute_walltime_hosts.py` | A27 part 2 (D5 sanity) | does the host (15 tissues, no times) read the recorded time of death? `Dataset.timed` in fast, ventilated, Hardy 3 and Hardy 4 donors, with the gate | sanity lines of `glucocorticoid/scripts/diag5_int_ext.py` |

The per-script deviation notes for these three scripts (A33 GC adjustment etc.) are in `organs/README.md` (deviations section), written before the split.
