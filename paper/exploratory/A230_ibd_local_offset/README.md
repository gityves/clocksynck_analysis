# A230: epithelial clock offset from the stromal host, inflamed vs non-inflamed intestine

Registration: `paper/docs/registered_experiments.md` (A230). GSE266616 blind primary, TAURUS non-blind secondary.

| step | command | output |
|---|---|---|
| data | `paper/common/datasets/ibd_gse266616/build.py` (`unpack()` once, then `build.py gse266616\|taurus <dir>`), run with `python -I` | `output/data/<cohort>/{pseudobulk,donors}.csv`; raw files and `PROVENANCE.md` in `$SCCC_DATA/gse266616/` |
| power (before scoring) | `RUN_CAP_GB=30 paper/common/run.sh paper/exploratory/A230_ibd_local_offset/scripts/compute_a230.py --stage power --cohort gse266616\|taurus --jobs 3 --n-rep 100` | `output/power_<cohort>.csv`, `gate_`, `cells_`, `eligible_samples_`, `swap_null_`, `coupling_mde_curve_` |
| genome-wide TAURUS pseudobulks (epithelial, fibroblast, endothelial; for the pseudo-clock sets) | `python -I paper/common/datasets/ibd_gse266616/build.py taurus x --genome` | `$SCCC_DATA/ibd/TAURUS/pseudobulk_genome/` (28 MB, genes with non-zero total) |
| tests (A230 design change: pseudo-clock check replaces the companion) | `RUN_CAP_GB=30 paper/common/run.sh paper/exploratory/A230_ibd_local_offset/scripts/compute_a230.py --stage tests --cohort taurus --jobs 3` (gse266616 refuses without `--force`) | `output/contrast_taurus.csv`, `pseudo.csv`, `tests.csv` |

Power stage: no observed contrast is computed. The within-patient swap null is taken from `Dataset.contrast` run on a random
swap of the labels (a null draw); only `Contrast.null` is kept, the contrast of that random labelling is discarded. The
tests stage has been smoke-tested on scratch output only (statistics not read).

Choices not fixed by the registration are listed in the docstring of `compute_a230.py` and `build.py`.
