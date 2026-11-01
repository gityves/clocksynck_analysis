# C4 fits: CoPhaser cell-cycle phase of malignant cells, remaining cohorts + seed stability (2026-11-01)

QC only, same STRICT rules as the pilot (no clock-gene statistic or plot; clock genes excluded from rhythmic, context and score sets; cell tables have no clock columns).
Script: `../pilot/scripts/fit_cellcycle.py` extended with `--out`, `--seeds`, `--extra`, `--load-only`, a Pal loader and per-seed theta saving (no copy). `scripts/run_all.sh` queue,
`scripts/extract_pal.py` (Pal tumours, `python -I zipfile` into the new folder `$SCCC_DATA/sc_breast_er/GSE161529/hb10x_tumours/`), `scripts/seed_stability.py`.
Same settings as the pilot (<= 3,000 cells/patient, >= 1,000 UMI, auto hyperparameters, 5 seeds x 200 epochs, kept = best non-vetoed seed). Pilot studies were refitted with the same seeds: kept seed, usable-seed count and cell set are identical to the pilot, and kept-seed theta equals the pilot's `cells_<study>.csv.gz` theta (circ. corr 1.000, median diff 0) for all six (the GPU fit is deterministic).

Outputs `output/`: `cells_<study>.csv.gz`, `validation_<study>.png`, `seeds_<study>.csv`, `seed_theta_<study>[_extra].npz` (aligned theta and cycling probability of every usable seed), `qc.csv` (pilot rows, `source=pilot`, plus new rows, `source=fits`; `qc_rerun_all.csv` = all rerun rows), `seed_stability.csv`, `patient_cycling*.csv`.

## Cohort selection
CNV-malignant cells, median UMI > 1,500, >= 3 patients with >= 300 cells (after the 3,000-cell cap, >= 1,000 UMI): ERP149257, GSE159115, GSE166555, GSE181919, GSE200997, GSE212966, GSE270680, GSE294300, GSE315534, GSE335070, SRP416180 (11 new). Not fitted: cohorts with < 3 such patients. Patient counts differ from the scout (h5ad patient labels). Pal GSE161529: 21 ER-labelled tumours, malignant = authors' `Tum` clusters (A170 definition), 43,693 cells, 19 patients >= 300.

## Validity (all 12 new figures inspected)
| study | cells / pts | seeds usable | cycling | S->G2M | judgement |
|---|---|---|---|---|---|
| GSE161529 Pal breast | 43,693 / 21 | 1/5 | 0.30 | holds (1.57) | valid: markers unimodal, S (0.35 pi) < histone peak (0.5 pi) < G2/M (0.8 pi), library 600 -> 1,700; G1 mode dense |
| ERP149257 CRC | 6,101 / 8 | 1/5 | 0.41 | holds | valid, low confidence: only 1 seed, wide bands, histone plateau 0.35-0.55 pi, library 1.0k -> 2.0k |
| GSE166555 CRC | 6,647 / 12 | 5/5 | 0.44 | holds | **suspect**: markers and histone peak fine, but library size peaks at -0.5 pi (G1) and falls through S/G2 (no rise); seeds unstable (circ. corr -0.82 to 0.43) |
| GSE181919 HNSC | 4,905 / 19 | 5/5 | 0.54 | holds | valid: tight bands, histone peak 0.5 pi, library 3.1k -> 7.0k late; 4 seeds corr 0.85-0.91 |
| GSE200997 CRC | 2,673 / 13 | 5/5 | 0.53 | holds | **weak**: markers ordered, histone peak 0.55 pi, but library size ragged/flat, MCM6 noisy, flat phase histogram; seeds unstable (-0.47 to 0.58) |
| GSE294300 CRC | 5,863 / 18 | 4/5 | 0.54 | holds | valid: clean markers, histone peak 0.45 pi, library 200 -> 1,100 |
| GSE315534 CRC | 9,563 / 5 | 5/5 | 0.60 | holds | **suspect**: markers and histone peak (0.35-0.45 pi) fine, but library size falls with phase (650 at -0.5 pi -> 300 at 0.8 pi), opposite to expectation |
| GSE335070 ampullary | 4,950 / 4 | 2/5 | 0.43 | holds | **weak**: histone profile has two peaks (0.45, 0.85 pi), library ragged without rise; 4 patients, 1 comparable seed (corr 0.45) |
| SRP416180 bladder | 8,795 / 7 | 4/5 | 0.56 | holds | valid: histone peak 0.55 pi, library 600 -> 1,200, seeds corr 0.54-0.86 |
| GSE159115 kidney | 6,108 / 7 | **0/5** | 0.05 | - | **invalid, exclude**: all seeds vetoed, 5 % cycling, multimodal phase histogram |
| GSE212966 pancreatic | 4,184 / 6 | **0/5** | 0.25 | - | **invalid, exclude**: all seeds vetoed (library ratio 6.8); markers look ordered but fit not accepted |
| GSE270680 gastric | 2,121 / 18 | **0/5** | 0.82 | - | **invalid, exclude**: all vetoed, 82 % cycling, G2/M peak wraps at +-pi, histone profile noisy |

Pilot six (unchanged, refit reproduced exactly; see `../pilot/PILOT.md`): GSE178341, GSE132465, GSE183904, SRP327758, GSE176078 valid; GSE131907 valid, narrow bands.

## Seed stability (`output/seed_stability.csv`)
Circular correlation (Jammalamadaka-SenGupta) and median |angular difference| (rad) between the kept seed and each other usable seed, over cells the kept seed calls cycling (p > 0.5), after the package orientation. Median difference is bounded by pi, a value near 0.2-0.3 with low correlation means the seeds agree on the G1 / cycling split but not on the order within the cycle.
- Seeds usable and agreeing: GSE181919 (0.85-0.91), GSE315534 (0.70-0.74), GSE178341 (0.93, then 0.00 and -0.03), GSE132465 (0.74, 0.88, 0.88, 0.02), SRP327758 (0.80-0.85), SRP416180 (0.54-0.86), GSE294300 (0.33-0.77).
- Unstable: GSE183904 (-0.41, -0.26, 0.33), GSE166555 (-0.82 to 0.43), GSE200997 (-0.47 to 0.58), GSE335070 (0.45, single comparison). Negative values mean a reversed or scrambled phase ordering in that seed.
- Extra seeds 5-9 (GSE131907, Wu GSE176078), compared with the kept seed (seed 1 / seed 4):
  - GSE131907: 3 of 5 usable (seeds 7, 8, 9); corr 0.82, 0.52, 0.24 (median diff 0.22-0.35). So 4 usable of 10 seeds; phase only moderately stable.
  - GSE176078: 3 of 5 usable (5, 7, 8); corr -0.47, 0.93, 0.91 (median diff 0.28, 0.14, 0.15). 4 usable of 10; two agree strongly, one is reversed.
- Pal and ERP149257 have a single usable seed (no comparison possible).
- The kept seed is not unique in several studies; phase from one seed should not be treated as definitive for the unstable studies.
