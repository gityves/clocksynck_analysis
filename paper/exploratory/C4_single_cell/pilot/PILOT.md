# C4 pilot, method QC: CoPhaser cell-cycle phase of CNV-malignant cells (2026-11-01)

QC only. No statistic of clock-gene expression across phase was computed or stored; clock genes (14 core + BHLHE40/41, NFIL3, CIPC, BMAL1 alias) are
removed from the rhythmic set (96-98 packaged cell-cycle genes), the context set (2,000 HVG + canonical histones, 2,013-2,046 genes, n_latent 10) and the S/G2M score lists.
CoPhaser commit 93861a8 with local diff (sha256 4795ee58...b7bf16, in qc.csv), torch 2.9.1, scanpy 1.11.5, anndata 0.12.7. Script `scripts/fit_cellcycle.py` (one model per study, CNV-malignant cells,
<= 3,000 cells/patient, >= 1,000 UMI, `auto_hyperparameters(cycle="cell_cycle")`, `fit_with_restarts` 5 seeds x 200 epochs, GPU; kept = best-ranked run that is finite and not vetoed). 6-10 min per study.
Orientation: the package's `orient_align_pseudotimes` (G1/S at pi/4). S/G2M scores: scanpy Tirosh scoring on the full gene matrix.

| study | cells | pts | seeds kept/run | rank | cycling (p>0.5) | S peak -> G2M peak (rad) | marker ratio | lib-size ratio | order S->G2M |
|---|---|---|---|---|---|---|---|---|---|
| GSE178341 CRC | 64,231 | 60 | 4/5 | 0.65 | 0.495 | 0.74 -> 2.14 (gap 1.40) | 4.9 | 2.0 | holds |
| GSE132465 CRC | 16,674 | 23 | 5/5 | 0.55 | 0.463 | 0.64 -> 1.78 (1.14) | 5.0 | 1.6 | holds |
| GSE131907 LUAD | 8,625 | 15 | 1/5 | 0.80 | 0.301 | 0.91 -> 2.06 (1.15) | 5.0 | 2.6 | holds |
| GSE183904 gastric | 9,904 | 25 | 4/5 | 0.65 | 0.448 | 0.72 -> 2.09 (1.37) | 3.7 | 1.9 | holds |
| SRP327758 oesoph. | 11,552 | 43 | 4/5 | 0.65 | 0.412 | 0.67 -> 1.97 (1.30) | 4.0 | 2.4 | holds |
| GSE176078 Wu breast | 22,099 | 20 | 1/5 | 0.50 | 0.261 | 0.73 -> 2.29 (1.55) | 8.5 | 2.7 | holds |

All numbers: `output/qc.csv` (plus seeds_<study>.csv for every restart: vetoed runs (score < 0) were rejected; the veto reason per seed is not stored). Per patient cycling fraction: `output/patient_cycling.csv`.

## Reading of the validation figures (`output/validation_<study>.png`)
Common pattern in all six: PCNA and MCM6 (S) are unimodal and peak at about 0.1-0.35 pi; TOP2A and MKI67 (G2/M) are unimodal and peak later, 0.5-0.8 pi, with a clean trough (about -0.7 to -0.3 pi) where all four are near zero (G1);
the histone fraction has one sharp peak (~0.45 pi) between the S and G2/M markers; library size rises through S/G2 and drops after the histone peak (ratios 1.6-2.7, below the 3.0 veto); the S -> G2M order holds in forward direction.
The phase histogram has a dense G1 mode near -0.6 pi (mostly non-cycling cells, high-posterior-width cells pile there) and a flat low-density tail across the cycle, as expected when 50-70 % of tumour cells are non-cycling.
Phase uniformity (cycling cells) is 0.99-1.00. The modelled Fourier profiles (panel B) are smooth and agree with the data profiles (A). Valid in each case:
- **GSE178341**: cleanest (tight bands, 64k cells); MCM6 broader than PCNA with a late shoulder. Marker resultant is low (0.17) because most cells sit in G1; not a problem.
- **GSE132465**: valid; PCNA has a small second bump near +pi, MKI67 plateau; 5/5 seeds usable, scores nearly tied.
- **GSE131907 (LUAD)**: valid but only 1 of 5 seeds passed vetoes (the other four were vetoed) and bands are wider (8.6k cells, 15 patients, median 325 cells/patient); lowest cycling fraction (0.30).
- **GSE183904 (gastric)**: now valid. Old fit: 100 % cycling and S->G2M order failed. New: 45 % cycling, S peak 0.72, G2M 2.09, histone peak between them. 4/5 seeds usable.
- **SRP327758 (oesophageal)**: now valid. Old fit: 100 % cycling (branch off). New: 41 % cycling, ordered S -> histone -> G2M, 4/5 seeds usable. Caveat: many patients have very few cells (median 105, 25th pct 16; only 14 patients >= 300), so per-patient cycling fractions are noisy.
- **Wu breast (Cancer Epithelial)**: valid; sharpest marker contrast (ratio 8.5) but the cycling fraction is low (0.26; per-patient 0.00-0.90, ER-dependent, median 0.20), only 1/5 seeds kept, and the narrow-posterior cells concentrate at the G1 mode.
Wu used raw counts from the mtx directly (no extra loading problems).

## Caveats
- Orientation is the package's own (RPE alignment), not the post hoc G2M-peak orientation of the archive; the S->G2M check is an independent confirmation, not a selection criterion.
- Seeds differ in score by only 0.05-0.15 ranks in most studies; the kept seed is not uniquely best. For GSE131907 and Wu only one seed was usable, so phase stability across seeds should be checked (circular correlation of theta between kept seeds) before registration; not done here.
- Cycling-probability calls depend on the Bernoulli branch (prior 0.1 chosen automatically because 62-77 % of cells score G1).
- Malignant labels come from CNV calls; patient counts differ from the scout (it counted patients in the legacy cycle tables).
