# method: GTEx benchmark with hidden times, and validations against true times

The method against cyclic-ordering methods on GTEx fast deaths with the time of death hidden, and the untimed host against known
sampling times (mouse colon cell types A31, baboon A36, living human skin A37). Estimates come from `clocksync` (`Dataset.fit`, `gate`,
`timed`, `gene_rhythms`, `paired.cross`, `waveform.fit_waveform`, `timed.agreement`); GTEx through `common/gtex.py` (aggregates only);
the comparison methods (COFE, CYCLOPS 2.0, CHIRAL, TimeTeller) are paper-only code in `common/benchmark/` (ported from
`legacy/src/sccc/benchmark/`, config imports changed). No plotting code. Outputs: `paper/results/figure_1/output/method/` (`.../subsample/` for
checks); the legacy runs are in `outputs/paper/fig1_method/`.

```
run.sh [--subsample N]       all scripts (full run: benchmark and CYCLOPS take hours)
_shared.py                   GTEx loader for the benchmark, local timed loading
compute_gtex_benchmark.py, compute_celltype_truth.py, compute_baboon_truth.py, compute_skin_truth.py,
compute_tcga_cofe.py, compute_cyclops_anchoring.py
```

## Analyses
| script | amendment | question and estimand | legacy script replaced |
|---|---|---|---|
| `compute_gtex_benchmark.py` | design of figure 1 (no amendment) | 8 GTEx tissues, draws of n = 10-80 fast deaths from the whole day or 08-18 h; tissue log2 amplitude and offset (reference free) of COFE (genome, clock genes), CYCLOPS, CHIRAL, TimeTeller (supervised), `clocksync` (waveform refitted without the draw) and the timed oracle, against the SVD of the recorded times of all fast deaths; donor-time error and skill of orderings; clocksync interval coverage and runtime; CCD | `paper_code/fig1_method/run_gtex_benchmark.py`, notebooks 2 and 6 |
| `compute_celltype_truth.py` | A31 (+ E1 window, post hoc) | per mouse colon cell type, untimed gene phases (`Dataset.gene_rhythms`) vs timed cosinor; pooled R of the 32 frozen followers vs host-permutation null; sensitivity arm `--arms C-DC,R-PC --tag sens_RPC`; `--window` | `validation/celltype_truth/scripts/{score,window}.py` |
| `compute_baboon_truth.py` | A36 | 12 animals x 63 tissues: host phase vs ZT (R, Fisher-Lee, Jammalamadaka-SenGupta), leave-one-animal-out, untimed tissue offsets vs timed gains, cross-species gain phases vs GTEx | `validation/baboon/scripts/{prep,score,human_gain}.py` |
| `compute_skin_truth.py` | A37 | dermis minus epidermis untimed vs timed offset (subject bootstrap), host phase vs clock time; 4 sensitivities | `validation/skin_paired/scripts/skin_paired.py` |
| `compute_tcga_cofe.py` | post hoc | COFE as published on 11 TCGA cancers (orderings, CV curve, validity flag) | `run_tcga_cofe.py` |
| `compute_cyclops_anchoring.py` | post hoc | can whole-day samples anchor CYCLOPS on hospital-hours samples (5 conditions, 3 tissues; breast case) | `run_cyclops_anchoring.py`, notebook 4 |

## Outputs
`gtex_benchmark/{tissues,samples,genes,coverage,runtime,ccd,summary,truth,truth_model,gene_truth}.csv` (per-run caches `runs/*.pkl`; delete to rerun),
`celltype_truth_{scores,genes,hosts}_<tag>.csv`, `celltype_truth_window_E1.csv`, `baboon_{readouts,animals,tissues,pairs,stability,gate}.csv`,
`skin_truth_{summary,boot,gate}.csv`, `tcga_cofe/`, `cyclops_anchoring/{scores,genes,breast_summary,breast_genes}.csv`. No per-donor
orderings or times are written (the legacy cached them locally). `--subsample N` means N donors per draw for the benchmark and CYCLOPS, N tumours for TCGA.

## Checks
| readout | here | legacy | note |
|---|---|---|---|
| benchmark truth (SVD on all fast deaths) | max difference 3e-7 | `gtex/truth.csv` | |
| clocksync n = 80 draw 0, full / window: tissue log2 amplitude, offset | 0.012, 2e-4 h / 0.005, 5e-4 h (max difference) | `runs/sccc_n80_*_0.csv` | same draws and cross-fitted waveform; 0 bootstraps |
| timed oracle, 2 draws | 3e-7 | | |
| CHIRAL n = 80, 16 tissue runs (2 draws) | 13 identical to 1e-5 h, 3 differ (Thyroid, Adipose) | `runs/chiral_*` | restarts are unstable to 1e-5 numerical differences |
| A31 pooled R of followers (`a1`), offset, n | 0.982, +0.15 h, 33; null 95 % 0.780, p 0.005 | 0.98, +0.2 h, 33; 0.82, 0.005 | all R/offset rows identical to 1e-12; nulls differ (within-arm permutation) |
| A31 sensitivity arm (C-DC + R-PC) | R 0.986, n 41, p 0.005 | 0.986, 41, 0.005 | |
| A31 E1 window, median pooled R followers | 3 ZTs 0.936, 2 ZTs 0.259 | 0.937, 0.271 | max difference 0.047 (100 permutations) |
| A36 readout 1: R, Fisher-Lee, median error | 0.987, 0.980 (p 0.002 at 500 permutations), 0.39 h | 0.987, 0.980, 0.39 h | JS r +0.676 here vs -0.874 (see Deviations) |
| A36 readout 2 (63 tissues): R, Fisher-Lee, median error | 0.972, 0.912, 0.46 h | 0.961, 0.859, 0.50 h | leave-one-out vs full-host offsets |
| A36 readout 3 (25 pairs): Fisher-Lee, p | 0.177, 0.198 | 0.177, 0.195 | |
| A37 untimed, timed, difference [95 % CI] | -1.69, -2.14, +0.44 [0.09, 0.75] h | -1.63, -2.16, +0.53 [0.16, 0.87] | FAIL on the CI clause in both |
| A37 host phase: JS rho, R, median error | 0.505, 0.706, 1.95 h | identical | |
| TimeTeller, CYCLOPS, CYCLOPS anchoring | run on `--subsample` (timeteller n = 24, 3 tissues; cyclops muscle, W = R = 30: skill window 0.19, pooled 0.12) | not compared | |

Not rerun: COFE genome-wide and CYCLOPS draws, TimeTeller at full n, benchmark bootstraps and coverage, CCD, TCGA COFE on the 11 cancers, the CYCLOPS
anchoring at n = 80 and the breast case.

## Deviations
- **Pooled null (A31).** Animals are permuted within arm (legacy: over all animals), so the null 95 % quantile is slightly lower (0.78 vs 0.82); p identical.
- **Timed truth (A31, A36, A37)** uses local least-squares fits (`timed_loading`, F test, gain cosinor), not the legacy weighted `factor.timed`.
- **A36 host phase.** `Dataset.timed` does not centre the host state; the legacy did. R, Fisher-Lee and error are unchanged, the (degenerate) JS correlation of the
  pre-registered primary changes sign. Permutation p at 500 (subsample) or 10,000. Readout 2 offsets come from the fit with a leave-one-out host (legacy: same,
  full-host differences in the amplitude-free offset). Human gains: 50 donor bootstraps, covariates as legacy.
- **A37 untimed offset** = phase of the paired cross-moment (`paired.cross`), which equals the free-host loading phase for two compartments; permutation of times is global (not within
  subject); host-phase intervals from 20-100 subject bootstraps.
- **Benchmark.** Duplicate samples keep the first (`common.gtex.load`); clocksync fits all tissues as anchors; no per-donor orderings are cached.

## Not ported
Notebooks 1, 3, 5, 7 of `fig1_method` (eigen-core clock, hospital hours, figure, estimand gap: scoring and plotting of the runs above), `ed_truth_validation.py` (plot; its panels
A40, A41 belong to `windows`, A48 to `species`), the TimeTeller environment setup.

## Local helpers (not in clocksync)
`timed_loading` (complex loading on clock time), `timed_cosinor` (F test), `timed_gain`, `human_gains`, `csvd_estimates`, `model_truth`, `consensus_times`, `cross_fit_waveform` (wraps `waveform.fit_waveform`).

## API gaps
`Dataset.timed` with centring and within-stratum permutation; `Dataset.timed_loadings`; `Dataset.fit` for two compartments (free host); a donor-bootstrap by subject; a public host-given coupling.
