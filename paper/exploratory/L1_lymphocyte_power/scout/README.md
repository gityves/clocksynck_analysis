# L1 scout: T-cell vs fibroblast timed clock amplitude, mouse colon chrono-atlas (not registered, no conclusions)

Script `amplitude_ratio_colon.py`; outputs in `output/` (`amplitude_ratios.csv` = per arm x compartment x gene set, amplitudes, depth and ratio vs fibroblast; `amplitudes.csv`, `ratios_vs_fibroblast.csv` the two halves).

- Data: `common.data.load("colon_chrono")` (legacy conversion; A273's table), arms C-DC, R-DC, R-PC; pseudobulks with n_cells >= 30 (A273 rule), >= 6 mice.
- Expression: log2(1 + count / depth * 1e6), depth = all-gene UMIs of the pseudobulk. The table holds only 14 clock genes: ARNTL CLOCK NPAS2 PER1-3 CRY1-2 NR1D1-2 DBP RORA TEF CIART (no HLF, BHLHE40/41, NFIL3). Sets: `all14`; `core13` = all14 minus RORA. NPAS2 is ~0 everywhere.
- Amplitude A: per gene OLS of expression on 1, cos(2 pi ZT/24), sin(...) over mice; squared amplitude a^2+b^2 minus its null expectation s2*(Cov_cc+Cov_ss) (s2 = residual variance, n-3 df; for balanced ZT this is 4 s2/n, i.e. both coefficients, not 2 s2/n); summed over genes, floored at 0, A = sqrt. Units: log2-CPM.
- Bootstrap: 2000 draws, seed 1, mice resampled within ZT (2 per ZT per arm, so each ZT keeps 2 draws with replacement); the debiasing is re-estimated per draw. Caveat: duplicated mice shrink the bootstrap residual variance, so the bootstrap distribution sits above the point estimate (boot median > point in the table); the upper bound is therefore if anything conservative. A (bootstrap) ratio uses A_fib floored at 1e-6.
- `clocksync_abs_amplitude` / `clocksync_coupling` = `Dataset.timed_loadings("zt")` (|a| in clock-plane units, coupling to known time), different units from A.
