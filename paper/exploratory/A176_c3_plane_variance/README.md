# A176: swamped or disordered clock? (C3)

Registered in `paper/docs/registered_experiments.md` (A176). For the seven registered coherence contrasts (SCAN-B H127.1/H127.2,
METABRIC H166.1, GSE72094, GSE68465, GSE62254, GSE15459) the script computes, low-coherence arm minus reference arm and averaged over the
registered strata, dR, dlog_plane, dlog_off (`clocksync.decompose.components`) and ds = dlog_plane - dlog_off.

    scripts/run_a176.py --stage calibration   # reproduces dR (assert 1e-6), nulls, pseudo-set regression, MDE; no observed ds
    scripts/run_a176.py --stage tests         # observed statistics, one-sided p, Holm over seven, H176.2 calls, descriptives

Run through `RUN_CAP_GB=20 paper/common/run.sh`. Output (git-ignored) in `output/`: `reproduction.csv`, `calibration.csv`,
`null_draws.csv`, `pseudo_stats.csv` (registered 200 pseudo sets per cohort, observed labels), `tests.csv`, `verdict.csv`,
`descriptive_without_per1_nr1d1.csv`, `descriptive_prolif_pool.csv` (+ pool sets / stats). One permutation draw gives all four statistics.
