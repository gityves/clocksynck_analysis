"""Cyclic-ordering methods compared with `clocksync` in figure 1 (paper-only, not part of the method package).

Ported unchanged from `legacy/src/sccc/benchmark/` (only the config imports changed): COFE (`cofe`), CYCLOPS 2.0 (`cyclops`,
Julia subprocess in `$SCCC_DATA/tools`), CHIRAL through pychiral (`chiral`), TimeTeller (`timeteller`, R subprocess), the scores
(`evaluate`) and bulk-expression helpers (`data`: GTEx fast deaths, TCGA recount3). They never touch `clocksync` internals; the
host-referenced and the oracle estimates of the benchmark are computed in `method/_shared.py` with the public API.
"""
