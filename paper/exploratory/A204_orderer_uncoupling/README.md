# A204: orderer outputs under working hours (COFE labels, tumour-normal differences; TimeTeller at REMAGUS size)

Registered in `paper/docs/registered_experiments.md` (A204, commit 70e3564).
- `scripts/compute_a204.py`: simulated fits, cached per replicate in `output/runs/` (run through `paper/common/run.sh`).
- `scripts/score_a204.py`: simulator gate, H204.1-4 (Holm), descriptives -> `output/{gate,tests,cofe_cells,auc,tt_cells,*_per_rep}.csv`.
- `paper/common/benchmark/cofe.py` gained `test=`: COFE's test rule (training means/SDs, projection on the loadings);
  checked on training samples used as test: median label difference 0.013 h (max 0.10 h).
