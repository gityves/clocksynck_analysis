# A126 menstrual-cycle natural experiment (registered in paper/docs/registered_experiments.md)
Order of scripts (python -I, sccc env; Kumar extraction in census_env via run_kumar.sh through paper/common/run.sh):
1. `h126_0_endometrium.py` (GSE266265 per-cell export, A24 quartile rule) -> `h126_0_test.py` (H126.0 + sensitivities + phase descriptive)
2. `make_panel.py` -> `run_kumar.sh` (Census, authors' broad_cell_type) ; `parse_murrow.py` (steps 1, 2; Murrow rds, Python only) ; `build_pb.py`
3. `score_validity.py` (luteal score, validity, feasibility; written before any clock statistic)
4. `run_h126.py` (H126.1, H126.3, pseudo-clock companion, Holm, MDE, descriptives)
Outputs in `output/` (see run_h126.log). Not run: breast cycling-vs-non-cycling descriptive, core-cue offset vs luteal score, GTEx breast.
