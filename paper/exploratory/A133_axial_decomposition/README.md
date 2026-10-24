# A133: circular / axial decomposition of the host covariance

Registered in `paper/docs/registered_experiments.md` (A133). This folder holds only the simulation gate V133.1.

- `scripts/axial.py`: estimator (rank-1 Takagi factor of the pseudo-covariance, iterative imputation for missing pairs,
  `C_circ = C - w w^H`); `tests/test_axial.py`: exact recovery in the noiseless case.
- `scripts/compute_v133_sims.py`: generative model (A109 `simulate_group` with host phase law and per-tissue directions),
  grid f x n x host R x direction, 100 reps per cell; gate criteria.
- `run.sh [--subsample 1] [--seed S]`: capped run; outputs `output/v133_grid.csv`, `output/v133_gate.csv`.

Choices beyond the registration: c = 0.45 (A109 C_target); amplitudes and Poisson variance from the A109 ventilated
clock-plane table; true phases fitted to the A109 fast-death pair offsets (SD 1.3 h); f = 0.5 cell added to the grid;
pair observed when >= 20 shared donors; supplementary presence 0.8 rows; gate values are the worst cell of the stated set.
H133.1 is not run.
