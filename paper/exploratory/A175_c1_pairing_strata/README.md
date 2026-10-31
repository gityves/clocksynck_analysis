# A175 (C1): pairing of tumour to own normal tissue by proliferation and ER stratum

Registered in `paper/docs/registered_experiments.md` (A175). Secondary tests on the H138.4 set (517 pairs, 10 gated cancers).

    RUN_CAP_GB=20 paper/common/run.sh paper/exploratory/A175_c1_pairing_strata/scripts/run_a175.py --stage power --jobs 4
    RUN_CAP_GB=20 paper/common/run.sh paper/exploratory/A175_c1_pairing_strata/scripts/run_a175.py --stage tests --jobs 4

- `power`: permutation null of D (labels only, observed D never evaluated) for H175.1 and H175.2 -> `output/power.csv`, `power_arms.csv`.
- `tests`: stop rule on H175.2 (MDE > 0.3 -> not run), observed D, p, Holm, pseudo-clock and responsive companions, descriptives
  -> `output/tests.csv`, `per_cancer.csv`, `descriptives.csv`, `h138_4_reproduction.csv`.

Library code is imported from `../A138_cancer_estimability/scripts` (`q3_lib`, `h138_8.Cancer`). The script first reproduces H138.4 = 0.3647.
