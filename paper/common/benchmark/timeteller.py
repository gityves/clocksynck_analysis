"""TimeTeller (Vlachou et al. 2024, PLoS Comput Biol 20:e1011779), a supervised clock-time and clock-dysfunction
estimator, run in its own R environment (`$SCCC_DATA/tools/timeteller_env`, R 4.3, built by `setup_timeteller.sh`).

TimeTeller models the clock-gene vector at each training time point (local PCA, multivariate normal, periodic
interpolation over the day) and predicts a sample's time as the maximum of its likelihood curve; Theta, the fraction
of the day where the likelihood ratio exceeds a cosine envelope, is its dysfunction score. It needs discrete,
replicated training times: recorded times are binned (`BIN_H`), and bins with fewer than `MIN_PER_BIN` samples are
left out of training. With 2 h bins some GTEx bins are sparse (few night-time deaths) and the per-bin Gaussians
unstable; 3 h bins (8 time points) gave the lowest error on a held-out split of skeletal muscle (2.6 h vs 4.3 h for
2 h bins; thyroid 2.5-2.7 h for any width). `cross_fit` gives every sample a prediction from a model not trained on it.
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd

from common.paths import DATA

RSCRIPT = DATA / "tools" / "timeteller_env" / "bin" / "Rscript"
DRIVER = Path(__file__).with_name("timeteller_run.R")
BIN_H, MIN_PER_BIN = 3.0, 5


def run_timeteller(train: pd.DataFrame, train_times: pd.Series, test: pd.DataFrame, out_dir: Path,
                   cores: int = 1, bin_h: float = BIN_H) -> pd.DataFrame:
    """`train`, `test`: samples x genes (log expression); `train_times`: h. Returns time_pred (h), theta, max_loglik
    per test sample."""
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    t = (np.floor(train_times.loc[train.index].to_numpy() % 24 / bin_h) + 0.5) * bin_h
    keep = pd.Series(t).map(pd.Series(t).value_counts()).to_numpy() >= MIN_PER_BIN
    tr = train.loc[keep].assign(time=t[keep])
    tr.rename_axis("sample").reset_index().to_csv(out_dir / "train.csv", index=False)
    test.rename_axis("sample").reset_index().to_csv(out_dir / "test.csv", index=False)
    env = {**os.environ, "R_ENVIRON_USER": os.devnull, "R_LIBS_USER": str(RSCRIPT.parents[1] / "lib" / "R" / "library"),
           "OMP_NUM_THREADS": "1", "OPENBLAS_NUM_THREADS": "1"}  # the environment's library only, not ~/R
    with open(out_dir / "timeteller.log", "w") as log:
        subprocess.run([str(RSCRIPT), str(DRIVER), str(out_dir / "train.csv"), str(out_dir / "test.csv"),
                        str(out_dir / "pred.csv"), str(cores)], env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
    res = pd.read_csv(out_dir / "pred.csv", index_col="sample")
    res.index = res.index.astype(test.index.dtype)
    return res.loc[test.index]


def cross_fit(logexpr: pd.DataFrame, times: pd.Series, out_dir: Path, k: int = 5, seed: int = 0,
              cores: int = 1) -> pd.DataFrame:
    """Out-of-fold TimeTeller predictions for every sample (k folds, random)."""
    fold = np.random.default_rng(seed).permutation(len(logexpr)) % k
    out = [run_timeteller(logexpr[fold != f], times, logexpr[fold == f], Path(out_dir) / f"fold{f}", cores)
           for f in range(k)]
    return pd.concat(out).loc[logexpr.index].assign(fold=fold)
