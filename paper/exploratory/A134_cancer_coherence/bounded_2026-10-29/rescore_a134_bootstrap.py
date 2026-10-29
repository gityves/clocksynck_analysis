"""A134 H134.1: what the permutation null permutes, and a bootstrap-based p. DESCRIPTIVE / POST HOC (2026-10-29, roadmap 5.4):
the registered verdict (p_perm one-sided, weak) is not changed.

Registered null (`score_h134.run`): whole-patient covariate rows (copy-number state of the clock-gene loci; one row per
patient, all genes together) are permuted across patients, malignant and normal independently, with expression and
coherence estimator fixed, so the statistic (gain of malignant minus gain of normal coherence after regressing on the
covariates) is recomputed on 79 intact patients. Its spread therefore contains only the chance fit of unrelated covariates,
not the sampling variability of the coherence estimator over patients that the patient bootstrap resamples.
Here, with the registered seeds and the same code (`score_h134.stat`), per replicate: the bootstrap statistic and its two
gains, the registered permutation null (reproduced), and a patient-resampled permutation null (resample patients, then
permute covariate rows: both sources of variation).
Bootstrap p's: share of replicates <= 0; bootstrap-t (stat / boot SD, normal); centred-bootstrap (share of boot - stat >=
stat); p from the patient-resampled permutation null.

    python rescore_a134_bootstrap.py [--reps 1000]
Outputs: output/{h134_1_bootstrap,h134_1_replicates}.csv (+ same for H134.2 joint, descriptive)
"""
from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[0] / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import norm  # noqa: E402

import score_h134 as H  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import data  # noqa: E402

OUT = HERE / "output"
PROG = H.PROG


def run(name, case, ref, Zm, Zn, reps, seed):
    rng = np.random.default_rng(seed)  # same stream order as score_h134.run: bootstrap, then permutation
    n = len(case)
    full = np.arange(n)
    obs = H.stat(case, ref, Zm, Zn, full)
    boot = np.array([H.stat(case, ref, Zm, Zn, rng.integers(0, n, n))[:3] for _ in range(reps)])
    perm = np.array([H.stat(case, ref, Zm, Zn, full, rng.permutation(n), rng.permutation(n))[:3] for _ in range(reps)])
    rng2 = np.random.default_rng([seed, 99])
    rp = np.array([H.stat(case, ref, Zm, Zn, (i := rng2.integers(0, n, n)), rng2.permutation(n), rng2.permutation(n))[:3]
                   for _ in range(reps)])
    s = obs[0]
    p_perm = (1 + (perm[:, 0] >= s).sum()) / (1 + reps)
    row = dict(test=name, patients=n, reps=reps, stat=s, p_perm_registered_check=p_perm,
               perm_null_sd=perm[:, 0].std(ddof=1), perm_null_mean=perm[:, 0].mean(),
               boot_sd=boot[:, 0].std(ddof=1), boot_mean=boot[:, 0].mean(), boot_lo=np.quantile(boot[:, 0], .025), boot_hi=np.quantile(boot[:, 0], .975),
               sd_ratio_boot_over_perm=boot[:, 0].std(ddof=1) / perm[:, 0].std(ddof=1),
               p_boot_le0=(1 + (boot[:, 0] <= 0).sum()) / (1 + reps),
               boot_t=s / boot[:, 0].std(ddof=1), p_boot_t_normal=1 - norm.cdf(s / boot[:, 0].std(ddof=1)),
               p_boot_centred=(1 + ((boot[:, 0] - s) >= s).sum()) / (1 + reps),
               resampled_perm_sd=rp[:, 0].std(ddof=1), resampled_perm_mean=rp[:, 0].mean(),
               p_resampled_perm=(1 + (rp[:, 0] >= s).sum()) / (1 + reps),
               gain_mal=obs[1], gain_norm=obs[2],
               boot_sd_gain_mal=boot[:, 1].std(ddof=1), boot_sd_gain_norm=boot[:, 2].std(ddof=1),
               boot_corr_gains=np.corrcoef(boot[:, 1], boot[:, 2])[0, 1],
               perm_sd_gain_mal=perm[:, 1].std(ddof=1), perm_sd_gain_norm=perm[:, 2].std(ddof=1),
               resampled_perm_sd_gain_mal=rp[:, 1].std(ddof=1), resampled_perm_sd_gain_norm=rp[:, 2].std(ddof=1))
    rep = pd.DataFrame(dict(boot=boot[:, 0], perm=perm[:, 0], resampled_perm=rp[:, 0]))
    return row, rep


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=1000)
    a = ap.parse_args()
    limit_threads(1)
    warnings.filterwarnings("ignore")
    O = H.OUT
    pts = pd.read_csv(O / "scout_patients79.csv")
    c = pd.read_csv(O / "scout_patient_covariates.csv", dtype={"patient": str})
    c["donor"] = c.study + "|" + c.patient
    pb, _ = data.load("cancer_sc")
    pb = pb[pb.donor.isin(pts.donor)]
    mal = pb[pb.compartment == H.MAL].drop_duplicates("donor").set_index("donor")
    nor = pb[pb.compartment == H.NORM].drop_duplicates("donor").set_index("donor")
    ok = sorted(set(mal.index) & set(nor.index))
    assert len(ok) == 79
    case, ref = mal.loc[ok], nor.loc[ok]
    cm = c[c.group == "malignant"].set_index("donor").loc[ok]
    cn_ = c[c.group == "normal"].set_index("donor").loc[ok]
    cols = [f"cn_{g}_mean" for g in H.GENES]
    Zm1, Zn1 = cm[cols].to_numpy(float)[:, :, None], cn_[cols].to_numpy(float)[:, :, None]
    row, rep = run("H134.1_copy_number", case, ref, Zm1, Zn1, a.reps, 1341)
    pd.DataFrame([row]).to_csv(OUT / "h134_1_bootstrap.csv", index=False)
    rep.to_csv(OUT / "h134_1_replicates.csv", index=False)
    print(pd.Series(row).to_string(), flush=True)


if __name__ == "__main__":
    main()
