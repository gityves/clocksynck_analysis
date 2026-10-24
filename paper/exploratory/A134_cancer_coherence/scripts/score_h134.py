"""A134 scoring: H134.1 (inferred copy number) and H134.2 (MYC, HYPOXIA, cycling; one joint regression) as registered in the
A134 note (registered_experiments.md). Statistic: coherence gain (adjusted minus raw, A87 coherence) of the malignant
epithelium minus the gain of the normal epithelium given the same regression on its own covariates; per gene OLS on the
stratum-centred covariates, A87 `components` on the residuals; patient bootstrap (1,000) interval; p from 1,000 permutations of
whole-patient covariate rows across patients (malignant and normal permuted independently), one-sided > 0.

    python score_h134.py [--reps 1000]
Outputs: output/a134_h1_h2.csv (statistics), a134_h1_h2_gene_effects.csv (per-gene slope on the malignant covariates).
"""
from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

PAPER = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(PAPER / "results" / "figure_5" / "scripts" / "cancer"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from clocksync.config import limit_threads  # noqa: E402
from clocksync.decompose import components, logfrac  # noqa: E402
from common import data  # noqa: E402

import _shared as sh  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "output"
GENES, WF, MAL, NORM = sh.GENES, sh.WF, sh.MAL, sh.NORM
PROG = ["myc", "hypoxia", "cycling"]


def centre(x, g):
    return x - pd.DataFrame(x).groupby(g).transform("mean").to_numpy()


def coh(y, pv, k):
    return components(y, pv, k, WF)["coherence"]


def residual(y, Z, g):
    """Per gene OLS residuals of y (n x G) on Z (n x G x p, stratum-centred inside)."""
    r = y.copy()
    for j in range(y.shape[1]):
        z = centre(Z[:, j, :], g)
        b = np.linalg.lstsq(z, y[:, j], rcond=None)[0]
        r[:, j] = y[:, j] - z @ b
    return r


def gain(y, pv, k, Z, g):
    raw = coh(y, pv, k)
    return coh(residual(y, Z, g), pv, k) - raw, raw


def stat(case, ref, Zm, Zn, idx, pm=None, pn=None):
    n = len(idx)
    out = []
    for d, Z, p in ((case, Zm, pm), (ref, Zn, pn)):
        x = d.iloc[idx]
        y, pv, k = logfrac(x, GENES, x.stratum.to_numpy())
        Zi = Z[idx]
        if p is not None:
            Zi = Zi[p]
        out.append(gain(y, pv, k, Zi, x.stratum.to_numpy()))
    return out[0][0] - out[1][0], out[0][0], out[1][0], out[0][1], out[1][1]


def run(name, case, ref, Zm, Zn, reps, seed):
    rng = np.random.default_rng(seed)
    n = len(case)
    obs = stat(case, ref, Zm, Zn, np.arange(n))
    boot = np.array([stat(case, ref, Zm, Zn, rng.integers(0, n, n))[0] for _ in range(reps)])
    perm = []
    for _ in range(reps):
        perm.append(stat(case, ref, Zm, Zn, np.arange(n), rng.permutation(n), rng.permutation(n))[0])
    perm = np.array(perm)
    return dict(test=name, patients=n, reps=reps, stat=obs[0], boot_lo=np.quantile(boot, .025), boot_hi=np.quantile(boot, .975),
                boot_sd=boot.std(ddof=1), p_perm_one_sided=(1 + (perm >= obs[0]).sum()) / (1 + reps),
                null_mean=perm.mean(), null_sd=perm.std(ddof=1), mde=2.49 * perm.std(ddof=1),
                gain_malignant=obs[1], gain_normal=obs[2], raw_coh_malignant=obs[3], raw_coh_normal=obs[4],
                adj_coh_malignant=obs[3] + obs[1], adj_coh_normal=obs[4] + obs[2])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=1000)
    a = ap.parse_args()
    limit_threads(2)
    warnings.filterwarnings("ignore")
    pts = pd.read_csv(OUT / "scout_patients79.csv")
    c = pd.read_csv(OUT / "scout_patient_covariates.csv", dtype={"patient": str})
    c["donor"] = c.study + "|" + c.patient
    pb, _ = data.load("cancer_sc")
    pb = pb[pb.donor.isin(pts.donor)]
    mal = pb[pb.compartment == MAL].drop_duplicates("donor").set_index("donor")
    nor = pb[pb.compartment == NORM].drop_duplicates("donor").set_index("donor")
    ok = sorted(set(mal.index) & set(nor.index))
    assert len(ok) == 79
    case, ref = mal.loc[ok], nor.loc[ok]
    cm = c[c.group == "malignant"].set_index("donor").loc[ok]
    cn_ = c[c.group == "normal"].set_index("donor").loc[ok]
    cols = [f"cn_{g}_mean" for g in GENES]
    Zm1 = cm[cols].to_numpy(float)[:, :, None]
    Zn1 = cn_[cols].to_numpy(float)[:, :, None]
    P = lambda d: np.repeat(d[PROG].to_numpy(float)[:, None, :], len(GENES), 1)  # noqa: E731
    cm = cm.rename(columns={"cycling": "cycling"})
    rows = [run("H134.1_copy_number", case, ref, Zm1, Zn1, a.reps, 1341),
            run("H134.2_programmes_joint", case, ref, P(cm), P(cn_), a.reps, 1342)]
    # secondary descriptives (not registered as tests): each programme alone
    for p in PROG:
        f = lambda d, p=p: np.repeat(d[[p]].to_numpy(float)[:, None, :], len(GENES), 1)  # noqa: E731
        rows.append(run(f"H134.2_single_{p}", case, ref, f(cm), f(cn_), a.reps, 1342))
    r = pd.DataFrame(rows)
    r.to_csv(OUT / "a134_h1_h2.csv", index=False)
    print(r.round(4).T.to_string())
    # per gene slopes (malignant), standardised
    y, pv, k = logfrac(case, GENES, case.stratum.to_numpy())
    g = case.stratum.to_numpy()
    ge = []
    for j, gene in enumerate(GENES):
        z = centre(Zm1[:, j, :], g)[:, 0]
        ge.append(dict(gene=gene, r_cn=np.corrcoef(z, y[:, j])[0, 1], r2_cn=np.corrcoef(z, y[:, j])[0, 1] ** 2))
    pd.DataFrame(ge).to_csv(OUT / "a134_h1_h2_gene_effects.csv", index=False)


if __name__ == "__main__":
    main()
