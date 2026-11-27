"""A146 addendum 1 (POST HOC, DESCRIPTIVE; registered in paper/docs/registered_experiments.md): is the internal-time null of
A146 part 1 explained by the error of the internal-time estimate?

  E1  noise-injection null: est_day = a + b x recorded + e; t_sim = a + b x recorded + permuted e (and, as sensitivity,
      recorded + permuted (est_day - recorded)); S(t_sim) on the tumour programs and on the normal-tissue omnibus (P2-P10).
  E2  Spearman of e with each tumour program score, stromal and immune content and subtype (Kruskal-Wallis).
                                                                     -> output/addendum1_{e1,e2}.csv
    python -I addendum1_noise.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

import score_part1 as P  # noqa: E402

SEED, N = 1461, 10_000
CONTENT = {"stromal": "COL1A1 COL1A2 DCN LUM", "immune": "PTPRC"}


def main():
    rng = np.random.default_rng(SEED)
    tr = pd.read_csv(P.TRUTH)
    tr = tr[tr.note.fillna("") == ""].set_index("patient_id")
    d = tr.join(pd.read_csv(P.V25).set_index("patient")[["est_day"]], how="inner")
    assert len(d) == 42
    pts = d.index.to_numpy()
    rec, est = d.resection_h.to_numpy(float), d.est_day.to_numpy(float)
    x = pd.read_csv(P.TPM).drop_duplicates("gene_name").set_index("gene_name")
    Y = {tis: P.scores(x, d[f"tpm_column_{tis}"].to_numpy(), pts, []) for tis in ("tumour", "normal")}
    cols = {"tumour": list(Y["tumour"].columns), "normal": [k for k in Y["normal"].columns if k != "P1_proliferation"]}
    obs = {tis: (P.omnibus(Y[tis][cols[tis]].to_numpy(), est, [P.H[k] for k in cols[tis]]),
                 P.omnibus(Y[tis][cols[tis]].to_numpy(), rec, [P.H[k] for k in cols[tis]])) for tis in cols}
    assert abs(obs["tumour"][0] - (-0.0560)) < 5e-4 and abs(obs["tumour"][1] - 1.3393) < 5e-4, obs  # A146 part 1 reproduced

    A = np.c_[np.ones(len(rec)), rec]
    coef = np.linalg.lstsq(A, est, rcond=None)[0]
    e = est - A @ coef
    perms = np.array([rng.permutation(len(rec)) for _ in range(N)])
    sims = {"calibrated (a + b x recorded + e)": A @ coef + e[perms], "uncompressed (recorded + error)": rec + (est - rec)[perms]}
    rows = []
    for tis in cols:
        Yk, hs = Y[tis][cols[tis]].to_numpy(), np.array([P.H[k] for k in cols[tis]])
        for lab, T in sims.items():
            S = P.tstat(Yk, P.design(T, hs)).mean(1)
            o = obs[tis][0]
            rows.append(dict(tissue=tis, null=lab, a=coef[0], b=coef[1], resid_sd=e.std(ddof=2), S_internal=o, S_recorded=obs[tis][1],
                             sim_mean=S.mean(), sim_q05=np.quantile(S, .05), sim_q95=np.quantile(S, .95),
                             percentile=100 * np.mean(S < o), p_one_sided_low=(1 + np.sum(S <= o)) / (N + 1)))
    E1 = pd.DataFrame(rows)
    E1.to_csv(P.OUT / "addendum1_e1.csv", index=False)

    rows = []
    feats = dict(Y["tumour"].items())
    for nm, g in CONTENT.items():
        found = [q for q in g.split() if q in x.index]
        feats[nm] = P.z(np.log2(x.loc[found, d.tpm_column_tumour].T.to_numpy() + 1)).mean(1)
    for nm, v in feats.items():
        v = np.asarray(v, float)
        rho = stats.spearmanr(e, v).statistic
        null = np.array([stats.spearmanr(e[p], v).statistic for p in perms[:N]])
        rows.append(dict(feature=nm, spearman=rho, p_perm=(1 + np.sum(np.abs(null) >= abs(rho))) / (N + 1)))
    sub = d.subtype_geo.astype(str).to_numpy()
    kw = stats.kruskal(*[e[sub == s] for s in sorted(set(sub))])
    rows.append(dict(feature="subtype (Kruskal-Wallis H)", spearman=kw.statistic, p_perm=kw.pvalue))
    E2 = pd.DataFrame(rows)
    E2.to_csv(P.OUT / "addendum1_e2.csv", index=False)
    pd.set_option("display.width", 250)
    print(E1.round(4).to_string(index=False))
    print(E2.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
