"""A146 note 2: origin of the recorded-time pattern of part 1 (GSE233242 tumours).

Registered in paper/docs/registered_experiments.md (A146 note 2). H146.4 glucocorticoid-response score vs recorded
time; H146.5 omnibus S with programs residualised on age and subtype; D146.2 immediate-early and hypoxia scores;
D146.3 S under shuffled predicted hours, and S by program group. Outputs: output/checks_*.csv.
    python -I checks_part1.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

import score_part1 as P  # noqa: E402

GC = "TSC22D3 FKBP5 DDIT4 ZBTB16 KLF9".split()
IEG = "FOS FOSB JUN EGR1 ATF3 IER2".split()
HYPOXIA = "VEGFA SLC2A1 PGK1 CA9 ADM BNIP3".split()
IMMUNE = ["P4_Tcell_activation", "P5_myeloid_AP", "P6_phagocytosis", "P7_Tcell_abundance", "P8_myeloid_abundance",
          "P10_PD1"]
N_SHUFFLE = 10_000


def gene_score(x, genes, cols):
    g = [q for q in genes if q in x.index]
    return P.z(np.log2(x.loc[g, cols].T.to_numpy() + 1)).mean(1), g


def residualise(Y, X):
    A = np.c_[np.ones(len(Y)), X]
    return Y - A @ np.linalg.lstsq(A, Y, rcond=None)[0]


def slope_test(y, t, perms):
    """t of the slope of y on t and its one-sided permutation p in both directions."""
    obs = P.tstat(y[:, None], t[:, None])[0]
    null = P.tstat(np.repeat(y[:, None], 1, 1), t[perms][..., None])[:, 0]
    return obs, (1 + (null <= obs).sum()) / (len(null) + 1), (1 + (null >= obs).sum()) / (len(null) + 1), null.std(ddof=1)


def omni_test(Yk, t, hs, perms):
    obs = P.tstat(Yk, P.design(t, hs)).mean()
    null = P.tstat(Yk, P.design(t[perms], hs)).mean(1)
    return obs, (1 + (null >= obs).sum()) / (len(null) + 1), null


def main():
    rng = np.random.default_rng(P.SEED)
    tr = pd.read_csv(P.TRUTH)
    tr = tr[tr.note.fillna("") == ""].set_index("patient_id")
    v = pd.read_csv(P.V25).set_index("patient")
    d = tr.join(v[["est_day"]], how="inner")
    assert len(d) == 42
    pts, t = d.index.to_numpy(), d.resection_h.to_numpy()
    x = pd.read_csv(P.TPM).drop_duplicates("gene_name").set_index("gene_name")
    perms = np.array([rng.permutation(len(pts)) for _ in range(P.NPERM)])  # the part 1 stream

    Y = P.scores(x, d.tpm_column_tumour.to_numpy(), pts, [])
    cols = list(Y.columns)
    hs = np.array([P.H[k] for k in cols])
    Yk = Y[cols].to_numpy()
    S0, p0, _ = omni_test(Yk, t, hs, perms)

    rows, desc = [], []
    # H146.4: glucocorticoid score falls with recorded time
    for tis in ("tumour", "normal"):
        g, found = gene_score(x, GC, d[f"tpm_column_{tis}"].to_numpy())
        obs, p_lo, p_hi, sd = slope_test(g, t, perms)
        r = dict(check="H146.4" if tis == "tumour" else "H146.4 normal (descriptive)", tissue=tis, genes=" ".join(found),
                 t_slope=obs, p_one_sided_negative=p_lo, null_sd=sd, spearman=stats.spearmanr(g, t).statistic)
        (rows if tis == "tumour" else desc).append(r)

    # H146.5: programs residualised on age and subtype
    sub = pd.get_dummies(d.subtype_geo.astype(str), drop_first=True).to_numpy(float)
    cov = np.c_[d.age.to_numpy(float), sub]
    S5, p5, n5 = omni_test(residualise(Yk, cov), t, hs, perms)
    rows.append(dict(check="H146.5", tissue="tumour", S_unadjusted=S0, S_adjusted=S5, p_one_sided=p5,
                     null_sd=n5.std(ddof=1), covariates="age + subtype (" + ",".join(sorted(d.subtype_geo.astype(str).unique())) + ")"))
    desc.append(dict(check="D146 scheduling: age", spearman=stats.spearmanr(d.age, t).statistic,
                     p=stats.spearmanr(d.age, t).pvalue))
    groups = [t[d.subtype_geo.astype(str).to_numpy() == s] for s in sorted(d.subtype_geo.astype(str).unique())]
    kw = stats.kruskal(*[g for g in groups if len(g) > 1])
    desc.append(dict(check="D146 scheduling: subtype (Kruskal-Wallis)", statistic=kw.statistic, p=kw.pvalue,
                     group_sizes=";".join(f"{s}={len(g)}" for s, g in zip(sorted(d.subtype_geo.astype(str).unique()), groups))))
    for k in ("P7_Tcell_abundance", "P8_myeloid_abundance"):
        s = stats.spearmanr(Y[k], t)
        desc.append(dict(check=f"D146 scheduling: {k}", spearman=s.statistic, p=s.pvalue))

    # D146.2: handling
    hd = []
    for name, genes in (("immediate_early", IEG), ("hypoxia", HYPOXIA)):
        g, found = gene_score(x, genes, d.tpm_column_tumour.to_numpy())
        obs, p_lo, p_hi, sd = slope_test(g, t, perms)
        desc.append(dict(check=f"D146.2 {name}", genes=" ".join(found), t_slope=obs, p_one_sided_negative=p_lo,
                         p_one_sided_positive=p_hi, spearman=stats.spearmanr(g, t).statistic))
        hd.append(g)
    S2, p2, _ = omni_test(residualise(Yk, np.c_[hd[0], hd[1]]), t, hs, perms)
    desc.append(dict(check="D146.2 S residualised on IEG + hypoxia", S_unadjusted=S0, S_adjusted=S2, p=p2))

    # D146.3: shuffled predicted hours; program groups
    T = P.tstat(Yk, P.design(t, hs))  # observed per-program t under the registered hours
    rs = np.random.default_rng([P.SEED, 3])
    shuf = np.array([P.tstat(Yk, P.design(t, hs[rs.permutation(len(hs))])).mean() for _ in range(N_SHUFFLE)])
    desc.append(dict(check="D146.3 shuffled hours", S_registered=S0, shuffled_mean=shuf.mean(), shuffled_q95=np.quantile(shuf, .95),
                     share_shuffled_ge_obs=float((shuf >= S0).mean()), n_shuffle=N_SHUFFLE))
    for name, grp in (("immune programs", [c for c in cols if c in IMMUNE]), ("other programs", [c for c in cols if c not in IMMUNE])):
        i = [cols.index(c) for c in grp]
        So, ps, _ = omni_test(Yk[:, i], t, hs[i], perms)
        desc.append(dict(check=f"D146.3 {name}", programs=" ".join(grp), S=So, p=ps,
                         per_program_t=";".join(f"{c}={T[j]:.2f}" for c, j in zip(grp, i))))

    R = pd.DataFrame(rows)
    p = np.array([R.loc[R.check == "H146.4", "p_one_sided_negative"].iloc[0], R.loc[R.check == "H146.5", "p_one_sided"].iloc[0]])
    R["p_holm"] = P.holm(p)
    ok = [R.t_slope.iloc[0] < 0, R.S_adjusted.iloc[1] > 0]
    R["verdict"] = np.where((R.p_holm < 0.05) & np.array(ok), "pass", "fail (not detected at this design)")
    out = P.OUT
    R.to_csv(out / "checks_tests.csv", index=False)
    pd.DataFrame(desc).to_csv(out / "checks_descriptive.csv", index=False)
    pd.set_option("display.width", 250, "display.max_colwidth", 120)
    print(R.round(4).T.to_string())
    print(pd.DataFrame(desc).round(4).T.to_string())


if __name__ == "__main__":
    main()
