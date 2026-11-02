"""A196 scoring (registered_experiments.md, A196, commit a0dc0c5): does arresting proliferation restore tumour clock coherence?
H196.1: NeoPalAna (GSE93204) T = coh(C1D15) - coh(C1D1), anastrozole + palbociclib vs anastrozole alone, within patient.
H196.2 (post hoc): reversal calibration; observed T against the rise P predicted by each series' own cross-sectional
proliferation-coherence relation (POETIC AI, GSE20181, GSE87411, GSE93204), pooled by inverse variance.
Coherence, designs, pseudo sets and scores are those of A168 (`a168_lib`).

Order: treatment check and power (--power-only stops there) -> H196.1, specificity, descriptives -> H196.2.
Run: python -I score_a196.py [--nperm 2000] [--npseudo 200] [--power-only]
"""
import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[2] / "A168_endocrine" / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

import a168_lib as L  # noqa: E402

OUT = HERE.parents[1] / "output"
A168_OUT = L.OUT
K, LIMIT, SEED = 2.49, 0.40, 20261103
say = lambda *a: print(pd.Timestamp.now().strftime("%F %T"), *a, flush=True)


def series(key):
    m, X = L.load(key, OUT if key == "neopal" else A168_OUT)
    cg = L.clock_genes(X)
    return m, X, L.Yof(X, cg), L.wf_of(cg)


def p2(null, t):
    c = null.mean()
    return float((1 + np.sum(np.abs(null - c) >= abs(t - c))) / (1 + len(null)))


def pairs(m, pre, post, arms=None):
    a, b, _ = L.pairs_from_meta(m, pre, post, arms)
    return L.Pairs("swap", "swap", a, b)


def calibration(name, m, Y, wf, d, nperm, rng):
    """Observed T and the rise P predicted from the pre-treatment proliferation tertiles (A196 H196.2)."""
    T = d.T(Y, wf); nT = d.null(Y, wf, nperm, rng)
    pre = d.pre
    pr = m.prolif.to_numpy()
    ter = L.A166.tertile(pr[pre])
    num = lambda t: L.coh(Y, pre[t == 0], wf) - L.coh(Y, pre[t == 2], wf)
    N = num(ter); nN = np.array([num(rng.permutation(ter)) for _ in range(nperm)])
    den = pr[pre][ter == 2].mean() - pr[pre][ter == 0].mean()
    dpro = float((pr[d.post] - pr[d.pre]).mean())
    scale = -dpro / den
    Nc, Tc = N - nN.mean(), T - nT.mean()
    P, sdP, sdT = Nc * scale, nN.std(ddof=1) * abs(scale), nT.std(ddof=1)
    return dict(series=name, n_pairs=d.n, T=T, T_centred=Tc, sd_T=sdT, numerator_low_minus_high=N, numerator_centred=Nc,
                sd_numerator=nN.std(ddof=1), p_numerator_one_sided=float((1 + np.sum(nN >= N)) / (1 + nperm)), prolif_high_minus_low=den,
                d_prolif=dpro, P=P, sd_P=sdP, T_minus_P=Tc - P, z=(Tc - P) / np.hypot(sdT, sdP))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--nperm", type=int, default=2000); ap.add_argument("--npseudo", type=int, default=200)
    ap.add_argument("--power-only", action="store_true")
    a = ap.parse_args()
    m, X, Y, wf = series("neopal")
    d = pairs(m, "C1D1", "C1D15")
    pr = m.prolif.to_numpy()
    dp = pr[d.post] - pr[d.pre]
    check = dict(n_pairs=d.n, d_prolif_mean=float(dp.mean()), d_prolif_sd=float(dp.std(ddof=1)),
                 wilcoxon_p=float(stats.wilcoxon(dp).pvalue), acted=bool(dp.mean() < 0 and stats.wilcoxon(dp).pvalue < .05))
    rng = np.random.default_rng(SEED)
    null = d.null(Y, wf, a.nperm, rng)
    power = dict(test="H196.1", null_mean=float(null.mean()), null_sd=float(null.std(ddof=1)), mde=K * float(null.std(ddof=1)), limit=LIMIT,
                 estimable=K * float(null.std(ddof=1)) <= LIMIT, n_perm=a.nperm)
    pd.DataFrame([{**check, **power}]).to_csv(OUT / "power.csv", index=False)
    say("treatment check", check); say("power", power)
    if a.power_only:
        return

    # ---- H196.1 ----
    T = d.T(Y, wf)
    used = X.loc[m.tp.isin(["C1D1", "C1D15"]).to_numpy()]
    cg = list(wf.genes)
    sets = L.pseudo_sets(used, L.pool_genes(used, L.estrogen_genes()), cg, a.npseudo, np.random.default_rng([SEED, 10]))
    pT = np.array([d.T(L.Yof(X, g), L.wf_pseudo(wf, g)) for g in sets])
    res = dict(test="H196.1", T=T, p_two_sided=p2(null, T), **{k: power[k] for k in ("null_mean", "null_sd", "mde", "estimable")},
               treatment_acted=check["acted"], pseudo_percentile=float((pT < T).mean() * 100), pseudo_q95=float(np.quantile(pT, .95)),
               pseudo_median=float(np.median(pT)), n_pairs=d.n)
    res["reading"] = ("uninformative (treatment check failed)" if not check["acted"] else
                      "not estimable: descriptive only" if not power["estimable"] else
                      "pass" if res["p_two_sided"] < .05 else "fail")
    pd.DataFrame([res]).to_csv(OUT / "tests.csv", index=False)
    say("H196.1", res)

    # ---- descriptives ----
    desc = []
    resid = L.residual_coh(m[["prolif", "purity"]].to_numpy())
    for lab, pre, post in (("BL -> C1D1 (anastrozole)", "BL", "C1D1"), ("BL -> C1D15", "BL", "C1D15"), ("C1D15 -> surgery", "C1D15", "Surg"),
                           ("C1D1 -> C1D15 (H196.1)", "C1D1", "C1D15")):
        dd = pairs(m, pre, post)
        nl = dd.null(Y, wf, 1000, np.random.default_rng([SEED, len(desc)]))
        t = dd.T(Y, wf)
        desc.append(dict(contrast=lab, n_pairs=dd.n, T=t, p_two_sided=p2(nl, t), null_sd=float(nl.std(ddof=1)),
                         T_resid_prolif_purity=dd.T(Y, wf, resid), d_prolif=float((pr[dd.post] - pr[dd.pre]).mean()),
                         d_estro=float((m.estro.to_numpy()[dd.post] - m.estro.to_numpy()[dd.pre]).mean())))
    for tp, g in m.groupby("tp"):
        desc.append(dict(contrast=f"coherence at {tp}", n_pairs=len(g), T=L.coh(Y, g.pos.to_numpy(), wf)))
    pd.DataFrame(desc).to_csv(OUT / "descriptive.csv", index=False)

    # ---- H196.2 (post hoc) ----
    cal = []
    for j, (name, key, pre, post, arms) in enumerate((("POETIC AI B->S", "poetic", "B", "S", ["AI"]), ("GSE20181 A->B", "letrozole", "A", "B", None),
                                                       ("GSE87411 base->w2", "z1031b", "base", "w2", None), ("GSE93204 C1D1->C1D15", "neopal", "C1D1", "C1D15", None))):
        mm, _, YY, ww = series(key)
        cal.append(calibration(name, mm, YY, ww, pairs(mm, pre, post, arms), a.nperm, np.random.default_rng([SEED, 20 + j])))
    cal = pd.DataFrame(cal)
    w = 1 / (cal.sd_T ** 2 + cal.sd_P ** 2)
    wn = 1 / cal.sd_numerator ** 2
    zn = float((cal.numerator_centred * wn).sum() / np.sqrt(wn.sum()))
    pooled = dict(series="pooled (inverse variance)", T_minus_P=float((cal.T_minus_P * w).sum() / w.sum()), sd_T=float(1 / np.sqrt(w.sum())),
                  z=float((cal.T_minus_P * w).sum() / np.sqrt(w.sum())), numerator_centred=float((cal.numerator_centred * wn).sum() / wn.sum()),
                  p_numerator_one_sided=float(stats.norm.sf(zn)))
    pooled["p_two_sided"] = float(2 * stats.norm.sf(abs(pooled["z"])))
    wt = 1 / cal.sd_T ** 2  # pooled observed change alone (descriptive)
    pooled.update(T_centred=float((cal.T_centred * wt).sum() / wt.sum()), sd_T_alone=float(1 / np.sqrt(wt.sum())), n_pairs=int(cal.n_pairs.sum()))
    pooled["informative"] = pooled["p_numerator_one_sided"] < .05
    cal = pd.concat([cal, pd.DataFrame([pooled])], ignore_index=True)
    cal.to_csv(OUT / "calibration_posthoc.csv", index=False)
    pd.set_option("display.width", 250)
    print(pd.DataFrame(desc).round(4).to_string(index=False))
    print(cal.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
