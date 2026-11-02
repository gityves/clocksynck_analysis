"""A192 scoring (registered in `paper/docs/registered_experiments.md`, commit a6d5922): is the prognostic value of clock genes in
breast cancer (Cadenas et al. 2014) carried by proliferation?

Mainz (GSE11121), Rotterdam (GSE2034) and Transbig (GSE7390): node-negative, systemically untreated, U133A, distant-metastasis-
free survival (`paper/common/datasets/breast_prognosis_geo/build.py`). Cox (PHReg, Efron), strata = cohort, readouts
standardised, within-cohort z-scores of log2 expression.
H192.1 (pipeline check): univariate HR of the 8 genes Cadenas reported.
H192.2 (primary): clock score S = mean z of CLOCK, PER1, PER2, PER3, CRY2, NPAS2, RORC; (a) unadjusted z against 1,000
proliferation-matched pseudo-scores (20 nearest pool genes in standardised mean, SD and correlation with G2M); (b) S adjusted
for G2M, ESR1, ERBB2.
H192.3 (primary): per-tumour in-plane log-odds (13 clock genes) adjusted for G2M, ESR1, ERBB2; 200 expression-matched pseudo
sets (A191 construction). Holm over H192.2(b) and H192.3. Power step first. Only summary statistics are written.

    python score_a192.py [--check] [--power-only]
"""
import argparse
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve()
ROOT = HERE.parents[4]
EXP = ROOT / "paper" / "exploratory"
sys.path[:0] = [str(EXP / "A166_metabric" / "scripts"), str(EXP / "C8_internal_time_outcomes" / "scripts")]
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import norm  # noqa: E402
from statsmodels.duration.hazard_regression import PHReg  # noqa: E402

import a166_lib as L  # noqa: E402  (GENES, WF, sh: hallmark_g2m, PURITY)
from score_a191 import cox, holm, mde, pseudo_sets_generic  # noqa: E402
from clocksync._core.controls import CLOCK_LIKE  # noqa: E402
from clocksync.decompose import sample_log_odds  # noqa: E402
from clocksync.genes import canonical  # noqa: E402
from clocksync.waveform import Waveform  # noqa: E402

OUT = HERE.parents[1] / "output"
DATA = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian/breast_prognosis_geo")
SEED, LN2 = 192, np.log(2)
PROTECTIVE = ["CLOCK", "PER1", "PER2", "PER3", "CRY2", "NPAS2", "RORC"]
CADENAS = {**{g: "protective" for g in PROTECTIVE}, "TIMELESS": "adverse"}
CADENAS_R_PROLIF = {"CRY2": -0.612, "PER2": -0.440, "PER3": -0.409, "TIMELESS": 0.674}  # their Table, proliferation metagene
T0 = time.time()


def log(*a):
    print(f"[{(time.time() - T0) / 60:5.1f} min]", *a, flush=True)


def load():
    E = pd.read_parquet(DATA / "expr.parquet")
    c = pd.read_csv(DATA / "clinical.csv", index_col=0).loc[E.index]
    Z = E.groupby(c.cohort).transform(lambda x: (x - x.mean()) / x.std())  # within-cohort z
    g2m = [g for g in L.sh.hallmark_g2m() if g in Z.columns]
    c["g2m"] = Z[g2m].mean(axis=1)
    for g in ("ESR1", "ERBB2", "MKI67", "AURKA"):
        c[g] = Z[g]
    return E, Z, c, g2m


def score(Z, genes):
    return Z[genes].mean(axis=1).to_numpy()


def prolif_matched_sets(E, Z, c, genes, pool, n, rng):
    """Each gene replaced by one of its 20 nearest pool genes in standardised (mean, SD, correlation with G2M); mean and SD
    are averaged over the per-cohort values (cohorts differ in scale: MAS5 vs RMA)."""
    mu = E.groupby(c.cohort).mean().mean()
    sd = E.groupby(c.cohort).std().mean()
    r = Z.corrwith(c.g2m)
    F = pd.DataFrame({"mu": mu, "sd": sd, "r": r})
    Fp = F.loc[pool]
    P = ((Fp - Fp.mean()) / Fp.std()).to_numpy()
    T = ((F.loc[genes] - Fp.mean()) / Fp.std()).to_numpy()
    near = [np.argsort(((P - T[k]) ** 2).sum(1))[:40] for k in range(len(genes))]
    sets = []
    for _ in range(n):
        used, s = set(), []
        for k in rng.permutation(len(genes)):
            cand = [j for j in near[k] if j not in used][:20]
            j = int(rng.choice(cand))
            used.add(j)
            s.append((k, pool[j]))
        sets.append([g for _, g in sorted(s)])
    match = pd.DataFrame({"clock_r_g2m": F.loc[genes, "r"].mean(), "pseudo_r_g2m_mean": np.mean([F.loc[s, "r"].mean() for s in sets]),
                          "clock_sd": F.loc[genes, "sd"].mean(), "pseudo_sd_mean": np.mean([F.loc[s, "sd"].mean() for s in sets])},
                         index=[0])
    return sets, match


def log_odds(E, c, genes, wf):
    Y = E[genes].to_numpy(float) * LN2
    Y = Y - pd.DataFrame(Y).groupby(c.cohort.to_numpy()).transform("mean").to_numpy()
    return sample_log_odds(Y, wf)


def fit(c, x, cov, tcol="time", ecol="event"):
    d = c.assign(x=x).dropna(subset=[tcol, ecol, "x"] + cov)
    Zc = d[cov].to_numpy(float) if cov else None
    r = cox(d[tcol].to_numpy(float), d[ecol].to_numpy(float), d.x.to_numpy(float), Zc, strata=d.cohort.to_numpy())
    return dict(n=len(d), events=int(d[ecol].sum()), **{k: v for k, v in r.items() if k != "llf"})


def logb(c, x, cov, idx):
    d = c.assign(x=x).iloc[idx].reset_index(drop=True)
    d = d.assign(cohort=c.cohort.to_numpy()[idx])
    x = (d.x - d.x.mean()) / d.x.std()
    X = np.column_stack([x, d[cov].to_numpy(float)]) if cov else x.to_numpy()[:, None]
    return PHReg(d.time.to_numpy(float), X, status=d.event.to_numpy(float), strata=d.cohort.to_numpy(), ties="efron").fit(disp=0).params[0]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--check", action="store_true", help="20 / 5 pseudo sets, 20 bootstraps, outputs suffixed _check")
    ap.add_argument("--power-only", action="store_true", help="stop after the power step (no outcome association)")
    a = ap.parse_args()
    n_pm, n_ps, n_boot = (20, 5, 20) if a.check else (1000, 200, 1000)
    sfx = "_check" if a.check else ""
    OUT.mkdir(exist_ok=True)
    rng = np.random.default_rng(SEED)
    E, Z, c, g2m = load()
    genes13 = [g for g in L.GENES if g in E.columns]
    wf13 = L.WF.subset(genes13)
    cov = ["g2m", "ESR1", "ERBB2"]
    c["S"] = score(Z, PROTECTIVE)
    c["s"] = log_odds(E, c, genes13, wf13)
    log("cohorts", c.groupby("cohort").agg(n=("event", "size"), events=("event", "sum")).to_dict(), "clock genes", len(genes13))

    # ---- power step (no outcome association) ----
    power = []
    for test, x in (("H192.2b", c.S), ("H192.3", c.s)):
        d = c.dropna(subset=["time", "event"] + cov)
        md, r2 = mde(x.loc[d.index].to_numpy(), np.column_stack([d[cov], pd.get_dummies(d.cohort, drop_first=True, dtype=float)]),
                     d.event.sum())
        power.append(dict(test=test, n=len(d), events=int(d.event.sum()), r2_on_covariates=r2, mde_hr_per_sd=md, estimable=md <= 1.5))
    power = pd.DataFrame(power)
    power.to_csv(OUT / f"a192_power{sfx}.csv", index=False)
    log("power\n" + power.round(4).to_string(index=False))
    if a.power_only:
        return

    # ---- H192.1 reproduction ----
    rep = []
    for g, d in CADENAS.items():
        r = fit(c, c[g] if g in c else Z[g], [])
        rep.append(dict(gene=g, published=d, r_g2m=float(np.corrcoef(Z[g], c.g2m)[0, 1]), r_prolif_cadenas=CADENAS_R_PROLIF.get(g, np.nan), **r))
    rep = pd.DataFrame(rep)
    ps = np.sort(rep.p.to_numpy())
    m = len(ps)
    bh = np.minimum.accumulate((m / np.arange(m, 0, -1)) * ps[::-1])[::-1]  # BH step-up on the sorted p
    rep["q_bh"] = rep.p.map(dict(zip(ps, np.minimum(bh, 1))))
    rep["reproduced"] = (rep.q_bh < .05) & np.where(rep.published == "protective", rep.hr < 1, rep.hr > 1)
    rep.to_csv(OUT / f"a192_reproduction{sfx}.csv", index=False)
    log("H192.1\n" + rep.round(4).to_string(index=False))

    # ---- H192.2 ----
    clock = set(CLOCK_LIKE) | set(g2m) | set(L.sh.PURITY) | set(L.GENES) | set(CADENAS)
    pool = [g for g in E.columns if canonical(g) not in clock and g not in clock]
    unadj = fit(c, c.S, [])
    sets, match = prolif_matched_sets(E, Z, c, PROTECTIVE, pool, n_pm, rng)
    zp = np.array([fit(c, score(Z, s), [])["z"] for s in sets])
    lo, hi = np.quantile(zp, [.025, .975])
    adj = fit(c, c.S, cov)
    zp_adj = np.array([fit(c, score(Z, s), cov)["z"] for s in sets])
    tests = [dict(test="H192.2a unadjusted S vs prolif-matched", **unadj, pseudo_lo=lo, pseudo_hi=hi,
                  pseudo_pct=100 * np.mean(zp < unadj["z"]), clock_specific=not (lo <= unadj["z"] <= hi)),
             dict(test="H192.2b S adjusted (G2M, ESR1, ERBB2)", **adj, pseudo_lo=np.quantile(zp_adj, .025),
                  pseudo_hi=np.quantile(zp_adj, .975), pseudo_pct=100 * np.mean(zp_adj < adj["z"]))]
    log("H192.2", tests)

    # ---- H192.3 ----
    r3 = fit(c, c.s, cov)
    allg = E
    pool3 = [g for g in E.columns if canonical(g) not in set(CLOCK_LIKE) | set(g2m) | set(L.sh.PURITY) and g not in genes13]
    zs = []
    for gs in pseudo_sets_generic(allg, genes13, pool3, n_ps, rng):
        sp = log_odds(E, c, gs, Waveform(tuple(gs), wf13.amp, wf13.psi, None, None, "set"))
        zs.append(fit(c, sp, cov)["z"])
    zs = np.abs(zs)
    tests.append(dict(test="H192.3 in-plane log-odds adjusted", **r3, pseudo_q95_abs_z=np.quantile(zs, .95),
                      pseudo_pct=100 * np.mean(zs < abs(r3["z"]))))
    tt = pd.DataFrame(tests)
    prim = tt.test.str.startswith(("H192.2b", "H192.3"))
    tt.loc[prim, "p_holm"] = holm(tt.loc[prim, "p"])
    a_spec, b = bool(tt.iloc[0].clock_specific), tt.iloc[1]
    tt.loc[0, "reading_H192.2"] = ("carried by proliferation" if (not a_spec and b.p >= .05) else
                                   "clock-specific prognosis" if (a_spec and b.p < .05 and b.hr < 1) else "mixed")
    tt.to_csv(OUT / f"a192_tests{sfx}.csv", index=False)
    match.to_csv(OUT / f"a192_prolif_match{sfx}.csv", index=False)

    # ---- secondary (not counted) ----
    sec = []
    for lab, cv in (("S + MKI67 only", ["MKI67"]), ("S + AURKA only", ["AURKA"]), ("S + G2M only", ["g2m"]),
                    ("S + ESR1 + ERBB2 (no proliferation)", ["ESR1", "ERBB2"])):
        sec.append(dict(test=lab, **fit(c, c.S, cv)))
    mt = c[c.cohort.isin(["Mainz", "Transbig"])]
    for g in ("PER3", "RORC", "TIMELESS"):
        for lab, cv in (("clinical", ["grade", "size", "ESR1", "ERBB2"]), ("clinical + G2M", ["grade", "size", "ESR1", "ERBB2", "g2m"])):
            sec.append(dict(test=f"Cadenas multivariate {g}, Mainz+Transbig, {lab}", **fit(mt, Z.loc[mt.index, g], cv)))
    sec.append(dict(test="TIMELESS + G2M", **fit(c, Z.TIMELESS, ["g2m"])))
    sec.append(dict(test="G2M alone (positive control)", **fit(c, c.g2m, [])))
    sec.append(dict(test="G2M + ESR1 + ERBB2 (positive control)", **fit(c, c.g2m, ["ESR1", "ERBB2"])))
    for coh in c.cohort.unique():
        cc = c[c.cohort == coh]
        sec.append(dict(test=f"S unadjusted, {coh}", **fit(cc, cc.S, [])))
        sec.append(dict(test=f"S adjusted, {coh}", **fit(cc, cc.S, cov)))
        sec.append(dict(test=f"in-plane adjusted, {coh}", **fit(cc, cc.s, cov)))
    erp = c[c.ESR1 > c.groupby("cohort").ESR1.transform("median")] if c.er.isna().any() else c[c.er == 1]
    sec.append(dict(test="S adjusted, ESR1 above cohort median (ER+ proxy)", **fit(erp, erp.S, cov)))
    sec.append(dict(test="in-plane unadjusted", **fit(c, c.s, [])))
    # share of the unadjusted log HR of S absorbed by G2M (bootstrap over tumours, within cohort)
    d = c.dropna(subset=["time", "event"] + cov)
    b0, b1 = logb(d, d.S, [], np.arange(len(d))), logb(d, d.S, ["g2m"], np.arange(len(d)))
    shares = []
    for _ in range(n_boot):
        idx = np.concatenate([rng.choice(np.flatnonzero(d.cohort.to_numpy() == k), (d.cohort == k).sum()) for k in d.cohort.unique()])
        u, v = logb(d, d.S, [], idx), logb(d, d.S, ["g2m"], idx)
        shares.append(1 - v / u)
    sec.append(dict(test="share of log HR(S) absorbed by G2M", estimate=1 - b1 / b0, lo=np.quantile(shares, .025),
                    hi=np.quantile(shares, .975)))
    # pooling with A191 H191.1 (inverse variance on log HR) and TOST, margin HR 1.15 per SD
    a191 = pd.read_csv(EXP / "C8_internal_time_outcomes" / "output" / "tests.csv").set_index("test").loc["H191.1"]
    est = []
    for r in (a191, tt.iloc[2]):
        bb = np.log(r.hr)
        se = (np.log(r.hr_hi) - np.log(r.hr_lo)) / (2 * 1.959964)
        est.append((bb, se))
    w = np.array([1 / s ** 2 for _, s in est])
    bp = float(np.sum(w * np.array([b_ for b_, _ in est])) / w.sum())
    sp = float(1 / np.sqrt(w.sum()))
    Q = float(np.sum(w * (np.array([b_ for b_, _ in est]) - bp) ** 2))
    marg = np.log(1.15)
    p_tost = max(norm.sf((bp + marg) / sp), norm.cdf((bp - marg) / sp))
    sec.append(dict(test="in-plane: A191 H191.1 + H192.3 pooled (fixed effect)", hr=np.exp(bp), hr_lo=np.exp(bp - 1.96 * sp),
                    hr_hi=np.exp(bp + 1.96 * sp), cochran_q=Q, p_tost_margin_1_15=p_tost))
    sec = pd.DataFrame(sec)
    sec.to_csv(OUT / f"a192_secondary{sfx}.csv", index=False)
    pd.set_option("display.width", 250)
    print(tt.round(4).to_string(index=False))
    print(match.round(3).to_string(index=False))
    print(sec.round(4).to_string(index=False))
    log("done")


if __name__ == "__main__":
    main()
