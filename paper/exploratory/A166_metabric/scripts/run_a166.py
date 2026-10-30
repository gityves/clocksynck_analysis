"""A166: METABRIC replication of the ER (H166.1, H166.2) and proliferation (H166.3) clock-coherence losses.
Order: (1) power.csv (null SDs and MDE = 2.49 x SD of the permutation nulls; no observed statistic is computed before it is written),
(2) observed T, one-sided permutation p, Holm over the three, (3) pseudo-clock specificity (200 variance-matched sets),
(4) sensitivities, gene floor, (5) TCGA BRCA secondary (H166.2, H166.3, raw p, same rules).
Seeds: primary nulls default_rng(166); pseudo sets default_rng(1660) METABRIC, default_rng(1661) TCGA; sensitivities default_rng([166, job]);
TCGA nulls default_rng(1662).
    python -I run_a166.py [--nperm 2000] [--npseudo 200] [--out DIR] [--jobs 4]"""
import argparse
import os
import sys
from pathlib import Path

for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "1")
sys.path.insert(0, str(Path(__file__).parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402
from scipy.stats import norm  # noqa: E402

import a166_lib as L  # noqa: E402
import clocksync as cs  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--nperm", type=int, default=2000); ap.add_argument("--npseudo", type=int, default=200)
ap.add_argument("--jobs", type=int, default=4); ap.add_argument("--out", default=str(L.OUT))
args = ap.parse_args()
OUT = Path(args.out); OUT.mkdir(parents=True, exist_ok=True)
LN2, K = L.LN2, norm.ppf(0.95) + norm.ppf(0.80)  # 2.49
SCANB = {"H166.1": -0.343, "H166.2": -0.268, "H166.3": -0.123}
NAMES = {"H166.1": "ER- minus ER+ within COHORT x purity tertile (= H127.2)",
         "H166.2": "ER- minus ER+ within PAM50 x proliferation tertile (= A127 post hoc)",
         "H166.3": "high minus low proliferation within PAM50 x purity tertile (= H127.1)"}
pl = lambda nl, o: cs.stats.perm_p(nl, o, "less", tie=0.0)
holm = lambda p: np.minimum(1, np.maximum.accumulate(np.sort(p) * (len(p) - np.arange(len(p)))))[np.argsort(np.argsort(p))]
Yof = lambda X, genes: LN2 * X[list(genes)].to_numpy(float)

m, X = L.load("metabric"); m["pos"] = np.arange(len(m)); assert (m.index == X.index).all()
Yc = Yof(X, L.GENES)
D = L.build(m)
for d in D: print(d.name, "groups", len(d.G), "tumours", d.n_tumours, flush=True)

# ---- (1) power, before any observed statistic -----------------------------------------------------------------------------------
rng = np.random.default_rng(166)
nulls = {d.name: d.null(Yc, L.WF, args.nperm, rng) for d in D}
pw = pd.DataFrame([dict(hypothesis=d.name, null_mean=nulls[d.name].mean(), null_sd=nulls[d.name].std(ddof=1),
                        MDE_80pct=K * nulls[d.name].std(ddof=1), scanb_effect=SCANB[d.name], n_groups=len(d.G), n_tumours=d.n_tumours,
                        n_perm=args.nperm) for d in D])
pw["scanb_effect_detectable"] = pw.MDE_80pct <= pw.scanb_effect.abs()
pw.to_csv(OUT / "power.csv", index=False); print(pw.round(4).to_string(), flush=True)

# ---- (2) observed, p, Holm ------------------------------------------------------------------------------------------------------
obs = {d.name: d.T(Yc, L.WF) for d in D}
p = np.array([pl(nulls[d.name], obs[d.name][0]) for d in D]); ph = holm(p)
pd.DataFrame({k: v for k, v in nulls.items()}).to_csv(OUT / "a166_null_draws.csv", index=False)

# ---- (3) specificity ------------------------------------------------------------------------------------------------------------
pool = L.pool_genes(X)
sets = L.pseudo_sets(X, pool, args.npseudo, np.random.default_rng(1660))
pd.DataFrame(sets).to_csv(OUT / "pseudo_sets_metabric.csv", index=False)
tp = {d.name: np.array([d.T(Yof(X, s), L.wf_for(s))[0] for s in sets]) for d in D}
pd.DataFrame(tp).to_csv(OUT / "a166_pseudo_T_draws.csv", index=False)
rows = []
for d, pp, hh in zip(D, p, ph):
    T, q05 = obs[d.name][0], np.quantile(tp[d.name], 0.05)
    verdict = "pass" if hh < 0.05 else "fail"
    spec = ("clock-specific" if T < q05 else "coherence loss not specific to the clock") if verdict == "pass" else \
           ("not applicable (test failed); T below pseudo q05" if T < q05 else "not applicable (test failed)")
    rows.append(dict(hypothesis=d.name, description=NAMES[d.name], T=T, p_one_lower=pp, p_holm=hh, null_sd=nulls[d.name].std(ddof=1),
                     MDE_80pct=K * nulls[d.name].std(ddof=1), scanb_effect=SCANB[d.name], n_groups=len(d.G), n_tumours=d.n_tumours,
                     pseudo_percentile=100 * np.mean(tp[d.name] <= T), pseudo_q05=q05, pseudo_median=np.median(tp[d.name]),
                     T_below_pseudo_q05=bool(T < q05), verdict=verdict, specificity=spec, n_perm=args.nperm, n_pseudo=len(sets)))
res = pd.DataFrame(rows); res.to_csv(OUT / "a166_tests.csv", index=False); print(res.drop(columns="description").round(4).to_string(), flush=True)

# ---- cells / strata tables ------------------------------------------------------------------------------------------------------
# the groups lose their key; recompute keys by re-running the grouping for the primary designs
a = m[m.er.notna()].copy(); a["pt"] = a.groupby("cohort").purity.transform(L.tertile)
b = a[a.pam50.notna()].copy(); b["ft"] = b.groupby("pam50").prolif.transform(L.tertile)
c = m[m.pam50.notna()].copy(); c["pt"] = c.groupby("pam50").purity.transform(L.tertile)
ce = []
for key, frame, d in ((["cohort", "pt"], a, D[0]), (["pam50", "ft"], b, D[1])):
    for k, g in frame.groupby(key, observed=True):
        i = g.pos.to_numpy(); e = d.er[i]
        if min((e == 0).sum(), (e == 1).sum()) >= L.MIN_ARM:
            ce.append(dict(hypothesis=d.name, key=str(k), n_er_neg=int((e == 0).sum()), n_er_pos=int((e == 1).sum()),
                           coh_er_neg=L.coh(Yc, i[e == 0], L.WF), coh_er_pos=L.coh(Yc, i[e == 1], L.WF),
                           prolif_mean_er_neg=g.prolif.to_numpy()[e == 0].mean(), prolif_mean_er_pos=g.prolif.to_numpy()[e == 1].mean()))
        else:
            ce.append(dict(hypothesis=d.name, key=str(k), n_er_neg=int((e == 0).sum()), n_er_pos=int((e == 1).sum()), dropped="< 10 per arm"))
ce = pd.DataFrame(ce); ce["diff"] = ce.coh_er_neg - ce.coh_er_pos; ce.to_csv(OUT / "cells.csv", index=False)
st = []
for k, g in c.groupby(["pam50", "pt"], observed=True):
    i = g.pos.to_numpy()
    if len(i) >= L.MIN_STRATUM:
        l = L.tertile(g.prolif.to_numpy())
        st.append(dict(key=str(k), n=len(i), n_low=int((l == 0).sum()), n_high=int((l == 2).sum()), coh_low=L.coh(Yc, i[l == 0], L.WF), coh_high=L.coh(Yc, i[l == 2], L.WF)))
    else:
        st.append(dict(key=str(k), n=len(i), dropped="< 30 tumours"))
st = pd.DataFrame(st); st["diff"] = st.coh_high - st.coh_low; st.to_csv(OUT / "strata.csv", index=False)
print(ce.round(3).to_string()); print(st.round(3).to_string(), flush=True)

# ---- (4) sensitivities ----------------------------------------------------------------------------------------------------------
def sens_job(label, mm, ercol, pamcol, by_cohort, which, seed, nperm):
    """one sensitivity: observed T, null (nperm), p, SD, MDE for the designs in `which`."""
    Y = Yc[mm.pos.to_numpy()]; mm = mm.copy(); mm["pos"] = np.arange(len(mm))  # subsets index their own rows
    ds = L.build(mm, ercol, pamcol, by_cohort); r = np.random.default_rng([166, seed]); out = []
    for d in ds:
        if d.name not in which: continue
        if len(d.G) == 0:
            out.append(dict(analysis=label, hypothesis=d.name, T=np.nan, n_groups=0, n_tumours=0)); continue
        o = d.T(Y, L.WF)[0]; nl = d.null(Y, L.WF, nperm, r)
        out.append(dict(analysis=label, hypothesis=d.name, T=o, p_one_lower=pl(nl, o), null_sd=nl.std(ddof=1), MDE_80pct=K * nl.std(ddof=1),
                        n_groups=len(d.G), n_tumours=d.n_tumours))
    return out
jobs = [("ER by ER_STATUS", m, "er_status_bin", "pam50", True, ("H166.1", "H166.2")),
        ("claudin-low as its own PAM50 stratum", m.assign(pam_cl=m.claudin.where(m.claudin.isin(L.SUBTYPES + ["claudin-low"]))), "er", "pam_cl", True, ("H166.2", "H166.3")),
        ("without COHORT stratification (cohort-wide purity tertile)", m, "er", "pam50", False, ("H166.1",))]
for c_ in sorted(m.cohort.unique()):
    jobs.append((f"COHORT {c_} only", m[m.cohort == c_], "er", "pam50", True, ("H166.1", "H166.2", "H166.3")))
res_s = Parallel(n_jobs=args.jobs)(delayed(sens_job)(*j[:6], 1 + k, args.nperm) for k, j in enumerate(jobs))
sens = pd.DataFrame([r for rs in res_s for r in rs]); sens.to_csv(OUT / "sensitivity.csv", index=False); print(sens.round(4).to_string(), flush=True)

# gene floor (clock genes, METABRIC)
g = X[L.GENES]; ern = m.er.to_numpy()
gf = pd.DataFrame({"gene": L.GENES, "mean": g.mean().values, "sd": g.std().values, "min": g.min().values, "q05": g.quantile(.05).values,
                   "mean_minus_floor": g.mean().values - L.FLOOR, "mean_minus_floor_in_sd": (g.mean().values - L.FLOOR) / g.std().values,
                   "frac_within_0.25_of_floor": (g < L.FLOOR + 0.25).mean().values,
                   "mean_er_neg": g[ern == 0].mean().values, "mean_er_pos": g[ern == 1].mean().values})
gf["array_floor"] = L.FLOOR; gf["gene_mean_sd_all_genes"] = X.mean().std()
gf.to_csv(OUT / "gene_floor.csv", index=False); print(gf.round(3).to_string(), flush=True)

# ---- (5) TCGA BRCA secondary ------------------------------------------------------------------------------------------------------
mt, Xt = L.load("tcga"); mt["pos"] = np.arange(len(mt)); Yt = Yof(Xt, L.GENES)
Dt = L.build(mt)[1:]  # H166.2, H166.3 (H166.1 is not part of the TCGA secondary)
rng = np.random.default_rng(1662)
setst = L.pseudo_sets(Xt, L.pool_genes(Xt), args.npseudo, np.random.default_rng(1661)); pd.DataFrame(setst).to_csv(OUT / "pseudo_sets_tcga.csv", index=False)
rows = []
for d in Dt:
    if len(d.G) == 0:
        rows.append(dict(hypothesis=d.name, n_groups=0, note="no cell with >= 10 tumours per ER arm" if d.kind == "er" else "no stratum >= 30")); continue
    o = d.T(Yt, L.WF); nl = d.null(Yt, L.WF, args.nperm, rng); t_ = np.array([d.T(Yof(Xt, s), L.wf_for(s))[0] for s in setst])
    rows.append(dict(hypothesis=d.name, T=o[0], p_one_lower=pl(nl, o[0]), null_sd=nl.std(ddof=1), MDE_80pct=K * nl.std(ddof=1), n_groups=len(d.G),
                     n_tumours=d.n_tumours, pseudo_percentile=100 * np.mean(t_ <= o[0]), pseudo_q05=np.quantile(t_, .05), T_below_pseudo_q05=bool(o[0] < np.quantile(t_, .05)),
                     per_group_diff=";".join(f"{x:.3f}" for x in o[1]), n_perm=args.nperm, n_pseudo=len(t_)))
pd.DataFrame(rows).to_csv(OUT / "tcga_secondary.csv", index=False); print(pd.DataFrame(rows).round(4).to_string(), flush=True)
# TCGA cell table (all PAM50 x proliferation-tertile cells with their arm sizes)
bt = mt[mt.er.notna() & mt.pam50.notna()].copy(); bt["ft"] = bt.groupby("pam50").prolif.transform(L.tertile)
tc = []
for k, gg in bt.groupby(["pam50", "ft"]):
    tc.append(dict(key=str(k), n_er_neg=int((gg.er == 0).sum()), n_er_pos=int((gg.er == 1).sum())))
ct = mt[mt.pam50.notna()].copy(); ct["pt"] = ct.groupby("pam50").purity.transform(L.tertile)
for k, gg in ct.groupby(["pam50", "pt"]): tc.append(dict(key="H166.3 " + str(k), n=len(gg)))
pd.DataFrame(tc).to_csv(OUT / "tcga_cells.csv", index=False); print(pd.DataFrame(tc).to_string())
print("done")
