"""A127: SCAN-B (GSE96058) clock-gene coherence vs proliferation (H127.1) and ER status (H127.2), pseudo-clock specificity
companion, MDE, Holm, exploratory (survival, coherence per PAM50, purity companion). python -I run_a127.py [--nperm 2000]
Coherence = clocksync.decompose.components (A97 definition); y = ln2 * log2(FPKM + 0.1), centred within the group; no Poisson term."""
import argparse, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
import numpy as np, pandas as pd
import a127_lib as L
import clocksync as cs
from clocksync._core.controls import CLOCK_LIKE
from clocksync.genes import canonical
from statsmodels.duration.hazard_regression import PHReg
from scipy.stats import norm

ap = argparse.ArgumentParser(); ap.add_argument("--nperm", type=int, default=2000); ap.add_argument("--nboot", type=int, default=1000)
ap.add_argument("--npseudo", type=int, default=200); args = ap.parse_args()
OUT, LN2 = L.OUT, np.log(2)
rng = np.random.default_rng(127)
m, X = L.load(); n = len(m)
Yc = LN2 * X[L.GENES].to_numpy(float)
g2m_all = set(L.sh.hallmark_g2m())

# ---- pseudo-clock sets: pool = detected in >=90% samples, not clock-like / G2M / purity; matched on mean and sd -------------------
pool = [g for g in X.columns if canonical(g) not in CLOCK_LIKE and g not in g2m_all and g not in L.sh.PURITY and g not in L.GENES]
raw = pd.read_csv(OUT / "prep_scale.csv", header=None, index_col=0)[1]
pm, ps = X[pool].mean(), X[pool].std()
zs = lambda a, ref: (a - ref.mean()) / ref.std()
P = np.column_stack([zs(pm.values, pm), zs(ps.values, ps)])
target = np.column_stack([zs(X[L.GENES].mean().values, pm), zs(X[L.GENES].std().values, ps)])
sets = []
for _ in range(args.npseudo):
    used, s = set(), []
    for t in rng.permutation(len(L.GENES)):
        near = [j for j in np.argsort(((P - target[t]) ** 2).sum(1))[:40] if j not in used][:20]
        j = int(rng.choice(near)); used.add(j); s.append((t, pool[j]))
    sets.append([g for _, g in sorted(s)])
pd.DataFrame(sets).to_csv(OUT / "pseudo_sets_variance_matched.csv", index=False)
mv = pd.DataFrame({"clock_mean": X[L.GENES].mean().values, "clock_sd": X[L.GENES].std().values,
                   "pseudo_mean_avg": [X[[s[i] for s in sets]].mean().mean() for i in range(14)],
                   "pseudo_sd_avg": [X[[s[i] for s in sets]].std().mean() for i in range(14)]}, index=L.GENES)
mv.to_csv(OUT / "pseudo_match_quality.csv"); print(mv.round(2))
wfs = {"clock": L.WF}
pseudo_wfs = [L.wf_for(s) for s in sets]
api = cs.paired.pseudo_waveforms(X[pool + L.GENES].mean(), args.npseudo, np.random.default_rng(1270), L.WF)  # mean-only (clocksync API), sensitivity
Ys = {"clock": Yc}
def Yof(genes): return LN2 * X[list(genes)].to_numpy(float)
Ypseudo = [Yof(s) for s in sets]; Yapi = [Yof(w.genes) for w in api]

# ---- strata and labels --------------------------------------------------------------------------------------------------------------
m["pur_t"] = m.groupby("pam50").purity.transform(lambda s: L.tertile(s))
m["pur_t_cohort"] = L.tertile(m.purity)
def strata_of(key):
    out = []
    for _, g in m.groupby(key):
        idx = g.index.to_numpy()
        if len(idx) >= 30: out.append(np.array([m.index.get_loc(i) for i in idx]))
    return out
S1 = strata_of(["pam50", "pur_t"])
S1c = strata_of(["pam50", "pur_t_cohort"])  # sensitivity: cohort-wide purity tertile
lab1 = [L.tertile(m.prolif.to_numpy()[i]) for i in S1]
lab1c = [L.tertile(m.prolif.to_numpy()[i]) for i in S1c]
def T1(Y, wf, S, labs, perms=None):
    d = []
    for k, (i, l) in enumerate(zip(S, labs)):
        l = l if perms is None else l[perms[k]]
        d.append(L.coh(Y, i[l == 2], wf) - L.coh(Y, i[l == 0], wf))
    return float(np.mean(d))
# H127.2
ptert = [np.flatnonzero(m.pur_t_cohort.to_numpy() == t) for t in range(3)]
er = m.er.to_numpy()
def T2(Y, wf, groups=ptert, lab=er, perms=None):
    d = []
    for k, i in enumerate(groups):
        l = lab[i] if perms is None else lab[i][perms[k]]
        d.append(L.coh(Y, i[l == 0], wf) - L.coh(Y, i[l == 1], wf))
    return float(np.mean(d))
ptert_prol = [np.flatnonzero(L.tertile(m.prolif.to_numpy()) == t) for t in range(3)]  # sensitivity: within proliferation tertile

obs1, obs2 = T1(Yc, L.WF, S1, lab1), T2(Yc, L.WF)
obs2s = T2(Yc, L.WF, ptert_prol)
print("obs", obs1, obs2, obs2s, flush=True)
null1, null2, null2s, nullpur = [], [], [], []
S_pam = [np.flatnonzero(m.pam50.to_numpy() == s) for s in L.SUBTYPES]
purlab = [L.tertile(m.purity.to_numpy()[i]) for i in S_pam]
def Tpur(Y, wf, perms=None):  # companion: low-purity tertile (high purity score) minus high-purity tertile, within PAM50
    d = []
    for k, (i, l) in enumerate(zip(S_pam, purlab)):
        l = l if perms is None else l[perms[k]]
        d.append(L.coh(Y, i[l == 2], wf) - L.coh(Y, i[l == 0], wf))
    return float(np.mean(d))
obspur = Tpur(Yc, L.WF)
for b in range(args.nperm):
    null1.append(T1(Yc, L.WF, S1, lab1, [rng.permutation(len(i)) for i in S1]))
    null2.append(T2(Yc, L.WF, perms=[rng.permutation(len(i)) for i in ptert]))
    null2s.append(T2(Yc, L.WF, ptert_prol, perms=[rng.permutation(len(i)) for i in ptert_prol]))
    nullpur.append(Tpur(Yc, L.WF, [rng.permutation(len(i)) for i in S_pam]))
    if b % 250 == 0: print("perm", b, flush=True)
null1, null2, null2s, nullpur = map(np.array, (null1, null2, null2s, nullpur))
pl = lambda nl, o: cs.stats.perm_p(nl, o, "less", tie=0.0)
p1, p2 = pl(null1, obs1), pl(null2, obs2)
# Holm over H127.1-2
ps_ = np.array([p1, p2]); o = np.argsort(ps_); holm = np.empty(2)
run = 0
for r, k in enumerate(o):
    run = max(run, (2 - r) * ps_[k]); holm[k] = min(1, run)
z_a, z_b = norm.ppf(0.95), norm.ppf(0.80)
mde = lambda nl: (z_a + z_b) * np.std(nl, ddof=1)
res = pd.DataFrame([
    dict(test="H127.1 high-low proliferation (PAM50 x purity strata)", T=obs1, null_mean=null1.mean(), null_sd=null1.std(ddof=1),
         null_q05=np.quantile(null1, .05), p_one_lower=p1, p_holm=holm[0], MDE_80pct=mde(null1), n_strata=len(S1)),
    dict(test="H127.2 ER- minus ER+ (within purity tertile)", T=obs2, null_mean=null2.mean(), null_sd=null2.std(ddof=1),
         null_q05=np.quantile(null2, .05), p_one_lower=p2, p_holm=holm[1], MDE_80pct=mde(null2), n_strata=3),
    dict(test="sensitivity ER- minus ER+ within proliferation tertile (not tested)", T=obs2s, null_mean=null2s.mean(),
         null_sd=null2s.std(ddof=1), null_q05=np.quantile(null2s, .05), p_one_lower=pl(null2s, obs2s), p_holm=np.nan,
         MDE_80pct=mde(null2s), n_strata=3),
    dict(test="companion: low-purity minus high-purity tertile (within PAM50)", T=obspur, null_mean=nullpur.mean(),
         null_sd=nullpur.std(ddof=1), null_q05=np.quantile(nullpur, .05), p_one_lower=pl(nullpur, obspur), p_holm=np.nan,
         MDE_80pct=mde(nullpur), n_strata=5)]).assign(n_perm=args.nperm, n_tumours=n)
res.to_csv(OUT / "a127_tests.csv", index=False); print(res.round(4).to_string(), flush=True)
pd.DataFrame(dict(null_h127_1=null1, null_h127_2=null2)).to_csv(OUT / "a127_null_draws.csv", index=False)

# ---- specificity companion ------------------------------------------------------------------------------------------------------------
def spec(name, obs, Tfun, sets_Y, wfs_):
    tp = np.array([Tfun(y, w) for y, w in zip(sets_Y, wfs_)])
    return dict(test=name, T_clock=obs, T_pseudo_mean=tp.mean(), T_pseudo_q05=np.quantile(tp, .05), T_pseudo_median=np.median(tp),
                pseudo_percentile_of_clock=100 * np.mean(tp <= obs), clock_below_q05=bool(obs < np.quantile(tp, .05)), n_sets=len(tp)), tp
rows, tps = [], {}
for nm, Ysets, wl in (("variance-matched (primary)", Ypseudo, pseudo_wfs), ("mean-matched (clocksync pseudo_waveforms)", Yapi, api)):
    for hn, f, ob in (("H127.1", lambda y, w: T1(y, w, S1, lab1), obs1), ("H127.2", lambda y, w: T2(y, w), obs2),
                      ("purity companion", lambda y, w: Tpur(y, w), obspur)):
        r, tp = spec(f"{hn} / {nm}", ob, f, Ysets, wl); rows.append(r); tps[f"{hn}|{nm}"] = tp
pd.DataFrame(rows).to_csv(OUT / "a127_specificity.csv", index=False); print(pd.DataFrame(rows).round(4).to_string(), flush=True)
pd.DataFrame({k: v for k, v in tps.items()}).to_csv(OUT / "a127_pseudo_T_draws.csv", index=False)

# ---- sensitivities: cohort-wide purity tertile strata; platform-centred; ----------------------------------------------------------------
sens = []
sens.append(dict(analysis="H127.1 strata with cohort-wide purity tertile", T=T1(Yc, L.WF, S1c, lab1c)))
nl = np.array([T1(Yc, L.WF, S1c, lab1c, [rng.permutation(len(i)) for i in S1c]) for _ in range(args.nperm)])
sens[-1].update(p_one_lower=pl(nl, sens[-1]["T"]), null_sd=nl.std(ddof=1), MDE_80pct=mde(nl), n_strata=len(S1c))
Yp = Yc.copy()
for pf in (0, 1): Yp[m.platform_i.to_numpy() == pf] -= Yp[m.platform_i.to_numpy() == pf].mean(0)
sens.append(dict(analysis="H127.1, expression centred per sequencing platform", T=T1(Yp, L.WF, S1, lab1)))
nl = np.array([T1(Yp, L.WF, S1, lab1, [rng.permutation(len(i)) for i in S1]) for _ in range(args.nperm // 2)])
sens[-1].update(p_one_lower=pl(nl, sens[-1]["T"]), null_sd=nl.std(ddof=1), MDE_80pct=mde(nl), n_strata=len(S1))
sens.append(dict(analysis="H127.2, expression centred per sequencing platform", T=T2(Yp, L.WF)))
nl = np.array([T2(Yp, L.WF, perms=[rng.permutation(len(i)) for i in ptert]) for _ in range(args.nperm // 2)])
sens[-1].update(p_one_lower=pl(nl, sens[-1]["T"]), null_sd=nl.std(ddof=1), MDE_80pct=mde(nl), n_strata=3)
# A97 check: proliferation within purity halves (within PAM50), tertile within each half
d = []
for half in (0, 1):
    sub = []
    for i in S_pam:
        h = L.tertile(m.purity.to_numpy()[i], 2) == half
        j = i[h]; sub.append((j, L.tertile(m.prolif.to_numpy()[j])))
    f = lambda perms=None: np.mean([L.coh(Yc, j[(l if perms is None else l[perms[k]]) == 2], L.WF) - L.coh(Yc, j[(l if perms is None else l[perms[k]]) == 0], L.WF) for k, (j, l) in enumerate(sub)])
    o_ = f(); nl = np.array([f([rng.permutation(len(j)) for j, _ in sub]) for _ in range(args.nperm // 2)])
    sens.append(dict(analysis=f"proliferation within {'high' if half == 0 else 'low'}-purity half (within PAM50)", T=o_, p_one_lower=pl(nl, o_),
                     null_sd=nl.std(ddof=1), MDE_80pct=mde(nl), n_strata=len(sub)))
pd.DataFrame(sens).to_csv(OUT / "a127_sensitivity.csv", index=False); print(pd.DataFrame(sens).round(4).to_string(), flush=True)

# ---- per-stratum descriptive ----------------------------------------------------------------------------------------------------------
rows = []
keys = [k for k, g in m.groupby(["pam50", "pur_t"])]
for (k, i, l) in zip(keys, S1, lab1):
    rows.append(dict(pam50=k[0], purity_tertile=k[1], n=len(i), n_per_tertile=int((l == 0).sum()), coh_low=L.coh(Yc, i[l == 0], L.WF),
                     coh_mid=L.coh(Yc, i[l == 1], L.WF), coh_high=L.coh(Yc, i[l == 2], L.WF)))
pd.DataFrame(rows).assign(d=lambda d: d.coh_high - d.coh_low).to_csv(OUT / "a127_strata_descriptive.csv", index=False)

# ---- exploratory: coherence per PAM50 (patient bootstrap) --------------------------------------------------------------------------------
rows = []
for s, i in zip(L.SUBTYPES, S_pam):
    b = [L.coh(Yc, rng.choice(i, len(i)), L.WF) for _ in range(args.nboot)]
    rows.append(dict(pam50=s, n=len(i), coherence=L.coh(Yc, i, L.WF), lo=np.quantile(b, .025), hi=np.quantile(b, .975)))
for lab, i in (("ER-", np.flatnonzero(er == 0)), ("ER+", np.flatnonzero(er == 1))):
    b = [L.coh(Yc, rng.choice(i, len(i)), L.WF) for _ in range(args.nboot)]
    rows.append(dict(pam50=lab, n=len(i), coherence=L.coh(Yc, i, L.WF), lo=np.quantile(b, .025), hi=np.quantile(b, .975)))
pd.DataFrame(rows).to_csv(OUT / "a127_coherence_per_pam50.csv", index=False); print(pd.DataFrame(rows).round(3).to_string(), flush=True)

# ---- exploratory: survival, off-plane score (A97 off_raw without normals: deviation from the cohort mean; no Poisson term) -------------
M_, Q_ = np.eye(14) - 1 / 14, np.eye(14) - L.PL
dev = (Yc - Yc.mean(0)) @ M_
off_raw = np.sum((dev @ Q_) ** 2, 1) / np.maximum(np.sum(dev ** 2, 1), 1e-12)
d = m.assign(off=off_raw)
d["time"] = pd.to_numeric(d["overall survival days"]); d["event"] = pd.to_numeric(d["overall survival event"]).astype(int)
d["age10"] = pd.to_numeric(d["age at diagnosis"]) / 10; d["size_cm"] = pd.to_numeric(d["tumor size"]) / 10
d["node"] = (d["lymph node status"] == "NodePositive").astype(float).where(d["lymph node status"].notna())
for s in ("LumB", "Her2", "Basal", "Normal"): d["pam_" + s] = (d.pam50 == s).astype(float)
pam = ["pam_LumB", "pam_Her2", "pam_Basal", "pam_Normal"]
def cox(df, cols):
    mod = PHReg(df.time.to_numpy(float), df[cols].astype(float).to_numpy(), status=df.event.to_numpy(int), ties="efron").fit(disp=False)
    return mod.params, mod.bse
rows = []
d["z_off"] = (d.off - d.off.mean()) / d.off.std(); d["z_prolif"] = (d.prolif - d.prolif.mean()) / d.prolif.std()
for nm, cols in (("off-plane + age, size, node, PAM50 (registered)", ["z_off", "age10", "size_cm", "node"] + pam),
                 ("off-plane + age, size, node, PAM50 + proliferation (post hoc)", ["z_off", "age10", "size_cm", "node", "z_prolif"] + pam),
                 ("proliferation + age, size, node, PAM50 (post hoc, reference)", ["z_prolif", "age10", "size_cm", "node"] + pam)):
    dd = d.dropna(subset=cols + ["time", "event"]); dd = dd[dd.time > 0]
    b, se = cox(dd, cols)
    rows.append(dict(model=nm, n=len(dd), events=int(dd.event.sum()), logHR_per_SD=b[0], se=se[0], HR=np.exp(b[0]),
                     lo=np.exp(b[0] - 1.96 * se[0]), hi=np.exp(b[0] + 1.96 * se[0]), p=float(2 * (1 - norm.cdf(abs(b[0] / se[0]))))))
pd.DataFrame(rows).to_csv(OUT / "a127_survival.csv", index=False); print(pd.DataFrame(rows).round(4).to_string())
pd.Series(dict(spearman_off_prolif=d.off.corr(d.prolif, method="spearman"), spearman_off_purity=d.off.corr(d.purity, method="spearman"),
               spearman_prolif_purity=d.prolif.corr(d.purity, method="spearman"), mean_off=d.off.mean())).to_csv(OUT / "a127_offplane_descriptive.csv", header=False)
print("done")
