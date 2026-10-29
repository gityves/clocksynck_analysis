"""A126 steps 3-4: H126.1 (luteal top vs bottom tertile, HR+ luminal coherence, n-weighted over studies, labels permuted within study),
H126.3 (BMI >= 30 vs < 30), pseudo-clock companion (200 variance-matched sets, A127 draw), Holm, MDE = 2.49 x null SD, descriptives.
python -I run_h126.py [--nperm 2000] [--npseudo 200]. Requires score_validity.py outputs (written before this runs)."""
import argparse, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
import numpy as np, pandas as pd
import a126_lib as L
import clocksync as cs
ap = argparse.ArgumentParser(); ap.add_argument("--nperm", type=int, default=2000); ap.add_argument("--npseudo", type=int, default=200); args = ap.parse_args()
OUT, G = L.OUT, L.GENES
rng = np.random.default_rng(1261)
pb = pd.read_parquet(OUT / "pb_all.parquet"); sc = pd.read_csv(OUT / "luteal_scores.csv", dtype={"donor": str})
dm = pd.read_csv(OUT / "donors_all.csv", dtype={"donor": str})
hrp = pb[(pb.comp == "hr")].set_index(["study", "donor"])
panel_cols = [c for c in pb.columns if c not in ("study", "donor", "comp", "n_cells", "depth")]
pool_all = pd.read_csv(OUT / "pseudo_pool.csv", header=None)[0].tolist()
pool = [g for g in pool_all if g in panel_cols]

def rows_of(df, comp="hr"):
    x = pb[pb.comp == comp].set_index(["study", "donor"])
    return x.reindex(pd.MultiIndex.from_frame(df[["study", "donor"]]))

# ---- pseudo-clock sets (mean and sd of log fractions over the HR+ pseudobulks, centred within study for the sd) ----------------------------------
hr_all = hrp.reset_index()
lf = np.log((hr_all[panel_cols].to_numpy(float) + .5) / hr_all.depth.to_numpy(float)[:, None])
lf = pd.DataFrame(lf, columns=panel_cols)
mean = lf.mean(); sd = (lf - lf.groupby(hr_all.study.to_numpy()).transform("mean")).std()
sets = L.pseudo_sets(mean, sd, pool, args.npseudo, np.random.default_rng(126))
pd.DataFrame(sets).to_csv(OUT / "pseudo_sets.csv", index=False)
mq = pd.DataFrame({"clock_mean": mean[G].values, "clock_sd": sd[G].values, "pseudo_mean_avg": [mean[[s[i] for s in sets]].mean() for i in range(14)],
                   "pseudo_sd_avg": [sd[[s[i] for s in sets]].mean() for i in range(14)]}, index=G)
mq.to_csv(OUT / "pseudo_match_quality.csv")
WFS = [L.wf_for(s) for s in sets]

# ---- generic two-group T ----------------------------------------------------------------------------------------------------------------
def make(df, comp="hr", genes=G):
    """df: study, donor, lab (int/NaN). Returns per-study list of (counts n x genes, depth, labels)."""
    r = rows_of(df, comp)
    ok = r.n_cells.notna().to_numpy() & (r.n_cells.fillna(0).to_numpy() >= 50)
    d, r = df[ok].reset_index(drop=True), r[ok]
    return d, r

def study_blocks(d, r, genes):
    out = []
    for s in d.study.unique():
        i = np.flatnonzero((d.study == s).to_numpy())
        out.append((s, i, r[genes].to_numpy(float)[i], r.depth.to_numpy(float)[i]))
    return out

def T(blocks, lab, hi, lo, wf, min_n=4, perm=None):
    num = den = 0.0
    for s, i, C, dep in blocks:
        l = lab[i] if perm is None else lab[i][perm[s]]
        a, b = l == hi, l == lo
        if a.sum() < min_n or b.sum() < min_n: continue
        t = L.coh_counts(C[a], dep[a], wf) - L.coh_counts(C[b], dep[b], wf)
        num += len(i) * t; den += len(i)
    return num / den if den else np.nan

def test(name, df, comp, hi, lo, tail, genes_sets=True, nperm=args.nperm, min_n=4, descriptive=False):
    d, r = make(df, comp)
    lab = d.lab.to_numpy()
    bl = study_blocks(d, r, G)
    obs = T(bl, lab, hi, lo, L.WF, min_n)
    null = np.empty(nperm)
    for k in range(nperm):
        perm = {s: rng.permutation(len(i)) for s, i, _, _ in bl}
        null[k] = T(bl, lab, hi, lo, L.WF, min_n, perm)
    null = null[np.isfinite(null)]
    p = cs.stats.perm_p(null, obs, "less", tie=0.0) if tail == "less" else (1 + np.sum(np.abs(null) >= abs(obs))) / (len(null) + 1)
    row = dict(test=name, comp=comp, T=obs, null_mean=null.mean(), null_sd=null.std(ddof=1), p=p, tail=tail, MDE=2.49 * null.std(ddof=1), n_perm=len(null),
               n_hi=int((lab == hi).sum()), n_lo=int((lab == lo).sum()), n_by_study=";".join(f"{s}:{len(i)}" for s, i, _, _ in bl))
    if genes_sets:
        tp = []
        for s_, w in zip(sets, WFS):
            blp = study_blocks(d, r, s_)
            tp.append(T(blp, lab, hi, lo, w, min_n))
        tp = np.array(tp)
        row.update(pseudo_percentile=100 * np.mean(tp <= obs), pseudo_q05=np.nanquantile(tp, .05), pseudo_q95=np.nanquantile(tp, .95), pseudo_median=np.nanmedian(tp),
                   clock_below_pseudo_q05=bool(obs < np.nanquantile(tp, .05)), clock_beyond_pseudo_5_95=bool(obs < np.nanquantile(tp, .05) or obs > np.nanquantile(tp, .95)))
        pd.DataFrame({"pseudo_T": tp}).to_csv(OUT / f"pseudoT_{name.split()[0]}_{comp}.csv", index=False)
    pd.DataFrame({"null": null}).to_csv(OUT / f"null_{name.split()[0]}_{comp}.csv", index=False)
    return row

pre = sc[sc.set == "pre"].copy(); pre["lab"] = pre.tert.astype(float)
# ---- H126.1 and H126.3 --------------------------------------------------------------------------------------------------------------------
r1 = test("H126.1 top vs bottom luteal tertile, HR+ luminal", pre, "hr", 2, 0, "less")
b = pre[pre.bmi.notna()].copy(); b["lab"] = (b.bmi >= 30).astype(float)
r3 = test("H126.3 BMI>=30 vs <30, HR+ luminal", b, "hr", 1, 0, "two")
main = pd.DataFrame([r1, r3]); main["p_holm"] = L.holm(main.p.to_numpy())
main.to_csv(OUT / "h126_main_tests.csv", index=False); print(main.round(4).T.to_string(), flush=True)

# ---- descriptives ---------------------------------------------------------------------------------------------------------------------------
desc = []
for comp in ("lumprog", "basal"):
    desc.append(test(f"desc top vs bottom luteal tertile (HR+ score), {comp}", pre, comp, 2, 0, "less", genes_sets=False, nperm=1000))
# tertile-wise coherence of HR+ (and lumprog, basal)
tw = []
for comp in ("hr", "lumprog", "basal"):
    d, r = make(pre, comp)
    for s in list(d.study.unique()) + ["pooled_groups"]:
        for t in range(3):
            m_ = ((d.lab == t) & ((d.study == s) if s != "pooled_groups" else True)).to_numpy()
            if m_.sum() >= 4:
                tw.append(dict(comp=comp, study=s, tertile=t, n=int(m_.sum()), coherence=L.coh_counts(r[G].to_numpy(float)[m_], r.depth.to_numpy(float)[m_])))
pd.DataFrame(tw).to_csv(OUT / "desc_coherence_by_tertile.csv", index=False)
# negative control: Kumar postmenopausal split by the same score
post = sc[(sc.set == "post") & (sc.study == "kumar") & (sc.tert >= 0)].copy(); post["lab"] = post.tert.astype(float)
if len(post) >= 12: desc.append(test("negctrl Kumar postmenopausal top vs bottom score tertile", post, "hr", 2, 0, "less", nperm=2000))
# Murrow without hormonal-contraception users: tertiles re-cut within Murrow non-users
hm = pre.merge(dm[["study", "donor", "hc"]].rename(columns={"hc": "hc2"}), on=["study", "donor"], how="left")
nh = hm[~((hm.study == "murrow") & (hm.hc2 != "none"))].copy()
for s, g in nh.groupby("study"): nh.loc[g.index, "lab"] = pd.qcut(g.cut.rank(method="first"), 3, labels=False)
desc.append(test("sens H126.1 Murrow hormonal-contraception users removed (tertiles re-cut)", nh, "hr", 2, 0, "less", nperm=args.nperm))
nh_n = nh[nh.study == "murrow"].shape[0]
# parity (Reed + Murrow nulliparous vs parous), BRCA (Reed)
pp = pre[pre.parity.notna() & pre.study.isin(["reed", "murrow"])].copy(); pp["lab"] = (pp.parity >= 1).astype(float)
desc.append(test("desc parity>=1 vs 0, HR+ luminal (Reed+Murrow)", pp, "hr", 1, 0, "two", genes_sets=False, nperm=1000))
br = pre[(pre.study == "reed") & pre.brca.isin(["BRCA1", "BRCA2", "WT", "assume_WT"])].copy(); br["lab"] = br.brca.isin(["BRCA1", "BRCA2"]).astype(float)
desc.append(test("desc BRCA carrier vs non-carrier, HR+ luminal (Reed, age<=45)", br, "hr", 1, 0, "two", genes_sets=False, nperm=1000))
D = pd.DataFrame(desc); D.to_csv(OUT / "h126_descriptive_tests.csv", index=False); print(D.round(4).T.to_string(), flush=True)

# ---- HR+ phase relative to the stromal host (clocksync Dataset), by luteal tertile -----------------------------------------------------
try:
    rows = []
    for s in pre.study.unique():
        lab = pre[pre.study == s].set_index("donor").tert
        for _, x in pb[(pb.study == s) & pb.donor.isin(lab.index)].iterrows():
            if x.comp in ("fibro", "endo") and x.n_cells >= 50: rows.append(dict(donor=f"{s}|{x.donor}", compartment=x.comp, stratum=s, depth=x.depth, n_cells=x.n_cells, **{g: x[g] for g in G}))
            if x.comp == "hr" and x.n_cells >= 50:
                for t in ("all", int(lab[x.donor])): rows.append(dict(donor=f"{s}|{x.donor}", compartment=f"hr|{t}", stratum=s, depth=x.depth, n_cells=x.n_cells, **{g: x[g] for g in G}))
    tab = pd.DataFrame(rows)
    ds = cs.Dataset(tab, anchors=["fibro", "endo"], name="a126")
    res = ds.fit(n_boot=0, gate=False, targets=[c for c in ds.compartments if c.startswith("hr")])
    t = res.table; t.to_csv(OUT / "desc_hr_phase_vs_host.csv"); print(t.round(3).to_string())
except Exception as e:
    (OUT / "desc_hr_phase_vs_host.txt").write_text(f"not run: {type(e).__name__}: {e}\n"); print("host fit failed:", e)
