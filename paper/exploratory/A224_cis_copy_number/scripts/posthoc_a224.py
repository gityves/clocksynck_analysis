"""A224 POST HOC descriptives (not registered): per-cell values, leave-one-gene adjustment, CN by subtype, trans component of beta.
Reuses compute_a224.py up to the stage branch (executed in a scratch output dir, so no registered output is touched or rerun).
    python -I posthoc_a224.py [--jobs 4]"""
import os
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT_REAL = HERE.parent / "output"
src = (HERE / "compute_a224.py").read_text(); src = src[: src.index('HYPS = ["H224.1"')]
tmp = tempfile.mkdtemp(prefix="a224_posthoc_")
sys.argv = ["compute_a224.py", "--stage", "power", "--out", tmp, "--jobs", "4"]
ns = {"__file__": str(HERE / "compute_a224.py"), "__name__": "a224"}
exec(compile(src, str(HERE / "compute_a224.py"), "exec"), ns)
globals().update({k: v for k, v in ns.items() if not k.startswith("__")})
import numpy as np, pandas as pd  # noqa: E402

# ---- 1. per cell -----------------------------------------------------------------------------------------------------------------
rows = []
def cellrows(h, coh, cells, Yr, Ya):
    vr, va = CellD(cells).values(Yr, L.WF), CellD(cells).values(Ya, L.WF)
    return [dict(hypothesis=h, cohort=coh, cell=c["key"], n=len(c["lab"]), n_arm1=int((c["lab"] == 1).sum()), n_arm0=int((c["lab"] == 0).sum()), unadjusted=a, adjusted=b, diff=a - b)
            for c, a, b in zip(cells, vr, va)]
rows += cellrows("H224.1", "METABRIC", c_mb, Ymb0, Ymb) + cellrows("H224.1", "TCGA_BRCA", c_br, Y0["BRCA"], Yt["BRCA"])
for p in PROJECTS:
    if cells_p[p]:
        rows += cellrows("H224.2", p, cells_p[p], Y0[p], Yt[p])
df = pd.DataFrame(rows)
mean_rows = []
for (h, coh), g in df.groupby(["hypothesis", "cohort"], sort=False):
    mean_rows.append(dict(hypothesis=h, cohort=coh, cell="MEAN_cells", n=g.n.sum(), unadjusted=g.unadjusted.mean(), adjusted=g.adjusted.mean(), diff=g["diff"].mean()))
for h in ("H224.1", "H224.2"):
    g = df[df.hypothesis == h]
    mean_rows.append(dict(hypothesis=h, cohort="POOLED", cell="MEAN_cells", n=g.n.sum(), unadjusted=g.unadjusted.mean(), adjusted=g.adjusted.mean(), diff=g["diff"].mean()))
pd.concat([df, pd.DataFrame(mean_rows)], ignore_index=True).to_csv(OUT_REAL / "posthoc_cells.csv", index=False)
print(pd.DataFrame(mean_rows).round(4).to_string(), flush=True)

# ---- 2. leave-one-gene -----------------------------------------------------------------------------------------------------------
coh_data = {"METABRIC": (c_mb, Ymb0, Ymb), "TCGA_BRCA": (c_br, Y0["BRCA"], Yt["BRCA"])}
def D(raw, adj, mask):
    """H224.1 per-cohort and pooled D with the genes in `mask` (bool, 14) adjusted, others raw."""
    v = {}
    for k, (cells, R, A) in coh_data.items():
        v[k] = CellD(cells).values(np.where(mask, A, R), L.WF)
    return {**{k: x.mean() for k, x in v.items()}, "POOLED": np.concatenate(list(v.values())).mean()}
lo = []
for mode in ("none", "all"):
    lo.append(dict(mode=mode, gene="-", **D(None, None, np.full(14, mode == "all"))))
for j, g in enumerate(GENES):
    m = np.zeros(14, bool); m[j] = True
    lo.append(dict(mode="only_this_gene_adjusted", gene=g, **D(None, None, m)))
    lo.append(dict(mode="all_but_this_gene_adjusted", gene=g, **D(None, None, ~m)))
lo = pd.DataFrame(lo)
for k in ("METABRIC", "TCGA_BRCA", "POOLED"):
    lo[k + "_attenuation_vs_none"] = (lo[k].iloc[0] - lo[k]) / lo[k].iloc[0]
lo.to_csv(OUT_REAL / "posthoc_leave_one_gene.csv", index=False); print(lo.round(4).to_string(), flush=True)

# ---- 3. CN by subtype within the H224.1 cells ------------------------------------------------------------------------------------
cs_rows = []
for k, (cells, _, _) in coh_data.items():
    idx = np.concatenate([c["idx"] for c in cells]); lab = np.concatenate([c["lab"] for c in cells])
    ids = (mb.patient.to_numpy() if k == "METABRIC" else tc["BRCA"]["sample"].to_numpy())[idx]; CN = CNm if k == "METABRIC" else CNt
    for g in GENES:
        cn = cn_of(CN, g, ids); r = dict(cohort=k, gene=g)
        for nm, m in (("LumA", lab == 1), ("others", lab == 0)):
            x = cn[m]; x = x[np.isfinite(x)]
            r.update({f"n_{nm}": len(x), f"mean_cn_{nm}": x.mean() if len(x) else np.nan, f"var_cn_{nm}": x.var(ddof=1) if len(x) > 1 else np.nan,
                      f"frac_nonzero_{nm}": (x != 0).mean() if len(x) else np.nan})
        r["var_ratio_others_over_LumA"] = r["var_cn_others"] / r["var_cn_LumA"]; cs_rows.append(r)
pd.DataFrame(cs_rows).to_csv(OUT_REAL / "posthoc_cn_by_subtype.csv", index=False); print(pd.DataFrame(cs_rows).round(3).to_string(), flush=True)

# ---- 4. beta with trans covariates ----------------------------------------------------------------------------------------------
fga_t = pd.read_csv(EXP / "A213_heterogeneity_coherence/output/fga_tcga.csv", index_col=0).fga
fga_m = pd.read_csv(EXP / "A213_heterogeneity_coherence/output/fga_metabric.csv", index_col=0).fga
def dummies(lbl):
    c = pd.Categorical(lbl.fillna("NA")); return pd.get_dummies(c, drop_first=True).to_numpy(float)
def refit(y, cn, ph, grp, covs):
    ok = np.isfinite(y) & np.isfinite(cn) & np.isfinite(ph) & np.all(np.isfinite(covs["_all"]), axis=1)
    g = np.unique(grp[ok], return_inverse=True)[1]; base = np.column_stack([cn, np.cos(ph), np.sin(ph)])[ok]
    out = {"n": int(ok.sum())}
    for name, keys in (("base", []), ("prolif", ["prolif"]), ("prolif_sub", ["prolif", "sub"]), ("prolif_sub_fga", ["prolif", "sub", "fga"])):
        X = np.column_stack([base] + [covs[k][ok] for k in keys]); b, se, _, _ = ols(y[ok], X, g)
        out[f"beta_{name}"], out[f"se_{name}"] = b[0] / LN2, se[0] / LN2
    return out
tr = []
setups = {"TCGA": (PROJECTS, [Y0[p] for p in PROJECTS], [tc[p] for p in PROJECTS], "sample", CNt, fga_t, PH["TCGA"]),
          "METABRIC": (["METABRIC"], [Ymb0], [mb], "patient", CNm, fga_m, PH["METABRIC"])}
for coh, (names, Ys, tabs, idcol, CN, fga, ph) in setups.items():
    grp = np.concatenate([np.full(len(y), i) for i, y in enumerate(Ys)]); Y = np.vstack(Ys)
    ids = np.concatenate([t[idcol].to_numpy() for t in tabs]); prol = np.concatenate([t.prolif.to_numpy(float) for t in tabs])
    sub = pd.concat([(t["sub"].astype(object).where(t["sub"].notna(), None).map(lambda s, n=n: None if s is None else f"{n}:{s}")) for n, t in zip(names, tabs)], ignore_index=True)
    sub = sub.fillna("NA")
    fg = fga.reindex(ids).to_numpy(float)
    covs = {"prolif": prol[:, None], "sub": dummies(sub), "fga": fg[:, None]}; covs["_all"] = np.column_stack([covs["prolif"], covs["fga"]])
    for j, g in enumerate(GENES):
        cn = cn_of(CN, g, ids)
        if not np.isfinite(cn).any() or np.nanstd(cn) == 0:
            tr.append(dict(cohort=coh, gene=g)); continue
        r = refit(Y[:, j], cn, PH[coh][:, j], grp, covs)
        r.update(cohort=coh, gene=g, share_beta_kept_prolif_sub_fga=r["beta_prolif_sub_fga"] / r["beta_base"]); tr.append(r)
tr = pd.DataFrame(tr)[["cohort", "gene", "n"] + [f"{a}_{b}" for b in ("base", "prolif", "prolif_sub", "prolif_sub_fga") for a in ("beta", "se")] + ["share_beta_kept_prolif_sub_fga"]]
tr.to_csv(OUT_REAL / "posthoc_beta_trans.csv", index=False); print(tr.round(3).to_string(), flush=True)
