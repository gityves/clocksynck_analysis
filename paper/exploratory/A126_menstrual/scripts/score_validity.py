"""A126 step 2b: luteal score, tertiles, score-validity checks and feasibility table. NO clock statistic is computed here.
Score = mean z (within study) of log2 CPM of PGR, TNFSF11, WNT4, AREG in the HR+ luminal pseudobulk, donors with >= 50 HR+ cells.
If |Spearman(score, HR+ cell number)| > 0.3 (pooled over studies) tertiles are cut on the residual after regressing the score on log cell
number within study. Outputs: score_validity.md/.csv, feasibility.csv, luteal_scores.csv (donor level, public data)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
import numpy as np, pandas as pd
from scipy.stats import spearmanr, mannwhitneyu
import a126_lib as L
pb = pd.read_parquet(L.OUT / "pb_all.parquet"); dm = pd.read_csv(L.OUT / "donors_all.csv", dtype={"donor": str})
hr = pb[(pb.comp == "hr") & (pb.n_cells >= 50)].merge(dm, on=["study", "donor"])
lc = np.log2(hr[L.LUTEAL].div(hr.depth, axis=0) * 1e6 + 1)
hr = hr.assign(**{f"lc_{g}": lc[g] for g in L.LUTEAL})
def zscore(df):
    z = df[[f"lc_{g}" for g in L.LUTEAL]]
    return ((z - z.mean()) / z.std(ddof=0)).mean(axis=1)
hr["group"] = np.where(hr.premeno, "pre", np.where(hr.post, "post", "unk"))
hr = hr[hr.group.isin(["pre", "post"])].copy()
hr["score"] = np.nan
for (s, gp), idx in hr.groupby(["study", "group"]).groups.items():
    hr.loc[idx, "score"] = zscore(hr.loc[idx])
hr["logn"] = np.log(hr.n_cells.astype(float))
pre = hr[hr.group == "pre"].copy()
rho_all, _ = spearmanr(pre.score, pre.logn)
rows = [dict(check="spearman score vs log HR+ cells, pooled premenopausal", value=rho_all, n=len(pre))]
for s, g in pre.groupby("study"):
    rows.append(dict(check=f"spearman score vs log HR+ cells, {s}", value=spearmanr(g.score, g.logn)[0], n=len(g)))
    rows.append(dict(check=f"score SD / IQR / range, {s}", value=g.score.std(), n=len(g), extra=f"IQR {g.score.quantile(.75) - g.score.quantile(.25):.2f}; range {g.score.min():.2f} to {g.score.max():.2f}"))
    ct = pd.cut(g.score.rank(method="first"), 3, labels=False)
    rows.append(dict(check=f"tertile mean score low/mid/high, {s}", value=g.score[ct == 2].mean() - g.score[ct == 0].mean(), n=len(g),
                     extra=" / ".join(f"{g.score[ct == k].mean():.2f}" for k in range(3))))
resid = abs(rho_all) > 0.3
pre["cut"] = pre.score
if resid:
    for s, g in pre.groupby("study"):
        b = np.polyfit(g.logn, g.score, 1); pre.loc[g.index, "cut"] = g.score - np.polyval(b, g.logn)
rows.append(dict(check="residual rule applied (|rho| > 0.3)", value=float(resid), n=len(pre)))
# age direction
for s in ("reed", "murrow"):
    g = pre[(pre.study == s) & pre.age.notna()]
    rows.append(dict(check=f"spearman score vs age, {s} (expected < 0)", value=spearmanr(g.score, g.age)[0], n=len(g)))
rows.append(dict(check="spearman score vs age, Reed+Murrow pooled (z within study)", value=spearmanr(pre[pre.study != "kumar"].score, pre[pre.study != "kumar"].age)[0], n=int((pre.study != "kumar").sum())))
# Kumar: premenopausal vs postmenopausal HR+ luminal marker level (raw log CPM mean z over a common reference: z computed within group, so use log CPM mean)
k = hr[hr.study == "kumar"].copy()
k["raw"] = k[[f"lc_{g}" for g in L.LUTEAL]].mean(axis=1)
u = mannwhitneyu(k.raw[k.group == "pre"], k.raw[k.group == "post"], alternative="greater")
rows.append(dict(check="Kumar mean log2CPM of 4 luteal genes, pre minus post (expected > 0); one-sided MWU p", value=k.raw[k.group == "pre"].mean() - k.raw[k.group == "post"].mean(), n=len(k), extra=f"p={u.pvalue:.3g}; n pre {int((k.group == 'pre').sum())}, post {int((k.group == 'post').sum())}"))
rows.append(dict(check="Reed: mean log2CPM, age<=45 minus age>=55 (expected > 0)", value=np.nan, n=0))
r_ = hr[hr.study == "reed"].assign(raw=lambda d: d[[f"lc_{g}" for g in L.LUTEAL]].mean(axis=1))
rows[-1]["value"] = r_.raw[r_.group == "pre"].mean() - r_.raw[r_.group == "post"].mean() if (r_.group == "post").any() else np.nan
rows[-1]["n"] = len(r_); rows[-1]["extra"] = f"n pre {int((r_.group=='pre').sum())}, post (>=55) {int((r_.group=='post').sum())}"
val = pd.DataFrame(rows); val.to_csv(L.OUT / "score_validity.csv", index=False)
# tertiles within study (cut score)
pre["tert"] = -1
for s, g in pre.groupby("study"): pre.loc[g.index, "tert"] = pd.qcut(g.cut.rank(method="first"), 3, labels=False)
# negative control: Kumar post
post = hr[(hr.group == "post")].copy(); post["tert"] = -1
for s, g in post.groupby("study"):
    if len(g) >= 9: post.loc[g.index, "tert"] = pd.qcut(g.score.rank(method="first"), 3, labels=False)
out = pd.concat([pre.assign(set="pre"), post.assign(set="post")])[["study", "donor", "set", "n_cells", "score", "cut", "tert", "age", "bmi", "parity", "brca", "hc", "site"]]
out.to_csv(L.OUT / "luteal_scores.csv", index=False)
fe = pre.groupby("study").tert.value_counts().unstack(fill_value=0); fe.columns = [f"tertile_{c}" for c in fe.columns]
fe.loc["pooled"] = fe.sum()
fe["n_eligible_hr50"] = pre.groupby("study").size().reindex(fe.index).fillna(len(pre)).astype(int)
bm = pre[pre.bmi.notna()]
fe["n_with_bmi"] = bm.groupby("study").size().reindex(fe.index).fillna(len(bm)).astype(int)
fe["n_bmi_ge30"] = bm[bm.bmi >= 30].groupby("study").size().reindex(fe.index).fillna(int((bm.bmi >= 30).sum())).astype(int)
fe["n_bmi_lt30"] = bm[bm.bmi < 30].groupby("study").size().reindex(fe.index).fillna(int((bm.bmi < 30).sum())).astype(int)
fe["feasible_ge10_per_tertile"] = fe[["tertile_0", "tertile_1", "tertile_2"]].min(axis=1) >= 10
fe.to_csv(L.OUT / "feasibility.csv")
md = ["# A126 score validity and feasibility (written before any clock statistic)\n", val.round(3).to_string(index=False), "\n", fe.to_string(),
      f"\nKumar postmenopausal negative-control donors with >= 50 HR+ cells: {int((post.study == 'kumar').sum())}; Reed age >= 55: {int((post.study == 'reed').sum())}\n"]
(L.OUT / "score_validity.md").write_text("\n".join(md)); print("\n".join(md))
