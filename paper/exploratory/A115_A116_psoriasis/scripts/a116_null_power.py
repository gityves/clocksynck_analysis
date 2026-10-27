"""A116 feasibility: patient counts per arm (metadata only) and the null of the A56 treated-minus-placebo change of the
lesional offset under RANDOM arm labels (shuffled on the pooled patients; the true labels are never applied to any clock
statistic). Only the null SD and quantiles are written. A56 estimator: per series and visit one projector for lesional +
uninvolved, ccorr(lesional, uninvolved) of the patient subset, pooled with weights 1 / se^2 of the baseline variance ratio
(all patients, patient bootstrap 500, arm-agnostic); change = baseline minus post; statistic = change(group A) - change(group B).
Secondary (PASI75 within treated): random labels with the responder / non-responder counts (counted, never linked to offsets).

    python -I a116_null_power.py
Outputs: output/a116_counts.csv, output/a116_null.csv, output/a116_min_p.csv
"""
import sys
from math import lgamma, log
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np, pandas as pd
import _lib as L

OUT = Path(__file__).resolve().parents[1] / "output"
N, Z = 2000, 1.645 + 0.842
PD = L.PD
FILES = {"GSE117468": "GSE117468.parquet", "GSE201827": "GSE201827.parquet"}


def load(f):
    d = L.read_series(L.ROOT / "processed_treatment" / f, L.CORE)
    counts = bool((d.value_type == "counts").all())
    t = {(v, s): d[(d.visit == v) & (d.status == s)].drop_duplicates("patient").set_index("patient") for v in ("baseline", "post") for s in ("nonlesional", "lesional")}
    pts = sorted(set.intersection(*[set(x.index) for x in t.values()]))
    return {k: x.loc[pts] for k, x in t.items()}, counts


def logC(n, k): return (lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)) / log(10)


rng = np.random.default_rng(116)
data, cnt = [], []
for s, f in FILES.items():
    t, counts = load(f)
    arm = t[("baseline", "lesional")].arm.astype(str)
    pl = arm.str.lower().str.contains("placebo").to_numpy()
    drug = arm.str.replace(r"\s+\d+.*", "", regex=True)
    wk = t[("post", "lesional")].weeks.dropna().unique() if "weeks" in t[("post", "lesional")] else []
    cnt.append(dict(series=s, visit_weeks=";".join(map(str, wk)), n_patients_4_samples=len(pl), n_placebo=int(pl.sum()), n_treated=int((~pl).sum()),
                    arms=";".join(f"{k}:{v}" for k, v in arm.value_counts().items()), pasi_available=bool(t[("post", "lesional")].pasi.notna().all()),
                    pasi_nonnull_patients=int(t[("post", "lesional")].pasi.notna().sum())))
    L2 = {k: PD.cohort_states(t[(v, "lesional")], t[(v, "nonlesional")], counts, L.WF) for k, v in (("b", "baseline"), ("p", "post"))}
    se = PD.boot_se(*L2["b"], n_boot=500, rng=rng)
    data.append(dict(series=s, b=L2["b"], p=L2["p"], w=1 / se ** 2, n=len(pl), n_pl=int(pl.sum()), pasi_b=t[("baseline", "lesional")].pasi.to_numpy(float),
                     pasi_p=t[("post", "lesional")].pasi.to_numpy(float), pl=pl))
cdf = pd.DataFrame(cnt)
# PASI75 counts within treated (counting only; no clock statistic)
for d, c in zip(data, cnt):
    imp = 1 - d["pasi_p"] / d["pasi_b"]
    ok = ~d["pl"] & np.isfinite(imp)
    c["treated_with_pasi_both_visits"] = int(ok.sum()); c["treated_pasi75"] = int((ok & (imp >= 0.75)).sum()); c["treated_not_pasi75"] = int((ok & (imp < 0.75)).sum())
    c["placebo_pasi75"] = int((d["pl"] & np.isfinite(imp) & (imp >= 0.75)).sum())
cdf = pd.DataFrame(cnt); cdf.to_csv(OUT / "a116_counts.csv", index=False)


def sub(b, m): return tuple(x[m] for x in b)


def change(ds, sel_per_series, w):
    ob = PD.pooled_offset([PD.ccorr(*sub(d["b"], m)) for d, m in zip(ds, sel_per_series)], w)
    op = PD.pooled_offset([PD.ccorr(*sub(d["p"], m)) for d, m in zip(ds, sel_per_series)], w)
    return L.wrap(ob - op)


def null(series_idx, sizes, tag, subset_masks=None):
    """random labels: in each series, `sizes[i]` of the eligible patients get label A (rest B)."""
    ds = [data[i] for i in series_idx]; w = np.array([d["w"] for d in ds])
    elig = [np.ones(d["n"], bool) if subset_masks is None else subset_masks[i] for i, d in zip(series_idx, ds)]
    out = np.empty(N)
    for r in range(N):
        A = []
        for d, e, k in zip(ds, elig, sizes):
            lab = np.zeros(d["n"], bool); idx = np.flatnonzero(e); lab[rng.permutation(idx)[:k]] = True; A.append(lab)
        selB = [e & ~a for e, a in zip(elig, A)]
        out[r] = L.wrap(change(ds, selB, w) - change(ds, A, w))  # group B (the larger group) minus group A
    return out


rows = []
def record(tag, series_idx, sizes, **kw):
    x = null(series_idx, sizes, tag, **kw)
    rows.append(dict(analysis=tag, series=";".join(data[i]["series"] for i in series_idx), group_sizes=";".join(map(str, sizes)),
                     other_group_sizes=";".join(str(int((np.ones(data[i]["n"], bool) if kw.get("subset_masks") is None else kw["subset_masks"][i]).sum()) - k) for i, k in zip(series_idx, sizes)),
                     n_perm=N, null_mean_h=x.mean(), null_sd_h=x.std(ddof=1), null_q05_h=np.quantile(x, .05), null_q95_h=np.quantile(x, .95),
                     mde_h=Z * x.std(ddof=1), min_p_monte_carlo=1 / (1 + N),
                     log10_distinct_labelings=sum(logC(int((np.ones(data[i]["n"], bool) if kw.get("subset_masks") is None else kw["subset_masks"][i]).sum()), k) for i, k in zip(series_idx, sizes))))
    print(rows[-1], flush=True)

# primary: pooled, group A = placebo-sized group (statistic = change(B) - change(A), B treated-sized: sign as H116.1 treated - placebo)
record("H116.1 pooled (GSE117468 + GSE201827)", [0, 1], [d["n_pl"] for d in data])
record("GSE117468 alone", [0], [data[0]["n_pl"]])
record("GSE201827 alone", [1], [data[1]["n_pl"]])
# secondary: PASI75 within treated (GSE117468 only: PASI not in GSE201827)
c0 = cdf.iloc[0]
masks = {0: ~data[0]["pl"] & np.isfinite(1 - data[0]["pasi_p"] / data[0]["pasi_b"])}
record("H116.2 PASI75 within treated (GSE117468)", [0], [int(c0.treated_pasi75)], subset_masks=masks)
pd.DataFrame(rows).to_csv(OUT / "a116_null.csv", index=False)
pd.set_option("display.width", 250); print(cdf.T.to_string()); print(pd.DataFrame(rows).round(3).T.to_string())
