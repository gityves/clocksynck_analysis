"""A115 null and power without healthy data. Per series the UNINVOLVED samples are split at random into two groups
(sizes proportional to the series' uninvolved and healthy counts, summing to the uninvolved count: availability cap) and the
within-group core-cue offset difference (group A minus group B) is computed; 2,000 splits, same splits pooled over series
with the A54 weights. Sensitivity: sampling with replacement at the full (uninvolved, healthy) sizes. Only uninvolved
samples are read. MDE = (1.645 + 0.842) * SD (one-sided alpha 0.05, 80 % power, normal approximation).

    python -I a115_null_power.py
Outputs: output/a115_null_series.csv, output/a115_null_pooled.csv, output/a115_min_p.csv
"""
import sys
from math import lgamma, log
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np, pandas as pd
import _lib as L

OUT = Path(__file__).resolve().parents[1] / "output"
N = 2000
Z = 1.645 + 0.842
sets = L.cue_sets(); genes = L.all_set_genes(sets); w54 = L.weights54()
# healthy counts: metadata only (output/raw_metadata_scan.csv, GSE121212 / GSE13355 / GSE201827 / GSE268958 lines)
H = {("GSE13355", "psoriasis"): 64, ("GSE121212", "psoriasis"): 38, ("GSE201827", "psoriasis"): 17, ("GSE268958", "psoriasis"): 19,
     ("GSE121212", "atopic_dermatitis"): 38}
FILE = {"GSE13355": "processed_heldout", "GSE121212": "processed", "GSE201827": "processed", "GSE268958": "processed"}


def prep(s, disease):
    d = L.read_series(L.ROOT / FILE[s] / f"{s}.parquet", genes)
    d = d[d.disease == disease]
    counts = bool((d.value_type == "counts").all())
    pat = sorted(set(d[d.status == "lesional"].patient) & set(d[d.status == "nonlesional"].patient))  # as A54: patients with a pair
    x = d[(d.status == "nonlesional") & d.patient.isin(pat)].drop_duplicates("patient")   # uninvolved only
    out = {}
    for k, wf in {"cue": sets["cue"], "systemic": sets["systemic"]}.items():
        bc, vc = L.states_for(x, L.WF, counts)[0]; bk, vk = L.states_for(x, wf, counts)[0]
        out[k] = (bc, vc, bk, vk)
    return len(x), out


def off(t, idx):
    return L.PD.S.ccorr(t[0][idx], t[2][idx], t[1][idx], t[3][idx]) if False else L.S.ccorr(t[0][idx], t[2][idx], t[1][idx], t[3][idx])


def logC(n, k):
    return (lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)) / log(10)


def run(items, tag, rng):
    """items: list of (series, disease, nU, dict). Returns per-series rows and pooled rows."""
    rows, pool = [], []
    for scheme in ("proportional_split", "with_replacement_full_size"):
        draws = {k: {it[0] + "|" + it[1]: np.empty((N, 2), complex) for it in items} for k in ("cue", "systemic")}
        for s, dis, n, st in items:
            u, h = n, H[(s, dis)]
            a = int(round(n * u / (u + h))); b = n - a
            for i in range(N):
                if scheme == "proportional_split":
                    p = rng.permutation(n); ia, ib = p[:a], p[a:a + b]
                else:
                    ia, ib = rng.integers(0, n, u), rng.integers(0, n, h)
                for k in ("cue", "systemic"):
                    draws[k][s + "|" + dis][i] = (off(st[k], ia), off(st[k], ib))
            for k in ("cue", "systemic"):
                dd = draws[k][s + "|" + dis]
                diff = np.array([L.wrap((np.angle(x[0]) - np.angle(x[1])) * L.HOURS) for x in dd])
                rows.append(dict(analysis=tag, scheme=scheme, set=k, series=s, disease=dis, n_uninvolved=n, n_healthy_available=h,
                                 size_A=a if scheme == "proportional_split" else u, size_B=b if scheme == "proportional_split" else h,
                                 null_mean_h=diff.mean(), null_sd_h=diff.std(ddof=1), mde_h=Z * diff.std(ddof=1),
                                 sd_scale_to_full=np.sqrt((1 / u + 1 / h) / (1 / a + 1 / b)) if scheme == "proportional_split" else 1.0))
        if len(items) > 1 or True:
            keys = [it[0] + "|" + it[1] for it in items]
            w = np.array([w54.get(it[0], 1.0) for it in items])
            for k in ("cue", "systemic"):
                pa = np.array([L.wrap(np.angle(np.sum(w * np.array([draws[k][key][i][0] for key in keys]))) * L.HOURS -
                                      np.angle(np.sum(w * np.array([draws[k][key][i][1] for key in keys]))) * L.HOURS) for i in range(N)])
                pool.append(dict(analysis=tag, scheme=scheme, set=k, n_series=len(items), series=";".join(keys), weights=";".join(f"{x:.1f}" for x in w),
                                 null_mean_h=pa.mean(), null_sd_h=pa.std(ddof=1), mde_h=Z * pa.std(ddof=1)))
    return rows, pool


rng = np.random.default_rng(115)
ps = [(s, "psoriasis") for s in ("GSE13355", "GSE121212", "GSE201827", "GSE268958")]
items = []
for s, dis in ps:
    n, st = prep(s, dis); items.append((s, dis, n, st)); print(s, dis, n, flush=True)
r1, p1 = run(items, "psoriasis_pooled", rng)
n, st = prep("GSE121212", "atopic_dermatitis"); print("AD", n)
r2, p2 = run([("GSE121212", "atopic_dermatitis", n, st)], "atopic_dermatitis", rng)
pd.DataFrame(r1 + r2).to_csv(OUT / "a115_null_series.csv", index=False)
pool = pd.DataFrame(p1 + p2); pool.to_csv(OUT / "a115_null_pooled.csv", index=False)
mp = []
for tag, its in (("psoriasis_pooled", items), ("atopic_dermatitis", [(("GSE121212"), "atopic_dermatitis", n, None)])):
    lg = [logC(it[2] + H[(it[0], it[1])], it[2]) for it in its]
    for n_perm in (2000, 10000):
        mp.append(dict(analysis=tag, series=";".join(f"{it[0]}(u={it[2]},h={H[(it[0], it[1])]})" for it in its),
                       log10_distinct_labelings_per_series=";".join(f"{x:.1f}" for x in lg), log10_joint=sum(lg), n_perm=n_perm,
                       min_p_monte_carlo=1 / (1 + n_perm), min_p_exact=min(1.0, 10 ** (-sum(lg)))))
pd.DataFrame(mp).to_csv(OUT / "a115_min_p.csv", index=False)
pd.set_option("display.width", 250)
print(pd.DataFrame(r1 + r2).round(3).to_string()); print(pool.round(3).to_string()); print(pd.DataFrame(mp).to_string())
