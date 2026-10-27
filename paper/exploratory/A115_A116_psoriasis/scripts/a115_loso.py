"""A115 post hoc sensitivity (2026-10-27, after scoring; labelled post hoc): leave one series out of the pooled D115.1.
Same estimator, weights and permutation scheme as score_a115.py (2,000 permutations, seed 1150).

    python a115_loso.py
Output: output/a115_loso.csv
"""
import numpy as np, pandas as pd
import score_a115 as A

SERIES = ("GSE13355", "GSE121212", "GSE201827", "GSE268958")
items = {s: A.prep(s, "psoriasis") for s in SERIES}
rng = np.random.default_rng(1150)
rows = []
for out in SERIES:
    keys = [s for s in SERIES if s != out]
    base = {k: (np.where(~items[k]["grp"])[0], np.where(items[k]["grp"])[0]) for k in keys}
    obs = A.stat(items, base, keys)[0]
    null = np.empty(2000)
    for i in range(2000):
        idx = {}
        for k in keys:
            p = rng.permutation(len(items[k]["grp"])); idx[k] = (p[:len(base[k][0])], p[len(base[k][0]):])
        null[i] = A.stat(items, idx, keys)[0]
    rows.append(dict(left_out=out, D115_h=obs, p_one_sided_le=(1 + np.sum(null <= obs)) / 2001, perm_null_sd=null.std(ddof=1)))
d = pd.DataFrame(rows); d.to_csv(A.OUT / "a115_loso.csv", index=False); print(d.round(3).to_string(index=False))
