"""A126 H126.0 test on the per-cell table from h126_0_endometrium.py. T0 = coh(cycling) - coh(non-cycling) over donors with >= 30
cycling epithelial cells; null: within-donor label swaps (10,000); one-sided p (loss). Sensitivities (not registered, labelled):
depth-matched non-cycling subsample, G1 comparator, proliferative-phase donors only. Descriptive: proliferative minus secretory."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
import numpy as np, pandas as pd
import a126_lib as L
G = L.GENES
cells = pd.read_parquet(L.OUT / "endometrium_epithelial_cells_clock.parquet")
elig = pd.read_csv(L.OUT / "h126_0_eligibility_counts.csv", index_col=0)
don = elig.index[elig.n_cyc >= 30]
print("eligible", len(don), "need >= 15")
rng = np.random.default_rng(1260)

def pb(df): return df[G].sum().to_numpy(float), df.depth.sum()
def pair_table(sel_a, sel_b, donors, match=False):
    A, B, da, db = [], [], [], []
    for d in donors:
        c = cells[cells.donor == d]
        a, b = c[sel_a(c)], c[sel_b(c)]
        if match: b = b.sample(min(len(a), len(b)), random_state=int(rng.integers(1 << 30)))
        x, y = pb(a), pb(b); A.append(x[0]); da.append(x[1]); B.append(y[0]); db.append(y[1])
    return np.array(A), np.array(da), np.array(B), np.array(db)

def run(name, A, da, B, db, nperm=10000, tail="less"):
    obs = L.coh_counts(A, da) - L.coh_counts(B, db)
    n = len(A); null = np.empty(nperm)
    for k in range(nperm):
        f = rng.random(n) < .5
        A2, B2 = np.where(f[:, None], B, A), np.where(f[:, None], A, B)
        d1, d2 = np.where(f, db, da), np.where(f, da, db)
        null[k] = L.coh_counts(A2, d1) - L.coh_counts(B2, d2)
    p = (1 + np.sum(null <= obs)) / (nperm + 1)
    return dict(test=name, n_donors=n, coh_cycling=L.coh_counts(A, da), coh_other=L.coh_counts(B, db), T0=obs, null_mean=null.mean(),
                null_sd=null.std(ddof=1), p_one_lower=p, MDE=2.49 * null.std(ddof=1), n_perm=nperm)
cyc, non, g1 = (lambda c: c.cyc), (lambda c: ~c.cyc), (lambda c: c.g1)
rows = [run("H126.0 cycling vs non-cycling (registered)", *pair_table(cyc, non, don))]
rows.append(run("sens: non-cycling subsampled to cycling cell number", *pair_table(cyc, non, don, True)))
rows.append(run("sens: cycling vs G1 (A24 note b)", *pair_table(cyc, g1, don)))
dp = [d for d in don if elig.phase[d] == "proliferative"]
rows.append(run("sens: proliferative-phase donors only", *pair_table(cyc, non, dp)) if len(dp) >= 5 else {})
res = pd.DataFrame(rows); res.to_csv(L.OUT / "h126_0_tests.csv", index=False); print(res.round(4).to_string())

# descriptive: phase contrast, all epithelial cells per donor
ph = elig.phase.to_dict(); P = {}
for d, c in cells.groupby("donor"):
    P[d] = pb(c)
D = pd.DataFrame({"ph": pd.Series(ph), "n": elig.n_epi})
out = []
for a, b in (("proliferative", "secretory"), ("proliferative", "periovulatory"), ("periovulatory", "secretory")):
    ia, ib = [d for d in P if ph[d] == a], [d for d in P if ph[d] == b]
    f = lambda ids: L.coh_counts(np.array([P[d][0] for d in ids]), np.array([P[d][1] for d in ids]))
    obs = f(ia) - f(ib); allid = ia + ib; null = []
    for _ in range(2000):
        s = rng.permutation(len(allid)); null.append(f([allid[i] for i in s[:len(ia)]]) - f([allid[i] for i in s[len(ia):]]))
    null = np.array(null)
    out.append(dict(contrast=f"{a} - {b}", n_a=len(ia), n_b=len(ib), coh_a=f(ia), coh_b=f(ib), diff=obs, null_sd=null.std(ddof=1),
                    p_two=(1 + np.sum(np.abs(null) >= abs(obs))) / 2001))
pd.DataFrame(out).to_csv(L.OUT / "h126_0_endometrium_phase_descriptive.csv", index=False); print(pd.DataFrame(out).round(4).to_string())
