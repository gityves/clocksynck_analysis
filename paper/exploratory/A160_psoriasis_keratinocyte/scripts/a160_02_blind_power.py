"""A160 step 2: feasibility counts and blind power. Uses PN and healthy libraries only (no PP statistic is computed;
PP libraries enter only the pair counts, i.e. cell numbers and the detection rule, which fix the number of pairs).

Estimand: within-group core-minus-cue offset arg E[c_core c_cue*] (one projector per set over both groups, as A115),
group A (PP) minus group B (PN). Null: PN and healthy libraries randomly assigned to two pseudo-groups of the same size
as the PP / PN groups ('boot': n libraries each drawn with replacement from the pool, the pool being smaller than 2n) or
as a disjoint random split of the pool ('split', smaller groups: its draws are rescaled to n libraries per group by
sqrt((1/n + 1/n) / (1/a + 1/b)), as A115). SD and MDE = (1.645 + 0.842) * SD. The gate uses the larger MDE of the two.
Variants: kc (keratinocyte-only) and bins (inverse-variance bootstrap weighted circular mean over the eligible A63 bins).

    python a160_02_blind_power.py [--n-null 2000]
Outputs: output/a160_feasibility.csv, a160_plan.json, a160_null_draws.csv, a160_blind_power.csv
"""
import argparse, json, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np, pandas as pd
import _lib as L

ap = argparse.ArgumentParser(); ap.add_argument("--n-null", type=int, default=2000); a = ap.parse_args()
X, obs = L.load()
st_of = lambda k: obs.loc[k, "status"]
rows, plan = [], {}
pr, pool, o = L.pairs_for(X, obs, "kc", L.MIN_LIB)
plan["kc"] = dict(pairs=pr, pool=pool, genes=L.cue_genes(X, pool))
rows.append(dict(variant="kc", unit="kc", pairs=len(pr), pool_PN=sum(st_of(k) == "nonlesional" for k in pool),
                 pool_healthy=sum(st_of(k) == "healthy" for k in pool), cue_genes_used=len(plan["kc"]["genes"]), eligible=len(pr) >= 1))
for b in L.BINS:
    pr_b, pool_b, _ = L.pairs_for(X, obs, b, L.MIN_BIN)
    ok = len(pr_b) >= L.MIN_PAIRS_BIN
    rows.append(dict(variant="bins", unit=b, pairs=len(pr_b), pool_PN=sum(st_of(k) == "nonlesional" for k in pool_b),
                     pool_healthy=sum(st_of(k) == "healthy" for k in pool_b),
                     cue_genes_used=len(L.cue_genes(X, pool_b)) if ok else np.nan, eligible=ok))
    if ok:
        plan[b] = dict(pairs=pr_b, pool=pool_b, genes=L.cue_genes(X, pool_b))
feas = pd.DataFrame(rows); feas.to_csv(L.OUT / "a160_feasibility.csv", index=False)
print(feas.to_string(index=False), flush=True)
json.dump({k: dict(pairs=v["pairs"].to_dict("list"), pool=v["pool"], genes=v["genes"]) for k, v in plan.items()},
          open(L.OUT / "a160_plan.json", "w"))  # anonymous public GEO libraries

rng = np.random.default_rng(160)
units = list(plan)
draws = []
for scheme in ("boot", "split"):
    for i in range(a.n_null):
        r = dict(scheme=scheme, draw=i)
        for u in units:
            p = plan[u]; n = len(p["pairs"])
            A, B = L.pseudo_groups(p["pool"], n, rng, scheme)
            st = L.states(X, A + B, p["genes"]); ia, ib = np.arange(len(A)), len(A) + np.arange(len(B))
            r[f"d_{u}"] = L.diff(st, ia, ib); r[f"se_{u}"] = L.boot_se(st, ia, ib, rng) if u != "kc" else np.nan
            r[f"sc_{u}"] = np.sqrt((1 / n + 1 / n) / (1 / len(A) + 1 / len(B)))
        draws.append(r)
        if i % 250 == 0:
            print(scheme, "draw", i, flush=True)
D = pd.DataFrame(draws); D.to_csv(L.OUT / "a160_null_draws.csv", index=False)
bins = [u for u in units if u != "kc"]
res = []
for scheme, g in D.groupby("scheme"):
    def scaled(u, col="d"):
        x = g[f"{col}_{u}"].to_numpy(); m = np.nanmean(x); return m + g[f"sc_{u}"].to_numpy() * (x - m) if scheme == "split" else x
    variants = {"kc": scaled("kc")}
    d = np.column_stack([scaled(u) for u in bins]); se = np.column_stack([g[f"se_{u}"].to_numpy() * (g[f"sc_{u}"].to_numpy() if scheme == "split" else 1) for u in bins])
    variants["bins"] = np.array([L.wmean(d[j][np.isfinite(se[j])], se[j][np.isfinite(se[j])]) for j in range(len(g))])
    for v, x in variants.items():
        x = x[np.isfinite(x)]
        res.append(dict(variant=v, scheme=scheme, n_draws=len(x), null_mean_h=float(x.mean()), null_sd_h=float(x.std(ddof=1)),
                        mde_h=L.Z * float(x.std(ddof=1))))
res = pd.DataFrame(res)
gate = res.groupby("variant").mde_h.max().rename("mde_gate_h")
res = res.merge(gate, on="variant"); res["confirmatory"] = res.mde_gate_h <= 1.7
res.to_csv(L.OUT / "a160_blind_power.csv", index=False); print(res.to_string(index=False))
