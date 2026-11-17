"""A302 (registered in `paper/docs/registered_experiments.md`): do the departures of meal-driven tissues share a direction
with the liver's within a person, beyond the body's time?

GTEx fast deaths, A104 periphery rows. For the liver and each tissue j: host from all periphery tissues except both;
each state regressed on it (strata intercepts); residual cross-moment minus the part induced by the common host error
(a_l conj(a_j) x mean trace of the host posterior covariance); corrected |correlation| r_j. T = median r_j (meal-driven
group G) - median r_j (other tissues N). Null (i): random groups of |G| from G + N (10,000); (ii): 100 matched sets.

    python score_a302.py [--n-sets 100] [--jobs 4] [--check]
Outputs (`../output/`): pairs.csv, sets.csv, tests.csv. Aggregates only.
"""
import argparse
import sys
import time
import warnings
from multiprocessing import Pool
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[3]))
sys.path.insert(0, str(HERE.parents[3] / "results" / "figure_2" / "scripts" / "backbone"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
import compute_oscillator_gtex as A104  # noqa: E402
from clocksync._core import host as HC  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.paired import pseudo_waveforms  # noqa: E402
from common import gtex  # noqa: E402

OUT = HERE.parents[1] / "output"
ANCHOR = "Liver"
G = ["Stomach", "Colon - Transverse", "Colon - Sigmoid", "Esophagus - Mucosa", "Adipose - Subcutaneous",
     "Adipose - Visceral (Omentum)"]
EXCLUDE = ["Whole Blood", "Testis"]
MIN_SHARED, N_DRAWS = 30, 10_000
META = ["donor", "compartment", "depth"]
S = {}


def residuals(d, k, h):
    pos = pd.Series(np.arange(len(h.resections)), index=[str(r) for r in h.resections])
    i = np.flatnonzero((d.compartment == k) & np.isin(d.resection.astype(str), pos.index))
    j = pos[d.resection[i].astype(str)].to_numpy()
    x = d.beta[i, 0] + 1j * d.beta[i, 1]
    f = h.f[j, 0] + 1j * h.f[j, 1]
    Gm = pd.get_dummies(pd.Series(d.group[i])).to_numpy(float)
    Z = np.column_stack([Gm, f])
    coef, *_ = np.linalg.lstsq(Z, x, rcond=None)
    e = x - Z @ coef
    tr = h.f_cov[j, 0, 0] + h.f_cov[j, 1, 1]
    donors = np.array([str(r).split("#")[0] for r in d.resection[i]])
    return pd.Series(e, index=donors), coef[-1], pd.Series(tr, index=donors)


def pair_r(d, comps, j):
    h = HC.host(d, [c for c in comps if c not in (ANCHOR, j)])
    el, al, tl = residuals(d, ANCHOR, h)
    ej, aj, tj = residuals(d, j, h)
    sh = el.index.intersection(ej.index)
    if len(sh) < MIN_SHARED:
        return dict(tissue=j, n=len(sh), r=np.nan, r_raw=np.nan)
    a, b, t = el[sh].to_numpy(), ej[sh].to_numpy(), tl[sh].to_numpy()
    m = np.mean(a * b.conj())
    L = al * np.conj(aj) * t.mean()
    vl = np.mean(np.abs(a) ** 2) - abs(al) ** 2 * t.mean()
    vj = np.mean(np.abs(b) ** 2) - abs(aj) ** 2 * t.mean()
    r = abs(m - L) / np.sqrt(vl * vj) if vl > 0 and vj > 0 else np.nan
    r_raw = abs(m) / np.sqrt(np.mean(np.abs(a) ** 2) * np.mean(np.abs(b) ** 2))
    return dict(tissue=j, n=len(sh), r=float(r), r_raw=float(r_raw))


def all_pairs(pb, wf, name):
    comps = sorted(pb.compartment.unique())
    ds = cs.Dataset(pb[META + list(wf.genes)], anchors=comps, waveform=wf, name=name)
    d = ds.data
    others = [c for c in comps if c != ANCHOR]
    return pd.DataFrame([pair_r(d, comps, j) for j in others])


def T_of(p):
    g = p[p.tissue.isin(G)].r.dropna()
    n = p[~p.tissue.isin(G)].r.dropna()
    return float(g.median() - n.median()) if len(g) and len(n) else np.nan


def set_task(i):
    warnings.filterwarnings("ignore")
    limit_threads(1)
    try:
        p = all_pairs(S["pb"], S["sets"][i], f"m{i}")
        return dict(set=i, T=T_of(p), median_r=float(p.r.median()))
    except Exception as e:  # noqa: BLE001
        return dict(set=i, error=f"{type(e).__name__}: {e}"[:200])


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--n-sets", type=int, default=100)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    warnings.filterwarnings("ignore")
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    import anndata as ad
    h5 = ad.read_h5ad(gtex.LABELLED, backed="r")
    genes = list(h5.var_names)
    h5.file.close()
    pb, _ = gtex.load(min_donors=A104.MIN_DONORS, genes=genes)
    limit_threads(1)
    pb = pb[~pb.compartment.str.startswith("Brain") & ~pb.compartment.isin(EXCLUDE)].reset_index(drop=True)
    wf = cs.gtex()
    p = all_pairs(pb, wf, "clock")
    p["group"] = np.where(p.tissue.isin(G), "G", "N")
    T = T_of(p)
    print(f"clock pairs done {time.time() - t0:.0f}s", flush=True)
    if a.check:
        print(p[["tissue", "n"]].to_string(index=False), f"{len(p)} pairs (r not printed)")
        return
    rng = np.random.default_rng(302)
    pool_t = p.dropna(subset=["r"]).tissue.to_numpy()
    ng = int(p[p.group == "G"].r.notna().sum())
    null = []
    for _ in range(N_DRAWS):
        gi = set(rng.choice(pool_t, ng, replace=False))
        rr = p.dropna(subset=["r"])
        null.append(rr[rr.tissue.isin(gi)].r.median() - rr[~rr.tissue.isin(gi)].r.median())
    null = np.array(null)
    p_lab = float((1 + np.sum(null >= T)) / (1 + N_DRAWS))
    p.to_csv(OUT / "pairs.csv", index=False)
    xg = pb[[c for c in pb.columns if c not in META]].to_numpy(np.float32)
    gn = [c for c in pb.columns if c not in META]
    lm = pd.Series(np.log((xg / pb.depth.to_numpy(float)[:, None]).mean(0) + 1e-12), index=gn)[xg.mean(0) >= 1]
    del xg
    sets = pseudo_waveforms(lm, a.n_sets, np.random.default_rng([302, 1]), wf, kind="matched")
    keep = sorted(set(g for w in sets for g in w.genes))
    S.update(pb=pb[META + keep].copy(), sets=sets)
    del pb
    with Pool(a.jobs) as pool:
        st = pd.DataFrame(pool.map(set_task, range(len(sets))))
    st.to_csv(OUT / "sets.csv", index=False)
    q95 = float(st["T"].quantile(0.95))
    tests = pd.DataFrame([dict(test="H302.1", T=T, p_labels=p_lab, matched_median_T=float(st["T"].median()), matched_q95_T=q95,
                               matched_n=int(st["T"].notna().sum()),
                               median_r_G=float(p[p.group == "G"].r.median()), median_r_N=float(p[p.group == "N"].r.median()),
                               passes=bool(T > 0 and p_lab < 0.05 and T > q95))])
    tests.to_csv(OUT / "tests.csv", index=False)
    print(p.round(3).to_string(index=False))
    print(tests.round(3).T.to_string())
    print(f"done {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
