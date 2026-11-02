"""A189 scoring (registered in `paper/docs/registered_experiments.md`, commit ded6066): fetal clock organisation within
one compartment (H189.1, L1) and within one organ across cell classes (H189.2 kidney, blind; H189.3 heart and lung,
A107b re-reads), Cao 2020.

H189.1: coherence R (`decompose.components`, GTEx waveform) against per-gene donor permutations (10,000), one-sided,
Holm over the 9 registered compartments; clock-specific if R > the 95th percentile of R over 200 expression-matched
pseudo-clock sets (`paired.pseudo_waveforms`) of the compartment; responsive sets reported. Leave-one-gene-out R.
H189.2/3: `Dataset.gate` (10,000 permutations), pass = p_pairing and p_acrophase < 0.05; clock-specific if p_pairing is
below the 5th percentile of the pseudo-set p_pairing (gate with each pseudo waveform, 1,000 permutations).
Pseudo sets come from the genome-wide pseudobulks (`build_genome.py`); log mean = log of the mean fraction over the
compartment's (or organ's classes') pseudobulks, genes with mean count >= 1.

    python score_a189.py [--jobs 4] [--n-perm 10000] [--n-pseudo 200] [--check]
"""
import argparse
import sys
import time
from multiprocessing import Pool
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[3]))
import anndata as ad  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.decompose import components, logfrac  # noqa: E402
from clocksync.paired import pseudo_waveforms  # noqa: E402
from common import data  # noqa: E402
from common.paths import DATA  # noqa: E402

OUT = HERE.parents[1] / "output"
GENOME = DATA / "natexp" / "cao2020" / "genome"
SEED = 189
WF = cs.gtex()
GENES = list(WF.genes)
L1 = {"Liver": "organ", "Lung": "organ", "Adrenal": "organ", "Heart": "organ", "Heart:stromal": "class",
      "Lung:stromal": "class", "Heart:endothelial": "class", "Lung:epithelial": "class", "Kidney:parenchyma": "class"}
L2 = {"H189.2 kidney": ["Kidney:parenchyma", "Kidney:stromal"],
      "H189.3 heart": ["Heart:parenchyma", "Heart:stromal", "Heart:endothelial"],
      "H189.3 lung": ["Lung:epithelial", "Lung:stromal", "Lung:endothelial"]}
MIN_CELLS = {"organ": 1000, "class": 200}
G = {}  # per-process: genome tables by level


def holm(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    adj = np.minimum(1, np.maximum.accumulate((len(p) - np.arange(len(p))) * p[o]))
    out = np.empty_like(adj)
    out[o] = adj
    return out


def genome(level):
    if level not in G:
        g = ad.read_h5ad(GENOME / f"pseudobulk_{level}.h5ad")
        g = g[g.obs.n_cells >= MIN_CELLS[level]]
        G[level] = g
    return G[level]


def table(level, comps, genes):
    """Counts of `genes` (columns) with donor, compartment, depth, stratum for the given compartments."""
    g = genome(level)
    g = g[g.obs.compartment.isin(comps)]
    X = g[:, genes].X.toarray()
    t = pd.DataFrame(X, columns=genes)
    for c in ("donor", "compartment", "n_cells", "depth"):
        t[c] = g.obs[c].to_numpy()
    return t.assign(stratum="one")


def log_mean(level, comps):
    g = genome(level)
    g = g[g.obs.compartment.isin(comps)]
    X = g.X.toarray()
    keep = X.mean(0) >= 1
    frac = X / g.obs.depth.to_numpy()[:, None]
    return pd.Series(np.log(frac.mean(0) + 1e-12), index=g.var_names)[keep]


def coherence(tab, wf):
    y, pv, k = logfrac(tab, list(wf.genes))
    return components(y, pv, k, wf)["coherence"], y, pv, k


def l1_job(args):
    comp, n_perm, n_pseudo = args
    limit_threads(1)
    level = L1[comp]
    tab = table(level, [comp], GENES)
    r, y, pv, k = coherence(tab, WF)
    rng = np.random.default_rng([SEED, 1, list(L1).index(comp)])
    null = np.array([components(np.column_stack([rng.permutation(c) for c in y.T]), pv, k, WF)["coherence"]
                     for _ in range(n_perm)])
    row = dict(compartment=comp, level=level, n_fetuses=tab.donor.nunique(), R=r,
               p=(1 + np.sum(null >= r)) / (1 + n_perm), null_mean=null.mean(), null_q95=np.quantile(null, .95))
    lm = log_mean(level, [comp])
    for kind, tag in (("matched", 2), ("responsive", 3)):
        ws = pseudo_waveforms(lm, n_pseudo, np.random.default_rng([SEED, tag, list(L1).index(comp)]), WF, kind)
        rs = np.array([coherence(table(level, [comp], list(w.genes)), w)[0] for w in ws])
        row[f"{kind}_q95"] = np.quantile(rs, .95)
        row[f"{kind}_pct"] = 100 * np.mean(rs < r)
    logo = {}
    for g in GENES:
        w = WF.drop([g])
        logo[g] = coherence(tab, w)[0]
    row["logo_min"], row["logo_max"] = min(logo.values()), max(logo.values())
    row["logo_min_gene"] = min(logo, key=logo.get)
    row["ga_days"] = ";".join(str(int(x)) for x in sorted(data.load("fetal_cao", level=level)[1]
                                                         .loc[tab.donor.unique(), "ga_days"]))
    return row


def gate(tab, wf, n_perm, seed):
    return cs.Dataset(tab, anchors=None, waveform=wf, name="a189").gate(n_perm=n_perm, seed=seed)


def l2_job(args):
    name, n_perm, n_pseudo = args
    limit_threads(1)
    comps = L2[name]
    tab = table("class", comps, GENES)
    full = tab.groupby("donor").compartment.nunique()
    tab = tab[tab.donor.isin(full.index[full == len(comps)])]
    g = gate(tab, WF, n_perm, SEED)
    row = dict(test=name, compartments="+".join(comps), n_fetuses=tab.donor.nunique(), p_pairing=g.p_pairing,
               p_acrophase=g.p_acrophase, clock_share=g.clock_share, passed=bool(g.p_pairing < .05 and
                                                                                   g.p_acrophase < .05))
    lm = log_mean("class", comps)
    for kind, tag in (("matched", 4), ("responsive", 5)):
        ws = pseudo_waveforms(lm, n_pseudo, np.random.default_rng([SEED, tag, list(L2).index(name)]), WF, kind)
        ps = []
        for j, w in enumerate(ws):
            t = table("class", comps, list(w.genes))
            ps.append(gate(t[t.donor.isin(tab.donor)], w, 1000, SEED + j).p_pairing)
        row[f"{kind}_q05_p_pairing"] = np.quantile(ps, .05)
        row[f"{kind}_pct"] = 100 * np.mean(np.array(ps) <= g.p_pairing)
    return row


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--n-perm", type=int, default=10000)
    ap.add_argument("--n-pseudo", type=int, default=200)
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    t0 = time.time()
    OUT.mkdir(exist_ok=True)
    limit_threads(1)
    n_perm, n_pseudo = (200, 5) if a.check else (a.n_perm, a.n_pseudo)
    sfx = "_check" if a.check else ""
    with Pool(a.jobs, initializer=limit_threads, initargs=(1,)) as pool:
        r2 = pool.map_async(l2_job, [(k, n_perm, n_pseudo) for k in L2])
        r1 = pool.map(l1_job, [(c, n_perm, n_pseudo) for c in L1])
        print(f"L1 done {time.time() - t0:.0f}s", flush=True)
        r2 = r2.get()
    t1 = pd.DataFrame(r1)
    t1["p_holm"] = holm(t1.p)
    t1["clock_specific"] = t1.R > t1.matched_q95
    t1["passed"] = (t1.p_holm < .05) & t1.clock_specific
    t2 = pd.DataFrame(r2)
    t2["clock_specific"] = t2.p_pairing < t2.matched_q05_p_pairing
    t1.to_csv(OUT / f"h189_1_coherence{sfx}.csv", index=False)
    t2.to_csv(OUT / f"h189_2_3_gates{sfx}.csv", index=False)
    pd.set_option("display.width", 250)
    print(t1.round(4).to_string(index=False))
    print(t2.round(4).to_string(index=False))
    print(f"runtime {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
