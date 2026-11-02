"""A190 scoring (registered in `paper/docs/registered_experiments.md`, commit 0e88acc): is the fetal adrenal clock shared
across its cell types? Cao 2020, genome-wide class pseudobulks from A189 (`build_genome.py`).

H190.1: `Dataset.gate` (GTEx waveform, 10,000 permutations) on Adrenal:parenchyma, stromal and endothelial in the fetuses
with all three; pass = p_pairing and p_acrophase < 0.05; clock-specific if p_pairing is below the 5th percentile of the
pseudo-set p_pairing (200 expression-matched sets, gate with each, 1,000 permutations); responsive sets descriptive.
Descriptive: L1 coherence R per class (clock genes below the mean-count floor of `log_mean` are added back so that the
matched sets can be drawn in the shallow classes) against per-gene donor permutations (10,000), matched pseudo 95th percentile.

    python score_a190.py [--jobs 4] [--n-perm 10000] [--n-pseudo 200] [--check]
"""
import argparse
import sys
import time
from multiprocessing import Pool
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path[:0] = [str(HERE.parents[3]), str(HERE.parents[2] / "A189_fetal_within" / "scripts")]
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from clocksync.config import limit_threads  # noqa: E402
from clocksync.decompose import components  # noqa: E402
from clocksync.paired import pseudo_waveforms  # noqa: E402
from score_a189 import GENES, WF, coherence, gate, genome, log_mean, table  # noqa: E402

OUT = HERE.parents[1] / "output"
SEED = 190
COMPS = ["Adrenal:parenchyma", "Adrenal:stromal", "Adrenal:endothelial"]


def fetuses():
    t = table("class", COMPS, GENES)
    full = t.groupby("donor").compartment.nunique()
    return set(full.index[full == len(COMPS)])


def log_mean_clock(comps):
    """`log_mean` (genes with mean count >= 1) plus any clock gene below that floor (shallow classes), so the
    matched sets can be drawn for every clock gene."""
    lm = log_mean("class", comps)
    miss = [x for x in GENES if x not in lm.index]
    if miss:
        g = genome("class")
        g = g[g.obs.compartment.isin(comps)]
        frac = g[:, miss].X.toarray() / g.obs.depth.to_numpy()[:, None]
        lm = pd.concat([lm, pd.Series(np.log(frac.mean(0) + 1e-12), index=miss)])
    return lm


def pseudo_gate(args):
    w, donors, j = args
    t = table("class", COMPS, list(w.genes))
    return gate(t[t.donor.isin(donors)], w, 1000, SEED + j).p_pairing


def l1(args):
    comp, donors, n_perm, n_pseudo, i = args
    tab = table("class", [comp], GENES)
    tab = tab[tab.donor.isin(donors)]
    r, y, pv, k = coherence(tab, WF)
    rng = np.random.default_rng([SEED, 1, i])
    null = np.array([components(np.column_stack([rng.permutation(c) for c in y.T]), pv, k, WF)["coherence"]
                     for _ in range(n_perm)])
    ws = pseudo_waveforms(log_mean_clock([comp]), n_pseudo, np.random.default_rng([SEED, 2, i]), WF, "matched")
    rs = []
    for w in ws:
        t = table("class", [comp], list(w.genes))
        rs.append(coherence(t[t.donor.isin(donors)], w)[0])
    return dict(compartment=comp, n_fetuses=tab.donor.nunique(), median_depth=tab.depth.median(), R=r,
                p=(1 + np.sum(null >= r)) / (1 + n_perm), null_mean=null.mean(), null_q95=np.quantile(null, .95),
                matched_q95=np.quantile(rs, .95), matched_pct=100 * np.mean(np.array(rs) < r))


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
    donors = fetuses()
    tab = table("class", COMPS, GENES)
    tab = tab[tab.donor.isin(donors)]
    lm = log_mean_clock(COMPS)
    with Pool(a.jobs, initializer=limit_threads, initargs=(1,)) as pool:
        r1 = pool.map_async(l1, [(c, donors, n_perm, n_pseudo, i) for i, c in enumerate(COMPS)])
        g = gate(tab, WF, n_perm, SEED)
        row = dict(test="H190.1 adrenal", compartments="+".join(COMPS), n_fetuses=len(donors), p_pairing=g.p_pairing,
                   p_acrophase=g.p_acrophase, clock_share=g.clock_share,
                   passed=bool(g.p_pairing < .05 and g.p_acrophase < .05))
        print(f"gate done {time.time() - t0:.0f}s: {row}", flush=True)
        for kind, tag in (("matched", 4), ("responsive", 5)):
            ws = pseudo_waveforms(lm, n_pseudo, np.random.default_rng([SEED, tag]), WF, kind)
            ps = np.array(pool.map(pseudo_gate, [(w, donors, j) for j, w in enumerate(ws)]))
            row[f"{kind}_q05_p_pairing"] = np.quantile(ps, .05)
            row[f"{kind}_pct"] = 100 * np.mean(ps <= g.p_pairing)
            print(f"{kind} done {time.time() - t0:.0f}s", flush=True)
        t1 = pd.DataFrame(r1.get())
    t2 = pd.DataFrame([row])
    t2["clock_specific"] = t2.p_pairing < t2.matched_q05_p_pairing
    t2.to_csv(OUT / f"h190_1_gate{sfx}.csv", index=False)
    t1.to_csv(OUT / f"l1_classes{sfx}.csv", index=False)
    pd.set_option("display.width", 250)
    print(t2.round(4).to_string(index=False))
    print(t1.round(4).to_string(index=False))
    print(f"runtime {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
