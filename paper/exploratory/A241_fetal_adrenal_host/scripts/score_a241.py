"""A241 scoring (registered in `paper/docs/registered_experiments.md`): do the organs of a human fetus share the clock
state of its adrenal? Cao 2020, the 7 fetuses with an adrenal pseudobulk (the A189 positive), whole organs with
>= 1,000 cells.

Host-restricted gate, `Dataset.gate(targets=organs, hosts=["Adrenal"], min_shared=3)`: only organ x adrenal
cross-covariances, pairing null = the fetuses of each organ permuted. p = the larger of p_pairing and p_acrophase.
Companions: the same gate with 200 expression-matched and 200 responsive pseudo-clock sets (1,000 permutations, same
seed), compared on `pairing_z`. Adrenal RNA in the other organs: per fetus and organ, the fraction a of adrenal marker
counts relative to the adrenal of the same fetus (an upper bound: the markers are not strictly adrenal; not valid in
the placenta, which expresses CYP11A1); the expected adrenal-derived clock-gene counts are subtracted (floored at 0,
rounded) and the gate is rerun (`<test>:corrected`).

Tests: `organs6` (primary, H241.1): lung, heart, liver, cerebellum, cerebrum, muscle. Descriptive: `organs6+placenta`,
one organ at a time (`only:<organ>`), and the lung on log fractions residualised on a linear age trend within organ
(`only:Lung:age`; adrenal 7 and lung 6 fetuses).

    python score_a241.py [--jobs 8] [--n-perm 10000] [--n-pseudo 200] [--design] [--check]

Outputs (`../output/`): design_fetal.csv, l2.csv, l2_pseudo.csv, contamination.csv, verdicts.csv (`_check` with
--check: 200 permutations, 10 sets).
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
from clocksync.paired import pseudo_waveforms  # noqa: E402
from common import data  # noqa: E402
from common import development as dev  # noqa: E402
from common.paths import DATA  # noqa: E402

OUT = HERE.parents[1] / "output"
GENOME = DATA / "natexp" / "cao2020" / "genome" / "pseudobulk_organ.h5ad"
SEED = 241
WF = cs.gtex()
GENES = list(WF.genes)
HOST, MIN_CELLS, MIN_SHARED = "Adrenal", 1000, 3
PRIMARY = ["Lung", "Heart", "Liver", "Cerebellum", "Cerebrum", "Muscle"]  # organs with an adult stand-in (power)
TESTS = {"organs6": PRIMARY, "organs6+placenta": PRIMARY + ["Placenta"],
         **{f"only:{o}": [o] for o in PRIMARY + ["Placenta"]}}
COMPANIONS = ("organs6", "organs6+placenta")  # tests scored with pseudo-clock sets
MARKERS = ["STAR", "CYP11A1", "CYP17A1", "CYP21A2", "CYP11B1", "NR5A1"]  # adrenal cortex (A190 post hoc)
FAMILIES = {"H241.1": ["organs6"]}
G = {}  # genome-wide pseudobulks of the adrenal fetuses (inherited by the workers)


def load():
    """Genome-wide organ pseudobulks (>= MIN_CELLS cells) of the fetuses that have an adrenal one; age in days."""
    g = ad.read_h5ad(GENOME)
    g = g[g.obs.n_cells >= MIN_CELLS]
    g = g[g.obs.donor.isin(g.obs.donor[g.obs.compartment == HOST])].copy()
    g.obs = g.obs.astype({"donor": str, "compartment": str})
    g.obs["age"] = data.load("fetal_cao", level="organ")[1].ga_days.reindex(g.obs.donor).to_numpy(float)
    G["g"] = g


def table(targets, genes):
    """Counts of `genes` for the adrenal and the target organs: donor, compartment, depth, age."""
    g = G["g"][G["g"].obs.compartment.isin([HOST, *targets])]
    t = g.obs[["donor", "compartment", "n_cells", "depth", "age"]].reset_index(drop=True)
    t[list(genes)] = g[:, list(genes)].X.toarray()
    return t


def gate(tab, wf, targets, n_perm, seed=SEED):
    """Adrenal-host gate over the target organs with >= MIN_SHARED adrenal fetuses; None when not estimable."""
    tg = [o for o in targets if (tab.compartment == o).sum() >= MIN_SHARED]
    try:
        return cs.Dataset(tab[tab.compartment.isin([HOST, *tg])], anchors=None, waveform=wf, min_donors=MIN_SHARED,
                          name="a241").gate(n_perm=n_perm, seed=seed, targets=tg, hosts=[HOST],
                                            min_shared=MIN_SHARED)
    except ValueError:
        return None


def log_mean(targets):
    g = G["g"][G["g"].obs.compartment.isin([HOST, *targets])]
    return dev.log_mean(pd.DataFrame(g.X.toarray(), columns=g.var_names), g.obs.depth)


def contamination(targets):
    """Adrenal RNA in the target organs: (summary per organ, table with the expected adrenal counts removed)."""
    t = table(targets, GENES + MARKERS)
    host = t[t.compartment == HOST].set_index("donor")
    rows, out = [], [t[t.compartment == HOST]]
    for o in targets:
        x = t[t.compartment == o].set_index("donor")
        h = host.loc[x.index]
        a = (x[MARKERS].sum(1) / x.depth) / (h[MARKERS].sum(1) / h.depth)
        exp = a.to_numpy()[:, None] * x.depth.to_numpy()[:, None] * (h[GENES] / h.depth.to_numpy()[:, None]).to_numpy()
        share = pd.Series(exp.sum(0) / x[GENES].sum(0).to_numpy(), index=GENES)
        rows.append(dict(organ=o, n_fetuses=len(x), a_median=a.median(), a_max=a.max(),
                         clock_share_median=share.median(), clock_share_max=share.max(),
                         clock_share_max_gene=share.idxmax()))
        x = x.copy()
        x[GENES] = np.round(np.clip(x[GENES].to_numpy(float) - exp, 0, None))
        out.append(x.reset_index())
    return pd.DataFrame(rows), pd.concat(out, ignore_index=True)[t.columns]


def age_adjusted(targets, genes):
    """Table with log fractions residualised on a linear age trend within each organ (counts at a common depth)."""
    return dev.adjust(table(targets, genes), list(genes), lambda g: dev.design(g, "age", 1))[0]


def _task(t):
    kind, name = t[:2]
    if kind == "gate":
        return t, gate(table(TESTS[name], GENES), WF, TESTS[name], G["n_perm"])
    if kind == "corrected":
        return t, gate(contamination(TESTS[name])[1], WF, TESTS[name], G["n_perm"])
    if kind == "age":
        return t, gate(age_adjusted(["Lung"], GENES), WF, ["Lung"], G["n_perm"])
    w = G["ws"][name][t[2]][t[3]]
    g = gate(table(TESTS[name], w.genes), w, TESTS[name], min(G["n_perm"], 1000))
    return t, dict(pairing_z=np.nan, p_pairing=np.nan) if g is None else dict(pairing_z=g.pairing_z,
                                                                             p_pairing=g.p_pairing)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--n-perm", type=int, default=10000)
    ap.add_argument("--n-pseudo", type=int, default=200)
    ap.add_argument("--design", action="store_true", help="write the design table only (no clock statistic)")
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    limit_threads(1)
    OUT.mkdir(exist_ok=True)
    pd.set_option("display.width", 250)
    load()
    o = G["g"].obs
    design = o.groupby("compartment").agg(n_fetuses=("donor", "nunique"), median_cells=("n_cells", "median"),
                                          median_depth=("depth", "median"), age_min=("age", "min"),
                                          age_max=("age", "max")).reset_index()
    design["primary"] = design.compartment.isin(PRIMARY)
    design.to_csv(OUT / "design_fetal.csv", index=False)
    print(design.to_string(index=False), flush=True)
    if a.design:
        return
    sfx = "_check" if a.check else ""
    G["n_perm"], n_pseudo = (200, 10) if a.check else (a.n_perm, a.n_pseudo)
    G["ws"] = {name: {kind: pseudo_waveforms(log_mean(TESTS[name]), n_pseudo, np.random.default_rng([SEED, j, i]), WF,
                                             kind) for j, kind in enumerate(("matched", "responsive"), 1)}
               for i, name in enumerate(COMPANIONS)}
    tasks = ([("gate", name) for name in TESTS] + [("corrected", name) for name in TESTS if name != "only:Placenta"]
             + [("age", "only:Lung:age")]
             + [("pseudo", name, kind, i) for name in COMPANIONS for kind in G["ws"][name] for i in range(n_pseudo)])
    t0, res = time.time(), {}
    with Pool(a.jobs, initializer=limit_threads, initargs=(1,)) as pool:
        for k, (t, r) in enumerate(pool.imap_unordered(_task, tasks)):
            res[t] = r
            if (k + 1) % 50 == 0 or k + 1 == len(tasks):
                print(f"{k + 1}/{len(tasks)} {time.time() - t0:.0f}s", flush=True)
    pseudo = pd.DataFrame([dict(test=t[1], kind=t[2], set=t[3], **r) for t, r in res.items() if t[0] == "pseudo"])
    pseudo = pseudo.sort_values(["test", "kind", "set"])
    rows = []
    for t, g in res.items():
        if t[0] == "pseudo":
            continue
        comp = {kind: pseudo[(pseudo.test == t[1]) & (pseudo.kind == kind)] for kind in ("matched", "responsive")}
        name = t[1] if t[0] in ("gate", "age") else f"{t[1]}:corrected"
        rows.append(dict(test=name, **dev.l2_row(g, comp if t[0] == "gate" and t[1] in COMPANIONS else {})))
    l2 = pd.DataFrame(rows)
    l2["clock_specific"] = l2.clock_specific.where(l2.test.isin(COMPANIONS))  # no companions: not assessed
    cont = contamination(TESTS["organs6+placenta"])[0]
    for t, name in ((l2, "l2"), (pseudo, "l2_pseudo"), (cont, "contamination")):
        t.to_csv(OUT / f"{name}{sfx}.csv", index=False)
    print(l2.round(4).to_string(index=False), cont.round(4).to_string(index=False), sep="\n", flush=True)
    v = dev.verdicts(pd.DataFrame(), l2, FAMILIES)
    c = l2.set_index("test").loc["organs6:corrected"]
    v["p_corrected"] = max(c.p_pairing, c.p_acrophase) if c.estimable else np.nan
    v.loc[(v.verdict == "PASS") & ~(v.p_corrected < .05), "verdict"] = "not shown (adrenal RNA)"
    v.to_csv(OUT / f"verdicts{sfx}.csv", index=False)
    print(v.round(4).to_string(index=False), f"runtime {time.time() - t0:.0f}s", sep="\n", flush=True)


if __name__ == "__main__":
    main()
