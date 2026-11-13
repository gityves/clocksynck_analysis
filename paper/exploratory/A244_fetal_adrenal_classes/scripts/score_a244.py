"""A244 scoring (registered in `paper/docs/registered_experiments.md`, with its design addendum): is the clock-shaped
covariation of the fetal adrenal carried by the cortex, by the medulla, or by both?

One pseudobulk per fetus and cell class (`prep_a244.py`: marker-defined cortex, medulla, stroma; >= 200 cells). L1 signed
coherence C across fetuses (`common.development`), pooled over cohorts with >= 2 fetuses: cohort indicator and a common
quadratic age trend, per-gene permutations within cohort; matched and responsive companions. As the `all` readout of A242.

Tests
- `cortex`, `medulla` (H244.1, Holm over the two).
- `<class>:corrected`: the expected clock-gene counts of the other class subtracted (`prep_a244.foreign`).
- descriptive: `<class>:corrected:ref200` (only fetuses whose reference pseudobulk has >= 200 cells), `stroma`, the windows
  `<class>:late` (>= 12 weeks, linear age) and `<class>:early` (< 12 weeks, quadratic), `<class>:-<cohort>` (leave one
  cohort out).

    python score_a244.py [--jobs 8] [--n-perm 10000] [--n-pseudo 200] [--design] [--check]

Outputs (`../output/`): design.csv, l1.csv, verdicts.csv; then `python score_a244.py --post` writes verdicts_a244.csv (the
registered PASS rule with the corrected readout) and position.csv (each C among the adult subsamples at the full clock).
Not run (registered descriptive): the depth-matched comparison of the two classes.
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[3]))
sys.path.insert(0, str(HERE.parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from common import development as dev  # noqa: E402
import prep_a244 as PR  # noqa: E402

OUT = HERE.parents[1] / "output"
SEED = 244
WF = cs.gtex()
LATE, MIN_FETUSES = 12.0, 2
FAMILIES = {"H244.1": ["cortex", "medulla"]}


def load():
    genes = list(WF.genes)
    pb = pd.read_parquet(PR.PB)
    sheet = pd.read_csv(PR.SHEET)
    sizes = pd.read_csv(OUT / "class_pseudobulk_all_sizes.csv")
    pb.attrs["cells"] = {(f, k): n for f, k, n in zip(sizes.fetus, sizes.cls, sizes.cells)}
    f, corr = PR.foreign(pb, sheet, genes)
    obs = pd.DataFrame(dict(donor=sheet.donor.to_numpy(), compartment=sheet.compartment.to_numpy(), cohort=sheet.cohort.to_numpy(),
                            block=sheet.cohort.to_numpy(), age=sheet.age.to_numpy(float), depth=sheet.depth.to_numpy(float), stratum="all"),
                       index=sheet.id.to_numpy())
    counts = pb.loc[obs.index].copy()
    ref200 = set(f.id[f.ref_ge_min])
    extra_c, extra_o = [], []
    for k in ("cortex", "medulla"):
        ids = obs.index[obs.compartment == k]
        c = counts.loc[ids].copy()
        c.loc[:, genes] = corr.loc[ids, genes].to_numpy()
        for tag, keep in ((":corrected", ids), (":corrected:ref200", [i for i in ids if i in ref200])):
            cc = c.loc[list(keep)].copy()
            cc.index = [f"{i}#{tag}" for i in keep]
            oo = obs.loc[list(keep)].copy()
            oo.index, oo["compartment"] = cc.index, k + tag
            extra_c.append(cc)
            extra_o.append(oo)
    return pd.concat([counts, *extra_c]), pd.concat([obs, *extra_o])


def cov(df: int):
    return lambda g: dev.design(g, "age" if df else None, df, ("block",))


def pooled(obs, rows, without=()):
    rows = rows & ~obs.cohort.isin(list(without))
    n = obs[rows].groupby("cohort").donor.transform("size").reindex(obs.index, fill_value=0)
    return rows & (n >= MIN_FETUSES)


def tests(obs):
    l1 = {}
    for k in ("cortex", "medulla"):
        cls = obs.compartment == k
        l1[k] = (pooled(obs, cls), cov(2), "block")
        for tag in (":corrected", ":corrected:ref200"):
            l1[k + tag] = (pooled(obs, obs.compartment == k + tag), cov(2), "block")
        l1[k + ":late"] = (pooled(obs, cls & (obs.age >= LATE)), cov(1), "block")
        l1[k + ":early"] = (pooled(obs, cls & (obs.age < LATE)), cov(2), "block")
        for c in sorted(obs.cohort[l1[k][0]].unique()):
            l1[f"{k}:-{c}"] = (pooled(obs, cls, [c]), cov(2), "block")
    l1["stroma"] = (pooled(obs, obs.compartment == "stroma"), cov(2), "block")
    return l1, {}


def post():
    l1 = pd.read_csv(OUT / "l1.csv").set_index("test")
    v = pd.read_csv(OUT / "verdicts.csv").set_index("test")
    rows = []
    for k in ("cortex", "medulla"):
        r, pc = v.loc[k], l1.loc[f"{k}:corrected"]
        if r.p_holm < 0.05 and r.verdict == "PASS":
            verdict = "PASS" if pc.p < 0.05 else "not shown (RNA of the other class)"
        elif r.p_holm < 0.05:
            verdict = "not clock-specific"
        elif r.p < 0.05:
            verdict = "suggestive, not established"
        else:
            verdict = "fail"
        rows.append(dict(test=k, n=r.n, C=r.stat, p=r.p, p_holm=r.p_holm, matched_pct=r.matched_pct, responsive_pct=r.responsive_pct,
                         C_corrected=pc.C, p_corrected=pc.p, verdict=verdict))
    pd.DataFrame(rows).to_csv(OUT / "verdicts_a244.csv", index=False)
    sub = pd.read_csv(OUT / "feasibility_subsamples.csv")
    pos = []
    for t in [x for x in l1.index if x.split(":")[0] in ("cortex", "medulla", "stroma") and x in set(sub.test)]:
        a = sub[(sub.test == t) & (sub.setting == "a1")].C
        pos.append(dict(test=t, C=l1.loc[t, "C"], adult_median=a.median(), percentile_among_adult=100 * float(np.mean(a < l1.loc[t, "C"]))))
    pd.DataFrame(pos).to_csv(OUT / "position.csv", index=False)
    pd.set_option("display.width", 250)
    print(pd.DataFrame(rows).round(4).to_string(index=False), pd.DataFrame(pos).round(3).to_string(index=False), sep="\n")


if __name__ == "__main__":
    if "--post" in sys.argv:
        post()
    else:
        dev.main(__doc__, load, tests, WF, SEED, OUT, FAMILIES)
