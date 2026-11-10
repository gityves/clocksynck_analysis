"""Feasibility of the fetal gonad design (Garcia-Alonso 2022) before any registration: adult nuclei cut to it.

No fetal expression value is read: only the pseudobulk sheet (donor, cell class, cells, depth, age, sex, assay) built
by `prep_fetal_gonad/scripts/prep_gonad.py`. Adult single-nucleus atlases stand for the fetal gonad: every fetal donor
is replaced by one adult donor having its cell classes, each pseudobulk thinned to the depth of the fetal one, and the
share of subsamples where the planned readout passes is its sensitivity for an adult-like clock (feasible if >= 0.8).

Readouts, as they would be registered (same functions as `feasibility_fetal.py`):
1. L1 per cell class: signed coherence C across donors, each gene residualised on age (quadratic), sex and assay, against
   per-gene permutations of the residuals; pass = one-sided p < 0.05.
2. L2 across the cell classes of one gonad: `Dataset.gate` on the residualised log fractions, free host, strata =
   sex x (before | from 13 weeks); pass = p_pairing < 0.05 and p_acrophase < 0.05.
The adult donor receives the age, sex and assay of the fetal donor it stands for (`unadjusted` rows: no adjustment).

Proxies. No adult gonad atlas is in the project: breast (10x 3' v3 stratum) and heart (Reichart 2022 stratum) stand
in, with fibroblasts for the gonadal stroma, the parenchymal class for the supporting cells and endothelium for
endothelium. The fetal surface epithelium has no adult counterpart: mural cells take its place in `classes4`, which
fixes the number of compartments and their depth, not the biology (`classes3` leaves it out). The sampling times of
the adult donors are unknown.

`--sex-trend`: the age trend is fitted per sex (testis and ovary diverge from 6-8 weeks), as registered in A238; the
first run (`passrates_gonad.csv`) used one trend for both sexes.

    python feasibility_gonad.py [--n-sub 200] [--jobs 4] [--sex-trend] [--subset all] [--check]
"""
import argparse
import sys
import time
from multiprocessing import Pool
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[4]))
sys.path.insert(0, str(HERE.parent))
import anndata as ad  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from clocksync.config import limit_threads  # noqa: E402
from common import data  # noqa: E402
from feasibility_fetal import DEPTH, GENES, NATEXP, OUT, POOL_FACTOR, SEED, basis, draw, gate, l1, logfrac, residual, run  # noqa: E402

MIN_CELLS, MIN_N, SPLIT = 200, 10, 13  # cells per fetal pseudobulk; fetal donors per job; weeks post fertilisation
CLASSES = ["mesenchymal_stromal", "supporting", "endothelial", "epithelial"]
PROXY = {"breast": ("organ_breast", "10x 3' v3", ["breast:fibroblast_stromal", "breast:epithelial",
                                                   "breast:endothelial", "breast:mural_smooth_muscle"]),
         "heart": ("heart_dcm_healthy", "reichart2022_dcm|10x 3' v3", ["heart:fibroblast", "heart:cardiomyocyte",
                                                                        "heart:endothelial", "heart:mural"])}
SUBSETS = {"all": lambda f: f.wpf > 0, "female": lambda f: f.sex == "female", "male": lambda f: f.sex == "male",
           "13+wpf": lambda f: f.wpf >= SPLIT}
JOBS, POOLS = [], {}  # (kind, label, proxy, fetal rows of the job); proxy -> {donor: table}
SEX_TREND = False  # `--sex-trend`


def fetal_sheet():
    """One row per fetal donor x class with >= MIN_CELLS cells (gonad libraries; metadata only)."""
    f = pd.concat(ad.read_h5ad(NATEXP / f"garcia_alonso2022/garcia_{s}_pseudobulk_class.h5ad", backed="r").obs
                  for s in ("female", "male"))
    f = f[f.gonad & (f.n_cells >= MIN_CELLS) & f.compartment.isin(CLASSES)]
    return f.astype({"donor": str, "compartment": str, "sex": str, "assay": str}).reset_index(drop=True)


def covariates(f, adjust):
    """Per-donor design: quadratic in age (one trend per sex with `--sex-trend`), sex, assay (constant columns
    dropped); empty when unadjusted."""
    d = f.drop_duplicates("donor").set_index("donor")
    age = basis(d.wpf, 2)
    if SEX_TREND:
        age = np.column_stack([age * (d.sex == s).to_numpy()[:, None] for s in ("female", "male")])
    Z = np.column_stack([age, d.sex == "female", d.assay == d.assay.iloc[0]]).astype(float)
    Z = Z[:, np.ptp(Z, 0) > 0] if adjust else Z[:, :0]
    return pd.DataFrame(Z, index=d.index)


def thin(tab, target, rng):
    """Binomial thinning of every pseudobulk to its target depth (never deepened)."""
    t = tab.copy()
    p = np.minimum(1.0, np.asarray(target, float) / t.depth.to_numpy(float))
    t[GENES] = rng.binomial(t[GENES].to_numpy(np.int64), p[:, None])
    t["depth"] = t.depth * p
    return t


def one(args):
    j, i = args
    kind, label, proxy, f = JOBS[j]
    pool, rng = POOLS[proxy], np.random.default_rng([SEED, j, i])
    sets = f.groupby("donor", sort=False).compartment.agg(set)
    pick = draw(list(sets), rng, pool)  # (fetal donor index, adult donor)
    fd = sets.index[[k for k, _ in pick]]
    tab = []
    for k, d in pick:  # the adult donor takes the name, classes and depths of its fetal donor
        t = pool[d][pool[d].compartment.isin(sets.iloc[k])].assign(donor=sets.index[k])
        tab.append(t.merge(f.loc[f.donor == sets.index[k], ["compartment", "depth"]], on="compartment",
                           suffixes=("", "_fetal")))
    tab = pd.concat(tab, ignore_index=True)
    tab = thin(tab, tab.depth_fetal, rng)
    f = f[f.donor.isin(fd)]
    rows = []
    for tag in ("adjusted", "unadjusted"):
        Z = covariates(f, tag == "adjusted")
        out = []
        for _, g in tab.groupby("compartment", sort=False):
            g = g.copy()
            g[GENES] = np.exp(residual(logfrac(g), Z.loc[g.donor].to_numpy())) * DEPTH
            g["depth"] = DEPTH
            out.append(g)
        t = pd.concat(out, ignore_index=True)
        if kind == "L1":
            r, p = l1(logfrac(t), rng)
            rows.append(dict(kind=kind, label=label, sub=i, adjust=tag, n=len(pick), stat=r, p=p, passed=p < .05))
            continue
        d = f.drop_duplicates("donor").set_index("donor")
        strat = d.sex + np.where(d.wpf >= SPLIT, "|late", "|early")
        t["stratum"] = t.donor.map(strat) if tag == "adjusted" else "one"
        try:
            share, pp, pa = gate(t, SEED + i)
        except ValueError:  # no class pair shares >= 5 donors: the gate is not estimable at this design
            share, pp, pa = np.nan, np.nan, np.nan
        rows.append(dict(kind=kind, label=label, sub=i, adjust=tag, n=len(pick), stat=share, p=pp, p_acrophase=pa,
                         passed=bool(pp < .05 and pa < .05)))
    return rows


def setup():
    f = fetal_sheet()
    design = []
    for proxy, (name, stratum, adult) in PROXY.items():
        a, _ = data.load(name)
        a = a[(a.stratum == stratum) & a.compartment.isin(adult)].copy()
        a["compartment"] = a.compartment.map(dict(zip(adult, CLASSES)))
        POOLS[proxy] = {d: g.reset_index(drop=True) for d, g in a.groupby("donor")}
        n_pool = a.groupby("compartment").donor.nunique()
        for sub, rule in SUBSETS.items():
            s = f[rule(f)]
            jobs = [("L1", c, s[s.compartment == c]) for c in CLASSES]
            for tag, use in (("classes4", CLASSES), ("classes3", CLASSES[:3]), ("pair", CLASSES[:2])):
                x = s[s.compartment.isin(use)]
                jobs.append(("L2", tag, x[x.groupby("donor").compartment.transform("nunique") >= 2]))
            for kind, tag, x in jobs:
                n = x.donor.nunique()
                n_adult = n_pool.get(tag, 0) if kind == "L1" else len(POOLS[proxy])
                ok = n >= MIN_N and n_adult >= POOL_FACTOR * n
                if ok:
                    JOBS.append((kind, f"{proxy}:{tag}:{sub}", proxy, x))
                design.append(dict(kind=kind, label=f"{proxy}:{tag}:{sub}", fetal_n=n, adult_pool=n_adult,
                                   fetal_depth=";".join(f"{c}={v:.3g}" for c, v in
                                                        x.groupby("compartment").depth.median().items()),
                                   adult_depth=";".join(f"{c}={v:.3g}" for c, v in a[a.compartment.isin(
                                       x.compartment.unique())].groupby("compartment").depth.median().items()),
                                   calibrated=ok))
    return pd.DataFrame(design)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--n-sub", type=int, default=200)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--sex-trend", action="store_true", help="one age trend per sex (output suffix _sextrend)")
    ap.add_argument("--subset", help="run the jobs of one subset only (all, female, male, 13+wpf)")
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    global SEX_TREND
    SEX_TREND = a.sex_trend
    OUT.mkdir(exist_ok=True)
    t0 = time.time()
    limit_threads(1)
    sfx = "_gonad" + ("_sextrend" if a.sex_trend else "") + ("_check" if a.check else "")
    design = setup()
    design.to_csv(OUT / f"design{sfx}.csv", index=False)
    print(design.to_string(index=False), flush=True)
    n = 3 if a.check else a.n_sub
    tasks = [(j, i) for j in range(len(JOBS)) for i in range(n) if a.subset in (None, JOBS[j][1].rsplit(":", 1)[1])]
    print(f"{len(tasks) // n} calibration jobs x {n} subsamples", flush=True)
    run(tasks, one, a.jobs, sfx, t0)

if __name__ == "__main__":
    main()
