"""A263 (registered in `paper/docs/registered_experiments.md`; NOT BLIND): is the lower coupling of brain to the body's
clock state specific to the clock genes, in human, sheep and baboon?

Per species, one table of all compartments (GTEx v10 fast deaths with tissues of >= 40 donors, as A142; SheepGTEx and
baboon genome-wide tables of `build_genomewide.py`, same rows as the A42 builders). For the 14 clock genes and for each of
200 expression-matched and 200 responsive pseudo-clock sets (`clocksync.paired.pseudo_waveforms`, matching variable =
mean log fraction over all rows of the species, genes with mean count >= 1):
  R2   coupling of every tissue to the peripheral-anchor host (`Dataset.fit(n_boot=0, gate=False)`, anchors = periphery,
       peripheral tissues read leave-one-out; A42 readout 2); L = log(median coupling brain / median coupling periphery).
  R1   (human and sheep) block agreement r* of brain vs periphery hosts (`Dataset.agreement`, n_boot=0; A42 readout 1).
Brain / periphery / excluded tissues as A42 (`compute_brain_periphery_species.groups`); human: brain = the "Brain -"
regions, excluded = Pituitary, periphery = the rest. Sheep with mixed strata (as A42).

    python score_a263.py --species human|sheep|baboon [--n-sets 200] [--jobs 8] [--check]
Outputs (`../output/`): sets_<species>.csv (one row per gene set: kind, L, medians, r*), tests_<species>.csv.
"""
import argparse
import sys
import time
from multiprocessing import Pool
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[3]))
sys.path.insert(0, str(HERE.parents[3] / "results" / "figure_2" / "scripts" / "organs"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
import compute_brain_periphery_species as SP  # noqa: E402
import _shared as SH  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.paired import pseudo_waveforms  # noqa: E402
from common.paths import DATA  # noqa: E402

OUT = HERE.parents[1] / "output"
SEED = 263
WF = cs.gtex()
S = {}


def groups(species, comps):
    if species == "human":
        brain = [c for c in comps if c.startswith("Brain -")]
        excl = [c for c in comps if c == "Pituitary"]
        return brain, [c for c in comps if c not in brain and c not in excl], excl
    return SP.groups(species, comps)


def table(species, genes=None):
    """All rows of a species with the given gene columns (None: every gene)."""
    if species == "human":
        pb, _ = SH.cohort("fast", genes="clock" if genes is None else list(genes))
        return pb
    t = pd.read_parquet(DATA / "derived" / "A263" / f"{species}_genomewide.parquet",
                        columns=None if genes is None else ["donor", "compartment", "depth",
                                                            *(["stratum"] if species == "sheep" else []), *genes])
    return t


def log_mean(species):
    """Mean log fraction of every gene (mean count >= 1) over all rows: the matching variable."""
    if species == "human":
        from common import gtex
        import anndata as ad
        a = ad.read_h5ad(gtex.LABELLED, backed="r")
        genes = list(a.var_names)
        a.file.close()
        pb = table("human", genes)
    else:
        pb = table(species)
        genes = [c for c in pb.columns if c not in ("donor", "compartment", "depth", "stratum")]
    x = pb[genes].to_numpy(np.float32)
    lm = pd.Series(np.log((x / pb.depth.to_numpy(float)[:, None]).mean(0) + 1e-12), index=genes)
    return lm[x.mean(0) >= 1]


def readouts(pb, wf, species):
    brain, per, excl = groups(species, sorted(pb.compartment.unique()))
    cols = [c for c in ("donor", "compartment", "depth", "stratum") if c in pb] + list(wf.genes)
    ds = cs.Dataset(pb[cols], anchors=per, waveform=wf, mixed_strata=species == "sheep", name=species, min_donors=5)
    tab = ds.fit(n_boot=0, gate=False, targets=brain + per + excl).table
    cb, cp = tab.coupling.reindex(brain).dropna(), tab.coupling.reindex(per).dropna()
    out = dict(median_brain=float(cb.median()), median_periph=float(cp.median()), n_brain=len(cb), n_periph=len(cp))
    out["L"] = float(np.log(out["median_brain"] / out["median_periph"])) if out["median_brain"] > 0 and out["median_periph"] > 0 else np.nan
    if species in ("human", "sheep"):
        try:
            a = ds.agreement({"brain": brain, "periphery": per}, n_boot=0, seed=SEED)
            out["r_star"] = float(a.r_star)
        except Exception as e:  # undefined reliability (no shared state in a block)
            out["r_star"], out["r_star_error"] = np.nan, type(e).__name__
    return out


def task(i):
    kind, w = S["sets"][i]
    pb = S["pb"]
    t = pb[[c for c in pb.columns if c in ("donor", "compartment", "depth", "stratum")] + list(w.genes)]
    try:
        r = readouts(t, w, S["species"])
    except Exception as e:
        r = dict(error=f"{type(e).__name__}: {e}"[:200])
    return dict(set=i, kind=kind, genes=";".join(w.genes), **r)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--species", choices=["human", "sheep", "baboon"], required=True)
    ap.add_argument("--n-sets", type=int, default=200)
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    limit_threads(1)
    OUT.mkdir(exist_ok=True)
    t0, n = time.time(), (3 if a.check else a.n_sets)
    lm = log_mean(a.species)
    print(f"{a.species}: {len(lm)} genes for matching {time.time() - t0:.0f}s", flush=True)
    sets = [("clock", WF)]
    for j, kind in enumerate(("matched", "responsive"), 1):
        sets += [(kind, w) for w in pseudo_waveforms(lm, n, np.random.default_rng([SEED, j]), WF, kind)]
    genes = sorted({g for _, w in sets for g in w.genes})
    pb = table(a.species, genes)
    S.update(sets=sets, pb=pb, species=a.species)
    rows = []
    with Pool(a.jobs, initializer=limit_threads, initargs=(1,)) as p:
        for k, r in enumerate(p.imap_unordered(task, range(len(sets)))):
            rows.append(r)
            if (k + 1) % 25 == 0 or k + 1 == len(sets):
                print(f"{k + 1}/{len(sets)} {time.time() - t0:.0f}s", flush=True)
    R = pd.DataFrame(rows).sort_values("set")
    sfx = "_check" if a.check else ""
    R.to_csv(OUT / f"sets_{a.species}{sfx}.csv", index=False)
    c = R[R.kind == "clock"].iloc[0]
    out = dict(species=a.species, L_clock=c.L, median_brain=c.median_brain, median_periph=c.median_periph,
               n_brain=c.n_brain, n_periph=c.n_periph)
    spec = True
    for kind in ("matched", "responsive"):
        x = R[R.kind == kind]
        out[f"{kind}_valid"] = int(x.L.notna().sum())
        out[f"{kind}_L_q05"] = float(np.nanquantile(x.L, 0.05))
        out[f"{kind}_L_median"] = float(np.nanmedian(x.L))
        out[f"{kind}_pct_below_clock"] = float(100 * np.mean(x.L.dropna() <= c.L))
        spec &= bool(c.L < out[f"{kind}_L_q05"])
    out["L_clock_specific"] = spec
    if "r_star" in R:
        out["r_star_clock"] = c.r_star
        spec_r = True
        for kind in ("matched", "responsive"):
            x = R[R.kind == kind].r_star
            out[f"{kind}_rstar_defined"] = int(x.notna().sum())
            out[f"{kind}_rstar_q05"] = float(np.nanquantile(x, 0.05)) if x.notna().any() else np.nan
            out[f"{kind}_rstar_median"] = float(np.nanmedian(x)) if x.notna().any() else np.nan
            spec_r &= bool(c.r_star < out[f"{kind}_rstar_q05"]) if x.notna().any() else False
        out["rstar_clock_specific"] = spec_r
    pd.DataFrame([out]).to_csv(OUT / f"tests_{a.species}{sfx}.csv", index=False)
    print(pd.Series(out).to_string(), f"runtime {time.time() - t0:.0f}s", sep="\n", flush=True)


if __name__ == "__main__":
    main()
