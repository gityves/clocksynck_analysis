"""A270 (registered in `paper/docs/registered_experiments.md`; NOT BLIND on the clock): do the responsive companion sets
carry internal time?

Per species, the clock and 200 responsive and 200 expression-matched pseudo-clock sets (`clocksync.paired.pseudo_waveforms`,
GTEx waveform amplitudes and phases; matching variable = mean log fraction over the rows used, genes with mean count >= 1).
Per set, the untimed host (every tissue an anchor, as A36 and `compute_walltime_hosts.py`) is read against recorded time
with `Dataset.timed` (R, permutation p, Fisher-Lee r); the reflected host is read against -t (descriptive).
  baboon  GSE98965, 12 animals, 63 tissues (A263 genome-wide table), ZT.
  human   GTEx fast deaths (Hardy 1-2), the 15 tissues of `V2_TISSUES`, depth = library size, time of death. Aggregates only.

    python score_a270.py --species baboon|human [--n-sets 200] [--n-perm 1000] [--jobs 6] [--check]
Outputs (`../output/`): sets_<species>.csv (one row per set: kind, R, p_R, r_fl, R_reflected), tests_<species>.csv.
"""
import argparse
import sys
import time
import warnings
from multiprocessing import Pool
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[3]))
sys.path.insert(0, str(HERE.parents[3] / "results" / "figure_2" / "scripts" / "organs"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.paired import pseudo_waveforms  # noqa: E402
from common import data  # noqa: E402
from common.paths import DATA  # noqa: E402

OUT = HERE.parents[1] / "output"
SEED = 270
WF = cs.gtex()
META = ["donor", "compartment", "depth"]
S = {}


def load(species):
    """(table of all genes, donor -> recorded hour)."""
    if species == "baboon":
        pb = pd.read_parquet(DATA / "derived" / "A263" / "baboon_genomewide.parquet")
        _, donors = data.load("baboon")
        return pb, donors.zt.astype(float)
    import anndata as ad
    from _shared import V2_TISSUES
    from common import gtex
    a = ad.read_h5ad(gtex.LABELLED, backed="r")
    genes = list(a.var_names)
    a.file.close()
    pb, donors = gtex.load(tissues=list(V2_TISSUES), hardy=(1, 2), genes=genes, depth="library_size")
    return pb, donors.tod.astype(float) % 24


def log_mean(pb):
    genes = [c for c in pb.columns if c not in META]
    x = pb[genes].to_numpy(np.float32)
    lm = pd.Series(np.log((x / pb.depth.to_numpy(float)[:, None]).mean(0) + 1e-12), index=genes)
    return lm[x.mean(0) >= 1]


def task(i):
    warnings.filterwarnings("ignore")
    kind, w = S["sets"][i]
    pb, hours, n_perm = S["pb"], S["hours"], S["n_perm"]
    out = dict(set=i, kind=kind, genes=";".join(w.genes))
    try:
        ds = cs.Dataset(pb[META + list(w.genes)], anchors=sorted(pb.compartment.unique()), waveform=w, name=f"set{i}")
        r = ds.timed(hours, n_perm=n_perm, seed=SEED)
        m = ds.timed((-hours) % 24, n_perm=1, seed=SEED)  # R only
        out.update(n=r.n, R=r.R, p_R=r.p_R, offset_h=r.offset_h, r_fl=r.r_fl, p_fl=r.p_fl, R_reflected=m.R)
    except Exception as e:  # noqa: BLE001
        out["error"] = f"{type(e).__name__}: {e}"[:200]
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--species", choices=["baboon", "human"], required=True)
    ap.add_argument("--n-sets", type=int, default=200)
    ap.add_argument("--n-perm", type=int, default=1000)
    ap.add_argument("--jobs", type=int, default=6)
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    t0, n = time.time(), (3 if a.check else a.n_sets)
    pb, hours = load(a.species)  # before limit_threads: pyarrow with OMP_NUM_THREADS=1 peaks far higher on wide tables
    limit_threads(1)
    OUT.mkdir(exist_ok=True)
    lm = log_mean(pb)
    print(f"{a.species}: {pb.donor.nunique()} donors, {pb.compartment.nunique()} tissues, {len(lm)} genes for matching "
          f"{time.time() - t0:.0f}s", flush=True)
    sets = [("clock", WF)]
    for j, kind in enumerate(("responsive", "matched"), 1):
        sets += [(kind, w) for w in pseudo_waveforms(lm, n, np.random.default_rng([SEED, j]), WF, kind)]
    genes = sorted({g for _, w in sets for g in w.genes})
    S.update(sets=sets, pb=pb[META + genes].copy(), hours=hours, n_perm=a.n_perm)
    del pb
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
    out = dict(species=a.species, n=c.n, R_clock=c.R, p_R_clock=c.p_R, r_fl_clock=c.r_fl)
    for kind in ("responsive", "matched"):
        x = R[R.kind == kind]
        ok = x.R.notna()
        out[f"{kind}_valid"] = int(ok.sum())
        out[f"{kind}_share_p05"] = float((x.p_R[ok] < 0.05).mean())
        out[f"{kind}_R_median"] = float(x.R[ok].median())
        out[f"{kind}_R_q95"] = float(x.R[ok].quantile(0.95))
        out[f"{kind}_maxR_median"] = float(np.maximum(x.R[ok], x.R_reflected[ok]).median())
        out[f"{kind}_pct_below_clock"] = float(100 * (x.R[ok] <= c.R).mean())
    out["carries_time"] = out["responsive_share_p05"] >= 0.5
    pd.DataFrame([out]).to_csv(OUT / f"tests_{a.species}{sfx}.csv", index=False)
    print(pd.Series(out).to_string(), f"runtime {time.time() - t0:.0f}s", sep="\n", flush=True)


if __name__ == "__main__":
    main()
