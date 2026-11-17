"""A300 (registered in `paper/docs/registered_experiments.md`): do the tissues facing the external environment change their
clock relation to the body with the season?

GTEx fast deaths, A104 periphery rows. Summer and Winter death seasons (`DTHSEASON`, protected; aggregates only). Per
arm the A104 periphery fit (every tissue an anchor, leave-one-out host); per tissue offset (h) and log z2.
Delta = offset_summer - offset_winter. T = median |Delta| (frontier) - median |Delta| (internal). Null: season labels
permuted over the summer and winter donors (500), both arms refitted.

    python score_a300.py [--n-perm 500] [--jobs 4]
Outputs (`../output/`): tissues.csv, tests.csv.
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
from clocksync.config import limit_threads  # noqa: E402
from common import gtex  # noqa: E402

OUT = HERE.parents[1] / "output"
FRONTIER = ["Lung", "Skin - Sun Exposed (Lower leg)", "Skin - Not Sun Exposed (Suprapubic)", "Esophagus - Mucosa",
            "Colon - Transverse", "Stomach"]
EXCLUDE = ["Whole Blood", "Testis"]
MIN_ARM = 20
S = {}


def wrap(h):
    return (np.asarray(h, float) + 12) % 24 - 12


def arm_fit(donors):
    pb = S["pb"]
    t = pb[pb.donor.isin(donors)]
    ds = cs.Dataset(t, anchors=S["comps"], name="arm")
    tab = ds.fit(n_boot=0, gate=False).table
    return tab.offset_h.reindex(S["comps"]), np.log(tab.z2.reindex(S["comps"]))


def stat(lab):
    """lab: boolean array over S['donors'] (True = summer)."""
    d = S["donors"]
    o1, z1 = arm_fit(d[lab])
    o0, z0 = arm_fit(d[~lab])
    delta = pd.Series(wrap(o1 - o0), index=S["comps"])
    dz = z1 - z0
    fr = delta.index.isin(FRONTIER)
    T = float(np.nanmedian(np.abs(delta[fr])) - np.nanmedian(np.abs(delta[~fr])))
    Tz = float(np.nanmedian(np.abs(dz[fr])) - np.nanmedian(np.abs(dz[~fr])))
    return T, Tz, delta, dz


def perm(seed):
    warnings.filterwarnings("ignore")
    limit_threads(1)
    lab = np.random.default_rng([300, seed]).permutation(S["lab"])
    try:
        T, Tz, _, _ = stat(lab)
        return T, Tz
    except Exception:  # noqa: BLE001
        return np.nan, np.nan


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--n-perm", type=int, default=500)
    ap.add_argument("--jobs", type=int, default=4)
    a = ap.parse_args()
    warnings.filterwarnings("ignore")
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    pb, donors = gtex.load(min_donors=A104.MIN_DONORS)
    limit_threads(1)
    season = gtex.phenotypes(["DTHSEASON"]).DTHSEASON.reindex(donors.index).astype(str)
    keep = season.index[season.isin(["Summer", "Winter"])]
    pb = pb[~pb.compartment.str.startswith("Brain") & ~pb.compartment.isin(EXCLUDE) & pb.donor.isin(keep)]
    n = pb.assign(s=pb.donor.map(season)).groupby(["compartment", "s"]).donor.nunique().unstack()
    comps = sorted(n.index[(n.min(axis=1) >= MIN_ARM)])
    pb = pb[pb.compartment.isin(comps)]
    d = np.array(sorted(pb.donor.unique()))
    S.update(pb=pb, comps=comps, donors=d, lab=(season.reindex(d) == "Summer").to_numpy())
    print(f"{len(comps)} tissues ({sum(c in FRONTIER for c in comps)} frontier), {S['lab'].sum()} summer / {(~S['lab']).sum()} winter donors "
          f"{time.time() - t0:.0f}s", flush=True)
    T, Tz, delta, dz = stat(S["lab"])
    with Pool(a.jobs) as pool:
        null = np.array(pool.map(perm, range(a.n_perm)))
    ok = np.isfinite(null[:, 0])
    p = (1 + np.sum(null[ok, 0] >= T)) / (1 + ok.sum())
    pz = (1 + np.sum(null[ok, 1] >= Tz)) / (1 + ok.sum())
    tis = pd.DataFrame(dict(frontier=[c in FRONTIER for c in comps], delta_offset_h=delta.to_numpy(), delta_log_z2=dz.to_numpy(),
                            n_summer=n.loc[comps, "Summer"].to_numpy(), n_winter=n.loc[comps, "Winter"].to_numpy()), index=comps)
    tis.to_csv(OUT / "tissues.csv")
    tests = pd.DataFrame([dict(test="H300.1", T_h=T, p_one_sided=p, n_perm=int(ok.sum()),
                               median_abs_delta_frontier_h=float(tis[tis.frontier].delta_offset_h.abs().median()),
                               median_abs_delta_internal_h=float(tis[~tis.frontier].delta_offset_h.abs().median()),
                               passes=bool(T > 0 and p < 0.05)),
                          dict(test="descriptive log z2", T_h=Tz, p_one_sided=pz, n_perm=int(ok.sum()))])
    tests.to_csv(OUT / "tests.csv", index=False)
    print(tis.round(3).to_string())
    print(tests.round(3).T.to_string())
    print(f"done {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
