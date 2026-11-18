"""A310 Stage 2a, block A input: donor-bootstrap intervals of each GTEx peripheral tissue's own variance 2 psi.

A104 rows and fit (GTEx fast deaths, tissues with >= 40 donors, periphery, every tissue an anchor, leave-one-out host);
200 donor resamples (`Dataset.resample`, host refitted). Replicate 0 is the full data and must reproduce A294b's clock
own variances. Aggregates only (GTEx controlled access).

    python boot_a310.py [--n-boot 200] [--jobs 4] [--check]
Output (`../output/`): boot.csv (replicate, tissue, S, own, coupling).
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
A294 = HERE.parents[2] / "A294_floor" / "output" / "tissues.csv"
SEED = 310
G = {}


def one(b):
    warnings.filterwarnings("ignore")
    limit_threads(1)
    ds = G["ds"] if b == 0 else G["ds"].resample(np.random.default_rng([SEED, b]))
    try:
        t = A104.per_tissue(ds)
    except (ValueError, np.linalg.LinAlgError) as e:
        return pd.DataFrame([dict(replicate=b, error=f"{type(e).__name__}: {e}"[:200])])
    return pd.DataFrame(dict(replicate=b, tissue=t.index, S=(t.z2 * t.V).to_numpy(), own=(2 * t.psi).to_numpy(),
                             coupling=t.coupling.to_numpy()))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--n-boot", type=int, default=200)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    t0 = time.time()
    pb, _ = gtex.load(min_donors=A104.MIN_DONORS)
    pb = pb[~pb.compartment.str.startswith("Brain")].reset_index(drop=True)
    comps = sorted(pb.compartment.unique())
    G["ds"] = cs.Dataset(pb, anchors=comps, name="gtex_periphery")
    full = one(0).set_index("tissue")
    ref = pd.read_csv(A294, index_col=0)
    dev = float(np.max(np.abs(full.own - ref.own.reindex(full.index))))
    print(f"{len(comps)} tissues; replicate 0 vs A294b own: max |diff| {dev:.2e} ({time.time() - t0:.0f}s)", flush=True)
    if dev > 1e-3:
        raise SystemExit("A294b own variances not reproduced: stop")
    n = 4 if a.check else a.n_boot
    OUT.mkdir(parents=True, exist_ok=True)
    res = [full.reset_index()]
    with Pool(a.jobs) as pool:
        for i, r in enumerate(pool.imap_unordered(one, range(1, n + 1)), 1):
            res.append(r)
            if i % 10 == 0:
                print(f"{i}/{n} replicates ({time.time() - t0:.0f}s)", flush=True)
    out = pd.concat(res, ignore_index=True)
    if a.check:
        print(out.groupby("replicate").size().to_string(), flush=True)
        return
    out.to_csv(OUT / "boot.csv", index=False)
    print(f"failed replicates {out.get('error', pd.Series(dtype=str)).notna().sum()}; done {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
