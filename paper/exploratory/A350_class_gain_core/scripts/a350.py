"""A350: A293 H293.4 (relative class gain conserved across organs) on the core loop.

Registered in paper/docs/registered_experiments.md (A350). Refits the A17 healthy atlas fits of A293 Part 3 (the atlases
not excluded there) with `compute_atlas_fits.fit_target` unchanged except the waveform (core loop = cs.gtex() without
PER1, PER2, NR1D1; learned anchors frozen from common/anchors_a17.csv; n_boot 500), then runs the A293 Part 3 code
unchanged on that table, and sets the 14-gene families beside the core ones (trio contribution, descriptive).

    python a350.py [--jobs 4] [--n-boot 500] [--only NAME ...]      (long: paper/common/run.sh, RUN_CAP_GB=30)
Outputs: ../output/{a17_healthy_fits_core.csv, gates_core.csv, part3/atlases.csv, part3/families.csv, compare.csv}
"""
from __future__ import annotations

import argparse
import contextlib
import io
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(PAPER / "results/figure_2/scripts/atlas_fits"))
sys.path.insert(0, str(PAPER / "exploratory/A293_organisation_levels/scripts"))
sys.path.insert(0, str(PAPER / "exploratory/A336_two_readouts/scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import _shared  # noqa: E402
import a336_common as C  # noqa: E402
import compute_atlas_fits as F  # noqa: E402
import part3_relative_gain as P3  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402

OUT = HERE.parents[1] / "output"
CORE = C.core()


def _core_dataset(self, anchors=None, **kw):
    import clocksync as cs
    kw.setdefault("waveform", CORE)
    return cs.Dataset(self.pb, anchors=anchors, donors=self.donors, name=self.name, **kw)


def fit_core(name, i, args, n_boot):
    _shared.Target.dataset = _core_dataset        # patched inside the worker (loky does not inherit the patch)
    F.Target.dataset = _core_dataset
    return F.fit_target(name, i, args, n_boot, F.N_PERM_GATE)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--n-boot", type=int, default=F.N_BOOT)
    ap.add_argument("--only", nargs="*", default=None)
    args = ap.parse_args()
    limit_threads(1)
    fa = argparse.Namespace(subsample=None, seed=0, anchor_readout="loo")   # compute_atlas_fits defaults (cli seed 0)
    reg = _shared.registry()
    names = [n for n in reg.name[reg.kind == "healthy"] if n not in P3.EXCLUDE and (args.only is None or n in args.only)]
    print(time.strftime("%F %T"), "atlases", len(names), flush=True)
    allh = list(reg.name[reg.kind == "healthy"])   # gate seed = index in the full healthy list, as the 14-gene run
    res = Parallel(n_jobs=args.jobs)(delayed(fit_core)(n, allh.index(n), fa, args.n_boot) for n in names)
    for r in res:
        if "error" in r:
            print("FAILED", r["error"], flush=True)
    res = [r for r in res if "error" not in r]
    OUT.mkdir(parents=True, exist_ok=True)
    fits = pd.concat([r["rows"] for r in res], ignore_index=True)
    fits.to_csv(OUT / "a17_healthy_fits_core.csv", index=False)
    pd.DataFrame([r["gate"] for r in res]).to_csv(OUT / "gates_core.csv", index=False)
    print(time.strftime("%F %T"), "fits done", flush=True)
    # A293 Part 3 unchanged on the core table
    P3.FITS, P3.OUT = OUT / "a17_healthy_fits_core.csv", OUT / "part3"
    P3.main()
    core = pd.read_csv(OUT / "part3/families.csv")
    f14 = pd.read_csv(PAPER / "exploratory/A293_organisation_levels/output/part3/families.csv")
    cmp_ = f14.merge(core, on="family", suffixes=("_14", "_core"))
    cmp_["trio_contribution_g"] = cmp_.pooled_14 - cmp_.pooled_core
    cmp_["se_ratio_core_14"] = cmp_.se_pooled_core / cmp_.se_pooled_14
    cmp_.to_csv(OUT / "compare.csv", index=False)
    print(cmp_.round(3).to_string(index=False))
    t = core[core.tested]
    print("H350.1", "PASS" if len(t) and t.within_bound.all() else "FAIL", f"({len(t)} families tested)", flush=True)


if __name__ == "__main__":
    main()
