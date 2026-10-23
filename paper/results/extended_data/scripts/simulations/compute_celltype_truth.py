"""V13: cell-type coupling against a timed truth in mouse (known ZT). Which estimator is right per cell type?

Method check of A19 part 4 (docs/PLAN_internal_desync.md (archive); "not a hypothesis test", reported as run). Mouse single-cell time courses in
which every animal's time is known; no sample is assigned a time by the method. Per compartment, the coupling to the known time
(the same estimator with the host replaced by f = (cos t, sin t), no posterior uncertainty: `species/_shared.known_host`,
local: not in clocksync) against
  free_corrected: every compartment on the host of all the others (`Dataset.fit` with every compartment an anchor; A17 corrections);
  anchored: host from the stromal-vascular compartments (fixed A17 rule).
(The A1-A16 free estimator, `Study.fit` without small-sample corrections, is not in the public API: not ported.)
Data: colon = the pre-publication chrono-atlas (8 animals x proximal / distal colon, ZT 1/7/13/19; it is the regenerating colon,
R-DC + R-PC, and is reported as such; `data/colon_prepub`); skin = GSE223109 dorsal skin, 18 animals (ZT 2..70 every 4 h),
`data/skin_dermis_zt`, fibroblasts of each animal split at random into `--parts` exchangeable compartments used as anchors.
Pseudobulks >= 30 cells, class in >= 75 % of the animals. Intervals: animal bootstrap (`Dataset.resample`, host refitted), 2.5 / 97.5 %.
Readout per estimator: error against the truth over compartments (mean, mean absolute, by family), Spearman.

    python compute_celltype_truth.py [colon skin] [--subsample N] [--n-boot 200] [--parts 3] [--jobs 4]
`--subsample N`: bootstraps = N (default 3), skin restricted to nothing else (animals are few).
Outputs: paper/results/extended_data/output/simulations/v13_{colon,skin}_{couplings,summary}.csv
Legacy: legacy/scripts/v13_celltype_truth.py (outputs/v13_celltype_truth/)
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import cli  # noqa: E402
from common.paths import out_dir  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))  # paper/
from common.paths import DATASETS, group_dir  # noqa: E402
sys.path.insert(0, str(group_dir("species")))
import _shared as S12  # noqa: E402  (local helpers of figure_S12: known_host)

FOLDER = "simulations"
W = S12.W
SV = ("fibroblast", "smooth_muscle", "pericyte", "endothelial", "lymphatic_EC")  # stromal-vascular classes (A17 rule)
FAMILY = {"fibroblast": "fibroblast_stromal", "smooth_muscle": "mural", "pericyte": "mural", "endothelial": "endothelial",
          "lymphatic_EC": "endothelial", "T": "lymphocyte", "B": "lymphocyte", "plasma": "lymphocyte", "macrophage": "myeloid",
          "DC": "myeloid", "colonocyte": "epithelial", "secretory": "epithelial", "stem": "epithelial", "immune": "immune"}
MIN_CELLS, MIN_PRESENT = 30, 0.75


def _load_builder(name):
    import importlib.util
    f = DATASETS / name / "build.py"
    spec = importlib.util.spec_from_file_location(f"build_{name}", f)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def load(name, parts, seed):
    """Dataset-ready (pseudobulk, zt per animal, anchors)."""
    if name == "colon":
        pb, dt = _load_builder("colon_prepub").build()
    else:
        pb, dt = _load_builder("skin_dermis_zt").build(parts=parts, seed=seed)
    pb = pb[pb.n_cells >= MIN_CELLS]
    present = pb.groupby("compartment").donor.nunique() / pb.donor.nunique()
    pb = pb[pb.compartment.isin(present.index[present >= MIN_PRESENT])].reset_index(drop=True)
    comps = sorted(pb.compartment.unique())
    cls = lambda k: k.split(":")[-1]  # noqa: E731
    anchors = [k for k in comps if cls(k) in SV] if name == "colon" else [k for k in comps if cls(k).startswith("fibroblast_")]
    return pb, dt.zt, anchors


def truth(ds, zt) -> pd.Series:
    """Coupling of each compartment to the known time, on the rows of a (possibly resampled) dataset."""
    d = ds.data
    hours = pd.Series(d.resection).str.split("#").str[0].map(zt).to_numpy(float)
    f = np.column_stack([np.cos(W * hours), np.sin(W * hours)])
    comp = np.asarray(d.compartment)
    out = {}
    for k in ds.compartments:
        i = np.flatnonzero(comp == k)
        if len(i) >= 5:
            out[k] = S12.known_host(d.beta[i], d.V[i], f[i], np.asarray(d.group)[i])["coupling"]
    return pd.Series(out, dtype=float)


def one_boot(pb, donors, anchors, zt, seed):
    limit_threads(1)
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            base = cs.Dataset(pb, anchors, donors=donors, mixed_strata=True)
            b = base.resample(np.random.default_rng(seed))
            return dict(truth=truth(b, zt), anchored=b.fit(n_boot=0, gate=False).table.coupling,
                        free_corrected=cs.Dataset(b.table, None, donors=donors, mixed_strata=True).fit(n_boot=0, gate=False).table.coupling)
    except (ValueError, np.linalg.LinAlgError) as e:
        print(f"replicate {seed}: {type(e).__name__}", flush=True)
        return None


def run_one(name, n_boot, parts, seed, jobs, out):
    pb, zt, anchors = load(name, parts, seed)
    donors = pd.DataFrame({"zt": zt})
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        ds = cs.Dataset(pb, anchors, donors=donors, mixed_strata=True)
        point = {"truth": truth(ds, zt), "anchored": ds.fit(n_boot=0, gate=False).table.coupling,
                 "free_corrected": cs.Dataset(pb, None, donors=donors, mixed_strata=True).fit(n_boot=0, gate=False).table.coupling}
    print(f"{name}: {pb.donor.nunique()} animals, {len(ds.compartments)} compartments, anchors {anchors}", flush=True)
    reps = [r for r in Parallel(n_jobs=jobs)(delayed(one_boot)(pb, donors, anchors, zt, seed + b) for b in range(n_boot)) if r]
    rows = []
    for est, vals in point.items():
        for k in ds.compartments:
            x = np.array([r[est].get(k, np.nan) for r in reps], float)
            x = x[np.isfinite(x)]
            rows.append(dict(compartment=k, family=FAMILY.get(k.split(":")[-1], "other"), anchor=k in anchors, estimator=est,
                             coupling=vals.get(k, np.nan), lo=np.quantile(x, .025) if len(x) else np.nan,
                             hi=np.quantile(x, .975) if len(x) else np.nan, n_boot=len(x),
                             n_animals=int((pb.compartment == k).sum())))
    tab = pd.DataFrame(rows)
    tab.to_csv(out / f"v13_{name}_couplings.csv", index=False)
    wide = tab.pivot(index="compartment", columns="estimator", values="coupling")
    wide["family"] = [FAMILY.get(k.split(":")[-1], "other") for k in wide.index]
    summ = []
    for est in ("free_corrected", "anchored"):
        err = wide[est] - wide["truth"]
        for fam, e in [("all", err)] + list(err.groupby(wide.family)):
            summ.append(dict(estimator=est, family=fam, n=int(e.notna().sum()), mean_error=e.mean(), mean_abs_error=e.abs().mean()))
        summ.append(dict(estimator=est, family="spearman_vs_truth", n=int(err.notna().sum()),
                         mean_error=wide[est].corr(wide["truth"], method="spearman"), mean_abs_error=np.nan))
    summ = pd.DataFrame(summ)
    summ.to_csv(out / f"v13_{name}_summary.csv", index=False)
    print(wide.round(2).to_string(), "\n", summ.round(3).to_string(), flush=True)


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("sets", nargs="*", choices=["colon", "skin"])
    ap.add_argument("--n-boot", type=int, default=200)
    ap.add_argument("--parts", type=int, default=3, help="skin: fibroblast parts used as anchors (legacy 2; the package needs >= 3 anchors)")
    args = ap.parse_args(argv)
    run = cli.setup(args)
    n_boot = args.subsample if run["quick"] else args.n_boot
    out = out_dir(FOLDER, run["quick"])
    for name in args.sets or ["colon", "skin"]:
        run_one(name, n_boot, args.parts, args.seed, args.jobs, out)


if __name__ == "__main__":
    main()
