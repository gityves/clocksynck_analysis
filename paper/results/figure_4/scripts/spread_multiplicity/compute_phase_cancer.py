"""A23, R23.4 (cancer arm): the phase of tumour epithelium relative to adjacent-normal epithelium on the stromal-vascular host.

Pre-registered in docs/PLAN_internal_desync.md (archive) (A23, R23.4). Per informative cancer study with learned A17 anchors: one host on
the anchors minus tumour:epithelial and adjacent_normal:epithelial (the registration says >= 2 left; `Dataset.fit` needs 3, see
the README); on it y4 = arg(z_tumour epithelium / z_adjacent-normal epithelium) x 24 / 2 pi (h, wrapped; > 0: the tumour
epithelium runs later). SE and interval by patient bootstrap (500, host refitted); studies pooled by DerSimonian-Laird (rows with
SE > 4 h excluded). A null is stated with its interval; shifts below the MDE (2.8 SE) are not "no shift".

    python compute_phase_cancer.py [--subsample N] [--jobs 4]
Outputs (paper/results/figure_4/output/spread_multiplicity/): a23_cancer.csv (per study and the pooled row)
Legacy: legacy/desync/common/a23_phase.py (`cancer`, `_cancer_one`).
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import cli, data  # noqa: E402
from common.anchors import TABLE, anchors  # noqa: E402
from common.paths import out_dir  # noqa: E402

FOLDER = "spread_multiplicity"
TE, NE = "tumour:epithelial", "adjacent_normal:epithelial"
SE_MAX = 4.0


def wrap(h):
    return (np.asarray(h, float) + 12.0) % 24.0 - 12.0


def y4(ds: cs.Dataset) -> float:
    """Tumour minus adjacent-normal epithelium phase (h) on the host of the anchors (offsets of one fit)."""
    t = ds.fit(n_boot=0, gate=False, targets=[TE, NE]).table
    ok = (t.z2 > 0).all()
    return float(wrap(t.offset_h[TE] - t.offset_h[NE])) if ok else np.nan


def one(study: str, seed: int, n_boot: int, subsample, sseed) -> dict:
    limit_threads(1)
    warnings.filterwarnings("ignore")
    name = f"cancer_{study}"
    pb, _ = data.load("cancer_sc", study=study, cnv=False)
    if subsample is not None:
        keep = np.random.default_rng(sseed).permutation(sorted(pb.donor.unique()))[:max(subsample, 30)]
        pb = pb[pb.donor.isin(keep)].reset_index(drop=True)
    comps = set(pb.compartment)
    anc = [a for a in anchors(name, "learned") if a not in (TE, NE) and a in comps]
    both = set(pb[pb.compartment == TE].donor) & set(pb[pb.compartment == NE].donor)
    if TE not in comps or NE not in comps or len(both) < 5:
        return dict(study=name, status="not estimable", n_patients=len(both), n_anchors=len(anc))
    if len(anc) < 3:
        return dict(study=name, status="not estimable (Dataset.fit needs 3 anchors)", n_patients=len(both), n_anchors=len(anc))
    ds = cs.Dataset(pb, anchors=anc, name=study)
    try:
        y = y4(ds)
    except (ValueError, np.linalg.LinAlgError):
        y = np.nan
    if not np.isfinite(y):
        return dict(study=name, status="not estimable", n_patients=len(both), n_anchors=len(anc))
    rng = np.random.default_rng(seed)
    B = []
    for _ in range(n_boot):
        try:
            B.append(y4(ds.resample(rng)))
        except (ValueError, np.linalg.LinAlgError):
            continue
    B = np.array(B, float)
    B = y + wrap(B[np.isfinite(B)] - y)
    return dict(study=name, status="ok", y4=y, se4=float(np.std(B, ddof=1)), lo=float(np.quantile(B, 0.025)),
                hi=float(np.quantile(B, 0.975)), n_boot=len(B), n_patients=len(both), n_anchors=len(anc))


def main(argv=None):
    args = cli.parser(__doc__).parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    n_boot = 20 if run["quick"] else 500
    reg = pd.read_csv(TABLE)
    studies = sorted(t.removeprefix("cancer_") for t in reg.target.unique() if t.startswith("cancer_"))
    if run["quick"]:
        studies = [s for s in studies if s in ("GSE132465", "GSE178341", "SRP273482")]
    res = Parallel(n_jobs=args.jobs)(delayed(one)(s, 8000 + i, n_boot, args.subsample, args.seed) for i, s in enumerate(studies))
    t = pd.DataFrame(res)
    ok = t[(t.status == "ok") & (t.se4 <= SE_MAX)] if "se4" in t else t.iloc[:0]
    if len(ok) >= 2:
        m = cs.stats.meta(ok.y4, ok.se4)
        t = pd.concat([t, pd.DataFrame([dict(study="pooled (DL)", status=f"{len(ok)} studies", y4=m["random"], se4=m["se_random"],
                                             lo=m["lo"], hi=m["hi"], mde_h=2.8 * m["se_random"])])], ignore_index=True)
    t.to_csv(out / "a23_cancer.csv", index=False)
    print(t.round(2).to_string(index=False))


if __name__ == "__main__":
    main()
