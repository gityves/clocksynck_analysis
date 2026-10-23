"""A17 step 2 for the 21 cancer targets: the anchors of every single-cell cancer cohort, learned leave one study out.

Pre-registered in docs/PLAN_internal_desync.md (archive) (A17 step 2; the 21 informative cancer studies of the 2026-10-03 pipeline). Per
cohort (main table: tumour, adjacent-normal and metastasis compartments, free host): the model-free clock-plane coherence of
every compartment pair over the patients it shares (centred within strata; excess over 200 shuffled pairings; label-free, no
bootstrap: the cancer tables enter the anchor rule only through their own pairs, as the contrast tables do). The anchor rule is
that of host_identity (`anchor_rule.select`): the learning set is the gate-passing healthy fits of
`host_identity/a17_{gates,coherence_classes}.csv` (run `compute_coherence.py` there first), one value per study x
family, anchors from the other studies; the cohort's own pair coherence can only drop an anchor; fewer than 3 anchors:
"no anchor". The fixed rule (stromal-vascular compartments) is the sensitivity. The paper uses the frozen table
(`common/anchors_a17.csv`, committed before any anchored coupling); this script re-derives the 21 cancer rows as a check.

    python compute_anchors_cancer.py [--subsample N] [--jobs 4]
Outputs (paper/results/figure_5/output/cancer/): anchors_cancer_relearned.csv (columns of common/anchors_a17.csv),
anchors_cancer_vs_frozen.csv (per cohort: compartments, anchors learned here / frozen, identical?; fixed rule likewise)
Legacy: legacy/desync/common/anchors.py (cancer targets of `targets`, `table_of`, `_one`; step select).
"""

from __future__ import annotations

import importlib.util
import sys
import warnings
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAPER = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(Path(__file__).resolve().parents[4]))  # paper/
from common.paths import group_dir, group_out  # noqa: E402
sys.path[:0] = [str(group_dir("host_identity")), str(group_dir("atlas_fits")), str(Path(__file__).resolve().parents[4])]
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import cli, data  # noqa: E402
from common.anchors import TABLE  # noqa: E402
from common.paths import out_dir  # noqa: E402

FOLDER = "cancer"
_spec = importlib.util.spec_from_file_location("s04_compute_anchors", group_dir("host_identity") / "compute_anchors.py")
S04 = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(S04)  # uses atlas_fits/_shared.py (registry of the healthy and contrast targets)
from _shared import registry  # noqa: E402  (that module)

QUICK = ("SRP273482", "GSE132465", "GSE181919")


def pairs_of(study: str, i: int, subsample, seed: int, n_perm: int) -> pd.DataFrame:
    """Pair coherence of the cohort's compartments (free host, `Dataset.coherence`, label-free)."""
    limit_threads(1)
    warnings.filterwarnings("ignore")
    pb, _ = data.load("cancer_sc", study=study, cnv=False)
    if subsample is not None:
        keep = np.random.default_rng(seed).permutation(sorted(pb.donor.unique()))[:max(subsample, 12)]
        pb = pb[pb.donor.isin(keep)].reset_index(drop=True)
    ds = cs.Dataset(pb, anchors=None, name=study)
    pairs, _ = ds.coherence(n_perm=n_perm, n_boot=0, seed=i)
    return pairs.assign(name=f"cancer_{study}")


def main(argv=None):
    args = cli.parser(__doc__).parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    s04 = group_out("host_identity") / ("subsample" if run["quick"] else "")
    gates, classes = pd.read_csv(s04 / "a17_gates.csv"), pd.read_csv(s04 / "a17_coherence_classes.csv")
    fz = pd.read_csv(TABLE)
    cz = fz[fz.target.str.startswith("cancer_")].groupby("target").first()
    studies = [t.removeprefix("cancer_") for t in cz.index]
    studies = [s for s in studies if s in QUICK] if run["quick"] else studies
    res = Parallel(n_jobs=args.jobs)(delayed(pairs_of)(s, i, args.subsample, args.seed, 50 if run["quick"] else 200)
                                     for i, s in enumerate(studies))
    pairs = pd.concat(res, ignore_index=True)
    reg = pd.concat([registry(), pd.DataFrame({"name": cz.index, "kind": "cancer", "study": cz.study, "tissue": cz.tissue},
                                              index=cz.index)])
    reg = reg[~reg.index.duplicated()]
    tab, _ = S04.learn(gates, classes, pairs, reg.loc[[n for n in reg.index if n in set(gates.name) or n in set(pairs.name)]])
    tab.to_csv(out / "anchors_cancer_relearned.csv", index=False)
    cmp = S04.compare(tab[tab.target.isin(cz.index)])
    cmp.to_csv(out / "anchors_cancer_vs_frozen.csv", index=False)
    print(cmp.to_string(index=False))
    print(cmp[["learned_identical", "fixed_identical"]].mean().round(3).to_dict(), len(cmp), "cohorts")


if __name__ == "__main__":
    main()
