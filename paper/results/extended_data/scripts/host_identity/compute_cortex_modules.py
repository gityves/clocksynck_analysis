"""A18 part 3.5: cortex modules. Does the neuronal host differ from the non-neuronal host?

Pre-registered as A18 part 3 item 5 (2026-10-12; docs/PLAN_internal_desync.md (archive); exploratory): "Cortex modules. Heatmap of the pair excess coherence; neuron-host vs
glia-host similarity per cortex dataset." The pair excess coherence is `a17_coherence_pairs.csv` of `compute_coherence.py`. Here the host similarity: per cortex dataset
(PsychAD MSSM and HBCC controls, SEA-AD DLPFC and MTG, Parkinson controls of the A9 regions), the host fitted on the neuron compartments and the host fitted on the
non-neuronal compartments (glia, myeloid, endothelial, mural; families of A17, `anchor_rule.family`), and the |complex correlation| of the two donor states
(centred within strata; `Dataset.host`); the within-module reference is the same similarity between two random halves of the same family (10 splits), which
shows how high the similarity of two hosts of one module is. Estimand: similarity(neuron host, non-neuronal host) against the within-module values; no null test.

    python compute_cortex_modules.py [--subsample N]
Outputs: paper/results/extended_data/output/host_identity/a18_modules.csv
Legacy: legacy/desync/common/followups.py (`modules`) -> common_output/a18_modules.csv
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))  # paper/
from common.paths import group_dir  # noqa: E402
sys.path.insert(0, str(group_dir("atlas_fits")))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from _shared import FOLDER_S04, Target, donor_strata, similarity  # noqa: E402
from anchor_rule import family  # noqa: E402
from common import cli  # noqa: E402
from common.paths import out_dir  # noqa: E402

TARGETS = ("psychad_mssm", "psychad_hbcc", "seaad_dlpfc", "seaad_mtg", "pd_control")
NON_NEURONAL = ("glia", "myeloid", "endothelial", "mural")


def sim(ds, a, b, strata):
    return similarity(ds.host(anchors=a)[["x", "y"]], ds.host(anchors=b)[["x", "y"]], strata)


def main(argv=None):
    ap = cli.parser(__doc__)
    args = ap.parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER_S04, run["quick"])
    warnings.filterwarnings("ignore")
    rows = []
    for name in TARGETS:
        T = Target(name, args.subsample, args.seed)
        ds = T.dataset(anchors=None)
        comps = sorted(set(T.pb.compartment))
        neu = [k for k in comps if family(k) == "neuron"]
        gli = [k for k in comps if family(k) in NON_NEURONAL]
        if len(neu) < 2 or len(gli) < 2:
            continue
        st = donor_strata(ds)
        rows.append(dict(target=name, comparison="neuron vs non-neuronal", neurons=len(neu), non_neuronal=len(gli),
                         donors=len(ds.donor_ids), similarity=sim(ds, neu, gli, st)))
        rng = np.random.default_rng(0)
        for fam, ks in (("neuron", neu), ("non-neuronal", gli)):
            if len(ks) < 4:
                continue
            for i in range(3 if run["quick"] else 10):
                p = rng.permutation(ks)
                rows.append(dict(target=name, comparison=f"{fam} half vs half", split=i, neurons=len(neu), non_neuronal=len(gli),
                                 donors=len(ds.donor_ids), similarity=sim(ds, sorted(p[: len(p) // 2]), sorted(p[len(p) // 2:]), st)))
        print(name, flush=True)
    R = pd.DataFrame(rows)
    R.to_csv(out / "a18_modules.csv", index=False)
    pd.set_option("display.width", 200)
    print(R.groupby(["target", "comparison"]).similarity.agg(["mean", "min", "max", "count"]).round(3))


if __name__ == "__main__":
    main()
