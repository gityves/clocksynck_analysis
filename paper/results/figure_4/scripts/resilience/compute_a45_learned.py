"""A45 sensitivity "the learned A17 anchor rule": inflamed vs non-inflamed gut (TAURUS) on anchors learned from other studies.

Pre-registered in docs/PLAN_internal_desync.md (archive) (A45, Sensitivities). Estimand and null as the A45 primary (compute_masld_gut.py,
`one_analysis`: T-cell d = c_inflamed - c_non-inflamed, companion d_T - d_fibroblast, labels permuted within patient, A22
own-d* detector, gates). Only the anchors differ: the A17 learned rule (leave one study out) applied to TAURUS, which is not
in the learning set. The rule (`host_identity/anchor_rule.py`, step `select`): per study x family the mean excess
of the healthy fits whose gate passes (tables a17_coherence_classes.csv and a17_gates.csv of host_identity, run
its compute_coherence.py first), pooled by tissue (gut) or over tissues; families that qualify (lower bound > 0) in the top half;
the target's own pair coherences (here `Dataset.coherence` on all samples) can only drop an anchor; < 3 anchors: no anchor.
The fixed rule's anchors (fibroblast, endothelial, mural) are the A45 primary.

    python compute_a45_learned.py [--subsample N] [--seed S]
Outputs (paper/results/figure_4/output/resilience/): a45_learned_{analyses,gates,contrast,companion,anchors,families}.csv
Legacy: legacy/desync/validation/ibd_taurus/scripts/a45_taurus.py (`a45_learned`, `learned_anchors`).
"""

from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(Path(__file__).resolve().parents[4]))  # paper/
from common.paths import group_dir  # noqa: E402
sys.path.insert(0, str(group_dir("host_identity")))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import anchor_rule as R  # noqa: E402  (local import of the S04 rule: the one definition of the learned rule)
import clocksync as cs  # noqa: E402
from common import cli, data  # noqa: E402
from common.paths import out_dir  # noqa: E402

import _shared as sh  # noqa: E402
import compute_masld_gut as M  # noqa: E402

FOLDER = "resilience"
S04 = "host_identity"


def learning_tables(quick: bool):
    """(family table per study, datasets) from the healthy fits of figure_S04 (gate passed, tissue not 'multi')."""
    d = out_dir(S04, False)
    d = d if (d / "a17_gates.csv").exists() else out_dir(S04, True)
    gates = pd.read_csv(d / "a17_gates.csv").set_index("name")
    classes = pd.read_csv(d / "a17_coherence_classes.csv")
    reg = sh.targets().set_index("name")
    names = [n for n in gates.index if gates.passed[n] and n in reg.index and reg.tissue[n] != "multi"]
    ft = R.family_table(classes[classes.name.isin(names)].rename(columns={"name": "dataset"}))
    ft["study"], ft["tissue"] = ft.dataset.map(reg.study), ft.dataset.map(reg.tissue)
    per = ft.groupby(["study", "tissue", "family"], as_index=False).agg(excess=("excess", "mean"), se=("se", "mean"))
    learn = per.rename(columns={"study": "dataset"})[["dataset", "family", "excess", "se"]]
    datasets = per[["study", "tissue"]].drop_duplicates().assign(dataset=lambda x: x.study, gate_pass=True)
    return learn, datasets


def main(argv=None):
    args = cli.parser(__doc__).parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    pb, sm = data.load("taurus")
    if args.subsample:
        pat = np.random.default_rng(args.seed).permutation(sorted(sm.patient.unique()))[:max(args.subsample // 4, 6)]
        sm = sm[sm.patient.isin(pat)]
        pb = pb[pb.donor.isin(sm.index)]
    el = M.eligible_samples(pb, sm, sm.index, "patient")
    lab, blk = sm.arm.reindex(el), sm.patient.reindex(el).astype(str)
    pbl = pb[pb.donor.isin(lab.index)].reset_index(drop=True)
    comps = sorted(set(pbl.compartment))
    # the target's own pair coherences (the legacy `coherence_confirm` on the fixed anchors; here all compartments)
    ds0 = cs.Dataset(pbl, anchors=None, donors=sm.loc[lab.index], name="a45_learned")
    pairs, _ = ds0.coherence(n_perm=run["n_perm"], n_boot=0, seed=args.seed)
    learn, datasets = learning_tables(run["quick"])
    anchors, status, fam = R.select(dict(study="TAURUS", tissue="gut"), comps, pairs, learn, datasets)
    fam.assign(analysis="a45_learned").to_csv(out / "a45_learned_families.csv", index=False)
    pd.DataFrame([dict(analysis="a45_learned", step="learned rule", anchors=";".join(anchors), status=status,
                       fixed_anchors=";".join(M.GUT_ANCHORS))]).to_csv(out / "a45_learned_anchors.csv", index=False)
    print("learned anchors:", anchors, status)
    if len(anchors) < 3:
        print("no anchor: contrast not run")
        return
    d = sm.assign(block=blk.reindex(sm.index))
    res = M.one_analysis("a45_learned", pb, d, lab, anchors, M.FIB, "block", run, args.seed + 500)
    M.merge([res], out, "a45_learned")
    for f in ("analyses", "contrast", "companion"):
        p = out / f"a45_learned_{f}.csv"
        if p.exists():
            print(f, "\n", pd.read_csv(p).round(3).to_string(index=False))


if __name__ == "__main__":
    main()
