"""A17 step 2: anchors learned leave one study out from the coherence of the healthy fits; the fixed stromal-vascular rule
as sensitivity. Reads the coherence tables of `compute_coherence.py` (same output folder) and writes the anchor table of
every healthy and contrast target; compares it with the frozen table `common/anchors_a17.csv` (the legacy run).

Pre-registered in docs/PLAN_internal_desync.md (archive) (A17 step 2). Learning set: healthy fits whose gate passes (not Tabula
Sapiens, not tissue "multi"); one value per study x family (mean over the study's fits, family value = mean excess of its
classes, SE = mean SE); anchors per target from the other studies; the target's own pair coherence can only drop an anchor;
< 3 anchors: "no anchor". The paper uses the frozen table (committed before any anchored coupling); this script
re-derives it as a check and is the code of the rule.

    python compute_anchors.py [--subsample N]      (run compute_coherence.py first)
Outputs: anchors_a17_relearned.csv (columns of common/anchors_a17.csv), a17_anchor_families.csv, anchors_vs_frozen.csv
(per target: compartments, anchors learned here / frozen, identical?; fixed rule likewise).
Legacy: legacy/desync/common/anchors.py (step select), src/sccc/atlas/anchor_rule.py.
"""

from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(Path(__file__).resolve().parents[4]))  # paper/
from common.paths import group_dir  # noqa: E402
sys.path.insert(0, str(group_dir("atlas_fits")))
import pandas as pd  # noqa: E402

import anchor_rule as R  # noqa: E402
from _shared import FOLDER_S04, frozen, registry  # noqa: E402
from common import cli  # noqa: E402
from common.paths import out_dir  # noqa: E402


def learn(gates: pd.DataFrame, classes: pd.DataFrame, pairs: pd.DataFrame, reg: pd.DataFrame):
    """Anchor table (and family tables) of every target of `reg` that has pairs, from the healthy coherence tables."""
    gates = gates.set_index("name")
    learn_names = [n for n in gates.index if gates.passed[n] and reg.tissue[n] != "multi"]
    ft = R.family_table(classes[classes.name.isin(learn_names)].rename(columns={"name": "dataset"}))
    ft["study"] = ft.dataset.map(reg.study)
    ft["tissue"] = ft.dataset.map(reg.tissue)
    per_study = ft.groupby(["study", "tissue", "family"], as_index=False).agg(excess=("excess", "mean"),
                                                                             se=("se", "mean"))
    learning = per_study.rename(columns={"study": "dataset"})[["dataset", "family", "excess", "se"]]
    datasets = per_study[["study", "tissue"]].drop_duplicates().assign(dataset=lambda d: d.study, gate_pass=True)
    rows, fams = [], []
    for name in reg.name:
        p_t = pairs[pairs.name == name]
        if p_t.empty:
            continue
        comps = sorted(set(p_t.k) | set(p_t.l))
        t = reg.loc[name]
        anchors, status, tab = R.select(dict(study=t.study, tissue=t.tissue), comps, p_t, learning, datasets)
        fx = R.fixed(dict(tissue=t.tissue), comps)
        fams.append(tab.assign(target=name))
        for k in comps:
            rows.append(dict(target=name, kind=t.kind, study=t.study, tissue=t.tissue, compartment=k,
                             family=R.family(k), anchor=k in anchors, anchor_fixed=k in fx, status=status))
    return pd.DataFrame(rows), pd.concat(fams)


def compare(new: pd.DataFrame) -> pd.DataFrame:
    old = frozen()
    out = []
    for name, g in new.groupby("target"):
        o = old[old.target == name].set_index("compartment")
        g = g.set_index("compartment")
        com = g.index.intersection(o.index)
        out.append(dict(target=name, kind=g.kind.iloc[0], compartments=len(g), compartments_frozen=len(o),
                        anchors=int(g.anchor.sum()), anchors_frozen=int(o.anchor.sum()),
                        learned_identical=bool(g.anchor[com].equals(o.anchor[com].astype(bool)) and len(com) == len(g)),
                        fixed_identical=bool(g.anchor_fixed[com].equals(o.anchor_fixed[com].astype(bool))),
                        status=g.status.iloc[0], status_frozen=o.status.iloc[0] if len(o) else ""))
    return pd.DataFrame(out)


def main(argv=None):
    args = cli.parser(__doc__).parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER_S04, run["quick"])
    gates, classes = pd.read_csv(out / "a17_gates.csv"), pd.read_csv(out / "a17_coherence_classes.csv")
    pairs = pd.read_csv(out / "a17_coherence_pairs.csv")
    tab, fams = learn(gates, classes, pairs, registry())
    tab.to_csv(out / "anchors_a17_relearned.csv", index=False)
    fams.to_csv(out / "a17_anchor_families.csv", index=False)
    cmp = compare(tab)
    cmp.to_csv(out / "anchors_vs_frozen.csv", index=False)
    print(cmp[["learned_identical", "fixed_identical"]].mean().round(3).to_dict(), len(cmp), "targets", flush=True)


if __name__ == "__main__":
    main()
