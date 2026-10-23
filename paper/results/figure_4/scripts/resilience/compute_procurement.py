"""A19 part 1: procurement check of the fibrotic-lung T-cell result (HLCA).

Pre-registered in docs/PLAN_internal_desync.md (archive) (A19 part 1). HLCA fibrotic cases are mostly explants and controls
mostly ventilated organ-donor lungs. On the A17 table and learned anchors of `hlca_fibrotic` (fixed rule as
sensitivity), one anchored host on all its donors, and the contrasts
  P1a: healthy living-patient / autopsy tissue (surgical resection, biopsy, autopsy) minus healthy donor lungs;
  P1b: fibrotic explants minus healthy living-patient / autopsy tissue (the disease contrast without donor-lung
       controls).
Per compartment d, d* (common host spread), labels permuted freely (as the original contrast; 1,000 permutations
shared by all compartments), donor-bootstrap interval of d; companion lung:fibroblast on the same permutations;
lymphocyte module (mean d over the lymphocyte compartments, A18 addendum 1 statistic 1). MDE of lung:T with --mde.
Reading (A19): P1b negative (p < 0.05) and P1a not -> the T-cell loss is not explained by procurement; P1a negative
with the sign of the original contrast -> procurement can produce it; both null -> inconclusive, with the MDE.

    python compute_procurement.py [--subsample N] [--seed S] [--mde]
Outputs (paper/results/figure_4/output/resilience/): a19_procurement.csv, a19_procurement_companion.csv,
a19_procurement_module.csv, a19_procurement_mde.csv (--mde)
Legacy: legacy/desync/common/a19_procurement.py.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from common import cli, data  # noqa: E402
from common.paths import out_dir  # noqa: E402

import _shared as sh  # noqa: E402

FOLDER = "resilience"
LIVING = ("surgical_resection", "biopsy", "autopsy")
sys.path.insert(0, str(Path(__file__).resolve().parents[4]))  # paper/
from common.paths import DATASETS  # noqa: E402
DONORS = DATASETS / "_atlas" / "inputs" / "lung_donors.csv"  # HLCA donor table of the builders
COMPANION = "lung:fibroblast"


def procurement(index) -> pd.Series:
    """Donor -> 'donor_lung', 'explant', 'living' or 'other' (HLCA `tissue_sampling_method_major`)."""
    m = pd.read_csv(DONORS).set_index("donor_id").tissue_sampling_method_major.astype(str)
    raw = pd.Series(m.reindex([i.split(":", 1)[-1] for i in index]).fillna("").to_numpy())
    out = np.where(raw.str.startswith("donor_lung"), "donor_lung",
                   np.where(raw.str.startswith("lung_explant"), "explant",
                            np.where(raw.str.replace("_cryopreserved", "").isin(LIVING), "living", "other")))
    return pd.Series(out, index=index)


def arms(lab: pd.Series) -> dict:
    p = procurement(lab.index)
    healthy, fibrotic = lab == "control", lab == "case"
    a = pd.concat([pd.Series("control", index=lab.index[healthy & (p == "donor_lung")]),
                   pd.Series("case", index=lab.index[healthy & (p == "living")])])
    b = pd.concat([pd.Series("control", index=lab.index[healthy & (p == "living")]),
                   pd.Series("case", index=lab.index[fibrotic & (p == "explant")])])
    return {"P1a_living_vs_donor_lung": a, "P1b_fibrotic_explant_vs_living": b}


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--mde", action="store_true", help="also the MDE of lung:T (slow)")
    args = ap.parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    fam = sh.family_map()
    pb, donors = data.load("hlca_fibrotic")
    lab = donors["hlca_fibrotic"].dropna()
    lab = lab[lab.isin(sh.GROUPS)]
    rows, comps, mods, mdes = [], [], [], []
    for rule in sh.RULES:
        anc = sh.anchor_set("hlca_fibrotic", rule)
        ds = cs.Dataset(pb, anchors=anc, donors=donors, name="hlca_fibrotic")
        for name, cl in arms(lab).items():
            cl = sh.subsample_donors(cl, args.subsample, args.seed)
            kw = dict(levels=sh.GROUPS, blocks="none", n_perm=run["n_perm"], n_boot=run["n_boot"], seed=args.seed)
            c = ds.contrast(cl, companion=COMPANION, **kw)
            t = c.table.reset_index().assign(contrast=name, rule=rule)
            lym = t[t.compartment.map(fam).eq("lymphocyte") & t.d.notna()]
            mods.append(dict(contrast=name, rule=rule, compartments=len(lym), module_d=lym.d.mean(),
                             module_d_star=lym.d_star.mean()))
            rows.append(t)
            if c.companion is not None:
                comps.append(c.companion.assign(contrast=name, rule=rule))
            if args.mde and rule == "learned":
                m = ds.contrast(cl, targets=["lung:T"], levels=sh.GROUPS, blocks="none", n_perm=0, n_boot=0,
                                seed=args.seed, mde=True, n_mde=10 if run["quick"] else 100).mde.iloc[0]
                mdes.append(dict(contrast=name, rule=rule, compartment="lung:T", mde=m["mde"], c_pooled=m["c_pooled"]))
            print(name, rule, (cl == "control").sum(), (cl == "case").sum(), flush=True)
    pd.concat(rows, ignore_index=True).to_csv(out / "a19_procurement.csv", index=False)
    pd.concat(comps, ignore_index=True).to_csv(out / "a19_procurement_companion.csv", index=False)
    pd.DataFrame(mods).to_csv(out / "a19_procurement_module.csv", index=False)
    if mdes:
        pd.DataFrame(mdes).to_csv(out / "a19_procurement_mde.csv", index=False)
    t = pd.concat(rows)
    print(t[t.compartment.isin(["lung:T", COMPANION])][["contrast", "rule", "compartment", "d", "p", "d_star",
                                                         "p_star"]].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
