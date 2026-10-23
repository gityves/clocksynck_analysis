"""A105 (H105.3, D2, D3) and A105c (post hoc): does a longer ventilation loosen the lung epithelium? PMI drift.

Pre-registered in docs/PLAN_internal_desync.md (archive) (A105: H105.3, D2, D3; A105c is POST HOC-INFORMED, no verdict).
Estimand: epithelial coupling c on the endothelial + fibroblast host fitted on the pooled donors of the contrast;
d = c(IMV >= 9 days) - c(IMV <= 1 day) and d* at a common host spread, with the spread-free log ratios of z2 and psi.
Null: arm labels permuted with the host fixed (2000 permutations; within study in A105c); one-sided p for d < 0 is
`p_less`, with a donor-bootstrap interval (host fixed). Specificity companion on the same permutations: myeloid.
  H105.3   Melms COVID donors, 9 long vs 8 short (supported if d < 0, p_less < 0.05 and d* < 0); MDE of the epithelium.
  D3       Delorey lung classes (>= 20 nuclei), same arms, without the donors that are also in Melms (D25-D32).
  A105c    Melms + Delorey lungs pooled (classes harmonised to endothelial, epithelial, fibroblast, myeloid), study as
           stratum, labels permuted within study, bootstrap interval of d (500).
  D2       within Melms COVID donors, the change of the epithelial and myeloid clock with PMI, days from symptom onset
           and IMV days (z-scored): `Dataset.covariate` (host-locked amplitude, drift in h per SD, coupling slope).

    python compute_covid_ventilation.py [--subsample N] [--seed S]
Outputs: paper/results/figure_4/output/natural_experiments/a105_imv_contrast.csv (+ _companion, _mde), a105_delorey_lung_imv.csv,
         a105c_pooled.csv (+ _companion), a105_covariates.csv
Legacy:  legacy/desync/validation/natural_experiments/scripts/a105_covid.py (contrast, D2, D3), a105c_pooled.py.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from common import cli, data  # noqa: E402
from common.paths import out_dir  # noqa: E402
from results.figure_2.scripts.embryo import _shared as S  # noqa: E402

FOLDER = "natural_experiments"
ANCH = ["Endothelial cells", "Fibroblasts"]
EPI, MYE = "Epithelial cells", "Myeloid"
OVERLAP = {f"D{i}" for i in range(25, 33)}
HARMONISED = {"Endothelial cells": "endothelial", "Epithelial cells": "epithelial", "Fibroblasts": "fibroblast",
              "Myeloid": "myeloid", "lung:endothelial": "endothelial", "lung:epithelial": "epithelial",
              "lung:fibroblast": "fibroblast", "lung:myeloid": "myeloid"}
LEVELS = ("short", "long")
COLS = ["donors_short", "donors_long", "c_short", "c_long", "d", "lo", "hi", "p", "p_less", "d_star", "p_star",
        "log_z2_ratio", "log_psi_ratio", "status"]


def arms(imv: pd.Series) -> pd.Series:
    """IMV >= 9 days 'long', <= 1 day 'short' (others dropped)."""
    return pd.Series(np.where(imv >= 9, "long", np.where(imv <= 1, "short", None)), index=imv.index).dropna()


def write(c, out, name):
    c.table.to_csv(out / f"{name}.csv")
    if c.companion is not None:
        c.companion.to_csv(out / f"{name}_companion.csv", index=False)
    print(name); print(c.table[[x for x in COLS if x in c.table]].round(3).to_string(), flush=True)
    if c.companion is not None:
        print(c.companion.round(3).to_string(index=False), flush=True)


def main(argv=None):
    args = cli.parser(__doc__).parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    n_perm, n_boot = (200 if run["quick"] else 2000), run["n_boot"]
    melms, ms = data.load("covid_melms", subsample=args.subsample, seed=args.seed)
    covid = ms.index[ms.group != "Control"]
    lung = S.eligible(melms, 20, donors=covid)

    # H105.3
    lab = arms(ms.intubation_days.reindex(covid))
    lab = lab[lab.index.isin(lung.donor)]
    # the Dataset holds the labelled donors only: states are centred per compartment over the donors of the dataset
    # (legacy: `Study` of the labelled donors)
    ds = S.dataset(lung[lung.donor.isin(lab.index)].reset_index(drop=True), ms, anchors=ANCH)
    c = ds.contrast(lab, levels=LEVELS, targets=[EPI, MYE], companion=MYE, blocks="none", n_perm=n_perm,
                    n_boot=n_boot, seed=7)
    write(c, out, "a105_imv_contrast")
    mde = ds.contrast(lab, levels=LEVELS, targets=[EPI], blocks="none", n_perm=0, n_boot=0, mde=True,
                      n_mde=3 if run["quick"] else 100, seed=70)
    mde.mde.rename_axis("compartment").to_csv(out / "a105_mde.csv")
    print(mde.mde.round(3).to_string(), flush=True)

    # D3: Delorey lung classes without the Melms-overlap donors
    cls, cl = data.load("covid_delorey", level="class", subsample=args.subsample, seed=args.seed)
    lc = S.eligible(cls, 20)
    lc = lc[lc.compartment.str.startswith("lung:") & ~lc.donor.isin(OVERLAP)].reset_index(drop=True)
    lab3 = arms(cl.imv_days)
    lab3 = lab3[lab3.index.isin(lc.donor)]
    print("D3 Delorey lung donors by arm:", lab3.value_counts().to_dict(), flush=True)
    if len(lab3.value_counts()) == 2 and lab3.value_counts().min() >= (2 if run["quick"] else 4):
        tg = [k for k in ("lung:epithelial", "lung:myeloid") if k in set(lc.compartment)]
        ds3 = S.dataset(lc[lc.donor.isin(lab3.index)].reset_index(drop=True), cl, anchors=["lung:endothelial",
                                                                                         "lung:fibroblast"])
        write(ds3.contrast(lab3, levels=LEVELS, targets=tg, blocks="none", n_perm=n_perm, n_boot=n_boot, seed=9),
              out, "a105_delorey_lung_imv")

    # A105c (post hoc): pooled lungs, study as stratum
    pm = S.eligible(melms, 20, donors=covid, compartments=[k for k in HARMONISED if k in set(melms.compartment)])
    pd_ = lc.copy()
    pd_ = pd_[pd_.compartment.isin(HARMONISED)]
    pool = pd.concat([pm.assign(stratum="melms"), pd_.assign(stratum="delorey")], ignore_index=True)
    pool["compartment"] = "lung:" + pool.compartment.map(HARMONISED)
    imv = pd.concat([ms.intubation_days.reindex(covid), cl.imv_days])
    labp = arms(imv)
    labp = labp[labp.index.isin(pool.donor)]
    pool = pool[pool.donor.isin(labp.index)].reset_index(drop=True)
    don = pd.DataFrame(index=pd.Index(sorted(pool.donor.unique()), name="donor"))
    print("A105c arms by study:", pd.crosstab(pool.drop_duplicates("donor").set_index("donor").stratum.reindex(labp.index),
                                              labp).to_dict(), flush=True)
    dsp = S.dataset(pool, don, anchors=["lung:endothelial", "lung:fibroblast"])
    write(dsp.contrast(labp, levels=LEVELS, targets=["lung:epithelial", "lung:myeloid"], companion="lung:myeloid",
                       n_perm=n_perm, n_boot=500 if not run["quick"] else n_boot, seed=13), out, "a105c_pooled")

    # D2: covariate model within COVID donors
    ds = S.dataset(lung, ms, anchors=ANCH)
    rows = []
    for var in ("pmi_h", "interval_death_symptoms_onset_days", "intubation_days"):
        x = ms[var].reindex(covid).dropna()
        x = (x - x.mean()) / x.std()
        r = ds.covariate(x, targets=[EPI, MYE], center=False, n_null=run["n_null"], seed=0)
        rows.append(r.reset_index().assign(variable=var))
    cov = pd.concat(rows, ignore_index=True)
    cov.to_csv(out / "a105_covariates.csv", index=False)
    print(cov.round(3).to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
