"""Standard tables from the pickles of the A17 re-analysis (`$SCCC_DATA/a17/tables/<target>.pkl`, written by the
legacy `desync/common/a17_tables.py`; our own files, not downloads).

Each pickle holds the exact input of one healthy fit or contrast: the donor x compartment table, the case/control
labels, the reference compartment, the permutation blocks. Targets that read the same atlas and the same samples are
one standard table here: the union of their rows (identical where they overlap; checked), with one donor-table
column per target: `<target>` = "case" / "control" / NaN, `healthy` = donor of the healthy fit, and
`block_<target>` = the permutation block when it is not the stratum.
"""

from __future__ import annotations

import pickle

import numpy as np
import pandas as pd

from census import ATLAS, GENES, DATA, standard  # noqa: F401

TABLES = DATA / "a17" / "tables"
HEAD = ["resection", "compartment", "stratum", "n_cells", "depth"]


def read(name: str) -> dict:
    with open(TABLES / f"{name}.pkl", "rb") as f:
        return pickle.load(f)


def from_pickles(names, healthy=None, prefix_strata: dict | None = None, covariates: pd.DataFrame | None = None):
    """(pb, donors) from the pickles `names` (union of rows). `healthy`: name of the healthy-fit pickle whose donors
    get `healthy = True` (default: every pickle of kind "healthy" among `names`)."""
    recs = {n: read(n) for n in names}
    parts, extra = [], {}
    for n, r in recs.items():
        t = r["table"].copy()
        t["stratum"] = t["stratum"].astype(str)
        if prefix_strata and n in prefix_strata:
            t["stratum"] = t.resection.str.split(":").str[0].map(prefix_strata[n]) + t.stratum
        parts.append(t[HEAD + GENES])
    pb = parts[0]
    for t in parts[1:]:  # union; rows present in both must agree
        a, b = pb.set_index(["resection", "compartment"]), t.set_index(["resection", "compartment"])
        common = a.index.intersection(b.index)
        num = ["n_cells", "depth"] + GENES
        if len(common) and np.abs(a.loc[common, num].to_numpy(float) - b.loc[common, num].to_numpy(float)).max() > 0:
            raise ValueError("pickles disagree on shared rows")
        if len(common) and (a.loc[common, "stratum"] != b.loc[common, "stratum"]).any():
            raise ValueError("pickles disagree on strata")
        pb = pd.concat([pb, t.set_index(["resection", "compartment"]).drop(common).reset_index()], ignore_index=True)
    pb = standard(pb)
    donors = pd.DataFrame(index=pd.Index(sorted(pb.donor.unique()), name="donor"))
    if covariates is not None:
        donors = donors.join(covariates)
    heal = [n for n, r in recs.items() if r["kind"] == "healthy"] if healthy is None else [healthy]
    if heal:
        members = set().union(*[set(recs[n]["table"].resection) for n in heal])
        donors["healthy"] = donors.index.isin(members)
    for n, r in recs.items():
        if r["kind"] == "healthy":
            continue
        lab = r["labels"]
        donors[n] = lab.reindex(donors.index).to_numpy() if lab is not None else np.nan
        w = r.get("within")
        if isinstance(w, pd.Series):
            strat = r["table"].groupby("resection").stratum.first().astype(str)
            if not w.reindex(strat.index).astype(str).equals(strat):  # a block other than the stratum
                donors[f"block_{n}"] = w.reindex(donors.index).to_numpy()
    return pb, donors


def meta(name: str) -> dict:
    """Metadata of one pickle (no table): kind, study, tissue, label, ref, groups, within, primary."""
    r = read(name)
    return {k: v for k, v in r.items() if k not in ("table", "labels")}
