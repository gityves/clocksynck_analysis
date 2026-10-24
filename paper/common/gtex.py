"""GTEx bulk tissues as compartments of one donor (v10 counts, v8 time of death).

Controlled access (dbGaP phs000424): the time of death (`tod`), exact age, ventilation and the other protected
phenotypes stay on this machine. Compute scripts may hold them in memory and write aggregates (per tissue, per
group, fitted parameters, figures); they never write per-donor rows.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

import clocksync as cs
from common.paths import GTEX_PHENOTYPES

LABELLED = Path("/home/maxine/Documents/paychere/scRNAseq_database/phase_inference/baselines/data/gtex/"
                "gtex_v10_labelled.h5ad")
PHENOTYPES = GTEX_PHENOTYPES / "phs000424.v8.pht002742.v8.p2.c1.GTEx_Subject_Phenotypes.GRU.txt"
FAST = (1, 2)  # Hardy scale: violent and fast deaths
VENTILATOR = (0,)


AGE = {"20-29": 25, "30-39": 35, "40-49": 45, "50-59": 55, "60-69": 65, "70-79": 75}


def load(tissues=None, hardy=FAST, genes="clock", min_donors: int = 0, depth: str = "sum",
         require_tod: bool = True) -> tuple[pd.DataFrame, pd.DataFrame]:
    """(pseudobulk, donors): one row per (donor, tissue), cultured cells excluded, first sample of a duplicate kept.

    depth: "sum" (counts summed over all genes, as 01_gtex_organs and most analyses) or "library_size" (the
    annotated library size, as v2_gtex and the A46/A47 runs). min_donors: tissues with fewer donors dropped.
    donors: Hardy class, time of death `tod` (h, controlled), sex, age (bracket midpoint) and ischaemic time (h).
    """
    import anndata as ad

    a = ad.read_h5ad(LABELLED, backed="r")
    try:
        o = a.obs
        keep = o.DTHHRDY.isin(list(hardy)).to_numpy() & ~o.tissue.astype(str).str.startswith("Cells -").to_numpy()
        if require_tod:
            keep &= o.tod.notna().to_numpy()
        if tissues is not None:
            keep &= o.tissue.astype(str).isin(list(tissues)).to_numpy()
        idx = np.flatnonzero(keep)
        names = list(cs.gtex().genes) if genes == "clock" else list(genes) if genes != "all" else list(a.var_names)
        from clocksync.genes import gene_columns
        cols = gene_columns(a.var_names, names, missing="raise")
        Xa = a.X[idx]
        tot = np.asarray(Xa.sum(1)).ravel()
        X = Xa[:, cols]
        X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
        ob = o.iloc[idx]
    finally:
        a.file.close()
    pb = pd.DataFrame(X, columns=names)
    pb.insert(0, "donor", ob.subject.astype(str).to_numpy())
    pb.insert(1, "compartment", ob.tissue.astype(str).to_numpy())
    pb.insert(2, "depth", tot if depth == "sum" else ob.library_size.to_numpy(float))
    pb = pb.drop_duplicates(["donor", "compartment"]).reset_index(drop=True)
    if min_donors:
        n = pb.groupby("compartment").donor.nunique()
        pb = pb[pb.compartment.isin(n.index[n >= min_donors])].reset_index(drop=True)
    d = ob.assign(donor=ob.subject.astype(str)).drop_duplicates("donor").set_index("donor")
    donors = pd.DataFrame({"hardy": d.DTHHRDY.astype(int), "tod": d.tod.astype(float) % 24,
                           "male": (d.sex.astype(str) == "male").astype(float),
                           "age": d.age_bracket.astype(str).map(AGE).astype(float),
                           "ischaemic_h": d.ischemic_time.astype(float) / 60})
    return pb, donors


def phenotypes(columns=None) -> pd.DataFrame:
    """Protected subject phenotypes (dbGaP pht002742), indexed by donor. In memory only."""
    p = pd.read_csv(PHENOTYPES, sep="\t", comment="#", low_memory=False).set_index("SUBJID")
    p.index = p.index.astype(str)
    return p if columns is None else p[list(columns)]
