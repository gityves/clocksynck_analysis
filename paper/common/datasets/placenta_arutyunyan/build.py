"""Arutyunyan et al. 2023 (Nature): first-trimester maternal-fetal interface (decidua, placenta, maternal blood).

Source:   CELLxGENE collection dataset 7a4dd01b (E-MTAB-12421), `$SCCC_DATA/natexp/arutyunyan2023/raw/7a4dd01b.h5ad`
          (raw counts in `raw/X`, CSR; gene symbols in `raw/var.feature_name`).
Replaces: legacy/desync/validation/natural_experiments/scripts/build_maternal_fetal.py (+ fetal_common.py), build_fetal.md.
Used by:  A106 (maternal gate, fetal x maternal cross gate, ambient-RNA check).

Choices (as the legacy builder)
- unit (donor) = pregnancy: `donor_id` without the `_fetus` / `_mother` suffix; the origin of a cell is the genotype
  call `origin_M_F` (M mother, F fetus); cells with another origin (`soc_didnt_work`, ambiguous or doublet) are dropped;
- compartment = origin x class (level "origin", default): `fet:<class>`, `mat:<class>` (decidua and placenta libraries
  summed), `matblood:<class>` (maternal blood libraries). Classes: the author `coarse_annot` (genotype consistent) and,
  to split trophoblast (VCT, SCT, EVT) and decidual stroma (perivascular PV1/PV2), the author fine label
  `celltype_predictions` inside a coarse class only; `<origin>:other` collects unmapped classes. Level "tissue":
  compartment `<tissue>:<origin>:<class>` (tissue of the library, decidua basalis merged into decidua), so maternal and
  fetal compartments can be taken from different libraries (ambient RNA);
- rows are pregnancy x compartment x assay (sc / sn / multiome); each pregnancy keeps the rows of its majority
  assay (most cells, over the classes other than `:other`), which becomes the stratum;
- extra columns in the table, not clock genes: `XIST` and `Y_total` (summed counts of XIST and of the Y-linked sex
  genes) for the ambient-RNA estimate of A106; every pseudobulk with at least one cell is kept.

Donor table (pregnancy): `ga_pcw`, `dev_age`, libraries, `n_fetal_cells`, `n_mat_cells`, `n_matblood_cells`,
`fetal_y_per_million` and the inferred fetal sex (Y counts per million fetal counts > 5: male; post hoc, see build_fetal.md).
"""

from __future__ import annotations

import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
sys.path[:0] = [str(HERE.parents[1]), str(Path(__file__).resolve().parents[3])]
import clocksync as cs  # noqa: E402
from common.paths import DATA  # noqa: E402
from natexp_common import GENES, compare, majority_stratum, read_legacy, standard  # noqa: E402

ROOT = DATA / "natexp" / "arutyunyan2023"
H5AD = ROOT / "raw" / "7a4dd01b.h5ad"
LEVELS = ("origin", "tissue")
Y_GENES = ["RPS4Y1", "DDX3Y", "UTY", "KDM5D", "EIF1AY", "ZFY", "USP9Y", "TXLNGY", "NLGN4Y", "TBL1Y"]
ASSAY = {"10X_3'_scRNA-seq": "10x 3' sc", "10X_3'_snRNA-seq": "10x 3' sn", "10X_multiome": "10x multiome"}
TROPH = {"1_VCT": "VCT", "9_VCT": "VCT", "20_SCT": "SCT", "13_EVT": "EVT", "16_EVT": "EVT"}
COARSE = {("F", "Trophoblast"): None, ("F", "fF"): "fibroblast", ("F", "Endothelial"): "endothelial",
          ("F", "Myeloid"): "Hofbauer",
          ("M", "dS_uSMC"): None, ("M", "Endothelial"): "endothelial", ("M", "Epithelial"): "epithelial",
          ("M", "NK"): "dNK", ("M", "Myeloid"): "myeloid", ("M", "T"): "T"}
BLOOD = {"Myeloid": "myeloid", "NK": "NK", "T": "T"}
GA = {"4-5_PCW": 4.5, "7-8_PCW": 7.5, "8-9_PCW": 8.5, "12-13_PCW": 12.5}


def compartment(origin, coarse, pred, tissue):
    if origin == "M" and tissue == "blood":
        c = BLOOD.get(coarse)
        return f"matblood:{c}" if c else None
    pre = "fet" if origin == "F" else "mat"
    if (origin, coarse) == ("F", "Trophoblast"):
        return f"fet:{TROPH.get(pred, 'other')}"
    if (origin, coarse) == ("M", "dS_uSMC"):
        return "mat:perivascular" if pred in ("14.1_PV1", "14.2_PV2") else "mat:dS"
    c = COARSE.get((origin, coarse))
    return f"{pre}:{c}" if c else f"{pre}:other"


def cell_labels() -> pd.DataFrame:
    """Per cell: pregnancy, compartment (`comp`), library tissue, assay and the key columns of the two levels."""
    a = ad.read_h5ad(H5AD, backed="r")
    obs = a.obs[["donor_id", "origin_M_F", "coarse_annot", "celltype_predictions", "tissue", "technology", "dev_age",
                 "sample", "dataset"]].astype(str).copy()
    a.file.close()
    o = obs
    o["tissue"] = o.tissue.replace({"decidua basalis": "decidua"})
    o["donor"] = o.donor_id.str.replace(r"_(fetus|mother)$", "", regex=True)
    o["comp"] = [compartment(a_, b, c, d) if a_ in ("M", "F") else None for a_, b, c, d in
                 zip(o.origin_M_F, o.coarse_annot, o.celltype_predictions, o.tissue)]
    o["assay"] = o.technology.map(ASSAY).astype(str)
    return o


def y_present() -> list:
    """Y-linked genes of the list that exist in the matrix (as the legacy builder: absent ones are skipped)."""
    import h5py
    from anndata.io import read_elem
    with h5py.File(H5AD, "r") as f:
        names = set(read_elem(f["raw/var"])["feature_name"].astype(str))
    return [g for g in Y_GENES if g in names]


def _pb(labels: pd.DataFrame, level: str, donors) -> pd.DataFrame:
    y = y_present()
    ok = labels.comp.notna()
    comp = labels.comp if level == "origin" else labels.tissue + ":" + labels.comp.astype(str)
    lab = pd.DataFrame({"donor": labels.donor, "compartment": (comp + "|" + labels.assay).where(ok)})
    if donors is not None:
        lab.loc[~lab.donor.isin(donors), "compartment"] = np.nan
    fine = cs.pseudobulk(H5AD, donor="donor", compartment="compartment", labels=lab, layer="raw",
                         gene_symbols="feature_name", genes=GENES + ["XIST"] + y, min_cells=1)
    fine[["compartment", "assay"]] = fine.compartment.str.rsplit("|", n=1, expand=True)
    fine["Y_total"] = fine[y].sum(axis=1)
    return fine.drop(columns=y)


def donor_table(labels: pd.DataFrame, fine: pd.DataFrame) -> pd.DataFrame:
    g = labels.groupby("donor")
    pg = pd.DataFrame({"dev_age": g.dev_age.first(), "samples": g["sample"].agg(lambda x: ";".join(sorted(set(x)))),
                       "dataset": g.dataset.agg(lambda x: ";".join(sorted(set(x)))),
                       "technologies": g.assay.agg(lambda x: ";".join(sorted(set(x[x != "nan"]))))})
    gnum = lambda s: GA.get(s, float(s.split("_")[0]) if "-" not in s else np.nan)  # noqa: E731
    pg["ga_pcw"] = pg.dev_age.map(gnum)
    origin = fine.compartment.str.split(":").map(lambda p: p[-2] if len(p) == 3 else p[0])
    fet = fine[origin == "fet"].groupby("donor")[["n_cells", "depth", "Y_total"]].sum()
    pg["n_fetal_cells"] = fet.n_cells.reindex(pg.index).fillna(0).astype(int)
    pg["fetal_y_per_million"] = (1e6 * fet.Y_total / fet.depth).reindex(pg.index)
    for name in ("mat", "matblood"):
        pg[f"n_{name}_cells"] = fine[origin == name].groupby("donor").n_cells.sum().reindex(pg.index).fillna(0).astype(int)
    pg["fetal_sex_inferred"] = np.where(pg.fetal_y_per_million.isna() | (pg.n_fetal_cells < 500), "unknown",
                                        np.where(pg.fetal_y_per_million > 5, "male", "female"))
    return pg.rename_axis("donor")


def build(level: str = "origin", donors: list | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    """(pseudobulk, donors). `donors`: restrict to these pregnancies (checks; the matrix is still read through)."""
    labels = cell_labels()
    fine = _pb(labels, level, donors)
    don = donor_table(labels, fine)
    pb = majority_stratum(fine, vote=~fine.compartment.str.endswith(":other"))
    return standard(pb, extra=("XIST", "Y_total")), (don if donors is None else don.loc[don.index.isin(donors)])


def from_legacy(level: str = "origin") -> tuple[pd.DataFrame, pd.DataFrame]:
    """Standard table from the legacy `pseudobulk[_tissue]/pseudobulk.h5ad` (+ `ambient_genes.csv`, `samples.csv`)."""
    d = ROOT / ("pseudobulk" if level == "origin" else "pseudobulk_tissue")
    pb = read_legacy(d / "pseudobulk.h5ad", extra=("assay",))
    amb = pd.read_csv(d / "ambient_genes.csv", dtype={"donor_id": str}).rename(
        columns={"donor_id": "donor", "cell_type": "compartment"})
    pb = pb.merge(amb[["donor", "compartment", "assay", "XIST", "Y_total"]], on=["donor", "compartment", "assay"])
    pb = majority_stratum(pb, vote=~pb.compartment.str.endswith(":other"))
    don = pd.read_csv(ROOT / "pseudobulk" / "samples.csv", dtype={"pregnancy": str}).set_index("pregnancy").rename_axis("donor")
    return standard(pb, extra=("XIST", "Y_total")), don


if __name__ == "__main__":
    sl = ["6", "10", "Hrv43"]  # Hrv43: more than one assay
    for lev in LEVELS:
        a, da = build(lev, sl)
        b, db = from_legacy(lev)
        m = a.merge(b, on=["donor", "compartment"])
        print("XIST/Y_total max diff", float(np.abs(m.XIST_x - m.XIST_y).max()), float(np.abs(m.Y_total_x - m.Y_total_y).max()),
              "strata", sorted(a.stratum.unique()))
        print(lev, compare(a, b, sl), "| donors table equal:",
              da[["ga_pcw", "n_fetal_cells", "fetal_sex_inferred"]].round(3).astype(str).equals(
                  db.loc[da.index, ["ga_pcw", "n_fetal_cells", "fetal_sex_inferred"]].round(3).astype(str)))
