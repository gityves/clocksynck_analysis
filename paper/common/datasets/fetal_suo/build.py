"""Suo et al. 2022 (Science): developing human immune atlas, 25 fetuses, 9 organs (908,046 cells; includes the
Popescu 2019, Park 2020, Jardine 2021, Elmentaite 2021 and Stewart 2019 samples).

Source:   CELLxGENE collection dataset 1c47f6db, `$SCCC_DATA/natexp/suo2022/raw/1c47f6db.h5ad` (raw counts in `raw/X`).
Replaces: legacy/desync/validation/natural_experiments/scripts/build_fetal_organs.py (+ fetal_common.py), build_fetal.md.
Used by:  A107 (organ gate, macrophages and fibroblast/stromal per organ, per-organ coupling).

Choices (as the legacy builder)
- unit (donor) = fetus (`donor_id`); organ from `tissue` (`ORGAN`: liver, skin, spleen, thymus, bone_marrow, gut,
  yolk_sac, kidney, mesenteric_LN);
- excluded cells everywhere: author labels DOUBLET*, LOW_Q*, HIGH_MITO and `predicted_doublets`; so "all cells of an
  organ" means all retained cells;
- level "organ_all" (default): compartment = organ, all sort fractions; "organ_unsorted": only `Sort_id == TOT`
  libraries; "class": `organ:class`, the author `celltype_annotation` mapped to broad classes by `CLASSMAP` (first
  matching regex; unmatched labels form `other`). CD45P and CD45N fractions were sequenced separately, so their sum does
  not reproduce the tissue composition;
- a row sums all assays (3' v2, 5' v1, 5' v2) of its fetus x compartment; its assay is the majority one by cells; each
  fetus keeps the rows of its majority assay (most cells; the class `other` does not vote), which becomes the stratum;
- every pseudobulk with at least one cell is kept.

Donor table (fetus): `ga_pcw` (weeks; Carnegie stages converted, approximate), `development_stage`, `sex`, `n_cells`.
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
sys.path[:0] = [str(HERE.parents[1]), str(Path(__file__).resolve().parents[3])]
import clocksync as cs  # noqa: E402
from common.paths import DATA  # noqa: E402
from natexp_common import compare, majority_stratum, read_legacy, standard  # noqa: E402

ROOT = DATA / "natexp" / "suo2022"
H5AD = ROOT / "raw" / "1c47f6db.h5ad"
LEVELS = ("organ_all", "organ_unsorted", "class")
ORGAN = {"liver": "liver", "skin of body": "skin", "spleen": "spleen", "thymus": "thymus", "bone marrow": "bone_marrow",
         "gut wall": "gut", "yolk sac": "yolk_sac", "kidney": "kidney", "mesenteric lymph node": "mesenteric_LN"}
JUNK = r"DOUBLET|^LOW_Q|^HIGH_MITO"
CLASSMAP = [  # (class, regex); first match wins; unmatched -> other
    ("macrophage", r"^MACROPHAGE|^OSTEOCLAST$"),
    ("monocyte", r"^(MONOCYTE|PROMONOCYTE)"),
    ("DC", r"^(CYCLING_)?(DC|PDC)|^(DC1|DC2|MIGRATORY_DC|AS_DC|PRE_DC2|DC_PROGENITOR|LANGERHANS)"),
    ("endothelial", r"^ENDOTHELIUM"),
    ("fibroblast_stromal", r"FIBROBLAST|^VSMC|^SMOOTH_MUSCLE|^MESOTHELIUM|^YS_STROMA|^MESENCHYMAL|^OSTEOBLAST|^CHONDROCYTE"),
    ("epithelial", r"EPITHELIUM|^KERATINOCYTE|^ENTEROENDOCRINE"),
    ("hepatocyte", r"^HEPATOCYTE"),
    ("erythroid", r"(^|_)ERY$|^(MID|EARLY|LATE|YS)_ERY|^CYCLING_YS_ERY"),
    ("megakaryocyte", r"_MK$"),
    ("progenitor_HSC", r"^(HSC_MPP|CYCLING_MPP|CMP|GMP|MEMP|CYCLING_MEMP|MEP|LMPP_MLP|MOP|PROMYELOCYTE|MYELOCYTE)$"),
    ("B", r"(^|_)B\d?$|^B1$|PRE_B|PRO_B|^PLASMA_B|^MATURE_B|^IMMATURE_B|^CYCLING_B"),
    ("T", r"(^|_)T$|^TREG|^CD4\+T|^CD8|^ABT|INNATE_T|^DN|^DP"),
    ("NK", r"^(CYCLING_)?NK$"),
    ("ILC", r"ILC"),
]
CARNEGIE = {"Carnegie stage 13": 4.6, "Carnegie stage 19": 6.9, "Carnegie stage 20": 7.3, "Carnegie stage 22": 7.7,
            "Carnegie stage 23": 8.1}  # CS -> post-conception days / 7 (approximate)


def classify(label: str):
    if re.search(JUNK, label):
        return None
    return next((c for c, p in CLASSMAP if re.search(p, label)), "other")


def ga_pcw(stage: str) -> float:
    m = re.match(r"(\d+)(st|nd|rd|th) week", stage)
    return float(m.group(1)) if m else CARNEGIE.get(stage, np.nan)


def cell_labels() -> pd.DataFrame:
    a = ad.read_h5ad(H5AD, backed="r")
    o = a.obs[["donor_id", "tissue", "celltype_annotation", "Sort_id", "assay", "sex", "development_stage",
               "predicted_doublets"]].astype(str).copy()
    a.file.close()
    cmap = {l: classify(l) for l in o.celltype_annotation.unique()}
    o["cls"] = o.celltype_annotation.map(cmap)
    o["organ"] = o.tissue.map(ORGAN)
    o["ok"] = o.cls.notna() & o.organ.notna() & (o.predicted_doublets != "True")
    return o


def build(level: str = "organ_all", donors: list | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    """(pseudobulk, donors). `donors`: restrict to these fetuses (checks; the matrix is still read through)."""
    o = cell_labels()
    use = o.ok & (o.donor_id.isin(donors) if donors is not None else True)
    if level == "organ_unsorted":
        use &= o.Sort_id == "TOT"
    comp = o.organ if level != "class" else o.organ + ":" + o.cls.astype(str)
    lab = pd.DataFrame({"donor": o.donor_id, "compartment": comp.where(use)})
    pb = cs.pseudobulk(H5AD, donor="donor", compartment="compartment", labels=lab, layer="raw", min_cells=1,
                       gene_symbols="feature_name")
    d = o[use].assign(compartment=comp[use])
    n = d.groupby(["donor_id", "compartment", "assay"]).size().rename("n").reset_index()  # row assay: majority by cells
    row_assay = n.sort_values("n").drop_duplicates(["donor_id", "compartment"], keep="last").rename(
        columns={"donor_id": "donor"})[["donor", "compartment", "assay"]]
    pb = pb.merge(row_assay, on=["donor", "compartment"])
    pb = majority_stratum(pb, vote=~pb.compartment.str.endswith(":other"))
    f = o.drop_duplicates("donor_id").set_index("donor_id")
    don = pd.DataFrame({"development_stage": f.development_stage, "sex": f.sex,
                        "ga_pcw": f.development_stage.map(ga_pcw), "n_cells": o[o.ok].groupby("donor_id").size()})
    don = don.rename_axis("donor")
    return standard(pb), (don if donors is None else don.loc[don.index.isin(donors)])


def from_legacy(level: str = "organ_all") -> tuple[pd.DataFrame, pd.DataFrame]:
    """Standard table from the legacy `pseudobulk/pseudobulk_<level>/pseudobulk.h5ad` (class level: `pseudobulk_class`)."""
    pb = read_legacy(ROOT / "pseudobulk" / f"pseudobulk_{level}" / "pseudobulk.h5ad", extra=("assay", "ga_pcw", "sex"))
    don = pb.drop_duplicates("donor").set_index("donor")[["ga_pcw", "sex"]]
    pb = majority_stratum(pb, vote=~pb.compartment.str.endswith(":other")).drop(columns=["ga_pcw", "sex"])
    return standard(pb), don


if __name__ == "__main__":
    sl = ["F45", "F37", "F21"]
    for lev in LEVELS:
        a, da = build(lev, sl)
        b, db = from_legacy(lev)
        print(lev, compare(a, b, sl), "| strata", sorted(a.stratum.unique()), "| ga equal:",
              bool(np.allclose(da.ga_pcw.reindex(sl).dropna(), db.ga_pcw.reindex(sl).dropna(), equal_nan=True)))
