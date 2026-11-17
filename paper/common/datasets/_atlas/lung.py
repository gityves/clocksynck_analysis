"""Lung atlases: HLCA full (healthy and IPF; fibrotic vs healthy with one T class) and Natri et al. 2024 ILD (fibrotic vs
control lung; non-cycling T control): shared builder of `lung_hlca`, `hlca_fibrotic`, `lung_natri`,
`lung_natri_control`. Amendments A4, A5, A6, A12.

Legacy: desync/06_lung/scripts/lung_atlas.py (`_groups`, `load`), desync/07_lung_ild/scripts/ild.py (`_groups`, `load`,
`labels`), desync/08_lung_t_residency/scripts/residency.py (`groups`, `load`, `donor_sets`), as captured by
desync/common/a17_tables.py: A17 targets `lung_hlca` (healthy), `lung_ipf_hlca`, `hlca_fibrotic`, `ild_natri`,
`lung_natri_control`.
Raw inputs (group pseudobulks): $SCCC_DATA/atlas/hlca_full, hlca_tcell, natri2024_ild, natri2024_tcell; donor labels that
are per cell in the Census (`lung_condition`, `study_major`, `tissue_sampling_method_major`) from `inputs/lung_donors.csv`
(desync/06_lung/scripts/lung_donors.py reads the source obs over S3).

Choices.
- HLCA (`lung_hlca`): lung-tissue sampling methods (not nose); classes from `ann_level_3` (`L3`; smooth muscle = all of
  `ann_level_2` "Smooth muscle"); majority assay per donor; donors with a single `lung_condition`; stratum = majority
  study | assay; pseudobulks >= 50 cells. Union of two legacy tables, each with its own donor rule and its own `>= 20
  donors` class rule: the healthy fit (conditions Healthy, Healthy (tumor adjacent); `healthy`) and IPF vs healthy in the
  four studies with both (`lung_ipf_hlca`, cases = IPF; permutation block = study, column `block_lung_ipf_hlca`).
- `hlca_fibrotic` (hlca_tcell build, one T class "lung:T"): the same rules with healthy donors (control) and donors with
  a fibrotic condition (`FIBROTIC`, case); all classes with >= 50-cell pseudobulks kept.
- Natri (`lung_natri`): primary donors; classes from `annotation_level_1` (`L1`; fibroblast also from the lineages of
  `FIB_LINEAGES`); all cells of a donor summed whatever the assay, stratum = majority processing site; pseudobulks >= 20
  cells; control = Status "Control", case = "Disease" except sarcoidosis. `lung_natri_control`: control donors only, from
  natri2024_tcell with the T class replaced by non-cycling T cells (`T_noncycling`), classes with >= 20 donors.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

import census as C
import legacy as L
import organ as O

HLCA_TISSUE = ("donor_lung", "donor_lung_cryopreserved", "lung_explant", "surgical_resection", "autopsy",
               "autopsy_cryopreserved", "biopsy", "biopsy_cryopreserved")
L3 = {"AT1": "AT1", "AT2": "AT2", "Multiciliated lineage": "multiciliated", "Secretory": "secretory", "Basal": "basal",
      "EC capillary": "EC_capillary", "EC venous": "EC_venous", "EC arterial": "EC_arterial",
      "Lymphatic EC mature": "lymphatic_EC", "Fibroblasts": "fibroblast", "Macrophages": "macrophage",
      "Monocytes": "monocyte", "Dendritic cells": "DC", "Mast cells": "mast", "T cell lineage": "T",
      "Innate lymphoid cell NK": "NK_ILC", "B cell lineage": "B"}
L1 = {"AT1": "AT1", "AT2": "AT2", "Ciliated": "ciliated", "Secretory": "secretory", "Basal": "basal",
      "Arteriole": "EC_arteriole", "Capillary": "EC_capillary", "Venous": "EC_venous", "Lymphatic": "lymphatic_EC",
      "Fibroblast": "fibroblast", "Macrophage": "macrophage", "Monocyte": "monocyte", "Dendritic Cells": "DC",
      "Dendritic cells": "DC", "Myeloid": "mast", "T cells": "T", "Lymphoid": "NK", "Bcells": "B"}
FIB_LINEAGES = ("PLIN2+ FB", "MyoFB", "MyoFB - Activated")
HEALTHY = ("Healthy", "Healthy (tumor adjacent)")
FIBROTIC = ("IPF", "Systemic sclerosis-associated ILD", "HP", "NSIP", "End-stage lung fibrosis, unknown etiology",
            "Myositis-associated ILD", "ILD")
IPF_STUDIES = ("Banovich_Kropski_2020", "Kaminski_2020", "Misharin_Budinger_2018", "Sheppard_2020")
LEGACY = {"lung_hlca": ["lung_hlca", "lung_ipf_hlca"], "hlca_fibrotic": ["hlca_fibrotic"], "lung_natri": ["ild_natri"],
          "lung_natri_control": ["lung_natri_control"]}
INPUTS = C.HERE / "inputs"


def hlca_donors() -> pd.DataFrame:
    d = pd.read_csv(INPUTS / "lung_donors.csv").set_index("donor_id")
    d.index = "hlca_full:" + d.index.astype(str)
    return d


def _hlca_cov() -> pd.DataFrame:
    """Census donors.csv of hlca_full plus the per-donor labels of lung_donors.csv (age filled from the HLCA field)."""
    c = C.census_donors(C.ATLAS / "hlca_full", "hlca_full")
    ld = hlca_donors()
    c["age"] = c.age.fillna(pd.to_numeric(ld.age_or_mean_of_age_range, errors="coerce").reindex(c.index))
    for k in ("lung_condition", "study_major", "tissue_sampling_method_major"):
        c[k] = ld[k].reindex(c.index)
    return c


def _hlca_frame(dataset: str, donors, genes, tcell: bool = False) -> pd.DataFrame:
    o, X = C.read(C.ATLAS / dataset, genes, donors)
    cls = o.ann_level_3.astype(str).map(L3)
    cls[o.ann_level_2.astype(str) == "Smooth muscle"] = "smooth_muscle"
    keep = o.tissue_sampling_method.astype(str).isin(HLCA_TISSUE).to_numpy() & (o.tissue.astype(str) != "nose").to_numpy()
    F = C.frame(o, X, "hlca_full", "lung:" + cls.astype(str), o.assay, keep=keep & cls.notna().to_numpy())
    return C.keep_majority(F)


def _natri_cls(o) -> pd.Series:
    cls = o.annotation_level_1.astype(str).map(L1)
    cls[o.lineage.astype(str).isin(FIB_LINEAGES)] = "fibroblast"
    return cls


def _natri_labels(d: pd.DataFrame) -> pd.Series:
    prim = d["is_primary_data_h5ad" if "is_primary_data_h5ad" in d else "is_primary_data"].astype(str) == "True"
    d = d[prim.to_numpy()]
    case = (d.Status == "Disease") & (d.Diagnosis != "Sarcoidosis")
    lab = pd.Series(np.where(d.Status == "Control", "control", np.where(case, "case", "")), index=d.index)
    return lab[lab != ""]


def _min_donors(pb: pd.DataFrame, keep, n: int) -> pd.DataFrame:
    sub = pb[pb.resection.isin(keep)]
    c = sub.groupby("compartment").resection.nunique()
    return sub[sub.compartment.isin(c[c >= n].index)]


def build(name: str, donors=None, genes="clock"):
    full = donors is None
    if name == "lung_hlca":
        F = _hlca_frame("hlca_full", donors, genes)
        ld = hlca_donors()
        ok = ld.index[~ld.lung_condition.str.contains(";")]
        pb = C.donor_table(F, ok, 50, stratum=None)
        pb["stratum"] = pb.resection.map(ld.study_major) + "|" + _assay_of(F, pb)
        h = ld.index[ld.lung_condition.isin(HEALTHY) & ~ld.lung_condition.str.contains(";")]
        i = ld.index[ld.lung_condition.isin(HEALTHY + ("IPF",)) & ld.study_major.isin(IPF_STUDIES)
                     & ~ld.lung_condition.str.contains(";")]
        th, ti = _min_donors(pb, h, 20 if full else 1), _min_donors(pb, i, 20 if full else 1)
        out = pd.concat([th, ti]).drop_duplicates(["resection", "compartment"]).reset_index(drop=True)
        lab = pd.Series(np.where(ld.lung_condition == "IPF", "case", "control"), index=ld.index).reindex(i)
        extra = {"lung_ipf_hlca": lab.rename("lung_ipf_hlca"), "healthy": O.flag(th.resection.unique()),
                 "block_lung_ipf_hlca": ld.study_major.reindex(i)}
        return O.finish(out, _hlca_cov(), extra, cols=("age", "sex", "disease", "lung_condition", "study_major",
                                                      "tissue_sampling_method_major"))
    if name == "hlca_fibrotic":
        F = _hlca_frame("hlca_tcell", donors, genes)  # T class is one compartment "lung:T"
        ld = hlca_donors()
        ld = ld[~ld.lung_condition.str.contains(";")]
        h, f = ld.index[ld.lung_condition.isin(HEALTHY)], ld.index[ld.lung_condition.isin(FIBROTIC)]
        lab = pd.concat([pd.Series("control", index=h), pd.Series("case", index=f)])
        pb = C.donor_table(F, lab.index, 50, stratum=None)
        pb["stratum"] = pb.resection.map(ld.study_major) + "|" + _assay_of(F, pb)
        return O.finish(pb, _hlca_cov(), {"hlca_fibrotic": lab.rename("hlca_fibrotic")},
                        cols=("age", "sex", "disease", "lung_condition", "study_major"))
    if name == "lung_natri":
        d = C.census_donors(C.ATLAS / "natri2024_ild", "natri2024_ild")
        lab = _natri_labels(d)
        root = C.ATLAS_WG if isinstance(genes, str) and genes == "all" else C.ATLAS  # all genes: the genome-wide build
        o, X = C.read(root / "natri2024_ild", genes, donors)
        cls = _natri_cls(o)
        F = C.frame(o, X, "natri2024_ild", "lung:" + cls.astype(str), o.Processing_site, keep=cls.notna().to_numpy())
        if donors is not None:
            lab = lab[lab.index.isin(F.resection.unique())]
        pb = C.donor_table(F, lab.index, 20, stratum=C.majority(F[F.resection.isin(lab.index)]))
        return O.finish(pb, d, {"ild_natri": lab.rename("ild_natri")}, cols=("age", "sex", "disease", "Diagnosis", "Status"))
    if name == "lung_natri_control":
        d = C.census_donors(C.ATLAS / "natri2024_ild", "natri2024_ild")
        lab = _natri_labels(d)
        ctrl = lab.index[lab == "control"]
        o, X = C.read(C.ATLAS / "natri2024_tcell", genes, donors)
        cls = _natri_cls(o)
        is_t = cls == "T"
        comp = cls.where(~is_t, "T_noncycling").where(~is_t | (o.t_cycling.astype(str) == "noncycling"))
        F = C.frame(o, X, "natri2024_ild", "lung:" + comp.astype(str), o.Processing_site, keep=comp.notna().to_numpy())
        stratum = C.majority(F)
        pb = C.donor_table(F, ctrl, 20, stratum=stratum, min_donors=None if not full else 20)
        return O.finish(pb, d, {"healthy": O.flag(pb.resection.unique())}, cols=("age", "sex", "disease", "Diagnosis", "Status"))
    raise KeyError(name)


def _assay_of(F: pd.DataFrame, pb: pd.DataFrame) -> pd.Series:
    """Majority assay of each pb row's donor (each donor has one assay after the majority step)."""
    return pb.resection.map(F.drop_duplicates("resection").set_index("resection").batch)


def from_legacy(name: str):
    if name == "lung_hlca":
        d = _hlca_cov()
        cols = ("age", "sex", "disease", "lung_condition", "study_major", "tissue_sampling_method_major")
    elif name == "hlca_fibrotic":
        d = _hlca_cov()
        cols = ("age", "sex", "disease", "lung_condition", "study_major")
    else:
        d = C.census_donors(C.ATLAS / "natri2024_ild", "natri2024_ild")
        cols = ("age", "sex", "disease", "Diagnosis", "Status")
    return L.from_pickles(LEGACY[name], covariates=C.donors_frame(d, d.index, cols=cols))


def check_legacy(name: str, pb: pd.DataFrame, donors: pd.DataFrame) -> dict:
    """Donor rules recomputed from the legacy table (they need the whole dataset, so not checkable on a slice)."""
    if name == "lung_hlca":
        ld = hlca_donors()
        h = ld.index[ld.lung_condition.isin(HEALTHY) & ~ld.lung_condition.str.contains(";")]
        i = ld.index[ld.lung_condition.isin(HEALTHY + ("IPF",)) & ld.study_major.isin(IPF_STUDIES)
                     & ~ld.lung_condition.str.contains(";")]
        th = O.healthy_members(pb, h, 20)
        ti = _min_donors(pb.rename(columns={"donor": "resection"}), i, 20).resection.unique()
        return dict(healthy_rule_ok=set(th) == set(donors.index[donors.healthy]),
                    ipf_donors_ok=set(ti) == set(donors.index[donors.lung_ipf_hlca.notna()]))
    return {}
