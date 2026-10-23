"""Cortex atlases (PsychAD x3 cohorts, SEA-AD x2 regions): shared builder of the five `<dataset>/build.py`.

Legacy: `desync/common/run_atlas.py::load` (classes `BRAIN`), `desync/03_cortex_conditions/scripts/brain_conditions.py`
(contrast definitions), `desync/common/a17_tables.py::cortex_healthy|cortex_contrast`.

Choices (as the legacy A1 / A13 runs): 10x groups only; classes of the Census cell_type by the first matching
`BRAIN` pattern (case-sensitive, unmatched dropped); compartment `brain:<class>`; each donor keeps its majority
assay; stratum `<dataset>|<assay>`; pseudobulks of at least 50 nuclei; every class kept (the contrasts apply their own
donor rule, the healthy fit kept classes with >= 20 healthy donors: the `healthy` column marks the donors of that
fit). Contrast labels follow `brain_conditions.CONTRASTS` (case/control per donor, with its age restriction;
tertiles of age among normal donors >= 18 y for the age contrasts).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

import census as C
import legacy as L

BRAIN = (
    ("excitatory", r"glutamatergic"), ("inhibitory_MGE", r"pvalb|sst|chandelier"),
    ("inhibitory_CGE", r"VIP|lamp5|sncg|caudal ganglionic"), ("inhibitory_other", r"GABAergic"),
    ("astrocyte", r"astrocyte"), ("oligodendrocyte", r"^oligodendrocyte$"), ("OPC", r"oligodendrocyte precursor"),
    ("microglia", r"microglial"), ("vascular", r"endothelial|pericyte|smooth muscle|leptomeningeal"),
)
MIN_CELLS, MIN_DONORS = 50, 20


def _is(col, *values):
    return lambda d: d[col].isin(values)


# contrast: (dataset, case, control, restriction); case/control are functions of the donor table or the age tertiles
CONTRASTS = {
    "AD_mssm": ("psychad_mssm", _is("AD_status", "Yes"), _is("disease", "normal"), lambda d: d.age_bc >= 65),
    "AD_radc": ("psychad_radc", _is("AD_status", "Yes"), _is("disease", "normal"), lambda d: d.age_bc >= 65),
    "AD_seaad_dlpfc": ("seaad_dlpfc", _is("Cognitive status", "Dementia"), _is("Cognitive status", "No dementia"), None),
    "AD_seaad_mtg": ("seaad_mtg", _is("Cognitive status", "Dementia"), _is("Cognitive status", "No dementia"), None),
    "ADNC_seaad_dlpfc": ("seaad_dlpfc", _is("ADNC", "High"), _is("ADNC", "Not AD", "Low", "Intermediate"), None),
    "ADNC_seaad_mtg": ("seaad_mtg", _is("ADNC", "High"), _is("ADNC", "Not AD", "Low", "Intermediate"), None),
    "SCZ_hbcc": ("psychad_hbcc", _is("Schizophrenia", "Yes"), _is("disease", "normal"), lambda d: d.age_bc >= 18),
    "SCZ_mssm": ("psychad_mssm", lambda d: (d.Schizophrenia == "Yes") & (d.AD_status == "No"),
                 _is("disease", "normal"), lambda d: d.age_bc >= 18),
    "BD_hbcc": ("psychad_hbcc", _is("Bipolar_Disorder", "Yes"), _is("disease", "normal"), lambda d: d.age_bc >= 18),
    "age_mssm": ("psychad_mssm", "top_tertile", "bottom_tertile",
                 lambda d: (d.disease == "normal") & (d.age_bc >= 18)),
    "age_hbcc": ("psychad_hbcc", "top_tertile", "bottom_tertile",
                 lambda d: (d.disease == "normal") & (d.age_bc >= 18)),
    "DLBD_mssm": ("psychad_mssm", _is("DLBD_status", "Yes"), _is("disease", "normal"), lambda d: d.age_bc >= 65),
    "DLBpure_mssm": ("psychad_mssm", lambda d: (d.DLBD_status == "Yes") & (d.AD_status == "No"),
                     _is("disease", "normal"), lambda d: d.age_bc >= 65),
    "PD_mssm": ("psychad_mssm", _is("Parkinson_disease", "Yes"), _is("disease", "normal"), lambda d: d.age_bc >= 65),
    "VaD_mssm": ("psychad_mssm", _is("Vascular_status", "Yes"), _is("disease", "normal"), lambda d: d.age_bc >= 65),
    "ASCVD_hbcc": ("psychad_hbcc", _is("ASCVD_status", "Yes"), _is("disease", "normal"), lambda d: d.age_bc >= 18),
    "DM_radc": ("psychad_radc", lambda d: d.disease.astype(str).str.contains("diabetes"), _is("disease", "normal"),
                lambda d: d.age_bc >= 18),
}
LEGACY = {  # dataset -> pickles of the A17 re-analysis that read its table
    "psychad_mssm": ["psychad_mssm", "AD_mssm", "SCZ_mssm", "DLBD_mssm", "DLBpure_mssm", "PD_mssm", "VaD_mssm",
                     "age_mssm"],
    "psychad_hbcc": ["psychad_hbcc", "SCZ_hbcc", "BD_hbcc", "age_hbcc", "ASCVD_hbcc"],
    "psychad_radc": ["AD_radc", "DM_radc"],
    "seaad_dlpfc": ["seaad_dlpfc", "AD_seaad_dlpfc", "ADNC_seaad_dlpfc"],
    "seaad_mtg": ["seaad_mtg", "AD_seaad_mtg", "ADNC_seaad_mtg"],
}
HEALTHY_FIT = {"psychad_mssm", "psychad_hbcc", "seaad_dlpfc", "seaad_mtg"}  # radc had no healthy fit


def _donors(ds: str) -> pd.DataFrame:
    d = C.census_donors(C.ATLAS / ds, ds)
    s = d.development_stage.astype(str)  # age of brain_conditions.donors (80+ coded 90), used for the restrictions
    d["age_bc"] = s.str.extract(r"(\d+)")[0].astype(float).where(~s.str.contains("80 year-old and over"), 90.0).to_numpy()
    return d


def labels(contrast: str, d: pd.DataFrame) -> pd.Series:
    """'case' / 'control' per donor (index of `d`) of one contrast; donors outside the contrast are absent."""
    _, case, control, keep = CONTRASTS[contrast]
    if keep is not None:
        d = d[keep(d).to_numpy()]
    if case == "top_tertile":
        lo, hi = d.age_bc.quantile([1 / 3, 2 / 3])
        lab = np.select([d.age_bc >= hi, d.age_bc <= lo], ["case", "control"], "")
    else:
        lab = np.select([case(d).to_numpy(), control(d).to_numpy()], ["case", "control"], "")
    lab = pd.Series(lab, index=d.index)
    return lab[lab != ""]


def healthy_donors(pb: pd.DataFrame, d: pd.DataFrame) -> pd.Index:
    """Donors of the healthy fit: normal donors, with the classes that have >= 20 of them (`pb`: resection column)."""
    col = "resection" if "resection" in pb else "donor"
    h = pb[pb[col].isin(d.index[d.disease.astype(str).eq("normal")])]
    n = h.groupby("compartment")[col].nunique()
    return pd.Index(h[h.compartment.isin(n[n >= MIN_DONORS].index)][col].unique())


def check_legacy(ds: str, pb: pd.DataFrame, donors: pd.DataFrame) -> dict:
    """Recompute the `healthy` column from the legacy table with the rule of `build` (a rule that needs the whole
    dataset and so cannot be checked on a slice)."""
    if "healthy" not in donors:
        return {}
    h = healthy_donors(pb, _donors(ds))
    return dict(healthy_rule_ok=bool(set(h) == set(donors.index[donors.healthy])))


def build(ds: str, donors=None, genes="clock"):
    """(pb, donors) of one cortex dataset from `<SCCC_DATA>/atlas/<ds>/pseudobulk.h5ad`."""
    o, X = C.read(C.ATLAS / ds, genes, donors)
    cell = C.classify(o.cell_type.astype(str), BRAIN, ignore_case=False)
    umi = o.assay.astype(str).str.startswith("10x").to_numpy()
    F = C.frame(o, X, ds, o.tissue_general.astype(str) + ":" + cell.astype(str), o.assay, keep=umi & cell.notna().to_numpy())
    F = C.keep_majority(F)
    pb = C.donor_table(F, None, MIN_CELLS, stratum="batch")
    pb["stratum"] = ds + "|" + pb.stratum
    d = _donors(ds)
    extra = {c: labels(c, d).rename(c) for c, v in CONTRASTS.items() if v[0] == ds}
    if donors is None and ds in HEALTHY_FIT:
        extra["healthy"] = pd.Series(True, index=healthy_donors(pb, d))
    pb = C.standard(pb)
    out = C.donors_frame(d, pb.donor.unique(), extra)
    if "healthy" in out:
        out["healthy"] = out.healthy.fillna(False).astype(bool)
    return pb, out


def from_legacy(ds: str):
    d = _donors(ds)
    pb, donors = L.from_pickles(LEGACY[ds], covariates=C.donors_frame(d, d.index))
    return pb, donors
