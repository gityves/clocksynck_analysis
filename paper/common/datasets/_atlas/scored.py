"""Cohort tables with the A51 score genes (Hallmark TNFA_SIGNALING_VIA_NFKB and INTERFERON_GAMMA_RESPONSE, minus the clock,
immediate-early, hypoxia and housekeeping genes of the panel; 359 genes, `inputs/score_all.txt`) added after the 14 clock
genes: the `<cohort>_scored` datasets, inputs of the immune-damping test (A51, `desync/validation/immune_damping`).

Rows are those of the A17 table of the same cohort (checked by the legacy `rebuild_check.csv`: resection, compartment,
stratum, n_cells, depth and the clock columns equal). `build` runs the builder of the base dataset on the rebuilt group
pseudobulks of `$SCCC_DATA/immune_damping/{atlas,geo}` with `genes="scored"` (the group files were rebuilt from the same
Census datasets / GEO matrices with the extended gene panel by `desync/validation/immune_damping/scripts/build_census.sh`,
`build_local.sh`); `from_legacy` reads the tables `$SCCC_DATA/immune_damping/tables/<cohort>.pkl` of the legacy run.
Donor columns: those of the base dataset.
"""

from __future__ import annotations

import importlib.util
import pickle
from pathlib import Path

import census as C

BASE = {  # dataset: (base dataset, immune_damping cohort pickle)
    "lung_natri_scored": ("lung_natri", "ild_natri"), "hlca_fibrotic_scored": ("hlca_fibrotic", "hlca"),
    "heart_dcm_scored": ("heart_dcm", "dcm"), "heart_myocarditis_scored": ("heart_myocarditis", "myocarditis"),
    "kidney_sc_scored": ("kidney_sc", "ckd_sc"), "ia_crohn_scored": ("ia_crohn", "ia_crohn"),
    "skin_ssc_scored": ("skin_ssc", "skin_ssc"), "skin_ad_scored": ("skin_ad", "skin_ad"),
    "endometrium_scored": ("endometrium", "endometriosis"),
}
TABLES = C.DATA / "immune_damping" / "tables"


def base_module(name: str):
    f = C.HERE.parent / name / "build.py"
    spec = importlib.util.spec_from_file_location(f"build_{name}", f)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def build(name: str, donors=None, genes="scored"):
    base, _ = BASE[name]
    pb, d = base_module(base).build(donors=donors, genes=genes)
    return pb, d


def from_legacy(name: str):
    base, cohort = BASE[name]
    with open(TABLES / f"{cohort}.pkl", "rb") as f:
        t = pickle.load(f)["table"]
    t["stratum"] = t["stratum"].astype(str)
    pb = C.standard(t)
    d = base_module(base).from_legacy()[1]
    return pb, d.reindex(pb.donor.unique())
