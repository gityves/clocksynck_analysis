"""Heart atlases (Census): healthy heart (Litviňuková/Kanemaru heart_global, nuclei, left ventricle), dilated
cardiomyopathy (Reichart 2022), myocarditis and the cross-study DCM vs healthy contrast. Amendment A8.
Shared builder of `heart_global`, `heart_dcm`, `heart_myocarditis`, `heart_dcm_healthy`.

Legacy: desync/10_heart/scripts/heart.py (`healthy_groups`, `disease_groups`, `dcm`, `dcm_healthy`, `myocarditis`).
Choices: left-ventricle groups only; classes `CLASSES` (case-insensitive, first match); compartment `heart:<class>`;
majority assay (heart_global: assay | suspension, nuclei only after the majority step); pseudobulks >= 20 cells;
heart_global keeps classes with >= 20 donors (the healthy fit; all donors are `healthy`); disease tables hold the
donors of the contrast (control = Census disease "normal", case = "dilated cardiomyopathy" / "myocarditis").
`heart_dcm_healthy`: DCM cases of Reichart 2022 against the heart_global donors, strata prefixed by study.
"""

from __future__ import annotations

import pandas as pd

import census as C
import legacy as L
import organ as O

CLASSES = (
    ("cardiomyocyte", r"ventricular cardiac myocyte|^cardiac muscle cell$"),
    ("atrial_cardiomyocyte", r"atrial cardiac myocyte"),
    ("fibroblast", r"fibroblast"),
    ("endothelial", r"endothelial|endocardial"),
    ("mural", r"pericyte|smooth muscle|mural"),
    ("myeloid", r"myeloid|macrophage|monocyte|dendritic"),
    ("lymphocyte", r"lymphocyte|(?-i:\bT cell)|natural killer"),
    ("neuron", r"neuron"),
    ("mast", r"mast cell"),
    ("adipocyte", r"adipocyte"),
    ("mesothelial", r"mesothelial"),
)
LV = "heart left ventricle"
DISEASE = {"heart_dcm": ("reichart2022_dcm", "dilated cardiomyopathy", "dcm"),
           "heart_myocarditis": ("heart_myocarditis", "myocarditis", "myocarditis")}
LEGACY = {"heart_global": ["heart"], "heart_dcm": ["dcm"], "heart_myocarditis": ["myocarditis"],
          "heart_dcm_healthy": ["dcm_healthy"]}


def _lv(o):
    return o.tissue.astype(str).to_numpy() == LV


def spec_global() -> O.Spec:
    return O.Spec(C.ATLAS / "heart_global", "heart_global", "heart", CLASSES, keep=_lv, suspension=True,
                  post=lambda F: F[F.batch.str.endswith("|nucleus")])


def spec_disease(ds: str) -> O.Spec:
    return O.Spec(C.ATLAS / ds, ds, "heart", CLASSES, keep=_lv)


def _global_table(donors, genes, min_donors):
    F = O.groups(spec_global(), donors, genes)
    return O.table(F, F.resection.unique(), 20, min_donors=min_donors)


def _disease_table(ds, case, donors, genes):
    d = C.census_donors(C.ATLAS / ds, ds)
    lab = O.labels_from(d, "disease", ["normal"], [case])
    if donors is not None:
        lab = lab[lab.index.isin([f"{ds}:{x}" for x in donors])]
    return d, lab, O.table(O.groups(spec_disease(ds), donors, genes), lab.index, 20)


def build(name: str, donors=None, genes="clock"):
    if name == "heart_global":
        pb = _global_table(donors, genes, None if donors is not None else 20)
        d = C.census_donors(C.ATLAS / "heart_global", "heart_global")
        return O.finish(pb, d, {"healthy": pd.Series(True, index=pb.resection.unique())})
    if name in DISEASE:
        ds, case, tag = DISEASE[name]
        d, lab, pb = _disease_table(ds, case, donors, genes)
        return O.finish(pb, d, {tag: lab.rename(tag)})
    if name == "heart_dcm_healthy":
        ds, case, _ = DISEASE["heart_dcm"]
        d, lab, pbd = _disease_table(ds, case, donors, genes)
        pbd = pbd[pbd.resection.isin(lab.index[lab == "case"])]
        pbd["stratum"] = ds + "|" + pbd.stratum
        pbh = _global_table(donors, genes, None)
        pbh["stratum"] = "heart_global|" + pbh.stratum
        pb = pd.concat([pbd, pbh], ignore_index=True)
        dh = C.census_donors(C.ATLAS / "heart_global", "heart_global")
        dd = pd.concat([d, dh])
        grp = pd.concat([pd.Series("case", index=pbd.resection.unique()), pd.Series("control", index=pbh.resection.unique())])
        return O.finish(pb, dd, {"dcm_healthy": grp.rename("dcm_healthy")})
    raise KeyError(name)


def from_legacy(name: str):
    ds = {"heart_global": "heart_global", "heart_dcm": "reichart2022_dcm", "heart_myocarditis": "heart_myocarditis"}
    if name == "heart_dcm_healthy":
        d = pd.concat([C.census_donors(C.ATLAS / "reichart2022_dcm", "reichart2022_dcm"),
                       C.census_donors(C.ATLAS / "heart_global", "heart_global")])
    else:
        d = C.census_donors(C.ATLAS / ds[name], ds[name])
    pb, donors = L.from_pickles(LEGACY[name], covariates=C.donors_frame(d, d.index))
    return pb, donors
