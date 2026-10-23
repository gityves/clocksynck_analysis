"""Immune Aging atlas (Census immune_aging_a7: Dominguez Conde et al. / Immune Aging 1b350d0a, 12 donors' tissues) and the
Gut Cell Atlas ileum (Crohn's vs normal): shared builder of `ia_spleen`, `ia_thoracic_ln`, `ia_mesenteric_ln`,
`ia_jejunum`, `ia_lung` (healthy fits, A15) and `ia_crohn` (A15 part 4).

Legacy: desync/17_immune_aging/scripts/immune_aging.py (`fit`, `crohn`), A17 targets of the same names.
Raw inputs $SCCC_DATA/atlas/immune_aging_a7 and $SCCC_DATA/atlas/gut_atlas_a7 (group pseudobulks).

Choices: normal groups of the tissue(s) (jejunum: lamina propria and epithelium); classes `COARSE`; compartment
`ia:<class>` (`ileum:<class>` for the gut atlas); majority assay; pseudobulks >= 20 cells; classes with >= `MIN_DONORS`
donors (20; jejunum and lung 15). `ia_crohn`: ileum groups of normal and "Crohn disease" donors (majority assay per
disease set), controls = normal donors, cases = Crohn donors, all classes kept.
"""

from __future__ import annotations

import pandas as pd

import census as C
import legacy as L
import organ as O

IA, GUT = "immune_aging_a7", "gut_atlas_a7"
TISSUES = {"ia_spleen": ("spleen",), "ia_thoracic_ln": ("thoracic lymph node",),
           "ia_mesenteric_ln": ("mesenteric lymph node",), "ia_jejunum": ("jejunum lamina propria", "jejunal epithelium"),
           "ia_lung": ("lung",)}
MIN_DONORS = {"ia_spleen": 20, "ia_thoracic_ln": 20, "ia_mesenteric_ln": 20, "ia_jejunum": 15, "ia_lung": 15}


def _spec(tissues) -> O.Spec:
    keep = lambda o: o.tissue.astype(str).isin(tissues).to_numpy() & (o.disease.astype(str) == "normal").to_numpy()  # noqa: E731,E501
    return O.Spec(C.ATLAS / IA, IA, "ia", O.COARSE, keep=keep)


def _gut_spec(dis: str) -> O.Spec:
    keep = lambda o: (o.tissue.astype(str) == "ileum").to_numpy() & (o.disease.astype(str) == dis).to_numpy()  # noqa: E731,E501
    return O.Spec(C.ATLAS / GUT, GUT, "ileum", O.COARSE, keep=keep)


def build(name: str, donors=None, genes="clock"):
    if name == "ia_crohn":
        sn, sc = _gut_spec("normal"), _gut_spec("Crohn disease")
        Fn, Fc = O.groups(sn, donors, genes), O.groups(sc, donors, genes)
        pb = O.table(pd.concat([Fn, Fc], ignore_index=True), None, 20)
        lab = pd.concat([pd.Series("control", index=Fn.resection.unique()).drop(Fc.resection.unique(), errors="ignore"),
                         pd.Series("case", index=Fc.resection.unique())])
        return O.finish(pb, C.census_donors(sn.path, GUT), {"ia_crohn": lab.rename("ia_crohn")})
    s = _spec(TISSUES[name])
    F = O.groups(s, donors, genes)
    pb = O.table(F, F.resection.unique(), 20, min_donors=None if donors is not None else MIN_DONORS[name])
    return O.finish(pb, C.census_donors(s.path, IA), {"healthy": O.flag(pb.resection.unique())})


def from_legacy(name: str):
    d = C.census_donors(C.ATLAS / (GUT if name == "ia_crohn" else IA), GUT if name == "ia_crohn" else IA)
    return L.from_pickles([name], covariates=C.donors_frame(d, d.index))
