"""Parkinson's disease brain regions (snRNA-seq of PFC, M1, V1, medial globus pallidus and the dorsal motor nucleus of the
vagus, controls and PD): donor x region x class pseudobulks.

Source: CELLxGENE Census dataset of the PD multi-region atlas (id d3cb449b..., "PD_d3cb449b"), raw input
$SCCC_DATA/atlas/parkinson_regions. Amendment A9. A17 targets: `pd_control` (healthy fit: control donors, compartments
with >= 15 of them; donor column `healthy`) and `pd` (donor column `pd`, case = PD).
Legacy builder replaced: desync/12_parkinson/scripts/parkinson.py (`X`, `labels`, `groups`, `contrast`).

Choices: classes `CLASSES` (case-insensitive, first match), compartment `brain:<region>_<class>` (regions of
`REGIONS`; other tissues dropped); majority assay; pseudobulks >= 20 cells; donors = Census disease "normal" or
"Parkinson disease"; one table for both (the control fit is its control part).
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_atlas"))
import census as C  # noqa: E402
import legacy as L  # noqa: E402
import organ as O  # noqa: E402

DATASET = DS = "parkinson_regions"
REGIONS = {"prefrontal cortex": "PFC", "primary motor cortex": "M1", "primary visual cortex": "V1",
           "medial globus pallidus": "GPi", "dorsal motor nucleus of vagus nerve": "DMV"}
CLASSES = (
    ("glutamatergic", r"glutamatergic"), ("GABAergic", r"GABAergic"), ("neuron_other", r"neuron"),
    ("OPC", r"oligodendrocyte precursor"), ("oligodendrocyte", r"^oligodendrocyte$"), ("astrocyte", r"astrocyte"),
    ("microglia", r"central nervous system macrophage|microglia"), ("vascular", r"endothelial|mural|pericyte"),
    ("ependymal", r"ependymal"),
)
SPEC = O.Spec(C.ATLAS / DS, DS, "brain", CLASSES, region=REGIONS)
REF = "brain:PFC_oligodendrocyte"


def _lab(d):
    return O.labels_from(d, "disease", ["normal"], ["Parkinson disease"])


def build(donors=None, genes="clock"):
    d = C.census_donors(C.ATLAS / DS, DS)
    lab = _lab(d)
    if donors is not None:
        lab = lab[lab.index.isin([f"{DS}:{x}" for x in donors])]
    pb = O.table(O.groups(SPEC, donors, genes), lab.index, 20)
    extra = {"pd": lab.rename("pd")}
    if donors is None:
        extra["healthy"] = O.flag(O.healthy_members(pb, lab.index[lab == "control"], 15))
    return O.finish(pb, d, extra)


def from_legacy():
    d = C.census_donors(C.ATLAS / DS, DS)
    return L.from_pickles(["pd_control", "pd"], covariates=C.donors_frame(d, d.index))


def check_legacy(ds, pb, donors):
    d = C.census_donors(C.ATLAS / DS, DS)
    lab = _lab(d)
    return dict(healthy_rule_ok=set(O.healthy_members(pb, lab.index[lab == "control"], 15)) ==
                set(donors.index[donors.healthy]))
