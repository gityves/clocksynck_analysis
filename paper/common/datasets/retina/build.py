"""Human retina snRNA-seq atlas (Li et al. 2023, 95 individual donors): donor x retinal class pseudobulks.

Source: CELLxGENE Census dataset d6505c89-c43d-4c28-8c4f-7351a5fd5528 (Li, Jiang, ... Chen; Cell 2023), raw input
`$SCCC_DATA/atlas/retina_li2023` (group pseudobulk of the Census build). Amendment A2.
Legacy builder replaced: desync/common/run_atlas.py::load (classes `RETINA`) + desync/04_retina/scripts/retina_checks.py
::table (A17 target `retina`).

Choices (as the legacy run): 10x groups; classes by the first matching `RETINA` pattern on the Census cell_type
(case-sensitive); compartment `eye:<class>`; majority assay per donor; pseudobulks >= 50 nuclei; stratum = the donor's
majority study (`donor_strata.csv`, from `desync/04_retina/scripts/retina_strata.py`, which reads the source h5ad
obs over S3; kept here as a small input file); the 9 pooled-sample donor ids are not in that file and are dropped;
classes with >= 20 donors.
"""

import sys
from pathlib import Path

import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[0] / "_atlas"))
import census as C  # noqa: E402
import legacy as L  # noqa: E402

DATASET, PREFIX, MIN_CELLS, MIN_DONORS = "retina", "retina_li2023", 50, 20
RETINA = (
    ("rod_bipolar", r"^rod bipolar"), ("rod", r"retinal rod cell"), ("cone", r"cone cell$"),
    ("cone_bipolar", r"bipolar"), ("horizontal", r"horizontal"), ("amacrine", r"amacrine"),
    ("ganglion", r"ganglion|OFFx"), ("Mueller", r"Mueller"), ("astrocyte", r"^astrocyte$"),
    ("microglia", r"microglial"),
)
REF = "eye:amacrine"  # reference compartment of the generation-1 fit (most donors)


def strata() -> pd.Series:
    s = pd.read_csv(HERE / "donor_strata.csv")
    return s.set_index(s.dataset + ":" + s.donor_id.astype(str)).stratum


def build(donors=None, genes="clock"):
    o, X = C.read(C.ATLAS / "retina_li2023", genes, donors)
    cell = C.classify(o.cell_type.astype(str), RETINA, ignore_case=False)
    umi = o.assay.astype(str).str.startswith("10x").to_numpy()
    F = C.keep_majority(C.frame(o, X, PREFIX, o.tissue_general.astype(str) + ":" + cell.astype(str), o.assay,
                                keep=umi & cell.notna().to_numpy()))
    pb = C.donor_table(F, None, MIN_CELLS, stratum="batch")
    pb["stratum"] = pb.resection.map(strata())
    pb = pb[pb.stratum.notna()]
    if donors is None:
        n = pb.groupby("compartment").resection.nunique()
        pb = pb[pb.compartment.isin(n[n >= MIN_DONORS].index)]
    pb = C.standard(pb)
    d = C.census_donors(C.ATLAS / "retina_li2023", PREFIX)
    out = C.donors_frame(d, pb.donor.unique(), {"study": strata().str.split("|").str[1]})
    return pb, out


def from_legacy():
    d = C.census_donors(C.ATLAS / "retina_li2023", PREFIX)
    pb, donors = L.from_pickles(["retina"], covariates=C.donors_frame(d, d.index, {"study": strata().str.split("|").str[1]}))
    return pb, donors
