"""IBDverse (Crohn's disease gut atlas, CELLxGENE collection 0e5c75da; ileum and rectum biopsies): shared builder of
`ibd_ileum` (healthy ileum + Crohn's inflamed / uninflamed ileum) and `ibd_rectum` (healthy rectum).

Legacy: desync/16_ibd/scripts/ibd.py (`X`, `healthy`, `contrast`); A17 targets `ibd_ileum`, `ibd_rectum` (healthy fits)
and `ibd_crohn`, `ibd_crohn_inflamed`, `ibd_crohn_uninflamed` (ileum contrasts, A14). Raw input (group level)
$SCCC_DATA/geo/IBDverse/pseudobulk, built by desync/16_ibd/scripts/build_ibd.py from the 12.7 GB h5ad
$SCCC_DATA/geo/IBDverse/raw (raw/X, author classes).

Choices: groups of the tissue|condition keys of the table (`healthy`, `crohn_inflamed`, `crohn_uninflamed`); classes
`NAMES` (exact); compartment `<tissue>:<class>`; assay "10x 3' v3" is the stratum; pseudobulks >= 20 cells; donor key
`ibdverse:<donor_id>`. Case donors of the Crohn contrasts are donors with inflamed and/or uninflamed ileal samples
(each has one), controls = healthy ileum donors; healthy fits keep classes with >= 20 donors. Age (decade midpoint,
`author_age_binned_10` of the source obs) from `ibd_ileum/ibdverse_age.csv`.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

import census as C
import legacy as L
import organ as O

PATH = C.GEO / "IBDverse" / "pseudobulk"
DS = "ibdverse"
NAMES = ("T", "B", "plasma", "myeloid", "enterocyte", "colonocyte", "secretory", "stem", "stromal", "endothelial")
CLASSES = tuple((c, f"^{c}$") for c in NAMES)
CONDS = {"ibd_ileum": ("ileum", ("healthy", "crohn_inflamed", "crohn_uninflamed")), "ibd_rectum": ("rectum", ("healthy",))}
CROHN = {"ibd_crohn": ("crohn_inflamed", "crohn_uninflamed"), "ibd_crohn_inflamed": ("crohn_inflamed",),
         "ibd_crohn_uninflamed": ("crohn_uninflamed",)}
LEGACY = {"ibd_ileum": ["ibd_ileum", "ibd_crohn", "ibd_crohn_inflamed", "ibd_crohn_uninflamed"], "ibd_rectum": ["ibd_rectum"]}
AGE = Path(__file__).resolve().parents[1] / "ibd_ileum" / "ibdverse_age.csv"


def _age() -> pd.Series:
    a = pd.read_csv(AGE, index_col=0).age_bin
    return (a + 4.5).rename("age").rename(lambda x: f"{DS}:{x}")


def _donors() -> pd.DataFrame:
    d = C.census_donors(PATH, DS)
    d["age"] = _age().reindex(d.index).to_numpy()
    return d


def _controls(name, pb, extra):
    if name == "ibd_ileum":
        return pd.Index(extra["ibd_crohn"].index[extra["ibd_crohn"] == "control"])
    return pd.Index(pb.resection.unique())


def check_legacy(ds, pb, donors):
    """The healthy fit's donors recomputed from the legacy table (needs the whole dataset: not checkable on a slice)."""
    ctrl = donors.index[donors.ibd_crohn == "control"] if "ibd_crohn" in donors else donors.index
    return dict(healthy_rule_ok=set(O.healthy_members(pb, ctrl, 20)) == set(donors.index[donors.healthy]))


def build(name: str, donors=None, genes="clock"):
    tissue, conds = CONDS[name]
    keys = {f"{tissue}|{c}" for c in conds}
    spec = O.Spec(PATH, DS, tissue, CLASSES, keep=lambda o: o.tissue.astype(str).isin(keys).to_numpy())
    o, X = C.read(PATH, genes, donors)
    F = O.groups(spec, read=(o, X))
    pb = O.table(F, F.resection.unique(), 20)
    d = _donors()
    extra = {}
    if name == "ibd_ileum":
        by = o.assign(r=f"{DS}:" + o.donor_id.astype(str), c=o.tissue.astype(str).str.split("|").str[1])
        sets = by.groupby("r").c.agg(set)
        sets = sets[sets.index.isin(pb.resection)]
        healthy = pd.Series(sets.index[sets.map(lambda s: "healthy" in s)])
        for n, cc in CROHN.items():
            case = sets.index[sets.map(lambda s: bool(s & set(cc)))]
            lab = pd.Series("control", index=sets.index[sets.map(lambda s: "healthy" in s)])
            lab = lab.drop(case, errors="ignore")
            lab = pd.concat([lab, pd.Series("case", index=case)])
            extra[n] = lab.rename(n)
    if donors is None:
        extra["healthy"] = O.flag(O.healthy_members(pb, _controls(name, pb, extra), 20))
    return O.finish(pb, d, extra, cols=("age", "disease"))


def from_legacy(name: str):
    d = _donors()
    return L.from_pickles(LEGACY[name], covariates=C.donors_frame(d, d.index, cols=("age", "disease")))
