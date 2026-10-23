"""Held-out cohorts of the immune-damping test (A51) and the gut / liver inflammation studies whose tables they reuse:
shared builder of `kong_cd`, `periodontitis`, `taurus`, `masld_pooled` and their `_scored` twins (clock genes followed by the 359
A51 score genes).

Legacy: desync/validation/immune_damping/scripts/ext_heldout.py (`census_cohort`, `taurus`, `masld`), with the loaders of the
original analyses: desync/validation/ibd_taurus/scripts/a45_taurus.py (`samples`, `table`; A45) and
desync/validation/masld/scripts/a44_pooled.py (`patients`, `table`; A43/A44 pool; A38 for GSE202379).
Raw inputs (group pseudobulks): kong_cd = $SCCC_DATA/atlas/{ileum_kong,kong_ti_str} (Kong 2023 immune and stromal fractions, Census
a37f857c / 0f4865d5); periodontitis = $SCCC_DATA/immune_damping/atlas/periodontitis (Census 9fcb0b73; only the extended build exists,
its clock columns are used for the clock-only table); taurus = $SCCC_DATA/ibd/TAURUS/pseudobulk (Thomas et al. 2024, Zenodo 13768607,
sample x class, from the 12.7 GB h5ad by desync/validation/ibd_taurus/scripts/build_taurus.py); masld_pooled =
$SCCC_DATA/masld/{GSE202379,GSE289173,GSE244832,GSE256398}/pseudobulk (patient x class, built by desync/validation/masld/scripts/
build_*.py with the annotations of each study).

Choices.
- kong_cd: ileal epithelium and lamina propria samples, normal and Crohn disease donors, classes `COARSE`, compartment `ileum:<class>`,
  majority assay per (file, disease) call, key `kong:<donor>`, pseudobulks >= 20 cells; control = normal, case = Crohn.
- periodontitis: normal, gingivitis and periodontitis donors (control = normal, case = the others), classes `COARSE`, compartment
  `mucosa:<class>`, key `perio:<donor>`, >= 20 cells.
- taurus: unit = biopsy sample (key `taurus:<sample>`, its own host state), classes `TAURUS` (exact), majority assay, >= 20 cells,
  stratum = gut site; permutation block = patient (donor column `patient`; A45 also uses `patient_site`), arm = inflamed (case) /
  non-inflamed (control) among treated patients, healthy subjects are `control` in the `taurus` label of the A51 table.
- masld_pooled: unit = patient (GSE202379: one sample per patient, the right lobe when there are two), classes `MASLD` plus plasma,
  >= 20 cells, stratum = study; D1 / D2 arms of A44 in the donor table (`arm_d1`, `arm_d2`); `living_donor` marks the GSE289173 SNUH
  controls. The A51 table kept the 126 patients with a D2 label; the rebuilt table keeps all 130.
Anchors used by the legacy runs are in `ANCHORS` (A17 fixed rule or learned rule per cohort; A38 for MASLD).
"""

from __future__ import annotations

import numpy as np
import pandas as pd

import census as C
import organ as O

BASES = ("kong_cd", "periodontitis", "taurus", "masld_pooled", "masld_pooled_a96")
ANCHORS = {"kong_cd": ["ileum:endothelial", "ileum:fibroblast_stromal", "ileum:mural_smooth_muscle"],
           "periodontitis": ["mucosa:endothelial", "mucosa:fibroblast_stromal", "mucosa:mural_smooth_muscle"],
           "taurus": ["gut:endothelial", "gut:fibroblast", "gut:mural"],
           "masld_pooled": ["liver:endothelial", "liver:stellate", "liver:macrophage"]}
REF = {"kong_cd": "ileum:CD8_T", "periodontitis": "mucosa:endothelial", "taurus": "gut:endothelial",
       "masld_pooled": "liver:endothelial"}  # reference compartments of the legacy fits
IMMUNE = C.DATA / "immune_damping"
TABLES = IMMUNE / "tables"
LEGACY_COHORT = {"kong_cd": "kong_cd", "periodontitis": "periodontitis", "taurus": "taurus", "masld_pooled": "masld"}
# ---- kong_cd, periodontitis ----------------------------------------------------------------------------------------------
KONG = (("ileum_kong", "kong_ti_str"), "ileum", "kong", ("ileal epithelium", "ileum lamina propria"),
        ("normal", "Crohn disease"), "normal")
PERIO = (("periodontitis",), "mucosa", "perio", None, ("normal", "gingivitis", "periodontitis"), "normal")


def _path(base: str, d: str):
    return (C.ATLAS / d) if base == "kong_cd" else IMMUNE / "atlas" / d


def _census(base: str, donors, genes):
    dirs, organ, ds, tissues, diseases, control = KONG if base == "kong_cd" else PERIO
    parts, dis = [], {}
    donor_tabs = []
    for d in dirs:
        donor_tabs.append(C.census_donors(_path(base, d), ds))
        for dd in diseases:
            keep = lambda o, dd=dd: ((o.disease.astype(str) == dd).to_numpy() if tissues is None else  # noqa: E731
                                     (o.tissue.astype(str).isin(tissues) & (o.disease.astype(str) == dd)).to_numpy())
            g = O.groups(O.Spec(_path(base, d), ds, organ, O.COARSE, keep=keep), donors, genes)
            parts.append(g)
            dis.update({r: dd for r in g.resection.unique()})
    F = pd.concat(parts, ignore_index=True)
    d = pd.concat(donor_tabs).loc[lambda t: ~t.index.duplicated()]
    pb = O.table(F, F.resection.unique(), 20)
    lab = pd.Series(dis).map(lambda v: "control" if v == control else "case")
    return pb, d, lab[lab.index.isin(pb.resection.unique())], pd.Series(dis)


# ---- taurus ---------------------------------------------------------------------------------------------------------------
TAURUS = C.DATA / "ibd" / "TAURUS" / "pseudobulk"
TAURUS_CLASSES = tuple((c, f"^{c}$") for c in ("fibroblast", "endothelial", "mural", "T", "B_plasma", "myeloid", "epithelial"))


def _taurus_samples() -> pd.DataFrame:
    s = pd.read_csv(TAURUS / "samples.csv")
    s["resection"] = "taurus:" + s["sample"].astype(str)
    s["arm"] = s.inflammation.map({"Inflamed": "case", "Non_Inflamed": "control"})
    return s.set_index("resection")


def _taurus(donors, genes):
    sm = _taurus_samples()
    o, X = C.read(TAURUS, genes, donors)  # genes="scored": the rebuild in immune_damping/ibd/TAURUS
    F = O.groups(O.Spec(TAURUS, "taurus", "gut", TAURUS_CLASSES), read=(o, X))
    pb = C.donor_table(F, sm.index, 20, stratum=None)
    pb["stratum"] = pb.resection.map(sm.site).astype(str)
    return pb, sm


# ---- masld_pooled ----------------------------------------------------------------------------------------------------------
MASLD_ROOT = C.DATA / "masld"
MASLD = tuple((c, f"^{c}$") for c in ("hepatocyte", "cholangiocyte", "endothelial", "stellate", "macrophage", "lymphocyte",
                                       "B", "neutrophil", "plasma"))
F244 = {"D11": "F4", "D12": "F4", "D15": "F4", "D18": "F4", "D14": "F3", "D16": "F3", "D10": "F2-3", "D13": "F2-3", "D17": "F2-3"}
STUDIES = (("GSE202379", "masld"), ("GSE289173", "gse289173"), ("GSE244832", "gse244832"), ("GSE256398", "gse256398"))


def _masld_patients(root=MASLD_ROOT) -> pd.DataFrame:
    rows = []
    s = pd.read_csv(root / "GSE202379" / "pseudobulk" / "samples.csv")  # A38: one sample per patient, right lobe first
    s["right"] = s.Lobe.eq("Right")
    s = s.sort_values(["right", "orig.ident"], ascending=[False, True]).drop_duplicates("Patient.ID")
    s["resection"] = "masld:P" + s["Patient.ID"].astype(str)
    m = s["SAF.Score"].str.extract(r"S(\d)A(\d)F(\d)").astype(float)
    s["F"] = m[2]
    p = s.set_index("resection")
    st = p["Disease.status"]
    d1 = np.where(st.isin(["Healthy control", "NAFLD"]), "control", "case")
    F = np.where(st.eq("end stage"), 4, np.where(p.F.isna() & st.eq("Healthy control"), 0, p.F))
    d2 = np.where(pd.isna(F), None, np.where(F.astype(float) >= 3, "case", "control"))
    rows.append(pd.DataFrame(dict(study="GSE202379", sample=p["orig.ident"].to_numpy(), diagnosis=st.to_numpy(),
                                  fibrosis=pd.Series(F).astype(str).to_numpy(), arm_d1=d1, arm_d2=d2, hospital="",
                                  age=np.nan), index=p.index))
    s = pd.read_csv(root / "GSE289173" / "pseudobulk" / "samples.csv")
    d1 = np.where(s.group.isin(["eMASH", "aMASH"]), "case", "control")
    d2 = np.where(s.group.eq("aMASH"), "case", "control")
    rows.append(pd.DataFrame(dict(study="GSE289173", sample=s["orig.ident"].to_numpy(), diagnosis=s.group.to_numpy(),
                                  fibrosis=np.where(s.group.eq("aMASH"), "F3", "F0-1"), arm_d1=d1, arm_d2=d2,
                                  hospital=s.hospital.to_numpy(), age=np.nan),
                             index="gse289173:" + s["orig.ident"].astype(str)))
    s = pd.read_csv(root / "GSE244832" / "pseudobulk" / "samples.csv")
    fib = s.donor.map(F244).fillna("F0-2 (not MASH)")
    d1 = np.where(s.condition.eq("MASH"), "case", "control")
    d2 = np.where(fib.isin(["F3", "F4"]), "case", np.where(fib.eq("F2-3"), None, "control"))
    rows.append(pd.DataFrame(dict(study="GSE244832", sample=s["orig.ident"].to_numpy(),
                                  diagnosis=(s.condition + " " + s.donor).to_numpy(), fibrosis=fib.to_numpy(), arm_d1=d1,
                                  arm_d2=d2, hospital="", age=np.nan), index="gse244832:" + s["orig.ident"].astype(str)))
    s = pd.read_csv(root / "GSE256398" / "pseudobulk" / "samples.csv")
    dg = s.diagnosis
    d1 = np.where(dg.str.startswith("MASH"), "case", "control")
    d2 = np.where(dg.eq("MASH Fibrosis F2-3"), None, np.where(dg.isin(["MASH Cirrhosis", "MASH Fibrosis F3"]), "case", "control"))
    rows.append(pd.DataFrame(dict(study="GSE256398", sample=s["orig.ident"].to_numpy(), diagnosis=dg.to_numpy(),
                                  fibrosis=dg.str.extract(r"(F[\d-]+|Cirrhosis)")[0].fillna("F0").to_numpy(), arm_d1=d1,
                                  arm_d2=d2, hospital="", age=s.age.str.extract(r"(\d+)")[0].astype(float).to_numpy()),
                             index="gse256398:" + s["orig.ident"].astype(str)))
    pt = pd.concat(rows)
    pt.index.name = "resection"
    pt["living_donor"] = pt.study.eq("GSE289173") & pt.hospital.eq("SNUH")
    return pt


def _masld(donors, genes, root=MASLD_ROOT):
    pt = _masld_patients(root)
    parts = []
    for study, ds in STUDIES:
        path = root / study / "pseudobulk"
        q = pt[pt.study == study]
        smp = pd.Series(q.index, index=f"{ds}:" + q["sample"].astype(str))
        sl = None
        if donors is not None:  # slice: the samples of the requested patients
            sl = [s for s, r in zip(q["sample"], q.index) if r.split(":", 1)[1] in {str(x) for x in donors}]
            if not sl:
                continue
        o, X = C.read(path, genes, sl)
        F = O.groups(O.Spec(path, ds, "liver", MASLD), read=(o, X))
        F = F[F.resection.isin(smp.index)].copy()
        F["resection"] = F.resection.map(smp).to_numpy()
        tb = C.donor_table(F, q.index, 20, stratum=None)
        tb["stratum"] = study
        parts.append(tb)
    return pd.concat(parts, ignore_index=True), pt


# ---- interface ------------------------------------------------------------------------------------------------------------
def _split(name: str):
    base = name.removesuffix("_scored")
    return base, ("scored" if name.endswith("_scored") else "clock")


def build(name: str, donors=None, genes=None):
    """`donors`: slice (ids after the dataset prefix: Census donor ids; TAURUS samples; MASLD: `<study>:...` keys' second part)."""
    base, g = _split(name)
    genes = g if genes is None else genes
    if base in ("kong_cd", "periodontitis"):
        pb, d, lab, dis = _census(base, donors, genes)
        return O.finish(pb, d, {base: lab.rename(base), "disease_label": dis.rename("disease_label")}, cols=("age", "sex", "disease"))
    if base == "taurus":
        pb, sm = _taurus(donors, genes)
        pb = C.standard(pb)
        cols = [c for c in sm.columns if c not in ("sample",)]
        lab = pd.Series(np.where(sm.inflammation == "Healthy", "control", "case"), index=sm.index)
        d = sm.loc[sm.index.isin(pb.donor.unique()), cols].copy()
        d["taurus"] = lab.reindex(d.index)
        d.index.name = "donor"
        return pb, d.reindex(pb.donor.unique())
    if base == "masld_pooled_a96":  # A96: panel of core + cue + locked genes, patients with a D1 arm
        pb, pt = _masld(donors, "all" if genes in (None, "clock", "scored") else genes, MASLD_ROOT / "a96")
        pb = C.standard(pb.rename(columns=lambda c: _canon(c) if c not in ("resection", "compartment", "stratum", "n_cells",
                                                                         "depth") else c))
        d = pt.reindex(pb.donor.unique()).rename_axis("donor")
        d["arm"] = d.arm_d1
        return pb, d
    pb, pt = _masld(donors, genes)
    pb = C.standard(pb)
    d = pt.reindex(pb.donor.unique()).rename_axis("donor")
    return pb, d


from common.paths import ARCHIVE  # noqa: E402
def _canon(g):
    from clocksync.genes import canonical
    return canonical(g)


def _a96_legacy(dataset: str):
    f = ARCHIVE / "legacy" / "desync" / "validation" / "masld" / "output" / "a96_counts.csv"
    t = pd.read_csv(f)
    t = t[t.dataset == dataset].dropna(axis=1, how="all").rename(columns={"resection": "donor"})
    arm = t.groupby("donor").arm.first()
    t = t.drop(columns=["dataset", "arm"])
    return C.standard(t.rename(columns={"donor": "resection"})), arm


def from_legacy(name: str):
    import pickle

    if name == "masld_pooled_a96":
        pb, arm = _a96_legacy("liver")
        pt = _masld_patients(MASLD_ROOT / "a96").reindex(pb.donor.unique()).rename_axis("donor")
        pt["arm"] = arm.reindex(pt.index)
        return pb, pt
    base, g = _split(name)
    with open(TABLES / f"{LEGACY_COHORT[base]}.pkl", "rb") as f:
        rec = pickle.load(f)
    t = rec["table"].copy()
    t["stratum"] = t["stratum"].astype(str)
    if g == "clock":
        t = t[["resection", "compartment", "stratum", "n_cells", "depth"] + list(C.GENES)]
    pb = C.standard(t)
    lab = rec["labels"]
    d = pd.DataFrame(index=pd.Index(pb.donor.unique(), name="donor"))
    if base in ("kong_cd", "periodontitis"):
        d[base] = lab.reindex(d.index)
        d["disease_label"] = rec["disease"].reindex(d.index)
        cov = C.census_donors(_path(base, rec_dirs(base)[0]), KONG[2] if base == "kong_cd" else PERIO[2])
        d = d.join(C.donors_frame(cov, cov.index, cols=("age", "sex", "disease")))
    elif base == "taurus":
        sm = rec["samples"]
        d = sm.reindex(d.index).drop(columns=["sample"], errors="ignore")
        d["taurus"] = lab.reindex(d.index)
    else:
        d = rec["patients"].reindex(d.index)
    return pb, d


def rec_dirs(base):
    return KONG[0] if base == "kong_cd" else PERIO[0]
