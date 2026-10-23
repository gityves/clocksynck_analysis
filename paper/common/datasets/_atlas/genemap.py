"""Pooled normal adult lung and DLPFC atlases of the gene-rhythm map (A50): shared builder of `gene_map_lung` and `gene_map_dlpfc`.

Legacy: desync/validation/gene_map/scripts/pool.py (`prep`, `read_groups`, `majority_assay`, `aggregate`; output `agg_<atlas>.pkl`)
with the A31b lung rules of desync/validation/gonogo_new/scripts/atlases.py::_lung_groups.
Raw inputs: `$SCCC_DATA/gene_map/<atlas>/<dataset id[:8]>/pseudobulk.h5ad` (normal, adult, primary, UMI cells of each selected Census
dataset, genome-wide, built by desync/common/build_pseudobulk.py from the queries of `inputs/gene_map_selection.csv`, i.e.
desync/validation/gene_map/scripts/select_datasets.py), and for the lung the HLCA healthy donors of `atlas_wg/lung_hlca` (as
`lung_hlca_genome`) with `inputs/lung_donors.csv`.

Choices: donor key `<collection id>:<donor id>`; lung classes: exact HLCA labels (Census cell types whose HLCA cells all share one A4 class)
plus the synonyms of `LUNG_SYN`; DLPFC classes: `BRAIN` (first matching pattern, case-insensitive); compartment `lung:<class>` /
`brain:<class>`; majority assay per donor over all datasets; groups summed per donor x compartment; pseudobulks >= 20 cells; stratum =
assay | block (block = collection; PsychAD's three sequencing cohorts are three blocks; HLCA donors: majority study). No minimum of donors
per class: the analysis sets of A50 (`others`, `pooled`, `*_nondisease` for lung; `pooled`, `psychad`, `other` for DLPFC) are cut from the
donor columns `source`, `group`, `disease_study` and apply their own rule (>= 20 donors per class). Genes: every gene (the legacy tables
were cut to the GTEx-significant genes: `from_legacy` has that panel only).
"""

from __future__ import annotations

import pickle
import re

import numpy as np
import pandas as pd

import census as C
import cortex as CX
import lung as LU
import organ as O

SELECTION = C.HERE / "inputs" / "gene_map_selection.csv"
ROOT = C.DATA / "gene_map"
from common.paths import ARCHIVE  # noqa: E402
LEGACY_OUT = ARCHIVE / "legacy" / "desync" / "validation" / "gene_map" / "output"
LUNG_SYN = (
    ("EC_capillary", r"^(pulmonary capillary|alveolar capillary type [12]|lung microvascular) endothelial cell$"),
    ("EC_venous", r"^vein endothelial cell of respiratory system$"),
    ("EC_arterial", r"^endothelial cell of (artery|arteriole)$"),
    ("lymphatic_EC", r"^endothelial cell of respiratory system lymphatic vessel$"),
    ("fibroblast", r"^(fibroblast of lung|adventitial cell|adventitial fibroblast|myofibroblast cell|"
                   r"secondary crest myofibroblast|pulmonary interstitial fibroblast|lung perichondrial fibroblast|"
                   r"pericyte)$"),
    ("smooth_muscle", r"^(smooth muscle cell|vascular associated smooth muscle cell|smooth muscle cell of the "
                      r"pulmonary artery|bronchial smooth muscle cell|blood vessel smooth muscle cell|"
                      r"bronchiolar smooth muscle cell)$"),
    ("macrophage", r"^(macrophage|alternatively activated macrophage|lung interstitial macrophage|"
                   r"metallothionein-positive alveolar macrophage)$"),
    ("monocyte", r"^(monocyte|intermediate monocyte|CD14-positive monocyte)$"),
    ("DC", r"dendritic cell"),
    ("mast", r"^mast cell$"),
    ("NK_ILC", r"natural killer cell|innate lymphoid cell"),
    ("T", r"\bT cell\b|thymocyte|helper T|T follicular helper"),
    ("B", r"\bB cell\b|plasma cell|plasmablast"),
    ("basal", r"^(basal cell|respiratory tract suprabasal cell|respiratory tract hillock cell)$"),
    ("multiciliated", r"^(ciliated cell|lung multiciliated epithelial cell|multiciliated epithelial cell)$"),
    ("secretory", r"^(club cell|lung secretory cell|secretory cell|goblet cell|lung goblet cell|respiratory tract goblet "
                  r"cell|bronchial goblet cell|nasal mucosa goblet cell|tracheobronchial goblet cell)$"),
)
PSYCHAD = "Population-scale cross-disorder"


def _hlca_label_map() -> dict:
    """Census cell_type -> A4 class for labels whose HLCA cells (lung tissue) all share one class."""
    import anndata as ad

    o = ad.read_h5ad(C.ATLAS_WG / "lung_hlca" / "pseudobulk.h5ad", backed="r").obs
    cls = o.ann_level_3.astype(str).map(LU.L3)
    cls[o.ann_level_2.astype(str) == "Smooth muscle"] = "smooth_muscle"
    keep = (o.tissue_sampling_method.astype(str).isin(LU.HLCA_TISSUE) & (o.tissue.astype(str) != "nose") & cls.notna()).to_numpy()
    x = pd.DataFrame(dict(ct=o.cell_type.astype(str).to_numpy(), cls=cls.to_numpy(), n=o.n_cells.to_numpy()))[keep]
    t = x.groupby(["ct", "cls"]).n.sum().unstack(fill_value=0)
    share = t.max(axis=1) / t.sum(axis=1)
    return t.idxmax(axis=1)[share >= 0.999].to_dict()


def _lung_class(labels, exact: dict) -> pd.Series:
    out = {}
    for lab in pd.unique(labels):
        c = exact.get(lab)
        if c is None:
            c = next((n for n, p in LUNG_SYN if re.search(p, lab)), None)
        out[lab] = c
    return pd.Series(labels).map(out)


def _sum(parts):
    """Majority assay over all frames, then per-dataset sparse sums, outer-joined and summed per donor x compartment x assay."""
    F = pd.concat([f.assign(_k=i) for i, (f, _, _) in enumerate(parts)], ignore_index=True)
    major = C.majority(F)
    out = []
    for f, X, names in parts:
        f = f[f.batch.to_numpy() == major.reindex(f.resection).to_numpy()]
        out.append(C.sum_groups_sparse(f, X, names, None, 1, "batch"))
    t = pd.concat(out, ignore_index=True, join="outer").fillna(0)
    return t.groupby(["resection", "compartment", "stratum"], as_index=False).sum(numeric_only=True)


def build(name: str, donors=None, genes="all"):
    atlas = name.removeprefix("gene_map_")
    sel = pd.read_csv(SELECTION)
    sel = sel[sel.atlas == atlas]
    dummy = lambda o: pd.DataFrame(index=range(len(o)))  # noqa: E731
    parts, metas, donor_tabs = [], [], []
    if atlas == "lung":
        exact = _hlca_label_map()
        classer = lambda lab: _lung_class(lab, exact)  # noqa: E731
        prefix = "lung"
        ld = LU.hlca_donors()
        o, X, names = C.read_sparse(C.ATLAS_WG / "lung_hlca", donors)
        cls = o.ann_level_3.astype(str).map(LU.L3)
        cls[o.ann_level_2.astype(str) == "Smooth muscle"] = "smooth_muscle"
        res = "hlca_full:" + o.donor_id.astype(str)
        healthy = ld.index[~ld.lung_condition.str.contains(";") & ld.lung_condition.isin(LU.HEALTHY)]
        keep = (o.tissue_sampling_method.astype(str).isin(LU.HLCA_TISSUE) & (o.tissue.astype(str) != "nose")).to_numpy() \
            & cls.notna().to_numpy() & res.isin(healthy).to_numpy()
        Fh = C.keep_majority(C.frame(o, dummy(o), "hlca_full", "lung:" + cls.astype(str), o.assay, keep=keep))
        hl = C.sum_groups_sparse(Fh, X, names, None, 1, "batch")
        metas.append(pd.DataFrame(dict(resection=Fh.resection.unique())).assign(
            collection_id="hlca_full", dataset_id="9f222629", disease_study=True, source="hlca", block=lambda t: t.resection.map(ld.study_major),
            group="hlca", collection="HLCA"))
        donor_tabs.append(C.census_donors(C.ATLAS_WG / "lung_hlca", "hlca_full"))
    else:
        classer = lambda lab: C.classify(lab, CX.BRAIN, ignore_case=True)  # noqa: E731
        prefix = "brain"
    for r in sel.itertuples():
        d = ROOT / atlas / r.dataset_id[:8]
        if not (d / "PROVENANCE.md").exists():
            continue
        o, X, names = C.read_sparse(d, donors)
        if not len(o):
            continue
        cls = pd.Series(np.asarray(classer(o.cell_type.astype(str))), index=o.index)
        F = C.frame(o, dummy(o), r.collection_id, prefix + ":" + cls.astype(str), o.assay, keep=cls.notna().to_numpy())
        if not len(F):
            continue
        parts.append((F, X, names))
        m = pd.DataFrame(dict(resection=F.resection.unique())).assign(
            collection_id=r.collection_id, dataset_id=r.dataset_id, collection=r.collection, disease_study=r.disease_study,
            source="other" if atlas == "lung" else "dlpfc")
        metas.append(m)
        t = C.census_donors(d, r.collection_id)
        donor_tabs.append(t)
    pb = _sum(parts)
    if atlas == "lung":
        pb = pd.concat([hl, pb], ignore_index=True, join="outer").fillna(0)
        pb = pb.groupby(["resection", "compartment", "stratum"], as_index=False).sum(numeric_only=True)
    meta = pd.concat(metas, ignore_index=True)
    meta = meta.drop_duplicates("resection", keep="first")
    if atlas == "dlpfc":
        psy = meta.collection.str.startswith(PSYCHAD)
        meta["block"] = np.where(psy, meta.dataset_id.str[:8], meta.collection_id)
        meta["group"] = np.where(psy, "PsychAD", "other")
    elif "block" not in meta or meta.block.isna().any():
        hb = meta.get("block")
        meta["block"] = np.where(meta.source == "other", meta.collection_id, meta.block if hb is not None else "")
        meta["group"] = np.where(meta.source == "hlca", "hlca", "other")
    meta = meta.set_index("resection")
    pb["stratum"] = pb.stratum.astype(str) + "|" + pb.resection.map(meta.block).astype(str)
    pb = pb[pb.n_cells >= 20]
    d = pd.concat(donor_tabs).loc[lambda t: ~t.index.duplicated()]
    pb = C.standard(pb)
    out = meta.reindex(pb.donor.unique()).join(d[["age", "sex"]], how="left")
    out["male"] = (out.sex.astype(str) == "male").astype(float)
    return pb, out.rename_axis("donor")


def from_legacy(name: str):
    atlas = name.removeprefix("gene_map_")
    with open(LEGACY_OUT / f"agg_{atlas}.pkl", "rb") as f:
        c = pickle.load(f)
    agg, meta = c["agg"], c["meta"]
    agg = agg.loc[:, ~agg.columns.duplicated()]
    pb = agg.rename(columns={"assay": "stratum"})
    pb["stratum"] = pb.stratum.astype(str) + "|" + pb.resection.map(meta.block).astype(str)
    pb = pb[pb.n_cells >= 20].rename(columns=lambda c_: C.canon(c_) if c_ not in ("resection", "compartment", "stratum", "n_cells", "depth") else c_)
    pb = pb.loc[:, ~pb.columns.duplicated()]
    pb = C.standard(pb)
    m = meta.copy()
    m["age"] = [C.stage_age(s) for s in m.development_stage]
    m["male"] = (m.sex.astype(str) == "male").astype(float)
    return pb, m.reindex(pb.donor.unique()).drop(columns=["development_stage"]).rename_axis("donor")
