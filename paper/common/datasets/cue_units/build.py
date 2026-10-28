"""Donor x tissue x coarse-class pseudobulks of the genome-wide Census atlases for the cue-following analysis (A85).

Source:   `$SCCC_DATA/atlas_wg/<atlas>/pseudobulk.h5ad` (Census group pseudobulks of all genes: donor x tissue x cell type x assay
          [x study]; the raw input of the `*_genome` builders in `_atlas/genome.py`); 9 atlases: lung_hlca, breast_kumar,
          breast_hbca, colon_crc, heart_global, liver_xu, skin_ganier, adipose_emont, muscle.
Table:    one row per unit (donor | tissue) and coarse class: counts of the requested genes (any list, e.g. the 14 core and the
          cue genes), `depth`, `n_cells`, `stratum` (assay, plus study where recorded). Healthy donors only (disease == normal
          where recorded; nuclei only where `suspension_type` is recorded), cell types mapped to `COARSE` (first match wins),
          majority assay per donor, pseudobulks of >= 20 cells, classes with >= 20 units.
Replaces: `table()` and `coarse()` of legacy/desync/validation/paired_damping/scripts/a85_cue_following.py (and `COARSE` of
          legacy/desync/common/run_atlas.py).
Used by:  A85 (cancer/compute_cue_following.py).
Check:    the legacy table is not stored; the A85 estimates computed from it equal the legacy output (see the README of cancer).
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from clocksync.genes import gene_columns  # noqa: E402
from common.paths import DATA  # noqa: E402

WG = DATA / "atlas_wg"
ATLASES = ["lung_hlca", "breast_kumar", "breast_hbca", "colon_crc", "heart_global", "liver_xu", "skin_ganier", "adipose_emont",
           "muscle"]
MIN_CELLS, MIN_UNITS = 20, 20
COARSE = (  # legacy/desync/common/run_atlas.py
    ("Treg", r"regulatory T"), ("CD4_T", r"CD4-positive"), ("CD8_T", r"CD8-positive"),
    ("T_other", r"(?-i:\bT cell)|thymocyte|T-helper|T follicular"), ("NK", r"natural killer"), ("ILC", r"innate lymphoid"),
    ("plasma", r"plasma cell|plasmablast"), ("B", r"(?-i:\bB cell)"), ("monocyte", r"monocyte"),
    ("macrophage_DC", r"macrophage|dendritic|phagocyte|Langerhans|microglial|myeloid cell|myeloid leukocyte"),
    ("mast", r"mast cell"), ("granulocyte", r"neutrophil|\bbasophil\b|granulocyte"),
    ("progenitor", r"progenitor cell$|stem cell$|precursor|myeloid progenitor"), ("erythroid", r"erythro"),
    ("endothelial", r"endothelial"),
    ("mural_smooth_muscle", r"pericyte|smooth muscle|mural|myometrial|myoid"),
    ("fibroblast_stromal", r"fibroblast|stromal|mesenchymal|stellate|adventitial|connective|tendon|keratocyte|trabecular meshwork"),
    ("muscle", r"muscle cell|myocyte|satellite"), ("adipocyte", r"adipocyte"), ("melanocyte", r"melanocyte"),
    ("glia_neuron", r"glial|Schwann|neuron|Mueller|astrocyte|oligodendrocyte"),
    ("epithelial", r"epitheli|enterocyte|goblet|basal cell|club|ciliated|alveolar|hepatocyte|cholangiocyte|acinar|"
                   r"duct|luminal|urothelial|keratinocyte|tuft|paneth|enteroendocrine|enterochromaffin|ionocyte|"
                   r"pancreatic|transit amplifying|crypt|secret|sebocyte|serous|mucus|salivary|mesothelial|colonocyte|mucous|foveolar"),
)


def coarse(cell_type: str):
    for name, pat in COARSE:
        if re.search(pat, cell_type, flags=re.I):
            return name
    return None


def table(name: str, genes) -> pd.DataFrame:
    """Unit x class counts of `genes` with depth, n_cells, stratum (genes absent from the atlas are dropped)."""
    genes = list(genes)
    a = ad.read_h5ad(WG / name / "pseudobulk.h5ad", backed="r")
    o = a.obs.copy()
    keep = np.ones(len(o), bool)
    if "disease" in o:
        keep &= (o.disease.astype(str) == "normal").to_numpy()
    if "suspension_type" in o:
        keep &= (o.suspension_type.astype(str) == "nucleus").to_numpy()
    o["cls"] = o.cell_type.astype(str).map(coarse)
    keep &= o.cls.notna().to_numpy()
    o = o[keep]
    maj = o.groupby("donor_id", observed=True).assay.agg(lambda x: x.value_counts().index[0])
    o = o[o.assay.astype(str).to_numpy() == o.donor_id.map(maj).astype(str).to_numpy()]
    names = a.var.feature_name.astype(str) if "feature_name" in a.var else a.var_names.astype(str)
    cols = gene_columns(names, genes)
    present = [g for g, c in zip(genes, cols) if c >= 0]
    cols = cols[cols >= 0]
    rows = np.flatnonzero(a.obs.index.isin(o.index))
    X = a.X[rows][:, cols]
    X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
    d = pd.DataFrame(X, columns=present, index=a.obs.index[a.obs.index.isin(o.index)])
    o = o.loc[d.index]
    d["unit"] = (o.donor_id.astype(str) + "|" + o.tissue.astype(str)).to_numpy()
    d["cls"] = o.cls.to_numpy()
    d["stratum"] = (o.assay.astype(str) + ("|" + o.study.astype(str) if "study" in o else "")).to_numpy()
    d["depth"], d["n_cells"] = o.depth.to_numpy(float), o.n_cells.to_numpy(float)
    g = d.groupby(["unit", "cls"])
    t = g[present + ["depth", "n_cells"]].sum()
    t["stratum"] = g.stratum.first()
    t = t[t.n_cells >= MIN_CELLS].reset_index()
    n = t.groupby("cls").unit.nunique()
    return t[t.cls.isin(n[n >= MIN_UNITS].index)]
