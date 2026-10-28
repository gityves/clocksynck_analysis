"""Organ atlases built with the legacy `organ_fit` steps (heart, kidney, Parkinson regions, healthy organs, islets,
IBDverse, Immune Aging, Kong, endometrium, skin, blood): `Spec` = how the groups of one Census/GEO group
pseudobulk become donor x compartment pseudobulks.

Legacy: `src/sccc/atlas/organ_fit.py::groups|donor_table` and the per-folder scripts `desync/10_heart ... 19_endometrium`
(class tuples, sample filters, region maps, blocks) as captured by `desync/common/a17_tables.py`.

`groups(spec)` reproduces `organ_fit.groups`: class = first matching pattern on the Census cell_type (case-insensitive;
unmatched dropped), compartment `<organ>:<class>` (with a region map `<organ>:<region>_<class>`), the donor's majority
assay (plus the suspension type if asked) over the kept groups. `table(...)` reproduces `organ_fit.donor_table`.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd

import census as C


# Census-wide cell classes (`desync/common/run_atlas.py::COARSE`; first matching pattern on the Census cell_type)
COARSE = (
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


@dataclass
class Spec:
    path: Path                          # folder with pseudobulk.h5ad and donors.csv
    prefix: str                         # donor key prefix ("<prefix>:<donor_id>")
    organ: str                          # compartment prefix
    classes: tuple                      # ((name, regex), ...)
    keep: Callable | None = None        # obs -> boolean array of groups to use
    region: dict | None = None          # Census tissue -> short name (compartment organ:region_class)
    suspension: bool = False            # assay | suspension_type
    post: Callable | None = None        # frame -> frame, applied after the majority-assay step
    ignore_case: bool = True


def groups(spec: Spec, donors=None, genes="clock", read=None) -> pd.DataFrame:
    """Groups of the dataset as rows (resection, compartment, batch=assay, n_cells, depth, counts), majority assay."""
    o, X = read if read is not None else C.read(spec.path, genes, donors)
    cls = C.classify(o.cell_type.astype(str), spec.classes, spec.ignore_case)
    ok = cls.notna().to_numpy() & (np.ones(len(o), bool) if spec.keep is None else np.asarray(spec.keep(o)))
    if spec.region is not None:
        reg = o.tissue.astype(str).map(spec.region)
        ok &= reg.notna().to_numpy()
        cls = reg + "_" + cls
    batch = o.assay.astype(str)
    if spec.suspension:
        batch = batch + "|" + o.suspension_type.astype(str)
    F = C.frame(o, X, spec.prefix, spec.organ + ":" + cls.astype(str), batch, keep=ok)
    F = C.keep_majority(F)
    return spec.post(F) if spec.post else F


def concat_groups(specs, donors=None, genes="clock") -> pd.DataFrame:
    """Several datasets under one donor key (Kong: epithelium, immune and stromal files)."""
    return pd.concat([groups(s, donors, genes) for s in specs], ignore_index=True)


def table(F: pd.DataFrame, keep=None, min_cells: int = 20, block: pd.Series | None = None,
          min_donors: int | None = None) -> pd.DataFrame:
    """`organ_fit.donor_table`: groups summed per donor x compartment x assay; stratum = assay (| block)."""
    pb = C.donor_table(F, keep, min_cells, stratum="batch", min_donors=min_donors)
    if block is not None:
        pb["stratum"] = pb.stratum + "|" + pb.resection.map(block).astype(str)
    return pb


def labels_from(d: pd.DataFrame, col: str, control, case) -> pd.Series:
    """'control' / 'case' from a donor column: values in `control` / in `case`; other donors absent."""
    v = d[col].astype(str)
    lab = pd.Series(np.where(v.isin(control), "control", np.where(v.isin(case), "case", "")), index=d.index)
    return lab[lab != ""]


def ref_most_donors(pb: pd.DataFrame) -> str:
    return pb.groupby("compartment").resection.nunique().sort_values(ascending=False, kind="stable").index[0]


def healthy_members(pb: pd.DataFrame, members, min_donors: int) -> pd.Index:
    """Donors of a healthy fit: `members` (donors) kept in the compartments that have >= `min_donors` of them."""
    col = "resection" if "resection" in pb else "donor"
    h = pb[pb[col].isin(members)]
    n = h.groupby("compartment")[col].nunique()
    return pd.Index(h[h.compartment.isin(n[n >= min_donors].index)][col].unique())


def finish(pb: pd.DataFrame, d: pd.DataFrame, extra: dict | None = None, cols=("age", "sex", "disease")):
    """Standard (pb, donors): `d` = census_donors table; `extra` = donor columns (Series); the `healthy` column is
    a boolean (False where absent)."""
    pb = C.standard(pb)
    out = C.donors_frame(d, pb.donor.unique(), extra, cols)
    if "healthy" in out:
        out["healthy"] = out.healthy.fillna(False).astype(bool)
    return pb, out


def flag(index) -> pd.Series:
    """A True per donor of `index` (donor column of a membership flag such as `healthy`)."""
    return pd.Series(True, index=pd.Index(index))
