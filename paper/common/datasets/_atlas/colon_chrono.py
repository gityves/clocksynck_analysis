"""Mouse colon chrono-atlas (Carmona-Alcocer et al., Nat Commun 2026; Zenodo 10.5281/zenodo.20048987; 24 samples, 16 mice, known
ZT 1/7/13/19, 10x): the timed single-cell atlas used for the known-time validations (A21, A31). Shared builder of
`colon_chrono` (14 clock genes, CellBender counts), `colon_chrono_raw` (14 clock genes, Cell Ranger raw counts of the same
cells; sensitivity) and `colon_chrono_genome` (every gene, CellBender counts; input of the cell-type gene-rhythm truth test).

Legacy: desync/21_colon_injury/scripts/colon_injury.py::build (A21; `pseudobulk_<counts>.csv`) and
desync/validation/celltype_truth/scripts/build.py (A31; `pseudobulk.parquet`), cell classes `COLON_CLASSES`
(scripts/v13_celltype_truth.py) plus the published names of A21.
Raw inputs: $SCCC_DATA/colon_chrono_atlas/export (the author's Seurat object streamed without R by
desync/21_colon_injury/scripts/stream_rds.py: CellBender counts as CSC arrays + meta.csv with the author cell types) and
$SCCC_DATA/colon_chrono_atlas/raw_h5 (Cell Ranger raw matrices).

Choices: arms from (condition, region): CTRL/DC = C-DC (undamaged distal colon), DSS/DC = R-DC (distal, day 17 after DSS), DSS/PC =
R-PC (proximal colon of the R-DC mice); author cell types mapped to classes (`CLASSES`; others dropped); donor = `<arm>:<mouse>`
(the compartments of one mouse and arm share its sampling time); depth = all counts of the cell; no stratum. `colon_chrono`:
pseudobulks >= 30 cells and classes in >= 6 of the 8 mice of an arm (A21). `colon_chrono_genome`: every group, genes as canonical
symbols (first column of a duplicated symbol), the cell and mouse rules are left to the analysis (A31: >= 30 cells, >= 3
animals at every ZT). Donor table: arm, mouse, condition, region, sex, zt (h).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp

import census as C
import legacy as L

ROOT = C.DATA / "colon_chrono_atlas"
EXPORT = ROOT / "export"
COLON_CLASSES = {
    "Interstitial": "fibroblast", "Trophocytes": "fibroblast", "Telocytes": "fibroblast", "Reticular": "fibroblast",
    "MM": "smooth_muscle", "MP": "smooth_muscle", "LPM": "smooth_muscle", "Pericytes": "pericyte",
    "Endothelial": "endothelial", "Endothelial_proliferative": "endothelial", "Lymphatic": "lymphatic_EC",
    "T": "T", "T_proliferative": "T", "B": "B", "B_proliferative": "B", "Plasma": "plasma",
    "Macrophages": "macrophage", "Dendritic": "DC",
    "Colonocytes": "colonocyte", "Colonocytes_TA": "colonocyte", "Goblet": "secretory", "Goblet_TA": "secretory",
    "ISC": "stem", "ISC_proliferative": "stem",
}
CLASSES = dict(COLON_CLASSES) | {"Trans_Colonocytes": "colonocyte", "Trans_Goblet": "secretory",
                                 "TA/ISC_Proliferative": "stem", "B_Proliferative": "B", "T_Proliferative": "T",
                                 "ILC3": "ILC"}
ARMS = {("CTRL", "DC"): "C-DC", ("DSS", "DC"): "R-DC", ("DSS", "PC"): "R-PC"}
MIN_CELLS, MIN_MICE = 30, 6
LEGACY_TABLE = {"colon_chrono": "desync/21_colon_injury/output/pseudobulk_cellbender.csv",
                "colon_chrono_raw": "desync/21_colon_injury/output/pseudobulk_raw.csv",
                "colon_chrono_genome": "desync/validation/celltype_truth/output/pseudobulk.parquet"}


def _meta() -> pd.DataFrame:
    m = pd.read_csv(EXPORT / "meta.csv", index_col=0, low_memory=False)
    m["arm"] = [ARMS.get((c, r)) for c, r in zip(m.condition, m.region)]
    m["cls"] = m.celltype.map(CLASSES)
    return m


def _cellbender(cells: np.ndarray):
    """CSC matrix (genes x selected cells) of the CellBender counts and the gene list (columns read from the memory-mapped export)."""
    p = np.load(EXPORT / "counts_p.npy")
    genes = open(EXPORT / "genes.txt").read().split("\n")[:-1]
    i = np.memmap(EXPORT / "counts_i.int32", dtype="<i4", mode="r")
    x = np.memmap(EXPORT / "counts_x.float32", dtype="<f4", mode="r")
    if len(cells) > 5000:  # most of the atlas: one pass over the whole matrix
        full = sp.csc_matrix((np.asarray(x, np.float64), np.asarray(i), p), shape=(len(genes), len(p) - 1))
        return full[:, cells], genes
    cols = [sp.csc_matrix((np.asarray(x[p[c]:p[c + 1]], np.float64), np.asarray(i[p[c]:p[c + 1]]),
                           np.array([0, p[c + 1] - p[c]])), shape=(len(genes), 1)) for c in cells]
    return sp.hstack(cols).tocsc(), genes


def _raw(m: pd.DataFrame, cells: np.ndarray, genes_wanted):
    """Cell Ranger raw counts of the selected cells (names `<sample>_raw_feature_bc_matrix_<barcode>`): panel counts, depth."""
    import scanpy as sc
    from clocksync.genes import gene_columns

    names = m.index[cells]
    samp = names.str.extract(r"^(.*)_raw_feature_bc_matrix_(.*)$")
    out, depth = np.zeros((len(cells), len(genes_wanted))), np.zeros(len(cells))
    for s in samp[0].unique():
        a = sc.read_10x_h5(ROOT / "raw_h5" / f"{s}_raw_feature_bc_matrix.h5")
        a.var_names_make_unique()
        idx = np.flatnonzero(samp[0].to_numpy() == s)
        pos = pd.Series(np.arange(a.n_obs), index=a.obs_names).reindex(samp[1].to_numpy()[idx])
        if pos.isna().any():
            raise KeyError(f"{s}: barcodes absent from the raw matrix")
        X = a.X[pos.to_numpy(int)]
        depth[idx] = np.asarray(X.sum(1)).ravel()
        out[idx] = X[:, gene_columns(a.var_names, genes_wanted)].toarray()
    return out, depth


def build(name: str, donors=None, genes=None):
    """`donors`: slice, mouse ids ("F1", ...) (the part of the donor id after `<arm>:`)."""
    m = _meta()
    ok = (m.arm.notna() & m.cls.notna()).to_numpy().copy()
    if donors is not None:
        ok &= m.mouse.astype(str).isin([str(d) for d in donors]).to_numpy()
    cells = np.flatnonzero(ok)
    mm = m.iloc[cells]
    codes, groups = pd.MultiIndex.from_arrays([mm.arm, mm.mouse, mm.cls]).factorize()
    Ind = sp.csr_matrix((np.ones(len(cells)), (np.arange(len(cells)), codes)), shape=(len(cells), len(groups)))
    pb = pd.DataFrame(groups.tolist(), columns=["arm", "mouse", "compartment"])
    pb["n_cells"] = np.bincount(codes)
    if name == "colon_chrono_raw":
        panel = list(C.GENES)
        out, depth = _raw(m, cells, panel)
        P = (Ind.T @ sp.csr_matrix(out)).toarray()
        pb["depth"] = np.bincount(codes, weights=depth)
        counts = pd.DataFrame(P, columns=panel)
    else:
        X, gl = _cellbender(cells)
        depth = np.asarray(X.sum(0)).ravel()
        pb["depth"] = np.bincount(codes, weights=depth)
        P = (X @ Ind).T.tocsr()  # groups x genes
        if name == "colon_chrono":
            from clocksync.genes import canonical, gene_columns
            cols = gene_columns([canonical(g) for g in gl], C.GENES)
            if (cols < 0).any():
                raise KeyError("clock genes missing")
            counts = pd.DataFrame(P[:, cols].toarray(), columns=C.GENES)
        else:
            from clocksync.genes import canonical
            sym = np.array([canonical(g) for g in gl])
            first = np.flatnonzero(~pd.Index(sym).duplicated())
            counts = pd.DataFrame(P[:, first].toarray(), columns=sym[first])
    pb = pd.concat([pb, counts], axis=1)
    pb["donor"] = pb.arm + ":" + pb.mouse
    if name != "colon_chrono_genome" and donors is None:
        pb = pb[pb.n_cells >= MIN_CELLS]
        ok_ = pb.groupby(["arm", "compartment"]).mouse.nunique() >= MIN_MICE
        pb = pb[ok_.reindex(pd.MultiIndex.from_frame(pb[["arm", "compartment"]])).to_numpy()]
    elif name != "colon_chrono_genome":
        pb = pb[pb.n_cells >= MIN_CELLS]
    pb = pb.drop(columns=["arm", "mouse"]).reset_index(drop=True)
    return _std(pb, m)


from common.paths import ARCHIVE  # noqa: E402
def _std(pb: pd.DataFrame, m: pd.DataFrame):
    head = ["donor", "compartment", "n_cells", "depth"]
    pb = pb[head + [c for c in pb.columns if c not in head]].sort_values(["donor", "compartment"]).reset_index(drop=True)
    pb["n_cells"] = pb.n_cells.astype(int)
    mm = m[m.arm.notna()]
    d = (mm.assign(donor=mm.arm + ":" + mm.mouse).groupby("donor")
         .agg(arm=("arm", "first"), mouse=("mouse", "first"), condition=("condition", "first"), region=("region", "first"),
              sex=("sex", "first"), zt=("ZT", "first")))
    return pb, d.reindex(pb.donor.unique()).astype({"zt": float})


def from_legacy(name: str):
    f = ARCHIVE / "legacy" / LEGACY_TABLE[name]
    t = pd.read_parquet(f) if f.suffix == ".parquet" else pd.read_csv(f)
    t["donor"] = t.arm + ":" + t.mouse
    t = t.drop(columns=[c for c in ("arm", "mouse", "resection", "stratum", "zt") if c in t])
    return _std(t, _meta())
