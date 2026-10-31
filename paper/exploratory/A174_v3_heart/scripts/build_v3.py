"""A174 (V3 port): donor x compartment clock-gene pseudobulks of GSE270258 (Jia 2024), human left-ventricle snRNA-seq of 16
patients operated at 06-10 h (AM) or 14-18 h (PM). No scoring here: gate, host fit and AM/PM tests are separate.

Source: `$SCCC_DATA/validation/GSE270258/GSE270258.h5ad` (raw counts, no author labels) and `samples.csv` (window, arm).
Legacy builder replaced: legacy/scripts/v3_heart.py (`assign`, pseudobulk step), same rule and parameters.

Cell classes. The COARSE rules match Census cell-type names, which this dataset lacks, so cells are called from markers,
frozen here: score of a class = mean log1p(CP10K) over its marker genes (CP10K from the genome-wide UMI count); a cell
takes the best class if best > MIN_SCORE and second < MAX_RATIO * best, else it is dropped. No clock gene is a marker.
Class names follow `common/datasets/_atlas/heart.py::CLASSES`; compartment `heart:<class>`. Pericyte and smooth muscle
are merged into `mural`, lymphoid is `lymphocyte` (the legacy kept them separate; they are merged in the standard heart
vocabulary). Pseudobulks of >= MIN_CELLS cells (repo standard 20; the legacy used 50: filter on `n_cells` to compare).
Donor = patient (one library each), prefixed `gse270258:`. The rule stamp (`r1-<hash>` of the constants and the
normalised source, as common/data.py) is written to the parquet metadata under common.data.RULES_KEY.

    python scripts/build_v3.py     # -> output/v3_heart.parquet, v3_heart_donors.parquet, cells_per_class.csv, labels.csv.gz
"""

from __future__ import annotations

import hashlib
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
import scipy.sparse as sp

PAPER = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PAPER))
import anndata as ad  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import data as D  # noqa: E402
from common.paths import DATA, exploratory_dir  # noqa: E402

SRC = DATA / "validation" / "GSE270258"
MARKERS = {
    "cardiomyocyte": ["TTN", "RYR2", "MYH7", "TNNT2", "SORBS2", "TECRL", "CTNNA3", "FHL2"],
    "fibroblast": ["DCN", "LUM", "PDGFRA", "ABCA8", "BICC1", "NEGR1", "C7", "FBLN1"],
    "endothelial": ["VWF", "PECAM1", "CDH5", "EGFL7", "FLT1", "EMCN"],
    "mural": ["RGS5", "PDGFRB", "KCNJ8", "ABCC9", "MYH11", "ACTA2", "TAGLN", "CNN1"],
    "myeloid": ["CD163", "F13A1", "MRC1", "CD14", "MS4A6A", "LYVE1"],
    "lymphocyte": ["PTPRC", "CD3E", "CD2", "IL7R", "NKG7", "SKAP1"],
    "adipocyte": ["ADIPOQ", "PLIN1", "GPAM", "PPARG"],
    "neuron": ["NRXN1", "NRXN3", "XKR4", "PLP1"],
}
MIN_SCORE, MAX_RATIO, MIN_CELLS, CHUNK = 0.3, 0.7, 20, 50_000


def assign(path: Path) -> pd.DataFrame:
    """Per-cell class (NaN if unassigned) by marker score; reads the h5ad backed in row chunks."""
    a = ad.read_h5ad(path, backed="r")
    genes = [g for gs in MARKERS.values() for g in gs]
    cols = cs.genes.gene_columns(a.var_names, genes, missing="raise")
    pos = {g: i for i, g in enumerate(genes)}
    names, out = list(MARKERS), []
    for s in range(0, a.n_obs, CHUNK):
        X = sp.csr_matrix(a.X[s:s + CHUNK])[:, cols].tocsr()
        X = X.multiply(1e4 / np.maximum(a.obs.total_umi.iloc[s:s + CHUNK].to_numpy(), 1)[:, None]).tocsr()
        X.data = np.log1p(X.data)
        S = np.column_stack([np.asarray(X[:, [pos[g] for g in MARKERS[n]]].mean(1)).ravel() for n in names])
        o = np.argsort(-S, 1)
        best, second = S[np.arange(len(S)), o[:, 0]], S[np.arange(len(S)), o[:, 1]]
        out.append(np.where((best > MIN_SCORE) & (second < MAX_RATIO * best), np.array(names)[o[:, 0]], None))
    lab = pd.DataFrame({"patient": a.obs.patient.astype(str).to_numpy(), "total_umi": a.obs.total_umi.to_numpy(),
                        "compartment": np.concatenate(out)}, index=a.obs_names.astype(str))
    a.file.close()
    lab["compartment"] = lab.compartment.map(lambda c: f"heart:{c}" if isinstance(c, str) else None)
    return lab


def build() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    lab = assign(SRC / "GSE270258.h5ad")
    pb = cs.pseudobulk(SRC / "GSE270258.h5ad", "patient", "compartment", labels=lab[lab.compartment.notna()],
                       depth="total_umi", min_cells=MIN_CELLS)
    pb["donor"] = "gse270258:" + pb.donor
    s = pd.read_csv(SRC / "samples.csv").set_index("patient")
    donors = s[["group", "window", "window_start_h", "window_end_h", "exact_crossclamp_time", "sex", "age"]].copy()
    donors.index = "gse270258:" + donors.index
    donors["window_mid_h"] = (donors.window_start_h + donors.window_end_h) / 2
    return pb, donors, lab


def stamp() -> str:
    """Rule fingerprint, same scheme as `common.data.rules_hash` (constants + normalised source of this file)."""
    h = hashlib.sha256()
    h.update(repr(sorted((k, repr(D._plain(v))) for k, v in globals().items() if k.isupper() and D._plain(v) is not None)).encode())
    h.update(D._normalised_source(Path(__file__)).encode())
    return f"{D._SCHEME}-{h.hexdigest()[:16]}"


def main():
    limit_threads(4)
    out = exploratory_dir("A174_v3_heart")
    pb, donors, lab = build()
    for df, f, idx in [(pb, "v3_heart.parquet", False), (donors.rename_axis("donor"), "v3_heart_donors.parquet", True)]:
        t = pa.Table.from_pandas(df, preserve_index=idx)
        pq.write_table(t.replace_schema_metadata({**(t.schema.metadata or {}), D.RULES_KEY: stamp().encode()}), out / f)
    lab.to_csv(out / "labels.csv.gz")
    n = lab.assign(compartment=lab.compartment.fillna("unassigned")).pivot_table(
        index="patient", columns="compartment", values="total_umi", aggfunc="size", fill_value=0)
    n.to_csv(out / "cells_per_class.csv")
    print(n.to_string(), "\n", f"stamp {stamp()}  rows {len(pb)}  ->", out)


if __name__ == "__main__":
    main()
