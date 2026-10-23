"""Cao et al. 2020 (Science): human fetal atlas by sci-RNA-seq3, 4,062,980 cells, 28 fetuses, 15 organs.

Source:   GEO GSE156793: `GSE156793_S3_gene_count.loom.gz` (HDF5 `matrix`, genes x cells, dense float32, chunks
          64 x 64; row attribute `gene_short_name`; column attributes sample, Fetus_id, Organ, Main_cluster_name, Assay,
          Sex, RT_group, Development_day) and the cell metadata `GSE156793_S1_metadata_cells.txt.gz`, in
          `$SCCC_DATA/natexp/cao2020/raw/`. Counts are exonic + intronic gene-assigned reads (about 0.83 x the author's
          All_reads); all values are checked to be non-negative integers.
Replaces: legacy/desync/validation/natural_experiments/scripts/build_cao.py (+ fetal_common.py), build_cao.md.
Used by:  A107b (fetal organ and within-organ gates).

HDF5 cannot be read through gzip. `build()` therefore decompresses the loom once (about 21 GB) into `scratch` (default
`$SCCC_DATA/natexp/cao2020/loom/`, deleted afterwards unless `keep_loom=True`); pass `loom=` to use a decompressed file.
The matrix is read in slabs of 6,400 cells (all genes) by up to 4 worker processes; each slab is summed to
(fetus, compartment) with `clocksync.pseudobulk` and the slab tables are added up.

Choices (as the legacy builder)
- unit (donor) = `Fetus_id`; no cell is excluded (the release is already filtered; no doublet class);
- level "organ" (default): compartment = Organ (all RT groups and clusters); level "class": `Organ:class`, the author
  `Main_cluster_name` mapped to broad classes by `CLASSMAP` (first matching regex; unmatched labels form `other`);
- stratum: the assay (sci-RNA-seq3, one value); every pseudobulk with at least one cell is kept.

Donor table (fetus): `ga_days` (Development_day), `ga_pcw` (days / 7), `sex`.
"""

from __future__ import annotations

import gzip
import os
import re
import shutil
import sys
from multiprocessing import Pool
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
import anndata as ad  # noqa: E402
import h5py  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scipy.sparse as sp  # noqa: E402

HERE = Path(__file__).resolve()
sys.path[:0] = [str(HERE.parents[1]), str(Path(__file__).resolve().parents[3])]
import clocksync as cs  # noqa: E402
from clocksync.genes import canonical, gene_columns  # noqa: E402
from common.paths import DATA  # noqa: E402
from natexp_common import GENES, compare, read_legacy, standard, sum_tables  # noqa: E402

ROOT = DATA / "natexp" / "cao2020"
LEVELS = ("organ", "class")
ASSAY = "sci-RNA-seq3"
SLAB = 6400
CLASSMAP = [  # (class, regex); first match wins; unmatched -> other
    ("endothelial", r"endothelial|^Endocardial"),
    ("stromal", r"^Stromal|^Smooth muscle|^Mesangial|^Stellate|^Epicardial fat|^Mesothelial"),
    ("immune_myeloid", r"^Myeloid|^Microglia|^Antigen presenting"),
    ("immune_lymphoid", r"^Lymphoid|^Thymocytes"),
    ("erythroid", r"^Erythroblasts"),
    ("megakaryocyte", r"^Megakaryocytes"),
    ("epithelial", r"epithelial|^Goblet|^Parietal and chief|^Ductal"),
    ("neural", r"neurons|^Astrocytes|^Oligodendrocytes|glia|^Schwann|^Granule|^Purkinje|^Unipolar brush|^Amacrine|"
               r"^Bipolar|^Ganglion|^Horizontal|^Photoreceptor|^Retinal progenitors|^Sympathoblasts"),
    ("parenchyma", r"^Adrenocortical|^Chromaffin|^Cardiomyocytes|^Hepatoblasts|^Metanephric|^Ureteric|"
                   r"^Skeletal muscle|^Satellite|^Acinar|^Islet|trophoblast|^Retinal pigment|^Lens|^Neuroendocrine"),
]


def classify(label: str) -> str:
    return next((c for c, p in CLASSMAP if re.search(p, label)), "other")


def cell_attributes(loom: Path, meta_file: Path | None = None) -> tuple[pd.DataFrame, np.ndarray, int]:
    """Column attributes of the loom as a table, the gene symbols, and the number of genes (cross-checked against the
    GEO cell metadata when given)."""
    with h5py.File(loom, "r") as f:
        ca = f["col_attrs"]
        dec = lambda k: np.char.decode(ca[k][:].astype("S"))  # noqa: E731
        obs = pd.DataFrame({k: dec(k) for k in ("sample", "Fetus_id", "Organ", "Main_cluster_name", "Sex")})
        obs["day"] = ca["Development_day"][:]
        sym = np.char.decode(f["row_attrs/gene_short_name"][:].astype("S"))
        n_genes, n_cells = f["matrix"].shape
    assert n_cells == len(obs)
    if meta_file is not None:
        m = pd.read_csv(meta_file, sep=None, engine="python", usecols=["sample", "Fetus_id", "Organ"])
        if not (len(m) == len(obs) and (m["sample"].to_numpy() == obs["sample"].to_numpy()).all()
                and (m.Fetus_id.to_numpy() == obs.Fetus_id.to_numpy()).all()):
            raise ValueError("loom and cell metadata differ in order")
    assert obs.groupby("Fetus_id").day.nunique().max() == 1
    return obs, sym, n_genes


def _slab(args):
    loom, c0, c1, obs, cols = args
    with h5py.File(loom, "r") as f:
        X = f["matrix"][:, c0:c1]
    if (X < 0).any() or (X != np.round(X)).any():
        raise ValueError(f"slab {c0}: counts are not non-negative integers")
    o = obs.iloc[c0:c1].copy()
    o["depth"] = X.sum(0, dtype=np.float64)
    o.index = o.index.astype(str)
    cells = ad.AnnData(sp.csr_matrix(X[cols, :].T.astype(np.float64)), obs=o,
                       var=pd.DataFrame(index=[canonical(g) for g in GENES]))
    return cs.pseudobulk(cells, donor="Fetus_id", compartment="compartment", depth="depth", min_cells=1)


def decompressed(gz: Path, scratch: Path) -> Path:
    out = scratch / gz.name.removesuffix(".gz")
    if not out.exists():
        scratch.mkdir(parents=True, exist_ok=True)
        with gzip.open(gz, "rb") as f, open(out, "wb") as g:
            shutil.copyfileobj(f, g, 1 << 24)
    return out


def build(level: str = "organ", loom: Path | None = None, n_jobs: int = 4, slab: int = SLAB,
          scratch: Path | None = None, keep_loom: bool = False, max_cells: int | None = None, meta_file: Path | None = None):
    """(pseudobulk, donors). `max_cells`: only the first cells of the matrix (checks); `meta_file`:
    GEO cell metadata to cross-check the loom order against (default for the downloaded loom)."""
    scratch = scratch or ROOT / "loom"
    made = loom is None
    if made:
        loom = decompressed(ROOT / "raw" / "GSE156793_S3_gene_count.loom.gz", scratch)
    try:
        obs, sym, _ = cell_attributes(loom, meta_file if meta_file is not None else
                                      ROOT / "raw" / "GSE156793_S1_metadata_cells.txt.gz" if made else None)
        obs["compartment"] = obs.Organ if level == "organ" else obs.Organ + ":" + obs.Main_cluster_name.map(classify)
        cols = gene_columns(sym, GENES)
        if (cols < 0).any():
            raise KeyError(f"clock genes missing: {np.array(GENES)[cols < 0]}")
        n = len(obs) if max_cells is None else min(max_cells, len(obs))
        tasks = [(loom, c0, min(c0 + slab, n), obs, cols) for c0 in range(0, n, slab)]
        with Pool(n_jobs) as pool:
            tables = list(pool.imap_unordered(_slab, tasks))
        pb = sum_tables(tables)
        pb["stratum"] = ASSAY
        first = obs.iloc[:n].drop_duplicates("Fetus_id").set_index("Fetus_id")
        don = pd.DataFrame({"ga_days": first.day, "ga_pcw": first.day / 7, "sex": first.Sex}).rename_axis("donor")
        return standard(pb), don.loc[sorted(set(pb.donor))]
    finally:
        if made and not keep_loom:
            Path(loom).unlink(missing_ok=True)


def from_legacy(level: str = "organ") -> tuple[pd.DataFrame, pd.DataFrame]:
    """Standard table from the legacy `pseudobulk_<level>/pseudobulk.h5ad`."""
    pb = read_legacy(ROOT / f"pseudobulk_{level}" / "pseudobulk.h5ad", extra=("ga_days", "ga_pcw", "sex"))
    don = pb.drop_duplicates("donor").set_index("donor")[["ga_days", "ga_pcw", "sex"]].rename_axis("donor")
    pb["stratum"] = ASSAY
    return standard(pb), don


if __name__ == "__main__":
    # class map against the legacy one, on the author labels of the GEO cell metadata
    lab = pd.read_csv(ROOT / "raw" / "GSE156793_S1_metadata_cells.txt.gz", sep=None, engine="python",
                      usecols=["Main_cluster_name"]).Main_cluster_name.unique()
    leg = pd.read_csv(ROOT / "class_map.csv").set_index("author_label")["class"]
    print("class map equals legacy for", len(lab), "labels:", bool((pd.Series({l: classify(l) for l in lab}) == leg.reindex(lab)).all()))
    # no decompressed loom on disk (21 GB): the reader is checked on a small synthetic loom with the same layout
    import tempfile

    rng = np.random.default_rng(0)
    n_cells, n_genes = 700, 120
    sym = np.array(GENES + [f"G{i}" for i in range(n_genes - len(GENES))])
    rng.shuffle(sym)
    labs = np.array(["Endothelial cells", "Stromal cells", "Myeloid cells", "Foo cells", "Hepatoblasts"])
    attrs = dict(sample=rng.choice(["s1", "s2"], n_cells), Fetus_id=rng.choice(["H1", "H2", "H3"], n_cells),
                 Organ=rng.choice(["Heart", "Lung"], n_cells), Main_cluster_name=rng.choice(labs, n_cells),
                 Sex=rng.choice(["M"], n_cells))
    M = rng.poisson(2.0, (n_genes, n_cells)).astype(np.float32)
    with tempfile.TemporaryDirectory() as t:
        p = Path(t) / "x.loom"
        with h5py.File(p, "w") as f:
            f["matrix"] = M
            for k, v in attrs.items():
                f[f"col_attrs/{k}"] = v.astype("S")
            f["col_attrs/Development_day"] = np.array([{"H1": 100, "H2": 110, "H3": 120}[x] for x in attrs["Fetus_id"]])
            f["row_attrs/gene_short_name"] = sym.astype("S")
        for lev in LEVELS:
            pb, don = build(lev, loom=p, n_jobs=2, slab=200)
            obs = pd.DataFrame(attrs)
            comp = obs.Organ if lev == "organ" else obs.Organ + ":" + obs.Main_cluster_name.map(classify)
            ref = pd.DataFrame(M.T, columns=sym)[GENES].assign(donor=obs.Fetus_id, compartment=comp, n_cells=1,
                                                             depth=M.sum(0)).groupby(["donor", "compartment"], as_index=False).sum()
            print(lev, "synthetic loom:", compare(pb.drop(columns="stratum"), ref[pb.drop(columns="stratum").columns]))
