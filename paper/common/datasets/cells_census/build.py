"""Per-cell clock-gene counts of CELLxGENE Census datasets (the cell-level input of A18 part 2: per-cell amplitude vs
synchrony).

Source:   Census 2025-11-08 (same version and dataset filter as the pseudobulk builders), streamed with `cellxgene_census` /
          `tiledbsoma`, which live in their own environment (`census_env`, not the `sccc` one): `build()` imports them lazily.
Output:   `$SCCC_DATA/cells/<name>/cells.parquet`: one row per cell, `observation_joinid`, Census `donor_id`, `cell_type`,
          `assay`, `tissue`, `raw_sum` (the cell's depth over all genes), the extra source columns of `DATASETS[name]` (joined on
          `observation_joinid` from the source h5ad obs over anonymous S3), and one integer column per clock gene. Every cell is
          kept (no primary-data filter, as the pseudobulks). `PROVENANCE.md` next to it.
Replaces: legacy/desync/common/extract_cells.py (+ `read_source_obs` of build_pseudobulk.py).
Used by:  A18 part 2 (`backbone/compute_per_cell.py`), A19 part 3 / V14 per-cell parts, A18 part 3 item 3 (PD).
Check:    the parquet files on disk (7 datasets) were made by the legacy script; `check_columns()` compares the schema of a
          fresh `build()` with them (needs network and `census_env`: not run here; the code is the legacy code, lines kept).

    python build.py NAME [--threads 2]        # one dataset of `DATASETS` into $SCCC_DATA/cells/NAME
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from clocksync import gtex  # noqa: E402
from common.paths import DATA  # noqa: E402

GENES = list(gtex().genes)
CENSUS_COLS = ["donor_id", "cell_type", "assay", "tissue"]
CENSUS_VERSION = "2025-11-08"
ALIAS = {"BMAL1": "ARNTL", "CHRONO": "CIART", "C1orf51": "CIART"}  # Census symbol -> panel symbol (as A152)
# name -> (Census dataset id, extra obs columns of the source h5ad)
DATASETS = {
    "hlca_full": ("9f222629-9e39-47d0-b83f-e08d610c7479",
                  ["ann_level_2", "ann_level_3", "tissue_sampling_method", "study", "lung_condition"]),
    "islets_t2d": ("e51bae9a-c747-4b64-904a-4da7cda218ab", ["Annotated_Clusters"]),
    "natri2024_ild": ("f14bc322-1322-4184-8d16-409557525ea5",
                      ["lineage", "annotation_level_1", "manual_annotation_1", "Processing_site", "Diagnosis", "Status",
                       "Sample_Type"]),
    "parkinson_regions": ("d3cb449b-c2b1-4b50-a7f1-21203535fe61", ["derived_class2", "path_braak_lb"]),
    "psychad_mssm": ("37a17b78-4864-4a42-b67b-31c00962795a", ["class", "subclass"]),
    "psychad_radc": ("5e57cd50-8e42-42d6-940d-5c1660d06864", []),
    "seaad_dlpfc": ("6f7fd0f1-a2ed-4ff1-80d3-33dde731cbc3", []),
}


def path(name: str) -> Path:
    return DATA / "cells" / name / "cells.parquet"


def load(name: str, columns=None) -> pd.DataFrame:
    """The per-cell table of a dataset (build it with `build` when missing)."""
    if not path(name).exists():
        build(name)
    return pd.read_parquet(path(name), columns=columns)


def read_source_obs(dataset_id: str, cols, version: str = CENSUS_VERSION) -> pd.DataFrame:
    """Columns of the source h5ad obs over anonymous S3 (only these arrays are fetched)."""
    import cellxgene_census as cc
    import fsspec
    import h5py
    uri = cc.get_source_h5ad_uri(dataset_id, census_version=version)["uri"]
    fs = fsspec.filesystem("s3", anon=True, default_block_size=2**22)
    out = {}
    with fs.open(uri, "rb") as f, h5py.File(f, "r") as h:
        for c in cols:
            g = h["obs"][c]
            if isinstance(g, h5py.Group):
                cats = np.array([x.decode() if isinstance(x, bytes) else x for x in g["categories"][()]], dtype=object)
                codes = g["codes"][()]
                out[c] = np.where(codes < 0, "NA", cats[np.clip(codes, 0, None)]).astype(str)
            else:
                v = g[()]
                out[c] = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in v], dtype=object) \
                    if v.dtype.kind in "SO" else v
    return pd.DataFrame(out)


def build(name: str, threads: int = 2, out: Path | None = None, version: str = CENSUS_VERSION) -> Path:
    """Stream the clock-gene counts of every cell of dataset `name` from the Census (run in `census_env`)."""
    import cellxgene_census as cc
    import tiledbsoma as soma
    dataset_id, source_cols = DATASETS[name]
    out = Path(out) if out else path(name).parent
    out.mkdir(parents=True, exist_ok=True)
    ctx = soma.SOMATileDBContext(tiledb_config={"soma.init_buffer_bytes": str(256 * 2**20),
                                                "sm.compute_concurrency_level": str(threads),
                                                "sm.io_concurrency_level": str(threads)})
    t0 = time.time()
    census = cc.open_soma(census_version=version, context=ctx)
    try:
        exp = census["census_data"]["homo_sapiens"]
        filt = f"dataset_id == '{dataset_id}'"
        obs = exp.obs.read(value_filter=filt, column_names=["soma_joinid", "observation_joinid", "raw_sum"]
                           + CENSUS_COLS).concat().to_pandas()
        if source_cols:
            src = read_source_obs(dataset_id, ["observation_joinid"] + source_cols, version)
            obs = obs.merge(src, on="observation_joinid", how="left", validate="one_to_one")
        # Census symbols of panel genes: ARNTL is BMAL1 in Census 2025-11-08 (the first version requested "ARNTL" only and
        # wrote an all-zero ARNTL column for the 7 datasets of DATASETS; fixed 2026-11-04, A203)
        names = [g for g in GENES] + list(ALIAS)
        vf = "feature_name in [%s]" % ", ".join(repr(g) for g in names)
        var = exp.ms["RNA"].var.read(value_filter=vf).concat().to_pandas()
        var["panel"] = var.feature_name.map(lambda g: ALIAS.get(g, g))
        var = var.drop_duplicates("panel")
        vlut = np.full(int(var.soma_joinid.max()) + 1, -1, dtype=np.int64)
        vlut[var.soma_joinid.to_numpy()] = [GENES.index(g) for g in var.panel]
        missing = sorted(set(GENES) - set(var.panel))
        if missing:
            raise ValueError(f"panel genes not found in the Census var: {missing}")
        lut = np.full(int(obs.soma_joinid.max()) + 1, -1, dtype=np.int64)
        lut[obs.soma_joinid.to_numpy()] = np.arange(len(obs))
        X = np.zeros((len(obs), len(GENES)), dtype=np.int32)
        q = exp.axis_query("RNA", obs_query=soma.AxisQuery(value_filter=filt), var_query=soma.AxisQuery(value_filter=vf))
        for tbl in q.X("raw").tables():
            r, c = lut[tbl["soma_dim_0"].to_numpy()], vlut[tbl["soma_dim_1"].to_numpy()]
            ok = (r >= 0) & (c >= 0)
            np.add.at(X, (r[ok], c[ok]), np.rint(tbl["soma_data"].to_numpy()[ok]).astype(np.int32))
        q.close()
    finally:
        census.close()
    df = pd.concat([obs.drop(columns=["soma_joinid"]).reset_index(drop=True), pd.DataFrame(X, columns=GENES)], axis=1)
    for c in df.columns:
        if df[c].dtype == object or str(df[c].dtype) == "category":
            df[c] = df[c].astype(str)
    df.to_parquet(out / "cells.parquet", index=False)
    (out / "PROVENANCE.md").write_text(
        f"# Per-cell clock-gene counts of Census dataset {dataset_id}\n- `paper/common/datasets/cells_census/build.py`, Census {version}, "
        f"{len(df)} cells, {df.donor_id.nunique()} donors, genes {' '.join(GENES)}; source columns {' '.join(source_cols)}; "
        f"{time.time() - t0:.0f}s.\n")
    return out / "cells.parquet"


def check_columns(name: str, fresh: Path, skip=("ARNTL",)) -> bool:
    """Same columns, cell ids and counts as the stored table (after a `build(name, out=fresh)`), genes in `skip` excepted
    (ARNTL: all zero in the first version)."""
    g = [x for x in GENES if x not in skip]
    a, b = pd.read_parquet(path(name), columns=["observation_joinid"] + g), pd.read_parquet(fresh / "cells.parquet")
    return bool(set(a.columns) <= set(b.columns) and a.observation_joinid.equals(b.observation_joinid)
                and (a[g].to_numpy() == b[g].to_numpy()).all())


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("name", choices=list(DATASETS))
    ap.add_argument("--threads", type=int, default=2)
    ap.add_argument("--out", type=Path, default=None, help="output folder (default: the stored dataset folder)")
    a = ap.parse_args()
    print(build(a.name, a.threads, out=a.out))
