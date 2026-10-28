"""A152: per-cell residency / clock counts of the lymphocyte cell types of a Census dataset (run in `census_env`).

Pulls, for the cells of the cell types in a JSON list (those of the A17 lymphocyte classes), the 7 residency markers and
the 14 clock genes (raw counts, X "raw"), `raw_sum` (depth over all genes) and the grouping columns of the atlas
pseudobulk (donor_id, tissue_general, tissue, cell_type, assay, suspension_type). Census 2025-11-08 as the pseudobulks.
Output $SCCC_DATA/cells/<name>/cells.parquet + PROVENANCE.md (public Census data; no per-donor controlled access).

    census_env/bin/python extract_cells.py NAME DATASET_ID CELLTYPES.json GENES.txt [--threads 4]
"""
import argparse
import json
import os
import time
from pathlib import Path

import numpy as np
import pandas as pd

VERSION = "2025-11-08"
COLS = ["observation_joinid", "donor_id", "tissue_general", "tissue", "cell_type", "assay", "suspension_type", "raw_sum"]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("name"); ap.add_argument("dataset_id"); ap.add_argument("celltypes"); ap.add_argument("genes")
    ap.add_argument("--threads", type=int, default=4)
    a = ap.parse_args()
    import cellxgene_census as cc
    import tiledbsoma as soma
    cts = json.loads(Path(a.celltypes).read_text())
    genes = Path(a.genes).read_text().split()
    out = Path(os.environ.get("SCCC_DATA", "/mnt/max8tb/data/CoPhaser_data/cancer_circadian")) / "cells" / a.name
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    ctx = soma.SOMATileDBContext(tiledb_config={"soma.init_buffer_bytes": str(256 * 2**20),
                                                "sm.compute_concurrency_level": str(a.threads),
                                                "sm.io_concurrency_level": str(a.threads)})
    census = cc.open_soma(census_version=VERSION, context=ctx)
    try:
        exp = census["census_data"]["homo_sapiens"]
        filt = f"dataset_id == '{a.dataset_id}' and cell_type in [{', '.join(repr(c) for c in cts)}]"
        obs = exp.obs.read(value_filter=filt, column_names=["soma_joinid"] + COLS).concat().to_pandas()
        print(len(obs), "cells", flush=True)
        alias = {"BMAL1": "ARNTL", "CHRONO": "CIART", "C1orf51": "CIART"}  # Census symbol -> panel symbol (clocksync.genes)
        vf = "feature_name in [%s]" % ", ".join(repr(g) for g in list(genes) + list(alias))
        var = exp.ms["RNA"].var.read(value_filter=vf).concat().to_pandas().drop_duplicates("feature_name")
        var["feature_name"] = var.feature_name.replace(alias)
        var = var.drop_duplicates("feature_name")
        missing = sorted(set(genes) - set(var.feature_name))
        assert not missing, missing
        order = {g: i for i, g in enumerate(genes)}
        vlut = np.full(int(var.soma_joinid.max()) + 1, -1, dtype=np.int64)
        vlut[var.soma_joinid.to_numpy()] = [order[g] for g in var.feature_name]
        lut = np.full(int(obs.soma_joinid.max()) + 1, -1, dtype=np.int64)
        lut[obs.soma_joinid.to_numpy()] = np.arange(len(obs))
        X = np.zeros((len(obs), len(genes)), dtype=np.int32)
        q = exp.axis_query("RNA", obs_query=soma.AxisQuery(value_filter=filt), var_query=soma.AxisQuery(value_filter=vf))
        for tbl in q.X("raw").tables():
            r, c = lut[tbl["soma_dim_0"].to_numpy()], vlut[tbl["soma_dim_1"].to_numpy()]
            ok = (r >= 0) & (c >= 0)
            np.add.at(X, (r[ok], c[ok]), np.rint(tbl["soma_data"].to_numpy()[ok]).astype(np.int32))
        q.close()
    finally:
        census.close()
    df = pd.concat([obs.drop(columns=["soma_joinid"]).reset_index(drop=True), pd.DataFrame(X, columns=genes)], axis=1)
    for c in df.columns:
        if df[c].dtype == object or str(df[c].dtype) in ("category", "str"):
            df[c] = df[c].astype(str)
    df.to_parquet(out / "cells.parquet", index=False)
    (out / "PROVENANCE.md").write_text(
        f"# A152 per-cell counts, Census dataset {a.dataset_id}\n- Census {VERSION}, X raw, public download via cellxgene_census "
        f"{cc.__version__} (`paper/exploratory/A152_lymphocyte_residency/scripts/extract_cells.py`), {time.strftime('%F')}.\n"
        f"- Cells {len(df)} of cell types {cts}; donors {df.donor_id.nunique()}; genes {' '.join(genes)}; "
        f"columns {' '.join(COLS)}; {time.time() - t0:.0f}s.\n- sha256 of cells.parquet: "
        f"{__import__('hashlib').sha256((out / 'cells.parquet').read_bytes()).hexdigest()}\n")
    print("DONE", round(time.time() - t0), "s", flush=True)


if __name__ == "__main__":
    main()
