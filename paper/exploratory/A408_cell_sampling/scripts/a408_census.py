"""A408: per-cell counts of the 50 companion sets' genes for the target-class cells of the A17 healthy HLCA donors (CELLxGENE
Census query; run in `census_env`, which has `cellxgene_census` / `tiledbsoma`).

Input  (`../output/census_request/`): cells.csv (`observation_joinid` of the cells, from `a408.py --stage prep`), genes.csv
       (gene symbols; canonical, BMAL1 -> ARNTL).
Output `$SCCC_DATA/cells/hlca_a408_companions/`: counts.npz (scipy CSR, cells x genes, int32), cells.csv (row order),
       genes.csv (column order), PROVENANCE.md (Census version, dataset id, query, date, sha256 of counts.npz).

    ~/miniconda3/envs/census_env/bin/python a408_census.py
"""
import hashlib
import sys
import time
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

HERE = Path(__file__).resolve()
REQ = HERE.parents[1] / "output" / "census_request"
DEST = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian/cells/hlca_a408_companions")
VERSION, DATASET = "2025-11-08", "9f222629-9e39-47d0-b83f-e08d610c7479"
ALIAS = {"BMAL1": "ARNTL", "CHRONO": "CIART", "C1orf51": "CIART"}


def main():
    import cellxgene_census as cc
    import tiledbsoma as soma
    t0 = time.time()
    want_cells = pd.read_csv(REQ / "cells.csv").observation_joinid.astype(str)
    genes = pd.read_csv(REQ / "genes.csv", header=None)[0].astype(str).tolist()
    ctx = soma.SOMATileDBContext(tiledb_config={"soma.init_buffer_bytes": str(256 * 2**20),
                                                "sm.compute_concurrency_level": "4", "sm.io_concurrency_level": "4"})
    census = cc.open_soma(census_version=VERSION, context=ctx)
    try:
        exp = census["census_data"]["homo_sapiens"]
        filt = f"dataset_id == '{DATASET}'"
        obs = exp.obs.read(value_filter=filt, column_names=["soma_joinid", "observation_joinid"]).concat().to_pandas()
        obs["observation_joinid"] = obs.observation_joinid.astype(str)
        obs = obs[obs.observation_joinid.isin(set(want_cells))].set_index("observation_joinid").loc[want_cells.tolist()]
        # genes: Census feature_name, mapped to the canonical symbol (as `clocksync.genes.canonical`)
        var = exp.ms["RNA"].var.read(column_names=["soma_joinid", "feature_name"]).concat().to_pandas()
        var["sym"] = var.feature_name.map(lambda g: ALIAS.get(g, g.upper() if g.upper() in set(genes) else g))
        var = var[var.sym.isin(genes)].drop_duplicates("sym")
        missing = sorted(set(genes) - set(var.sym))
        print(f"cells {len(obs)}/{len(want_cells)}, genes {len(var)}/{len(genes)}, missing {missing}", flush=True)
        col = {s: i for i, s in enumerate(genes)}
        vlut = np.full(int(var.soma_joinid.max()) + 1, -1, dtype=np.int64)
        vlut[var.soma_joinid.to_numpy()] = [col[s] for s in var.sym]
        sj = obs.soma_joinid.to_numpy()
        lut = np.full(int(sj.max()) + 1, -1, dtype=np.int64)
        lut[sj] = np.arange(len(sj))
        # contiguous soma_joinid slices (a query on scattered coordinates is ~40x slower); rows of unrequested cells dropped
        rows, cols, vals = [], [], []
        vco = np.sort(var.soma_joinid.to_numpy())
        lo, hi, step = int(sj.min()), int(sj.max()), 100_000
        for a0 in range(lo, hi + 1, step):
            sub = sj[(sj >= a0) & (sj < a0 + step)]
            if not len(sub):
                continue
            q = exp.axis_query("RNA", obs_query=soma.AxisQuery(coords=(slice(a0, min(a0 + step, hi + 1) - 1),)),
                               var_query=soma.AxisQuery(coords=(vco,)))
            for tbl in q.X("raw").tables():
                r, c = lut[tbl["soma_dim_0"].to_numpy()], vlut[tbl["soma_dim_1"].to_numpy()]
                ok = (r >= 0) & (c >= 0)
                rows.append(r[ok]), cols.append(c[ok]), vals.append(np.rint(tbl["soma_data"].to_numpy()[ok]).astype(np.int32))
            q.close()
            print(f"slice from {a0} {time.time() - t0:.0f}s", flush=True)
    finally:
        census.close()
    X = sp.csr_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                      shape=(len(obs), len(genes)), dtype=np.int32)
    DEST.mkdir(parents=True, exist_ok=True)
    sp.save_npz(DEST / "counts.npz", X)
    pd.Series(obs.index).to_csv(DEST / "cells.csv", index=False, header=["observation_joinid"])
    pd.Series(genes).to_csv(DEST / "genes.csv", index=False, header=False)
    sha = hashlib.sha256((DEST / "counts.npz").read_bytes()).hexdigest()
    size = (DEST / "counts.npz").stat().st_size
    (DEST / "PROVENANCE.md").write_text(
        f"# A408 companion per-cell counts (HLCA)\n- CELLxGENE Census {VERSION} (homo_sapiens, X \"raw\"), dataset {DATASET} "
        f"(HLCA core and extended); downloaded {date.today()} via cellxgene_census {cc.__version__}.\n"
        f"- Query: obs coords = the {len(obs)} cells of `{REQ / 'cells.csv'}` (target-class cells of the A17 healthy donors with "
        f">= 400 cells of the class, from `paper/exploratory/A408_cell_sampling/scripts/a408.py --stage prep`); var = the "
        f"{len(var)} genes of `genes.csv` (50 companion sets of 14 genes, A263 matching, seed [408, 1]; Census feature_name, "
        f"BMAL1 -> ARNTL). Missing in Census: {missing or 'none'}.\n"
        f"- Script `paper/exploratory/A408_cell_sampling/scripts/a408_census.py`; {time.time() - t0:.0f} s; counts.npz "
        f"{size / 1e6:.1f} MB, {X.nnz} non-zeros.\n- sha256 of counts.npz: {sha}\n- Untrusted data: parse only.\n")
    print(f"saved {size / 1e6:.1f} MB nnz {X.nnz} sha {sha} {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
