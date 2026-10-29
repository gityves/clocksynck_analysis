"""A126: Kumar HBCA (Census 842c6f5d, 2025-11-08) pseudobulks by donor x authors' `broad_cell_type`, panel genes only (raw counts).
Run in census_env: census_env/bin/python extract_kumar_pb.py [--threads 4]. Output OUT/kumar_pb_counts.parquet (groups x genes) +
kumar_pb_groups.csv (donor, broad, n_cells, depth = sum of raw_sum, sample_source, assay). Public data."""
import argparse, sys, time
from pathlib import Path
import numpy as np, pandas as pd
OUT = Path(__file__).resolve().parents[1] / "output"
DS = "842c6f5d-4a94-4eef-8510-8c792d1124bc"; VERSION = "2025-11-08"
ap = argparse.ArgumentParser(); ap.add_argument("--threads", type=int, default=4); a = ap.parse_args()
import cellxgene_census as cc, tiledbsoma as soma


def read_source_obs(dataset_id, cols, version=VERSION):
    """Columns of the source h5ad obs over anonymous S3 (as paper/common/datasets/cells_census/build.py)."""
    import fsspec, h5py
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
                out[c] = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in v], dtype=object) if v.dtype.kind in "SO" else v
    return pd.DataFrame(out)
panel = (OUT / "panel_genes.txt").read_text().split()
alias = {"ARNTL": "BMAL1", "CIART": "CHRONO", "BMAL1": "ARNTL", "CHRONO": "CIART"}
want = list(panel) + [alias[g] for g in panel if g in alias]
t0 = time.time()
ctx = soma.SOMATileDBContext(tiledb_config={"soma.init_buffer_bytes": str(256 * 2**20), "sm.compute_concurrency_level": str(a.threads),
                                            "sm.io_concurrency_level": str(a.threads)})
census = cc.open_soma(census_version=VERSION, context=ctx)
exp = census["census_data"]["homo_sapiens"]
filt = f"dataset_id == '{DS}'"
obs = exp.obs.read(value_filter=filt, column_names=["soma_joinid", "observation_joinid", "donor_id", "raw_sum", "assay"]).concat().to_pandas()
src = read_source_obs(DS, ["observation_joinid", "broad_cell_type", "author_cell_type"])
obs = obs.merge(src, on="observation_joinid", how="left", validate="one_to_one")
print(len(obs), "cells; broad types:", obs.broad_cell_type.value_counts().to_dict(), flush=True)
site = pd.read_csv(Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian/atlas_wg/breast_kumar/donors.csv"), usecols=["donor_id", "sample_source"]).set_index("donor_id").sample_source
obs["grp"] = obs.donor_id.astype(str) + "|" + obs.broad_cell_type.astype(str)
g = obs.grp.astype("category"); gi = g.cat.codes.to_numpy()
var = exp.ms["RNA"].var.read(value_filter="feature_name in [%s]" % ", ".join(repr(x) for x in sorted(set(want)))).concat().to_pandas()
var["name"] = [x if x in set(panel) else alias.get(x, x) for x in var.feature_name]
var = var[var.name.isin(panel)].drop_duplicates("name")
print("panel genes found", len(var), "of", len(panel), flush=True)
vl = np.full(int(var.soma_joinid.max()) + 1, -1); vl[var.soma_joinid.to_numpy()] = np.arange(len(var))
lut = np.full(int(obs.soma_joinid.max()) + 1, -1); lut[obs.soma_joinid.to_numpy()] = np.arange(len(obs))
M = np.zeros((len(g.cat.categories), len(var)), np.float64)
q = exp.axis_query("RNA", obs_query=soma.AxisQuery(value_filter=filt), var_query=soma.AxisQuery(value_filter="feature_name in [%s]" % ", ".join(repr(x) for x in var.feature_name)))
nb = 0
for tbl in q.X("raw").tables():
    r, c = lut[tbl["soma_dim_0"].to_numpy()], vl[tbl["soma_dim_1"].to_numpy()]
    ok = (r >= 0) & (c >= 0)
    np.add.at(M, (gi[r[ok]], c[ok]), np.rint(tbl["soma_data"].to_numpy()[ok]))
    nb += 1
    if nb % 50 == 0: print("batches", nb, round(time.time() - t0), "s", flush=True)
q.close(); census.close()
meta = obs.groupby("grp", observed=True).agg(donor=("donor_id", "first"), broad=("broad_cell_type", "first"), n_cells=("raw_sum", "size"),
                                              depth=("raw_sum", "sum"), assay=("assay", "first")).reindex(g.cat.categories)
meta["sample_source"] = meta.donor.map(site)
meta.to_csv(OUT / "kumar_pb_groups.csv")
pd.DataFrame(M, index=g.cat.categories, columns=var.name.to_numpy()).astype(np.float32).to_parquet(OUT / "kumar_pb_counts.parquet")
print("DONE", round(time.time() - t0), "s", flush=True)
