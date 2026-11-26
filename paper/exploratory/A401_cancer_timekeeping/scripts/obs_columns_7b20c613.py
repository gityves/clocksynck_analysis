"""A401: all obs columns of the source h5ad of Census dataset 7b20c613 (column, n unique, examples), saved with the cells' Census ids."""
import sys
import numpy as np, pandas as pd, fsspec, h5py, cellxgene_census as cc
DS = "7b20c613-9add-43d1-87e9-defd3d9b9f8c"
OUT = "/mnt/max8tb/data/CoPhaser_data/cancer_circadian/cellxgene_tumour_a401/"
uri = cc.get_source_h5ad_uri(DS, census_version="2025-11-08")["uri"]
fs = fsspec.filesystem("s3", anon=True, default_block_size=2**22)
cols = {}
with fs.open(uri, "rb") as f, h5py.File(f, "r") as h:
    o = h["obs"]
    for c in o.keys():
        try:
            g = o[c]
            if isinstance(g, h5py.Group) and "categories" in g:
                cats = np.array([x.decode() if isinstance(x, bytes) else x for x in g["categories"][()]], dtype=object)
                codes = g["codes"][()]; cols[c] = np.where(codes < 0, "NA", cats[np.clip(codes, 0, None)]).astype(str)
            elif isinstance(g, h5py.Dataset):
                v = g[()]; cols[c] = np.array([x.decode() if isinstance(x, bytes) else str(x) for x in v], dtype=object) if v.dtype.kind in "SO" else v.astype(str)
        except Exception as e:
            print("skip", c, e)
df = pd.DataFrame(cols)
rows = [dict(column=c, n_unique=df[c].nunique(), examples=" | ".join(map(str, df[c].value_counts().index[:6]))[:300]) for c in df.columns]
pd.DataFrame(rows).to_csv(OUT + "integratedcancerce_7b20c613_obs_columns.csv", index=False)
df.to_parquet("/tmp/claude-1000/obs_7b20c613_full.parquet")
print(pd.DataFrame(rows).to_string())
