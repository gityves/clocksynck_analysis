"""A394 feasibility: nnz (stored non-zeros) of the tumour cells of datasets with qualifying donors, to size a counts download."""
from pathlib import Path
import cellxgene_census as cc, pandas as pd, tiledbsoma as soma
OUT = Path(__file__).resolve().parents[1] / "output" / "feasibility"
r = pd.read_csv(OUT / "census_counts.csv")
ids = sorted(r[r.donors_qual1 > 0].dataset_id.unique())
dis = pd.read_csv(OUT / "census_disease_labels_kept.csv").disease.tolist()
ctx = soma.SOMATileDBContext(tiledb_config={"sm.compute_concurrency_level": "2", "sm.io_concurrency_level": "2"})
c = cc.open_soma(census_version="2025-11-08", context=ctx)
exp = c["census_data"]["homo_sapiens"]
f = "is_primary_data == True and dataset_id in [%s] and disease in [%s]" % (", ".join(map(repr, ids)), ", ".join(map(repr, dis)))
acc = []
for t in exp.obs.read(value_filter=f, column_names=["dataset_id", "nnz"]):
    acc.append(t.to_pandas().groupby("dataset_id").nnz.agg(["size", "sum"]))
s = pd.concat(acc).groupby(level=0).sum().rename(columns={"size": "tumour_cells", "sum": "nnz"})
s.to_csv(OUT / "census_nnz_by_dataset.csv"); print(s.sum())
