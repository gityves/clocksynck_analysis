"""A226 scout 1b: CELLxGENE Census (2025-11-08 LTS) metadata only. Per dataset x donor x disease class x cell class counts of
fibroblast / endothelial cells, primary data. Writes scout/census_donor_counts.csv."""
import re, sys
import cellxgene_census as cc
import pandas as pd
OUT = sys.argv[1]
FIB = re.compile(r"fibroblast|stromal cell|myofibroblast|pericyte|smooth muscle", re.I)
END = re.compile(r"endothelial", re.I)
with cc.open_soma(census_version="2025-11-08") as C:
    obs = C["census_data"]["homo_sapiens"].obs
    ds = C["census_info"]["datasets"].read().concat().to_pandas().set_index("dataset_id")
    cols = ["dataset_id", "donor_id", "disease", "tissue_general", "cell_type", "tissue_type", "suspension_type", "assay"]
    # datasets containing any non-normal disease cell, primary
    parts = []
    it = obs.read(column_names=cols, value_filter="is_primary_data == True and disease != 'normal'").concat().to_pandas()
    dsets = sorted(it.dataset_id.unique())
    print(len(it), "non-normal primary cells;", len(dsets), "datasets", flush=True)
    del it
    full = obs.read(column_names=cols, value_filter="is_primary_data == True and dataset_id in %r" % (dsets,)).concat().to_pandas()
full["cls"] = "other"
full.loc[full.cell_type.map(lambda x: bool(FIB.search(str(x)))), "cls"] = "fibroblast"
full.loc[full.cell_type.map(lambda x: bool(END.search(str(x)))), "cls"] = "endothelial"
full["dstatus"] = (full.disease == "normal").map({True: "normal", False: "disease"})
g = full[full.cls != "other"].groupby(["dataset_id", "donor_id", "dstatus", "cls", "disease", "tissue_general", "tissue_type"], observed=True).size().rename("n").reset_index()
g["title"] = g.dataset_id.map(ds.dataset_title)
g["collection"] = g.dataset_id.map(ds.collection_name)
g["doi"] = g.dataset_id.map(ds.collection_doi)
g.to_csv(OUT, index=False)
allc = full.groupby(["dataset_id", "donor_id", "dstatus"], observed=True).size().rename("n_all").reset_index()
allc.to_csv(OUT.replace(".csv", "_allcells.csv"), index=False)
print(len(g), flush=True)
