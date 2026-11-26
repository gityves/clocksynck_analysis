"""A394 feasibility: Census obs-only scout of primary human tumour single-cell data with stroma (metadata only, no expression).

Census 2025-11-08 (as the project's pseudobulk builders), run in `census_env`. Writes
output/feasibility/census_counts.csv (dataset x disease x tissue_general) and census_summary.md.
"""
import re
import sys
from pathlib import Path

import cellxgene_census as cc
import numpy as np
import pandas as pd
import pyarrow as pa
import tiledbsoma as soma

VERSION = "2025-11-08"
OUT = Path(__file__).resolve().parents[1] / "output" / "feasibility"
OUT.mkdir(parents=True, exist_ok=True)
DATA = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian")
MIN = 50

TUMOUR = re.compile(r"neoplasm|carcinoma|cancer|tumou?r|melanoma|glioma|glioblastoma|sarcoma|blastoma|adenoma|mesothelioma", re.I)
EXCLUDE = re.compile(r"benign|leuk|lymphoma|myeloma|non-?neoplastic|hyperplasia", re.I)
MALIG = re.compile(r"malignant|neoplastic|tumou?r cell|cancer cell|glioblastoma|melanoma cell", re.I)
EPI = re.compile(r"epithelial|luminal|basal cell|hepatocyte|enterocyte|keratinocyte|club cell|goblet|pneumocyte|"
                 r"secretory|ductal|acinar|cholangiocyte|colonocyte|paneth|urothelial|prostate|enteroendocrine|tuft", re.I)
FIB = re.compile(r"fibroblast|stromal cell|stromal", re.I)
ENDO_LYMPH = re.compile(r"lymphatic", re.I)
ENDO = re.compile(r"endothelial", re.I)
MURAL = re.compile(r"pericyte|smooth muscle|mural", re.I)
GROUP = ["dataset_id", "disease", "tissue_general", "tissue", "donor_id", "cell_type", "assay", "suspension_type", "tissue_type"]

TCGA = [  # (code, disease regex, tissue regex)
    ("LUAD", r"lung adenocarcinoma|adenocarcinoma.*lung|lung.*adenocarcinoma", r"lung"),
    ("LUSC", r"squamous.*lung|lung.*squamous", r"lung"),
    ("BRCA", r"breast|ductal carcinoma|lobular", r"breast"),
    ("COAD/READ", r"colorectal|colon|rectal|rectum", r"colon|intestine|rectum|large"),
    ("KIRC", r"clear cell|renal cell carcinoma", r"kidney"),
    ("KIRP", r"papillary renal", r"kidney"),
    ("LIHC", r"hepatocellular|liver|cholangio", r"liver"),
    ("PRAD", r"prostate", r"prostate"),
    ("PAAD", r"pancrea", r"pancrea"),
    ("STAD", r"gastric|stomach", r"stomach"),
    ("ESCA", r"esophag|oesophag", r"esophag"),
    ("HNSC", r"head and neck|oral|laryn|oropharyn|pharyn|tongue", r"mouth|oral|head|tongue|pharyn|laryn"),
    ("BLCA", r"bladder|urothelial", r"bladder"),
    ("THCA", r"thyroid", r"thyroid"),
    ("UCEC", r"endometri|uterine", r"uterus|endometri"),
    ("OV", r"ovar|serous", r"ovar|fallopian"),
    ("CESC", r"cervi", r"cervi"),
    ("GBM/LGG", r"glioblastoma|glioma|astrocytoma|oligodendro", r"brain|cerebr|spinal|nervous"),
    ("SKCM", r"melanoma", r"skin"),
]


def tcga(disease, tissue):
    d, t = disease.lower(), (tissue or "").lower()
    for code, dr, tr in TCGA:
        if re.search(dr, d) and (re.search(tr, t) or code in ("GBM/LGG", "SKCM", "BRCA", "COAD/READ", "KIRC")):
            return code
    for code, dr, tr in TCGA:
        if re.search(tr, t) and re.search(r"carcinoma|cancer|neoplasm", d):
            return code + "?"  # tissue match only (generic disease label)
    return "other"


def batches(exp, filt, cols):
    ctx_cols = cols
    for tbl in exp.obs.read(value_filter=filt, column_names=ctx_cols):
        yield tbl.to_pandas()


def main():
    ctx = soma.SOMATileDBContext(tiledb_config={"soma.init_buffer_bytes": str(128 * 2**20),
                                                "sm.compute_concurrency_level": "2", "sm.io_concurrency_level": "2"})
    census = cc.open_soma(census_version=VERSION, context=ctx)
    try:
        exp = census["census_data"]["homo_sapiens"]
        dsets = census["census_info"]["datasets"].read().concat().to_pandas().set_index("dataset_id")
        # 1) distinct disease labels
        cnt = {}
        for b in batches(exp, "is_primary_data == True", ["disease"]):
            for k, v in b.disease.value_counts().items():
                cnt[k] = cnt.get(k, 0) + v
        dis = pd.Series(cnt)
        keep = sorted(d for d in dis.index if TUMOUR.search(d) and not EXCLUDE.search(d))
        print(len(dis), "distinct diseases;", len(keep), "kept", flush=True)
        pd.DataFrame({"disease": keep, "n_cells": dis[keep].values}).to_csv(OUT / "census_disease_labels_kept.csv", index=False)
        # 2) datasets with tumour cells; all their cells (incl. normal), aggregated
        flt = "is_primary_data == True and disease in [%s]" % ", ".join(repr(d) for d in keep)
        ds_ids = set()
        for b in batches(exp, flt, ["dataset_id"]):
            ds_ids |= set(b.dataset_id.unique())
        print(len(ds_ids), "datasets with tumour cells", flush=True)
        agg = []
        ids = sorted(ds_ids)
        for i in range(0, len(ids), 40):
            f = "is_primary_data == True and dataset_id in [%s]" % ", ".join(repr(x) for x in ids[i:i + 40])
            for b in batches(exp, f, GROUP):
                agg.append(b.groupby(GROUP, observed=True, dropna=False).size().rename("n").reset_index())
            print(" chunk", i, flush=True)
    finally:
        census.close()
    a = pd.concat(agg).groupby(GROUP, dropna=False, observed=True).n.sum().reset_index()
    for c in GROUP:  # categoricals -> str (else groupby counts every unused category as a donor)
        a[c] = a[c].astype(str)
    a = a[a.tissue_type == "tissue"]  # drop organoids / cell cultures
    a["is_tumour"] = a.disease.map(lambda d: bool(TUMOUR.search(d)) and not EXCLUDE.search(d))
    ct = a.cell_type
    a["malig"] = ct.str.contains(MALIG)
    a["epi"] = ct.str.contains(EPI) & ~a.malig
    a["fib"] = ct.str.contains(FIB) & ~a.malig
    a["lymph"] = ct.str.contains(ENDO_LYMPH)
    a["endo"] = ct.str.contains(ENDO) & ~a.lymph
    a["mural"] = ct.str.contains(MURAL)
    a["fib"] &= ~a.mural
    t = a[a.is_tumour]
    # per dataset: does it label malignant cells at all (else tumour epithelium stands in)
    has_mal = t.groupby("dataset_id").malig.any()
    t = t.assign(tumour_cell=np.where(t.dataset_id.map(has_mal), t.malig, t.epi))
    normal = a[(a.disease == "normal")]
    normal_donors = normal.groupby("dataset_id").donor_id.apply(set)
    rows = []
    for (ds, dis_, tg), g in t.groupby(["dataset_id", "disease", "tissue_general"]):
        dn = g.groupby("donor_id").apply(lambda x: pd.Series({
            "tum": x.n[x.tumour_cell].sum(), "fib": x.n[x.fib].sum(), "endo": x.n[x.endo].sum(),
            "mural": x.n[x.mural].sum(), "lymph": x.n[x.lymph].sum()}), include_groups=False)
        nstr = (dn[["fib", "endo", "mural"]] >= MIN).sum(axis=1)
        tumok = dn.tum >= MIN
        nd = set(normal_donors.get(ds, set()))
        meta = dsets.loc[ds]
        rows.append(dict(
            dataset_id=ds, collection_id=meta.collection_id, collection_name=meta.collection_name,
            dataset_title=meta.dataset_title, dataset_cells=int(meta.dataset_total_cell_count),
            disease=dis_, tissue_general=tg, tissues=";".join(sorted(g.tissue.unique()))[:120],
            tcga=tcga(dis_, tg), malignant_label=bool(has_mal[ds]), cells_in_group=int(g.n.sum()),
            assays=";".join(sorted(g.assay.unique()))[:80], suspension=";".join(sorted(g.suspension_type.unique())),
            n_donors=len(dn), donors_tumour50=int(tumok.sum()),
            donors_qual2=int((tumok & (nstr >= 2)).sum()), donors_qual1=int((tumok & (nstr >= 1)).sum()),
            donors_fib50=int((dn.fib >= MIN).sum()), donors_endo50=int((dn.endo >= MIN).sum()),
            donors_mural50=int((dn.mural >= MIN).sum()), donors_lymph_endo50=int((dn.lymph >= MIN).sum()),
            has_normal_in_dataset=ds in normal_donors.index,
            normal_same_donors=len(set(dn.index) & nd),
            qual1_cells=int(g.n.sum()),
        ))
    r = pd.DataFrame(rows)
    # donor_id overlap between datasets of one collection
    don = t.groupby("dataset_id").donor_id.apply(set)
    col = dsets.collection_id
    shared = []
    for ds in r.dataset_id.unique():
        others = [o for o in don.index if o != ds and col[o] == col[ds]]
        shared.append((ds, ";".join(o[:8] for o in others if don[ds] & don[o])))
    r["donor_overlap_datasets"] = r.dataset_id.map(dict(shared))
    # on disk already?
    folders = [p.name for p in DATA.iterdir()]
    prov = {}
    for p in DATA.glob("*/PROVENANCE.md"):
        try:
            prov[p.parent.name] = p.read_text(errors="ignore")
        except OSError:
            pass
    def ondisk(row):
        hits = [f for f, txt in prov.items() if row.dataset_id in txt or row.dataset_id[:8] in f]
        hits += [f for f in folders if row.dataset_id[:8] in f and f not in hits]
        return ";".join(hits)
    r["on_disk"] = r.apply(ondisk, axis=1)
    r = r.sort_values(["tcga", "donors_qual2", "dataset_cells"], ascending=[True, False, False])
    r.to_csv(OUT / "census_counts.csv", index=False)
    print(r.shape)


if __name__ == "__main__":
    main()
