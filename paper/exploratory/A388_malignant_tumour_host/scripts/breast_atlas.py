"""A388 feasibility, breast atlases (counts only; no clock statistic). Wu 2021 (GSE176078): authors' cell types, ER and
clinical subtype per patient. Pal 2021 (GSE161529): authors' tumour clusters (malignant), stromal classes from the
`cancer_common` marker rule (the authors give none), ER / clinical subtype from the GEO titles (A170 scout).
Rules as `cancer_sc`: >= 50 cells per pseudobulk, sample median UMI >= 1,500.  Writes output/feasibility/breast_patients.csv"""
import gzip, io, sys, warnings, zipfile
from pathlib import Path
import numpy as np, pandas as pd, scipy.sparse as sp
warnings.filterwarnings("ignore")
HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[2]
A170 = PAPER / "exploratory/A170_sc_breast_er"
sys.path.insert(0, str(A170 / "scripts")); sys.path.insert(0, str(PAPER / "common" / "datasets"))
import prep_a170 as P  # noqa: E402  (functions only)
from cancer_common.compartments import marker_compartments  # noqa: E402
OUT = HERE.parent / "output/feasibility"
MINC, MINUMI = 50, 1500


def wu():
    m = pd.read_csv(P.WU / "metadata.csv", index_col=0)
    rec = pd.read_csv(P.SCOUT / "wu_reconciled.csv").set_index("orig.ident")
    rows = []
    for p, g in m.groupby("orig.ident"):
        c = g.celltype_major.value_counts()
        rows.append(dict(cohort="Wu", patient=p, er=rec.final_ER[p], subtype=rec.final_label[p], median_umi=g.nCount_RNA.median(),
                         malignant=c.get("Cancer Epithelial", 0), fibroblast=c.get("CAFs", 0),
                         endothelial=c.get("Endothelial", 0), mural=c.get("PVL", 0)))
    return pd.DataFrame(rows)


def pal():
    z = zipfile.ZipFile(P.PAL_ZIP)
    lab = pd.read_csv(PAPER / "exploratory/scouts_2026-10-31/sc_breast_er/pooled_patient_table.csv").query("cohort=='Pal'").set_index("patient")
    rows = []
    for obj, (er, tum) in P.PAL_OBJ.items():
        md = P.R.read_rds(P.D / f"GSE161529/figshare17058077/SeuratObject_{obj}.rds")["meta.data"]
        vals = dict(zip(list(md["names"].data), md.data))
        cells = pd.DataFrame({"group": P.fac(vals["group"]), "cluster": P.fac(vals["seurat_clusters"]).astype(int)},
                             index=list(md["row.names"].data))
        for g, d in cells.groupby("group"):
            bc, sym, gi, ci, v = P.read10x(z, g.replace("_", "-"))
            pos = pd.Series(np.arange(len(bc)), index=[f"{g}_{b}" for b in bc])[d.index].to_numpy()
            X = sp.csr_matrix((v, (ci, gi)), shape=(len(bc), len(sym)))[pos]
            umi = np.asarray(X.sum(1)).ravel()
            comp = np.asarray(marker_compartments(X, sym, umi))
            tumc = d.cluster.isin(tum).to_numpy()
            c = pd.Series(comp[~tumc]).value_counts()
            rows.append(dict(cohort="Pal", patient=g, er=lab.ER_arm[g], subtype=lab.label[g], median_umi=float(np.median(umi)),
                             malignant=int(tumc.sum()), fibroblast=int(c.get("fibroblast", 0)),
                             endothelial=int(c.get("endothelial", 0)), mural=int(c.get("mural", 0))))
            print("Pal", g, rows[-1], flush=True)
    return pd.DataFrame(rows)


if __name__ == "__main__":
    t = pd.concat([wu(), pal()], ignore_index=True)
    ok = t.median_umi >= MINUMI
    for k in ("malignant", "fibroblast", "endothelial", "mural"):
        t[f"has_{k}"] = (t[k] >= MINC) & ok
    t["n_anchors"] = t[["has_fibroblast", "has_endothelial", "has_mural"]].sum(1)
    t["mal_ge2"] = t.has_malignant & (t.n_anchors >= 2)
    t["mal_ge3"] = t.has_malignant & (t.n_anchors >= 3)
    t.to_csv(OUT / "breast_patients_private_ids.csv", index=False)
