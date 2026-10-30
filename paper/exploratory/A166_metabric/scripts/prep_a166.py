"""A166 data preparation (blind to the clock statistic): parse METABRIC (cBioPortal brca_metabric) and TCGA BRCA (recount3) into
tidy per-cohort tables cached under ../output/: {cohort}_meta.parquet (tumour metadata + proliferation and purity scores),
{cohort}_expr.parquet (samples x genes: clock, G2M, purity and pseudo-pool genes, float32).
Run: python -I prep_a166.py [metabric|tcga ...]   (reads downloaded text files only; nothing downloaded is executed)"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import a166_lib as L  # noqa: E402

D = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian")
L.OUT.mkdir(exist_ok=True)


def zscores(X, genes, groups=None):
    """A97 `_zmean` (mean over genes of gene z-scores), gene z within `groups` when given; the score is then z-scored within groups."""
    def one(x): return L.sh._zmean(x, genes)
    if groups is None:
        s = one(X); return (s - s.mean()) / s.std()
    s = X.groupby(groups.to_numpy(), group_keys=False).apply(one)
    return s.groupby(groups.to_numpy()).transform(lambda v: (v - v.mean()) / v.std())


def metabric():
    P = D / "metabric/brca_metabric"
    ex = pd.read_csv(P / "data_mrna_illumina_microarray.txt", sep="\t", index_col=0).drop(columns="Entrez_Gene_Id")
    ex = ex.groupby(level=0).mean()  # 218 symbols have >1 probe row: averaged
    floor = float(np.nanmin(ex.to_numpy()))
    X = ex.T.astype(np.float32)
    nan_cells = int(X.isna().sum().sum())
    X = X.fillna(X.mean())  # 16 NaN cells: gene mean
    pt = pd.read_csv(P / "data_clinical_patient.txt", sep="\t", comment="#")
    sm = pd.read_csv(P / "data_clinical_sample.txt", sep="\t", comment="#")
    cl = sm.merge(pt, on="PATIENT_ID").set_index("SAMPLE_ID").loc[X.index]
    m = pd.DataFrame({"patient": cl.PATIENT_ID, "cohort": cl.COHORT.astype(int), "claudin": cl.CLAUDIN_SUBTYPE,
                      "er_ihc": cl.ER_IHC, "er_status": cl.ER_STATUS}, index=X.index)
    m["er"] = m.er_ihc.map({"Positve": 1, "Positive": 1, "Negative": 0})  # the cBioPortal file spells the level "Positve"
    m["er_status_bin"] = m.er_status.map({"Positive": 1, "Negative": 0})
    m["pam50"] = m.claudin.where(m.claudin.isin(L.SUBTYPES))
    g2m = [g for g in L.sh.hallmark_g2m() if g in X.columns]
    m["prolif"] = zscores(X, g2m, m.cohort).to_numpy()
    m["purity"] = zscores(X, L.sh.PURITY, m.cohort).to_numpy()
    pool = L.pool_genes(X)
    keep = sorted(set(pool) | set(L.GENES) | set(g2m) | set(L.sh.PURITY))
    X[keep].to_parquet(L.OUT / "metabric_expr.parquet"); m.to_parquet(L.OUT / "metabric_meta.parquet")
    pd.Series(dict(floor_observed=floor, n_tumours=len(X), n_genes=X.shape[1], n_kept=len(keep), n_pool=len(pool), nan_cells_imputed=nan_cells,
                   n_g2m=len(g2m))).to_csv(L.OUT / "metabric_prep.csv", header=False)
    print("metabric", len(X), "floor", floor, "NaN cells", nan_cells, "G2M", len(g2m)); print(m.er_ihc.value_counts(dropna=False).to_dict())


def tcga():
    L.sh.FOLDER  # noqa: B018
    tum, _ = L.sh.load_tumours("BRCA")  # one primary per patient (lowest barcode), recount3 AUC-scaled counts, cached by the repo
    b = L.sh.tcga_builder()
    s = b.samples_of("BRCA")
    t = s[s.sample_type == "Primary Tumor"].sort_values("barcode").drop_duplicates("patient")
    assert list(t.barcode) == list(tum.barcode)
    ct = b.counts_of("BRCA", t)  # genome-wide scaled counts, index barcode
    depth = ct.sum(axis=1)
    assert np.allclose(depth.to_numpy(), tum.depth.to_numpy())
    lc = np.log2(ct.div(depth, axis=0) * 1e6 + 1)  # repo log-CPM convention (as load_tumours)
    g2m = [g for g in L.sh.hallmark_g2m() if g in lc.columns]
    md = b._md("BRCA", "tcga").drop_duplicates("tcga_barcode").set_index("tcga_barcode").reindex(t.barcode)
    st = pd.read_csv(D / "tcga/pancan_subtypes/brca_tcga_pan_can_atlas_2018/data_clinical_patient.txt", sep="\t", comment="#")
    st = st.drop_duplicates("PATIENT_ID").set_index("PATIENT_ID").SUBTYPE
    m = pd.DataFrame({"patient": t.patient.to_numpy(), "er_xml": md["xml_breast_carcinoma_estrogen_receptor_status"].to_numpy()}, index=t.barcode.to_numpy())
    m["er"] = m.er_xml.map({"Positive": 1, "Negative": 0})  # Indeterminate and missing out
    m["pam50"] = m.patient.map(st).str.replace("BRCA_", "", regex=False).where(lambda x: x.isin(L.SUBTYPES))
    m["cohort"] = 1
    m["prolif"] = zscores(lc, g2m).to_numpy(); m["purity"] = zscores(lc, L.sh.PURITY).to_numpy()
    assert np.corrcoef(m.prolif, tum.prolif)[0, 1] > 0.9999 and np.corrcoef(m.purity, tum.purity)[0, 1] > 0.9999  # = load_tumours up to scale
    detected = (ct > 0).mean()
    pool = [g for g in L.pool_genes(lc) if detected[g] >= 0.9]  # pool: counts > 0 in >= 90 % of tumours (as A127)
    keep = sorted(set(pool) | set(L.GENES) | set(g2m) | set(L.sh.PURITY))
    lc[keep].astype(np.float32).to_parquet(L.OUT / "tcga_expr.parquet"); m.to_parquet(L.OUT / "tcga_meta.parquet")
    pd.Series(dict(n_tumours=len(m), n_kept=len(keep), n_pool=len(pool), n_g2m=len(g2m))).to_csv(L.OUT / "tcga_prep.csv", header=False)
    print("tcga", len(m), m.er_xml.value_counts(dropna=False).to_dict(), m.pam50.value_counts(dropna=False).to_dict())


if __name__ == "__main__":
    for c in (sys.argv[1:] or ["metabric", "tcga"]):
        {"metabric": metabric, "tcga": tcga}[c]()
