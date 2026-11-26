"""GSE103322 (Puram 2017 HNSCC, Smart-seq2, 18 patients, primary + lymph-node metastasis): author cancer / non-cancer labels.
Deposited values are log2(TPM/10+1); back-transformed to TPM/10 and summed (units column)."""
from geobuild import *
G = "GSE103322"
f = OUT / G / "raw" / "GSE103322_HNSCC_all_data.txt.gz"
df = pd.read_csv(f, sep="\t", index_col=0, low_memory=False)
h = df.iloc[:5].T
h.columns = ["maxima", "ln", "cancer", "noncancer", "ctype"]
E = df.iloc[5:].astype(np.float32)
sym = np.array([s.strip("'") for s in E.index])
pat = h.index.str.extract(r"^(?:HNSCC_?|HN)(\d+)", expand=False).to_numpy()
nopat = pd.isna(pat)  # 109 "combo" plate cells carry no patient id: dropped
ctype = h.ctype.astype(str).str.lstrip("-")
m = {"Fibroblast": "fibroblast", "Endothelial": "endothelial", "T cell": "T_NK", "B cell": "B_plasma",
     "Macrophage": "myeloid", "Dendritic": "myeloid"}
comp = np.where(h.cancer.astype(int) == 1, "malignant", ctype.map(m).fillna("drop").to_numpy())
tissue = np.where(h.ln.astype(int) == 1, "metastasis", "tumour")
lab = np.array([f"{t}:{c}" if c != "drop" else "drop" for t, c in zip(tissue, comp)])
lab[nopat] = "drop"
pat = np.where(nopat, "NA", pat)
X = sp.csr_matrix(np.expm1(E.to_numpy().T * np.log(2)))  # 2^x - 1 = TPM/10
o = pd.DataFrame({"donor": [f"{G}|HN{p}" for p in pat], "sample": [f"HN{p}_{'LN' if l else 'primary'}" for p, l in zip(pat, h.ln.astype(int))],
                  "compartment": lab, "stratum": "Smart-seq2"})
pb = pseudobulk(X, sym, o, "TPM/10 summed (back-transformed from deposited log2(TPM/10+1))")
donors = (o.groupby(["donor", "sample"]).size().rename("cells").reset_index().assign(
    patient=lambda d: d.donor.str.split("|").str[1], subtype="HNSCC", treatment="not stated in GEO", response="", time_fields="none",
    tissue=lambda d: np.where(d["sample"].str.endswith("LN"), "lymph-node metastasis", "primary tumour")))
cm = pd.crosstab(pd.Series(ctype.to_numpy()).where(h.cancer.astype(int).to_numpy() == 0, "(cancer cell)"), lab).stack().reset_index()
cm.columns = ["author_label", "compartment", "n_cells"]; cm = cm[cm.n_cells > 0]
cm["note"] = "author cancer-cell call -> malignant; Mast, myocyte and unlabelled non-cancer dropped; LN cells use prefix metastasis:"
finish(G, pb, donors, cm)
