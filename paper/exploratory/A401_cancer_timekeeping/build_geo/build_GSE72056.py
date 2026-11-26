"""GSE72056 (Tirosh 2016 melanoma, Smart-seq2, 19 tumours, log2(TPM/10+1)): author malignant / non-malignant calls."""
from geobuild import *
G = "GSE72056"
f = OUT / G / "raw" / "GSE72056_melanoma_single_cell_revised_v2.txt.gz"
df = pd.read_csv(f, sep="\t", index_col=0, low_memory=False)
tum = df.iloc[0].astype(int).to_numpy(); mal = df.iloc[1].astype(int).to_numpy(); nm = df.iloc[2].astype(int).to_numpy()
E = df.iloc[3:].astype(np.float32)
sym = np.array(E.index.astype(str))
nmmap = {1: "T_NK", 2: "B_plasma", 3: "myeloid", 4: "endothelial", 5: "fibroblast", 6: "T_NK"}
comp = [("malignant" if m == 2 else nmmap.get(n, "drop") if m == 1 else "drop") for m, n in zip(mal, nm)]
X = sp.csr_matrix(np.expm1(E.to_numpy().T * np.log(2)))
o = pd.DataFrame({"donor": [f"{G}|Mel{t}" for t in tum], "sample": [f"Mel{t}" for t in tum],
                  "compartment": ["tumour:" + c if c != "drop" else "drop" for c in comp], "stratum": "Smart-seq2"})
pb = pseudobulk(X, sym, o, "TPM/10 summed (back-transformed from deposited log2(TPM/10+1))")
donors = pd.DataFrame({"patient": [f"Mel{t}" for t in sorted(set(tum))]}).assign(sample=lambda d: d.patient, subtype="melanoma",
    treatment="not in GEO file", response="", time_fields="none",
    note="same Tirosh 2016 patients re-appear in GSE115978 (cohort Tirosh): do not pool both")
donors["cells"] = [int((tum == int(p[3:])).sum()) for p in donors.patient]
cm = pd.DataFrame({"author_code": ["malignant=2", "non-malignant 1", "2", "3", "4", "5", "6", "malignant=0 / non-malignant 0"],
                   "compartment": ["tumour:malignant", "T_NK", "B_plasma", "myeloid", "endothelial", "fibroblast", "T_NK", "drop"]})
finish(G, pb, donors, cm)
