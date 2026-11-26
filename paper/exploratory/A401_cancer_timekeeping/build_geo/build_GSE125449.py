"""GSE125449 (Ma 2019 liver cancer, 10x, 12+7 tumours HCC / ICC): author cell types."""
import gzip
from geobuild import *
G = "GSE125449"
raw = OUT / G / "raw"
SERIES = {  # GEO series matrices (sample title -> cancer type); Set1 = GPL18573 (NextSeq 500), Set2 = GPL20301 (HiSeq 4000)
    "Set1": {"S16_P10_LCP18": "HCC", "S02_P01_LCP21": "HCC", "S10_P05_LCP23": "HCC", "S09_P04_LCP25": "ICC", "S08_P03_LCP26": "ICC",
             "S07_P02_LCP28": "HCC", "S11_P06_LCP29": "ICC", "S12_P07_LCP30": "HCC", "S20_P12_LCP35": "ICC", "S21_P13_LCP37": "HCC",
             "S15_P09_LCP38": "HCC", "S19_P11_LCP39": "ICC"},
    "Set2": {"S351_P10_LCP34": "HCC", "S355_P13_LCP42": "ICC", "S358_P16_LCP46": "ICC", "S305_P06_LCP56": "ICC", "S300_P02_LCP60": "ICC",
             "S364_P21_LCP65": "HCC", "S365_P22_LCP66": "ICC"}}
MAP = {"Malignant cell": "malignant", "CAF": "fibroblast", "TEC": "endothelial", "TAM": "myeloid", "T cell": "T_NK",
       "B cell": "B_plasma", "HPC-like": "HPC_like", "unclassified": "drop"}
items, meta, cm = [], [], []
for s in ("Set1", "Set2"):
    X, sym, bc = read_10x(raw / f"GSE125449_{s}_")
    smp = pd.read_csv(raw / f"GSE125449_{s}_samples.txt.gz", sep="\t")
    assert len(smp) == X.shape[0] and (smp["Cell Barcode"].to_numpy() == bc).all()
    for name in smp.Sample.unique():
        m = (smp.Sample == name).to_numpy()
        auth = np.array([MAP[t] for t in smp.Type[m]], dtype=object)
        pat = name.split("_")[1]
        items.append((name, pat, X[np.flatnonzero(m)], sym, bc[m], auth))
        meta.append(dict(patient=pat, sample=name, set=s, subtype=SERIES[s][name], stratum="10x_NextSeq500" if s == "Set1" else "10x_HiSeq4000"))
    cm.append(pd.crosstab(smp.Type.rename("a"), smp.Type.map(MAP).rename("b")).stack().rename("n").reset_index())
strat = {m["sample"]: m["stratum"] for m in meta}
X, sym, obs = by_sample(items, G, lambda n: "tumour", lambda n: strat[n])
pb = pseudobulk(X, sym, obs, "counts")
donors = pd.DataFrame(meta).assign(treatment="not stated (GEO: treatment-naive resections per paper not verified)", response="", time_fields="none")
donors["tissue"] = "tumour"
cm = pd.concat(cm).drop_duplicates().rename(columns={"a": "author_type", "b": "compartment"}); cm = cm[cm.n > 0]
cm["compartment"] = "tumour:" + cm.compartment; cm["note"] = "HPC-like kept apart (malignant sub-state in the paper); unclassified dropped"
finish(G, pb, donors.drop(columns="stratum"), cm)
