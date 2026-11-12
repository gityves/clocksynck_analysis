"""Scout (design only): GSE174554 (Wang et al. 2022) snRNA-seq: per sample tumour vs normal nucleus counts, clock-gene UMI sum, MGMT UMI.

The only released nucleus annotation is Tumor / Normal (GSE174554_Tumor_normal_metadata.txt.gz); no finer cell types.
Parse only: tar members are read in memory, nothing is extracted or executed.
"""
import gzip
import io
import re
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

from clocksync.config import limit_threads
from clocksync.genes import CLOCK_GENES, gene_columns

limit_threads(4)
D = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian/geo_gse174554_wang2022")
OUT = Path(__file__).resolve().parents[1] / "output"
meta = pd.read_csv(D / "GSE174554_Tumor_normal_metadata.txt.gz", sep=r"\s+")
meta.columns = ["key", "tn"]
meta["key"] = meta.key.str.replace(r"[-.]1$", "", regex=True)  # barcode suffix is -1 in the matrices, .1 or none in the metadata
meta = meta.drop_duplicates("key").set_index("key").tn  # 8347 keys repeat with the same label



def read_mtx(raw):
    """Minimal MatrixMarket coordinate reader (pandas), robust to a malformed file (the C readers segfault on some)."""
    txt = gzip.decompress(raw)
    lines = io.BytesIO(txt)
    while True:
        pos = lines.tell()
        line = lines.readline()
        if not line.startswith(b"%"):
            break
    nr, nc, nnz = map(int, line.split())
    df = pd.read_csv(lines, sep=r"\s+", header=None, names=["i", "j", "v"], dtype={"i": np.int64, "j": np.int64, "v": np.float64})
    if len(df) != nnz:
        print("  nnz mismatch", len(df), nnz, flush=True)
    return sp.csr_matrix((df.v.values, (df.i.values - 1, df.j.values - 1)), shape=(nr, nc))


tf = tarfile.open(D / "GSE174554_RAW.tar")
names = [m.name for m in tf.getmembers() if m.name.endswith("_matrix.mtx.gz")]
rows = []
for nm in sorted(names):
    stem = re.sub(r"^GSM\d+_", "", nm).replace("_matrix.mtx.gz", "")
    if "snATAC" in stem or "Transcriptomics" in nm:
        continue
    sample = stem.replace("_batch2", "v2")
    rd = lambda suf: tf.extractfile(nm.replace("matrix.mtx.gz", suf)).read()
    bc = pd.read_csv(io.BytesIO(rd("barcodes.tsv.gz")), header=None, compression="gzip")[0].str.replace(r"[-.]1$", "", regex=True)
    ft = pd.read_csv(io.BytesIO(rd("features.tsv.gz")), header=None, sep="\t", compression="gzip")
    sym = ft.iloc[:, 1] if ft.shape[1] > 1 else ft.iloc[:, 0]
    try:
        X = read_mtx(tf.extractfile(nm).read())  # genes x cells (or transposed)
    except Exception as e:  # malformed member: skip and report
        print("SKIPPED", nm, repr(e)[:150], flush=True)
        continue
    if X.shape[0] != len(sym):
        X = X.T.tocsr()
    clock = gene_columns(sym.astype(str).tolist(), CLOCK_GENES)
    mg = gene_columns(sym.astype(str).tolist(), ["MGMT"])[0]
    tot = np.asarray(X.sum(axis=0)).ravel()
    clk = np.asarray(X[clock[clock >= 0]].sum(axis=0)).ravel()
    mgv = np.asarray(X[mg].todense()).ravel() if mg >= 0 else np.zeros(X.shape[1])
    tn = pd.Series([f"{sample}_{b}" for b in bc]).map(meta)
    for lab in ["Tumor", "Normal", None]:
        m = tn.isna().values if lab is None else (tn == lab).values
        rows.append(dict(sample=sample, file_stem=stem, cls=lab or "unmatched", cells=int(m.sum()), clock_umi=clk[m].sum(),
                         mgmt_umi=mgv[m].sum(), total_umi=tot[m].sum(), n_clock_genes=int((clock >= 0).sum())))
    print(sample, X.shape, flush=True)
df = pd.DataFrame(rows)
df["patient"] = df["sample"].str.replace(r"v2$", "", regex=True)  # inferred: 'v2' = second batch / sample of the same SF id
df.to_csv(OUT / "wang2022_sample_counts.csv", index=False)

pt = df.groupby(["patient", "cls"]).agg(cells=("cells", "sum"), clock=("clock_umi", "sum"), mgmt=("mgmt_umi", "sum")).reset_index()
res = []
mal = set(pt[(pt.cls == "Tumor") & (pt.cells >= 100)].patient)
for c in ["Tumor", "Normal"]:
    x = pt[pt.cls == c]
    big = x[x.cells >= 100]
    r = dict(study="GSE174554", cls=c, n_patients_any=int((x.cells > 0).sum()), n_patients_ge100=len(big),
             n_patients_ge100_with_tumor=len(set(big.patient) & mal), median_cells_per_patient=x.cells.median(),
             median_clock_umi_pseudobulk=big.clock.median(), median_clock_umi_per_cell=(big.clock / big.cells).median())
    if c == "Tumor":
        r.update(median_mgmt_umi_pseudobulk=big.mgmt.median(), frac_patients_mgmt_ge10=(big.mgmt >= 10).mean())
    res.append(r)
res = pd.DataFrame(res)
res.to_csv(OUT / "wang2022_design_table.csv", index=False)
print(res.T.to_string())
print("matched fraction of nuclei:", 1 - df[df.cls == "unmatched"].cells.sum() / df.groupby("sample").cells.sum().sum() * 1.0 if False else df.groupby("cls").cells.sum().to_dict())
