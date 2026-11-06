"""A226 step 2: build per-cohort h5ad (raw integer counts; obs sample, patient, sample_type, author_cell_type, sorting,
treatment, total_umi, chemistry) for the new paired tumour / adjacent-normal cohorts, in `$SCCC_DATA/a226/<study>/<study>.h5ad`.
Same format and legacy helpers as `A149_paired_tumour_normal/scripts/build_geo_a149_part1.py`; files are parsed, never executed.

    python build_a226.py GSE308103 GSE183219 GSE155109 cxg_b6b5ea88 cxg_5af90777
Rules as A149 note 1: only patients with a tumour and an adjacent-normal sample are kept; integer counts asserted.
"""
import gzip, io, re, sys, tarfile
from pathlib import Path

sys.path.insert(0, "/home/maxine/Documents/paychere/scCircadianCancer/legacy/scripts/geo")
import anndata as ad  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scipy.sparse as sp  # noqa: E402
import build_geo as B  # noqa: E402
from scipy.io import mmread  # noqa: E402

A226 = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian/a226")
B.ROOT = A226                       # fetch / soft_samples / finalize write and read under a226/
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "A149_paired_tumour_normal/scripts"))


def _pair_stats(obs):
    g = obs.assign(sample_type=obs.sample_type.astype(str), patient=obs.patient.astype(str)).groupby("patient").sample_type.agg(set)
    return int(sum(("tumour" in s) and ("adjacent_normal" in s) for s in g)), int(len(g))


B.pair_stats = _pair_stats


def keep_paired(obs):
    t = obs.groupby("patient").sample_type.agg(set)
    ok = {p for p, s in t.items() if {"tumour", "adjacent_normal"} <= s}
    return obs.patient.isin(ok).to_numpy(), sorted(set(t.index) - ok)


def assemble(parts):
    mats, gl, ob = [], [], []
    for M, gn, bc, spec in parts:
        o = pd.DataFrame(index=[f"{spec['sample']}_{b}" for b in bc])
        for k, v in spec.items():
            o[k] = v
        mats.append(M); gl.append(gn); ob.append(o)
        print(f"  {spec['sample']} ({spec['patient']}, {spec['sample_type']}) {M.shape[0]} cells", flush=True)
    X, genes, _ = B.stack_samples(mats, gl)
    return X, genes, pd.concat(ob)


def read_dense_member(tf, member, tmp):
    with open(tmp, "wb") as f:
        f.write(tf.extractfile(member).read())
    try:
        return B.read_dense_gene_by_cell(tmp)
    finally:
        tmp.unlink()


def GSE308103():
    gse = "GSE308103"
    recs = B.fetch_all(gse, [f"{gse}_RAW.tar"])
    soft = B.soft_samples(gse)
    tf = tarfile.open(A226 / gse / f"{gse}_RAW.tar")
    mem = {re.match(r"(GSM\d+)_", m.name).group(1): m for m in tf.getmembers()}
    parts, dropped_lesions = [], []
    for gsm, v in soft.items():
        lib = [x for x in v["other"] if "Library name" in x][0].split("Library name:")[1].strip()
        m = re.fullmatch(r"(P\d+)_(Normal|LUAD|AIS|AAH|MIA)[-]?(\d*)", lib)
        assert m, lib
        if m.group(2) not in ("Normal", "LUAD"):
            dropped_lesions.append(lib)
            continue
        X, genes, cells = read_dense_member(tf, mem[gsm], A226 / gse / "_tmp_member.txt.gz")
        spec = dict(sample=lib, patient=m.group(1), sample_type="tumour" if m.group(2) == "LUAD" else "adjacent_normal",
                    chemistry="snRNA FFPE 10x", lesion=m.group(2))
        parts.append((X, genes, cells, spec))
    X, genes, obs = assemble(parts)
    obs["author_cell_type"], obs["sorting"], obs["treatment"] = "", "", "naive"
    k, dropped = keep_paired(obs)
    X, obs = X[np.flatnonzero(k)], obs[k]
    notes = ["Only the libraries `<patient>_Normal[1]` (adjacent_normal) and `<patient>_LUAD[1]` (tumour: invasive LUAD) are used; precursor "
             f"lesion libraries (AAH, AIS, MIA: {len(dropped_lesions)}) are dropped before any statistic (A226 choice; A149 kept all stages).",
             "patient = library prefix P<n> (GSM title 'patient n'); treatment: naive (series summary: treatment-naive lungs); FFPE single nuclei; "
             "raw counts per library as deposited (dense genes x cells), no author cell-type labels.",
             f"Patients without both a Normal and a LUAD library dropped: {dropped or 'none'}."]
    return B.finalize(gse, X, genes, obs, B._prov(recs), notes, title="Early LUAD snRNA-seq: normal lung and invasive LUAD of the same patients")


def call_cells(M, n_expected=5000, min_genes=200):
    """Cell Ranger v2-style ordmag rule, as `build_geo_a149_part1.call_cells`."""
    umi = np.asarray(M.sum(axis=1)).ravel()
    top = np.sort(umi)[::-1][:n_expected]
    cut = max(top[int(0.01 * len(top))] / 10.0, 1.0)
    return (umi >= cut) & (np.diff(M.indptr) >= min_genes), float(cut)


def GSE183219():
    gse = "GSE183219"
    recs = B.fetch_all(gse, [f"{gse}_RAW.tar"])
    soft = B.soft_samples(gse)
    tf = tarfile.open(A226 / gse / f"{gse}_RAW.tar")
    mem = {re.match(r"(GSM\d+)_", m.name).group(1): m for m in tf.getmembers() if m.name.endswith(".tar.gz")}
    parts = []
    for gsm, v in soft.items():
        m = re.fullmatch(r"scRNAseq_patient_(\d+)_(adjacent|tumor)", v["title"])
        assert m, v["title"]
        inner = tarfile.open(fileobj=io.BytesIO(tf.extractfile(mem[gsm]).read()), mode="r:gz")
        names = {Path(n).name: n for n in inner.getnames()}
        bc = pd.read_csv(inner.extractfile(names["barcodes.tsv"]), header=None)[0].astype(str).tolist()
        ft = pd.read_csv(inner.extractfile(names["features.tsv"]), header=None, sep="\t", dtype=str)
        M = B._mm_to_cells_by_genes(inner.extractfile(names["matrix.mtx"]))
        assert M.shape == (len(bc), len(ft))
        if ft.shape[1] > 2:
            keep = np.flatnonzero((ft[2] == "Gene Expression").to_numpy()); M, ft = M[:, keep], ft.iloc[keep]
        mask, cut = call_cells(M)       # raw (unfiltered) droplet matrices: label-free cell calling per library
        print(f"  {v['title']}: {M.shape[0]} barcodes, {int(mask.sum())} called (UMI >= {cut:.0f})", flush=True)
        M, bc = M[mask], [b for b, k in zip(bc, mask) if k]
        spec = dict(sample=v["title"], patient="p" + m.group(1), sample_type="tumour" if m.group(2) == "tumor" else "adjacent_normal",
                    chemistry="10x 3'")
        parts.append((M, ft[1].tolist(), bc, spec))
    X, genes, obs = assemble(parts)
    obs["author_cell_type"], obs["sorting"], obs["treatment"] = "", "stromal (CD45-depleted, sorted stromal cells)", "unknown"
    k, dropped = keep_paired(obs)
    X, obs = X[np.flatnonzero(k)], obs[k]
    notes = ["GSM title scRNAseq_patient_<n>_(adjacent|tumor); `tissue type`/`cell type: sorted stromal cells` from the SOFT file. Cell Ranger MEX per GSM.",
             "Stromal-sorted fractions (no epithelium: unusable for A149, fine for stromal pairing). sorting is recorded as stromal, not CD45+, so the A149 "
             "CD45+ exclusion does not remove it. treatment not stated in the record -> unknown.", f"Unpaired patients dropped: {dropped or 'none'}."]
    return B.finalize(gse, X, genes, obs, B._prov(recs), notes, title="NSCLC sorted stromal cells, tumour and adjacent")


def GSE155109():
    gse = "GSE155109"
    fn = [f"{gse}_bcTax_{x}.csv.gz" for x in ("EC_enriched_metadata", "EC_enriched_raw_matrix", "stroma_metadata", "stroma_raw_matrix")]
    recs = B.fetch_all(gse, fn)
    parts = []
    for frac in ("EC_enriched", "stroma"):
        meta = pd.read_csv(A226 / gse / f"{gse}_bcTax_{frac}_metadata.csv.gz", index_col=0, sep=";")
        X, genes, cells = B.read_dense_gene_by_cell(A226 / gse / f"{gse}_bcTax_{frac}_raw_matrix.csv.gz", sep=",")
        assert meta.index.is_unique
        miss = [c for c in cells if c not in meta.index]
        assert not miss, (frac, miss[:3])
        meta = meta.reindex(cells)
        for (pt, cond), idx in meta.reset_index(drop=True).groupby(["patient", "condition"]).groups.items():
            typ = "tumour" if cond.lower().startswith("tumo") else ("adjacent_normal" if "peri" in cond.lower() else None)
            if typ is None:
                continue
            idx = np.asarray(idx)
            parts.append((X[idx], genes, [cells[i] for i in idx], dict(sample=f"{pt}_{cond}_{frac}", patient=pt, sample_type=typ,
                                                                    chemistry="10x 3'", fraction=frac)))
    X, genes, obs = assemble(parts)
    obs["author_cell_type"], obs["sorting"], obs["treatment"] = "", "FACS fraction (EC-enriched or stroma; CD45-depleted)", "unknown"
    k, dropped = keep_paired(obs)
    X, obs = X[np.flatnonzero(k)], obs[k]
    notes = ["Two fractions per condition (EC-enriched and stroma) from the deposited raw matrices; condition tumor -> tumour, 'peritumoral' / 'enriched pEC' -> adjacent_normal.",
             f"Unpaired patients dropped: {dropped or 'none'}. treatment not stated -> unknown."]
    return B.finalize(gse, X, genes, obs, B._prov(recs), notes, title="Breast cancer tumour and peritumoural endothelial / stromal cells")


def census(pre, patient, sample, sample_type, chem, treatment, extra=()):
    """Census source h5ad (raw counts in `raw.X`, symbols in raw.var.feature_name)."""
    d = A226 / f"cellxgene_{pre}"
    a = ad.read_h5ad(d / f"{pre}.h5ad")
    X = sp.csr_matrix(a.raw.X)
    genes = a.raw.var["feature_name"].astype(str).tolist()
    obs = pd.DataFrame({"sample": a.obs[sample].astype(str).to_numpy(), "patient": a.obs[patient].astype(str).to_numpy(),
                        "sample_type": sample_type(a.obs), "chemistry": chem}, index=a.obs.index.astype(str))
    for c in extra:
        obs[c] = a.obs[c].astype(str).to_numpy()
    obs["author_cell_type"] = a.obs.get("cell_type", pd.Series("", index=a.obs.index)).astype(str).to_numpy()
    obs["sorting"], obs["treatment"] = "", treatment
    obs = obs[obs.sample_type.isin(["tumour", "adjacent_normal"])]
    X = X[a.obs.index.get_indexer(obs.index)]
    k, dropped = keep_paired(obs)
    return X[np.flatnonzero(k)], genes, obs[k], dropped, d


def cxg_b6b5ea88():
    pre = "b6b5ea88"
    X, genes, obs, dropped, d = census(pre, "donor_id", "sample_id",
                                       lambda o: np.where(o["source"].astype(str) == "tumour", "tumour", "adjacent_normal"), "10x 3' v3", "unknown",
                                       extra=("hpv_stat",))
    # tumour libraries (TC, TL, T) and adjacent libraries of the same donor; counts as deposited
    notes = ["CELLxGENE Census dataset b6b5ea88 (OPSCC, doi 10.1186/s12943-024-02191-9): `source` tumour / normal, donor_id, sample_id; raw counts from raw.X; "
             "treatment not stated -> unknown; author cell labels (Labels_broad) not used, marker rules as A149.", f"Unpaired donors dropped: {dropped or 'none'}."]
    prov = dict(files=[dict(file=f"{pre}.h5ad", url="CELLxGENE Census 2025-11-08 download_source_h5ad (see PROVENANCE.md of the download)", bytes=(d / f"{pre}.h5ad").stat().st_size, sha256=B.sha256(d / f"{pre}.h5ad"))], filtering="None.")
    out = A226 / "cxg_b6b5ea88"
    B.ROOT = A226
    (A226 / "OPSCC_b6b5ea88").mkdir(exist_ok=True)
    return B.finalize("OPSCC_b6b5ea88", X, genes, obs, prov, notes, title="OPSCC tumour and normal oropharynx (CELLxGENE)")


def cxg_5af90777():
    pre = "5af90777"
    def st(o):
        s = o["summaryDescription"].astype(str)
        return np.where(s == "Tumour", "tumour", np.where(s == "Normal kidney", "adjacent_normal", "other"))
    X, genes, obs, dropped, d = census(pre, "donor_id", "orig.ident", st, "10x 5'", "unknown", extra=("summaryDescription", "region"))
    notes = ["CELLxGENE Census dataset 5af90777 (Li 2022 Cancer Cell, kidney tumour territories): tumour = summaryDescription 'Tumour' (RCC, oncocytoma "
             "and benign tumour libraries included as tumour libraries of the donor, disease label per cell), adjacent_normal = 'Normal kidney'; "
             "'Tumour-normal', blood, adrenal, fat, thrombus and metastasis libraries are not used.", f"Unpaired donors dropped: {dropped or 'none'}."]
    prov = dict(files=[dict(file=f"{pre}.h5ad", url="CELLxGENE Census 2025-11-08 download_source_h5ad", bytes=(d / f"{pre}.h5ad").stat().st_size, sha256=B.sha256(d / f"{pre}.h5ad"))], filtering="None.")
    (A226 / "RCC_5af90777").mkdir(exist_ok=True)
    return B.finalize("RCC_5af90777", X, genes, obs, prov, notes, title="Kidney tumour and normal kidney (CELLxGENE, Li 2022)")


def _qian(study, prefix, tissue):
    """Qian et al. 2020 (Cell Res) donors of the CELLxGENE multi-tissue TME atlas b617ee1b (donor_id `Qian<tissue>_...`): disease
    'normal' = normal tissue of the donor, otherwise tumour; CD45+-sorted fractions (harm_cd45pos == yes) flagged `CD45+`."""
    pre = "b617ee1b"
    d = A226 / f"cellxgene_{pre}"
    a = ad.read_h5ad(d / f"{pre}.h5ad", backed="r")
    sel = np.flatnonzero(a.obs.donor_id.astype(str).str.startswith(prefix).to_numpy())
    o = a.obs.iloc[sel]
    X = sp.csr_matrix(a.raw.X[sel])
    genes = a.raw.var["feature_name"].astype(str).tolist()
    obs = pd.DataFrame({"sample": (o.donor_id.astype(str) + "|" + o["harm_sample.type"].astype(str) + "|" + o["harm_cd45pos"].astype(str)).to_numpy(),
                        "patient": o.donor_id.astype(str).to_numpy(),
                        "sample_type": np.where(o.disease.astype(str) == "normal", "adjacent_normal", "tumour"),
                        "chemistry": o.assay.astype(str).to_numpy(), "author_cell_type": o.cell_type.astype(str).to_numpy(),
                        "sorting": np.where(o["harm_cd45pos"].astype(str) == "yes", "CD45+", ""), "treatment": "unknown"}, index=o.index.astype(str))
    k, dropped = keep_paired(obs)
    X, obs = X[np.flatnonzero(k)], obs[k]
    notes = [f"Donors `{prefix}*` of CELLxGENE Census dataset b617ee1b (Qian et al. 2020, Cell Res; harmonised multi-tissue TME atlas): disease normal -> adjacent_normal, "
             "else tumour; raw counts from raw.X; sorting CD45+ where harm_cd45pos is yes. Samples are donor x sample type x CD45 fraction (the per-sample "
             "library is not in the Census object). Treatment not stated in the atlas -> unknown.", f"Unpaired donors dropped: {dropped or 'none'}."]
    prov = dict(files=[dict(file=f"{pre}.h5ad", url="CELLxGENE Census 2025-11-08 download_source_h5ad", bytes=(d / f"{pre}.h5ad").stat().st_size, sha256=B.sha256(d / f"{pre}.h5ad"))], filtering="Donor subset.")
    (A226 / study).mkdir(exist_ok=True)
    return B.finalize(study, X, genes, obs, prov, notes, title=f"Qian 2020 {tissue}: tumour and normal tissue of the same donors (CELLxGENE atlas)")


def Qian2020_CRC():
    return _qian("Qian2020_CRC", "QianCRC", "colorectal")


def Qian2020_LUAD():
    return _qian("Qian2020_LUAD", "QianLUAD", "lung")


if __name__ == "__main__":
    for s in sys.argv[1:]:
        globals()[s]()
