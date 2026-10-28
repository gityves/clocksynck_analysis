"""A149 part 2: GEO processed count matrices -> one h5ad per cohort, in the format of the legacy GEO builder.

Usage:  python build_geo_a149_part2.py GSE206785 GSE150290 ...

Helpers (fetch, fetch_all, soft_samples, assert_integer, finalize, write_provenance, stack_samples, _tar_10x,
_mm_to_cells_by_genes) are imported from the legacy builder (read-only), so output layout, obs columns and PROVENANCE.md are
identical: <ROOT>/<GSE>/<GSE>.h5ad + PROVENANCE.md. X = raw integer UMI counts (CSR int32, cells x genes, symbols).
obs: sample, patient, sample_type (tumour | adjacent_normal | metastasis | other), author_cell_type, sorting ("" unless sorted),
treatment (legacy vocabulary naive | treated | unknown; finalize() asserts it), total_umi, chemistry.

A cohort whose supplied matrix fails the integer check is DROPPED (a DROPPED.md is written in its folder), never reconstructed.
Downloads are only parsed (gzip/tar/mtx/csv readers); nothing downloaded is executed; .rds files are never read.
"""
import gzip
import io
import itertools
import re
import sys
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

sys.path.insert(0, "/home/maxine/Documents/paychere/scCircadianCancer/legacy/scripts/geo")
import build_geo as bg  # noqa: E402
from build_geo import (ROOT, _mm_to_cells_by_genes, _prov, _tar_10x, assert_integer, fetch, fetch_all,  # noqa: E402,F401
                       finalize, sha256, soft_samples, stack_samples)


def _pair_stats(obs):
    """(patients with tumour and adjacent_normal, patients); replaces the legacy version, whose groupby.agg(set) fails in current pandas."""
    pts = {p: set(g) for p, g in obs.groupby("patient", observed=True).sample_type}
    return int(sum(("tumour" in s) and ("adjacent_normal" in s) for s in pts.values())), len(pts)


bg.pair_stats = _pair_stats


class Dropped(Exception):
    """Supplied matrix is not raw integer counts."""


# ----------------------------------------------------------------------------- local helpers
def read_dense(src, sep="\t", skip_genes=(), label=""):
    """Streaming reader for a dense genes x cells text matrix (path or binary file object, gz). One gene row at a time.
    Same as legacy read_dense_gene_by_cell, plus `skip_genes` (annotation rows embedded in the matrix) and file objects.
    Header may or may not carry a leading empty label. Returns CSR cells x genes (int32), genes, cells; asserts integers."""
    if isinstance(src, (str, Path)):
        f = gzip.open(src, "rt") if str(src).endswith(".gz") else open(src)
    else:
        f = io.TextIOWrapper(gzip.GzipFile(fileobj=src), errors="replace")
    genes, indptr, buf_i, buf_d, nnz = [], [0], [], [], 0
    with f:
        header = f.readline().rstrip("\n").replace('"', "").split(sep)
        first = f.readline()
        ncell = first.count(sep)
        cells = header[1:] if len(header) == ncell + 1 else header
        assert len(cells) == ncell, (len(cells), ncell)
        for k, line in enumerate(itertools.chain([first], f), 1):
            line = line.rstrip("\n")
            if not line:
                continue
            i = line.index(sep)
            g = line[:i].strip('"')
            if g in skip_genes:
                continue
            v = np.fromstring(line[i + 1:], dtype=np.float32, sep=sep)
            assert v.size == ncell, (g, v.size, ncell)
            nz = np.flatnonzero(v)
            d = v[nz]
            if not (np.all(d == np.round(d)) and d.min(initial=0) >= 0):
                raise Dropped(f"non-integer/negative value in row {g} (e.g. {d[d != np.round(d)][:3].tolist()})")
            genes.append(g)
            buf_i.append(nz.astype(np.int32))
            buf_d.append(d.astype(np.int32))
            nnz += nz.size
            indptr.append(nnz)
            if k % 4000 == 0:
                print(f"  {label} read {k} genes", flush=True)
    M = sp.csr_matrix((np.concatenate(buf_d), np.concatenate(buf_i), np.asarray(indptr, dtype=np.int64)),
                      shape=(len(genes), ncell))
    return M.T.tocsr(), genes, cells


def probe_dense_integer(path, sep, nrows=400, skip=0):
    """Read the first `nrows` gene rows of a dense gz matrix (header + rows) and return the fraction of non-integer values."""
    with gzip.open(path, "rt") as f:
        f.readline()
        vals = []
        for k, line in enumerate(f):
            if k < skip:
                continue
            if k >= nrows + skip:
                break
            vals.append(np.fromstring(line.split(sep, 1)[1], dtype=np.float64, sep=sep))
    v = np.concatenate(vals)
    nz = v[v != 0]
    return float(np.mean(nz != np.round(nz))) if nz.size else 0.0, nz[:5].tolist()


def soft_field(soft, gsm, key):
    """Value of a '!Sample_<key>' line of a GSM from soft_samples()['other'] (or title)."""
    for l in soft[gsm].get("other", []):
        if l.startswith(f"!Sample_{key} "):
            return l.split("=", 1)[1].strip()
    return ""


def tar_extract(gse, name):
    """Copy one tar member of <gse>_RAW.tar to <ROOT>/<gse>/_members/<basename> (parse-only data; basename only)."""
    dest = ROOT / gse / "_members" / Path(name).name
    if not dest.exists():
        dest.parent.mkdir(exist_ok=True)
        with tarfile.open(ROOT / gse / f"{gse}_RAW.tar") as tf, open(dest, "wb") as out:
            out.write(tf.extractfile(name).read())
    return dest


def raw_tar_members(gse):
    with tarfile.open(ROOT / gse / f"{gse}_RAW.tar") as tf:
        return sorted(m.name for m in tf.getmembers() if m.isfile())


def _obs(cells, **cols):
    o = pd.DataFrame(index=[str(c) for c in cells])
    for k, v in cols.items():
        o[k] = v
    return o


def _common(obs, **defaults):
    for k, v in dict(author_cell_type="", sorting="", treatment="unknown", chemistry="").items():
        if k not in obs:
            obs[k] = defaults.get(k, v)
    return obs


def _dropped(gse, reason, recs=None):
    (ROOT / gse / "DROPPED.md").write_text(f"# {gse} dropped (A149 part 2)\n\n{reason}\n")
    raise Dropped(reason)


# ============================================================================ cohorts
def GSE206785():
    """Gastric, 24 patients: scgex.txt.gz (cells x genes, comma). Suspected normalised -> integer check only."""
    gse = "GSE206785"
    fetch_all(gse, [f"{gse}_scgex.txt.gz", f"{gse}_metadata.txt.gz"])
    # file is cells x genes: rows are cells, the first column the cell id
    with gzip.open(ROOT / gse / f"{gse}_scgex.txt.gz", "rt") as f:
        f.readline()
        vals = np.concatenate([np.fromstring(f.readline().split(",", 1)[1], dtype=np.float64, sep=",") for _ in range(300)])
    nz = vals[vals != 0]
    frac = float(np.mean(nz != np.round(nz)))
    if frac > 0:
        _dropped(gse, f"scgex.txt.gz fails the integer check: {frac:.1%} of nonzero values in the first 300 cells are non-integer "
                      f"(smallest nonzero values {np.unique(nz)[:3].tolist()}; ln 2 = 0.6931 suggests log1p of counts). "
                      "GEO text: SCTransform-normalised. Counts not reconstructed.")
    raise RuntimeError("GSE206785 unexpectedly integer: implement the reader")


def GSE290298():
    """HCC N1-N4 / T1-T4: single All_processed_data.csv.gz (genes x cells, comma). Integer check on the first 400 genes."""
    gse = "GSE290298"
    fetch_all(gse, [f"{gse}_All_processed_data.csv.gz"])
    frac, ex = probe_dense_integer(ROOT / gse / f"{gse}_All_processed_data.csv.gz", ",", nrows=400, skip=0)
    if frac > 0:
        _dropped(gse, f"All_processed_data.csv.gz fails the integer check: {frac:.1%} of nonzero values in the first 400 gene rows "
                      f"are non-integer (e.g. {ex}); the file is normalised/log-scaled. Counts not reconstructed.")
    raise RuntimeError("GSE290298 unexpectedly integer: implement the reader")


def GSE150290():
    """Gastric (intestinal/diffuse): per-sample genes x cells UMI tables `<GSM>_PatNN-{A,B}.txt.gz` in RAW.tar (authors' QC'd cells).
    The raw_gene_bc_matrices tarballs (unfiltered droplets) are not used."""
    gse = "GSE150290"
    recs = fetch_all(gse, [f"{gse}_RAW.tar"])
    soft = soft_samples(gse)
    names = [n for n in raw_tar_members(gse) if re.fullmatch(r"GSM\d+_Pat\d+-[AB]\.txt\.gz", n)]
    mats, gl, obs_parts = [], [], []
    for n in names:
        gsm, pat, ab = re.fullmatch(r"(GSM\d+)_(Pat\d+)-([AB])\.txt\.gz", n).groups()
        src = soft_field(soft, gsm, "source_name_ch1")
        assert soft[gsm]["title"] == f"{pat}-{ab}", (gsm, soft[gsm]["title"])
        assert src in {"Adjacent normal cells", "Cancer cells", "Control"}, src
        if src == "Control":
            st = "other"
            assert int(pat[3:]) >= 25
        else:
            st = "adjacent_normal" if src == "Adjacent normal cells" else "tumour"
            assert (ab == "A") == (st == "adjacent_normal") and int(pat[3:]) <= 24, (n, src)
        X, genes, cells = read_dense(tar_extract(gse, n), label=pat + ab)
        mats.append(X)
        gl.append(genes)
        o = _obs(cells, sample=f"{pat}-{ab}", patient=pat, sample_type=st, gsm=gsm, source_name=src)
        lau = soft[gsm]["chars"].get("lauren classification", "")
        o["lauren"] = lau
        obs_parts.append(o)
        print(f"  {pat}-{ab} {X.shape[0]} cells x {len(genes)} genes ({st})", flush=True)
    X, genes, _ = stack_samples(mats, gl)
    obs = _common(pd.concat(obs_parts), chemistry="10x 3' v2")
    notes = [
        "Files: per-sample tab-delimited UMI tables (genes x cells) `<GSM>_PatNN-{A,B}.txt.gz` from RAW.tar: the authors' quality-controlled cells "
        "(series text: 13,022 cells; thresholds nUMI 100-20000, mito/Hb < 10%). The `raw_gene_bc_matrices.tar.gz` members (unfiltered "
        "droplets, 112k cells in the paper) are not used. Gene lists differ per sample (zero-count genes dropped by the authors) and are unioned.",
        "patient: Pat<NN> from the GSM title `PatNN-A/-B` (asserted equal to file name).",
        "sample_type: GSM source_name_ch1 `Adjacent normal cells` (-A) -> adjacent_normal, `Cancer cells` (-B) -> tumour, `Control` "
        "(Pat25-Pat29, non-cancer controls, A only) -> other; asserted consistent with the A/B suffix. Pat21 has an adjacent normal only.",
        "chemistry: 10x 3' v2 (GSM extract protocol `10x signle cell 3' v2 kit`). treatment: unknown (not stated). sorting: '' (none stated).",
        "author_cell_type: none deposited -> ''. `lauren` column = GSM `lauren classification` (Intestinal/Diffuse).",
    ]
    return finalize(gse, X, genes, obs, _prov(recs, filtering="Authors' QC'd cells (per-sample UMI tables); raw droplet matrices not used."),
                    notes, title="Gastric cancer, adjacent normal and non-cancer controls")


def GSE167297():
    """Diffuse-type gastric: per-sample `<GSM>_PtN_{Normal,Superficial,Deep}_CountMatrix.txt.gz` (genes x cells)."""
    gse = "GSE167297"
    recs = fetch_all(gse, [f"{gse}_RAW.tar"])
    soft = soft_samples(gse)
    mats, gl, obs_parts = [], [], []
    for n in raw_tar_members(gse):
        gsm, pt, layer = re.fullmatch(r"(GSM\d+)_(Pt\d+)_(Normal|Superficial|Deep)_CountMatrix\.txt\.gz", n).groups()
        X, genes, cells = read_dense(tar_extract(gse, n), label=f"{pt}_{layer}")
        mats.append(X)
        gl.append(genes)
        o = _obs(cells, sample=f"{pt}_{layer}", patient=pt, sample_type="adjacent_normal" if layer == "Normal" else "tumour",
                 gsm=gsm, gsm_title=soft[gsm]["title"], tumour_layer="" if layer == "Normal" else layer.lower())
        obs_parts.append(o)
        print(f"  {pt}_{layer} {X.shape[0]} cells", flush=True)
    X, genes, _ = stack_samples(mats, gl)
    obs = pd.concat(obs_parts)
    obs.index = obs["sample"] + "_" + pd.Index(obs.index)  # barcodes repeat across samples
    obs = _common(obs, chemistry="10x 3' v3")
    notes = [
        "Per-sample dense count tables (genes x cells) from RAW.tar; barcodes are made unique by prefixing the sample id.",
        "patient: Pt<N> in the file name (GSM titles `PtN_normal` / `PtN_superficial_cancer` / `PtN_deep_cancer`). Pt2 has no normal.",
        "sample_type: Normal -> adjacent_normal (`paired normal gastric tissue`); Superficial and Deep layers of the cancer are both tumour "
        "(distinct sample ids `PtN_Superficial` / `PtN_Deep`; `tumour_layer` column keeps the layer).",
        "chemistry: 10x 3' v3 (Cell Ranger 3.1.0 pipeline; v3 per candidate sheet). treatment: unknown (not stated). sorting: ''. "
        "author_cell_type: none deposited -> ''.",
    ]
    return finalize(gse, X, genes, obs, _prov(recs), notes, title="Diffuse-type gastric cancer (superficial/deep layers) with matched normal")


def GSE199654():
    """ESCC epithelial cells only: UMI matrix (genes x cells), metadata with sampleType N/D/T."""
    gse = "GSE199654"
    pre = f"{gse}_scTDN"
    files = [f"{pre}_UMI_matrix_epithelial_cells.txt.gz", f"{pre}_METADATA_epithelial_cells.txt.gz", f"{pre}_barcode_epithelial_cells.txt.gz",
             f"{pre}_feature_epithelial_cells.txt.gz"]
    recs = fetch_all(gse, files)
    X, genes, cells = read_dense(ROOT / gse / files[0], label="ESCC")
    cells = [re.sub(r"\.(\d+)$", r"-\1", c) for c in cells]
    meta = pd.read_csv(ROOT / gse / files[1], sep="\t", index_col="cellID")
    bc = pd.read_csv(ROOT / gse / files[3 - 1], sep="\t", index_col=0)["x"].tolist()
    assert cells == bc == list(meta.index), "matrix header, barcode file and metadata disagree"
    m = meta.reindex(cells)
    obs = _obs(cells, sample=m.sampleID.values, patient=m.patientID.values,
               sample_type=m.sampleType.map({"N": "adjacent_normal", "T": "tumour", "D": "other"}).values)
    assert obs.sample_type.notna().all()
    obs["stage"] = m.stage.values
    obs = _common(obs, chemistry="10x 5' v1", treatment="naive")
    notes = [
        "Epithelial-cell matrix only (authors deposit epithelium only, for patient privacy). Matrix header uses `.1` barcode suffix; converted to `-1` "
        "to match the barcode and metadata files (all three asserted identical).",
        "patient: metadata `patientID` (LZE<n>). sample: metadata `sampleID` (LZE<n>N/D/T, with replicate numbers for some).",
        "sample_type: metadata `sampleType` N (non-tumour oesophagus, >= 5 cm from tumour) -> adjacent_normal, T -> tumour, D (low/high-grade "
        "intraepithelial neoplasia, <= 2 cm from tumour) -> other. `stage` = metadata stage column.",
        "chemistry: Chromium Single Cell 5' v1 (GSM extract protocol). treatment: naive (series: 14 treatment-naive ESCC patients). sorting: ''. "
        "author_cell_type: none (all epithelial) -> ''.",
        "Only the UMI-count matrix is used (the normalised file is ignored).",
    ]
    return finalize(gse, X, genes, obs, _prov(recs, filtering="Epithelial cells as deposited (no stroma/immune)."), notes,
                    title="ESCC epithelial cells (non-tumour / dysplasia / tumour)")


def GSE160269():
    """ESCC: CD45- and CD45+ dense UMI tables (genes x cells, space separated) with cell/sample/annotated_type files."""
    gse = "GSE160269"
    recs_f = []
    for pop in ("neg", "pos"):
        recs_f += [f"{gse}_CD45{pop}_UMIs.txt.gz", f"{gse}_CD45{pop}_cells.txt.gz"]
    recs = fetch_all(gse, recs_f)
    soft = soft_samples(gse)
    mats, gl, obs_parts = [], [], []
    for pop, lab in (("neg", "CD45-"), ("pos", "CD45+")):
        X, genes, cells = read_dense(ROOT / gse / f"{gse}_CD45{pop}_UMIs.txt.gz", sep=" ", label=f"CD45{pop}")
        ann = pd.read_csv(ROOT / gse / f"{gse}_CD45{pop}_cells.txt.gz", sep=" ", index_col="cell")
        assert list(ann.index) == cells, f"CD45{pop}: matrix header differs from cells file"
        pt = ann["sample"].str.extract(r"^(P\d+)([NT])$")
        assert pt.notna().all().all()
        o = _obs(cells, sample=(ann["sample"] + "_" + lab).values, patient=pt[0].values,
                 sample_type=pt[1].map({"N": "adjacent_normal", "T": "tumour"}).values,
                 sorting=lab, author_cell_type=ann.annotated_type.values)
        mats.append(X)
        gl.append(genes)
        obs_parts.append(o)
        print(f"  {lab}: {X.shape}", flush=True)
    X, genes, _ = stack_samples(mats, gl)
    obs = _common(pd.concat(obs_parts), chemistry="10x 5'")
    # cross-check sample/type against GSM info lines
    info = {}
    for d in soft.values():
        txt = next((r.split(":", 1)[1].strip() for r in d["raw"] if r.startswith("info:")), "")
        m = re.match(r"10X 5' scRNA-Seq for (CD45[+-]) cells from (P\d+) (TUMOR|Normal Tissue)", txt)
        if m:
            info.setdefault((m.group(2), m.group(1)), set()).add(m.group(3))
    tn = obs.groupby(["patient", "sorting"]).sample_type.agg(lambda s: set(s))
    for (p, so), s in tn.items():
        want = {"adjacent_normal" if x == "Normal Tissue" else "tumour" for x in info.get((p, so), set())}
        assert want and s <= want, (p, so, s, want)
    notes = [
        "Two matrices: `CD45neg_UMIs` (CD45- sorted, 5' 10x) and `CD45pos_UMIs` (CD45+ sorted), genes x cells, space separated; both kept, "
        "`sorting` = CD45- / CD45+. Cell barcodes carry the sample/sort prefix (`P1T-E-...` / `-I-`), unique across both files. `sample` = "
        "`<sample>_CD45-` / `_CD45+`.",
        "patient: sample prefix `P<n>` of the authors' cells file (`P126N`, `P126T`...). sample_type: suffix N (`Tumor Normal Adjacent tissue removed "
        "from esophageal carcinoma`) -> adjacent_normal, T -> tumour; cross-checked against the GSM info lines (`... from P126 Normal Tissue / TUMOR`). "
        "Only P126, P127, P128, P130 have a normal (4 pairs); 60 tumours in all.",
        "author_cell_type: authors' `annotated_type` from the cells files. chemistry: 10x 5' (GSM info lines). treatment: unknown (not stated). "
        "Raw data: GSA HRA000195.",
    ]
    return finalize(gse, X, genes, obs, _prov(recs), notes, title="ESCC tumours and adjacent normals, CD45- and CD45+ fractions")


def GSE178481():
    """ccRCC: per-sample `<GSM>_RCC-<pt>-<Normal|PTumor[k]>.count.csv.gz` (genes x cells, comma)."""
    gse = "GSE178481"
    recs = fetch_all(gse, [f"{gse}_RAW.tar"])
    soft = soft_samples(gse)
    mats, gl, obs_parts = [], [], []
    for n in raw_tar_members(gse):
        gsm, pt, kind = re.fullmatch(r"(GSM\d+)_RCC-([A-Z]+\d+)-(Normal|PTumor\d*)\.count\.csv\.gz", n).groups()
        ch = soft[gsm]["chars"]
        assert soft[gsm]["title"] == f"RCC-{pt}-{kind}"
        assert (ch["condition"] == "Normal") == (kind == "Normal")
        thr = "thrombus" in ch["treatment"].lower()
        st = "adjacent_normal" if kind == "Normal" else ("other" if thr else "tumour")
        tr = "treated" if "chemotherapy" in ch["treatment"].lower() else "naive"
        X, genes, cells = read_dense(tar_extract(gse, n), sep=",", label=f"{pt}-{kind}")
        mats.append(X)
        gl.append(genes)
        o = _obs(cells, sample=f"{pt}-{kind}", patient=pt, sample_type=st, treatment=tr, gsm=gsm, histology=ch["Stage"],
                 metastatic_sites=ch["metastasic sites"], sorting="live CD235a-")
        obs_parts.append(o)
        print(f"  {pt}-{kind} {X.shape[0]} cells ({st}, {tr})", flush=True)
    X, genes, _ = stack_samples(mats, gl)
    obs = _common(pd.concat(obs_parts), chemistry="10x 3' v2")
    notes = [
        "Per-sample dense count tables (genes x cells, csv) from RAW.tar; cell ids already carry the sample name.",
        "patient: PR1-PR9 and BM1-BM2 from the GSM title `RCC-<patient>-Normal` / `-PTumor[k]`. sample_type: Normal -> adjacent_normal; PTumor, "
        "PTumor1-3 (multi-region) -> tumour. BM1/BM2 are primary kidney tumours of patients with metastatic disease (GSM source `kidney`, "
        "`metastasic sites` Bone / Bone, lung, adrenal gland, kidney): the matrix has NO brain-metastasis sample, so none is labelled metastasis. "
        "RCC-PR9-PTumor2 is the primary-tumour thrombus (GSM `treatment: Primary RCC tumor Thrombus`) -> other. PR7 has a tumour only, PR8 a normal only.",
        "treatment: GSM `treatment` field: `No` -> naive; the one `chemotherapy for collicular lymphoma` sample -> treated. "
        "sorting: `live CD235a-` (flow sorted for live non-erythroid cells, GSM extract protocol). chemistry: 10x 3' v2 (Cell Ranger 2.2).",
        "author_cell_type: none deposited -> ''. `histology` = GSM Stage field (ccRCC / pRCC).",
    ]
    return finalize(gse, X, genes, obs, _prov(recs), notes, title="Renal cell carcinoma, tumour regions and adjacent normal kidney")


def GSE242299():
    """ccRCC, inDrops-type (Atrandi/1CellBio): the authors' filtered all-cells h5ad (50,236 cells x 33,538 genes)."""
    gse = "GSE242299"
    fn = f"{gse}_all_cells_50236_33538.h5ad.gz"
    recs = fetch_all(gse, [fn])
    h5 = ROOT / gse / fn[:-3]
    if not h5.exists():
        with gzip.open(ROOT / gse / fn, "rb") as f, open(h5, "wb") as out:
            while True:
                b = f.read(1 << 24)
                if not b:
                    break
                out.write(b)
    import anndata as ad
    a = ad.read_h5ad(h5)
    o = a.obs
    d = a.X.data[:5_000_000]
    assert np.all(d == np.round(d)) and np.allclose(np.asarray(a.X.sum(1)).ravel(), o.total_counts.astype(float).to_numpy()), "X is not raw counts"
    assert set(o.tissue.astype(str)) == {"Tumor", "Healthy"}
    pt = o.patient.astype(str)
    assert (o["sample"].astype(str).str[1:] == pt.str[1:]).all()
    obs = _obs(o["library_name"].astype(str) + "_" + o.index.astype(str), sample=o["sample"].astype(str).values, patient=pt.values,
               sample_type=o.tissue.astype(str).map({"Tumor": "tumour", "Healthy": "adjacent_normal"}).values,
               author_cell_type=o.cell_type.astype(str).values)
    obs["emulsion"] = o.library_name.astype(str).values
    obs["author_broad_cell_type"] = o.broad_cell_type.astype(str).values
    obs["beads"] = o.beads.astype(str).values
    obs["pt_stage"] = o["pT stage"].astype(str).values
    obs = _common(obs, chemistry="inDrops-type (Atrandi/1CellBio)", treatment="naive")
    notes = [
        "Authors' filtered all-cells h5ad (50,236 cells x 33,538 genes, raw counts in .X: integer-valued and row sums equal the authors' `total_counts`, both asserted); "
        "decompressed next to the .gz and read with anndata. The per-emulsion RAW.tar MEX files hold 147,456 unfiltered barcodes per emulsion "
        "and two are damaged (T2_1 matrix truncated, N3_1 barcodes/features empty), so they are not used. Cell ids = emulsion name + the h5ad row label.",
        "patient: authors' `patient` (P1-P9; = GSM title `patient N ...`). sample = authors' `sample` (T<n> tumour, N<n> healthy-adjacent kidney): the "
        "emulsions (`library_name`, kept as `emulsion`) of one tissue are merged per sample. Patient 1 has a normal only.",
        "sample_type: `tissue` Tumor -> tumour, Healthy (GSM `healthy-adjacent kidney`) -> adjacent_normal. author_cell_type: authors' `cell_type` "
        "(`broad_cell_type` kept). chemistry: inDrops-type droplet barcoding (Atrandi/1CellBio beads, STARsolo), not 10x. "
        "treatment: naive (series: no patient had received prior systemic therapy). sorting: ''.",
    ]
    return finalize(gse, a.X, list(a.var_names), obs, _prov(recs, filtering="Authors' filtered all-cells h5ad; per-emulsion raw MEX files not used."),
                    notes, title="ccRCC and healthy-adjacent kidney (inDrops-type)")


def GSE282701():
    """HBV+/HBV- HCC, SeekOne droplet: per-sample 10x-style MEX (tumor T / para-cancerous P) in RAW.tar."""
    gse = "GSE282701"
    recs = fetch_all(gse, [f"{gse}_RAW.tar"])
    soft = soft_samples(gse)
    tf = tarfile.open(ROOT / gse / f"{gse}_RAW.tar")
    names = sorted(m.name for m in tf.getmembers())
    mats, gl, il, obs_parts = [], [], [], []
    for gsm, d in sorted(soft.items()):
        pt, kind, hbv = d["title"].split(",")
        assert kind in ("tumor", "adjacent")
        stem = [n for n in names if n.startswith(gsm + "_") and n.endswith("_matrix.mtx.gz")][0][:-len("_matrix.mtx.gz")]
        sid = stem.split("_", 1)[1]
        assert sid.endswith("T" if kind == "tumor" else "P"), (sid, kind)
        M, genes, gids, bc = _tar_10x(tf, names, stem)
        mats.append(M)
        gl.append(genes)
        il.append(gids)
        o = _obs(bc, sample=f"{pt}_{kind}", patient=pt, sample_type="tumour" if kind == "tumor" else "adjacent_normal", gsm=gsm,
                 gsm_title=d["title"], hbv_status="HBV+" if hbv == "HBV-po" else "HBV-")
        o.index = o["sample"] + "_" + pd.Index(o.index)
        obs_parts.append(o)
        print(f"  {pt}_{kind} {M.shape[0]} barcodes x {M.shape[1]} genes", flush=True)
    X, genes, gids = stack_samples(mats, gl, il)
    obs = _common(pd.concat(obs_parts), chemistry="SeekOne DD 3'", treatment="naive")
    notes = [
        "Per-sample MEX (barcodes/features/matrix) from RAW.tar as deposited (cell-called matrices; barcode counts per sample in the table below).",
        "patient: GSM title `P<n>,<tumor|adjacent>,<HBV-po|HBV-ne>` -> P<n> (P10-P15). sample_type: tumor -> tumour; adjacent (para-cancerous "
        "tissue, file suffix P) -> adjacent_normal (asserted against the T/P suffix of the file stem). `hbv_status` from the title.",
        "chemistry: SeekOne Digital Droplet Single Cell 3' (SeekGene). treatment: naive (series: no anti-tumour treatment before surgery). "
        "sorting: ''. author_cell_type: none -> ''. Adjacent liver is para-cancerous and may be cirrhotic/hepatitic.",
    ]
    return finalize(gse, X, genes, obs, _prov(recs), notes, gene_ids=gids, title="HBV+/- hepatocellular carcinoma and adjacent liver (SeekOne)")


def GSE242889():
    """HCC (MVI), Singleron CeleScope: RAW.tar of per-sample plain tar archives (matrix.mtx, barcodes.tsv, genes.tsv)."""
    gse = "GSE242889"
    recs = fetch_all(gse, [f"{gse}_RAW.tar"])
    soft = soft_samples(gse)
    mats, gl, il, obs_parts = [], [], [], []
    with tarfile.open(ROOT / gse / f"{gse}_RAW.tar") as tf:
        for m in sorted(tf.getmembers(), key=lambda m: m.name):
            gsm, tag, ident = re.fullmatch(r"(GSM\d+)_(\dN?T)_([PC]\d+)\.tar\.gz", m.name).groups()
            title = soft[gsm]["title"]
            assert title == tag, (title, tag)
            num = tag[0]
            nested = tarfile.open(fileobj=tf.extractfile(m), mode="r:")  # plain tar despite the .tar.gz name
            mem = {Path(x.name).name: x for x in nested.getmembers() if x.isfile()}
            assert set(mem) == {"matrix.mtx", "barcodes.tsv", "genes.tsv"}, mem.keys()
            M = _mm_to_cells_by_genes(nested.extractfile(mem["matrix.mtx"]))
            bc = pd.read_csv(nested.extractfile(mem["barcodes.tsv"]), header=None, dtype=str)[0].tolist()
            ft = pd.read_csv(nested.extractfile(mem["genes.tsv"]), header=None, sep="\t", dtype=str)
            assert M.shape == (len(bc), len(ft)), (m.name, M.shape, len(bc), len(ft))
            mats.append(M)
            gl.append(ft[1].tolist())
            il.append(ft[0].tolist())
            o = _obs(bc, sample=f"{num}{'NT' if 'N' in tag else 'T'}", patient=f"HCC{num}",
                     sample_type="adjacent_normal" if "N" in tag else "tumour", gsm=gsm, gsm_title=title, library_id=ident,
                     mvi=soft_field(soft, gsm, "source_name_ch1"))
            o.index = o["sample"] + "_" + pd.Index(o.index)
            obs_parts.append(o)
            print(f"  {tag} {M.shape[0]} cells x {M.shape[1]} genes", flush=True)
    X, genes, gids = stack_samples(mats, gl, il)
    obs = _common(pd.concat(obs_parts), chemistry="Singleron CeleScope", treatment="naive")
    notes = [
        "RAW.tar members are plain tar archives named `.tar.gz` (not compressed) holding `matrix.mtx` (genes x cells), `barcodes.tsv`, `genes.tsv` "
        "(Ensembl id, symbol); read member by member, no extraction to disk. Cells as deposited.",
        "patient: leading number of the GSM title (`1NT`, `1T` -> HCC1); the library ids differ in letter only (`1NT_P21` normal, `1T_C21` cancer, "
        "same number 21), kept as `library_id`. sample_type: `NT` (GSM source `adjacent tissues`) -> adjacent_normal, `T` -> tumour; `mvi` = GSM source name.",
        "chemistry: Singleron Matrix / CeleScope (STAR, featureCounts, GRCh38 Ensembl 92). treatment: naive (series: no prior therapeutic "
        "intervention). sorting: ''. author_cell_type: none -> ''.",
    ]
    return finalize(gse, X, genes, obs, _prov(recs), notes, gene_ids=gids, title="HCC (MVI) tumour and adjacent tissue, Singleron")


def GSE144240():
    """Cutaneous SCC, P1-P10 tumour and patient-matched normal skin. SuperSeries: the scRNA matrix is in sub-series GSE144236
    (`GSE144236_cSCC_counts.txt.gz`, genes x cells, with two annotation rows) + `patient_metadata_new.txt.gz`."""
    gse = "GSE144240"
    sub = "GSE144236"
    base = bg._base(sub)
    import subprocess
    recs = []
    for fname, kind in [(f"{sub}_cSCC_counts.txt.gz", "suppl"), (f"{sub}_patient_metadata_new.txt.gz", "suppl"), (f"{sub}_family.soft.gz", "soft")]:
        dest = ROOT / gse / fname
        url = f"{base}/{kind}/{fname}"
        if not dest.exists():
            subprocess.run(["curl", "-sS", "-f", "-C", "-", "--retry", "5", "-o", str(dest), url], check=True)
        recs.append(dict(file=fname, url=url, bytes=dest.stat().st_size, sha256=sha256(dest)))
    f0 = fetch_all(gse, [], soft=True)  # SuperSeries SOFT (GSE144240_family.soft.gz)
    recs += f0
    X, genes, cells = read_dense(ROOT / gse / f"{sub}_cSCC_counts.txt.gz", skip_genes={"Patient", "Tissue: 0=Normal, 1=Tumor"}, label="cSCC")
    meta = pd.read_csv(ROOT / gse / f"{sub}_patient_metadata_new.txt.gz", sep="\t", index_col=0)
    m = meta.reindex(cells)
    assert m.patient.notna().all() and len(meta) == len(cells)
    pre = pd.Series(cells).str.extract(r"^(P\d+)_(Tumor|Normal)_")
    assert (pre[0].values == m.patient.values).all() and (pre[1].values == m["tum.norm"].values).all()
    obs = _obs(cells, sample=(pre[0] + "_" + pre[1]).values, patient=pre[0].values,
               sample_type=pre[1].map({"Tumor": "tumour", "Normal": "adjacent_normal"}).values,
               author_cell_type=m.level2_celltype.fillna("").astype(str).values)
    obs["author_cell_type_l1"] = m.level1_celltype.values
    obs["author_cell_type_l3"] = m.level3_celltype.values
    obs["author_ncount"] = m.nCount_RNA.values
    obs = _common(obs, chemistry="10x 3' v2")
    notes = [
        "GSE144240 is a SuperSeries; its RAW.tar holds spatial (ST) and WES files only. The scRNA-seq matrix is `GSE144236_cSCC_counts.txt.gz` of "
        "sub-series GSE144236 (48,164 cells, hg19 symbols; the two annotation rows `Patient` and `Tissue: 0=Normal, 1=Tumor` at the top of the "
        "matrix are skipped) with `GSE144236_patient_metadata_new.txt.gz` (author cell types). Data stored in the GSE144240 folder.",
        "patient: barcode prefix P<n> (asserted equal to metadata `patient`). sample_type: `Tumor` -> tumour, `Normal` -> adjacent_normal (GSM `tissue: Normal Skin`; "
        "patient-matched normal skin adjacent to the cSCC). Replicate libraries (e.g. P1 scRNA_1/_2) cannot be separated and form one sample per patient and tissue.",
        "author_cell_type: metadata `level2_celltype` (level 1 and 3 kept as extra columns). chemistry: 10x 3' v2 (series design). "
        "treatment: unknown (primary cSCC, not stated). sorting: ''.",
    ]
    return finalize(gse, X, genes, obs, _prov(recs), notes, title="Cutaneous SCC and patient-matched normal skin (scRNA, via GSE144236)")


COHORTS = {k: v for k, v in list(globals().items()) if re.fullmatch(r"GSE\d+", k)}

if __name__ == "__main__":
    for g in sys.argv[1:]:
        print(f"=== {g}", flush=True)
        try:
            COHORTS[g]()
        except Dropped as e:
            print(f"[{g}] DROPPED: {e}", flush=True)
            (ROOT / g / "DROPPED.md").write_text(f"# {g} dropped (A149 part 2)\n\n{e}\n")
