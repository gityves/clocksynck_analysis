"""A149 part 1: paired tumour / adjacent-normal single-cell cohorts (lung, prostate, bladder) built as h5ad.

Usage:  python build_geo_a149_part1.py GSE164789 [GSE229253 ...]

Same output format as the legacy GEO builder (imported unchanged for fetch / finalize / readers):
<ROOT>/<GSE>/<GSE>.h5ad (X = raw integer UMI, CSR cells x genes, symbols as var_names) and PROVENANCE.md.
obs: sample, patient, sample_type (tumour / adjacent_normal / other), author_cell_type, sorting ('' unless sorted),
treatment (naive / treated / unknown, the legacy vocabulary), total_umi, chemistry (+ cohort-specific columns).
Rules (A149 note 1): only patients with both a tumour and an adjacent-normal sample are kept (as the legacy paired builds);
a cohort failing `assert_integer` is dropped (finalize raises, nothing is written); TCR/BCR libraries are skipped.
Only GEO supplementary files and SOFT metadata are downloaded; they are parsed, never executed.
"""
import re
import sys
import tarfile
from pathlib import Path

sys.path.insert(0, "/home/maxine/Documents/paychere/scCircadianCancer/legacy/scripts/geo")
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scipy.sparse as sp  # noqa: E402
import build_geo as B  # noqa: E402  (legacy helpers, unchanged)
from build_geo import ROOT, _prov, _tar_10x, fetch_all, finalize, read_dense_gene_by_cell, soft_samples, stack_samples  # noqa: E402

ST = {"T": "tumour", "N": "adjacent_normal"}


def _pair_stats(obs):
    """Robust replacement of the legacy `pair_stats` (fails with the pandas / anndata of the sccc env, where obs columns are
    categorical after `write_h5ad`): (patients with tumour and adjacent_normal, patients)."""
    g = obs.assign(sample_type=obs.sample_type.astype(str), patient=obs.patient.astype(str)).groupby("patient").sample_type.agg(lambda x: set(x))
    return int(sum(("tumour" in s) and ("adjacent_normal" in s) for s in g)), int(len(g))


B.pair_stats = _pair_stats      # write_provenance looks pair_stats up in the legacy module namespace


# ----------------------------------------------------------------------------- helpers
def _stem(names, gsm):
    """Common stem of the per-GSM 10x files inside the tar (everything before [._]matrix.mtx.gz)."""
    m = [n for n in names if n.startswith(gsm + "_") and re.search(r"matrix\.mtx\.gz$", n)]
    assert len(m) == 1, (gsm, m)
    return re.sub(r"[._]matrix\.mtx\.gz$", "", m[0])


def _keep_paired(specs):
    """Keep samples of patients that have both tumour and adjacent_normal; return (kept, dropped patient list)."""
    by = {}
    for s in specs:
        by.setdefault(s["patient"], set()).add(s["sample_type"])
    ok = {p for p, t in by.items() if {"tumour", "adjacent_normal"} <= t}
    return [s for s in specs if s["patient"] in ok], sorted(set(by) - ok)


def call_cells(M, n_expected=5000, min_genes=200):
    """Cell calling for unfiltered (raw-droplet) matrices, label-free and per library (Cell Ranger v2-style ordmag rule): cutoff =
    UMI of the barcode at the 99th percentile of the `n_expected` top barcodes, divided by 10; also >= `min_genes` genes.
    Returns the boolean mask of called barcodes and the cutoff."""
    umi = np.asarray(M.sum(axis=1)).ravel()
    top = np.sort(umi)[::-1][:n_expected]
    cut = max(top[int(0.01 * len(top))] / 10.0, 1.0)
    genes = np.diff(M.indptr)
    return (umi >= cut) & (genes >= min_genes), float(cut)


def _assemble(parts, defaults):
    """parts: list of (M, genes, ids, barcodes, spec). Returns X, genes, gene_ids, obs."""
    mats, gl, il, ob = [], [], [], []
    for M, gn, ids, bc, s in parts:
        o = pd.DataFrame(index=[f"{s['sample']}_{b}" for b in bc])
        for k, v in s.items():
            o[k] = v
        mats.append(M); gl.append(gn); il.append(ids); ob.append(o)
        print(f"  {s['sample']} ({s['patient']}, {s['sample_type']}) {M.shape[0]} cells", flush=True)
    X, genes, gids = stack_samples(mats, gl, il if all(i is not None for i in il) else None)
    obs = pd.concat(ob)
    for k, v in defaults.items():
        if k not in obs:
            obs[k] = v
    return X, genes, gids, obs


def _read_mex_specs(gse, tf, names, specs, feat="features"):
    parts = []
    for s in specs:
        M, gn, ids, bc = _tar_10x(tf, names, _stem(names, s["gsm"]), feat=feat)
        parts.append((M, gn, ids, bc, s))
    return parts


def _dense_member(gse, tf, member):
    """Extract one dense genes x cells member to a temp file in the cohort folder (parsed only), read, delete."""
    tmp = ROOT / gse / "_tmp_member.txt.gz"
    with open(tmp, "wb") as f:
        f.write(tf.extractfile(member).read())
    try:
        X, genes, cells = read_dense_gene_by_cell(tmp)
    finally:
        tmp.unlink()
    return X, genes, cells


def _pair_note(kept_specs, dropped, what="patients"):
    return (f"Only {what} with both a tumour and an adjacent_normal sample are kept ({len(set(s['patient'] for s in kept_specs))}); "
            f"dropped (single type): {', '.join(dropped) if dropped else 'none'}.")


# ============================================================================ cohorts
def GSE164789():
    """Precursor and invasive LUAD (AIS/MIA/IAC): 3' v3 and 5' GEX libraries per GSM in RAW.tar; TCR libraries skipped."""
    gse = "GSE164789"
    recs = fetch_all(gse, [f"{gse}_RAW.tar"])
    soft = soft_samples(gse)
    tf = tarfile.open(ROOT / gse / f"{gse}_RAW.tar")
    names = sorted(m.name for m in tf.getmembers())
    specs = []
    for gsm, v in soft.items():
        desc = " ".join(x for x in v.get("other", []) if x.startswith("!Sample_description"))
        if "TCR" in desc or v["title"].endswith("-TCR"):
            continue
        m = re.fullmatch(r"([A-Z]+)-(K?)([AT])(\d?)", v["title"])
        assert m, v["title"]
        chem = "10x 5'" if "5'" in desc else "10x 3' v3"
        assert "GEX" in desc
        specs.append(dict(gsm=gsm, sample=v["title"], patient=m.group(1), sample_type=ST["T" if m.group(3) == "T" else "N"],
                          chemistry=chem, lesion=v["chars"]["histology"], sex=v["chars"].get("Sex", ""), age=v["chars"].get("age", "")))
    specs, dropped = _keep_paired(specs)
    X, genes, gids, obs = _assemble(_read_mex_specs(gse, tf, names, specs, feat="genes"), {})
    obs["author_cell_type"], obs["sorting"], obs["treatment"] = "", "", "naive"
    npat = len({s["patient"] for s in specs})
    notes = [
        _pair_note(specs, dropped) + " TCR libraries (16 GSMs, titles *-TCR) skipped: gene-expression libraries only.",
        "patient: GSM title prefix before '-' (e.g. `CXL-A` -> `CXL`; `LRH-KA`, `LRH-KT` -> `LRH`). LY (38 y F) and LYJ (60 y M) are different "
        f"patients by age/sex, so {npat} patients. Patients with two samples per type (HYY A1/A2, T1/T2) keep distinct sample ids.",
        "sample_type: title suffix A (adjacent, `-A`, `-A1`, `-KA`) -> adjacent_normal; T -> tumour (`-T`, `-T1`, `-KT`). All lesion stages are 'tumour'; "
        "`lesion` records the GSM characteristic `histology` (AIS / MIA / IAC), which is per patient and also given for the adjacent sample.",
        "chemistry: GSM description '10x Genomics 3' GEX library' -> `10x 3' v3` (extract protocol: Chromium Single Cell 3' V3, 46 samples); "
        "'5' GEX library' -> `10x 5'` (16 samples). Chemistry is mixed within the cohort, one chemistry per patient.",
        "treatment: naive (GSM treatment_protocol: 'treatment naive'). sorting: '' (no sorting stated). author_cell_type: none deposited.",
        "Cell Ranger 3.1.0 filtered matrices as deposited (GRCh38, genes.tsv, 2 columns); the authors state cells were filtered (>200 genes, mito < 20%).",
    ]
    return finalize(gse, X, genes, obs, _prov(recs, filtering="Cell Ranger matrices as deposited; GEX libraries of patients with both tumour and adjacent sample."),
                    notes, gene_ids=gids, title="Precursor and invasive lung adenocarcinoma, tumour and adjacent lung")


def GSE229253():
    """Early-stage LUAD snRNA-seq (STARsolo): per-GSM MEX in RAW.tar."""
    gse = "GSE229253"
    recs = fetch_all(gse, [f"{gse}_RAW.tar"])
    soft = soft_samples(gse)
    tf = tarfile.open(ROOT / gse / f"{gse}_RAW.tar")
    names = sorted(m.name for m in tf.getmembers())
    specs = []
    for gsm, v in soft.items():
        ch = v["chars"]
        t = {"Tumor": "T", "Lung adjacent normal": "N"}[ch["tissue"]]
        assert v["title"].startswith(ch["patient"] + "-" + t), v["title"]
        specs.append(dict(gsm=gsm, sample=v["title"].split()[0], patient=ch["patient"], sample_type=ST[t], chemistry="snRNA 10x 3'",
                          recurrence=ch.get("recurrence status", "")))
    specs, dropped = _keep_paired(specs)
    X, genes, gids, obs = _assemble(_read_mex_specs(gse, tf, names, specs), {})
    obs["author_cell_type"], obs["sorting"], obs["treatment"] = "", "", "unknown"
    notes = [
        _pair_note(specs, dropped),
        "patient: GSM characteristic `patient` (NYU<id>) = title prefix (asserted). sample_type: GSM `tissue` Tumor -> tumour, 'Lung adjacent normal' -> "
        "adjacent_normal; agrees with title suffix -T / -N (asserted).",
        "chemistry: single-nucleus 10x Genomics (Celsius Therapeutics, STARsolo); the 3' kit version is not stated in the record -> `snRNA 10x 3'`.",
        "treatment: unknown (stage I resections, prospectively banked 2005-2015; no statement on prior therapy). sorting: '' (nuclei, no sorting stated). "
        "`recurrence` = GSM `recurrence status`.",
        "Nuclei filtered by the authors (>500 genes, >1000 UMIs, <10% mito, scDblFinder doublets removed); counts as deposited (raw UMI, `counts_matrix`).",
    ]
    return finalize(gse, X, genes, obs, _prov(recs, filtering="STARsolo matrices as deposited (author-filtered nuclei); patients with tumour and adjacent normal only."),
                    notes, gene_ids=gids, title="Early-stage lung adenocarcinoma snRNA-seq, tumour and adjacent normal")


def GSE333596():
    """LUAD stage I / III scRNA-seq: per-GSM MEX in RAW.tar; matched tumour / adjacent normal libraries."""
    gse = "GSE333596"
    recs = fetch_all(gse, [f"{gse}_RAW.tar"])
    soft = soft_samples(gse)
    tf = tarfile.open(ROOT / gse / f"{gse}_RAW.tar")
    names = sorted(m.name for m in tf.getmembers())
    specs = []
    for gsm, v in soft.items():
        lib = re.sub(r"^GSM\d+_", "", _stem(names, gsm))          # library name = file stem (equals the GSM 'Library name')
        m = re.fullmatch(r"([TN])(\d[AB])|S(\d)([TN])", lib)
        assert m, lib
        if m.group(4):      # S1T/S1N, S3T/S3N: inconsistent libraries, not used (see notes)
            continue
        t, pid = m.group(1), m.group(2)
        specs.append(dict(gsm=gsm, sample=lib, patient=pid, sample_type=ST[t], chemistry="10x 3' v3.1", stage=v["chars"]["Stage"]))
    specs, dropped = _keep_paired(specs)
    parts, calls = [], []
    for s in specs:
        M, gn, ids, bc = _tar_10x(tf, names, _stem(names, s["gsm"]))
        mask, cut = call_cells(M)
        print(f"  {s['sample']}: {M.shape[0]} barcodes supplied, {int(mask.sum())} called (UMI >= {cut:.0f})", flush=True)
        calls.append(f"{s['sample']} {M.shape[0]} -> {int(mask.sum())} (cutoff {cut:.0f} UMI)")
        s["n_barcodes_supplied"] = int(M.shape[0])
        parts.append((M[mask], gn, ids, [b for b, k in zip(bc, mask) if k], s))
    X, genes, gids, obs = _assemble(parts, {})
    obs["author_cell_type"], obs["sorting"], obs["treatment"] = "", "", "naive"
    notes = [
        _pair_note(specs, dropped),
        "DROPPED LIBRARIES S1T/S1N and S3T/S3N (pairs S1, S3): the `*N` matrices hold 1,474,560 barcodes = two stacked 737,280-barcode raw libraries (suffix -1 and -2) "
        "and the `*T` matrices (16,046 and 61,481 barcodes) are a barcode subset of them with identical UMI counts (checked for S1: all 16,046 barcodes, equal UMI "
        "vector), i.e. `*T` is the filtered version of an aggregate that contains both libraries, so tumour and normal cannot be separated and the `*N` file is not a normal-only library. "
        "The four A/B pairs (T1A/N1A, T1B/N1B, T2A/N2A, T2B/N2B), separate raw libraries with 737,280 barcodes each, are used.",
        "CELL CALLING (rule added by this builder, not the authors'): the deposited matrices, although described as filtered, hold every barcode of the whitelist "
        "(737,280 barcodes per library for T1A-N2B). Cells were called per library by the Cell Ranger v2-style rule "
        "(cutoff = UMI of the barcode at the 99th percentile of the top 5,000 barcodes / 10, and >= 200 genes). Barcodes supplied -> called: " + "; ".join(calls) + ".",
        "patient: matched pairs by library name: T1A/N1A -> `1A`, T1B/N1B -> `1B`, T2A/N2A -> `2A`, T2B/N2B -> `2B` (GSM title: 'Adjacent normal tissue N1A "
        "(matched to T1A)'); (S1T/S1N and S3T/S3N not used.) Whether 1A/1B (2A/2B) are different patients is not stated in the record; "
        "the series design says stage I (n = 3) and stage III (n = 3), i.e. six patients in the record, of which four (1A, 1B, 2A, 2B) are usable here (S1, S3 dropped).",
        "sample_type: library name T / N (GSM `tissue`: lung tumor / lung). chemistry: Chromium Next GEM 3' v3.1 (GSM extract protocol). `stage` = GSM Stage.",
        "treatment: naive (GSM treatment_protocol: 'No treatment was applied prior to sample collection'). sorting: '' (viable cells, no sorting stated).",
    ]
    return finalize(gse, X, genes, obs, _prov(recs, filtering="Cell Ranger 7.0 raw-droplet matrices; cells called per library (see notes)."),
                    notes, gene_ids=gids, title="Lung adenocarcinoma stage I / III, tumour and matched adjacent normal")


def GSE153935():
    """NSCLC fibroblast study (Drop-seq): one merged genes x cells UMI table of all cells; samples from the cell-id prefix D<n>T / D<n>N."""
    gse = "GSE153935"
    f1 = f"{gse}_TLDS_AllCells.txt.gz"
    recs = fetch_all(gse, [f1])
    soft = soft_samples(gse)
    X, genes, cells = read_dense_gene_by_cell(ROOT / gse / f1)
    samp = pd.Series([c.split("_")[0] for c in cells])
    titles = {v["chars"]["donor number"] and re.fullmatch(r"Donor (\d+) (tumour|non-involved lung)", v["title"]).groups(): g for g, v in soft.items()}
    assert set(titles) >= {(m.group(1), "tumour" if m.group(2) == "T" else "non-involved lung")
                           for m in map(lambda x: re.fullmatch(r"D(\d+)([TN])", x), samp.unique())}, "sample prefixes not in GSM titles"
    obs = pd.DataFrame(index=cells)
    obs["sample"] = samp.values
    obs["patient"] = samp.str.extract(r"^(D\d+)")[0].values
    obs["sample_type"] = np.where(samp.str.endswith("T"), "tumour", "adjacent_normal")
    obs["chemistry"] = "Drop-seq"
    obs["author_cell_type"], obs["sorting"], obs["treatment"] = "", "", "unknown"
    keep_p = {p for p, g in obs.groupby("patient") if set(g.sample_type) == {"tumour", "adjacent_normal"}}
    dropped = sorted(set(obs.patient) - keep_p)
    m = obs.patient.isin(keep_p).to_numpy()
    X, obs = X[m], obs[m]
    notes = [
        f"Only donors with both a tumour (D<n>T) and a non-involved lung (D<n>N) sample are kept ({len(keep_p)}: {', '.join(sorted(keep_p))}); dropped (single type): {', '.join(dropped)}.",
        "patient / sample: prefix of the cell id (D2T_<barcode>) = GSM title 'Donor 2 tumour' / 'Donor 2 non-involved lung' (asserted present). The GSM characteristic "
        "`donor number` reads 'Donor 1' for donors 10-12 (typo); the title and the cell-id prefix are used. sample_type: T -> tumour, N (non-involved lung) -> adjacent_normal.",
        "Single file `TLDS_AllCells.txt.gz`: digital expression (UMI) matrix of 12 donors, hg19, merged and QC-filtered by the authors (empty droplets, high mito removed). "
        "The `Merged_StromalCells` file (third-party KU Leuven stromal cells) is not used.",
        "chemistry: Drop-seq (extract protocol), 12-mer barcodes. treatment: unknown (not stated). sorting: '' (the study is fibroblast-focused but no sorting is stated for these cells).",
    ]
    return finalize(gse, X, genes, obs, _prov(recs, filtering="Author-filtered all-cell UMI table; donors with tumour and non-involved lung only."),
                    notes, title="NSCLC and non-involved lung (Drop-seq), fibroblast study")


def GSE310916():
    """NSCLC, matched nuclei (snRNA 3' v3, 5 patients) and cells (5' v2, 4 patients): per-GSM MEX in RAW.tar (Seurat rds not downloaded)."""
    gse = "GSE310916"
    recs = fetch_all(gse, [f"{gse}_RAW.tar"])
    soft = soft_samples(gse)
    tf = tarfile.open(ROOT / gse / f"{gse}_RAW.tar")
    names = sorted(m.name for m in tf.getmembers())
    specs = []
    for gsm, v in soft.items():
        m = re.fullmatch(r"LUNG_([TN])(\d+)_(NUC|CELL)", v["title"])
        assert m, v["title"]
        t, pid, mod = m.groups()
        assert v["chars"]["tissue"] == {"T": "tumor tissue", "N": "normal tissue"}[t]
        specs.append(dict(gsm=gsm, sample=v["title"], patient=f"P{pid}", sample_type=ST[t],
                          chemistry="snRNA 10x 3' v3" if mod == "NUC" else "scRNA 10x 5' v2", modality="nuclei" if mod == "NUC" else "cells"))
    specs, dropped = _keep_paired(specs)
    X, genes, gids, obs = _assemble(_read_mex_specs(gse, tf, names, specs), {})
    obs["author_cell_type"], obs["sorting"], obs["treatment"] = "", "", "naive"
    notes = [
        _pair_note(specs, dropped),
        "Nuclei (`_NUC`, snRNA 3' v3, patients 35 36 38 39 44) and cells (`_CELL`, scRNA 5' v2, patients 35 36 39 44; the library file names end `_5PADT`) are "
        "separate samples of the same patient; `chemistry` records which and `modality` = nuclei / cells. Patient 38 has no cell library.",
        "If the 5' libraries carry an antibody-capture block, only 'Gene Expression' rows are kept (legacy `_tar_10x`).",
        "patient: `P<number>` from the GSM title. sample_type: T / N in the title = GSM `tissue` tumor tissue / normal tissue ('distant normal tissue'), asserted.",
        "treatment: naive (series design: none of the patients received neoadjuvant treatment). sorting: '' (cells: viable cells by Ficoll, cryopreserved; nuclei: DAPI-stained, no lineage sorting).",
        "The integrated Seurat rds (2 GB) is not read or downloaded.",
    ]
    return finalize(gse, X, genes, obs, _prov(recs, filtering="Cell Ranger 7.0.1 matrices as deposited (RAW.tar only)."),
                    notes, gene_ids=gids, title="NSCLC matched nuclei and cells, tumour and distant normal")


def GSE200972():
    """Multiple primary lung cancers: per-lesion tumour libraries and matched normal; MEX (P1-P3) and dense txt (P4) in RAW.tar."""
    gse = "GSE200972"
    recs = fetch_all(gse, [f"{gse}_RAW.tar"])
    soft = soft_samples(gse)
    tf = tarfile.open(ROOT / gse / f"{gse}_RAW.tar")
    members = {m.name: m for m in tf.getmembers()}
    names = sorted(members)
    specs = []
    for gsm, v in soft.items():
        m = re.fullmatch(r"(P\d)-([TN])([A-Z])(\d?)-scRNAseq", v["title"])
        assert m, v["title"]
        pid, t, site, rep = m.groups()
        specs.append(dict(gsm=gsm, sample=v["title"].replace("-scRNAseq", ""), patient=pid, sample_type=ST[t], chemistry="10x 3'", lesion=site + rep))
    specs, dropped = _keep_paired(specs)
    parts = []
    for s in specs:
        dense = [n for n in names if n.startswith(s["gsm"] + "_") and n.endswith("_matrix.tsv.gz")]
        if dense:
            M, gn, cells = _dense_member(gse, tf, members[dense[0]])
            parts.append((M, gn, None, cells, s))
        else:
            M, gn, ids, bc = _tar_10x(tf, names, _stem(names, s["gsm"]))
            parts.append((M, gn, ids, bc, s))
    X, genes, gids, obs = _assemble(parts, {})
    obs["author_cell_type"], obs["sorting"], obs["treatment"] = "", "", "unknown"
    notes = [
        _pair_note(specs, dropped),
        "patient: GSM title `P<n>-...`; the file names of P3 say `P8_*` and those of P4 `P4-2T1` etc.; the GSM titles (P3, P4-TS1/TS2/NS/TI/NI) are used and "
        "the GSM accession links title and file. sample_type: first letter after the patient (T / N). Every tumour lesion of a patient is 'tumour' and every "
        "normal sample 'adjacent_normal' (normal adjacent to each lesion); sample ids stay distinct (e.g. P1-TI, P1-TM); `lesion` = lesion code of the title "
        "(I, M, S plus replicate digit).",
        "P1-P3: Cell Ranger 4.0 filtered MEX. P4: dense genes x cells UMI txt per sample (no gene ids); the GSM files for P4 have no barcodes/features files.",
        "chemistry: `10x 3'` (the record names only 10x Genomics, Cell Ranger 4.0; kit version not stated). treatment: unknown (not stated). sorting: ''.",
    ]
    return finalize(gse, X, genes, obs, _prov(recs, filtering="Cell Ranger matrices as deposited."),
                    notes, gene_ids=gids if gids is not None else None, title="Multiple primary lung cancers, lesions and matched normal")


def GSE117570():
    """Early-stage NSCLC (P1-P4): per-sample dense genes x cells UMI txt in RAW.tar."""
    gse = "GSE117570"
    recs = fetch_all(gse, [f"{gse}_RAW.tar"])
    soft = soft_samples(gse)
    tf = tarfile.open(ROOT / gse / f"{gse}_RAW.tar")
    members = {m.name: m for m in tf.getmembers()}
    parts = []
    for gsm, v in soft.items():
        m = re.fullmatch(r"(P\d)_(Tumor|Normal)", v["title"])
        assert m, v["title"]
        assert v["chars"]["tissue"].lower() in {"tumor", "adjacent normal"}
        f = [n for n in members if n.startswith(gsm + "_")]
        assert len(f) == 1
        M, gn, cells = _dense_member(gse, tf, members[f[0]])
        s = dict(gsm=gsm, sample=v["title"], patient=m.group(1), sample_type="tumour" if m.group(2) == "Tumor" else "adjacent_normal",
                 chemistry="10x 3' v2", diagnosis=v["chars"]["diagnosis"], tnm=v["chars"].get("tnm status", ""))
        parts.append((M, gn, None, cells, s))
    X, genes, gids, obs = _assemble(parts, {})
    obs["author_cell_type"], obs["sorting"], obs["treatment"] = "", "", "naive"
    notes = [
        "patient: `P<n>` from the GSM title; sample_type: Tumor -> tumour, Normal ('tissue: Adjacent normal') -> adjacent_normal. `diagnosis` and `tnm` from the GSM "
        "characteristics (P2 is squamous, P1/P3/P4 adenocarcinoma).",
        "Files `*_processed_data.txt.gz` are dense genes x cells UMI tables (no gene ids; header has no leading label); the integer check of the legacy reader is "
        "applied to every value.",
        "chemistry: 10x Chromium Single Cell 3' (P/N 120237 = v2; Read1 26 bp) -> `10x 3' v2`. treatment: naive (GSM 'treatment-naive'). sorting: ''.",
    ]
    return finalize(gse, X, genes, obs, _prov(recs, filtering="cellranger-count UMI tables as deposited."),
                    notes, title="Early-stage NSCLC, tumour and adjacent normal")


def GSE193337():
    """Cribriform prostate adenocarcinoma: benign and tumour scRNA per patient, per-GSM MEX in RAW.tar (bulk RNA-seq files not read)."""
    gse = "GSE193337"
    recs = fetch_all(gse, [f"{gse}_RAW.tar"])
    soft = soft_samples(gse)
    tf = tarfile.open(ROOT / gse / f"{gse}_RAW.tar")
    names = sorted(m.name for m in tf.getmembers())
    specs = []
    for gsm, v in soft.items():
        m = re.fullmatch(r"scRNA-seq of patient (\d) (benign|tumor) sample", v["title"])
        if not m:
            continue
        specs.append(dict(gsm=gsm, sample=f"P{m.group(1)}{'n' if m.group(2) == 'benign' else 't'}", patient=f"P{m.group(1)}",
                          sample_type="adjacent_normal" if m.group(2) == "benign" else "tumour", chemistry="10x 3' v3.1",
                          pathology=v["chars"]["disease state"]))
    specs, dropped = _keep_paired(specs)
    X, genes, gids, obs = _assemble(_read_mex_specs(gse, tf, names, specs), {})
    obs["author_cell_type"], obs["sorting"], obs["treatment"] = "", "", "unknown"
    notes = [
        _pair_note(specs, dropped) + " The bulk RNA-seq (NEC / TEC) GSMs are not used.",
        "patient: GSM title 'scRNA-seq of patient <n> benign / tumor sample' -> `P<n>`; sample_type: benign -> adjacent_normal (paired benign 1x1 cm biopsy of the same "
        "prostate), tumor -> tumour. `pathology` = GSM disease state (benign / pT stage).",
        "chemistry: Chromium Next GEM 3' v3.1 (GSM extract protocol). treatment: unknown (radical prostatectomy specimens, no statement on prior therapy). "
        "sorting: '' (the protocol says viable cells were sorted, a live-cell sort without lineage selection; not recorded as sorting).",
        "Cell Ranger 3.1 matrices as deposited.",
    ]
    return finalize(gse, X, genes, obs, _prov(recs, filtering="Cell Ranger matrices as deposited."),
                    notes, gene_ids=gids, title="Cribriform prostate cancer, tumour and benign (scRNA part)")


def _seqwell_cohort(gse, specs, tf, members, sources_note, notes, title, filtering):
    parts = []
    for s in specs:
        f = [n for n in members if n.startswith(s["gsm"] + "_")]
        assert len(f) == 1, (s["gsm"], f)
        M, gn, cells = _dense_member(gse, tf, members[f[0]])
        parts.append((M, gn, None, cells, s))
    X, genes, gids, obs = _assemble(parts, {})
    return X, genes, obs


def GSE176031():
    """Prostate Seq-Well S3: radical prostatectomy N1/N2/T1/T2 libraries of 4 patients (PR5249, PR5251, PR5254, PR5261); dense dge in RAW.tar."""
    gse = "GSE176031"
    recs = fetch_all(gse, [f"{gse}_RAW.tar"])
    soft = soft_samples(gse)
    tf = tarfile.open(ROOT / gse / f"{gse}_RAW.tar")
    members = {m.name: m for m in tf.getmembers()}
    specs = []
    for gsm, v in soft.items():
        m = re.fullmatch(r"PA_(PR\d+)_([NT])(\d)_S\d+_L\d+", v["title"])
        if not m:        # biopsies, tumour-only prostatectomies without N/T label, organoids
            continue
        pid, t, rep = m.groups()
        specs.append(dict(gsm=gsm, sample=f"{pid}_{t}{rep}", patient=pid, sample_type=ST[t], chemistry="Seq-Well S3", geo_patient=v["chars"]["patient"]))
    specs, dropped = _keep_paired(specs)
    X, genes, obs = _seqwell_cohort(gse, specs, tf, members, "", [], "", "")
    obs["author_cell_type"], obs["sorting"], obs["treatment"] = "", "", "unknown"
    notes = [
        _pair_note(specs, dropped) + " Not used: the 6 biopsy libraries (Patients 1-3), the tumour-only prostatectomies (PR5186, PR5196, PR5199) and the unlabelled "
        "PR5269_1-4 libraries (no N/T in the title), and all organoid libraries.",
        "patient: PR id in the GSM title (GSM `patient` = Patient 7 to 10 is kept as `geo_patient`). sample_type: N<k> -> adjacent_normal, T<k> -> tumour; N1/N2 and "
        "T1/T2 are replicate Seq-Well arrays of the same sample (distinct sample ids).",
        "chemistry: Seq-Well S3 (GEO overview). Files `*_dge.txt.gz` are dense genes x cells raw gene counts (hg19 symbols; header `GENE`).",
        "treatment: unknown (not stated in the GEO record). sorting: ''.",
    ]
    return finalize(gse, X, genes, obs, _prov(recs, filtering="dge tables as deposited; radical-prostatectomy N / T libraries of 4 paired patients only."),
                    notes, title="Prostate cancer Seq-Well, paired normal and tumour radical prostatectomies")


def GSE293189():
    """Bladder histologic variants, Seq-Well: per-GSM dense dge in RAW.tar; the 2 paraganglioma patients' samples are 'other'."""
    gse = "GSE293189"
    recs = fetch_all(gse, [f"{gse}_RAW.tar"])
    soft = soft_samples(gse)
    tf = tarfile.open(ROOT / gse / f"{gse}_RAW.tar")
    members = {m.name: m for m in tf.getmembers()}
    specs = []
    for gsm, v in soft.items():
        m = re.fullmatch(r"Patient (\S+) \((.+?)\) (Paired Normal|Tumor)(\d?)(?: Replica(\d))?", v["title"])
        assert m, v["title"]
        pid, hist, kind, k, rep = m.groups()
        f = [n for n in members if n.startswith(gsm + "_")]
        assert len(f) == 1
        lib = re.fullmatch(rf"{gsm}_(.+)_dge\.txt\.gz", f[0]).group(1)          # library name = file stem
        if pid == "21226" and "HY_21262" in lib:        # title typo: library HY_21262 is a different (tumour-only) patient
            pid = "21262"
        stype = "adjacent_normal" if kind == "Paired Normal" else "tumour"
        if pid.startswith("PG_"):
            stype = "other"
        specs.append(dict(gsm=gsm, sample=lib, patient=pid, sample_type=stype, chemistry="Seq-Well", histology=hist,
                          tissue_label=kind + k))
    # pairing check on the underlying tissue labels, before 'other' is applied to paragangliomas
    pre = {}
    for s in specs:
        pre.setdefault(s["patient"], set()).add("normal" if s["tissue_label"].startswith("Paired") else "tumour")
    pairs = sorted(p for p, t in pre.items() if t == {"normal", "tumour"})
    keep = set(pairs)
    specs = [s for s in specs if s["patient"] in keep]
    X, genes, obs = _seqwell_cohort(gse, specs, tf, members, "", [], "", "")
    obs["author_cell_type"], obs["sorting"], obs["treatment"] = "", "", "unknown"
    notes = [
        f"Only patients with a Paired Normal and a Tumor sample in the record are kept ({len(pairs)}: {', '.join(pairs)}); tumour-only patients (11734, 12041, 12049, HG, "
        "PS, SG, 21262) are dropped.",
        "patient: GSM title `Patient <id> (<histology>)`. Patient FG has CIS and Nested tumours under one id and one paired normal (histology kept per sample in `histology`). "
        "Two titles say `Patient 21226 (Small large cell)` but their library names are HY_21262_*: treated as patient 21262 (tumour only, dropped).",
        "sample_type: `Paired Normal` -> adjacent_normal, `Tumor` -> tumour; all samples of the two paraganglioma patients (PG_M, PG_S) -> `other` (no epithelium of origin; "
        "A149 note 1). Replica and Tumor1/2/3 libraries keep distinct sample ids (library name).",
        "chemistry: Seq-Well. Files `*_dge.txt.gz` are dense genes x cells raw gene counts (hg19, header `GENE`; 12-mer barcodes, so cell ids are prefixed by the library). "
        "treatment: unknown (not stated). sorting: ''.",
    ]
    return finalize(gse, X, genes, obs, _prov(recs, filtering="dge tables as deposited; patients with paired normal and tumour only."),
                    notes, title="Bladder cancer histologic variants (Seq-Well), tumour and paired normal")


COHORTS = {k: v for k, v in list(globals().items()) if re.fullmatch(r"GSE\d+", k)}

if __name__ == "__main__":
    for g in sys.argv[1:]:
        print(f"=== {g}", flush=True)
        COHORTS[g]()
