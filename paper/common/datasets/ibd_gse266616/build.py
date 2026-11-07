"""Crohn's disease terminal ileum / ascending colon biopsies, inflamed vs adjacent / non-inflamed (GSE266616; Zhou, Xie,
Basu et al., bioRxiv 10.1101/2025.03.26.645552) and, with the same cell rules, TAURUS (GSE282122 / Zenodo 13768607):
sample x class pseudobulks of the 14 clock genes for amendment A230 (epithelial offset against the stromal host).

Source:   `$SCCC_DATA/gse266616/raw/` (GEO series matrices + `GSE266616_RAW.tar`; see `gse266616/PROVENANCE.md`),
          unpacked by `unpack()` to `$SCCC_DATA/gse266616/samples/<GSM>_<title>/` (10x matrix.mtx / features / barcodes);
          TAURUS: `$SCCC_DATA/ibd/TAURUS/TAURUS_raw_counts_annotated_final.h5ad` (raw counts, read in row chunks).
Used by:  `paper/exploratory/A230_ibd_local_offset/scripts/compute_a230.py`.

Choices (A230, fixed before any count was read)
- no author annotation: cells classed by `cancer_common.compartments.marker_compartments` (the `cancer_sc` rule: epithelial,
  fibroblast, endothelial, mural, myeloid, T_NK, B_plasma, mast; cells whose second-best class score is >= 0.7 x the best, or
  whose best score is <= 0.3, are `unassigned` = doublet-like / unclear and dropped);
- unit (donor) = sample (biopsy); compartment = `gut:<class>`; pseudobulks of < 50 cells dropped; samples with median UMI per
  cell < 1,500 dropped (median over all cells of the sample, as `cancer_common.studies`); stratum = site;
- donor table: patient, site, region, condition, arm (`inflamed` / `reference`), reference type (`adjacent` / `non`),
  `patient_site`, median UMI, cells. Arms (`arms`): per patient and site an inflamed sample against the same patient's
  adjacent-to-inflamed (preferred, if both exist) or non-inflamed sample at that site; controls and unpaired samples have
  no arm. TAURUS: arm = inflamed (`inflammation == Inflamed`) vs `reference` (Non_Inflamed), every sample with an arm.
"""

from __future__ import annotations

import gzip
import io
import re
import sys
import tarfile
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

HERE = Path(__file__).resolve()
sys.path[:0] = [str(HERE.parents[1]), str(HERE.parents[3])]
import clocksync as cs  # noqa: E402
from common.paths import DATA  # noqa: E402
from natexp_common import clock_adata, standard  # noqa: E402

from cancer_common.compartments import marker_compartments  # noqa: E402

ROOT = DATA / "gse266616"
RAW, SAMPLES = ROOT / "raw", ROOT / "samples"
TAURUS_H5 = DATA / "ibd" / "TAURUS" / "TAURUS_raw_counts_annotated_final.h5ad"
MIN_CELLS, MIN_UMI = 50, 1500


# ---- sample table ------------------------------------------------------------------------------------------------
def sample_table() -> pd.DataFrame:
    """One row per GSE266616 sample (index = GEO title, e.g. `HA48TI_inf`): gsm, patient (from the title: the patient
    field of the last three samples is malformed), site (TI / AC), region (inf / adjacent / non / control), condition,
    treatment, age."""
    rows = []
    for f in sorted(RAW.glob("*series_matrix.txt.gz")):
        with gzip.open(f, "rt") as h:
            lines = [ln.rstrip("\n").split("\t") for ln in h if ln.startswith("!Sample_")]
        d = {}
        for ln in lines:
            vals = [v.strip('"') for v in ln[1:]]
            key = ln[0]
            if key == "!Sample_characteristics_ch1":
                name = vals[0].split(":")[0]
                key = f"char_{name}"
                vals = [v.split(": ", 1)[1] if ": " in v else v for v in vals]
            d[key] = vals
        t = pd.DataFrame({"title": d["!Sample_title"], "gsm": d["!Sample_geo_accession"],
                          "tissue": d["char_tissue"], "condition": d["char_condition"], "age": d["char_age"],
                          "treatment": d["char_treatment (biological)"], "region": d["char_sample region"]})
        rows.append(t)
    t = pd.concat(rows, ignore_index=True)
    m = t.title.str.extract(r"^([A-Z]+\d+)(AC|TI)")
    t["patient"], t["site"] = m[0], m[1]
    assert t.site.notna().all() and (t.site.map({"AC": "ascending_colon", "TI": "terminal_ileum"}) == t.tissue).all()
    return t.set_index("title")


def arms(st: pd.DataFrame) -> pd.DataFrame:
    """Inflamed vs reference samples: per CD patient and site, the inflamed sample against the adjacent-to-inflamed
    sample (non-inflamed if there is none). Columns: patient, site, sample, arm, ref_type, `patient_site`. Raises if a
    patient-site has several inflamed or several candidate reference samples of the chosen type."""
    out = []
    for (p, s), g in st[st.condition == "CD"].groupby(["patient", "site"]):
        inf, ad, non = (g.index[g.region == r].tolist() for r in ("inf", "adjacent", "non"))
        if not inf or not (ad or non):
            continue
        ref, typ = (ad, "adjacent") if ad else (non, "non")
        if len(inf) != 1 or len(ref) != 1:
            raise ValueError(f"{p} {s}: {len(inf)} inflamed, {len(ref)} {typ} samples")
        out += [dict(patient=p, site=s, sample=inf[0], arm="inflamed", ref_type=typ, patient_site=f"{p}|{s}"),
                dict(patient=p, site=s, sample=ref[0], arm="reference", ref_type=typ, patient_site=f"{p}|{s}")]
    return pd.DataFrame(out).set_index("sample")


# ---- raw files (untrusted: parsed, never executed) ---------------------------------------------------------------
def unpack(force: bool = False) -> None:
    """Unpack `GSE266616_RAW.tar` and the per-sample `.tar.gz` into `samples/<GSM>_<title>/` (data-only members)."""
    SAMPLES.mkdir(parents=True, exist_ok=True)
    with tarfile.open(RAW / "GSE266616_RAW.tar") as outer:
        for m in outer:
            if not m.isfile():
                continue
            dest = SAMPLES / re.sub(r"\.tar\.gz$", "", m.name)
            if dest.exists() and not force:
                continue
            dest.mkdir(parents=True, exist_ok=True)
            with tarfile.open(fileobj=io.BytesIO(outer.extractfile(m).read()), mode="r:gz") as inner:
                for im in inner:
                    if im.isfile():
                        (dest / Path(im.name).name).write_bytes(inner.extractfile(im).read())


def _folder(gsm: str) -> Path:
    f = sorted(SAMPLES.glob(f"{gsm}_*"))
    if len(f) != 1:
        raise FileNotFoundError(f"{gsm}: {len(f)} folders in {SAMPLES}")
    return f[0]


def read_sample(gsm: str):
    """(X cells x genes CSR, symbols, barcodes) of one unpacked 10x sample (mtx, features / genes, barcodes)."""
    d = _folder(gsm)
    pick = lambda stem: next(p for p in sorted(d.iterdir()) if stem in p.name)  # noqa: E731
    op = lambda p: gzip.open(p, "rt") if p.suffix == ".gz" else open(p)  # noqa: E731
    with op(pick("matrix")) as h:
        M = scipy.io.mmread(h).tocsc()
    with op(pick("features") if any("features" in p.name for p in d.iterdir()) else pick("genes")) as h:
        feat = pd.read_csv(h, sep="\t", header=None)
    with op(pick("barcodes")) as h:
        bc = pd.read_csv(h, sep="\t", header=None)[0].astype(str).to_numpy()
    sym = feat[1].astype(str).to_numpy() if feat.shape[1] > 1 else feat[0].astype(str).to_numpy()
    X = sp.csr_matrix(M.T)
    if X.shape != (len(bc), len(sym)):
        raise ValueError(f"{gsm}: matrix {X.shape} vs {len(bc)} barcodes, {len(sym)} features")
    return X, sym, bc


def _pseudobulk(X, sym, umi, bc, patient: str) -> tuple[pd.DataFrame, dict]:
    comp = marker_compartments(X, sym, umi)
    obs = pd.DataFrame({"donor": patient, "compartment": ["gut:" + c for c in comp]}, index=pd.Index(bc).astype(str))
    t = cs.pseudobulk(clock_adata(X, sym, umi, obs), donor="donor", compartment="compartment", depth="depth", min_cells=1)
    info = dict(n_cells=len(umi), median_umi=float(np.median(umi)), n_unassigned=int((comp == "unassigned").sum()))
    return t, info


def _finish(parts, infos, donors):
    pb = pd.concat(parts, ignore_index=True)
    pb = pb[(pb.n_cells >= MIN_CELLS) & (pb.compartment != "gut:unassigned")]
    keep = set(infos.index[infos.median_umi >= MIN_UMI])
    pb = pb[pb.donor.isin(keep)].reset_index(drop=True)
    pb["stratum"] = pb.donor.map(donors.site)
    donors = donors.join(infos)
    donors["kept"] = donors.index.isin(keep)
    return standard(pb), donors


def build(which: str = "gse266616", samples: list | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    """(pseudobulk, donors) of `gse266616` (all 85 samples; `samples`: restrict to these titles) or `taurus`. All classes with
    >= 50 cells in samples with median UMI >= 1,500; donor = sample; donor table in the module docstring."""
    return _gse266616(samples) if which == "gse266616" else _taurus(samples)


def _gse266616(samples):
    st = sample_table()
    ar = arms(st)
    st = st.join(ar[["arm", "ref_type", "patient_site"]])
    if samples is not None:
        st = st.loc[samples]
    parts, infos = [], {}
    for title, r in st.iterrows():
        X, sym, bc = read_sample(r.gsm)
        umi = np.asarray(X.sum(1)).ravel()
        t, infos[title] = _pseudobulk(X, sym, umi, bc, title)
        parts.append(t)
        print(f"{title}: {infos[title]['n_cells']} cells, median UMI {infos[title]['median_umi']:.0f}", flush=True)
    donors = st.rename_axis("donor")
    return _finish(parts, pd.DataFrame(infos).T.astype({"n_cells": int, "n_unassigned": int, "median_umi": float})
                   .rename_axis("donor"), donors.assign(site=donors.site))


def _taurus(samples):
    import h5py
    from anndata.experimental import read_elem

    sm = pd.read_csv(DATA / "ibd" / "TAURUS" / "pseudobulk" / "samples.csv").set_index("sample")
    sm["arm"] = sm.inflammation.map({"Inflamed": "inflamed", "Non_Inflamed": "reference"})
    sm["ref_type"] = np.where(sm.arm == "reference", "non", "")
    sm = sm.rename(columns={"site": "site"})
    if samples is not None:
        sm = sm.loc[samples]
    parts, infos = [], {}
    with h5py.File(TAURUS_H5, "r") as f:
        obs = read_elem(f["obs"])
        var = read_elem(f["var"])
        sym = (var["gene_symbol"] if "gene_symbol" in var
               else pd.Series(var.index, index=var.index)).astype(str).to_numpy()
        smp = obs["sample_id"].astype(str)
        indptr = f["X/indptr"][:]
        for s in sm.index:
            rows = np.flatnonzero(smp.to_numpy() == s)
            if not len(rows):
                continue
            if not (np.diff(rows) == 1).all():  # non-contiguous: read row by row blocks
                blocks = np.split(rows, np.flatnonzero(np.diff(rows) > 1) + 1)
            else:
                blocks = [rows]
            mats = []
            for b in blocks:
                a, z = indptr[b[0]], indptr[b[-1] + 1]
                mats.append(sp.csr_matrix((f["X/data"][a:z], f["X/indices"][a:z], indptr[b[0]:b[-1] + 2] - a),
                                          shape=(len(b), len(sym))))
            X = sp.vstack(mats).tocsr()
            umi = np.asarray(X.sum(1)).ravel()
            t, infos[s] = _pseudobulk(X, sym, umi, obs.index[rows].astype(str), s)
            parts.append(t)
            print(f"{s}: {infos[s]['n_cells']} cells, median UMI {infos[s]['median_umi']:.0f}", flush=True)
    donors = sm.loc[list(infos)].rename_axis("donor")
    return _finish(parts, pd.DataFrame(infos).T.astype({"n_cells": int, "n_unassigned": int, "median_umi": float})
                   .rename_axis("donor"), donors)


if __name__ == "__main__":
    import argparse

    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("which", choices=["gse266616", "taurus"])
    ap.add_argument("out", help="output directory (pseudobulk.csv, donors.csv)")
    ap.add_argument("--unpack", action="store_true", help="unpack the GEO tar first (gse266616)")
    a = ap.parse_args()
    if a.unpack and a.which == "gse266616":
        unpack()
    pb, donors = build(a.which)
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    pb.to_csv(out / "pseudobulk.csv", index=False)
    donors.to_csv(out / "donors.csv")
    print(f"{a.which}: {pb.donor.nunique()} samples kept of {len(donors)}, {len(pb)} pseudobulks")
