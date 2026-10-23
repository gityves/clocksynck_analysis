"""Copy-number calls for epithelial cells, in the manner of inferCNV (Patel 2014, Tirosh 2016).

Profile of a cell: log-normalised expression of expressed autosomal genes minus a baseline, clipped, smoothed by
a running mean over `window` genes in genome order within each chromosome, centred on the cell's median.
Score: mean squared profile; call: `clonality` (correlation with the patient's own tumour consensus minus that with
other patients'), which separates clonal copy number from generic epithelial patterns and depth noise.

The baseline matters. Against immune and stromal cells (classic inferCNV), normal epithelium already looks
aberrant because epithelial genes cluster in the genome (keratin, surfactant loci): in adjacent-normal lung 11-49 %
of epithelial cells were called malignant. The baseline of an epithelial cell is therefore the mean of the study's
normal epithelium (adjacent-normal or healthy samples) of the *other* patients, and the null distribution of the
statistic is that of normal epithelial cells scored the same way (leave-one-patient-out). A cell is `malignant` above the
95th percentile of that null, `normal` below its 80th, otherwise `ambiguous`. Studies without normal epithelium get
no calls. On adjacent-normal vs tumour epithelium (2026-10-03): colorectal ERP149257 0-11 % vs 19-93 % malignant;
lung ERP153272 0-6 % vs 21-78 % in 4 of 5 patients (patient 21: 56 % of its normal, unresolved).

Gene positions: GENCODE v46 basic (`DATA/reference/gene_positions_gencode46.tsv`, docs/DATA.md).

Ported unchanged from legacy/src/sccc/prep/cnv.py (paper-only code: copy-number calls of the malignant epithelium of
the single-cell cancer cohorts, `study_calls`; the calls of the 29 cohorts are in `$SCCC_DATA/cnv/<study>.csv.gz`).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import scipy.sparse as sp

import sys
from pathlib import Path

sys.path[:0] = [str(Path(__file__).resolve().parents[3])]
from common.paths import DATA  # noqa: E402

POSITIONS = DATA / "reference" / "gene_positions_gencode46.tsv"


def positions(path=POSITIONS) -> pd.DataFrame:
    p = pd.read_csv(path, sep="\t", header=None, names=["gene_id", "symbol", "chrom", "start", "end", "type"])
    p = p[p.chrom.isin([f"chr{i}" for i in range(1, 23)])].drop_duplicates("gene_id")
    p["chr"] = p.chrom.str[3:].astype(int)
    return p.set_index("gene_id")[["symbol", "chr", "start"]]


def genome_order(genes, pos: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    """Columns of autosomal genes with a position, in genome order, and their chromosome. `genes`: Ensembl ids
    (version suffix ignored) or, if most are not, symbols."""
    ids = pd.Index([str(g).split(".")[0] for g in genes])
    if np.mean(ids.str.startswith("ENSG")) < 0.5:
        pos = pos.reset_index().drop_duplicates("symbol").set_index("symbol")
        ids = pd.Index([str(g).upper() for g in genes])
    cols = np.flatnonzero(ids.isin(pos.index) & ~ids.duplicated())
    p = pos.loc[ids[cols]]
    o = np.lexsort((p.start.to_numpy(), p.chr.to_numpy()))
    return cols[o], p.chr.to_numpy()[o]


def log_norm(X, total_umi, cols) -> np.ndarray:
    """log1p(counts per 10k of the cell's genome-wide UMI) for columns `cols`, dense float32."""
    X = sp.csr_matrix(X)[:, cols]
    E = X.multiply(1e4 / np.maximum(np.asarray(total_umi, float), 1.0)[:, None]).tocsr()
    return np.log1p(E.toarray().astype(np.float32))


def _running_mean(M: np.ndarray, chrom: np.ndarray, window: int) -> np.ndarray:
    out = np.empty_like(M)
    for c in np.unique(chrom):
        j = np.flatnonzero(chrom == c)
        h = min(window, len(j)) // 2
        cs = np.cumsum(np.pad(M[:, j], ((0, 0), (1, 0))), axis=1)
        lo = np.clip(np.arange(len(j)) - h, 0, len(j))
        hi = np.clip(np.arange(len(j)) + h + 1, 0, len(j))
        out[:, j] = (cs[:, hi] - cs[:, lo]) / (hi - lo)
    return out


def profiles(E: np.ndarray, baseline: np.ndarray, chrom: np.ndarray, expressed: np.ndarray, window: int = 101,
             clip: float = 3.0) -> np.ndarray:
    """Smoothed profiles of the cells in `E` (cells x genes in genome order) against `baseline` (genes)."""
    D = np.clip(E[:, expressed] - baseline[expressed], -clip, clip)
    P = _running_mean(D, chrom[expressed], window)
    return P - np.median(P, 1, keepdims=True)


def score(P: np.ndarray) -> np.ndarray:
    return np.mean(P ** 2, 1)


def corr(P: np.ndarray, c: np.ndarray) -> np.ndarray:
    Pc, cc = P - P.mean(1, keepdims=True), c - c.mean()
    return Pc @ cc / np.maximum(np.linalg.norm(Pc, axis=1) * np.linalg.norm(cc), 1e-12)


def clonality(P: dict, kind: dict, patient: dict) -> dict:
    """Per sample: correlation of each cell with its own patient's tumour consensus minus the mean correlation with
    the other patients' consensuses. Consensus = mean profile of the upper half (by score) of a patient's tumour
    epithelial cells. Generic epithelial patterns correlate with every consensus and cancel; clonal copy number
    is patient-specific."""
    tum = {}
    for s, M in P.items():
        if kind[s] == "tumour":
            tum.setdefault(patient[s], []).append(M)
    cons = {}
    for p, Ms in tum.items():
        M = np.vstack(Ms)
        sc = score(M)
        cons[p] = M[sc >= np.median(sc)].mean(0)
    out = {}
    for s, M in P.items():
        own = patient[s]
        if own not in cons:
            continue
        other = [corr(M, c) for q, c in cons.items() if q != own]
        out[s] = corr(M, cons[own]) - np.mean(other, 0)
    return out


def call(stat: np.ndarray, null: np.ndarray) -> np.ndarray:
    """`malignant` above the 95th percentile of the null (normal epithelium), `normal` below its 80th."""
    lo, hi = np.quantile(null, 0.80), np.quantile(null, 0.95)
    return np.where(stat > hi, "malignant", np.where(stat < lo, "normal", "ambiguous"))


def study_calls(samples, pos: pd.DataFrame | None = None, min_mean: float = 0.1, min_null: int = 200) -> pd.DataFrame | None:
    """Calls for the epithelial cells of a study, in two passes over `samples`.

    `samples`: list of (sample, patient, sample_type, load) where `load()` returns (X, genes, total_umi, compartment,
    barcodes) and compartment comes from `cancer_common.compartments`. Normal epithelium = epithelial cells of adjacent_normal or
    healthy samples. Profiles are taken against the normal epithelium of the other patients; the call uses
    `clonality`, with the normal epithelium of patients that also have a tumour sample as null. Returns one row per
    epithelial cell of patients with a tumour sample (sample, barcode, patient, sample_type, cnv_score, clonality,
    cnv_call), or None when fewer than `min_null` null cells or fewer than 2 tumour patients exist.
    """
    pos = positions() if pos is None else pos
    normal = ("adjacent_normal", "healthy")
    cols = chrom = None
    s_norm, n_norm, s_counts, n_epi = {}, {}, 0.0, 0
    for sample, patient, kind, load in samples:  # pass 1: normal-epithelium sums per patient, expressed genes
        X, genes, umi, comp, _ = load()
        if cols is None:
            cols, chrom = genome_order(genes, pos)
        e = comp == "epithelial"
        if not e.any():
            continue
        s_counts = s_counts + np.asarray(sp.csr_matrix(X)[e][:, cols].sum(0)).ravel()
        n_epi += int(e.sum())
        if kind in normal:
            E = log_norm(X[e], umi[e], cols)
            s_norm[patient] = s_norm.get(patient, 0) + E.sum(0)
            n_norm[patient] = n_norm.get(patient, 0) + len(E)
    if len(n_norm) < 2:
        return None
    expressed = s_counts / n_epi >= min_mean
    tot_s, tot_n = sum(s_norm.values()), sum(n_norm.values())
    P, kinds, pats, bcs = {}, {}, {}, {}
    for sample, patient, kind, load in samples:  # pass 2: profiles against the other patients' normal epithelium
        X, genes, umi, comp, barcodes = load()
        e = comp == "epithelial"
        if e.sum() < 5:
            continue
        base = (tot_s - s_norm.get(patient, 0)) / (tot_n - n_norm.get(patient, 0))
        P[sample] = profiles(log_norm(X[e], umi[e], cols), base, chrom, expressed)
        kinds[sample], pats[sample], bcs[sample] = kind, patient, np.asarray(barcodes)[e]
    if len({pats[s] for s in P if kinds[s] == "tumour"}) < 2:
        return None
    clon = clonality(P, kinds, pats)
    out = pd.concat([pd.DataFrame({"sample": s, "barcode": bcs[s], "patient": pats[s], "sample_type": kinds[s],
                                   "cnv_score": score(P[s]), "clonality": c}) for s, c in clon.items()],
                    ignore_index=True)
    null = out.clonality[out.sample_type.isin(normal)].to_numpy()
    if len(null) < min_null:
        return None
    out["cnv_call"] = call(out.clonality.to_numpy(), null)
    return out
