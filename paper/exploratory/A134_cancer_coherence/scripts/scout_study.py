"""A134 scout, per cohort: whole-transcriptome per-patient covariates of the malignant and adjacent-normal epithelium of
the A87 patients (feasibility only; no regression, no coherence).

Per (patient, group in {malignant, normal}) of the patients with malignant and adjacent-normal epithelium and a stromal host
(the 79 of A87; `cancer_sc.paired_epithelium`):
  cells; mean copy-number profile value at each of the 14 clock-gene loci (inferCNV-like profile of
  `cancer_common.cnv`: log1p CP10k minus the mean of the other patients' normal epithelium, clipped at +-3, running mean
  of 101 expressed genes in genome order, centred on the cell median; own gene removed from its window; a gene that is
  not expressed takes the window of the nearest expressed gene on its chromosome); fraction of cycling cells (Seurat
  phase call on the Tirosh S / G2M lists, as A24) and mean MYC_TARGETS_V1 / HYPOXIA scores (log1p CP10k of the set minus
  expression-matched controls, as `cells_cancer.scores`).
Malignant = stored copy-number call `malignant` of the tumour epithelium; normal = adjacent-normal epithelium.

    python scout_study.py STUDY      -> output/scout_cells_<STUDY>.csv, output/scout_geneset_<STUDY>.csv
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

PAPER = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(PAPER / "common" / "datasets"))
sys.path.insert(0, str(PAPER / "results" / "figure_5" / "scripts" / "cancer"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scipy.sparse as sp  # noqa: E402

from clocksync.config import limit_threads  # noqa: E402
from clocksync.genes import CLOCK_GENES, canonical, gene_columns  # noqa: E402
from common import data  # noqa: E402
from common.paths import DATA  # noqa: E402

from cancer_common import cnv  # noqa: E402
from cancer_common import studies as S  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "output"
GMT = DATA / "reference" / "h.all.v2024.1.Hs.symbols.gmt"
GENES = list(CLOCK_GENES)
POS_NAME = {"ARNTL": "BMAL1"}  # GENCODE v46 symbol
SETS = ("HALLMARK_MYC_TARGETS_V1", "HALLMARK_HYPOXIA")
CHUNK = 4000


def hallmark(gmt=GMT) -> dict:
    out = {}
    for line in Path(gmt).read_text().splitlines():
        f = line.split("\t")
        if f[0] in SETS:
            out[f[0]] = [g for g in f[2:] if canonical(g) not in {canonical(x) for x in GENES}]
    return out


def locus_table(pos: pd.DataFrame) -> pd.DataFrame:
    p = pos.reset_index().drop_duplicates("symbol").set_index("symbol")
    return pd.DataFrame({g: p.loc[POS_NAME.get(g, g), ["chr", "start"]] for g in GENES}).T


def ordered(var, sym, pos):
    """Positioned autosomal genes in genome order: columns, chromosome, start, symbol."""
    cols, chrom = cnv.genome_order(var, pos)
    ids = pd.Index([str(g).split(".")[0] for g in var])
    if np.mean(ids.str.startswith("ENSG")) < 0.5:
        q = pos.reset_index().drop_duplicates("symbol").set_index("symbol")
        start = q.loc[[str(var[c]).upper() for c in cols], "start"].to_numpy()
    else:
        start = pos.loc[ids[cols], "start"].to_numpy()
    return cols, chrom, start


def main(study: str):
    limit_threads(4)
    warnings.filterwarnings("ignore")
    import importlib.util
    spec = importlib.util.spec_from_file_location("cb", PAPER / "common" / "datasets" / "cells_cancer" / "build.py")
    cb = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cb)
    pos = cnv.positions()
    loci = locus_table(pos)
    pts = pd.read_csv(OUT / "scout_patients79.csv")
    keep = set(pts.donor[pts.study == study].str.split("|", n=1).str[1])
    us = S.majority_chemistry(S.units(study)[0])
    calls = S.read_calls(study)
    hm = hallmark()

    # pass 1: expressed genes, normal-epithelium baselines, mean expression of all genes (controls)
    cols = chrom = start = sym = None
    s_counts, n_epi, s_norm, n_norm, lsum = 0.0, 0, {}, {}, 0.0
    for sample, patient, kind, load, chem, _ in us:
        X, var, umi, comp, bc, sym = load()
        if cols is None:
            cols, chrom, start = ordered(var, sym, pos)
        e = comp == "epithelial"
        if not e.any():
            continue
        Xe = sp.csr_matrix(X)[e]
        s_counts = s_counts + np.asarray(Xe[:, cols].sum(0)).ravel()
        n_epi += int(e.sum())
        M = Xe.multiply(1e4 / np.maximum(umi[e], 1)[:, None]).tocsr()
        M.data = np.log1p(M.data)
        lsum = lsum + np.asarray(M.sum(0)).ravel()
        if kind in ("adjacent_normal", "healthy"):
            E = cnv.log_norm(Xe, umi[e], cols)
            s_norm[patient] = s_norm.get(patient, 0) + E.sum(0)
            n_norm[patient] = n_norm.get(patient, 0) + len(E)
    expressed = s_counts / n_epi >= 0.1
    ex = np.flatnonzero(expressed)
    tot_s, tot_n = sum(s_norm.values()), sum(n_norm.values())
    log_mean = lsum / n_epi

    # gene-set coverage
    symc = pd.Index([canonical(s) for s in sym])
    cov = []
    for name, g in hm.items():
        pres = [x for x in g if canonical(x) in symc]
        fin = [x for x in pres if np.isfinite(log_mean[symc.get_loc(canonical(x))]) and log_mean[symc.get_loc(canonical(x))] > 0]
        cov.append(dict(study=study, set=name, genes_in_set=len(g), present=len(pres), expressed=len(fin),
                        matrix_genes=len(sym), positioned_genes=len(cols), expressed_positioned=len(ex)))
    pd.DataFrame(cov).to_csv(OUT / f"scout_geneset_{study}.csv", index=False)

    # clock-gene loci in the expressed, ordered list
    chrom_e, start_e = chrom[ex], start[ex]
    sym_o = np.array([canonical(sym[c]) for c in cols])[ex]
    lp = []
    for g in GENES:
        c0, s0 = int(loci.loc[g, "chr"]), int(loci.loc[g, "start"])
        j = np.flatnonzero(chrom_e == c0)
        own = [k for k in j if sym_o[k] == canonical(g)]
        if own:
            k, dist = own[0], 0
        else:
            k = j[np.argmin(np.abs(start_e[j] - s0))]
            dist = int(abs(start_e[k] - s0))
        jj = j  # genes of the chromosome, in order
        i = int(np.flatnonzero(jj == k)[0])
        h = min(101, len(jj)) // 2
        lo, hi = max(i - h, 0), min(i + h + 1, len(jj))
        lp.append(dict(gene=g, k=k, own=bool(own), dist_bp=dist, win=hi - lo))
    lp = pd.DataFrame(lp)

    sg = cb.control_sets(sym, log_mean, np.random.default_rng(0), {"S": cb.S_GENES, "G2M": cb.G2M_GENES})
    hs = cb.control_sets(sym, log_mean, np.random.default_rng(0), hm)

    rows = []
    for sample, patient, kind, load, chem, _ in us:
        if patient not in keep or kind not in ("tumour", "adjacent_normal"):
            continue
        X, var, umi, comp, bc, _ = load()
        e = np.flatnonzero(comp == "epithelial")
        if kind == "tumour":
            call = calls.reindex(pd.MultiIndex.from_arrays([[sample] * len(e), np.asarray(bc).astype(str)[e]])).to_numpy()
            e = e[call == "malignant"]
            grp = "malignant"
        else:
            grp = "normal"
        if len(e) == 0:
            continue
        base = (tot_s - s_norm.get(patient, 0)) / (tot_n - n_norm.get(patient, 0))
        Xs = sp.csr_matrix(X)[e]
        for a in range(0, len(e), CHUNK):
            Xc, uc = Xs[a:a + CHUNK], umi[e][a:a + CHUNK]
            D = np.clip(cnv.log_norm(Xc, uc, cols)[:, ex] - base[ex], -3.0, 3.0)
            Sm = cnv._running_mean(D, chrom_e, 101)
            med = np.median(Sm, 1)
            L = np.empty((len(D), len(GENES)))
            for q, r in lp.iterrows():
                k = int(r.k)
                if r.own:
                    j = np.flatnonzero(chrom_e == chrom_e[k])
                    i = int(np.flatnonzero(j == k)[0])
                    h = min(101, len(j)) // 2
                    cnt = min(i + h + 1, len(j)) - max(i - h, 0)
                    L[:, q] = (Sm[:, k] * cnt - D[:, k]) / (cnt - 1) - med
                else:
                    L[:, q] = Sm[:, k] - med
            sc = cb.scores(Xc, uc, {**sg, **hs})
            ph = np.where((sc.S > sc.G2M) & (sc.S > 0), "S", np.where((sc.G2M >= sc.S) & (sc.G2M > 0), "G2M", "G1"))
            rows.append(pd.DataFrame(L, columns=[f"cn_{g}" for g in GENES]).assign(
                patient=patient, group=grp, sample=sample, cycling=(ph != "G1").astype(float), S=sc.S.to_numpy(),
                G2M=sc.G2M.to_numpy(), myc=sc[SETS[0]].to_numpy(), hypoxia=sc[SETS[1]].to_numpy()))
    c = pd.concat(rows, ignore_index=True)
    g = c.drop(columns="sample").groupby(["patient", "group"])
    cn = [f"cn_{x}" for x in GENES]
    agg = pd.concat([g[cn].mean().add_suffix("_mean"), g[cn].median().add_suffix("_median"),
                     g[["cycling", "S", "G2M", "myc", "hypoxia"]].mean()], axis=1)
    n = c.groupby(["patient", "group"]).size().rename("cells")
    out = pd.concat([n, agg], axis=1).reset_index().assign(study=study)
    out.to_csv(OUT / f"scout_cells_{study}.csv", index=False)
    lp.drop(columns="k").assign(study=study).to_csv(OUT / f"scout_loci_{study}.csv", index=False)
    print(study, "patients", c.patient.nunique(), "rows", len(out), flush=True)


if __name__ == "__main__":
    main(sys.argv[1])
