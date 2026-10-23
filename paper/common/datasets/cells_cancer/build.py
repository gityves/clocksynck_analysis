"""Per-cell tables of the single-cell cancer cohorts (cell-level inputs of A18 part 2 for cancer, A24 and the cycling
follow-ups T1-T3, C1-C4 of the 2026-10-04 amendments).

Source:   the same cohorts, samples, chemistry rule and cell classes as `cancer_sc` (`cancer_common.studies.units`,
          `majority_chemistry`); cell by cell instead of summed.
Tables (one pair per cohort; `$SCCC_DATA/cells/cancer/<study>_{dispersion,epithelial,host}.csv.gz`):
- `dispersion`: every classified cell of the eligible samples: the 14 clock-gene counts, `resection` (patient), `sample`,
  `barcode`, `compartment` (`<sample type>:<class>`), `stratum` (chemistry), `depth` (genome-wide UMI) and, for T/NK cells,
  `subtype` (CD8 / CD4 / T / NK from CD3, CD8A/B, CD4, NKG7/GNLY counts). Input of the dispersion D (`cells.vonmises`).
- `epithelial`: every epithelial cell: the 14 clock counts (`ARNTL` ...), the counts of the first expression-matched pseudo-clock
  gene set (`p0_<GENE>`, T2), Tirosh S and G2M module scores (Seurat `CellCycleScoring` style: control genes of the same
  expression bin, S = first 43 genes of the cycling list, G2M the rest), the combined `cycling` score and the other program
  scores (hypoxia, stress, interferon, pEMT), `erbb2` (log CP10k), `resection`, `kind` (tumour / adjacent_normal / ...),
  `stratum`, `depth`.
- `host`: pseudobulks of every non-epithelial compartment (>= 50 cells) of the same samples, clock counts, `depth`, `n_cells`.
`load(study, table)` reads the stored table of `$SCCC_DATA/cells/cancer`, else the legacy file of the run
(`outputs/celldisp/<study>_cells.csv.gz`, `outputs/states/cells/<study>[_host].csv.gz`), else builds it.
Replaces: legacy/scripts/x4_cell_dispersion.py (`extract`), s2_cellcycle.py (`collect`), s1_states.py (`collect`),
          legacy/src/sccc/prep/states.py (programs, `control_sets`, `scores`; copied below, "local: not in clocksync").
Check:    `python build.py check STUDY`: rebuilds the tables of a (small) cohort and compares them with the legacy files
          (clock counts, S, G2M, erbb2, host; see `INDEX_natexp_cancer.md`).

    python build.py build STUDY [--tables dispersion epithelial]       # needs the cohort's raw files (see cancer_common.studies)
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp

HERE = Path(__file__).resolve()
sys.path[:0] = [str(HERE.parents[1]), str(Path(__file__).resolve().parents[3])]
import clocksync as cs  # noqa: E402
from clocksync.genes import CLOCK_GENES, canonical, gene_columns  # noqa: E402
from clocksync.paired import pseudo_waveforms  # noqa: E402
from common.paths import DATA  # noqa: E402

from cancer_common import studies as S  # noqa: E402

OUT = DATA / "cells" / "cancer"
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))  # paper/
from common.paths import ARCHIVE  # noqa: E402
LEGACY = ARCHIVE / "outputs"
GENES = list(CLOCK_GENES)
SUB = ("CD3D", "CD3E", "CD8A", "CD8B", "CD4", "NKG7", "GNLY")
MIN_HOST = 50

# ---- local: not in clocksync (legacy/src/sccc/prep/states.py, verbatim) ------------------------------------------------
CLOCK_LIKE = set(GENES) | {
    "ARNTL2", "BHLHE40", "BHLHE41", "NFIL3", "HLF", "RORB", "RORC", "CIPC", "TIMELESS", "CSNK1D", "CSNK1E",
    "USP2", "TSPAN4", "WEE1", "TSC22D3", "KLF9", "KLF15", "ZBTB16", "FKBP5", "PDK4", "TXNIP", "SGK1", "DDIT4",
}
PROGRAMS = {
    "cycling": ["MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2", "MCM6", "CDCA7", "DTL", "PRIM1",
                "UHRF1", "HELLS", "RFC2", "RPA2", "NASP", "RAD51AP1", "GMNN", "WDR76", "SLBP", "CCNE2", "UBR7", "POLD3",
                "MSH2", "ATAD2", "RAD51", "RRM2", "CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM", "CASP8AP2", "USP1",
                "CLSPN", "POLA1", "CHAF1B", "BRIP1", "E2F8", "HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A",
                "NDC80", "CKS2", "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "SMC4", "CCNB2", "CKAP2L", "CKAP2",
                "AURKB", "BUB1", "KIF11", "ANP32E", "TUBB4B", "GTSE1", "KIF20B", "HJURP", "CDCA3", "CDC20", "TTK",
                "CDC25C", "KIF2C", "RANGAP1", "NCAPD2", "DLGAP5", "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR", "AURKA",
                "PSRC1", "ANLN", "LBR", "CKAP5", "CENPE", "CTCF", "NEK2", "G2E3", "GAS2L3", "CBX5", "CENPA"],
    "hypoxia": ["VEGFA", "SLC2A1", "PGK1", "CA9", "BNIP3", "NDRG1", "ADM", "ANKRD37", "P4HA1", "LDHA", "ALDOA", "MXI1",
                "PFKFB3", "EGLN3", "ANGPTL4", "ERO1A", "FAM162A", "PFKFB4", "BNIP3L", "ENO2"],
    "stress": ["FOS", "FOSB", "JUN", "JUNB", "EGR1", "ATF3", "IER2", "IER3", "DUSP1", "ZFP36", "KLF6", "NR4A1", "HSPA1A",
               "HSPA1B", "DNAJB1", "HSPH1", "PPP1R15A"],
    "interferon": ["ISG15", "IFIT1", "IFIT2", "IFIT3", "MX1", "OAS1", "OAS2", "IFI6", "IFI44L", "IFITM1", "STAT1",
                   "RSAD2", "XAF1", "IFI27", "BST2"],
    "pEMT": ["LAMC2", "LAMB3", "LAMA3", "TGFBI", "ITGA5", "PDPN", "MMP10", "PTHLH", "VIM", "TNC", "INHBA", "ITGA6",
             "P4HA2", "CDH13", "MMP1", "SERPINE1", "COL17A1", "LGALS1"],
}
PROGRAMS = {k: [g for g in v if canonical(g) not in CLOCK_LIKE] for k, v in PROGRAMS.items()}
_CC = PROGRAMS["cycling"]
S_GENES, G2M_GENES = _CC[: _CC.index("E2F8") + 1], _CC[_CC.index("E2F8") + 1:]  # Tirosh 2016 S and G2M lists


def control_sets(symbols, log_mean: np.ndarray, rng: np.random.Generator, programs=None, n_bins: int = 24,
                 n_ctrl: int = 5) -> dict:
    """Per program: (program gene columns, control gene columns), controls matched by expression bin."""
    programs = PROGRAMS if programs is None else programs
    sym = pd.Index([canonical(s) for s in symbols])
    ok = np.isfinite(log_mean) & ~sym.isin(list(CLOCK_LIKE))
    edges = np.quantile(log_mean[ok], np.linspace(0, 1, n_bins + 1))
    bins = np.clip(np.searchsorted(edges, log_mean, side="right") - 1, 0, n_bins - 1)
    in_any = sym.isin([canonical(g) for v in programs.values() for g in v])
    pool = {b: np.flatnonzero(ok & ~in_any & (bins == b)) for b in range(n_bins)}
    out = {}
    for name, genes in programs.items():
        cols = np.array([i for i in (sym.get_loc(canonical(g)) if canonical(g) in sym else -1 for g in genes)
                         if isinstance(i, (int, np.integer)) and i >= 0 and ok[i]], int)
        ctrl = np.concatenate([rng.choice(pool[bins[c]], n_ctrl, replace=False) for c in cols]) if len(cols) else np.array([], int)
        out[name] = (cols, ctrl)
    return out


def scores(X, total_umi, sets: dict) -> pd.DataFrame:
    """Program scores of the cells of X (cells x genes counts): mean log1p(CP10k) of the program minus its controls."""
    X = sp.csr_matrix(X)
    scale = 1e4 / np.maximum(np.asarray(total_umi, float), 1.0)
    out = {}
    for name, (cols, ctrl) in sets.items():
        if len(cols) == 0:
            out[name] = np.full(X.shape[0], np.nan)
            continue

        def mean_log(c):
            M = X[:, c].multiply(scale[:, None]).tocsr()
            M.data = np.log1p(M.data)
            return np.asarray(M.mean(1)).ravel()
        out[name] = mean_log(cols) - mean_log(ctrl)
    return pd.DataFrame(out)


# ---- builders -----------------------------------------------------------------------------------------------------------

def dispersion_cells(study: str) -> pd.DataFrame:
    """Per-cell clock counts of every classified cell of the cohort's eligible samples (legacy x4 `extract`)."""
    parts = []
    for sample, patient, kind, ld, chem, _ in S.majority_chemistry(S.units(study)[0]):
        X, _, umi, comp, bc, sym = ld()
        cols, sc = gene_columns(sym, GENES), gene_columns(sym, SUB)
        keep = comp != "unassigned"
        Xk = sp.csr_matrix(X)[np.flatnonzero(keep)]
        C, M = Xk[:, cols].toarray(), Xk[:, np.where(sc >= 0, sc, cols[0])].toarray() * (sc >= 0)
        d = pd.DataFrame(np.asarray(C, np.int32), columns=GENES)
        cd3, cd8, cd4, nk = M[:, 0] + M[:, 1], M[:, 2] + M[:, 3], M[:, 4], M[:, 5] + M[:, 6]
        sub = np.where(cd3 > 0, np.where(cd8 > 0, "CD8", np.where(cd4 > 0, "CD4", "T")), np.where(nk > 0, "NK", ""))
        ck = comp[keep]
        parts.append(d.assign(resection=patient, sample=sample, barcode=np.asarray(bc).astype(str)[keep],
                              compartment=[f"{kind}:{c}" for c in ck], stratum=chem, depth=umi[keep],
                              subtype=np.where(ck == "T_NK", sub, "")))
    return pd.concat(parts, ignore_index=True)


def epithelial_cells(study: str, pseudo: int | None = 0) -> tuple[pd.DataFrame, pd.DataFrame]:
    """(epithelial cells, host pseudobulks): clock and pseudo-clock counts, S / G2M and program scores, ERBB2; host =
    pseudobulks of the other compartments (>= 50 cells). Legacy s2 `collect` (S, G2M, erbb2, host) and s1_states
    `collect` (programs, pseudo-clock set `pseudo`)."""
    us = S.majority_chemistry(S.units(study)[0])
    sums, n = 0.0, 0
    for _, _, _, ld, _, _ in us[:8]:  # pass 1: expression of epithelial cells (matches the controls)
        X, _, umi, comp, _, sym = ld()
        e = comp == "epithelial"
        if e.any():
            M = sp.csr_matrix(X)[e].multiply(1e4 / np.maximum(umi[e], 1)[:, None]).tocsr()
            M.data = np.log1p(M.data)
            sums, n = sums + np.asarray(M.sum(0)).ravel(), n + int(e.sum())
    log_mean = sums / max(n, 1)
    sets_sg = control_sets(sym, log_mean, np.random.default_rng(0), {"S": S_GENES, "G2M": G2M_GENES})
    sets_all = control_sets(sym, log_mean, np.random.default_rng(0))
    erbb2 = gene_columns(sym, ["ERBB2"])[0]
    cols = gene_columns(sym, GENES)
    pcols = None
    if pseudo is not None:
        lm = pd.Series(log_mean, index=[str(x).upper() for x in sym])
        lm = lm[~lm.index.duplicated()]
        names = pseudo_waveforms(lm, pseudo + 1, np.random.default_rng(1000), cs.gtex())[pseudo].genes
        pcols = gene_columns(sym, names)
    host, cells = [], []
    for _, pat, kind, ld, chem, _ in us:  # pass 2
        X, _, umi, comp, _, sym = ld()
        Xs = sp.csr_matrix(X)
        C = Xs[:, cols].toarray().astype(float)
        lab = np.array([f"{kind}:{c}" for c in comp], dtype=object)
        h = pd.DataFrame(C, columns=GENES).assign(resection=pat, compartment=lab, stratum=chem, depth=umi, n_cells=1)
        host.append(h[~np.isin(comp, ["epithelial", "unassigned"])].groupby(
            ["resection", "compartment", "stratum"], as_index=False).sum(numeric_only=True))
        e = comp == "epithelial"
        if e.any():
            Xe = Xs[e]
            sg = scores(Xe, umi[e], sets_sg)
            pr = scores(Xe, umi[e], sets_all).rename(columns={"cycling": "cycling"})
            er = np.log1p(np.asarray(Xe[:, erbb2].todense()).ravel() * 1e4 / umi[e])
            df = pd.concat([pd.DataFrame(C[e], columns=GENES), sg, pr.drop(columns=[]).add_prefix("prog_")], axis=1)
            if pcols is not None:
                df = pd.concat([df, pd.DataFrame(Xs[e][:, pcols].toarray().astype(float),
                                                 columns=[f"p0_{g}" for g in GENES])], axis=1)
            cells.append(df.assign(resection=pat, kind=kind, stratum=chem, depth=umi[e], erbb2=er))
    host = pd.concat(host).groupby(["resection", "compartment", "stratum"], as_index=False).sum(numeric_only=True)
    return pd.concat(cells, ignore_index=True), host[host.n_cells >= MIN_HOST].reset_index(drop=True)


def build(study: str, tables=("dispersion", "epithelial")) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    if "dispersion" in tables:
        dispersion_cells(study).to_csv(OUT / f"{study}_dispersion.csv.gz", index=False)
    if "epithelial" in tables:
        c, h = epithelial_cells(study)
        c.to_csv(OUT / f"{study}_epithelial.csv.gz", index=False)
        h.to_csv(OUT / f"{study}_host.csv.gz", index=False)


def load(study: str, table: str = "dispersion", legacy: bool = True) -> pd.DataFrame:
    """Stored table (`dispersion`, `epithelial`, `host`); the legacy run's file when `legacy` and no stored one."""
    dt = dict(resection=str, stratum=str, sample=str, barcode=str)
    f = OUT / f"{study}_{table}.csv.gz"
    if not f.exists() and legacy:
        f = {"dispersion": LEGACY / "celldisp" / f"{study}_cells.csv.gz",
             "epithelial": LEGACY / "states" / "cells" / f"{study}.csv.gz",
             "host": LEGACY / "states" / "cells" / f"{study}_host.csv.gz"}[table]
    if not f.exists():
        build(study, ("dispersion",) if table == "dispersion" else ("epithelial",))
        f = OUT / f"{study}_{table}.csv.gz"
    return pd.read_csv(f, dtype=dt, keep_default_na=False)


def check(study: str) -> None:
    """Rebuild one cohort and compare with the legacy files."""
    new = dispersion_cells(study)
    old = pd.read_csv(LEGACY / "celldisp" / f"{study}_cells.csv.gz", dtype=dict(resection=str, sample=str, barcode=str),
                      keep_default_na=False)
    print("dispersion: rows", len(new), len(old), "clock counts equal:",
          bool((new[GENES].to_numpy() == old[GENES].to_numpy()).all()) if len(new) == len(old) else "n/a",
          "| barcodes equal:", (bool((new.barcode.to_numpy() == old.barcode.to_numpy()).all()) if len(new) == len(old)
                              else "n/a") if "barcode" in old else "(legacy table has no barcodes)")
    c, h = epithelial_cells(study)
    o = pd.read_csv(LEGACY / "states" / "cells" / f"{study}.csv.gz", dtype=dict(resection=str, stratum=str), keep_default_na=False)
    oh = pd.read_csv(LEGACY / "states" / "cells" / f"{study}_host.csv.gz", dtype=dict(resection=str, stratum=str), keep_default_na=False)
    print("epithelial: rows", len(c), len(o), "max |diff| S, G2M, erbb2:",
          [float(np.abs(c[k].to_numpy() - o[k].to_numpy()).max()) if len(c) == len(o) else "n/a" for k in ("S", "G2M", "erbb2")])
    print("host: rows", len(h), len(oh), "counts equal:", bool((h[GENES].to_numpy() == oh[GENES].to_numpy()).all())
          if len(h) == len(oh) else "n/a")


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    b = sub.add_parser("build")
    b.add_argument("study")
    b.add_argument("--tables", nargs="+", default=["dispersion", "epithelial"])
    sub.add_parser("check").add_argument("study")
    a = ap.parse_args()
    build(a.study, a.tables) if a.cmd == "build" else check(a.study)
