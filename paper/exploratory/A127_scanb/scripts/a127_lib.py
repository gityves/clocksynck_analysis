"""A127 shared loader and definitions (A97 definitions applied to the SCAN-B log2(FPKM + 0.1) matrix)."""
import sys
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "paper"))
sys.path.insert(0, str(ROOT / "paper/results/figure_5/scripts/cancer"))
import _shared as sh  # noqa: E402  (PURITY, hallmark_g2m, _zmean: A97 definitions)
from clocksync.decompose import components  # noqa: E402
from clocksync.waveform import Waveform  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "output"
GENES, WF = sh.GENES, sh.WF
PL = sh.plane()
SUBTYPES = ["LumA", "LumB", "Her2", "Basal", "Normal"]


def tertile(x, k=3):
    return pd.qcut(pd.Series(x).rank(method="first"), k, labels=False).to_numpy()


def load():
    """(tumours: one row per patient with PAM50 and ER recorded, expression: samples x genes float32), scores z-scored within
    this cohort. Replicates dropped (first listed kept)."""
    meta = pd.read_csv(OUT / "scanb_samples_parsed.csv", index_col=0)
    E = pd.read_parquet(OUT / "scanb_genes_kept.parquet")
    m = meta[~meta.is_repl.astype(bool)].sort_values("col_order").drop_duplicates("patient")
    m = m[m["pam50 subtype"].isin(SUBTYPES) & m["er status"].isin([0, 1, "0", "1"])].copy()
    m["er"] = m["er status"].astype(int)
    X = E.loc[m.index]
    g2m = [g for g in sh.hallmark_g2m() if g in X.columns]
    m["prolif"] = sh._zmean(X, g2m).to_numpy()
    m["purity"] = sh._zmean(X, sh.PURITY).to_numpy()
    m["pam50"] = m["pam50 subtype"]
    m["platform_i"] = (m.platform == "NextSeq").astype(int)
    return m, X


def wf_for(genes) -> Waveform:
    """Waveform with the clock's amplitudes / phases on `genes` (the clock genes themselves or a pseudo set)."""
    return Waveform(tuple(genes), WF.amp, WF.psi, None, None, "set")


def coh(Y, idx, wf) -> float:
    """A97 coherence (decompose.components): Y = ln2 * log2 FPKM of the set's genes (samples x 14); no Poisson term (FPKM)."""
    y = Y[idx] - Y[idx].mean(0)
    return components(y, np.zeros(Y.shape[1]), 1, wf)["coherence"]
