"""Post hoc, descriptive (2026-11-03; user observation stated before running): is the time-of-death signal of GTEx Whole Blood
carried by blood cell composition rather than by a coherent core clock?

GTEx v10 fast deaths (Hardy 1-2) with recorded time of death, Whole Blood and five solid tissues for comparison. Per tissue, on
log2-normalised counts of expressed genes (mean log2 >= 3):
- composition scores: mean within-tissue z of a priori marker genes (neutrophil, T, B, NK, monocyte, eosinophil, erythroid);
- per-gene cosinor R^2 on time of death, and the partial R^2 of the time terms after the composition scores;
- top-100 time genes (clock genes excluded) chosen in one random donor half and evaluated in the other (both directions,
  averaged), against 100 random expressed genes;
- the 14 clock genes: median R^2 before / after composition, and phase concordance with the GTEx waveform
  |sum_g w_g exp(i(phase_g - psi_g))| / sum_g w_g, w = cosinor amplitude (1: genes keep their expected relative phases).
Controlled access: only per-tissue aggregates are written; no donor time or value leaves memory.

Extension (post hoc, 2026-11-25, user request for slide 3.6): the same readouts on the other low-gain tissues, testis and the
13 brain regions (`--low-gain`, written to `blood_composition_low_gain.csv`; the six tissues above are not rerun). The
composition scores stay the blood-cell markers (blood content of the tissue).

    python blood_composition.py              # a few minutes
    python blood_composition.py --low-gain
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
sys.path[:0] = [str(PAPER), str(PAPER / "results/figure_1/scripts/method")]
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from _shared import load  # noqa: E402
from clocksync import paired as P  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common.benchmark.data import log2_normalised  # noqa: E402

OUT = HERE.parents[1] / "output"
TISSUES = ("Whole Blood", "Muscle - Skeletal", "Lung", "Thyroid", "Skin - Sun Exposed (Lower leg)", "Adipose - Subcutaneous")
MARKERS = {"neutrophil": ["FCGR3B", "CSF3R", "CXCR2", "FPR1", "MNDA", "S100A12", "CMTM2", "MMP25"],
           "T": ["CD3E", "CD3D", "CD3G", "CD2", "CD247", "IL7R", "LCK", "CD5"],
           "B": ["MS4A1", "CD79A", "CD79B", "CD19", "BANK1", "CD22"],
           "NK": ["GNLY", "NKG7", "KLRD1", "KLRF1", "PRF1", "GZMB"],
           "monocyte": ["CD14", "FCN1", "VCAN", "CSF1R", "LYZ"],
           "eosinophil": ["CLC", "IL5RA", "CCR3", "SIGLEC8"],
           "erythroid": ["HBB", "HBA1", "HBA2", "ALAS2", "CA1"]}
LOW_GAIN = ("Testis", "Brain - Amygdala", "Brain - Anterior cingulate cortex (BA24)", "Brain - Caudate (basal ganglia)",
            "Brain - Cerebellar Hemisphere", "Brain - Cerebellum", "Brain - Cortex", "Brain - Frontal Cortex (BA9)",
            "Brain - Hippocampus", "Brain - Hypothalamus", "Brain - Nucleus accumbens (basal ganglia)",
            "Brain - Putamen (basal ganglia)", "Brain - Spinal cord (cervical c-1)", "Brain - Substantia nigra")
WF = P.gtex()
CLOCK = list(WF.genes)
N_TOP, SEED = 100, 20261103


def z(x):
    return (x - x.mean(0)) / x.std(0)


def resid(Y, C):
    C = np.column_stack([np.ones(len(Y)), C])
    return Y - C @ np.linalg.lstsq(C, Y, rcond=None)[0]


def cosinor(Y, h, C=None):
    """Per-gene R^2 of the time terms (partial, after C if given), amplitude and acrophase (h)."""
    w = 2 * np.pi * h / 24
    T = np.column_stack([np.cos(w), np.sin(w)])
    base = np.zeros((len(Y), 0)) if C is None else C
    r0 = resid(Y, base)
    full = np.column_stack([base, T])
    X = np.column_stack([np.ones(len(Y)), full])
    beta = np.linalg.lstsq(X, Y, rcond=None)[0]
    r1 = Y - X @ beta
    r2 = 1 - (r1 ** 2).sum(0) / (r0 ** 2).sum(0)
    a, b = beta[-2], beta[-1]
    return r2, np.hypot(a, b), (np.arctan2(b, a) * 24 / (2 * np.pi)) % 24


def concordance(amp, phase_h, genes):
    psi = pd.Series(WF.psi, index=CLOCK).loc[genes].to_numpy()
    ph = 2 * np.pi * phase_h / 24
    return float(abs(np.sum(amp * np.exp(1j * (ph - psi)))) / amp.sum())


def tissue(t, counts, tod, rng):
    tod = tod.dropna()
    L = log2_normalised(counts.loc[tod.index])
    L = L.loc[:, L.mean() >= 3]
    h = tod.to_numpy(float)
    comp = pd.DataFrame({k: z(L[[g for g in v if g in L.columns]]).mean(1) for k, v in MARKERS.items()
                         if sum(g in L.columns for g in v) >= 2})
    C = comp.to_numpy()
    rows = dict(tissue=t, donors=len(L), genes=L.shape[1], composition_scores=";".join(comp.columns))
    # composition scores' own rhythm
    cr2, _, cph = cosinor(z(C), h)
    for k, r, p in zip(comp.columns, cr2, cph):
        rows[f"R2_{k}"], rows[f"peak_h_{k}"] = float(r), float(p)
    # top time genes: select in one half, evaluate in the other
    pool = [g for g in L.columns if g not in CLOCK]
    Y = L[pool].to_numpy(float)
    half = rng.permutation(len(L)) < len(L) // 2
    ev = []
    for sel, tst in ((half, ~half), (~half, half)):
        r2s, _, _ = cosinor(Y[sel], h[sel])
        top = np.argsort(r2s)[::-1][:N_TOP]
        rnd = rng.choice(Y.shape[1], N_TOP, replace=False)
        r2t, _, _ = cosinor(Y[tst][:, top], h[tst])
        r2c, _, _ = cosinor(Y[tst][:, top], h[tst], C[tst])
        r2r, _, _ = cosinor(Y[tst][:, rnd], h[tst])
        cc = np.corrcoef(np.column_stack([Y[tst][:, top], C[tst]]).T)[:N_TOP, N_TOP:]
        ev.append([np.median(r2t), np.median(r2c), np.median(r2r), np.median(np.abs(cc).max(1)),
                   np.mean(np.abs(cc).max(1) > 0.5)])
    ev = np.mean(ev, 0)
    rows.update(top100_R2=ev[0], top100_R2_after_composition=ev[1], random100_R2=ev[2],
                top100_share_removed=1 - ev[1] / ev[0], top100_median_max_abs_r_composition=ev[3],
                top100_share_abs_r_gt_0_5=ev[4])
    # clock genes
    cg = [g for g in CLOCK if g in L.columns]
    Yc = L[cg].to_numpy(float)
    r2, amp, ph = cosinor(Yc, h)
    r2a, ampa, pha = cosinor(Yc, h, C)
    allr2, _, _ = cosinor(L.to_numpy(float), h)
    rows.update(clock_genes=len(cg), clock_R2_median=float(np.median(r2)), clock_R2_after_composition=float(np.median(r2a)),
                clock_R2_percentile_median=float(np.median([(allr2 < x).mean() * 100 for x in r2])),
                clock_phase_concordance=concordance(amp, ph, cg), clock_phase_concordance_after_composition=concordance(ampa, pha, cg))
    return rows


def main():
    limit_threads()
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    low = "--low-gain" in sys.argv
    tissues, name = (LOW_GAIN, "blood_composition_low_gain.csv") if low else (TISSUES, "blood_composition.csv")
    data = load(list(tissues))
    res = pd.DataFrame([tissue(t, *data[t], rng) for t in tissues if t in data])
    res.to_csv(OUT / name, index=False)
    pd.set_option("display.width", 250)
    print(res.T.to_string())


if __name__ == "__main__":
    main()
