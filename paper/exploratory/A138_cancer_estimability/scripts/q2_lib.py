"""A138 Q2: simulator (p = 2,000 genes, TCGA BRCA non-clock structure, GTEx clock waveforms) and the ordering methods.

Registered in paper/docs/registered_experiments.md (A138, Q2). Choices not in the registration are listed in
`q2_gate.py` ("Choices not in the registration").

Simulator. Panel: the 15 CoFE clock genes, the CYCLOPS default seed genes present in the BRCA table, 15 % of p rhythmic
followers (GTEx: p < 0.01 in >= 3 tissues; GTEx pooled amplitude and phase) and random other BRCA genes. Log2 expression
of sample j: BRCA gene mean + the leading K non-clock principal components of the 782 BRCA tumours (score rows of
random tumours, so the variances and the subtype / purity structure are the real ones) + independent Gaussian residual
with the BRCA residual SD of each gene + the clock term amp_g cos(w (t_j + e_j - acrophase_g)) (log2 units, GTEx).
"""
from __future__ import annotations

import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[2]
sys.path.insert(0, str(PAPER))

from common.benchmark import evaluate as ev  # noqa: E402
from common.benchmark.tcga import COFE_CLOCK  # noqa: E402

OUT = HERE.parent / "output" / "q2_gate"
GTEX_COS = PAPER / "results" / "figure_2" / "output" / "gene_rhythms" / "gtex_cosinor.csv.gz"
P_GENES, K_PCS, FOLLOW = 2000, 10, 0.15
W = 2 * np.pi / 24
LN2 = np.log(2.0)
CLOCK = list(COFE_CLOCK)
N_MODELS = 80          # CYCLOPS template
SEED_MTH_GENE = 1500   # CYCLOPS expression-rank cutoff: 75 % of p (template: 10,000 of the ~13,000 TCGA genes)
PANEL_SEED = 138


def brca_table() -> pd.DataFrame:
    """BRCA tumours as COFE selects them, log2 normalised (cached in /tmp for repeated builds)."""
    f = Path(tempfile.gettempdir()) / "brca_a138.pkl"
    if f.exists():
        return pd.read_pickle(f)
    from common.benchmark.data import cofe_filter, tcga_counts, tcga_tumours
    X = cofe_filter(tcga_counts("BRCA", tcga_tumours("BRCA")))
    X.to_pickle(f)
    return X


def build_panel() -> dict:
    import clocksync as cs
    from common.benchmark.cyclops import seeds as cyclops_seeds
    X = brca_table()
    rng = np.random.default_rng(PANEL_SEED)
    g = pd.read_csv(GTEX_COS)
    sig = g[g.p < 0.01]
    z = sig.assign(c=sig.b_cos + 1j * sig.b_sin).groupby("gene").c.mean()
    amp_all = sig.groupby("gene").amp.mean()
    ntis = sig.groupby("gene").tissue.nunique()
    # waveform of the 15 clock genes: clocksync's GTEx waveform where it has the gene, pooled GTEx cosinor otherwise
    w0 = cs.gtex()
    amp, psi = {}, {}
    gc = g.assign(c=g.b_cos + 1j * g.b_sin).groupby("gene").agg(c=("c", "mean"), amp=("amp", "mean"))
    for x in CLOCK:
        if x in w0.genes:
            j = w0.genes.index(x)
            amp[x], psi[x] = float(w0.amp[j]), float(w0.psi[j])
        else:
            amp[x], psi[x] = float(gc.amp[x]), float(np.angle(gc.c[x]) % (2 * np.pi))
    n_f = int(round(FOLLOW * P_GENES))
    cand = [x for x in ntis.index[ntis >= 3] if x in X.columns and x not in CLOCK]
    if len(cand) < n_f:
        cand = [x for x in ntis.index[ntis >= 2] if x in X.columns and x not in CLOCK]
    fol = sorted(rng.choice(sorted(cand), n_f, replace=False))
    for x in fol:
        amp[x], psi[x] = float(amp_all[x]), float(np.angle(z[x]) % (2 * np.pi))
    seeds = [x for x in cyclops_seeds() if x in X.columns and x not in CLOCK and x not in fol]
    pool = sorted(set(X.columns) - set(CLOCK) - set(fol) - set(seeds))
    others = list(seeds) + list(rng.choice(pool, P_GENES - len(CLOCK) - n_f - len(seeds), replace=False))
    panel = CLOCK + list(fol) + others
    Y = X[panel]
    mesor = Y.mean().to_numpy()
    Yc = (Y - Y.mean()).to_numpy()
    U, s, Vt = np.linalg.svd(Yc, full_matrices=False)
    scores = U[:, :K_PCS] * s[:K_PCS]
    V = Vt[:K_PCS].T
    resid_sd = (Yc - scores @ V.T).std(0)
    lam = s ** 2 / (len(Yc) - 1)
    a = np.array([amp.get(x, 0.0) for x in panel]) / LN2  # natural log -> log2 units
    ps = np.array([psi.get(x, 0.0) for x in panel])
    OUT.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"gene": panel, "amp_log2": a, "acro_h": ps / W, "mesor": mesor, "resid_sd": resid_sd,
                  "kind": ["clock"] * len(CLOCK) + ["follower"] * n_f + ["other"] * len(others)}).to_csv(
        OUT / "panel_genes.csv", index=False)
    pd.DataFrame({"pc": np.arange(1, K_PCS + 1), "lambda": lam[:K_PCS], "share": lam[:K_PCS] / lam.sum()}).to_csv(
        OUT / "panel_pcs.csv", index=False)
    return dict(genes=np.array(panel, dtype=str), amp=a, psi=ps, mesor=mesor, scores=scores, V=V, resid_sd=resid_sd,
                n_seeds_in=len(seeds))


def load_panel() -> dict:
    f = OUT / "panel.npz"
    if not f.exists():
        P = build_panel()
        np.savez(f, **P)
    z = np.load(f, allow_pickle=False)
    P = {k: z[k] for k in z.files}
    P["genes"] = [str(x) for x in P["genes"]]
    return P


def simulate(P: dict, n: int, seed: int, arm: str = "full", sd_h: float = 0.0, window=(8.0, 18.0)):
    """Samples x genes log2 expression and true clock time (h). Arms: full (24 h), window, null (no clock)."""
    rng = np.random.default_rng(seed)
    t = rng.uniform(*(window if arm == "window" else (0.0, 24.0)), n)
    idx = rng.choice(len(P["scores"]), n, replace=False)
    amp = np.zeros_like(P["amp"]) if arm == "null" else P["amp"]
    phi = W * (t + rng.normal(0, sd_h, n))[:, None]  # clock state of the tumour (SD 0: its own sampling time)
    Y = (P["mesor"] + P["scores"][idx] @ P["V"].T + rng.normal(size=(n, len(P["genes"]))) * P["resid_sd"]
         + amp * np.cos(phi - P["psi"]))
    return pd.DataFrame(Y, columns=P["genes"], index=[f"s{i}" for i in range(n)]), t


def run_method(method: str, Y: pd.DataFrame, seed: int, tmp: Path, threads: int = 1):
    """Phases (h, arbitrary origin and direction) and a metadata dict."""
    if method in ("cofe_genome", "cofe_clock"):
        from common.benchmark.cofe import run_cofe
        X = Y if method == "cofe_genome" else Y[CLOCK]
        r = run_cofe(X, seed=seed, ncores=threads)
        return r["phase"].rename(None), dict(best_s=r["best_s"], interior_min=bool(r["interior_min"]), n_genes=X.shape[1])
    if method == "cyclops":
        from common.benchmark.cyclops import run_cyclops, seeds
        sg = [x for x in seeds() if x in Y.columns]
        expr = (2 ** Y - 1).T
        r = run_cyclops(expr, sg, tmp, seed=seed, threads=threads, seed_mth_Gene=SEED_MTH_GENE, train_n_models=N_MODELS)
        return (r["phase_raw"] * 24 / (2 * np.pi)) % 24, dict(n_seeds=len(sg))
    raise ValueError(method)


def score(phase: pd.Series, Y: pd.DataFrame, t: np.ndarray, P: dict) -> tuple[dict, pd.DataFrame]:
    """Readouts as the GTEx benchmark (`compute_gtex_benchmark.score_run`): ordering aligned to the truth by the best
    rotation / reflection (oracle), clock-gene acrophases on the aligned ordering against the cosinor of the true times
    (`arrangement_error`, benchmark metric) and against the true waveform acrophases."""
    p = phase.to_numpy(float)
    aligned, mae, _ = ev.align(p, t)
    est = ev.cosinor(Y[CLOCK], aligned, CLOCK)
    truth = ev.cosinor(Y[CLOCK], t, CLOCK)
    wf = pd.DataFrame({"amplitude": P["amp"][: len(CLOCK)], "phase": (P["psi"][: len(CLOCK)] / W) % 24}, index=CLOCK)
    a1 = ev.arrangement_error(est, truth)
    a2 = ev.arrangement_error(est, wf)
    genes = a1.rename(columns={"phase_err": "err_cosinor_h"}).join(a2.phase_err.rename("err_waveform_h"))
    m = dict(mae_h=mae, circ_corr=abs(ev.circ_corr(p, t)), skill=ev.skill(mae, t),
             clock_err_med_h=float(a1.phase_err.abs().median()), clock_err_wf_med_h=float(a2.phase_err.abs().median()),
             clock_err_mean_h=float(a1.phase_err.abs().mean()))
    return m, genes
