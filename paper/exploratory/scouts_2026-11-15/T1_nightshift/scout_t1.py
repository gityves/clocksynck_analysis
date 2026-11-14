"""Scout for roadmap item T1 (night-shift muscle, genome-wide): can programs of genes be assigned to the time they follow?
Simulation on the real design; no statistic of the later amendment is computed on the real labels.

Design (GSE106800, Wefers 2018): vastus lateralis at 07:00 and 19:00 after an aligned and after a 12 h-misaligned
protocol, 11 subjects with all four arrays. For a gene with reference amplitude A and acrophase phi (clock time), the
evening-minus-morning difference is 2 A cos(w (19 - phi - s)) under a shift s, times an excursion ratio r. Hypotheses:
  central time  s = 0 h,    r = 1      (the gene keeps its aligned difference)
  tissue clock  s = 2.68 h, r = 0.233  (the muscle clock, A257 point estimates)
  behaviour     s = 12 h,   r = 1      (the difference reverses)
Estimator for a program of G genes: regress the subject-mean differences of each condition on (c_g, s_g) = 2 A_g (cos, sin)
(w (19 - phi_g)), giving complex z_aligned and z_misaligned; shift = arg(z_mis / z_al) / w, ratio = |z_mis / z_al|;
90 % jackknife-t intervals over subjects (log ratio and angle). A hypothesis
(s, r) is excluded when the Hotelling T2 of the per-subject z_mis - r e^{iws} z_al rejects at 0.10 (F(2, n - 2)). The amplitude scale between reference and arrays cancels in the ratio.

Inputs used: reference phases and amplitudes from the timed muscle cohort (GSE108539, 10 subjects x 6 biopsies,
per-gene cosinor with subject intercepts, exon layer); residual vectors of the GSE106800 differences around their
condition means (design noise with the real gene-gene correlation; no hypothesis contrast is computed on them).
Truth per gene = reference estimate + its sampling error + a between-cohort phase difference (sd 1 h).

    python scout_t1.py [--reps 200] [--jobs 4]
Outputs (output/): reference_genes.csv (aggregate per-gene cosinor of the public reference), noise_summary.csv,
runs.csv, separability.csv.
"""

from __future__ import annotations

import argparse
import sys
from multiprocessing import Pool
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[2]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(PAPER / "exploratory/A257_nightshift_muscle/scripts"))
from clocksync.config import limit_threads  # noqa: E402
from common.paths import DATA  # noqa: E402

OUT = HERE / "output"
W = 2 * np.pi / 24
SEED = 1115
HYP = {"central": (0.0, 1.0), "tissue_clock": (2.68, 0.233), "behaviour": (12.0, 1.0)}
SIZES = [5, 10, 20, 50, 100]
PHASE_SD_COHORT = 1.0


def reference() -> pd.DataFrame:
    """Per-gene cosinor (log2 RPKM, exon layer) with subject intercepts on GSE108539."""
    from scipy import stats
    import common.datasets.muscle_timed_gse108539.build as B
    x = pd.read_csv(DATA / "muscle_timed/GSE108539/GSE108539_norm_reads_count.txt.gz", sep="\t", index_col=0)
    x.index = x.index.str.split(".").str[0]
    x = x.loc[:, [c.startswith("Exon") for c in x.columns]]
    sym = B._symbols().reindex(x.index)
    x = x[sym.notna().to_numpy()].groupby(sym.dropna().to_numpy()).mean()
    meta = x.columns.str.extract(r"Exon_(\d+)_S\d+_(\d+)\.\d+$")
    subj, t = meta[0].to_numpy(), meta[1].astype(float).to_numpy()
    keep = (x.mean(axis=1) > 0).to_numpy()  # log2 RPKM > 0
    Y = x[keep].to_numpy().T
    Dm = pd.get_dummies(subj).to_numpy(float)
    D = np.c_[Dm, np.cos(W * t), np.sin(W * t)]
    B_, *_ = np.linalg.lstsq(D, Y, rcond=None)
    res = Y - D @ B_
    dfr = len(Y) - D.shape[1]
    s2 = (res ** 2).sum(0) / dfr
    XtXi = np.linalg.pinv(D.T @ D)
    a, b = B_[-2], B_[-1]
    rss0 = ((Y - Dm @ np.linalg.lstsq(Dm, Y, rcond=None)[0]) ** 2).sum(0)
    F = ((rss0 - (res ** 2).sum(0)) / 2) / s2
    p = stats.f.sf(F, 2, dfr)
    se = np.sqrt(s2 * (XtXi[-2, -2] + XtXi[-1, -1]) / 2)
    return pd.DataFrame(dict(gene=x.index[keep], amp=np.hypot(a, b), acro_h=np.mod(np.arctan2(b, a), 2 * np.pi) / W,
                             se=se, p=p))


def arrays():
    """Subject-level 19:00 minus 07:00 differences (log2) per condition: dict cond -> DataFrame subjects x genes."""
    from compute_a257 import load
    smp, y = load()
    y = y / np.log(2.0)  # back to log2
    d = {}
    for c in ("ctl", "mis"):
        rows = {}
        for s, g in smp[smp.cond == c].groupby("subjectid"):
            if set(g.pm) == {0, 1}:
                rows[s] = y.iloc[g.index[g.pm == 1][0]] - y.iloc[g.index[g.pm == 0][0]]
        d[c] = pd.DataFrame(rows).T
    both = sorted(set(d["ctl"].index) & set(d["mis"].index))
    return {c: d[c].loc[both] for c in d}


def fit_z(D, C, S):
    """Complex program loading from subject-mean differences D (genes) on columns C, S (genes)."""
    X = np.c_[C, S]
    b, *_ = np.linalg.lstsq(X, D, rcond=None)
    return b[0] + 1j * b[1]


def one(task):
    rep, G, kind, hyp, ref, res_al, res_mis = task
    limit_threads(1)
    rng = np.random.default_rng([SEED, rep, G, len(kind), len(hyp)])
    if kind == "random":
        idx = rng.choice(len(ref), G, replace=False)
    else:  # phase-clustered: genes within +-3 h of a random centre
        c = rng.uniform(0, 24)
        pool = np.flatnonzero(np.abs((ref.acro_h.to_numpy() - c + 12) % 24 - 12) <= 3)
        idx = rng.choice(pool, min(G, len(pool)), replace=False)
    r_ = ref.iloc[idx]
    s_true, r_true = HYP[hyp]
    phi_t = r_.acro_h.to_numpy() + rng.normal(0, 1, len(r_)) * r_.se.to_numpy() / r_.amp.to_numpy() / W \
        + rng.normal(0, PHASE_SD_COHORT, len(r_))
    a_t = np.maximum(r_.amp.to_numpy() + rng.normal(0, 1, len(r_)) * r_.se.to_numpy(), 0)
    sig_al = 2 * a_t * np.cos(W * (19 - phi_t))
    sig_mis = r_true * 2 * a_t * np.cos(W * (19 - phi_t - s_true))
    n = res_al.shape[0]
    gcols = r_.col.to_numpy()
    pick = rng.integers(0, n, n)  # subjects drawn jointly for both conditions (within-subject correlation kept)
    da = sig_al[None, :] + res_al[pick][:, gcols]
    dm = sig_mis[None, :] + res_mis[pick][:, gcols]
    C = 2 * r_.amp.to_numpy() * np.cos(W * (19 - r_.acro_h.to_numpy()))
    S = 2 * r_.amp.to_numpy() * np.sin(W * (19 - r_.acro_h.to_numpy()))

    def est(ia):
        za, zm = fit_z(da[ia].mean(0), C, S), fit_z(dm[ia].mean(0), C, S)
        return np.log(zm / za)  # log ratio + i * angle

    from scipy import stats
    q0 = est(np.arange(n))
    jk = np.array([est(np.delete(np.arange(n), i)) for i in range(n)])
    jk = jk.real + 1j * (q0.imag + np.angle(np.exp(1j * (jk.imag - q0.imag))))
    se = lambda v: np.sqrt((n - 1) / n * ((v - v.mean()) ** 2).sum())
    t = stats.t.ppf(0.95, n - 1)
    s_hat, r_hat = q0.imag / W, float(np.exp(q0.real))
    half_s, half_r = t * se(jk.imag) / W, t * se(jk.real)
    lo_s, hi_s = s_hat - half_s, s_hat + half_s
    lo_r, hi_r = r_hat * np.exp(-half_r), r_hat * np.exp(half_r)
    inside = lambda h: abs((h - s_hat + 12) % 24 - 12) <= half_s
    out = dict(rep=rep, G=len(r_), kind=kind, truth=hyp, shift_h=s_hat, ratio=r_hat, lo_s=lo_s, hi_s=hi_s,
               lo_r=lo_r, hi_r=hi_r)
    zA = np.array([fit_z(da[i], C, S) for i in range(n)])  # per-subject program loadings
    zM = np.array([fit_z(dm[i], C, S) for i in range(n)])
    for h, (s, r) in HYP.items():
        # Hotelling T2 of the per-subject residual z_mis - r e^{iws} z_al (mean 0 under h, whatever the amplitude scale)
        dd = zM - r * np.exp(1j * W * s) * zA
        X = np.c_[dd.real, dd.imag]
        m, Sg = X.mean(0), np.cov(X.T)
        T2 = n * m @ np.linalg.solve(Sg, m)
        p_h = stats.f.sf((n - 2) / (2 * (n - 1)) * T2, 2, n - 2)
        out[f"p_{h}"] = p_h
        out[f"excl_{h}"] = bool(p_h < 0.10)
        out[f"excl_shift_{h}"] = bool(not inside(s))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=200)
    ap.add_argument("--jobs", type=int, default=4)
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    ref = reference()
    ref.to_csv(OUT / "reference_genes.csv", index=False)
    d = arrays()
    genes = sorted(set(ref.gene) & set(d["ctl"].columns))
    rhy = ref[ref.gene.isin(genes) & (ref.p < 0.01)].reset_index(drop=True)
    al, mis = d["ctl"][rhy.gene], d["mis"][rhy.gene]
    res_al, res_mis = (al - al.mean()).to_numpy(), (mis - mis.mean()).to_numpy()
    rhy["col"] = np.arange(len(rhy))
    pd.DataFrame([dict(subjects=len(al), reference_genes=len(ref), reference_rhythmic_p01=int((ref.p < 0.01).sum()),
                       on_array=len(rhy), median_amp_log2=rhy.amp.median(),
                       median_resid_sd_aligned=float(np.median(res_al.std(0, ddof=1))),
                       median_resid_sd_misaligned=float(np.median(res_mis.std(0, ddof=1))),
                       median_signal_2A=float(np.median(2 * rhy.amp)))]).to_csv(OUT / "noise_summary.csv", index=False)
    tasks = [(r, G, k, h, rhy, res_al, res_mis) for G in SIZES for k in ("random", "clustered") for h in HYP
             for r in range(a.reps)]
    with Pool(a.jobs) as pool:
        R = pd.DataFrame(pool.map(one, tasks, chunksize=20))
    R.to_csv(OUT / "runs.csv", index=False)
    S = R.groupby(["kind", "G", "truth"]).agg(
        shift_median=("shift_h", "median"), ratio_median=("ratio", "median"),
        **{f"excl_{h}": (f"excl_{h}", "mean") for h in HYP},
        **{f"excl_shift_{h}": (f"excl_shift_{h}", "mean") for h in HYP}).reset_index()
    S.to_csv(OUT / "separability.csv", index=False)
    pd.set_option("display.width", 250)
    print(pd.read_csv(OUT / "noise_summary.csv").T, "\n", S.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
