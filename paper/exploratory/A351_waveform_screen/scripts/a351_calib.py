"""A351 calibration (V351, stop rule): can the gene-level screen (Q_a, Q_phi) separate uniform damping (D) from shared (W) and
per-tumour (H) rewiring on the A212 plasmode?  Bases, generators, matched strengths (A212 output/matched.csv) and seeds are
A212's own (compute_a212.py is exec'd up to its main block; replicate k, r of the 'clouds' stage uses the same SeedSequence, so
the same case halves and mechanism draws are reproduced).

Estimator (group design, per replicate and readout = 14 genes or core loop = a336_common.core):
  each stratum's case and control half is centred; host score s = y X_c (X_c: gene-centred waveform design of the readout's
  waveform, the same fixed projector for both halves, so uniform damping rescales amplitudes but not loadings);
  per gene g: OLS of y_g on s over the stacked strata -> complex coefficient c_g (b_cos + i b_sin);
  amplitude A_g = sd of the fitted component; a_g = ln(A_case / A_ctrl); phi_g = arg(c_case / c_ctrl).
  Variances: unit bootstrap (N_BOOT = 50), tumours resampled within stratum x half; var(a_g) over boots, var(phi_g) of the wrapped
  deviation from the point estimate.  Q_a = sum w_g (a_g - weighted mean)^2; Q_phi = sum w_g wrap(phi_g - mu)^2 with mu the
  weighted circular mean; w = 1/var.
AUC: 5-fold cross-validated LDA on (ln Q_a, ln Q_phi) (joint) and the raw statistic with rewiring scoring higher (each alone).
    python a351_calib.py [--jobs 4] [--reps 200]"""
import argparse
import sys
import zlib
from pathlib import Path

HERE = Path(__file__).resolve().parent
EXP = HERE.parents[1]
ap = argparse.ArgumentParser(); ap.add_argument("--jobs", type=int, default=4); ap.add_argument("--reps", type=int, default=200)
A = ap.parse_args()
src = (EXP / "A212_coherence_mechanisms/scripts/compute_a212.py").read_text().split("# ---- main")[0]
sys.argv = ["compute_a212.py", "--stage", "clouds", "--jobs", "1"]
G = {"__file__": str(EXP / "A212_coherence_mechanisms/scripts/compute_a212.py"), "__name__": "a212"}
exec(compile(src, "compute_a212.py", "exec"), G)
sys.path.insert(0, str(EXP / "A336_two_readouts/scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis  # noqa: E402
from sklearn.metrics import roc_auc_score  # noqa: E402
from sklearn.model_selection import StratifiedKFold, cross_val_predict  # noqa: E402

import a336_common as C336  # noqa: E402

OUT = HERE.parent / "output" / "calibration"; OUT.mkdir(parents=True, exist_ok=True)
N_BOOT, SEED, split, apply_, design = 50, G["SEED"], G["split"], G["apply"], G["design"]
M = pd.read_csv(EXP / "A212_coherence_mechanisms/output/matched.csv").set_index(["base", "mechanism"])
say = G["say"]
wrap = lambda x: (x + np.pi) % (2 * np.pi) - np.pi


def halves(parts, wf, mech, s, seed, sd_e):
    """A212 replicate's random halves with the mechanism applied: [(case, control)] per stratum (rng order as A212)."""
    rng = np.random.default_rng(seed); out = []
    if mech == "W":
        Xw = design(wf, rng.normal(0.0, s, len(wf.genes)))
    for C, E, beta in parts:
        i = rng.permutation(len(C)); a, b = i[: len(i) // 2], i[len(i) // 2:]
        ya = beta[a] @ Xw.T + E[a] if mech == "W" else apply_(mech, s, C[a], E[a], beta[a], wf, rng, sd_e)
        out.append((ya, C[b] + E[b]))
    return out


def fit(ys, idx, cols, Xc):
    """ys: per-stratum arrays; idx: per-stratum row indices. -> complex c_g and amplitude A_g over `cols`."""
    y = np.vstack([v[i] - v[i].mean(0) for v, i in zip(ys, idx)])[:, cols]
    s = y @ Xc
    coef = np.linalg.solve(s.T @ s, s.T @ y)
    return coef[0] + 1j * coef[1], np.sqrt(np.mean((s @ coef) ** 2, 0))


def screen(case, ctrl, cols, Xc, rng):
    def est(ic, ib):
        c1, A1 = fit(case, ic, cols, Xc); c0, A0 = fit(ctrl, ib, cols, Xc)
        return np.log(A1 / A0), np.angle(c1 / c0)
    full = lambda ys: [np.arange(len(v)) for v in ys]
    a, p = est(full(case), full(ctrl))
    da, dp = [], []
    for _ in range(N_BOOT):
        ic = [rng.integers(0, len(v), len(v)) for v in case]; ib = [rng.integers(0, len(v), len(v)) for v in ctrl]
        ab, pb = est(ic, ib); da.append(ab - a); dp.append(wrap(pb - p))
    wa, wp = 1 / np.var(da, 0), 1 / np.var(dp, 0)
    qa = np.sum(wa * (a - np.sum(wa * a) / wa.sum()) ** 2)
    mu = np.angle(np.sum(wp * np.exp(1j * p)))
    return qa, np.sum(wp * wrap(p - mu) ** 2)


def replicate(parts, wf, mech, s, seed, sd_e, readouts):
    hv = halves(parts, wf, mech, s, seed, sd_e)
    case, ctrl = [h[0] for h in hv], [h[1] for h in hv]
    rng = np.random.default_rng(np.random.SeedSequence([SEED, 351, *seed.spawn_key[-3:], zlib.crc32(mech.encode())]))
    res = {}
    for name, w in readouts.items():
        cols = [wf.genes.index(g) for g in w.genes]
        Xc = w.X - w.X.mean(0)
        res[name] = screen(case, ctrl, cols, Xc, rng)
    return res


rows = []
for bi, b in enumerate(["B1", "B2"]):
    strata, wf = dict(B1=G["base_b1"], B2=G["base_b2"])[b]()
    parts = [split(Y, wf) for Y in strata]
    sd_e = float(np.median(np.concatenate([E for _, E, _ in parts]).std(0)))
    readouts = {"g14": wf, "core": C336.core(wf)}
    for mi, mech in enumerate(G["MECHS"]):
        if mech == "S":
            continue
        s = M.loc[(b, mech), "matched_strength"]
        t0 = pd.Timestamp.now()
        seeds = [np.random.SeedSequence([SEED, bi, mi, 1, 0, r]) for r in range(A.reps)]
        out = Parallel(n_jobs=A.jobs)(delayed(replicate)(parts, wf, mech, s, sd, sd_e, readouts) for sd in seeds)
        for r, o in enumerate(out):
            for name, (qa, qp) in o.items():
                rows.append(dict(base=b, mechanism=mech, strength=s, rep=r, readout=name, Qa=qa, Qphi=qp))
        say(b, mech, "done", round((pd.Timestamp.now() - t0).total_seconds()), "s"); pd.DataFrame(rows).to_csv(OUT / "replicates.csv", index=False)
R = pd.DataFrame(rows); R.to_csv(OUT / "replicates.csv", index=False)

cv = StratifiedKFold(5, shuffle=True, random_state=SEED)
res = []
for b in ["B1", "B2"]:
    for ro in ["g14", "core"]:
        for m1 in ["W", "H"]:
            d = R[(R.base == b) & (R.readout == ro) & R.mechanism.isin([m1, "D"])]
            y = (d.mechanism == m1).to_numpy().astype(int)
            x = np.log(d[["Qa", "Qphi"]].to_numpy())
            p = cross_val_predict(LinearDiscriminantAnalysis(), x, y, cv=cv, method="predict_proba")[:, 1]
            res.append(dict(base=b, readout=ro, pair=f"{m1} vs D", auc_joint=roc_auc_score(y, p), auc_Qa=roc_auc_score(y, d.Qa),
                            auc_Qphi=roc_auc_score(y, d.Qphi), median_Qa_mech=d[y == 1].Qa.median(), median_Qa_D=d[y == 0].Qa.median(),
                            median_Qphi_mech=d[y == 1].Qphi.median(), median_Qphi_D=d[y == 0].Qphi.median()))
Au = pd.DataFrame(res)
g = Au[Au.readout == "g14"]
Au["V351_pass_g14"] = bool((g.auc_joint >= 0.9).all()); Au.to_csv(OUT / "auc.csv", index=False)
print(Au.round(3).to_string()); say("done")
