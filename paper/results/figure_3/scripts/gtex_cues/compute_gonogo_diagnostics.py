"""A27 part 1, diagnostics D1, D2 and D4: the untimed per-gene step on GTEx bulk (post hoc, exploratory).

Post hoc, exploratory (docs/PLAN_internal_desync.md (archive), A27 part 1; no hypothesis ids). Estimand: per target tissue, every expressed
gene (median CPM >= 5) is regressed on the host state alone (log CPM ~ 1 + D2 covariates + host f, 2 columns; complex loading
c_est = beta_cos... of the host columns, peak phase in the host frame) and compared with the timed cosinor loading c_true of the
same tissue (recorded time of death, same covariates); score on the genes with timed p < 0.01 outside the clock frame:
R = |mean exp(i (arg c_est - arg c_true))| (concentration of the phase error), common rotation, median error after the rotation,
and, for windows, the up/down agreement along the window (slope at the window centre). Hosts: `clocksync`
(`Dataset.fit(...).host`) from the clock genes of the other `TISSUES`, target left out, no times used.
  d1  the step on the same donors: designs all fast deaths and the 08-18 h window; null = host permuted within the design (100)
  d2  degradation: host from random estimate donors (half, 45, 25), 14 or 3 other tissues, with or without the window; truth on
      the disjoint donors; a timed oracle on the same estimate donors; 5 replicates
  d4  ventilator cases (Hardy 0) as the "atlas": untimed host from the ventilated donors, genes of the target tissue regressed
      on it, scored against the timed cosinor of fast deaths; control: fast-death donors of the same n against the others;
      also the timed fit of the ventilated donors on recorded time; clock-frame genes reported separately
D3 (the atlas pseudobulk instead of the cell-type recombination) needs the cell-type atlases and is in compute_gonogo_atlas.py
(not here). Aggregates only (controlled access: per-donor values stay in memory).

    python compute_gonogo_diagnostics.py [--step d1 d2 d4] [--subsample N] [--jobs 4]
Outputs (paper/results/figure_3/output/gtex_cues/): diag1_summary.csv, diag1_genes.csv.gz, diag2_reps.csv, diag4_reps.csv
Legacy: legacy/desync/validation/gonogo/scripts/{diag1_gtex,diag2_degrade,diag4_ventilator}.py. Local: the per-gene regressions
(`_shared.partial_f`, `timed_cosinor`; not in clocksync).
"""

from __future__ import annotations

import zlib
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402
from scipy import stats  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import cli  # noqa: E402
from common.paths import out_dir  # noqa: E402

import _genes as G  # noqa: E402
import _shared as S  # noqa: E402
from _shared import W  # noqa: E402

FOLDER = G.FOLDER
WINDOW = (8.0, 18.0)
FRAME = set(cs.gtex().genes) | {"BMAL1"}
TARGETS_D1 = ["Heart - Atrial Appendage", "Muscle - Skeletal", "Lung", "Thyroid", "Artery - Tibial",
              "Adipose - Visceral (Omentum)"]
TARGETS_D2 = ["Heart - Atrial Appendage", "Muscle - Skeletal", "Lung", "Adipose - Visceral (Omentum)"]
CONFIGS = {  # name: (n estimate donors or "half", host tissues, window)
    "reference (half, 14, day)": ("half", 14, False),
    "45 donors": (45, 14, False),
    "25 donors": (25, 14, False),
    "45 donors, 3-tissue host": (45, 3, False),
    "45 donors, window": (45, 14, True),
    "atlas-like (45, 3, window)": (45, 3, True),
}


def in_window(h):
    return (h >= WINDOW[0]) & (h < WINDOW[1])


def load_tissue(tis, hardy, donors=None):
    """log(CPM + 0.5) of all genes, gene names and annotations of the rows of one tissue (time of death known)."""
    Y, genes, obs = S.expression(tis, hardy=hardy, tod=True, median_cpm=None)
    if donors is not None:
        m = obs.subject.isin(donors).to_numpy()
        Y, obs = Y[m], obs[m].reset_index(drop=True)
    return Y, genes, obs


def expressed(Y, cols=None):
    """Genes with median CPM >= 5 over the given rows."""
    return np.median(np.exp(Y if cols is None else Y[cols]) - 0.5, 0) >= 5


def loadings(Y, f, Z):
    D = np.c_[np.ones(len(Y)), Z, f]
    B, *_ = np.linalg.lstsq(D, Y, rcond=None)
    return B[-2] + 1j * B[-1]


def timed_truth(Y, t, Z):
    """Complex timed loading and 2-df F p of every gene (cosinor on t with covariates)."""
    return S.timed_cosinor(Y, t, Z)


def score(c_est, c_true):
    R = np.exp(1j * (np.angle(c_est) - np.angle(c_true))).mean()
    err = np.abs(np.angle(np.exp(1j * (np.angle(c_est) - np.angle(c_true) - np.angle(R))))) / W
    return abs(R), float(np.angle(R) / W % 24), float(np.median(err))


def host_of(pb, host_t, donors):
    """Host of the donors from the clock genes of the listed tissues that have rows among them (e.g. the ventilated donors
    have no pituitary or brain rows)."""
    sub = pb[pb.donor.isin(set(donors))]
    return S.host_states(pb, [t for t in host_t if t in set(sub.compartment)], donors=donors)


# ---- D1 -------------------------------------------------------------------------------------------------------------

def d1_one(design, tis, pb, seed, donors, n_perm):
    limit_threads(1)
    rng = np.random.default_rng(seed)
    Y, genes, ob = load_tissue(tis, (1, 2), donors)
    t = ob.tod.to_numpy(float)
    if design == "window":
        m = in_window(t)
        Y, ob, t = Y[m], ob[m].reset_index(drop=True), t[m]
        pbw = pb[pb.donor.isin(set(ob.subject))]
    else:
        pbw = pb
    host_t = [x for x in S.TISSUES if x != tis and x in set(pbw.compartment)]
    h = host_of(pbw, host_t, set(ob.subject))
    m = ob.subject.isin(h.index).to_numpy()
    Y, ob, t = Y[m], ob[m].reset_index(drop=True), t[m]
    ex = expressed(Y)
    Y, g = Y[:, ex], genes[ex]
    Z = S.d2_covariates(ob)
    f = h[["x", "y"]].loc[ob.subject].to_numpy()
    c_true, p = timed_truth(Y, t, Z)
    c_est = loadings(Y, f, Z)
    frame = np.isin(g, list(FRAME))
    sel = (p < 0.01) & ~frame
    nulls = [abs(np.exp(1j * (np.angle(loadings(Y[:, sel], f[rng.permutation(len(f))], Z)) - np.angle(c_true[sel]))).mean())
             for _ in range(n_perm)]
    R, rot, med = score(c_est[sel], c_true[sel])
    out = dict(design=design, tissue=tis, donors=len(t), host_tissues=len(host_t), n_rhythmic=int(sel.sum()), R=R, rot_h=rot,
               med_err_h=med, null_R_95=float(np.quantile(nulls, 0.95)), p_R=(1 + np.sum(np.array(nulls) >= R)) / (n_perm + 1))
    if design == "window":
        tc = np.angle(np.exp(1j * W * t).mean()) / W
        e = np.exp(1j * W * tc)
        se_, st_ = np.imag(np.conj(c_est[sel]) * e), np.imag(np.conj(c_true[sel]) * e)
        out.update(slope_rho=stats.spearmanr(se_, st_)[0], slope_sign=float(np.mean(np.sign(se_) == np.sign(st_))))
    fr = score(c_est[frame], c_true[frame])
    out.update(frame_R=fr[0], frame_rot_h=fr[1])
    return out, pd.DataFrame(dict(design=design, tissue=tis, gene=g, p=p, c_true=c_true, c_est=c_est))


# ---- D2 -------------------------------------------------------------------------------------------------------------

def d2_one(tis, cfg, rep, pb, data, quick):
    limit_threads(1)
    n, nh, window = CONFIGS[cfg]
    Y, genes, ob = data
    rng = np.random.default_rng(1000 * rep + zlib.crc32(cfg.encode()) % 997)
    subj = ob.subject.to_numpy()
    t = ob.tod.to_numpy(float)
    pool = np.flatnonzero(in_window(t)) if window else np.arange(len(t))
    k = len(t) // 2 if n == "half" else n
    if len(pool) < 25 or k > len(pool) - 10:
        return dict(tissue=tis, config=cfg, rep=rep, status="too few donors")
    est = rng.choice(pool, size=min(k, len(pool) - 10), replace=False)
    tru = np.setdiff1d(np.arange(len(t)), est)
    c_true, p = timed_truth(Y[tru], t[tru], S.d2_covariates(ob.iloc[tru]))
    others = [x for x in S.TISSUES if x != tis]
    host_t = list(rng.choice(others, nh, replace=False)) if nh < len(others) else others
    try:
        h = host_of(pb, host_t, set(subj[est]))
    except Exception as ex:  # noqa: BLE001
        return dict(tissue=tis, config=cfg, rep=rep, status=f"host error {ex}")
    e2 = est[np.isin(subj[est], h.index)]
    f = h[["x", "y"]].loc[subj[e2]].to_numpy()
    Ze = S.d2_covariates(ob.iloc[e2])
    c_est = loadings(Y[e2], f, Ze)
    c_or, _ = timed_truth(Y[e2], t[e2], Ze)
    frame = np.isin(genes, list(FRAME))
    sel = (p < 0.01) & ~frame
    R, rot, med = score(c_est[sel], c_true[sel])
    rot_f = np.angle(np.exp(1j * (np.angle(c_est[frame]) - np.angle(c_true[frame]))).mean())
    tc = np.angle(np.exp(1j * W * t[e2]).mean()) / W
    e = np.exp(1j * W * tc)
    ce = c_est[sel] * np.exp(-1j * rot_f)
    se_, st_ = np.imag(np.conj(ce) * e), np.imag(np.conj(c_true[sel]) * e)
    so_ = np.imag(np.conj(c_or[sel]) * e)
    big = np.abs(st_) > np.quantile(np.abs(st_), 0.5)
    Ro, _, medo = score(c_or[sel], c_true[sel])
    return dict(tissue=tis, config=cfg, rep=rep, status="ok", oracle_R=Ro, oracle_med_err_h=medo,
                oracle_slope_sign=float(np.mean(np.sign(so_) == np.sign(st_))),
                oracle_slope_sign_clear=float(np.mean(np.sign(so_[big]) == np.sign(st_[big]))), n_est=len(e2),
                n_truth=len(tru), host_tissues=len(host_t), genes=int(sel.sum()), R=R, med_err_h=med, rot_h=rot,
                frame_rot_h=float(rot_f / W), slope_sign=float(np.mean(np.sign(se_) == np.sign(st_))),
                slope_sign_clear=float(np.mean(np.sign(se_[big]) == np.sign(st_[big]))),
                slope_rho=stats.spearmanr(se_, st_)[0])


# ---- D4 -------------------------------------------------------------------------------------------------------------

def d4_one(tis, arm, rep, pb_f, pb_v, data):
    limit_threads(1)
    (Yf, obf), (Yv, obv), genes, n = data
    rng = np.random.default_rng(97 * rep + (arm == "ventilator"))
    tf = obf.tod.to_numpy(float)
    if arm == "ventilator":
        est_Y, est_ob, pbx = Yv, obv, pb_v
        ie, tru = rng.choice(len(Yv), size=n, replace=False), np.arange(len(Yf))
    else:
        est_Y, est_ob, pbx = Yf, obf, pb_f
        perm = rng.permutation(len(Yf))
        ie, tru = perm[:n], perm[n:]
    c_true, p = timed_truth(Yf[tru], tf[tru], S.d2_covariates(obf.iloc[tru]))
    subj = est_ob.subject.to_numpy()[ie]
    host_t = [x for x in S.TISSUES if x != tis]
    h = host_of(pbx, host_t, set(subj))
    m = np.isin(subj, h.index)
    ie, subj = ie[m], subj[m]
    f = h[["x", "y"]].loc[subj].to_numpy()
    Ze = S.d2_covariates(est_ob.iloc[ie])
    c_est = loadings(est_Y[ie], f, Ze)
    c_timed, _ = timed_truth(est_Y[ie], est_ob.tod.to_numpy(float)[ie], Ze)
    frame = np.isin(genes, list(FRAME))
    sel = (p < 0.01) & ~frame
    out = dict(tissue=tis, arm=arm, rep=rep, n=len(ie), genes=int(sel.sum()))
    for lab, c in (("untimed", c_est), ("timed", c_timed)):
        out[f"{lab}_R"], out[f"{lab}_rot_h"], _ = score(c[sel], c_true[sel])
        out[f"{lab}_frame_R"], out[f"{lab}_frame_rot_h"], _ = score(c[frame], c_true[frame])
    return out


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--step", nargs="+", default=["d1", "d2", "d4"], choices=["d1", "d2", "d4"])
    args = ap.parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    o = G.obs_table()
    sub = G.pick_donors(o, args.subsample, args.seed)
    n_rep = 2 if run["quick"] else 5
    n_perm = 20 if run["quick"] else 100
    pb_f, _ = S.clock_table(hardy=(1, 2), require_tod=True)
    pd.set_option("display.width", 250)
    if "d1" in args.step:
        tasks = [(d, t) for d in ("all", "window") for t in TARGETS_D1]
        res = Parallel(n_jobs=args.jobs)(delayed(d1_one)(d, t, pb_f[pb_f.donor.isin(sub)] if sub else pb_f, i, sub, n_perm)
                                         for i, (d, t) in enumerate(tasks))
        df = pd.DataFrame([r[0] for r in res])
        df.to_csv(out / "diag1_summary.csv", index=False)
        pd.concat([r[1] for r in res]).to_csv(out / "diag1_genes.csv.gz", index=False)
        print(df.round(3).to_string(index=False), flush=True)
    if "d2" in args.step:
        pbs = pb_f[pb_f.donor.isin(sub)] if sub else pb_f
        data = {}
        for t in TARGETS_D2:
            Y, genes, ob = load_tissue(t, (1, 2), sub)
            ex = expressed(Y)
            data[t] = (Y[:, ex], genes[ex], ob)
        tasks = [(t, c, r) for t in TARGETS_D2 for c in CONFIGS for r in range(n_rep)]
        rows = Parallel(n_jobs=args.jobs)(delayed(d2_one)(t, c, r, pbs, data[t], run["quick"]) for t, c, r in tasks)
        df = pd.DataFrame(rows)
        df.to_csv(out / "diag2_reps.csv", index=False)
        ok = df[df.status == "ok"]
        print(ok.groupby("config", sort=False)[["R", "oracle_R", "med_err_h", "oracle_med_err_h", "slope_sign",
                                               "oracle_slope_sign", "n_est", "genes"]].mean().round(2).to_string(), flush=True)
        print(df.status.value_counts().to_string())
    if "d4" in args.step:
        pb_v, _ = S.clock_table(hardy=(0,), require_tod=True)
        if sub:
            pb_f, pb_v = pb_f[pb_f.donor.isin(sub)], pb_v[pb_v.donor.isin(sub)]
        data, tasks = {}, []
        for t in TARGETS_D1:
            (Yf, gf, obf), (Yv, gv, obv) = load_tissue(t, (1, 2), sub), load_tissue(t, (0,), sub)
            if len(Yv) < 20 or len(Yf) < 40:
                print("D4 skipped (too few donors):", t, len(Yf), len(Yv))
                continue
            ex = expressed(Yf) & expressed(Yv)
            data[t] = ((Yf[:, ex], obf), (Yv[:, ex], obv), gf[ex], min(len(Yv), len(Yf) // 2))
            tasks += [(t, a, r) for a in ("fast-death control", "ventilator") for r in range(n_rep)]
        rows = Parallel(n_jobs=args.jobs)(delayed(d4_one)(t, a, r, pb_f, pb_v, data[t]) for t, a, r in tasks)
        df = pd.DataFrame(rows)
        df.to_csv(out / "diag4_reps.csv", index=False)
        cm = lambda x: float(np.angle(np.exp(1j * W * x).mean()) / W % 24)  # noqa: E731
        agg = {c: (cm if c.endswith("rot_h") else "mean") for c in df.columns if c not in ("tissue", "arm", "rep")}
        print(df.groupby("arm").agg(agg).round(2).to_string(), flush=True)


if __name__ == "__main__":
    main()
