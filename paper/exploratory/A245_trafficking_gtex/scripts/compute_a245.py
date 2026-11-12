"""A245: do endothelial adhesion genes and the T-cell content of healthy human tissues follow the hours implied by A229?

Registered in paper/docs/registered_experiments.md (A245). GTEx v10 fast deaths (Hardy 1-2) with a time of death; target
tissues = non-cultured tissues except brain regions, whole blood and spleen, with >= 40 donors. Two time axes per donor:
recorded time of death, and internal time = phase of the untimed body host of A234 (refitted without the target tissue when
it is one of the 12 host tissues). Readouts (log2(CPM + 1)): VCAM1, ICAM1, SELE as residuals on an endothelial content
score (PECAM1, CDH5, VWF) within tissue; T-cell content (CD3D, CD3E, CD2, TRAC); PTPRC.

Statistic: per tissue the OLS t of y on cos(w (time - h)), h the predicted peak, with Hardy x sex indicators, exact age and
the sample's ischaemic time; S = mean of t over tissues. Phase-free companion: R = |mean over tissues of (t_cos + i t_sin)|
of the two-term cosinor with the same covariates. Null: times permuted across donors within Hardy x sex, one permutation
for all tissues, readouts and both statistics. Stages:
  power  tissues, donors, null SD of S and R (permuted times only; no observed statistic is computed or printed), MDE_80 =
         2.49 x null SD, and the log2 amplitude the MDE of S corresponds to in the median tissue
                                                              -> output/{tissues,direction,power}.csv
  score  (after the power note is dated) observed S, R, p, Holm, matched-gene and pseudo-clock companions, bootstrap of
         S(internal) - S(recorded), descriptives   -> output/score_{tests,matched,pseudo,per_tissue,sweep,descriptive}.csv
  `--placebo`: recorded and internal times shuffled once over donors (code test; output/placebo/).

    python compute_a245.py --stage power|score [--subsample N] [--jobs 4]
Controlled access: time of death, exact age and ischaemic time stay in memory; only aggregates are written.

Choices not fixed by the registration
- Donors: fast deaths with a time of death and >= 5 of the 12 body tissues (the A234 body host); the same donors on both
  axes. Leave-one-out hosts are fitted on these donors without the target tissue. Each host is oriented by the sign of its
  Fisher-Lee correlation with time of death (D234.0).
- Strata with < 3 donors are merged into the other sex of the same Hardy class (A234). Tissue ischaemic time is the
  sample's; a missing value takes the tissue median.
- Matched companions: genes with mean count >= 1, outside the clock, readout and endothelial-score genes; for a gene readout
  200 genes drawn from the 400 closest in mean log2(CPM + 1) over the target tissues, residualised on the endothelial score
  like the readout; for the T-cell score 200 sets, one gene drawn for each of the four from its 400 closest.
- Pseudo-clock hosts (internal time only): 50 expression-matched sets (`compute_a208.pseudo_sets`), same leave-one-out
  scheme and orientation.
- Bootstrap of S(internal) - S(recorded): donors resampled with replacement, hosts not refitted.
- Ventilated donors (descriptive): Hardy 0, own body host fitted with the same rules, internal time only.
- Holm: primaries H245.1-2; secondaries H245.3 (ICAM1, SELE) and H245.4 (VCAM1, T-cell content on internal time) together.
"""

from __future__ import annotations

import sys
import time
import zlib
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "exploratory" / "A208_sex_coupling" / "scripts"))
sys.path.insert(0, str(ROOT / "exploratory" / "A234_internal_time_age" / "scripts"))
import clocksync as cs  # noqa: E402
from clocksync import stats as CS  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import cli, gtex  # noqa: E402

import compute_a234 as H  # noqa: E402  (body host, orientation, strata)

OUT = HERE.parent / "output"
W = 2 * np.pi / 24
ENDO = ["PECAM1", "CDH5", "VWF"]
TCELL = ["CD3D", "CD3E", "CD2", "TRAC"]
# readout -> (genes, residualised on the endothelial score, predicted peak h of A229 or None)
READOUTS = {"VCAM1": (["VCAM1"], True, 20.5), "T_cell": (TCELL, False, 4.8), "ICAM1": (["ICAM1"], True, 18.5),
            "SELE": (["SELE"], True, 19.0), "PTPRC": (["PTPRC"], False, None), "T_cell_resid_HBB": (TCELL, False, 4.8)}
TESTS = [("H245.1", "primary", "VCAM1", "recorded"), ("H245.2", "primary", "T_cell", "recorded"),
         ("H245.3", "secondary", "ICAM1", "recorded"), ("H245.3", "secondary", "SELE", "recorded"),
         ("H245.4", "secondary", "VCAM1", "internal"), ("H245.4", "secondary", "T_cell", "internal")]
AXES = ("recorded", "internal")
MIN_DONORS, MIN_TISSUES, MIN_HOST = 40, 10, 5
N_PERM, N_MATCH, N_POOL, N_PSEUDO, N_BOOT, MDE_K = 10_000, 200, 400, 50, 2000, 2.49


def say(*a):
    print(time.strftime("%F %T"), *a, flush=True)


def seed_of(*parts) -> int:
    return zlib.crc32("|".join(map(str, parts)).encode()) % (2 ** 31)


# ---- data -------------------------------------------------------------------------------------------------------------

def target(t: str) -> bool:
    return not (t.startswith("Cells -") or t.startswith("Brain") or t in ("Whole Blood", "Spleen"))


def load_lcpm(hardy):
    """(log2(CPM + 1) samples x genes as float32, gene names, sample frame: donor, tissue, ischaemic_h) of target tissues."""
    import anndata as ad
    a = ad.read_h5ad(gtex.LABELLED, backed="r")
    try:
        o = a.obs
        tis = o.tissue.astype(str)
        idx = np.flatnonzero(o.DTHHRDY.isin(list(hardy)).to_numpy() & o.tod.notna().to_numpy() & tis.map(target).to_numpy())
        X = a.X[idx]
        X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
        genes = np.asarray(a.var_names)
        ob = o.iloc[idx]
    finally:
        a.file.close()
    mean_count = X.mean(0)
    L = np.log2(X / X.sum(1, keepdims=True) * 1e6 + 1).astype(np.float32)
    sm = pd.DataFrame({"donor": ob.subject.astype(str).to_numpy(), "tissue": ob.tissue.astype(str).to_numpy(),
                       "ischaemic_h": ob.ischemic_time.astype(float).to_numpy() / 60})
    first = ~sm.duplicated(["donor", "tissue"]).to_numpy()
    return L[first], genes, sm[first].reset_index(drop=True), mean_count


def hosts(hardy, tissues, waveform=None, donors_keep=None):
    """Donor frame (hardy, tod, male, exact age) and {tissue: oriented internal time (h) per donor}, key '' = full host."""
    pb, donors = gtex.load(tissues=list(H.BODY), hardy=hardy, genes="clock" if waveform is None else list(waveform.genes))
    n = pb.groupby("donor").compartment.nunique()
    keep = n.index[n >= MIN_HOST] if donors_keep is None else donors_keep
    pb, donors = pb[pb.donor.isin(keep)].reset_index(drop=True), donors.loc[sorted(keep)].copy()
    donors["age_exact"] = gtex.phenotypes(["AGE"]).AGE.reindex(donors.index).astype(float)
    out, fl = {}, {}
    for t in [""] + [t for t in tissues if t in H.BODY]:
        p = pb[pb.compartment != t]
        kw = {} if waveform is None else {"waveform": waveform}
        ph = cs.Dataset(p, anchors=[a for a in H.BODY if a in set(p.compartment)], **kw).host().phase_h
        ph, fl[t] = H.orient(ph.reindex(donors.index), donors.tod)
        out[t] = ph % 24
    return donors, out, fl


def design(L, genes, sm, donors, times):
    """Per tissue: sample rows, donor positions, covariate matrix C (intercept, strata, age, ischaemic time) and its Q."""
    col = {g: i for i, g in enumerate(genes)}
    sm = sm[sm.donor.isin(donors.index)]
    n = sm.groupby("tissue").donor.nunique()
    codes = pd.Series(H.strata(donors, donors.index), index=donors.index)
    units = {}
    for t in sorted(n.index[n >= MIN_DONORS]):
        s = sm[sm.tissue == t]
        d = s.donor.to_numpy()
        isc = s.ischaemic_h.fillna(s.ischaemic_h.median()).to_numpy()
        st = pd.get_dummies(codes.loc[d].to_numpy(), drop_first=True).to_numpy(float)
        C = np.column_stack([np.ones(len(d)), st, donors.age_exact.loc[d].to_numpy(), isc])
        units[t] = dict(rows=s.index.to_numpy(), pos=donors.index.get_indexer(d), C=C, Q=np.linalg.qr(C)[0])
    return units, col, codes.to_numpy()


def resid(Q, Y):
    return Y - Q @ (Q.T @ Y)


def readout(L, col, u, name=None, genes_=None, endo=None):
    """y of one tissue: mean log2(CPM + 1) over the genes, residualised on the endothelial score (or HBB) if asked."""
    g, on_endo, _ = READOUTS[name] if name else (genes_, endo, None)
    y = L[np.ix_(u["rows"], [col[x] for x in g])].astype(float).mean(1)
    cov = ENDO if on_endo else ["HBB"] if name == "T_cell_resid_HBB" else None
    if cov:
        e = L[np.ix_(u["rows"], [col[x] for x in cov])].astype(float).mean(1)
        A = np.column_stack([np.ones(len(e)), e])
        y = y - A @ np.linalg.lstsq(A, y, rcond=None)[0]
    return y


# ---- statistics -------------------------------------------------------------------------------------------------------

def t_cos(Q, Yr, x):
    """t of the slope on x (1-D or perms x n) for every column of Yr (already residualised on the covariates)."""
    xr = np.atleast_2d(x) - (np.atleast_2d(x) @ Q) @ Q.T
    sxx = (xr ** 2).sum(1)[:, None]
    b = (xr @ Yr) / sxx
    df = Yr.shape[0] - Q.shape[1] - 1
    rss = (Yr ** 2).sum(0)[None, :] - b ** 2 * sxx
    return b / np.sqrt(rss / df / sxx)


def t_cosinor(Q, Yr, tm):
    """t_cos + i t_sin of the two-term cosinor (tm: perms x n hours), plus amplitude and peak hour of the fit."""
    c, s = np.cos(W * np.atleast_2d(tm)), np.sin(W * np.atleast_2d(tm))
    c, s = c - (c @ Q) @ Q.T, s - (s @ Q) @ Q.T
    cc, ss, sc = (c * c).sum(1)[:, None], (s * s).sum(1)[:, None], (c * s).sum(1)[:, None]
    det = cc * ss - sc ** 2
    cy, sy = c @ Yr, s @ Yr
    bc, bs = (ss * cy - sc * sy) / det, (cc * sy - sc * cy) / det
    df = Yr.shape[0] - Q.shape[1] - 2
    sig2 = ((Yr ** 2).sum(0)[None, :] - bc * cy - bs * sy) / df
    z = bc / np.sqrt(sig2 * ss / det) + 1j * bs / np.sqrt(sig2 * cc / det)
    return z, np.hypot(bc, bs), (np.arctan2(bs, bc) / W) % 24


def permutations(codes, n, seed):
    """n permutations of the donor positions within strata (rows)."""
    rng = np.random.default_rng(seed)
    idx = [np.flatnonzero(codes == c) for c in np.unique(codes)]
    P = np.tile(np.arange(len(codes)), (n, 1))
    for i in idx:
        P[:, i] = i[np.argsort(rng.random((n, len(i))), axis=1)]
    return P


def tissue_stats(u, Y, hs, donor_time, P=None):
    """Per axis: t (perms x readouts) at the predicted hours and complex cosinor t. P None = observed (identity)."""
    Yr = resid(u["Q"], Y)
    out = {}
    for ax in AXES:
        full = donor_time[ax] if ax == "recorded" else u["time_full"]
        tm = full[u["pos"]][None, :] if P is None else full[P][:, u["pos"]]
        t = np.column_stack([t_cos(u["Q"], Yr[:, [j]], np.cos(W * (tm - h)))[:, 0] if h is not None
                             else np.full(len(tm), np.nan) for j, h in enumerate(hs)])
        out[ax] = (t, t_cosinor(u["Q"], Yr, tm)[0])
    return out


def attach_full_time(units, donors, times):
    """Internal time over ALL donors for each tissue (its own leave-one-out host), so a permuted donor carries a time."""
    for t, u in units.items():
        u["time_full"] = times.get(t, times[""]).reindex(donors.index).to_numpy(float)


def pooled(units, L, col, names, donors, P=None):
    """S (mean t over tissues) and R (|mean complex t|) per readout and axis: arrays (perms or 1) x readouts."""
    hs = [READOUTS[k][2] for k in names]
    dt = {"recorded": donors.tod.to_numpy(float)}
    acc = {ax: [0, 0] for ax in AXES}
    for u in units.values():
        Y = np.column_stack([readout(L, col, u, k) for k in names])
        r = tissue_stats(u, Y, hs, dt, P)
        for ax in AXES:
            acc[ax][0] = acc[ax][0] + r[ax][0]
            acc[ax][1] = acc[ax][1] + r[ax][1]
    k = len(units)
    return {ax: (acc[ax][0] / k, np.abs(acc[ax][1] / k)) for ax in AXES}


def _null_chunk(units, L, col, names, donors, codes, n, seed):
    limit_threads(1)
    return pooled(units, L, col, names, donors, permutations(codes, n, seed))


def null_of(units, L, col, names, donors, codes, n, jobs, seed):
    parts = Parallel(n_jobs=jobs)(delayed(_null_chunk)(units, L, col, names, donors, codes, len(c), seed + j)
                                  for j, c in enumerate(np.array_split(np.arange(n), max(jobs, 1) * 4)))
    return {ax: (np.vstack([p[ax][0] for p in parts]), np.vstack([p[ax][1] for p in parts])) for ax in AXES}


# ---- stages -----------------------------------------------------------------------------------------------------------

def prepare(hardy, args, quick):
    L, genes, sm, mean_count = load_lcpm(hardy)
    tissues = sorted(sm.tissue.unique())
    donors, times, fl = hosts(hardy, tissues)
    if quick:
        donors = donors.loc[np.random.default_rng(args.seed).permutation(sorted(donors.index))[:args.subsample]]
    units, col, codes = design(L, genes, sm, donors, times)
    if args.placebo:
        perm = np.random.default_rng(seed_of("A245placebo", args.seed)).permutation(len(donors))
        donors = donors.assign(tod=donors.tod.to_numpy()[perm])
        times = {k: pd.Series(v.reindex(donors.index).to_numpy()[perm], index=donors.index) for k, v in times.items()}
    attach_full_time(units, donors, times)
    return L, genes, sm, mean_count, donors, times, fl, units, col, codes


def power(args, quick):
    L, genes, sm, _, donors, times, fl, units, col, codes = prepare(gtex.FAST, args, quick)
    pd.DataFrame([dict(tissue=t, donors=len(u["rows"]), host="leave-one-out" if t in times else "full") for t, u in units.items()]
                 ).to_csv(OUT / "tissues.csv", index=False)
    pd.DataFrame([dict(host=k or "full", fisher_lee_tod=v, negated=v < 0, donors=len(donors)) for k, v in fl.items()]
                 ).to_csv(OUT / "direction.csv", index=False)
    names = list(READOUTS)
    nul = null_of(units, L, col, names, donors, codes, 200 if quick else N_PERM, args.jobs, seed_of("A245null", args.seed))
    rows = []
    for ax in AXES:
        for j, k in enumerate(names):
            h = READOUTS[k][2]
            se = []                                           # SE of the slope per tissue: y and time never paired
            for u in units.values():
                yr = resid(u["Q"], readout(L, col, u, k))
                tm = (donors.tod.to_numpy(float) if ax == "recorded" else u["time_full"])[u["pos"]]
                x = np.cos(W * (tm - (h or 0)))
                xr = x - u["Q"] @ (u["Q"].T @ x)
                se.append(np.sqrt((yr ** 2).sum() / (len(yr) - u["Q"].shape[1] - 1) / (xr ** 2).sum()))
            sd_s, sd_r = (float(nul[ax][0][:, j].std()) if h is not None else np.nan), float(nul[ax][1][:, j].std())
            rows.append(dict(axis=ax, readout=k, predicted_peak_h=h, tissues=len(units), donors=len(donors), null_sd_S=sd_s,
                             mde_S=MDE_K * sd_s, mde_amplitude_log2_median_tissue=MDE_K * sd_s * float(np.median(se)),
                             null_mean_R=float(nul[ax][1][:, j].mean()), null_sd_R=sd_r, null_q95_R=float(np.quantile(nul[ax][1][:, j], .95)),
                             scored=len(units) >= MIN_TISSUES))
    P = pd.DataFrame(rows)
    P.to_csv(OUT / "power.csv", index=False)
    print(P.round(4).to_string(), flush=True)


def matched_sets(name, genes, col, mean_l, ok, rng):
    """N_MATCH gene lists expression-matched to the readout's genes (one draw per gene from its N_POOL closest)."""
    cand = np.flatnonzero(ok)
    picks = []
    for g in READOUTS[name][0]:
        near = cand[np.argsort(np.abs(mean_l[cand] - mean_l[col[g]]))[:N_POOL]]
        picks.append(rng.choice(near, N_MATCH, replace=len(near) < N_MATCH))
    return [[genes[p[i]] for p in picks] for i in range(N_MATCH)]


def _pseudo_one(w, hardy, tissues, donors_index, L, genes, sm, names):
    limit_threads(1)
    donors, times, _ = hosts(hardy, tissues, waveform=w, donors_keep=donors_index)
    units, col, _ = design(L, genes, sm, donors, times)
    attach_full_time(units, donors, times)
    return pooled(units, L, col, names, donors)["internal"][0][0]


def score(args, quick):
    import compute_a208 as B
    L, genes, sm, mean_count, donors, times, fl, units, col, codes = prepare(gtex.FAST, args, quick)
    names = list(READOUTS)
    n_perm = 200 if quick else N_PERM
    obs = pooled(units, L, col, names, donors)
    nul = null_of(units, L, col, names, donors, codes, n_perm, args.jobs, seed_of("A245null", args.seed))
    say("observed and null done")

    # matched-gene companion (both axes), same tissues and hours
    excl = set(cs.gtex().genes) | set(ENDO) | {"HBB"} | {g for v in READOUTS.values() for g in v[0]}
    ok = (mean_count >= 1) & ~np.isin(genes, list(excl))
    mean_l = np.mean([L[u["rows"]].mean(0) for u in units.values()], axis=0)
    dt = {"recorded": donors.tod.to_numpy(float)}
    mrows, mq = [], {}
    for k in ("VCAM1", "T_cell", "ICAM1", "SELE"):
        sets = matched_sets(k, genes, col, mean_l, ok, np.random.default_rng(seed_of("A245match", k)))
        acc = {ax: 0 for ax in AXES}
        for u in units.values():
            Y = np.column_stack([readout(L, col, u, genes_=s, endo=READOUTS[k][1]) for s in sets])
            r = tissue_stats(u, Y, [READOUTS[k][2]] * len(sets), dt)
            for ax in AXES:
                acc[ax] = acc[ax] + r[ax][0][0]
        for ax in AXES:
            v = acc[ax] / len(units)
            mq[(k, ax)] = float(np.quantile(v, .95))
            mrows += [dict(readout=k, axis=ax, set=i, genes=";".join(s), S=float(x)) for i, (s, x) in enumerate(zip(sets, v))]
    pd.DataFrame(mrows).to_csv(OUT / "score_matched.csv", index=False)
    say("matched companions done")

    # pseudo-clock hosts (internal time)
    ws = B.pseudo_sets(gtex.FAST, 4 if quick else N_PSEUDO, seed_of("A245pseudo", args.seed))
    ps = np.array(Parallel(n_jobs=args.jobs)(delayed(_pseudo_one)(w, gtex.FAST, sorted(units), donors.index, L, genes, sm, names)
                                             for w in ws))
    pd.DataFrame(ps, columns=names).to_csv(OUT / "score_pseudo.csv", index=False)
    say("pseudo-clock hosts done")

    # bootstrap of S(internal) - S(recorded), hosts fixed
    rng = np.random.default_rng(seed_of("A245boot", args.seed))
    D = []
    ids = np.arange(len(donors))
    for _ in range(50 if quick else N_BOOT):
        draw = rng.choice(ids, len(ids))
        cnt = np.bincount(draw, minlength=len(ids))
        acc = {ax: 0 for ax in AXES}
        k = 0
        for u in units.values():
            wts = cnt[u["pos"]]
            if (wts > 0).sum() < u["C"].shape[1] + 5:
                continue
            rep = np.repeat(np.arange(len(wts)), wts)
            Y = np.column_stack([readout(L, col, u, n)[rep] for n in names[:2]])
            Q = np.linalg.qr(u["C"][rep])[0]
            Yr = resid(Q, Y)
            for ax in AXES:
                tm = (dt["recorded"] if ax == "recorded" else u["time_full"])[u["pos"]][rep]
                acc[ax] = acc[ax] + np.array([t_cos(Q, Yr[:, [j]], np.cos(W * (tm - READOUTS[n][2])))[0, 0]
                                              for j, n in enumerate(names[:2])])
            k += 1
        D.append((acc["internal"] - acc["recorded"]) / k)
    D = np.array(D)
    say("bootstrap done")

    rows = []
    for hyp, role, k, ax in TESTS:
        j = names.index(k)
        s, n0 = float(obs[ax][0][0, j]), nul[ax][0][:, j]
        r, n1 = float(obs[ax][1][0, j]), nul[ax][1][:, j]
        row = dict(hypothesis=hyp, role=role, readout=k, axis=ax, predicted_peak_h=READOUTS[k][2], tissues=len(units),
                   donors=len(donors), S=s, null_sd_S=float(n0.std()), p_one_sided=(1 + np.sum(n0 >= s)) / (1 + len(n0)),
                   matched_q95=mq[(k, ax)], above_matched=s > mq[(k, ax)], R=r, null_q95_R=float(np.quantile(n1, .95)),
                   p_R=(1 + np.sum(n1 >= r)) / (1 + len(n1)), mde_S=MDE_K * float(n0.std()))
        if ax == "internal":
            row.update(pseudo_q95=float(np.quantile(ps[:, j], .95)), above_pseudo=s > float(np.quantile(ps[:, j], .95)),
                       D_internal_minus_recorded=s - float(obs["recorded"][0][0, j]),
                       D_lo90=float(np.quantile(D[:, j], .05)), D_hi90=float(np.quantile(D[:, j], .95)))
        rows.append(row)
    T = pd.DataFrame(rows)
    for role in ("primary", "secondary"):
        m = (T.role == role).to_numpy()
        p = T.loc[m, "p_one_sided"].to_numpy()
        o = np.argsort(p)
        T.loc[m, "p_holm"] = np.maximum.accumulate(np.minimum(1, p[o] * (len(p) - np.arange(len(p)))))[np.argsort(o)]
    T["pass"] = (T.p_holm < 0.05) & T.above_matched & (T.get("above_pseudo", True).fillna(True).astype(bool))
    T.to_csv(OUT / "score_tests.csv", index=False)
    print(T.round(4).to_string(), flush=True)

    # descriptives: every readout and axis (S, R), per-tissue cosinor, offset sweep, ventilated donors
    desc = [dict(cohort="fast", axis=ax, readout=k, predicted_peak_h=READOUTS[k][2], S=float(obs[ax][0][0, j]),
                 p_S=(1 + np.sum(nul[ax][0][:, j] >= obs[ax][0][0, j])) / (1 + n_perm) if READOUTS[k][2] is not None else np.nan,
                 R=float(obs[ax][1][0, j]), p_R=(1 + np.sum(nul[ax][1][:, j] >= obs[ax][1][0, j])) / (1 + n_perm))
            for ax in AXES for j, k in enumerate(names)]
    pt, sw = [], []
    for t, u in units.items():
        Y = np.column_stack([readout(L, col, u, k) for k in names])
        Yr = resid(u["Q"], Y)
        for ax in AXES:
            tm = (dt["recorded"] if ax == "recorded" else u["time_full"])[u["pos"]]
            z, amp, pk = t_cosinor(u["Q"], Yr, tm)
            for j, k in enumerate(names):
                pt.append(dict(tissue=t, axis=ax, readout=k, donors=len(tm), amplitude_log2=float(amp[0, j]), peak_h=float(pk[0, j]),
                               t_cos=float(z[0, j].real), t_sin=float(z[0, j].imag)))
    pd.DataFrame(pt).to_csv(OUT / "score_per_tissue.csv", index=False)
    for off in (-6, -4, -2, 0, 2, 4, 6):
        acc = {ax: 0 for ax in AXES}
        for u in units.values():
            Y = np.column_stack([readout(L, col, u, k) for k in names])
            r = tissue_stats(u, Y, [None if READOUTS[k][2] is None else READOUTS[k][2] + off for k in names], dt)
            for ax in AXES:
                acc[ax] = acc[ax] + r[ax][0][0]
        sw += [dict(offset_h=off, axis=ax, readout=k, S=float(acc[ax][j] / len(units))) for ax in AXES for j, k in enumerate(names)
               if READOUTS[k][2] is not None]
    pd.DataFrame(sw).to_csv(OUT / "score_sweep.csv", index=False)
    if not quick and not args.placebo:
        Lv, gv, smv, _ = load_lcpm(gtex.VENTILATOR)
        dv, tv, _ = hosts(gtex.VENTILATOR, sorted(smv.tissue.unique()))
        uv, cv, codv = design(Lv, gv, smv, dv, tv)
        attach_full_time(uv, dv, tv)
        ov = pooled(uv, Lv, cv, names, dv)["internal"]
        nv = null_of(uv, Lv, cv, names, dv, codv, 2000, args.jobs, seed_of("A245vent", args.seed))["internal"]
        desc += [dict(cohort="ventilated", axis="internal", readout=k, predicted_peak_h=READOUTS[k][2], tissues=len(uv),
                      donors=len(dv), S=float(ov[0][0, j]),
                      p_S=(1 + np.sum(nv[0][:, j] >= ov[0][0, j])) / 2001 if READOUTS[k][2] is not None else np.nan,
                      R=float(ov[1][0, j]), p_R=(1 + np.sum(nv[1][:, j] >= ov[1][0, j])) / 2001) for j, k in enumerate(names)]
    pd.DataFrame(desc).to_csv(OUT / "score_descriptive.csv", index=False)
    say("done")


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--stage", choices=["power", "score"], required=True)
    ap.add_argument("--placebo", action="store_true", help="times shuffled once over donors (code test; output/placebo/)")
    args = ap.parse_args(argv)
    run = cli.setup(args)
    global OUT
    OUT = OUT / "placebo" if args.placebo else OUT / "quick" if run["quick"] else OUT
    OUT.mkdir(parents=True, exist_ok=True)
    {"power": power, "score": score}[args.stage](args, run["quick"])


if __name__ == "__main__":
    main()
