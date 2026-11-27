"""A409 (registered in paper/docs/registered_experiments.md, "A409"): does the patient's clock loudness reach the tumour?

A387 pooled table (14 genes). Radius covariance (`clocksync.elasticity.radius_covariance`) between a target's clock-plane
state and a host state of the same patients, radii residualised on the phase harmonics of both.
  power  null of the H409.1 leverage from within-stratum shuffles of the tumour-fibroblast states (identity pairing not
         evaluated); MDE                                                         -> output/a409/power.csv
  score  H409.1 (tumour fibroblast vs normal-stroma host), references, companions (200 matched / 200 responsive pseudo
         waveforms, the same for target and host), D409 myeloid slopes             -> output/a409/{tests,companions,d409}.csv
    python a409.py --stage power|score [--jobs 8] [--genome-jobs 3]           (long: paper/common/run.sh)
"""
from __future__ import annotations

import argparse
import time
import warnings

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

import a387 as R
import a387_score as SC
from a387 import M, cs, limit_threads
from clocksync import pairwise
from clocksync.elasticity import radius_covariance

OUT = R.OUT / "a409"
BASE = R.BASE
NF, NE, NM = "adjacent_normal:fibroblast", "adjacent_normal:endothelial", "adjacent_normal:mural"
TF, TE, MY = "tumour:fibroblast", "tumour:endothelial", "tumour:myeloid"
NORMAL_ANCH, TUMOUR_ANCH = [NF, NE, NM], [TF, TE]
PAIRS = {"H409.1 tumour fibroblast ~ normal-stroma host": (TF, NORMAL_ANCH),
         "ref: normal fibroblast ~ normal endothelium + mural host": (NF, [NE, NM]),
         "desc: tumour endothelium ~ normal-stroma host": (TE, NORMAL_ANCH)}
N_SHUF, N_SET, N_PERM, SEED, GATE = 2000, 200, 10_000, 409, 0.30
T0 = time.time()


def log(*a):
    print(f"[{(time.time() - T0) / 60:6.1f} min]", *a, flush=True)


def centred_by(z, stratum):
    s = pd.Series(z)
    return (s - s.groupby(np.asarray(stratum)).transform("mean")).to_numpy()


def pair_states(table, target, anchors, wf):
    """(donors, target state, host state, stratum) of the patients with both; states centred within stratum."""
    limit_threads(1)
    warnings.filterwarnings("ignore")
    h = cs.Dataset(table[table.compartment.isin(anchors)], anchors=anchors, waveform=wf, mixed_strata=True).host()
    rows = table[table.compartment == target].drop_duplicates("donor")
    zt = pd.Series(pairwise.states(rows, wf), index=rows.donor.to_numpy())
    ids = zt.index.intersection(h.index)
    st = rows.set_index("donor").stratum.reindex(ids).to_numpy()
    a = centred_by(zt.reindex(ids).to_numpy(), st)
    b = centred_by((h.x + 1j * h.y).reindex(ids).to_numpy(), st)
    return np.asarray(ids), a, b, st


def shuffles(a, st, n, rng):
    """n within-stratum permutations of a (index arrays)."""
    P = np.tile(np.arange(len(a)), (n, 1))
    for s in np.unique(st):
        i = np.flatnonzero(st == s)
        P[:, i] = i[np.argsort(rng.random((n, len(i))), axis=1)]
    return P


def null_lev(a, b, st, n, seed):
    P = shuffles(a, st, n, np.random.default_rng(seed))
    return np.array([radius_covariance(a[p], b, groups=st)["leverage"] for p in P])


def stage_power(args):
    table = R.pooled_table(args.jobs, False)
    ids, a, b, st = pair_states(table, TF, NORMAL_ANCH, BASE)
    nul = null_lev(a, b, st, N_SHUF, [SEED, 0])
    T = pd.DataFrame([dict(pair="H409.1", patients=len(ids), strata=len(set(st)), null_mean=float(nul.mean()), null_sd=float(nul.std()),
                           mde=2.49 * float(nul.std()), gate=GATE, scored=bool(2.49 * nul.std() <= GATE))])
    T.to_csv(OUT / "power.csv", index=False)
    log("power\n", T.round(4).to_string(index=False))


def comp_one(table, wfs):
    out = []
    for fam, j, w in wfs:
        r = dict(family=fam, set=j)
        for lab, (tgt, anch) in PAIRS.items():
            try:
                _, a, b, st = pair_states(table, tgt, anch, w)
                rc = radius_covariance(a, b, groups=st)
                r[lab] = rc["leverage"]
            except Exception:  # noqa: BLE001  (a pseudo waveform without a fit: no value)
                r[lab] = np.nan
        out.append(r)
    return out


def shares(table_all):
    """Per patient: log2 tumour myeloid share of all tumour cells (>= 100), as A246 P8 / A402b."""
    t = table_all[table_all.compartment.str.startswith("tumour:")]
    N = t.groupby("donor").n_cells.sum()
    n = t[t.compartment == MY].groupby("donor").n_cells.sum()
    ok = n.index[(N.reindex(n.index) >= 100).to_numpy()]
    return np.log2((n.reindex(ok) + .5) / (N.reindex(ok) + 1))


def full_tables(studies, jobs):
    def one(s):
        limit_threads(1)
        warnings.filterwarnings("ignore")
        try:
            pb = M.table(s)
        except Exception:  # noqa: BLE001
            return None
        pb = pb[["donor", "compartment", "n_cells"]].copy()
        pb["donor"] = s + "|" + pb.donor.astype(str)
        return pb
    return pd.concat([t for t in Parallel(n_jobs=jobs)(delayed(one)(s) for s in studies) if t is not None], ignore_index=True)


def slope_test(y, x, cohort, stratum, n_perm, seed):
    """Pooled within-cohort slope of z-scored y on z-scored x (stratum-centred), patient weights; permutation null."""
    d = pd.DataFrame(dict(y=y, x=x, c=cohort, s=stratum)).dropna()
    units = []
    for c, g in d.groupby("c"):
        if len(g) < 4:
            continue
        yy, xx = centred_by(g.y.to_numpy(), g.s), centred_by(g.x.to_numpy(), g.s)
        if yy.std() > 0 and xx.std() > 0:
            units.append((len(g), yy / yy.std(), xx / xx.std()))
    if not units:
        return dict(patients=0)
    w = np.array([u[0] for u in units], float)
    obs = float(sum(n * np.mean(yy * xx) for n, yy, xx in units) / w.sum())
    rng = np.random.default_rng(seed)
    nul = np.zeros(n_perm)
    for n, yy, xx in units:
        P = np.argsort(rng.random((n_perm, n)), axis=1)
        nul += n * (yy[P] * xx[None, :]).mean(1)
    nul /= w.sum()
    return dict(patients=int(w.sum()), cohorts=len(units), slope=obs, null_sd=float(nul.std()),
                p_two_sided=float((np.sum(np.abs(nul) >= abs(obs) - 1e-12) + 1) / (n_perm + 1)))


def stage_score(args):
    pw = pd.read_csv(OUT / "power.csv")
    scored = bool(pw.scored.iloc[0])
    table = R.pooled_table(args.jobs, False)
    rows = []
    for k, (lab, (tgt, anch)) in enumerate(PAIRS.items()):
        ids, a, b, st = pair_states(table, tgt, anch, BASE)
        rc = radius_covariance(a, b, groups=st)
        nul = null_lev(a, b, st, N_SHUF, [SEED, k])
        rows.append(dict(pair=lab, patients=len(ids), leverage=rc["leverage"], cv2=rc["cv2"], null_mean=float(nul.mean()),
                         null_sd=float(nul.std()), p_one_sided=float((np.sum(nul >= rc["leverage"]) + 1) / (N_SHUF + 1)),
                         hypothesis=("H409.1" if scored else "descriptive (power)") if k == 0 else "descriptive"))
    log("observed\n", pd.DataFrame(rows).round(4).to_string(index=False))
    # companions on genome-wide pseudobulks (A387 part 2 construction)
    keep = set(zip(table.donor, table.compartment))
    studies = sorted({d.split("|", 1)[0] for d in table.donor})
    Gs = Parallel(n_jobs=args.genome_jobs)(delayed(SC.genome_cohort)(s, {x for x in keep if x[0].startswith(s + "|")}) for s in studies)
    G = pd.concat([g for g in Gs if g is not None and len(g)], ignore_index=True)
    G = G.drop(columns=["stratum"]).merge(table[["donor", "compartment", "stratum"]], on=["donor", "compartment"])
    genes = [c for c in G.columns if c not in ("donor", "compartment", "stratum", "depth", "n_cells")]
    genes = [g for g in genes if G[g].notna().all()]
    ref = G[G.compartment.isin(NORMAL_ANCH)]
    frac = ref[genes].div(ref.depth, axis=0)
    lm = np.log(frac.mean() + 1e-12)[(ref[genes].mean() >= 1) | frac.columns.isin(list(BASE.genes))]
    from clocksync.paired import pseudo_waveforms
    wfs = [(f, j, w) for code, f in enumerate(("matched", "responsive"))
           for j, w in enumerate(pseudo_waveforms(lm, N_SET, np.random.default_rng([SEED, code]), BASE, f))]
    cols = ["donor", "compartment", "stratum", "depth", "n_cells"]
    chunks = [wfs[i::args.jobs] for i in range(args.jobs)]
    res = Parallel(n_jobs=args.jobs)(delayed(comp_one)(G[cols + sorted({g for _, _, w in ch for g in w.genes} | set(BASE.genes))], ch)
                                     for ch in chunks)
    Cp = pd.DataFrame([r for x in res for r in x])
    Cp.to_csv(OUT / "companions.csv", index=False)
    for r in rows:
        for fam in ("matched", "responsive"):
            v = Cp[Cp.family == fam][r["pair"]].dropna().to_numpy()
            r.update({f"{fam}_n": len(v), f"{fam}_median": float(np.median(v)), f"{fam}_q95": float(np.quantile(v, .95)),
                      f"{fam}_pct": float(100 * np.mean(v < r["leverage"]))})
    T = pd.DataFrame(rows)
    T["pass"] = ((T.p_one_sided < 0.05) & (T.leverage > T.matched_q95)).where(T.hypothesis == "H409.1")
    T.to_csv(OUT / "tests.csv", index=False)
    log("tests\n", T.round(4).to_string(index=False))
    # D409: myeloid share on normal-stroma vs tumour-stroma loudness (same patients)
    full = full_tables(studies, args.jobs)
    y = shares(full)
    hn = cs.Dataset(table[table.compartment.isin(NORMAL_ANCH)], anchors=NORMAL_ANCH, waveform=BASE, mixed_strata=True).host()
    ht = cs.Dataset(table[table.compartment.isin(TUMOUR_ANCH)], anchors=TUMOUR_ANCH, waveform=BASE, mixed_strata=True).host()
    st = table.drop_duplicates("donor").set_index("donor").stratum
    ids = y.index.intersection(hn.index).intersection(ht.index)
    s_ = st.reindex(ids).to_numpy()
    ln = centred_by(np.log(hn.radius.reindex(ids).to_numpy()), s_)
    lt = centred_by(np.log(ht.radius.reindex(ids).to_numpy()), s_)
    coh = np.array([d.split("|")[0] for d in ids])
    D = []
    for lab, x in (("normal-stroma loudness", ln), ("tumour-stroma loudness", lt)):
        D.append(dict(loudness=lab, **slope_test(y.reindex(ids).to_numpy(), x, coh, s_, N_PERM, [SEED, 9, len(lab)])))
    D.append(dict(loudness="corr(normal, tumour loudness), within stratum", slope=float(np.corrcoef(ln, lt)[0, 1]), patients=len(ids)))
    D = pd.DataFrame(D)
    D.to_csv(OUT / "d409.csv", index=False)
    log("D409\n", D.round(4).to_string(index=False))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stage", choices=("power", "score"), required=True)
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--genome-jobs", type=int, default=3)
    args = ap.parse_args()
    limit_threads(1)
    OUT.mkdir(parents=True, exist_ok=True)
    (stage_power if args.stage == "power" else stage_score)(args)
    log("done")


if __name__ == "__main__":
    main()
