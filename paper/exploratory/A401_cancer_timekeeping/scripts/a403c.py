"""A403c (POST HOC, DESCRIPTIVE CHECK; registered in paper/docs/registered_experiments.md, "A403c"): is the colder-louder lead
of A403 a phase effect seen through the host radius?

Same pool, hosts, loudness, readouts and covariates as A403 / A403b.
  C1  per cancer R^2 of stratum-centred loudness on cos / sin of phi and 2 phi (phi = host phase); the same for the 200
      matched and 200 responsive pseudo-hosts of each cancer (A403 sets)            -> output/a403c/c1.csv
  C2  A403b statistics (pooled and per cancer) with loudness replaced by its C1 residual (pseudo-hosts likewise)
                                                                                     -> output/a403c/c2.csv
    python a403c.py [--jobs 12] [--genome-jobs 4]                       (long: paper/common/run.sh)
"""
from __future__ import annotations

import argparse
import time
import warnings

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

import a403_a405_score as S
import a403b as B
from a402_a405 import BASE, META, U, cs, log, trim, centred, cancer_hosts, loudness_frame, load_slim, readouts
from clocksync.config import limit_threads

OUT = U.OUT / "a403c"
N_SET = 200


def harmonics(phi, stratum):
    X = np.column_stack([np.cos(phi), np.sin(phi), np.cos(2 * phi), np.sin(2 * phi)])
    return np.column_stack([centred(X[:, j], stratum) for j in range(4)])


def phase_free(loud, phi, stratum, cancer):
    """Per cancer: residual of the stratum-centred loudness on the centred phase harmonics, and R^2 per cancer."""
    res = pd.Series(np.nan, index=loud.index)
    r2 = {}
    for c in pd.unique(cancer):
        m = (cancer == c) & np.isfinite(loud.to_numpy()) & np.isfinite(phi)
        if m.sum() < 10:
            continue
        st = pd.Series(stratum[m])
        y = centred(loud.to_numpy()[m], st)
        X = harmonics(phi[m], st)
        b = np.linalg.lstsq(X, y, rcond=None)[0]
        e = y - X @ b
        res.iloc[np.flatnonzero(m)] = e
        r2[c] = (float(1 - np.sum(e ** 2) / np.sum(y ** 2)) if np.sum(y ** 2) > 0 else np.nan, int(m.sum()))
    return res, r2


def pseudo_states(tab, classes, w):
    limit_threads(1)
    warnings.filterwarnings("ignore")
    try:
        h = cs.Dataset(tab, anchors=classes, waveform=w, mixed_strata=True).host()
        return pd.DataFrame({"z": (h.x + 1j * h.y).to_numpy(), "radius": h.radius.to_numpy()}, index=h.index)
    except Exception:  # noqa: BLE001  (a pseudo-host that cannot be fitted: no value for that set)
        return pd.DataFrame(columns=["z", "radius"])


def companion_states(hosts, args):
    """{family: [DataFrame donor -> pseudo-host z, radius]}; the A403 sets (seeds [401, 0] / [401, 1], per cancer)."""
    from clocksync.paired import pseudo_waveforms
    out = {"matched": [[] for _ in range(N_SET)], "responsive": [[] for _ in range(N_SET)]}
    for c, h in hosts.items():
        Gc = h["Gc"]
        genes = [g for g in Gc.columns if g not in META]
        frac = Gc[genes].div(Gc.depth, axis=0)
        lm = np.log(frac.mean() + 1e-12)[(Gc[genes].mean() >= 1) | frac.columns.isin(list(BASE.genes))]
        for code, fam in enumerate(("matched", "responsive")):
            ws = pseudo_waveforms(lm, N_SET, np.random.default_rng([U.SEED, code]), BASE, fam)
            res = Parallel(n_jobs=args.jobs)(delayed(pseudo_states)(Gc[META + sorted(set(w.genes) | set(BASE.genes))], h["classes"], w)
                                             for w in ws)
            for j, r in enumerate(res):
                out[fam][j].append(r)
        log(f"companion states {c}")
    return {f: [pd.concat(v) if v else pd.DataFrame(columns=["z", "radius"]) for v in sets] for f, sets in out.items()}


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--jobs", type=int, default=12)
    ap.add_argument("--genome-jobs", type=int, default=4)
    args = ap.parse_args()
    args.smoke = False
    limit_threads(1)
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    frames, donors, f388 = U.pool_sources(args, OUT)
    from joblib.externals.loky import get_reusable_executor
    get_reusable_executor().shutdown(wait=True)
    p401 = pd.read_csv(U.OUT / "power" / "p401.csv")
    cancers = [c for c in p401[(p401.readout.astype(str) == "14") & p401.p401_pass].cancer if c in set(donors.cancer) and c != "unspecified"]
    donors = donors[donors.cancer.isin(cancers)].reset_index(drop=True)
    hosts = cancer_hosts(frames, donors, cancers)
    sources, f388 = list(frames), set(f388)
    frames.clear()
    trim()
    loud = loudness_frame(hosts)
    assert len(hosts) == 12 and len(loud) == 574, "pool differs from A403"
    slim, errs = load_slim(sources, donors, f388, args)
    assert len(slim) == 51 and errs.empty
    ro = readouts(slim, loud)
    df = loud.merge(ro.drop(columns="cohort"), on="donor").merge(S.covariates(hosts), on="donor", how="left")
    df["cohort"] = df.donor.str.split("|").str[0]
    rep = B.slope(df, "immune_hot", B.BASECOV)[0]
    a403 = pd.read_csv(U.OUT / "score_a403_a405" / "a403_tests.csv").query("adjusted and readout == 'immune_hot'").slope.iloc[0]
    assert abs(rep - a403) < 1e-9, (rep, a403)
    log(f"A403 reproduced ({rep:.4f}); table {len(df)} patients")

    # C1: observed
    phi = np.angle(df.z.to_numpy())
    df["loud_pf"], r2 = phase_free(df.loud, phi, df.stratum.to_numpy(), df.cancer.to_numpy())
    Z = companion_states(hosts, args)
    get_reusable_executor().shutdown(wait=True)
    st_of = loud.set_index("donor").stratum
    ca_of = loud.set_index("donor").cancer
    comp_raw, comp_pf, r2c = {}, {}, {}
    for fam, sets in Z.items():
        comp_raw[fam], comp_pf[fam], r2c[fam] = [], [], []
        for zd in sets:
            zd = zd[np.isfinite(zd.radius.astype(float).to_numpy()) & (zd.radius.astype(float).to_numpy() > 0)]
            zs = zd.z
            if zd.empty:
                comp_raw[fam].append(pd.Series(dtype=float)), comp_pf[fam].append(pd.Series(dtype=float)), r2c[fam].append({})
                continue
            st = st_of.reindex(zs.index)
            lr = pd.Series(centred(np.log(zd.radius.astype(float).to_numpy()), st), index=zd.index)
            pf, rr = phase_free(lr, np.angle(zs.to_numpy().astype(complex)), st.to_numpy(), ca_of.reindex(zs.index).to_numpy())
            comp_raw[fam].append(lr)
            comp_pf[fam].append(pf.dropna())
            r2c[fam].append(rr)
    rows = []
    cancers_ = sorted(r2)
    for c in cancers_ + ["pooled"]:
        if c == "pooled":
            obs = sum(r2[k][0] * r2[k][1] for k in r2) / sum(r2[k][1] for k in r2)
            n = sum(r2[k][1] for k in r2)
        else:
            obs, n = r2[c]
        row = dict(cancer=c, patients=n, r2=obs)
        for fam, lst in r2c.items():
            if c == "pooled":
                v = np.array([sum(d[k][0] * d[k][1] for k in d) / sum(d[k][1] for k in d) for d in lst if d])
            else:
                v = np.array([d[c][0] for d in lst if c in d])
            v = v[np.isfinite(v)]
            row.update({f"{fam}_n": len(v), f"{fam}_median": float(np.median(v)) if len(v) else np.nan,
                        f"{fam}_q95": float(np.quantile(v, .95)) if len(v) else np.nan,
                        f"{fam}_percentile": float(100 * np.mean(v < obs)) if len(v) else np.nan})
        rows.append(row)
    C1 = pd.DataFrame(rows)
    C1.to_csv(OUT / "c1.csv", index=False)
    pd.set_option("display.width", 300, "display.max_columns", 40)
    log("C1\n" + C1.round(4).to_string(index=False))

    # C2: A403b statistics on phase-free loudness (and on raw loudness, same code, as the reference)
    C2 = []
    for lab, col, comp in (("raw", "loud", comp_raw), ("phase_free", "loud_pf", comp_pf)):
        for rd in S.ROUT:
            for r in B.per_cancer(df, rd, B.BASECOV, comp, f"c2{lab}", loud_col=col):
                C2.append(dict(r, loudness=lab))
        log(f"C2 {lab} done")
    C2 = pd.DataFrame(C2)
    C2.to_csv(OUT / "c2.csv", index=False)
    log("C2 pooled\n" + C2[C2.cancer == "pooled"].round(4).to_string(index=False))
    log(f"done {(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    main()
