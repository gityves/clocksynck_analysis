"""A146 part 2 (pan-cancer) scoring: H146.3 omnibus S, within-cohort permutation null, Holm over H146.1-3, D146.2.

Reuses part2pc_prep.py unchanged (units_of, omnibus_null, UNITS/EXTRA, HUMAN_H, keeping, SEED) and its caches. Pool = gate-passing
cohorts with >= 6 timed patients (as in the prep). Statistic per program: slope t of stacked OLS of the within-cohort z score on the
cohort-centred cos(2 pi (t - h_k)/24), df = n - cohorts - 1; S = mean t_k. Null: times permuted within cohort (10,000, seed 146),
one permutation per cohort and replicate shared by all programs; one-sided p = (1 + #{S0 >= S}) / (1 + n_perm).
D146.2: per program and per cohort t (df n - 2), omnibus11, keeping / non-keeping / breast subsets, offset sweep (+-2, 4, 6 h).
Subset and sweep p-values are descriptive (same permutation machinery, same seed).

    python score_part2.py [--n-perm 10000]
Outputs: output/part2_{tests,per_program,per_cohort,subsets,offset_sweep}.csv
"""
from __future__ import annotations
import sys, argparse
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import t as tdist

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import part2pc_prep as P  # noqa: E402
q = P.q
OUT = P.BASE
W = 2 * np.pi / 24
BASE_H = dict(P.HUMAN_H)


def build_pools():
    H = pd.read_csv(P.CACHE / "host_cohorts.csv")
    T = pd.read_csv(P.CACHE / "host_times.csv", dtype={"patient": str})
    pools, meta = {}, {}
    for r in H.itertuples():
        if not (r.gate_passes is True or r.gate_passes == True):  # noqa: E712
            continue
        timed = T[T.study == r.study]
        tmap = dict(zip(timed.patient, timed.internal_h))
        if len(tmap) < P.MIN_UNIT:
            continue
        pb, mal = q.study_table(r.study)
        tum = pb[pb.compartment.str.startswith("tumour:")]
        c, tot, comp = P.cohort_inputs(r.study, None, mal)
        pbc = {"mal": set(tum[tum.compartment == mal].donor.astype(str)),
               "T_NK": set(tum[tum.compartment == "tumour:T_NK"].donor.astype(str)),
               "myeloid": set(tum[tum.compartment == "tumour:myeloid"].donor.astype(str))}
        ids = sorted(tmap)
        pools[r.study] = (np.array([tmap[p] for p in ids]), ids, P.units_of(r.study, c, tot, comp, pbc, ids))
        cancer = c.cancer.iloc[0]
        tcga, conf, cls = P.keeping(cancer)
        meta[r.study] = dict(cancer=cancer, tcga=tcga, cls=cls, breast=("breast" in cancer) or r.study == "SRP329496")
    return pools, meta


def to_co(pools, studies, progs):
    co = {}
    for s in studies:
        t, ids, un = pools[s]
        ix = {p: i for i, p in enumerate(ids)}
        co[s] = (t, {pr: (np.array([ix[p] for p in sc.index]), sc.to_numpy()) for pr, sc in un.items() if pr in progs})
    return co


def hof(pr):
    return BASE_H["P9"] if pr == "P9m" else BASE_H[pr]


def tstat(co, pr, shift=0.0):
    """Observed stacked slope t for program pr; returns (t, n, n_cohorts, r_partial)."""
    sxy = sxx = syy = 0.0; n = nc = 0
    for s, (t, u) in co.items():
        if pr not in u:
            continue
        idx, y = u[pr]
        cs = np.cos(W * (t[idx] - hof(pr) - shift)); cs = cs - cs.mean()
        sxy += cs @ y; sxx += (cs ** 2).sum(); syy += (y ** 2).sum(); n += len(y); nc += 1
    if n - nc - 1 < 3:
        return np.nan, n, nc, np.nan
    b = sxy / sxx
    s2 = (syy - b * sxy) / (n - nc - 1)
    tt = b / np.sqrt(s2 / sxx)
    return float(tt), n, nc, float(tt / np.sqrt(tt ** 2 + n - nc - 1))


def cohort_t(t, idx, y, pr):
    cs = np.cos(W * (t[idx] - hof(pr))); cs = cs - cs.mean()
    sxx = (cs ** 2).sum(); b = cs @ y / sxx
    yc = y - y.mean()
    s2 = ((yc - b * cs) ** 2).sum() / (len(y) - 2)
    return float(b / np.sqrt(s2 / sxx))


def run_S(co, progs, n_perm, shift=0.0):
    used = [p for p in progs if any(p in v[1] for v in co.values())]
    obs = [tstat(co, p, shift)[0] for p in used]
    keep = [p for p, o in zip(used, obs) if np.isfinite(o)]
    S0 = float(np.mean([o for o in obs if np.isfinite(o)]))
    for k in P.HUMAN_H:
        P.HUMAN_H[k] = BASE_H[k] + shift
    try:
        null = P.omnibus_null(co, keep, n_perm, P.SEED)
    finally:
        P.HUMAN_H.update(BASE_H)
    p = (1 + np.sum(null >= S0)) / (1 + n_perm)
    return dict(S=S0, p=float(p), null_mean=float(null.mean()), null_sd=float(null.std(ddof=1)),
                null_q95=float(np.quantile(null, .95)), MDE=float(np.quantile(null, .95) + .84 * null.std(ddof=1)),
                n_programs=len(keep), n_cohorts=len(co), n_patients_scored=len({(s, i) for s, v in co.items() for p in keep if p in v[1] for i in v[1][p][0]}))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("--n-perm", type=int, default=10000); a = ap.parse_args()
    NP = a.n_perm
    pools, meta = build_pools()
    allc = list(pools)
    prog10 = [u[0] for u in P.UNITS]; prog11 = prog10 + ["P9m"]
    print("pool:", allc, "timed patients:", sum(len(v[1]) for v in pools.values()), flush=True)

    # H146.3 and tests
    co = to_co(pools, allc, prog11)
    main = run_S(co, prog10, NP)
    print("H146.3", main, flush=True)
    p1 = pd.read_csv(OUT / "part1_tests.csv").set_index("hypothesis")
    raw = {"H146.1": float(p1.loc["H146.1", "p"]), "H146.2": float(p1.loc["H146.2", "p"]), "H146.3": main["p"]}
    order = sorted(raw, key=raw.get); holm = {}; run = 0.0
    for i, h in enumerate(order):
        run = max(run, min(1.0, (len(order) - i) * raw[h])); holm[h] = run
    rows = []
    for h in ("H146.1", "H146.2"):
        r = p1.loc[h]
        rows.append(dict(hypothesis=h, time=r.time, n_patients=r.n_patients, n_programs=r.n_programs, n_cohorts=np.nan, S=r.S, p=r.p,
                         null_mean=r.null_mean, null_sd=r.null_sd, null_q95=r.null_q95, MDE=r.MDE, p_holm=holm[h],
                         verdict=("pass" if holm[h] < .05 and r.S > 0 else "fail"), source="part1_tests.csv (raw p; Holm recomputed over 3)"))
    rows.append(dict(hypothesis="H146.3", time="internal", n_patients=main["n_patients_scored"], n_programs=main["n_programs"],
                     n_cohorts=main["n_cohorts"], S=main["S"], p=main["p"], null_mean=main["null_mean"], null_sd=main["null_sd"],
                     null_q95=main["null_q95"], MDE=main["MDE"], p_holm=holm["H146.3"],
                     verdict=("pass" if holm["H146.3"] < .05 and main["S"] > 0 else "fail (not detected at this design)" if main["S"] < main["MDE"] else "fail"),
                     source="score_part2.py", n_perm=NP, seed=P.SEED))
    pd.DataFrame(rows).to_csv(OUT / "part2_tests.csv", index=False)

    # per program, per cohort
    pp = []
    for pr in prog11:
        tt, n, nc, r = tstat(co, pr)
        pp.append(dict(program=pr, compartment=dict((u[0], u[1]) for u in P.UNITS + P.EXTRA)[pr], predicted_h=hof(pr), t=tt, partial_r=r,
                       n_patients=n, n_cohorts=nc, in_omnibus10=pr in prog10))
    pd.DataFrame(pp).to_csv(OUT / "part2_per_program.csv", index=False)
    pc = []
    for s, (t, ids, un) in pools.items():
        for pr in prog11:
            if pr in un:
                idx = np.array([ids.index(p) for p in un[pr].index]); y = un[pr].to_numpy()
                pc.append(dict(study=s, cancer=meta[s]["cancer"], clock_class=meta[s]["cls"], breast=meta[s]["breast"], program=pr,
                               n=len(y), t=cohort_t(t, idx, y, pr)))
    pd.DataFrame(pc).to_csv(OUT / "part2_per_cohort.csv", index=False)

    # subsets (descriptive)
    sub_rows = []
    sets = {"omnibus10": (allc, prog10), "omnibus11": (allc, prog11),
            "keeping": ([s for s in allc if meta[s]["cls"] == "keeping"], prog10),
            "not_keeping": ([s for s in allc if meta[s]["cls"] == "not_keeping"], prog10),
            "breast": ([s for s in allc if meta[s]["breast"] and s in ("ERP130387", "SRP273482")], prog10)}
    for name, (ss, pg) in sets.items():
        if not ss:
            continue
        r = run_S(to_co(pools, ss, pg), pg, NP)
        sub_rows.append(dict(subset=name, cohorts=";".join(ss), n_timed_patients=sum(len(pools[s][1]) for s in ss), **r))
        print(name, r, flush=True)
    pd.DataFrame(sub_rows).to_csv(OUT / "part2_subsets.csv", index=False)

    # offset sweep
    sw = []
    for sh in (-6, -4, -2, 0, 2, 4, 6):
        r = run_S(co, prog10, NP, shift=float(sh)); sw.append(dict(offset_h=sh, **r)); print("offset", sh, r["S"], r["p"], flush=True)
    pd.DataFrame(sw).to_csv(OUT / "part2_offset_sweep.csv", index=False)
