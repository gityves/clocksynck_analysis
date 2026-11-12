"""A248: glucocorticoid-driven T-cell rhythm in GTEx blood and tissues (follow-up of A245).

Registered in paper/docs/registered_experiments.md (A248). Donors, covariates, strata, time axes and statistics of A245
(`compute_a245`); whole blood is added as its own compartment. Readouts (mean log2(CPM + 1)): T-cell content, glucocorticoid
response GC, naive and cytotoxic effector scores, a neutrophil score for one sensitivity.

Tests (recorded time; S = mean over the compartments of the OLS t at the predicted hour, permutation of time within Hardy x sex):
  H248.1 blood T-cell content at 02:10 h      H248.2 GC of the 26 tissues at 10:40 h      H248.4 tissue T-cell content at 14:00 h
  H248.3 blood: partial correlation of T-cell content with GC given the covariates (Fisher z), < 0; null: GC permuted across
         donors within Hardy x sex; companion: GC replaced by matched 5-gene sets.
Stages:
  power  null SD and MDE_80 (2.49 x SD) of every test from permuted times / permuted GC only; no observed statistic
                                                                            -> output/{power,compartments}.csv
  score  (after the power note is dated) observed statistics, p, Holm, matched companions, descriptives
                                                 -> output/score_{tests,matched,descriptive,per_compartment}.csv
  `--placebo`: times and GC order shuffled once over donors (code test; output/placebo/).

    python compute_a248.py --stage power|score [--subsample N] [--jobs 4]
Controlled access: time of death, exact age and ischaemic time stay in memory; only aggregates are written.

Choices not fixed by the registration
- Donors: the 233 donors of A245 (fast deaths, time of death, body host); blood enters with those that have a blood sample.
- Hours: 02:10 = 2.167 h, 10:40 = 10.667 h; naive score against 02:25 h (the mean of 02:09 and 02:41), effector 15:34 h.
- GC permutation: one random key per donor and iteration, ordered within stratum among the donors of the compartment (the
  same keys for all compartments).
- Matched companions as A245 (one gene per readout gene from its 400 closest in mean log2(CPM + 1) over the compartments of
  the test; clock, readout and endothelial-score genes excluded).
- Descriptive D = S(internal) - S(recorded) is given without an interval; ventilated donors (Hardy 0, own host) are read
  against recorded time with 2,000 permutations.
- Holm: primaries H248.1-2; secondaries H248.3-4.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2]))
sys.path.insert(0, str(HERE.parents[1] / "A245_trafficking_gtex" / "scripts"))
import clocksync as cs  # noqa: E402
from common import cli, gtex  # noqa: E402

import compute_a245 as M  # noqa: E402

OUT = HERE.parent / "output"
BLOOD = "Whole Blood"
SETS = {"T_cell": M.TCELL, "GC": ["TSC22D3", "FKBP5", "DDIT4", "ZBTB16", "KLF9"], "naive": ["CCR7", "SELL", "LEF1", "TCF7"],
        "effector": ["GZMB", "PRF1", "CX3CR1", "FGFBP2"], "neutrophil": ["CSF3R", "FCGR3B", "CXCR2"]}
# test -> (hypothesis, role, readout, compartments, predicted hour)
TIME_TESTS = {"blood_T": ("H248.1", "primary", "T_cell", "blood", 2 + 10 / 60), "tissue_GC": ("H248.2", "primary", "GC", "tissues", 10 + 40 / 60),
              "tissue_T": ("H248.4", "secondary, post hoc to A245", "T_cell", "tissues", 14.0)}
DESC_TIME = {"blood_GC": ("GC", "blood", 10 + 40 / 60), "blood_naive": ("naive", "blood", 2 + 25 / 60),
             "blood_effector": ("effector", "blood", 15 + 34 / 60)}
N_PERM, N_MATCH, N_POOL, MDE_K = 10_000, 200, 400, 2.49
W = M.W


def keep(t):
    return M.target(t) or t == BLOOD


def prepare(hardy, args, quick, placebo):
    L, genes, sm, mean_count = M.load_lcpm(hardy, keep)
    donors, times, _ = M.hosts(hardy, sorted(sm.tissue.unique()))
    if quick:
        donors = donors.loc[np.random.default_rng(args.seed).permutation(sorted(donors.index))[:args.subsample]]
    units, col, codes = M.design(L, genes, sm, donors, times)
    if placebo:
        perm = np.random.default_rng(M.seed_of("A248placebo", args.seed)).permutation(len(donors))
        donors = donors.assign(tod=donors.tod.to_numpy()[perm])
        times = {k: pd.Series(v.reindex(donors.index).to_numpy()[perm], index=donors.index) for k, v in times.items()}
    M.attach_full_time(units, donors, times)
    groups = {"blood": {k: u for k, u in units.items() if k == BLOOD}, "tissues": {k: u for k, u in units.items() if k != BLOOD}}
    return L, genes, mean_count, donors, groups, col, codes


def score_of(L, col, u, genes_, resid_on=None):
    y = L[np.ix_(u["rows"], [col[g] for g in genes_])].astype(float).mean(1)
    if resid_on:
        e = L[np.ix_(u["rows"], [col[g] for g in resid_on])].astype(float).mean(1)
        A = np.column_stack([np.ones(len(e)), e])
        y = y - A @ np.linalg.lstsq(A, y, rcond=None)[0]
    return y


def S_time(units, Y, h, tod, axis="recorded", P=None):
    """Mean over compartments of the t at hour h. Y: {compartment: samples x readouts}; returns perms (or 1) x readouts."""
    acc = 0
    for k, u in units.items():
        full = tod if axis == "recorded" else u["time_full"]
        tm = full[u["pos"]][None, :] if P is None else full[P][:, u["pos"]]
        acc = acc + M.t_cos(u["Q"], M.resid(u["Q"], Y[k]), np.cos(W * (tm - h)))
    return acc / len(units)


def cosinor(units, Y, tod, axis="recorded"):
    rows = []
    for k, u in units.items():
        tm = (tod if axis == "recorded" else u["time_full"])[u["pos"]]
        z, amp, pk = M.t_cosinor(u["Q"], M.resid(u["Q"], Y[k]), tm)
        rows.append((k, len(tm), amp[0], pk[0], z[0]))
    return rows


def partial_z(units, Ya, Yb, order=None, with_time=None):
    """Mean over compartments of Fisher z of the partial correlation of Ya (1 column) with each column of Yb given the covariates.
    order: perms x donors random keys -> Ya... is fixed and Yb permuted within stratum among the compartment's donors."""
    acc = 0
    for k, u in units.items():
        Q = u["Q"]
        if with_time is not None:
            tm = with_time[u["pos"]]
            Q = np.linalg.qr(np.column_stack([u["C"], np.cos(W * tm), np.sin(W * tm)]))[0]
        a = M.resid(Q, Ya[k])[:, 0]
        a = a / np.sqrt((a ** 2).sum())
        if order is None:
            b = M.resid(Q, Yb[k])
            r = a @ (b / np.sqrt((b ** 2).sum(0)))
        else:
            st, key = u["codes"], order[:, u["pos"]]
            idx = np.tile(np.arange(len(a)), (len(key), 1))
            for c in np.unique(st):
                i = np.flatnonzero(st == c)
                idx[:, i] = i[np.argsort(key[:, i], axis=1)]
            b = Yb[k][:, 0][idx]                                   # perms x samples
            b = b - (b @ Q) @ Q.T
            r = (b @ a) / np.sqrt((b ** 2).sum(1))
        acc = acc + np.arctanh(np.clip(r, -0.999999, 0.999999))
    return acc / len(units)


def matched(genes, col, mean_l, ok, names, rng):
    cand = np.flatnonzero(ok)
    picks = [rng.choice(cand[np.argsort(np.abs(mean_l[cand] - mean_l[col[g]]))[:N_POOL]], N_MATCH, replace=False) for g in names]
    return [[genes[p[i]] for p in picks] for i in range(N_MATCH)]


def run(args, quick, stage):
    L, genes, mean_count, donors, groups, col, codes = prepare(gtex.FAST, args, quick, args.placebo)
    for g in groups.values():
        for u in g.values():
            u["codes"] = codes[u["pos"]]
    tod = donors.tod.to_numpy(float)
    n_perm = 200 if quick else N_PERM
    P = M.permutations(codes, n_perm, M.seed_of("A248null", args.seed))
    keys = np.random.default_rng(M.seed_of("A248gc", args.seed)).random((n_perm, len(donors)))
    Y = {name: {grp: {k: score_of(L, col, u, SETS[name])[:, None] for k, u in g.items()} for grp, g in groups.items()} for name in SETS}
    if args.placebo:                                               # GC order shuffled once (the time shuffle is in prepare)
        rng = np.random.default_rng(M.seed_of("A248placebo_gc", args.seed))
        for grp in Y["GC"]:
            Y["GC"][grp] = {k: v[rng.permutation(len(v))] for k, v in Y["GC"][grp].items()}
    pd.DataFrame([dict(group=grp, compartment=k, donors=len(u["rows"])) for grp, g in groups.items() for k, u in g.items()]
                 ).to_csv(OUT / "compartments.csv", index=False)

    nul_t = {t: S_time(groups[grp], Y[r][grp], h, tod, P=P)[:, 0] for t, (_, _, r, grp, h) in TIME_TESTS.items()}
    nul_gc = partial_z(groups["blood"], Y["T_cell"]["blood"], Y["GC"]["blood"], order=keys)
    rows = []
    for t, (hyp, role, r, grp, h) in TIME_TESTS.items():
        se = []
        for k, u in groups[grp].items():
            yr = M.resid(u["Q"], Y[r][grp][k])[:, 0]
            x = np.cos(W * (tod[u["pos"]] - h))
            xr = x - u["Q"] @ (u["Q"].T @ x)
            se.append(np.sqrt((yr ** 2).sum() / (len(yr) - u["Q"].shape[1] - 1) / (xr ** 2).sum()))
        sd = float(nul_t[t].std())
        rows.append(dict(test=t, hypothesis=hyp, role=role, readout=r, group=grp, compartments=len(groups[grp]),
                         donors=int(max(len(u["rows"]) for u in groups[grp].values())), predicted_peak_h=h, null_sd=sd, mde=MDE_K * sd,
                         mde_amplitude_log2_median=MDE_K * sd * float(np.median(se))))
    sd = float(nul_gc.std())
    rows.append(dict(test="blood_T_vs_GC", hypothesis="H248.3", role="secondary", readout="T_cell ~ GC", group="blood", compartments=1,
                     donors=len(groups["blood"][BLOOD]["rows"]), null_sd=sd, mde=MDE_K * sd, mde_partial_r=float(np.tanh(MDE_K * sd))))
    pw = pd.DataFrame(rows)
    if stage == "power":
        pw.to_csv(OUT / "power.csv", index=False)
        print(pw.round(4).to_string(), flush=True)
        return

    excl = set(cs.gtex().genes) | set(M.ENDO) | {"HBB"} | {g for v in SETS.values() for g in v}
    ok = (mean_count >= 1) & ~np.isin(genes, list(excl))
    tests, mrows = [], []
    for t, (hyp, role, r, grp, h) in TIME_TESTS.items():
        g = groups[grp]
        s = float(S_time(g, Y[r][grp], h, tod)[0, 0])
        mean_l = np.mean([L[u["rows"]].mean(0) for u in g.values()], axis=0)
        sets = matched(genes, col, mean_l, ok, SETS[r], np.random.default_rng(M.seed_of("A248match", t)))
        Ym = {k: np.column_stack([score_of(L, col, u, x) for x in sets]) for k, u in g.items()}
        sm_ = S_time(g, Ym, h, tod)[0]
        zc = np.mean([c[4][0] for c in cosinor(g, Y[r][grp], tod)])
        Rn = np.abs(np.mean([M.t_cosinor(u["Q"], M.resid(u["Q"], Y[r][grp][k]), tod[P][:, u["pos"]])[0][:, 0] for k, u in g.items()], axis=0))
        tests.append(dict(test=t, hypothesis=hyp, role=role, readout=r, group=grp, compartments=len(g), predicted_peak_h=h, statistic=s,
                          null_sd=float(nul_t[t].std()), p_one_sided=(1 + np.sum(nul_t[t] >= s)) / (1 + n_perm), matched_q95=float(np.quantile(sm_, .95)),
                          above_matched=bool(s > np.quantile(sm_, .95)), fitted_peak_h=float((np.angle(zc) / W) % 24), R=float(abs(zc)),
                          p_R=(1 + np.sum(Rn >= abs(zc))) / (1 + n_perm), S_internal=float(S_time(g, Y[r][grp], h, tod, "internal")[0, 0])))
        mrows += [dict(test=t, set=i, genes=";".join(x), statistic=float(v)) for i, (x, v) in enumerate(zip(sets, sm_))]
    g = groups["blood"]
    z = float(partial_z(g, Y["T_cell"]["blood"], Y["GC"]["blood"])[0])
    mean_l = L[g[BLOOD]["rows"]].mean(0)
    sets = matched(genes, col, mean_l, ok, SETS["GC"], np.random.default_rng(M.seed_of("A248match", "gc")))
    zm = partial_z(g, Y["T_cell"]["blood"], {BLOOD: np.column_stack([score_of(L, col, g[BLOOD], x) for x in sets])})
    tests.append(dict(test="blood_T_vs_GC", hypothesis="H248.3", role="secondary", readout="T_cell ~ GC", group="blood", compartments=1,
                      statistic=z, partial_r=float(np.tanh(z)), null_sd=float(nul_gc.std()), p_one_sided=(1 + np.sum(nul_gc <= z)) / (1 + n_perm),
                      matched_q95=float(np.quantile(zm, .05)), above_matched=bool(z < np.quantile(zm, .05))))
    mrows += [dict(test="blood_T_vs_GC", set=i, genes=";".join(x), statistic=float(v)) for i, (x, v) in enumerate(zip(sets, zm))]
    T = pd.DataFrame(tests)
    T["D_internal_minus_recorded"] = T.S_internal - T.statistic
    for fam in (["H248.1", "H248.2"], ["H248.3", "H248.4"]):
        m = T.hypothesis.isin(fam).to_numpy()
        p = T.loc[m, "p_one_sided"].to_numpy()
        o = np.argsort(p)
        T.loc[m, "p_holm"] = np.maximum.accumulate(np.minimum(1, p[o] * (len(p) - np.arange(len(p)))))[np.argsort(o)]
    T["pass"] = (T.p_holm < 0.05) & T.above_matched
    T.to_csv(OUT / "score_tests.csv", index=False)
    pd.DataFrame(mrows).to_csv(OUT / "score_matched.csv", index=False)
    print(T.round(4).T.to_string(), flush=True)

    # descriptives
    desc = []
    for name, (r, grp, h) in DESC_TIME.items():
        s, n0 = float(S_time(groups[grp], Y[r][grp], h, tod)[0, 0]), S_time(groups[grp], Y[r][grp], h, tod, P=P)[:, 0]
        c = cosinor(groups[grp], Y[r][grp], tod)[0]
        desc.append(dict(item=name, predicted_peak_h=h, statistic=s, p_one_sided=(1 + np.sum(n0 >= s)) / (1 + n_perm),
                         fitted_peak_h=float(c[3][0]), amplitude_log2=float(c[2][0])))
    u = g[BLOOD]
    yn = {BLOOD: score_of(L, col, u, SETS["T_cell"], resid_on=SETS["neutrophil"])[:, None]}
    h = TIME_TESTS["blood_T"][4]
    s, n0 = float(S_time(g, yn, h, tod)[0, 0]), S_time(g, yn, h, tod, P=P)[:, 0]
    desc.append(dict(item="blood_T_resid_neutrophil", predicted_peak_h=h, statistic=s, p_one_sided=(1 + np.sum(n0 >= s)) / (1 + n_perm)))
    gt = groups["tissues"]
    for label, wt in (("tissues_T_vs_GC", None), ("tissues_T_vs_GC_given_time", tod)):
        z = float(partial_z(gt, Y["T_cell"]["tissues"], Y["GC"]["tissues"], with_time=wt)[0])
        n0 = partial_z(gt, Y["T_cell"]["tissues"], Y["GC"]["tissues"], order=keys[:2000], with_time=wt)
        desc.append(dict(item=label, statistic=z, partial_r=float(np.tanh(z)), null_sd=float(n0.std()),
                         p_two_sided=(1 + np.sum(np.abs(n0) >= abs(z))) / (1 + len(n0))))
    pc = []
    for r in ("T_cell", "GC", "naive", "effector"):
        for grp in groups:
            for ax in M.AXES:
                pc += [dict(readout=r, group=grp, axis=ax, compartment=k, donors=n, amplitude_log2=float(a[0]), peak_h=float(p[0]),
                            t_cos=float(zz[0].real), t_sin=float(zz[0].imag)) for k, n, a, p, zz in cosinor(groups[grp], Y[r][grp], tod, ax)]
    pd.DataFrame(pc).to_csv(OUT / "score_per_compartment.csv", index=False)
    if not quick and not args.placebo:                              # ventilated donors, recorded time
        Lv, _, _, dv, gv, cv, codv = prepare(gtex.VENTILATOR, args, False, False)
        tv, Pv = dv.tod.to_numpy(float), M.permutations(codv, 2000, M.seed_of("A248vent", args.seed))
        for t, (_, _, r, grp, h) in TIME_TESTS.items():
            Yv = {k: score_of(Lv, cv, u, SETS[r])[:, None] for k, u in gv[grp].items()}
            s, n0 = float(S_time(gv[grp], Yv, h, tv)[0, 0]), S_time(gv[grp], Yv, h, tv, P=Pv)[:, 0]
            c = np.mean([x[4][0] for x in cosinor(gv[grp], Yv, tv)])
            desc.append(dict(item=f"ventilated_{t}", predicted_peak_h=h, statistic=s, p_one_sided=(1 + np.sum(n0 >= s)) / 2001,
                             fitted_peak_h=float((np.angle(c) / W) % 24), donors=len(dv), compartments=len(gv[grp])))
    pd.DataFrame(desc).to_csv(OUT / "score_descriptive.csv", index=False)
    print(pd.DataFrame(desc).round(4).to_string(), flush=True)


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--stage", choices=["power", "score"], required=True)
    ap.add_argument("--placebo", action="store_true", help="times and GC order shuffled once (code test; output/placebo/)")
    args = ap.parse_args(argv)
    r = cli.setup(args)
    global OUT
    OUT = OUT / "placebo" if args.placebo else OUT / "quick" if r["quick"] else OUT
    OUT.mkdir(parents=True, exist_ok=True)
    run(args, r["quick"], args.stage)


if __name__ == "__main__":
    main()
