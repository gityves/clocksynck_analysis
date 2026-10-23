"""Aim 1 (GTEx, timed): the shared departure of the organs from recorded time. Rank, gains, brain vs periphery, covariates.

Pre-registered in docs/PLAN_internal_desync.md (archive) (Aim 1; the pre-A1 analysis of `desync/01_gtex_organs`, message 2 of the
paper). Every tissue of a donor is fitted on the recorded time of death in the clock plane,
beta_rk = mu_k + a_k e^{i w t_r} + e_rk; the residual in the tissue's rotating frame, u = e / (a e^{iwt}), splits into a
phase departure (Im u / w, h; > 0: ahead of recorded time) and a log-amplitude departure (Re u). Per cohort (fast deaths,
Hardy 1-2; ventilator cases, Hardy 0):
  1. rank of the shared departure: factor analysis with missing data of the phase, amplitude and joint departures, number
     of factors q = 0..4 by 10-fold cross-validated held-out log-density against a null with donors permuted within
     tissue, for all tissues and for the tissues outside the brain; two-factor varimax loadings;
  2. gains on the shared departure (one-factor loadings scaled to a mean of 1), share of each tissue's departure that is
     shared, absolute timed amplitude, ischaemic-time slope of the phase departure (h per h);
  3. the untimed host-referenced fit of the same tissues (`Dataset.fit`: coupling, relative amplitude, offset);
  4. brain vs periphery: correlation over donors of the mean standardised departure of the brain regions and of the other
     tissues, disattenuated by split-half reliability (r*, message 2), donor-bootstrap interval of r;
  5. covariates (age, sex, ischaemic time, Hardy 2): complex least squares beta = mu + d x + (a + c x) e^{it} per tissue,
     log2 amplitude change and phase advance per unit, pooled over tissues (fixed effect, Q, I2).
Intervals: donor bootstrap (host refitted for the untimed columns). Aggregates only (GTEx time of death is controlled
access). Seeds fixed.

    python compute_organs_timed.py [--subsample N] [--jobs 1] [--cohorts fast ventilator]
Outputs: paper/results/figure_2/output/organs/{rank_cv,tissues,covariates,covariates_pooled,phase_loadings_q2,brain_periphery,
phase_corr_<cohort>,phase_resid_corr_<cohort>}.csv
Legacy: legacy/desync/01_gtex_organs/scripts/gtex_organs.py (outputs legacy/desync/01_gtex_organs/output/).
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from common import cli  # noqa: E402
from common.paths import out_dir  # noqa: E402

import _shared as S  # noqa: E402

FOLDER = "organs"
W = S.W


def attach(ds: cs.Dataset, pb, donors, isch) -> None:
    """Row covariates onto `ds.table` (rows of `ds.table` follow `pb`), so that bootstrap replicates carry them."""
    d = donors.reindex(pb.donor)
    ds.table = ds.table.assign(tod=d.tod.to_numpy(), age=d.age.to_numpy(), male=d.male.to_numpy(),
                               hardy2=(d.hardy.to_numpy() == 2).astype(float), isch=isch)


def untimed(ds: cs.Dataset) -> pd.DataFrame:
    """Coupling, log2 relative amplitude and offset of every tissue on the host of the others (`Dataset.fit`; the
    legacy took the coupling from the all-tissue host, here an anchor is read on its leave-one-out host)."""
    t = ds.fit(n_boot=0, gate=False).table
    return pd.DataFrame({"coupling": t.coupling, "untimed_log2_amp": t.log_amp_rel / np.log(2),
                         "untimed_offset_h": t.offset_h})


def tissue_summary(ds: cs.Dataset) -> pd.DataFrame:
    beta, V, tb = ds.data.beta, ds.data.V, ds.table
    donor, tissue = ds.data.resection, ds.data.compartment
    t_h, isch = tb.tod.to_numpy(), tb.isch.to_numpy()
    lab, _, a, _, e = S.timed_fit(beta, V, t_h, tissue)
    ph, am = S.departures(e, a[np.searchsorted(lab, tissue)], t_h)
    Wp, Wa = S.wide(ph, donor, tissue), S.wide(am, donor, tissue)
    out = pd.DataFrame({"donors": Wp.notna().sum(), "abs_amp": np.abs(a), "clock_offset_h": np.angle(a) / W},
                       index=pd.Index(lab, name="tissue"))
    for name, Wd in (("phase", Wp), ("amp", Wa)):
        Y = (Wd - Wd.mean()).to_numpy()
        B, psi_d, _, _ = S.fa(Y, 1)
        sign = np.sign(B[:, 0].sum())
        out[f"{name}_dep_sd"] = np.sqrt(np.nanvar(Y, 0))
        out[f"{name}_gain"] = sign * B[:, 0] / np.mean(sign * B[:, 0])
        out[f"{name}_shared_frac"] = B[:, 0] ** 2 / (B[:, 0] ** 2 + psi_d)
    d = pd.DataFrame({"donor": donor, "tissue": tissue, "ph": ph, "isch": isch}).dropna()
    d[["ph_w", "isch_w"]] = d[["ph", "isch"]] - d.groupby("donor")[["ph", "isch"]].transform("mean")
    out["phase_per_isch_h"] = d.groupby("tissue").apply(
        lambda g: np.sum(g.ph_w * g.isch_w) / np.sum(g.isch_w ** 2), include_groups=False)
    out = out.join(untimed(ds))
    out["timed_log2_amp"] = np.log2(out.abs_amp) - np.log2(out.abs_amp).mean()
    return out


def covariate_effects(beta, t_h, tissue, X: pd.DataFrame) -> pd.DataFrame:
    """Per tissue, beta = mu + d x + (a + c x) e^{it} by complex least squares (x centred within tissue): log2
    amplitude change Re(c/a)/ln 2 and phase advance Im(c/a)/w (h) per unit of each column of X."""
    y = beta[:, 0] + 1j * beta[:, 1]
    rows = []
    for k in np.unique(tissue):
        i = tissue == k
        x = X[i].to_numpy(float)
        x = x - x.mean(0)
        ok = x.std(0) > 0
        x = x[:, ok]
        e = np.exp(1j * t_h[i] * W)
        D = np.column_stack([np.ones(i.sum()), x, e, x * e[:, None]])
        coef = np.linalg.lstsq(D, y[i], rcond=None)[0]
        a, c = coef[1 + x.shape[1]], coef[2 + x.shape[1]:]
        for name, cj in zip(X.columns[ok], c):
            rows.append(dict(tissue=k, covariate=name, log2_amp=(cj / a).real / np.log(2),
                             phase_advance_h=(cj / a).imag / W))
    return pd.DataFrame(rows)


def cov_table(ds: cs.Dataset, cols) -> pd.DataFrame:
    tb = ds.table
    return covariate_effects(ds.data.beta, tb.tod.to_numpy(), ds.data.compartment, tb[cols].reset_index(drop=True))


def boot_ci(point: pd.DataFrame, reps: pd.DataFrame, cols) -> pd.DataFrame:
    q = reps.groupby("tissue")[cols].quantile([0.025, 0.975]).unstack()
    q.columns = [f"{c}_{'lo' if p == 0.025 else 'hi'}" for c, p in q.columns]
    return point.join(q)


def pooled(ce: pd.DataFrame, name: str) -> list:
    """Fixed-effect pooling over tissues (`cs.stats.meta`) with Cochran's Q and I2 of each covariate effect."""
    rows = []
    for cov, g in ce.groupby("covariate"):
        for est in ("log2_amp", "phase_advance_h"):
            m = cs.stats.meta(g[est], g[f"{est}_se"])
            w = 1 / g[f"{est}_se"] ** 2
            Q = float(np.sum(w * (g[est] - m["fixed"]) ** 2))
            rows.append(dict(cohort=name, covariate=cov, estimand=est, pooled=m["fixed"], se=m["se_fixed"], Q=Q,
                             df=len(g) - 1, I2=max(0.0, (Q - len(g) + 1) / Q) if Q > 0 else 0.0, n_tissues=len(g)))
    return rows


def run_cohort(name: str, args, run, out) -> dict:
    pb, donors = S.cohort(name, genes="clock", subsample=args.subsample, seed=args.seed)
    isch = S.sample_ischaemic_h(pb)
    ds = S.dataset(pb)
    attach(ds, pb, donors, isch)
    tb = ds.table
    beta, V, donor, tissue = ds.data.beta, ds.data.V, ds.data.resection, ds.data.compartment
    t_h = tb.tod.to_numpy()
    print(f"{name}: {pb.donor.nunique()} donors, {pb.compartment.nunique()} tissues, {len(pb)} samples", flush=True)
    qs = (0, 1, 2) if run["quick"] else (0, 1, 2, 3, 4)
    folds = 5 if run["quick"] else 10

    brain = pd.Series(tissue).str.startswith("Brain").to_numpy()
    cv = []
    for subset, m in (("all", np.ones(len(tb), bool)), ("periphery", ~brain)):
        for p in (False, True):
            cv.append(S.rank_cv(beta[m], V[m], t_h[m], donor[m], tissue[m], qs=qs, n_folds=folds, permute=p,
                                seed=args.seed).assign(subset=subset))
    cv = pd.concat(cv).assign(cohort=name)

    covs = ["age", "male", "isch"] + (["hardy2"] if name == "fast" else [])
    ts = tissue_summary(ds)
    ce = cov_table(ds, covs)

    rng = np.random.default_rng(args.seed)
    rep_ts, rep_ce = [], []
    for b in range(run["n_boot"]):
        rb = ds.resample(rng)
        try:
            rep_ts.append(tissue_summary(rb).reset_index().assign(rep=b))
            rep_ce.append(cov_table(rb, covs).assign(rep=b))
        except (np.linalg.LinAlgError, ValueError):
            continue
    num = [c for c in ts.columns if c != "donors"]
    ts = boot_ci(ts, pd.concat(rep_ts), num).assign(cohort=name) if rep_ts else ts.assign(cohort=name)
    ce = ce.rename(columns={"covariate": "covariate"})
    if rep_ce:
        rce = pd.concat(rep_ce)
        q = rce.groupby(["tissue", "covariate"])[["log2_amp", "phase_advance_h"]].agg(
            [lambda s: s.quantile(0.025), lambda s: s.quantile(0.975), "std"])
        q.columns = [f"{c}_{s}" for c in ("log2_amp", "phase_advance_h") for s in ("lo", "hi", "se")]
        ce = ce.join(q, on=["tissue", "covariate"])
    ce = ce.assign(cohort=name).replace({"covariate": {"isch": "ischaemic_h"}})
    meta = pooled(ce, name) if rep_ce else []

    # shared factor of the phase departure, its residual correlation, varimax loadings, brain vs periphery
    lab, _, a, _, e = S.timed_fit(beta, V, t_h, tissue)
    ph, am = S.departures(e, a[np.searchsorted(lab, tissue)], t_h)
    Wp = S.wide(ph, donor, tissue)
    Y = (Wp - Wp.mean()).to_numpy()
    B, psi, m, _ = S.fa(Y, 1)
    resid = pd.DataFrame(Y - m @ B.T, Wp.index, Wp.columns)
    Z = ((Wp - Wp.mean()) / Wp.std()).to_numpy()
    B2, psi2, _, _ = S.fa(Z, 2)
    L = S.varimax(B2)
    L *= np.sign(L.sum(0))
    loadings = pd.DataFrame(L, Wp.columns, ["factor_a", "factor_b"]).assign(psi=psi2, cohort=name)
    rows = []
    n_split = 50
    for dname, val in (("phase", ph), ("amplitude", am)):
        Zw = S.zwide(val, donor, tissue)
        if Zw.columns.str.startswith("Brain").sum() < 4:
            break
        st = S.block_stat(Zw, rng, n_split)
        lo, hi = S.block_boot(Zw, rng, run["n_boot"])
        isb = Zw.columns.str.startswith("Brain")
        both = pd.concat([S._score(Zw.loc[:, ~isb]), S._score(Zw.loc[:, isb])], axis=1).dropna().to_numpy()
        bs = [np.corrcoef(both[i].T)[0, 1] for i in (rng.integers(0, len(both), len(both))
                                                   for _ in range(run["n_boot"]))]
        rows.append(dict(departure=dname, donors=st["donors"], r=st["r"], r_lo=np.quantile(bs, .025),
                         r_hi=np.quantile(bs, .975), reliability_periphery=st["rel_periph"],
                         reliability_brain=st["rel_brain"], r_disattenuated=st["r_star"], r_star_lo=lo,
                         r_star_hi=hi, cohort=name))
    Wp.corr(min_periods=20).to_csv(out / f"phase_corr_{name}.csv")
    resid.corr(min_periods=20).to_csv(out / f"phase_resid_corr_{name}.csv")
    return dict(rank_cv=cv, tissues=ts, covariates=ce, covariates_pooled=pd.DataFrame(meta),
                phase_loadings_q2=loadings, brain_periphery=pd.DataFrame(rows))


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--cohorts", nargs="+", default=list(S.COHORTS))
    ap.add_argument("--n-boot", type=int, default=None, help="donor bootstrap replicates (default 200; 20 with --subsample)")
    args = ap.parse_args(argv)
    run = cli.setup(args)
    if args.n_boot is not None:
        run["n_boot"] = args.n_boot
    out = out_dir(FOLDER, run["quick"])
    outs = [run_cohort(c, args, run, out) for c in args.cohorts]
    for name in ("rank_cv", "tissues", "covariates", "covariates_pooled", "phase_loadings_q2", "brain_periphery"):
        t = pd.concat([o[name] for o in outs])
        t.to_csv(out / f"{name}.csv", index=name in ("tissues", "phase_loadings_q2"))
    print(pd.concat([o["brain_periphery"] for o in outs]).round(3).to_string())


if __name__ == "__main__":
    main()
