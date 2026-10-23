"""Where in the internal day were the samples of a cohort taken? (A25)

Pre-registered as A25 (2026-10-14; docs/PLAN_internal_desync.md (archive), run note appended: V25.2 passes, V25.3 and V25.4 fail
as specified, V25.1 fails at 6 h windows, T1 and T2 pass, T3 not supported). The estimand is the time read by the host's
anchor clocks in the GTEx clock frame (clock time plus the cohort's mean clock offset, not identified without recorded
times). Estimator (`_internal_time.py`, local): the host covariance gives the window centre modulo 12 h (major axis),
its width w (ratio r of radial to tangential variance) and amplitude; the 12 h flip is resolved by curvature (`c_curv`)
or by the hospital-hours prior (`c_day`, an assumption); donor bootstrap (host refitted); the axis is interpretable when
the upper bound of r is below 0.6.
Parts (`--parts`):
  sim      V25.1 simulation (`clocksync.simulate`): 4 compartments, n 20/40/80, windows 6/9/12/18/24 h at centres 10 and 14 h,
           coupling 0.5/0.8; acceptance at n = 40, coupling 0.8, w <= 12 h: median axis error <= 1.5 h and flip correct >= 80 %.
  gtex     V25.2 GTEx fast deaths, tissues with >= 30 donors; windows centred 0, 3, ..., 21 h of width 6/9/12/18 h plus all
           donors (>= 20 donors, tissues with >= 8): median |axis error| <= 2 h and flip correct >= 75 % over r < 0.6.
  colon    V25.4 mouse colon chrono-atlas (arm = stratum; host on fibroblast, smooth muscle, endothelium): 2, 3 or 4 consecutive
           ZTs; every subset with r < 0.6 within 2 h of the common offset.
  cohorts  R25: every cohort with an anchored host (living: procedure in working hours; deceased: any hour; other: reported
           only), T1 (axis centres of interpretable living cohorts concentrated towards 13 h mod 12, V test, mean within +-3 h),
           T2 (more curvature centres in 07-19 h than chance, binomial), T3 (r larger in deceased than living, Mann-Whitney).
V25.3 (Chen 2016 prefrontal cortex, microarray) is not ported: `Dataset` takes counts, not array states.

    python compute_internal_time.py [--subsample N] [--seed S] [--jobs 4] [--parts ...] [--reps N]
Outputs (paper/results/figure_1/output/windows/): internal_time_{sim,sim_summary,gtex,colon,cohorts,tests}.csv
GTEx rows are windows (aggregates), never donors.
Legacy: legacy/desync/common/a25_internal_time.py -> common_output/a25_{sim,sim_summary,gtex,colon,cohorts,tests}.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402
from scipy import stats as sps  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import cli, data  # noqa: E402
from common.paths import out_dir  # noqa: E402

import _internal_time as T  # noqa: E402
import _shared as sh  # noqa: E402

FOLDER = "windows"
PARTS = ("sim", "gtex", "colon", "cohorts")
LIVING = ["kidney_sc", "ckd_sc", "endometriosis", "skin_ad", "skin_ssc", "kong_ileum", "blood_onek1k", "blood_aida",
          "lupus", "lung_hlca|hospital", "dcm|case"]
DECEASED = ["psychad_mssm", "psychad_hbcc", "seaad_dlpfc", "pd", "AD_radc", "retina", "ia_spleen", "tabula_sapiens",
            "heart", "lung_natri_control", "lung_hlca|anytime"]
OTHER = ["islets", "organ_breast", "organ_adipose", "organ_muscle", "organ_bone_marrow", "organ_colon", "ia_crohn"]
HOSPITAL = {"surgical_resection", "biopsy", "brush", "biopsy;brush"}
ANYTIME = {"donor_lung", "donor_lung_cryopreserved", "lung_explant", "autopsy"}
COLON_ANCHORS = ["fibroblast", "smooth_muscle", "endothelial"]  # lymphatic EC decoupled in R-DC (A21)


# ---- V25.1 ----------------------------------------------------------------------------------------------------

def sim_one(n, w_h, c_h, cp, rep, K=4):
    limit_threads(1)
    rng = np.random.default_rng(10_000 * int(w_h) + 100 * n + rep)
    z = np.exp(1j * T.W * rng.normal(0, 1.0, K)) * np.exp(rng.normal(0, 0.2, K))
    z /= z[0] / abs(z[0])
    psi = 0.5 * (1 - cp) / cp  # coupling at full-day spread, lam = 1
    wf = cs.gtex()
    comps = {f"k{j}": cs.Compartment(z=complex(z[j]), psi=psi) for j in range(K)}
    tab, t = cs.simulate(wf, comps, n=n, window_h=(c_h - w_h / 2, c_h + w_h / 2), jitter_h=0.0, seed=int(rng.integers(1 << 30)))
    shift = np.angle((z / abs(z)).mean())
    try:
        ds = cs.Dataset(tab, anchors=list(comps), name="sim")
        e, ft = T.host_estimate(ds)
    except Exception as ex:
        return dict(n=n, w_set=w_h, c_set=c_h, rep=rep, coupling=cp, status=f"error: {type(ex).__name__}")
    tw = T.true_window(t + shift * T.HOURS)
    truth = pd.Series(t, index=sorted(set(tab.resection))).reindex(e["donor_h"].index)
    return dict(n=n, w_set=w_h, c_set=c_h, rep=rep, coupling=cp, status="ok",
                **{k: v for k, v in e.items() if k not in ("donor_h", "tangent", "donor_day_h", "tangent_day", "sigma_f")}, **tw,
                donor_err=float(np.median(np.abs(T.cdiff(e["donor_h"].to_numpy(), T.wrap24(truth.to_numpy() + shift * T.HOURS))))))


def sim(args, run, out):
    reps = args.reps if args.reps is not None else (2 if run["quick"] else 50)
    ns, ws = ((40,), (9, 18)) if run["quick"] else ((20, 40, 80), (6, 9, 12, 18, 24))
    grid = [(n, w, c, cp, r) for n in ns for w in ws for c in (10.0, 14.0) for cp in (0.5, 0.8) for r in range(reps)]
    d = pd.DataFrame(Parallel(n_jobs=args.jobs)(delayed(sim_one)(*g) for g in grid))
    d = d[d.status == "ok"].copy()
    d["err_axis12"] = T.cdiff(2 * d.c_axis12, 2 * np.mod(d.c_true, 12.0)) / 2
    d["err_curv"] = T.cdiff(d.c_curv, d.c_true)
    d["flip_ok"] = np.abs(d.err_curv) < 6
    d.to_csv(out / "internal_time_sim.csv", index=False)
    s = d.groupby(["n", "w_set", "coupling"]).agg(
        axis_mae=("err_axis12", lambda x: np.median(np.abs(x))), flip_ok=("flip_ok", "mean"), w_est=("w", "median"),
        w_true=("w_true", "median"), r=("r", "median"), donor_mae=("donor_err", "median")).reset_index()
    s.to_csv(out / "internal_time_sim_summary.csv", index=False)
    acc = s[(s.n == 40) & (s.coupling == 0.8) & (s.w_set <= 12)]
    print(s.round(2).to_string(index=False))
    print("V25.1 acceptance (n = 40, coupling 0.8, w <= 12 h): axis error <= 1.5 h and flip >= 80 %:",
          bool(len(acc) and (acc.axis_mae <= 1.5).all() and (acc.flip_ok >= 0.8).all()))


# ---- V25.2 / V25.4 ------------------------------------------------------------------------------------------

def window_row(ds_all, keep, t_true, row, anchors, min_donors=20, min_comp_donors=8):
    """Internal-time window of the donors `keep` against their true times `t_true` (h)."""
    dd = ds_all.subset(donors=keep)
    cnt = dd.table.groupby("compartment").resection.nunique()
    anc = [a for a in anchors if cnt.get(a, 0) >= min_comp_donors]
    row = dict(row, n=len(dd.donor_ids), n_anchors=len(anc))
    if row["n"] < min_donors or len(anc) < 2:
        return row | dict(status="too few")
    try:
        dd = ds_all.subset(donors=keep, compartments=anc, anchors=anc)
        e, ft = T.host_estimate(dd, anc)
    except Exception as ex:
        return row | dict(status=f"error: {type(ex).__name__}: {ex}")
    t = pd.Series(t_true).reindex(e["donor_h"].index)
    tw = T.true_window(t.to_numpy())
    est = e["donor_h"].to_numpy()
    sa, sb = np.sin(T.W * est - np.angle(np.exp(1j * T.W * est).mean())), np.sin(T.W * t.to_numpy() - np.angle(np.exp(1j * T.W * t.to_numpy()).mean()))
    return row | {k: v for k, v in e.items() if k not in ("donor_h", "tangent", "donor_day_h", "tangent_day", "sigma_f")} | tw | dict(
        status="ok", err_axis12=T.cdiff(2 * e["c_axis12"], 2 * np.mod(tw["c_true"], 12.0)) / 2,
        err_curv=T.cdiff(e["c_curv"], tw["c_true"]), err_day=T.cdiff(e["c_day"], tw["c_true"]),
        donor_mae=float(np.median(np.abs(T.cdiff(est, t.to_numpy())))),
        donor_ccorr=float((sa * sb).sum() / np.sqrt((sa ** 2).sum() * (sb ** 2).sum())),
        flip_ok=bool(abs(T.cdiff(e["c_curv"], tw["c_true"])) < 6))


def gtex_windows(args, run, out):
    pb, donors = sh.gtex_fast(30, args.subsample, args.seed)
    ds = cs.Dataset(pb, anchors=None, donors=donors, name="gtex_fast")
    tod = donors.tod.dropna() % 24
    anchors = list(ds.compartments)
    sets = [("all", np.nan, np.nan, tod.index.to_numpy())]
    centres, widths = (range(0, 24, 6), (9, 18)) if run["quick"] else (range(0, 24, 3), (6, 9, 12, 18))
    for c in centres:
        for w in widths:
            sets.append((f"{c:g}h+-{w / 2:g}", c, w, tod.index[np.abs(T.cdiff(tod.to_numpy(), c)) <= w / 2].to_numpy()))
    rows = Parallel(n_jobs=args.jobs)(delayed(window_row)(ds, keep, tod, dict(data="gtex_fast", window=nm, c_set=c, w_set=w), anchors)
                                      for nm, c, w, keep in sets)
    d = pd.DataFrame(rows)
    d.to_csv(out / "internal_time_gtex.csv", index=False)
    ok = d[(d.status == "ok")]
    good = ok[(ok.r < T.R_MAX) & ok.window.ne("all")] if len(ok) else ok
    print(ok[["window", "n", "n_anchors", "c_true", "c_axis12", "c_curv", "err_axis12", "err_curv", "r", "r_true", "w",
              "w_true"]].round(2).to_string(index=False))
    if len(good):
        print(f"V25.2: median |axis error| {good.err_axis12.abs().median():.2f} h (<= 2), flip correct "
              f"{good.flip_ok.mean():.2f} (>= 0.75) over {len(good)} windows with r < 0.6")


def colon_windows(args, run, out):
    pb, donors = data.load("colon_chrono", subsample=args.subsample, seed=args.seed)
    ds = cs.Dataset(pb, anchors=COLON_ANCHORS, donors=donors, strata="arm", name="colon_chrono")
    zt = donors.zt.astype(float)
    rows, zts = [], (1.0, 7.0, 13.0, 19.0)
    for m in (2, 3, 4):
        for i in range(4 if m < 4 else 1):
            pick = [zts[(i + j) % 4] for j in range(m)]
            keep = zt.index[zt.isin(pick)].to_numpy()
            r = window_row(ds, keep, zt, dict(zts="/".join(f"{p:g}" for p in pick)), COLON_ANCHORS, min_donors=6, min_comp_donors=4)
            rows.append(r)
    d = pd.DataFrame(rows)
    ok = d.status == "ok"
    if ok.any():
        off = np.angle(np.exp(1j * T.W * d.loc[ok, "err_curv"]).mean()) * T.HOURS
        d["offset_common"] = off
        d["err_relative"] = T.cdiff(d.err_curv, off)
    d.to_csv(out / "internal_time_colon.csv", index=False)
    print(d.drop(columns=[c for c in ("status",) if False]).round(2).to_string(index=False))


# ---- R25 ------------------------------------------------------------------------------------------------------

def cohort(entry: str, subsample, seed):
    """(pseudobulk, donors, anchors, rule, name) of one entry (`name|part`: case / control / hospital / anytime)."""
    name, _, part = entry.partition("|")
    if name.startswith("cancer_"):
        pb, donors = data.load("cancer_sc", study=name.removeprefix("cancer_"))
        lab = None
    else:
        row = sh.targets().loc[name]
        pb, donors = data.load(row.dataset)
        lab = donors[name] if (row.kind == "contrast" and name in donors) else None
        if lab is not None:
            lab = lab[lab.isin(["control", "case"])]
            pb = pb[pb.donor.isin(lab.index)]
        elif name == "lung_hlca":
            pb = pb[pb.donor.isin(donors.index[donors.healthy.astype(bool)])]
    if part in ("case", "control"):
        pb = pb[pb.donor.isin(lab.index[lab == part])]
    elif part:
        m = donors.tissue_sampling_method_major.astype(str).reindex(pb.donor.unique())
        pb = pb[pb.donor.isin(m.index[m.isin(HOSPITAL if part == "hospital" else ANYTIME)])]
    if subsample is not None:
        keep = np.random.default_rng(seed).permutation(sorted(pb.donor.unique()))[:subsample]
        pb = pb[pb.donor.isin(keep)]
    pb = pb.reset_index(drop=True)
    for rule in ("learned", "fixed"):
        try:
            anc = sh.anchor_set(name, rule)
        except KeyError:
            anc = None
        anc = [a for a in (anc or []) if a in set(pb.compartment)]
        if len(anc) >= 2:
            return pb, donors, anc, rule
    return pb, donors, None, None


def cohort_row(entry, klass, run, subsample, seed):
    limit_threads(1)
    try:
        pb, donors, anc, rule = cohort(entry, subsample, seed)
        if anc is None:
            return dict(entry=entry, klass=klass, status="no anchors")
        ds = cs.Dataset(pb, anchors=anc, name=entry)
        e, ft = T.host_estimate(ds, anc)
        reps = T.boot(ds, anc, e, run["n_boot"], seed)
        return dict(entry=entry, klass=klass, rule=rule, n=len(ds.donor_ids), n_anchors=len(anc), anchors=";".join(anc),
                    status="ok", **T.summarise(e, reps))
    except Exception as ex:
        return dict(entry=entry, klass=klass, status=f"error: {type(ex).__name__}: {ex}")


def tests(out_df: pd.DataFrame) -> pd.DataFrame:
    """T1-T3 of A25."""
    ok = out_df[out_df.status == "ok"]
    liv, dec = ok[ok.klass == "living"], ok[ok.klass == "deceased"]
    li = liv[liv.interpretable]
    res = {}
    a = 2 * np.pi * li.c_axis12.to_numpy() / 12
    n = len(a)
    mu0 = 2 * np.pi * 1.0 / 12  # 13 h mod 12
    V = np.sum(np.cos(a - mu0)) if n else np.nan
    u = V * np.sqrt(2 / n) if n else np.nan
    mean12 = np.mod(np.angle(np.exp(1j * a).mean()) * 12 / (2 * np.pi), 12) if n else np.nan
    p1 = 1 - sps.norm.cdf(u) if n else np.nan
    res["T1"] = dict(n=n, mean_axis12=mean12, R=abs(np.exp(1j * a).mean()) if n else np.nan, u=u, p=p1,
                     accept=bool(n and p1 < 0.05 and abs(T.cdiff(2 * mean12, 2.0) / 2) <= 3))
    inday = int(((li.c_curv >= 7) & (li.c_curv < 19)).sum())
    res["T2"] = dict(n=n, in_day=inday, p=sps.binomtest(inday, n, 0.5, alternative="greater").pvalue if n else np.nan)
    res["T3"] = dict(n_living=len(liv), n_deceased=len(dec), r_living=liv.r.median(), r_deceased=dec.r.median(),
                     interpretable_living=liv.interpretable.mean(), interpretable_deceased=dec.interpretable.mean(),
                     p=sps.mannwhitneyu(dec.r, liv.r, alternative="greater").pvalue if len(liv) and len(dec) else np.nan)
    return pd.DataFrame(res).T


def cohorts(args, run, out):
    cancer = sorted(t for t in pd.read_csv(sh.AN.TABLE).target.unique() if str(t).startswith("cancer_"))
    todo = [(e, "living") for e in LIVING + cancer] + [(e, "deceased") for e in DECEASED] + [(e, "other") for e in OTHER]
    if run["quick"]:
        todo = [t for t in todo if t[0] in ("skin_ad", "kidney_sc", "heart", "lung_natri_control", "islets", "ckd_sc",
                                            "psychad_mssm", "lung_hlca|hospital")]
    rows = Parallel(n_jobs=args.jobs)(delayed(cohort_row)(e, k, run, args.subsample, args.seed) for e, k in todo)
    d = pd.DataFrame(rows)
    d.to_csv(out / "internal_time_cohorts.csv", index=False)
    if (d.status == "ok").any():
        t = tests(d)
        t.to_csv(out / "internal_time_tests.csv")
        print(t.to_string())


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--parts", nargs="+", default=list(PARTS), choices=PARTS)
    ap.add_argument("--reps", type=int, default=None, help="V25.1 replicates per cell (default 50; 2 with --subsample)")
    args = ap.parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    pd.set_option("display.width", 220)
    for p in args.parts:
        {"sim": sim, "gtex": gtex_windows, "colon": colon_windows, "cohorts": cohorts}[p](args, run, out)


if __name__ == "__main__":
    main()
