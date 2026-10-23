"""A21: known-time validation of disease contrasts in injured mouse colon (chrono-atlas, ZT 1, 7, 13, 19).

Pre-registered (A21, written before running; the paper's findings and the pre-publication file were known, the C-DC arm
not). Arms C-DC (undamaged distal colon), R-DC (distal, day 17 after DSS), R-PC (proximal colon of the R-DC mice), 8 mice each,
2 mice per ZT; pseudobulks >= 30 cells, a class enters an arm if present in >= 6 of its 8 mice. Anchors: the fixed
stromal-vascular rule (fibroblast, smooth_muscle, pericyte, endothelial, lymphatic_EC present). Truths free of the anchored
estimator: the host replaced by the known time (loadings a_k and own noise psi_k of the clock-plane states on cos / sin of ZT,
Poisson-weighted fit; coupling to time c = |a|^2 / (|a|^2 + 2 psi), host variance 1).
- R21.1 (per arm): anchored coupling (`Dataset.fit`) against the coupling to known time across compartments; accept Spearman >= 0.7
  and mean error within +-0.15 (healthy check on C-DC).
- R21.2 (R-DC minus C-DC): untimed d*, log z2 ratio, log psi ratio (`Dataset.contrast`, one pooled anchored host, arm labels permuted
  within ZT blocks) against the timed contrasts (animal bootstrap within arm); accept Spearman >= 0.5 for each and the same sign in
  >= 80 % of the compartments whose timed interval excludes 0.
- R21.3 (immune families): verdict damping (slope truth < 0 with interval excluding 0 and |delta c| <= 0.15), decoupling (delta c < -0.15
  with interval excluding 0) or indeterminate, from the timed and the untimed estimates; agreement reported.
- R21.4: phase of each non-anchor compartment relative to the anchors per arm (untimed: own anchored host per arm; timed: arg a_k),
  change R-DC minus C-DC; epithelial mean (colonocyte, secretory, stem); accept untimed within 3 h of timed and both intervals excluding 0.
- Power: `Dataset.contrast(mde=True)` for fibroblast, colonocyte, T at this design.

    python compute_colon_injury.py [--subsample N] [--counts cellbender|raw] [--n-boot B] [--jobs J]
    (the atlas has 24 mice: `--subsample` only reduces permutations, bootstraps and the power replicates)
Outputs: paper/results/figure_1/output/species/a21_{r211_couplings,r211_summary,r212_contrast,r212_recovery,r213_immune,r214_phase,mde}[_raw].csv
Legacy: 21_colon_injury/scripts/colon_injury.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

import _shared as S  # noqa: E402
import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import cli, data  # noqa: E402
from common.paths import out_dir  # noqa: E402

FOLDER = "species"
W, HOURS = S.W, S.HOURS
STROMAL_VASCULAR = ("fibroblast", "smooth_muscle", "pericyte", "endothelial", "lymphatic_EC")
EPITHELIAL = ("colonocyte", "secretory", "stem")
IMMUNE = {"lymphocyte": ("T", "B", "plasma", "ILC"), "myeloid": ("macrophage", "DC")}
FAMILY = {"fibroblast": "fibroblast_stromal", "smooth_muscle": "mural_smooth_muscle", "pericyte": "mural_smooth_muscle",
          "endothelial": "endothelial", "lymphatic_EC": "endothelial", "colonocyte": "epithelial", "secretory": "epithelial",
          "stem": "epithelial", "T": "lymphocyte", "B": "lymphocyte", "plasma": "lymphocyte", "ILC": "lymphocyte",
          "macrophage": "myeloid", "DC": "myeloid"}
MDE_COMPS = ("fibroblast", "colonocyte", "T")


def wrap(a):
    return (np.asarray(a) + np.pi) % (2 * np.pi) - np.pi


def cmean(a):
    return np.angle(np.mean(np.exp(1j * np.asarray(a))))


def arm_dataset(pb, donors, arm):
    p = pb[pb.donor.isin(donors.index[donors.arm == arm])]
    comps = sorted(p.compartment.unique())
    anchors = [k for k in comps if k in STROMAL_VASCULAR]
    return cs.Dataset(p, anchors, donors=donors, name=arm), anchors


def timed_loadings(ds, zt):
    """Known-time truth of every compartment of a (possibly resampled) dataset: the host is the known time (local)."""
    d = ds.data
    hours = pd.Series(d.resection).str.split("#").str[0].map(zt).to_numpy(float)
    f = np.column_stack([np.cos(W * hours), np.sin(W * hours)])
    comp, rows = np.asarray(d.compartment), {}
    for k in ds.compartments:
        i = np.flatnonzero(comp == k)
        if len(i) >= 5:
            rows[k] = S.known_host(d.beta[i], d.V[i], f[i], np.asarray(d.group)[i])
    return pd.DataFrame(rows).T.astype({"z2": float, "psi": float, "V": float, "coupling": float})


def timed_coupling(t):
    return t.coupling


def anchored(ds):
    t = ds.fit(n_boot=0, gate=False).table
    return t


def _boot_r211(ds, zt, seed):
    limit_threads(1)
    try:
        b = ds.resample(np.random.default_rng(seed))
        return anchored(b).coupling, timed_coupling(timed_loadings(b, zt))
    except (ValueError, np.linalg.LinAlgError):
        return None


def r211(pb, donors, zt, n_boot, jobs):
    rows, summ = [], []
    for arm in ("C-DC", "R-DC", "R-PC"):
        ds, anchors = arm_dataset(pb, donors, arm)
        est, tru = anchored(ds).coupling, timed_coupling(timed_loadings(ds, zt))
        reps = [r for r in Parallel(n_jobs=jobs)(delayed(_boot_r211)(ds, zt, 10_000 + b) for b in range(n_boot)) if r is not None]
        for src, point, j in (("anchored", est, 0), ("truth", tru, 1)):
            for k in ds.compartments:
                x = np.array([r[j].get(k, np.nan) for r in reps], float)
                x = x[np.isfinite(x)]
                rows.append(dict(arm=arm, compartment=k, family=FAMILY.get(k, "other"), anchor=k in anchors, estimator=src,
                                 coupling=point.get(k, np.nan), lo=np.quantile(x, .025) if len(x) else np.nan,
                                 hi=np.quantile(x, .975) if len(x) else np.nan, n_boot=len(x)))
        w = pd.DataFrame(rows)
        w = w[w.arm == arm].pivot(index="compartment", columns="estimator", values="coupling").dropna()
        summ.append(dict(arm=arm, estimator="anchored", n=len(w), spearman=spearmanr(w.anchored, w.truth).statistic,
                         mean_error=(w.anchored - w.truth).mean(), mean_abs_error=(w.anchored - w.truth).abs().mean()))
    return pd.DataFrame(rows), pd.DataFrame(summ)


def timed_contrast(t_ctrl, t_case, V0=None):
    """Timed truths of R-DC minus C-DC at one host spread V0 (default 1): d*, log z2 ratio, log psi ratio per compartment."""
    k = t_ctrl.index.intersection(t_case.index)
    V0 = pd.Series(1.0, index=k) if V0 is None else V0.reindex(k)
    cstar = lambda t: t.z2[k] * V0 / (t.z2[k] * V0 + 2 * t.psi[k])  # noqa: E731
    with np.errstate(divide="ignore", invalid="ignore"):
        out = pd.DataFrame(dict(d_star=cstar(t_case) - cstar(t_ctrl), log_z2_ratio=np.log(t_case.z2[k] / t_ctrl.z2[k]),
                                log_psi_ratio=np.log(t_case.psi[k] / t_ctrl.psi[k])))
    return out.replace([np.inf, -np.inf], np.nan)


def _boot_timed(pb, donors, zt, comps, V0, seed):
    limit_threads(1)
    rng = np.random.default_rng(seed)
    out = {}
    try:
        for arm in ("C-DC", "R-DC"):
            ds, _ = arm_dataset(pb[pb.compartment.isin(comps)], donors, arm)
            out[arm] = timed_loadings(ds.resample(rng), zt)
        return timed_contrast(out["C-DC"], out["R-DC"], V0)
    except (ValueError, np.linalg.LinAlgError):
        return None


def r212(pb, donors, zt, n_boot, n_perm, jobs):
    p = pb[pb.donor.str.startswith(("C-DC", "R-DC"))]
    both = p.groupby("compartment").donor.apply(lambda s: s.str.split(":").str[0].nunique()) == 2
    comps = list(both.index[both])
    p = p[p.compartment.isin(comps)]
    anchors = [k for k in comps if k in STROMAL_VASCULAR]
    ds = cs.Dataset(p, anchors, donors=donors, name="colon_DSS")
    con = ds.contrast("arm", levels=("C-DC", "R-DC"), blocks="zt", n_perm=n_perm, n_boot=n_boot).table
    t_ctrl = timed_loadings(arm_dataset(p, donors, "C-DC")[0], zt)
    t_case = timed_loadings(arm_dataset(p, donors, "R-DC")[0], zt)
    V0 = timed_loadings(ds, zt).V
    pt = timed_contrast(t_ctrl, t_case, V0)
    reps = [r for r in Parallel(n_jobs=jobs)(delayed(_boot_timed)(p, donors, zt, comps, V0, 20_000 + b) for b in range(n_boot)) if r is not None]
    rows = []
    for k in comps:
        r = con.loc[k]
        row = dict(compartment=k, family=FAMILY.get(k, "other"), anchor=k in anchors, n_control=r["donors_C-DC"],
                   n_case=r["donors_R-DC"], p_perm=r.p_star, status=r.status)
        for nm in ("d_star", "log_z2_ratio", "log_psi_ratio"):
            row[f"{nm}_untimed"] = r[nm]
            row[f"{nm}_timed"] = pt[nm].get(k, np.nan)
            x = np.array([q[nm].get(k, np.nan) for q in reps], float)
            x = x[np.isfinite(x)]
            row[f"{nm}_timed_lo"] = np.quantile(x, .025) if len(x) > 20 else np.nan
            row[f"{nm}_timed_hi"] = np.quantile(x, .975) if len(x) > 20 else np.nan
        row["d_untimed_lo"], row["d_untimed_hi"] = r.lo, r.hi
        rows.append(row)
    tab = pd.DataFrame(rows)
    rec = []
    for nm in ("d_star", "log_z2_ratio", "log_psi_ratio"):
        ok = tab[[f"{nm}_untimed", f"{nm}_timed"]].dropna()
        sig = tab[(tab[f"{nm}_timed_lo"] > 0) | (tab[f"{nm}_timed_hi"] < 0)]
        same = np.sign(sig[f"{nm}_untimed"]) == np.sign(sig[f"{nm}_timed"])
        rec.append(dict(estimand=nm, n=len(ok), spearman=spearmanr(ok.iloc[:, 0], ok.iloc[:, 1]).statistic, n_timed_excl0=len(sig),
                        same_sign=float(same.mean()) if len(sig) else np.nan))
    # R21.3: immune families. Timed: bootstrap draws (same replicate across compartments). Untimed: point estimates, interval =
    # mean of the members' bounds (the contrast does not expose its draws; approximation)
    fam = []
    for fname, members in IMMUNE.items():
        ks = [k for k in members if k in pt.index]
        if not ks:
            continue
        D = np.nanmean(np.stack([np.stack([q.reindex(ks)[nm] for nm in ("d_star", "log_z2_ratio", "log_psi_ratio")], 1)
                                 for q in reps]), axis=1)  # draws x 3
        for src in ("timed", "untimed"):
            pe = tab.set_index("compartment").loc[ks, [f"{nm}_{src}" for nm in ("d_star", "log_z2_ratio", "log_psi_ratio")]].mean(0).to_numpy()
            if src == "timed":
                lo, hi = np.nanquantile(D, .025, 0), np.nanquantile(D, .975, 0)
            else:
                lo = np.array([con.loc[ks, "lo"].mean(), np.nan, np.nan])
                hi = np.array([con.loc[ks, "hi"].mean(), np.nan, np.nan])
            slope_hi = hi[1] if src == "timed" else np.nan
            if src == "timed":
                verdict = ("damping" if slope_hi < 0 and abs(pe[0]) <= 0.15 else "decoupling" if hi[0] < 0 and pe[0] < -0.15 else "indeterminate")
            else:  # untimed rule on d* with the members' mean bounds; slope interval unavailable: point estimate < 0
                verdict = ("damping" if pe[1] < 0 and abs(pe[0]) <= 0.15 else "decoupling" if hi[0] < 0 and pe[0] < -0.15 else "indeterminate")
            fam.append(dict(family=fname, source=src, compartments=",".join(ks), d_star=pe[0], d_star_lo=lo[0], d_star_hi=hi[0],
                            log_z2_ratio=pe[1], log_z2_lo=lo[1], log_z2_hi=hi[1], log_psi_ratio=pe[2], verdict=verdict))
    return tab, pd.DataFrame(rec), pd.DataFrame(fam), ds


def _rel_untimed(ds, anchors):
    t = ds.fit(n_boot=0, gate=False).table
    off = t.offset_h / HOURS
    anc = cmean(off[anchors])
    return {k: off[k] - anc for k in t.index if k not in anchors}


def _rel_timed(ds, anchors, zt):
    t = timed_loadings(ds, zt)
    ang = np.angle(t.z.astype(complex))
    anc = cmean(ang[t.index.isin(anchors)])
    return {k: ang[i] - anc for i, k in enumerate(t.index) if k not in anchors}


def _phase_once(pb, donors, zt, seed):
    limit_threads(1)
    rng = None if seed is None else np.random.default_rng(seed)
    out = {}
    try:
        for arm in ("C-DC", "R-DC"):
            ds, anchors = arm_dataset(pb, donors, arm)
            if rng is not None:
                ds = ds.resample(rng)
            out[arm] = (_rel_untimed(ds, anchors), _rel_timed(ds, anchors, zt))
        return out
    except (ValueError, np.linalg.LinAlgError):
        return None


def r214(pb, donors, zt, n_boot, jobs):
    p = pb[pb.donor.str.startswith(("C-DC", "R-DC"))]
    both = p.groupby("compartment").donor.apply(lambda s: s.str.split(":").str[0].nunique()) == 2
    p = p[p.compartment.isin(both.index[both])]

    def shifts(o):
        s = {}
        for src, j in (("untimed", 0), ("timed", 1)):
            c, r = o["C-DC"][j], o["R-DC"][j]
            ks = sorted(set(c) & set(r))
            s[src] = {k: wrap(r[k] - c[k]) for k in ks}
            ep = [s[src][k] for k in ks if k in EPITHELIAL]
            if ep:
                s[src]["epithelial_mean"] = cmean(ep)
        return s

    pt = shifts(_phase_once(p, donors, zt, None))
    reps = [shifts(o) for o in Parallel(n_jobs=jobs)(delayed(_phase_once)(p, donors, zt, 30_000 + b) for b in range(n_boot)) if o is not None]
    rows = []
    for k in pt["untimed"]:
        row = dict(compartment=k)
        for src in ("untimed", "timed"):
            v = pt[src].get(k, np.nan)
            x = np.array([r[src].get(k, np.nan) for r in reps], float)
            x = v + wrap(x[np.isfinite(x)] - v)
            row[f"shift_h_{src}"] = v * HOURS
            row[f"shift_h_{src}_lo"] = np.quantile(x, .025) * HOURS if len(x) else np.nan
            row[f"shift_h_{src}_hi"] = np.quantile(x, .975) * HOURS if len(x) else np.nan
        row["untimed_minus_timed_h"] = wrap(pt["untimed"].get(k, np.nan) - pt["timed"].get(k, np.nan)) * HOURS
        row["n_boot"] = len(reps)
        rows.append(row)
    return pd.DataFrame(rows)


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--counts", default="cellbender", choices=["cellbender", "raw"])
    ap.add_argument("--n-boot", type=int, default=None, help="animal bootstrap replicates (default 500; 20 for a check)")
    args = ap.parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    tag = "" if args.counts == "cellbender" else "_raw"
    pb, donors = data.load("colon_chrono" if args.counts == "cellbender" else "colon_chrono_raw")
    zt = donors.zt
    n_boot = args.n_boot or (20 if run["quick"] else 500)
    n_perm = 100 if run["quick"] else 1000
    pd.set_option("display.width", 250)
    cpl, summ = r211(pb, donors, zt, min(n_boot, 200), args.jobs)
    cpl.to_csv(out / f"a21_r211_couplings{tag}.csv", index=False)
    summ.to_csv(out / f"a21_r211_summary{tag}.csv", index=False)
    print(summ.round(3).to_string(index=False), flush=True)
    tab, rec, fam, ds = r212(pb, donors, zt, n_boot, n_perm, args.jobs)
    tab.to_csv(out / f"a21_r212_contrast{tag}.csv", index=False)
    rec.to_csv(out / f"a21_r212_recovery{tag}.csv", index=False)
    fam.to_csv(out / f"a21_r213_immune{tag}.csv", index=False)
    print(rec.round(3).to_string(index=False))
    print(fam.round(3).to_string(index=False), flush=True)
    ph = r214(pb, donors, zt, n_boot, args.jobs)
    ph.to_csv(out / f"a21_r214_phase{tag}.csv", index=False)
    print(ph.round(2).to_string(index=False), flush=True)
    if args.counts == "cellbender":
        mde = ds.contrast("arm", levels=("C-DC", "R-DC"), blocks="zt", targets=list(MDE_COMPS), n_perm=50 if run["quick"] else 200,
                          n_boot=0, mde=True, n_mde=5 if run["quick"] else 100)
        mde.mde.to_csv(out / "a21_mde.csv") if mde.mde is not None else None
        print(mde.mde)


if __name__ == "__main__":
    main()
