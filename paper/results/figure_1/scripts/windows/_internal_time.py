"""local: not in clocksync. Internal-time window of a cohort from the spread of its host states (A25 estimator,
legacy `desync/common/a25_internal_time.py`), written on the public `Dataset.fit` output.

Over a collection window of width w centred on internal time c, the host covariance Sigma_f (noise removed) has
tangential variance lam^2 (1/2 - sin w / 2w) and radial variance lam^2 (1/2 + sin w / 2w - R^2), R = sin(w/2)/(w/2).
The major axis gives c modulo 12 h; r = radial / tangential gives the uniform-equivalent width w; the 12 h flip is
resolved by the curvature kappa = cov(rho, t^2) of the donors' host means (kappa < 0 at the true centre) or, as the
hospital-hours assumption, by the candidate nearer 13:00 (`c_day`).

`Dataset.fit` returns the donors' posterior host means (`Fit.host`) but not the host covariance Sigma_f; `host_cov`
rebuilds it from the law of total variance, Sigma_f = Cov(posterior means) + E[posterior covariance], with the
posterior covariance of donor s from the loadings of the anchors it has, (Sigma_f^-1 + sum_k z2_k / psi_k I)^-1
(API gap: expose the host covariance).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.optimize import brentq

import clocksync as cs

W = 2 * np.pi / 24
HOURS = 1 / W
R_MAX = 0.6  # axis interpretable when the bootstrap 95 % upper bound of r is below this (pre-registered, A25)


def wrap24(h):
    return np.mod(np.asarray(h, float), 24.0)


def cdiff(a, b):
    """Circular difference a - b in hours, in (-12, 12]."""
    return -(np.mod(-(np.asarray(a, float) - np.asarray(b, float)) + 12.0, 24.0) - 12.0)


def window_r(w):
    """Radial / tangential variance ratio of phases uniform on a window of width w (radians)."""
    w = np.asarray(w, float)
    s = np.sin(w) / (2 * w)
    R = np.sin(w / 2) / (w / 2)
    return (0.5 + s - R ** 2) / (0.5 - s)


def width_of(r: float) -> float:
    """Uniform-equivalent window width (h) for a variance ratio r; 24 h at r >= 1."""
    if not np.isfinite(r) or r <= 0:
        return 0.0
    if r >= 1:
        return 24.0
    lo, hi = 0.05, 2 * np.pi - 1e-6
    if r <= window_r(lo):
        return lo * HOURS
    if r >= window_r(hi):
        return 24.0
    return brentq(lambda w: window_r(w) - r, lo, hi) * HOURS


def true_window(hours) -> dict:
    """Circular mean, resultant length, uniform-equivalent width and r of recorded times."""
    e = np.exp(1j * W * np.asarray(hours, float))
    c = np.angle(e.mean()) * HOURS
    x = np.conj(e) * np.exp(1j * c * W)
    tang, rad = np.var(x.imag), np.var(x.real)
    r = rad / tang if tang > 0 else np.nan
    return dict(c_true=wrap24(c), R_true=abs(e.mean()), r_true=r, w_true=width_of(r))


def centre_by(f, strata):
    f = np.array(f, float, copy=True)
    if strata is None:
        return f - f.mean(0)
    strata = np.asarray(strata)
    for q in np.unique(strata):
        i = strata == q
        f[i] -= f[i].mean(0)
    return f


def host_cov(f, strata, precision) -> np.ndarray:
    """Sigma_f (2 x 2): covariance of the posterior host means `f` (centred within strata) plus the mean posterior
    covariance, `precision` = per-donor information sum_k z2_k / psi_k (isotropic), by fixed-point iteration."""
    fc = centre_by(f, strata)
    C = fc.T @ fc / (len(fc) - 1)
    S = C.copy()
    for _ in range(200):
        inv = np.linalg.inv(S + 1e-12 * np.eye(2))
        post = np.mean([np.linalg.inv(inv + p * np.eye(2)) for p in precision], axis=0)
        new = C + post
        if np.abs(new - S).max() < 1e-10:
            break
        S = new
    return (S + S.T) / 2


def estimate(f, S, shift: float, strata=None, donors=None) -> dict:
    """Window of the donors with host means `f` (n x 2) and covariance `S`; `shift` (rad): reference frame ->
    mean-anchor frame (circular mean of arg z over the host compartments)."""
    ev, evec = np.linalg.eigh((S + S.T) / 2)
    l2, l1 = max(ev[0], 0.0), ev[1]
    axis = np.arctan2(evec[1, 1], evec[0, 1])  # major axis (tangent), mod pi
    r = l2 / l1 if l1 > 0 else np.nan
    w_h = width_of(r)
    w = w_h * W
    lam = np.sqrt(l1 / (0.5 - np.sin(w) / (2 * w))) if w > 0 else np.nan
    fc = centre_by(f, strata)
    cand = axis + np.pi / 2
    e_c, e_t = np.array([np.cos(cand), np.sin(cand)]), np.array([-np.sin(cand), np.cos(cand)])
    rho, t = fc @ e_c, fc @ e_t
    kappa = float(np.cov(rho, t ** 2)[0, 1]) if len(fc) > 2 else np.nan
    c_curv = cand if kappa < 0 else cand + np.pi
    sgn = 1.0 if kappa < 0 else -1.0
    R = np.sin(w / 2) / (w / 2) if w > 0 else 1.0
    phase = c_curv + np.arctan2(sgn * t, sgn * rho + lam * R)
    c_axis = wrap24((cand + shift) * HOURS)
    a, b = wrap24(c_axis), wrap24(c_axis + 12.0)
    c_day = a if abs(cdiff(a, 13.0)) <= abs(cdiff(b, 13.0)) else b
    out = dict(sigma_f=S, c_axis12=np.mod(c_axis, 12.0), c_curv=wrap24((c_curv + shift) * HOURS), c_day=c_day, r=r, w=w_h,
               lam=lam, kappa=kappa, l1=l1, l2=l2, shift_h=shift * HOURS)
    out["donor_h"] = pd.Series(wrap24((phase + shift) * HOURS), index=donors)
    out["tangent"] = pd.Series(sgn * t, index=donors)  # tangent coordinate along the resolved direction
    # the same with the 12 h flip resolved by the hospital prior (centre nearer 13:00)
    use = cand if abs(cdiff(c_axis, 13.0)) <= abs(cdiff(c_axis + 12.0, 13.0)) else cand + np.pi
    u_c, u_t = np.array([np.cos(use), np.sin(use)]), np.array([-np.sin(use), np.cos(use)])
    rho_d, t_d = fc @ u_c, fc @ u_t
    out["donor_day_h"] = pd.Series(wrap24(wrap24((use + shift) * HOURS) + np.arctan2(t_d, rho_d + lam * R) * HOURS),
                                   index=donors)
    out["tangent_day"] = pd.Series(t_d, index=donors)
    return out


def host_estimate(ds: cs.Dataset, anchors=None, seed: int = 0) -> tuple[dict, cs.Fit]:
    """Anchored host of `ds` (`Dataset.fit`, no gate, no bootstrap) and its internal-time window."""
    anchors = list(anchors if anchors is not None else ds.anchors)
    sub = ds if anchors == list(ds.anchors) else ds.subset(anchors=anchors, compartments=anchors)
    ft = sub.fit(n_boot=0, gate=False, targets=anchors, seed=seed)
    tb = ft.table.loc[anchors]
    ref = max(anchors, key=lambda a: (sub.table.compartment == a).sum())
    shift = -np.deg2rad(tb.offset_h[ref] * 15)  # circular mean of arg z in the frame of the reference anchor
    t = sub.table
    strata = t.groupby("resection").stratum.first().reindex(ft.host.index).to_numpy() if "stratum" in t else None
    info = (tb.z2 / np.maximum(tb.psi, 1e-3)).reindex(anchors)
    prec = (t[t.compartment.isin(anchors)].assign(i=lambda x: x.compartment.map(info))
            .groupby("resection").i.sum().reindex(ft.host.index).fillna(0).to_numpy())
    f = ft.host[["x", "y"]].to_numpy()
    S = host_cov(f, strata, prec)
    return estimate(f, S, shift, strata, ft.host.index), ft


def boot(ds: cs.Dataset, anchors, point: dict, n_boot: int, seed: int) -> pd.DataFrame:
    """Donor bootstrap of the estimate (host refitted on `Dataset.resample` replicates is not possible with a gated fit,
    so each replicate rebuilds a Dataset from resampled pseudobulk rows); angles aligned to the point estimate."""
    rng = np.random.default_rng(seed)
    t = ds.table.rename(columns={"resection": "donor"})
    ids = np.array(sorted(set(t.donor)))
    by = {d: x for d, x in t.groupby("donor")}
    out = []
    for b in range(n_boot):
        pick = rng.choice(ids, len(ids))
        tb = pd.concat([by[d].assign(donor=f"{d}#{j}") for j, d in enumerate(pick)], ignore_index=True)
        try:
            e, _ = host_estimate(cs.Dataset(tb, anchors=list(anchors), name=ds.name), anchors)
        except Exception:
            continue
        out.append(dict(rep=b, r=e["r"], w=e["w"], lam=e["lam"], kappa=e["kappa"],
                        c_axis12=point["c_axis12"] + cdiff(2 * e["c_axis12"], 2 * point["c_axis12"]) / 2,
                        c_curv=point["c_curv"] + cdiff(e["c_curv"], point["c_curv"]),
                        flip_same=float(abs(cdiff(e["c_curv"], point["c_curv"])) < 6)))
    return pd.DataFrame(out)


def summarise(point: dict, reps: pd.DataFrame) -> dict:
    q = lambda c: reps[c].quantile([0.025, 0.975]).to_numpy() if len(reps) else (np.nan, np.nan)
    row = {k: v for k, v in point.items() if k not in ("donor_h", "tangent", "donor_day_h", "tangent_day", "sigma_f")}
    for c in ("r", "w", "c_axis12", "c_curv"):
        row[f"{c}_lo"], row[f"{c}_hi"] = q(c)
    row["c_axis12_se"] = reps.c_axis12.std() if len(reps) else np.nan
    row["flip_support"] = reps.flip_same.mean() if len(reps) else np.nan
    row["n_boot"] = len(reps)
    row["interpretable"] = bool(np.isfinite(row["r_hi"]) and row["r_hi"] < R_MAX)
    return row
