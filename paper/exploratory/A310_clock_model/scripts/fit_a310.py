"""A310 Stage 2a: fit of the mean-field population model to blocks A (GTEx personal departures), B (cell-class gains
and phases) and C (liver feeding drive), the registered adequacy criterion, identifiability profiles and the
alternatives. Not blind, descriptive (registration in `paper/docs/registered_experiments.md`).

Locked state with nu = 0 and a real drive: R solves -gamma R + (K R + F)(1 - R^2) / 2 = 0, so a tissue's input R_o
fixes the varying constant in closed form (gamma_o, K_o or F_o). Departure covariance from `model.departure` logic
(Lyapunov at the locked state).

    python fit_a310.py
Inputs: ../output/boot.csv (boot_a310.py), A294b tissues.csv, A293 Part 2-3 atlases.csv, A297 tests.csv.
Outputs (`../output/`): blockA_tissues.csv, blockA_variants.csv, blockB_gains.csv, blockB_phases.csv, blockC.csv,
profiles.csv, composition_bound.csv, verdict.csv.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.linalg import solve_continuous_lyapunov
from scipy.optimize import fsolve, minimize
from scipy.stats import theilslopes

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parent))
import model  # noqa: E402

EXP = HERE.parents[2]
OUT = HERE.parents[1] / "output"
SLOPE_OBS = (0.04437239247573796, 0.33155911954731887)  # A294b H294b.2 90 % interval (A294_floor/output/tests.csv)
Z90 = 1.6448536269514722
RNG = np.random.default_rng(310)


# ----------------------------------------------------------------------------------------------- block A
def departure_at(R, variant, f):
    """Total departure variance (q = 1), radial and tangential parts, for a tissue at synchrony R (array).
    Variant: the constant that varies between tissues; the shared ones are 1 (time unit) and f."""
    R = np.asarray(R, float)
    if variant == "dispersion":
        K, F = 1.0, f
        g = (K * R + F) * (1 - R ** 2) / (2 * R)
        pars = [(gi, K, F) for gi in g]
    elif variant == "coupling":
        g, F = 1.0, f
        K = (2 * g * R - F * (1 - R ** 2)) / (R * (1 - R ** 2))
        pars = [(g, Ki, F) for Ki in K]
    else:  # drive
        g, K = 1.0, f
        F = (2 * g * R - K * R * (1 - R ** 2)) / (1 - R ** 2)
        pars = [(g, K, Fi) for Fi in F]
    out = []
    for (g_, K_, F_), r in zip(pars, R):
        if g_ <= 0 or K_ < 0 or F_ <= 0 or not 0 < r < 1:
            out.append((np.nan, np.nan))
            continue
        J = model.jacobian(complex(r, 0.0), g_, K_, F_)
        if np.max(np.linalg.eigvals(J).real) >= 0:
            out.append((np.nan, np.nan))
            continue
        C = solve_continuous_lyapunov(J, -np.eye(2))
        out.append((C[0, 0], C[1, 1]))  # w* real: x = radial, y = tangential
    return np.array(out)


def block_a_predict(theta, variant, S):
    logQ, Rmax_logit, logf = theta
    Rmax = 1 / (1 + np.exp(-Rmax_logit))
    R = Rmax * np.sqrt(S / S.max())
    d = departure_at(R, variant, np.exp(logf))
    # bulk units: S_o = c R_o^2 with c = S_max / Rmax^2; E_o = c q (rad + tan) -> one scale Q
    return np.exp(logQ) * d.sum(1), d, R


def block_a_sse(theta, variant, S, logE, se):
    E, _, _ = block_a_predict(theta, variant, S)
    if not np.all(np.isfinite(E)) or np.any(E <= 0):
        return 1e9
    return float(np.sum(((np.log(E) - logE) / se) ** 2))


def fit_block_a(S, logE, se, variant):
    best = None
    for Rl in (-2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 6.0):
        for lf in (-4.0, -3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0):
            E0 = block_a_predict([0.0, Rl, lf], variant, S)[0]
            if not np.all(np.isfinite(E0)):
                continue
            x0 = [float(np.mean(logE) - np.mean(np.log(E0))), Rl, lf]
            r = minimize(block_a_sse, x0, args=(variant, S, logE, se), method="Nelder-Mead",
                         options=dict(xatol=1e-6, fatol=1e-8, maxiter=20000))
            if best is None or r.fun < best.fun:
                best = r
    return best


def profile(fun, theta, i, lo, hi, n=41):
    """Profile of parameter i over [lo, hi] (other parameters re-optimised); returns grid and minimum SSE."""
    grid = np.linspace(lo, hi, n)
    vals = []
    for v in grid:
        def f(rest):
            t = np.insert(rest, i, v)
            return fun(t)
        rest0 = np.delete(theta, i)
        r = minimize(f, rest0, method="Nelder-Mead", options=dict(xatol=1e-5, fatol=1e-7, maxiter=4000))
        vals.append(r.fun)
    return grid, np.array(vals)


def identifiable(grid, vals, est, fmin):
    left, right = vals[grid < est], vals[grid > est]
    return bool(len(left) and len(right) and left.max() - fmin >= 2.71 and right.max() - fmin >= 2.71)


# ----------------------------------------------------------------------------------------------- block C
def liver_state(F, gamma, K, w0=None):
    """Locked state for a constant complex drive (nu = 0). The equation is rotation-covariant, so w* = R(|F|) e^{i arg F}
    with R the unique positive root of -K R^3 / 2 - |F| R^2 / 2 + (K / 2 - gamma) R + |F| / 2 = 0 (w0 unused)."""
    a = abs(F)
    if not np.all(np.isfinite([a, gamma, K])) or K <= 0 or gamma <= 0:
        return complex(np.nan, np.nan)
    roots = np.roots([-K / 2, -a / 2, K / 2 - gamma, a / 2])
    R = roots[(abs(roots.imag) < 1e-9) & (roots.real > 0)].real
    return complex(R.max() * np.exp(1j * np.angle(F))) if len(R) else complex(np.nan, np.nan)


def block_c_predict(theta):
    Ff, alpha, gamma, K = np.exp(theta[0]), theta[1], np.exp(theta[2]), np.exp(theta[3])
    al = liver_state(1 + Ff * np.exp(1j * alpha), gamma, K, 0.5 + 0j)
    if not np.isfinite(al):
        return None
    inv = liver_state(1 + Ff * np.exp(1j * (alpha + np.pi)), gamma, K, al)
    rem = liver_state(1 + 0j, gamma, K, al)
    if not (np.isfinite(inv) and np.isfinite(rem)):
        return None

    def contrast(w):  # amplitude ratio, acrophase shift in h (later = positive, mod 24)
        return abs(w) / abs(al), (-(np.angle(w) - np.angle(al)) * 24 / (2 * np.pi)) % 24
    return dict(inv=contrast(inv), rem=contrast(rem), nr=(1.0, 0.0))


def linear_c_predict(theta):
    Ff, alpha = np.exp(theta[0]), theta[1]
    al = 1 + Ff * np.exp(1j * alpha)
    inv, rem = 1 + Ff * np.exp(1j * (alpha + np.pi)), 1 + 0j

    def contrast(w):
        return abs(w) / abs(al), (-(np.angle(w) - np.angle(al)) * 24 / (2 * np.pi)) % 24
    return dict(inv=contrast(inv), rem=contrast(rem), nr=(1.0, 0.0))


def c_table():
    t = pd.read_csv(EXP / "A297_feeding_model" / "output" / "tests.csv")
    keep = {"Manella": "inv", "Xin": "inv", "observed Greenwell": "rem", "control Greenwell": "nr"}
    rows = []
    for k, arm in keep.items():
        r = t[t.contrast.str.startswith(k)].iloc[0]
        rows.append(dict(contrast=r.contrast, arm=arm, amp=r.amp_ratio, amp_lo=r.amp_lo90, amp_hi=r.amp_hi90,
                         shift=r.shift_h, shift_lo=r.shift_lo90, shift_hi=r.shift_hi90))
    return pd.DataFrame(rows)


def wrap(d):
    return (d + 12) % 24 - 12


def c_sse(theta, tab, pred_fun):
    p = pred_fun(theta)
    if p is None:
        return 1e9
    s = 0.0
    for r in tab.itertuples():
        a, h = p[r.arm]
        se_a = (np.log(r.amp_hi) - np.log(r.amp_lo)) / (2 * Z90)
        se_h = (r.shift_hi - r.shift_lo) / (2 * Z90)
        s += ((np.log(a) - np.log(r.amp)) / se_a) ** 2 + (wrap(h - r.shift) / se_h) ** 2
    return float(s)


def c_hits(theta, tab, pred_fun):
    p = pred_fun(theta)
    rows = []
    for r in tab.itertuples():
        a, h = p[r.arm]
        rows.append(dict(contrast=r.contrast, amp_obs=r.amp, amp_pred=a, amp_hit=bool(r.amp_lo <= a <= r.amp_hi),
                         shift_obs=r.shift, shift_pred=h,
                         shift_hit=bool(wrap(h - r.shift_lo) >= 0 and wrap(r.shift_hi - h) >= 0)))
    return pd.DataFrame(rows)


def fit_c(tab, pred_fun, k):
    best = None
    starts = RNG.uniform([-2, -np.pi, -3, -3][:k], [2, np.pi, 2, 2][:k], size=(30, k))
    for x0 in starts:
        r = minimize(c_sse, x0, args=(tab, pred_fun), method="Nelder-Mead",
                     options=dict(xatol=1e-6, fatol=1e-8, maxiter=8000))
        if best is None or r.fun < best.fun:
            best = r
    return best


# ----------------------------------------------------------------------------------------------- main
def main():
    OUT.mkdir(parents=True, exist_ok=True)
    verdict = {}
    # ---------------- block A
    floor = pd.read_csv(EXP / "A294_floor" / "output" / "tissues.csv", index_col=0)
    boot = pd.read_csv(OUT / "boot.csv").dropna(subset=["own"])
    boot["E"] = boot.own - boot.tissue.map(floor.own_matched_median)
    full = boot[boot.replicate == 0].set_index("tissue")
    reps = boot[boot.replicate > 0]
    ta = pd.DataFrame(dict(S=full.S, E=full.E))
    ta["E_lo"] = reps.groupby("tissue").E.quantile(0.05)
    ta["E_hi"] = reps.groupby("tissue").E.quantile(0.95)
    ta["se_log"] = reps.assign(l=np.log(reps.E.clip(lower=1e-3))).groupby("tissue").l.std()
    ta["boot_E_nonpos"] = reps.groupby("tissue").E.apply(lambda e: float((e <= 0).mean()))
    print(f"block A: {len(ta)} tissues, {reps.replicate.nunique()} replicates", flush=True)
    S, logE, se = ta.S.to_numpy(), np.log(ta.E.to_numpy()), ta.se_log.to_numpy()
    obs_slope = theilslopes(logE, np.log(S))[0]
    vrows, fits = [], {}
    for v in ("dispersion", "coupling", "drive"):
        r = fit_block_a(S, logE, se, v)
        if r is None:
            vrows.append(dict(variant=v, sse=np.inf, hit_share=0.0, model_slope=np.nan, slope_in_interval=False))
            continue
        E, d, R = block_a_predict(r.x, v, S)
        ok = np.isfinite(E)
        hit = (E >= ta.E_lo.to_numpy()) & (E <= ta.E_hi.to_numpy())
        slope = theilslopes(np.log(E[ok]), np.log(S[ok]))[0] if ok.sum() > 2 else np.nan
        vrows.append(dict(variant=v, sse=r.fun, hit_share=float(hit.mean()), model_slope=slope,
                          slope_in_interval=bool(SLOPE_OBS[0] <= slope <= SLOPE_OBS[1]), logQ=r.x[0],
                          Rmax=1 / (1 + np.exp(-r.x[1])), f=np.exp(r.x[2]), R_min=float(np.nanmin(R)),
                          share_median=float(np.nanmedian((d[:, 1] - d[:, 0]) / d.sum(1)))))
        fits[v] = (r, E, d, R, hit)
    va = pd.DataFrame(vrows)
    va["meets"] = (va.hit_share >= 0.8) & va.slope_in_interval
    va.to_csv(OUT / "blockA_variants.csv", index=False)
    print(va.round(3).to_string(), flush=True)
    print(f"observed Theil-Sen slope of log E on log S (bootstrap-input E): {obs_slope:.3f}", flush=True)
    chosen = va[va.meets].sort_values("sse").variant.iloc[0] if va.meets.any() else None
    verdict["A_variant"] = chosen
    verdict["A_pass"] = chosen is not None
    best_v = chosen or va.sort_values("sse").variant.iloc[0]
    r, E, d, R, hit = fits[best_v]
    ta = ta.assign(R=R, E_pred=E, hit=hit, var_rad=d[:, 0], var_tan=d[:, 1],
                   share_pred=(d[:, 1] - d[:, 0]) / d.sum(1))
    ta.to_csv(OUT / "blockA_tissues.csv")
    prof = []
    names = ["logQ", "Rmax_logit", "logf"]
    fun = lambda t: block_a_sse(t, best_v, S, logE, se)  # noqa: E731
    for i, nm in enumerate(names):
        est = r.x[i]
        span = np.log(100) if nm != "Rmax_logit" else 6.0
        g, vals = profile(fun, r.x, i, est - span, est + span)
        prof.append(dict(block="A", variant=best_v, parameter=nm, estimate=est, sse_min=r.fun,
                         identifiable=identifiable(g, vals, est, r.fun),
                         lo90=float(g[(vals - r.fun) < 2.71].min()), hi90=float(g[(vals - r.fun) < 2.71].max())))
    # ---------------- block B
    gb = pd.read_csv(EXP / "A293_organisation_levels" / "output" / "part3" / "atlases.csv")
    fam = pd.read_csv(EXP / "A293_organisation_levels" / "output" / "part3" / "families.csv").set_index("family")
    w = 1 / gb.se ** 2
    theta = (gb.g * w).groupby(gb.family).sum() / w.groupby(gb.family).sum()
    gb["fit"] = gb.family.map(theta)
    gb["within"] = (gb.g - gb.fit).abs() <= Z90 * gb.se
    gb["within_re_pooled"] = (gb.g - gb.family.map(fam.pooled)).abs() <= Z90 * gb.se
    gb.to_csv(OUT / "blockB_gains.csv", index=False)
    ph = pd.read_csv(EXP / "A293_organisation_levels" / "output" / "part2" / "atlases.csv")
    lead_rows = ph.class_pair.str.startswith("endothelial")
    lam = float(ph.offset_h[lead_rows].mean())
    ph["pred"] = np.where(lead_rows, lam, 0.0)
    ph["resid"] = ph.offset_h - ph.pred
    rms = float(np.sqrt(np.mean(ph.resid ** 2)))
    ph.to_csv(OUT / "blockB_phases.csv", index=False)
    verdict.update(B_gain_within=float(gb.within.mean()), B_gain_within_re=float(gb.within_re_pooled.mean()),
                   B_lead_h=lam, B_phase_rms_h=rms, B_phase_rms_no_lead_h=float(np.sqrt(np.mean(ph.offset_h ** 2))))
    verdict["B_pass"] = bool(verdict["B_gain_within"] >= 0.8 and rms <= 0.6)
    print(f"block B: gains within {verdict['B_gain_within']:.3f} (RE pooled {verdict['B_gain_within_re']:.3f}); "
          f"lead {lam:.2f} h, RMS {rms:.2f} h", flush=True)
    # ---------------- block C
    tab = c_table()
    rc = fit_c(tab, block_c_predict, 4)
    hc = c_hits(rc.x, tab, block_c_predict).assign(model="population")
    rl = fit_c(tab, linear_c_predict, 2)
    hl = c_hits(rl.x, tab, linear_c_predict).assign(model="linear")
    hits_c = int(hc.amp_hit.sum() + hc.shift_hit.sum())
    hits_l = int(hl.amp_hit.sum() + hl.shift_hit.sum())
    pd.concat([hc, hl]).to_csv(OUT / "blockC.csv", index=False)
    pc = dict(zip(["log_Ff", "alpha", "log_gamma", "log_K"], rc.x))
    al = liver_state(1 + np.exp(pc["log_Ff"]) * np.exp(1j * pc["alpha"]), np.exp(pc["log_gamma"]), np.exp(pc["log_K"]), 0.5)
    verdict.update(C_hits=hits_c, C_hits_linear=hits_l, C_sse=rc.fun, C_sse_linear=rl.fun, C_pass=bool(hits_c >= 6),
                   C_R_adlib=float(abs(al)), C_K_over_2gamma=float(np.exp(pc["log_K"]) / (2 * np.exp(pc["log_gamma"]))),
                   **{f"C_{k}": v for k, v in pc.items()})
    print(pd.concat([hc, hl]).round(3).to_string(), flush=True)
    funC = lambda t: c_sse(t, tab, block_c_predict)  # noqa: E731
    for i, nm in enumerate(pc):
        est = rc.x[i]
        span = np.pi if nm == "alpha" else np.log(100)
        g, vals = profile(funC, rc.x, i, est - span, est + span, n=25)
        prof.append(dict(block="C", variant="liver", parameter=nm, estimate=est, sse_min=rc.fun,
                         identifiable=identifiable(g, vals, est, rc.fun),
                         lo90=float(g[(vals - rc.fun) < 2.71].min()), hi90=float(g[(vals - rc.fun) < 2.71].max())))
    pd.DataFrame(prof).to_csv(OUT / "profiles.csv", index=False)
    print(pd.DataFrame(prof).round(3).to_string(), flush=True)
    # ---------------- alternative (b): composition alone
    ratio = float(np.exp(fam.pooled.max() - min(fam.pooled.min(), 0.0)))
    cb = pd.DataFrame(dict(S=ta.S, bound=np.sqrt(ta.S / ta.S.max()) * ratio)).sort_values("bound")
    cb.to_csv(OUT / "composition_bound.csv")
    verdict["b_class_gain_ratio"] = ratio
    verdict["b_tissues_bound_below_half"] = int((cb.bound < 0.5).sum())
    # ---------------- global
    n_par, n_obs = 4 + 6 + 4, 28 + 59 + 8
    verdict.update(n_par=n_par, n_obs=n_obs, global_pass=bool(n_par < n_obs / 3))
    verdict["verdict"] = ("ADEQUATE" if verdict["A_pass"] and verdict["B_pass"] and verdict["C_pass"] and verdict["global_pass"]
                          else "PARTIAL" if verdict["A_pass"] and verdict["B_pass"] else "INADEQUATE")
    pd.Series(verdict).to_csv(OUT / "verdict.csv", header=["value"])
    print(pd.Series(verdict).to_string(), flush=True)


if __name__ == "__main__":
    main()
