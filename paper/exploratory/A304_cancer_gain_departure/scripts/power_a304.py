"""A304 power stage (A2): what a pure gain loss and a pure departure rise look like on the paired TCGA host, and the
power of each readout at the real design of each contrast.

Registered in paper/docs/registered_experiments.md (A304, "Power stage"). Computes no observed contrast: the real data
enter only as arm sizes, per-sample Poisson variances and the all-pairs moments of each cancer (V_N, V_T, |C|: the A138
data already seen), from which S = rel V_N, g = |C| / S, 2 psi_N = V_N - S, 2 psi_T = V_T - g^2 S (floor 0.05 V_T), with
rel the GTEx coupling of the matching tissue. Per replicate, complex Gaussian states h ~ CN(0, S); N = h + e_N + Poisson,
T = g_i h + e_T,i + Poisson; the case arm gets g x k_G (gain loss) or 2 psi_T x k_D (departure rise).

    python power_a304.py [--reps 500] [--perms 500] [--jobs 4] [--smoke]
Outputs (output/power/): calibration.csv (per cancer S, g, psi, rel), design.csv (cells and arm sizes per contrast),
power.csv (per contrast x scenario: rejection of each readout at 0.05 and 0.025, mean Q_G / Q_R / Q_N), mde.csv.
"""

from __future__ import annotations

import argparse
import time

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

import lib_a304 as L

SCEN = [("null", 1.0, 1.0)] + [(f"gain x{k}", k, 1.0) for k in (0.85, 0.7, 0.5, 0.35)] + \
       [(f"departure x{k}", 1.0, k) for k in (1.25, 1.5, 2.0, 3.0)]
T0 = time.time()


def log(*a):
    print(f"[{(time.time() - T0) / 60:6.1f} min]", *a, flush=True)


def cn(rng, var, size):
    """Complex Gaussian with E|x|^2 = var (array or scalar)."""
    s = np.sqrt(np.asarray(var, float) / 2)
    return s * (rng.normal(size=size) + 1j * rng.normal(size=size))


def calibrate(states: dict, rel: dict) -> pd.DataFrame:
    rows = []
    for c, u in states.items():
        bt, vt, bn, vn = u["raw"]
        zt, zn = L.q.S.complex_state(bt), L.q.S.complex_state(bn)
        zt, zn = zt - zt.mean(), zn - zn.mean()
        n = len(zt)
        VN = float(np.sum(np.abs(zn) ** 2) / (n - 1) - vn.mean())
        VT = float(np.sum(np.abs(zt) ** 2) / (n - 1) - vt.mean())
        C = abs(np.sum(zt * np.conj(zn)) / (n - 1))
        S = rel[c] * VN
        g = C / S
        psiT2 = max(VT - g ** 2 * S, 0.05 * VT)
        rows.append(dict(cancer=c, n=n, rel=rel[c], V_N=VN, V_T=VT, abs_C=C, S=S, g=g, psi_N2=VN - S, psi_T2=psiT2,
                         psi_T2_floored=VT - g ** 2 * S < 0.05 * VT))
    return pd.DataFrame(rows).set_index("cancer")


def one_rep(cells, cal, kG, kD, n_perm, seed):
    rng = np.random.default_rng(seed)
    labs, zts, zns = [], [], []
    for cell in cells:
        k = cal.loc[cell["cancer"]]
        n, lab = cell["n"], cell["lab"]
        h = cn(rng, k.S, n)
        zn = h + cn(rng, k.psi_N2, n) + cn(rng, cell["vn"], n)
        g = np.where(lab == 1, k.g * kG, k.g)
        p2 = np.where(lab == 1, k.psi_T2 * kD, k.psi_T2)
        zt = g * h + cn(rng, p2, n) + cn(rng, cell["vt"], n)
        P = np.vstack([lab, L.perm_labels(lab, n_perm, rng)])
        labs.append(P); zts.append(zt); zns.append(zn)
    vals = [L.cell_values(c, P, zt, zn) for c, P, zt, zn in zip(cells, labs, zts, zns)]
    Q = L.pooled(cells, vals)
    out = {}
    for key in ("Q_G", "Q_R", "Q_N"):
        out[key] = Q[key][0]
    out["p_gain"] = L.p_less(Q["Q_G"][1:], Q["Q_G"][0])
    out["p_dep"] = L.p_greater(Q["Q_R"][1:], Q["Q_R"][0])
    out["p_host"] = min(1.0, 2 * min(L.p_less(Q["Q_N"][1:], Q["Q_N"][0]), L.p_greater(Q["Q_N"][1:], Q["Q_N"][0])))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=500)
    ap.add_argument("--perms", type=int, default=500)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--smoke", action="store_true")
    a = ap.parse_args()
    from clocksync.config import limit_threads
    limit_threads(1)
    if a.smoke:
        a.reps, a.perms = 10, 50
    out = L.OUT / ("power_smoke" if a.smoke else "power")
    out.mkdir(parents=True, exist_ok=True)
    samples = L.sh.load_tcga()
    states = L.q.tcga_states()
    assert all(u["gated"] for u in states.values())
    ann = L.annotations(samples)
    cal = calibrate(states, L.reliability())
    cal.to_csv(out / "calibration.csv")
    log("calibration\n", cal.round(3).to_string())
    drows, prow = [], []
    for ci, (h, name) in enumerate(L.CONTRASTS.items()):
        cells = L.design(states, ann, h)
        for c in cells:
            drows.append(dict(contrast=h, cancer=c["cancer"], n=c["n"], n_case=int((c["lab"] == 1).sum()),
                              n_ref=int((c["lab"] == 0).sum())))
        if not cells:
            continue
        for si, (sname, kG, kD) in enumerate(SCEN):
            res = Parallel(n_jobs=a.jobs)(delayed(one_rep)(cells, cal, kG, kD, a.perms, [304, 2, ci, si, r])
                                          for r in range(a.reps))
            r = pd.DataFrame(res)
            prow.append(dict(contrast=h, scenario=sname, k_G=kG, k_D=kD, reps=len(r),
                             power_gain_05=float((r.p_gain < 0.05).mean()), power_gain_025=float((r.p_gain < 0.025).mean()),
                             power_dep_05=float((r.p_dep < 0.05).mean()), power_dep_025=float((r.p_dep < 0.025).mean()),
                             host_flag_05=float((r.p_host < 0.05).mean()),
                             mean_Q_G=float(r.Q_G.mean()), sd_Q_G=float(r.Q_G.std()), mean_Q_R=float(r.Q_R.mean()),
                             sd_Q_R=float(r.Q_R.std()), mean_Q_N=float(r.Q_N.mean()), nan_Q_G=float(r.Q_G.isna().mean()),
                             nan_Q_R=float(r.Q_R.isna().mean())))
            log(h, sname, {k: round(v, 3) for k, v in prow[-1].items() if k.startswith(("power", "mean_Q"))})
    pd.DataFrame(drows).to_csv(out / "design.csv", index=False)
    pw = pd.DataFrame(prow)
    pw.to_csv(out / "power.csv", index=False)
    # MDE: smallest effect on the grid with power >= 0.8 at 0.025 (registered rule points: k_G 0.5, k_D 2)
    m = []
    for h, g in pw.groupby("contrast"):
        gg = g[g.scenario.str.startswith("gain")].sort_values("k_G", ascending=False)
        dd = g[g.scenario.str.startswith("departure")].sort_values("k_D")
        okg = gg[gg.power_gain_025 >= 0.8]
        okd = dd[dd.power_dep_025 >= 0.8]
        m.append(dict(contrast=h, mde_k_G=okg.k_G.max() if len(okg) else np.nan,
                      mde_k_D=okd.k_D.min() if len(okd) else np.nan,
                      power_gain_at_0_5=float(gg.loc[gg.k_G == 0.5, "power_gain_025"].iloc[0]),
                      power_dep_at_2=float(dd.loc[dd.k_D == 2.0, "power_dep_025"].iloc[0]),
                      typeI_gain=float(g.loc[g.scenario == "null", "power_gain_05"].iloc[0]),
                      typeI_dep=float(g.loc[g.scenario == "null", "power_dep_05"].iloc[0])))
    m = pd.DataFrame(m)
    m["gain_scored"] = m.power_gain_at_0_5 >= 0.8
    m["departure_scored"] = m.power_dep_at_2 >= 0.8
    m.to_csv(out / "mde.csv", index=False)
    log("MDE\n", m.round(3).to_string())


if __name__ == "__main__":
    main()
