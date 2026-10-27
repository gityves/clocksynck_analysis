"""Score A139 from output/cache (or cache_quick with --quick): v139.csv (V139.1, V139.2 on the primary cell),
d139_3.csv (whole grid), failures.csv, gene_level.csv (false-flag rate and power per gene and flag type)."""
from __future__ import annotations

import argparse
import pickle
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import beta as beta_dist

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import lib  # noqa: E402

PRIMARY = ("08-18h", 80)
FLAGS = {"a_boot": "a_flag_boot", "a_lr": "a_flag_lr", "b_boot": "b_flag_boot"}


def ub(k, n, a=0.05):
    """One-sided 95 % Clopper-Pearson upper bound."""
    return 1.0 if k >= n else float(beta_dist.ppf(1 - a, k + 1, n - k))


def load(cache):
    G, S, F = [], [], []
    for f in sorted(cache.glob("*.pkl")):
        r = pickle.load(open(f, "rb"))
        if isinstance(r[0], str):
            F.append(r[1])
        else:
            G.append(r[0]); S.append(r[1])
    return (pd.concat(G, ignore_index=True) if G else pd.DataFrame(), pd.DataFrame(S), pd.DataFrame(F))


def gene_flags(G, S):
    """Per (window, n, scenario, gene, flag type): flags, replicates, rate, upper bound. LR is replicate-level."""
    rows = []
    keys = ["window", "n", "scenario"]
    for (w, n, sc), d in G.groupby(keys):
        for g, dg in d.groupby("gene", sort=False):
            for ft, col in FLAGS.items():
                if col not in dg or dg[col].isna().all():
                    continue
                k, m = int(dg[col].sum()), len(dg)
                rows.append(dict(window=w, n=n, scenario=sc, gene=g, shifted=bool(dg.shifted.iloc[0]), flag=ft,
                                 flags=k, reps=m, rate=k / m, upper95=ub(k, m)))
    return pd.DataFrame(rows)


def errors(G, scen, est, col_truth="truth_h"):
    d = G[G.scenario == scen]
    out = {}
    for lab, sel in (("shifted", d.shifted), ("unshifted", ~d.shifted)):
        x = d[sel]
        if len(x) and f"{est}_dev" in x:
            out[lab] = float(np.median(np.abs(lib.wrap(x[f"{est}_dev"] - x[col_truth]))))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    a = ap.parse_args()
    cache = lib.OUT / ("cache_quick" if a.quick else "cache")
    sfx = "_quick" if a.quick else ""
    G, S, F = load(cache)
    F.to_csv(lib.OUT / f"failures{sfx}.csv", index=False)
    if G.empty:
        return
    gf = gene_flags(G, S)
    gf.to_csv(lib.OUT / f"gene_level{sfx}.csv", index=False)

    # ---- D139.3: whole grid ----
    rows = []
    for (w, n, sc), d in G.groupby(["window", "n", "scenario"]):
        s = S[(S.window == w) & (S.n == n) & (S.scenario == sc)]
        r = dict(window=w, n=n, scenario=sc, reps=len(s))
        for est in ("a", "b"):
            if f"{est}_dev" not in d:
                continue
            d_ = d.dropna(subset=[f"{est}_dev"])
            if d_.empty:
                continue
            for lab, sel in (("shifted", d_.shifted), ("unshifted", ~d_.shifted)):
                if sel.any():
                    x = d_[sel]
                    r[f"{est}_mae_{lab}_h"] = float(np.median(np.abs(lib.wrap(x[f"{est}_dev"] - x.truth_h))))
                    r[f"{est}_mae_centred_{lab}_h"] = float(np.median(np.abs(lib.wrap(x[f"{est}_dev"] - x.truth_centred_h))))
            for ft in ("boot",) + (("lr",) if est == "a" else ()):
                col = f"{est}_flag_{ft}"
                if col in d_:
                    sh = gf[(gf.window == w) & (gf.n == n) & (gf.scenario == sc) & (gf.flag == f"{est}_{ft}")]
                    if ft == "lr":
                        r["a_lr_flag_rate"] = float(s.a_lr_p.lt(0.05).mean())
                    else:
                        r[f"{est}_boot_power_mean"] = float(sh[sh.shifted].rate.mean()) if sh.shifted.any() else np.nan
                        r[f"{est}_boot_power_min"] = float(sh[sh.shifted].rate.min()) if sh.shifted.any() else np.nan
                        r[f"{est}_boot_ffr_mean"] = float(sh[~sh.shifted].rate.mean())
                        r[f"{est}_boot_ffr_max"] = float(sh[~sh.shifted].rate.max())
        if "b_converged" in s:
            r["b_converged_share"] = float(s.b_converged.mean()); r["b_rounds_median"] = float(s.b_rounds.median())
        if sc in ("S3", "SW"):
            r["gate_pass_share"] = float(s.gate_pass.astype(float).mean()); r["gate_clock_share_median"] = float(s.gate_share.median())
        r["z_off_bias_h_E4"] = float(np.median(lib.wrap(s.a_z_off_h - 3.0)))
        r["z_log_amp_bias_E4"] = float(np.log(s.a_z_amp / 1.0).median())
        r["fit_off_bias_h"] = float(np.median(lib.wrap(s.fit_off_h - 3.0)))
        r["fit_log_amp_bias"] = float(s.fit_logamp.median())
        r["n_boot_ok_min"] = int(min(s.a_n_boot_ok.min(), s.get("b_n_boot_ok", s.a_n_boot_ok).min()))
        rows.append(r)
    grid = pd.DataFrame(rows)
    grid.to_csv(lib.OUT / f"d139_3{sfx}.csv", index=False)

    # ---- verdicts on the primary cell ----
    P = G[(G.window == PRIMARY[0]) & (G.n == PRIMARY[1])]
    Sp = S[(S.window == PRIMARY[0]) & (S.n == PRIMARY[1])]
    v = []
    est_of = {"S1": "a", "S3": "a", "SW": "b"}
    ok1 = []
    for sc, est in est_of.items():
        if not (P.scenario == sc).any():
            continue
        e = errors(P[P.scenario == sc], sc, est)
        # descriptive: the other estimator and the centred truth
        e_alt = errors(P[P.scenario == sc], sc, "a" if est == "b" else "b") if est == "b" else {}
        ec = errors(P[P.scenario == sc], sc, est, "truth_centred_h")
        p = e.get("shifted", np.nan) <= 1.5 and e.get("unshifted", np.nan) <= 1.0
        ok1.append(p)
        v.append(dict(test="V139.1", scenario=sc, estimator=f"E-{est}", mae_shifted_h=e.get("shifted"),
                      mae_unshifted_h=e.get("unshifted"), threshold_shifted=1.5, threshold_unshifted=1.0,
                      mae_shifted_centred_truth_h=ec.get("shifted"), mae_unshifted_centred_truth_h=ec.get("unshifted"),
                      mae_shifted_E_a_h=e_alt and None, reps=int((Sp.scenario == sc).sum()), pass_=bool(p)))
    # E-a in SW and E-b in S0 (descriptive)
    for sc, est in (("SW", "a"), ("S0", "a"), ("S0", "b")):
        if (P.scenario == sc).any() and f"{est}_dev" in P:
            e = errors(P[P.scenario == sc], sc, est)
            v.append(dict(test="V139.1-descriptive", scenario=sc, estimator=f"E-{est}", mae_shifted_h=e.get("shifted"),
                          mae_unshifted_h=e.get("unshifted"), reps=int((Sp.scenario == sc).sum())))
    v.append(dict(test="V139.1", scenario="all", estimator="E-a (S1,S3), E-b (SW)", pass_=bool(len(ok1) == 3 and all(ok1))))

    pf = gf[(gf.window == PRIMARY[0]) & (gf.n == PRIMARY[1])]
    ok2 = []
    for ft in ("a_boot", "b_boot"):
        d = pf[(pf.flag == ft) & (pf.scenario == "S0")]
        if len(d):
            p = bool((d.upper95 <= 0.05).all())
            ok2.append(p)
            v.append(dict(test="V139.2-false-flag", scenario="S0", estimator=ft, max_rate=d.rate.max(), max_upper95=d.upper95.max(),
                          genes_over=int((d.upper95 > 0.05).sum()), n_genes=len(d), reps=int(d.reps.iloc[0]), pass_=p,
                          detail="; ".join(f"{r.gene}:{r.flags}/{r.reps} ub={r.upper95:.3f}" for r in d.itertuples())))
    for sc, ft in (("S1", "a_boot"), ("S3", "a_boot"), ("SW", "b_boot")):
        d = pf[(pf.flag == ft) & (pf.scenario == sc) & pf.shifted]
        if len(d):
            p = bool((d.rate >= 0.8).all())
            ok2.append(p)
            v.append(dict(test="V139.2-power", scenario=sc, estimator=ft, power_min=d.rate.min(), power_mean=d.rate.mean(),
                          threshold=0.8, reps=int(d.reps.iloc[0]), pass_=p,
                          detail="; ".join(f"{r.gene}:{r.rate:.2f}" for r in d.itertuples())))
    for sc in ("S1", "S3", "SW"):  # descriptive: the other flag types
        for ft in ("a_boot", "b_boot"):
            d = pf[(pf.flag == ft) & (pf.scenario == sc) & pf.shifted]
            if len(d) and not ((sc, ft) in (("S1", "a_boot"), ("S3", "a_boot"), ("SW", "b_boot"))):
                v.append(dict(test="V139.2-power-descriptive", scenario=sc, estimator=ft, power_min=d.rate.min(), power_mean=d.rate.mean(),
                              reps=int(d.reps.iloc[0])))
    for sc in ("S0", "S1", "S3", "SW"):  # E4 likelihood-ratio flag (replicate level)
        s = Sp[Sp.scenario == sc]
        if len(s):
            k = int((s.a_lr_p < 0.05).sum())
            v.append(dict(test="E4-LR-flag", scenario=sc, estimator="a_lr", rate=k / len(s), upper95=ub(k, len(s)), reps=len(s),
                          detail="false-flag rate" if sc == "S0" else "power (any shift in the target)"))
    n_exp = 2 + 3  # S0 a_boot, S0 b_boot, three power rows
    v.append(dict(test="V139.2", scenario="all", estimator="a_boot (S0,S1,S3), b_boot (S0,SW)", pass_=bool(len(ok2) == n_exp and all(ok2))))
    pd.DataFrame(v).drop(columns=["mae_shifted_E_a_h"], errors="ignore").to_csv(lib.OUT / f"v139{sfx}.csv", index=False)
    print(pd.DataFrame(v).drop(columns=["detail"], errors="ignore").to_string())


if __name__ == "__main__":
    main()
