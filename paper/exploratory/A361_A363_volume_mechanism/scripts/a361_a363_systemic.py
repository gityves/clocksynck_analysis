"""A361 (A360 on the core loop) and A363 part 1 (inflammatory state): which systemic states lower the person's volume?

GTEx fast deaths, A360 rows and covariates. Person volume = log host radius from all 25 tissues (A360), on the 14 genes
(A360, not blind) and on the core loop (PER1, PER2, NR1D1 removed; blind). Companions: 200 expression-matched and 200
responsive pseudo-clock waveforms (seeds [361, 1] / [361, 2]), each read on its 14 positions and its core positions,
with the same person covariates (the coefficients of a generic donor scale on the same persons).

    python a361_a363_systemic.py --stage power   # A360 reproduced; companions; analytic power; no core clock value
    python a361_a363_systemic.py --stage score   # registered tests (after the power note is committed)
Outputs (`../output/`): systemic_companions.csv, systemic_power.csv | systemic_tests.csv. Aggregates only.
"""
import argparse
import time

import numpy as np
import pandas as pd
from scipy.stats import norm

import vm_common as V
from clocksync.config import limit_threads
from clocksync.paired import pseudo_waveforms

N_SETS, N_PERM, N_BOOT = 200, 10_000, 2000
A360_MODEL = V.BASE + V.TECH + V.SYS
A363_MODEL = A360_MODEL + ["inf"]


def partial(per, log_r, block, given, seed):
    df = per.assign(log_r=log_r.reindex(per.index)).dropna(subset=["log_r"] + block + given)
    return V.a360.partial_r2(df, block, given, np.random.default_rng(seed), N_PERM)


def boot_coefs(per, log_rs, terms, seed):
    """Donor bootstrap of the standardised coefficients, shared draws across the volumes in `log_rs`."""
    df = per.dropna(subset=terms).copy()
    for k, v in log_rs.items():
        df[k] = v.reindex(df.index)
    df = df.dropna(subset=list(log_rs))
    rng = np.random.default_rng(seed)
    out = {k: [] for k in log_rs}
    for _ in range(N_BOOT):
        s = df.iloc[rng.integers(0, len(df), len(df))]
        for k in log_rs:
            out[k].append(V.a360.std_coefs(s.assign(log_r=s[k]), terms))
    return {k: np.array(v) for k, v in out.items()}, len(df)


def main():
    V.quiet()
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=["power", "score"], required=True)
    a = ap.parse_args()
    limit_threads(4)
    t0 = time.time()
    V.OUT.mkdir(parents=True, exist_ok=True)
    pb, _, gc, ieg = V.a359.load()
    inf = V.inflammatory_genes(pb.columns)
    comps = sorted(pb.compartment.unique())
    per = V.persons(pb, comps, gc, ieg, inf)
    print(f"{len(comps)} tissues, {len(per)} persons, inflammatory genes {len(inf)} ({time.time() - t0:.0f}s)", flush=True)

    # A360 reproduced first (not blind)
    v14 = V.volume(pb, comps, V.WF)
    ref = pd.read_csv(V.EXP / "A358_A360_gain_sources/output/a360_coefficients.csv")
    ref = ref[ref.model == "registered"].set_index("term").std_coef
    got = V.coefs(per, v14, A360_MODEL)
    if not np.allclose(got[ref.index].to_numpy(), ref.to_numpy(), atol=1e-9):
        raise SystemExit(f"A360 not reproduced: {got.round(4).to_dict()} vs {ref.round(4).to_dict()}")
    print("A360 reproduced (14 genes, registered model)", flush=True)

    if a.stage == "power":
        lm = V.lm_of(pb)
        rows = []
        for kind, tag in (("matched", 1), ("responsive", 2)):
            for j, w in enumerate(pseudo_waveforms(lm, N_SETS, np.random.default_rng([361, tag]), V.WF, kind=kind)):
                for gs, ww in (("14 genes", w), ("core", V.core_of(w, f"core{j}"))):
                    lr = V.volume(pb, comps, ww)
                    c60, c63 = V.coefs(per, lr, A360_MODEL), V.coefs(per, lr, A363_MODEL)
                    rows.append(dict(kind=kind, set=j, gene_set=gs, gc=c60.gc, ieg=c60.ieg, inf=c63.inf))
            print(f"{kind} done ({time.time() - t0:.0f}s)", flush=True)
        comp = pd.DataFrame(rows)
        comp.to_csv(V.OUT / "systemic_companions.csv", index=False)
        # analytic power (standardised scale)
        a360c = pd.read_csv(V.EXP / "A358_A360_gain_sources/output/a360_coefficients.csv")
        a360c = a360c[a360c.model == "registered"].set_index("term")
        se_gc = float((a360c.loc["gc", "hi90"] - a360c.loc["gc", "lo90"]) / (2 * 1.645))  # A360 bootstrap 90 % width
        df = per.dropna(subset=A363_MODEL)
        X = df[A360_MODEL].to_numpy(float)
        X = np.column_stack([np.ones(len(X)), (X - X.mean(0)) / X.std(0, ddof=1)])
        y = (df.inf - df.inf.mean()) / df.inf.std(ddof=1)
        r2_inf = 1 - float(((y - X @ np.linalg.lstsq(X, y, rcond=None)[0]) ** 2).sum()) / float((y ** 2).sum())
        r2_full = float(pd.read_csv(V.EXP / "A358_A360_gain_sources/output/a360_tests.csv").set_index("test").loc["full model R2", "partial_r2"])
        se_inf = float(np.sqrt((1 - r2_full) / (len(df) * (1 - r2_inf))))
        out = {}
        for gs in ("core", "14 genes"):
            c = comp[comp.gene_set == gs]
            for kind in ("matched", "responsive"):
                for term in ("gc", "inf"):
                    v = c[c.kind == kind][term]
                    out[f"{gs}|{kind}|{term}|median"] = float(v.median())
                    out[f"{gs}|{kind}|{term}|q05"] = float(v.quantile(0.05))
        out.update(persons=len(df), se_gc_a360=se_gc, se_inf_analytic=se_inf, r2_inf_given_covariates=r2_inf)
        thr_gc = min(out["core|matched|gc|q05"], -1.96 * se_gc)
        thr_inf = min(out["core|matched|inf|q05"], -1.96 * se_inf)
        for lab, alt in (("full", -0.235), ("half", -0.1175)):
            out[f"A361_power_{lab}"] = float(norm.cdf((thr_gc - alt) / se_gc))
        b_inf = -np.sqrt(0.03 * (1 - r2_full) / (1 - r2_inf))   # partial R^2 0.03 in standardised units
        for lab, alt in (("pr2_0.03", b_inf), ("pr2_0.015", b_inf / np.sqrt(2))):
            out[f"A363.1_power_{lab}"] = float(norm.cdf((thr_inf - alt) / se_inf))
        out.update(A361_threshold=thr_gc, A363_1_threshold=thr_inf, A363_1_alt_coef=float(b_inf))
        out["A361_gate"] = bool(out["A361_power_full"] >= 0.8)
        out["A363_1_gate"] = bool(out["A363.1_power_pr2_0.03"] >= 0.8)
        pd.Series(out).to_csv(V.OUT / "systemic_power.csv", header=["value"])
        print(pd.Series(out).to_string(), f"\ndone {time.time() - t0:.0f}s", flush=True)
        return

    # ---- score ---------------------------------------------------------------------------------------------------------
    pw = pd.read_csv(V.OUT / "systemic_power.csv", index_col=0)["value"]
    comp = pd.read_csv(V.OUT / "systemic_companions.csv")
    vc = V.volume(pb, comps, V.CORE)
    rows = []
    if str(pw["A361_gate"]) == "True":
        r2, p = partial(per, vc, V.SYS, V.BASE + V.TECH, 361)
        c = V.coefs(per, vc, A360_MODEL)
        bt, n = boot_coefs(per, {"core": vc, "14": v14}, A360_MODEL, [361, 3])
        k = A360_MODEL.index("gc")
        m = comp[(comp.gene_set == "core") & (comp.kind == "matched")].gc
        r = comp[(comp.gene_set == "core") & (comp.kind == "responsive")].gc
        row = dict(test="H361.1", persons=n, partial_r2=r2, p_perm=p, gc_coef=c.gc, gc_lo90=np.quantile(bt["core"][:, k], .05),
                   gc_hi90=np.quantile(bt["core"][:, k], .95), ieg_coef=c.ieg, matched_q05=m.quantile(.05),
                   pct_matched=float((m < c.gc).mean() * 100), responsive_q05=r.quantile(.05),
                   pct_responsive=float((r < c.gc).mean() * 100), gc_coef_14=V.coefs(per, v14, A360_MODEL).gc)
        d = bt["14"][:, k] - bt["core"][:, k]
        row.update(diff_14_minus_core=row["gc_coef_14"] - c.gc, diff_lo90=np.quantile(d, .05), diff_hi90=np.quantile(d, .95))
        row["passes"] = bool(p < 0.05 and c.gc < 0 and c.gc < row["matched_q05"])
        row["reading"] = ("core: systemic and beyond a generic donor scale" if row["passes"] else
                          "core: not beyond matched genes" if p < 0.05 and c.gc < 0 else "core: not detected")
        rows.append(row)
    if str(pw["A363_1_gate"]) == "True":
        for gs, v in (("core", vc), ("14 genes", v14)):
            r2, p = partial(per, v, ["inf"], A360_MODEL, 363)
            c = V.coefs(per, v, A363_MODEL)
            bt, n = boot_coefs(per, {gs: v}, A363_MODEL, [363, 3])
            m = comp[(comp.gene_set == gs) & (comp.kind == "matched")].inf
            r = comp[(comp.gene_set == gs) & (comp.kind == "responsive")].inf
            row = dict(test="H363.1" if gs == "core" else "H363.1 (14 genes, reported)", persons=n, partial_r2=r2,
                       p_perm=p, inf_coef=c.inf, inf_lo90=np.quantile(bt[gs][:, -1], .05), inf_hi90=np.quantile(bt[gs][:, -1], .95),
                       gc_coef=c.gc, matched_q05=m.quantile(.05), pct_matched=float((m < c.inf).mean() * 100),
                       responsive_q05=r.quantile(.05), pct_responsive=float((r < c.inf).mean() * 100))
            row["passes"] = bool(p < 0.05 and c.inf < 0 and c.inf < row["matched_q05"])
            rows.append(row)
    t = pd.DataFrame(rows)
    t.to_csv(V.OUT / "systemic_tests.csv", index=False)
    print(t.round(4).T.to_string(), f"\ndone {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
