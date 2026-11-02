"""A193 scoring (registered in `paper/docs/registered_experiments.md`, commit a6d5922): does high CYCLOPS magnitude predict
worse survival in luminal A breast cancer (Li et al., PNAS 2024)?

Readout: `output/cyclops_magnitude.csv` (`cyclops_magnitude.py`; mean of within-fit z-scored magnitudes over three fits).
Reproduction (TCGA, not counted): the authors' "died within 5 years" (death at <= 60 months = 1, every other patient 0, as in
their script), top tertile vs the rest, logistic with age and metastasis at diagnosis, then + MKI67; Cox OS (top tertile and
per SD) adjusted for age, T, N and G2M.
H193.1: METABRIC + SCAN-B luminal A, OS, top tertile (within cohort) vs the rest, Cox stratified by cohort, age + G2M.
H193.2: magnitude per SD, cohort-specific clinical covariates + G2M. Holm over two; replicated if Holm p < 0.05 and HR > 1.
Power step first. Only summary statistics are written.

    python score_a193.py [--power-only]
"""
import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import statsmodels.api as sm  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

import outcomes as O  # noqa: E402
from score_a191 import Z80, cox, holm, metabric_readout  # noqa: E402
from clocksync.decompose import sample_log_odds  # noqa: E402

OUT = HERE.parents[1] / "output"
LN2 = np.log(2)
CLIN = {"SCAN-B": ["age", "log_size", "node_pos", "grade3", "grade1"], "METABRIC": ["age", "log_size", "nodes", "grade3", "grade1"]}


def readout(cohort):
    m = pd.read_csv(OUT / "cyclops_magnitude.csv")
    m = m[(m.cohort == cohort) & (m.is_tumour.astype(bool))].set_index("sample")
    return m


def cohort_table(k):
    c, X = O.LOADERS[k]()
    r = readout(k)
    c = c.join(r[["magnitude"]], how="inner").assign(cohort=k)
    c = c[c.pam50 == "LumA"].copy()
    c["mag_z"] = (c.magnitude - c.magnitude.mean()) / c.magnitude.std()
    q = pd.qcut(c.magnitude.rank(method="first"), 3, labels=False)
    c["mag_top"] = (q == 2).astype(float)
    return c, X, r


def run(parts, cov, x, extra=(), tcol="time", ecol="event"):
    d = pd.concat([p.dropna(subset=[tcol, ecol, x] + cov.get(p.cohort.iloc[0], []) + list(extra)) for p in parts])
    Z = O.block_design(d, cov, extra)
    r = cox(d[tcol].to_numpy(float), d[ecol].to_numpy(float), d[x].to_numpy(float), Z, strata=d.cohort.to_numpy())
    if set(d[x].dropna().unique()) <= {0.0, 1.0}:  # report the binary contrast, not per SD of the indicator
        sd = d[x].std()
        for k in ("hr", "hr_lo", "hr_hi"):
            r[k] = r[k] ** (1 / sd)
    return d, Z, dict(n=len(d), events=int(d[ecol].sum()), **{k: v for k, v in r.items() if k != "llf"})


def mde_unit(x, Z, events):
    """Schoenfeld MDE of the HR per unit of x (binary: the contrast), with 1 - R^2 of x on Z."""
    r2 = sm.OLS(x, sm.add_constant(Z)).fit().rsquared if Z is not None else 0.0
    return float(np.exp(Z80 / np.sqrt(events * np.var(x) * (1 - r2)))), float(r2)


def logistic(d, cols):
    y = d.died5.to_numpy(float)
    X = sm.add_constant(d[cols].to_numpy(float))
    r = sm.Logit(y, X).fit(disp=0)
    return dict(n=len(d), deaths_5y=int(y.sum()), or_top=float(np.exp(r.params[1])), p=float(r.pvalues[1]))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--power-only", action="store_true")
    a = ap.parse_args()
    tabs = {k: cohort_table(k) for k in ("METABRIC", "SCAN-B", "TCGA-BRCA")}
    prim = [tabs["METABRIC"][0], tabs["SCAN-B"][0]]
    cov1 = {k: ["age"] for k in ("METABRIC", "SCAN-B")}
    cov2 = CLIN
    power = []
    for test, cov, x, lim in (("H193.1", cov1, "mag_top", 2.0), ("H193.2", cov2, "mag_z", 1.5)):
        d = pd.concat([p.dropna(subset=["time", "event", x, "g2m"] + cov[p.cohort.iloc[0]]) for p in prim])
        Z = np.column_stack([O.block_design(d, cov, ("g2m",)), (d.cohort == "SCAN-B").to_numpy(float)])
        md, r2 = mde_unit(d[x].to_numpy(float), Z, d.event.sum())
        power.append(dict(test=test, n=len(d), events=int(d.event.sum()), r2_on_covariates=r2, mde_hr=md, limit=lim,
                          estimable=md <= lim))
    power = pd.DataFrame(power)
    power.to_csv(OUT / "a193_power.csv", index=False)
    print(power.round(4).to_string(index=False), flush=True)
    if a.power_only:
        return

    tests = []
    for test, cov, x in (("H193.1 top tertile, age + G2M", cov1, "mag_top"), ("H193.2 per SD, clinical + G2M", cov2, "mag_z")):
        tests.append(dict(test=test, **run(prim, cov, x, ("g2m",))[2]))
    tt = pd.DataFrame(tests)
    tt["p_holm"] = holm(tt.p)
    tt["replicated"] = (tt.p_holm < .05) & (tt.hr > 1)
    tt.to_csv(OUT / "a193_tests.csv", index=False)

    # ---- reproduction in TCGA (the authors' cohort; not counted) ----
    t = tabs["TCGA-BRCA"][0]
    t["died5"] = ((t.event == 1) & (t.time <= 60)).astype(float)
    rep = []
    tl = t.dropna(subset=["age", "time"])
    rep.append(dict(test="TCGA died <= 5 y ~ top tertile + age + M1 (their model)", **logistic(tl, ["mag_top", "age", "m1"])))
    rep.append(dict(test="TCGA died <= 5 y ~ top tertile + age + M1 + MKI67", **logistic(tl, ["mag_top", "age", "m1", "MKI67"])))
    rep.append(dict(test="TCGA died <= 5 y ~ top tertile", **logistic(tl, ["mag_top"])))
    tc = {"TCGA-BRCA": ["age", "t2", "t34", "node_pos"]}
    for lab, x, extra in (("top tertile", "mag_top", ()), ("top tertile + G2M", "mag_top", ("g2m",)),
                          ("per SD", "mag_z", ()), ("per SD + G2M", "mag_z", ("g2m",))):
        rep.append(dict(test=f"TCGA Cox OS, {lab}, age + T + N", **run([t], tc, x, extra)[2]))
    rep.append(dict(test="TCGA Cox OS, top tertile, unadjusted", **run([t], {"TCGA-BRCA": []}, "mag_top")[2]))
    rep = pd.DataFrame(rep)
    rep.to_csv(OUT / "a193_reproduction.csv", index=False)

    # ---- secondary ----
    sec = []
    for p in prim:
        k = p.cohort.iloc[0]
        sec.append(dict(test=f"{k} top tertile, age + G2M", **run([p], {k: ["age"]}, "mag_top", ("g2m",))[2]))
        sec.append(dict(test=f"{k} per SD, clinical + G2M", **run([p], {k: cov2[k]}, "mag_z", ("g2m",))[2]))
        sec.append(dict(test=f"{k} top tertile, unadjusted", **run([p], {k: []}, "mag_top")[2]))
    sec.append(dict(test="pooled top tertile, age + MKI67", **run(prim, cov1, "mag_top", ("MKI67",))[2]))
    sec.append(dict(test="pooled top tertile, age only", **run(prim, cov1, "mag_top")[2]))
    mb = prim[0]
    sec.append(dict(test="METABRIC RFS top tertile, age + G2M", **run([mb], {"METABRIC": ["age"]}, "mag_top", ("g2m",), "rfs_t", "rfs_e")[2]))
    sec.append(dict(test="METABRIC RFS per SD, clinical + G2M", **run([mb], {"METABRIC": cov2["METABRIC"]}, "mag_z", ("g2m",), "rfs_t", "rfs_e")[2]))
    for p in prim:  # their 5-year logistic endpoint
        k = p.cohort.iloc[0]
        q = p.dropna(subset=["age", "time"]).assign(died5=lambda d: ((d.event == 1) & (d.time <= 60)).astype(float))
        sec.append(dict(test=f"{k} died <= 5 y ~ top tertile + age (logistic)", **logistic(q, ["mag_top", "age"])))
    # correlations: magnitude vs G2M, MKI67, in-plane log-odds; seed concordance
    m_c, m_X = O.S.metabric()
    inpl = {"METABRIC": pd.Series(metabric_readout(m_X, m_c, O.L.GENES), index=m_c.index)}
    c_s, X_s = O.scanb()
    Ys = X_s[O.A.GENES].to_numpy(float) * LN2
    inpl["SCAN-B"] = pd.Series(sample_log_odds(Ys - Ys.mean(0), O.A.WF), index=X_s.index)
    for k, (p, _, r) in tabs.items():
        row = dict(test=f"Spearman magnitude, {k} LumA", rho_g2m=spearmanr(p.magnitude, p.g2m)[0],
                   rho_mki67=spearmanr(p.magnitude, p.MKI67)[0])
        if k in inpl:
            row["rho_inplane"] = spearmanr(p.magnitude, inpl[k].reindex(p.index), nan_policy="omit")[0]
        s = r.loc[r.index.isin(p.index), [c for c in r.columns if c.startswith("magnitude_s")]]
        row["seed_concordance_min_rho"] = float(s.corr(method="spearman").min().min())
        sec.append(row)
    sec = pd.DataFrame(sec)
    sec.to_csv(OUT / "a193_secondary.csv", index=False)
    pd.set_option("display.width", 250)
    print(tt.round(4).to_string(index=False))
    print(rep.round(4).to_string(index=False))
    print(sec.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
