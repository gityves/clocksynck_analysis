"""A194 scoring (registered in `paper/docs/registered_experiments.md`, commit a6d5922, addendum 1dfaa8e): is TimeTeller's
clock-dysfunction score Theta an independent prognostic factor in breast cancer, and does it survive a proliferation term?

Readout: `output/timeteller_theta.csv`, variant "aligned" (addendum). SCAN-B + TCGA BRCA, overall survival, Cox (PHReg,
Efron) stratified by cohort with cohort-specific clinical covariates (`outcomes.COVARIATES`), Theta standardised within
cohort. H194.1: clinical model; H194.2: + G2M. Holm over two. Power step first. Only summary statistics are written.

    python score_a194.py [--power-only]
"""
import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

import outcomes as O  # noqa: E402
from score_a191 import cox, holm, mde, metabric_readout  # noqa: E402
from clocksync.decompose import sample_log_odds  # noqa: E402

OUT = HERE.parents[1] / "output"
LN2 = np.log(2)


def readout(cohort):
    t = pd.read_csv(OUT / "timeteller_theta.csv")
    t = t[(t.cohort == cohort) & (t.variant == "aligned")].set_index("sample").theta
    return t


def pooled(cohorts):
    parts = []
    for k in cohorts:
        c, _ = O.LOADERS[k]()
        c = c.join(readout(k).rename("theta"), how="inner").assign(cohort=k)
        c["theta_z"] = (c.theta - c.theta.mean()) / c.theta.std()
        parts.append(c)
    return parts


def run(d, cov, extra=(), tcol="time", ecol="event", x="theta_z"):
    need = [tcol, ecol, x] + sorted({v for k in cov for v in cov[k]}) + list(extra)
    d = pd.concat([p.dropna(subset=[tcol, ecol, x] + cov.get(p.cohort.iloc[0], []) + list(extra)) for p in d])
    del need
    Z = O.block_design(d, cov, extra)
    r = cox(d[tcol].to_numpy(float), d[ecol].to_numpy(float), d[x].to_numpy(float), Z, strata=d.cohort.to_numpy())
    return d, Z, dict(n=len(d), events=int(d[ecol].sum()), **{k: v for k, v in r.items() if k != "llf"})


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--power-only", action="store_true")
    a = ap.parse_args()
    prim = pooled(["SCAN-B", "TCGA-BRCA"])
    cov = {k: O.COVARIATES[k] for k in ("SCAN-B", "TCGA-BRCA")}
    power = []
    for test, extra in (("H194.1", ()), ("H194.2", ("g2m",))):
        d = pd.concat([p.dropna(subset=["time", "event", "theta_z"] + cov[p.cohort.iloc[0]] + list(extra)) for p in prim])
        Z = O.block_design(d, cov, extra)
        md, r2 = mde(d.theta_z.to_numpy(), np.column_stack([Z, (d.cohort == "SCAN-B").to_numpy(float)]), d.event.sum())
        power.append(dict(test=test, n=len(d), events=int(d.event.sum()), r2_on_covariates=r2, mde_hr_per_sd=md, estimable=md <= 1.5))
    power = pd.DataFrame(power)
    power.to_csv(OUT / "a194_power.csv", index=False)
    print(power.round(4).to_string(index=False), flush=True)
    if a.power_only:
        return

    tests = []
    for test, extra in (("H194.1 clinical", ()), ("H194.2 clinical + G2M", ("g2m",))):
        _, _, r = run(prim, cov, extra)
        tests.append(dict(test=test, **r))
    tt = pd.DataFrame(tests)
    tt["p_holm"] = holm(tt.p)
    tt["passed_published_direction"] = (tt.p_holm < .05) & (tt.hr > 1)
    tt.to_csv(OUT / "a194_tests.csv", index=False)

    sec = []
    for extra, lab in (((), "unadjusted"), (("MKI67",), "clinical + MKI67")):
        cv = cov if lab != "unadjusted" else {k: [] for k in cov}
        sec.append(dict(test=f"SCAN-B + TCGA, {lab}", **run(prim, cv, extra)[2]))
    for p in prim:
        k = p.cohort.iloc[0]
        sec.append(dict(test=f"{k} clinical", **run([p], {k: cov[k]})[2]))
        sec.append(dict(test=f"{k} clinical + G2M", **run([p], {k: cov[k]}, ("g2m",))[2]))
    g12 = [p[p.grade12 == 1] for p in prim if "grade12" in p and p.grade12.notna().any()]
    sec.append(dict(test="SCAN-B grade 1-2, clinical (preprint subgroup)", **run(g12, {"SCAN-B": cov["SCAN-B"]})[2]))
    for p in prim:  # tertiles: top vs bottom third (worse vs better clock)
        k = p.cohort.iloc[0]
        q = pd.qcut(p.theta.rank(method="first"), 3, labels=False)
        p["theta_top"], p["theta_mid"] = (q == 2).astype(float), (q == 1).astype(float)
    sec.append(dict(test="SCAN-B + TCGA, top vs bottom Theta tertile, clinical (mid tertile in model)",
                    **run(prim, cov, ("theta_mid",), x="theta_top")[2]))
    mb = pooled(["METABRIC"])[0]
    mcov = {"METABRIC": O.COVARIATES["METABRIC"]}
    for lab, extra, tc, ec in (("OS clinical", (), "time", "event"), ("OS clinical + G2M", ("g2m",), "time", "event"),
                               ("RFS clinical", (), "rfs_t", "rfs_e"), ("RFS clinical + G2M", ("g2m",), "rfs_t", "rfs_e")):
        sec.append(dict(test=f"METABRIC (arrays, secondary) {lab}", **run([mb.assign(cohort="METABRIC")], mcov, extra, tc, ec)[2]))
    # correlations with proliferation and with the A191 in-plane log-odds
    c_s, X_s = O.scanb()
    sp = pd.Series(sample_log_odds((X_s[O.A.GENES].to_numpy(float) * LN2) - (X_s[O.A.GENES].to_numpy(float) * LN2).mean(0),
                                   O.A.WF), index=X_s.index)
    m_c, m_X = O.S.metabric()
    mp = pd.Series(metabric_readout(m_X, m_c, O.L.GENES), index=m_c.index)
    for p, inplane in ((prim[0], sp), (prim[1], None), (mb, mp)):
        k = p.cohort.iloc[0]
        row = dict(test=f"Spearman Theta vs G2M, {k}", rho=spearmanr(p.theta, p.g2m, nan_policy="omit")[0])
        sec.append(row)
        if inplane is not None:
            sec.append(dict(test=f"Spearman Theta vs in-plane log-odds, {k}",
                            rho=spearmanr(p.theta, inplane.reindex(p.index), nan_policy="omit")[0]))
    sec = pd.DataFrame(sec)
    sec.to_csv(OUT / "a194_secondary.csv", index=False)
    pd.set_option("display.width", 250)
    print(tt.round(4).to_string(index=False))
    print(sec.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
