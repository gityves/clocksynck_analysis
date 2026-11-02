"""A191 (C8) scoring (registered in `paper/docs/registered_experiments.md`, commit ff1c1b4): tumour clock state and outcome.

H191.1 (L1): METABRIC, ER+ (IHC) with hormone therapy, relapse-free survival: per-tumour in-plane log-odds
(`clocksync.decompose.sample_log_odds`, 14 genes, log2 x ln 2, gene-centred within METABRIC cohort), Cox (PHReg, strata =
cohort) adjusted for age, size, grade, log(1 + positive nodes), PAM50 and the G2M score.
H191.2 (L3): GSE14520 HCC, recurrence: per-patient folded tumour-to-own-liver phase distance relative to the cohort offset
(`paired.cohort_states`, arrays, 13 genes), Cox adjusted for age, sex, BCLC (B/C), cirrhosis, AFP high, size large,
multinodular.
Power step first (MDE of the HR per SD from events and 1 - R^2 of the readout on the covariates; > 1.5: not run). Specificity:
200 expression-matched pseudo-clock sets (A166/A127 construction), clock-specific if |z| > q95 |z_pseudo|. Holm over two.
Secondary: OS, all METABRIC, unadjusted HRs, GSE14520 internal time at surgery (liver phase, cos/sin, LR test 2 df).
Only summary statistics are written.

    python score_a191.py [--n-pseudo 200] [--check]
"""
import argparse
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve()
ROOT = HERE.parents[4]
EXP = ROOT / "paper" / "exploratory"
sys.path[:0] = [str(EXP / "A166_metabric" / "scripts"), str(EXP / "A177_c1_replication" / "scripts")]
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import statsmodels.api as sm  # noqa: E402
from scipy.stats import chi2, norm  # noqa: E402
from statsmodels.duration.hazard_regression import PHReg  # noqa: E402

import a166_lib as L  # noqa: E402  (GENES, WF, wf_for, pseudo_sets, pool_genes, load, sh)
import power_a177 as Q  # noqa: E402  (load: GSE14520 pairs)
from clocksync.decompose import sample_log_odds  # noqa: E402
from clocksync.paired import cohort_states  # noqa: E402
from clocksync.stats import complex_state  # noqa: E402
from clocksync._core.controls import CLOCK_LIKE  # noqa: E402
from clocksync.genes import canonical  # noqa: E402
from clocksync.waveform import Waveform  # noqa: E402

OUT = HERE.parents[1] / "output"
SCCC = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian")
SEED, LN2, HOURS, Z80 = 191, np.log(2), 12 / np.pi, 1.959964 + 0.841621
T0 = time.time()


def log(*a):
    print(f"[{(time.time() - T0) / 60:5.1f} min]", *a, flush=True)


def holm(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    adj = np.minimum(1, np.maximum.accumulate((len(p) - np.arange(len(p))) * p[o]))
    out = np.empty_like(adj)
    out[o] = adj
    return out


def design(df, cols, cat):
    """Covariate matrix: numeric columns as is, categorical ones as dummies (first level dropped)."""
    X = df[cols].copy()
    for c in cat:
        X = pd.concat([X.drop(columns=c), pd.get_dummies(X[c].astype(str), prefix=c, drop_first=True, dtype=float)], axis=1)
    return X.astype(float)


def cox(time_, event, x, Z, strata=None):
    """Wald z, HR per SD and its interval for the standardised readout x, adjusted for Z; log-likelihood."""
    x = (x - x.mean()) / x.std()
    X = np.column_stack([x, Z]) if Z is not None and Z.shape[1] else x[:, None]
    r = PHReg(time_, X, status=event, strata=strata, ties="efron").fit(disp=0)
    b, se = r.params[0], r.bse[0]
    return dict(z=b / se, hr=np.exp(b), hr_lo=np.exp(b - 1.96 * se), hr_hi=np.exp(b + 1.96 * se),
                p=float(2 * norm.sf(abs(b / se))), llf=r.llf)


def mde(x, Z, events):
    """Schoenfeld MDE of the HR per SD with 1 - R^2 of x on the covariates."""
    x = (x - x.mean()) / x.std()
    r2 = sm.OLS(x, sm.add_constant(Z)).fit().rsquared if Z is not None and Z.shape[1] else 0.0
    return float(np.exp(Z80 / np.sqrt(events * (1 - r2)))), float(r2)


# ---- H191.1 METABRIC ------------------------------------------------------------------------------------------------------

def metabric():
    m, X = L.load("metabric")
    assert X.index.equals(m.index)
    d = SCCC / "metabric" / "brca_metabric"
    pt = pd.read_csv(d / "data_clinical_patient.txt", sep="\t", comment="#").set_index("PATIENT_ID")
    smp = pd.read_csv(d / "data_clinical_sample.txt", sep="\t", comment="#").set_index("SAMPLE_ID")
    c = m.join(pt[["AGE_AT_DIAGNOSIS", "LYMPH_NODES_EXAMINED_POSITIVE", "HORMONE_THERAPY", "RFS_MONTHS", "RFS_STATUS",
                   "OS_MONTHS", "OS_STATUS"]], on="patient").join(smp[["TUMOR_SIZE", "GRADE"]])
    c["nodes"] = np.log1p(c.LYMPH_NODES_EXAMINED_POSITIVE)
    c["rfs_t"], c["rfs_e"] = c.RFS_MONTHS, c.RFS_STATUS.astype(str).str.startswith("1").astype(float)
    c["os_t"], c["os_e"] = c.OS_MONTHS, c.OS_STATUS.astype(str).str.startswith("1").astype(float)
    c["subtype"] = c.claudin
    return c, X


def metabric_readout(X, c, genes):
    Y = X[genes].to_numpy(float) * LN2
    Y = Y - pd.DataFrame(Y).groupby(c.cohort.to_numpy()).transform("mean").to_numpy()
    return sample_log_odds(Y, L.wf_for(genes))


# ---- H191.2 GSE14520 ------------------------------------------------------------------------------------------------------

def gse14520():
    t, n, p, wf = Q.load("gse14520")
    s = pd.read_csv(SCCC / "scout_luad_lihc" / "GSE14520" / "GSE14520_Extra_Supplement.txt.gz", sep="\t")
    s = s[s["Tissue Type"] == "Tumor"].drop_duplicates("ID").set_index("ID")
    c = p[["patient"]].join(s, on="patient")
    c["bclc_bc"] = c["BCLC staging"].isin(["B", "C"]).astype(float).where(c["BCLC staging"].isin(["0", "A", "B", "C"]))
    c["cirrhosis"] = c.Cirrhosis.map({"Y": 1.0, "N": 0.0})
    c["afp_high"] = c["AFP (>/<=300ng/ml)"].map({"high": 1.0, "low": 0.0})
    c["size_large"] = c["Main Tumor Size (>/<=5 cm)"].map({"large": 1.0, "small": 0.0})
    c["multinod"] = c.Multinodular.map({"Y": 1.0, "N": 0.0})
    c["male"] = c.Gender.map({"M": 1.0, "F": 0.0})
    c["age"] = pd.to_numeric(c.Age, errors="coerce")
    c["rec_t"], c["rec_e"] = pd.to_numeric(c["Recurr months"], errors="coerce"), pd.to_numeric(c["Recurr status"], errors="coerce")
    c["os_t"], c["os_e"] = pd.to_numeric(c["Survival months"], errors="coerce"), pd.to_numeric(c["Survival status"], errors="coerce")
    return t, n, c.reset_index(drop=True), wf


def phase_readouts(t, n, wf):
    """Folded tumour-to-liver distance (h) relative to the cohort offset, and the liver phase relative to the cohort (rad)."""
    g = list(wf.genes)
    bT, _, bN, _ = cohort_states(t[g], n[g], counts=False, waveform=wf)
    zT, zN = complex_state(bT), complex_state(bN)
    cc = np.sum(zT * np.conj(zN))
    d = np.abs(np.angle(zT * np.conj(zN) * np.conj(cc) / abs(cc))) * HOURS
    return d, np.angle(zN)


def pseudo_sets_generic(X, genes, pool, n, rng):
    """A166 `pseudo_sets` for an arbitrary clock gene list (mean and SD matching, 20 nearest pool genes)."""
    pm, ps = X[pool].mean(), X[pool].std()
    zs = lambda a, ref: (a - ref.mean()) / ref.std()  # noqa: E731
    P = np.column_stack([zs(pm.values, pm), zs(ps.values, ps)])
    target = np.column_stack([zs(X[genes].mean().values, pm), zs(X[genes].std().values, ps)])
    sets = []
    for _ in range(n):
        used, s = set(), []
        for k in rng.permutation(len(genes)):
            near = [j for j in np.argsort(((P - target[k]) ** 2).sum(1))[:40] if j not in used][:20]
            j = int(rng.choice(near))
            used.add(j)
            s.append((k, pool[j]))
        sets.append([g for _, g in sorted(s)])
    return sets


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--n-pseudo", type=int, default=200)
    ap.add_argument("--check", action="store_true", help="5 pseudo sets, outputs suffixed _check")
    ap.add_argument("--power-only", action="store_true", help="stop after the power step (no outcome association)")
    a = ap.parse_args()
    n_ps = 5 if a.check else a.n_pseudo
    sfx = "_check" if a.check else ""
    OUT.mkdir(exist_ok=True)
    power, tests, sec = [], [], []

    # ---- METABRIC ----
    c, X = metabric()
    cov1, cat1 = ["AGE_AT_DIAGNOSIS", "TUMOR_SIZE", "GRADE", "nodes", "subtype", "prolif"], ["GRADE", "subtype"]
    c["s"] = metabric_readout(X, c, L.GENES)
    prim = c[(c.er == 1) & (c.HORMONE_THERAPY == "YES")].dropna(subset=cov1 + ["rfs_t", "rfs_e", "s"])
    Z1 = design(prim, cov1, cat1).to_numpy()
    md, r2 = mde(prim.s.to_numpy(), np.column_stack([Z1, pd.get_dummies(prim.cohort, drop_first=True, dtype=float)]),
                 prim.rfs_e.sum())
    power.append(dict(test="H191.1", n=len(prim), events=int(prim.rfs_e.sum()), r2_on_covariates=r2, mde_hr_per_sd=md,
                      estimable=md <= 1.5))
    log("H191.1 power", power[-1])

    # ---- GSE14520 ----
    t, n, g, wf = gse14520()
    g["d"], g["phi"] = phase_readouts(t, n, wf)
    cov2 = ["age", "male", "bclc_bc", "cirrhosis", "afp_high", "size_large", "multinod"]
    prim2 = g.dropna(subset=cov2 + ["rec_t", "rec_e"])
    Z2 = prim2[cov2].to_numpy(float)
    md2, r22 = mde(prim2.d.to_numpy(), Z2, prim2.rec_e.sum())
    power.append(dict(test="H191.2", n=len(prim2), events=int(prim2.rec_e.sum()), r2_on_covariates=r22, mde_hr_per_sd=md2,
                      estimable=md2 <= 1.5))
    log("H191.2 power", power[-1])
    pd.DataFrame(power).to_csv(OUT / f"power{sfx}.csv", index=False)
    if a.power_only:
        return

    # ---- primary tests (only the estimable ones) ----
    rng = np.random.default_rng(SEED)
    if power[0]["estimable"]:
        r = cox(prim.rfs_t.to_numpy(float), prim.rfs_e.to_numpy(float), prim.s.to_numpy(), Z1, strata=prim.cohort.to_numpy())
        pool = L.pool_genes(X)
        zs = []
        for gs in L.pseudo_sets(X, pool, n_ps, rng):
            sp = pd.Series(metabric_readout(X, c, gs), index=c.index).loc[prim.index].to_numpy()
            zs.append(cox(prim.rfs_t.to_numpy(float), prim.rfs_e.to_numpy(float), sp, Z1, strata=prim.cohort.to_numpy())["z"])
        zs = np.abs(zs)
        tests.append(dict(test="H191.1", **{k: v for k, v in r.items() if k != "llf"}, pseudo_q95_abs_z=np.quantile(zs, .95),
                          pseudo_pct=100 * np.mean(zs < abs(r["z"]))))
        log("H191.1", tests[-1])
    if power[1]["estimable"]:
        r = cox(prim2.rec_t.to_numpy(float), prim2.rec_e.to_numpy(float), prim2.d.to_numpy(), Z2)
        allg = pd.concat([t, n])
        clock = set(CLOCK_LIKE) | set(L.sh.hallmark_g2m()) | set(L.sh.PURITY)
        pool2 = [x for x in allg.columns if canonical(x) not in clock and x not in wf.genes]
        zs = []
        for gs in pseudo_sets_generic(allg, list(wf.genes), pool2, n_ps, rng):
            dp, _ = phase_readouts(t, n, Waveform(tuple(gs), wf.amp, wf.psi, None, None, "set"))
            dp = pd.Series(dp).loc[prim2.index].to_numpy()
            zs.append(cox(prim2.rec_t.to_numpy(float), prim2.rec_e.to_numpy(float), dp, Z2)["z"])
        zs = np.abs(zs)
        tests.append(dict(test="H191.2", **{k: v for k, v in r.items() if k != "llf"}, pseudo_q95_abs_z=np.quantile(zs, .95),
                          pseudo_pct=100 * np.mean(zs < abs(r["z"]))))
        log("H191.2", tests[-1])
    tt = pd.DataFrame(tests)
    if len(tt):
        tt["p_holm"] = holm(tt.p)
        tt["clock_specific"] = tt.z.abs() > tt.pseudo_q95_abs_z
        tt["passed"] = (tt.p_holm < .05) & tt.clock_specific
    tt.to_csv(OUT / f"tests{sfx}.csv", index=False)

    # ---- secondary (not counted) ----
    os1 = c[(c.er == 1) & (c.HORMONE_THERAPY == "YES")].dropna(subset=cov1 + ["os_t", "os_e", "s"])
    sec.append(dict(test="H191.1 OS", n=len(os1), events=int(os1.os_e.sum()),
                    **cox(os1.os_t.to_numpy(float), os1.os_e.to_numpy(float), os1.s.to_numpy(), design(os1, cov1, cat1).to_numpy(),
                          strata=os1.cohort.to_numpy())))
    allm = c.dropna(subset=cov1 + ["er", "rfs_t", "rfs_e", "s"])
    sec.append(dict(test="H191.1 all METABRIC RFS (+ER)", n=len(allm), events=int(allm.rfs_e.sum()),
                    **cox(allm.rfs_t.to_numpy(float), allm.rfs_e.to_numpy(float), allm.s.to_numpy(),
                          design(allm, cov1 + ["er"], cat1).to_numpy(), strata=allm.cohort.to_numpy())))
    sec.append(dict(test="H191.1 unadjusted RFS", n=len(prim), events=int(prim.rfs_e.sum()),
                    **cox(prim.rfs_t.to_numpy(float), prim.rfs_e.to_numpy(float), prim.s.to_numpy(), None,
                          strata=prim.cohort.to_numpy())))
    o2 = g.dropna(subset=cov2 + ["os_t", "os_e"])
    sec.append(dict(test="H191.2 OS", n=len(o2), events=int(o2.os_e.sum()),
                    **cox(o2.os_t.to_numpy(float), o2.os_e.to_numpy(float), o2.d.to_numpy(), o2[cov2].to_numpy(float))))
    sec.append(dict(test="H191.2 unadjusted recurrence", n=len(prim2), events=int(prim2.rec_e.sum()),
                    **cox(prim2.rec_t.to_numpy(float), prim2.rec_e.to_numpy(float), prim2.d.to_numpy(), None)))
    # internal time at surgery: liver phase (cos, sin) added to the covariate model, LR test 2 df
    for lab, tc, ec in (("recurrence", "rec_t", "rec_e"), ("OS", "os_t", "os_e")):
        q = g.dropna(subset=cov2 + [tc, ec])
        base = PHReg(q[tc].to_numpy(float), q[cov2].to_numpy(float), status=q[ec].to_numpy(float), ties="efron").fit(disp=0)
        full = PHReg(q[tc].to_numpy(float), np.column_stack([q[cov2].to_numpy(float), np.cos(q.phi), np.sin(q.phi)]),
                     status=q[ec].to_numpy(float), ties="efron").fit(disp=0)
        lr = 2 * (full.llf - base.llf)
        sec.append(dict(test=f"GSE14520 internal time at surgery, {lab} (cos/sin, LR 2 df)", n=len(q), events=int(q[ec].sum()),
                        lr=lr, p=float(chi2.sf(lr, 2))))
    pd.DataFrame(sec).drop(columns="llf", errors="ignore").to_csv(OUT / f"secondary{sfx}.csv", index=False)
    pd.set_option("display.width", 250)
    print(pd.DataFrame(power).round(4).to_string(index=False))
    print(tt.round(4).to_string(index=False))
    print(pd.DataFrame(sec).drop(columns="llf", errors="ignore").round(4).to_string(index=False))
    log("done")


if __name__ == "__main__":
    main()
