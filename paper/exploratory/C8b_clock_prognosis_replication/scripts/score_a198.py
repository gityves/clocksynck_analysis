"""A198 scoring (registered in `paper/docs/registered_experiments.md`, commit 7d4fe8c): does the tumour clock signal modify the
benefit of chemotherapy or endocrine therapy in breast cancer, beyond proliferation?

Readout: per-tumour in-plane log-odds s (A191), genes centred within METABRIC cohort or SCAN-B platform, s standardised within
cohort. Cox OS, METABRIC + SCAN-B, strata = cohort; terms s x T, s, T, G2M x T, G2M and cohort-specific clinical covariates;
statistic = Wald z of s x T. H198.1 T = chemotherapy (all tumours); H198.2 T = endocrine therapy (ER+). Holm over two.
Power step first (MDE of the interaction HR ratio per SD; above 1.5: not estimable). Specificity: 200 expression-matched
pseudo-clock sets. Only summary statistics are written.

    python score_a198.py [--power-only] [--npseudo 200]
"""
import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import statsmodels.api as sm  # noqa: E402
from scipy.stats import norm  # noqa: E402
from statsmodels.duration.hazard_regression import PHReg  # noqa: E402

import outcomes as O  # noqa: E402
from score_a191 import SCCC, Z80, holm  # noqa: E402
from clocksync.decompose import sample_log_odds  # noqa: E402

OUT = HERE.parents[1] / "output"
LN2, LIMIT, SEED = np.log(2), 1.5, 20261103
CLIN = {"METABRIC": ["age", "log_size", "nodes", "grade3", "grade1", "erbb2"],
        "SCAN-B": ["age", "log_size", "node_pos", "grade3", "grade1", "her2"]}
TESTS = {"H198.1": ("chemo", False), "H198.2": ("endo", True)}  # (treatment, ER+ only)


def log_odds(X, genes, groups):
    Y = X[genes].to_numpy(float) * LN2
    Y = Y - pd.DataFrame(Y).groupby(np.asarray(groups)).transform("mean").to_numpy()
    s = sample_log_odds(Y, O.L.wf_for(genes))
    return (s - s.mean()) / s.std()


def cohorts():
    """{cohort: (table, expression, centring groups)} with chemo, endo, subtype, s."""
    c, X = O.metabric()
    m, _ = O.S.metabric()
    pt = pd.read_csv(SCCC / "metabric" / "brca_metabric" / "data_clinical_patient.txt", sep="\t", comment="#").set_index("PATIENT_ID")
    yn = {"YES": 1.0, "NO": 0.0}
    c["chemo"] = pt.CHEMOTHERAPY.map(yn).reindex(m.patient).to_numpy()
    c["endo"] = m.HORMONE_THERAPY.map(yn)
    c["subtype"] = m.claudin
    mb = (c, X, m.cohort.to_numpy())
    c, X = O.scanb()
    a, _ = O.A.load()
    c["chemo"], c["endo"] = a["chemo treated"].astype(float), a["endocrine treated"].astype(float)
    c["subtype"] = c.pam50
    sb = (c, X, a.platform.to_numpy())
    out = {}
    for k, (c, X, grp) in (("METABRIC", mb), ("SCAN-B", sb)):
        c = c.assign(cohort=k, s=log_odds(X, O.L.GENES, grp))
        c["g2m"] = (c.g2m - c.g2m.mean()) / c.g2m.std()
        for lev in sorted(c.subtype.dropna().unique()):
            if lev != "LumA":
                c[f"sub_{lev}"] = (c.subtype == lev).astype(float).where(c.subtype.notna())
        out[k] = (c, X, grp)
    return out


def frame(data, treat, er_only, x="s", extra_cov=(), tcol="time", ecol="event"):
    """Complete-case pooled table (time / event taken from tcol / ecol) and its design; returns (table, x x T, other terms)."""
    parts, cov = [], {}
    for k, (c, _, _) in data.items():
        cv = CLIN[k] + [v for v in c.columns if v.startswith("sub_")] + ([] if er_only else ["er"]) + list(extra_cov)
        d = c[c.er == 1] if er_only else c
        d = d.dropna(subset=[tcol, ecol, x, "g2m", treat] + cv).assign(time=lambda q: q[tcol], event=lambda q: q[ecol])
        parts.append(d); cov[k] = cv
    d = pd.concat(parts)
    T = d[treat].to_numpy(float)
    xs, g = d[x].to_numpy(float), d.g2m.to_numpy(float)
    Z = np.column_stack([xs, T, g * T, g, O.block_design(d, cov)])
    return d, xs * T, Z


def fit(d, inter, Z):
    r = PHReg(d.time.to_numpy(float), np.column_stack([inter, Z]), status=d.event.to_numpy(float),
              strata=d.cohort.to_numpy(), ties="efron").fit(disp=0)
    b, se = r.params[0], r.bse[0]
    return dict(n=len(d), events=int(d.event.sum()), z=b / se, ratio=np.exp(b), ratio_lo=np.exp(b - 1.96 * se),
                ratio_hi=np.exp(b + 1.96 * se), p=float(2 * norm.sf(abs(b / se))), hr_s_untreated=float(np.exp(r.params[1])))


def mde(d, inter, Z):
    C = np.column_stack([Z, (d.cohort == "SCAN-B").to_numpy(float)])
    res = inter - sm.OLS(inter, sm.add_constant(C)).fit().fittedvalues
    return float(np.exp(Z80 / np.sqrt(d.event.sum() * res.var())))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--power-only", action="store_true"); ap.add_argument("--npseudo", type=int, default=200)
    a = ap.parse_args()
    data = cohorts()
    power = []
    for h, (t, er) in TESTS.items():
        d, i, Z = frame(data, t, er)
        md = mde(d, i, Z)
        power.append(dict(test=h, treatment=t, er_only=er, n=len(d), events=int(d.event.sum()),
                          treated=int(d[t].sum()), mde_ratio=md, limit=LIMIT, estimable=md <= LIMIT))
    power = pd.DataFrame(power)
    power.to_csv(OUT / "a198_power.csv", index=False)
    print(power.round(4).to_string(index=False), flush=True)
    if a.power_only:
        return

    # pseudo-clock readouts, one matched set per cohort per draw
    pseudo = {}
    for j, (k, (c, X, grp)) in enumerate(data.items()):
        pool = O.L.pool_genes(X)
        sets = O.L.pseudo_sets(X[pool + O.L.GENES], pool, a.npseudo, np.random.default_rng([SEED, j]))
        pseudo[k] = [log_odds(X, s, grp) for s in sets]
    rows = []
    for h, (t, er) in TESTS.items():
        if not power.set_index("test").estimable[h]:
            rows.append(dict(test=h, treatment=t, p=np.nan)); continue
        d, i, Z = frame(data, t, er)
        r = fit(d, i, Z)
        zp = []
        for b in range(a.npseudo):
            for k, (c, _, _) in data.items():
                c["s_pseudo"] = pseudo[k][b]
            dp, ip, Zp = frame(data, t, er, "s_pseudo")
            zp.append(fit(dp, ip, Zp)["z"])
        zp = np.abs(zp)
        rows.append(dict(test=h, treatment=t, **r, pseudo_abs_z_q95=float(np.quantile(zp, .95)),
                         pseudo_percentile=float(np.mean(zp <= abs(r["z"])) * 100), clock_specific=abs(r["z"]) > np.quantile(zp, .95)))
    tt = pd.DataFrame(rows)
    ok = tt.p.notna()
    tt.loc[ok, "p_holm"] = holm(tt.loc[ok, "p"])
    tt["verdict"] = np.select([~ok, (tt.p_holm < .05) & tt.clock_specific.fillna(False).astype(bool)], ["not estimable", "pass"], "fail")
    tt.to_csv(OUT / "a198_tests.csv", index=False)

    # secondary
    sec = []
    for h, (t, er) in TESTS.items():
        d, i, Z = frame(data, t, er)
        sec.append(dict(test=f"{h} without G2M x T", **fit(d, i, np.delete(Z, 2, axis=1))))
        for k in data:
            sub = (d.cohort == k).to_numpy()
            sec.append(dict(test=f"{h} {k} only", **fit(d[sub], i[sub], Z[sub][:, Z[sub].std(0) > 0])))
        dr, ir, Zr = frame({"METABRIC": data["METABRIC"]}, t, er, tcol="rfs_t", ecol="rfs_e")
        sec.append(dict(test=f"{h} METABRIC RFS", **fit(dr, ir, Zr[:, Zr.std(0) > 0])))
        for k in data:  # MKI67 in place of G2M
            data[k][0]["g2m_saved"] = data[k][0].g2m
            data[k][0]["g2m"] = data[k][0].MKI67
        dm_, im_, Zm_ = frame(data, t, er)
        sec.append(dict(test=f"{h} MKI67 in place of G2M", **fit(dm_, im_, Zm_)))
        for k in data:
            data[k][0]["g2m"] = data[k][0].g2m_saved
    d, i, Z = frame(data, "chemo", True)
    sec.append(dict(test="H198.1 within ER+", **fit(d, i, Z)))
    sec = pd.DataFrame(sec)
    sec.to_csv(OUT / "a198_secondary.csv", index=False)
    pd.set_option("display.width", 250)
    print(tt.round(4).to_string(index=False)); print(sec.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
