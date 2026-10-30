"""A153 scoring: primaries, H153.1 (Holm effect), H153.2 (TOST equivalence), H153.3 (pooled family), readout tables.

Reads output/ of compute_a153.py and output/family_eligibility.csv. Scope `all5` is the registered eligibility (A153
note of 2026-10-30); `rowcheck3` is the descriptive sensitivity (fibroblast_stromal and endothelial both eligible).
Writes, per scope: readout_<scope>.csv (target x compartment, class and reason), summary_<scope>.csv (class per
disease), h1531_<scope>.csv (primaries with Holm), h1533_<scope>.csv (pooled families), a17_comparison_<scope>.csv
(changed verdicts first); and mde_table.csv, target_summary.csv.

    python score_a153.py [--subsample]
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[2]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(PAPER / "results/extended_data/scripts/host_identity"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import t as student  # noqa: E402

import clocksync as cs  # noqa: E402
from common.paths import PSEUDOBULK  # noqa: E402
from anchor_rule import family  # noqa: E402

OUT = HERE.parent / "output"
A17 = PAPER / "results/figure_4/output/resilience/a17_contrast.csv"
A159 = PAPER / "exploratory/A159_specific_coupling/output/d159_compartments.csv"
MDE_MAX, EQ, Z90 = 0.2, 0.15, 1.645
ATLAS = {"lung": "lung_hlca", "heart": "heart"}  # D159.1 healthy atlas of the same tissue (descriptive flag)
DISEASE = dict(lung_ipf_hlca="lung fibrosis", hlca_fibrotic="lung fibrosis", ild_natri="lung fibrosis", lupus="lupus",
               dcm="DCM", dcm_healthy="DCM", myocarditis="myocarditis", ckd_sc="CKD", skin_ad="atopic dermatitis",
               skin_ad_nonlesional="atopic dermatitis", skin_ssc="systemic sclerosis", islets_t2d="T2D islets",
               ibd_crohn="Crohn", ibd_crohn_inflamed="Crohn", ibd_crohn_uninflamed="Crohn", ia_crohn="Crohn",
               endometriosis="endometriosis", endometriosis_severe="endometriosis")


def load_mde() -> pd.DataFrame:
    rows = [json.loads(f.read_text()) for f in sorted((OUT / "mde").glob("*.json"))]
    m = pd.DataFrame(rows)
    if m.empty:
        return pd.DataFrame(columns=["target", "compartment", "mde_dstar"])
    m["mde_dstar"] = pd.to_numeric(m.mde_dstar, errors="coerce").where(m.mde_dstar != "inf", np.inf)
    return m.drop(columns=["curve"], errors="ignore")


def compartments_of(dataset: str) -> list:
    return sorted(pd.read_parquet(PSEUDOBULK / f"{dataset}.parquet", columns=["compartment"]).compartment.unique())


def classify(scope: str, elig: set, status: pd.DataFrame, tab: pd.DataFrame, mde: pd.DataFrame, brain: pd.DataFrame):
    rows = []
    for r in status.itertuples():
        base = dict(target=r.target, disease=DISEASE.get(r.target, r.target))
        if r.status != "ok":
            for k in compartments_of(r.dataset):
                rows.append(dict(base, compartment=k, family=family(k), cls="not estimable", reason=r.reason))
    t = tab.merge(mde[["target", "compartment", "mde_dstar", "c0", "s_star", "error"]], on=["target", "compartment"], how="left")
    for r in t.itertuples():
        base = dict(target=r.target, disease=DISEASE.get(r.target, r.target), compartment=r.compartment, family=r.family)
        if str(r.status) != "ok" or not np.isfinite(r.p_star) or not np.isfinite(r.se_d_star):
            reason = ("both arms at the boundary" if "boundary" in str(r.status) else
                      "too few donors or failed fit" if str(r.status) == "ok" else str(r.status))
            rows.append(dict(base, cls="not estimable", reason=reason))
        elif r.family not in elig:
            rows.append(dict(base, cls="not eligible", reason=f"family {r.family} not eligible ({scope})"))
        elif pd.isna(r.mde_dstar):
            rows.append(dict(base, cls="not eligible", reason="MDE not computed" if pd.isna(r.error) else f"MDE failed: {r.error}"))
        elif r.mde_dstar > MDE_MAX:
            rows.append(dict(base, cls="not eligible", reason="MDE > max (power 0.8 not reached)" if np.isinf(r.mde_dstar)
                             else f"MDE {r.mde_dstar:.3f} > {MDE_MAX}"))
        else:
            rows.append(dict(base, cls="primary", reason=f"MDE {r.mde_dstar:.3f} <= {MDE_MAX}"))
    for r in brain.itertuples():
        rows.append(dict(target=r.target, disease="brain: " + r.target.split("_")[0], compartment="(all)", family="", cls="not eligible",
                         reason=r.reason))
    ro = pd.DataFrame(rows)
    ro = ro.merge(tab.drop(columns=["family", "status"]), on=["target", "compartment"], how="left")
    ro = ro.merge(mde[["target", "compartment", "mde_dstar", "c0", "s_star"]], on=["target", "compartment"], how="left")
    prim = ro.cls == "primary"
    ro["p_holm"] = np.nan
    if prim.any():
        ro.loc[prim, "p_holm"] = cs.stats.holm(ro.loc[prim, "p_star"].to_numpy())
    ro["ci90_lo"], ro["ci90_hi"] = ro.d_star - Z90 * ro.se_d_star, ro.d_star + Z90 * ro.se_d_star
    eff = prim & (ro.p_holm < 0.05)
    eq = prim & ~eff & (ro.ci90_lo > -EQ) & (ro.ci90_hi < EQ)
    ro.loc[prim, "cls"] = "inconclusive"
    ro.loc[eq, "cls"] = "equivalent"
    ro.loc[eff, "cls"] = "effect"
    ro["reason"] = ro.reason.where(~prim, ro.reason + np.where(eff, "; Holm p* < 0.05", np.where(eq, "; 90% CI within +-0.15",
                                                                                                    "; no effect, 90% CI beyond +-0.15")))
    ro["specificity"] = np.where(eff, np.where(ro.p_companion_d_star < 0.05, "target-specific (companion p < 0.05)",
                                               "not target-specific (companion p >= 0.05)"), "")
    ro["pseudo_clock"] = "not available (no genome-wide table for the disease control arm)"
    ro["primary"] = prim
    # D159.1 flag of the healthy atlas compartment of the same name (descriptive)
    d159 = pd.read_csv(A159)
    flag = d159.set_index(["atlas", "compartment"]).clock_specific
    tis = ro.target.map(status.set_index("target").tissue)
    ro["d159_clock_specific"] = [flag.get((ATLAS.get(ti), k), np.nan) for ti, k in zip(tis, ro.compartment)]
    return ro


def h1533(ro: pd.DataFrame, fam: str) -> dict:
    p = ro[ro.primary & (ro.family == fam)]
    per = []
    for tg, g in p.groupby("target"):
        w = 1 / g.se_d_star ** 2
        per.append(dict(target=tg, n_comp=len(g), y=float((w * g.d_star).sum() / w.sum()), se=float(g.se_d_star.mean())))
    per = pd.DataFrame(per)
    out = dict(family=fam, n_targets=len(per), targets=";".join(per.target) if len(per) else "")
    if len(per) < 3:
        return dict(out, verdict="not estimable (< 3 targets)")
    m = cs.stats.meta(per.y, per.se, method="hk")
    k, mu, se = m["k"], m["random"], m["se_random"]
    h = student.ppf(0.95, k - 1) * se
    verdict = ("effect (95% interval excludes 0)" if (m["lo"] > 0 or m["hi"] < 0) else
               "equivalent (90% interval within +-0.15)" if (mu - h > -EQ and mu + h < EQ) else "inconclusive")
    return dict(out, mean=mu, se=se, lo95=m["lo"], hi95=m["hi"], lo90=mu - h, hi90=mu + h, tau2=m["tau2"], k=k, verdict=verdict)


def a17_compare(ro: pd.DataFrame) -> pd.DataFrame:
    a = pd.read_csv(A17)
    a = a[a.rule == "learned"][["target", "compartment", "status", "p", "p_star", "d_star", "primary", "p_holm_primary"]]
    def verdict(r):
        if r.status != "ok" or pd.isna(r.p):
            return "not estimable"
        if bool(r.primary) and r.p_holm_primary < 0.05:
            return "effect (A17 primary, Holm)"
        if bool(r.primary):
            return "primary, no effect"
        return "non-primary, nominal p < 0.05" if r.p < 0.05 else "non-primary, n.s."
    a["a17_verdict"] = a.apply(verdict, axis=1)
    a = a.rename(columns=dict(status="a17_status", p="a17_p", p_star="a17_p_star", d_star="a17_d_star", primary="a17_primary",
                              p_holm_primary="a17_p_holm_primary"))
    c = ro[ro.compartment != "(all)"].merge(a, on=["target", "compartment"], how="left")
    cat = lambda v: "not estimable" if v == "not estimable" else "effect" if str(v).startswith("effect") else "estimable"
    c["a17_cat"] = c.a17_verdict.map(lambda v: cat(v) if isinstance(v, str) else "absent in A17")
    c["a153_cat"] = c.cls.map(lambda v: "not estimable" if v == "not estimable" else "effect" if v == "effect" else "estimable")
    c["changed"] = (c.a17_cat != c.a153_cat) & (c.a17_cat != "absent in A17")
    c["change"] = np.where(~c.changed, "", c.a17_cat + " -> " + c.a153_cat)
    cols = ["target", "disease", "compartment", "family", "a17_verdict", "a17_p", "a17_d_star", "cls", "d_star", "p_star", "changed", "change"]
    return c.sort_values(["changed", "target", "compartment"], ascending=[False, True, True])[cols]


def main():
    global OUT
    if "--subsample" in sys.argv:
        OUT = OUT / "subsample"
    status = pd.read_csv(OUT / "target_status.csv").fillna({"reason": ""})
    brain = pd.read_csv(OUT / "brain_targets.csv")
    tab = pd.concat([pd.read_csv(f) for f in sorted((OUT / "targets").glob("*_table.csv"))], ignore_index=True)
    mde = load_mde()
    mde.to_csv(OUT / "mde_table.csv", index=False)
    status.to_csv(OUT / "target_summary.csv", index=False)
    fe = pd.read_csv(HERE.parent / "output" / "family_eligibility.csv")
    for scope in ("all5", "rowcheck3"):
        elig = set(fe[(fe.scope == scope) & fe.eligible].family)
        ro = classify(scope, elig, status, tab, mde, brain)
        ro.to_csv(OUT / f"readout_{scope}.csv", index=False)
        summ = ro.groupby(["disease", "cls"]).size().unstack(fill_value=0)
        summ["targets"] = ro.groupby("disease").target.nunique()
        summ.to_csv(OUT / f"summary_{scope}.csv")
        prim = ro[ro.primary].sort_values("p_star")
        prim[["target", "disease", "compartment", "family", "d_star", "se_d_star", "p_star", "p_holm", "ci90_lo", "ci90_hi", "cls",
              "mde_dstar", "p_companion_d_star", "d_star_minus_companion", "specificity", "d159_clock_specific"]] \
            .to_csv(OUT / f"h1531_{scope}.csv", index=False)
        pd.DataFrame([h1533(ro, f) for f in ("fibroblast_stromal", "endothelial")]).to_csv(OUT / f"h1533_{scope}.csv", index=False)
        a17_compare(ro).to_csv(OUT / f"a17_comparison_{scope}.csv", index=False)
        print(f"== {scope}: eligible families {sorted(elig)}")
        print(ro.cls.value_counts().to_string())
        print(summ.to_string())
        print(pd.read_csv(OUT / f"h1533_{scope}.csv").to_string())
        print(prim[["target", "compartment", "d_star", "se_d_star", "p_star", "p_holm", "cls", "mde_dstar"]].round(3).to_string())


if __name__ == "__main__":
    main()
