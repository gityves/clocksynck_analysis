"""A165 scoring (registered in paper/docs/registered_experiments.md, A165): H165.1 per-compartment phase, amplitude and
noise effects (Holm over all primary rows x readouts), H165.2 pooled directional damping and noise per family (HK
meta-analysis, Holm over 4), approximate MDEs from the bootstrap SE, and the descriptive tables.

    python score_a165.py [--out DIR]
Writes (DIR): readout.csv (one row per target x compartment x readout), h1651.csv (primary rows), h1652.csv,
gate_share.csv (gate clock share per arm), summary.csv.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[2]
sys.path.insert(0, str(PAPER / "results/extended_data/scripts/host_identity"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402
from statsmodels.stats.multitest import multipletests  # noqa: E402

import clocksync as cs  # noqa: E402
from anchor_rule import family  # noqa: E402

OUT = HERE.parent / "output"
FAMILIES = ("fibroblast_stromal", "endothelial")
K_MDE = stats.norm.ppf(0.975) + stats.norm.ppf(0.80)  # 2.80: 80 % power, two-sided 0.05
# readout: (estimate, p, bootstrap SE, companion p)
READOUTS = {"phase": ("offset_diff_h", "p_offset", "se_offset_diff_h", "p_companion_offset_diff_h"),
            "amplitude": ("log_z2_ratio", "p_log_z2", "se_log_z2_ratio", "p_companion_log_z2_ratio"),
            "noise": ("log_psi_ratio", "p_log_psi", "se_log_psi_ratio", "p_companion_log_psi_ratio")}


def excluded(r, readout: str) -> str:
    """Registered per-readout exclusion: phase / amplitude need z2 > 0 and no truncated z2 in either arm; noise needs
    neither arm at the psi floor."""
    arms = (r.boundary_type_control, r.boundary_type_case)
    if readout in ("phase", "amplitude"):
        if "z2_truncated" in arms or not (r.z2_control > 0 and r.z2_case > 0):
            return "z2 truncated or 0 in an arm"
    elif "psi_floor" in arms:
        return "psi floor in an arm"
    return ""


def long_table(tab: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for r in tab.itertuples(index=False):
        for name, (est, p, se, pc) in READOUTS.items():
            x = dict(target=r.target, compartment=r.compartment, family=r.family, is_anchor=r.is_anchor, readout=name,
                     estimate=getattr(r, est), p=getattr(r, p), se=getattr(r, se), p_companion=getattr(r, pc),
                     donors_control=r.donors_control, donors_case=r.donors_case, log_V_ratio=r.log_V_ratio,
                     reason=excluded(r, name))
            if not x["reason"] and not (np.isfinite(x["estimate"]) and np.isfinite(x["p"])):
                x["reason"] = "not estimable (non-finite estimate or p)"
            mde = K_MDE * x["se"]
            x.update(mde=mde, mde_unit={"phase": "h", "amplitude": "fold amplitude", "noise": "fold psi"}[name],
                     mde_natural=mde if name == "phase" else np.exp(mde / 2) if name == "amplitude" else np.exp(mde))
            x["primary"] = r.family in FAMILIES and not x["reason"]
            rows.append(x)
    return pd.DataFrame(rows)


def h1652(prim: pd.DataFrame) -> pd.DataFrame:
    out = []
    for fam in FAMILIES:
        for name, sign in (("amplitude", -1), ("noise", +1)):
            g = prim[(prim.family == fam) & (prim.readout == name)]
            per = g.groupby("target").apply(lambda d: pd.Series(dict(
                y=np.sum(d.estimate / d.se ** 2) / np.sum(1 / d.se ** 2), se=d.se.mean(), n=len(d))), include_groups=False)
            row = dict(family=fam, readout=name, direction="< 0 (damping)" if sign < 0 else "> 0 (more noise)",
                       n_targets=len(per), targets=";".join(per.index))
            if len(per) < 3:
                out.append(dict(row, verdict="not estimable (< 3 targets)"))
                continue
            m = cs.stats.meta(per.y, per.se, method="hk")
            t = m["random"] / m["se_random"]
            p = stats.t.sf(sign * t, m["k"] - 1)
            out.append(dict(row, pooled=m["random"], se=m["se_random"], lo95=m["lo"], hi95=m["hi"], tau2=m["tau2"],
                            p_one_sided=p))
    res = pd.DataFrame(out).reindex(columns=[*pd.DataFrame(out).columns.union(["p_one_sided", "verdict"], sort=False)])
    ok = res.p_one_sided.notna()
    res["p_holm"], res["verdict"] = np.nan, res.verdict.astype(object)
    if ok.any():
        res.loc[ok, "p_holm"] = multipletests(res.loc[ok, "p_one_sided"], method="holm")[1]
    res.loc[ok, "verdict"] = np.where(res.loc[ok, "p_holm"] < 0.05, "effect", "not detected")
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", type=Path, default=OUT)
    args = ap.parse_args()
    st = pd.read_csv(args.out / "target_status.csv")
    ok = st[st.status == "ok"].target
    tab = pd.concat([pd.read_csv(args.out / "targets" / f"{t}_table.csv") for t in ok], ignore_index=True)
    tab = tab[tab.status == "ok"].copy()
    tab["family"] = tab.compartment.map(family)
    lt = long_table(tab)
    prim = lt[lt.primary].copy()
    prim["p_holm"] = multipletests(prim.p, method="holm")[1] if len(prim) else []
    prim["class"] = np.where(prim.p_holm < 0.05, "effect", "not detected")
    prim["specificity"] = np.where(prim.p_companion < 0.05, "target-specific",
                                   np.where(prim.p_companion.notna(), "not target-specific", "no companion"))
    # registered: a phase effect also requires the shift not to be shared by the companion set
    ph = (prim.readout == "phase") & (prim["class"] == "effect") & ~(prim.p_companion < 0.05)
    prim.loc[ph, "class"] = "effect shared with the companion (not internal desynchrony)"
    lt = lt.merge(prim[["target", "compartment", "readout", "p_holm", "class", "specificity"]],
                  on=["target", "compartment", "readout"], how="left")
    other = pd.Series(np.where(lt.reason != "", "excluded: " + lt.reason, "descriptive (non-primary family)"), index=lt.index)
    lt["class"] = lt["class"].fillna(other)
    lt.to_csv(args.out / "readout.csv", index=False)
    prim.to_csv(args.out / "h1651.csv", index=False)
    h2 = h1652(prim)
    h2.to_csv(args.out / "h1652.csv", index=False)
    arms = [c for c in st.columns if c.startswith("arm_") and c.endswith("_clock_share")]
    gs = st.loc[st.status == "ok", ["target", *arms]].copy()
    gs["case_minus_control"] = gs["arm_case_clock_share"] - gs["arm_control_clock_share"]
    gs.to_csv(args.out / "gate_share.csv", index=False)
    summ = [dict(item="targets estimable", value=int(len(ok))),
            dict(item="primary rows x readouts", value=int(len(prim))),
            *[dict(item=f"primaries {r}", value=int((prim.readout == r).sum())) for r in READOUTS],
            dict(item="H165.1 effects", value=int((prim["class"] == "effect").sum())),
            *[dict(item=f"median MDE {r} ({u})", value=float(prim[prim.readout == r].mde_natural.median()))
              for r, u in (("phase", "h"), ("amplitude", "fold"), ("noise", "fold"))],
            *[dict(item=f"H165.2 {r.family} {r.readout}", value=r.verdict) for r in h2.itertuples()]]
    pd.DataFrame(summ).to_csv(args.out / "summary.csv", index=False)
    print(pd.DataFrame(summ).to_string(index=False))


if __name__ == "__main__":
    main()
