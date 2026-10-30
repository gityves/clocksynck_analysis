"""A167 scoring (registered in paper/docs/registered_experiments.md, A167 and amendment 1): H167.1 (d*, Holm over the
3 primaries), H167.2 (phase, amplitude, noise, Holm over 9), the liver-wide specificity qualifier, MDEs, and the
descriptive secondaries (D1, leave-one-study-out, other compartments).

    python score_a167.py
Writes output/: h1671.csv, h1672.csv, secondary.csv, summary.csv.
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "A165_disease_readouts" / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from statsmodels.stats.multitest import multipletests  # noqa: E402

from score_a165 import K_MDE, READOUTS, excluded  # noqa: E402

OUT = HERE.parent / "output"
PRIMARIES = ["liver:stellate", "liver:endothelial", "liver:hepatocyte"]


def specific(p_comp) -> str:
    return "compartment-specific" if p_comp < 0.05 else "liver-wide, not compartment-specific" \
        if np.isfinite(p_comp) else "no companion"


def h1671(t: pd.DataFrame, power: pd.DataFrame) -> pd.DataFrame:
    x = t[t.compartment.isin(PRIMARIES)][["compartment", "d_star", "p_star", "se_d_star", "d_star_minus_companion",
                                           "p_companion_d_star", "log_V_ratio"]].copy()
    x["lo95_d_star"], x["hi95_d_star"] = x.d_star - 1.96 * x.se_d_star, x.d_star + 1.96 * x.se_d_star  # bootstrap SE
    x = x.merge(power[["compartment", "mde_dstar"]], on="compartment", how="left")
    ok = x.p_star.notna()
    x["p_holm"] = np.nan
    x.loc[ok, "p_holm"] = multipletests(x.loc[ok, "p_star"], method="holm")[1]
    x["class"] = np.where(~ok, "not estimable", np.where(x.p_holm < 0.05, "effect", "not detected"))
    x["specificity"] = np.where(x["class"] == "effect", x.p_companion_d_star.map(specific), "")  # qualifies effects only
    return x


def h1672(t: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for r in t[t.compartment.isin(PRIMARIES)].itertuples(index=False):
        for name, (est, p, se, pc) in READOUTS.items():
            e, s = getattr(r, est), getattr(r, se)
            mde = K_MDE * s
            rows.append(dict(compartment=r.compartment, readout=name, estimate=e, se=s, p=getattr(r, p),
                             p_companion=getattr(r, pc), reason=excluded(r, name), mde=mde,
                             mde_natural=mde if name == "phase" else np.exp(mde / 2) if name == "amplitude" else np.exp(mde)))
    x = pd.DataFrame(rows)
    ok = (x.reason == "") & x.p.notna()
    x["p_holm"] = np.nan
    x.loc[ok, "p_holm"] = multipletests(x.loc[ok, "p"], method="holm")[1]
    x["class"] = np.where(~ok, "excluded: " + x.reason.where(x.reason != "", "not estimable"),
                          np.where(x.p_holm < 0.05, "effect", "not detected"))
    x["specificity"] = np.where(x["class"] == "effect", x.p_companion.map(specific), "")
    # registered: a phase effect also requires the companion phase p < 0.05
    ph = (x.readout == "phase") & (x["class"] == "effect") & ~(x.p_companion < 0.05)
    x.loc[ph, "class"] = "effect shared with the liver (not internal desynchrony)"
    return x


def main():
    an = pd.read_csv(OUT / "analyses.csv").set_index("analysis")
    power = pd.read_csv(OUT / "power.csv")
    summ = [dict(item=f"{a} status", value=f"{r.status} {r.reason if isinstance(r.reason, str) else ''}".strip())
            for a, r in an.iterrows()]
    if an.loc["d2", "status"] == "ok":
        t = pd.read_csv(OUT / "contrast_d2.csv")
        a, b = h1671(t, power), h1672(t)
        a.to_csv(OUT / "h1671.csv", index=False)
        b.to_csv(OUT / "h1672.csv", index=False)
        summ += [dict(item=f"H167.1 {r.compartment}", value=f"{r['class']} (d* {r.d_star:+.3f}, Holm {r.p_holm:.3g}, "
                                                           f"MDE {r.mde_dstar:.3f}{'; ' + r.specificity if r.specificity else ''})")
                 for _, r in a.iterrows()]
        summ.append(dict(item="H167.2 effects", value=int((b["class"] == "effect").sum())))
    sec = []
    for f in sorted(OUT.glob("contrast_*.csv")):
        t = pd.read_csv(f)
        cols = ["analysis", "compartment", "d_star", "p_star", "offset_diff_h", "p_offset", "log_z2_ratio", "p_log_z2",
                "log_psi_ratio", "p_log_psi", "p_companion_d_star", "log_V_ratio", "is_primary"]
        sec.append(t[[c for c in cols if c in t]])
    pd.concat(sec, ignore_index=True).to_csv(OUT / "secondary.csv", index=False)
    s = pd.DataFrame(summ)
    s.to_csv(OUT / "summary.csv", index=False)
    print(s.to_string(index=False))


if __name__ == "__main__":
    main()
