"""Slide figure (results, drives, 4.4): a systemic glucocorticoid-like state displaces some organs' clocks from the body's
time.

GTEx fast deaths; aggregates only. One panel (B below; half slide). `volume` (A) is drawn on slide 3.5, panel D
(user, 2026-11-25: the volume link belongs to 3.5, with its agonal caveat, A367).
A. Person clock volume (log host radius from all tissues, 237 persons) regressed on the person's systemic states and on
   technical terms (A360, registered; number of tissues and time of death in the model, not shown): standardised
   coefficients with 90 % intervals. Glucocorticoid-target state on the 14 clock genes (A360) and on the core loop
   (A361); grey: the same coefficient for 200 expression-matched (light) and 200 stress-responsive (dark) gene sets,
   median and 5-95 % range (A361).
B. Per tissue (25, organ-system colours), the share of its clock-state departure from the body's clock explained by the
   person's signal factor read from the person's other tissues (R^2, A357 H357.2); grey tick: the tissue's 95th
   percentile over 100 expression-matched sets. Differences between tissues are descriptive (H357.1 not scored).

Status: A360 H360.1 pass, A361 H361.1 pass (core loop blind), A357 H357.2 pass; display choices post hoc (user,
2026-11-25).
Data: exploratory/A358_A360_gain_sources/output/a360_coefficients.csv; exploratory/A361_A363_volume_mechanism/output/
{systemic_tests,systemic_companions}.csv; exploratory/A357_signal_departure/output/score_tissues.csv.
Output: presentation/figures/systemic_state/systemic_state.{svg,png} (half).

    python scripts/systemic_state.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D

import paperfig as pf
from _systems import COLOR, gtex_system

PAPER = Path(__file__).resolve().parents[2]
EXP = PAPER / "exploratory"
OUT = PAPER / "presentation" / "figures" / "systemic_state"
CS = sns.color_palette().as_hex()[0]
MATCHED, RESPONSIVE = "#BBBBBB", "#6E6E6E"
SHORT = {"Adipose - Subcutaneous": "fat, subcutaneous", "Adipose - Visceral (Omentum)": "fat, visceral",
         "Adrenal Gland": "adrenal", "Artery - Aorta": "aorta", "Artery - Coronary": "coronary artery",
         "Artery - Tibial": "tibial artery", "Breast - Mammary Tissue": "breast", "Colon - Sigmoid": "colon, sigmoid",
         "Colon - Transverse": "colon, transverse", "Esophagus - Gastroesophageal Junction": "oes.-gastric junction",
         "Esophagus - Mucosa": "oes. mucosa", "Esophagus - Muscularis": "oes. muscularis",
         "Heart - Atrial Appendage": "heart, atrium", "Heart - Left Ventricle": "heart, ventricle",
         "Kidney - Cortex": "kidney", "Liver": "liver", "Lung": "lung", "Muscle - Skeletal": "skeletal muscle",
         "Nerve - Tibial": "tibial nerve", "Pituitary": "pituitary", "Prostate": "prostate",
         "Skin - Not Sun Exposed (Suprapubic)": "skin, not exposed", "Skin - Sun Exposed (Lower leg)": "skin, exposed",
         "Stomach": "stomach", "Thyroid": "thyroid"}


def volume(ax, coef, tests, comp, xlabel="effect on the person's shared\nloudness [standardised coefficient]"):
    """Standardised effects on the person's clock volume (log host radius; A360, A361) with companions."""
    c = coef[coef.model == "registered"].set_index("term")
    t = tests.set_index("test").loc["H361.1"]
    rows = [("glucocorticoid state\n14 clock genes", c.loc["gc", "std_coef"], c.loc["gc", "lo90"], c.loc["gc", "hi90"],
             "14 genes"),
            ("glucocorticoid state\ncore loop", t.gc_coef, t.gc_lo90, t.gc_hi90, "core"),
            ("immediate-early state", c.loc["ieg", "std_coef"], c.loc["ieg", "lo90"], c.loc["ieg", "hi90"], None),
            ("RNA quality (RIN)", c.loc["rin", "std_coef"], c.loc["rin", "lo90"], c.loc["rin", "hi90"], None),
            ("ischaemic time", c.loc["ischaemic_h", "std_coef"], c.loc["ischaemic_h", "lo90"],
             c.loc["ischaemic_h", "hi90"], None)]
    ax.axvline(0, color=pf.GREY, lw=0.8, ls="--", zorder=0)
    out = []
    for i, (lab, x, lo, hi, gs) in enumerate(rows):
        if gs:
            for kind, col, dy in (("matched", MATCHED, -0.22), ("responsive", RESPONSIVE, 0.22)):
                v = comp[(comp.kind == kind) & (comp.gene_set == gs)].gc
                ax.plot(v.quantile([0.05, 0.95]), [i + dy] * 2, color=col, lw=2.2, solid_capstyle="butt", zorder=1)
                ax.scatter([v.median()], [i + dy], s=18, color=col, zorder=2)
                out.append(dict(row=lab, family=kind, median=v.median(), q05=v.quantile(0.05)))
        ax.plot([lo, hi], [i, i], color=CS, lw=1.6, zorder=3)
        ax.scatter([x], [i], s=40, color=CS, zorder=4)
        out.append(dict(row=lab, family="clock", median=x, q05=lo))
    ax.set_yticks(range(len(rows)), [r[0] for r in rows])
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(len(rows) - 0.5, -0.6)
    ax.set_xlabel(xlabel)
    return pd.DataFrame(out)


def departure(ax, t):
    t = t.sort_values("r2_sig").reset_index(drop=True)
    for i, r in t.iterrows():
        ax.plot([r.r2_matched_q95] * 2, [i - 0.35, i + 0.35], color=MATCHED, lw=1.6, zorder=1)
        ax.scatter([r.r2_sig], [i], s=26, color=COLOR[gtex_system(r.tissue)], zorder=3)
    ax.set_yticks(range(len(t)), [SHORT[x] for x in t.tissue], fontsize=8)
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(-0.7, len(t) - 0.3)
    ax.set_xlim(-0.004, None)
    ax.set_xlabel("share of the organ's departure explained\nby the person's state [R²]")


def main():
    pf.use_style("slide")
    OUT.mkdir(parents=True, exist_ok=True)
    tis = pd.read_csv(EXP / "A357_signal_departure/output/score_tissues.csv")
    lay = pf.Layout.slide("B", kind="half")
    departure(lay.axes["B"], tis)
    present = {gtex_system(x) for x in tis.tissue}
    h = [Line2D([], [], color=MATCHED, lw=2.2, label="matched genes, 95th pct")]
    h += [Line2D([], [], ls="", marker="o", ms=7, color=c, label=s) for s, c in COLOR.items() if s in present]
    lay.shared_legend("B", h, [x.get_label() for x in h], side="bottom", ncol=3, fontsize=9.5, handletextpad=0.3,
                      columnspacing=1.0)
    fig = lay.finalize()
    pd.set_option("display.width", 200)
    print(tis[["tissue", "r2_sig", "r2_matched_q95"]].sort_values("r2_sig", ascending=False).round(3).to_string(index=False))
    print("median R2", round(tis.r2_sig.median(), 3), "median matched median", round(tis.r2_matched_median.median(), 3),
          "above own q95:", int((tis.r2_sig > tis.r2_matched_q95).sum()), "of", len(tis))
    print(pf.collisions(fig))
    pf.save(fig, OUT / "systemic_state")


if __name__ == "__main__":
    main()
