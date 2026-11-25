"""Slide figure (results, 5.1): chronic disease leaves the clock of its own organs in place: the desynchrony predicted
from mouse is not seen, at stated detection limits.

GTEx (fast and ventilated deaths pooled), cases matched to up to 3 controls within cohort (A164 rule). Rows: type 2
diabetes and obesity (module: liver, both adiposes, skeletal muscle), hypertension (heart, arteries, kidney cortex),
ischaemic heart disease (heart, coronary artery, aorta), COPD and smoking (lung), depression (brain regions).
A. Coupling of the disease's own organs to the body's clock state, case minus control at a common host spread (module
   d*), 90 % interval; dashed: the registered equivalence bound +-0.15 (A164 H164.1-2 for diabetes and obesity; A233
   H233.c.3, not an equivalence test, for the others).
B. Peak time of the disease's own organs, case minus control (module offset, drawn + = later: the clocksync offset
   negated), 90 % interval; dashed: +-1 h. The mouse high-fat diet advances the liver by hours (A235: excluded).
   Diabetes and obesity: A235 (A164 host, pooled); the others: A233.
C. Loudness of the disease's own organs relative to the companion organs (tibial nerve, both skins): module minus
   companion difference of the log z2 ratio (case over control), interval +-1.645 SD of the permutation null; grey
   tick: the MDE_80 on each side (A306; T2D, hypertension, smoking; smoking not powered, open). Other conditions not
   tested.

Status: A164, A233, A235, A306 registered; display post hoc (user, 2026-11-25). GTEx aggregates only.
Data: exploratory/A164_gtex_metabolic_matched/output/a164_tests.csv; exploratory/A233_chronic_conditions/output/
score_pooled.csv; exploratory/A235_metabolic_offset/output/score_pooled.csv; exploratory/A306_disease_gain_departure/
output/score_tests.csv.
Output: presentation/figures/chronic_disease/chronic_disease.{svg,png} (wide).

    python scripts/chronic_disease.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D

import paperfig as pf

PAPER = Path(__file__).resolve().parents[2]
EXP = PAPER / "exploratory"
OUT = PAPER / "presentation" / "figures" / "chronic_disease"
CS = sns.color_palette().as_hex()[0]
MUTED = pf.GREY
ROWS = {"t2d": "type 2 diabetes", "obese": "obesity", "htn": "hypertension", "ihd": "ischaemic heart disease",
        "copd": "COPD", "depression": "depression", "smoking": "smoking"}


def tables():
    a164 = pd.read_csv(EXP / "A164_gtex_metabolic_matched/output/a164_tests.csv")
    a164 = a164[(a164.statistic == "module") & (a164.set == "pooled")].set_index("design")
    a233 = pd.read_csv(EXP / "A233_chronic_conditions/output/score_pooled.csv")
    a233 = a233[(a233.selection == "module") & (a233.set == "pooled")]
    a235 = pd.read_csv(EXP / "A235_metabolic_offset/output/score_pooled.csv")
    a235 = a235[(a235.host == "A164") & (a235.readout == "offset_h") & (a235.set == "pooled")].set_index("design")
    rows = []
    for k in ROWS:
        if k in a164.index:
            d = a164.loc[k]
            dstar = (d.D, d.lo90, d.hi90, d.mde_80)
            o = a235.loc[k]
        else:
            d = a233[(a233.design == k) & (a233.readout == "d_star")].iloc[0]
            dstar = (d.estimate, d.lo90, d.hi90, d.mde_80)
            o = a233[(a233.design == k) & (a233.readout == "offset_h")].iloc[0]
        rows.append(dict(condition=k, d_star=dstar[0], d_lo=dstar[1], d_hi=dstar[2], d_mde=dstar[3],
                         peak_h=-o.estimate, peak_lo=-o.hi90, peak_hi=-o.lo90, peak_mde=o.mde_80))
    g = pd.read_csv(EXP / "A306_disease_gain_departure/output/score_tests.csv")
    g = g[g.readout == "gain"].set_index("condition")
    return pd.DataFrame(rows).set_index("condition"), g


def forest(ax, t, col, lo, hi, bound, xlabel, xlim):
    y = np.arange(len(t))
    ax.hlines(y, t[lo], t[hi], color=CS, lw=1.6, zorder=2)
    ax.scatter(t[col], y, s=32, color=CS, zorder=3)
    for b in (-bound, bound):
        ax.axvline(b, color=MUTED, lw=0.8, ls="--", zorder=0)
    ax.axvline(0, color=MUTED, lw=0.6, zorder=0)
    ax.set_xlim(*xlim)
    ax.set_ylim(len(t) - 0.5, -0.5)
    ax.set_xlabel(xlabel)
    ax.tick_params(axis="y", length=0)


def loudness(ax, g, order):
    for i, k in enumerate(order):
        if k not in g.index:
            continue
        r = g.loc[k]
        h = 1.645 * r.null_sd
        open_ = r.p_holm != r.p_holm  # not scored under Holm: not powered
        ax.plot([r.delta - h, r.delta + h], [i, i], color=CS, lw=1.6, zorder=2)
        ax.scatter([r.delta], [i], s=32, color="white" if open_ else CS, edgecolor=CS, lw=1.2, zorder=3)
        for m in (-r.mde_80, r.mde_80):
            ax.plot([m, m], [i - 0.25, i + 0.25], color="#BBBBBB", lw=2, zorder=1)
    ax.axvline(0, color=MUTED, lw=0.6, zorder=0)
    ax.set_ylim(len(order) - 0.5, -0.5)
    ax.set_xlim(-1.1, 0.8)
    ax.set_xlabel("loudness, own organs\nvs other organs [log ratio]")
    ax.tick_params(axis="y", length=0)


def main():
    pf.use_style("slide")
    OUT.mkdir(parents=True, exist_ok=True)
    t, g = tables()
    lay = pf.Layout.slide("ABC", kind="wide", col_ratios=[1, 1, 1])
    forest(lay.axes["A"], t, "d_star", "d_lo", "d_hi", 0.15, "coupling of own organs,\ncase - control [d*]",
           (-0.45, 0.45))
    lay.axes["A"].set_yticks(range(len(t)), [ROWS[k] for k in t.index])
    forest(lay.axes["B"], t, "peak_h", "peak_lo", "peak_hi", 1.0, "peak time of own organs,\ncase - control [h]",
           (-4.6, 2.0))
    loudness(lay.axes["C"], g, list(t.index))
    for k in "BC":
        lay.axes[k].set_yticks(range(len(t)), [""] * len(t))
    h = [Line2D([], [], color=MUTED, lw=0.8, ls="--", label="equivalence bound (±0.15, ±1 h)"),
         Line2D([], [], color="#BBBBBB", lw=2, label="detection limit (80 % power)"),
         Line2D([], [], ls="", marker="o", ms=6, mfc="white", mec=CS, label="not powered")]
    lay.shared_legend("ABC", h, [x.get_label() for x in h], side="bottom", ncol=3, fontsize=10, handletextpad=0.4)
    fig = lay.finalize()
    pd.set_option("display.width", 200)
    print(t.round(3).to_string())
    print(g[["delta", "null_sd", "mde_80", "p", "p_holm", "verdict"]].round(3).to_string())
    print(pf.collisions(fig))
    pf.save(fig, OUT / "chronic_disease")


if __name__ == "__main__":
    main()
