"""Slide figure (results, local disease, 6.1): tumours keep their patient's time, carried by their fibroblasts, and
express less of it as they proliferate.

A. Tumour-normal pairing (TCGA, 517 pairs, 10 cancers; A138 H138.4, registered): Poisson-corrected |correlation| of the
   tumour's and its own patient's normal clock state, per cancer (dots; open grey: the cancer's pairing-null mean) and
   pooled (14 genes, filled; core loop, open, A336, with its 90 % bootstrap interval); dark tick: pooled pairing-null
   q95; light ticks: the expression-matched set q95 (14 genes and core loop).
B. Who carries it (A346, registered; single-cell atlases): pairing of each compartment of the tumour with the same
   compartment of the patient's normal tissue (cancer cells: with normal epithelium), 14 genes (filled) and core loop
   (open); grey ticks: the expression-matched set q95 of each readout.
C. Proliferation (A304, registered): high- against low-proliferation tumours, pooled log ratio of gain (the tumour's
   expression of its patient's time) and of departure (what it does on its own), interval +-1.645 null SD; grey: the
   5-95 % range of 200 expression-matched sets.

Status: A138 H138.4 pass (core loop at the line, A336); A346 H346.3 pass on 14 genes, H346.1 descriptive (fibroblasts
meet every condition, endothelium falls below matched genes); A304 H304.1 gain pass, departure fail (powered).
Data: exploratory/A138_cancer_estimability/output/score/{tests,per_unit}.csv; exploratory/A336_two_readouts/output/
part1/a138.csv; exploratory/A346_tumour_stroma_core/output/score/tests.csv; exploratory/A304_cancer_gain_departure/
output/score/tests.csv. Output: presentation/figures/tumours/tumours.{svg,png} (wide).

    python scripts/tumours.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D

import paperfig as pf

PAPER = Path(__file__).resolve().parents[2]
EXP = PAPER / "exploratory"
OUT = PAPER / "presentation" / "figures" / "tumours"
CS = sns.color_palette().as_hex()[0]
MUTED, MATCHED = pf.GREY, "#BBBBBB"
CANCER = {"BLCA": "bladder", "BRCA": "breast", "COAD": "colon", "KIRC": "kidney, clear cell", "KIRP": "kidney, papillary",
          "LIHC": "liver", "LUAD": "lung", "PRAD": "prostate", "THCA": "thyroid", "UCEC": "endometrium"}
COMP = {"fibroblast": "fibroblasts", "endothelial": "endothelium", "epithelium": "cancer cells"}


def tick(ax, x, y, h=0.3, color=MATCHED):
    ax.plot([x, x], [y - h, y + h], color=color, lw=2, zorder=1)


def pairing(ax):
    t = pd.read_csv(EXP / "A138_cancer_estimability/output/score/tests.csv").set_index("hypothesis").loc["H138.4"]
    u = pd.read_csv(EXP / "A138_cancer_estimability/output/score/per_unit.csv")
    u = u[u.hypothesis == "H138.4"].sort_values("statistic", ascending=False)
    core = pd.read_csv(EXP / "A336_two_readouts/output/part1/a138.csv").set_index("id").loc["A138_H138.4"]
    for i, r in enumerate(u.itertuples()):
        ax.scatter([r.null_mean], [i], s=22, facecolor="white", edgecolor=MUTED, zorder=2)
        ax.scatter([r.statistic], [i], s=26, color=CS, zorder=3)
    y = len(u) + 0.6
    ax.scatter([t.statistic], [y], s=44, color=CS, zorder=3)
    lo, hi = (float(x) for x in str(core.notes).split("core 90 % bootstrap ")[1].split(";")[0].split(" to "))
    ax.plot([lo, hi], [y + 0.7] * 2, color=CS, lw=1.4, zorder=2)
    ax.scatter([core.value_core], [y + 0.7], s=44, facecolor="white", edgecolor=CS, lw=1.3, zorder=3)
    tick(ax, t.null_q95, y, color=MUTED)
    tick(ax, t.companion_q95, y)
    tick(ax, float(str(core.notes).split("matched q95 ")[1].split(",")[0]), y + 0.7)
    ax.set_yticks(list(range(len(u))) + [y, y + 0.7],
                  [CANCER[x] for x in u.unit] + ["pooled", "pooled, core loop"])
    ax.set_ylim(y + 1.3, -0.7)
    ax.set_xlim(0, 0.62)
    ax.set_xlabel("tumour follows own normal\n[pairing, |r|]")
    ax.tick_params(axis="y", length=0)
    return u[["unit", "statistic", "null_mean"]], t, core


def compartments(ax):
    t = pd.read_csv(EXP / "A346_tumour_stroma_core/output/score/tests.csv")
    t["readout"] = t.readout.astype(str)
    for i, (k, lab) in enumerate(COMP.items()):
        for dy, rd, open_ in ((-0.15, "14", False), (0.15, "core", True)):
            r = t[(t.kind == k) & (t.readout == rd)].iloc[0]
            tick(ax, r.matched_q95, i + dy, 0.12)
            ax.scatter([r.clock], [i + dy], s=40, facecolor="white" if open_ else CS, edgecolor=CS, lw=1.3, zorder=3)
    ax.set_yticks(range(len(COMP)), list(COMP.values()))
    ax.set_ylim(len(COMP) - 0.5, -0.5)
    ax.set_xlim(0.25, 0.8)
    ax.set_xlabel("tumour compartment follows\nthe same in normal tissue [|r|]")
    ax.tick_params(axis="y", length=0)
    return t[["kind", "readout", "clock", "matched_q95"]]


def proliferation(ax):
    t = pd.read_csv(EXP / "A304_cancer_gain_departure/output/score/tests.csv")
    t = t[t.contrast == "H304.1"].set_index("readout")
    for i, (k, lab) in enumerate((("gain", "expression of the\npatient's time"), ("departure", "own departure"))):
        r = t.loc[k]
        ax.plot([r.matched_q05, r.matched_q95], [i + 0.22] * 2, color=MATCHED, lw=2.2, zorder=1)
        h = 1.645 * r.null_sd
        ax.plot([r.statistic - h, r.statistic + h], [i, i], color=CS, lw=1.6, zorder=2)
        ax.scatter([r.statistic], [i], s=40, color=CS, zorder=3)
    ax.axvline(0, color=MUTED, lw=0.8, ls="--", zorder=0)
    ax.set_yticks([0, 1], ["expression of the\npatient's time", "own departure"])
    ax.set_ylim(1.6, -0.6)
    ax.set_xlim(-0.85, 0.5)
    ax.set_xlabel("high / low proliferation\n[log ratio]")
    ax.tick_params(axis="y", length=0)
    return t[["statistic", "null_sd", "p", "matched_q05", "matched_q95"]]


def main():
    pf.use_style("slide")
    OUT.mkdir(parents=True, exist_ok=True)
    lay = pf.Layout.slide("AxBC", kind="wide", col_ratios=[1.05, 0.12, 0.9, 0.9], blank=("x",))
    a, t, core = pairing(lay.axes["A"])
    b = compartments(lay.axes["B"])
    c = proliferation(lay.axes["C"])
    h = [Line2D([], [], ls="", marker="o", ms=7, color=CS, label="clock genes"),
         Line2D([], [], ls="", marker="o", ms=7, mfc="white", mec=CS, mew=1.3, label="core loop"),
         Line2D([], [], ls="", marker="o", ms=6, mfc="white", mec=MUTED, label="pairing null mean"),
         Line2D([], [], color=MUTED, lw=2, label="pairing null, 95th pct"),
         Line2D([], [], color=MATCHED, lw=2.2, label="matched genes (95th pct; C: 5-95 %)")]
    lay.shared_legend("AxBC", h, [x.get_label() for x in h], side="bottom", ncol=5, fontsize=9.5, handletextpad=0.3)
    fig = lay.finalize()
    pd.set_option("display.width", 200)
    print(a.round(3).to_string(index=False))
    print("pooled", round(t.statistic, 3), "null q95", round(t.null_q95, 3), "matched q95", round(t.companion_q95, 3),
          "core", round(core.value_core, 3))
    print(b.round(3).to_string(index=False))
    print(c.round(3).to_string())
    print(pf.collisions(fig))
    pf.save(fig, OUT / "tumours")


if __name__ == "__main__":
    main()
