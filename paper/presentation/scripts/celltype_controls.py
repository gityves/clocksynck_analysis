"""Slide figure (results, validation, cell-type level): untimed offsets between cell types of one sample match known time in
mouse tissues, also when the cell classes move apart; without a clock, pairs stay correlated but offsets no longer track.

A. Phase offset between two cell types of the same sample [h], clocksync (untimed pairwise, `clocksync.pairwise`) against
   timed (each compartment's phase against ZT): mouse colon chrono-atlas (A273: undamaged, 28 pairs; regenerating after
   injury, 15 pairs; pairs with both compartments timed-rhythmic, timed coupling >= 0.5; 8 mice per arm, ZT 1/7/13/19)
   and mouse aorta (A279, Auerbach 2026, aligned arm, 6 pairs of EC, SMC, fibroblast, macrophage; 14 libraries). The
   untimed value is drawn on the turn nearest the timed one (timed + error); dashed: identity.
B. Aorta, |r| between the clock-plane states of two cell types over libraries, per pair: clock intact (aligned arm,
   filled; upper row) and Bmal1 knockout (open; lower row; timed coupling 0.00-0.20 in every class, no clock), each
   beside the 5-95 % range (bar) and median (tick) of 200 expression-matched (light grey) and 200 stress-responsive
   (dark grey) gene sets of the same arm (A354; H354.2 PASS, H354.1 FAIL).

Status: A273 registered calibration (V273.1, V273.2 PASS; not blind for the per-class timed phases); A279 registered
(V279.1 PASS; not blind for A201's per-cell-type shifts); knockout descriptive in A279; companions A354
registered (not blind for the clock |r|, blind for every companion value).
Data: exploratory/A273_pairwise_offset_truth/output/pairs.csv; exploratory/A279_aorta_pairwise_truth/output/pairs.csv;
exploratory/A279_aorta_pairwise_truth/output/a354/pairs.csv.
Output: presentation/figures/celltype_controls/celltype_controls.{svg,png} (half).

    python scripts/celltype_controls.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D

import paperfig as pf

PAPER = Path(__file__).resolve().parents[2]
EXP = PAPER / "exploratory"
OUT = PAPER / "presentation" / "figures" / "celltype_controls"
TAB = sns.color_palette().as_hex()
COLON, AORTA, INK = TAB[2], TAB[3], "#222222"  # organ systems: digestive green, cardiovascular red
MATCHED, RESPONSIVE = "#BBBBBB", "#6E6E6E"
SHORT = {"EC": "EC", "Fibroblast": "fib", "SMC": "SMC", "Macrophage": "mac"}


def offsets(ax):
    c = pd.read_csv(EXP / "A273_pairwise_offset_truth/output/pairs.csv")
    c = c[c.arm.isin(["C-DC", "R-DC"]) & c.rhythmic]
    a = pd.read_csv(EXP / "A279_aorta_pairwise_truth/output/pairs.csv").query("arm == 'aligned'")
    ax.plot([-12, 12], [-12, 12], color=pf.GREY, lw=1, ls="--", zorder=1)
    for d, col, face, lab in ((c[c.arm == "C-DC"], COLON, COLON, "colon"),
                              (c[c.arm == "R-DC"], COLON, "white", "colon, regenerating"),
                              (a, AORTA, AORTA, "aorta")):
        ax.scatter(d.timed_h, d.timed_h + d.error_h, s=34, facecolor=face, edgecolor=col, lw=1.2, zorder=3, label=lab)
    ax.set_xlim(-12, 12)
    ax.set_ylim(-12, 12)
    ax.set_xticks([-12, -6, 0, 6, 12])
    ax.set_yticks([-12, -6, 0, 6, 12])
    ax.set_aspect("equal")
    ax.set_xlabel("offset between cell types,\ntimed [h]")
    ax.set_ylabel("offset between cell types,\nclocksync [h]")
    ax.legend(loc="upper left", fontsize=10, frameon=False, handletextpad=0.1, borderaxespad=0.1)
    return c, a


def knockout(ax):
    """Per pair, clock intact (top) and Bmal1 knockout (bottom): clock |r| beside the 5-95 % range and median of the
    expression-matched (light grey) and stress-responsive (dark grey) gene sets of the same arm (A354)."""
    d = pd.read_csv(EXP / "A279_aorta_pairwise_truth/output/a354/pairs.csv")
    d = d[d.arm.isin(["aligned", "bmal1_ko"])].assign(pair=lambda x: x.a.map(SHORT) + "-" + x.b.map(SHORT))
    order = (d[d.arm == "bmal1_ko"].sort_values(["stromal", "clock_r"]).pair.tolist())
    ticks, labels = [], []
    for i, pair in enumerate(order):
        for arm, dy in (("aligned", -0.2), ("bmal1_ko", 0.2)):
            r = d[(d.pair == pair) & (d.arm == arm)].iloc[0]
            y = i + dy
            for fam, col, off in (("matched", MATCHED, -0.07), ("responsive", RESPONSIVE, 0.07)):
                ax.plot([r[f"{fam}_q05"], r[f"{fam}_q95"]], [y + off] * 2, color=col, lw=3, solid_capstyle="butt",
                        zorder=1)
                ax.plot([r[f"{fam}_median"]] * 2, [y + off - 0.06, y + off + 0.06], color=INK, lw=0.8, zorder=2)
            ax.scatter([r.clock_r], [y], s=40, zorder=3, facecolor=AORTA if arm == "aligned" else "white",
                       edgecolor=AORTA, lw=1.2)
        ticks.append(i)
        labels.append(pair)
    ax.set_yticks(ticks, labels)
    ax.tick_params(axis="y", length=0)
    ax.set_xlim(0, 1.03)
    ax.set_ylim(len(order) - 0.45, -0.55)
    ax.set_xlabel("correlation between\ncell types, |r|")
    h = [Line2D([], [], ls="", marker="o", ms=7, color=AORTA, label="clock genes, clock intact"),
         Line2D([], [], ls="", marker="o", ms=7, mfc="white", mec=AORTA, label="clock genes, Bmal1 knockout"),
         Line2D([], [], color=MATCHED, lw=3, label="matched genes, 5-95 %"),
         Line2D([], [], color=RESPONSIVE, lw=3, label="stress-responsive genes, 5-95 %")]
    return d, h


def main():
    pf.use_style("slide")
    OUT.mkdir(parents=True, exist_ok=True)
    lay = pf.Layout.slide("AB", kind="half", height=120, col_ratios=[1, 1])
    c, a = offsets(lay.axes["A"])
    w, h = knockout(lay.axes["B"])
    lay.shared_legend("AB", h, [x.get_label() for x in h], side="bottom", ncol=2, fontsize=10, columnspacing=1.0)
    fig = lay.finalize()
    for name, d in (("colon", c), ("aorta", a)):
        print(name, len(d), "pairs; median |error|", round(d.error_h.abs().median(), 2))
    print(w[["arm", "pair", "clock_r", "matched_q95", "responsive_q95"]].round(3).to_string(index=False))
    print(pf.collisions(fig))
    pf.save(fig, OUT / "celltype_controls")


if __name__ == "__main__":
    main()
