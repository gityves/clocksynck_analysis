"""Slide figure (results, R2.2): the organ sets the time; each cell class sets how strongly it expresses it.
Lymphocytes are the zero-gain end of that scale: no readable clock of their own.

A. Gain of each cell class relative to the endothelium of the same atlas (A293 Part 3, registered, H293.4 pass): ratio
   of the shared-rhythm amplitudes, log scale; coloured dots: one per atlas (mean over the family's compartments);
   black: random-effects pooled value with 95 % interval (pooled +- 1.96 SE). Lymphocytes: no gain is estimated; the
   arrow marks them at the zero-gain end (A281, A287, A294: no readable clock on the full clock-gene set).
B. Lead of the first class within human organs (A389, normal tissue: healthy atlases and tumour-adjacent normal samples,
   10 organs / 608 donors for endothelium-fibroblast; pairwise offsets summed per organ): dots, organs (14 genes);
   black, random-effects pooled with 90 % interval (14 genes); open blue, core loop; open diamonds: mouse colon and aorta
   against known time (A273, A279). Positive = the first class earlier. H389.1 FAIL on the registered companion
   condition; standardised, beyond 93.5-98 % of matched sets (A389 note 2, post hoc).
Lymphocytes: no gain on this scale; they do not share the tissue's clock state (pairwise, A272-A288; coupling median
0.13 in the A17 atlas fits). Whether they keep a clock of their own is open (A281, A321, bound A326: speaker notes).
A also shows the core-loop pooled gains (open blue, A350, blind; mural tau above the bound: "not conserved").

Data: exploratory/A293_organisation_levels/output/part3/{atlases,families}.csv;
exploratory/A293_organisation_levels/output/part2/tests.csv (mouse); exploratory/A389_class_phase_order/output/score/{organs,pooled}.csv; exploratory/A281_lymphocyte_clock/output/pairs.csv;
exploratory/A294_lymphocyte_hormonal/output/pairs.csv;
exploratory/A350_class_gain_core/output/compare.csv; exploratory/L1_lymphocyte_power/output/design_power.csv;
exploratory/L1_lymphocyte_power/scout/output/ratios_vs_fibroblast.csv. Output: presentation/figures/celltypes_health/celltypes_health.{svg,png}.

    python scripts/celltypes_health.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D

import paperfig as pf

PAPER = Path(__file__).resolve().parents[2]
EXP = PAPER / "exploratory"
OUT = PAPER / "presentation" / "figures" / "celltypes_health"
INK, MUTED, MATCHED = "#222222", pf.GREY, "#BBBBBB"
TAB = sns.color_palette().as_hex()
CORE = TAB[0]
CELL = {"stromal-vascular": TAB[5], "macrophage": TAB[1], "lymphocyte": TAB[9]}
FAMILY = {  # pair-table family: (label, colour)
    "fibroblast_stromal": ("fibroblast", TAB[5]), "endothelial": ("endothelial", TAB[5]), "mural": ("mural", TAB[5]),
    "epithelial": ("epithelial", TAB[2]), "parenchyma": ("parenchyma", TAB[4]), "other": ("glia, melanocyte,\nprogenitor", TAB[8]),
    "myeloid": ("myeloid", TAB[1]), "lymphocyte": ("lymphocyte", TAB[9])}
ATLAS = {"lung_hlca": "lung", "heart": "heart", "organ_breast": "breast", "organ_colon": "colon (Pelka)",
         "organ_muscle": "muscle", "kidney_sc": "kidney", "kong_ileum": "ileum (Kong)", "kong_colon": "colon (Kong)", "organ_eye_anterior": "anterior eye",
         "organ_optic_nerve": "optic nerve", "ibd_ileum": "ileum (IBD)", "ibd_rectum": "rectum (IBD)"}


GAIN_ORDER = ("fibroblast_stromal", "endothelial", "mural", "parenchyma", "epithelial", "myeloid", "lymphocyte")
PAIRS = {"endothelial-fibroblast_stromal": "endothelial\nvs fibroblast", "endothelial-mural": "endothelial\nvs mural",
         "mural-fibroblast_stromal": "mural\nvs fibroblast"}


def gain(ax, at, fam, core):
    """Relative gain to endothelium per family: atlases (dots), pooled value and 95 % interval (black, 14 genes); open
    blue: the same on the core loop (A350, blind)."""
    fam, core = fam.set_index("family"), core.set_index("family")
    lo_edge = 0.2
    for i, f in enumerate(GAIN_ORDER):
        col = FAMILY[f][1]
        if f == "endothelial":
            ax.scatter([1], [i], s=40, marker="D", color=col, edgecolor="white", lw=0.5, zorder=3)
            continue
        if f == "lymphocyte":  # no gain on this scale: they do not share the tissue's clock state (A272-A288)
            ax.text(0.5, i, "not coupled to the tissue's clock", color=col, va="center", ha="center", fontsize=9)
            continue
        a = at[at.family == f]
        y = i + np.random.default_rng(i).uniform(-0.2, 0.2, len(a))
        ax.scatter(np.exp(a.g), y, s=18, color=col, alpha=0.6, edgecolor="none", zorder=1)
        r = fam.loc[f]
        ax.plot(np.exp([r.pooled - 1.96 * r.se_pooled, r.pooled + 1.96 * r.se_pooled]), [i, i], color=INK, lw=1.6,
                zorder=2)
        ax.scatter([np.exp(r.pooled)], [i], s=34, color=INK, zorder=3)
        if f in core.index:
            c = core.loc[f]
            yc = i + 0.32
            ax.plot(np.exp([c.pooled_core - 1.96 * c.se_pooled_core, c.pooled_core + 1.96 * c.se_pooled_core]), [yc, yc],
                    color=CORE, lw=1.4, zorder=2)
            ax.scatter([np.exp(c.pooled_core)], [yc], s=30, facecolor="white", edgecolor=CORE, lw=1.3, zorder=3)
            if not c.within_bound_core:  # A350 H350.1: tau above the registered bound on the core loop
                ax.annotate("not conserved", (np.exp(c.pooled_core + 1.96 * c.se_pooled_core), yc), xytext=(4, 0),
                            textcoords="offset points", ha="left", va="center", fontsize=8, color=CORE)
    ax.axvline(1, color=MUTED, lw=0.8, ls="--", zorder=0)
    ax.set_xscale("log")
    ax.set_xlim(lo_edge, 2.5)
    ax.set_xticks([0.25, 0.5, 1, 2], ["0.25", "0.5", "1", "2"])
    ax.minorticks_off()
    ax.set_yticks(range(len(GAIN_ORDER)), [FAMILY[f][0] for f in GAIN_ORDER])
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(len(GAIN_ORDER) - 0.5, -0.5)
    ax.set_xlabel("gain relative to endothelium")


A389 = {"endothelial-fibroblast_stromal": "endo-fib", "endothelial-mural": "endo-mural", "mural-fibroblast_stromal": "mural-fib"}


def phase(ax, organs, pooled, tests):
    """Lead of the first class within organs (A389, normal tissue): organs, random-effects pooled (14 genes, core loop),
    mouse against known time. Lead = minus the A389 peak-time offset (positive = first class earlier)."""
    tests = tests.set_index("class_pair")
    for i, (p, label) in enumerate(PAIRS.items()):
        k = A389[p]
        o = organs[(organs.pair == k) & (organs.readout == "14")]
        y = i - 0.12 + np.random.default_rng(i).uniform(-0.12, 0.12, len(o))
        ax.scatter(-o.peak_h, y, s=18, color=CELL["stromal-vascular"], alpha=0.6, edgecolor="none", zorder=1)
        for dy, rd, open_ in ((-0.12, "14", False), (0.12, "core", True)):
            r = pooled[(pooled.pair == k) & (pooled.readout == rd)].iloc[0]
            ax.plot([-r.hi90, -r.lo90], [i + dy] * 2, color=CORE if open_ else INK, lw=1.6, zorder=2)
            ax.scatter([-r.pooled_peak_h], [i + dy], s=34, facecolor="white" if open_ else INK,
                       edgecolor=CORE if open_ else INK, lw=1.3, zorder=3)
        m = tests.loc[p]
        ax.scatter([m.mouse_colon_timed_h, m.mouse_aorta_timed_h], [i + 0.34] * 2, s=26, marker="D", facecolor="white",
                   edgecolor=INK, lw=0.9, zorder=3)
    ax.axvline(0, color=MUTED, lw=0.8, ls="--", zorder=0)
    ax.set_xlim(-2.0, 3.0)
    ax.set_ylim(len(PAIRS) - 0.4, -0.5)
    ax.set_yticks(range(len(PAIRS)), list(PAIRS.values()))
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("lead of the first class [h]")
    h = [Line2D([], [], ls="", marker="o", ms=5, color=CELL["stromal-vascular"], alpha=0.6, label="human organ"),
         Line2D([], [], ls="", marker="D", ms=5, mfc="white", mec=INK, label="mouse, known time")]
    ax.legend(handles=h, loc="upper right", frameon=False, fontsize=8, handletextpad=0.1)
    return pooled[["pair", "readout", "organs", "donors", "pooled_peak_h", "lo90", "hi90", "Q_p", "tau_h"]]


def main():
    pf.use_style("slide")
    OUT.mkdir(parents=True, exist_ok=True)
    lay = pf.Layout.slide("AB", kind="wide", col_ratios=[1.2, 1])
    p3, p2 = EXP / "A293_organisation_levels/output/part3", EXP / "A293_organisation_levels/output/part2"
    gain(lay.axes["A"], pd.read_csv(p3 / "atlases.csv"), pd.read_csv(p3 / "families.csv"),
         pd.read_csv(EXP / "A350_class_gain_core/output/compare.csv"))
    a389 = EXP / "A389_class_phase_order/output/score"
    org, pool = pd.read_csv(a389 / "organs.csv"), pd.read_csv(a389 / "pooled.csv")
    b = phase(lay.axes["B"], org[org["sample"] == "normal"], pool[pool["sample"] == "normal"], pd.read_csv(p2 / "tests.csv"))
    h = [Line2D([], [], ls="", marker="o", ms=7, color=INK, label="pooled, 14 genes"),
         Line2D([], [], ls="", marker="o", ms=7, mfc="white", mec=CORE, label="pooled, core loop")]
    lay.shared_legend("A", h, [x.get_label() for x in h], side="bottom", ncol=2, fontsize=10, handletextpad=0.3,
                      columnspacing=1.0)
    fig = lay.finalize()
    pd.set_option("display.width", 200)
    print(b.round(3).to_string(index=False))
    print(pf.collisions(fig))
    pf.save(fig, OUT / "celltypes_health")


if __name__ == "__main__":
    main()
