"""Slide figure (results, R2.2): the organ sets the time; each cell class sets how strongly it expresses it.
Lymphocytes are the zero-gain end of that scale: no readable clock of their own.

A. Gain of each cell class relative to the endothelium of the same atlas (A293 Part 3, registered, H293.4 pass): ratio
   of the shared-rhythm amplitudes, log scale; coloured dots: one per atlas (mean over the family's compartments);
   black: random-effects pooled value with 95 % interval (pooled +- 1.96 SE). Lymphocytes: no gain is estimated; the
   arrow marks them at the zero-gain end (A281, A287, A294: no readable clock on the full clock-gene set).
B. Phase offsets between stromal-vascular classes of the same atlas (A293 Part 2; H293.3 FAIL, the order is not shown):
   coloured dots, human atlases (untimed pairwise offsets); black: circular mean over atlases with atlas bootstrap
   90 % interval; black open diamonds: mouse colon and aorta against known time (A273, A279). Positive = the first
   class earlier.
C. Do lymphocytes keep a clock of their own? Signed clock-gene coherence between donors (A281, `decompose.signed_coherence`,
   per-sample offset fitted out): lymphocytes, fibroblast-stromal cells of the same donors thinned to the lymphocytes'
   clock-gene counts (median of 20 thinnings), and the 5-95 % range of the lymphocytes' 200 expression-matched gene sets
   (grey bar; A294 recompute of A281's sets and seeds, reproducing its median and percentile), nine
   atlases (A281 H281.1 fail: 2 of 9 lymphocyte compartments coherent; H281.2 pass).

Data: exploratory/A293_organisation_levels/output/part3/{atlases,families}.csv;
exploratory/A293_organisation_levels/output/part2/{atlases,tests}.csv; exploratory/A281_lymphocyte_clock/output/pairs.csv;
exploratory/A294_lymphocyte_hormonal/output/pairs.csv. Output: presentation/figures/celltypes_health/celltypes_health.{svg,png}.

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
INK, MUTED = "#222222", pf.GREY
TAB = sns.color_palette().as_hex()
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


def gain(ax, at, fam):
    """Relative gain to endothelium per family: atlases (dots), pooled value and 95 % interval (black)."""
    fam = fam.set_index("family")
    lo_edge = 0.2
    for i, f in enumerate(GAIN_ORDER):
        col = FAMILY[f][1]
        if f == "endothelial":
            ax.scatter([1], [i], s=40, marker="D", color=col, edgecolor="white", lw=0.5, zorder=3)
            continue
        if f == "lymphocyte":
            ax.annotate("", xy=(lo_edge * 1.05, i), xytext=(0.5, i),
                        arrowprops=dict(arrowstyle="-|>", color=col, lw=1.6, mutation_scale=10))
            ax.text(lo_edge * 1.05, i - 0.3, "no readable clock", color=col, va="center", ha="left", fontsize=8)
            continue
        a = at[at.family == f]
        y = i + np.random.default_rng(i).uniform(-0.2, 0.2, len(a))
        ax.scatter(np.exp(a.g), y, s=18, color=col, alpha=0.6, edgecolor="none", zorder=1)
        r = fam.loc[f]
        ax.plot(np.exp([r.pooled - 1.96 * r.se_pooled, r.pooled + 1.96 * r.se_pooled]), [i, i], color=INK, lw=1.6,
                zorder=2)
        ax.scatter([np.exp(r.pooled)], [i], s=34, color=INK, zorder=3)
    ax.axvline(1, color=MUTED, lw=0.8, ls="--", zorder=0)
    ax.set_xscale("log")
    ax.set_xlim(lo_edge, 2.5)
    ax.set_xticks([0.25, 0.5, 1, 2], ["0.25", "0.5", "1", "2"])
    ax.minorticks_off()
    ax.set_yticks(range(len(GAIN_ORDER)), [FAMILY[f][0] for f in GAIN_ORDER])
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(len(GAIN_ORDER) - 0.5, -0.5)
    ax.set_xlabel("gain relative to endothelium")


def phase(ax, at, tests):
    """Untimed offsets between stromal-vascular classes: human atlases, their circular mean, mouse against known time."""
    tests = tests.set_index("class_pair")
    for i, (p, label) in enumerate(PAIRS.items()):
        a = at[at.class_pair == p]
        y = i + np.random.default_rng(i).uniform(-0.18, 0.18, len(a))
        ax.scatter(a.offset_h, y, s=18, color=CELL["stromal-vascular"], alpha=0.6, edgecolor="none", zorder=1)
        r = tests.loc[p]
        ax.plot([r.boot_lo90, r.boot_hi90], [i, i], color=INK, lw=1.6, zorder=2)
        ax.scatter([r.human_circular_mean_h], [i], s=34, color=INK, zorder=3)
        ax.scatter([r.mouse_colon_timed_h, r.mouse_aorta_timed_h], [i + 0.3] * 2, s=26, marker="D", facecolor="white",
                   edgecolor=INK, lw=0.9, zorder=3)
    ax.axvline(0, color=MUTED, lw=0.8, ls="--", zorder=0)
    ax.set_xlim(-1.5, 3.5)
    ax.set_ylim(len(PAIRS) - 0.4, -0.5)
    ax.set_yticks(range(len(PAIRS)), list(PAIRS.values()))
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("lead of the first class [h]")
    h = [Line2D([], [], ls="", marker="o", ms=5, color=CELL["stromal-vascular"], alpha=0.6, label="human atlas"),
         Line2D([], [], ls="", marker="D", ms=5, mfc="white", mec=INK, label="mouse, known time")]
    ax.legend(handles=h, loc="upper right", frameon=False, fontsize=8, handletextpad=0.1)


def own_clock(ax, p, rng_):
    p = p.merge(rng_[["atlas", "matched_q05", "matched_q95"]], on="atlas").sort_values("C_stromal_thinned_median",
                                                                                         ascending=False)
    for i, r in enumerate(p.itertuples()):
        ax.plot([r.matched_q05, r.matched_q95], [i, i], color=MUTED, lw=6, alpha=0.3, solid_capstyle="butt", zorder=1)
        ax.scatter([r.C_stromal_thinned_median], [i], s=34, color=CELL["stromal-vascular"], edgecolor="white", lw=0.5,
                   zorder=3)
        ax.scatter([r.C_lymph], [i], s=34, color=CELL["lymphocyte"], edgecolor="white", lw=0.5, zorder=3)
    ax.set_yticks(range(len(p)), [ATLAS[a] for a in p.atlas])
    ax.tick_params(axis="y", length=0)
    ax.invert_yaxis()
    ax.axvline(0, color=MUTED, lw=0.8, ls="--")
    ax.set_xlim(-0.8, 1)
    ax.set_xlabel("own clock: coherence\nbetween donors")
    h = [Line2D([], [], ls="", marker="o", ms=6, color=CELL["stromal-vascular"], label="stromal-vascular"),
         Line2D([], [], ls="", marker="o", ms=6, color=CELL["lymphocyte"], label="lymphocyte"),
         Line2D([], [], color=MUTED, lw=6, alpha=0.3, label="random genes")]
    ax.legend(handles=h, loc="upper center", bbox_to_anchor=(0.45, -0.24), ncol=3, frameon=False, fontsize=8, columnspacing=0.8,
              handletextpad=0.4)


def main():
    pf.use_style("slide")
    OUT.mkdir(parents=True, exist_ok=True)
    lay = pf.Layout.slide("ABC", kind="wide", col_ratios=[1.1, 1, 1])
    p3, p2 = EXP / "A293_organisation_levels/output/part3", EXP / "A293_organisation_levels/output/part2"
    gain(lay.axes["A"], pd.read_csv(p3 / "atlases.csv"), pd.read_csv(p3 / "families.csv"))
    phase(lay.axes["B"], pd.read_csv(p2 / "atlases.csv"), pd.read_csv(p2 / "tests.csv"))
    own_clock(lay.axes["C"], pd.read_csv(EXP / "A281_lymphocyte_clock/output/pairs.csv"),
              pd.read_csv(EXP / "A294_lymphocyte_hormonal/output/pairs.csv"))
    fig = lay.finalize()
    print(pf.collisions(fig))
    pf.save(fig, OUT / "celltypes_health")


if __name__ == "__main__":
    main()
