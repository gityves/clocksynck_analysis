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
C. Do lymphocytes keep a clock of their own? Per atlas, signed clock-gene coherence between donors (A281,
   `decompose.signed_coherence`, per-sample offset fitted out) of the lymphocytes (cyan) and of the same donors'
   fibroblast-stromal cells thinned to the lymphocytes' clock-gene counts (brown, median of 20 thinnings), joined by a
   line; grey tick: 95th percentile of the lymphocytes' 200 expression-matched gene sets (A294 recompute of A281's sets
   and seeds). Nine atlases (A281 H281.1 fail: 2 of 9 lymphocyte compartments coherent; H281.2 pass).
D. What "no readable clock" excludes (A326, registered; not blind for the human results): chance that the designs read
   a lymphocyte clock scaled to r times the fibroblasts' (A281, 14 genes, nine atlases; A321, core loop, four fresh
   atlases, blind); dotted: 0.8; cyan band: mouse colon T-cell / fibroblast amplitude ratio, point 0.246 to q95 0.404
   (L1 scout, undamaged arm).
A also shows the core-loop pooled gains (open blue, A350, blind; mural tau above the bound: "not conserved").

Data: exploratory/A293_organisation_levels/output/part3/{atlases,families}.csv;
exploratory/A293_organisation_levels/output/part2/{atlases,tests}.csv; exploratory/A281_lymphocyte_clock/output/pairs.csv;
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
    """Per atlas: coherence between donors of the lymphocytes' clock genes (cyan) and of the same donors' fibroblast-
    stromal cells thinned to the lymphocytes' counts (brown), joined by a line; grey tick: 95th percentile of the
    lymphocytes' 200 expression-matched gene sets (A294)."""
    p = p.merge(rng_[["atlas", "matched_q95"]], on="atlas").sort_values("C_stromal_thinned_median", ascending=False)
    for i, r in enumerate(p.itertuples()):
        ax.plot([r.C_lymph, r.C_stromal_thinned_median], [i, i], color="#DDDDDD", lw=1.2, zorder=1)
        ax.plot([r.matched_q95] * 2, [i - 0.32, i + 0.32], color=MATCHED, lw=2.2, solid_capstyle="butt", zorder=2)
        ax.scatter([r.C_stromal_thinned_median], [i], s=34, color=CELL["stromal-vascular"], edgecolor="white", lw=0.5,
                   zorder=3)
        ax.scatter([r.C_lymph], [i], s=34, color=CELL["lymphocyte"], edgecolor="white", lw=0.5, zorder=3)
    ax.set_yticks(range(len(p)), [ATLAS[a] for a in p.atlas])
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(len(p) - 0.5, -0.5)
    ax.axvline(0, color=MUTED, lw=0.8, ls="--", zorder=0)
    ax.set_xlim(-0.6, 1)
    ax.set_xlabel("clock-gene coherence\nbetween donors")
    return p


def bound(ax, dp, mouse):
    """What 'no readable clock' excludes (A326): chance that the lymphocyte designs would have read a clock of size r
    relative to the fibroblasts' (A281, 14 genes, nine atlases; A321, core loop, fresh atlases, blind); band: the
    mouse colon T-cell / fibroblast amplitude ratio, point to q95 (L1 scout)."""
    ax.axvspan(mouse.ratio_point, mouse.ratio_q95, color=CELL["lymphocyte"], alpha=0.15, lw=0, zorder=0)
    ax.text((mouse.ratio_point + mouse.ratio_q95) / 2, 0.45, "mouse\ncolon T", ha="center", va="center", fontsize=9,
            color=CELL["lymphocyte"])
    for d, ls, lab in (("A281", "-", "14 genes"), ("A321", "--", "core loop")):
        x = dp[dp.design == d].sort_values("r")
        ax.plot(x.r, x.design_power, color=INK, lw=1.4, ls=ls, marker="o", ms=4, label=lab, zorder=2)
    ax.axhline(0.8, color=MUTED, lw=0.6, ls=":", zorder=0)
    ax.set_xlim(0, 1.05)
    ax.set_ylim(-0.03, 1.03)
    ax.set_xlabel("lymphocyte clock size\n[fraction of fibroblasts']")
    ax.set_ylabel("chance it is read")
    ax.legend(loc="upper left", frameon=False, fontsize=9, handlelength=1.6)


def main():
    pf.use_style("slide")
    OUT.mkdir(parents=True, exist_ok=True)
    lay = pf.Layout.slide("ABCD", kind="wide", col_ratios=[1.1, 0.95, 0.95, 0.8])
    p3, p2 = EXP / "A293_organisation_levels/output/part3", EXP / "A293_organisation_levels/output/part2"
    gain(lay.axes["A"], pd.read_csv(p3 / "atlases.csv"), pd.read_csv(p3 / "families.csv"),
         pd.read_csv(EXP / "A350_class_gain_core/output/compare.csv"))
    phase(lay.axes["B"], pd.read_csv(p2 / "atlases.csv"), pd.read_csv(p2 / "tests.csv"))
    oc = own_clock(lay.axes["C"], pd.read_csv(EXP / "A281_lymphocyte_clock/output/pairs.csv"),
                   pd.read_csv(EXP / "A294_lymphocyte_hormonal/output/pairs.csv"))
    L1 = EXP / "L1_lymphocyte_power"
    mouse = pd.read_csv(L1 / "scout/output/ratios_vs_fibroblast.csv").query(
        "arm == 'C-DC' and compartment == 'T' and gene_set == 'all14'").iloc[0]
    dp = pd.read_csv(L1 / "output/design_power.csv")
    bound(lay.axes["D"], dp, mouse)
    h = [Line2D([], [], ls="", marker="o", ms=7, color=INK, label="pooled, 14 genes"),
         Line2D([], [], ls="", marker="o", ms=7, mfc="white", mec=CORE, label="pooled, core loop"),
         Line2D([], [], ls="", marker="o", ms=7, color=CELL["stromal-vascular"], label="fibroblasts, thinned"),
         Line2D([], [], ls="", marker="o", ms=7, color=CELL["lymphocyte"], label="lymphocytes"),
         Line2D([], [], color=MATCHED, lw=2.2, label="matched genes, 95th pct")]
    lay.shared_legend("ABCD", h, [x.get_label() for x in h], side="bottom", ncol=5, fontsize=10, handletextpad=0.3,
                      columnspacing=1.0)
    fig = lay.finalize()
    pd.set_option("display.width", 200)
    print(oc[["atlas", "C_lymph", "C_stromal_thinned_median", "matched_q95"]].round(3).to_string(index=False))
    print("mouse T/fibroblast ratio", round(mouse.ratio_point, 3), "q95", round(mouse.ratio_q95, 3))
    print(pf.collisions(fig))
    pf.save(fig, OUT / "celltypes_health")


if __name__ == "__main__":
    main()
