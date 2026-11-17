"""Slide figure (results, R2.2): inside tissues, stromal-vascular cells and macrophages share the tissue's clock state;
lymphocytes do not, and carry no readable clock of their own.

Readout (anchor-free, `clocksync.pairwise`): for two cell classes of the same donors, |r| = the complex correlation over
donors of their clock-plane states; a class's sharing = the median |r| over its pairs with the non-immune classes of
the tissue (X). Companions: 200 expression-matched and 200 stress-responsive gene sets per atlas; calibrated p from
5000 placebo draws (A283 form).
A. Every cell class of the ten healthy body atlases of A272 and A274 (one dot per compartment): clock state shared with
   the tissue's non-immune cells, median |r| over its pairs with them, by cell family; grey open circle: the same
   compartment's expression-matched gene sets (median over its pairs). Post hoc, descriptive (from the scored pair
   tables; the registered statistics are in B).
B. Do lymphocytes keep a clock of their own? Signed clock-gene coherence between donors (A281, `decompose.signed_coherence`,
   per-sample offset fitted out): lymphocytes, fibroblast-stromal cells of the same donors thinned to the lymphocytes'
   clock-gene counts (median of 20 thinnings), and the lymphocytes' expression-matched gene sets (grey, median), nine
   atlases (A281 H281.1 fail: 2 of 9 lymphocyte compartments coherent; H281.2 pass).

Data: exploratory/A272_backbone_pairwise/output/pairs.csv; exploratory/A274_celltype_map/output/pairs.csv;
exploratory/A281_lymphocyte_clock/output/pairs.csv. Output: presentation/figures/celltypes_health/celltypes_health.{svg,png}.

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
BODY = ("lung_hlca", "heart", "organ_breast", "organ_muscle", "organ_colon", "organ_optic_nerve", "organ_eye_anterior",
        "kidney_sc", "kong_ileum", "kong_colon")
ATLAS = {"lung_hlca": "lung", "heart": "heart", "organ_breast": "breast", "organ_colon": "colon (Pelka)",
         "organ_muscle": "muscle", "kidney_sc": "kidney", "kong_ileum": "ileum (Kong)", "kong_colon": "colon (Kong)", "organ_eye_anterior": "anterior eye",
         "organ_optic_nerve": "optic nerve", "ibd_ileum": "ileum (IBD)", "ibd_rectum": "rectum (IBD)"}


def compartments():
    """Per compartment of the body atlases: median clock |r| and median matched-set |r| over its pairs with non-immune
    classes."""
    P = pd.concat([pd.read_csv(EXP / f"{d}/output/pairs.csv") for d in ("A272_backbone_pairwise", "A274_celltype_map")])
    P = P[P.atlas.isin(BODY)]
    cols = ["atlas", "c", "f", "g", "r_abs", "matched_median"]
    L = pd.concat([P.rename(columns={"a": "c", "family_a": "f", "family_b": "g"})[cols],
                   P.rename(columns={"b": "c", "family_b": "f", "family_a": "g"})[cols]])
    L = L[~L.g.isin(["lymphocyte", "myeloid"]) & L.f.isin(FAMILY)]
    return L.groupby(["atlas", "c", "f"]).agg(clock=("r_abs", "median"), matched=("matched_median", "median")).reset_index()


def families(ax, C):
    order = C.groupby("f").clock.median().sort_values(ascending=False).index
    for i, f in enumerate(order):
        c = C[C.f == f]
        y = i + np.random.default_rng(i).uniform(-0.22, 0.22, len(c))
        ax.scatter(c.matched, y, s=14, facecolor="white", edgecolor=MUTED, lw=0.8, zorder=1)
        ax.scatter(c.clock, y, s=26, color=FAMILY[f][1], edgecolor="white", lw=0.4, zorder=2)
        ax.plot([c.clock.median()] * 2, [i - 0.32, i + 0.32], color=INK, lw=1.8, zorder=3)
    ax.set_yticks(range(len(order)), [FAMILY[f][0] for f in order])
    ax.tick_params(axis="y", length=0)
    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.set_xlabel("clock state shared with\nthe tissue's other cells")
    h = [Line2D([], [], ls="", marker="o", ms=5, mfc="white", mec=MUTED, label="random genes")]
    ax.legend(handles=h, loc="lower right", frameon=False, fontsize=8, handletextpad=0.1)


def own_clock(ax, p):
    p = p.sort_values("C_stromal_thinned_median", ascending=False)
    for i, r in enumerate(p.itertuples()):
        ax.plot([r.C_lymph, r.C_stromal_thinned_median], [i, i], color=MUTED, lw=0.8, zorder=1)
        ax.scatter([r.matched_C_median], [i], s=22, facecolor="white", edgecolor=MUTED, lw=1.2, zorder=2)
        ax.scatter([r.C_stromal_thinned_median], [i], s=34, color=CELL["stromal-vascular"], edgecolor="white", lw=0.5,
                   zorder=3)
        ax.scatter([r.C_lymph], [i], s=34, color=CELL["lymphocyte"], edgecolor="white", lw=0.5, zorder=3)
    ax.set_yticks(range(len(p)), [ATLAS[a] for a in p.atlas])
    ax.tick_params(axis="y", length=0)
    ax.invert_yaxis()
    ax.axvline(0, color=MUTED, lw=0.8, ls="--")
    ax.set_xlim(-0.6, 1)
    ax.set_xlabel("own clock: coherence\nbetween donors")
    h = [Line2D([], [], ls="", marker="o", ms=6, color=CELL["stromal-vascular"], label="stromal-vascular,\nsame depth"),
         Line2D([], [], ls="", marker="o", ms=6, color=CELL["lymphocyte"], label="lymphocyte"),
         Line2D([], [], ls="", marker="o", ms=5, mfc="white", mec=MUTED, label="random genes")]
    ax.legend(handles=h, loc="lower right", frameon=False, fontsize=8, handletextpad=0.1)


def main():
    pf.use_style("slide")
    OUT.mkdir(parents=True, exist_ok=True)
    lay = pf.Layout.slide("AB", kind="wide", col_ratios=[1.1, 1])
    families(lay.axes["A"], compartments())
    own_clock(lay.axes["B"], pd.read_csv(EXP / "A281_lymphocyte_clock/output/pairs.csv"))
    fig = lay.finalize()
    print(pf.collisions(fig))
    pf.save(fig, OUT / "celltypes_health")


if __name__ == "__main__":
    main()
