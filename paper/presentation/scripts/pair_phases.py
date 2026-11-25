"""Slide figure (results, drives, 4.2): organ pairs keep their phase differences in ventilated donors with a living
brain, not in brain death.

One dot per pair of the 12 body tissues (65 pairs, A317): phase difference between the two organs' clocks in the group
against fast deaths (+ = first organ peaks later). Left: ventilated post-mortem donors with a living brain (79); right:
brain-dead organ donors (circular mean over 20 draws of 79 donors, as the registered test). Dashed: identity
(arrangement kept). Bars: 95 % donor-bootstrap intervals of a group-sized estimate (200 draws), in both directions. Open
circles: pairs with the lung.

Status: display of the registered statistic (A293 H293.2 pass; A317 H317.1 pass); the display itself is post hoc (user,
2026-11-25). Pairs share organs, so they are not independent points. GTEx aggregates only (one row per tissue pair).
Data: presentation/output/pair_phases/pairs.csv (`pair_phases_compute.py`).
Output: presentation/figures/pair_phases/pair_phases.{svg,png} (half).

    python scripts/pair_phases.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D
from scipy.stats import pearsonr

import paperfig as pf

PAPER = Path(__file__).resolve().parents[2]
SRC = PAPER / "presentation" / "output" / "pair_phases" / "pairs.csv"
OUT = PAPER / "presentation" / "figures" / "pair_phases"
CS = sns.color_palette().as_hex()[0]
GROUPS = {"ventilated": "ventilated, living brain", "brain_death": "brain death"}
LIM = (-3.6, 4.2)


def panel(ax, t, g):
    lung = t.tissue_a.eq("Lung") | t.tissue_b.eq("Lung")
    ax.plot(LIM, LIM, color=pf.GREY, lw=0.8, ls="--", zorder=0)
    ax.axhline(0, color="#e6e6e6", lw=0.6, zorder=0)
    ax.axvline(0, color="#e6e6e6", lw=0.6, zorder=0)
    for r in t.itertuples():
        ax.plot([r.fast] * 2, [getattr(r, f"{g}_lo"), getattr(r, f"{g}_hi")], color=CS, lw=0.6, alpha=0.18, zorder=1)
        ax.plot([r.fast_lo, r.fast_hi], [getattr(r, g)] * 2, color=CS, lw=0.6, alpha=0.18, zorder=1)
    ax.scatter(t.fast[~lung], t[g][~lung], s=22, color=CS, edgecolor="white", lw=0.4, zorder=3)
    ax.scatter(t.fast[lung], t[g][lung], s=22, facecolor="white", edgecolor=CS, lw=1.1, zorder=3)
    ax.set(xlim=LIM, ylim=LIM, aspect="equal", title=GROUPS[g])
    ax.set_xlabel("phase difference, fast deaths [h]")
    return dict(group=g, pairs=len(t), pearson=pearsonr(t.fast, t[g])[0], slope=np.polyfit(t.fast, t[g], 1)[0],
                pearson_no_lung=pearsonr(t.fast[~lung], t[g][~lung])[0], sd=t[g].std(), sd_fast=t.fast.std(),
                median_halfwidth=float(np.median((t[f"{g}_hi"] - t[f"{g}_lo"]) / 2)))


def main():
    pf.use_style("slide")
    OUT.mkdir(parents=True, exist_ok=True)
    t = pd.read_csv(SRC)
    lay = pf.Layout.slide("AB", kind="half")
    rows = [panel(lay.axes[k], t, g) for g, k in zip(GROUPS, "AB")]
    lay.axes["A"].set_ylabel("phase difference in the group [h]")
    lay.axes["B"].tick_params(labelleft=False)
    h = [Line2D([], [], ls="", marker="o", ms=6, color=CS, mec="white", label="organ pair"),
         Line2D([], [], ls="", marker="o", ms=6, mfc="white", mec=CS, mew=1.1, label="pair with the lung"),
         Line2D([], [], color=pf.GREY, lw=0.8, ls="--", label="arrangement kept")]
    lay.shared_legend("AB", h, [x.get_label() for x in h], side="bottom", ncol=3, fontsize=10, handletextpad=0.3,
                      columnspacing=1.2)
    fig = lay.finalize()
    print(pd.DataFrame(rows).round(2).to_string(index=False))
    print("fast median half-width", round(float(np.median((t.fast_hi - t.fast_lo) / 2)), 2))
    print(pf.collisions(fig))
    pf.save(fig, OUT / "pair_phases")


if __name__ == "__main__":
    main()
