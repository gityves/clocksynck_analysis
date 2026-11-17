"""Slide figure (results, R2.1): the organs of a person share one clock state; the brain is less coupled to it.

A. Coupling of every GTEx tissue to the body's clock state (peripheral-anchor host, peripheral tissues read leave-one-out;
   A263 R2 fit, registered readout, here with a donor bootstrap: `organs_compute.py`): one dot per tissue with its 95 %
   interval, sorted within periphery (26 tissues) and brain (13 regions), coloured by organ system; black: group medians;
   grey dashed: 95th percentile of the same group median over 200 stress-responsive gene sets (A263 sets, post hoc
   comparison). Pituitary, in neither registered group, is not shown (coupling in the caption).
D. How far a tissue's clock departs from its person's body time (A259, registered): SD across donors [h] with its 90 %
   interval, 11 body tissues. An upper bound: the departure is not above that of the companion gene sets (H259.2).

Data: presentation/output/organs/couplings.csv; exploratory/A263_companions/output/sets_human.csv;
exploratory/A259_tissue_deviation/output/score/targets.csv. GTEx aggregates only (one value per tissue).
Output: presentation/figures/organs_health/organs_health.{svg,png}.

    python scripts/organs_health.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import paperfig as pf

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _systems import COLOR, gtex_system, legend_in  # noqa: E402

PAPER = Path(__file__).resolve().parents[2]
EXP = PAPER / "exploratory"
OUT = PAPER / "presentation" / "figures" / "organs_health"
INK, MUTED = "#222222", pf.GREY
SHORT = {"Adipose - Visceral (Omentum)": "adipose, visceral", "Adipose - Subcutaneous": "adipose, subcut.",
         "Skin - Not Sun Exposed (Suprapubic)": "skin, not sun-exp.", "Skin - Sun Exposed (Lower leg)": "skin, sun-exposed",
         "Heart - Atrial Appendage": "heart, atrium", "Esophagus - Mucosa": "oesophagus mucosa",
         "Muscle - Skeletal": "skeletal muscle", "Artery - Tibial": "tibial artery", "Nerve - Tibial": "tibial nerve",
         "Thyroid": "thyroid", "Lung": "lung"}


def couplings(ax, t, sets):
    x0 = 0
    ticks = []
    for g, col in (("periphery", "median_periph"), ("brain", "median_brain")):
        s = t[t.group == g].sort_values("coupling", ascending=False)
        x = x0 + np.arange(len(s))
        c = [COLOR[gtex_system(k)] for k in s.index]
        ax.vlines(x, s.coupling_lo, s.coupling_hi, colors=c, lw=1, alpha=0.6, zorder=1)
        ax.scatter(x, s.coupling, s=22, c=c, edgecolor="white", lw=0.4, zorder=2)
        ax.plot([x[0] - 0.4, x[-1] + 0.4], [s.coupling.median()] * 2, color=INK, lw=1.6, zorder=3)
        q = sets.loc[sets.kind == "responsive", col].quantile(0.95)
        ax.plot([x[0] - 0.4, x[-1] + 0.4], [q] * 2, color=MUTED, lw=1.4, ls="--", zorder=0)
        ticks.append((x.mean(), g))
        x0 += len(s) + 2
    ax.annotate("other gene sets\n(95th percentile)", (x[0] - 0.4, q), xytext=(0, -4), textcoords="offset points",
                ha="left", va="top", fontsize=8, color=MUTED)
    ax.set_xticks([p for p, _ in ticks], [lab for _, lab in ticks])
    ax.tick_params(axis="x", length=0)
    ax.set_xlim(-1.5, x0 - 1)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("coupling to the\nbody's clock state")
    legend_in(ax, fontsize=8, loc="lower left")


def departure(ax, d):
    d = d.sort_values("hours")
    y = np.arange(len(d))
    c = [COLOR[gtex_system(k)] for k in d.target]
    ax.hlines(y, d.hours_lo, d.hours_hi, colors=c, lw=1.6)
    ax.scatter(d.hours, y, s=24, c=c, zorder=3)
    ax.set_yticks(y, [SHORT.get(k, k) for k in d.target], fontsize=9)
    ax.tick_params(axis="y", length=0)
    ax.set_xlim(0, 4)
    ax.set_xlabel("departure from\nbody time, SD [h]")


def main():
    pf.use_style("slide")
    OUT.mkdir(parents=True, exist_ok=True)
    t = pd.read_csv(PAPER / "presentation/output/organs/couplings.csv", index_col="tissue")
    sets = pd.read_csv(EXP / "A263_companions/output/sets_human.csv")
    dep = pd.read_csv(EXP / "A259_tissue_deviation/output/score/targets.csv")
    lay = pf.Layout.slide("AD", kind="wide", col_ratios=[1.7, 1])
    couplings(lay.axes["A"], t, sets)
    departure(lay.axes["D"], dep)
    fig = lay.finalize()
    print(t.groupby("group").coupling.median().round(3).to_dict())
    print(pf.collisions(fig))
    pf.save(fig, OUT / "organs_health")


if __name__ == "__main__":
    main()
