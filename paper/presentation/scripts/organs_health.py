"""Slide figure (results, R2.1): the organs of a person share one clock state; the brain is less coupled to it.

A. Coupling of every GTEx tissue to the body's clock state (peripheral-anchor host, peripheral tissues read leave-one-out;
   A263 R2 fit, registered readout, here with a donor bootstrap: `organs_compute.py`): one dot per tissue with its 95 %
   interval, sorted within periphery and brain, coloured by organ system (pituitary, excluded from both groups, shown apart).
B. Agreement of a block of tissues with the rest of the body (r*, disattenuated host correlation over donors, A142
   registered, 200 donor bootstraps): brain (13 regions) and digestive (7 tissues) against random peripheral blocks of the
   same size (grey: median and 5-95 % of 1000 blocks).
C. Is the brain's lower coupling specific to the clock genes? log(median coupling brain / periphery) for the clock genes
   (blue line) and for 200 expression-matched and 200 stress-responsive gene sets (A263 H263.1, registered: fail).
D. How far a tissue's clock departs from its person's body time (A259, registered): SD across donors [h] with its 90 %
   interval, 11 body tissues. An upper bound: the departure is not above that of the companion gene sets (H259.2).

Data: presentation/output/organs/couplings.csv; exploratory/A142_organ_blocks/output/{tests,reference_summary}.csv;
exploratory/A263_companions/output/sets_human.csv; exploratory/A259_tissue_deviation/output/score/targets.csv.
GTEx aggregates only (one value per tissue or block). Output: presentation/figures/organs_health/organs_health.{svg,png}.

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
CS = "#1f77b4"
SHORT = {"Adipose - Visceral (Omentum)": "adipose, visceral", "Adipose - Subcutaneous": "adipose, subcut.",
         "Skin - Not Sun Exposed (Suprapubic)": "skin, not sun-exp.", "Skin - Sun Exposed (Lower leg)": "skin, sun-exposed",
         "Heart - Atrial Appendage": "heart, atrium", "Esophagus - Mucosa": "oesophagus mucosa",
         "Muscle - Skeletal": "skeletal muscle", "Artery - Tibial": "tibial artery", "Nerve - Tibial": "tibial nerve",
         "Thyroid": "thyroid", "Lung": "lung"}


def couplings(ax, t):
    x0 = 0
    ticks = []
    for g, lab in (("periphery", "periphery"), ("excluded", "pituitary"), ("brain", "brain")):
        s = t[t.group == g].sort_values("coupling", ascending=False)
        x = x0 + np.arange(len(s))
        c = [COLOR[gtex_system(k)] for k in s.index]
        ax.vlines(x, s.coupling_lo, s.coupling_hi, colors=c, lw=1, alpha=0.6, zorder=1)
        ax.scatter(x, s.coupling, s=22, c=c, edgecolor="white", lw=0.4, zorder=2)
        if len(s) > 1:
            ax.plot([x[0] - 0.4, x[-1] + 0.4], [s.coupling.median()] * 2, color=INK, lw=1.6, zorder=3)
        ticks.append((x.mean(), lab))
        x0 += len(s) + 2
    ax.set_xticks([p for p, _ in ticks], [lab for _, lab in ticks])
    ax.tick_params(axis="x", length=0)
    ax.set_xlim(-1.5, x0 - 1)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("coupling to the\nbody's clock state")
    legend_in(ax, fontsize=8, loc="lower left")


def blocks(ax, tests, ref):
    for j, (b, k, lab) in enumerate((("brain", 13, "brain"), ("digestive", 7, "digestive"))):
        r = tests.set_index("block").loc[b]
        q = ref.set_index("k").loc[k]
        ax.plot([j - 0.12] * 2, [q.q05, q.q95], color=MUTED, lw=5, alpha=0.35, solid_capstyle="butt")
        ax.scatter([j - 0.12], [q.q50], marker="_", s=120, color=MUTED, lw=2)
        c = COLOR["brain, eye"] if b == "brain" else COLOR["digestive"]
        ax.errorbar([j + 0.12], [r.r_star], yerr=[[r.r_star - r.r_star_lo], [r.r_star_hi - r.r_star]], fmt="o", color=c,
                    ms=6, lw=1.6, capsize=0)
    ax.set_xticks([0, 1], ["brain", "digestive"])
    ax.tick_params(axis="x", length=0)
    ax.set_xlim(-0.5, 1.5)
    ax.set_ylim(0.4, 1.05)
    ax.set_ylabel("agreement with\nthe rest of the body")
    ax.annotate("random\nperipheral\nblocks", (0 - 0.12, ref.set_index("k").loc[13].q05), xytext=(0, -6),
                textcoords="offset points", ha="center", va="top", fontsize=8, color=MUTED)


def specificity(ax, sets):
    clock = sets.loc[sets.kind == "clock", "L"].iloc[0]
    for j, (k, lab) in enumerate((("matched", "expression-\nmatched"), ("responsive", "stress-\nresponsive"))):
        v = sets.loc[sets.kind == k, "L"].dropna().to_numpy()
        y = j + np.random.default_rng(j).uniform(-0.18, 0.18, len(v))
        ax.scatter(v, y, s=6, color=MUTED, alpha=0.5, lw=0)
    ax.axvline(clock, color=CS, lw=2.2)
    ax.annotate("clock genes", (clock, 1.45), xytext=(-4, 0), textcoords="offset points", fontsize=9, color=CS, va="center",
                ha="right")
    ax.set_yticks([0, 1], ["expression-\nmatched", "stress-\nresponsive"])
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(-0.5, 1.7)
    ax.set_xlim(-4, 2)
    ax.axvline(0, color=MUTED, lw=0.8, ls="--")
    ax.set_xlabel("brain / periphery\ncoupling [log ratio]")


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
    tests = pd.read_csv(EXP / "A142_organ_blocks/output/tests.csv")
    ref = pd.read_csv(EXP / "A142_organ_blocks/output/reference_summary.csv").query("plane == 'base'")
    sets = pd.read_csv(EXP / "A263_companions/output/sets_human.csv")
    dep = pd.read_csv(EXP / "A259_tissue_deviation/output/score/targets.csv")
    lay = pf.Layout.slide("ABCD", kind="wide", col_ratios=[1.5, 0.6, 0.8, 0.9])
    couplings(lay.axes["A"], t)
    blocks(lay.axes["B"], tests, ref)
    specificity(lay.axes["C"], sets)
    departure(lay.axes["D"], dep)
    fig = lay.finalize()
    print(t.groupby("group").coupling.median().round(3).to_dict())
    print(pf.collisions(fig))
    pf.save(fig, OUT / "organs_health")


if __name__ == "__main__":
    main()
