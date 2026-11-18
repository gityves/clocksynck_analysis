"""Slide figure (results, R3.1): what holds the organisation. The organ arrangement survives the loss of daily cues but
not the loss of the brain; a peripheral clock keeps its amplitude when its main cue is removed.

A. GTEx: does a death group keep the organ arrangement of fast deaths? Fisher-Lee correlation of the organ-pair phase
   differences of each group with those of fast deaths (12 body tissues, pairs with >= 20 shared donors), for the clock
   genes (blue) and 200 expression-matched and 200 stress-responsive gene sets (grey dots, same donors and pairs)
   (A293 Part 1; slow deaths H293.1 and ventilated living-brain donors H293.2 registered, Hardy 3 and brain-dead organ
   donors descriptive).
B. Mouse liver: amplitude of the clock state relative to ad libitum feeding (A297, registered, H297.1 fail above the
   interval), 90 % library bootstrap interval. Feeding inverted (Manella day feeding; Xin day- vs night-restricted) is
   the calibration; grey open: the amplitude predicted for removal of the feeding rhythm if cues add, (1 + rho) / 2;
   black: observed removal (Greenwell arrhythmic feeding) and the night-restricted control. Right: phase shift [h].

Data: exploratory/A293_organisation_levels/output/part1/{groups,sets}.csv; exploratory/A297_feeding_model/output/tests.csv.
GTEx aggregates only. Output: presentation/figures/holds_organisation/holds_organisation.{svg,png}.

    python scripts/holds_organisation.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns

import paperfig as pf

PAPER = Path(__file__).resolve().parents[2]
EXP = PAPER / "exploratory"
OUT = PAPER / "presentation" / "figures" / "holds_organisation"
INK, MUTED = "#222222", pf.GREY
CS = sns.color_palette().as_hex()[0]
GROUPS = {"slow_hardy4": "slow deaths", "ventilated_postmortem": "ventilated,\nbrain alive",
          "hardy3": "intermediate\ndeaths", "organ_donor": "brain-dead\norgan donors"}
CONTRASTS = [  # (label in tests.csv starts with, display label, kind)
    ("Manella DF", "inverted\n(Manella)", "obs"), ("Xin DRF", "inverted\n(Xin)", "obs"),
    ("prediction from Manella", "removed, additive\nprediction (Manella)", "pred"),
    ("prediction from Xin", "removed, additive\nprediction (Xin)", "pred"),
    ("observed Greenwell arrhythmic", "removed,\nobserved", "obs"),
    ("control Greenwell night", "night-restricted\n(control)", "obs")]


def arrangement(ax, g, sets):
    g = g.set_index("group")
    for i, k in enumerate(GROUPS):
        for j, fam in enumerate(("matched", "responsive")):
            v = sets[(sets.group == k) & (sets.family == fam)].rho.to_numpy()
            y = i + (j - 0.5) * 0.3 + np.random.default_rng(10 * i + j).uniform(-0.08, 0.08, len(v))
            ax.scatter(v, y, s=5, color=MUTED, alpha=0.45, lw=0, zorder=1)
        ax.scatter([g.loc[k, "rho_clock"]], [i], s=48, color=CS, zorder=3)
    ax.axvline(0, color=MUTED, lw=0.8, ls="--", zorder=0)
    ax.set_xlim(-0.3, 1)
    ax.set_yticks(range(len(GROUPS)), list(GROUPS.values()))
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(len(GROUPS) - 0.5, -0.5)
    ax.set_xlabel("organ arrangement shared\nwith fast deaths")
    ax.text(0.03, len(GROUPS) - 0.62, "grey, upper: random genes; lower: stress-responsive genes", fontsize=8, color=MUTED, va="bottom")


def liver(ax, t):
    for i, (key, label, kind) in enumerate(CONTRASTS):
        r = t[t.contrast.str.startswith(key)].iloc[0]
        col = MUTED if kind == "pred" else INK
        ax.plot([r.amp_lo90, r.amp_hi90], [i, i], color=col, lw=1.6, zorder=1)
        ax.scatter([r.amp_ratio], [i], s=36, color="white" if kind == "pred" else col, edgecolor=col, lw=1.2, zorder=2)
        ax.text(1.32, i, f"{r.shift_h:+.1f} h", va="center", ha="left", fontsize=8.5, color=col)
    ax.text(1.32, -0.75, "shift", va="center", ha="left", fontsize=8.5, color=MUTED)
    ax.axvline(1, color=MUTED, lw=0.8, ls="--", zorder=0)
    ax.set_xlim(0, 1.3)
    ax.set_yticks(range(len(CONTRASTS)), [c[1] for c in CONTRASTS])
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(len(CONTRASTS) - 0.5, -1)
    ax.set_xlabel("liver clock amplitude,\nrelative to ad libitum")


def main():
    pf.use_style("slide")
    OUT.mkdir(parents=True, exist_ok=True)
    p1 = EXP / "A293_organisation_levels/output/part1"
    lay = pf.Layout.slide("AB", kind="wide", col_ratios=[1, 1])
    arrangement(lay.axes["A"], pd.read_csv(p1 / "groups.csv"), pd.read_csv(p1 / "sets.csv"))
    liver(lay.axes["B"], pd.read_csv(EXP / "A297_feeding_model/output/tests.csv"))
    fig = lay.finalize()
    print(pf.collisions(fig))
    pf.save(fig, OUT / "holds_organisation")


if __name__ == "__main__":
    main()
