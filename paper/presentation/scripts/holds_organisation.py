"""Slide figures (results, R3.1 and R3.2): what holds the organisation.

R3.1 (`holds_organisation`), GTEx death types:
A. Does a death group keep the organ arrangement of fast deaths? Fisher-Lee correlation of the organ-pair phase
   differences of each group with those of fast deaths (12 body tissues, pairs with >= 20 shared donors): bars for the
   clock genes (blue) and the medians of 200 expression-matched (light grey) and 200 stress-responsive gene sets (dark
   grey), whisker to their 95th percentile, same donors and pairs; fast deaths are the reference (1 by construction)
   (A293 Part 1; slow deaths H293.1 and ventilated living-brain donors H293.2 registered, Hardy 3 and brain-dead organ
   donors descriptive). Rows in the same order in A and B.
B. Does the body's clock still read the external day? Alignment of the host phase (clock genes of 15 reference tissues,
   fitted without time) to the recorded hour of death, R = |mean exp(i (host phase - w t_death))|, per death type, with
   its permutation p (times shuffled over donors, 10,000): fast, slow and intermediate deaths from A27 (post hoc
   readout), ventilated donors with a living brain and brain-dead organ donors from A314 (descriptive, post hoc split).
R3.2 (`cue_removal`), designed removal of a cue in one organ:
A. Mouse liver (A297, registered, H297.1 fail above the interval): core-clock amplitude (blue, clock genes) with the
   feeding rhythm removed (Greenwell arrhythmic feeding), relative to ad libitum, 90 % library bootstrap interval; grey
   open: the amplitude expected if the clock state were the sum of its cues, (1 + rho) / 2, rho calibrated on an
   inverted feeding rhythm in two datasets (Manella, Xin).
B. Human muscle after three days of 12 h-shifted sleep, light, activity and meals (Wefers 2018, 11 men, biopsies at
   07:00 and 19:00): shift of the 07:00-to-19:00 change against the aligned condition, 90 % interval, for the core
   clock (14 clock genes, A257 H257.2) and the rhythmic outputs (1,587 reference-rhythmic genes, A261 H261.3); right:
   amplitude of the change relative to aligned (A257 H257.3, A261).

Data: exploratory/A293_organisation_levels/output/part1/{groups,sets}.csv; results/figure_3/output/cues_lost/walltime_hosts.csv;
exploratory/A314_vent_walltime/output/groups.csv;
exploratory/A297_feeding_model/output/tests.csv; exploratory/A257_nightshift_muscle/output/score_tests.csv;
exploratory/A261_nightshift_programs/output/programs.csv. GTEx aggregates only.
Output: presentation/figures/holds_organisation/{holds_organisation,cue_removal}.{svg,png}.

    python scripts/holds_organisation.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

import paperfig as pf

PAPER = Path(__file__).resolve().parents[2]
EXP = PAPER / "exploratory"
OUT = PAPER / "presentation" / "figures" / "holds_organisation"
INK, MUTED = "#222222", pf.GREY
CS = sns.color_palette().as_hex()[0]
FAMILY = {"matched": ("random genes", "#BBBBBB"), "responsive": ("stress-responsive genes", "#6E6E6E")}
ROWS = ["fast deaths", "slow deaths", "intermediate\ndeaths", "ventilated,\nbrain alive", "brain-dead\norgan donors"]
GROUPS = {"slow_hardy4": 1, "hardy3": 2, "ventilated_postmortem": 3, "organ_donor": 4}  # A293 group -> row
H = 0.26  # bar height in panel A


def rows(ax):
    ax.set_yticks(range(len(ROWS)), ROWS)
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(len(ROWS) - 0.45, -0.55)


def arrangement(ax, g):
    """Organ arrangement shared with fast deaths: clock genes and the medians of the two companion families, whisker to
    their 95th percentile; fast deaths are the reference (1 by construction)."""
    g = g.set_index("group")
    ax.barh(0, 1, height=H, facecolor="white", edgecolor=CS, lw=1.2)
    ax.text(0.5, 0, "reference", ha="center", va="center", fontsize=8, color=CS)
    for k, i in GROUPS.items():
        r = g.loc[k]
        ax.barh(i - H, r.rho_clock, height=H, color=CS)
        for j, (fam, (_, col)) in enumerate(FAMILY.items()):
            y = i + j * H
            ax.barh(y, r[f"{fam}_median"], height=H, color=col)
            ax.plot([r[f"{fam}_median"], r[f"{fam}_q95"]], [y, y], color=col, lw=1.2)
            ax.plot([r[f"{fam}_q95"]] * 2, [y - H / 3, y + H / 3], color=col, lw=1.2)
    ax.axvline(0, color=MUTED, lw=0.8, zorder=0)
    ax.set_xlim(-0.05, 1.05)
    rows(ax)
    ax.set_xlabel("organ arrangement shared\nwith fast deaths")
    ax.set_title("organ arrangement", fontsize=13, color=INK)
    h = [Patch(color=CS, label="clock genes")] + [Patch(color=c, label=lab) for lab, c in FAMILY.values()]
    ax.legend(handles=h, loc="upper center", bbox_to_anchor=(0.45, -0.24), ncol=3, frameon=False, fontsize=8,
              handlelength=1, columnspacing=0.8)


def external(ax, w, v):
    """Alignment of the body's clock to the recorded hour of death, one bar per row: fast, slow and intermediate deaths
    from A27, the two ventilated groups from A314 (all donors of each)."""
    w, v = w.set_index("cohort"), v.set_index("group")
    R = [w.loc["fast", "R"], w.loc["hardy4", "R"], w.loc["hardy3", "R"], v.loc["ventilated_postmortem", "R"],
         v.loc["organ_donor", "R"]]
    ax.barh(range(len(ROWS)), R, height=0.55, color=CS)
    ax.set_xlim(0, 1.05)
    rows(ax)
    ax.set_yticklabels([])
    ax.set_xlabel("body clock aligned to\nthe hour of death, R")
    ax.set_title("external time", fontsize=13, color=INK)


def liver(ax, t):
    obs = t[t.contrast.str.startswith("observed Greenwell")].iloc[0]
    ax.scatter([1], [0], s=40, color=CS, zorder=2)
    ax.plot([obs.amp_lo90, obs.amp_hi90], [1, 1], color=CS, lw=1.6, zorder=1)
    ax.scatter([obs.amp_ratio], [1], s=40, color=CS, zorder=2)
    ax.annotate(f"{obs.shift_h:+.1f} h", (obs.amp_hi90, 1), xytext=(6, 0), textcoords="offset points", va="center",
                fontsize=8.5, color=CS)
    for j, key in enumerate(("prediction from Manella", "prediction from Xin")):
        r = t[t.contrast.str.startswith(key)].iloc[0]
        y = 2 + (j - 0.5) * 0.25
        ax.plot([r.amp_lo90, r.amp_hi90], [y, y], color=MUTED, lw=1.6, zorder=1)
        ax.scatter([r.amp_ratio], [y], s=36, facecolor="white", edgecolor=MUTED, lw=1.2, zorder=2)
    ax.axvline(1, color=MUTED, lw=0.8, ls="--", zorder=0)
    ax.set_xlim(0, 1.2)
    ax.set_yticks([0, 1, 2], ["fed ad libitum", "feeding rhythm\nremoved", "expected if\ncues add"])
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(2.5, -0.5)
    ax.set_xlabel("liver clock amplitude\n(ad libitum = 1)")
    ax.set_title("mouse liver", fontsize=13, color=INK)


def muscle(ax, core, outputs):
    rows = [("core clock", core.s, core.lo90, core.hi90, core.R, CS),
            ("rhythmic\noutputs", outputs.shift_h, outputs.shift_lo, outputs.shift_hi, outputs.ratio, INK)]
    for i, (_, s, lo, hi, ratio, col) in enumerate(rows):
        ax.plot([lo, hi], [i, i], color=col, lw=1.6, zorder=1)
        ax.scatter([s], [i], s=40, color=col, zorder=2)
        ax.text(14.2, i, f"×{ratio:.2f}", va="center", ha="left", fontsize=8.5, color=col)
    ax.text(14.2, -0.3, "amplitude", va="center", ha="left", fontsize=8.5, color=MUTED)
    for x, lab in ((0, "central time"), (12, "behaviour")):
        ax.axvline(x, color=MUTED, lw=0.8, ls="--", zorder=0)
        ax.text(x + 0.3, -0.65, lab, ha="left", va="center", fontsize=8.5, color=MUTED)
    ax.set_xlim(-1, 14)
    ax.set_xticks([0, 3, 6, 9, 12])
    ax.set_yticks([0, 1], [r[0] for r in rows])
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(1.9, -0.9)
    ax.set_xlabel("shift after 12 h shift\nof behaviour [h]")
    ax.set_title("human muscle", fontsize=13, color=INK)


def main():
    pf.use_style("slide")
    OUT.mkdir(parents=True, exist_ok=True)
    p1 = EXP / "A293_organisation_levels/output/part1"
    lay = pf.Layout.slide("AB", kind="wide", col_ratios=[1.2, 1], align="all")
    arrangement(lay.axes["A"], pd.read_csv(p1 / "groups.csv"))
    external(lay.axes["B"], pd.read_csv(PAPER / "results/figure_3/output/cues_lost/walltime_hosts.csv"),
             pd.read_csv(EXP / "A314_vent_walltime/output/groups.csv"))
    fig = lay.finalize()
    print(pf.collisions(fig))
    pf.save(fig, OUT / "holds_organisation")

    lay = pf.Layout.slide("AB", kind="wide", col_ratios=[1, 1])
    liver(lay.axes["A"], pd.read_csv(EXP / "A297_feeding_model/output/tests.csv"))
    t257 = pd.read_csv(EXP / "A257_nightshift_muscle/output/score_tests.csv")
    core = t257[t257.test.str.startswith("H257.2")].iloc[0].copy()
    core["R"] = t257[t257.test.str.startswith("H257.3")].iloc[0].R
    prog = pd.read_csv(EXP / "A261_nightshift_programs/output/programs.csv").set_index("program")
    muscle(lay.axes["B"], core, prog.loc["ALL_RHYTHMIC"])
    fig = lay.finalize()
    print(pf.collisions(fig))
    pf.save(fig, OUT / "cue_removal")


if __name__ == "__main__":
    main()
