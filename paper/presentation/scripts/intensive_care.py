"""Slide figure (results, drives, 4.2): the organ arrangement across death types, with intensive care and brain death;
same panels as slide 3.last (`organ_arrangement.py`), two rows added.

Rows: GTEx fast deaths (Hardy 1-2; reference), slow deaths (Hardy 4), intermediate deaths (Hardy 3), ventilated
post-mortem donors with a living brain, brain-dead organ donors; 12 body tissues.
A. Correlation of each group's organ-pair phase differences with those of fast deaths (Fisher-Lee rho; A293 part 1):
   clock genes (blue bar, filled when above both lines); lines: 95th percentile of 200 expression-matched (light grey)
   and 200 stress-responsive (dark grey) gene sets. Diamond and thin line: median and 5-95 % expected from sampling
   noise alone, n fast deaths against the others (A374 post hoc note; brain death at 79 donors, as its rho).
B. Alignment of the body's clock to the recorded hour of death, R, all donors of each group (fast, slow, intermediate:
   A27 readout; ventilated with a living brain and brain death: A314, descriptive split).
C. 90 % compatibility interval of the organ scatter given the group's arrangement correlation (A374, 14 genes; arrow:
   upper bound beyond the simulated grid); dashed: the size of the arrangement.

Status: A293 part 1 registered (H293.1 slow, H293.2 ventilated: pass; intermediate and brain death descriptive); A374
registered; panel B post hoc readout (A27) and descriptive (A314). GTEx aggregates only (one value per group).
Data: exploratory/A293_organisation_levels/output/part1/groups.csv; results/figure_3/output/cues_lost/walltime_hosts.csv;
exploratory/A314_vent_walltime/output/groups.csv; exploratory/A374_order_interval/output/{intervals,noise_reference}.csv.
Output: presentation/figures/intensive_care/intensive_care.{svg,png} (wide).

    python scripts/intensive_care.py
"""
from pathlib import Path

import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

import paperfig as pf
from organ_arrangement import CS, EDGE, MATCHED, RESPONSIVE, ROWS, WALL, arrangement, external, interval

PAPER = Path(__file__).resolve().parents[2]
EXP = PAPER / "exploratory"
OUT = PAPER / "presentation" / "figures" / "intensive_care"
R = {**ROWS, "ventilated_postmortem": "ventilated,\nliving brain", "organ_donor": "brain\ndeath"}
W = {**WALL, "ventilated_postmortem": "ventilated_postmortem", "organ_donor": "organ_donor"}


def wall():
    """A27 groups and the two A314 groups (all donors) in one table keyed by cohort."""
    a = pd.read_csv(PAPER / "results/figure_3/output/cues_lost/walltime_hosts.csv")
    b = pd.read_csv(EXP / "A314_vent_walltime/output/groups.csv").rename(columns={"group": "cohort"})
    return pd.concat([a[["cohort", "donors", "R", "p_R"]], b[["cohort", "donors", "R", "p_R"]]])


def main():
    pf.use_style("slide")
    OUT.mkdir(parents=True, exist_ok=True)
    lay = pf.Layout.slide("ABC", kind="wide", col_ratios=[1, 0.8, 0.8], align="all")
    a = arrangement(lay.axes["A"], pd.read_csv(EXP / "A293_organisation_levels/output/part1/groups.csv"),
                    pd.read_csv(EXP / "A374_order_interval/output/noise_reference.csv"), R)
    b = external(lay.axes["B"], wall(), R, W)
    c = interval(lay.axes["C"], pd.read_csv(EXP / "A374_order_interval/output/intervals.csv"), R, xmax=6.6)
    h = [Patch(facecolor=CS, edgecolor=CS, label="clock genes"),
         Line2D([], [], color=MATCHED, lw=2.2, label="matched genes, 95th pct"),
         Line2D([], [], color=RESPONSIVE, lw=2.2, path_effects=EDGE, label="stress-responsive genes, 95th pct"),
         Line2D([], [], color=CS, lw=1.2, marker="D", ms=6, mfc="white", mec=CS, label="clock genes: expected if unchanged")]
    lay.shared_legend("ABC", h, [x.get_label() for x in h], side="bottom", ncol=2, fontsize=10, handletextpad=0.4,
                      columnspacing=1.0)
    fig = lay.finalize()
    pd.set_option("display.width", 200)
    print(a.round(3).to_string(index=False), b.round(4).to_string(index=False), c.round(3).to_string(index=False), sep="\n")
    print(pf.collisions(fig))
    pf.save(fig, OUT / "intensive_care")


if __name__ == "__main__":
    main()
