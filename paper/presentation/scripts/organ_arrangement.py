"""Slide figure (results, adult health, last slide of section 3): the organ arrangement is internal: it holds when the body's clock no longer
reads the hour of death.

Rows: GTEx fast deaths (Hardy 1-2; reference), slow deaths (Hardy 4) and intermediate deaths (Hardy 3).
A. Correlation of each group's organ-pair phase differences with those of fast deaths (Fisher-Lee rho, 12 body tissues;
   A293 part 1): clock genes (blue bar); lines: 95th percentile of 200 expression-matched (light grey) and 200
   stress-responsive (dark grey) gene sets, same donors and pairs. Fast deaths are the reference (1 by construction,
   outlined). Filled: above both lines. Open diamond and thin line: median and 5-95 % of the value expected from
   sampling noise alone, n fast deaths (the group's donor count) against the other fast deaths (post hoc, 200 draws;
   `A374_order_interval/scripts/posthoc_noise_reference.py`).
B. Alignment of the body's clock (host phase from the clock genes, fitted without time) to the recorded hour of death,
   R = |mean exp(i (host phase - w t_death))|, all donors of each group (A27 readout).
C. How far the organs can have moved: 90 % compatibility interval of the organ scatter (SD of each organ's shift, h)
   given the group's arrangement correlation (A374, 14 genes; core loop in the speaker notes); dashed: the size of the
   arrangement itself (median organ-pair offset in fast deaths).

Status: A293 part 1 registered, NOT BLIND (clock rho known from A109 / A132; companions new): H293.1 (slow) pass;
intermediate deaths descriptive. Panel B post hoc readout (A27). Panel C A374 registered (intervals blind). Ventilated and brain-dead groups are shown in sections
4 and 7 (`holds_organisation.py`). GTEx aggregates only (one value per group).
Data: exploratory/A293_organisation_levels/output/part1/groups.csv; results/figure_3/output/cues_lost/walltime_hosts.csv;
exploratory/A374_order_interval/output/{intervals,noise_reference}.csv.
Output: presentation/figures/organ_arrangement/organ_arrangement.{svg,png} (half).

    python scripts/organ_arrangement.py
"""
from pathlib import Path

import pandas as pd
import seaborn as sns
from matplotlib import patheffects
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

import paperfig as pf

PAPER = Path(__file__).resolve().parents[2]
EXP = PAPER / "exploratory"
OUT = PAPER / "presentation" / "figures" / "organ_arrangement"
CS = sns.color_palette().as_hex()[0]
MATCHED, RESPONSIVE = "#BBBBBB", "#6E6E6E"
EDGE = [patheffects.Stroke(linewidth=3.6, foreground="white"), patheffects.Normal()]
ROWS = {"fast": "fast\ndeaths", "slow_hardy4": "slow\ndeaths", "hardy3": "intermediate\ndeaths"}
WALL = {"fast": "fast", "slow_hardy4": "hardy4", "hardy3": "hardy3"}
# group -> row of the A374 interval table; noise reference: the row of the same donor count
IVL = {"slow_hardy4": "slow deaths (Hardy 4)", "hardy3": "Hardy 3 (descriptive)",
       "ventilated_postmortem": "ventilated, living brain", "organ_donor": "brain-dead organ donors"}
NOISE = {**IVL, "organ_donor": "ventilated, living brain"}  # brain-dead rho is read at 79 donors (A293)
H = 0.55


def rows(ax, R=ROWS):
    ax.set_yticks(range(len(R)), list(R.values()))
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(len(R) - 0.5, -0.5)


def qline(ax, y, x, col):
    ax.plot([x, x], [y - H / 2, y + H / 2], color=col, lw=2.2, solid_capstyle="butt", zorder=3,
            path_effects=EDGE if col == RESPONSIVE else None)


def arrangement(ax, g, nz, R=ROWS):
    g = g.set_index("group")
    nz = nz[nz.gene_set == "14 genes"].set_index("group")
    ax.barh(0, 1, height=H, facecolor="white", edgecolor=CS, lw=1.4, ls="--")
    ax.text(0.5, 0, "reference", ha="center", va="center", fontsize=10, color=CS)
    out = []
    for i, k in enumerate(R):
        if k == "fast":
            continue
        r = g.loc[k]
        ok = r.rho_clock > max(r.matched_q95, r.responsive_q95)
        ax.barh(i, r.rho_clock, height=H, color=CS if ok else "white", edgecolor=CS, lw=1.4, zorder=2)
        qline(ax, i, r.matched_q95, MATCHED)
        qline(ax, i, r.responsive_q95, RESPONSIVE)
        e = nz.loc[NOISE[k]]
        ax.plot([e.rho_q05, e.rho_q95], [i + H / 2 + 0.12] * 2, color=CS, lw=1.2, zorder=4)
        ax.plot(e.rho_median, i + H / 2 + 0.12, ls="", marker="D", ms=6, mfc="white", mec=CS, mew=1.4, zorder=5)
        out.append(dict(group=k, n=r.n, rho=r.rho_clock, matched_q95=r.matched_q95, responsive_q95=r.responsive_q95,
                        noise_median=e.rho_median, noise_q05=e.rho_q05, noise_q95=e.rho_q95, test=r.test))
    ax.set_xlim(0, 1.05)
    rows(ax, R)
    ax.set_xlabel("organ arrangement\nshared with fast deaths")
    return pd.DataFrame(out)


def external(ax, w, R=ROWS, wall=WALL):
    w = w.set_index("cohort")
    out = []
    for i, k in enumerate(R):
        r = w.loc[wall[k]]
        ax.barh(i, r.R, height=H, color=CS)
        out.append(dict(group=k, donors=r.donors, R=r.R, p_R=r.p_R))
    ax.set_xlim(0, 1.05)
    rows(ax, R)
    ax.set_yticklabels([])
    ax.set_xlabel("body clock aligned to\nthe hour of death, R")
    return pd.DataFrame(out)


def interval(ax, iv, R=ROWS, xmax=2.5):
    iv = iv[iv.gene_set == "14 genes"].set_index("group")
    scale = float(iv.median_abs_pair_offset_fast_h.iloc[0])
    ax.axvline(scale, color=pf.GREY, lw=1, ls="--", zorder=0)
    ax.text(scale, -0.42, " size of the\n arrangement", ha="left", va="top", fontsize=9, color=pf.GREY)
    out = []
    for i, k in enumerate(R):
        if k == "fast":
            continue
        r = iv.loc[IVL[k]]
        ax.plot([r.lower_h, r.upper_h], [i, i], color=CS, lw=4, solid_capstyle="butt")
        if r.upper_at_grid_edge:  # upper bound not reached on the simulated grid
            ax.annotate("", (r.upper_h + 0.5, i), xytext=(r.upper_h - 0.05, i),
                        arrowprops=dict(arrowstyle="-|>", color=CS, lw=1.5, mutation_scale=12))
        out.append(dict(group=k, lower_h=r.lower_h, upper_h=r.upper_h, scale_h=scale))
    ax.set_xlim(0, xmax)
    rows(ax, R)
    ax.set_yticklabels([])
    ax.set_xlabel("how far the organs\ncan have moved [h]")
    return pd.DataFrame(out)


def main():
    pf.use_style("slide")
    OUT.mkdir(parents=True, exist_ok=True)
    lay = pf.Layout.slide("ABC", kind="half", height=100, col_ratios=[1, 0.8, 0.8], align="all")
    a = arrangement(lay.axes["A"], pd.read_csv(EXP / "A293_organisation_levels/output/part1/groups.csv"),
                    pd.read_csv(EXP / "A374_order_interval/output/noise_reference.csv"))
    b = external(lay.axes["B"], pd.read_csv(PAPER / "results/figure_3/output/cues_lost/walltime_hosts.csv"))
    c = interval(lay.axes["C"], pd.read_csv(EXP / "A374_order_interval/output/intervals.csv"))
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
    pf.save(fig, OUT / "organ_arrangement")


if __name__ == "__main__":
    main()
