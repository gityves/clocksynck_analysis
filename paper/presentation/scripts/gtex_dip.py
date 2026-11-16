"""Slide figure (results, R1.2): CHIRAL's donor internal phase (DIP) against clocksync's untimed internal time, GTEx fast deaths.

Reads the aggregates written by `gtex_dip_compute.py` (no per-donor value: controlled access). Two panels, clocksync fitted on
A. all fast deaths (full day) and B. the fast deaths who died 08-18 h only (hospital hours). x: CHIRAL's DIP (learnt on all
GTEx donors, whole day); y: clocksync's internal time, aligned to the DIP by the best rotation and reflection. Dots and bars:
circular mean and circular SD of the internal time per 3 h bin of DIP (bins with < 5 donors omitted); dashed: identity.
C. median absolute error of each phase against the recorded time of death on the same donors (each aligned to it), with donor
bootstrap 95 % intervals. Colours: seaborn default palette (clocksync first, blue).

Output: paper/presentation/figures/gtex_dip/gtex_dip.{svg,png} (wide slide figure).

    python scripts/gtex_dip.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns

import paperfig as pf

SRC = Path(__file__).resolve().parents[1] / "output" / "gtex_dip"
OUT = Path(__file__).resolve().parents[1] / "figures" / "gtex_dip"
INK, MUTED = "#222222", pf.GREY
CS, CHIRAL = sns.color_palette().as_hex()[:2]
TITLE = {"full": "full day", "window": "hospital hours"}


def panel(ax, design, bins):
    b = bins[(bins.design == design)].dropna(subset=["it_mean_h"])
    x = (b.dip_lo + b.dip_hi) / 2
    y = b.it_mean_h + 24 * np.round((x - b.it_mean_h) / 24)  # unwrap towards the diagonal
    ax.plot([0, 24], [0, 24], color=MUTED, lw=1, ls="--", zorder=1)
    ax.errorbar(x, y, yerr=b.it_circ_sd_h, fmt="o", color=CS, ms=6, lw=1.6, capsize=0, zorder=3)
    ax.set_xlim(0, 24)
    ax.set_ylim(-6, 30)
    ax.set_xticks(range(0, 25, 6))
    ax.set_yticks(range(0, 25, 6))
    ax.set_aspect("equal")
    ax.set_xlabel("CHIRAL donor internal phase [h]")
    ax.set_ylabel("clocksync internal time [h]")
    ax.set_title(TITLE[design], fontsize=13, color=INK)


def tod_panel(ax, summ):
    """Median absolute error of each internal phase against the recorded time of death, same donors, bootstrap 95 % CI."""
    y = {"full": 1, "window": 0}
    for design, yy in y.items():
        s = summ.loc[design]
        for k, c, dy in (("dip", CHIRAL, 0.06), ("it", CS, -0.06)):
            m = s[f"{k}_mae_tod_h"]
            ax.errorbar([m], [yy + dy], xerr=[[m - s[f"{k}_mae_tod_lo"]], [s[f"{k}_mae_tod_hi"] - m]], fmt="o", ms=8,
                        color=c, lw=1.6, capsize=0, zorder=3)
    ax.set_yticks(list(y.values()), [TITLE[k] for k in y])
    ax.set_ylim(-0.6, 1.8)
    ax.set_xlim(0, 4)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("median absolute error vs\nrecorded time of death [h]")
    s = summ.loc["full"]
    ax.annotate("CHIRAL DIP", (s.dip_mae_tod_h, 1.06), xytext=(0, 12), textcoords="offset points", ha="center", fontsize=10,
                color=CHIRAL)
    ax.annotate("clocksync", (s.it_mae_tod_h, 0.94), xytext=(0, -12), textcoords="offset points", ha="center", va="top",
                fontsize=10, color=CS)


def main():
    pf.use_style("slide")
    OUT.mkdir(parents=True, exist_ok=True)
    summ = pd.read_csv(SRC / "summary.csv").set_index("design")
    bins = pd.read_csv(SRC / "bins.csv")
    lay = pf.Layout.slide("ABC", kind="wide")
    for key, design in zip("AB", ("full", "window")):
        panel(lay.axes[key], design, bins)
    tod_panel(lay.axes["C"], summ)
    fig = lay.finalize()
    print(summ.round(3).T.to_string())
    print(pf.collisions(fig))
    pf.save(fig, OUT / "gtex_dip")


if __name__ == "__main__":
    main()
