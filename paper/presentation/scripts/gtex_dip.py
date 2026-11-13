"""Slide figure (results, R1.2): CHIRAL's donor internal phase (DIP) against clocksync's untimed internal time, GTEx fast deaths.

Reads the aggregates written by `gtex_dip_compute.py` (no per-donor value: controlled access). Two panels, clocksync fitted on
A. all fast deaths (whole day) and B. the fast deaths who died 08-18 h only. x: CHIRAL's DIP (learnt on all GTEx donors, whole
day); y: clocksync's internal time, aligned to the DIP by the best rotation and reflection. Grey cells: number of donors per
3 h x 3 h cell (cells with < 5 donors left blank); dots and bars: circular mean and circular SD of the internal time per 3 h bin
of DIP (bins with < 5 donors omitted); dashed: identity. C. |circular r| of each internal phase with the recorded time of death,
on the same donors (summary statistics).

Output: paper/presentation/figures/gtex_dip/gtex_dip.{svg,png} (wide slide figure).

    python scripts/gtex_dip.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

import paperfig as pf

SRC = Path(__file__).resolve().parents[1] / "output" / "gtex_dip"
OUT = Path(__file__).resolve().parents[1] / "figures" / "gtex_dip"
INK, MUTED = "#222222", pf.GREY
COLOR = pf.OKABE_ITO[3]
CMAP = LinearSegmentedColormap.from_list("cells", ["#f4f4f4", "#7a7a7a"])
TITLE = {"full": "clocksync on the whole day", "window": "clocksync on 08-18 h only"}


def panel(ax, design, s, bins, hist, vmax):
    h = hist[hist.design == design]
    for r in h.dropna(subset=["n"]).itertuples():
        ax.add_patch(pf_rect(r.dip_lo, r.it_lo, CMAP(r.n / vmax)))
    b = bins[(bins.design == design)].dropna(subset=["it_mean_h"])
    x = (b.dip_lo + b.dip_hi) / 2
    y = b.it_mean_h + 24 * np.round((x - b.it_mean_h) / 24)  # unwrap towards the diagonal
    ax.plot([0, 24], [0, 24], color=MUTED, lw=1, ls="--", zorder=1)
    ax.errorbar(x, y, yerr=b.it_circ_sd_h, fmt="o", color=COLOR, ms=6, lw=1.6, capsize=0, zorder=3)
    ax.set_xlim(0, 24)
    ax.set_ylim(-6, 30)
    ax.set_xticks(range(0, 25, 6))
    ax.set_yticks(range(0, 25, 6))
    ax.set_aspect("equal")
    ax.set_xlabel("CHIRAL donor internal phase (h)")
    ax.set_ylabel("clocksync internal time (h, aligned)")
    ax.set_title(f"{TITLE[design]} ({s.n_donors} donors)\ncircular r = {s.circ_corr_dip:.2f}, "
                 f"median |difference| {s.median_abs_diff_h:.1f} h", fontsize=12, color=INK)


def tod_panel(ax, summ):
    """Agreement of each internal phase with the recorded time of death, same donors (|circular r|)."""
    y = {"full": 1, "window": 0}
    for design, yy in y.items():
        s = summ.loc[design]
        ax.plot([s.circ_corr_dip_tod, s.circ_corr_it_tod], [yy, yy], color="#dddddd", lw=2, zorder=1)
        ax.scatter([s.circ_corr_dip_tod], [yy], s=70, color=MUTED, zorder=2)
        ax.scatter([s.circ_corr_it_tod], [yy], s=70, color=COLOR, zorder=3)
    ax.set_yticks(list(y.values()), ["all fast deaths\n(whole day)", "died 08-18 h"])
    ax.set_ylim(-0.6, 1.8)
    ax.set_xlim(0, 1)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("|circular r| with recorded time of death")
    ax.annotate("CHIRAL DIP", (summ.loc["full", "circ_corr_dip_tod"], 1), xytext=(0, 10), textcoords="offset points",
                ha="center", fontsize=10, color=MUTED)
    ax.annotate("clocksync", (summ.loc["full", "circ_corr_it_tod"], 1), xytext=(0, -12), textcoords="offset points",
                ha="center", va="top", fontsize=10, color=COLOR)
    ax.set_title("which phase follows the clock?", fontsize=12, color=INK)


def pf_rect(x, y, c, w=3.0):
    from matplotlib.patches import Rectangle
    return Rectangle((x, y), w, w, fc=c, ec="white", lw=0.8, zorder=0)


def main():
    pf.use_style("slide")
    OUT.mkdir(parents=True, exist_ok=True)
    summ = pd.read_csv(SRC / "summary.csv").set_index("design")
    bins, hist = pd.read_csv(SRC / "bins.csv"), pd.read_csv(SRC / "hist2d.csv")
    vmax = float(hist.n.max())
    lay = pf.Layout.slide("ABC", kind="wide")
    for key, design in zip("AB", ("full", "window")):
        panel(lay.axes[key], design, summ.loc[design], bins, hist, vmax)
    tod_panel(lay.axes["C"], summ)
    fig = lay.finalize()
    print(summ.round(3).T.to_string())
    print(pf.collisions(fig))
    pf.save(fig, OUT / "gtex_dip")


if __name__ == "__main__":
    main()
