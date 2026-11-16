"""Slide figure (results, R1.1): the GTEx organ clock map, as in the CHIRAL paper, from deaths between 08 and 18 h.

One dot per tissue on a clock face: angle = clock phase offset of the tissue (h, centred on the circular mean over tissues,
later clockwise), radius = relative clock amplitude (linear, relative to the geometric mean over tissues).
A. Timed reference, truth (i): per-tissue cosinor of the 14 clock genes on the recorded time of death of all fast deaths
   (Hardy 1-2, whole day), rank-1 complex SVD.
B. clocksync on the 111 fast deaths who died 08-18 h, times hidden (A123, waveform cross-fitted outside the window).
C. CHIRAL (pychiral) on the same donors, per-donor phases combined across tissues, then cosinor + SVD as A.
In B and C the reference (grey ring) is A and a thin line joins each tissue to it.

Data: exploratory/A123_working_hours/output/organ_map.csv (per-tissue aggregates; A123 registered, H123.1-2 scored
2026-10-24) and tests.csv (the statistics in the titles). 40 tissues; Whole Blood has no truth (ii) but has truth (i).

Output: paper/presentation/figures/gtex_organ_map/gtex_organ_map.{svg,png} (wide slide figure).

    python scripts/gtex_organ_map.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.patches import Arc

import paperfig as pf
from _systems import COLOR, gtex_system, legend_row

PAPER = Path(__file__).resolve().parents[2]
SRC = PAPER / "exploratory/A123_working_hours/output"
OUT = Path(__file__).resolve().parents[1] / "figures" / "gtex_organ_map"
HOURS = 24 / (2 * np.pi)
INK, MUTED, GHOST = "#222222", pf.GREY, "#c8c8c8"
RINGS = (0.5, 1.0, 2.0)


def xy(offset_h, log2_amp):
    """Clock face: 0 h on top, clockwise; radius the linear relative amplitude."""
    th, r = np.asarray(offset_h) / HOURS, 2.0 ** np.asarray(log2_amp)
    return np.column_stack([r * np.sin(th), r * np.cos(th)])


def face(ax, rmax):
    """Upper half of the clock face (-6 h to +6 h around the mean): the tissue offsets span a few hours."""
    for r in RINGS:
        ax.add_patch(Arc((0, 0), 2 * r, 2 * r, theta1=0, theta2=180, ec="#e6e6e6", lw=0.8, zorder=0))
        a = -5.0 / HOURS  # ring labels on the -5 h spoke, away from the data
        ax.annotate(f"{r:g}x", (r * np.sin(a), r * np.cos(a)), ha="center", va="center", fontsize=8, color=MUTED, zorder=0)
    for h in range(-6, 7):
        a = h / HOURS
        ax.plot([0, rmax * np.sin(a)], [0, rmax * np.cos(a)], color="#f0f0f0" if h % 3 else "#e0e0e0", lw=0.6, zorder=0)
    for h in (-6, -3, 0, 3, 6):
        a = h / HOURS
        lab = {0: "0 h", -3: "-3 h", 3: "+3 h", -6: "-6 h\nearlier", 6: "+6 h\nlater"}[h]
        ax.annotate(lab, (1.1 * rmax * np.sin(a), 1.1 * rmax * np.cos(a) + (0.12 * rmax if abs(h) == 6 else 0)),
                    ha="center", va="center", fontsize=9, color=MUTED)
    ax.annotate("radius: relative amplitude", (0, -0.12 * rmax), ha="center", va="top", fontsize=9, color=MUTED)
    ax.set_xlim(-1.25 * rmax, 1.25 * rmax)
    ax.set_ylim(-0.3 * rmax, 1.2 * rmax)
    ax.set_aspect("equal")
    ax.axis("off")


def panel(ax, m, est, ref, title, rmax):
    face(ax, rmax)
    p = xy(m[f"{est}_offset_h"], m[f"{est}_log2_amp"])
    col = [COLOR[gtex_system(t)] for t in m.index]
    if ref is not None:
        q = xy(m[f"{ref}_offset_h"], m[f"{ref}_log2_amp"])
        for a, b in zip(q, p):
            ax.plot([a[0], b[0]], [a[1], b[1]], color=GHOST, lw=0.7, zorder=1)
        ax.scatter(*q.T, s=22, facecolor="none", edgecolor=GHOST, lw=0.8, zorder=2)
    ax.scatter(*p.T, s=34, c=col, edgecolor="white", lw=0.5, zorder=3)
    ax.set_title(title, fontsize=12, color=INK)


def main():
    pf.use_style("slide")
    OUT.mkdir(parents=True, exist_ok=True)
    m = pd.read_csv(SRC / "organ_map.csv", index_col="tissue").dropna(subset=["est_offset_h", "chiral_offset_h", "truth_i_offset_h"])
    tests = pd.read_csv(SRC / "tests.csv").set_index(["method", "truth"])
    rmax = 1.05 * float(2 ** m[["est_log2_amp", "chiral_log2_amp", "truth_i_log2_amp"]].max().max())
    lay = pf.Layout.slide("ABC;LLL", kind="wide", row_heights=[1, 0.1])
    panel(lay.axes["A"], m, "truth_i", None, "timed, full day", rmax)
    panel(lay.axes["B"], m, "est", "truth_i", "clocksync, hospital hours", rmax)
    panel(lay.axes["C"], m, "chiral", "truth_i", "CHIRAL, hospital hours", rmax)
    legend_row(lay.axes["L"])
    fig = lay.finalize()
    print(f"{len(m)} tissues; rmax {rmax:.2f}")
    print(tests.loc[[("clocksync", "truth_i"), ("chiral", "truth_i")], ["circ_corr_offset", "median_abs_offset_err_h",
                                                                       "pearson_log2_amp"]].round(3).to_string())
    print(pf.collisions(fig))
    pf.save(fig, OUT / "gtex_organ_map")


if __name__ == "__main__":
    main()
