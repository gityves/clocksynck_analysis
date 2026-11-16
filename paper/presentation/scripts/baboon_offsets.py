"""Slide figure (results, R1.4): a primate sampled around the clock (baboon atlas, Mure 2018; A36, archive amendment re-run on
current code, descriptive).

A. Organ map as the GTEx figure: per tissue, angle = untimed offset (later clockwise), radius = untimed relative amplitude,
   coloured by organ system; grey rings joined by a line: the timed offset and amplitude (cosinor gain) of the same tissue.
B. 63 tissues: untimed offset (clocksync, every tissue an anchor, GTEx waveform) against the timed offset (phase of the
   tissue's 12-ZT cosinor gain on the waveform), both centred on the timed circular mean after one common shift; coloured
   by organ system as the GTEx organ map (`_systems.py`); dashed: identity.
C. Host phase of each animal (one animal per ZT) against its ZT, after one common shift; dashed: identity.

Data: results/figure_1/output/species/a36_tissues.csv; presentation/output/baboon/host_phase.csv (`baboon_compute.py`).
Output: presentation/figures/baboon_offsets/baboon_offsets.{svg,png} (wide slide figure).

    python scripts/baboon_offsets.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns

import paperfig as pf
from _systems import COLOR, baboon_system, legend_in
from gtex_organ_map import GHOST, face, xy

PAPER = Path(__file__).resolve().parents[2]
OUT = PAPER / "presentation" / "figures" / "baboon_offsets"
INK, MUTED = "#222222", pf.GREY
TAB = sns.color_palette().as_hex()
CS = TAB[0]
W = 2 * np.pi / 24


def wrap(h):
    return (np.asarray(h, float) + 12) % 24 - 12


def cmean(h):
    return np.angle(np.exp(1j * W * np.asarray(h, float)).mean()) / W


def clock_map(ax, t):
    """Organ map as the GTEx one: angle = untimed offset, radius = relative amplitude (both clocksync, untimed); grey rings and
    lines: the timed values of the same tissue."""
    timed = t.timed_phase / W
    untimed = t.offset_h - cmean(t.offset_h - timed)
    c = cmean(timed)
    la = np.log2(t.amplitude) - np.log2(t.amplitude).mean()
    lg = np.log2(t.g_abs) - np.log2(t.g_abs).mean()
    p, q = xy(wrap(untimed - c), la), xy(wrap(timed - c), lg)
    rmax = 1.05 * float(2 ** max(la.max(), lg.max()))
    face(ax, rmax, label_h=-5.5, labelled=(1.0, 2.0))  # the 0.5x ring sits in the brain cluster
    for a, b in zip(q, p):
        ax.plot([a[0], b[0]], [a[1], b[1]], color=GHOST, lw=0.7, zorder=1)
    ax.scatter(*q.T, s=18, facecolor="none", edgecolor=GHOST, lw=0.8, zorder=2)
    ax.scatter(*p.T, s=26, c=[COLOR[baboon_system(k)] for k in t.tissue], edgecolor="white", lw=0.5, zorder=3)
    ax.set_title("organ map", fontsize=13, color=INK)


def tissues(ax, t):
    timed = t.timed_phase / W
    untimed = t.offset_h - cmean(t.offset_h - timed)
    x = wrap(timed - cmean(timed))
    y = x + wrap(untimed - timed)
    ax.plot([-7, 7], [-7, 7], color=MUTED, lw=1, ls="--", zorder=1)
    ax.scatter(x, y, s=26, c=[COLOR[baboon_system(k)] for k in t.tissue], edgecolor="white", lw=0.5, zorder=3)
    ax.set_xlim(-7, 7)
    ax.set_ylim(-7, 7)
    ax.set_xticks(range(-6, 7, 3))
    ax.set_yticks(range(-6, 7, 3))
    ax.set_aspect("equal")
    ax.set_xlabel("timed offset [h]")
    ax.set_ylabel("untimed offset [h]")
    ax.set_title("tissues", fontsize=13, color=INK)


def animals(ax, h):
    y = h.zt_h + wrap(h.host_h - cmean(h.host_h - h.zt_h) - h.zt_h)
    ax.plot([-2, 24], [-2, 24], color=MUTED, lw=1, ls="--", zorder=1)
    ax.scatter(h.zt_h, y, s=40, color=CS, edgecolor="white", lw=0.5, zorder=3)
    ax.set_xlim(-2, 24)
    ax.set_ylim(-2, 24)
    ax.set_xticks(range(0, 25, 6))
    ax.set_yticks(range(0, 25, 6))
    ax.set_aspect("equal")
    ax.set_xlabel("zeitgeber time [h]")
    ax.set_ylabel("internal time [h]")
    ax.set_title("animals", fontsize=13, color=INK)


def main():
    pf.use_style("slide")
    OUT.mkdir(parents=True, exist_ok=True)
    t = pd.read_csv(PAPER / "results/figure_1/output/species/a36_tissues.csv")
    h = pd.read_csv(PAPER / "presentation/output/baboon/host_phase.csv")
    lay = pf.Layout.slide("ABC", kind="wide", col_ratios=[1.5, 1, 1])
    clock_map(lay.axes["A"], t)
    tissues(lay.axes["B"], t)
    animals(lay.axes["C"], h)
    legend_in(lay.axes["B"])
    fig = lay.finalize()
    print(len(t), "tissues;", len(h), "animals")
    print(pf.collisions(fig))
    pf.save(fig, OUT / "baboon_offsets")


if __name__ == "__main__":
    main()
