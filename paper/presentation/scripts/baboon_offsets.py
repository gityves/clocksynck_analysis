"""Slide figure (results, R1.4): a primate sampled around the clock (baboon atlas, Mure 2018; A36, archive amendment re-run on
current code, descriptive).

A. 63 tissues: untimed offset (clocksync, every tissue an anchor, GTEx waveform) against the timed offset (phase of the
   tissue's 12-ZT cosinor gain on the waveform), both centred on the timed circular mean after one common shift; peripheral
   and neural / ocular tissues in two colours; dashed: identity.
B. Host phase of each animal (one animal per ZT) against its ZT, after one common shift; dashed: identity.

Data: results/figure_1/output/species/a36_tissues.csv; presentation/output/baboon/host_phase.csv (`baboon_compute.py`).
Output: presentation/figures/baboon_offsets/baboon_offsets.{svg,png} (half-slide figure, beside the bullets).

    python scripts/baboon_offsets.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns

import paperfig as pf

PAPER = Path(__file__).resolve().parents[2]
OUT = PAPER / "presentation" / "figures" / "baboon_offsets"
INK, MUTED = "#222222", pf.GREY
TAB = sns.color_palette().as_hex()
CS, PERIPHERAL, NEURAL = TAB[0], TAB[8], TAB[6]
W = 2 * np.pi / 24


def wrap(h):
    return (np.asarray(h, float) + 12) % 24 - 12


def cmean(h):
    return np.angle(np.exp(1j * W * np.asarray(h, float)).mean()) / W


def tissues(ax, t):
    timed = t.timed_phase / W
    untimed = t.offset_h - cmean(t.offset_h - timed)
    x = wrap(timed - cmean(timed))
    y = x + wrap(untimed - timed)
    ax.plot([-7, 7], [-7, 7], color=MUTED, lw=1, ls="--", zorder=1)
    for m, c, lab in ((~t.neural, PERIPHERAL, "peripheral"), (t.neural, NEURAL, "neural, ocular")):
        ax.scatter(x[m], y[m], s=26, color=c, edgecolor="white", lw=0.5, zorder=3, label=lab)
    ax.set_xlim(-7, 7)
    ax.set_ylim(-7, 7)
    ax.set_xticks(range(-6, 7, 3))
    ax.set_yticks(range(-6, 7, 3))
    ax.set_aspect("equal")
    ax.set_xlabel("timed offset [h]")
    ax.set_ylabel("untimed offset [h]")
    ax.legend(loc="upper left", frameon=False, fontsize=10, handletextpad=0.2)
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
    lay = pf.Layout.slide("AB", kind="half", height=100)
    tissues(lay.axes["A"], t)
    animals(lay.axes["B"], h)
    fig = lay.finalize()
    print(len(t), "tissues;", len(h), "animals")
    print(pf.collisions(fig))
    pf.save(fig, OUT / "baboon_offsets")


if __name__ == "__main__":
    main()
