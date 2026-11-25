"""Slide figure (results, drives, 4.2c): the reference map of slide 3.1 for fast deaths, slow deaths, ventilated donors
with a living brain and brain death.

One clock face per group, as slide 3.1: angle = peak time of the tissue's clock genes relative to the group's body time
(0 h; host from the peripheral tissues, read leave-one-out; later peaks clockwise, clocksync offset negated), radius =
relative amplitude (gain on the group's own body time). Bars: 95 % donor-bootstrap intervals (200 draws). Tissues with
>= 30 donors in the group; organ-system colours; same radius scale on every face. Each face is relative to its own
group's body clock: it shows arrangement and relative loudness, not the body-wide loss of loudness (slide 4.2b).

Status: descriptive, post hoc display (user, 2026-11-25); readout as slide 3.1 (A263 R2 fit). GTEx aggregates only
(one row per group and tissue).
Data: presentation/output/organs/couplings_{fast,slow,ventilated,brain_death}.csv (`organs_groups_compute.py`).
Output: presentation/figures/organs_groups/organs_groups.{svg,png} (wide).

    python scripts/organs_groups.py
"""
from pathlib import Path

import numpy as np
import pandas as pd

from matplotlib import patheffects
from matplotlib.lines import Line2D

import paperfig as pf

from _systems import COLOR, gtex_system
from gtex_organ_map import face
from reference_map import xy

PAPER = Path(__file__).resolve().parents[2]
SRC = PAPER / "presentation" / "output" / "organs"
OUT = PAPER / "presentation" / "figures" / "organs_groups"
GROUPS = {"fast": "fast deaths", "slow": "slow deaths", "ventilated": "ventilated,\nliving brain",
          "brain_death": "brain death"}
LABEL = {"Adrenal Gland": ("adrenal", (-6, 0), "right"), "Lung": ("lung", (7, -5), "left"),
         "Thyroid": ("thyroid", (-6, 5), "right")}  # label, offset [pt], alignment


def draw(ax, d, rmax, labels=tuple(LABEL)):
    face(ax, rmax, label_h=5.0, labelled=(1.0, 2.0))
    for r in d.itertuples():
        col = COLOR[gtex_system(r.tissue)]
        lo, hi = np.exp(r.log_amp_lo), np.exp(r.log_amp_hi)
        ax.plot(*xy([r.offset_h] * 2, [lo, hi]), color=col, lw=1.0, alpha=0.6, zorder=2)
        u = np.linspace(r.offset_lo, r.offset_hi, 30)
        ax.plot(*xy(u, np.full_like(u, r.rel_amplitude)), color=col, lw=1.0, alpha=0.6, zorder=2)
    x, y = xy(d.offset_h, d.rel_amplitude)
    ax.scatter(x, y, s=24, c=[COLOR[gtex_system(t)] for t in d.tissue], edgecolor="white", lw=0.5, zorder=3)
    for t, (lab, off, ha) in LABEL.items():
        if t in labels and t in set(d.tissue):
            r = d.set_index("tissue").loc[t]
            ax.annotate(lab, xy(r.offset_h, r.rel_amplitude), xytext=off, textcoords="offset points", ha=ha,
                        va="center", fontsize=9, color=COLOR[gtex_system(t)],
                        path_effects=[patheffects.withStroke(linewidth=2.5, foreground="white")])


def main():
    pf.use_style("slide")
    OUT.mkdir(parents=True, exist_ok=True)
    d = {g: pd.read_csv(SRC / f"couplings_{g}.csv") for g in GROUPS}
    rmax = 1.05 * max(float(np.exp(x.log_amp_hi).max()) for x in d.values())
    lay = pf.Layout.slide("ABCD", kind="wide", height=70)  # half faces: a shorter frame than the preset
    for (g, title), k in zip(GROUPS.items(), "ABCD"):
        draw(lay.axes[k], d[g], rmax, tuple(LABEL) if g == "fast" else ("Adrenal Gland", "Thyroid"))  # lung: crowded
        lay.axes[k].set_title(title, fontsize=11)
    h = [Line2D([], [], ls="", marker="o", ms=7, color=c, label=s) for s, c in COLOR.items()]
    lay.shared_legend("ABCD", h, [x.get_label() for x in h], side="bottom", ncol=6, fontsize=10, handletextpad=0.2,
                      columnspacing=1.2)
    fig = lay.finalize()
    pd.set_option("display.width", 200)
    for g, x in d.items():
        print(g, f"{x.donors_group.iloc[0]} donors, {len(x)} tissues, donors per tissue {x.n.min()}-{x.n.max()}")
        print(x[["tissue", "n", "offset_h", "offset_lo", "offset_hi", "rel_amplitude", "coupling"]].round(2)
              .to_string(index=False))
    print(f"rmax {rmax:.2f}")
    print(pf.collisions(fig))
    pf.save(fig, OUT / "organs_groups")


if __name__ == "__main__":
    main()
