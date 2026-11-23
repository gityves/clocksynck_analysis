"""Slide figure (results, adult health, 3.1): the reference organisation of the body's clock, GTEx fast deaths over the
whole day, read untimed by clocksync.

One dot per tissue on the clock face of slide 1.1: angle = peak time of the tissue's clock genes relative to the person's
body time (0 h, host from the peripheral tissues; peripheral tissues read leave-one-out; later peaks clockwise, the
clocksync offset negated), radius = relative amplitude
(gain on the body's time, linear). Bars: 95 % donor-bootstrap intervals (200 draws) of the offset (arc) and of the
amplitude (radial). Colours: organ systems; 13 brain regions read on the peripheral host. Direct labels on the
tissues discussed in the deck.

Status: A263 R2 fit (registered readout), per-tissue intervals added by `organs_compute.py`; display descriptive.
GTEx aggregates only (one row per tissue).
Data: presentation/output/organs/couplings.csv. Output: presentation/figures/reference_map/reference_map.{svg,png} (half).

    python scripts/reference_map.py
"""
from pathlib import Path

import numpy as np
import pandas as pd

import paperfig as pf
from _systems import COLOR, gtex_system, legend_in
from gtex_organ_map import HOURS, face

PAPER = Path(__file__).resolve().parents[2]
OUT = PAPER / "presentation" / "figures" / "reference_map"
LABEL = {"Adrenal Gland": ("adrenal", (-6, 4), "right"), "Lung": ("lung", (-6, 0), "right"),
         "Skin - Not Sun Exposed (Suprapubic)": ("skin", (-6, -2), "right"),
         "Whole Blood": ("blood", (-6, -6), "right"), "Testis": ("testis", (22, -26), "left"),
         "Thyroid": ("thyroid", (-8, 2), "right"), "Brain - Hippocampus": ("brain regions", (10, 8), "left")}


def xy(off, amp):
    th = -np.asarray(off) / HOURS  # peak time: clocksync offset (> 0 advanced) negated, as CHIRAL Fig. 1B
    return np.asarray(amp) * np.sin(th), np.asarray(amp) * np.cos(th)


def draw(ax, d, rmax):
    face(ax, rmax, label_h=5.0, labelled=(1.0, 2.0))
    for r in d.itertuples():
        col = COLOR[gtex_system(r.tissue)]
        lo, hi = np.exp(r.log_amp_lo), np.exp(r.log_amp_hi)
        ax.plot(*xy([r.offset_h] * 2, [lo, hi]), color=col, lw=1.0, alpha=0.6, zorder=2)
        u = np.linspace(r.offset_lo, r.offset_hi, 30)
        ax.plot(*xy(u, np.full_like(u, r.rel_amplitude)), color=col, lw=1.0, alpha=0.6, zorder=2)
    x, y = xy(d.offset_h, d.rel_amplitude)
    ax.scatter(x, y, s=34, c=[COLOR[gtex_system(t)] for t in d.tissue], edgecolor="white", lw=0.5, zorder=3)
    for t, (lab, off, ha) in LABEL.items():
        r = d.set_index("tissue").loc[t]
        ax.annotate(lab, xy(r.offset_h, r.rel_amplitude), xytext=off, textcoords="offset points", ha=ha, va="center",
                    fontsize=9, color=COLOR[gtex_system(t)],
                    arrowprops=dict(arrowstyle="-", color=COLOR[gtex_system(t)], lw=0.6, shrinkA=1, shrinkB=4)
                    if abs(off[0]) > 12 else None)


def main():
    pf.use_style("slide")
    OUT.mkdir(parents=True, exist_ok=True)
    d = pd.read_csv(PAPER / "presentation/output/organs/couplings.csv")
    rmax = 1.05 * float(np.exp(d.log_amp_hi).max())
    lay = pf.Layout.slide("A", kind="half")
    draw(lay.axes["A"], d, rmax)
    legend_in(lay.axes["A"], fontsize=10, loc="lower right")
    fig = lay.finalize()
    pd.set_option("display.width", 200)
    print(d[["tissue", "n", "offset_h", "offset_lo", "offset_hi", "rel_amplitude"]].round(2).to_string(index=False))
    print(f"{len(d)} tissues; donors per tissue {d.n.min()}-{d.n.max()}; rmax {rmax:.2f}")
    print(pf.collisions(fig))
    pf.save(fig, OUT / "reference_map")


if __name__ == "__main__":
    main()
