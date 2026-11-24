"""Slide figure (results, 3.4b): why louder organs are more coupled. Each organ's departure from the body's clock state
has a nearly fixed size, not a fixed angle, so a louder organ's time carries a smaller phase error.

GTEx fast deaths, 27 peripheral tissues (A104 / A294b table). S: clock-plane variance shared with the body's clock state
(sqrt S, the shared amplitude); D = 2 psi: the tissue's departure at a fixed body state, net of counting noise.
A. Schematic drawn to scale with real values: the shared state of a quiet tissue (sun-exposed skin) and of a loud one
   (thyroid) on circles of radius sqrt S, each with its departure (disc of 1 SD per axis, sqrt(D / 2)); the wedge is the
   phase scatter it implies (+- 1 SD). Post hoc display.
B. D against S per tissue (log axes); grey lines through the geometric-mean S and median D: fixed size (slope 0), fixed
   angle (slope 1, a phase jitter of constant size), individual clock periods (slope 2). Observed slope 0.07 (Theil-Sen,
   A294b descriptive).
C. Phase scatter implied by the departure, sqrt(D / (2 S)) in hours (dots; small-angle approximation, post hoc display),
   against the shared amplitude; black diamonds: departure SD in hours measured by A259 (registered; 90 % interval), a
   different estimator on the same donors. Grey: prediction of a fixed-size departure at the median D (falls as
   1 / amplitude) and of a fixed angle (flat). Whole blood and testis omitted (no shown clock; small-angle fails).
D. Direction of the departure (A311, registered): fraction of D along the phase direction, per tissue, for the clock
   genes (system colours) and the median of 100 expression-matched sets (grey); black: medians. 0.5: as much amplitude
   as phase.

Data: presentation/output/organs_departure/{tissues,summary}.csv (`organs_departure_compute.py`).
Output: presentation/figures/organs_departure/organs_departure.{svg,png} (wide). GTEx aggregates only.

    python scripts/organs_departure.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Arc, Circle, Wedge

import paperfig as pf
from _systems import COLOR, gtex_system, legend_in

PAPER = Path(__file__).resolve().parents[2]
DATA = PAPER / "presentation" / "output" / "organs_departure"
OUT = PAPER / "presentation" / "figures" / "organs_departure"
INK, MUTED = "#222222", pf.GREY
H = 24 / (2 * np.pi)
QUIET, LOUD = "Skin - Sun Exposed (Lower leg)", "Thyroid"
NAMES = {QUIET: "skin (sun-exposed)", LOUD: "thyroid"}
NOCLOCK = ["Whole Blood", "Testis"]


def schematic(ax, t):
    for k, ang, pos in ((QUIET, 135, (-2.3, -0.6)), (LOUD, 55, (2.45, 3.5))):
        r = t.loc[k]
        c = COLOR[gtex_system(k)]
        th = np.deg2rad(ang)
        x, y = r.amp * np.cos(th), r.amp * np.sin(th)
        half = np.rad2deg(r.phase_h / H)
        ax.add_patch(Arc((0, 0), 2 * r.amp, 2 * r.amp, theta1=0, theta2=180, color=MUTED, lw=0.8, ls="--"))
        ax.add_patch(Wedge((0, 0), r.amp + 0.9, ang - half, ang + half, color=c, alpha=0.12, lw=0))
        ax.add_patch(Circle((x, y), np.sqrt(r.D / 2), color=c, alpha=0.35, lw=0))
        ax.plot([0, x], [0, y], color=c, lw=1.2)
        ax.scatter([x], [y], s=20, color=c, zorder=3)
        ax.text(*pos, f"{NAMES[k]}\n\u00b1{r.phase_h:.1f} h", ha="center", va="center",
                fontsize=9, color=c)
    ax.scatter([0], [0], s=14, color=INK, zorder=3)
    ax.text(0, -0.15, "body's\nclock state", ha="center", va="top", fontsize=8, color=INK)
    ax.text(2.4, -0.6, "disc: departure (1 SD)\nwedge: phase scatter", ha="center", va="center", fontsize=8,
            color=MUTED)
    ax.set_xlim(-3.6, 3.6)
    ax.set_ylim(-1.8, 4.9)
    ax.set_aspect("equal")
    ax.axis("off")


def log_ticks(ax, ticks, axis):
    getattr(ax, f"set_{axis}scale")("log")
    getattr(ax, f"set_{axis}ticks")(ticks, [f"{v:g}" for v in ticks])
    ax.minorticks_off()


def models(ax, t, s):
    ax.scatter(t.S, t.D, s=24, c=[COLOR[gtex_system(k)] for k in t.index], edgecolor="white", lw=0.4, zorder=3)
    s0, d0 = s.geomean_S, s.median_D
    xs = np.geomspace(0.1, 13, 50)
    for k, ls, lab in ((0, "-", "fixed size"), (1, "--", "fixed angle"), (2, ":", "own clock\nperiods")):
        ys = d0 * (xs / s0) ** k
        ax.plot(xs, ys, color=MUTED, lw=1.1, ls=ls, zorder=1)
        i = np.searchsorted(ys, 40) - 1 if k else len(xs) - 1
        xl = xs[min(i, len(xs) - 1)] if k else 13
        ax.text(xl, min(ys[min(i, len(xs) - 1)], 35) * (1.08 if k else 0.78), lab, fontsize=8, color=MUTED,
                ha="right" if k else "right", va="bottom" if k else "top")
    ax.text(0.97, 0.03, f"observed slope {s.slope:.2f}", transform=ax.transAxes, ha="right", va="bottom", fontsize=8.5,
            color=INK)
    ax.set_xlim(0.1, 13)
    ax.set_ylim(0.1, 40)
    log_ticks(ax, [0.1, 0.3, 1, 3, 10], "x")
    log_ticks(ax, [0.1, 0.3, 1, 3, 10, 30], "y")
    ax.set_xlabel("shared with the body, S")
    ax.set_ylabel("own departure, D")


def hours(ax, t, s):
    b = t.drop(index=NOCLOCK)
    xs = np.geomspace(0.9, 3.5, 50)
    ax.plot(xs, np.sqrt(s.median_D / 2) / xs * H, color=MUTED, lw=1.1, zorder=1)
    ax.axhline(b.phase_h.median(), color=MUTED, lw=1.1, ls="--", zorder=1)
    ax.text(3.45, b.phase_h.median() + 0.08, "fixed angle", ha="right", va="bottom", fontsize=8, color=MUTED)
    ax.text(1.0, np.sqrt(s.median_D / 2) / 1.0 * H, "fixed size", ha="left", va="bottom", fontsize=8, color=MUTED)
    ax.scatter(b.amp, b.phase_h, s=24, c=[COLOR[gtex_system(k)] for k in b.index], edgecolor="white", lw=0.4, zorder=3)
    a = b.dropna(subset=["a259_hours"])
    ax.errorbar(a.amp * 1.03, a.a259_hours, yerr=[a.a259_hours - a.a259_hours_lo, a.a259_hours_hi - a.a259_hours],
                fmt="D", ms=4.5, mfc="white", mec=INK, ecolor=INK, elinewidth=0.8, capsize=0, zorder=2)
    log_ticks(ax, [1, 2, 3], "x")
    ax.set_xlim(0.9, 3.5)
    ax.set_ylim(0, 4.2)
    ax.set_xlabel("amplitude shared\nwith the body, √S")
    ax.set_ylabel("phase scatter [h]")
    h = [Line2D([], [], ls="", marker="o", ms=5, color=MUTED, label="implied by D"),
         Line2D([], [], ls="", marker="D", ms=5, mfc="white", mec=INK, label="measured (A259)")]
    ax.legend(handles=h, loc="upper right", frameon=False, fontsize=8, handletextpad=0.1)


def direction(ax, t):
    rows = [("ftan_clock", "clock genes"), ("ftan_matched", "matched\ngenes")]
    for i, (col, _) in enumerate(rows):
        y = i + np.random.default_rng(i).uniform(-0.15, 0.15, len(t))
        c = [COLOR[gtex_system(k)] for k in t.index] if i == 0 else MUTED
        ax.scatter(t[col], y, s=18, c=c, alpha=0.7, edgecolor="none", zorder=2)
        ax.plot([t[col].median()] * 2, [i - 0.3, i + 0.3], color=INK, lw=2, zorder=3)
    ax.axvline(0.5, color=MUTED, lw=0.8, ls="--", zorder=0)
    ax.set_xlim(0, 1)
    ax.set_xticks([0, 0.5, 1], ["amplitude\nonly", "equal", "phase\nonly"])
    ax.set_yticks(range(len(rows)), [l for _, l in rows])
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(len(rows) - 0.4, -0.6)
    ax.set_xlabel("direction of the departure")


def main():
    pf.use_style("slide")
    OUT.mkdir(parents=True, exist_ok=True)
    t = pd.read_csv(DATA / "tissues.csv", index_col="tissue")
    s = pd.read_csv(DATA / "summary.csv").iloc[0]
    lay = pf.Layout.slide("ABCD", kind="wide", col_ratios=[1.05, 1, 1, 0.85])
    schematic(lay.axes["A"], t)
    models(lay.axes["B"], t, s)
    legend_in(lay.axes["B"], fontsize=7.5, loc="upper left")
    hours(lay.axes["C"], t, s)
    direction(lay.axes["D"], t)
    fig = lay.finalize()
    print(s.round(3).to_string())
    print(pf.collisions(fig))
    pf.save(fig, OUT / "organs_departure")


if __name__ == "__main__":
    main()
