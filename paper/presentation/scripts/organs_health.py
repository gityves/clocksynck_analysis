"""Slide figure (results, 3.4a; formerly R2.1): louder organs are more coupled to the body's clock state, in human and
in baboon, brain included. Why is slide 3.4b (`organs_departure.py`).

A. GTEx fast deaths, 27 peripheral tissues (A104 / A294b table, leave-one-out peripheral host): coupling against the
   amplitude shared with the body's clock state, sqrt(S) (S = z2 V_h, gain squared times host spread), one dot per
   tissue coloured by organ system; grey: c = S / (S + D) at the median departure D (the definition of c drawn at a
   constant D, not a fit). Spearman rho with amplitude +0.734 (A104; A291 registered: a perfectly locked tissue reads
   0.97 at every design).
B. Baboon (Mure 2018, 12 animals): rhythm amplitude against known time (A36 timed gain |g|, estimated without the host)
   against coupling (A42 readout 2), 35 peripheral tissues and 22 brain regions; grey: the same curve as A,
   c = a^2 / (a^2 + D), with D fitted on the peripheral tissues (least squares, one parameter; amplitude units differ
   from A, so D is not comparable). Post hoc panel (A289 H289.3 registered the untimed amplitude: brain on the line). Pineal, pituitary and ocular
   tissues not shown.

Right axes: the phase SD a coupling implies, (24 / 2 pi) sqrt((1 - c) / (2 c)) hours, for a departure of equal
amplitude and phase parts and small angles (as slide 3.4b, panel C); approximate beyond about 3 h.

Data: exploratory/A294_floor/output/tissues.csv; results/figure_2/output/organs/r2_couplings.csv;
results/figure_1/output/species/a36_tissues.csv. GTEx aggregates only (one value per tissue).
Output: presentation/figures/organs_health/organs_health.{svg,png} (half).

    python scripts/organs_health.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar

import paperfig as pf

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _systems import COLOR, baboon_system, gtex_system, legend_in  # noqa: E402

PAPER = Path(__file__).resolve().parents[2]
EXP = PAPER / "exploratory"
OUT = PAPER / "presentation" / "figures" / "organs_health"
INK, MUTED = "#222222", pf.GREY


H = 24 / (2 * np.pi)


def hours_axis(ax, label=True):
    """Right axis: the phase SD a coupling implies, SD = H sqrt((1 - c) / (2 c)) (isotropic departure, small
    angle; as slide 3.4b)."""
    fwd = lambda c: H * np.sqrt((1 - np.clip(c, 1e-3, 1)) / (2 * np.clip(c, 1e-3, 1)))  # noqa: E731
    inv = lambda h: 1 / (1 + 2 * (np.asarray(h) / H) ** 2)  # noqa: E731
    sec = ax.secondary_yaxis("right", functions=(fwd, inv))
    sec.set_yticks([1, 1.5, 2, 3, 4, 6], ["1", "1.5", "2", "3", "4", "6"])
    if label:
        sec.set_ylabel("phase SD [h]")
    return sec


def scatter(ax, t, system, title, x="own_amplitude", xlabel="clock-gene variation"):
    """Amplitude column `x` against coupling, with c = a^2 / (a^2 + D), D fitted on the peripheral tissues (least
    squares, one parameter); returns D and the median residual of the peripheral and brain tissues."""
    t = t[t.group.isin(["periphery", "brain"])].rename(columns={x: "amp"})
    per = t[t.group == "periphery"]
    curve = lambda a, d: a ** 2 / (a ** 2 + d)  # noqa: E731
    fit = minimize_scalar(lambda ld: ((per.coupling - curve(per.amp, np.exp(ld))) ** 2).sum(), bounds=(-6, 6),
                          method="bounded")
    d = float(np.exp(fit.x))
    xs = np.geomspace(t.amp.min() * 0.8, t.amp.max() * 1.2, 80)
    ax.plot(xs, curve(xs, d), color=MUTED, lw=1.2, zorder=1)
    ax.text(0.97, 0.03, f"a\u00b2 / (a\u00b2 + {d:.2f})", transform=ax.transAxes, ha="right", va="bottom", fontsize=8.5,
            color=MUTED)
    ax.scatter(t.amp, t.coupling, s=24, c=[COLOR[system(k)] for k in t.index], edgecolor="white", lw=0.4,
               zorder=2)
    ax.set_xscale("log")
    ax.minorticks_off()
    ticks = [v for v in (0.25, 0.5, 1, 2, 4) if t.amp.min() * 0.8 <= v <= t.amp.max() * 1.2]
    ax.set_xticks(ticks, [f"{v:g}" for v in ticks])
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlabel(xlabel)
    ax.set_title(title, fontsize=13, color=INK)
    res = t.coupling - curve(t.amp, d)
    return d, float(res[t.group == "periphery"].median()), float(res[t.group == "brain"].median())


def coupling_curve(ax, t):
    """Coupling against the shared amplitude sqrt(S), with S / (S + D) at the median D."""
    ax.scatter(np.sqrt(t.S), t.coupling, s=24, c=[COLOR[gtex_system(k)] for k in t.index], edgecolor="white", lw=0.4,
               zorder=2)
    xs = np.geomspace(0.3, 4, 80)
    d = t.own.median()
    ax.plot(xs, xs**2 / (xs**2 + d), color=MUTED, lw=1.2, zorder=1)
    ax.text(0.97, 0.03, f"a\u00b2 / (a\u00b2 + {d:.2f})", transform=ax.transAxes, ha="right", va="bottom", fontsize=8.5,
            color=MUTED)
    ax.set_xscale("log")
    ax.set_xlim(0.3, 4)
    ax.set_xticks([0.5, 1, 2, 4], ["0.5", "1", "2", "4"])
    ax.minorticks_off()
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlabel("amplitude shared\nwith the body, \u221aS")
    ax.set_ylabel("coupling")


def main():
    pf.use_style("slide")
    OUT.mkdir(parents=True, exist_ok=True)
    t = pd.read_csv(EXP / "A294_floor/output/tissues.csv", index_col="tissue")
    bab = pd.read_csv(PAPER / "results/figure_2/output/organs/r2_couplings.csv").query("species == 'baboon'")
    lay = pf.Layout.slide("AxB", kind="half", col_ratios=[1, 0.13, 1], blank=("x",))
    coupling_curve(lay.axes["A"], t)
    legend_in(lay.axes["A"], fontsize=7.5, loc="upper left")
    lay.axes["A"].set_title("human", fontsize=13, color=INK)
    hours_axis(lay.axes["A"])
    timed = pd.read_csv(PAPER / "results/figure_1/output/species/a36_tissues.csv", index_col="tissue").g_abs
    bab = bab.set_index("compartment").assign(timed_amp=timed)
    d, r_per, r_brain = scatter(lay.axes["B"], bab, baboon_system, "baboon", x="timed_amp",
                                xlabel="amplitude against\nknown time")
    print(f"baboon fitted D {d:.3f}; median residual periphery {r_per:.3f}, brain {r_brain:.3f}")
    lay.axes["B"].set_ylabel("coupling")
    hours_axis(lay.axes["B"])
    fig = lay.finalize()
    print(pf.collisions(fig))
    pf.save(fig, OUT / "organs_health")


if __name__ == "__main__":
    main()
