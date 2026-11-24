"""Slide figure (results, 3.4a; formerly R2.1): louder organs are more coupled to the body's clock state, in human and
in baboon, brain included. Why is slide 3.4b (`organs_departure.py`).

A. GTEx fast deaths, 27 peripheral tissues (A104 / A294b table, leave-one-out peripheral host): coupling against the
   amplitude shared with the body's clock state, sqrt(S) (S = z2 V_h, gain squared times host spread), one dot per
   tissue coloured by organ system; grey: c = S / (S + D) at the median departure D (the definition of c drawn at a
   constant D, not a fit). Spearman rho with amplitude +0.734 (A104; A291 registered: a perfectly locked tissue reads
   0.97 at every design).
B. Baboon (Mure 2018, 12 animals): rhythm amplitude against known time (A36 timed gain |g|, estimated without the host)
   against coupling (A42 readout 2), 35 peripheral tissues and 22 brain regions; grey: the peripheral least-squares
   line. Post hoc panel (A289 H289.3 registered the untimed amplitude: brain on the line). Pineal, pituitary and ocular
   tissues not shown.

Data: exploratory/A294_floor/output/tissues.csv; results/figure_2/output/organs/r2_couplings.csv;
results/figure_1/output/species/a36_tissues.csv. GTEx aggregates only (one value per tissue).
Output: presentation/figures/organs_health/organs_health.{svg,png} (half).

    python scripts/organs_health.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

import paperfig as pf

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _systems import COLOR, baboon_system, gtex_system, legend_in  # noqa: E402

PAPER = Path(__file__).resolve().parents[2]
EXP = PAPER / "exploratory"
OUT = PAPER / "presentation" / "figures" / "organs_health"
INK, MUTED = "#222222", pf.GREY


def scatter(ax, t, system, title, x="own_amplitude", xlabel="clock-gene variation"):
    """Amplitude column `x` against coupling, peripheral line; `t`: tissue index, columns group, coupling, `x`."""
    t = t[t.group.isin(["periphery", "brain"])].rename(columns={x: "amp"})
    per = t[t.group == "periphery"]
    b, a = np.polyfit(np.log(per.amp), per.coupling, 1)
    xs = np.geomspace(t.amp.min() * 0.9, t.amp.max() * 1.1, 50)
    ax.plot(xs, a + b * np.log(xs), color=MUTED, lw=1.2, zorder=1)
    ax.scatter(t.amp, t.coupling, s=24, c=[COLOR[system(k)] for k in t.index], edgecolor="white", lw=0.4,
               zorder=2)
    ax.set_xscale("log")
    ax.minorticks_off()
    ticks = [v for v in (0.25, 0.5, 1, 2, 4) if t.amp.min() * 0.9 <= v <= t.amp.max() * 1.1]
    ax.set_xticks(ticks, [f"{v:g}" for v in ticks])
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlabel(xlabel)
    ax.set_title(title, fontsize=13, color=INK)


def coupling_curve(ax, t):
    """Coupling against the shared amplitude sqrt(S), with S / (S + D) at the median D."""
    ax.scatter(np.sqrt(t.S), t.coupling, s=24, c=[COLOR[gtex_system(k)] for k in t.index], edgecolor="white", lw=0.4,
               zorder=2)
    xs = np.geomspace(0.3, 4, 80)
    d = t.own.median()
    ax.plot(xs, xs**2 / (xs**2 + d), color=MUTED, lw=1.2, zorder=1)
    ax.text(2.2, 0.12, f"S / (S + {d:.1f})", ha="center", fontsize=8.5, color=MUTED)
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
    lay = pf.Layout.slide("AB", kind="half", col_ratios=[1, 1])
    coupling_curve(lay.axes["A"], t)
    legend_in(lay.axes["A"], fontsize=7.5, loc="upper left")
    lay.axes["A"].set_title("human", fontsize=13, color=INK)
    timed = pd.read_csv(PAPER / "results/figure_1/output/species/a36_tissues.csv", index_col="tissue").g_abs
    bab = bab.set_index("compartment").assign(timed_amp=timed)
    scatter(lay.axes["B"], bab, baboon_system, "baboon", x="timed_amp", xlabel="amplitude against\nknown time")
    lay.axes["B"].set_ylabel("coupling")
    fig = lay.finalize()
    print(pf.collisions(fig))
    pf.save(fig, OUT / "organs_health")


if __name__ == "__main__":
    main()
