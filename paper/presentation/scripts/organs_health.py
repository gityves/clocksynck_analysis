"""Slide figure (results, R2.1): each organ expresses the body's time with its own gain over a clock-specific departure
of nearly constant size; coupling follows, the brain included.

Notation (A104 table, GTEx fast deaths, peripheral tissues read leave-one-out on the peripheral host): S = z2 V, the
clock-plane variance a tissue shares with the body's clock state; D = 2 psi, its own variance; coupling c = S / (S + D).
A. Shared against own clock variation, one dot per tissue (log axes, same range), coloured by organ system; grey: Theil-Sen
   line of log D on log S (A294b descriptive slope). Post hoc display (A291 post hoc table).
B. Coupling against S; grey: c = S / (S + D) at the median D over tissues. The curve is the definition of c drawn at a
   constant D, not a fit: the finding is that D is nearly constant (A).
C. Is the own variation specific to the clock genes? Own variance of the clock genes against the median of 100
   expression-matched gene sets in the same tissue (A294b, registered, H294b.1 pass); grey: identity.
D. Baboon (Mure 2018, 12 animals): rhythm amplitude against known time (A36 timed gain |g|, estimated without the host)
   against coupling (A42 readout 2), 35 peripheral tissues and 22 brain regions; grey: the peripheral least-squares line.
   Post hoc panel (A289 H289.3 registered the untimed amplitude: brain on the line). Pineal, pituitary and ocular tissues
   not shown.

Data: exploratory/A294_floor/output/{tissues,tests}.csv; results/figure_2/output/organs/r2_couplings.csv;
results/figure_1/output/species/a36_tissues.csv. GTEx aggregates only (one value per tissue).
Output: presentation/figures/organs_health/organs_health.{svg,png}.

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


def log_ticks(ax, lo, hi, axis="x"):
    ticks = [v for v in (0.1, 0.3, 1, 3, 10) if lo <= v <= hi]
    getattr(ax, f"set_{axis}scale")("log")
    getattr(ax, f"set_{axis}lim")(lo, hi)
    getattr(ax, f"set_{axis}ticks")(ticks, [f"{v:g}" for v in ticks])
    ax.minorticks_off()


def shared_own(ax, t, slope):
    """S against D per tissue, with the Theil-Sen line of log D on log S."""
    c = [COLOR[gtex_system(k)] for k in t.index]
    ax.scatter(t.S, t.own, s=24, c=c, edgecolor="white", lw=0.4, zorder=2)
    xs = np.geomspace(t.S.min(), t.S.max(), 50)
    icpt = np.median(np.log(t.own) - slope * np.log(t.S))
    ax.plot(xs, np.exp(icpt + slope * np.log(xs)), color=MUTED, lw=1.2, zorder=1)
    log_ticks(ax, 0.08, 15)
    log_ticks(ax, 0.08, 15, "y")
    ax.set_xlabel("shared with the body, S")
    ax.set_ylabel("own departure, D")


def coupling_curve(ax, t):
    """Coupling against S, with S / (S + D) at the median D."""
    ax.scatter(t.S, t.coupling, s=24, c=[COLOR[gtex_system(k)] for k in t.index], edgecolor="white", lw=0.4, zorder=2)
    xs = np.geomspace(0.08, 15, 80)
    d = t.own.median()
    ax.plot(xs, xs / (xs + d), color=MUTED, lw=1.2, zorder=1)
    ax.text(3, 0.12, f"S / (S + {d:.1f})", ha="center", fontsize=8.5, color=MUTED)
    log_ticks(ax, 0.08, 15)
    ax.set_ylim(-0.05, 1.05)
    ax.set_xlabel("shared with the body, S")
    ax.set_ylabel("coupling to the\nbody's clock state")


def specific(ax, t, ratio):
    """Own variance of the clock genes against that of expression-matched genes, per tissue; `ratio`: median clock /
    matched over tissues."""
    ax.scatter(t.own_matched_median, t.own, s=24, c=[COLOR[gtex_system(k)] for k in t.index], edgecolor="white", lw=0.4,
               zorder=2)
    ax.plot([0.2, 6], [0.2, 6], color=MUTED, lw=0.8, ls="--", zorder=1)
    ax.plot([0.2, 6 / ratio], [0.2 * ratio, 6], color=MUTED, lw=0.8, ls=":", zorder=1)
    ax.text(2.2, 2.2 * 0.8, "1:1", fontsize=8.5, color=MUTED, ha="left", va="top")
    ax.text(1.6, 1.6 * ratio, f"×{ratio:.1f}", fontsize=8.5, color=MUTED, ha="left", va="top")
    log_ticks(ax, 0.2, 6)
    log_ticks(ax, 0.2, 6, "y")
    ax.set_xticks([0.3, 1, 3], ["0.3", "1", "3"])
    ax.set_yticks([0.3, 1, 3], ["0.3", "1", "3"])
    ax.set_xlabel("own departure,\nrandom genes")
    ax.set_ylabel("own departure, clock genes")


def main():
    pf.use_style("slide")
    OUT.mkdir(parents=True, exist_ok=True)
    t = pd.read_csv(EXP / "A294_floor/output/tissues.csv", index_col="tissue")
    slope = pd.read_csv(EXP / "A294_floor/output/tests.csv").set_index("test").loc[
        "descriptive", "slope_log_own_clock_on_log_S"]
    bab = pd.read_csv(PAPER / "results/figure_2/output/organs/r2_couplings.csv").query("species == 'baboon'")
    lay = pf.Layout.slide("ABCD", kind="wide", col_ratios=[1, 1, 0.85, 1])
    shared_own(lay.axes["A"], t, slope)
    legend_in(lay.axes["A"], fontsize=7.5, loc="upper left")
    coupling_curve(lay.axes["B"], t)
    tests = pd.read_csv(EXP / "A294_floor/output/tests.csv").set_index("test")
    specific(lay.axes["C"], t, float(np.exp(tests.loc["H294.1", "median_log_ratio"])))
    for k in "ABC":
        lay.axes[k].set_title("human", fontsize=13, color=INK)
    timed = pd.read_csv(PAPER / "results/figure_1/output/species/a36_tissues.csv", index_col="tissue").g_abs
    bab = bab.set_index("compartment").assign(timed_amp=timed)
    scatter(lay.axes["D"], bab, baboon_system, "baboon", x="timed_amp", xlabel="rhythm amplitude\nagainst known time")
    lay.axes["D"].set_ylabel("coupling")
    fig = lay.finalize()
    print(pf.collisions(fig))
    pf.save(fig, OUT / "organs_health")


if __name__ == "__main__":
    main()
