"""Slide figure (results, 3.4; formerly R2.1): coupling mostly reads loudness. Each organ expresses the body's time with
its own gain over a departure of nearly constant size; coupling follows, the brain included.

Notation (A104 table, GTEx fast deaths, peripheral tissues read leave-one-out on the peripheral host): S = z2 V_h, the
clock-plane variance a tissue shares with the body's clock state (gain squared times host spread); D = 2 psi, its clock-plane variance
at a fixed body state, net of counting noise (own departure in amplitude and phase); coupling c = S / (S + D).
Order: the observation (A, B: louder clocks are more coupled, in human and baboon), then why (C: the departure is
nearly constant), then the departure's specificity (D).
A. Coupling against the shared amplitude sqrt(S), one dot per tissue, coloured by organ system; grey: c = S / (S + D)
   at the median D over tissues (the definition of c drawn at a constant D, not a fit; explained by C).
B. Baboon (Mure 2018, 12 animals): rhythm amplitude against known time (A36 timed gain |g|, estimated without the host)
   against coupling (A42 readout 2), 35 peripheral tissues and 22 brain regions; grey: the peripheral least-squares line.
   Post hoc panel (A289 H289.3 registered the untimed amplitude: brain on the line). Pineal, pituitary and ocular tissues
   not shown.
C. Shared against own clock variation (log axes, same range); grey: Theil-Sen line of log D on log S (A294b descriptive
   slope). Post hoc display (A291 post hoc table).
D. Is the own variation specific to the clock genes? Own variance of the clock genes against the median of 100
   expression-matched gene sets in the same tissue (A294b, registered, H294b.1 pass); dashed: identity; dotted: median
   ratio. The core-loop values (A329) are not shown: the 11-gene projector has 3.1 times the noise gain, so its ratio is
   not comparable (A329 erratum).

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


def specific(ax, t, ratio):
    """Own variance of the clock genes against that of expression-matched genes, per tissue; `ratio`: median clock /
    matched over tissues."""
    ax.scatter(t.own_matched_median, t.own, s=24, c=[COLOR[gtex_system(k)] for k in t.index], edgecolor="white", lw=0.4,
               zorder=2)
    lo, hi = 0.25, 8
    ax.plot([lo, hi], [lo, hi], color=MUTED, lw=0.8, ls="--", zorder=1)
    ax.plot([lo, hi / ratio], [lo * ratio, hi], color=MUTED, lw=0.8, ls=":", zorder=1)
    ax.text(4, 4 * 0.85, "1:1", fontsize=8.5, color=MUTED, ha="left", va="top")
    ax.text(hi / ratio * 1.05, hi * 0.95, f"\u00d7{ratio:.1f}", fontsize=8.5, color=MUTED, ha="left", va="top")
    log_ticks(ax, lo, hi)
    log_ticks(ax, lo, hi, "y")
    ax.set_xticks([0.3, 1, 3], ["0.3", "1", "3"])
    ax.set_yticks([0.3, 1, 3], ["0.3", "1", "3"])
    ax.set_xlabel("own departure,\nexpression-matched genes")
    ax.set_ylabel("own departure, clock genes")


def main():
    pf.use_style("slide")
    OUT.mkdir(parents=True, exist_ok=True)
    t = pd.read_csv(EXP / "A294_floor/output/tissues.csv", index_col="tissue")
    slope = pd.read_csv(EXP / "A294_floor/output/tests.csv").set_index("test").loc[
        "descriptive", "slope_log_own_clock_on_log_S"]
    bab = pd.read_csv(PAPER / "results/figure_2/output/organs/r2_couplings.csv").query("species == 'baboon'")
    lay = pf.Layout.slide("ABCD", kind="wide", col_ratios=[1, 1, 1, 0.85])
    coupling_curve(lay.axes["A"], t)
    legend_in(lay.axes["A"], fontsize=7.5, loc="upper left")
    timed = pd.read_csv(PAPER / "results/figure_1/output/species/a36_tissues.csv", index_col="tissue").g_abs
    bab = bab.set_index("compartment").assign(timed_amp=timed)
    scatter(lay.axes["B"], bab, baboon_system, "baboon", x="timed_amp", xlabel="rhythm amplitude\nagainst known time")
    lay.axes["B"].set_ylabel("coupling")
    lay.axes["B"].set_xlabel("amplitude against\nknown time")
    shared_own(lay.axes["C"], t, slope)
    tests = pd.read_csv(EXP / "A294_floor/output/tests.csv").set_index("test")
    specific(lay.axes["D"], t, float(np.exp(tests.loc["H294.1", "median_log_ratio"])))
    print("S fold", round(t.S.max() / t.S.min(), 1), "gain fold", round(np.sqrt(t.S.max() / t.S.min()), 2),
          "D fold", round(t.own.max() / t.own.min(), 2), "median D", round(t.own.median(), 2))
    for k in "ACD":
        lay.axes[k].set_title("human", fontsize=13, color=INK)
    fig = lay.finalize()
    print(pf.collisions(fig))
    pf.save(fig, OUT / "organs_health")


if __name__ == "__main__":
    main()
