"""Slide figure (results, R2.1, merged with R2.2): organs share one clock state; weak clocks are the loosely coupled ones,
the brain included.

A. GTEx fast deaths: clock-gene variation (A104 own amplitude, sqrt of shared + own clock-plane variance) against coupling to the body's clock state (peripheral-anchor host, peripheral
   tissues read leave-one-out; A263 R2 fit, `organs_compute.py`), one dot per tissue, coloured by organ system (brain
   regions pink); grey: least-squares line coupling = a + b log(amplitude) over the peripheral tissues. Pituitary (in
   neither A263 group) not shown. The registered statistic is A104's (H104.1 not supported; observed rho +0.734).
B. Baboon (Mure 2018, 12 animals): rhythm amplitude against known time (A36 timed gain |g|, estimated without the host)
   against coupling to the body's clock state (A42 readout 2), 35 peripheral tissues and 22 brain regions; grey: the
   peripheral line. Post hoc panel (A289 H289.3 registered the untimed amplitude: brain on the line). Pineal, pituitary
   and ocular tissues not shown.
C. Is the relation specific to the clock genes? rho(amplitude, coupling) over the GTEx peripheral tissues for the clock
   genes (blue) and 200 expression-matched and 200 stress-responsive gene sets (A271 Part B, H271.2; calibrated: A290).
D. How far a tissue's clock departs from its person's body time (A259, registered): SD across donors [h] with its 90 %
   interval, 11 body tissues. An upper bound: the departure is not above that of the companion gene sets (H259.2).

Data: presentation/output/organs/couplings.csv; results/figure_2/output/organs/r2_couplings.csv;
results/figure_1/output/species/a36_tissues.csv;
exploratory/A271_responsive_companions/output/a104_sets.csv; exploratory/A259_tissue_deviation/output/score/targets.csv.
GTEx aggregates only (one value per tissue). Output: presentation/figures/organs_health/organs_health.{svg,png}.

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
CS = "#1f77b4"
SHORT = {"Adipose - Visceral (Omentum)": "adipose, visceral", "Adipose - Subcutaneous": "adipose, subcut.",
         "Skin - Not Sun Exposed (Suprapubic)": "skin, not sun-exp.", "Skin - Sun Exposed (Lower leg)": "skin, sun-exposed",
         "Heart - Atrial Appendage": "heart, atrium", "Esophagus - Mucosa": "oesophagus mucosa",
         "Muscle - Skeletal": "skeletal muscle", "Artery - Tibial": "tibial artery", "Nerve - Tibial": "tibial nerve",
         "Thyroid": "thyroid", "Lung": "lung"}


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


def specificity(ax, sets):
    s = sets[sets.block == "periphery"]
    clock = s.loc[s.kind == "clock", "rho"].iloc[0]
    for j, k in enumerate(("matched", "responsive")):
        v = s.loc[s.kind == k, "rho"].dropna().to_numpy()
        ax.scatter(v, j + np.random.default_rng(j).uniform(-0.18, 0.18, len(v)), s=6, color=MUTED, alpha=0.5, lw=0)
    ax.axvline(clock, color=CS, lw=2.2)
    ax.annotate("clock\ngenes", (clock, 1.5), xytext=(-4, 0), textcoords="offset points", ha="right", va="center",
                fontsize=9, color=CS)
    ax.axvline(0, color=MUTED, lw=0.8, ls="--")
    ax.set_yticks([0, 1], ["expression-\nmatched", "stress-\nresponsive"])
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(-0.5, 1.9)
    ax.set_xlim(-1, 1)
    ax.set_xlabel("variation-coupling\ncorrelation")


def departure(ax, d):
    d = d.sort_values("hours")
    y = np.arange(len(d))
    c = [COLOR[gtex_system(k)] for k in d.target]
    ax.hlines(y, d.hours_lo, d.hours_hi, colors=c, lw=1.6)
    ax.scatter(d.hours, y, s=24, c=c, zorder=3)
    ax.set_yticks(y, [SHORT.get(k, k) for k in d.target], fontsize=9)
    ax.tick_params(axis="y", length=0)
    ax.set_xlim(0, 4)
    ax.set_xlabel("departure from\nbody time, SD [h]")


def main():
    pf.use_style("slide")
    OUT.mkdir(parents=True, exist_ok=True)
    t = pd.read_csv(PAPER / "presentation/output/organs/couplings.csv", index_col="tissue")
    bab = pd.read_csv(PAPER / "results/figure_2/output/organs/r2_couplings.csv").query("species == 'baboon'")
    sets = pd.read_csv(EXP / "A271_responsive_companions/output/a104_sets.csv")
    dep = pd.read_csv(EXP / "A259_tissue_deviation/output/score/targets.csv")
    lay = pf.Layout.slide("ABCD", kind="wide", col_ratios=[1, 1, 0.7, 0.85])
    scatter(lay.axes["A"], t, gtex_system, "human")
    lay.axes["A"].set_ylabel("coupling to the\nbody's clock state")
    legend_in(lay.axes["A"], fontsize=8, loc="upper left")
    timed = pd.read_csv(PAPER / "results/figure_1/output/species/a36_tissues.csv", index_col="tissue").g_abs
    bab = bab.set_index("compartment").assign(timed_amp=timed)
    scatter(lay.axes["B"], bab, baboon_system, "baboon", x="timed_amp", xlabel="rhythm amplitude\nagainst known time")
    specificity(lay.axes["C"], sets)
    departure(lay.axes["D"], dep)
    fig = lay.finalize()
    print(pf.collisions(fig))
    pf.save(fig, OUT / "organs_health")


if __name__ == "__main__":
    main()
