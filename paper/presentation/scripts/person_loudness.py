"""Slide figure (results, adult health, 3.5): each person has a clock volume that every organ follows.

A. Loudness shared across a person's tissues (leverage: correlation of the clock-plane radii of a person's tissues, or
   tissue halves, after the phase harmonics are removed): GTEx fast deaths (A336 part 2, 244 donors, 27 tissues, 14
   genes; registered) and living women (A331, MuTHER fat and skin from one biopsy, 680 women; registered, H331.1 pass).
   Blue: clock genes; dots: 200 (GTEx) or 200 (MuTHER) expression-matched (light grey) and stress-responsive (dark
   grey) gene sets on the same persons.
B. Elasticity of each tissue's amplitude to the person's shared amplitude, relative to the other tissues' composite
   (A312, registered), with 90 % interval, against the tissue's gain sqrt(S). Grey band: 10-90 % of the estimates
   when every organ scales exactly in proportion (post hoc calibration, `A312_elasticity/scripts/calib_level.py`):
   the estimator reads about 1.1 for a true 1.
C. What tracks the person volume (A360, registered; H360.1 pass, H360.2 fail): standardised coefficients of the log
   host radius with donor-bootstrap 90 % intervals, 237 fast-death donors.

Data: exploratory/A336_two_readouts/output/{person_gtex,person_gtex_companions}.csv;
exploratory/P1_person_loudness/output/{tests,companions}.csv; exploratory/A312_elasticity/output/tissues.csv;
exploratory/A358_A360_gain_sources/output/a360_coefficients.csv; exploratory/A312_elasticity/output/calib_level.csv. GTEx
aggregates only.
Output: presentation/figures/person_loudness/person_loudness.{svg,png} (wide).

    python scripts/person_loudness.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D

import paperfig as pf
from _systems import COLOR, gtex_system

PAPER = Path(__file__).resolve().parents[2]
EXP = PAPER / "exploratory"
OUT = PAPER / "presentation" / "figures" / "person_loudness"
INK, MUTED = "#222222", pf.GREY
CS = sns.color_palette().as_hex()[0]
MATCHED, RESPONSIVE = "#BBBBBB", "#6E6E6E"
TERMS = {"gc": "glucocorticoid-target\nstate", "ieg": "immediate-early\nstate", "rin": "RIN", "ischaemic_h": "ischaemic time",
         "n_tissues": "tissues sampled", "tod_cos": "time of death (cos)", "tod_sin": "time of death (sin)"}


def shared(ax):
    g = pd.read_csv(EXP / "A336_two_readouts/output/person_gtex.csv").set_index("gene_set").loc["14 genes"]
    gc = pd.read_csv(EXP / "A336_two_readouts/output/person_gtex_companions.csv").query("gene_set == '14 genes'")
    m = pd.read_csv(EXP / "P1_person_loudness/output/tests.csv").set_index("test").loc["H331.1"]
    mc = pd.read_csv(EXP / "P1_person_loudness/output/companions.csv")
    rows = [("GTEx deaths\n(27 tissues)", g.leverage, gc), ("living women\n(fat, skin)", m.leverage, mc)]
    for i, (_, clock, comp) in enumerate(rows):
        for kind, col, dy in (("matched", MATCHED, -0.13), ("responsive", RESPONSIVE, 0.13)):
            v = comp[comp.kind == kind].leverage.to_numpy()
            y = i + dy + np.random.default_rng(i).uniform(-0.07, 0.07, len(v))
            ax.scatter(v, y, s=6, color=col, alpha=0.6, edgecolor="none", zorder=1)
        ax.plot([clock] * 2, [i - 0.32, i + 0.32], color=CS, lw=3, zorder=3)
    ax.set_yticks(range(len(rows)), [r[0] for r in rows])
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(len(rows) - 0.5, -0.5)
    ax.set_xlim(-0.05, 1)
    ax.set_xlabel("loudness shared across\na person's tissues [r]")
    h = [Line2D([], [], color=CS, lw=3, label="clock genes"),
         Line2D([], [], ls="", marker="o", ms=4, color=MATCHED, label="matched genes"),
         Line2D([], [], ls="", marker="o", ms=4, color=RESPONSIVE, label="stress-responsive genes")]
    ax.legend(handles=h, loc="lower right", frameon=False, fontsize=8, handletextpad=0.5, handlelength=1.2)
    return g, m


def elasticity(ax):
    t = pd.read_csv(EXP / "A312_elasticity/output/tissues.csv")
    x = np.sqrt(t.S)
    c = [COLOR[gtex_system(k)] for k in t.tissue]
    ax.errorbar(x, t.elasticity, yerr=[t.elasticity - t.lo90, t.hi90 - t.elasticity], fmt="none", ecolor=MUTED,
                elinewidth=0.7, zorder=1)
    ax.scatter(x, t.elasticity, s=24, c=c, edgecolor="white", lw=0.4, zorder=2)
    ax.axhline(1, color=MUTED, lw=0.8, ls="--", zorder=0)
    cal = pd.read_csv(EXP / "A312_elasticity/output/calib_level.csv").elasticity
    ax.axhspan(cal.quantile(0.1), cal.quantile(0.9), color=MUTED, alpha=0.12, lw=0, zorder=0)
    ax.text(1.0, 0.56, "every organ in proportion\n(simulated, as estimated)", ha="center", va="top",
            fontsize=8, color=MUTED)
    ax.set_xscale("log")
    ax.set_xticks([0.5, 1, 2, 3], ["0.5", "1", "2", "3"])
    ax.minorticks_off()
    ax.set_xlabel("organ gain, √S")
    ax.set_ylabel("follows the person's volume\n[elasticity, relative]")
    return t


def tracks(ax):
    co = pd.read_csv(EXP / "A358_A360_gain_sources/output/a360_coefficients.csv").query("model == 'registered'")
    co = co.set_index("term").loc[list(TERMS)]
    y = np.arange(len(co))
    sig = (co.lo90 > 0) | (co.hi90 < 0)
    ax.hlines(y, co.lo90, co.hi90, color=[INK if s else MUTED for s in sig], lw=1.6)
    ax.scatter(co.std_coef, y, s=30, color=[INK if s else MUTED for s in sig], zorder=3)
    ax.axvline(0, color=MUTED, lw=0.8, ls="--", zorder=0)
    ax.set_yticks(y, [TERMS[k] for k in co.index])
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(len(co) - 0.5, -0.5)
    ax.set_xlim(-0.45, 0.35)
    ax.set_xlabel("effect on the person's\nvolume [SD per SD]")
    return co


def main():
    pf.use_style("slide")
    OUT.mkdir(parents=True, exist_ok=True)
    lay = pf.Layout.slide("ABC", kind="wide", col_ratios=[1.1, 1, 1])
    g, m = shared(lay.axes["A"])
    t = elasticity(lay.axes["B"])
    co = tracks(lay.axes["C"])
    fig = lay.finalize()
    print("GTEx", round(g.leverage, 3), "matched q95", round(g.matched_q95, 3), "responsive q95", round(g.responsive_q95, 3))
    print("MuTHER", round(m.leverage, 3), "matched q95", round(m.matched_q95, 3), "responsive q95", round(m.responsive_q95, 3))
    print("elasticity median", round(t.elasticity.median(), 3), "range", round(t.elasticity.min(), 2), round(t.elasticity.max(), 2))
    print(co.round(3).to_string())
    print(pf.collisions(fig))
    pf.save(fig, OUT / "person_loudness")


if __name__ == "__main__":
    main()
