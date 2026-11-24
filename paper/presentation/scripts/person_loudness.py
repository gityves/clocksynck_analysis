"""Slide figure (results, adult health, 3.5): persons differ in a body-wide expression scale; in deaths the clock
follows it most reliably (revised 2026-11-24 after A361-A371; living panel removed, A331 withdrawn by A366).

A. How reliably loudness is shared across a person's tissues (leverage: correlation of the clock-plane radii of two
   halves of a person's tissues after the phase harmonics are removed), core loop: all GTEx fast deaths (A336 part 2,
   244 donors, 27 tissues; registered) and violent, instantaneous deaths (Hardy 1, A367, 34 donors; registered, scored
   at power 0.775 by the user's decision). Blue: clock genes; dots: 200 expression-matched (light grey) and 200
   stress-responsive (dark grey) gene sets on the same donors.
B. How large the shared loudness is (cv2: shared squared coefficient of variation of the radii), same sets, all fast
   deaths, core loop (post hoc size note): the clock is not above random genes.
C. Elasticity of each tissue's amplitude to the person's shared amplitude, relative to the other tissues' composite
   (A312, registered), with 90 % interval, against the tissue's gain sqrt(S). Grey band: 10-90 % of the estimates
   when every organ scales exactly in proportion (post hoc calibration, `A312_elasticity/scripts/calib_level.py`).
D. What tracks the person's scale (A360, registered; H360.1 pass, H360.2 fail): standardised coefficients of the log
   host radius (14 genes) with donor-bootstrap 90 % intervals, 237 fast-death donors.

Data: exploratory/A336_two_readouts/output/{person_gtex,person_gtex_companions}.csv;
exploratory/A367_hardy_agonal/output/{tests,companions}.csv; exploratory/A312_elasticity/output/{tissues,calib_level}.csv;
exploratory/A358_A360_gain_sources/output/a360_coefficients.csv. GTEx aggregates only.
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


def _dots(ax, i, comp, col_name):
    for kind, col, dy in (("matched", MATCHED, -0.13), ("responsive", RESPONSIVE, 0.13)):
        v = comp[comp.kind == kind][col_name].to_numpy()
        y = i + dy + np.random.default_rng(i).uniform(-0.07, 0.07, len(v))
        ax.scatter(v, y, s=6, color=col, alpha=0.6, edgecolor="none", zorder=1)


def _rows(ax, labels):
    ax.set_yticks(range(len(labels)), labels)
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(len(labels) - 0.5, -0.5)


def shared(ax):
    g = pd.read_csv(EXP / "A336_two_readouts/output/person_gtex.csv").set_index("gene_set").loc["core"]
    gc = pd.read_csv(EXP / "A336_two_readouts/output/person_gtex_companions.csv").query("gene_set == 'core'")
    h = pd.read_csv(EXP / "A367_hardy_agonal/output/tests.csv").set_index("gene_set").loc["core"]
    hc = pd.read_csv(EXP / "A367_hardy_agonal/output/companions.csv")
    _dots(ax, 0, gc, "leverage")
    _dots(ax, 1, hc, "leverage_h1")
    for i, v in enumerate((g.leverage, h.lev_h1)):
        ax.plot([v] * 2, [i - 0.32, i + 0.32], color=CS, lw=3, zorder=3)
    _rows(ax, ["all fast deaths", "violent deaths"])
    ax.set_xlim(-0.05, 1)
    ax.set_xlabel("how reliably loudness is shared\nacross a person's tissues [r]")
    hd = [Line2D([], [], color=CS, lw=3, label="clock genes"),
          Line2D([], [], ls="", marker="o", ms=4, color=MATCHED, label="matched genes"),
          Line2D([], [], ls="", marker="o", ms=4, color=RESPONSIVE, label="stress-responsive genes")]
    ax.legend(handles=hd, loc="lower left", frameon=False, fontsize=8, handletextpad=0.5, handlelength=1.2)
    return g, h


def size(ax):
    g = pd.read_csv(EXP / "A336_two_readouts/output/person_gtex.csv").set_index("gene_set").loc["core"]
    gc = pd.read_csv(EXP / "A336_two_readouts/output/person_gtex_companions.csv").query("gene_set == 'core'")
    _dots(ax, 0, gc, "cv2")
    ax.plot([g.cv2] * 2, [-0.32, 0.32], color=CS, lw=3, zorder=3)
    _rows(ax, [""])
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("how large the shared\nloudness is [cv²]")
    return g, gc


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
    ax.set_ylabel("follows the person's scale\n[elasticity, relative]")
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
    ax.set_xlabel("effect on the person's\nscale [SD per SD]")
    return co


def main():
    pf.use_style("slide")
    OUT.mkdir(parents=True, exist_ok=True)
    lay = pf.Layout.slide("ABCD", kind="wide", col_ratios=[1.1, 0.7, 1, 1])
    g, h = shared(lay.axes["A"])
    gs, gcs = size(lay.axes["B"])
    t = elasticity(lay.axes["C"])
    co = tracks(lay.axes["D"])
    fig = lay.finalize()
    print("GTEx core", round(g.leverage, 3), "matched q95", round(g.matched_q95, 3), "responsive q95", round(g.responsive_q95, 3))
    print("Hardy 1 core", round(h.lev_h1, 3), "matched q95", round(h.matched_q95_h1, 3), "responsive q95", round(h.responsive_q95_h1, 3))
    print("cv2 core", round(gs.cv2, 3), gcs.groupby("kind").cv2.median().round(3).to_dict())
    print("elasticity median", round(t.elasticity.median(), 3), "range", round(t.elasticity.min(), 2), round(t.elasticity.max(), 2))
    print(co.round(3).to_string())
    print(pf.collisions(fig))
    pf.save(fig, OUT / "person_loudness")


if __name__ == "__main__":
    main()
