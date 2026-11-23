"""Slide figure (results, before birth, 2.2): is the fetal clock-shaped state a clock? The core loop covaries on its own,
the states are not one slow factor, and the pre-mRNA lead excludes fast pulsatile input.

A. Clock-gene covariation between individuals C, 14 clock genes (ink) against the core loop alone (blue; 11 genes, without
   PER1, PER2, NR1D1), per fetal unit: whole adrenal, del Valle 2023 (A328 part a, BLIND; H328.1, H328.2), adrenal cortex
   cells (A244 / A266), liver and lung (A243 / A269). Filled: significant and above the 95th percentile of 200
   expression-matched and 200 stress-responsive gene sets; open: not. Lines across each bar: 95th percentile of the
   matched (light grey) and stress-responsive (dark grey, white edge) gene sets.
B. Elongation of the fetal clock-plane states, iota (0: a circle, the states spread evenly around the clock; 1: a line,
   one factor moving individuals back and forth), point and 90 % bootstrap interval within cohort (A340, BLIND, H340.1),
   for the three scored units (adrenal: whole adrenal, del Valle). References from the power stage, noise-matched to each unit at 0.75 of the adult clock:
   blue band, a daily clock sampled in hospital hours (08-18 h, 5-95 %); grey band, no clock (genes permuted, 5-95 %);
   dark line, a single slow factor (5th percentile of the line reference; above it = compatible with one line).
C. Pre-mRNA lead over mRNA [h] (unspliced over spliced clock-plane state, matched estimator, point and 90 % bootstrap
   interval; A328 part b, BLIND), for the same units and the running postnatal liver clock (grey; 12 libraries). A
   first-order filter (mRNA follows pre-mRNA with lifetime tau) driven at period T gives the lead
   (24 / 2 pi) arctan(2 pi tau / T) on the 24 h clock: with tau = 1 h (A339) 0.98 h for a daily input (dashed) and
   >= 3.09 h for an input of period <= 6 h (red band: fast, pulsatile). Minimum period T_min = 2 pi tau / tan(phi_hi)
   from the interval's upper bound. POST HOC: a combination of registered elements (A338 fast-input threshold, A339
   lifetime; the T_min companion of A339 is NOT BLIND). The matched estimator reads leads about 2-fold low (Greenwell,
   A337), so T_min may be overestimated up to 2-fold. Bounds written to presentation/output/fetal_clock/lead_bounds.csv.

Status: A243, A244 registered, blind; A266, A269 registered, not blind for the earlier fetal results; A328 part a and
A340 registered, BLIND. Aggregates only.
Data: exploratory/A328_fetal_running_loop/output/score_a/l1.csv; exploratory/A244_fetal_adrenal_classes/output/{l1,verdicts_a244}.csv;
exploratory/A266_fetal_adrenal_drive/output/l1.csv; exploratory/A243_fetal_organs_pooled/output/l1.csv;
exploratory/A269_fetal_organs_drive/output/l1.csv; exploratory/A340_fetal_geometry/output/score/verdicts.csv;
exploratory/A328_fetal_running_loop/output/score_b/verdicts.csv.
Output: presentation/figures/fetal_clock/fetal_clock.{svg,png} (wide); presentation/output/fetal_clock/lead_bounds.csv.

    python scripts/fetal_clock.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

import paperfig as pf
from matplotlib import patheffects

PAPER = Path(__file__).resolve().parents[2]
EXP = PAPER / "exploratory"
OUT = PAPER / "presentation" / "figures" / "fetal_clock"
TAB_OUT = PAPER / "presentation" / "output" / "fetal_clock"
TAU, FAST = 1.0, 6.0  # mRNA lifetime [h] (A339); fast-input period threshold [h] (A338)
CORE, ALL, INK = sns.color_palette().as_hex()[0], "#222222", "#222222"
MATCHED, RESPONSIVE, NOCLOCK = "#BBBBBB", "#6E6E6E", "#DDDDDD"
FASTC = sns.color_palette().as_hex()[3]  # fast, pulsatile input
GENES = {"all": ("14 clock genes", ALL, -0.19), "core": ("core loop alone", CORE, 0.19)}
BAR = 0.33
EDGE = [patheffects.Stroke(linewidth=3.6, foreground="white"), patheffects.Normal()]


def covariation():
    """Rows: unit, genes (all / core), C, ok (clock-specific), matched and responsive q95."""
    def row(unit, genes, r, ok=None):
        return dict(unit=unit, genes=genes, C=r.C, q95=r.matched_q95, rq95=r.responsive_q95,
                    ok=bool(r.clock_specific if ok is None else ok))
    dv = pd.read_csv(EXP / "A328_fetal_running_loop/output/score_a/l1.csv").set_index("test")
    a244 = EXP / "A244_fetal_adrenal_classes/output"
    cx = pd.read_csv(a244 / "l1.csv").set_index("test").loc["cortex"]
    cx_ok = pd.read_csv(a244 / "verdicts_a244.csv").set_index("test").loc["cortex", "verdict"] == "PASS"
    a266 = pd.read_csv(EXP / "A266_fetal_adrenal_drive/output/l1.csv").set_index("test")
    a243 = pd.read_csv(EXP / "A243_fetal_organs_pooled/output/l1.csv").set_index("test")
    a269 = pd.read_csv(EXP / "A269_fetal_organs_drive/output/l1.csv").set_index("test")
    rows = [row("whole\nadrenal,\nblind", "all", dv.loc["dv:all"]), row("whole\nadrenal,\nblind", "core", dv.loc["dv:core"]),
            row("adrenal\ncortex", "all", cx, cx_ok), row("adrenal\ncortex", "core", a266.loc["cortex:core"])]
    for o in ("liver", "lung"):
        rows += [row(o, "all", a243.loc[f"{o}:all"]), row(o, "core", a269.loc[f"{o}:core"])]
    return pd.DataFrame(rows)


def panel_a(ax, d):
    units = list(dict.fromkeys(d.unit))
    for _, r in d.iterrows():
        _, c, dx = GENES[r.genes]
        x = units.index(r.unit) + dx
        ax.bar([x], [r.C], width=BAR, color=c if r.ok else "white", edgecolor=c, lw=1.4, zorder=2)
        ax.plot([x - BAR / 2, x + BAR / 2], [r.q95] * 2, color=MATCHED, lw=2.2, solid_capstyle="butt", zorder=3)
        ax.plot([x - BAR / 2, x + BAR / 2], [r.rq95] * 2, color=RESPONSIVE, lw=2.2, solid_capstyle="butt", zorder=3,
                path_effects=EDGE)
    ax.axhline(0, color=pf.GREY, lw=0.6, zorder=1)
    ax.set_xticks(range(len(units)), units)
    ax.tick_params(axis="x", length=0)
    ax.set_xlim(-0.6, len(units) - 0.4)
    ax.set_ylim(0, 1)
    ax.set_ylabel("clock-gene covariation\nbetween individuals, C")


def panel_b(ax):
    v = pd.read_csv(EXP / "A340_fetal_geometry/output/score/verdicts.csv")
    lab = {("delvalle", "all"): "adrenal,\n14 genes", ("delvalle", "core"): "adrenal,\ncore loop",
           ("liver", "all"): "liver,\n14 genes"}
    for i, r in enumerate(v.itertuples()):
        ax.fill_between([i - 0.32, i + 0.32], r.daily_q05, r.daily_q95, color=CORE, alpha=0.18, lw=0, zorder=1)
        ax.fill_between([i - 0.32, i + 0.32], r.noclock_q05, r.noclock_q95, color=NOCLOCK, alpha=0.8, lw=0, zorder=0)
        ax.plot([i - 0.32, i + 0.32], [r.line_q05] * 2, color="#555555", lw=2, zorder=2)
        ax.plot([i, i], [r.iota_lo, r.iota_hi], color=INK, lw=1.8, zorder=3)
        ax.scatter([i], [r.iota], s=46, color=INK, zorder=4)
    ax.set_xticks(range(len(v)), [lab[u, w] for u, w in zip(v.unit, v.waveform)])
    ax.tick_params(axis="x", length=0)
    ax.set_xlim(-0.6, len(v) - 0.4)
    ax.set_ylim(0, 1.02)
    ax.set_ylabel("elongation of the states\n(0 circle, 1 line)")
    return v


def lead_h(T, tau=TAU):
    """Lead of pre-mRNA over mRNA on the 24 h clock [h] for an input of period T [h] (first-order filter)."""
    return 24 / (2 * np.pi) * np.arctan(2 * np.pi * tau / T)


def lead_bounds():
    """Fetal leads (matched estimator) with the minimum input period implied by each upper bound; written to a table."""
    b = pd.read_csv(EXP / "A328_fetal_running_loop/output/score_b/verdicts.csv").set_index(["unit", "waveform"])
    keys = [("delvalle", "all"), ("delvalle", "core"), ("liver", "all"), ("postnatal_liver", "all")]
    rows = []
    for u, w in keys:
        r = b.loc[(u, w)]
        phi = 2 * np.pi * r.hi_matched / 24
        rows.append(dict(unit=u, waveform=w, lead_h=r.delta_matched, lo_h=r.lo_matched, hi_h=r.hi_matched,
                         gate_p=r.gate_p_matched, T_min_h=2 * np.pi * TAU / np.tan(phi) if phi > 0 else np.inf,
                         lead_daily_h=lead_h(24), lead_fast_h=lead_h(FAST)))
    t = pd.DataFrame(rows)
    TAB_OUT.mkdir(parents=True, exist_ok=True)
    t.to_csv(TAB_OUT / "lead_bounds.csv", index=False)
    return t


def panel_c(ax, t):
    lab = {("delvalle", "all"): "adrenal,\n14 genes", ("delvalle", "core"): "adrenal,\ncore loop",
           ("liver", "all"): "liver,\n14 genes", ("postnatal_liver", "all"): "liver after\nbirth"}
    top = 4.5
    ax.fill_between([-0.6, len(t) - 0.4], t.lead_fast_h.iloc[0], top, color=FASTC, alpha=0.15, lw=0, zorder=0)
    ax.axhline(t.lead_daily_h.iloc[0], color="#555555", lw=1.2, ls="--", zorder=1)
    ax.axhline(0, color=pf.GREY, lw=0.6, zorder=1)
    for i, r in enumerate(t.itertuples()):
        c = pf.GREY if r.unit.startswith("postnatal") else INK
        ax.plot([i, i], [r.lo_h, r.hi_h], color=c, lw=1.8, zorder=3)
        ax.scatter([i], [r.lead_h], s=46, color=c, zorder=4)
    ax.set_xticks(range(len(t)), [lab[u, w] for u, w in zip(t.unit, t.waveform)])
    ax.tick_params(axis="x", length=0)
    ax.set_xlim(-0.6, len(t) - 0.4)
    ax.set_ylim(-1.5, top)
    ax.set_ylabel("pre-mRNA lead over mRNA [h]")


def main():
    pf.use_style("slide")
    OUT.mkdir(parents=True, exist_ok=True)
    lay = pf.Layout.slide("ABC", kind="wide", col_ratios=[1.25, 1, 1.1])
    d = covariation()
    panel_a(lay.axes["A"], d)
    v = panel_b(lay.axes["B"])
    t = lead_bounds()
    panel_c(lay.axes["C"], t)
    h = [Patch(facecolor=ALL, edgecolor=ALL, label="14 clock genes"),
         Patch(facecolor=CORE, edgecolor=CORE, label="core loop alone"),
         Patch(facecolor="white", edgecolor=INK, label="not clock-specific"),
         Line2D([], [], color=MATCHED, lw=2.2, label="matched genes, 95th pct"),
         Line2D([], [], color=RESPONSIVE, lw=2.2, path_effects=EDGE, label="stress-responsive genes, 95th pct"),
         Line2D([], [], ls="", marker="o", ms=7, color=INK, label="fetal states, 90 % CI"),
         Patch(facecolor=CORE, alpha=0.18, lw=0, label="daily clock, hospital hours"),
         Patch(facecolor=NOCLOCK, alpha=0.8, lw=0, label="no clock"),
         Line2D([], [], color="#555555", lw=2, label="one slow factor (lower bound)"),
         Line2D([], [], color="#555555", lw=1.2, ls="--", label="expected lead, daily input"),
         Patch(facecolor=FASTC, alpha=0.15, lw=0, label="expected lead, input period ≤ 6 h"),
         Line2D([], [], ls="", marker="o", ms=7, color=pf.GREY, label="running clock after birth")]
    lay.shared_legend("ABC", h, [x.get_label() for x in h], side="bottom", ncol=4, fontsize=10, handletextpad=0.4,
                      columnspacing=1.0)
    fig = lay.finalize()
    print(d.round(3).to_string(index=False))
    print(v[["unit", "waveform", "iota", "iota_lo", "iota_hi", "line_q05", "daily_q05", "daily_q95", "noclock_q05",
             "noclock_q95", "holm", "H340_1"]].round(3).to_string(index=False))
    print(t.round(2).to_string(index=False))
    print(pf.collisions(fig))
    pf.save(fig, OUT / "fetal_clock")


if __name__ == "__main__":
    main()
