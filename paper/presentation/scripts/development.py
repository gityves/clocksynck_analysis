"""Slide figure (results, development): before birth, some organs carry a clock-shaped state shared between individuals.

A. Clock-gene covariation between individuals (C, signed coherence of the 14 clock genes against the adult waveform,
   cohort / batch and age adjusted) per organ, human fetus (adrenal cortex A244; lung, liver, heart, kidney A243, pooled
   open cohorts, 4-20 wpc) and late mouse embryo (A309, Qiu 2024, E15-E18.75: adrenocortical cells, hepatocytes, lung
   progenitors, ventricular cardiomyocytes, nephron progenitors). Filled: significant (Holm) and above the 95th percentile
   of 200 expression-matched and 200 stress-responsive gene sets; open: not. Open diamond: median C of adult tissue
   subsampled to the same design (full-size clock). Grey tick: 95th percentile of the expression-matched sets.
B-E. Pairs of a hormone-induced clock gene (PER1, PER2, NR1D1) with a core clock gene (33 pairs): correlation between
   individuals against cos(adult phase difference), least-squares line through the origin; human adrenal cortex, lung,
   liver (A307 rows), mouse hepatocytes (A309). Filled = W (mean of r x expected) significant and above the matched q95.

Status: A243, A244 registered; A307, A309 registered, NOT BLIND for the earlier fetal results. Aggregates only.
Data: exploratory/A243_fetal_organs_pooled/output/{l1,verdicts,position_in_adult}.csv;
exploratory/A244_fetal_adrenal_classes/output/{l1,verdicts_a244,position}.csv;
exploratory/A309_mouse_embryo_organs/output/{l1,power_passrates,gc_pattern}.csv; exploratory/A307_fetal_gc_pattern/output/tests.csv;
presentation/output/development/pairs.csv (scripts/development_compute.py).
Output: presentation/figures/development/development.{svg,png}.

    python scripts/development.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D

import paperfig as pf

PAPER = Path(__file__).resolve().parents[2]
EXP = PAPER / "exploratory"
OUT = PAPER / "presentation" / "figures" / "development"
TAB = sns.color_palette().as_hex()
SPECIES = {"human": TAB[0], "mouse": TAB[1]}  # human fetus, mouse embryo
MATCHED = "#BBBBBB"
ORGANS = ["adrenal cortex", "liver", "lung", "heart", "kidney"]
DY = {"human": -0.17, "mouse": 0.17}


def covariation():
    """One row per organ x species: C, pass, matched q95, adult median C."""
    a243 = EXP / "A243_fetal_organs_pooled/output"
    l1 = pd.read_csv(a243 / "l1.csv").set_index("test")
    v = pd.read_csv(a243 / "verdicts.csv").set_index("test")
    adult = pd.read_csv(a243 / "position_in_adult.csv").query("setting == 'a1'").set_index("test")
    rows = [dict(organ=o, species="human", C=l1.loc[f"{o}:all", "C"], ok=v.loc[f"{o}:all", "verdict"] == "PASS",
                 q95=l1.loc[f"{o}:all", "matched_q95"], adult=adult.loc[f"{o}:all", "adult_median_C"])
            for o in ("liver", "lung", "heart", "kidney")]
    a244 = EXP / "A244_fetal_adrenal_classes/output"
    c = pd.read_csv(a244 / "l1.csv").set_index("test").loc["cortex"]
    rows.append(dict(organ="adrenal cortex", species="human", C=c.C, q95=c.matched_q95,
                     ok=pd.read_csv(a244 / "verdicts_a244.csv").set_index("test").loc["cortex", "verdict"] == "PASS",
                     adult=pd.read_csv(a244 / "position.csv").set_index("test").loc["cortex", "adult_median"]))
    a309 = EXP / "A309_mouse_embryo_organs/output"
    m = pd.read_csv(a309 / "l1.csv").set_index("readout")
    pw = pd.read_csv(a309 / "power_passrates.csv").query("setting == 'a1'").set_index("readout")
    for r, o in (("adrenal", "adrenal cortex"), ("liver", "liver"), ("lung", "lung"), ("heart", "heart"),
                 ("kidney", "kidney")):
        rows.append(dict(organ=o, species="mouse", C=m.loc[r, "C"], ok=m.loc[r, "verdict"] == "PASS",
                         q95=m.loc[r, "matched_q95"], adult=pw.loc[r, "median_C"]))
    return pd.DataFrame(rows)


def pattern():
    """Pair-level W: (signal-induced, core) pairs, observed r and the adult expectation; pass flags from A307 / A309."""
    d = pd.read_csv(PAPER / "presentation/output/development/pairs.csv")
    h = pd.read_csv(EXP / "A307_fetal_gc_pattern/output/tests.csv").query("statistic == 'W'").set_index("organ")["pass"]
    m = pd.read_csv(EXP / "A309_mouse_embryo_organs/output/gc_pattern.csv").query("statistic == 'W'").set_index("readout")["pass"]
    ok = {("human", o): bool(v) for o, v in h.items()} | {("mouse", o): bool(v) for o, v in m.items()}
    return d.assign(ok=[ok[s, o] for s, o in zip(d.species, d.organ)])


def mark(ax, x, y, species, ok, q95, adult=None):
    c = SPECIES[species]
    ax.plot([q95, q95], [y - 0.13, y + 0.13], color=MATCHED, lw=2.2, solid_capstyle="butt", zorder=1)
    if adult is not None:
        ax.scatter([adult], [y], s=34, marker="D", facecolor="white", edgecolor=c, lw=1.1, zorder=2)
    ax.scatter([x], [y], s=46, color=c if ok else "white", edgecolor=c, lw=1.4, zorder=3)


def rows(ax, n):
    ax.axvline(0, color=pf.GREY, lw=0.6, zorder=0)
    ax.set_yticks(range(n), ORGANS[:n])
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(n - 0.5, -0.5)


def panel_a(ax, d):
    for _, r in d.iterrows():
        mark(ax, r.C, ORGANS.index(r.organ) + DY[r.species], r.species, r.ok, r.q95, r.adult)
    rows(ax, len(ORGANS))
    ax.set_xlim(-0.6, 1.0)
    ax.set_xlabel("clock-gene covariation\nbetween individuals, C")


def panel_pairs(ax, d, species, organ, title, xlab, ylab):
    x = d[(d.species == species) & (d.organ == organ)]
    c, ok = SPECIES[species], bool(x.ok.iloc[0])
    ax.axhline(0, color=pf.GREY, lw=0.6, zorder=0)
    ax.axvline(0, color=pf.GREY, lw=0.6, zorder=0)
    ax.scatter(x.expected, x.r, s=18, color=c if ok else "white", edgecolor=c, lw=0.9, zorder=2)
    b = np.sum(x.expected * x.r) / np.sum(x.expected ** 2)
    ax.plot([-1, 1], [-b, b], color=c, lw=1.6, zorder=3)
    ax.set_xlim(-1.1, 1.1)
    ax.set_ylim(-1, 1)
    ax.set_xticks([-1, 0, 1])
    ax.set_yticks([-0.5, 0, 0.5])
    ax.set_title(title, fontsize=13)
    ax.set_xlabel(xlab)
    ax.set_ylabel(ylab)
    if not xlab:
        ax.set_xticklabels([])
    if not ylab:
        ax.set_yticklabels([])


def main():
    pf.use_style("slide")
    OUT.mkdir(parents=True, exist_ok=True)
    lay = pf.Layout.slide([["A", "B", "C"], ["A", "D", "E"]], kind="wide", col_ratios=[1.5, 1, 1])
    panel_a(lay.axes["A"], covariation())
    p = pattern()
    for k, (sp, org, title, xl, yl) in zip("BCDE", (
            ("human", "cortex", "adrenal cortex", "", "correlation\nbetween fetuses"),
            ("human", "lung", "lung", "", ""),
            ("human", "liver", "liver", "expected from\nthe adult clock", "correlation\nbetween fetuses"),
            ("mouse", "liver", "liver (mouse)", "expected from\nthe adult clock", ""))):
        panel_pairs(lay.axes[k], p, sp, org, title, xl, yl)
    h = [Line2D([], [], ls="", marker="o", ms=8, color=SPECIES["human"], label="human fetus"),
         Line2D([], [], ls="", marker="o", ms=8, color=SPECIES["mouse"], label="mouse embryo"),
         Line2D([], [], ls="", marker="o", ms=8, mfc="white", mec="#222222", label="not clock-specific"),
         Line2D([], [], ls="", marker="D", ms=7, mfc="white", mec="#222222", label="expected at adult size"),
         Line2D([], [], ls="", marker="|", ms=12, mew=2.2, color=MATCHED, label="matched genes, 95th pct")]
    lay.shared_legend("ABCDE", h, [x.get_label() for x in h], side="bottom", ncol=5, fontsize=11,
                      handletextpad=0.3, columnspacing=1.2)
    fig = lay.finalize()
    print(pf.collisions(fig))
    pf.save(fig, OUT / "development")


if __name__ == "__main__":
    main()
