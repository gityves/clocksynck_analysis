"""Slide figure (results, development): before birth, some organs carry a clock-shaped state shared between individuals.

A. Clock-gene covariation between individuals (C, signed coherence of the 14 clock genes against the adult waveform,
   cohort / batch and age adjusted) per organ, human fetus (adrenal cortex A244; lung, liver, heart, kidney A243, pooled
   open cohorts, 4-20 wpc) and late mouse embryo (A309, Qiu 2024, E15-E18.75: adrenocortical cells, hepatocytes, lung
   progenitors, ventricular cardiomyocytes, nephron progenitors). Filled: significant (Holm) and above the 95th percentile
   of 200 expression-matched and 200 stress-responsive gene sets; open: not. Open diamond: median C of adult tissue
   subsampled to the same design (full-size clock). Grey tick: 95th percentile of the expression-matched sets.
B. Link of the signal-induced clock genes (PER1, PER2, NR1D1) to glucocorticoid-target activity (L, mean r^2 with the
   FKBP5 / TSC22D3 / ZBTB16 / KLF9 score; A307 human, A309 mouse liver).
C. Agreement of the signal-induced genes with the core genes in the clock's phase pattern (W, mean r cos(adult phase
   difference) over the 33 pairs; A307, A309). B and C: filled = above the matched 95th percentile and significant.

Status: A243, A244 registered; A307, A309 registered, NOT BLIND for the earlier fetal results. Aggregates only.
Data: exploratory/A243_fetal_organs_pooled/output/{l1,verdicts,position_in_adult}.csv;
exploratory/A244_fetal_adrenal_classes/output/{l1,verdicts_a244,position}.csv;
exploratory/A309_mouse_embryo_organs/output/{l1,power_passrates,gc_pattern}.csv; exploratory/A307_fetal_gc_pattern/output/tests.csv.
Output: presentation/figures/development/development.{svg,png}.

    python scripts/development.py
"""
from pathlib import Path

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
    """L and W of the signal-induced genes: human cortex, lung, liver (A307), mouse liver (A309)."""
    h = pd.read_csv(EXP / "A307_fetal_gc_pattern/output/tests.csv").assign(species="human")
    m = pd.read_csv(EXP / "A309_mouse_embryo_organs/output/gc_pattern.csv").rename(columns={"readout": "organ"})
    d = pd.concat([h, m.assign(species="mouse")])
    d["organ"] = d.organ.replace({"cortex": "adrenal cortex"})
    d["ok"] = d["pass"].astype(bool)
    return d[["organ", "species", "statistic", "value", "matched_q95", "ok"]]


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


def panel_bc(ax, d, stat, xlabel, xlim):
    for _, r in d[d.statistic == stat].iterrows():
        mark(ax, r.value, ORGANS.index(r.organ) + DY[r.species], r.species, r.ok, r.matched_q95)
    rows(ax, len(ORGANS))
    ax.set_xlim(*xlim)
    ax.set_xlabel(xlabel)


def main():
    pf.use_style("slide")
    OUT.mkdir(parents=True, exist_ok=True)
    lay = pf.Layout.slide("ABC", kind="wide", col_ratios=[1.4, 1, 1])
    panel_a(lay.axes["A"], covariation())
    p = pattern()
    panel_bc(lay.axes["B"], p, "L", "link to glucocorticoid\ntargets, L", (-0.02, 0.55))
    panel_bc(lay.axes["C"], p, "W", "agreement with the clock's\nphase pattern, W", (-0.02, 0.13))
    for k in "BC":
        lay.axes[k].set_yticklabels([])
    lay.axes["C"].set_xticks([0, 0.05, 0.10])
    h = [Line2D([], [], ls="", marker="o", ms=8, color=SPECIES["human"], label="human fetus"),
         Line2D([], [], ls="", marker="o", ms=8, color=SPECIES["mouse"], label="mouse embryo"),
         Line2D([], [], ls="", marker="o", ms=8, mfc="white", mec="#222222", label="not clock-specific"),
         Line2D([], [], ls="", marker="D", ms=7, mfc="white", mec="#222222", label="expected at adult size"),
         Line2D([], [], ls="", marker="|", ms=12, mew=2.2, color=MATCHED, label="matched genes, 95th pct")]
    lay.shared_legend("ABC", h, [x.get_label() for x in h], side="bottom", ncol=5, fontsize=11,
                      handletextpad=0.3, columnspacing=1.2)
    fig = lay.finalize()
    print(pf.collisions(fig))
    pf.save(fig, OUT / "development")


if __name__ == "__main__":
    main()
