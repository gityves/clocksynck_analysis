"""Slide figure (results, before birth, 2.3): the fetal liver and lung covariation is not current glucocorticoid activity
under another name, and glucocorticoid activity does not rank the fetal organs.

A. 14 clock genes, clock-gene covariation between individuals C, before (filled, A243) and after the glucocorticoid score
   is regressed out (hatched; A269 H269.2, pass in both). The score: per fetus, the mean of the standardised log
   fractions of four direct glucocorticoid-receptor targets that are not clock genes (FKBP5, TSC22D3, ZBTB16, KLF9);
   each clock gene's log expression is regressed on cohort, quadratic age and the score (own slope per gene) and C is
   computed on the residuals; the companion gene sets get the same adjustment. Grey line: matched genes, 95th
   percentile; both readouts are also above the stress-responsive 95th percentile (clock-specific). This removes a
   single shared glucocorticoid axis, linear in the activity at sampling; it does not remove a rhythmic input acting
   with delays.
B. Seven fetal organs: C of the 14 genes (A243, A244, A239) against the between-fetus variance of the
   glucocorticoid-target score (counting-noise corrected; brain: three genes, uncorrected, open marker) (A341, descriptive;
   Spearman 0.25, one-sided p 0.30). Filled: clock-specific C. Colours: organ systems (`_systems.py`).

Status: A243 registered, blind; A269 registered, not blind for the earlier fetal results; A341 descriptive by its power
rule. Aggregates only.
Data: exploratory/A269_fetal_organs_drive/output/l1.csv; exploratory/A243_fetal_organs_pooled/output/l1.csv;
exploratory/A341_gr_ranking/output/{measures,spearman}.csv.
Output: presentation/figures/fetal_signal/fetal_signal.{svg,png} (half).

    python scripts/fetal_signal.py
"""
from pathlib import Path

import pandas as pd
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

import paperfig as pf
from _systems import COLOR

PAPER = Path(__file__).resolve().parents[2]
EXP = PAPER / "exploratory"
OUT = PAPER / "presentation" / "figures" / "fetal_signal"
INK, MATCHED = "#222222", "#BBBBBB"
ORGANS = ("liver", "lung")
SYSTEM = {"adrenal cortex": "endocrine", "adrenal medulla": "endocrine", "liver": "digestive", "lung": "other",
          "kidney": "other", "heart": "cardiovascular", "brain": "brain, eye"}
CLOCK_SPECIFIC = {"adrenal cortex", "liver", "lung"}  # 14-gene C clock-specific (A243, A244; slide 2.1)


def adjusted(ax):
    base = pd.read_csv(EXP / "A243_fetal_organs_pooled/output/l1.csv").set_index("test")
    adj = pd.read_csv(EXP / "A269_fetal_organs_drive/output/l1.csv").set_index("test")
    rows = []
    for i, o in enumerate(ORGANS):
        for dx, r, lab in ((-0.19, base.loc[f"{o}:all"], "14 genes"), (0.19, adj.loc[f"{o}:gc"], "adjusted")):
            ok = bool(r.clock_specific)
            hatch = None if lab == "14 genes" else "///"
            fill = (INK if hatch is None else "white") if ok else "white"
            ax.bar([i + dx], [r.C], width=0.33, color=fill, edgecolor=INK, hatch=hatch, lw=1.4, zorder=2)
            ax.plot([i + dx - 0.165, i + dx + 0.165], [r.matched_q95] * 2, color=MATCHED, lw=2.2,
                    solid_capstyle="butt", zorder=3)
            rows.append(dict(organ=o, readout=lab, C=r.C, p=r.p, matched_q95=r.matched_q95, clock_specific=ok))
    ax.set_xticks(range(len(ORGANS)), ORGANS)
    ax.tick_params(axis="x", length=0)
    ax.set_xlim(-0.6, len(ORGANS) - 0.4)
    ax.set_ylim(0, 0.8)
    ax.set_ylabel("14 clock genes, C")
    return pd.DataFrame(rows)


def ranking(ax):
    m = pd.read_csv(EXP / "A341_gr_ranking/output/measures.csv")
    nudge = {"adrenal cortex": (6, 4, "left"), "adrenal medulla": (-6, 6, "right"), "heart": (-6, 6, "right"),
             "brain": (6, 0, "left"), "kidney": (6, -8, "left"), "lung": (6, 2, "left"), "liver": (-6, -6, "right")}
    for r in m.itertuples():
        c = COLOR[SYSTEM[r.organ]]
        exact = pd.isna(r.gc_genes)
        ax.scatter([r.gc_var], [r.C], s=50, color=c if (r.organ in CLOCK_SPECIFIC and exact) else "white",
                   edgecolor=c, lw=1.4, zorder=3)
        dx, dy, ha = nudge[r.organ]
        ax.annotate(r.organ, (r.gc_var, r.C), xytext=(dx, dy), textcoords="offset points", fontsize=10, color=c, ha=ha,
                    va="center")
    ax.axhline(0, color=pf.GREY, lw=0.6, zorder=0)
    ax.set_xlim(-0.01, 0.2)
    ax.set_ylim(-0.3, 0.8)
    ax.set_xlabel("glucocorticoid-target score,\nvariance between fetuses")
    ax.set_ylabel("14 clock genes, C")
    return m


def main():
    pf.use_style("slide")
    OUT.mkdir(parents=True, exist_ok=True)
    lay = pf.Layout.slide("AB", kind="half", height=110, col_ratios=[0.8, 1.3])
    b = adjusted(lay.axes["A"])
    m = ranking(lay.axes["B"])
    h = [Patch(facecolor=INK, edgecolor=INK, label="14 genes"),
         Patch(facecolor="white", edgecolor=INK, hatch="///", label="glucocorticoid score regressed out"),
         Line2D([], [], color=MATCHED, lw=2.2, label="matched genes, 95th pct")]
    lay.shared_legend("AB", h, [x.get_label() for x in h], side="bottom", ncol=3, fontsize=10, handletextpad=0.4,
                      columnspacing=1.0)
    fig = lay.finalize()
    pd.set_option("display.width", 200)
    print(b.round(3).to_string(index=False),
          m[["organ", "n", "gc_var", "C"]].round(3).to_string(index=False),
          pd.read_csv(EXP / "A341_gr_ranking/output/spearman.csv").round(3).to_string(index=False), sep="\n")
    print(pf.collisions(fig))
    pf.save(fig, OUT / "fetal_signal")


if __name__ == "__main__":
    main()
