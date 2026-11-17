"""Slide figure (results, R1.2): error of the GTEx organ map against the number of donors, full day vs hospital hours.

GTEx benchmark of Figure 1 (`results/figure_1/scripts/method/compute_gtex_benchmark.py --tissue-set working_hours`, post hoc
figure design, run 2026-11-16 on current code): the 40 tissues with >= 20 fast deaths (Hardy 1-2) between 08 and 18 h
(Muscle - Skeletal the reference). A draw is n donors from all fast deaths (`full`, full day) or from those who died 08-18 h
(`window`); 10 draws per n (CYCLOPS: 3, for cost). Each method gives per-tissue offsets and relative amplitudes, scored against the timed truth:
per-tissue cosinor of the clock genes on the recorded time of death of all fast deaths, rank-1 complex SVD.
An untimed method fixes neither the direction of its cycle nor its origin, so every method (clocksync and the timed
oracle included) is scored after the best alignment of its offsets to the truth, per draw: reflection and rotation
(0.1 h grid) minimising the median |offset error| over tissues (user, 2026-11-17: the fair comparison).
Rows, per draw: median over tissues of |offset error|; share of tissues whose offset is off by > 2 h (a median over
tissues hides them); median over tissues of |log2 amplitude error|. Plotted: median over draws (line, dot) and
interquartile range over draws (band).
Methods: clocksync (untimed); CHIRAL, CoFE on the clock genes, CYCLOPS 2.0 on the clock genes as both expression and
seeds (untimed orderers); timed oracle (cosinor on the drawn donors' recorded times).

Output: paper/presentation/figures/gtex_sweep/gtex_sweep.{svg,png} (wide slide figure).

    python scripts/gtex_sweep.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D

import paperfig as pf

PAPER = Path(__file__).resolve().parents[2]
SRC = PAPER / "results/figure_1/output/method/gtex_benchmark_40/tissues.csv"
OUT = Path(__file__).resolve().parents[1] / "figures" / "gtex_sweep"
INK = "#222222"
TAB = sns.color_palette().as_hex()
METHODS = {  # key: (label, colour, line style)
    "clocksync": ("clocksync", TAB[0], "-"),
    "chiral": ("CHIRAL", TAB[1], "-"),
    "cofe_clock": ("CoFE, clock genes", TAB[2], "-"),
    "cyclops_clock": ("CYCLOPS 2.0, clock genes", TAB[4], "-"),
    "timed": ("timed oracle", INK, ":"),
}
DESIGNS = {"full": "full day", "window": "hospital hours"}
MISS_H = 2.0
METRICS = {"offset": "offset\nerror [h]", "miss": f"tissues off\nby > {MISS_H:g} h [%]", "amp": "amplitude\nerror [|log2|]"}


def wrap(h):
    return (np.asarray(h, float) + 12) % 24 - 12


def aligned_errors(g: pd.DataFrame) -> pd.Series:
    """Offset errors after the reflection and rotation of the estimates that minimise the median |error| over tissues."""
    est, true = (g.offset_h_true + g.offset_err_h).to_numpy(), g.offset_h_true.to_numpy()
    rot = np.arange(0, 24, 0.1)[:, None]
    best = min((wrap(d * est + rot - true) for d in (1, -1)), key=lambda e: np.nanmedian(np.abs(e), axis=1).min())
    return pd.Series(best[np.nanmedian(np.abs(best), axis=1).argmin()], index=g.index)


def per_draw(t: pd.DataFrame) -> pd.DataFrame:
    t = t.dropna(subset=["offset_h_true"])
    err = t.groupby(["method", "design", "n", "draw"], group_keys=False).apply(aligned_errors).abs()
    t = t.assign(offset=err, miss=100 * (err > MISS_H), amp=t.log2_amp_err.abs())
    return t.groupby(["method", "design", "n", "draw"]).agg(offset=("offset", "median"), miss=("miss", "mean"),
                                                            amp=("amp", "median")).reset_index()


def main():
    pf.use_style("slide")
    OUT.mkdir(parents=True, exist_ok=True)
    d = per_draw(pd.read_csv(SRC))
    present = [m for m in METHODS if m in set(d.method)]
    g = d.groupby(["method", "design", "n"])[list(METRICS)]
    med, lo, hi = g.median(), g.quantile(0.25), g.quantile(0.75)
    lay = pf.Layout.slide("AB;CD;EF;LL", kind="wide", row_heights=[1, 1, 1, 0.12], align="all")
    grid = dict(zip(((r, c) for r in METRICS for c in DESIGNS), "ABCDEF"))
    for (metric, design), key in grid.items():
        ax = lay.axes[key]
        for m in present:
            lab, c, ls = METHODS[m]
            y, a, b = (x.loc[(m, design), metric] for x in (med, lo, hi))
            n = y.index.to_numpy()
            ax.fill_between(n, a, b, color=c, alpha=0.1, lw=0)
            top = m == "clocksync"
            ax.plot(n, y, color=c, ls=ls, lw=2.4 if top else 1.4, marker="o", ms=4, zorder=3 if top else 2)
        ax.set_xscale("log", base=2)
        ax.set_xticks([10, 20, 40, 80], ["10", "20", "40", "80"])
        ax.minorticks_off()
        ax.set_ylim(0, None)
        if design == "full":
            ax.set_ylabel(METRICS[metric])
        if metric == "offset":
            ax.set_title(DESIGNS[design], fontsize=13, color=INK)
        if metric == "amp":
            ax.set_xlabel("donors")
        else:
            ax.tick_params(labelbottom=False)
    for left, right in ("AB", "CD", "EF"):
        lay.axes[right].sharey(lay.axes[left])
    h = [Line2D([], [], color=METHODS[m][1], ls=METHODS[m][2], lw=2, label=METHODS[m][0]) for m in present]
    L = lay.axes["L"]
    L.legend(handles=h, loc="center", ncol=len(h), frameon=False, fontsize=10, handlelength=1.6, columnspacing=1.5)
    L.axis("off")
    fig = lay.finalize()
    pd.set_option("display.width", 200)
    print(med.round(2).unstack("design").to_string())
    print(pf.collisions(fig))
    pf.save(fig, OUT / "gtex_sweep")


if __name__ == "__main__":
    main()
