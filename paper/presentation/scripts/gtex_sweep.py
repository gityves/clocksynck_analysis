"""Slide figure (results, R1.2): error of the GTEx organ map against the number of donors, whole day vs hospital hours.

GTEx benchmark of Figure 1 (`results/figure_1/scripts/method/compute_gtex_benchmark.py`, post hoc figure design; clocksync
refitted 2026-11-14 on current code, other methods run 2026-10-25): a draw is n donors from all fast deaths (Hardy 1-2; `full`, whole day) or from those who died 08-18 h (`window`);
8 tissues (Muscle - Skeletal the reference). Each method gives per-tissue offsets and relative amplitudes, scored against the
timed truth: per-tissue cosinor of the clock genes on the recorded time of death of all fast deaths, rank-1 complex SVD.
Per draw: median over tissues of |offset error| (h) and |log2 amplitude error|; plotted: median over draws (line, dot) and
interquartile range over draws (band). 10 draws per n; CYCLOPS 2.0 and CoFE genome-wide (cost) only at n = 80, 3 draws.
Methods: clocksync (untimed); CHIRAL, CoFE (clock genes, genome-wide), CYCLOPS 2.0 (untimed orderers); TimeTeller (trained on
recorded times, supervised reference); timed oracle (cosinor on the drawn donors' recorded times).

Output: paper/presentation/figures/gtex_sweep/gtex_sweep.{svg,png} (wide slide figure).

    python scripts/gtex_sweep.py
"""
from pathlib import Path

import numpy as np
import pandas as pd

import paperfig as pf

PAPER = Path(__file__).resolve().parents[2]
SRC = PAPER / "results/figure_1/output/method/gtex_benchmark/tissues.csv"
OUT = Path(__file__).resolve().parents[1] / "figures" / "gtex_sweep"
INK, MUTED = "#222222", pf.GREY
METHODS = {  # key: (label, colour, line style)
    "clocksync": ("clocksync", pf.OKABE_ITO[3], "-"),
    "chiral": ("CHIRAL", pf.OKABE_ITO[0], "-"),
    "cofe_clock": ("CoFE, clock genes", pf.OKABE_ITO[2], "-"),
    "cofe_genome": ("CoFE, genome-wide", pf.OKABE_ITO[4], "-"),
    "cyclops": ("CYCLOPS 2.0", pf.OKABE_ITO[1], "-"),
    "timeteller": ("TimeTeller (supervised)", MUTED, "--"),
    "timed": ("timed oracle", INK, ":"),
}
DESIGNS = {"full": "whole day", "window": "hospital hours (08-18 h)"}
METRICS = {"offset_err_h": "offset error (h)", "log2_amp_err": "amplitude error (|log2|)"}


def per_draw(t: pd.DataFrame) -> pd.DataFrame:
    return (t.assign(**{m: t[m].abs() for m in METRICS}).groupby(["method", "design", "n", "draw"])[list(METRICS)].median()
            .reset_index())


def main():
    pf.use_style("slide")
    OUT.mkdir(parents=True, exist_ok=True)
    d = per_draw(pd.read_csv(SRC))
    g = d.groupby(["method", "design", "n"])[list(METRICS)]
    med, lo, hi = g.median(), g.quantile(0.25), g.quantile(0.75)
    lay = pf.Layout.slide("AB;CD;LL", kind="wide", row_heights=[1, 1, 0.1])
    axes = {("offset_err_h", "full"): "A", ("offset_err_h", "window"): "B",
            ("log2_amp_err", "full"): "C", ("log2_amp_err", "window"): "D"}
    for (metric, design), key in axes.items():
        ax = lay.axes[key]
        for m, (lab, c, ls) in METHODS.items():
            if (m, design) not in med.index.droplevel("n"):
                continue
            y, a, b = (x.loc[(m, design), metric] for x in (med, lo, hi))
            n = y.index.to_numpy()
            if len(n) > 1:
                ax.fill_between(n, a, b, color=c, alpha=0.07, lw=0)
            ax.plot(n, y, color=c, ls=ls, lw=2.2 if m == "clocksync" else 1.4, marker="o", ms=4, zorder=3 if m == "clocksync" else 2)
            if len(n) == 1:
                ax.errorbar(n, y, yerr=[y - a, b - y], color=c, lw=1.2, capsize=0)
        ax.set_xscale("log", base=2)
        ax.set_xticks([10, 20, 40, 80], ["10", "20", "40", "80"])
        ax.set_ylim(0, None)
        ax.set_ylabel(METRICS[metric])
        if key in "AB":
            ax.set_title(DESIGNS[design], fontsize=13, color=INK)
        else:
            ax.set_xlabel("donors")
    for k in "CD":
        lay.axes[k].sharey(lay.axes["C"])
    lay.axes["B"].sharey(lay.axes["A"])
    L = lay.axes["L"]
    for j, (m, (lab, c, ls)) in enumerate(METHODS.items()):
        x, y = (0.0, 0.11, 0.22, 0.39, 0.56, 0.68, 0.87)[j], 0.5
        L.plot([x, x + 0.02], [y, y], color=c, ls=ls, lw=2, transform=L.transAxes, clip_on=False)
        L.annotate(lab, (x + 0.02, y), xycoords="axes fraction", xytext=(6, 0), textcoords="offset points", va="center",
                   fontsize=10, color=INK)
    L.axis("off")
    fig = lay.finalize()
    pd.set_option("display.width", 200)
    print(med.round(3).unstack("design").to_string())
    print(pf.collisions(fig))
    pf.save(fig, OUT / "gtex_sweep")


if __name__ == "__main__":
    main()
