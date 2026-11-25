"""Slide figure (results, adult health, 3.6): whole blood is the exception; its core clock genes keep no better time
than random genes.

GTEx v10 fast deaths with a recorded time of death; whole blood (212 donors) against five solid tissues (179-220).
A. Rhythmicity of the 14 core clock genes against the hour of death: median cosinor R^2 of the clock genes, as a
   percentile of all expressed genes of the tissue. Dashed: 50th percentile (a typical gene).
B. Phase agreement of the clock genes with the GTEx waveform, as measured (filled) and after regressing out the cell
   composition scores (open; neutrophil, T, B, NK, monocyte, eosinophil, erythroid markers).

Status: post hoc, descriptive (user observation stated before running, 2026-11-03; `results_preliminary.md`, "GTEx
Whole Blood has no core-clock time signal"). GTEx aggregates only (one row per tissue).
Data: exploratory/blood_composition_posthoc/output/blood_composition.csv.
Output: presentation/figures/blood_exception/blood_exception.{svg,png} (half).

    python scripts/blood_exception.py
"""
from pathlib import Path

import pandas as pd
from matplotlib.lines import Line2D

import paperfig as pf
from _systems import COLOR, gtex_system

PAPER = Path(__file__).resolve().parents[2]
OUT = PAPER / "presentation" / "figures" / "blood_exception"
INK, MUTED = "#222222", pf.GREY
SHORT = {"Whole Blood": "whole blood", "Muscle - Skeletal": "skeletal muscle", "Lung": "lung", "Thyroid": "thyroid",
         "Skin - Sun Exposed (Lower leg)": "skin", "Adipose - Subcutaneous": "subcutaneous fat"}


def rows():
    d = pd.read_csv(PAPER / "exploratory/blood_composition_posthoc/output/blood_composition.csv")
    blood = d[d.tissue == "Whole Blood"]
    solid = d[d.tissue != "Whole Blood"].sort_values("clock_phase_concordance", ascending=False)
    return pd.concat([blood, solid]).reset_index(drop=True)


def percentile(ax, d):
    for i, r in d.iterrows():
        ax.scatter([r.clock_R2_percentile_median], [i], s=40, color=COLOR[gtex_system(r.tissue)], zorder=3)
    ax.axvline(50, color=MUTED, lw=0.8, ls="--", zorder=0)
    ax.set_xlim(0, 105)
    ax.set_xticks([0, 50, 100])
    ax.set_xlabel("clock-gene rhythmicity\n[percentile of all genes]")
    ax.set_yticks(range(len(d)), [SHORT[t] for t in d.tissue])
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(len(d) - 0.5, -0.5)


def concordance(ax, d):
    for i, r in d.iterrows():
        c = COLOR[gtex_system(r.tissue)]
        x0, x1 = r.clock_phase_concordance, r.clock_phase_concordance_after_composition
        ax.plot([x0, x1], [i, i], color=c, lw=1, zorder=2)
        ax.scatter([x0], [i], s=40, color=c, zorder=3)
        ax.scatter([x1], [i], s=40, facecolor="white", edgecolor=c, lw=1.2, zorder=4)
    ax.set_xlim(0, 1)
    ax.set_xticks([0, 0.5, 1])
    ax.set_xlabel("clock-gene phase agreement\nwith the body waveform")
    ax.set_yticks(range(len(d)), [])
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(len(d) - 0.5, -0.5)
    h = [Line2D([], [], ls="", marker="o", ms=6, color=INK, label="as measured"),
         Line2D([], [], ls="", marker="o", ms=6, mfc="white", mec=INK, label="composition removed")]
    ax.legend(handles=h, loc="lower left", frameon=False, fontsize=9, handletextpad=0.1, borderaxespad=0.2)


def main():
    pf.use_style("slide")
    OUT.mkdir(parents=True, exist_ok=True)
    d = rows()
    lay = pf.Layout.slide("AB", kind="half")
    percentile(lay.axes["A"], d)
    concordance(lay.axes["B"], d)
    fig = lay.finalize()
    print(d[["tissue", "donors", "clock_R2_median", "clock_R2_percentile_median", "clock_phase_concordance",
             "clock_phase_concordance_after_composition", "top100_R2", "random100_R2",
             "top100_share_removed"]].round(4).to_string())
    print(pf.collisions(fig))
    pf.save(fig, OUT / "blood_exception")


if __name__ == "__main__":
    main()
