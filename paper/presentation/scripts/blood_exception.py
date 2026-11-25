"""Slide figure (results, adult health, 3.6): whole blood is the exception; its core clock genes keep no better time
than random genes.

GTEx v10 fast deaths with a recorded time of death; whole blood (212 donors) against five solid tissues (179-220).
A. Rhythmicity of the 14 core clock genes against the hour of death: median cosinor R^2 of the clock genes, as a
   percentile of all expressed genes of the tissue. Dashed: 50th percentile (a typical gene).
B. Phase agreement of the clock genes with the GTEx waveform, as measured (filled) and after regressing out the cell
   composition scores (open bars; neutrophil, T, B, NK, monocyte, eosinophil, erythroid markers).

Status: post hoc, descriptive (user observation stated before running, 2026-11-03; `results_preliminary.md`, "GTEx
Whole Blood has no core-clock time signal"). GTEx aggregates only (one row per tissue).
Data: exploratory/blood_composition_posthoc/output/blood_composition.csv.
Output: presentation/figures/blood_exception/blood_exception.{svg,png} (wide).

    python scripts/blood_exception.py
"""
from pathlib import Path

import pandas as pd
from matplotlib.patches import Patch

import paperfig as pf
from _systems import COLOR, gtex_system

PAPER = Path(__file__).resolve().parents[2]
OUT = PAPER / "presentation" / "figures" / "blood_exception"
INK, MUTED = "#222222", pf.GREY
SHORT = {"Whole Blood": "whole blood", "Muscle - Skeletal": "muscle", "Lung": "lung", "Thyroid": "thyroid",
         "Skin - Sun Exposed (Lower leg)": "skin", "Adipose - Subcutaneous": "fat"}


def rows():
    d = pd.read_csv(PAPER / "exploratory/blood_composition_posthoc/output/blood_composition.csv")
    blood = d[d.tissue == "Whole Blood"]
    solid = d[d.tissue != "Whole Blood"].sort_values("clock_phase_concordance", ascending=False)
    return pd.concat([blood, solid]).reset_index(drop=True)


def tissue_axis(ax, d):
    ax.set_xticks(range(len(d)), [SHORT[t] for t in d.tissue], rotation=45, ha="right", rotation_mode="anchor")
    ax.tick_params(axis="x", length=0, pad=6)
    ax.set_xlim(-0.6, len(d) - 0.4)


def percentile(ax, d):
    ax.bar(range(len(d)), d.clock_R2_percentile_median, width=0.7,
           color=[COLOR[gtex_system(t)] for t in d.tissue], zorder=2)
    ax.axhline(50, color=MUTED, lw=0.8, ls="--", zorder=3)
    ax.annotate("typical\ngene", (len(d) - 0.4, 50), xytext=(4, 0), textcoords="offset points", ha="left",
                va="center", fontsize=9, color=MUTED, annotation_clip=False)
    ax.set_ylim(0, 100)
    ax.set_yticks([0, 50, 100])
    ax.set_ylabel("clock-gene rhythmicity\n[percentile of all genes]")
    tissue_axis(ax, d)


def concordance(ax, d):
    w = 0.36
    for i, r in d.iterrows():
        c = COLOR[gtex_system(r.tissue)]
        ax.bar(i - w / 2, r.clock_phase_concordance, width=w, color=c, zorder=2)
        ax.bar(i + w / 2, r.clock_phase_concordance_after_composition, width=w, facecolor="white", edgecolor=c,
               lw=1.2, zorder=2)
    ax.set_ylim(0, 1)
    ax.set_yticks([0, 0.5, 1])
    ax.set_ylabel("clock-gene phase agreement\nwith the body waveform")
    tissue_axis(ax, d)
    h = [Patch(facecolor=INK, label="as measured"),
         Patch(facecolor="white", edgecolor=INK, lw=1.2, label="composition removed")]
    ax.legend(handles=h, loc="upper center", bbox_to_anchor=(0.5, 1.12), ncol=2, frameon=False, fontsize=9,
              handlelength=1, columnspacing=1)


def main():
    pf.use_style("slide")
    OUT.mkdir(parents=True, exist_ok=True)
    d = rows()
    lay = pf.Layout.slide("AxB", kind="wide", col_ratios=[1, 0.12, 1], blank=("x",))
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
