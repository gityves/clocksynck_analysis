"""Slide figure (results, adult health, second slide after the reference map): whole blood is the exception; its core
clock genes keep no better time than random genes, while the other quiet tissues (testis, brain) keep it.

GTEx v10 fast deaths with a recorded time of death. Whole blood (212 donors), the other low-gain tissues (testis, 135;
13 brain regions, 98-169) and five solid tissues (179-220). Brain: one bar at the median over the 13 regions, one dot per
region.
A. Rhythmicity of the 14 core clock genes against the hour of death: median cosinor R^2 of the clock genes, as a
   percentile of all expressed genes of the tissue. Dashed: 50th percentile (a typical gene).
B. Phase agreement of the clock genes with the GTEx waveform, as measured (filled bars) and after regressing out the
   blood-cell composition scores (open bars; neutrophil, T, B, NK, monocyte, eosinophil, erythroid markers).

Status: post hoc, descriptive (user observation stated before running, 2026-11-03; `results_preliminary.md`, "GTEx
Whole Blood has no core-clock time signal"); testis and brain added post hoc (2026-11-25) as low-gain comparisons.
GTEx aggregates only (one row per tissue).
Data: exploratory/blood_composition_posthoc/output/{blood_composition,blood_composition_low_gain}.csv.
Output: presentation/figures/blood_exception/blood_exception.{svg,png} (wide).

    python scripts/blood_exception.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.patches import Patch

import paperfig as pf
from _systems import COLOR, gtex_system

PAPER = Path(__file__).resolve().parents[2]
SRC = PAPER / "exploratory/blood_composition_posthoc/output"
OUT = PAPER / "presentation" / "figures" / "blood_exception"
INK, MUTED = "#222222", pf.GREY
SHORT = {"Whole Blood": "whole blood", "Testis": "testis", "Brain": "brain", "Muscle - Skeletal": "muscle",
         "Lung": "lung", "Thyroid": "thyroid", "Skin - Sun Exposed (Lower leg)": "skin", "Adipose - Subcutaneous": "fat"}
VALUES = ("clock_R2_percentile_median", "clock_phase_concordance", "clock_phase_concordance_after_composition")


def rows():
    """One row per bar (brain: median over regions) and the brain regions themselves."""
    d = pd.read_csv(SRC / "blood_composition.csv")
    low = pd.read_csv(SRC / "blood_composition_low_gain.csv")
    brain = low[low.tissue.str.startswith("Brain")]
    solid = d[d.tissue != "Whole Blood"].sort_values("clock_phase_concordance", ascending=False)
    bars = pd.concat([d[d.tissue == "Whole Blood"], low[low.tissue == "Testis"],
                      pd.DataFrame([dict(tissue="Brain", **brain[list(VALUES)].median())]), solid])
    return bars.reset_index(drop=True), brain


def colour(t):
    return COLOR[gtex_system(t)]


def tissue_axis(ax, d):
    ax.set_xticks(range(len(d)), [SHORT[t] for t in d.tissue], rotation=45, ha="right", rotation_mode="anchor")
    ax.tick_params(axis="x", length=0, pad=6)
    ax.set_xlim(-0.6, len(d) - 0.4)


def regions(ax, x, y):
    jitter = np.random.default_rng(0).uniform(-0.1, 0.1, len(y))
    ax.scatter(x + jitter, y, s=10, facecolor="white", edgecolor=INK, lw=0.8, zorder=4, clip_on=False)


def percentile(ax, d, brain):
    ax.bar(range(len(d)), d.clock_R2_percentile_median, width=0.7, color=[colour(t) for t in d.tissue], zorder=2)
    regions(ax, d.index[d.tissue == "Brain"][0], brain.clock_R2_percentile_median)
    ax.axhline(50, color=MUTED, lw=0.8, ls="--", zorder=3)
    ax.annotate("typical\ngene", (len(d) - 0.4, 50), xytext=(4, 0), textcoords="offset points", ha="left",
                va="center", fontsize=9, color=MUTED, annotation_clip=False)
    ax.set_ylim(0, 100)
    ax.set_yticks([0, 50, 100])
    ax.set_ylabel("clock-gene rhythmicity\n[percentile of all genes]")
    tissue_axis(ax, d)


def concordance(ax, d, brain):
    w = 0.36
    for i, r in d.iterrows():
        c = colour(r.tissue)
        ax.bar(i - w / 2, r.clock_phase_concordance, width=w, color=c, zorder=2)
        ax.bar(i + w / 2, r.clock_phase_concordance_after_composition, width=w, facecolor="white", edgecolor=c,
               lw=1.2, zorder=2)
    i = d.index[d.tissue == "Brain"][0]
    regions(ax, i - w / 2, brain.clock_phase_concordance)
    regions(ax, i + w / 2, brain.clock_phase_concordance_after_composition)
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
    d, brain = rows()
    lay = pf.Layout.slide("AxB", kind="wide", col_ratios=[1, 0.12, 1], blank=("x",))
    percentile(lay.axes["A"], d, brain)
    concordance(lay.axes["B"], d, brain)
    fig = lay.finalize()
    pd.set_option("display.width", 200)
    print(d[["tissue", *VALUES]].round(4).to_string())
    print(brain[list(VALUES)].agg(["min", "max"]).round(4).to_string())
    print(pf.collisions(fig))
    pf.save(fig, OUT / "blood_exception")


if __name__ == "__main__":
    main()
