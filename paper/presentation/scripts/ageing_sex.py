"""Slide figure (results, 3.6; formerly R2.4): age and sex change neither the organisation nor the loudness of the core
clock.

A. CHIRAL's output-rhythm statistic, re-read (WP0 item 5, confirmation of archive A101b / A102b): total log ratio of the
   number of rhythmic genes, old (60-79 y) over young (< 50 y) and female over male, with CHIRAL's normalisation (TMM),
   for CHIRAL's donor subsample (re-run here), all donors at equal n, donors matched on death type (Hardy class) and on death
   type x RIN tertile; tissue bootstrap 95 % interval.
B, C. Core-clock organisation (clocksync, GTEx fast and ventilated donors matched within cohort): pooled difference of
   coupling at a common host spread (d*) and of organ offsets, old - young (A231) and female - male (A208), with 90 %
   intervals; dashed: the registered equivalence bounds (+-0.15, +-1 h). Both registered (H231.1, H208.1 pass; offsets
   equivalent within +-1 h).
D. Loudness, old (60-70 y) over young (21-49 y), donors matched on death type (A231 pairs): amplitude ratio of the
   core clock (14 genes; A373 H373.2, registered, equivalence pass) and of each tissue's rhythmic output programme
   (A373 H373.3), 90 % intervals; dashed: the registered equivalence margin (+-0.21 on the log variance, amplitude
   0.90-1.11). Replaces the timing panel (user, 2026-11-25: body-wide shifts are out of scope); `timing()` (A234,
   internal time advance per decade) is kept for a backup.

Data: exploratory/WP0_integrity/output/chiral_tmm/{count_comparison_tmm_summary,sex_comparison_matched_summary}.csv;
exploratory/A231_age_coupling/output/score_tests.csv; exploratory/A208_sex_coupling/output/tests.csv;
exploratory/A373_age_loudness/output/{tests,out_tests}.csv (A234_internal_time_age/output/score_tests.csv for
`timing`). GTEx aggregates only.
Output: presentation/figures/ageing_sex/ageing_sex.{svg,png}.

    python scripts/ageing_sex.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D

import paperfig as pf

PAPER = Path(__file__).resolve().parents[2]
EXP = PAPER / "exploratory"
OUT = PAPER / "presentation" / "figures" / "ageing_sex"
INK, MUTED = "#222222", pf.GREY
TAB = sns.color_palette().as_hex()
CS, CHIRAL = TAB[0], TAB[1]
VERSIONS = {"R1": "CHIRAL's\nsubsample", "pooled": "all donors", "matched": "matched on\ndeath type",
            "matched_rin": "+ RNA quality"}
COMPARISON = {"age_old_vs_young": ("old / young", "o"), "sex_female_vs_male": ("female / male", "D")}


def chiral(ax):
    W = EXP / "WP0_integrity/output/chiral_tmm"
    d = pd.concat([pd.read_csv(W / "count_comparison_tmm_summary.csv"), pd.read_csv(W / "sex_comparison_matched_summary.csv")])
    d = d[(d.normalisation == "TMM") & (d.statistic == "total_logratio")]
    for j, (cmp, (label, m)) in enumerate(COMPARISON.items()):
        for i, v in enumerate(VERSIONS):
            r = d[(d.comparison == cmp) & (d.version == v)].iloc[0]
            y = i + (j - 0.5) * 0.3
            ax.plot([r.lo, r.hi], [y, y], color=CHIRAL, lw=1.4, zorder=1)
            ax.scatter([r.estimate], [y], s=32, marker=m, color=CHIRAL if j == 0 else "white", edgecolor=CHIRAL,
                       lw=1.2, zorder=2)
    ax.axvline(0, color=MUTED, lw=0.8, ls="--", zorder=0)
    ax.set_yticks(range(len(VERSIONS)), list(VERSIONS.values()))
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(len(VERSIONS) - 0.5, -0.5)
    ax.set_xlabel("rhythmic output genes,\nlog ratio")
    ax.set_title("CHIRAL", fontsize=13, color=INK)
    h = [Line2D([], [], ls="", marker=m, ms=6, mfc=CHIRAL if j == 0 else "white", mec=CHIRAL, label=lab)
         for j, (lab, m) in enumerate(COMPARISON.values())]
    ax.legend(handles=h, loc="lower left", frameon=False, fontsize=8, handletextpad=0.2)


def organisation(ax, readout, bound, xlabel, title):
    rows = [("old - young", pd.read_csv(EXP / "A231_age_coupling/output/score_tests.csv").set_index("readout")),
            ("female - male", pd.read_csv(EXP / "A208_sex_coupling/output/tests.csv").set_index("statistic"))]
    for i, (_, t) in enumerate(rows):
        r = t.loc[readout]
        ax.plot([r.lo90, r.hi90], [i, i], color=CS, lw=1.6, zorder=1)
        ax.scatter([r.estimate], [i], s=34, color=CS, zorder=2)
    for b in (-bound, bound):
        ax.axvline(b, color=MUTED, lw=0.8, ls="--", zorder=0)
    ax.axvline(0, color=MUTED, lw=0.6, zorder=0)
    ax.set_xlim(-1.4 * bound, 1.4 * bound)
    ax.set_yticks(range(len(rows)), [r[0] for r in rows])
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(len(rows) - 0.4, -0.6)
    ax.set_xlabel(xlabel)
    ax.set_title(title, fontsize=13, color=INK)


def timing(ax):
    t = pd.read_csv(EXP / "A234_internal_time_age/output/score_tests.csv").set_index("test")
    for i, (k, lab) in enumerate((("H234.1", "body"), ("H234.3", "brain"))):
        r = t.loc[k]
        ax.plot([r.lo90, r.hi90], [i, i], color=CS, lw=1.6, zorder=1)
        ax.scatter([r.beta_h_per_decade], [i], s=34, color=CS, zorder=2)
    ax.axvline(0, color=MUTED, lw=0.8, ls="--", zorder=0)
    ax.set_xlim(-0.5, 2)
    ax.set_yticks([0, 1], ["body", "brain"])
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(1.6, -0.6)
    ax.set_xlabel("internal time advance\n[h per decade]")
    ax.set_title("timing", fontsize=13, color=INK)


def loudness(ax):
    c = pd.read_csv(EXP / "A373_age_loudness/output/tests.csv").set_index("gene_set").loc["14 genes"]
    o = pd.read_csv(EXP / "A373_age_loudness/output/out_tests.csv").set_index("test").loc["H373.3"]
    rows = [("core clock", c, CS), ("rhythmic\noutputs", o, INK)]
    for i, (_, r, col) in enumerate(rows):
        ax.plot([r.amp_ratio_lo90, r.amp_ratio_hi90], [i, i], color=col, lw=1.6, zorder=1)
        ax.scatter([r.amp_ratio], [i], s=34, color=col, zorder=2)
    for b in (np.exp(-0.105), np.exp(0.105)):
        ax.axvline(b, color=MUTED, lw=0.8, ls="--", zorder=0)
    ax.axvline(1, color=MUTED, lw=0.6, zorder=0)
    ax.set_xlim(0.84, 1.18)
    ax.set_yticks(range(len(rows)), [r[0] for r in rows])
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(len(rows) - 0.4, -0.6)
    ax.set_xlabel("amplitude, old / young")
    ax.set_title("loudness", fontsize=13, color=INK)
    return pd.DataFrame([dict(row=lab, ratio=r.amp_ratio, lo=r.amp_ratio_lo90, hi=r.amp_ratio_hi90) for lab, r, _ in rows])


def main():
    pf.use_style("slide")
    OUT.mkdir(parents=True, exist_ok=True)
    lay = pf.Layout.slide("ABCD", kind="wide", col_ratios=[1.3, 1, 1, 1])
    chiral(lay.axes["A"])
    organisation(lay.axes["B"], "d_star", 0.15, "coupling difference, d*", "coupling")
    organisation(lay.axes["C"], "offset_h", 1.0, "organ offset difference [h]", "organ offsets")
    lo = loudness(lay.axes["D"])
    fig = lay.finalize()
    print(lo.round(3).to_string(index=False))
    print(pf.collisions(fig))
    pf.save(fig, OUT / "ageing_sex")


if __name__ == "__main__":
    main()
