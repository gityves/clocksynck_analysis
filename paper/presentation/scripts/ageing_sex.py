"""Slide figure (results, R2.4): ageing shifts the body's clock as a whole; neither age nor sex changes its organisation.

A. CHIRAL's output-rhythm statistic, re-read (WP0 item 5, confirmation of archive A101b / A102b): total log ratio of the
   number of rhythmic genes, old (60-79 y) over young (< 50 y) and female over male, with CHIRAL's normalisation (TMM),
   for CHIRAL's donor subsample (re-run here), all donors at equal n, donors matched on death type (Hardy class) and on death
   type x RIN tertile; tissue bootstrap 95 % interval.
B, C. Core-clock organisation (clocksync, GTEx fast and ventilated donors matched within cohort): pooled difference of
   coupling at a common host spread (d*) and of organ offsets, old - young (A231) and female - male (A208), with 90 %
   intervals; dashed: the registered equivalence bounds (+-0.15, +-1 h). Both registered (H231.1, H208.1 pass; offsets
   equivalent within +-1 h).
D. Internal time (untimed host phase) minus recorded time of death, change per decade of age, body and brain hosts
   (A234, registered, H234.1 and H234.3 pass), 90 % interval; fast deaths, strata Hardy x sex. A timed readout.

Data: exploratory/WP0_integrity/output/chiral_tmm/{count_comparison_tmm_summary,sex_comparison_matched_summary}.csv;
exploratory/A231_age_coupling/output/score_tests.csv; exploratory/A208_sex_coupling/output/tests.csv;
exploratory/A234_internal_time_age/output/score_tests.csv. GTEx aggregates only.
Output: presentation/figures/ageing_sex/ageing_sex.{svg,png}.

    python scripts/ageing_sex.py
"""
from pathlib import Path

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


def main():
    pf.use_style("slide")
    OUT.mkdir(parents=True, exist_ok=True)
    lay = pf.Layout.slide("ABCD", kind="wide", col_ratios=[1.3, 1, 1, 1])
    chiral(lay.axes["A"])
    organisation(lay.axes["B"], "d_star", 0.15, "coupling difference, d*", "coupling")
    organisation(lay.axes["C"], "offset_h", 1.0, "organ offset difference [h]", "organ offsets")
    timing(lay.axes["D"])
    fig = lay.finalize()
    print(pf.collisions(fig))
    pf.save(fig, OUT / "ageing_sex")


if __name__ == "__main__":
    main()
