"""Slide figure (results, drives, 4.1): intensive care withdraws the drives: tissue clocks get quieter, most where they
were loudest, and keep their organ order; brain death moves the order and keeps the amplitude.

Rows (A-C): GTEx fast deaths (Hardy 1-2; reference), ventilated post-mortem donors with a living brain, brain-dead
organ donors; 12 body tissues.
A. Correlation of each group's organ-pair phase differences with those of fast deaths (Fisher-Lee rho; A293 part 1):
   clock genes (blue bar, filled when above both lines); lines: 95th percentile of 200 expression-matched (light grey)
   and 200 stress-responsive (dark grey) gene sets, same donors and pairs. Open diamond and thin line: median and 5-95 %
   expected from sampling noise alone, 79 fast deaths against the others (A374 post hoc note; both groups are read at
   79 donors, the brain-dead value being the median over equal-n draws).
B. Clock-specific clock amplitude against the comparison group: ventilated against fast deaths, brain-dead against
   ventilated. exp(Q* / 2), Q* = median over tissues of the log ratio of the debiased clock-plane variance minus the
   median of 200 matched sets (A315, A335 part 1d, pooled donors). Lines: matched (light grey) and responsive (dark grey)
   95th percentiles on the same scale (ventilated row). Brain-dead row: 90 % donor-bootstrap interval and the
   registered equivalence margin (+-0.357 in Q*, shaded; A335 part 3).
C. How far the organs can have moved: 90 % compatibility interval of the organ scatter (SD of each organ's shift, h)
   given the group's arrangement correlation (A374, 14 genes; brain-dead upper bound at the grid edge, arrow); dashed:
   the size of the arrangement (median organ-pair offset in fast deaths).
D. Per tissue, loudness lost in intensive care (host-referenced log z2, fast over ventilated, at a common host spread;
   90 % donor-bootstrap interval, 500 draws) against the tissue's clock strength in fast deaths (z2; A349 D349.2-3,
   descriptive). Organ-system colours.

Status: A293 part 1 registered (H293.2 ventilated pass; brain-dead descriptive); A315 descriptive contrasts (registered)
and H315.1 FAIL (powered); A335 registered (core loop and age / sex matching in the speaker notes); A374 registered;
A349 descriptive. All 14 genes. GTEx aggregates only (one value per group or tissue).
Data: exploratory/A293_organisation_levels/output/part1/groups.csv; exploratory/A374_order_interval/output/
{intervals,noise_reference}.csv; exploratory/A335_review_batch/output/{a315_core,equivalence}.csv;
exploratory/A349_icu_limits/output/part2_tissues.csv.
Output: presentation/figures/intensive_care/intensive_care.{svg,png} (wide).

    python scripts/intensive_care.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import patheffects
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy.stats import spearmanr

import paperfig as pf
from _systems import COLOR, gtex_system

PAPER = Path(__file__).resolve().parents[2]
EXP = PAPER / "exploratory"
OUT = PAPER / "presentation" / "figures" / "intensive_care"
CS = sns.color_palette().as_hex()[0]
INK, MUTED = "#222222", pf.GREY
MATCHED, RESPONSIVE = "#BBBBBB", "#6E6E6E"
EDGE = [patheffects.Stroke(linewidth=3.6, foreground="white"), patheffects.Normal()]
ROWS = {"fast": "fast\ndeaths", "ventilated_postmortem": "ventilated,\nliving brain", "organ_donor": "brain\ndeath"}
IVL = {"ventilated_postmortem": "ventilated, living brain", "organ_donor": "brain-dead organ donors"}
SHORT = {"Thyroid": "thyroid", "Whole Blood": "blood", "Esophagus - Mucosa": "oes. mucosa",
         "Heart - Atrial Appendage": "heart", "Muscle - Skeletal": "muscle",
         "Skin - Not Sun Exposed (Suprapubic)": "skin"}
H = 0.55


def rows(ax, labels=True):
    ax.set_yticks(range(len(ROWS)), list(ROWS.values()) if labels else [])
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(len(ROWS) - 0.5, -0.5)


def qline(ax, y, x, col):
    ax.plot([x, x], [y - H / 2, y + H / 2], color=col, lw=2.2, solid_capstyle="butt", zorder=3,
            path_effects=EDGE if col == RESPONSIVE else None)


def arrangement(ax, g, nz):
    g = g.set_index("group")
    e = nz[(nz.gene_set == "14 genes") & (nz.group == "ventilated, living brain")].iloc[0]
    ax.barh(0, 1, height=H, facecolor="white", edgecolor=CS, lw=1.4, ls="--")
    ax.text(0.5, 0, "reference", ha="center", va="center", fontsize=10, color=CS)
    out = []
    for i, k in enumerate(ROWS):
        if k == "fast":
            continue
        r = g.loc[k]
        ok = r.rho_clock > max(r.matched_q95, r.responsive_q95)
        ax.barh(i, r.rho_clock, height=H, color=CS if ok else "white", edgecolor=CS, lw=1.4, zorder=2)
        qline(ax, i, r.matched_q95, MATCHED)
        qline(ax, i, r.responsive_q95, RESPONSIVE)
        ax.plot([e.rho_q05, e.rho_q95], [i + H / 2 + 0.12] * 2, color=CS, lw=1.2, zorder=4)
        ax.plot(e.rho_median, i + H / 2 + 0.12, ls="", marker="D", ms=6, mfc="white", mec=CS, mew=1.4, zorder=5)
        out.append(dict(group=k, n=r.n, rho=r.rho_clock, matched_q95=r.matched_q95, responsive_q95=r.responsive_q95,
                        noise_median=e.rho_median, noise_q05=e.rho_q05, noise_q95=e.rho_q95, test=r.test))
    ax.set_xlim(0, 1.05)
    rows(ax)
    ax.set_xlabel("organ arrangement\nshared with fast deaths")
    return pd.DataFrame(out)


def amplitude(ax, q, eq):
    """Clock-specific amplitude ratio exp(Q*/2) against the comparison group (variance ratio -> amplitude)."""
    q = q[q.gene_set == "14 genes"].set_index("test")
    amp = lambda x: np.exp(np.asarray(x) / 2)
    ax.axvline(1, color=MUTED, lw=0.8, ls="--", zorder=0)
    ax.text(1, 0, "reference", ha="center", va="center", fontsize=10, color=CS,
            bbox=dict(facecolor="white", edgecolor="none", pad=1), zorder=5)
    v = q.loc["ICU fast vs ventilated PM"]  # Q = log(fast / ventilated): ventilated / fast is -Q*
    x = amp(-v.Q_star)
    ax.barh(1, x - 1, left=1, height=H, color=CS, zorder=2)
    qline(ax, 1, amp(-(v.matched_q95 - v.matched_median)), MATCHED)
    qline(ax, 1, amp(-(v.resp_q95 - v.matched_median)), RESPONSIVE)
    b = eq[(eq.gene_set == "14 genes") & (eq.donors == "pooled")].iloc[0]  # Q = log(brain dead / ventilated)
    m = float(b.margin)
    ax.fill_betweenx([2 - H / 2, 2 + H / 2], amp(-m), amp(m), color=CS, alpha=0.12, lw=0, zorder=1)
    ax.plot(amp([b.lo90, b.hi90]), [2, 2], color=CS, lw=1.8, zorder=3)
    ax.plot(amp(b.Q_star), 2, ls="", marker="o", ms=7, color=CS, zorder=4)
    for y, lab in ((1, "vs fast deaths"), (2, "vs ventilated")):
        ax.text(0.8, y + H / 2 + 0.14, lab, ha="left", va="center", fontsize=9, color=MUTED)
    ax.set_xscale("log")
    ax.set_xlim(0.78, 1.25)
    ax.set_xticks([0.8, 0.9, 1, 1.1, 1.2], ["0.8", "0.9", "1", "1.1", "1.2"])
    ax.minorticks_off()
    rows(ax, labels=False)
    ax.set_xlabel("clock amplitude, clock-specific\n[ratio to comparison]")
    return pd.DataFrame([dict(row="ventilated vs fast", Q_star=v.Q_star, amp=x,
                              matched_q95_amp=amp(-(v.matched_q95 - v.matched_median)),
                              responsive_q95_amp=amp(-(v.resp_q95 - v.matched_median))),
                         dict(row="brain dead vs ventilated", Q_star=b.Q_star, amp=amp(b.Q_star),
                              lo=amp(b.lo90), hi=amp(b.hi90), margin_lo=amp(-m), margin_hi=amp(m))])


def interval(ax, iv):
    iv = iv[iv.gene_set == "14 genes"].set_index("group")
    scale = float(iv.median_abs_pair_offset_fast_h.iloc[0])
    ax.axvline(scale, color=MUTED, lw=1, ls="--", zorder=0)
    ax.text(scale, -0.42, " size of the\n arrangement", ha="left", va="top", fontsize=9, color=MUTED)
    out = []
    for i, k in enumerate(ROWS):
        if k == "fast":
            continue
        r = iv.loc[IVL[k]]
        ax.plot([r.lower_h, r.upper_h], [i, i], color=CS, lw=4, solid_capstyle="butt")
        if r.upper_at_grid_edge:
            ax.annotate("", (r.upper_h + 0.5, i), xytext=(r.upper_h - 0.05, i),
                        arrowprops=dict(arrowstyle="-|>", color=CS, lw=1.5, mutation_scale=12))
        out.append(dict(group=k, lower_h=r.lower_h, upper_h=r.upper_h, at_edge=r.upper_at_grid_edge, scale_h=scale))
    ax.set_xlim(0, 6.6)
    ax.set_xticks([0, 2, 4, 6])
    rows(ax, labels=False)
    ax.set_xlabel("how far the organs\ncan have moved [h]")
    return pd.DataFrame(out)


def per_tissue(ax, t):
    t = t[(t.gene_set == "14 genes") & (t.design == "pooled")]
    for _, r in t.iterrows():
        c = COLOR[gtex_system(r.tissue)]
        ax.plot([r.z2_fast] * 2, [r.lz2_lo90, r.lz2_hi90], color=c, lw=1, alpha=0.6, zorder=2)
        ax.scatter([r.z2_fast], [r.log_z2_ratio], s=30, color=c, zorder=3)
        if r.tissue in SHORT:
            left = r.tissue.startswith("Skin")
            ax.annotate(SHORT[r.tissue], (r.z2_fast, r.log_z2_ratio), xytext=(-5 if left else 5, 0),
                        textcoords="offset points", ha="right" if left else "left", va="center", fontsize=9, color=c)
    ax.axhline(0, color=MUTED, lw=0.8, ls="--", zorder=0)
    ax.set_xscale("log")
    ax.set_xticks([0.03, 0.1, 0.3, 1, 3, 10], ["0.03", "0.1", "0.3", "1", "3", "10"])
    ax.minorticks_off()
    ax.set_xlabel("clock strength in fast deaths [z2]")
    ax.set_ylabel("loudness lost in intensive care\n[log, fast / ventilated]")
    rho = spearmanr(t.z2_fast, t.log_z2_ratio).statistic
    return t[["tissue", "z2_fast", "log_z2_ratio", "lz2_lo90", "lz2_hi90"]], rho


def main():
    pf.use_style("slide")
    OUT.mkdir(parents=True, exist_ok=True)
    lay = pf.Layout.slide("ABCxD", kind="wide", col_ratios=[1, 0.8, 0.8, 0.25, 1.1], blank=("x",))
    a = arrangement(lay.axes["A"], pd.read_csv(EXP / "A293_organisation_levels/output/part1/groups.csv"),
                    pd.read_csv(EXP / "A374_order_interval/output/noise_reference.csv"))
    b = amplitude(lay.axes["B"], pd.read_csv(EXP / "A335_review_batch/output/a315_core.csv"),
                  pd.read_csv(EXP / "A335_review_batch/output/equivalence.csv"))
    c = interval(lay.axes["C"], pd.read_csv(EXP / "A374_order_interval/output/intervals.csv"))
    d, rho = per_tissue(lay.axes["D"], pd.read_csv(EXP / "A349_icu_limits/output/part2_tissues.csv"))
    h = [Patch(facecolor=CS, edgecolor=CS, label="clock genes"),
         Line2D([], [], color=MATCHED, lw=2.2, label="matched genes, 95th pct"),
         Line2D([], [], color=RESPONSIVE, lw=2.2, path_effects=EDGE, label="stress-responsive genes, 95th pct"),
         Line2D([], [], color=CS, lw=1.2, marker="D", ms=6, mfc="white", mec=CS, label="expected if unchanged"),
         Patch(facecolor=CS, alpha=0.12, label="equivalence margin")]
    lay.shared_legend("ABC", h, [x.get_label() for x in h], side="bottom", ncol=3, fontsize=10, handletextpad=0.4,
                      columnspacing=1.0)
    fig = lay.finalize()
    pd.set_option("display.width", 200)
    for x in (a, b, c, d):
        print(x.round(3).to_string(index=False))
    print("Spearman loss vs fast z2:", round(rho, 3))
    print(pf.collisions(fig))
    pf.save(fig, OUT / "intensive_care")


if __name__ == "__main__":
    main()
