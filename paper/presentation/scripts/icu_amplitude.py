"""Slide figure (results, drives, 4.1b): intensive care quiets the tissue clocks, most where they were loudest; brain
death keeps their amplitude.

A. Clock-specific clock amplitude, one row per comparison: ventilated post-mortem donors (living brain) against fast
   deaths; brain-dead organ donors against ventilated donors; organ donors with the longest against the shortest
   ventilation (tertiles). exp(Q* / 2), Q* = median over 12 body tissues of the log ratio of the debiased clock-plane
   variance minus the median of 200 expression-matched gene sets (A315; A335 part 1d, pooled donors). Lines: matched
   (light grey) and stress-responsive (dark grey) 95th percentiles on the same scale (first row; not stored for the
   tertile contrast). Brain-death row: 90 % donor-bootstrap interval and the registered equivalence margin (+-0.357 in
   Q*, shaded; A335 part 3).
B. Per tissue, loudness lost in intensive care (host-referenced log z2, fast over ventilated, at a common host spread;
   90 % donor-bootstrap interval, 500 draws; A349 D349.2, descriptive) against the tissue's loudness in fast deaths, its
   gain on the reference map of slide 3.1 (`organs_compute.py`). Organ-system colours. Post hoc display choices (user,
   2026-11-25): whole blood omitted (no readable clock, slide 3.1b); x = reference gain rather than the fast-death z2 of
   A349 D349.3, because z2 is read in each host's own scale and an anchor's leave-one-out host can differ (muscle: host
   spread 0.18 against about 3.5 for the other tissues, z2 inflated to 19.5); Spearman recomputed on the gain.
C. Per tissue, coupling to the body's clock lost in intensive care at a common host spread (d* = coupling fast minus
   ventilated at V0; 90 % donor-bootstrap interval; A349 D349.2), same x axis. Raw couplings are not compared: the host
   spread itself falls under ventilation.

Status: A315 descriptive contrasts (registered) and H315.1 FAIL (powered); A335 registered; A349 descriptive; the
amplitude display exp(Q* / 2) is post hoc. 14 genes (core loop in the speaker notes). GTEx aggregates only.
Data: exploratory/A335_review_batch/output/{a315_core,equivalence}.csv; exploratory/A315_braindeath_gain/output/
descriptive.csv; exploratory/A349_icu_limits/output/part2_tissues.csv; presentation/output/organs/couplings.csv.
Output: presentation/figures/icu_amplitude/icu_amplitude.{svg,png} (wide).

    python scripts/icu_amplitude.py
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
OUT = PAPER / "presentation" / "figures" / "icu_amplitude"
CS = sns.color_palette().as_hex()[0]
MUTED = pf.GREY
MATCHED, RESPONSIVE = "#BBBBBB", "#6E6E6E"
EDGE = [patheffects.Stroke(linewidth=3.6, foreground="white"), patheffects.Normal()]
ROWS = ["ventilated\nvs fast deaths", "brain death\nvs ventilated", "long vs short\nventilation"]
SHORT = {"Thyroid": "thyroid", "Esophagus - Mucosa": "oes. mucosa",
         "Heart - Atrial Appendage": "heart", "Muscle - Skeletal": "muscle",
         "Skin - Not Sun Exposed (Suprapubic)": "skin"}
# label offsets [pt] and alignment, placed clear of the neighbouring error bars
PLACE = {"Esophagus - Mucosa": ((16, -16), "left"), "Heart - Atrial Appendage": ((-7, 0), "right"),
         "Skin - Not Sun Exposed (Suprapubic)": ((-6, 0), "right")}
H = 0.55


def amp(q):
    """Clock-plane variance ratio (log) -> amplitude ratio."""
    return np.exp(np.asarray(q, float) / 2)


def qline(ax, y, x, col):
    ax.plot([x, x], [y - H / 2, y + H / 2], color=col, lw=2.2, solid_capstyle="butt", zorder=3,
            path_effects=EDGE if col == RESPONSIVE else None)


def amplitude(ax, core, eq, desc):
    q = core[core.gene_set == "14 genes"].set_index("test")
    v = q.loc["ICU fast vs ventilated PM"]  # Q = log(fast / ventilated): ventilated / fast is -Q*
    b = eq[(eq.gene_set == "14 genes") & (eq.donors == "pooled")].iloc[0]  # Q = log(brain dead / ventilated)
    t = desc.set_index("contrast").loc["od_top_vs_bottom_tertile"]  # Q = log(longest / shortest tertile)
    x = [amp(-v.Q_star), amp(b.Q_star), amp(t.Q_star)]
    ax.axvline(1, color=MUTED, lw=0.8, ls="--", zorder=0)
    m = float(b.margin)
    ax.fill_betweenx([1 - H / 2, 1 + H / 2], amp(-m), amp(m), color=CS, alpha=0.12, lw=0, zorder=1)
    for i, xi in enumerate(x):
        ax.barh(i, xi - 1, left=1, height=H, color=CS, zorder=2)
    qline(ax, 0, amp(-(v.matched_q95 - v.matched_median)), MATCHED)
    qline(ax, 0, amp(-(v.resp_q95 - v.matched_median)), RESPONSIVE)
    ax.plot(amp([b.lo90, b.hi90]), [1, 1], color="#222222", lw=1.6, zorder=4)
    ax.set_xscale("log")
    ax.set_xlim(0.78, 1.25)
    ax.set_xticks([0.8, 0.9, 1, 1.1, 1.2], ["0.8", "0.9", "1", "1.1", "1.2"])
    ax.minorticks_off()
    ax.set_yticks(range(len(ROWS)), ROWS)
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(len(ROWS) - 0.5, -0.5)
    ax.set_xlabel("clock amplitude, clock-specific\n[ratio to comparison]")
    return pd.DataFrame([dict(row="ventilated vs fast", Q_star=v.Q_star, amp=x[0],
                              matched_q95=amp(-(v.matched_q95 - v.matched_median)),
                              responsive_q95=amp(-(v.resp_q95 - v.matched_median))),
                         dict(row="brain death vs ventilated", Q_star=b.Q_star, amp=x[1], lo=amp(b.lo90),
                              hi=amp(b.hi90), margin_lo=amp(-m), margin_hi=amp(m)),
                         dict(row="long vs short ventilation", Q_star=t.Q_star, amp=x[2],
                              pct_matched=t.percentile_vs_matched)])


def tissues(t, gain):
    t = t[(t.gene_set == "14 genes") & (t.design == "pooled") & (t.tissue != "Whole Blood")]  # no readable clock (3.1b)
    return t.assign(gain=t.tissue.map(gain.rel_amplitude))


def per_tissue(ax, t, y, lo, hi, labels=SHORT):
    """One dot per tissue: y with its 90 % interval against the reference gain."""
    for _, r in t.iterrows():
        c = COLOR[gtex_system(r.tissue)]
        ax.plot([r.gain] * 2, [r[lo], r[hi]], color=c, lw=1, alpha=0.6, zorder=2)
        ax.scatter([r.gain], [r[y]], s=30, color=c, zorder=3)
        if r.tissue in labels:
            off, ha = PLACE.get(r.tissue, ((7, 0), "left"))
            ax.annotate(SHORT[r.tissue], (r.gain, r[y]), xytext=off, textcoords="offset points", ha=ha, va="center",
                        fontsize=9, color=c, arrowprops=dict(arrowstyle="-", color=c, lw=0.6, shrinkA=0, shrinkB=4)
                        if abs(off[1]) > 0 else None)
    ax.axhline(0, color=MUTED, lw=0.8, ls="--", zorder=0)
    ax.set_xlabel("loudness in fast deaths\n[gain, reference map]")
    return spearmanr(t.gain, t[y]).statistic


def main():
    pf.use_style("slide")
    OUT.mkdir(parents=True, exist_ok=True)
    lay = pf.Layout.slide("AxByC", kind="wide", col_ratios=[1, 0.12, 0.9, 0.12, 0.9], blank=("x", "y"))
    a = amplitude(lay.axes["A"], pd.read_csv(EXP / "A335_review_batch/output/a315_core.csv"),
                  pd.read_csv(EXP / "A335_review_batch/output/equivalence.csv"),
                  pd.read_csv(EXP / "A315_braindeath_gain/output/descriptive.csv"))
    d = tissues(pd.read_csv(EXP / "A349_icu_limits/output/part2_tissues.csv"),
                pd.read_csv(PAPER / "presentation/output/organs/couplings.csv", index_col="tissue"))
    rho = per_tissue(lay.axes["B"], d, "log_z2_ratio", "lz2_lo90", "lz2_hi90")
    lay.axes["B"].set_ylabel("loudness lost in intensive care\n[log, fast / ventilated]")
    rho_c = per_tissue(lay.axes["C"], d, "d_star", "dstar_lo90", "dstar_hi90",
                       labels=("Muscle - Skeletal", "Thyroid"))
    lay.axes["C"].set_ylabel("coupling lost in intensive care\n[fast minus ventilated]")
    h = [Patch(facecolor=CS, edgecolor=CS, label="clock genes"),
         Line2D([], [], color=MATCHED, lw=2.2, label="matched genes, 95th pct"),
         Line2D([], [], color=RESPONSIVE, lw=2.2, path_effects=EDGE, label="stress-responsive genes, 95th pct"),
         Line2D([], [], color="#222222", lw=1.6, label="90 % interval"),
         Patch(facecolor=CS, alpha=0.12, label="equivalence margin")]
    lay.shared_legend("A", h, [x.get_label() for x in h], side="bottom", ncol=2, fontsize=10, handletextpad=0.4,
                      columnspacing=1.0)
    fig = lay.finalize()
    pd.set_option("display.width", 200)
    print(a.round(3).to_string(index=False))
    print(d[["tissue", "gain", "log_z2_ratio", "lz2_lo90", "lz2_hi90", "d_star", "dstar_lo90", "dstar_hi90"]]
          .round(3).to_string(index=False))
    print("Spearman loudness loss vs gain:", round(rho, 3), " coupling loss vs gain:", round(rho_c, 3),
          " median d*:", round(d.d_star.median(), 3))
    print(pf.collisions(fig))
    pf.save(fig, OUT / "icu_amplitude")


if __name__ == "__main__":
    main()
