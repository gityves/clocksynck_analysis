"""Slide figures (results, local disease, 6.1a bulk and 6.1b single cells): tumours keep their patient's time, carried
by their fibroblasts, and express less of it as they proliferate.

6.1a (bulk, `tumours_bulk`)
A. Tumour-normal pairing (TCGA, 517 pairs, 10 cancers; A138 H138.4, registered): Poisson-corrected |correlation| of the
   tumour's and its own patient's normal clock state, bars per cancer (number of pairs in the label; dark ticks: the
   cancer's pairing-null mean) and pooled (the pair-weighted mean; 14 genes filled, core loop open with its 90 %
   bootstrap interval, A336; dark tick: pooled pairing-null q95; light ticks: expression-matched set q95).
B. Proliferation (A304, registered): high- against low-proliferation tumours, pooled log ratio of gain (the tumour's
   expression of its patient's time) and of departure (what it does on its own), interval +-1.645 null SD; grey: the
   5-95 % range of 200 expression-matched sets.

6.1b (single cells, full model, `tumours_sc`; A387, registered, blind): one model over 34 tumour / adjacent-normal
cohorts, host from the patient's normal-sample stroma (fibroblasts, endothelium, mural), cohort x chemistry strata;
rows fibroblasts, endothelium, cancer cells (patients with a host in the label); 14 genes filled, core loop open;
90 % donor-bootstrap intervals (grey).
A. Coupling c of each tumour compartment; light ticks: matched-set q95 (none for cancer cells); grey dots: per cancer
   on the pooled host (colorectal, lung, liver, prostate; A387 note 1, post hoc implementation, descriptive).
B. Peak time of each tumour compartment minus its normal counterpart [h] (+ = tumour later); dashed: the registered
   +-1.5 h bound of "keeps the patient's time".
C. Loudness ratio, tumour / normal counterpart (log axis).

Status: A138 H138.4 pass (core loop at the line, A336); A304 H304.1 gain pass, departure fail (powered); A387 G387 pass,
H387.3 fibroblasts pass on both readouts, endothelium pass on 14 genes and fail on the core loop, H387.4 fibroblasts pass;
cancer cells descriptive.
Data: exploratory/A138_cancer_estimability/output/score/{tests,per_unit}.csv; exploratory/A336_two_readouts/output/
part1/a138.csv; exploratory/A304_cancer_gain_departure/output/score/tests.csv; exploratory/A387_tumour_full_model/
output/score/{part2_tests,per_cancer_tests}.csv. Output: presentation/figures/tumours/tumours_{bulk,sc}.{svg,png}.

    python scripts/tumours.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D

import paperfig as pf

PAPER = Path(__file__).resolve().parents[2]
EXP = PAPER / "exploratory"
OUT = PAPER / "presentation" / "figures" / "tumours"
CS = sns.color_palette().as_hex()[0]
MUTED, MATCHED, INK = pf.GREY, "#BBBBBB", "#222222"
CANCER = {"BLCA": "bladder", "BRCA": "breast", "COAD": "colon", "KIRC": "kidney, clear cell", "KIRP": "kidney, papillary",
          "LIHC": "liver", "LUAD": "lung", "PRAD": "prostate", "THCA": "thyroid", "UCEC": "endometrium"}
COMP = {"fibroblast": "fibroblasts", "endothelial": "endothelium", "epithelium": "cancer cells"}
BOUND_H = 1.5


def tick(ax, x, y, h=0.3, color=MATCHED, zorder=1):
    ax.plot([x, x], [y - h, y + h], color=color, lw=2, zorder=zorder)


def pairing(ax):
    t = pd.read_csv(EXP / "A138_cancer_estimability/output/score/tests.csv").set_index("hypothesis").loc["H138.4"]
    u = pd.read_csv(EXP / "A138_cancer_estimability/output/score/per_unit.csv")
    u = u[u.hypothesis == "H138.4"].sort_values("statistic", ascending=False)
    core = pd.read_csv(EXP / "A336_two_readouts/output/part1/a138.csv").set_index("id").loc["A138_H138.4"]
    hb = 0.7
    ax.barh(range(len(u)), u.statistic, height=hb, color=CS, zorder=2)
    for i, r in enumerate(u.itertuples()):
        tick(ax, r.null_mean, i, hb / 2, color=INK, zorder=4)
    y = len(u) + 0.6
    ax.barh([y], [t.statistic], height=hb, color=CS, zorder=2)
    lo, hi = (float(x) for x in str(core.notes).split("core 90 % bootstrap ")[1].split(";")[0].split(" to "))
    ax.barh([y + 0.85], [core.value_core], height=hb, facecolor="white", edgecolor=CS, lw=1.3, zorder=2)
    ax.plot([lo, hi], [y + 0.85] * 2, color=CS, lw=1.4, zorder=3)
    tick(ax, t.null_q95, y, hb / 2, color=INK, zorder=4)
    tick(ax, t.companion_q95, y, hb / 2, zorder=4)
    tick(ax, float(str(core.notes).split("matched q95 ")[1].split(",")[0]), y + 0.85, hb / 2)
    ax.set_yticks(list(range(len(u))) + [y, y + 0.85],
                  [f"{CANCER[x]} ({n})" for x, n in zip(u.unit, u.n)]
                  + [f"pooled ({int(t.n_units)})", "pooled, core loop"])
    ax.set_ylim(y + 1.4, -0.6)
    ax.set_xlim(0, 0.62)
    ax.set_xlabel("tumour follows own normal\n[pairing, |r|]")
    ax.tick_params(axis="y", length=0)
    return u[["unit", "statistic", "null_mean"]], t, core


def sc_rows():
    t = pd.read_csv(EXP / "A387_tumour_full_model/output/score/part2_tests.csv")
    t["readout"] = t.readout.astype(str)
    return t


ROWS = ((-0.15, "14", False), (0.15, "core", True))


def sc_axis(ax, t, labels):
    ax.set_ylim(len(COMP) - 0.5, -0.5)
    if labels:
        ax.set_yticks(range(len(COMP)), [f"{lab} ({int(t[(t.kind == k) & (t.readout == '14')].n.iloc[0])})"
                                         for k, lab in COMP.items()])
    else:
        ax.set_yticks(range(len(COMP)), [""] * len(COMP))
    ax.tick_params(axis="y", length=0)


def coupling(ax, t):
    pc = pd.read_csv(EXP / "A387_tumour_full_model/output/score/per_cancer_tests.csv")
    pc = pc[pc.cancer != "other"].assign(readout=lambda d: d.readout.astype(str))
    hb = 0.28
    for i, k in enumerate(COMP):
        for dy, rd, open_ in ROWS:
            r = t[(t.kind == k) & (t.readout == rd)].iloc[0]
            ax.barh([i + dy], [r.coupling], height=hb, facecolor="white" if open_ else CS, edgecolor=CS, lw=1.3, zorder=2)
            ax.plot([r.coupling_lo90, r.coupling_hi90], [i + dy] * 2, color=MUTED, lw=1.4, zorder=3)
            if np.isfinite(r.matched_q95):
                tick(ax, r.matched_q95, i + dy, hb / 2, zorder=4)
            d = pc[(pc.kind == k) & (pc.readout == rd)]
            ax.scatter(d.coupling, np.full(len(d), i + dy), s=16, color=INK, alpha=0.55, lw=0, zorder=5)
    sc_axis(ax, t, True)
    ax.set_xlim(0, 0.9)
    ax.set_xlabel("coupled to the patient's\nnormal stroma [c]")
    return pc[["kind", "readout", "cancer", "n", "coupling"]]


def phase(ax, t):
    for i, k in enumerate(COMP):
        for dy, rd, open_ in ROWS:
            r = t[(t.kind == k) & (t.readout == rd)].iloc[0]
            ax.plot([r.peak_lo90, r.peak_hi90], [i + dy] * 2, color=MUTED, lw=1.4, zorder=2)
            ax.scatter([r.peak_h], [i + dy], s=40, facecolor="white" if open_ else CS, edgecolor=CS, lw=1.3, zorder=3)
    for x in (-BOUND_H, BOUND_H):
        ax.axvline(x, color=MUTED, lw=0.8, ls="--", zorder=0)
    ax.axvline(0, color=MATCHED, lw=0.8, zorder=0)
    sc_axis(ax, t, False)
    ax.set_xlim(-3.6, 3.6)
    ax.set_xlabel("peak time, tumour minus\nnormal counterpart [h]")


def loudness(ax, t):
    for i, k in enumerate(COMP):
        for dy, rd, open_ in ROWS:
            r = t[(t.kind == k) & (t.readout == rd)].iloc[0]
            ax.plot(np.exp([r.log_amp_lo90, r.log_amp_hi90]), [i + dy] * 2, color=MUTED, lw=1.4, zorder=2)
            ax.scatter([np.exp(r.log_amp_ratio)], [i + dy], s=40, facecolor="white" if open_ else CS, edgecolor=CS,
                       lw=1.3, zorder=3)
    ax.axvline(1, color=MATCHED, lw=0.8, zorder=0)
    ax.set_xscale("log")
    ax.set_xlim(0.12, 1.6)
    ax.set_xticks([0.125, 0.25, 0.5, 1], ["1/8", "1/4", "1/2", "1"])
    ax.minorticks_off()
    sc_axis(ax, t, False)
    ax.set_xlabel("loudness, tumour / normal\ncounterpart [ratio]")


def proliferation(ax):
    t = pd.read_csv(EXP / "A304_cancer_gain_departure/output/score/tests.csv")
    t = t[t.contrast == "H304.1"].set_index("readout")
    for i, (k, lab) in enumerate((("gain", "expression of the\npatient's time"), ("departure", "own departure"))):
        r = t.loc[k]
        ax.plot([r.matched_q05, r.matched_q95], [i + 0.22] * 2, color=MATCHED, lw=2.2, zorder=1)
        h = 1.645 * r.null_sd
        ax.plot([r.statistic - h, r.statistic + h], [i, i], color=CS, lw=1.6, zorder=2)
        ax.scatter([r.statistic], [i], s=40, color=CS, zorder=3)
    ax.axvline(0, color=MUTED, lw=0.8, ls="--", zorder=0)
    ax.set_yticks([0, 1], ["expression of the\npatient's time", "own departure"])
    ax.set_ylim(1.6, -0.6)
    ax.set_xlim(-0.85, 0.5)
    ax.set_xlabel("high / low proliferation\n[log ratio]")
    ax.tick_params(axis="y", length=0)
    return t[["statistic", "null_sd", "p", "matched_q05", "matched_q95"]]


def main():
    pf.use_style("slide")
    OUT.mkdir(parents=True, exist_ok=True)
    # 6.1a bulk
    lay = pf.Layout.slide("AxB", kind="wide", col_ratios=[1.05, 0.15, 0.9], blank=("x",))
    a, t, core = pairing(lay.axes["A"])
    c = proliferation(lay.axes["B"])
    h = [Line2D([], [], ls="", marker="o", ms=7, color=CS, label="clock genes"),
         Line2D([], [], ls="", marker="o", ms=7, mfc="white", mec=CS, mew=1.3, label="core loop"),
         Line2D([], [], color=INK, lw=2, label="pairing null (pooled: 95th pct)"),
         Line2D([], [], color=MATCHED, lw=2.2, label="matched genes (95th pct; B: 5-95 %)")]
    lay.shared_legend("AxB", h, [x.get_label() for x in h], side="bottom", ncol=4, fontsize=9.5, handletextpad=0.3)
    fig = lay.finalize()
    pd.set_option("display.width", 200)
    print(a.round(3).to_string(index=False))
    print("pooled", round(t.statistic, 3), "null q95", round(t.null_q95, 3), "matched q95", round(t.companion_q95, 3),
          "core", round(core.value_core, 3))
    print(c.round(3).to_string())
    print("bulk", pf.collisions(fig))
    pf.save(fig, OUT / "tumours_bulk")
    # 6.1b single cells
    lay = pf.Layout.slide("ABC", kind="wide", col_ratios=[1.15, 0.9, 0.9])
    st = sc_rows()
    pc = coupling(lay.axes["A"], st)
    phase(lay.axes["B"], st)
    loudness(lay.axes["C"], st)
    h = [Line2D([], [], ls="", marker="o", ms=7, color=CS, label="clock genes"),
         Line2D([], [], ls="", marker="o", ms=7, mfc="white", mec=CS, mew=1.3, label="core loop"),
         Line2D([], [], color=MUTED, lw=1.4, label="90 % interval"),
         Line2D([], [], color=MATCHED, lw=2.2, label="matched genes, 95th pct"),
         Line2D([], [], ls="", marker="o", ms=5, color=INK, alpha=0.55, label="per cancer"),
         Line2D([], [], color=MUTED, lw=0.8, ls="--", label="\u00b11.5 h")]
    lay.shared_legend("ABC", h, [x.get_label() for x in h], side="bottom", ncol=6, fontsize=9.5, handletextpad=0.3)
    fig = lay.finalize()
    print(st[["kind", "readout", "n", "coupling", "peak_h", "peak_lo90", "peak_hi90", "log_amp_ratio"]].round(3)
          .to_string(index=False))
    print(pc.round(3).to_string(index=False))
    print("sc", pf.collisions(fig))
    pf.save(fig, OUT / "tumours_sc")


if __name__ == "__main__":
    main()
