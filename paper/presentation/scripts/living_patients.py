"""Slide figure (results, R1.3): in living patients, the untimed internal time follows the recorded hour.

A. Breast surgery (GSE233242, 42 patients, tumour + adjacent normal; A200 registered, pass): one dot per patient, internal
   time with the 13:00 rule (the estimates' circular mean set to 13:00, no recorded time used) against the recorded resection
   hour; dashed: identity. B. Spearman of internal time with the resection hour for the clock genes (line) and for the 200
   expression-matched non-clock gene sets (histogram, the H200.1 null).
C. Heart surgery (GSE270258, 16 patients, 8 operated 06-10 h and 8 at 14-18 h; A174, confirmation): each patient's shared
   clock state (host) in the clock plane (0 h on top, later clockwise), group means joined by an arrow; dashed: the direction
   forward time predicts (18 h).
D. Subcutaneous adipose (GSE104674, 12 men, 4 biopsies each; A113 registered, pass): single-biopsy internal time against
   zeitgeber time, rotated by the common offset; per time point circular mean and circular SD; dashed: identity.

Data: presentation/output/living/breast_{patients,null}.csv (`living_breast_compute.py`), exploratory/A174_v3_heart/output/
v3_host.csv, exploratory/A113_living_time/output/a113_phases.csv. Output: presentation/figures/living_patients/
living_patients.{svg,png} (wide slide figure).

    python scripts/living_patients.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.patches import Circle

import paperfig as pf

PAPER = Path(__file__).resolve().parents[2]
SRC = PAPER / "presentation" / "output" / "living"
OUT = PAPER / "presentation" / "figures" / "living_patients"
INK, MUTED = "#222222", pf.GREY
TAB = sns.color_palette().as_hex()
CS, AM, PM = TAB[0], TAB[9], TAB[5]
W = 2 * np.pi / 24


def wrap(h):
    return (np.asarray(h, float) + 12) % 24 - 12


def breast(ax, p):
    ax.plot([7, 19], [7, 19], color=MUTED, lw=1, ls="--", zorder=1)
    y = p.resection_h + wrap(p.internal_h - p.resection_h)  # nearest turn to the recorded hour
    ax.scatter(p.resection_h, y, s=26, color=CS, edgecolor="white", lw=0.5, zorder=3)
    ax.set_xlim(7, 19)
    ax.set_xticks(range(8, 19, 2))
    ax.set_yticks(range(6, 23, 4))
    ax.set_xlabel("resection hour [h]")
    ax.set_ylabel("internal time [h]")
    ax.set_title("breast surgery", fontsize=13, color=INK)


def breast_null(ax, null, obs):
    ax.hist(null.dropna(), bins=np.arange(-0.6, 0.65, 0.05), color="#cfcfcf", edgecolor="white", lw=0.5)
    ax.axvline(obs, color=CS, lw=2)
    ax.annotate("clock\ngenes", (obs, ax.get_ylim()[1]), xytext=(-4, -2), textcoords="offset points", ha="right", va="top",
                fontsize=10, color=CS)
    ax.annotate("matched\nnon-clock\ngenes", (np.nanmedian(null), ax.get_ylim()[1]), xytext=(0, -2),
                textcoords="offset points", ha="center", va="top", fontsize=10, color=MUTED)
    ax.set_xlabel("Spearman with\nresection hour")
    ax.set_ylabel("gene sets")
    ax.set_xlim(-0.6, 0.7)


def heart(ax, h):
    a = np.arctan2(h.y, h.x)  # host angle; drawn with 0 h on top, later clockwise
    r = np.hypot(h.x, h.y)
    X, Y = r * np.sin(a), r * np.cos(a)
    rmax = 1.15 * r.max()
    for rr in np.linspace(0, rmax, 4)[1:]:
        ax.add_patch(Circle((0, 0), rr, fc="none", ec="#eeeeee", lw=0.8, zorder=0))
    for hh, lab in ((0, "0 h"), (6, "6 h"), (12, "12 h"), (18, "18 h")):
        t = hh * W
        ax.annotate(lab, (1.08 * rmax * np.sin(t), 1.08 * rmax * np.cos(t)), ha="center", va="center", fontsize=9, color=MUTED)
    means = {}
    for g, c, lab in (("AM", AM, "06-10 h"), ("PM", PM, "14-18 h")):
        i = (h.window == g).to_numpy()
        ax.scatter(X[i], Y[i], s=28, color=c, edgecolor="white", lw=0.5, zorder=3)
        means[g] = np.array([X[i].mean(), Y[i].mean()])
        ax.scatter(*means[g], s=110, marker="D", color=c, edgecolor="white", lw=1, zorder=4)
        ax.annotate(lab, means[g], xytext=(0, 12), textcoords="offset points", ha="center", fontsize=10, color=c)
    ax.annotate("", means["PM"], xytext=means["AM"], arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.6), zorder=5)
    d = np.linalg.norm(means["PM"] - means["AM"])
    t = 18 * W
    ax.annotate("", means["AM"] + d * np.array([np.sin(t), np.cos(t)]), xytext=means["AM"],
                arrowprops=dict(arrowstyle="-|>", color=MUTED, lw=1.2, ls="--"), zorder=2)
    ax.annotate("forward time", means["AM"] + d * np.array([np.sin(t), np.cos(t)]), xytext=(0, -8), textcoords="offset points",
                ha="center", va="top", fontsize=9, color=MUTED)
    ax.set_xlim(-1.25 * rmax, 1.25 * rmax)
    ax.set_ylim(-1.25 * rmax, 1.25 * rmax)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.set_anchor("N")
    ax.set_title("heart surgery", fontsize=13, color=INK)


def adipose(ax, p):
    off = np.angle(np.exp(1j * W * (p.phase_h - p.time_h)).mean()) / W
    p = p.assign(internal=p.time_h + wrap(p.phase_h - off - p.time_h))
    ax.plot([-3, 20], [-3, 20], color=MUTED, lw=1, ls="--", zorder=1)
    jit = np.random.default_rng(0).uniform(-0.35, 0.35, len(p))
    ax.scatter(p.time_h + jit, p.internal, s=14, color=CS, alpha=0.35, lw=0, zorder=2)
    for t, g in p.groupby("time_h"):
        z = np.exp(1j * W * g.internal).mean()
        m = t + wrap(np.angle(z) / W - t)
        ax.errorbar([t], [m], yerr=np.sqrt(-2 * np.log(abs(z))) / W, fmt="o", color=CS, ms=7, lw=1.8, capsize=0, zorder=3)
    ax.set_xlim(-2, 18)
    ax.set_xticks([0, 6, 12, 18])
    ax.set_yticks([0, 6, 12, 18])
    ax.set_xlabel("zeitgeber time [h]")
    ax.set_ylabel("internal time [h]")
    ax.set_title("adipose biopsies", fontsize=13, color=INK)


def main():
    pf.use_style("slide")
    OUT.mkdir(parents=True, exist_ok=True)
    bp = pd.read_csv(SRC / "breast_patients.csv")
    bn = pd.read_csv(SRC / "breast_null.csv").spearman
    bs = pd.read_csv(SRC / "breast_summary.csv").iloc[0]
    hh = pd.read_csv(PAPER / "exploratory/A174_v3_heart/output/v3_host.csv")
    ad = pd.read_csv(PAPER / "exploratory/A113_living_time/output/a113_phases.csv").query("dataset == 'adipose GSE104674'")
    lay = pf.Layout.slide("ABCD", kind="wide", col_ratios=[1.1, 0.8, 1.1, 1.1])
    breast(lay.axes["A"], bp)
    breast_null(lay.axes["B"], bn, bs.spearman)
    heart(lay.axes["C"], hh)
    adipose(lay.axes["D"], ad)
    fig = lay.finalize()
    print(bs.round(3).to_string(), f"\nheart n {len(hh)}, adipose n {len(ad)}")
    print(pf.collisions(fig))
    pf.save(fig, OUT / "living_patients")


if __name__ == "__main__":
    main()
