"""Slide figure (results, R1.3): in living patients, the untimed internal time follows the recorded hour.

A. Breast surgery (GSE233242, 42 patients, bulk tumour + adjacent normal = two compartments; A200 registered, pass): one dot
   per patient, internal time (host phase; the estimates' circular mean set to 13:00, no recorded time used) against the
   recorded resection hour, equal axes; dashed: identity.
B. Heart surgery (GSE270258, 16 patients, snRNA-seq cardiomyocyte, fibroblast, endothelial, neuron; 8 operated 06-10 h and 8
   at 14-18 h; A174, confirmation): each patient's internal time (host phase), by surgery window, with the circular mean.
C. Per-sample absolute error after removing one common offset (best rotation against the recorded times, both cohorts):
   breast (08-18 h, host of two compartments) vs subcutaneous adipose (GSE104674, 12 men x 4 biopsies over the day, one
   compartment, single-biopsy phase on the GTEx clock plane; A113 registered, pass). Lines: medians.

Data: presentation/output/living/breast_{patients,null}.csv (`living_breast_compute.py`), exploratory/A174_v3_heart/output/
v3_host.csv, exploratory/A113_living_time/output/a113_phases.csv. Output: presentation/figures/living_patients/
living_patients.{svg,png} (wide slide figure).

    python scripts/living_patients.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns

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
    ax.plot([8, 19], [8, 19], color=MUTED, lw=1, ls="--", zorder=1)
    y = p.resection_h + wrap(p.internal_h - p.resection_h)  # nearest turn to the recorded hour
    ax.scatter(p.resection_h, y, s=26, color=CS, edgecolor="white", lw=0.5, zorder=3)
    ax.set_xlim(8, 19)
    ax.set_ylim(8, 19)
    ax.set_xticks(range(8, 19, 2))
    ax.set_yticks(range(8, 19, 2))
    ax.set_aspect("equal")
    ax.set_xlabel("resection hour [h]")
    ax.set_ylabel("internal time [h]")
    ax.set_title("breast surgery", fontsize=13, color=INK)


def heart(ax, h):
    """Internal time per patient by surgery window, on 0-24 h (the two groups do not cross midnight)."""
    for k, (g, c, lab) in enumerate((("AM", AM, "06-10 h"), ("PM", PM, "14-18 h"))):
        x = h.phase_h[h.window == g].to_numpy()
        y = 1 - k + np.random.default_rng(k).uniform(-0.12, 0.12, len(x))
        ax.scatter(x, y, s=30, color=c, edgecolor="white", lw=0.5, zorder=3)
        m = np.mod(np.angle(np.exp(1j * W * x).mean()) / W, 24)
        ax.plot([m, m], [1 - k - 0.3, 1 - k + 0.3], color=c, lw=2.5, zorder=2)
        ax.annotate(lab, (0, 1 - k), xytext=(-6, 0), textcoords="offset points", ha="right", va="center", fontsize=11,
                    color=c, annotation_clip=False)
    ax.set_xlim(0, 24)
    ax.set_xticks(range(0, 25, 6))
    ax.set_ylim(-0.6, 1.6)
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("internal time [h]")
    ax.set_title("heart surgery", fontsize=13, color=INK)


def errors(ax, err):
    for k, (lab, e) in enumerate(err.items()):
        e = np.abs(np.asarray(e))
        x = k + np.random.default_rng(k).uniform(-0.15, 0.15, len(e))
        ax.scatter(x, e, s=14, color=CS, alpha=0.45, lw=0, zorder=2)
        ax.plot([k - 0.28, k + 0.28], [np.median(e)] * 2, color=INK, lw=2.2, zorder=3)
    ax.set_xticks(range(len(err)), list(err))
    ax.set_xlim(-0.6, len(err) - 0.4)
    ax.set_ylim(0, 6)
    ax.set_ylabel("absolute error [h]")


def offset_errors(est, true):
    d = wrap(np.asarray(est) - np.asarray(true))
    return wrap(d - np.angle(np.exp(1j * W * d).mean()) / W)


def main():
    pf.use_style("slide")
    OUT.mkdir(parents=True, exist_ok=True)
    bp = pd.read_csv(SRC / "breast_patients.csv")
    bs = pd.read_csv(SRC / "breast_summary.csv").iloc[0]
    hh = pd.read_csv(PAPER / "exploratory/A174_v3_heart/output/v3_host.csv")
    ad = pd.read_csv(PAPER / "exploratory/A113_living_time/output/a113_phases.csv").query("dataset == 'adipose GSE104674'")
    err = {"breast\n08-18 h": offset_errors(bp.internal_h, bp.resection_h),
           "adipose\nfull day": offset_errors(ad.phase_h, ad.time_h)}
    lay = pf.Layout.slide("ABC", kind="wide", col_ratios=[1, 1.2, 0.7])
    breast(lay.axes["A"], bp)
    heart(lay.axes["B"], hh)
    errors(lay.axes["C"], err)
    fig = lay.finalize()
    print(bs.round(3).to_string(), f"\nheart n {len(hh)}, adipose n {len(ad)}")
    print({k: round(float(np.median(np.abs(v))), 2) for k, v in err.items()})
    print(pf.collisions(fig))
    pf.save(fig, OUT / "living_patients")


if __name__ == "__main__":
    main()
