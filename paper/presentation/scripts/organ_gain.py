"""Slide figure (results, adult health, 3.3): organs differ in loudness, and the difference is neither their cell
composition nor post-mortem attenuation.

A. Tissues made mostly of one cell class (post hoc choice: squamous epithelium, adipocytes, smooth muscle): each tissue's
   gain on the body's time (GTEx fast deaths, `couplings.csv`), one labelled dot per tissue; right: highest over lowest.
B. Eight GTEx tissues of the single-nucleus atlas (A323, registered; H323.1 FAIL): gain predicted from composition and
   the conserved class gains (left) against the observed gain (right), each relative to its geometric mean over the
   eight; lines join a tissue. Fold ranges 1.29 (predicted) and 2.89 (observed).
C. Post-mortem attenuation (A356, registered; H356.1, H356.2 FAIL): per tissue (dots, 25 peripheral tissues), gain of
   the donors with the longer ischaemic time over the shorter (6.7 h apart) and with the higher RIN over the lower;
   black: pooled with 90 % permutation-calibrated interval. Dashed: no change.

Status: A323 registered, not blind for the GTEx and class gains; A356 registered, blind for every split gain; panel A
post hoc display. GTEx aggregates only (one value per tissue).
Data: presentation/output/organs/couplings.csv; exploratory/A323_gain_decomposition/output/{tissues,tests}.csv;
exploratory/A356_gain_attenuation/output/{score_tissues,score_pooled}.csv.
Output: presentation/figures/organ_gain/organ_gain.{svg,png} (wide).

    python scripts/organ_gain.py
"""
from pathlib import Path

import numpy as np
import pandas as pd

import paperfig as pf
from _systems import COLOR, gtex_system
from organs_composition import CLASSES, same_class

PAPER = Path(__file__).resolve().parents[2]
EXP = PAPER / "exploratory"
OUT = PAPER / "presentation" / "figures" / "organ_gain"
INK, MUTED = "#222222", pf.GREY
SHORT = {"Heart - Left Ventricle": "heart", "Breast - Mammary Tissue": "breast",
         "Esophagus - Muscularis": "oes. muscularis", "Muscle - Skeletal": "skeletal muscle", "Lung": "lung",
         "Esophagus - Mucosa": "oes. mucosa", "Prostate": "prostate", "Skin - Sun Exposed (Lower leg)": "skin"}


def class_rows():
    org = pd.read_csv(PAPER / "presentation/output/organs/couplings.csv").set_index("tissue")
    R = [dict(row=cls, label=t, value=float(np.exp(org.log_amp_rel[t])), key=gtex_system(t))
         for cls, ts in CLASSES.items() for t in ts]
    return pd.DataFrame(R)


def decomposition(ax):
    t = pd.read_csv(EXP / "A323_gain_decomposition/output/tissues.csv")
    pred = np.log(t.G) - np.log(t.G).mean()
    obs = np.log(t.sqrtS) - np.log(t.sqrtS).mean()
    for p, o, name in zip(pred, obs, t.tissue):
        c = COLOR[gtex_system(name)]
        ax.plot([0, 1], [np.exp(p), np.exp(o)], color=c, lw=1.2, alpha=0.8, zorder=2)
        ax.scatter([0, 1], [np.exp(p), np.exp(o)], s=26, color=c, zorder=3)
    order = np.argsort(obs.to_numpy())
    ly, prev = [], -np.inf
    for i in order:  # labels at the right, each at least GAP (log) above the previous one
        prev = max(obs.to_numpy()[i], prev + 0.075)
        ly.append(prev)
    for y, i in zip(ly, order):
        ax.annotate(SHORT[t.tissue[i]], (1, np.exp(obs.to_numpy()[i])), xytext=(1.1, np.exp(y)), textcoords="data",
                    ha="left", va="center", fontsize=8, color=COLOR[gtex_system(t.tissue[i])],
                    arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.5, shrinkA=0, shrinkB=3))
    ax.set_yscale("log")
    ax.set_xlim(-0.25, 1.75)
    ax.set_ylim(0.45, 1.8)
    ax.set_yticks([0.5, 0.75, 1, 1.5, 2], ["0.5", "0.75", "1", "1.5", "2"])
    ax.minorticks_off()
    ax.set_xticks([0, 1], ["predicted\nfrom cells", "observed"])
    ax.tick_params(axis="x", length=0)
    ax.axhline(1, color=MUTED, lw=0.8, ls="--", zorder=0)
    ax.set_ylabel("organ gain [fold of mean]")
    return t


def attenuation(ax):
    s = pd.read_csv(EXP / "A356_gain_attenuation/output/score_tissues.csv")
    p = pd.read_csv(EXP / "A356_gain_attenuation/output/score_pooled.csv").set_index("covariate")
    rows = [("isch", "ischaemia\nlong / short"), ("rin", "RIN\nhigh / low")]
    for i, (cov, lab) in enumerate(rows):
        d = s[s.covariate == cov]
        y = i + np.random.default_rng(i).uniform(-0.18, 0.18, len(d))
        ax.scatter(np.exp(d.delta), y, s=16, c=[COLOR[gtex_system(t)] for t in d.tissue], alpha=0.6, edgecolor="none",
                   zorder=1)
        r = p.loc[cov]
        ax.plot(np.exp([r.lo90, r.hi90]), [i, i], color=INK, lw=1.8, zorder=2)
        ax.scatter([np.exp(r.pooled_delta)], [i], s=40, color=INK, zorder=3)
    ax.axvline(1, color=MUTED, lw=0.8, ls="--", zorder=0)
    ax.set_xscale("log")
    ax.set_xlim(0.55, 1.6)
    ax.set_xticks([0.6, 0.8, 1, 1.25, 1.5], ["0.6", "0.8", "1", "1.25", "1.5"])
    ax.minorticks_off()
    ax.set_yticks(range(len(rows)), [l for _, l in rows])
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(len(rows) - 0.5, -0.5)
    ax.set_xlabel("gain ratio within tissue [fold]")
    return p


def main():
    pf.use_style("slide")
    OUT.mkdir(parents=True, exist_ok=True)
    lay = pf.Layout.slide("ABC", kind="wide", col_ratios=[1.15, 0.8, 0.95])
    a = lay.axes["A"]
    same_class(a, class_rows())
    a.axvline(1, color=MUTED, lw=0.8, ls="--", zorder=0)
    a.set_xscale("log")
    a.set_xlim(0.6, 3.2)
    a.minorticks_off()
    a.set_xticks([0.75, 1, 1.5, 2, 3], ["0.75", "1", "1.5", "2", "3"])
    a.set_xlabel("organ gain [fold]")
    t = decomposition(lay.axes["B"])
    p = attenuation(lay.axes["C"])
    fig = lay.finalize()
    pd.set_option("display.width", 200)
    print(pd.read_csv(EXP / "A323_gain_decomposition/output/tests.csv").iloc[0][["rho", "p_exact", "range_pred",
                                                                                  "range_obs"]].to_dict())
    print(p[["pooled_delta", "lo90", "hi90", "p_one_sided", "fold"]].round(3).to_string())
    print(pf.collisions(fig))
    pf.save(fig, OUT / "organ_gain")


if __name__ == "__main__":
    main()
