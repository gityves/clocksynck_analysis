"""Slide figure (results, R2.3): back to the organs. A sample's cell composition does not set its organ's phase, and it
explains only part of the organ's gain.

A. Phase. Top row: offset of each GTEx peripheral tissue from its person's body time (organ map, A104 table; one dot per
   tissue, coloured by organ system). Below: phase change of a tissue's clock state per SD of a marker score of its
   samples (A295, registered, H295.2 pass): stromal-vascular and lymphocyte scores, one dot per tissue, black: pooled
   value with 90 % interval; adrenal cortex score (A299, descriptive specificity arm), estimate and 90 % interval.
B. Gain, fold on a log scale. Top row: each tissue's gain on the body's time (relative amplitude, A104 table). Second
   row: tissues made mostly of the same cell class, higher over lower gain (post hoc pairs: squamous, oesophagus
   mucosa over skin; adipocyte, visceral over subcutaneous adipose and breast; smooth muscle, oesophagus muscularis
   over sigmoid colon). Below: gain change per SD of the marker scores (A295 descriptive for stromal-vascular, H295.1
   fail for lymphocytes; A299 cortex), as fold = 1 + change.
Status: display post hoc; every value is from a scored table. Writes the plotted rows to
presentation/output/organs_composition/rows.csv. GTEx aggregates only (one value per tissue).

Data: presentation/output/organs/couplings.csv; exploratory/A295_composition/output/{tissues,tests}.csv;
exploratory/A299_adrenal_medulla/output/tests.csv. Output: presentation/figures/organs_composition/organs_composition.{svg,png}.

    python scripts/organs_composition.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns

import paperfig as pf

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _systems import COLOR, gtex_system, legend_in  # noqa: E402

PAPER = Path(__file__).resolve().parents[2]
EXP = PAPER / "exploratory"
OUT = PAPER / "presentation" / "figures" / "organs_composition"
TAB_OUT = PAPER / "presentation" / "output" / "organs_composition"
INK, MUTED = "#222222", pf.GREY
TAB = sns.color_palette().as_hex()
SCORE = {"stromal_vascular": ("+1 SD stromal-\nvascular cells", TAB[5]), "lymphocyte": ("+1 SD lymphocytes", TAB[9])}
CORTEX = ("adrenal: +1 SD\ncortex cells", COLOR["endocrine"])
PAIRS = [("squamous", "Esophagus - Mucosa", "Skin - Sun Exposed (Lower leg)"),
         ("squamous", "Esophagus - Mucosa", "Skin - Not Sun Exposed (Suprapubic)"),
         ("adipocyte", "Adipose - Visceral (Omentum)", "Adipose - Subcutaneous"),
         ("adipocyte", "Adipose - Visceral (Omentum)", "Breast - Mammary Tissue"),
         ("smooth muscle", "Esophagus - Muscularis", "Colon - Sigmoid")]


def rows():
    """Long table of the plotted values: panel, row, label, value, interval, colour key."""
    org = pd.read_csv(PAPER / "presentation/output/organs/couplings.csv")
    org = org[org.group == "periphery"].set_index("tissue")
    tis = pd.read_csv(EXP / "A295_composition/output/tissues.csv")
    tests = pd.read_csv(EXP / "A295_composition/output/tests.csv").set_index("score")
    cx = pd.read_csv(EXP / "A299_adrenal_medulla/output/tests.csv")
    cx = cx[(cx.score == "cortex") & cx.donors.str.startswith("Hardy")].iloc[0]
    R = []
    for t, r in org.iterrows():
        R.append(dict(panel="phase", row="organ", label=t, value=r.offset_h, lo=r.offset_lo, hi=r.offset_hi,
                      key=gtex_system(t)))
        R.append(dict(panel="gain", row="organ", label=t, value=np.exp(r.log_amp_rel), lo=np.exp(r.log_amp_lo),
                      hi=np.exp(r.log_amp_hi), key=gtex_system(t)))
    for cls, a, b in PAIRS:
        R.append(dict(panel="gain", row="pair", label=f"{cls}: {a} / {b}",
                      value=np.exp(org.log_amp_rel[a] - org.log_amp_rel[b]), lo=np.nan, hi=np.nan, key=cls))
    for s in SCORE:
        t = tis[tis.score == s]
        for r in t.itertuples():
            R.append(dict(panel="phase", row=s, label=r.tissue, value=r.phase_change_h, lo=np.nan, hi=np.nan, key=s))
            R.append(dict(panel="gain", row=s, label=r.tissue, value=1 + r.gain_change, lo=np.nan, hi=np.nan, key=s))
        p = tests.loc[s]
        R.append(dict(panel="phase", row=s, label="pooled", value=p.pooled_phase_change_h, lo=p.phase_lo90,
                      hi=p.phase_hi90, key=s))
        R.append(dict(panel="gain", row=s, label="pooled", value=1 + p.pooled_gain_change,
                      lo=1 + p.pooled_gain_change - 1.645 * p.gain_se, hi=1 + p.pooled_gain_change + 1.645 * p.gain_se,
                      key=s))
    R.append(dict(panel="phase", row="cortex", label="Adrenal Gland", value=cx.phase_change_h, lo=cx.phase_lo90,
                  hi=cx.phase_hi90, key="cortex"))
    R.append(dict(panel="gain", row="cortex", label="Adrenal Gland", value=1 + cx.gain_change, lo=1 + cx.gain_lo90,
                  hi=1 + cx.gain_hi90, key="cortex"))
    return pd.DataFrame(R)


def draw(ax, R, order, labels):
    for i, row in enumerate(order):
        d = R[R.row == row]
        if row == "organ":
            y = i + np.random.default_rng(i).uniform(-0.2, 0.2, len(d))
            ax.scatter(d.value, y, s=22, c=[COLOR[k] for k in d.key], edgecolor="white", lw=0.3, zorder=2)
            continue
        if row == "pair":
            ax.scatter(d.value, [i] * len(d), s=22, color=INK, zorder=2)
            for j, (k, g) in enumerate(d.groupby("key", sort=False)):
                ax.annotate(k, (g.value.median(), i), xytext=(0, 8 if j % 2 == 0 else -9), textcoords="offset points",
                            ha="center", va="center", fontsize=7.5, color=MUTED)
            continue
        if row == "cortex":
            r = d.iloc[0]
            ax.plot([r.lo, r.hi], [i, i], color=CORTEX[1], lw=1.6, zorder=2)
            ax.scatter([r.value], [i], s=34, color=CORTEX[1], zorder=3)
            continue
        t, p = d[d.label != "pooled"], d[d.label == "pooled"].iloc[0]
        y = i + np.random.default_rng(i).uniform(-0.2, 0.2, len(t))
        ax.scatter(t.value, y, s=16, color=SCORE[row][1], alpha=0.6, edgecolor="none", zorder=1)
        ax.plot([p.lo, p.hi], [i, i], color=INK, lw=1.6, zorder=2)
        ax.scatter([p.value], [i], s=34, color=INK, zorder=3)
    ax.set_yticks(range(len(order)), labels)
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(len(order) - 0.5, -0.5)


def main():
    pf.use_style("slide")
    OUT.mkdir(parents=True, exist_ok=True)
    TAB_OUT.mkdir(parents=True, exist_ok=True)
    R = rows()
    R.to_csv(TAB_OUT / "rows.csv", index=False)
    lay = pf.Layout.slide("AB", kind="wide", col_ratios=[1, 1])
    a, b = lay.axes["A"], lay.axes["B"]
    comp = [SCORE["stromal_vascular"][0], SCORE["lymphocyte"][0], CORTEX[0]]
    draw(a, R[R.panel == "phase"], ["organ", "stromal_vascular", "lymphocyte", "cortex"], ["organ offset"] + comp)
    a.axvline(0, color=MUTED, lw=0.8, ls="--", zorder=0)
    a.set_xlim(-2.5, 5.5)
    a.set_xlabel("phase relative to the body [h]")
    ad = R[(R.panel == "phase") & (R.row == "organ") & (R.label == "Adrenal Gland")].iloc[0]
    y_ad = a.collections[0].get_offsets()[list(R[(R.panel == "phase") & (R.row == "organ")].label).index(
        "Adrenal Gland")][1]
    a.annotate("adrenal", (ad.value, y_ad), xytext=(-6, 0), textcoords="offset points", ha="right", va="center",
               fontsize=7.5, color=COLOR["endocrine"])
    legend_in(a, fontsize=7.5, loc="lower right")
    draw(b, R[R.panel == "gain"], ["organ", "pair", "stromal_vascular", "lymphocyte", "cortex"],
         ["organ gain", "same cell class,\ntwo organs"] + comp)
    b.axvline(1, color=MUTED, lw=0.8, ls="--", zorder=0)
    b.set_xscale("log")
    b.set_xlim(0.2, 3.2)
    b.set_xticks([0.25, 0.5, 1, 2], ["0.25", "0.5", "1", "2"])
    b.minorticks_off()
    b.set_xlabel("gain [fold]")
    fig = lay.finalize()
    print(pf.collisions(fig))
    pf.save(fig, OUT / "organs_composition")


if __name__ == "__main__":
    main()
