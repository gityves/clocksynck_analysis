"""Slide figure (results, R1.5): the readout moves when the clock is broken or shifted (mouse positive controls, A111/A129).

A. Hepatocyte Bmal1 knockout vs AlbCre control (Manella 2021, ad libitum, 24 vs 24 mice; H111.3 registered, pass): change in
   coupling of each tissue at the control's host spread (d*), point and 95 % interval from the permutations.
B. Night to day restricted feeding (Xin 2021, 4 tissues; H111.1 registered, pass): change of each tissue's untimed offset
   against the change of its timed offset; dashed: identity.

Data: exploratory/A111_feeding_controls/output/a111_h111_3.csv (primary, ad libitum rows), a111_h111_1.csv (xin2021).
Output: presentation/figures/positive_controls/positive_controls.{svg,png} (half-slide figure).

    python scripts/positive_controls.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns

import paperfig as pf

PAPER = Path(__file__).resolve().parents[2]
SRC = PAPER / "exploratory/A111_feeding_controls/output"
OUT = PAPER / "presentation" / "figures" / "positive_controls"
INK, MUTED = "#222222", pf.GREY
CS = sns.color_palette().as_hex()[0]
NAMES = {"Liver": "liver", "WAT": "white fat", "Lung": "lung", "liver": "liver", "VAT": "visceral fat", "heart": "heart",
         "kidney": "kidney"}


def knockout(ax, t):
    t = t.iloc[::-1].reset_index(drop=True)
    y = np.arange(len(t))
    ax.axvline(0, color=MUTED, lw=1, ls="--", zorder=1)
    ax.errorbar(t.d_star, y, xerr=[t.d_star - t.d_star_lo, t.d_star_hi - t.d_star], fmt="o", color=CS, ms=8, lw=1.8,
                capsize=0, zorder=3)
    ax.set_yticks(y, [NAMES[k] for k in t.target])
    ax.set_ylim(-0.6, len(t) - 0.4)
    ax.set_xlim(-0.5, 0.2)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlabel("change in coupling,\nknockout - control")
    ax.set_title("liver Bmal1 knockout", fontsize=13, color=INK)


def feeding(ax, r):
    tis = ["liver", "VAT", "heart", "kidney"]
    x = np.array([r[f"timed_{k}"] for k in tis], float)
    y = np.array([r[f"untimed_{k}"] for k in tis], float)
    ax.plot([-4, 11], [-4, 11], color=MUTED, lw=1, ls="--", zorder=1)
    ax.scatter(x, y, s=50, color=CS, edgecolor="white", lw=0.5, zorder=3)
    for k, a, b in zip(tis, x, y):
        ax.annotate(NAMES[k], (a, b), xytext=(7, -2), textcoords="offset points", fontsize=10, color=INK, va="top")
    ax.set_xlim(-4, 11)
    ax.set_ylim(-4, 11)
    ax.set_xticks([-3, 0, 3, 6, 9])
    ax.set_yticks([-3, 0, 3, 6, 9])
    ax.set_aspect("equal")
    ax.set_xlabel("timed offset change [h]")
    ax.set_ylabel("untimed offset change [h]")
    ax.set_title("day feeding", fontsize=13, color=INK)


def main():
    pf.use_style("slide")
    OUT.mkdir(parents=True, exist_ok=True)
    k = pd.read_csv(SRC / "a111_h111_3.csv")
    k = k[(k.dataset == "manella2021") & (k.variant == "primary") & k.contrast.str.startswith("AL:")]
    f = pd.read_csv(SRC / "a111_h111_1.csv").query("dataset == 'xin2021' and variant == 'primary'").iloc[0]
    lay = pf.Layout.slide("AB", kind="half", height=100)
    knockout(lay.axes["A"], k)
    feeding(lay.axes["B"], f)
    fig = lay.finalize()
    print(k[["target", "d_star", "d_star_lo", "d_star_hi"]].round(3).to_string(index=False))
    print(f[["r_fl", "p_exact"]].to_string())
    print(pf.collisions(fig))
    pf.save(fig, OUT / "positive_controls")


if __name__ == "__main__":
    main()
