"""Slide figure (results, validation, 1.6, negative control): tissues of one living person, sampled once without a recorded time, share
one clock state; the person's own cell line, which carries her genome but not her body time, does not.

State: per sample, the 14 clock genes' log expression projected on the clock plane of the GTEx waveform (unweighted
projector; no model fit), centred within tissue. Angle = direction of the centred state on the clock dial [h]: the direction of a person's deviation from the cohort mean, not her time of day.
A. MuTHER (E-TABM-1140), 680 women: pointwise mutual information log2 p(x,y) / (p(x) p(y)) [bits] of the clock phase
   in subcutaneous fat (x) and skin (y; one punch biopsy), from a product von Mises kernel density on the 24 h torus
   (kernel SD 1.5 h); red: phase pairs seen together more often than their own frequencies predict, blue: less often.
   Dividing by the marginals removes the clinic-hour clustering of the phases. Grey: fewer than 0.5 persons expected
   per 1 h x 1 h cell. Dashed: identity.
B. The same for fat against the women's own lymphoblastoid cell line (794 women). A and B are descriptive
   illustrations of A255 H255.1 and H255.2 (`living_tissues_compute.py` reproduces the scored r exactly).
C. Pairing r = |sum za conj(zb)| / sqrt(sum |za|^2 sum |zb|^2) over persons: 1 when every person's state in tissue b is
   her state in tissue a turned by one common angle (the offset), 0 when unrelated. 14 clock genes (ink; A251 skin,
   A255 fat, muscle, skin, cell line) and core loop (blue; A348, 11 genes). Lines: 95th percentile of 200
   expression-matched (light grey) and 200 stress-responsive (dark grey) gene sets (14 genes: A251 / A255 sets; core
   loop: A348 calibrated broad responsive sets). The cell line is scored on 14 genes only. Filled: above both; open: not.

Status: A251 and A255 registered (not blind for the earlier GTEx pairing); A348 blind for every core-loop value; A and B
descriptive. Open data, aggregates and per-person phases of open cohorts only.
Data: presentation/output/living_tissues/muther_phases.csv; exploratory/A251_skin_sites/output/score_tests.csv;
exploratory/A255_living_tissues/output/score_tests.csv; exploratory/A285_companion_calibration/output/a348/{a251,a255}/
tests.csv.
Output: presentation/figures/living_tissues/living_tissues.{svg,png} (wide).

    python scripts/living_tissues_compute.py && python scripts/living_tissues.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import patheffects
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

import paperfig as pf

PAPER = Path(__file__).resolve().parents[2]
EXP = PAPER / "exploratory"
SRC = PAPER / "presentation" / "output" / "living_tissues"
OUT = PAPER / "presentation" / "figures" / "living_tissues"
TAB = sns.color_palette().as_hex()
INK, CORE = "#222222", TAB[0]
MATCHED, RESPONSIVE = "#BBBBBB", "#6E6E6E"
EDGE = [patheffects.Stroke(linewidth=3.6, foreground="white"), patheffects.Normal()]
BAR = 0.33
GRID, SD_H, VMAX, MIN_EXP = 48, 1.5, 1.5, 0.5  # mesh, kernel SD [h], PMI range [bits], min expected persons per 1 h cell
W = 2 * np.pi / 24
NAME = {"fat": "fat", "skin": "skin", "lcl": "own cell line"}


def pmi(d, x, y):
    """Pointwise mutual information log2 p(x,y) / (p(x) p(y)) of the two phases, from a product von Mises kernel
    density on the 24 h torus (kernel SD about SD_H hours), on a GRID x GRID mesh. Also: expected persons per 1 h x 1 h
    cell under independence, used to grey out regions with too few persons to read."""
    d = d.dropna(subset=[f"{x}_phase_h", f"{y}_phase_h"])
    a, b = (d[f"{k}_phase_h"].to_numpy() * W for k in (x, y))
    g = (np.arange(GRID) + 0.5) * 24 / GRID * W
    kappa = 1 / (SD_H * W) ** 2
    ka = np.exp(kappa * (np.cos(g[:, None] - a[None, :]) - 1))  # grid x persons
    kb = np.exp(kappa * (np.cos(g[:, None] - b[None, :]) - 1))
    pxy = ka @ kb.T
    pxy /= pxy.sum()
    px, py = pxy.sum(1), pxy.sum(0)
    v = np.log2(pxy / np.outer(px, py))
    expected = np.outer(px, py) * len(d) * (GRID / 24) ** 2
    return np.clip(v, -VMAX, VMAX), expected, len(d)


def phases(ax, d, x, y):
    v, e, k = pmi(d, x, y)
    edges = np.linspace(0, 24, GRID + 1)
    im = ax.pcolormesh(edges, edges, np.ma.masked_where(e < MIN_EXP, v).T, cmap="RdBu_r", vmin=-VMAX, vmax=VMAX,
                       rasterized=True)
    ax.set_facecolor("#E6E6E6")
    ax.plot([0, 24], [0, 24], color=pf.GREY, lw=0.8, ls="--", zorder=3)
    ax.set_xlim(0, 24)
    ax.set_ylim(0, 24)
    ax.set_xticks([0, 6, 12, 18, 24])
    ax.set_yticks([0, 6, 12, 18, 24])
    ax.set_aspect("equal")
    ax.set_xlabel(f"clock-state angle, {NAME[x]} [h]")
    ax.set_ylabel(f"clock-state angle, {NAME[y]} [h]")
    ok = e >= MIN_EXP
    return im, dict(pair=f"{x}-{y}", persons=k, shown=float(ok.mean()), diag_pmi=float(np.mean(np.diag(v)[np.diag(ok)])),
                    off_pmi=float(np.median(v[ok])))


def qline(ax, x, y, col):
    ax.plot([x - BAR / 2, x + BAR / 2], [y] * 2, color=col, lw=2.2, solid_capstyle="butt", zorder=3,
            path_effects=EDGE if col == RESPONSIVE else None)


def pairing():
    """One row per tissue pair x gene set: r, matched q95, responsive q95."""
    s = pd.read_csv(EXP / "A251_skin_sites/output/score_tests.csv")
    s.index = s.family.str.split().str[0]
    t = pd.read_csv(EXP / "A255_living_tissues/output/score_tests.csv").set_index("test")
    c = pd.concat([pd.read_csv(EXP / f"A285_companion_calibration/output/a348/{p}/tests.csv") for p in ("a251", "a255")])
    c = c.pivot(index="test", columns="family", values=["clock", "companion_q95"])
    units = [("skin,\nthree\nsites", s, "H251.1", "H251.1"), ("skin\nlayers,\nsame site", s, "H251.2a", "H251.2a"),
             ("skin\nlayers,\nother site", s, "H251.2b", "H251.2b"),
             ("surgery,\ntwo\ntissues", t, "surgical cohorts pooled", "surgical cohorts pooled"),
             ("fat-skin,\none\npunch", t, "muther fat-skin", "muther fat-skin"),
             ("fat-own\ncell line", t, "muther fat-lcl", None), ("skin-own\ncell line", t, "muther skin-lcl", None)]
    rows = []
    for lab, tab, k14, kc in units:
        r = tab.loc[k14]
        rows.append(dict(unit=lab, genes="14 genes", r=r.r, mq95=r.matched_q95, rq95=r.responsive_q95))
        if kc:
            rows.append(dict(unit=lab, genes="core loop", r=c.loc[kc, ("clock", "matched")],
                             mq95=c.loc[kc, ("companion_q95", "matched")],
                             rq95=c.loc[kc, ("companion_q95", "responsive_broad")]))
    d = pd.DataFrame(rows)
    d["ok"] = (d.r > d.mq95) & (d.r > d.rq95)
    return d


def shared(ax, d):
    order = list(dict.fromkeys(d.unit))
    for r in d.itertuples():
        x = order.index(r.unit) + (-0.19 if r.genes == "14 genes" else 0.19)
        col = INK if r.genes == "14 genes" else CORE
        ax.bar([x], [r.r], width=BAR, color=col if r.ok else "white", edgecolor=col, lw=1.4, zorder=2)
        qline(ax, x, r.mq95, MATCHED)
        qline(ax, x, r.rq95, RESPONSIVE)
    ax.set_xticks(range(len(order)), order)
    ax.tick_params(axis="x", length=0)
    ax.set_xlim(-0.6, len(order) - 0.4)
    ax.set_ylim(0, 0.9)
    ax.set_ylabel("shared clock state\nbetween two tissues, r")


def main():
    pf.use_style("slide")
    OUT.mkdir(parents=True, exist_ok=True)
    lay = pf.Layout.slide("ABC", kind="wide", col_ratios=[1, 1, 2.3])
    ph = pd.read_csv(SRC / "muther_phases.csv")
    (im, ia), (_, ib) = phases(lay.axes["A"], ph, "fat", "skin"), phases(lay.axes["B"], ph, "fat", "lcl")
    cb = pf.colorbar(lay.axes["B"], im, label="[bits]", length=0.8)
    cb.set_ticks([-1.5, 0, 1.5])
    n = pd.DataFrame([ia, ib])
    d = pairing()
    shared(lay.axes["C"], d)
    h = [Patch(facecolor=INK, edgecolor=INK, label="14 clock genes"),
         Patch(facecolor=CORE, edgecolor=CORE, label="core loop"),
         Patch(facecolor="white", edgecolor=INK, label="not clock-specific"),
         Line2D([], [], color=MATCHED, lw=2.2, label="matched genes, 95th pct"),
         Line2D([], [], color=RESPONSIVE, lw=2.2, path_effects=EDGE, label="stress-responsive genes, 95th pct")]
    lay.shared_legend("AB", [Patch(facecolor="#E6E6E6", edgecolor="none")], ["too few persons to read"],
                      side="bottom", ncol=1, fontsize=10, handletextpad=0.4)
    lay.shared_legend("C", h, [x.get_label() for x in h], side="bottom", ncol=3, fontsize=10, handletextpad=0.4,
                      columnspacing=1.0)
    fig = lay.finalize()
    pd.set_option("display.width", 200)
    print(n.round(3).to_string(index=False))
    print(d.round(3).to_string(index=False))
    print(pf.collisions(fig))
    pf.save(fig, OUT / "living_tissues")


if __name__ == "__main__":
    main()
