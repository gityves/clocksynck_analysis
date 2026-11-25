"""Slide figures (results, R3.1 and R3.2): what holds the organisation.

R3.1 (`holds_organisation`), GTEx death types:
A. Does a death group keep the organ arrangement of fast deaths? Fisher-Lee correlation of the organ-pair phase
   differences of each group with those of fast deaths (12 body tissues, pairs with >= 20 shared donors): bars for the
   clock genes (blue) and the medians of 200 expression-matched (light grey) and 200 stress-responsive gene sets (dark
   grey), whisker to their 95th percentile, same donors and pairs; fast deaths are the reference (1 by construction)
   (A293 Part 1; slow deaths H293.1 and ventilated living-brain donors H293.2 registered, Hardy 3 and brain-dead organ
   donors descriptive). Rows in the same order in A and B.
B. Does the body's clock still read the external day? Alignment of the host phase (clock genes of 15 reference tissues,
   fitted without time) to the recorded hour of death, R = |mean exp(i (host phase - w t_death))|, per death type, with
   its permutation p (times shuffled over donors, 10,000): fast, slow and intermediate deaths from A27 (post hoc
   readout), ventilated donors with a living brain and brain-dead organ donors from A314 (descriptive, post hoc split).
R3.2 (`cue_removal`), designed removal or shift of a cue in one organ, as clock faces: angle = phase shift against
the reference (0 h at the top, clockwise), radius = amplitude relative to it; 90 % intervals as a radial segment and an
arc.
A. Mouse liver core clock (A297, registered, H297.1 fail above the interval) against ad libitum: feeding rhythm removed
   (Greenwell arrhythmic) and the night-restricted control (open), feeding inverted (Manella day feeding, Xin day- vs
   night-restricted; the calibration), and the state expected for removal if cues add, (1 + rho) / 2 (open diamonds, colour of removal). Colours are conditions (reference grey).
B. Mouse organs under inverted feeding (slide 4.1, introduction of section 4): Manella 2021 day-fed against ad libitum (liver, white fat, lung;
   A322, circles) and Xin 2021 day- against night-restricted (liver, heart, kidney, visceral fat; A324, squares):
   amplitude ratio (radius, 90 % library-bootstrap interval) and shift (angle) of each organ's clock-gene state.
   Dashed: the pure mix of the two drives, (1 - b) + b exp(i Delta), a straight chord from the reference to the fully
   shifted state (Delta fitted per dataset, 11.28 and 11.42 h); a pure mix puts every organ on it. Organ-system colours.
(`muscle`, human muscle under 12 h-shifted behaviour, A257 / A261, is kept for slide 4.3.)

Data: exploratory/A293_organisation_levels/output/part1/{groups,sets}.csv; results/figure_3/output/cues_lost/walltime_hosts.csv;
exploratory/A314_vent_walltime/output/groups.csv;
exploratory/A297_feeding_model/output/tests.csv; exploratory/A322_conflicting_drives/output/{ratios,fit}.csv and
xin/{ratios,fit}.csv. GTEx aggregates only.
Output: presentation/figures/holds_organisation/{holds_organisation,cue_removal}.{svg,png}.

    python scripts/holds_organisation.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib import patheffects
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

import paperfig as pf
from _systems import COLOR, gtex_system

PAPER = Path(__file__).resolve().parents[2]
EXP = PAPER / "exploratory"
OUT = PAPER / "presentation" / "figures" / "holds_organisation"
INK, MUTED = "#222222", pf.GREY
TAB = sns.color_palette().as_hex()
CS = TAB[0]
REF = "#6E6E6E"  # reference condition of a clock face (ad libitum, aligned)
FEED = {"night": TAB[2], "removed": TAB[4], "inverted": TAB[3]}  # mouse feeding arms; expected removal: open diamond
OUTPUTS = TAB[1]
FAMILY = {"matched": ("random genes", "#BBBBBB"), "responsive": ("stress-responsive genes", "#6E6E6E")}
ROWS = ["fast deaths", "slow deaths", "intermediate\ndeaths", "ventilated,\nbrain alive", "brain-dead\norgan donors"]
GROUPS = {"slow_hardy4": 1, "hardy3": 2, "ventilated_postmortem": 3, "organ_donor": 4}  # A293 group -> row
H = 0.26  # bar height in panel A


def rows(ax):
    ax.set_yticks(range(len(ROWS)), ROWS)
    ax.tick_params(axis="y", length=0)
    ax.set_ylim(len(ROWS) - 0.45, -0.55)


def arrangement(ax, g):
    """Organ arrangement shared with fast deaths: clock genes and the medians of the two companion families, whisker to
    their 95th percentile; fast deaths are the reference (1 by construction)."""
    g = g.set_index("group")
    ax.barh(0, 1, height=H, facecolor="white", edgecolor=CS, lw=1.2)
    ax.text(0.5, 0, "reference", ha="center", va="center", fontsize=8, color=CS)
    for k, i in GROUPS.items():
        r = g.loc[k]
        ax.barh(i - H, r.rho_clock, height=H, color=CS)
        for j, (fam, (_, col)) in enumerate(FAMILY.items()):
            y = i + j * H
            ax.barh(y, r[f"{fam}_median"], height=H, color=col)
            ax.plot([r[f"{fam}_median"], r[f"{fam}_q95"]], [y, y], color=col, lw=1.2)
            ax.plot([r[f"{fam}_q95"]] * 2, [y - H / 3, y + H / 3], color=col, lw=1.2)
    ax.axvline(0, color=MUTED, lw=0.8, zorder=0)
    ax.set_xlim(-0.05, 1.05)
    rows(ax)
    ax.set_xlabel("organ arrangement shared\nwith fast deaths")
    ax.set_title("organ arrangement", fontsize=13, color=INK)
    h = [Patch(color=CS, label="clock genes")] + [Patch(color=c, label=lab) for lab, c in FAMILY.values()]
    ax.legend(handles=h, loc="upper center", bbox_to_anchor=(0.45, -0.24), ncol=3, frameon=False, fontsize=8,
              handlelength=1, columnspacing=0.8)


def external(ax, w, v):
    """Alignment of the body's clock to the recorded hour of death, one bar per row: fast, slow and intermediate deaths
    from A27, the two ventilated groups from A314 (all donors of each)."""
    w, v = w.set_index("cohort"), v.set_index("group")
    R = [w.loc["fast", "R"], w.loc["hardy4", "R"], w.loc["hardy3", "R"], v.loc["ventilated_postmortem", "R"],
         v.loc["organ_donor", "R"]]
    ax.barh(range(len(ROWS)), R, height=0.55, color=CS)
    ax.set_xlim(0, 1.05)
    rows(ax)
    ax.set_yticklabels([])
    ax.set_xlabel("body clock aligned to\nthe hour of death, R")
    ax.set_title("external time", fontsize=13, color=INK)


def clockface(ax, rmax, rings, title, names=None):
    """Polar frame on a cartesian axis: angle = phase shift (0 h at the top, clockwise, 24 h a turn), radius =
    amplitude relative to the reference; grey rings at `rings`, spokes every 6 h."""
    th = np.linspace(0, 2 * np.pi, 200)
    for r in rings:
        ax.plot(r * np.sin(th), r * np.cos(th), color=MUTED, lw=0.6 if r != 1 else 0.9, ls="-" if r == 1 else ":",
                zorder=0)
        ax.text(r * np.sin(-0.6), r * np.cos(-0.6), f"{r:g}", fontsize=10, color=MUTED, ha="right", va="bottom")
    for h in (0, 6, 12, 18):
        t = 2 * np.pi * h / 24
        ax.plot([0, rmax * np.sin(t)], [0, rmax * np.cos(t)], color=MUTED, lw=0.5, zorder=0)
        lab = f"{h} h" + (f"\n{names[h]}" if names and h in names else "")
        ax.text(1.12 * rmax * np.sin(t), 1.12 * rmax * np.cos(t), lab, fontsize=11, color=MUTED, ha="center",
                va="center", linespacing=0.9)
    ax.set_aspect("equal")
    ax.set_xlim(-1.25 * rmax, 1.25 * rmax)
    ax.set_ylim(-1.25 * rmax, 1.25 * rmax)
    ax.axis("off")
    ax.set_title(title, fontsize=16, color=INK)


def vector(ax, amp, shift, color, label=None, amp_ci=None, shift_ci=None, open_=False, dx=8, dy=0, ha="left", size=40,
           marker="o"):
    """One condition: point at (amp, shift), 90 % intervals as a radial segment and an arc, direct label."""
    t = 2 * np.pi * shift / 24
    if amp_ci is not None:
        ax.plot([amp_ci[0] * np.sin(t), amp_ci[1] * np.sin(t)], [amp_ci[0] * np.cos(t), amp_ci[1] * np.cos(t)],
                color=color, lw=1.4, zorder=1)
    if shift_ci is not None:
        u = 2 * np.pi * np.linspace(shift_ci[0], shift_ci[1], 30) / 24
        ax.plot(amp * np.sin(u), amp * np.cos(u), color=color, lw=1.4, zorder=1)
    ax.scatter([amp * np.sin(t)], [amp * np.cos(t)], s=size, marker=marker, color="white" if open_ else color,
               edgecolor=color, lw=1.2, zorder=3)
    if label:
        ax.annotate(label, (amp * np.sin(t), amp * np.cos(t)), xytext=(dx, dy), textcoords="offset points",
                    ha=ha, va="center", fontsize=8, color=color, linespacing=0.9)


def legend_side(ax, entries):
    """One-column legend to the right of a clock face: (label, colour, marker, open)."""
    h = [Line2D([], [], ls="", marker=m, ms=8, mfc="white" if o else c, mec=c, mew=1.2, label=lab)
         for lab, c, m, o in entries]
    ax.legend(handles=h, loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False, fontsize=11.5,
              handletextpad=0.3, labelspacing=0.8)


def liver(ax, t):
    """Mouse liver clock state relative to ad libitum feeding (A297): amplitude ratio and phase shift."""
    clockface(ax, 1.2, (0.5, 1), "mouse liver")
    r = lambda key: t[t.contrast.str.startswith(key)].iloc[0]
    vector(ax, 1, 0, REF)
    for key in ("Manella DF", "Xin DRF"):
        x = r(key)
        vector(ax, x.amp_ratio, x.shift_h, FEED["inverted"], None, (x.amp_lo90, x.amp_hi90),
               (x.shift_lo90, x.shift_hi90), size=30)
    for key in ("prediction from Manella", "prediction from Xin"):
        x = r(key)
        vector(ax, x.amp_ratio, x.shift_h, FEED["removed"], open_=True, size=32, marker="D")
    x = r("control Greenwell night")
    vector(ax, x.amp_ratio, x.shift_h, FEED["night"], None, (x.amp_lo90, x.amp_hi90), (x.shift_lo90, x.shift_hi90))
    x = r("observed Greenwell arrhythmic")
    vector(ax, x.amp_ratio, x.shift_h, FEED["removed"], None, (x.amp_lo90, x.amp_hi90), (x.shift_lo90, x.shift_hi90))
    legend_side(ax, [("ad libitum (reference)", REF, "o", False), ("night-restricted", FEED["night"], "o", False),
                      ("feeding inverted", FEED["inverted"], "o", False),
                      ("feeding rhythm removed", FEED["removed"], "o", False),
                      ("removed: expected\nif cues add", FEED["removed"], "D", True)])


def muscle(ax, core, outputs):
    """Human muscle, misaligned against aligned (A257, A261): size and phase of the 07:00-to-19:00 change."""
    clockface(ax, 1.5, (0.5, 1, 1.5), "human muscle", {0: "central time", 12: "behaviour"})
    vector(ax, 1, 0, REF)
    vector(ax, core.R, core.s, CS, None, (core.R_lo, core.R_hi), (core.lo90, core.hi90))
    vector(ax, outputs.ratio, outputs.shift_h, OUTPUTS, None, (outputs.ratio_lo, outputs.ratio_hi),
           (outputs.shift_lo, outputs.shift_hi))
    legend_side(ax, [("aligned (reference)", REF, "o", False), ("core clock", CS, "o", False),
                      ("rhythmic outputs", OUTPUTS, "o", False)])


ORGAN = {"Liver": ("liver", "Liver"), "WAT": ("white fat", "Adipose"), "Lung": ("lung", "Lung"),
         "liver": ("liver", "Liver"), "heart": ("heart", "Heart"), "kidney": ("kidney", "Kidney"),
         "VAT": ("visceral fat", "Adipose")}  # tissue -> (label, GTEx-style name for the system colour)
DATASET = {"manella": ("o", "day-fed vs ad libitum"), "xin": ("s", "day- vs night-restricted")}
# label offsets [pt] and alignment per (dataset, tissue), clear of the neighbouring points and intervals
PLACE = {("manella", "Liver"): ((-9, 0), "right"), ("xin", "liver"): ((-9, 0), "right"),
         ("manella", "WAT"): ((9, 0), "left"), ("xin", "VAT"): ((-9, 0), "right"),
         ("manella", "Lung"): ((6, -12), "left"), ("xin", "heart"): ((0, 11), "center"),
         ("xin", "kidney"): ((9, 0), "left")}


def chord(ax, delta_h):
    """Pure mix of the two drives, (1 - b) + b exp(i Delta): a straight chord from the reference to the shifted state."""
    t = 2 * np.pi * delta_h / 24
    ax.plot([0, np.sin(t)], [1, np.cos(t)], color=MUTED, lw=1.0, ls="--", zorder=1)


def organs(ax, tabs):
    """Mouse organs under inverted feeding (A322 Manella, A324 Xin): amplitude ratio and shift against the reference
    arm, with the pure-mix chord of each dataset (an organ's expected state is the chord point at its weight)."""
    clockface(ax, 1.2, (0.5, 1), "mouse organs")
    vector(ax, 1, 0, REF)
    for name, (t, fit) in tabs.items():
        chord(ax, float(fit.delta_h.iloc[0]))
        mk = DATASET[name][0]
        for r in t.itertuples():
            lab, g = ORGAN[r.tissue]
            c = COLOR[gtex_system(g)]
            vector(ax, r.amp_ratio, r.shift_h, c, None, (r.amp_lo90, r.amp_hi90), None, marker=mk, size=40)
            off, ha = PLACE[(name, r.tissue)]
            th = 2 * np.pi * r.shift_h / 24
            ax.annotate(lab, (r.amp_ratio * np.sin(th), r.amp_ratio * np.cos(th)), xytext=off,
                        textcoords="offset points", ha=ha, va="center", fontsize=10, color=c,
                        path_effects=[patheffects.withStroke(linewidth=3, foreground="white")])
    h = [Line2D([], [], ls="", marker="o", ms=8, color=REF, label="reference arm")]
    h += [Line2D([], [], ls="", marker=m, ms=8, mfc=INK, mec=INK, label=lab) for m, lab in DATASET.values()]
    h += [Line2D([], [], color=MUTED, lw=1.0, ls="--", label="pure mix of drives")]
    ax.legend(handles=h, loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False, fontsize=11.5,
              handletextpad=0.3, labelspacing=0.8)


def main():
    pf.use_style("slide")
    OUT.mkdir(parents=True, exist_ok=True)
    p1 = EXP / "A293_organisation_levels/output/part1"
    lay = pf.Layout.slide("AB", kind="wide", col_ratios=[1.2, 1], align="all")
    arrangement(lay.axes["A"], pd.read_csv(p1 / "groups.csv"))
    external(lay.axes["B"], pd.read_csv(PAPER / "results/figure_3/output/cues_lost/walltime_hosts.csv"),
             pd.read_csv(EXP / "A314_vent_walltime/output/groups.csv"))
    fig = lay.finalize()
    print(pf.collisions(fig))
    pf.save(fig, OUT / "holds_organisation")

    lay = pf.Layout.slide("AB", kind="wide", col_ratios=[1, 1])
    liver(lay.axes["A"], pd.read_csv(EXP / "A297_feeding_model/output/tests.csv"))
    a322 = EXP / "A322_conflicting_drives/output"
    tabs = {"manella": (pd.read_csv(a322 / "ratios.csv"), pd.read_csv(a322 / "fit.csv")),
            "xin": (pd.read_csv(a322 / "xin/ratios.csv"), pd.read_csv(a322 / "xin/fit.csv"))}
    organs(lay.axes["B"], tabs)
    fig = lay.finalize()
    for name, (t, fit) in tabs.items():
        print(name, "Delta", fit.delta_h.iloc[0])
        print(t[["tissue", "amp_ratio", "amp_lo90", "amp_hi90", "shift_h", "b", "model_amp", "model_shift_h",
                 "mahalanobis"]].round(3).to_string(index=False))
    print(pf.collisions(fig))
    pf.save(fig, OUT / "cue_removal")


if __name__ == "__main__":
    main()
