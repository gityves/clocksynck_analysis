"""Slide figure (results, validation): in timed living human skin and muscle, the untimed readouts equal the timed ones.

A. Phase offset between two compartments of the same person [h], clocksync (untimed, host-referenced) against the timed
   truth (cosinor on the recorded sampling time):
   - skin, dermis minus epidermis at the forearm (A251 H251.4, registered, pass): clocksync on 148 untimed women
     (Wu 2020, clinic biopsies), timed on 20 men sampled over the day (same study);
   - skin, exposed forearm minus protected buttock (A342 H342.3, 20 volunteers x 4 times, GSE326530): the offsets agree,
     but the pairing is not clock-specific against companions (H342.3 FAIL on that condition only).
B. Pre-mRNA lead over mRNA [h] (unspliced over spliced clock state, direct introns), untimed lead statistic (matched
   estimator) against the timed lead: skin protected and exposed sites (A342 H342.1-2, registered, pass), muscle (Perrin
   2018 vastus lateralis, 10 subjects x 6 times; A337 Perrin H337.1-2, BLIND for unspliced counts; core loop H337.2 FAIL).
Filled: 14 clock genes; open: core loop (11 genes, without PER1, PER2, NR1D1). Point and 90 % interval (bootstraps over
women / volunteers / subjects).
C. Hospital hours: the same cohorts cut to subsets of clock times fixed before scoring (skin, A344: full day 00, 06, 12,
   18 h; three times 06, 12, 18 h, primary; two times 06, 12 h. Muscle, A353: full day 00-20 h every 4 h; three times
   08, 12, 16 h, primary; two times 08, 12 h): untimed estimate minus the full-cycle timed truth [h], 14 clock genes, 90 %
   paired bootstrap interval (volunteers / subjects). Skin lead at each site (H344.2 pass on three times; H344.3:
   protected lead fails on two), muscle lead (H353.1-2 pass), skin exposed-minus-protected offset (H344.1 undecidable:
   point kept, interval too wide). Grey bands: registered bounds (lead +-1.5 h, offset +-1 h). A344, A353 registered,
   NOT BLIND for the full cycle, BLIND for the subsets.

Data: exploratory/A251_skin_sites/output/{score_offsets,score_descriptive}.csv;
exploratory/A342_skin_timed/output/score_direct/{pairing,leads}.csv; exploratory/A337_timed_lead/output/perrin/verdicts.csv; exploratory/A342_skin_timed/output/score_a344/subsets.csv;
exploratory/A337_timed_lead/output/a353/subsets.csv.
Output: presentation/figures/living_timed/living_timed.{svg,png} (wide).

    python scripts/living_timed.py
"""
from pathlib import Path

import numpy as np
import pandas as pd
import seaborn as sns
from matplotlib.lines import Line2D

import paperfig as pf

PAPER = Path(__file__).resolve().parents[2]
EXP = PAPER / "exploratory"
OUT = PAPER / "presentation" / "figures" / "living_timed"
CS, TIMED, INK = sns.color_palette().as_hex()[0], "#222222", "#222222"
DY = 0.16
SUBSETS = {"full": "full day", "D3": "three day times", "D2": "two day times"}
SHADE = dict(zip(SUBSETS, sns.color_palette("Blues_r", 4).as_hex()[:3]))
READOUTS = [("lead_P_diff", "lead,\nskin P", 1.5), ("lead_E_diff", "lead,\nskin E", 1.5),
            ("lead_M_diff", "lead,\nmuscle", 1.5), ("off_diff", "offset,\nskin E-P", 1.0)]


def offsets():
    """Rows: label, readout, clocksync (est, lo, hi), timed (est, lo, hi)."""
    a = EXP / "A251_skin_sites/output"
    o = pd.read_csv(a / "score_offsets.csv").query("a == 'epidermis_forearm' and b == 'dermis_forearm'").iloc[0]
    t = pd.read_csv(a / "score_descriptive.csv").set_index("item").loc[
        "timed men forearm dermis minus epidermis offset (h)"]
    rows = [dict(label="skin: dermis\n- epidermis", core=False, u=(-o.offset_h, -o.hi90, -o.lo90),
                 t=(t.value, t.lo90, t.hi90))]
    p = pd.read_csv(EXP / "A342_skin_timed/output/score_direct/pairing.csv").query("variant == 'main'")
    for _, r in p.iterrows():
        rows.append(dict(label="skin: exposed\n- protected site", core=r.waveform == "core",
                         u=(r.off_u, r.off_u_lo, r.off_u_hi), t=(r.off_t, r.off_t_lo, r.off_t_hi)))
    return rows


def leads():
    rows = []
    s = pd.read_csv(EXP / "A342_skin_timed/output/score_direct/leads.csv").query("variant == 'main'")
    for site, name in (("P", "skin, protected"), ("E", "skin, exposed")):
        for _, r in s[s.site == site].iterrows():
            rows.append(dict(label=name, core=r.waveform == "core", u=(r.matched_delta, r.matched_lo, r.matched_hi),
                             t=(r.dT, r.dT_lo, r.dT_hi)))
    m = pd.read_csv(EXP / "A337_timed_lead/output/perrin/verdicts.csv").query("unspliced == 'direct'")
    for _, r in m.iterrows():
        rows.append(dict(label="muscle", core=r.waveform == "core", u=(r.delta_U, r.dU_lo, r.dU_hi),
                         t=(r.delta_T, r.dT_lo, r.dT_hi)))
    return rows


def forest(ax, rows, xlim, xlabel):
    labels, ys, y = [], [], 0.0
    for i, r in enumerate(rows):
        if i and r["label"] != rows[i - 1]["label"]:
            y += 0.35
        for (est, lo, hi), c, dy in ((r["u"], CS, -DY), (r["t"], TIMED, DY)):
            ax.plot([lo, hi], [y + dy] * 2, color=c, lw=1.8, solid_capstyle="butt", zorder=2)
            ax.scatter([est], [y + dy], s=48, marker="o" if c == CS else "D", lw=1.4, zorder=3,
                       facecolor="white" if r["core"] else c, edgecolor=c)
        if not i or r["label"] != rows[i - 1]["label"]:
            labels.append(r["label"])
            ys.append(y + (0.5 if i + 1 < len(rows) and rows[i + 1]["label"] == r["label"] else 0))
        y += 1
    ax.axvline(0, color=pf.GREY, lw=0.6, zorder=0)
    ax.set_yticks(ys, labels)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_ylim(y - 0.5, -0.6)
    ax.set_xlim(*xlim)
    ax.set_xlabel(xlabel)


def hospital(ax, d):
    """Untimed minus full-cycle timed truth per readout, one point per subset of clock times."""
    for i, (col, _, bound) in enumerate(READOUTS):
        ax.fill_between([i - 0.42, i + 0.42], -bound, bound, color="#eeeeee", lw=0, zorder=0)
        for s, dx in zip(SUBSETS, (-0.22, 0, 0.22)):
            r = d.loc[s]
            ax.plot([i + dx] * 2, [r[f"{col}_lo"], r[f"{col}_hi"]], color=SHADE[s], lw=2, zorder=2)
            ax.scatter([i + dx], [r[col]], s=40, color=SHADE[s], zorder=3)
    ax.axhline(0, color=pf.GREY, lw=0.6, zorder=1)
    ax.set_xticks(range(len(READOUTS)), [r[1] for r in READOUTS])
    ax.tick_params(axis="x", length=0)
    ax.set_xlim(-0.55, len(READOUTS) - 0.45)
    ax.set_ylim(-2.3, 2.3)
    ax.set_ylabel("untimed - timed truth [h]")
    h = [Line2D([], [], ls="", marker="o", ms=7, color=SHADE[s], label=v) for s, v in SUBSETS.items()]
    ax.legend(h, list(SUBSETS.values()), loc="lower left", fontsize=10, frameon=False, handletextpad=0.1,
              borderaxespad=0.1)


def main():
    pf.use_style("slide")
    OUT.mkdir(parents=True, exist_ok=True)
    lay = pf.Layout.slide("ABC", kind="wide", col_ratios=[1, 1.05, 1.15])
    a, b = offsets(), leads()
    forest(lay.axes["A"], a, (-3.5, 3.5), "phase offset between\ncompartments [h]")
    forest(lay.axes["B"], b, (-0.5, 4.2), "pre-mRNA lead over mRNA [h]")
    lay.axes["A"].set_title("offset", fontsize=13, color=INK)
    lay.axes["B"].set_title("pre-mRNA lead", fontsize=13, color=INK)
    hh = pd.read_csv(EXP / "A342_skin_timed/output/score_a344/subsets.csv").query("waveform == 'all'").set_index("subset")
    mu = pd.read_csv(EXP / "A337_timed_lead/output/a353/subsets.csv").query("waveform == 'all'").set_index("subset")
    for k, c in (("", "diff"), ("_lo", "diff_lo"), ("_hi", "diff_hi")):
        hh[f"lead_M_diff{k}"] = mu[c]
    hospital(lay.axes["C"], hh)
    lay.axes["C"].set_title("hospital hours", fontsize=13, color=INK)
    h = [Line2D([], [], ls="", marker="o", ms=8, color=CS, label="clocksync (untimed)"),
         Line2D([], [], ls="", marker="D", ms=7, color=TIMED, label="timed"),
         Line2D([], [], ls="", marker="o", ms=8, color="#222222", label="14 clock genes"),
         Line2D([], [], ls="", marker="o", ms=8, mfc="white", mec="#222222", label="core loop")]
    lay.shared_legend("AB", h, [x.get_label() for x in h], side="bottom", ncol=4, fontsize=11)
    fig = lay.finalize()
    for r in a + b:
        print(r["label"].replace("\n", " "), "core" if r["core"] else "14", np.round(r["u"], 2), np.round(r["t"], 2))
    print(hh[[c for r in READOUTS for c in (r[0], f"{r[0]}_lo", f"{r[0]}_hi")]].round(2).to_string())
    print(pf.collisions(fig))
    pf.save(fig, OUT / "living_timed")


if __name__ == "__main__":
    main()
