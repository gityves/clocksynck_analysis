"""Slide schematic (introduction, slide 4.4): a hospital-hours window forced onto a full circle.

Patients are sampled once between 10 and 18 h; tumour and normal tissue of one patient share the sampling time, so the
true tumour - normal phase difference is 0 h for every patient. An orderer that assumes the samples cover the 24 h cycle
spreads each tissue's samples over the whole circle (here: rank -> uniform phase, each tissue ordered separately, as
CoFE orders tumour and normal). The 8 h window is stretched threefold, and with it the estimation noise, which is larger
in the weaker tumour clock: the inferred tumour - normal differences spread over the cycle and read as "uncoupling".
Illustrative simulation, not fitted to any data.

Output: paper/presentation/figures/hospital_hours_ordering/hospital_hours_ordering.{svg,png} (320 x 110 mm, "wide").

    python scripts/hospital_hours_ordering.py
"""
from pathlib import Path

import numpy as np

import paperfig as pf

OUT = Path(__file__).resolve().parents[1] / "figures" / "hospital_hours_ordering"
N, WINDOW = 40, (10.0, 18.0)  # patients, sampling window (h)
SD = {"normal": 0.5, "tumour": 1.5}  # error of the phase each tissue's expression carries (h); weaker tumour clock
R = {"normal": 1.0, "tumour": 0.72}  # ring of each tissue on the dials
MARKER = {"normal": "o", "tumour": "^"}
SHOW = [3, 17, 31]  # patients followed across panels (indices in sampling order)
INK, MUTED, BAND = "#222222", pf.GREY, "#f2f2f2"


def simulate(seed=1):
    rng = np.random.default_rng(seed)
    t = np.sort(rng.uniform(*WINDOW, N))
    observed = {k: t + rng.normal(0, sd, N) for k, sd in SD.items()}
    mid = np.mean(WINDOW)  # the orderer's origin is aligned on the cohort (as with a reference acrophase)
    inferred = {k: (mid + 24 * ((np.argsort(np.argsort(o)) + 0.5) / N - 0.5)) % 24 for k, o in observed.items()}
    diff = (inferred["tumour"] - inferred["normal"] + 12) % 24 - 12
    return t, inferred, diff


def rad(h):
    return 2 * np.pi * np.asarray(h) / 24


def heading(ax, title, sub):
    """Title and grey sub-title at fixed offsets above the axes top (the dials and the histogram share it)."""
    for text, dy, size, color in ((title, 50, 13, INK), (sub, 30, 11, MUTED)):
        ax.annotate(text, (0.5, 1), xytext=(0, dy), xycoords="axes fraction", textcoords="offset points", ha="center",
                    va="bottom", fontsize=size, color=color, annotation_clip=False)


def dial(ax, title, sub):
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_ylim(0, 1.18)
    ax.set_yticks([])
    ax.set_xticks(rad([0, 6, 12, 18]), ["00 h", "06", "12", "18"], color=MUTED)
    ax.grid(False)
    ax.spines["polar"].set_color(MUTED)
    ax.spines["polar"].set_linewidth(0.8)
    heading(ax, title, sub)


def tissues(ax, hours):
    for k in ("normal", "tumour"):
        ax.scatter(rad(hours[k]), np.full(N, R[k]), s=34, marker=MARKER[k], facecolor="white", edgecolor=MUTED,
                   linewidth=0.9, zorder=3)
    for j, i in enumerate(SHOW):
        c = pf.OKABE_ITO[j]
        th = rad([hours["normal"][i], hours["tumour"][i]])
        ax.plot(th, [R["normal"], R["tumour"]], color=c, lw=1.6, zorder=4)
        for k, th_k in zip(("normal", "tumour"), th):
            ax.scatter(th_k, R[k], s=60, marker=MARKER[k], color=c, edgecolor="white", linewidth=0.8, zorder=5)


def main():
    pf.use_style("slide")
    t, inferred, diff = simulate()
    lay = pf.Layout.slide("ABC", "wide", polar="AB", col_ratios=[1, 1, 1.25], wgap=14)
    A, B, C = (lay.axes[k] for k in "ABC")

    dial(A, "Sampled in hospital hours", "tumour and normal of a patient: same time")
    A.bar(rad(np.mean(WINDOW)), 1.18, width=rad(np.ptp(WINDOW)), color=BAND, zorder=0)
    tissues(A, {"normal": t, "tumour": t})
    A.text(rad(14), 0.38, "10-18 h", ha="center", va="center", fontsize=11, color=MUTED)
    handles = [A.scatter([], [], s=40, marker=MARKER[k], facecolor="white", edgecolor=INK, linewidth=0.9) for k in R]
    A.legend(handles, ["normal", "tumour"], loc="center", bbox_to_anchor=(0.72, 0.74), fontsize=11, handletextpad=0.3)

    dial(B, "Ordered as if they covered 24 h", "each tissue on its own: 8 h stretched to 24 h")
    tissues(B, inferred)

    bins = np.arange(-12, 12.1, 2)
    C.hist(diff, bins=bins, color="#bdbdbd", edgecolor="white", linewidth=2, zorder=2)
    top = np.histogram(diff, bins)[0].max()
    C.axvline(0, color=INK, lw=1.6, ls=(0, (4, 3)), zorder=3)
    C.text(0.4, top * 1.32, "truth: 0 h for every patient", ha="left", va="center", fontsize=11, color=INK)
    C.text(-11.6, top * 1.12, "3 patients:", ha="left", va="center", fontsize=11, color=MUTED)
    for j, i in enumerate(SHOW):
        C.scatter(diff[i], top * 1.12, s=60, marker="v", color=pf.OKABE_ITO[j], edgecolor="white", linewidth=0.8, zorder=4)
    C.set_xlim(-12, 12)
    C.set_ylim(0, top * 1.45)
    C.set_xticks([-12, -6, 0, 6, 12])
    C.yaxis.get_major_locator().set_params(integer=True)
    C.set_xlabel("inferred tumour − normal (h)")
    C.set_ylabel("patients")
    heading(C, "Tumour − normal difference", "spread by several hours: reads as uncoupling")

    fig = lay.finalize()
    OUT.mkdir(parents=True, exist_ok=True)
    pf.save(fig, OUT / "hospital_hours_ordering")
    sd = np.sqrt(-2 * np.log(np.abs(np.mean(np.exp(1j * rad(diff))))))
    print(f"circular SD of inferred differences: {24 * sd / (2 * np.pi):.1f} h (truth 0 h); |diff| > 3 h in "
          f"{np.mean(np.abs(diff) > 3):.0%} of patients")


if __name__ == "__main__":
    main()
