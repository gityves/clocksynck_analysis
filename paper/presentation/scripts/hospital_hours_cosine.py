"""Slide schematic (introduction, slide 4.4): a hospital-hours fragment of the clock is fitted by a full, smaller cosine.

Patients are sampled once between 10 and 18 h, so each clock gene shows only an 8 h fragment of its 24 h cosine
(panel A, true time). An orderer that assumes the samples cover the cycle spreads them over a full period: here the
angle in the plane of the first two eigengenes (each scaled to unit SD, as in CYCLOPS), mapped by rank to a uniform
phase, with the origin anchored on the ARNTL acrophase. Against these assigned phases a full cosine with a smaller
amplitude fits the same samples about as well as the true fragment does (panel B), so a goodness of fit cannot tell a
full cycle from a stretched fragment.

Illustrative simulation, not fitted to data: unit amplitudes, Gaussian noise, acrophases ARNTL 20 h, NR1D1 8 h
(ARNTL antiphasic to NR1D1, ruben2018database) and PER2 14 h (after NR1D1, balsalobre1998serum); absolute times are
illustrative.

Output: paper/presentation/figures/hospital_hours_cosine/hospital_hours_cosine.{svg,png} (320 x 110 mm, "wide").

    python scripts/hospital_hours_cosine.py
"""
from pathlib import Path

import numpy as np

import paperfig as pf

OUT = Path(__file__).resolve().parents[1] / "figures" / "hospital_hours_cosine"
N, WINDOW, NOISE, SEED = 60, (10.0, 18.0), 0.15, 1
GENES = {"ARNTL": 20.0, "NR1D1": 8.0, "PER2": 14.0}  # acrophase (h)
COLOR = dict(zip(GENES, pf.OKABE_ITO))
INK, MUTED, BAND = "#222222", pf.GREY, "#ececec"
W = 2 * np.pi / 24


def simulate():
    rng = np.random.default_rng(SEED)
    t = rng.uniform(*WINDOW, N)
    y = np.column_stack([np.cos(W * (t - p)) for p in GENES.values()]) + rng.normal(0, NOISE, (N, len(GENES)))
    return t, y


def cosinor(h, y):
    """Least-squares y = m + a cos(W (h - p)); returns m, a, p (h) and R^2."""
    X = np.column_stack([np.ones_like(h), np.cos(W * h), np.sin(W * h)])
    b = np.linalg.lstsq(X, y, rcond=None)[0]
    return b[0], np.hypot(b[1], b[2]), (np.arctan2(b[2], b[1]) / W) % 24, 1 - np.var(y - X @ b) / np.var(y)


def order(y):
    """Assigned phase (h): eigengene angle -> rank -> uniform over 24 h; origin and direction set by the gene order."""
    u = np.linalg.svd(y - y.mean(0), full_matrices=False)[0][:, :2]
    u = u / u.std(0)
    rank = np.argsort(np.argsort(np.angle(u[:, 0] + 1j * u[:, 1])))
    h = 24 * (rank + 0.5) / len(y)
    best = None
    for s in (1, -1):  # direction: the one that keeps the true acrophase order
        hs = (s * h) % 24
        p = np.array([cosinor(hs, y[:, g])[2] for g in range(y.shape[1])])
        shift = GENES["ARNTL"] - p[0]  # anchor: ARNTL at its true acrophase
        err = np.abs((p + shift - np.array(list(GENES.values())) + 12) % 24 - 12).sum()
        if best is None or err < best[0]:
            best = (err, (hs + shift) % 24)
    return best[1]


def heading(ax, title, sub):
    for text, dy, size, color in ((title, 30, 13, INK), (sub, 10, 11, MUTED)):
        ax.annotate(text, (0.5, 1), xytext=(0, dy), xycoords="axes fraction", textcoords="offset points", ha="center",
                    va="bottom", fontsize=size, color=color, annotation_clip=False)


def main():
    pf.use_style("slide")
    t, y = simulate()
    h = order(y)
    true = {g: (*cosinor(t, y[:, i])[:2], GENES[g], cosinor(t, y[:, i])[3]) for i, g in enumerate(GENES)}
    fit = {g: cosinor(h, y[:, i]) for i, g in enumerate(GENES)}

    lay = pf.Layout.slide("AB", "wide", wgap=16)
    A, B = lay.axes["A"], lay.axes["B"]
    A.axvspan(*WINDOW, color=BAND, zorder=0, lw=0)
    A.text(np.mean(WINDOW), 1.6, "hospital hours", ha="center", va="center", fontsize=11, color=MUTED)
    grid = np.linspace(0, 24, 400)
    for i, g in enumerate(GENES):
        A.plot(grid, np.cos(W * (grid - GENES[g])), color=COLOR[g], lw=2, zorder=3)
        A.scatter(t, y[:, i], s=22, color=COLOR[g], edgecolor="white", linewidth=0.6, alpha=0.85, zorder=4)
    for i, g in enumerate(GENES):
        B.plot(grid, fit[g][0] + fit[g][1] * np.cos(W * (grid - fit[g][2])), color=COLOR[g], lw=2, zorder=3)
        B.scatter(h, y[:, i], s=22, color=COLOR[g], edgecolor="white", linewidth=0.6, alpha=0.85, zorder=4)
    for ax, xlabel in ((A, "time of day (h)"), (B, "phase assigned by the orderer (h)")):
        ax.set_xlim(0, 24)
        ax.set_xticks([0, 6, 12, 18, 24])
        ax.set_ylim(-1.45, 1.8)
        ax.set_yticks([-1, 0, 1])
        ax.set_xlabel(xlabel)
    A.set_ylabel("expression (a.u.)")
    lay.share("AB", y=True)

    handles = [A.plot([], [], color=COLOR[g], lw=2, marker="o", ms=6, mec="white")[0] for g in GENES]
    A.legend(handles, list(GENES), loc="upper left", bbox_to_anchor=(0.0, 1.0), ncol=1, fontsize=11, handlelength=1.6)
    key = lambda d: "   ".join(f"{g} {d[g][3]:.2f}" for g in GENES)
    A.text(0.5, -0.27, f"R² of the true cosine on the samples:  {key(true)}", transform=A.transAxes, ha="center",
           va="top", fontsize=10, color=MUTED)
    B.text(0.5, -0.27, f"R² of the fitted full cosine:  {key(fit)}", transform=B.transAxes, ha="center", va="top",
           fontsize=10, color=MUTED)

    heading(A, "True time: 8 h of each cosine", "samples taken between 10 and 18 h")
    heading(B, "Orderer: the same samples over a full cycle", "a full cosine with a smaller amplitude also fits")

    fig = lay.finalize()
    OUT.mkdir(parents=True, exist_ok=True)
    pf.save(fig, OUT / "hospital_hours_cosine")
    for g in GENES:
        print(f"{g}: true R2 {true[g][3]:.2f} amp {true[g][1]:.2f} | stretched R2 {fit[g][3]:.2f} amp {fit[g][1]:.2f} "
              f"acrophase {fit[g][2]:.1f} h (true {GENES[g]:.0f} h)")


if __name__ == "__main__":
    main()
