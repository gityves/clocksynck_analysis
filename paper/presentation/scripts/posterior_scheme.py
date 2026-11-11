"""Slide scheme (method, slide M4a, notes or backup): why the host posterior S_i is an ellipse under hospital hours.

The model's prior on the host is a Gaussian, f_i ~ N(m, Sigma_f), fitted to the cohort's spread. Each donor's
compartments add a round data term (conformal z_k: precision |z_k|^2 / (sigma_k^2 + v), the same in every direction).
The posterior is N(f_hat_i, S_i), S_i = (Sigma_f^-1 + tau I)^-1, f_hat_i = S_i (Sigma_f^-1 m + tau g_i), with g_i what the
compartments say alone.

A. Hospital hours (deaths 08-18 h): the donors lie on an arc; the Gaussian that fits it is long along the chord and
   thin across. The prior pins the donor across the arc (amplitude) and leaves it free along it (hour): S_i is an
   ellipse along the arc; a guess off across the arc is pulled back strongly, one off along it hardly moves. (The pull
   is towards the cohort's mean, slightly inside the arc: a Gaussian approximates the curve by its chord.)
B. Around the clock (0-24 h): the donors lie on a ring, but a Gaussian cannot say "at distance R": its fit is a round
   disc centred where no donor is. S_i is round, and f_hat_i is pulled towards the centre (lower amplitude). The readout
   puts the lost spread back with V = tr(cov f_hat + mean S).

Drawing convention (as the videos): phase 0 on top, turning clockwise. Schematic: one data precision for all donors;
two example donors per panel, whose compartments' guess g_i is the true host shifted by 0.3 across the arc (first)
or along it (second). Prior, data noise and posterior drawn at 1 SD.

Output: paper/presentation/figures/posterior_scheme/posterior_scheme.{svg,png} (wide slide figure).

    python scripts/posterior_scheme.py
"""
from pathlib import Path

import numpy as np
from matplotlib.patches import Circle, Ellipse

import paperfig as pf

OUT = Path(__file__).resolve().parents[1] / "figures" / "posterior_scheme"
HOURS = 24 / (2 * np.pi)
DATA_VAR = 0.12  # variance of the compartments' combined guess, per coordinate
CASES = {  # panel: (title, hour window, example donors (hour, guess offset radial, tangential, note))
    "A": ("hospital hours (08-18 h)", (8.0, 18.0), [(11.0, 0.3, 0.0, "off across the arc:\npulled back"),
                                                    (15.0, 0.0, 0.3, "off along the arc:\nbarely moves")]),
    "B": ("around the clock (0-24 h)", (0.0, 24.0), [(11.0, 0.3, 0.0, "pulled towards\nthe centre"),
                                                     (15.0, 0.0, 0.3, "pulled towards\nthe centre")]),
}
INK, MUTED, PRIOR, POST = "#222222", pf.GREY, "#bdbdbd", pf.OKABE_ITO[4]


def face(h):
    """Clock-face position of hour h on the unit circle: 0 h on top, clockwise."""
    th = np.asarray(h) / HOURS
    return np.column_stack([np.sin(th), np.cos(th)])


def radial(h):
    th = h / HOURS
    return np.array([np.sin(th), np.cos(th)]), np.array([np.cos(th), -np.sin(th)])


def ellipse(ax, centre, cov, n_sd, **kw):
    w, v = np.linalg.eigh(cov)
    ax.add_patch(Ellipse(centre, 2 * n_sd * np.sqrt(w[1]), 2 * n_sd * np.sqrt(w[0]),
                         angle=np.degrees(np.arctan2(v[1, 1], v[0, 1])), **kw))


def panel(ax, title, window, examples):
    hours = np.linspace(*window, 200, endpoint=window[1] - window[0] < 24)
    f = face(hours)
    m, sig = f.mean(0), np.cov(f.T)
    isig = np.linalg.inv(sig)
    S = np.linalg.inv(isig + np.eye(2) / DATA_VAR)

    ax.add_patch(Circle((0, 0), 1, fill=False, ec="#eeeeee", lw=1, zorder=0))
    for h in (0, 6, 18):
        ax.annotate(f"{h} h", 1.13 * face(h)[0], ha="center", va="center", fontsize=9, color=MUTED)
    ax.plot(*f.T, color=MUTED, lw=4, alpha=0.35, solid_capstyle="round", zorder=1)
    ellipse(ax, m, sig, 1, fill=True, fc=PRIOR, ec="none", alpha=0.25, zorder=0)
    ellipse(ax, m, sig, 1, fill=False, ec=MUTED, lw=1, ls="--", zorder=1)

    for (h, dr, dt, note), side in zip(examples, (1, -1)):
        r, t = radial(h)
        true = face(h)[0]
        g = true + dr * r + dt * t
        fhat = S @ (isig @ m + g / DATA_VAR)
        ax.add_patch(Circle(g, np.sqrt(DATA_VAR), fill=False, ec=INK, lw=0.8, ls=":", zorder=2))
        ellipse(ax, fhat, S, 1, fill=True, fc=POST, ec=POST, alpha=0.35, lw=1, zorder=3)
        ax.scatter(*true, s=30, color=INK, zorder=5)
        ax.scatter(*g, s=40, marker="x", color=INK, lw=1.2, zorder=5)
        ax.annotate("", fhat, g, arrowprops=dict(arrowstyle="-|>", color=POST, lw=1.6, shrinkA=3, shrinkB=2),
                    zorder=4)
        ax.scatter(*fhat, s=30, color=POST, edgecolor="white", lw=0.6, zorder=6)
        ax.annotate(note, g, xytext=(side * 52, 0), textcoords="offset points",
                    ha="left" if side > 0 else "right", va="center", fontsize=10, color=INK, zorder=7)

    ax.set_title(title, fontsize=13, color=INK)
    ax.set_xlim(-1.3, 1.3)
    ax.set_ylim(-1.6, 1.2)
    ax.set_aspect("equal")
    ax.axis("off")
    return sig, S


def legend(ax):
    """Key under the panels, text only beside its glyph."""
    items = [("o", INK, "true host (unknown)"), ("x", INK, r"$g_i$: the donor's compartments alone (dotted: their noise)"),
             ("o", POST, r"$\hat f_i$, $S_i$: after the prior"), ("s", PRIOR, r"prior $\Sigma_f$ (dashed): the cohort's spread")]
    for j, (mk, c, txt) in enumerate(items):
        x, y = 0.02 + 0.5 * (j % 2), 0.7 - 0.5 * (j // 2)
        ax.scatter([x], [y], s=40, marker=mk, color=c, transform=ax.transAxes, clip_on=False, lw=1.2)
        ax.annotate(txt, (x, y), xycoords="axes fraction", xytext=(10, 0), textcoords="offset points", va="center",
                    fontsize=10, color=INK)
    ax.axis("off")


def main():
    pf.use_style("slide")
    OUT.mkdir(parents=True, exist_ok=True)
    lay = pf.Layout.slide("AB;LL", kind="wide", row_heights=[1, 0.16])
    for key, (title, window, ex) in CASES.items():
        sig, S = panel(lay.axes[key], title, window, ex)
        print(f"{key}: Sigma_f eigenvalues {np.linalg.eigvalsh(sig).round(3)}, S eigenvalues {np.linalg.eigvalsh(S).round(3)}")
    legend(lay.axes["L"])
    fig = lay.finalize()
    print(pf.collisions(fig))
    pf.save(fig, OUT / "posterior_scheme")


if __name__ == "__main__":
    main()
