"""Slide figure (method, slide M1): clock-gene pseudobulks projected onto the clock plane.

A. The reference: the 14 genes of the GTEx waveform (clocksync.gtex()), angle = acrophase, length = relative amplitude.
   Fixed, never re-estimated. Its zero is a reference phase, not clock time.
B. One pseudobulk (one donor's fibroblasts): its gene-centred log fractions, minus the compartment's mean over donors
   (bars, genes in acrophase order), and the fitted plane X_g (beta - mean beta) (dots). The gene mesors are left out:
   they only add a constant to beta, absorbed by the compartment mean mu_k.
C. The clock plane: beta of every donor of three compartments, with 2 SD Poisson ellipses (V). Each cloud is an arc
   (08-18 h of its cycle) around its compartment's baseline; a constant per-gene offset b (platform, cell-type
   baseline) moves the whole cloud by L b (hepatocytes: ghost -> solid). Baselines are set for readability (BASE).

Simulated with clocksync.simulate (30 donors sampled 08-18 h, loadings as in model_animation.py: endothelial 0.75 at
+3 h, hepatocyte 1.15 at -4 h); beta and V from cs.Dataset (the package projection). Panel C uses the clock-face
convention of the video: phase 0 on top, turning clockwise. Illustrative, not fitted to data.

Output: paper/presentation/figures/projection/projection.{svg,png} (320 x 110 mm, "wide").

    python scripts/projection.py
"""
from pathlib import Path

import numpy as np
from matplotlib.patches import Ellipse, FancyArrowPatch

import clocksync as cs
import paperfig as pf

OUT = Path(__file__).resolve().parents[1] / "figures" / "projection"
N, WINDOW, DEPTH, SEED = 30, (8.0, 18.0), 1e7, 7
BASE = {  # panel C: where each compartment's cloud sits (display coordinates, plane units)
    "fibroblast": np.array([0.0, 0.0]),
    "endothelial": np.array([-2.7, 2.0]),
    "hepatocyte": np.array([2.9, 0.4]),
    "hepatocyte (before)": np.array([1.2, 2.9]),
}
HOURS = 24 / (2 * np.pi)
COMPS = {  # name: (z, marker), as in model_animation.py
    "fibroblast": (1.0 + 0j, "o"),
    "endothelial": (0.75 * np.exp(1j * 3 / HOURS), "^"),
    "hepatocyte": (1.15 * np.exp(-1j * 4 / HOURS), "s"),
}
COLOR = dict(zip(COMPS, pf.OKABE_ITO))
DONOR = 4  # the fibroblast pseudobulk shown in panel B
INK, MUTED = "#222222", pf.GREY


def face(b):
    """Clock-face display of plane states (x, y): phase 0 on top, turning clockwise (as model_animation.D)."""
    b = np.atleast_2d(b)
    return np.column_stack([b[:, 1], b[:, 0]])


def simulate():
    wf = cs.gtex()
    comps = {k: cs.Compartment(z, psi=0.004, mesor_sd=0.0) for k, (z, _) in COMPS.items()}
    df, _ = cs.simulate(wf, comps, n=N, window_h=WINDOW, jitter_h=0.5, depth=DEPTH, tau=0.02, shift_sd=0.2, seed=SEED)
    ds = cs.Dataset(df.rename(columns={"resection": "donor"}), anchors=["fibroblast"], min_donors=5)
    return wf, ds.study


def heading(ax, title, sub, dy=(30, 10)):
    for text, d, size, color in ((title, dy[0], 13, INK), (sub, dy[1], 11, MUTED)):
        ax.annotate(text, (0.5, 1), xytext=(0, d), xycoords="axes fraction", textcoords="offset points", ha="center",
                    va="bottom", fontsize=size, color=color, annotation_clip=False)


def panel_reference(ax, wf):
    psi, amp = np.mod(wf.psi, 2 * np.pi), np.asarray(wf.amp)
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    rmax = amp.max() * 1.12
    ax.set_ylim(0, rmax)
    ax.set_yticks([])
    ax.set_xticks(np.radians([0, 90, 180, 270]), ["0", "6", "12", "18"], color=MUTED)
    ax.grid(False)
    ax.spines["polar"].set_color(MUTED)
    ax.spines["polar"].set_linewidth(0.8)
    for p, a in zip(psi, amp):
        ax.annotate("", xy=(p, a), xytext=(p, 0), arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.3, mutation_scale=9))
    tips = np.column_stack([psi, amp])
    segs = np.array([[[p, 0], [p, a]] for p, a in zip(psi, amp)])
    pf.place_labels(ax, psi, amp, list(wf.genes), segments=segs, obstacles=tips, prefer=list(np.pi / 2 - psi),
                    fontsize=9, color=INK, inside=False)
    heading(ax, "Reference waveform (GTEx)", "acrophase and relative amplitude of 14 genes: fixed", dy=(48, 28))


def panel_profile(ax, wf, st):
    i = np.flatnonzero((st.compartment == "fibroblast"))
    y = np.log(np.maximum(st.counts[i], 0.5) / st.depth[i, None])
    y = y - y.mean(1, keepdims=True)  # centred over genes
    dev = y[DONOR] - y.mean(0)  # minus the compartment mean over donors
    order = np.argsort(np.mod(wf.psi, 2 * np.pi))
    Xc = wf.X - wf.X.mean(0)
    fitted = Xc @ (st.beta[i][DONOR] - st.beta[i].mean(0))
    fitted = fitted - fitted.mean()
    x = np.arange(len(order))
    ax.bar(x, dev[order], width=0.7, color="#c9c9c9", zorder=2)
    ax.scatter(x, fitted[order], s=40, color=COLOR["fibroblast"], edgecolor="white", linewidth=0.8, zorder=3)
    ax.axhline(0, color=MUTED, lw=0.8, zorder=1)
    ax.set_xticks(x, [wf.genes[g] for g in order], rotation=90, fontsize=9)
    ax.set_xlim(-0.7, len(x) - 0.3)
    ax.set_ylabel("log fraction vs compartment mean")
    b = st.beta[i][DONOR] - st.beta[i].mean(0)
    ax.text(0.03, 0.05, f"β − mean: λ {np.hypot(*b):.2f}, φ {np.mod(np.arctan2(b[1], b[0]), 2 * np.pi) * HOURS:.1f} h",
            transform=ax.transAxes, ha="left", va="bottom", fontsize=10, color=INK)
    heading(ax, "One pseudobulk → one point β", "bars: data; dots: weighted least-squares fit on the plane", dy=(48, 28))


def ellipse(ax, centre, cov, color, alpha=1.0):
    c = cov[::-1, ::-1]  # clock-face display swaps the axes
    w, v = np.linalg.eigh(c)
    ang = np.degrees(np.arctan2(v[1, 1], v[0, 1]))
    ax.add_patch(Ellipse(centre, 4 * np.sqrt(w[1]), 4 * np.sqrt(w[0]), angle=ang, fill=False, ec=color, lw=0.5,
                         alpha=alpha, zorder=2))


def panel_plane(ax, wf, st):
    """Each compartment's cloud is placed at its baseline BASE (display coordinates): its gene-mesor term L_k m is
    replaced by L_k b_k for a per-gene baseline b_k, exact because the projection is linear."""
    pts = []
    for k, (_, marker) in COMPS.items():
        i = st.compartment == k
        b = st.beta[i] - st.L[k] @ wf.mesor  # the arc z_k e^(i t_s) + noise, around 0
        if k == "hepatocyte":
            ghost = face(b) + BASE["hepatocyte (before)"]
            ax.scatter(*ghost.T, s=34, marker=marker, facecolor="none", edgecolor=COLOR[k], linewidth=0.7, alpha=0.5,
                       zorder=2)
            pts.append(ghost)
            a = ghost.mean(0)
        p = face(b) + BASE[k]
        pts.append(p)
        for c, v in zip(p, st.V[i]):
            ellipse(ax, c, v, COLOR[k], alpha=0.45)
        ax.scatter(*p.T, s=40, marker=marker, color=COLOR[k], edgecolor="white", linewidth=0.7, zorder=3, label=k)
    z = a + BASE["hepatocyte"] - BASE["hepatocyte (before)"]  # cloud mean before -> after
    ax.add_patch(FancyArrowPatch(a, z, arrowstyle="-|>", mutation_scale=14, color=INK, lw=1.4, zorder=4, shrinkA=4,
                                 shrinkB=4))
    pf.place_labels(ax, [(a[0] + z[0]) / 2], [(a[1] + z[1]) / 2], ["constant offset\n(platform, baseline)"],
                    obstacles=np.vstack(pts), segments=np.array([[a, z]]), fontsize=10, color=INK)
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xticks([])
    ax.set_yticks([])
    for sp in ("left", "bottom"):
        ax.spines[sp].set_visible(False)
    ax.legend(loc="lower left", fontsize=10, handletextpad=0.3, borderaxespad=0.2)
    heading(ax, "Clock plane: one point per donor × compartment", "ellipses: 2 SD Poisson error; offsets move a cloud",
            dy=(48, 28))


def main():
    pf.use_style("slide")
    wf, st = simulate()
    lay = pf.Layout.slide("ABC", "wide", polar="A", col_ratios=[1, 1.15, 1.15], wgap=18)
    panel_reference(lay.axes["A"], wf)
    panel_profile(lay.axes["B"], wf, st)
    panel_plane(lay.axes["C"], wf, st)
    fig = lay.finalize()
    OUT.mkdir(parents=True, exist_ok=True)
    pf.save(fig, OUT / "projection")
    failed = pf.place_all(fig)
    if failed:
        print("labels not placed:", list(failed))


if __name__ == "__main__":
    main()
