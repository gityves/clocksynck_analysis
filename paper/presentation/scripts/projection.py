"""Slide figures (method, slide M1): clock-gene pseudobulks projected onto the clock plane. Three one-panel figures,
shown in order beside the bullets.

A. The reference: the 14 genes of the GTEx waveform (clocksync.gtex()), one point per gene at angle = acrophase psi_g,
   radius = relative amplitude amp_g. Fixed, never re-estimated. Its zero is a reference phase, not clock time.
B. One pseudobulk (one donor's fibroblasts) against its compartment: each gene's log fraction minus the compartment's
   mean over donors (computable on real data), divided by the gene's fixed reference amplitude amp_g, against its
   acrophase psi_g. The model says these fall on one cosine, lam cos(phi - psi); the curve is the package's beta minus
   the compartment's mean beta (exact by linearity: L (y - mean y) = beta - mean beta), not refitted here. Its height
   lam and peak phi are the arrow from the cloud's mean to the ringed donor in C. Point area is the gene's weight in
   the fit (w_g amp_g^2 in this scaling): low-amplitude genes scatter most and count least. Gene centring (which removes
   the per-sample common shift s_ik) is not drawn: s_ik = 0 in this simulation (shift_sd 0).
C. The clock plane: beta of every donor of three compartments, 2 SD Poisson ellipses (V). Each cloud is an arc (08-18 h
   of its cycle) around its compartment's baseline; a constant per-gene offset b (platform, cell-type baseline) moves
   the whole cloud by L b (hepatocytes: ghost -> solid). Baselines are set for readability (BASE): the gene-mesor term
   L_k m is replaced by a chosen constant, exact because the projection is linear. The ringed fibroblast is panel B's.

Drawing convention (as the video): phase 0 on top, turning clockwise. With beta = lam (sin phi, cos phi) and
X_g = amp_g (sin psi_g, cos psi_g), X_g . beta = lam amp_g cos(phi - psi_g): beta_1 is horizontal, beta_2 vertical. The
package stores (cos, sin); the order of the two components is a drawing convention.

Simulated with clocksync.simulate (30 donors sampled 08-18 h, loadings as in model_animation.py: endothelial 0.75 at
+3 h, hepatocyte 1.15 at -4 h); beta and V from cs.Dataset (the package projection). Illustrative, not fitted to data.

Output: paper/presentation/figures/projection/projection_{A,B,C}.{svg,png} (150 x 120 mm each).

    python scripts/projection.py
"""
from pathlib import Path

import numpy as np
from matplotlib.patches import Ellipse, FancyArrowPatch

import clocksync as cs
import paperfig as pf

OUT = Path(__file__).resolve().parents[1] / "figures" / "projection"
SIZE = (150.0, 120.0)  # mm, one figure beside the bullets
N, WINDOW, DEPTH, SEED = 30, (8.0, 18.0), 1e7, 7
HOURS = 24 / (2 * np.pi)
COMPS = {  # name: (z, marker), as in model_animation.py
    "fibroblast": (1.0 + 0j, "o"),
    "endothelial": (0.75 * np.exp(1j * 3 / HOURS), "^"),
    "hepatocyte": (1.15 * np.exp(-1j * 4 / HOURS), "s"),
}
BASE = {  # panel C: where each compartment's cloud sits (display coordinates, plane units)
    "fibroblast": np.array([0.0, 0.0]),
    "endothelial": np.array([-2.7, 2.0]),
    "hepatocyte": np.array([2.9, 0.4]),
    "hepatocyte (before)": np.array([1.2, 2.9]),
}
COLOR = dict(zip(COMPS, pf.OKABE_ITO))
INK, MUTED = "#222222", pf.GREY


def donor(st):
    """Panel B's fibroblast pseudobulk (ringed in C): the one farthest from its cloud's mean, so the arrow reads."""
    b = st.beta[st.compartment == "fibroblast"]
    return int(np.argmax(np.hypot(*(b - b.mean(0)).T)))


def face(b):
    """Clock-face display of plane states stored as (cos, sin): (sin, cos) = (beta_1, beta_2), phase 0 on top."""
    b = np.atleast_2d(b)
    return np.column_stack([b[:, 1], b[:, 0]])


def simulate():
    wf = cs.gtex()
    comps = {k: cs.Compartment(z, psi=0.004, mesor_sd=0.0) for k, (z, _) in COMPS.items()}
    df, _ = cs.simulate(wf, comps, n=N, window_h=WINDOW, jitter_h=0.5, depth=DEPTH, tau=0.02, shift_sd=0.0, seed=SEED)
    ds = cs.Dataset(df.rename(columns={"resection": "donor"}), anchors=["fibroblast"], min_donors=5)
    return wf, ds.study


def title(ax, text, dy=14):
    ax.annotate(text, (0.5, 1), xytext=(0, dy), xycoords="axes fraction", textcoords="offset points", ha="center",
                va="bottom", fontsize=13, color=INK, annotation_clip=False)


def panel_reference(ax, wf):
    psi, amp = np.mod(wf.psi, 2 * np.pi), np.asarray(wf.amp)
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    ax.set_ylim(0, amp.max() * 1.15)
    ax.set_yticks([0.25, 0.5], ["0.25", "0.5"], color=MUTED, fontsize=9)
    ax.set_rlabel_position(225)
    ax.set_xticks(np.radians([0, 90, 180, 270]), ["0 h", "6", "12", "18"], color=MUTED)
    ax.grid(True, color="#e6e6e6", lw=0.7)
    ax.spines["polar"].set_color(MUTED)
    ax.spines["polar"].set_linewidth(0.8)
    ax.scatter(psi, amp, s=60, color=INK, edgecolor="white", linewidth=0.8, zorder=3)
    pf.place_labels(ax, psi, amp, list(wf.genes), prefer=list(np.pi / 2 - psi), fontsize=10, color=INK, inside=False,
                    clearance=1.2)
    title(ax, "Reference waveform (GTEx): fixed", dy=30)


def panel_profile(ax, wf, st):
    i, d = np.flatnonzero(st.compartment == "fibroblast"), donor(st)
    logf = np.log(np.maximum(st.counts[i], 0.5) / st.depth[i][:, None])
    y = logf[d] - logf.mean(0)  # minus the compartment's mean over donors
    amp, psi = np.asarray(wf.amp), np.mod(wf.psi, 2 * np.pi)
    b = st.beta[i][d] - st.beta[i].mean(0)  # the package's beta minus the cloud's mean
    lam, phi = np.hypot(*b), np.arctan2(b[1], b[0])
    w = 1.0 / (1.0 / np.maximum(np.median(st.counts[i], 0), 0.5) + st.tau2)  # the package's gene weights
    grid = np.linspace(0, 2 * np.pi, 300)
    ax.plot(grid * HOURS, lam * np.cos(phi - grid), color=COLOR["fibroblast"], lw=2.2, zorder=2)
    size = 220 * (w * amp ** 2) / (w * amp ** 2).max()
    ax.scatter(psi * HOURS, y / amp, s=np.maximum(size, 12), color="#9a9a9a", edgecolor="white", linewidth=0.8, zorder=3)
    ax.axhline(0, color="#dddddd", lw=0.8, zorder=0)
    peak = np.mod(phi, 2 * np.pi) * HOURS
    ax.annotate("", xy=(peak, lam), xytext=(peak, 0), arrowprops=dict(arrowstyle="<->", color=INK, lw=1.2,
                                                                       mutation_scale=10, shrinkA=0, shrinkB=0))
    ax.annotate("λ", (peak, lam / 2), xytext=(-6, 0), textcoords="offset points", ha="right", va="center", fontsize=12)
    ax.annotate("φ", (peak, 0), xytext=(0, -6), textcoords="offset points", ha="center", va="top", fontsize=12)
    ax.set_xlim(0, 24)
    ax.set_xticks([0, 6, 12, 18, 24])
    lim = 1.6 * lam
    ax.set_ylim(-lim, lim)
    ax.set_xlabel("gene acrophase ψ (h)")
    ax.set_ylabel("(log fraction − compartment mean)\n/ reference amplitude")
    pf.place_labels(ax, psi * HOURS, np.clip(y / amp, -lim, lim), list(wf.genes), fontsize=9, color=MUTED,
                    clearance=1.5)
    title(ax, "One pseudobulk → one β (vs its compartment)")


def ellipse(ax, centre, cov, color, alpha=1.0):
    c = cov[::-1, ::-1]  # the display swaps the stored components
    w, v = np.linalg.eigh(c)
    ang = np.degrees(np.arctan2(v[1, 1], v[0, 1]))
    ax.add_patch(Ellipse(centre, 4 * np.sqrt(w[1]), 4 * np.sqrt(w[0]), angle=ang, fill=False, ec=color, lw=0.5,
                         alpha=alpha, zorder=2))


def panel_plane(ax, wf, st):
    pts, d = [], donor(st)
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
        if k == "fibroblast":  # panel B's donor, and the arrow from the cloud's mean that B's cosine encodes
            m = p.mean(0)
            ax.scatter(*m, s=90, marker="+", color=INK, linewidth=1.6, zorder=4)
            ax.add_patch(FancyArrowPatch(m, p[d], arrowstyle="-|>", mutation_scale=11, color=INK, lw=1.1, zorder=4,
                                         shrinkA=2, shrinkB=7))
            ax.scatter(*p[d], s=170, facecolor="none", edgecolor=INK, linewidth=1.3, zorder=4)
            mean_at, mean_seg = m, [m, p[d]]
    z = a + BASE["hepatocyte"] - BASE["hepatocyte (before)"]  # cloud mean before -> after
    ax.add_patch(FancyArrowPatch(a, z, arrowstyle="-|>", mutation_scale=14, color=INK, lw=1.4, zorder=4, shrinkA=4,
                                 shrinkB=4))
    pf.place_labels(ax, [(a[0] + z[0]) / 2], [(a[1] + z[1]) / 2], ["constant offset\n(platform, baseline)"],
                    obstacles=np.vstack(pts), segments=np.array([[a, z]]), fontsize=10, color=INK)
    pf.place_labels(ax, [mean_at[0]], [mean_at[1]], ["mean"], obstacles=np.vstack(pts),
                    segments=np.array([mean_seg]), prefer=[np.pi], fontsize=10, color=INK)
    ax.set_aspect("equal", adjustable="datalim")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_xlabel(r"$\beta_1$", fontsize=13)
    ax.set_ylabel(r"$\beta_2$", fontsize=13, rotation=0, labelpad=12)
    for sp in ("left", "bottom"):
        ax.spines[sp].set_color(MUTED)
    ax.legend(loc="lower left", fontsize=10, handletextpad=0.3, borderaxespad=0.2)
    title(ax, "Clock plane: one point per donor × compartment")


def main():
    pf.use_style("slide")
    wf, st = simulate()
    OUT.mkdir(parents=True, exist_ok=True)
    for key, draw, polar in (("A", panel_reference, "A"), ("B", panel_profile, ""), ("C", panel_plane, "")):
        lay = pf.Layout.slide("A", width=SIZE[0], height=SIZE[1], polar=polar)
        draw(lay.axes["A"], wf, st) if key != "A" else draw(lay.axes["A"], wf)
        fig = lay.finalize()
        pf.save(fig, OUT / f"projection_{key}")
        failed = pf.place_all(fig)
        if failed:
            print(f"{key}: labels not placed:", list(failed))
    i = st.compartment == "fibroblast"
    b = st.beta[i][donor(st)] - st.beta[i].mean(0)
    print(f"panel B donor: lam {np.hypot(*b):.2f}, phi {np.mod(np.arctan2(b[1], b[0]), 2 * np.pi) * HOURS:.1f} h")


if __name__ == "__main__":
    main()
