"""Slide scheme (method, slide M4b): every output of the fit, drawn on one simulated fit. Three panels.

A. Per donor: the host state f_hat_i (dot) and its posterior covariance S_i (1 SD ellipse), in the host plane (units of
   the reference anchor, centred over donors). The dashed ellipse is cov(f_hat) + mean S at 1 SD; V is its trace, the
   spread of the host over the donors.
B. Per compartment (here macrophage), centred on its mean mu_k: its states beta_ik (dots); the host pattern (grey, what
   a compartment with z = 1, the reference, would show) and the fitted z_k f_hat (colour): the host pattern turned by
   arg z_k (the phase offset) and scaled by |z_k| (relative amplitude). Both patterns are drawn as curves: donors
   ordered along the host's long axis (no sampling time used), averaged in 7 equal-count bins. The dotted
   circle (1 SD) is sigma_k, its own noise. The arrows are z = 1 (grey) and z_k, at 1.8 host SD times |z|.
C. Coupling, per compartment: the clock variance it shares with the host, |z_k|^2 V (colour), and its own, 2 sigma_k^2
   (grey); c is the shared share. Anchors read on their leave-one-out host.

Drawing convention (as the videos): phase 0 on top, turning clockwise; the package stores (cos, sin). Simulation and
fit as scripts/fit.py (seed 6), point estimates only (n_boot 0). Illustrative, not fitted to data.

Output: paper/presentation/figures/outputs/outputs.{svg,png} (wide slide figure).

    python scripts/outputs.py
"""
import warnings
from pathlib import Path

import numpy as np
from matplotlib.patches import Arc, Circle, Ellipse

import clocksync as cs
import paperfig as pf
from fit import COMPS, HOURS, TARGETS, face, simulate

OUT = Path(__file__).resolve().parents[1] / "figures" / "outputs"
TARGET = "macrophage"
INK, MUTED, GHOST = "#222222", pf.GREY, "#c8c8c8"


def ellipse(ax, centre, cov, n_sd, **kw):
    """Ellipse of a model-space covariance on the clock face (display swaps the components)."""
    w, v = np.linalg.eigh(cov[::-1, ::-1])
    ax.add_patch(Ellipse(centre, 2 * n_sd * np.sqrt(w[1]), 2 * n_sd * np.sqrt(w[0]),
                         angle=np.degrees(np.arctan2(v[1, 1], v[0, 1])), **kw))


def clean(ax, lim):
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    ax.axis("off")


def panel_host(ax, fit):
    f = fit.host[["x", "y"]].to_numpy()
    f = f - f.mean(0)
    S = fit.host_cov.loc[fit.host.index, ["xx", "xy", "yy"]].to_numpy()
    covs = np.array([[[xx, xy], [xy, yy]] for xx, xy, yy in S])
    for c, cov in zip(face(f), covs):
        ellipse(ax, c, cov, 1, fill=False, ec=GHOST, lw=0.6, zorder=1)
    ax.scatter(*face(f).T, s=22, color=INK, lw=0, zorder=3)
    tot = np.cov(f.T) + covs.mean(0)
    ellipse(ax, (0, 0), tot, 1, fill=False, ec=MUTED, lw=1.2, ls="--", zorder=2)
    j = int(np.argmax(face(f)[:, 0]))
    ellipse(ax, face(f)[j], covs[j], 1, fill=True, fc=pf.OKABE_ITO[4], ec=pf.OKABE_ITO[4], alpha=0.35, zorder=2)
    ax.annotate(r"$\hat f_i$, $S_i$" + "\none per donor", face(f)[j], xytext=(14, -18), textcoords="offset points",
                fontsize=11, color=INK, va="top")
    lim = 1.3 * np.abs(face(f)).max()
    ax.annotate(rf"spread $V$ = tr(cov $\hat f$ + mean $S$) = {np.trace(tot):.2f}", (0, -0.92 * lim), ha="center",
                fontsize=11, color=MUTED)
    ax.set_title("per donor: host state", fontsize=13, color=INK)
    clean(ax, lim)
    return lim


def arc(f, n_bins=7):
    """The host pattern as a curve: donors ordered along the host's long axis (no sampling time used), averaged in
    n_bins equal-count bins."""
    order = np.argsort(f @ np.linalg.eigh(np.cov(f.T))[1][:, -1])
    return np.array([f[b].mean(0) for b in np.array_split(order, n_bins)])


def panel_compartment(ax, st, fit, lim):
    t = fit.table.loc[TARGET]
    z = t.rel_amplitude * np.exp(1j * t.offset_h / HOURS)
    host = fit.host[["x", "y"]]
    i = np.flatnonzero((st.compartment == TARGET) & np.isin(st.resection, host.index))
    f = host.loc[st.resection[i]].to_numpy()
    f = f - f.mean(0)
    b = st.beta[i] - st.beta[i].mean(0)
    Mz = cs._core.factor.conformal(z)
    curve = arc(f)
    color = TARGETS[TARGET]
    ax.scatter(*face(b).T, s=14, color=color, alpha=0.3, lw=0, zorder=1)
    ax.plot(*face(curve).T, color=GHOST, lw=5, solid_capstyle="round", zorder=2)
    ax.plot(*face(curve @ Mz.T).T, color=color, lw=3, solid_capstyle="round", zorder=3)
    sd = np.sqrt(t.V / 2)
    L = 1.8 * sd
    for zz, c, w in ((1.0 + 0j, MUTED, 1.4), (z, color, 2.2)):
        ax.annotate("", (L * zz.imag, L * zz.real), (0, 0),
                    arrowprops=dict(arrowstyle="-|>", color=c, lw=w, shrinkA=0, shrinkB=0), zorder=5)
    ax.scatter(0, 0, s=80, marker="+", color=INK, lw=1.8, zorder=6)
    r = 0.6 * L
    ax.add_patch(Arc((0, 0), 2 * r, 2 * r, theta1=90 - np.degrees(np.angle(z)), theta2=90, color=INK, lw=1, zorder=5))
    ax.annotate(f"offset arg $z_k$\n= {t.offset_h:+.1f} h", (0, L), xytext=(0, 6), textcoords="offset points",
                ha="center", va="bottom", fontsize=11, color=INK)
    ax.annotate(rf"$|z_k|$ = {t.rel_amplitude:.2f}", (L * z.imag, L * z.real), xytext=(6, 0),
                textcoords="offset points", fontsize=11, color=color, va="center")
    ends = face(curve @ Mz.T)[[0, -1]]
    end = ends[np.argmax(ends[:, 0])]
    ax.add_patch(Circle(end, np.sqrt(t.psi), fill=False, ec=INK, lw=1, ls=":", zorder=4))
    ax.annotate(r"$\sigma_k$: own noise", end, xytext=(0, -np.sqrt(t.psi) / lim * 160 - 6), textcoords="offset points",
                ha="center", va="top", fontsize=11, color=INK)
    ax.annotate("+: $\\mu_k$; grey: host pattern (z = 1)\ncolour: $z_k \\hat f$, turned and scaled; dots: data", (0, -0.95 * lim),
                ha="center", fontsize=10, color=MUTED)
    ax.set_title(f"per compartment ({TARGET})", fontsize=13, color=INK)
    clean(ax, lim)


def panel_coupling(ax, fit):
    tab = fit.table.loc[list(COMPS)]
    sig, own = tab.z2 * tab.V, 2 * tab.psi
    assert np.allclose(sig / (sig + own), tab.coupling, atol=1e-6)
    y = np.arange(len(tab))
    for yi, (k, s), o in zip(y, sig.items(), own):
        color = TARGETS.get(k, INK)
        ax.barh(yi, s, color=color, height=0.6, lw=0)
        ax.barh(yi, o, left=s, color=GHOST, height=0.6, lw=0)
        ax.annotate(f"c = {tab.loc[k, 'coupling']:.2f}", (s + o, yi), xytext=(5, 0), textcoords="offset points",
                    va="center", fontsize=10, color=INK)
    ax.set_yticks(y, [f"{k} (anchor)" if COMPS[k][1] else k for k in tab.index])
    ax.set_ylim(len(tab) - 0.5, -0.6)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlim(0, 1.3 * float((sig + own).max()))
    ax.set_xlabel(r"clock variance: shared $|z_k|^2 V$ (colour) + own $2\sigma_k^2$ (grey)")
    ax.set_title(r"coupling $c_k = |z_k|^2 V \,/\, (|z_k|^2 V + 2\sigma_k^2)$", fontsize=13, color=INK)


def main():
    pf.use_style("slide")
    OUT.mkdir(parents=True, exist_ok=True)
    ds, _ = simulate()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fit = ds.fit(n_boot=0, reference="fibroblast")
    print(fit.table[["coupling", "rel_amplitude", "offset_h", "z2", "psi", "V"]].round(3).to_string())
    lay = pf.Layout.slide("ABC", kind="wide")
    lim = panel_host(lay.axes["A"], fit)
    panel_compartment(lay.axes["B"], ds.study, fit, lim)
    panel_coupling(lay.axes["C"], fit)
    fig = lay.finalize()
    print(pf.collisions(fig))
    pf.save(fig, OUT / "outputs")


if __name__ == "__main__":
    main()
