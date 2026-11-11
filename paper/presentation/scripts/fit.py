"""Slide figures (method, slide M4): the fit, in two steps. Three one-panel figures, shown in order beside the bullets.

A. Step 1, the host: EM on the anchor compartments only (fibroblast, endothelial, smooth muscle) gives each donor's
   host state f_i: posterior mean (point) and 2 SD posterior ellipse, in the clock plane (fibroblast units, centred
   over donors). Points are coloured by the true sampling hour, which the method never sees: the host orders the
   donors along the arc of their hours.
B. Step 2, the readout: each compartment regressed on the host, beta_ik = mu_k + z_k f_i + e_ik. Observed state
   (beta_ik minus the compartment's mean) against the fitted z_k f_i (minus its mean), both plane coordinates pooled.
   A compartment that follows the host lies on the diagonal; its scatter around it is its unshared variance sigma_k^2 plus
   measurement noise. A compartment with its own, unshared clock has z near 0: all its variation is off the line.
C. The outputs: coupling c with its 95 % donor-bootstrap interval (host refitted in each replicate) against the true
   coupling of the simulation, z2 V / (z2 V + 2 psi_total) with V the spread of the true sampling phases and psi_total
   the compartment's noise plus its gene-level noise projected on the plane.

Drawing convention (as the video): phase 0 on top, turning clockwise; the package stores (cos, sin).

Simulation: 40 donors x 6 compartments, deaths 08-18 h (jitter 1 h), depth 1e7; anchors fibroblast (z 1), endothelial
(0.75 at +3 h), smooth muscle (0.9 at +1 h); targets hepatocyte (1.15 at -4 h), macrophage (0.6 at +2 h, higher noise),
lymphocyte (own clock, not shared). Seed 6 is typical: over seeds 0-11 the macrophage coupling has median 0.646 against
a true median 0.632 (n_boot 0); the anchors read low, as expected with three anchors (A154, A162).
`Dataset.fit(n_boot=500, reference="fibroblast")`. Illustrative, not fitted to data.

Output: paper/presentation/figures/fit/fit_{A,B,C}.{svg,png} (150 x 120 mm each).

    python scripts/fit.py
"""
import warnings
from pathlib import Path

import numpy as np
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Ellipse

import clocksync as cs
import paperfig as pf

OUT = Path(__file__).resolve().parents[1] / "figures" / "fit"
SIZE = (150.0, 120.0)  # mm, one figure beside the bullets
N, WINDOW, DEPTH, TAU, SEED, N_BOOT = 40, (8.0, 18.0), 1e7, 0.02, 6, 500
HOURS = 24 / (2 * np.pi)
COMPS = {  # name: (simulator compartment, anchor)
    "fibroblast": (cs.Compartment(1.0, psi=0.01, mesor_sd=0.3), True),
    "endothelial": (cs.Compartment(0.75 * np.exp(1j * 3 / HOURS), psi=0.01, mesor_sd=0.3), True),
    "smooth muscle": (cs.Compartment(0.9 * np.exp(1j * 1 / HOURS), psi=0.01, mesor_sd=0.3), True),
    "hepatocyte": (cs.Compartment(1.15 * np.exp(-1j * 4 / HOURS), psi=0.02, mesor_sd=0.3), False),
    "macrophage": (cs.Compartment(0.6 * np.exp(1j * 2 / HOURS), psi=0.05, mesor_sd=0.3), False),
    "lymphocyte": (cs.Compartment(0.8, psi=0.03, coupled=False, mesor_sd=0.3), False),
}
ANCHORS = [k for k, (_, a) in COMPS.items() if a]
TARGETS = {"hepatocyte": pf.OKABE_ITO[2], "macrophage": pf.OKABE_ITO[0], "lymphocyte": pf.OKABE_ITO[5]}
INK, MUTED = "#222222", pf.GREY
HOUR_CMAP = LinearSegmentedColormap.from_list("hour", ["#c6dbef", "#08306b"])


def simulate():
    df, t = cs.simulate(cs.gtex(), {k: c for k, (c, _) in COMPS.items()}, n=N, window_h=WINDOW, depth=DEPTH, tau=TAU,
                        shift_sd=0.1, seed=SEED)
    ds = cs.Dataset(df.rename(columns={"resection": "donor"}), anchors=ANCHORS, min_donors=5)
    hours = dict(zip([f"r{i:03d}" for i in range(N)], t))
    return ds, hours


def true_coupling(st, hours):
    th = np.array(list(hours.values())) / HOURS
    V = np.trace(np.cov(np.c_[np.cos(th), np.sin(th)].T))  # spread of the true host (fibroblast units)
    out = {}
    for k, (c, _) in COMPS.items():
        psi = c.psi + TAU ** 2 * np.trace(st.L[k] @ st.L[k].T) / 2  # unshared variance: compartment + gene noise in the plane
        z2 = abs(c.z) ** 2 if c.coupled else 0.0
        out[k] = z2 * V / (z2 * V + 2 * psi)
    return out


def face(b):
    """Clock-face display of plane states stored as (cos, sin): (sin, cos), phase 0 on top."""
    b = np.atleast_2d(b)
    return np.column_stack([b[:, 1], b[:, 0]])


def panel_host(ax, fit, hours):
    h, hc = fit.host, fit.host_cov.loc[fit.host.index]
    f = h[["x", "y"]].to_numpy()
    f = f - f.mean(0)
    hr = np.array([hours[d] for d in h.index])
    norm = Normalize(7, 19)
    for (x, y), (xx, xy, yy), t in zip(face(f), hc[["xx", "xy", "yy"]].to_numpy(), hr):
        w, v = np.linalg.eigh(np.array([[yy, xy], [xy, xx]]))  # display swaps the components
        ax.add_patch(Ellipse((x, y), 4 * np.sqrt(w[1]), 4 * np.sqrt(w[0]),
                             angle=np.degrees(np.arctan2(v[1, 1], v[0, 1])), fill=False, lw=0.5,
                             ec=HOUR_CMAP(norm(t)), alpha=0.6, zorder=1))
    sc = ax.scatter(*face(f).T, c=hr, cmap=HOUR_CMAP, norm=norm, s=40, edgecolor="white", linewidth=0.6, zorder=3)
    cb = ax.figure.colorbar(sc, cax=ax.inset_axes([1.04, 0.3, 0.035, 0.4]), ticks=[8, 13, 18])
    cb.outline.set_visible(False)
    cb.set_label("true sampling hour\n(not used)", color=MUTED)
    lim = 1.25 * np.abs(f).max()
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    ax.axhline(0, color="#eeeeee", lw=0.8, zorder=0)
    ax.axvline(0, color="#eeeeee", lw=0.8, zorder=0)
    ax.set_xlabel(r"host $f_1$")
    ax.set_ylabel(r"host $f_2$")


def readout(st, fit, k):
    """Observed (beta minus its mean) and fitted (z f minus its mean) plane states of compartment k."""
    t = fit.table.loc[k]
    z = t.rel_amplitude * np.exp(1j * t.offset_h / HOURS)
    host = fit.host[["x", "y"]]
    i = np.flatnonzero((st.compartment == k) & np.isin(st.resection, host.index))
    f = host.loc[st.resection[i]].to_numpy()
    pred = f @ cs._core.factor.conformal(z).T
    obs = st.beta[i]
    return (obs - obs.mean(0)).ravel(), (pred - pred.mean(0)).ravel()


def panel_readout(ax, st, fit):
    lim = 0
    for k, color in TARGETS.items():
        obs, pred = readout(st, fit, k)
        ax.scatter(pred, obs, s=16, color=color, alpha=0.7, lw=0, zorder=2)
        lim = max(lim, np.abs(obs).max(), np.abs(pred).max())
    lim *= 1.1
    ax.plot([-lim, lim], [-lim, lim], color=MUTED, lw=1, ls="--", zorder=1)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal")
    for j, (k, color) in enumerate(TARGETS.items()):
        ax.annotate(f"{k}   c = {fit.table.loc[k, 'coupling']:.2f}", (0.03, 0.97 - 0.07 * j), xycoords="axes fraction",
                    va="top", fontsize=11, color=color)
    ax.set_xlabel(r"fitted $z_k f_i$ (both coordinates)")
    ax.set_ylabel(r"observed $\beta_{ik}$ (centred)")


def panel_estimates(ax, fit, truth):
    tab = fit.table.loc[list(COMPS)]
    y = np.arange(len(tab))
    for yi, (k, r) in zip(y, tab.iterrows()):
        color = TARGETS.get(k, INK)
        ax.hlines(yi, r.coupling_lo, r.coupling_hi, color=color, lw=2.5, zorder=1)
        ax.scatter([r.coupling], [yi], s=60, color=color, zorder=3)
        ax.scatter([truth[k]], [yi], s=110, marker="|", color=INK, linewidth=2, zorder=4)
    ax.set_yticks(y, [f"{k} (anchor)" if COMPS[k][1] else k for k in tab.index])
    ax.set_ylim(len(tab) - 0.5, -0.5)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlim(-0.02, 1.02)
    ax.set_xlabel("coupling c (dot, 95 % bootstrap interval; bar: true)")


def main():
    pf.use_style("slide")
    OUT.mkdir(parents=True, exist_ok=True)
    ds, hours = simulate()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        fit = ds.fit(n_boot=N_BOOT, reference="fibroblast")
    truth = true_coupling(ds.study, hours)
    tab = fit.table.assign(true_c=truth)
    tab["true_amp"] = [abs(c.z) if c.coupled else np.nan for c, _ in (COMPS[k] for k in tab.index)]
    tab["true_offset_h"] = [np.angle(c.z) * HOURS if c.coupled else np.nan for c, _ in (COMPS[k] for k in tab.index)]
    print(f"gate: p_pairing {fit.gate.p_pairing:.3f}, p_acrophase {fit.gate.p_acrophase:.3f}")
    print(tab[["coupling", "coupling_lo", "coupling_hi", "true_c", "rel_amplitude", "true_amp", "offset_h",
               "offset_lo", "offset_hi", "true_offset_h"]].round(3).to_string())
    draws = (("A", lambda ax: panel_host(ax, fit, hours)), ("B", lambda ax: panel_readout(ax, ds.study, fit)),
             ("C", lambda ax: panel_estimates(ax, fit, truth)))
    for key, draw in draws:
        lay = pf.Layout.slide("A", width=SIZE[0], height=SIZE[1])
        draw(lay.axes["A"])
        fig = lay.finalize()
        pf.save(fig, OUT / f"fit_{key}")


if __name__ == "__main__":
    main()
