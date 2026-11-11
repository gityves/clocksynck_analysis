"""Slide figures (method, slide M3): the gate. Three one-panel figures, shown in order beside the bullets.

A. What the gate looks at: the cross-covariance C_kl (genes x genes, across donors) between two compartments of the
   same donors (fibroblast x hepatocyte), genes ordered by acrophase; drawn as a correlation (C_kl divided by the
   genes' SDs) so that low-amplitude genes show. Measurement noise is independent between compartments, so C_kl
   carries only what the two share. A shared clock gives C_kl = X z_k Sigma_f z_l^* X^T: rank 2, in the clock plane,
   sign set by the acrophases (genes peaking together co-vary, opposite genes anti-co-vary). With deaths at 08-18 h
   the host spread Sigma_f is elongated, so the pattern is blocks rather than bands parallel to the diagonal.
B. Null 1, shuffled donor pairing (is there shared signal in the clock plane?): clock_cov = sum_kl ||P C_kl P||^2,
   P the orthogonal projector on the clock plane, against donors permuted within each compartment.
C. Null 2, permuted acrophases (is the shared signal clock-shaped?): clock_share = clock_cov / sum_kl ||C_kl||^2
   against planes built from acrophases shuffled across genes.

Three simulated cohorts (30 donors x 3 compartments, deaths 08-18 h, loadings as in model_animation.py):
- "shared clock": every compartment follows the donor's clock: passes both;
- "clocks not shared": every compartment has a clock, at a phase independent of the donor's other compartments: its
  cross-covariance is clock-shaped noise, so the share is high, but pairing finds nothing: fails null 1;
- "shared non-clock factor": no shared clock, one donor factor with a random gene loading: shared signal, mostly
  outside the plane: fails null 2.
The seeds are chosen so that each cohort shows its case (illustrative; the calibrated pass rates are A154 and A181).
The nulls are recomputed here exactly as `Study.gate` does (same random stream); the p values are checked against
`Dataset.gate`.

Output: paper/presentation/figures/gate/gate_{A,B,C}.{svg,png} (150 x 120 mm each).

    python scripts/gate.py
"""
from pathlib import Path

import numpy as np
from matplotlib.colors import TwoSlopeNorm

import clocksync as cs
import paperfig as pf
from clocksync._core.controls import shuffle_pairing
from clocksync._core.study import _cross_covariances, centred_log_fraction, cross_covariances, plane

OUT = Path(__file__).resolve().parents[1] / "figures" / "gate"
SIZE = (150.0, 120.0)  # mm, one figure beside the bullets
N, WINDOW, DEPTH, N_PERM = 30, (8.0, 18.0), 4e6, 1000
HOURS = 24 / (2 * np.pi)
Z = {"fibroblast": 1.0 + 0j, "endothelial": 0.75 * np.exp(1j * 3 / HOURS), "hepatocyte": 1.15 * np.exp(-1j * 4 / HOURS)}
COHORTS = {  # label: (shared clock, shared factor loading sd, seed)
    "shared clock": (True, 0.0, 0),
    "clocks not shared": (False, 0.0, 0),
    "shared non-clock factor": (False, 0.4, 2),
}
INK, MUTED, NULL = "#222222", pf.GREY, "#c8c8c8"


def simulate(coupled, load_sd, seed):
    wf = cs.gtex()
    comps = {k: cs.Compartment(z, psi=0.05, coupled=coupled, mesor_sd=0.3) for k, z in Z.items()}
    load = np.random.default_rng(100 + seed).normal(0.0, load_sd, len(wf.genes)) if load_sd else None
    df, _ = cs.simulate(wf, comps, n=N, window_h=WINDOW, depth=DEPTH, tau=0.05, shared_load=load, seed=seed)
    return cs.Dataset(df.rename(columns={"resection": "donor"}), anchors=["fibroblast"], min_donors=5)


def nulls(st, seed=0):
    """Observed statistics and both null distributions, as `Study.gate` computes them (same random stream)."""
    rng = np.random.default_rng(seed)
    y = centred_log_fraction(st.counts, st.depth)
    cc = lambda Cs, P: sum(np.sum((P @ C @ P) ** 2) for C in Cs)
    Cs = _cross_covariances(y, st.resection, st.compartment, st.stratum)
    P0 = plane(st.waveform)
    cl, tot = cc(Cs, P0), sum(np.sum(C ** 2) for C in Cs)
    null_cl = np.array([cc(_cross_covariances(y, shuffle_pairing(st.resection, st.group, rng), st.compartment,
                                              st.stratum), P0) for _ in range(N_PERM)])
    null_sh = np.array([cc(Cs, plane(st.waveform.permuted(rng))) / tot for _ in range(N_PERM)])
    p = lambda null, obs: (1 + np.sum(null >= obs)) / (1 + len(null))
    return dict(cl=cl, null_cl=null_cl, p_pair=p(null_cl, cl), share=cl / tot, null_sh=null_sh,
                p_acro=p(null_sh, cl / tot))


def panel_matrix(ax, st):
    wf = st.waveform
    y = centred_log_fraction(st.counts, st.depth)
    C = {(k, l): Cv for k, l, _, Cv in cross_covariances(y, st.resection, st.compartment, st.stratum)}
    sd = {k: y[st.compartment == k].std(0, ddof=1) for k in ("fibroblast", "hepatocyte")}
    order = np.argsort(np.mod(wf.psi, 2 * np.pi))
    M = (C[("fibroblast", "hepatocyte")] / np.outer(sd["fibroblast"], sd["hepatocyte"]))[np.ix_(order, order)]
    lim = 1.0
    im = ax.imshow(M, cmap="RdBu_r", norm=TwoSlopeNorm(0, -lim, lim), interpolation="nearest")
    cb = ax.figure.colorbar(im, cax=ax.inset_axes([1.04, 0.3, 0.035, 0.4]), ticks=[-1, 0, 1])
    cb.outline.set_visible(False)
    cb.set_label("correlation", color=MUTED)
    genes = [wf.genes[g] for g in order]
    ax.set_xticks(range(len(genes)), genes, rotation=90, fontsize=8)
    ax.set_yticks(range(len(genes)), genes, fontsize=8)
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_xlabel("hepatocyte gene (by acrophase)")
    ax.set_ylabel("fibroblast gene (by acrophase)")


def panel_null(ax, res, key):
    """One row per cohort: the null distribution (grey dots) and the observed statistic (filled: p < 0.05)."""
    rng = np.random.default_rng(1)
    for row, (label, r) in enumerate(res.items()):
        if key == "pair":  # relative to the null mean, so the cohorts share one axis
            null, obs, p = r["null_cl"] / r["null_cl"].mean(), r["cl"] / r["null_cl"].mean(), r["p_pair"]
        else:
            null, obs, p = r["null_sh"], r["share"], r["p_acro"]
        ax.scatter(null, row + rng.uniform(-0.18, 0.18, len(null)), s=4, color=NULL, lw=0, zorder=1)
        ax.scatter([obs], [row], s=110, marker="D", zorder=3, color=INK if p < 0.05 else "white", edgecolor=INK,
                   linewidth=1.4)
        ax.annotate(f"p = {p:.3f}" if p >= 0.001 else "p < 0.001", (1, row), xycoords=("axes fraction", "data"),
                    xytext=(6, 0), textcoords="offset points", va="center", fontsize=10,
                    color=INK if p < 0.05 else MUTED, annotation_clip=False)
    ax.set_yticks(range(len(res)), list(res))
    ax.set_ylim(len(res) - 0.5, -0.5)
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    if key == "pair":
        ax.set_xscale("log")
        ax.axvline(1, color="#e6e6e6", lw=0.8, zorder=0)
        ax.set_xlabel("shared clock-plane covariance\n(relative to shuffled pairing)")
    else:
        ax.set_xlim(0, 1)
        ax.set_xlabel("clock share of the shared covariance\n(grey: permuted acrophases)")


def main():
    pf.use_style("slide")
    OUT.mkdir(parents=True, exist_ok=True)
    res, studies = {}, {}
    for label, (coupled, load_sd, seed) in COHORTS.items():
        ds = simulate(coupled, load_sd, seed)
        studies[label], res[label] = ds.study, nulls(ds.study)
        g = ds.gate(n_perm=N_PERM, seed=0)
        assert np.isclose(g.p_pairing, res[label]["p_pair"]) and np.isclose(g.p_acrophase, res[label]["p_acro"])
        print(f"{label}: p_pairing {g.p_pairing:.3f}, p_acrophase {g.p_acrophase:.3f}, share {g.clock_share:.2f} "
              f"(null mean {g.share_null_mean:.2f})")
    draws = (("A", lambda ax: panel_matrix(ax, studies["shared clock"])),
             ("B", lambda ax: panel_null(ax, res, "pair")),
             ("C", lambda ax: panel_null(ax, res, "acro")))
    for key, draw in draws:
        lay = pf.Layout.slide("A", width=SIZE[0], height=SIZE[1])
        draw(lay.axes["A"])
        fig = lay.finalize()
        pf.save(fig, OUT / f"gate_{key}")


if __name__ == "__main__":
    main()
