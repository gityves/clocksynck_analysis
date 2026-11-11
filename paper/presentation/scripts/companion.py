"""Slide figures (method, after the readout): the pseudo-clock companion and why it is needed. Three one-panel figures.

A. Why (simulated): coupling is the share of a compartment's clock-plane variation that it shares with the host.
   Anything a donor's compartments share moves the clock genes too: genotype (eQTL), sex, age, agonal state, ambient
   RNA. Here a compartment with its own clock, unshared with the donor ("lymphocyte", coupled=False), is read against
   a host built from two anchors that do share a clock, while donor effects shared by all compartments (one value per
   donor and gene, sd on the x axis) grow. Its coupling rises from 0 with no shared clock; a truly coupled compartment
   ("hepatocyte") for reference. The clock-gene readout alone cannot tell the two apart.
B. The registered test (A150, real data, class-fixed rerun): mean coupling over 53 compartments of 5 single-cell
   atlases, clock genes against 200 expression-matched pseudo-clock sets (non-clock genes given the clock genes'
   amplitudes and acrophases; `paired.pseudo_waveforms`), each run through the same pipeline.
C. Per family (A150): median coupling with clock genes against the median and 95th percentile over pseudo sets.

Simulation (A): 30 donors x 4 compartments, deaths 08-18 h, anchors fibroblast and endothelial; 5 seeds per value of
the donor effect sd; coupling from `Dataset.fit(gate=False)` (point estimates). The gate passes in all 5 runs at sd 0
and 4 of 5 at sd 0.3 (checked separately, 1,000 permutations). Illustrative, not fitted to data.

Output: paper/presentation/figures/companion/companion_{A,B,C}.{svg,png} (150 x 120 mm each).

    python scripts/companion.py
"""
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

import clocksync as cs
import paperfig as pf

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parents[1] / "figures" / "companion"
A150 = ROOT / "paper" / "exploratory" / "A150_atlas_specificity" / "output"
SIZE = (150.0, 120.0)  # mm, one figure beside the bullets
HOURS = 24 / (2 * np.pi)
SD_GRID, SEEDS = [0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3], range(5)
COMPS = {
    "fibroblast": cs.Compartment(1.0, psi=0.03, mesor_sd=0.3),
    "endothelial": cs.Compartment(0.75 * np.exp(1j * 3 / HOURS), psi=0.03, mesor_sd=0.3),
    "hepatocyte": cs.Compartment(1.15 * np.exp(-1j * 4 / HOURS), psi=0.15, mesor_sd=0.3),
    "lymphocyte": cs.Compartment(0.8, psi=0.03, coupled=False, mesor_sd=0.3),
}
SHOW = {"hepatocyte": ("shared clock", pf.OKABE_ITO[2]), "lymphocyte": ("own clock, not shared", pf.OKABE_ITO[5])}
FAMILY = {"fibroblast_stromal": "fibroblast", "endothelial": "endothelial", "epithelial": "epithelial",
          "mural": "mural", "myeloid": "myeloid", "lymphocyte": "lymphocyte", "other": "other",
          "parenchyma": "parenchyma", "neuron": "neuron"}
INK, MUTED, NULL = "#222222", pf.GREY, "#c8c8c8"


def simulate_why():
    rows = []
    for sd in SD_GRID:
        for seed in SEEDS:
            df, _ = cs.simulate(cs.gtex(), COMPS, n=30, window_h=(8.0, 18.0), depth=4e6, tau=0.05, shared_sd=sd,
                                seed=seed)
            ds = cs.Dataset(df.rename(columns={"resection": "donor"}), anchors=["fibroblast", "endothelial"],
                            min_donors=5)
            t = ds.fit(n_boot=0, gate=False).table
            rows += [dict(sd=sd, seed=seed, compartment=k, coupling=t.loc[k, "coupling"]) for k in SHOW]
    return pd.DataFrame(rows)


def panel_why(ax, sim):
    rng = np.random.default_rng(0)
    for k, (label, color) in SHOW.items():
        d = sim[sim.compartment == k]
        ax.scatter(d.sd + rng.uniform(-0.006, 0.006, len(d)), d.coupling, s=18, color=color, alpha=0.5, lw=0)
        med = d.groupby("sd").coupling.median()
        ax.plot(med.index, med.values, color=color, lw=2.2)
        ax.annotate(label, (med.index[-1], med.values[-1]), xytext=(6, 0), textcoords="offset points", va="center",
                    fontsize=11, color=color, annotation_clip=False)
    ax.set_xlim(-0.01, 0.31)
    ax.set_ylim(-0.02, 1)
    ax.set_xlabel("donor effects (sex, genotype, ...) shared by all\ncompartments (sd, log expression)")
    ax.set_ylabel("coupling c (clock genes)")


def panel_pooled(ax):
    h = pd.read_csv(A150 / "h150_1.csv").iloc[0]
    null = pd.read_csv(A150 / "pseudo_pooled.csv").pooled_mean.to_numpy()
    bins = np.linspace(0.1, 0.45, 36)
    ax.hist(null, bins=bins, color=NULL, lw=0)
    ax.axvline(h.pseudo_q95, color=MUTED, lw=1, ls="--")
    ax.annotate("95th pct", (h.pseudo_q95, 1), xycoords=("data", "axes fraction"), xytext=(-4, -2),
                textcoords="offset points", ha="right", va="top", fontsize=10, color=MUTED)
    ax.axvline(h.clock_pooled_mean, color=INK, lw=2.2)
    ax.annotate(f"clock genes\n{h.clock_pooled_mean:.3f}\np = {h.p:.3f}", (h.clock_pooled_mean, 0.6),
                xycoords=("data", "axes fraction"), xytext=(6, 0), textcoords="offset points", va="center",
                fontsize=11, color=INK)
    ax.annotate(f"{int(h.n_sets)} pseudo-clock sets", (np.median(null), 0.3), xycoords=("data", "axes fraction"),
                ha="center", fontsize=10, color=MUTED)
    ax.set_xlabel(f"mean coupling over {int(h.n_compartments)} compartments ({int(h.n_atlases)} atlases)")
    ax.set_ylabel("pseudo-clock sets")
    ax.set_yticks([])
    ax.spines["left"].set_visible(False)


def panel_family(ax):
    f = pd.read_csv(A150 / "per_family.csv").set_index("family").loc[list(FAMILY)]
    f = f.sort_values("c_clock_median")
    y = np.arange(len(f))
    ax.hlines(y, f.c_pseudo_median, f.c_pseudo_q95, color=NULL, lw=6, zorder=1)
    ax.scatter(f.c_pseudo_median, y, s=30, color=MUTED, zorder=2)
    above = f.c_clock_median > f.c_pseudo_q95
    ax.scatter(f.c_clock_median, y, s=90, marker="D", zorder=3, color=np.where(above, INK, "white"), edgecolor=INK,
               linewidth=1.4)
    ax.set_yticks(y, [f"{FAMILY[k]} ({n})" for k, n in f.n_compartments.items()])
    ax.tick_params(axis="y", length=0)
    ax.spines["left"].set_visible(False)
    ax.set_xlim(0, 1)
    ax.set_xlabel("median coupling\n(diamond: clock genes; grey: pseudo-clock median to 95th pct)")


def main():
    pf.use_style("slide")
    OUT.mkdir(parents=True, exist_ok=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sim = simulate_why()
    print(sim.groupby(["compartment", "sd"]).coupling.median().unstack().round(2))
    for key, draw in (("A", lambda ax: panel_why(ax, sim)), ("B", panel_pooled), ("C", panel_family)):
        lay = pf.Layout.slide("A", width=SIZE[0], height=SIZE[1])
        draw(lay.axes["A"])
        fig = lay.finalize()
        pf.save(fig, OUT / f"companion_{key}")


if __name__ == "__main__":
    main()
