"""A145 step 1: the shared_sd values of arms S and N (no method is run).

Registered: shared_sd "the A140 sizes on the log scale, matched by the share of clock-plane variance it takes". The A140
non-clock factor (q1_lib: patient score N(0, lambda) on a unit-norm loading over a random 10-30 % of 2,000 genes, lambda =
78.7 at scale 1; scale 0.5 for the small size) enters the log rates, as the simulator's shared program does (sd shared_sd
per gene, independent over genes). Its clock-plane variance is the trace of Cov(beta) of the 14-gene state it induces:
A140: lambda * |L_p l|^2 (L_p the projector at the A140 depth and mesor; l the loading on the 14 clock genes), E over loadings
= lambda * |L_p|_F^2 / G (G = 2,000); A145: shared_sd^2 * |L_p|_F^2 (L_p at the A145 depth and mesor).
Share = factor / (factor + clock), clock = z^2 V + 2 psi of an anchor (A140: 1.122, V = 1, psi = 0.061; A145 arm F: V / 0.9).
shared_sd is solved so that the A145 share equals the A140 share; the match is then verified on simulated data
(single compartment, shared program only, Poisson-corrected trace of Cov(beta) at shared_sd vs 0).
Output: output/simulator_check.csv (quantity, value, note).
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[2]))
sys.path.insert(0, str(HERE.parents[1] / "A140_dip_cyclops"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
import lib as A  # noqa: E402
import q1_lib as Q  # noqa: E402
from clocksync._core.plane import gene_weights, projector  # noqa: E402
from clocksync.simulate import Compartment  # noqa: E402

G_A140 = 2000


def plane_projector(depth: float, mesor: np.ndarray) -> np.ndarray:
    """L_p as the Dataset builds it: weights from the median counts depth * exp(mesor)."""
    return projector(A.WF, gene_weights(depth * np.exp(mesor)[None, :]))


def main():
    rows = []
    add = lambda q, v, note="": rows.append(dict(quantity=q, value=float(v), note=note))
    P = Q.load_panel()
    lam1 = float(P["lam"][0])
    L140 = plane_projector(P["depth"], P["mesor"][:14])
    fro140 = float(np.sum(L140 ** 2))
    clock140 = 1.0 + 2 * Q.psi_for(Q.COUPLING, Q.ARMS["A"]["window"])  # arm D: whole day (V = 1), psi of the 08-18 h spread
    # Monte Carlo over loadings (A140's own generator) for the spread
    draws = []
    for s in range(4000):
        l = Q.loadings(P, 1, s)[:14, 0]
        draws.append(lam1 * float(np.sum((L140 @ l) ** 2)))
    draws = np.array(draws)
    add("a140_lambda_scale1", lam1, "q1/panel.npz lam[0]")
    add("a140_clock_plane_var_per_compartment", clock140, "z^2 V + 2 psi, V = 1, psi = 0.061")
    add("a140_factor_plane_var_expected_scale1", lam1 * fro140 / G_A140, "lambda |L_p|_F^2 / G")
    add("a140_factor_plane_var_mc_mean_scale1", draws.mean(), "4000 loading draws")
    add("a140_factor_plane_var_mc_median_scale1", np.median(draws), "")
    F1 = lam1 * fro140 / G_A140
    shares = {}
    for key, scale in (("large", 1.0), ("small", 0.5)):
        F = F1 * scale ** 2
        shares[key] = F / (F + clock140)
        add(f"a140_share_{key}", shares[key], f"A140 factor scale {scale}: F / (F + clock)")
    # A145 projector (GTEx mesor, default depth)
    L145 = plane_projector(4e6, np.asarray(A.WF.mesor))
    fro145 = float(np.sum(L145 ** 2))
    clock145 = 1.0 / A.C_ANCHOR  # z = 1, V = 1 (full day), psi = V (1 - c) / (2 c): z^2 V + 2 psi = V / c
    add("a145_clock_plane_var_anchor_full_day", clock145, "V / 0.9")
    add("a145_projector_frobenius_sq", fro145, "")
    for key in ("large", "small"):
        sh = shares[key]
        sd = float(np.sqrt(sh / (1 - sh) * clock145 / fro145))
        add(f"shared_sd_{key}", sd, "solved: sd^2 |L_p|_F^2 / (sd^2 |L_p|_F^2 + clock) = A140 share")
    # verification on simulated data: factor-only compartment, Poisson-corrected trace Cov(beta), shared_sd vs 0
    wf0 = A.WF
    for key in ("large", "small"):
        sd = [r for r in rows if r["quantity"] == f"shared_sd_{key}"][0]["value"]
        tr = {}
        for s in (0.0, sd):
            tab, _ = cs.simulate(wf0, {"X": Compartment(z=0.0, psi=0.0, coupled=False)}, n=3000, window_h=(0, 24), jitter_h=0.0,
                                 amplitude=0.0, shared_sd=s, seed=1)
            b, v, _ = A.states(tab.rename(columns={"resection": "donor"}), "X")
            tr[s] = np.sum((b - b.mean(0)) ** 2) / (len(b) - 1) - np.mean(v)
        F_sim = tr[sd] - tr[0.0]
        add(f"verify_factor_plane_var_{key}", F_sim, f"simulated, shared_sd {sd:.4f}: trace Cov(beta) minus the same at 0")
        add(f"verify_share_{key}", F_sim / (F_sim + clock145), f"target {shares[key]:.4f}")
    pd.DataFrame(rows).to_csv(A.OUT / "simulator_check.csv", index=False)
    print(pd.DataFrame(rows).to_string())


if __name__ == "__main__":
    main()
