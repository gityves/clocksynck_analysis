"""A377 (registered 2026-11-25): is the brain-death loss of the organ arrangement a displacement along tissue-specific
axes rather than a change of the relations between organ clocks?

Axial factor: rank-1 Takagi fit of the pseudo-covariance P = E[z_a z_b] on the observed off-diagonal pairs (diagonal
imputed with the off-diagonal pairs, as registered: tissue-own non-circular noise does not enter). Axis alignment
A = |sum_a w_a^2| / sum_a |w_a|^2 (1: one common direction); axial size: pooled pseudo-covariance ratio (A109
`circ_ratio`).

Stages:
  power  V377.1 gate on the A133 simulator (ventilated design, fast-death offsets, n 79, 200 replicates per cell, f 0.2
         and 0.44, common vs tissue-specific directions); SD of D_A = A(post-mortem) - A(organ donors) and of the
         circularity-ratio difference under the A317 relabelling scheme (2000 relabellings of Hardy 0 donors into two
         groups of 79). No group-split observed value is computed.
  score  (written after the power note, before any group-split value) observed A and circularity ratio per group
         (post-mortem: all 79; organ donors: median over the 20 A317 draws); D_A with its relabelling p (same 2000
         relabellings as the power stage); H377.1 criteria; companions (200 matched, 200 responsive sets, same draws);
         per-tissue axis direction arg(w_a) / 2 [h, mod 12] and |w_a| (organ donors: fit on all organ donors,
         descriptive); host Rayleigh and axial R (host frame from the leading eigenvector of the cross-covariance).
         H377.2: the A133 generative model with each group's fitted axial factor (direction w_a / |w_a|, relative
         strength |w_a| / a_a scaled to mean square 1) in place of a common direction; f calibrated on a grid (0.05 to
         0.95, 50 replicates per value) to the group's observed circularity ratio; 200 replicates at that f, n 79;
         reference an independent f 0 draw at n 244; Fisher-Lee rho of pair phases with the reference.

    python a377.py --stage power | score
Outputs (`../output/`): gate.csv (cells), gate_reps.csv, power.csv | tests.csv, groups.csv, companions.csv, axes.csv,
sufficiency.csv, sufficiency_reps.csv. Aggregates only.
"""
import argparse
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve()
EXP = HERE.parents[2]
sys.path.insert(0, str(EXP / "A317_braindeath_confirm" / "scripts"))
sys.path.insert(0, str(EXP / "A133_axial_decomposition" / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from sklearn.metrics import roc_auc_score  # noqa: E402

import compute_v133_sims as S  # noqa: E402
import score_a317 as A  # noqa: E402
from axial import takagi_rank1  # noqa: E402
from score_a317 import B, H, gtex, raw  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402

OUT = HERE.parents[1] / "output"
N_REP, N_PERM, SEED = 200, 2000, 377


def fit_offdiag(P, observed, n_iter=500, tol=1e-10):
    """Rank-1 Takagi factor w (P ~ w w^T) fitted on the observed off-diagonal entries; every other entry, the diagonal
    included, imputed from the current fit (A133 `fit_axial` with the diagonal treated as unobserved)."""
    obs = observed | observed.T
    np.fill_diagonal(obs, False)
    P = (P + P.T) / 2
    w = takagi_rank1(np.where(obs, P, 0))
    for _ in range(n_iter):
        w_new = takagi_rank1(np.where(obs, P, np.outer(w, w)))
        done = np.linalg.norm(np.outer(w_new, w_new) - np.outer(w, w)) <= tol * max(np.linalg.norm(np.outer(w, w)), 1e-300)
        w = w_new
        if done:
            break
    return w


def alignment(w):
    return float(abs(np.sum(w ** 2)) / max(np.sum(abs(w) ** 2), 1e-300))


def circ_ratio(C, P, observed):
    iu = np.triu_indices(len(C), 1)
    o = observed[iu]
    return float(np.mean(abs(P[iu][o])) / np.mean(abs(C[iu][o])))


def readout(Z, M, observed):
    C, P, _ = H.cross_moments(Z, M)
    return alignment(fit_offdiag(P, observed)), circ_ratio(C, P, observed)


def centred(Z, M, idx):
    m = M[idx]
    mu = (Z[idx] * m).sum(0) / np.maximum(m.sum(0), 1)
    return (Z[idx] - mu) * m, m


def gate():
    a, V, phi = S.design()
    rng = np.random.default_rng(SEED)
    rows = []
    for f in (0.2, 0.44):
        for u in ("common", "pertissue"):
            for r in range(N_REP):
                Z, M = S.simulate(79, f, 0.0, u, a, V, phi, rng)
                _, _, N = H.cross_moments(Z, M)
                al, cr = readout(Z, M, N >= S.MIN_N)
                rows.append(dict(f=f, directions=u, rep=r, A=al, circ_ratio=cr))
    reps = pd.DataFrame(rows)
    cells = reps.groupby(["f", "directions"]).agg(A_median=("A", "median"), A_q05=("A", lambda x: x.quantile(0.05)),
                                                 A_q95=("A", lambda x: x.quantile(0.95)),
                                                 circ_ratio_median=("circ_ratio", "median")).reset_index()
    auc = {f: roc_auc_score((reps[reps.f == f].directions == "common").astype(int), reps[reps.f == f].A)
           for f in (0.2, 0.44)}
    cells["auc"] = cells.f.map(auc)
    c44 = cells[cells.f == 0.44].set_index("directions")
    passes = bool(c44.A_median["common"] >= 0.8 and c44.A_median["pertissue"] <= 0.5 and auc[0.44] >= 0.9)
    cells["V377_1_passes"] = passes
    return reps, cells


def setup(sets=False):
    T = H.Tissues(200, 0)
    T.build_states(np.random.default_rng(1), 0, real=True)
    if sets:
        A.RS.add_sets(T, 200, 0)
    coh = gtex.phenotypes(B.PH_COLS).reindex(T.donors).COHORT.to_numpy()
    fast, h0 = T.group_donors("fast"), T.group_donors("vent")
    od, pm = h0[coh[h0] == B.OPO], h0[coh[h0] == B.PM]
    rng = np.random.default_rng(A.SEED)
    od_draws = [rng.choice(od, len(pm), replace=False) for _ in range(A.N_DRAWS)]
    Z, M = raw(T, "clock", 0, T.R)
    mask = B.pair_mask(M, fast, pm, *od_draws[:1]) & B.pair_mask(M, od)  # the A317 pairs
    obs = np.zeros((len(H.BODY), len(H.BODY)), bool)
    for (a, b), m in zip(H.PAIRS, mask):
        obs[a, b] = obs[b, a] = m
    return T, Z, M, mask, obs, h0, od, pm, od_draws


def relabel_null(Z, M, obs, h0, n):
    rng = np.random.default_rng(SEED)
    d_a, d_c = [], []
    for _ in range(N_PERM):  # relabelled Hardy 0 donors only: no group-split value
        p = rng.permutation(h0)
        ra, rc = readout(*centred(Z, M, p[:n]), obs)
        oa, oc = readout(*centred(Z, M, p[n:2 * n]), obs)
        d_a.append(ra - oa)
        d_c.append(rc - oc)
    return np.array(d_a), np.array(d_c)


def power():
    limit_threads(4)
    T, Z, M, mask, obs, h0, od, pm, od_draws = setup()
    d_a, d_c = relabel_null(Z, M, obs, h0, len(pm))
    return pd.Series(dict(pairs=int(mask.sum()), n_group=len(pm), sd_D_A=float(np.std(d_a, ddof=1)),
                          mde_D_A=float(2.49 * np.std(d_a, ddof=1)), q95_D_A=float(np.quantile(d_a, 0.95)),
                          sd_D_circ=float(np.std(d_c, ddof=1)), mde_D_circ=float(2.49 * np.std(d_c, ddof=1))))


def group(Z, M, obs, pm, od_draws):
    """(A, circularity ratio) of the post-mortem donors and the organ donors (median over draws)."""
    a_pm, c_pm = readout(*centred(Z, M, pm), obs)
    r = np.array([readout(*centred(Z, M, d), obs) for d in od_draws])
    return a_pm, c_pm, float(np.median(r[:, 0])), float(np.median(r[:, 1]))


def host_R(Z, M, idx):
    Zc, m = centred(Z, M, idx)
    C, _, _ = H.cross_moments(Zc, m)
    u, _ = H.host_frame(C)
    phi = np.angle(H.host_states(Zc, m, u))
    return float(abs(np.mean(np.exp(1j * phi)))), float(abs(np.mean(np.exp(2j * phi))))


def simulate_fitted(n, f, w, a, V, phi, rng):
    """A133 `simulate` (uniform host phases) with the axial direction and relative strength of the fitted w."""
    T = len(a)
    rel = abs(w) / a
    rel = rel / np.sqrt(np.mean(rel ** 2))
    d = w / np.maximum(abs(w), 1e-300) * rel
    theta = rng.uniform(0, 2 * np.pi, n)[:, None]
    s_ = rng.standard_normal(n)[:, None]
    xi = (rng.standard_normal((n, T)) + 1j * rng.standard_normal((n, T))) / np.sqrt(2)
    z = a * (np.sqrt(S.C_COUP * f) * s_ * d + np.sqrt(S.C_COUP * (1 - f)) * np.exp(1j * (theta + phi))
             + np.sqrt(1 - S.C_COUP) * xi)
    z = z + np.sqrt(V / 2) * (rng.standard_normal((n, T)) + 1j * rng.standard_normal((n, T)))
    return z - z.mean(0), np.ones((n, T))


def sufficiency(w, observed_rho, target_ratio, label, rng):
    a, V, phi = S.design()
    full = np.ones((len(a), len(a)), bool)
    grid = []
    for f in np.round(np.arange(0.05, 0.96, 0.05), 2):
        cr = [circ_ratio(*H.cross_moments(*simulate_fitted(79, f, w, a, V, phi, rng))[:2], full) for _ in range(50)]
        grid.append((f, float(np.median(cr))))
    f = min(grid, key=lambda x: abs(x[1] - target_ratio))[0]
    reps = []
    for r in range(N_REP):
        Cg = H.cross_moments(*simulate_fitted(79, f, w, a, V, phi, rng))[0]
        Cr = H.cross_moments(*S.simulate(244, 0.0, 0.0, "common", a, V, phi, rng))[0]
        reps.append(dict(group=label, f=f, rep=r, rho=float(H.fisher_lee_fast(S.pair_phases(Cg), S.pair_phases(Cr)))))
    reps = pd.DataFrame(reps)
    q05, q50, q95 = reps.rho.quantile([0.05, 0.5, 0.95])
    row = dict(group=label, target_ratio=target_ratio, f=f, sim_ratio=dict(grid)[f], rho_q05=q05, rho_median=q50,
               rho_q95=q95, observed_rho=observed_rho,
               verdict="inside" if q05 <= observed_rho <= q95 else ("below" if observed_rho < q05 else "above"))
    return row, reps


def score():
    limit_threads(4)
    T, Z, M, mask, obs, h0, od, pm, od_draws = setup(sets=True)
    a_pm, c_pm, a_od, c_od = group(Z, M, obs, pm, od_draws)
    d_a, d_c = relabel_null(Z, M, obs, h0, len(pm))
    D_A, D_C = a_pm - a_od, c_pm - c_od
    gate_q95 = float(pd.read_csv(OUT / "gate.csv").set_index(["f", "directions"]).loc[(0.44, "pertissue"), "A_q95"])
    p_a = float((1 + np.sum(d_a >= D_A)) / (1 + len(d_a)))
    tests = pd.Series(dict(A_pm=a_pm, A_od=a_od, D_A=D_A, p_D_A=p_a, circ_pm=c_pm, circ_od=c_od, D_circ=D_C,
                           p_D_circ_two_sided=float((1 + np.sum(abs(d_c) >= abs(D_C))) / (1 + len(d_c))),
                           pertissue_A_q95=gate_q95,
                           H377_1_passes=bool(D_A > 0 and p_a < 0.05 and a_od < gate_q95)))
    comp = []
    for kind, lab in (("pseudo", "matched"), ("resp", "responsive")):
        for k in range(200):
            Zk, Mk = raw(T, kind, k, T.R)
            ap, cp, ao, co = group(Zk, Mk, obs, pm, od_draws)
            comp.append(dict(family=lab, set=k, A_pm=ap, A_od=ao, D_A=ap - ao, circ_pm=cp, circ_od=co))
    comp = pd.DataFrame(comp)
    axes, groups = [], []
    for g, idx in (("post-mortem", pm), ("organ donors (all)", od)):
        Zc, m = centred(Z, M, idx)
        C, P, _ = H.cross_moments(Zc, m)
        w = fit_offdiag(P, obs)
        axes += [dict(group=g, tissue=t, axis_h=float((np.angle(w[j]) * H.HOURS) % 12), strength=float(abs(w[j])))
                 for j, t in enumerate(H.BODY)]
        rR, aR = host_R(Z, M, idx)
        groups.append(dict(group=g, n=len(idx), A=alignment(w), circ_ratio=circ_ratio(C, P, obs), host_rayleigh_R=rR,
                           host_axial_R=aR))
    rng = np.random.default_rng(SEED + 1)
    w_od = fit_offdiag(H.cross_moments(*centred(Z, M, od))[1], obs)
    w_pm = fit_offdiag(H.cross_moments(*centred(Z, M, pm))[1], obs)
    rho = pd.read_csv(A.OUT / "posthoc_loto.csv").set_index("dropped").loc["none"]  # A317 reproduction, 65 pairs
    s_od, r_od = sufficiency(w_od, float(rho.rho_od), c_od, "organ donors", rng)
    s_pm, r_pm = sufficiency(w_pm, float(rho.rho_pm), c_pm, "post-mortem (model check)", rng)
    return tests, comp, pd.DataFrame(axes), pd.DataFrame(groups), pd.DataFrame([s_od, s_pm]), pd.concat([r_od, r_pm])


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=["power", "score"], required=True)
    stage = ap.parse_args().stage
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    if stage == "score":
        tests, comp, axes, groups, suff, reps = score()
        tests.to_csv(OUT / "tests.csv", header=["value"])
        comp.to_csv(OUT / "companions.csv", index=False)
        axes.to_csv(OUT / "axes.csv", index=False)
        groups.to_csv(OUT / "groups.csv", index=False)
        suff.to_csv(OUT / "sufficiency.csv", index=False)
        reps.to_csv(OUT / "sufficiency_reps.csv", index=False)
        print(tests.to_string())
        print(comp.groupby("family")[["A_pm", "A_od", "D_A", "circ_pm", "circ_od"]].quantile([0.5, 0.95]).round(3).to_string())
        print(groups.round(3).to_string(index=False))
        print(axes.round(2).to_string(index=False))
        print(suff.round(3).to_string(index=False), f"\n({time.time() - t0:.0f}s)")
        return
    reps, cells = gate()
    reps.to_csv(OUT / "gate_reps.csv", index=False)
    cells.to_csv(OUT / "gate.csv", index=False)
    print(cells.round(3).to_string(index=False), f"\n({time.time() - t0:.0f}s)", flush=True)
    pw = power()
    pw.to_csv(OUT / "power.csv", header=["value"])
    print(pw.round(4).to_string(), f"\n({time.time() - t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
