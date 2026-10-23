"""A109c D109c.3 / D109c.4 (descriptive): is the 0 / 12 h push of ventilated GTEx donors found in fatal COVID-19?

Pre-registered in paper/docs/registered_experiments.md, section A109c.
  D109c.3  Delorey 2021 whole-organ nuclei (heart, liver, kidney, lung; >= 50 nuclei), 6 pairs: pooled circularity
           ratio |E[z_a z_b]| / |E[z_a conj z_b]|, median folded pair phase distance to 0 / 12 h, Fisher-Lee with the
           GTEx fast-death relations of the matching tissues (A92 list: heart -> Heart - Atrial Appendage, lung ->
           Lung; liver and kidney have no match, so one pair is matched and Fisher-Lee is undefined), against 200
           simulations each (rotating, 50 % mixture, pure displacement) at the Delorey design and the GTEx fast /
           ventilated values of the matching pair and of all 66 pairs.
  D109c.4  Melms 2021 lung classes (>= 20 cells in >= 10 COVID donors, then >= 20 cells in the controls): the same
           ratio and folded distance, COVID donors vs the 7 surgical controls, donor-label permutation p (1,000) of
           the ratio difference (reported, not tested).
States: `clocksync.Dataset` (GTEx waveform, Poisson variance), centred within compartment over the donors of the arm;
pairs with fewer than MIN_SHARED shared donors are dropped. Statistics and simulator: A109 `_helpers.py`.

    python compute_covid_push.py [--subsample N] [--seed S]
Outputs (exploratory/A109_coupling_nulls/output/, .../subsample/ for checks): a109c_covid_delorey.csv,
a109c_covid_melms.csv, a109c_covid_sims.csv (dataset: delorey, melms_covid, melms_control).
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import _helpers as H  # noqa: E402
from clocksync import stats as ST  # noqa: E402
from common import cli, data  # noqa: E402
from common.paths import exploratory_dir  # noqa: E402
from results.figure_2.scripts.embryo import _shared as S  # noqa: E402

FOLDER = "A109_coupling_nulls"
MIN_SHARED = 5
GTEX_OF = {"heart": "Heart - Atrial Appendage", "lung": "Lung"}  # A92 tissues of the Delorey organs
MODELS = (("rotating", 0.0), ("mix50", 0.5), ("displacement", 1.0))


def matrices(pb: pd.DataFrame, donors_df, anchors=None):
    """Raw donor x compartment complex states, mask, Poisson variance (trace), donor and compartment labels."""
    ds = S.dataset(pb, donors_df, anchors=anchors)
    d = ds.data
    don = sorted(set(d.resection))
    comp = sorted(set(d.compartment))
    di, ci = {x: i for i, x in enumerate(don)}, {x: i for i, x in enumerate(comp)}
    Z, M, V = np.zeros((len(don), len(comp)), complex), np.zeros((len(don), len(comp))), np.zeros((len(don), len(comp)))
    for r in range(len(d.compartment)):
        i, j = di[d.resection[r]], ci[d.compartment[r]]
        Z[i, j] = d.beta[r, 0] + 1j * d.beta[r, 1]
        M[i, j] = 1
        V[i, j] = np.trace(d.V[r])
    return Z, M, V, np.array(don), np.array(comp)


def centre(Z, M, idx):
    """Rows idx, centred within compartment over the donors present."""
    m = M[idx]
    z = Z[idx] * m
    mu = z.sum(0) / np.maximum(m.sum(0), 1)
    return (z - mu * m) * m, m


def stats_of(Zc, m):
    """(ratio, folded distance h, n pairs, pair table) over the pairs with >= MIN_SHARED shared donors."""
    Hm, Pc, N = H.cross_moments(Zc, m)
    iu = np.triu_indices(len(Hm), 1)
    ok = N[iu] >= MIN_SHARED
    if ok.sum() == 0:
        return np.nan, np.nan, 0, None
    h, p = Hm[iu][ok], Pc[iu][ok]
    return float(abs(p).mean() / abs(h).mean()), H.folded_distance(np.angle(h)), int(ok.sum()), (iu, ok, h, p, N[iu][ok])


def coupling(Zc, m, V):
    """Mean over pairs of |noise-corrected correlation| (the C of A109)."""
    Hm, _, N = H.hermitian(Zc, m, V)
    d = np.sqrt(np.maximum(Hm.diagonal().real, 1e-12))
    c = abs(Hm) / np.outer(d, d)
    iu = np.triu_indices(len(Hm), 1)
    ok = N[iu] >= MIN_SHARED
    return float(np.clip(c[iu][ok].mean(), 0.01, 0.99))


def delorey(args, run, rng, n_sim):
    organ, cl = data.load("covid_delorey", level="organ", subsample=args.subsample, seed=args.seed)
    t = S.eligible(organ, 50)
    Z, M, V, don, comp = matrices(t, cl)
    Zc, m = centre(Z, M, np.arange(len(don)))
    ratio, dist, npair, tab = stats_of(Zc, m)
    iu, ok, h, p, n = tab
    names = [(comp[a], comp[b]) for a, b in zip(*iu)]
    names = [x for x, k in zip(names, ok) if k]
    rel = pd.read_csv(H.PAPER / "exploratory" / FOLDER / "output" / "pair_relations.csv")
    rows = []
    for (a, b), hh, pp, nn in zip(names, h, p, n):
        r = dict(level="pair", organ_a=a, organ_b=b, n_shared=int(nn), offset_h=float(np.angle(hh) * H.HOURS),
                 folded_h=H.folded_distance(np.array([np.angle(hh)])), ratio=float(abs(pp) / abs(hh)))
        if a in GTEX_OF and b in GTEX_OF:
            for g in ("fast", "vent"):
                q = rel[(rel.group == g) & (rel.tissue_a == GTEX_OF[a]) & (rel.tissue_b == GTEX_OF[b])]
                r[f"gtex_{g}_offset_h"] = float(q.offset_h.iloc[0])
                r[f"gtex_{g}_folded_h"] = H.folded_distance(np.array([q.offset_h.iloc[0] / H.HOURS]))
                r[f"gtex_{g}_ratio"] = float(q.pseudo_cov_ratio.iloc[0])
        rows.append(r)
    # Fisher-Lee over the matched pairs (needs >= 3)
    mt = [(i, a, b) for i, (a, b) in enumerate(names) if a in GTEX_OF and b in GTEX_OF]
    fl = np.nan
    if len(mt) >= 3:
        f = rel[rel.group == "fast"].set_index(["tissue_a", "tissue_b"]).offset_h
        fl = ST.fisher_lee(np.angle(h)[[i for i, _, _ in mt]], [f[(GTEX_OF[a], GTEX_OF[b])] / H.HOURS for _, a, b in mt])
    summ = dict(level="pooled", n_donors=len(don), n_organs=len(comp), n_pairs=npair, ratio=ratio, folded_h=dist,
                C=coupling(Zc, m, V), n_matched_pairs=len(mt), fisher_lee_vs_gtex_fast=fl,
                organs_without_gtex_match=";".join(c for c in comp if c not in GTEX_OF))
    for g in ("fast", "vent"):  # all 66 GTEx pairs (reference scale)
        q = rel[rel.group == g]
        summ[f"gtex_{g}_all_ratio"] = float((q.cov_modulus * q.pseudo_cov_ratio).sum() / q.cov_modulus.sum())
        summ[f"gtex_{g}_all_folded_h"] = H.folded_distance(q.offset_h.to_numpy() / H.HOURS)
    rows.append(summ)
    # simulations at the Delorey design (donors, organs, depths, observed spread and coupling)
    Cobs = summ["C"]
    sims = [dict(model="observed", rep=-1, ratio=ratio, folded_h=dist, C_target=Cobs)]
    for name, f in MODELS:
        for r in range(n_sim):
            zs = H.simulate_group(Zc, m, V, Cobs, f, H.GC_H, rng)
            ra, di, _, _ = stats_of(zs, m)
            sims.append(dict(model=name, rep=r, ratio=ra, folded_h=di, C_target=Cobs))
    sims = pd.DataFrame(sims)
    sims.insert(0, "dataset", "delorey")
    for name, _ in MODELS:
        s = sims[sims.model == name]
        summ[f"sim_{name}_ratio_median"] = float(s.ratio.median())
        summ[f"sim_{name}_ratio_q025"] = float(s.ratio.quantile(0.025))
        summ[f"sim_{name}_ratio_q975"] = float(s.ratio.quantile(0.975))
        summ[f"sim_{name}_folded_median"] = float(s.folded_h.median())
        summ[f"sim_{name}_frac_ratio_ge_obs"] = float((s.ratio >= ratio).mean())
    return pd.DataFrame(rows), sims


def simulate(Zc, m, V, n_sim, rng, label):
    """Observed row and n_sim simulations per model at the design (donors, classes, depths) of Zc."""
    ratio, dist, _, _ = stats_of(Zc, m)
    C = coupling(Zc, m, V)
    rows = [dict(dataset=label, model="observed", rep=-1, ratio=ratio, folded_h=dist, C_target=C)]
    for name, f in MODELS:
        for r in range(n_sim):
            ra, di, _, _ = stats_of(H.simulate_group(Zc, m, V, C, f, H.GC_H, rng), m)
            rows.append(dict(dataset=label, model=name, rep=r, ratio=ra, folded_h=di, C_target=C))
    return pd.DataFrame(rows)


def melms(args, run, rng, n_sim):
    pb, ms = data.load("covid_melms")  # always all 26 donors (the 10-donor eligibility needs them); --subsample only cuts permutations
    covid, ctrl = ms.index[ms.group != "Control"], ms.index[ms.group == "Control"]
    t = S.eligible(pb, 20, 10, donors=covid)  # classes passing in the COVID arm
    cls = sorted(t.compartment.unique())
    t = S.eligible(pb, 20, compartments=cls)
    Z, M, V, don, comp = matrices(t, ms)
    iscov = np.isin(don, covid)
    ic, ik = np.flatnonzero(iscov), np.flatnonzero(~iscov)

    def arm(idx):
        Zc, m = centre(Z, M, idx)
        return stats_of(Zc, m)[:3]

    rc, dc, nc = arm(ic)
    rk, dk, nk = arm(ik)
    n_perm = run["n_perm"]
    null = np.full(n_perm, np.nan)
    nulld = np.full(n_perm, np.nan)
    for i in range(n_perm):
        pm = rng.permutation(len(don))
        a, b = arm(pm[:len(ic)]), arm(pm[len(ic):])
        null[i], nulld[i] = a[0] - b[0], a[1] - b[1]
    d = rc - rk
    p = (1 + np.nansum(abs(null) >= abs(d))) / (1 + np.isfinite(null).sum())
    pd_ = (1 + np.nansum(abs(nulld) >= abs(dc - dk))) / (1 + np.isfinite(nulld).sum())
    sims, extra = [], {}
    for lab, idx in (("covid", ic), ("control", ik)):
        Zc, m = centre(Z, M, idx)
        sm = simulate(Zc, m, V[idx], n_sim, rng, f"melms_{lab}")
        sims.append(sm)
        obs = sm.ratio.iloc[0]
        extra[f"C_{lab}"] = float(sm.C_target.iloc[0])
        for name, _ in MODELS:
            x = sm[sm.model == name].ratio
            extra.update({f"sim_{lab}_{name}_ratio_median": float(x.median()),
                          f"sim_{lab}_{name}_ratio_q025": float(x.quantile(0.025)),
                          f"sim_{lab}_{name}_ratio_q975": float(x.quantile(0.975)),
                          f"sim_{lab}_{name}_frac_ratio_ge_obs": float((x >= obs).mean())})
    return pd.DataFrame([dict(
        classes=";".join(comp), n_classes=len(comp), donors_covid=len(ic), donors_control=len(ik),
        pairs_covid=nc, pairs_control=nk, ratio_covid=rc, ratio_control=rk, ratio_diff=d, p_perm_ratio=p,
        folded_h_covid=dc, folded_h_control=dk, folded_diff=dc - dk, p_perm_folded=pd_, n_perm=n_perm, **extra)]), pd.concat(sims)


def main(argv=None):
    args = cli.parser(__doc__).parse_args(argv)
    run = cli.setup(args)
    out = exploratory_dir(FOLDER, run["quick"])
    rng = np.random.default_rng(1090 + args.seed)
    n_sim = 20 if run["quick"] else 200
    d, sims = delorey(args, run, rng, n_sim)
    d.to_csv(out / "a109c_covid_delorey.csv", index=False)
    sims.to_csv(out / "a109c_covid_sims.csv", index=False)
    print(d.round(3).T.to_string(), flush=True)
    m, msims = melms(args, run, rng, n_sim)
    m.to_csv(out / "a109c_covid_melms.csv", index=False)
    pd.concat([sims, msims]).to_csv(out / "a109c_covid_sims.csv", index=False)
    print(m.round(3).T.to_string(), flush=True)


if __name__ == "__main__":
    main()
