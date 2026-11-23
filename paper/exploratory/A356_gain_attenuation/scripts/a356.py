"""A356 (check G2; registered in `paper/docs/registered_experiments.md`): is the GTEx organ gain a post-mortem attenuation?

GTEx v10 fast deaths (A263 R2 table), peripheral tissues of the A263 fit except whole blood. Per tissue t: leave-one-out
host h_d (`Dataset.host`, anchors = the other peripheral tissues); tissue states z_td (log fractions of the 14 genes,
centred within t, one projector); gain of a donor group G: |sum z conj h| / sum |h|^2. Donors split at the tissue median
of the sample's ischaemic time or RIN; Delta = log g(long / high) - log g(short / low); null: labels permuted among the
tissue's donors. Aggregates only (GTEx controlled access): per-tissue rows, no donor values.

  --stage power   permutation nulls (no split gain computed): SD, MDE, decidability -> output/power.csv
  --stage score   (power note dated) H356.1-3, D356.4 -> output/score_{tissues,pooled,across,descriptive}.csv

    python a356.py --stage power|score
"""
from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(PAPER / "exploratory/A263_companions/scripts"))
sys.path.insert(0, str(PAPER / "exploratory/A143_brain_artefact/scripts"))
sys.path.insert(0, str(PAPER / "results/figure_2/scripts/organs"))
import clocksync as cs  # noqa: E402
from clocksync import paired as P  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402

warnings.filterwarnings("ignore")
OUT = HERE.parents[1] / "output"
SEED, N_PERM, N_PERM_T, MDE_K, BOUND = 356, 1000, 10_000, 2.49, 0.10
COVS = ("isch", "rin")
SIGN = {"isch": -1, "rin": 1}  # technical-attenuation prediction for Delta (high minus low covariate)


def inputs(core=False):
    """Per peripheral tissue (no whole blood): donor ids, complex states, host states (leave-one-out), covariates."""
    import _shared as SH
    import compute_brain_artefact_checks as A143
    import score_a263 as A

    wf = A.WF.drop(["PER1", "PER2", "NR1D1"]) if core else A.WF
    pb = A.table("human")
    pb = pb[["donor", "compartment", "depth", *A.WF.genes]].reset_index(drop=True)
    brain, per, excl = A.groups("human", sorted(pb.compartment.unique()))
    att = A143.sample_attributes(pb)
    pb["isch"], pb["rin"] = SH.sample_ischaemic_h(pb), att.SMRIN.to_numpy(float)
    ds = cs.Dataset(pb[["donor", "compartment", "depth", *A.WF.genes]], anchors=per, waveform=wf, min_donors=5)
    L = P.projector(pb[list(wf.genes)], wf, counts=True)
    tissues = [t for t in per if t != "Whole Blood"]
    out = {}
    for t in tissues:
        h = ds.host(anchors=[p for p in per if p != t])
        r = pb[pb.compartment == t].set_index("donor")
        r = r[r.index.isin(h.index)]
        beta, _ = P.states(r, L, wf, depth=r.depth.to_numpy(float))
        beta = beta - beta.mean(0)
        hh = h.loc[r.index]
        out[t] = dict(donor=r.index.to_numpy(), z=beta[:, 0] + 1j * beta[:, 1], h=(hh.x + 1j * hh.y).to_numpy(),
                      isch=r.isch.to_numpy(float), rin=r.rin.to_numpy(float))
    return out


def gain(z, h):
    return float(abs(np.sum(z * np.conj(h))) / np.sum(abs(h) ** 2))


def split(x):
    """High half = above the tissue median (ties to the low half); NaN covariate: excluded."""
    ok = np.isfinite(x)
    return ok, ok & (x > np.nanmedian(x))


def delta(d, cov, hi=None):
    ok, high = split(d[cov]) if hi is None else (np.isfinite(d[cov]), hi)
    lo = ok & ~high
    return float(np.log(gain(d["z"][high], d["h"][high])) - np.log(gain(d["z"][lo], d["h"][lo])))


def perm_null(d, cov, k, ti, n=N_PERM):
    ok, high = split(d[cov])
    idx = np.flatnonzero(ok)
    out = np.empty(n)
    for j in range(n):
        rng = np.random.default_rng([SEED, k, ti, j])
        hi = np.zeros(len(high), bool)
        hi[rng.permutation(idx)[:high.sum()]] = True
        out[j] = delta(d, cov, hi)
    return out


def power():
    D = inputs()
    rows, nulls = [], {}
    for k, cov in enumerate(COVS):
        for ti, (t, d) in enumerate(D.items()):
            nulls[(cov, t)] = perm_null(d, cov, k, ti)
            ok, high = split(d[cov])
            rows.append(dict(covariate=cov, tissue=t, donors=int(ok.sum()), high=int(high.sum()),
                             null_sd=float(np.std(nulls[(cov, t)])),
                             median_high=float(np.median(d[cov][high])), median_low=float(np.median(d[cov][ok & ~high]))))
    t = pd.DataFrame(rows)
    pooled = []
    for cov in COVS:
        s = t[t.covariate == cov]
        w = 1 / s.null_sd.to_numpy() ** 2
        T0 = sum(wi * nulls[(cov, ti)] for wi, ti in zip(w, s.tissue)) / w.sum()
        sd = float(np.std(T0))
        pooled.append(dict(covariate=cov, tissues=len(s), null_sd=sd, mde=MDE_K * sd, bound=BOUND,
                           decidable=MDE_K * sd <= BOUND))
    np.savez(OUT / "perm_nulls.npz", **{f"{c}|{ti}": v for (c, ti), v in nulls.items()})
    # across tissues: power of a one-sided Spearman permutation test at rho -0.5 (bivariate normal ranks, 2000 sims)
    n, rng, hits = len(D), np.random.default_rng([SEED, 9]), 0
    from scipy.stats import spearmanr
    crit = np.quantile([spearmanr(rng.normal(size=n), rng.normal(size=n))[0] for _ in range(5000)], 0.05)
    for _ in range(2000):
        x = rng.normal(size=n)
        y = -0.5 * x + np.sqrt(0.75) * rng.normal(size=n)
        hits += spearmanr(x, y)[0] <= crit
    pooled.append(dict(covariate="across tissues (Spearman, rho -0.5)", tissues=n, null_sd=np.nan, mde=np.nan,
                       bound=np.nan, decidable=np.nan, power=hits / 2000))
    t.to_csv(OUT / "power_tissues.csv", index=False)
    pd.DataFrame(pooled).to_csv(OUT / "power.csv", index=False)
    pd.set_option("display.width", 200)
    print(pd.DataFrame(pooled).round(3).to_string(index=False))
    print(t.groupby("covariate")[["donors", "null_sd", "median_high", "median_low"]].median().round(3))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", choices=["power", "score"], required=True)
    a = ap.parse_args(argv)
    limit_threads(4)
    OUT.mkdir(parents=True, exist_ok=True)
    if a.stage == "power":
        power()
    else:
        raise SystemExit("score stage: written after the power note")


if __name__ == "__main__":
    main()
