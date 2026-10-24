"""Helpers of A111 / A112 / A129 (feeding and genetic controls): expression QC, timed truth, exact tissue-permutation
statistic and a contrast engine for hosts that `Dataset.contrast` does not cover (a single tissue as host).

Local: not in clocksync. `Dataset.contrast` needs >= 2 anchors, so a companion tissue read on one other tissue (WAT on
lung, lung on WAT; liver on heart) is computed here with the same internals (`clocksync._core.host.on_host`) and the
same permutation sequence as the public contrast; `check_against_public` asserts that both give the same null for the
liver.
"""

from __future__ import annotations

import itertools

import numpy as np
import pandas as pd

import clocksync as cs
from clocksync import stats as S
from clocksync._core import factor
from clocksync._core import host as A


HOURS = 24.0 / (2 * np.pi)
W = 2 * np.pi / 24
EXPR_MIN = 5  # median count of a clock gene in a tissue (A48 / WP3 note)
WF = cs.gtex()
GENES = list(WF.genes)


# ---- expression QC ---------------------------------------------------------------------------------------------------

def gene_qc(pb: pd.DataFrame, dataset: str) -> pd.DataFrame:
    """Median count of each clock gene per compartment; `kept` when >= 5 (the registered rule)."""
    x = pb[GENES]
    med = x.fillna(0).groupby(pb.compartment).median()  # a gene absent from a library (removed by the authors) counts 0
    q = med.stack().rename("median_count").reset_index()
    q.columns = ["compartment", "gene", "median_count"]
    q["n_missing"] = x.isna().groupby(pb.compartment).sum().stack().to_numpy()
    q["dataset"], q["kept"] = dataset, q.median_count >= EXPR_MIN
    return q[["dataset", "compartment", "gene", "median_count", "n_missing", "kept"]]


def waveform_without(genes) -> cs.Waveform:
    return WF.subset([g for g in GENES if g not in set(genes)])


# ---- circular helpers ------------------------------------------------------------------------------------------------

def circ_centre(h) -> np.ndarray:
    """Hours (wrapped) minus their circular mean."""
    h = np.asarray(h, float)
    return S.wrap(h - S.circmean(h))


def rel_to_rest(delta: pd.Series, k: str) -> float:
    """Shift of tissue k (h) relative to the circular mean of the other tissues' shifts."""
    others = delta.drop(k)
    return float(S.wrap(delta[k] - S.circmean(others.to_numpy())))


def exact_circ_corr(a_h, b_h) -> dict:
    """Fisher-Lee circular correlation of two series over tissues (hours in) and its exact one-sided permutation p over
    all n! tissue permutations (the observed pairing is one of them)."""
    a, b = np.asarray(a_h, float) * W, np.asarray(b_h, float) * W
    obs = S.fisher_lee(a, b)
    null = np.array([S.fisher_lee(a, b[list(p)]) for p in itertools.permutations(range(len(b)))])
    hit = int(np.sum(null >= obs - 1e-9))
    return dict(n=len(a), r_fl=obs, p_exact=hit / len(null), n_perm_total=len(null), n_ge=hit, p_min=1 / len(null))


# ---- timed truth ----------------------------------------------------------------------------------------------------

def timed_truth(ds: cs.Dataset, hours: str, n_boot: int, seed: int = 0) -> pd.DataFrame:
    """`Dataset.timed_loadings` plus donor-bootstrap draws of the clock offset (n_boot rows per compartment, long)."""
    t = ds.timed_loadings(hours)
    return t


def boot_timed(ds: cs.Dataset, hours: str, n_boot: int, seed: int) -> list:
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(n_boot):
        try:
            out.append(ds.resample(rng).timed_loadings(hours))
        except (ValueError, np.linalg.LinAlgError):
            continue
    return out


# ---- contrast engine -------------------------------------------------------------------------------------------------

def pseudo_host(d: A.Data, tissue: str) -> factor.FactorFit:
    """A single tissue's own clock-plane states as the host (measurement covariance = its Poisson V; the tissue's own
    biological noise stays in the host, so couplings are attenuated alike in every group)."""
    i = np.flatnonzero(d.compartment == tissue)
    return factor.FactorFit(np.array([tissue]), d.resection[i], tissue, np.array([1 + 0j]), np.zeros((1, 2)),
                            np.array([0.0]), np.zeros((2, 2)), d.beta[i], d.V[i], np.nan, 0, True, None)


def multi_host(d: A.Data, anchors, target: str) -> factor.FactorFit:
    """Host of the anchors (as `Dataset.contrast` does for a non-anchor target)."""
    return A._hosts(d, [a for a in anchors if np.any(d.compartment == a)], [target])[target]


def _grp(d, k, h, V0, dn, min_donors):
    if len(dn) < min_donors:
        return dict(c=np.nan, cs=np.nan, z=np.nan + 0j, boundary=False)
    e = A.on_host(d, k, h, dn)
    s = e["z2"] * V0
    return dict(c=e["coupling"], cs=s / (s + 2 * e["psi"]) if np.isfinite(s) else np.nan, z=e["z"],
                boundary=e["boundary"])


def contrast_engine(d: A.Data, keys: dict, lab: pd.Series, blk, levels: tuple, n_perm: int, n_boot: int, seed: int = 0,
                    min_donors: int = 5) -> dict:
    """Coupling contrast (case minus reference) of several (target, host) pairs on one set of permutations and
    bootstrap draws (identical to `Dataset.contrast`'s sequence: `rng = default_rng(seed)`, `permute_within(lab, blk)`).

    keys : name -> (target tissue, host FactorFit). d : Data restricted to the labelled donors.
    Returns dict with `obs` (DataFrame per key: c_ref, c_case, d, d_star, offset_diff_h, n_ref, n_case), `null`
    (name -> DataFrame of draws: d, d_star, off), `boot` (name -> DataFrame of replicate d, d_star).
    """
    prep = {}
    for name, (k, h) in keys.items():
        idx = pd.Index(np.intersect1d(d.donors(k), h.resections))
        V0 = A.on_host(d, k, h, idx.to_numpy())["V"]
        prep[name] = (k, h, idx, V0)

    def diffs(labels):
        out = {}
        for name, (k, h, idx, V0) in prep.items():
            lk = labels.reindex(idx).dropna()
            r = [_grp(d, k, h, V0, lk.index[lk == g].to_numpy(), min_donors) for g in levels]
            both = r[0]["boundary"] and r[1]["boundary"]
            off = S.wrap(np.angle(r[1]["z"] / r[0]["z"]) * HOURS) if abs(r[0]["z"]) > 0 else np.nan
            out[name] = (np.nan if both else r[1]["c"] - r[0]["c"], np.nan if both else r[1]["cs"] - r[0]["cs"],
                         float(off), r)
        return out

    obs = diffs(lab)
    rng = np.random.default_rng(seed)
    null = {n: [] for n in prep}
    for _ in range(n_perm):
        r = diffs(S.permute_within(lab, blk, rng))
        for n in prep:
            null[n].append(r[n][:3])
    null = {n: pd.DataFrame(v, columns=["d", "d_star", "off"], dtype=float) for n, v in null.items()}
    # bootstrap within groups, one draw shared by all keys (paired intervals of differences between keys)
    brng = np.random.default_rng(seed + 1)
    rows = {n: {r: np.flatnonzero((d.compartment == k) & (d.resection == r)) for r in prep[n][2]}
            for n, (k, *_rest) in prep.items()}
    boot = {n: [] for n in prep}
    common = prep[next(iter(prep))][2]
    for n_ in prep:
        common = common.intersection(prep[n_][2])
    don = {g: np.array([x for x in lab.index[lab == g] if x in common]) for g in levels}
    for _ in range(n_boot):
        draw = {g: brng.choice(don[g], len(don[g])) for g in levels}
        for n, (k, h, idx, V0) in prep.items():
            e = [A.on_host_drawn(d, k, h, draw[g], rows[n]) for g in levels]
            s = [x["z2"] * V0 for x in e]
            cs_ = [si / (si + 2 * x["psi"]) if np.isfinite(si) else np.nan for si, x in zip(s, e)]
            boot[n].append((e[1]["coupling"] - e[0]["coupling"], cs_[1] - cs_[0]))
    boot = {n: pd.DataFrame(v, columns=["d", "d_star"], dtype=float) for n, v in boot.items()}
    tab = pd.DataFrame({n: dict(target=prep[n][0], n_ref=int((lab.reindex(prep[n][2]) == levels[0]).sum()),
                                n_case=int((lab.reindex(prep[n][2]) == levels[1]).sum()),
                                c_ref=obs[n][3][0]["c"], c_case=obs[n][3][1]["c"], d=obs[n][0], d_star=obs[n][1],
                                offset_diff_h=obs[n][2], V0=prep[n][3]) for n in prep}).T
    return dict(obs=tab, null=null, boot=boot, levels=levels)


def summarise_engine(res: dict, companions: dict | None = None, target: str | None = None) -> pd.DataFrame:
    """One row per key: d and d* with bootstrap 95 % intervals, one-sided (less) and two-sided permutation p; plus, for
    `target`, the difference to each companion key (d* and d) on the same permutations and draws."""
    rows = []
    obs, null, boot = res["obs"], res["null"], res["boot"]
    for n in obs.index:
        r = dict(key=n, **{c: obs.loc[n, c] for c in obs.columns})
        for col in ("d", "d_star"):
            b = boot[n][col].dropna()
            r[f"{col}_lo"], r[f"{col}_hi"] = (np.quantile(b, 0.025), np.quantile(b, 0.975)) if len(b) > 10 else (np.nan,) * 2
            r[f"p_less_{col}"] = S.perm_p(null[n][col], float(obs.loc[n, col]), "less")
            r[f"p_two_{col}"] = S.perm_p(null[n][col], float(obs.loc[n, col]))
        r["p_offset"] = S.perm_p(null[n]["off"], float(obs.loc[n, "offset_diff_h"]), tie=0.0)
        rows.append(r)
    out = pd.DataFrame(rows)
    if target is not None:
        for comp in (companions or {}):
            for col in ("d", "d_star"):
                diff = float(obs.loc[target, col] - obs.loc[comp, col])
                nd = null[target][col] - null[comp][col]
                bd = (boot[target][col] - boot[comp][col]).dropna()
                out.loc[out.key == target, f"{col}_minus_{comp}"] = diff
                out.loc[out.key == target, f"{col}_minus_{comp}_lo"] = np.quantile(bd, 0.025) if len(bd) > 10 else np.nan
                out.loc[out.key == target, f"{col}_minus_{comp}_hi"] = np.quantile(bd, 0.975) if len(bd) > 10 else np.nan
                out.loc[out.key == target, f"p_less_{col}_minus_{comp}"] = S.perm_p(nd, diff, "less")
    return out


def pair_corr(d: A.Data, a: str, b: str, donors) -> complex:
    """Poisson-corrected complex correlation of two tissues' states over `donors` (arg: a minus b)."""
    ia = pd.Series(np.flatnonzero(d.compartment == a), index=d.resection[d.compartment == a])
    ib = pd.Series(np.flatnonzero(d.compartment == b), index=d.resection[d.compartment == b])
    sh = ia.index.intersection(ib.index).intersection(pd.Index(donors))
    if len(sh) < 5:
        return np.nan + 0j
    x, y = ia[sh].to_numpy(), ib[sh].to_numpy()
    return S.ccorr(d.beta[x], d.beta[y], np.trace(d.V[x], axis1=1, axis2=2), np.trace(d.V[y], axis1=1, axis2=2))
