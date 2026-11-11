"""Clock readouts in a developmental series (A237 onwards).

In a developmental series age moves every gene, and any programme shared by the compartments of an individual passes
a test built for the clock (A236). The readouts are therefore computed on log fractions from which age (and other
covariates) has been removed gene by gene within each compartment, and each is read beside the same statistic on
expression-matched and on responsive pseudo-clock gene sets (`clocksync.paired.pseudo_waveforms`).

- L1, `coherence`: signed clock-gene coherence C across the individuals of one compartment
  (`clocksync.decompose.signed_coherence`: per-sample offset fitted out, sign kept, so that a shift common to the
  clock genes is not read as a clock), against per-gene permutations of the residuals. The unsigned R of
  `decompose.components` (the A189 statistic) is reported beside it.
- L2, `shared`: `Dataset.gate` across the compartments of an individual (free host); a `stratum` column centres the
  cross-covariances and restricts the pairing permutations within strata (age bins). `pairing_z` (shared clock-plane
  covariance standardised by its pairing null) is the statistic compared with the pseudo sets.

Input: `counts`, samples x genes (symbols; counts, or 2^log2 intensity for arrays), and `obs` with one row per
sample: donor, compartment, depth, optional stratum, covariates. `score` runs a set of L1 and L2 tests in parallel;
`verdicts` applies the registered families (Holm within a family, both companions required).
"""
import argparse
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd

import clocksync as cs
from clocksync.decompose import components, signed_coherence
from clocksync.paired import pseudo_waveforms

DEPTH = 1e9  # common depth of the adjusted tables


def basis(age, df: int) -> np.ndarray:
    """Centred polynomial of age (df columns); empty when df = 0 or age is constant."""
    a = np.asarray(age, float)
    if df == 0 or np.ptp(a) == 0:
        return np.empty((len(a), 0))
    a = (a - a.mean()) / a.std()
    return np.column_stack([a ** k for k in range(1, df + 1)])


def design(obs: pd.DataFrame, age: str | None = "age", df: int = 2, factors=()) -> np.ndarray:
    """Covariate matrix of the rows of `obs`: polynomial of age and one indicator per extra level of each factor;
    constant columns are dropped."""
    cols = [basis(obs[age], df)] if age else []
    for f in factors:
        cols.append(pd.get_dummies(obs[f].astype(str), drop_first=True).to_numpy(float).reshape(len(obs), -1))
    Z = np.column_stack(cols) if cols else np.empty((len(obs), 0))
    return Z[:, np.ptp(Z, 0) > 0] if Z.size else Z


def residual(y, Z) -> np.ndarray:
    """y minus its regression on Z (intercept included; the mean is kept)."""
    M = np.column_stack([np.ones(len(y)), Z])
    return y - M @ np.linalg.lstsq(M, y, rcond=None)[0] + y.mean(0)


def table(counts: pd.DataFrame, obs: pd.DataFrame, genes) -> pd.DataFrame:
    """Long table of `genes` for the rows of `obs` (the columns of `obs` first)."""
    t = obs.reset_index(drop=True).copy()
    t[list(genes)] = counts.loc[obs.index, list(genes)].to_numpy(float)
    return t


def adjust(tab: pd.DataFrame, genes, cov) -> tuple[pd.DataFrame, int]:
    """Table whose log fractions are residualised within each compartment on `cov(rows)` (a design matrix), written
    back as counts at a common depth; and the largest number of fitted parameters per gene (intercept included)."""
    genes, out, k = list(genes), [], 1
    for _, g in tab.groupby("compartment", sort=False):
        g = g.copy()
        Z = cov(g)
        y = np.log((g[genes].to_numpy(float) + 0.5) / g.depth.to_numpy(float)[:, None])
        g[genes] = np.exp(residual(y, Z)) * DEPTH
        g["depth"] = DEPTH
        out.append(g)
        k = max(k, 1 + Z.shape[1])
    return pd.concat(out, ignore_index=True), k


def _y(tab, wf):
    y = np.log(tab[list(wf.genes)].to_numpy(float) / tab.depth.to_numpy(float)[:, None])
    return y - y.mean(0)


def _permuted(y, groups, rng):
    """Every column of y permuted on its own, within the row groups (index arrays)."""
    out = y.copy()
    for i in groups:
        for j in range(y.shape[1]):
            out[i, j] = y[rng.permutation(i), j]
    return out


def coherence(tab: pd.DataFrame, wf, n_perm: int = 0, rng=None, strata: str | None = None) -> tuple[float, float, float]:
    """Signed coherence C of one adjusted compartment, its one-sided p against per-gene permutations of the rows (nan
    when n_perm = 0) and the 95th percentile of that null. `strata`: column of `tab` within whose levels the rows are
    permuted (cohorts of a pooled compartment, which are not exchangeable with each other)."""
    y = _y(tab, wf)
    obs = signed_coherence(y, wf)
    if not n_perm:
        return obs, np.nan, np.nan
    groups = [np.arange(len(y))] if strata is None else list(tab.groupby(strata, sort=False).indices.values())
    null = np.array([signed_coherence(_permuted(y, groups, rng), wf) for _ in range(n_perm)])
    return obs, (1 + np.sum(null >= obs)) / (1 + n_perm), float(np.quantile(null, .95))


def unsigned(tab: pd.DataFrame, wf) -> float:
    """Unsigned coherence R of `decompose.components` (no offset; the A189 statistic)."""
    return components(_y(tab, wf), np.zeros(len(wf.genes)), 1, wf)["coherence"]


def shared(tab: pd.DataFrame, wf, n_perm: int, seed: int, min_donors: int = 5):
    """Gate across the compartments (>= `min_donors` individuals each) of an adjusted table; None when no
    compartment pair shares >= 5 individuals (not estimable)."""
    keep = tab.groupby("compartment").donor.transform("nunique") >= min_donors
    try:
        return cs.Dataset(tab[keep], anchors=None, waveform=wf, min_donors=min_donors, name="dev").gate(
            n_perm=n_perm, seed=seed)
    except ValueError:
        return None


def log_mean(counts: pd.DataFrame, depth, min_count: float = 1.0) -> pd.Series:
    """Log of the mean fraction of every gene with mean count >= `min_count` (the matching variable)."""
    x = counts.to_numpy(float)
    lm = pd.Series(np.log((x / np.asarray(depth, float)[:, None]).mean(0) + 1e-12), index=counts.columns)
    return lm[x.mean(0) >= min_count]


def surrogates(counts: pd.DataFrame, depth, Z=None, n: int = 1, exclude=(), min_count: float = 50.0) -> np.ndarray:
    """Surrogates of an unrecorded covariate (developmental age, composition): the first `n` principal components
    (samples x n, unit variance) of the log fractions of the genes with mean count >= `min_count`, `exclude` (the
    clock genes) left out, after the regression of every gene on Z."""
    x = counts.drop(columns=[g for g in exclude if g in counts]).to_numpy(float)
    x = x[:, np.nan_to_num(x.mean(0)) >= min_count]
    y = np.log((x + 0.5) / np.asarray(depth, float)[:, None])
    y = residual(y, np.empty((len(y), 0)) if Z is None else Z)
    u = np.linalg.svd(y - y.mean(0), full_matrices=False)[0][:, :n]
    return u / u.std(0)


def l1_test(counts, obs, cov, wf, n_perm: int, n_pseudo: int, seed, strata: str | None = None) -> dict:
    """L1 of one compartment (the rows of `obs`) with its companions and the leave-one-gene-out range.
    `seed`: list of ints (the three streams are [.., 0] permutations, [.., 1] matched, [.., 2] responsive);
    `strata`: column of `obs` within whose levels the permutations are done (`coherence`)."""
    genes = list(wf.genes)
    tab, k = adjust(table(counts, obs, genes), genes, cov)
    r, p, q95 = coherence(tab, wf, n_perm, np.random.default_rng([*seed, 0]), strata)
    row = dict(n=len(obs), k=k, C=r, p=p, null_q95=q95, R_unsigned=unsigned(tab, wf))
    lm = log_mean(counts.loc[obs.index], obs.depth)
    for j, kind in enumerate(("matched", "responsive"), 1):
        rs = []
        for w in pseudo_waveforms(lm, n_pseudo, np.random.default_rng([*seed, j]), wf, kind):
            t, _ = adjust(table(counts, obs, w.genes), w.genes, cov)
            rs.append(coherence(t, w)[0])
        row[f"{kind}_q95"], row[f"{kind}_pct"] = np.quantile(rs, .95), 100 * np.mean(np.array(rs) < r)
    logo = {}
    for g in genes:
        w = wf.drop([g])
        logo[g] = coherence(adjust(table(counts, obs, w.genes), w.genes, cov)[0], w)[0]
    row["logo_min"], row["logo_max"], row["logo_min_gene"] = min(logo.values()), max(logo.values()), min(logo, key=logo.get)
    row["clock_specific"] = bool(r > row["matched_q95"] and r > row["responsive_q95"])
    return row


def l2_prepare(counts, obs, cov, wf, n_pseudo: int, seed) -> tuple[pd.DataFrame, dict]:
    """Adjusted clock table and the pseudo waveforms of an L2 test (matched, responsive); `seed` as in `l1_test`."""
    genes = list(wf.genes)
    tab, _ = adjust(table(counts, obs, genes), genes, cov)
    lm = log_mean(counts.loc[obs.index], obs.depth)
    ws = {kind: pseudo_waveforms(lm, n_pseudo, np.random.default_rng([*seed, j]), wf, kind)
          for j, kind in enumerate(("matched", "responsive"), 1)}
    return tab, ws


def l2_pseudo(counts, obs, cov, w, n_perm: int, seed: int) -> dict:
    """Gate of one pseudo waveform on the same individuals, covariates, strata and pairing permutations (same seed)."""
    t, _ = adjust(table(counts, obs, w.genes), list(w.genes), cov)
    g = shared(t, w, n_perm, seed)
    return dict(pairing_z=np.nan, p_pairing=np.nan) if g is None else dict(pairing_z=g.pairing_z, p_pairing=g.p_pairing)


def l2_row(g, pseudo: dict) -> dict:
    """Result row of an L2 test from its gate and the pseudo-set gates ({kind: DataFrame of `l2_pseudo` rows})."""
    if g is None:
        return dict(estimable=False)
    row = dict(estimable=True, n=g.n_donors, n_compartments=g.n_compartments, n_pairs=g.n_pairs,
               clock_share=g.clock_share, share_null_mean=g.share_null_mean, p_pairing=g.p_pairing,
               p_acrophase=g.p_acrophase, pairing_z=g.pairing_z, gate=bool(g.passes))
    for kind, t in pseudo.items():
        row[f"{kind}_z_q95"] = t.pairing_z.quantile(.95)
        row[f"{kind}_pct"] = 100 * (t.pairing_z < g.pairing_z).mean()
        row[f"{kind}_share_p05"] = (t.p_pairing < .05).mean()
    row["clock_specific"] = bool(all(g.pairing_z > row[f"{kind}_z_q95"] for kind in pseudo))
    return row


_S = {}  # state of `score`, inherited by its workers (fork)


def _task(t):
    kind, name = t[:2]
    counts, obs, wf, seed = _S["counts"], _S["obs"], _S["wf"], _S["seed"]
    if kind == "L1":
        rows, cov, *strata = _S["l1"][name]
        return t, l1_test(counts, obs[rows], cov, wf, _S["n_perm"], _S["n_pseudo"],
                          [seed, 1, list(_S["l1"]).index(name)], *strata)
    rows, cov = _S["l2"][name]
    if kind == "L2":
        return t, shared(_S["tab"][name], wf, _S["n_perm"], seed)
    return t, l2_pseudo(counts, obs[rows], cov, _S["ws"][name][t[2]][t[3]], _S["n_perm_pseudo"], seed)


def score(counts, obs, l1: dict, l2: dict, wf, seed: int, n_perm: int = 10000, n_pseudo: int = 200,
          n_perm_pseudo: int = 1000, jobs: int = 4) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Run L1 tests ({name: (rows, cov)} or (rows, cov, strata), rows a boolean mask of `obs` selecting one
    compartment, strata a column of `obs` within which the permutations are done) and L2 tests
    ({name: (rows, cov)}, several compartments; `obs.stratum` gives the strata) with their companions.
    Returns the L1 table (`l1_test` rows), the L2 table (`l2_row` rows) and the pseudo-set gates of the L2 tests."""
    from clocksync.config import limit_threads
    _S.update(counts=counts, obs=obs, l1=l1, l2=l2, wf=wf, seed=seed, n_perm=n_perm, n_pseudo=n_pseudo,
              n_perm_pseudo=n_perm_pseudo, tab={}, ws={})
    tasks = [("L1", name) for name in l1]
    for j, (name, (rows, cov)) in enumerate(l2.items()):
        _S["tab"][name], _S["ws"][name] = l2_prepare(counts, obs[rows], cov, wf, n_pseudo, [seed, 2, j])
        tasks += [("L2", name)] + [("pseudo", name, kind, i) for kind in _S["ws"][name] for i in range(n_pseudo)]
    t0, res = time.time(), {}
    with Pool(jobs, initializer=limit_threads, initargs=(1,)) as pool:
        for k, (t, r) in enumerate(pool.imap_unordered(_task, tasks)):
            res[t] = r
            if (k + 1) % 50 == 0 or k + 1 == len(tasks):
                print(f"{k + 1}/{len(tasks)} {time.time() - t0:.0f}s", flush=True)
    t1 = pd.DataFrame([dict(test=name, **res[("L1", name)]) for name in l1])
    pseudo = pd.DataFrame([dict(test=t[1], kind=t[2], set=t[3], **r) for t, r in res.items() if t[0] == "pseudo"])
    t2 = pd.DataFrame([dict(test=name, **l2_row(res[("L2", name)], {
        kind: pseudo[(pseudo.test == name) & (pseudo.kind == kind)] for kind in _S["ws"][name]}
        if res[("L2", name)] is not None else {})) for name in l2])
    return t1, t2, pseudo.sort_values(["test", "kind", "set"]) if len(pseudo) else pseudo


def main(doc: str, load, tests, wf, seed: int, out, families: dict | None = None) -> None:
    """Command line of a scoring script: `load() -> (counts, obs)`, `tests(obs) -> (l1, l2)` as in `score`. Writes
    design.csv (always) and l1.csv, l2.csv, l2_pseudo.csv and, with `families`, verdicts.csv to `out` (`_check` with
    --check: 200 permutations, 10 sets)."""
    from clocksync.config import limit_threads
    ap = argparse.ArgumentParser(description=doc.split("\n")[0])
    ap.add_argument("--jobs", type=int, default=12)
    ap.add_argument("--n-perm", type=int, default=10000)
    ap.add_argument("--n-pseudo", type=int, default=200)
    ap.add_argument("--design", action="store_true", help="write the design table only (no clock statistic)")
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    limit_threads(1)
    out.mkdir(exist_ok=True)
    counts, obs = load()
    l1, l2 = tests(obs)
    design = pd.DataFrame([dict(kind=k, test=name, n_samples=int(rows.sum()), n_individuals=obs[rows].donor.nunique(),
                                compartments=obs[rows].compartment.nunique(), n_covariates=cov(obs[rows]).shape[1],
                                strata=obs[rows].stratum.nunique() if k == "L2" else 1)
                           for k, d in (("L1", l1), ("L2", l2)) for name, (rows, cov, *_) in d.items()])
    design.to_csv(out / "design.csv", index=False)
    print(design.to_string(index=False), flush=True)
    if a.design:
        return
    sfx = "_check" if a.check else ""
    n_perm, n_pseudo = (200, 10) if a.check else (a.n_perm, a.n_pseudo)
    t1, t2, pseudo = score(counts, obs, l1, l2, wf, seed, n_perm, n_pseudo, min(n_perm, 1000), a.jobs)
    for t, name in ((t1, "l1"), (t2, "l2"), (pseudo, "l2_pseudo")):
        t.to_csv(out / f"{name}{sfx}.csv", index=False)
    pd.set_option("display.width", 250)
    print(t1.round(4).to_string(index=False), t2.round(4).to_string(index=False), sep="\n", flush=True)
    if families:
        v = verdicts(t1, t2, families)
        v.to_csv(out / f"verdicts{sfx}.csv", index=False)
        print(v.round(4).to_string(index=False), flush=True)


def verdicts(l1: pd.DataFrame, l2: pd.DataFrame, families: dict) -> pd.DataFrame:
    """Registered reading of the tests of `score`: `families` = {hypothesis: [test names]}, Holm within a family.
    The p of an L2 test is the larger of its pairing and acrophase p. PASS: Holm p < 0.05 and clock-specific (above
    the q95 of both pseudo sets); a gate that is not estimable has no verdict."""
    cols = ["test", "level", "n", "stat", "p", "matched_pct", "responsive_pct", "clock_specific"]
    a = l1.assign(level="L1", stat=l1.C) if len(l1) else pd.DataFrame(columns=cols)
    b = l2.reindex(columns=[*l2.columns, *(c for c in ("test", "p_pairing", "p_acrophase", "pairing_z", "n",
                                                        "matched_pct", "responsive_pct", "clock_specific")
                                           if c not in l2)])  # an empty table (no L2 test) has no column
    b = b.assign(level="L2", stat=b.pairing_z, p=b[["p_pairing", "p_acrophase"]].astype(float).max(axis=1, skipna=False))
    t = pd.concat([a[cols], b[cols]]).set_index("test")
    rows = []
    for h, names in families.items():
        x = t.loc[names]
        ph = holm(x.p.astype(float).fillna(1.0))
        for (name, r), q in zip(x.iterrows(), ph):
            ok = bool(r.clock_specific) if pd.notna(r.clock_specific) else False
            verdict = ("not estimable" if pd.isna(r.p) else "PASS" if q < .05 and ok else
                       "not clock-specific" if q < .05 else "fail")
            rows.append(dict(hypothesis=h, test=name, level=r.level, n=r.n, stat=r.stat, p=r.p, p_holm=q,
                             matched_pct=r.matched_pct, responsive_pct=r.responsive_pct, verdict=verdict))
    return pd.DataFrame(rows)


def holm(p) -> np.ndarray:
    """Holm-adjusted p values."""
    p = np.asarray(p, float)
    o = np.argsort(p)
    adj = np.minimum(1, np.maximum.accumulate((len(p) - np.arange(len(p))) * p[o]))
    out = np.empty_like(adj)
    out[o] = adj
    return out
