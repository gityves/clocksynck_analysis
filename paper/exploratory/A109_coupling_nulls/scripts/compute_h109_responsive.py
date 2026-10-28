"""A109 note (2026-10-28, roadmap 5.1): H109.1-R and the PER1/PER2-free sensitivity of H109.1 (GTEx v10, 12 tissues).

Registered in paper/docs/registered_experiments.md, section "A109 / A138 note (2026-10-28 ...)". Secondary readouts: they
do not change the registered H109.1 verdict.
  H109.1-R   ventilated donors: C_clock exceeds the median of 200 responsive-gene sets (glucocorticoid targets and
             immediate-early genes, expression-matched, carrying the GTEx waveform loadings;
             `clocksync.paired.pseudo_waveforms(kind="responsive")`). One-sided add-one p < 0.05.
  Sensitivity  H109.1 with PER1 and PER2 dropped from the waveform (`Waveform.drop`): (a) against the registered
             matched pseudo-clock null (200 sets of 14 genes, as stored), (b) against 200 matched sets of 12 genes
             carrying the reduced waveform (same size as the statistic).
  Reported   mean Jaccard overlap of the responsive sets (per tissue and overall).
The registered statistic is reproduced first (C_clock, matched pseudo median and q95 against output/h109_1_null.csv).

    python compute_h109_responsive.py [--subsample N] [--seed S]
Outputs: paper/exploratory/A109_coupling_nulls/output/responsive/ (subsample/ for checks).
Aggregates only (GTEx controlled access): no per-donor values are written.
"""

from __future__ import annotations

import itertools
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import _helpers as H  # noqa: E402
import compute_coupling_nulls as base  # noqa: E402
from clocksync._core.controls import RESPONSIVE_GENES  # noqa: E402
from common import cli, gtex  # noqa: E402
from common.paths import exploratory_dir  # noqa: E402

FOLDER = "A109_coupling_nulls"
DROP = ("PER1", "PER2")
BODY, GROUPS, WF, P = H.BODY, H.GROUPS, H.WF, H.P
WF12 = WF.drop(DROP)


def reread_rows(T: H.Tissues):
    """Yield (tissue index, tissue, DataFrame of counts of all genes of its rows, aligned with `T.rows[t]`).

    The row selection of `H.Tissues` is repeated (same filters, same order) and checked against its donors and depth."""
    import anndata as ad

    a = ad.read_h5ad(gtex.LABELLED, backed="r")
    try:
        o = a.obs
        keep = (o.DTHHRDY.isin([0, 1, 2, 3, 4]) & o.tissue.astype(str).isin(BODY) & o.tod.notna()).to_numpy()
        sub = o.iloc[np.flatnonzero(keep)]
        subj = sub.subject.astype(str).to_numpy()
        ok = ~pd.Series(list(zip(subj, sub.tissue.astype(str)))).duplicated().to_numpy()
        sub, subj = sub[ok], subj[ok]
        names = np.asarray(a.var_names)
        for k, t in enumerate(BODY):
            m = (sub.tissue.astype(str) == t).to_numpy() & np.isin(subj, T.donors)
            ids = sub.index[m]
            pos = a.obs_names.get_indexer(ids)
            o_ = np.argsort(pos)
            pos, sd = pos[o_], subj[m][o_]
            lib = a.obs.library_size.to_numpy(float)[pos]
            assert (sd == T.row_donor[t]).all() and np.allclose(lib, T.rows[t].depth.to_numpy(float)), t
            X = a.X[pos]
            X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
            lm = pd.Series(np.log(X / lib[:, None] * 1e6 + 0.5).mean(0), index=names)
            yield k, t, X, lm, names
    finally:
        a.file.close()


def add_sets(T: H.Tissues, n_sets: int, seed: int) -> pd.DataFrame:
    """States of the responsive sets (`T.Z["resp"]`), of 12-gene matched sets (`"pm12"`) and of the real clock genes
    without PER1/PER2 (`"clock12"`). Returns the Jaccard overlap and pool size per tissue."""
    for d in (T.Z, T.V):
        d.update(resp={}, pm12={}, clock12={})
    jac = []
    for k, t, X, lm_all, names in reread_rows(T):
        dep = T.rows[t].depth.to_numpy(float)
        lm_r = lm_all[~lm_all.index.duplicated()]
        lm_m = lm_all[~lm_all.index.isin(H.CLOCK_OUTPUT - set(WF.genes))]
        lm_m = lm_m[~lm_m.index.duplicated()]
        col = {g: i for i, g in enumerate(names)}
        sets = {}
        sets["resp"] = P.pseudo_waveforms(lm_r, n_sets, np.random.default_rng([109, 7, seed, k]), WF, kind="responsive")
        sets["pm12"] = P.pseudo_waveforms(lm_m, n_sets, np.random.default_rng([109, 8, seed, k]), WF12)
        for kind, ws in sets.items():
            u = sorted({g for w in ws for g in w.genes})
            x = pd.DataFrame(X[:, [col[g] for g in u]].astype(float), columns=u)
            z, v = H._state_set(x, ws, dep)
            T.Z[kind][t], T.V[kind][t] = z, v
        z, v = H._state_set(T.rows[t], [WF12], dep)
        T.Z["clock12"][t], T.V["clock12"][t] = z, v
        gs = [set(w.genes) for w in sets["resp"]]
        pool = [g for g in RESPONSIVE_GENES if g in col]
        j = [len(a & b) / len(a | b) for a, b in itertools.combinations(gs, 2)]
        jac.append(dict(tissue=t, n_pool_present=len(pool), set_size=len(gs[0]), n_sets=len(gs),
                        mean_jaccard=float(np.mean(j)), min_jaccard=float(np.min(j)), max_jaccard=float(np.max(j)),
                        genes_used=len(set().union(*gs))))
        del X
    return pd.DataFrame(jac)


def row(g, test, kind_stat, kind_null, draws, T, **kw):
    c = H.coherence(T, kind_stat, draws[g])[0]
    nul = H.coherence(T, kind_null, draws[g])
    return dict(group=g, test=test, C_stat=c, null_median=float(np.median(nul)), null_q95=float(np.quantile(nul, 0.95)),
                null_sd=float(nul.std(ddof=1)), excess=float(c - np.median(nul)), p_one_sided=base.pval(c, nul),
                n_sets=len(nul), n_donors=len(T.group_donors(g)), n_pairs=len(H.PAIRS), **kw)


def reproduce(T, draws, out: Path, quick: bool) -> pd.DataFrame:
    """Registered C_clock and matched pseudo-clock null of H109.1 against output/h109_1_null.csv."""
    ref = pd.read_csv(out.parents[0] / ("subsample/" if quick else "") / "h109_1_null.csv")
    ref = ref[ref.null_type == "pseudo_clock"].set_index("group")
    rows = []
    for g in GROUPS:
        c = H.coherence(T, "clock", draws[g])[0]
        nul = H.coherence(T, "pseudo", draws[g])
        new = dict(C_clock=c, null_median=float(np.median(nul)), null_q95=float(np.quantile(nul, 0.95)),
                   p_one_sided=base.pval(c, nul))
        for key, v in new.items():
            rows.append(dict(group=g, quantity=key, reproduced=v, stored=ref.loc[g, key],
                             match=bool(np.isclose(v, ref.loc[g, key], rtol=1e-9, atol=1e-12))))
    r = pd.DataFrame(rows)
    if not quick:
        assert r.match.all(), f"registered H109.1 not reproduced:\n{r[~r.match]}"
    return r


def main(argv=None):
    ap = cli.parser(__doc__)
    args = ap.parse_args(argv)
    run = cli.setup(args)
    quick = run["quick"]
    sz = base.sizes(quick)
    out = exploratory_dir(FOLDER).joinpath("responsive")
    out = out / "subsample" if quick else out
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(109 + args.seed)
    T = H.Tissues(sz["n_sets"], args.seed, args.subsample)
    T.build_states(np.random.default_rng(2109 + args.seed), 1, real=True)  # acrophase sets are not used here
    n_eq = len(T.group_donors("H4"))
    draws = {g: H.draws_equal_n(T.group_donors(g), min(n_eq, len(T.group_donors(g))), sz["n_sub"], rng) for g in GROUPS}
    print({g: len(T.group_donors(g)) for g in GROUPS}, "equal n:", n_eq, flush=True)
    rep = reproduce(T, draws, out.parent if quick else out, quick)
    rep.to_csv(out / "reproduction_check.csv", index=False)
    print(rep.round(6).to_string(index=False), flush=True)
    n_sets = 20 if quick else 200
    jac = add_sets(T, n_sets, args.seed)
    jac.loc[len(jac)] = dict(tissue="mean over tissues", **jac.drop(columns="tissue").mean().to_dict())
    jac.to_csv(out / "jaccard_responsive_sets.csv", index=False)
    print(jac.round(3).to_string(index=False), flush=True)
    rows = []
    for g in GROUPS:
        rows.append(row(g, "H109.1-R: C_clock vs responsive sets", "clock", "resp", draws, T))
        rows.append(row(g, "PER-free C_clock vs registered matched null (14 genes)", "clock12", "pseudo", draws, T))
        rows.append(row(g, "PER-free C_clock vs matched 12-gene sets", "clock12", "pm12", draws, T))
        rows.append(row(g, "registered C_clock vs matched 12-gene sets (companion)", "clock", "pm12", draws, T))
    res = pd.DataFrame(rows)
    res.to_csv(out / "h109_1_responsive.csv", index=False)
    print(res.round(4).to_string(index=False), flush=True)
    print("done", flush=True)


if __name__ == "__main__":
    main()
