"""Loaders and small helpers shared by the scripts of lung_breaks.

Method calls are clocksync's. Everything marked "local: not in clocksync" is a short helper for something the package
has no public call for (listed in the folder README).
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import anchors as AN  # noqa: E402
from common import data  # noqa: E402
from common.paths import DATA  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))  # paper/
from common.paths import DATASETS  # noqa: E402
TARGETS = DATASETS / "_atlas" / "a17_targets.csv"
GROUPS = ("control", "case")


def targets() -> pd.DataFrame:
    """The A17 registry (`data/_atlas/a17_targets.csv`), indexed by name."""
    return pd.read_csv(TARGETS).set_index("name", drop=False)


def anchor_set(name: str, rule: str = "learned"):
    """Anchor compartments of an A17 target (None when the rule gives fewer than 3)."""
    t = pd.read_csv(AN.TABLE)
    t = t[t.target == name]
    if t.empty:
        return None
    anc = AN.anchors(name, rule)
    if len(anc) < 3 or (rule == "learned" and str(t.status.iloc[0]).startswith("no anchor")):
        return None
    return anc


def subsample_labels(lab: pd.Series, n: int | None, seed: int) -> pd.Series:
    """`n` random labelled donors in total, every group kept (at least 5 per group when available)."""
    if n is None or n >= len(lab):
        return lab
    rng = np.random.default_rng(seed)
    keep = []
    for g in lab.unique():
        idx = lab.index[lab == g].to_numpy()
        keep += list(rng.permutation(idx)[: max(n // lab.nunique(), 5)])
    return lab.loc[keep]


def braak_stage(index) -> pd.Series:
    """local: not in clocksync. Braak Lewy-body stage groups "0", "1-4", "5-6" of the Parkinson donors (donor ids
    `parkinson_regions:<id>`), from the public donor metadata of the raw atlas (`path_braak_lb`), as legacy `followups._stages`."""
    don = pd.read_csv(DATA / "atlas" / "parkinson_regions" / "donors.csv").set_index("donor_id")
    b = pd.to_numeric(don["path_braak_lb"], errors="coerce")
    b.index = "parkinson_regions:" + b.index.astype(str)
    return pd.cut(b.reindex(index), [-0.5, 0.5, 4.5, 6.5], labels=["0", "1-4", "5-6"]).astype(str)


def load_target(name: str, subsample: int | None = None, seed: int = 0, labels: pd.Series | None = None):
    """(pseudobulk, donor table, labels, blocks, row) of an A17 contrast target. Labels: the donor-table column named
    like the target ("control" / "case") or `labels` (donor -> control / case, or a function of the donor ids); blocks follow the registry's `within`
    ("stratum" -> None, empty -> "none", "block:<col>" -> that column)."""
    row = targets().loc[name]
    pb, donors = data.load(row.dataset)
    lab = (donors[name] if labels is None else labels(donors.index) if callable(labels) else labels).dropna()
    lab = lab[lab.isin(GROUPS)]
    lab = subsample_labels(lab, subsample, seed)
    pb = pb[pb.donor.isin(lab.index)].reset_index(drop=True)
    w = row.within
    blocks = None if w == "stratum" else "none" if not isinstance(w, str) else w.split(":", 1)[1]
    return pb, donors.loc[donors.index.isin(lab.index)], lab, blocks, row


def load_built(name: str, subsample: int | None = None, seed: int = 0, **kw):
    """Standard tables of a builder without a legacy table (e.g. `hlca_tstates`, mode=...): the stored parquet when
    present, else `build(**kw)` of the builder, stored with `data.save` (its options are part of the file name)."""
    from common.paths import PSEUDOBULK
    f = PSEUDOBULK / f"{data._key(name, kw)}.parquet"
    if not f.exists():
        pb, donors = data._builder(name).build(**kw)
        data.save(name, pb, donors, **kw)
    return data.load(name, subsample=subsample, seed=seed, **kw)


def dataset(pb, donors, anchors, lab=None, name="", **kw) -> cs.Dataset:
    """`cs.Dataset` with the labels (donor -> control / case) as donor-table column `lab`."""
    d = donors if lab is None else donors.assign(lab=lab)
    return cs.Dataset(pb, anchors=anchors, donors=d, name=name, **kw)


def block_series(pb: pd.DataFrame, blocks, donors: pd.DataFrame | None, lab: pd.Series):
    """local: not in clocksync. Permutation blocks of `Dataset.contrast` as a Series (None: one global permutation)
    for the label permutations of `module_test`: None -> the donor's stratum when there are several; "none" -> no
    blocks; a donor-table column name -> that column."""
    if isinstance(blocks, str) and blocks == "none":
        return None
    if blocks is None:
        st = pb.groupby("donor").stratum.first() if "stratum" in pb else None
        return None if st is None or st.nunique() < 2 else st.reindex(lab.index).fillna("")
    return donors[blocks].reindex(lab.index).astype(str).fillna("")


def module_test(ds: cs.Dataset, lab: pd.Series, blocks, blk: pd.Series | None, stats: dict, comps, n_perm: int,
                seed: int = 0, jobs: int = 1, scales=("d", "d_star")) -> pd.DataFrame:
    """local: not in clocksync (the package returns a p per compartment, not for functions of several).

    Permutation test of module statistics of the case-minus-control coupling differences. `stats`: name -> function of
    the dict {compartment: d} returning one number. One label permutation is shared by all compartments, within
    `blk` (None: global); the anchored host is that of `Dataset.contrast` on the labelled donors (it does not depend on
    the labels). Two-sided add-one p (|null| >= |obs|); non-finite d are left to the statistic. `scales`: "d" (group
    spread) and "d_star" (common host spread). Returns one row per (scale, statistic) with value, p and the d of every
    compartment."""
    comps = list(comps)

    def table(labels):
        t = ds.contrast(labels, levels=GROUPS, targets=comps, blocks=blocks if blocks is not None else "none",
                        n_perm=0, n_boot=0).table
        return {s: {k: (float(t.loc[k, s]) if k in t.index and t.loc[k].get("status") == "ok" else np.nan)
                    for k in comps} for s in scales}

    obs = table(lab)

    def one(i):
        limit_threads(1)
        rng = np.random.default_rng([seed, i])
        t = table(cs.stats.permute_within(lab, blk, rng))
        return {s: {n: f(t[s]) for n, f in stats.items()} for s in scales}

    null = Parallel(n_jobs=jobs)(delayed(one)(i) for i in range(n_perm))
    rows = []
    for s in scales:
        for n, f in stats.items():
            v = f(obs[s])
            x = np.array([r[s][n] for r in null], float)
            rows.append(dict(scale=s, statistic=n, value=v, p=cs.stats.perm_p(x, v, tie=1e-12),
                             n_perm_ok=int(np.isfinite(x).sum()), **{f"d[{k}]": obs[s][k] for k in comps}))
    return pd.DataFrame(rows)


def mean_of(names):
    """Statistic: mean of the finite entries of `names` (NaN when none)."""
    def f(dk):
        x = [dk[k] for k in names if np.isfinite(dk.get(k, np.nan))]
        return float(np.mean(x)) if x else np.nan
    return f


def difference(a, b):
    """Statistic: a(dk) - b(dk)."""
    return lambda dk: a(dk) - b(dk)
