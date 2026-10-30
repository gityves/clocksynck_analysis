"""A166 shared definitions: loaders of the cached tidy tables and the three test statistics (A127 definitions)."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "paper"))
sys.path.insert(0, str(ROOT / "paper/exploratory/A127_scanb/scripts"))
import a127_lib as A  # noqa: E402  (coh, wf_for, tertile, GENES, WF, sh)
from clocksync._core.controls import CLOCK_LIKE  # noqa: E402
from clocksync.genes import canonical  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "output"
GENES, WF, coh, wf_for, tertile, sh = A.GENES, A.WF, A.coh, A.wf_for, A.tertile, A.sh
SUBTYPES = A.SUBTYPES
LN2, MIN_ARM, MIN_STRATUM, FLOOR = np.log(2), 10, 30, 4.56


def load(cohort: str):
    """(meta, expression samples x genes float32) of 'metabric' or 'tcga'. tcga expression is log2 CPM, metabric log2 intensity."""
    return pd.read_parquet(OUT / f"{cohort}_meta.parquet"), pd.read_parquet(OUT / f"{cohort}_expr.parquet")


def pool_genes(X: pd.DataFrame):
    g2m = set(sh.hallmark_g2m())
    return [g for g in X.columns if canonical(g) not in CLOCK_LIKE and g not in g2m and g not in sh.PURITY and g not in GENES]


def pseudo_sets(X, pool, n, rng):
    """n pseudo-clock sets of 14 genes, each clock gene replaced by one of its 20 nearest pool genes in standardised (mean, SD)
    (A127 construction)."""
    pm, ps = X[pool].mean(), X[pool].std()
    zs = lambda a, ref: (a - ref.mean()) / ref.std()
    P = np.column_stack([zs(pm.values, pm), zs(ps.values, ps)])
    target = np.column_stack([zs(X[GENES].mean().values, pm), zs(X[GENES].std().values, ps)])
    sets = []
    for _ in range(n):
        used, s = set(), []
        for t in rng.permutation(len(GENES)):
            near = [j for j in np.argsort(((P - target[t]) ** 2).sum(1))[:40] if j not in used][:20]
            j = int(rng.choice(near)); used.add(j); s.append((t, pool[j]))
        sets.append([g for _, g in sorted(s)])
    return sets


# ---- groups (lists of integer index arrays into the full cohort table; m["pos"] = row position) --------------------------------
def er_groups(m, keys, ercol="er", min_arm=MIN_ARM):
    """Groups of tumours sharing `keys`, with >= min_arm tumours per ER arm (`ercol` 0/1, NaN rows must be removed beforehand)."""
    out = []
    for _, g in m.groupby(keys, observed=True):
        e = g[ercol].to_numpy()
        if min((e == 0).sum(), (e == 1).sum()) >= min_arm:
            out.append(g.pos.to_numpy())
    return out


def prolif_strata(m, keys, min_n=MIN_STRATUM):
    """(index arrays, proliferation tertile labels by rank within the stratum) of strata with >= min_n tumours."""
    S, lab = [], []
    for _, g in m.groupby(keys, observed=True):
        if len(g) >= min_n:
            S.append(g.pos.to_numpy()); lab.append(tertile(g.prolif.to_numpy()))
    return S, lab


def T_er(Y, wf, G, er, perms=None):
    """mean over groups of coherence(ER-) - coherence(ER+); perms: one permutation of each group's ER labels. NaN if no group."""
    d = []
    for k, i in enumerate(G):
        l = er[i] if perms is None else er[i][perms[k]]
        d.append(coh(Y, i[l == 0], wf) - coh(Y, i[l == 1], wf))
    return (float(np.mean(d)) if d else np.nan), d


def T_prolif(Y, wf, S, labs, perms=None):
    """mean over strata of coherence(high proliferation tertile) - coherence(low); perms permute the tertile labels within stratum."""
    d = []
    for k, (i, l) in enumerate(zip(S, labs)):
        l = l if perms is None else l[perms[k]]
        d.append(coh(Y, i[l == 2], wf) - coh(Y, i[l == 0], wf))
    return (float(np.mean(d)) if d else np.nan), d


class Design:
    """One test statistic: kind 'er' (groups G, labels er) or 'prolif' (strata G, tertile labels lab)."""
    def __init__(self, name, kind, G, er=None, lab=None):
        self.name, self.kind, self.G, self.er, self.lab = name, kind, G, er, lab

    def T(self, Y, wf, perms=None):
        return T_er(Y, wf, self.G, self.er, perms) if self.kind == "er" else T_prolif(Y, wf, self.G, self.lab, perms)

    def null(self, Y, wf, n, rng):
        return np.array([self.T(Y, wf, [rng.permutation(len(i)) for i in self.G])[0] for _ in range(n)])

    @property
    def n_tumours(self):
        return int(sum(len(i) for i in self.G))


def build(m, ercol="er", pamcol="pam50", by_cohort=True):
    """The three A166 designs on tumour table m (columns cohort, pam50-like, er-like, prolif, purity, pos).
    H1: strata COHORT x purity tertile (tertile within COHORT; cohort-wide when by_cohort is False), >= 10 per ER arm.
    H2: cells PAM50 x proliferation tertile (tertile within PAM50), >= 10 per ER arm.
    H3: strata PAM50 x purity tertile (within PAM50), >= 30 tumours; proliferation tertile within stratum."""
    m = m.copy()
    a = m[m[ercol].notna()].copy()
    a["pt"] = (a.groupby("cohort").purity.transform(tertile) if by_cohort else tertile(a.purity.to_numpy()))
    h1 = Design("H166.1", "er", er_groups(a, ["cohort", "pt"], ercol), er=m[ercol].to_numpy(float))
    b = a[a[pamcol].notna()].copy()
    b["ft"] = b.groupby(pamcol).prolif.transform(tertile)
    h2 = Design("H166.2", "er", er_groups(b, [pamcol, "ft"], ercol), er=m[ercol].to_numpy(float))
    c = m[m[pamcol].notna()].copy()
    c["pt"] = c.groupby(pamcol).purity.transform(tertile)
    S, lab = prolif_strata(c, [pamcol, "pt"])
    return [h1, h2, Design("H166.3", "prolif", S, lab=lab)]
