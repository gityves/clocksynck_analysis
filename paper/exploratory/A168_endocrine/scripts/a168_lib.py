"""A168 shared definitions: loaders of the cached series tables, the paired designs and their permutation nulls, pseudo-clock sets
for a variable gene list. Coherence, Waveform and the A97 scores are those of A127 / A166 (a127_lib, a166_lib)."""
import sys
import zlib
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "paper"))
sys.path.insert(0, str(ROOT / "paper/exploratory/A166_metabric/scripts"))
import a166_lib as A166  # noqa: E402  (coh, wf_for, GENES, WF, sh, pseudo_sets, pool_genes)
from clocksync._core.controls import CLOCK_LIKE  # noqa: E402
from clocksync.genes import canonical  # noqa: E402
from clocksync.waveform import Waveform  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "output"
GENES, WF, coh, sh = A166.GENES, A166.WF, A166.coh, A166.sh
LN2 = np.log(2)
SEED = zlib.crc32(b"A168") % (2 ** 31)
SERIES = {"poetic": "GSE105777", "letrozole": "GSE20181", "z1031b": "GSE87411"}


def hallmark(name: str) -> list[str]:
    """A Hallmark gene set from the same MSigDB file as `sh.hallmark_g2m` (clock genes removed)."""
    for line in sh.GMT.read_text().splitlines():
        f = line.split("\t")
        if f[0] == name:
            return [g for g in f[2:] if g not in GENES]
    raise KeyError(name)


def estrogen_genes() -> list[str]:
    """HALLMARK_ESTROGEN_RESPONSE_EARLY and _LATE (union)."""
    return sorted(set(hallmark("HALLMARK_ESTROGEN_RESPONSE_EARLY")) | set(hallmark("HALLMARK_ESTROGEN_RESPONSE_LATE")))


def load(series: str, out: Path = OUT):
    """(meta, expression samples x genes float32) of 'poetic', 'letrozole' or 'z1031b'. Expression is log2 (ratio for z1031b)."""
    return pd.read_parquet(out / f"{series}_meta.parquet"), pd.read_parquet(out / f"{series}_expr.parquet")


def clock_genes(X: pd.DataFrame) -> list[str]:
    return [g for g in GENES if g in X.columns]


def Yof(X, genes):
    """Y = ln 2 x log2 expression of `genes` (samples x genes)."""
    return LN2 * X[list(genes)].to_numpy(float)


def wf_of(genes) -> Waveform:
    """Waveform of the clock genes present (a subset of the 14 A127 genes, in A127 order)."""
    return WF.subset(genes) if len(genes) < len(GENES) else WF


def wf_pseudo(base: Waveform, genes) -> Waveform:
    """`base` amplitudes / phases carried by the pseudo set `genes` (gene i of the set stands for clock gene i)."""
    return Waveform(tuple(genes), base.amp, base.psi, None, None, "set")


def pool_genes(X, extra_exclude=()):
    """A166 pool (no clock-like, G2M, purity-proxy or clock genes) minus `extra_exclude` (here the oestrogen-response sets)."""
    ex = set(extra_exclude)
    return [g for g in A166.pool_genes(X) if g not in ex]


def pseudo_sets(X, pool, genes, n, rng):
    """A166 / A127 construction for an arbitrary list of clock genes `genes` (the 14 reproduce `a166_lib.pseudo_sets` exactly): each
    gene replaced by one of its 20 nearest pool genes in standardised (mean, SD), no pool gene used twice within a set."""
    pm, ps = X[pool].mean(), X[pool].std()
    zs = lambda a, ref: (a - ref.mean()) / ref.std()
    P = np.column_stack([zs(pm.values, pm), zs(ps.values, ps)])
    target = np.column_stack([zs(X[list(genes)].mean().values, pm), zs(X[list(genes)].std().values, ps)])
    sets = []
    for _ in range(n):
        used, s = set(), []
        for t in rng.permutation(len(genes)):
            near = [j for j in np.argsort(((P - target[t]) ** 2).sum(1))[:40] if j not in used][:20]
            j = int(rng.choice(near)); used.add(j); s.append((t, pool[j]))
        sets.append([g for _, g in sorted(s)])
    return sets


def residual_coh(C):
    """Coherence of Y residualised, within the group, on the covariate matrix C (rows aligned with Y; intercept added)."""
    def f(Y, idx, wf):
        c = np.column_stack([np.ones(len(idx)), C[idx]])
        y = Y[idx]
        return coh(y - c @ np.linalg.lstsq(c, y, rcond=None)[0], np.arange(len(idx)), wf)
    return f


class Pairs:
    """A paired design on integer row positions: one pre and one post sample per patient.
    kind 'did':  T = [coh(post | arm 1) - coh(pre | arm 1)] - [coh(post | arm 0) - coh(pre | arm 0)]; null: arm labels permuted over patients.
    kind 'swap': T = coh(post) - coh(pre); null: the two time labels swapped with probability 1/2 within each patient."""
    def __init__(self, name, kind, pre, post, arm=None):
        self.name, self.kind, self.pre, self.post = name, kind, np.asarray(pre), np.asarray(post)
        self.arm = None if arm is None else np.asarray(arm)

    @property
    def n(self):
        return len(self.pre)

    def draw(self, rng):
        return rng.permutation(self.arm) if self.kind == "did" else rng.random(self.n) < 0.5

    def T(self, Y, wf, cohf=coh, lab=None):
        if self.kind == "did":
            a = self.arm if lab is None else lab
            d = lambda k: cohf(Y, self.post[a == k], wf) - cohf(Y, self.pre[a == k], wf)
            return d(1) - d(0)
        s = np.zeros(self.n, bool) if lab is None else lab
        return cohf(Y, np.where(s, self.pre, self.post), wf) - cohf(Y, np.where(s, self.post, self.pre), wf)

    def null(self, Y, wf, n, rng, cohf=coh):
        return np.array([self.T(Y, wf, cohf, self.draw(rng)) for _ in range(n)])

    def counts(self):
        return dict(n_patients=self.n, **({} if self.arm is None else {"n_arm1": int((self.arm == 1).sum()), "n_arm0": int((self.arm == 0).sum())}))


def pairs_from_meta(m, pre, post, arms=None):
    """One (pre, post) row pair per patient of `m` (columns patient, tp, arm; m['pos'] = row position). arms: arm values kept.
    Returns the position arrays (pre, post, arm values) of patients having exactly one row at each time label."""
    d = m if arms is None else m[m.arm.isin(arms)]
    a = d[d.tp == pre].set_index("patient"); b = d[d.tp == post].set_index("patient")
    assert a.index.is_unique and b.index.is_unique, "more than one row per patient and time label"
    both = a.index.intersection(b.index).sort_values()
    return a.loc[both, "pos"].to_numpy(), b.loc[both, "pos"].to_numpy(), a.loc[both, "arm"].to_numpy()


def build(series, m):
    """The registered design(s) of a series on its meta table m (m['pos'] = row position): dict name -> Pairs."""
    if series == "poetic":
        pre, post, arm = [], [], []
        for (p0, p1, a) in (("B", "S", "AI"), ("D", "S", "UNT")):
            x, y, _ = pairs_from_meta(m[m.arm == a], p0, p1)
            pre.append(x); post.append(y); arm.append(np.full(len(x), int(a == "AI")))
        pre, post, arm = map(np.concatenate, (pre, post, arm))
        ai = arm == 1
        return {"H168.1": Pairs("H168.1", "did", pre, post, arm), "H168.2": Pairs("H168.2", "swap", pre[ai], post[ai])}
    if series == "letrozole":
        x, y, _ = pairs_from_meta(m, "A", "B")
        return {"H168.3": Pairs("H168.3", "swap", x, y)}
    raise KeyError(series)
