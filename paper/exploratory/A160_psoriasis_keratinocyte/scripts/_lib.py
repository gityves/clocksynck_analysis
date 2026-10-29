"""A160 helpers: pseudobulk tables, data rules, the A78 / A81b paired core-minus-cue estimator on `clocksync.paired`
building blocks, the bin-weighted combination, null by random pseudo-pairs and exact label swaps."""
from __future__ import annotations
import sys
from itertools import product
from pathlib import Path
import numpy as np, pandas as pd, anndata as ad

PAPER = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PAPER.parent)); sys.path.insert(0, str(PAPER))
import clocksync as cs
from clocksync import stats as S
from clocksync import paired as PD
from common.paths import DATA, group_out

HOURS = S.HOURS
OUT = Path(__file__).resolve().parents[1] / "output"
H5 = DATA / "psoriasis_sc" / "a160_pseudobulk.h5ad"
WF = cs.gtex()
CORE = list(WF.genes)
BINS = ["kc_cycling", "kc_d1", "kc_d2", "kc_d3", "kc_d4"]
MIN_LIB, MIN_BIN, MIN_PAIRS_BIN, MIN_CUE_FRAC = 200, 50, 8, 0.6
Z = 1.645 + 0.842  # one-sided 0.05, 80 % power
N_BOOT = 100


def cue_waveform() -> cs.Waveform:
    ref = pd.read_csv(group_out("gtex_cues") / "a81b_cue_reference.csv", index_col=0)
    x = ref[ref.set == "cue"]
    return cs.Waveform(tuple(x.index), x.amp.to_numpy(), x.acro_h.to_numpy() * 2 * np.pi / 24, None, None, "cue")


CUE = cue_waveform()


def restrict(wf, genes):
    idx = [wf.genes.index(g) for g in genes]
    return cs.Waveform(tuple(genes), np.asarray(wf.amp)[idx], np.asarray(wf.psi)[idx], None, None, wf.name)


def load():
    """(counts DataFrame rows `lib|group` x genes (core + cue + depth), obs)."""
    A = ad.read_h5ad(H5)
    A.var_names = [{"BMAL1": "ARNTL"}.get(g, g) for g in A.var_names]
    genes = [g for g in dict.fromkeys(CORE + list(CUE.genes)) if g in A.var_names]
    X = pd.DataFrame(A[:, genes].X.toarray(), index=A.obs.index, columns=genes)
    X["depth"] = A.obs.depth.to_numpy(float)
    return X, A.obs.copy()


def detected(x: pd.DataFrame) -> pd.Series:
    """Library rule: all 14 core genes with >= 1 count and >= 60 % of the 40-gene cue set detected."""
    core = (x[CORE] >= 1).all(1)
    cue = (x[[g for g in CUE.genes if g in x.columns]] >= 1).sum(1) / len(CUE.genes) >= MIN_CUE_FRAC
    return core & cue


def pairs_for(X, obs, group, min_cells):
    """Eligible PN/PP pairs of a group: patient-visit with both libraries >= min_cells and passing the detection rule.
    Returns (table with pair id, case_key, ref_key), pool of PN + healthy keys eligible in the group."""
    o = obs[(obs.group.astype(str) == group)].copy()
    for c in ("patient", "visit", "status"):
        o[c] = o[c].astype(str)
    o["ok"] = (o.n_cells >= min_cells) & detected(X.loc[o.index]).to_numpy()
    o["pair"] = o.patient + "|" + o.visit
    ok = o[o.ok]
    case = ok[ok.status == "lesional"].set_index("pair"); ref = ok[ok.status == "nonlesional"].set_index("pair")
    u = sorted(set(case.index) & set(ref.index))
    pr = pd.DataFrame({"pair": u, "case": [case.loc[p, "lib"] + "|" + group for p in u], "ref": [ref.loc[p, "lib"] + "|" + group for p in u]})
    pool = ok[ok.status.isin(["nonlesional", "healthy"])].index.tolist()
    return pr, pool, o


def cue_genes(X, keys) -> list[str]:
    """Cue genes with >= 1 count in every one of the libraries `keys` (blind: callers pass PN + healthy)."""
    g = [x for x in CUE.genes if x in X.columns]
    return [x for x in g if (X.loc[keys, x] >= 1).all()]


def states(X, keys, cuegenes, drop_core=None):
    """{'core': (beta, v), 'cue': (beta, v)} of the rows `keys` (duplicates allowed), one projector per set over all rows."""
    x = X.loc[list(keys)]
    out = {}
    for k, wf in (("core", WF.drop(drop_core) if drop_core else WF), ("cue", restrict(CUE, cuegenes))):
        g = list(wf.genes)
        Lp = PD.projector(x[g], wf, counts=True)
        out[k] = tuple(np.asarray(v) for v in PD.states(x, Lp, wf, depth=x.depth.to_numpy(float)))
    return out


def group_offset(st, idx) -> float:
    """Within-group core-minus-cue offset arg E[c_core c_cue*] (h) of the rows `idx` (states centred in the group)."""
    (bc, vc), (bk, vk) = st["core"], st["cue"]
    return float(np.angle(S.ccorr(bc[idx], bk[idx], vc[idx], vk[idx])) * HOURS)


def diff(st, ia, ib) -> float:
    """Group A minus group B core-cue offset (h, wrapped)."""
    return float(S.wrap(group_offset(st, ia) - group_offset(st, ib)))


def paired_core(st_core, ia, ib) -> float:
    """Paired offset of the core states of A (rows ia) relative to B (rows ib, aligned): arg ccorr (A minus B), h."""
    (b, v) = st_core
    return float(np.angle(PD.ccorr(b[ia], v[ia], b[ib], v[ib])) * HOURS)


def boot_se(st, ia, ib, rng, n_boot=N_BOOT):
    """Bootstrap SD (h) of the A-minus-B offset: libraries resampled within each group (circular spread)."""
    d0 = diff(st, ia, ib)
    d = np.array([diff(st, rng.choice(ia, len(ia)), rng.choice(ib, len(ib))) for _ in range(n_boot)])
    d = S.wrap(d - d0); d = d[np.isfinite(d)]
    return float(np.std(d, ddof=1)) if len(d) > 2 else np.nan


def wmean(d, se):
    """Inverse-variance weighted circular mean (h)."""
    d, se = np.asarray(d), np.asarray(se)
    w = 1 / np.maximum(se, 1e-3) ** 2
    return float(np.angle(np.sum(w * np.exp(1j * d / HOURS))) * HOURS)


def pseudo_groups(pool, n, rng, scheme):
    """Two groups of PN / healthy libraries. 'boot': n libraries each, drawn with replacement, independently (full pair
    count; the pool is smaller than 2n). 'split': the pool split at random into two disjoint groups (sizes P // 2, P - P // 2)."""
    pool = np.asarray(pool)
    if scheme == "boot":
        return list(rng.choice(pool, n)), list(rng.choice(pool, n))
    p = rng.permutation(len(pool)); h = len(pool) // 2
    return list(pool[p[:h]]), list(pool[p[h:]])


def all_flips(n):
    return np.array(list(product([False, True], repeat=n)))
