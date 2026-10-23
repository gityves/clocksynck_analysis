"""Mixing model of A60 / A61 / A69: which classes could carry the bulk psoriatic lesion delay?

beta_bulk = sum_k L diag(w_k) X beta_k: each class enters the bulk clock-plane state through the complex gain of its 2 x 2
map (w_kg = class share of clock gene g in the bulk, L the bulk projector, X the waveform design). The predicted bulk delay
is arg(z_bulk,L / z_bulk,N); composition-only (lesional shares, uninvolved class states), within-class-only (uninvolved
shares, lesional states) and per-class contributions. Local: not in clocksync (the model is specific to this analysis).
"""

from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(Path(__file__).resolve().parents[4]))  # paper/
from common.paths import group_dir  # noqa: E402
sys.path.insert(0, str(group_dir("psoriasis")))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.optimize import nnls  # noqa: E402

import _shared as H  # noqa: E402
import _singlecell as C  # noqa: E402

GENES, WF, HOURS = H.GENES, H.WF, H.HOURS
N_MARK, MIN_CPM = 30, 20.0
BULK = [("processed", "GSE197056"), ("processed", "GSE268958"), ("processed", "GSE121212"),
        ("processed_heldout", "GSE283313")]


def S_allgenes(name: str):
    """All-gene count sums of the psoriasis single-cell builders (`pseudobulk_all`, `pseudobulk_kc_states_all`)."""
    import importlib
    return H.builder("psoriasis_sc").S.allgenes(name)


def fit_z(t: pd.DataFrame, ks) -> dict:
    """Class state z = rel. amplitude x exp(i offset) (complex, reference-free frame) and coupling per compartment."""
    tab = C.fit_states(t[t.compartment.str[2:].isin(ks)].reset_index(drop=True))
    return {c: (tab.rel_amplitude[c] * np.exp(1j * tab.offset_h[c] / HOURS), tab.coupling[c]) for c in tab.index}


def class_cpm(a, groups: dict) -> pd.DataFrame:
    """CPM (genes x groups) of summed counts of the rows of AnnData `a` selected by boolean masks."""
    x = pd.DataFrame(np.asarray(a.X.todense() if hasattr(a.X, "todense") else a.X), columns=a.var_names)
    x = x.rename(columns={"BMAL1": "ARNTL"}).T.groupby(level=0).sum().T  # the studies name ARNTL differently
    dep = a.obs.depth.to_numpy()
    return pd.DataFrame({k: x[m].sum() / dep[m].sum() * 1e6 for k, m in groups.items()})


def signatures(a, ks, whole=("GSE173706", "GSE228421")):
    """Whole-skin class CPM per status ('L', 'N') and pooled; class mRNA share of the classed cells per status."""
    a = a[a.obs.study.isin(whole) & a.obs.compartment.isin(ks)]
    ob = a.obs.reset_index(drop=True)
    comp, stat = ob.compartment.to_numpy(), ob.status.to_numpy()
    sig = {s: class_cpm(a, {k: (comp == k) & (stat == st) for k in ks}) for st, s in C.STATUS.items()}
    sig["pooled"] = class_cpm(a, {k: comp == k for k in ks})
    prop = {s: ob[ob.status == st].groupby("compartment").depth.sum().reindex(ks).pipe(lambda v: v / v.sum())
            for st, s in C.STATUS.items()}
    return sig, prop


def markers(sig: pd.DataFrame, genes) -> list[str]:
    s = sig.loc[sig.index.intersection(genes)]
    out = []
    for k in s.columns:
        other = s.drop(columns=k).max(1)
        lr = np.log((s[k] + 1) / (other + 1))[s[k] >= MIN_CPM]
        out += list(lr.nlargest(N_MARK).index)
    return sorted(set(out))


def bulk() -> pd.DataFrame:
    """Bulk psoriasis count series (GSE197056, GSE268958, GSE121212, GSE283313): `status` and CPM of every gene
    common to the four series."""
    parts = []
    for folder, g in BULK:
        d = pd.read_parquet(H.DATA / "skin_lesional" / folder / f"{g}.parquet")
        parts.append(d[d.disease == "psoriasis"].rename(columns={"BMAL1": "ARNTL"}))
    meta = {"patient", "status", "sample", "dataset", "disease", "platform", "value_type", "depth", "n_replicates"}
    genes = sorted(set.intersection(*[set(p.columns) for p in parts]) - meta)
    return pd.concat([pd.concat([d[["status"]], d[genes].astype(float).div(d.depth.astype(float), axis=0) * 1e6], axis=1)
                      for d in parts], ignore_index=True)


def bulk_projector(b: pd.DataFrame) -> np.ndarray:
    """One projector of the bulk lesional series: gene weights from the median counts at 30 M reads."""
    return H.PD.projector(pd.DataFrame(b[GENES].to_numpy() * 30.0, columns=GENES), WF, counts=True)


def fractions(b: pd.DataFrame, sig: pd.DataFrame, mk, scaled=False) -> dict:
    """Mean mRNA fractions per status by NNLS (scaled, post hoc: markers divided by their largest signature CPM)."""
    S_ = sig.loc[mk].to_numpy()
    sc = S_.max(1) if scaled else np.ones(len(mk))
    S_ = S_ / sc[:, None]
    out = {}
    for st, s in C.STATUS.items():
        y = b.loc[b.status == st, mk].to_numpy() / sc
        m = np.array([nnls(S_, v)[0] for v in y])
        m = m / m.sum(1, keepdims=True)
        out[s] = pd.Series(m.mean(0), index=sig.columns)
    return out


def shares(m: pd.Series, c: pd.DataFrame) -> pd.DataFrame:
    """w_kg (clock genes x classes)."""
    num = c.loc[GENES] * m
    return num.div(num.sum(1), axis=0)


def gains(w: pd.DataFrame, L: np.ndarray) -> pd.Series:
    """Complex gain of each class's map L diag(w_k) X (rotation-invariant part)."""
    X = WF.X
    out = {}
    for k in w.columns:
        M = L @ (w[k].to_numpy()[:, None] * X)
        out[k] = 0.5 * (M[0, 0] + M[1, 1]) + 0.5j * (M[1, 0] - M[0, 1])
    return pd.Series(out)


def predict(z: dict, g: dict, ks) -> dict:
    """Predicted bulk offsets (h): full, composition only, within-class only, per-class contributions."""
    zs = lambda s, k: z.get(f"{s}:{k}", (np.nan, 0))[0]
    ok = [k for k in ks if np.isfinite(zs("L", k)) and np.isfinite(zs("N", k))]
    bulk_ = lambda gs, ss: sum(g[gs][k] * zs(ss[k], k) for k in ok)
    allL, allN = {k: "L" for k in ok}, {k: "N" for k in ok}
    off = lambda num, den: float(np.angle(num / den) * HOURS)
    zN = bulk_("N", allN)
    out = dict(full=off(bulk_("L", allL), zN), composition=off(bulk_("L", allN), zN), within=off(bulk_("N", allL), zN))
    for k in ok:
        keep = dict(allL, **{k: "N"})
        out[f"contrib_{k}"] = float(C.wrap(out["full"] - off(bulk_("L", keep), zN)))
    return out


def summary(z, ks) -> pd.DataFrame:
    rows = []
    for k in ks:
        for s in "LN":
            c = f"{s}:{k}"
            if c in z:
                rows.append(dict(cls=k, status=s, amp=abs(z[c][0]), offset_h=np.angle(z[c][0]) * HOURS, coupling=z[c][1]))
    return pd.DataFrame(rows)
