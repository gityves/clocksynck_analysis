"""Loaders and small helpers shared by the scripts of windows.

Method calls are clocksync's. Everything marked "local: not in clocksync" is for something the package has no public
call for (listed in the folder README); the A25 internal-time estimator is in `_internal_time.py`.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import clocksync as cs  # noqa: E402
import clocksync.timed  # noqa: E402,F401
from common import anchors as AN  # noqa: E402
from common import data, gtex  # noqa: E402


sys.path.insert(0, str(Path(__file__).resolve().parents[4]))  # paper/
from common.paths import DATASETS  # noqa: E402
TARGETS = DATASETS / "_atlas" / "a17_targets.csv"
W = 2 * np.pi / 24


def targets() -> pd.DataFrame:
    return pd.read_csv(TARGETS).set_index("name", drop=False)


def anchor_set(name: str, rule: str):
    """A17 anchors of a target (None when the rule gives fewer than 2 usable)."""
    t = pd.read_csv(AN.TABLE)
    t = t[t.target == name]
    if t.empty:
        return None
    anc = AN.anchors(name, rule)
    if len(anc) < 2 or (rule == "learned" and str(t.status.iloc[0]).startswith("no anchor")):
        return None
    return anc


def window_factor(t_h) -> float:
    """1 - R_W^2 of collection times (h): the share of between-donor variance of the coherent clock part that survives
    a sampling window (R_W = |mean exp(i omega t)|)."""
    return 1.0 - abs(np.mean(np.exp(1j * W * np.asarray(t_h, float)))) ** 2


def gtex_fast(min_donors: int = 30, subsample: int | None = None, seed: int = 0):
    """GTEx fast deaths (Hardy 1-2), tissues with >= `min_donors` donors, 14 clock genes. Controlled access: the time of
    death stays in memory; scripts write aggregates only. Returns (pseudobulk, donors with `tod`)."""
    pb, donors = gtex.load(min_donors=min_donors)
    if subsample is not None:
        keep = np.random.default_rng(seed).permutation(sorted(pb.donor.unique()))[:subsample]
        pb = pb[pb.donor.isin(keep)].reset_index(drop=True)
        donors = donors.loc[donors.index.isin(keep)]
    return pb, donors


def built(name: str, subsample: int | None = None, seed: int = 0, **kw):
    """Standard tables of a builder without a legacy table (`data/<name>/build.py`): the stored parquet when present,
    else `build(**kw)`, stored with `data.save`."""
    from common.paths import PSEUDOBULK
    if not (PSEUDOBULK / f"{data._key(name, kw)}.parquet").exists():
        pb, donors = data._builder(name).build(**kw)
        data.save(name, pb, donors, **kw)
    return data.load(name, subsample=subsample, seed=seed, **kw)


def baboon(subsample: int | None = None, seed: int = 0):
    """(pseudobulk, donors) of the baboon time course (12 animals = ZT0..ZT22, 63 kept tissues; `data/baboon`).
    Donors: `zt` (h)."""
    return built("baboon", subsample, seed)


def timed_tissue_gains(pb: pd.DataFrame, donors: pd.DataFrame, wf=None, min_cpm: float = 5.0) -> pd.DataFrame:
    """local: not in clocksync. Timed counterpart of a tissue's untimed offset (A36): per tissue, complex cosinor of
    log(CPM + 0.5) on ZT of the clock genes with median CPM >= `min_cpm` (`timed.cosinor`), and the complex gain
    g = conj(x) c / conj(x) x against the waveform's loadings x = amp e^{i psi}. The tissue's clock-plane phase (advance)
    is -arg g; `g_phi_h` is arg g in hours (lag)."""
    wf = wf or cs.gtex()
    genes = list(wf.genes)
    zt = donors.zt.astype(float)
    rows = []
    for k, sub in pb.groupby("compartment"):
        cpm = sub[genes].to_numpy(float) / sub.depth.to_numpy(float)[:, None] * 1e6
        ex = np.median(cpm, 0) >= min_cpm
        r = cs.timed.cosinor(np.log(cpm[:, ex] + 0.5), zt.reindex(sub.donor).to_numpy())
        c = r["amplitude"] * np.exp(1j * r["acrophase_h"] * W)
        x = (wf.amp * np.exp(1j * wf.psi))[ex]
        g = (np.conj(x) @ c) / (np.conj(x) @ x)
        rows.append(dict(tissue=k, g_abs=abs(g), g_arg=np.angle(g), g_genes=int(ex.sum())))
    t = pd.DataFrame(rows).set_index("tissue")
    t["timed_phase"] = -t.g_arg
    t["g_phi_h"] = (t.g_arg / W) % 24
    return t


def chen_load():
    """local: not in clocksync. Chen et al. 2016 prefrontal cortex (GSE71620; BA11 and BA47, time of death known, log2
    microarray): samples (donor, region, `tod_raw` = ZT of death), log2 expression of the waveform genes, the waveform.
    The BA11 / BA47 donor pairing is rebuilt from (age, sex, race, pmi, ph, tod) as in organs/compute_chen_age.py."""
    from clocksync.genes import canonical
    from common.paths import DATA
    raw = DATA / "brain_tod" / "GSE71620"
    key = ["age", "sex", "race", "pmi", "ph", "tod_raw"]
    s = pd.read_csv(raw / "samples.csv")
    s = s[s.tod_raw.notna()].copy()
    s["donor"] = s[key].astype(str).agg("|".join, axis=1)
    s["donor"] = s.donor.map({d: f"d{i:03d}" for i, d in enumerate(sorted(s.donor.unique()))})
    e = pd.read_csv(raw / "expression_clock.csv").set_index("gsm").loc[s.gsm]
    e.columns = [canonical(c.split("|")[0]) for c in e.columns]
    wf = cs.gtex()
    genes = [g for g in wf.genes if g in e.columns]
    y = e[genes].T.groupby(level=0).mean().T[genes]
    return s.reset_index(drop=True), y.reset_index(drop=True), wf.subset(genes)
