"""Shared pieces of the single-cell skin analyses (A58, A60, A61, A63-A66, A69): status-split tables, the free-host fit
with offsets and amplitudes of every compartment, unit bootstrap and status swaps, wrapped percentile intervals.

Imported by `compute_sc_status.py` here and by the single-cell scripts of `psoriasis_details`.
Estimator: `clocksync.Dataset(table, anchors=None).fit(n_boot=0, gate=False)`: the host is fitted on all compartments,
offsets and log amplitudes are read on that host and are reference-free (centred over compartments); all differences used
here (lesional minus uninvolved, bin minus bin) are independent of the centring. The legacy `Study.fit` took the loadings of
the joint factor model; `Dataset.fit` regresses each compartment on the posterior host state, so values differ in the
second digit (see the folder README).
"""

from __future__ import annotations

import importlib.util
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import clocksync as cs  # noqa: E402
from clocksync import stats as S  # noqa: E402


warnings.filterwarnings("ignore")
HOURS = S.HOURS
STATUS = {"lesional": "L", "nonlesional": "N"}
MIN_CELLS, MIN_UNITS = 20, 5
DERMAL = ["fibroblast", "endothelial", "mural"]
BINS = ["kc_d1", "kc_d2", "kc_d3", "kc_d4"]


def wrap(h):
    return S.wrap(h)


def cmean(h) -> float:
    """Circular mean of offsets in hours."""
    return float(S.circmean(np.asarray(h, float)))


sys.path.insert(0, str(Path(__file__).resolve().parents[4]))  # paper/
from common.paths import DATASETS  # noqa: E402
def load(name: str, subsample: int | None = None, seed: int = 0):
    """Standard table of a skin single-cell builder; `subsample`: keep that many random units (`donor`)."""
    f = DATASETS / name / "build.py"
    spec = importlib.util.spec_from_file_location(f"build_{name}", f)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    pb, donors = mod.build()
    if subsample is not None:
        keep = np.random.default_rng(seed).permutation(sorted(pb.donor.unique()))[:subsample]
        pb, donors = pb[pb.donor.isin(keep)], donors.loc[donors.index.isin(keep)]
    return pb.reset_index(drop=True), donors


def status_split(pb: pd.DataFrame, comps=None, min_cells: int = MIN_CELLS) -> pd.DataFrame:
    """Compartments `L:<class>` / `N:<class>`; units with both statuses only (A58 rule)."""
    d = pb[pb.n_cells >= min_cells].copy()
    if comps is not None:
        d = d[d.compartment.isin(comps)]
    d["compartment"] = d.status.map(STATUS) + ":" + d.compartment.astype(str)
    both = d.groupby("donor").status.nunique()
    return d[d.donor.isin(both.index[both == 2])].reset_index(drop=True)


def classes(t: pd.DataFrame, order=None, min_units: int = MIN_UNITS) -> list[str]:
    """Classes present in both statuses in >= `min_units` units."""
    n = t.groupby("compartment").donor.nunique()
    ks = order or sorted({c[2:] for c in n.index})
    return [c for c in ks if n.get(f"L:{c}", 0) >= min_units and n.get(f"N:{c}", 0) >= min_units]


def dataset(t: pd.DataFrame, min_donors: int = MIN_UNITS) -> cs.Dataset:
    return cs.Dataset(t, anchors=None, min_donors=min_donors)


def fit_states(t: pd.DataFrame, min_donors: int = MIN_UNITS) -> pd.DataFrame:
    """Offset (h), log amplitude, rel. amplitude and coupling of every compartment (free host, reference-free)."""
    tab = dataset(t, min_donors).fit(n_boot=0, gate=False).table
    return tab[["offset_h", "log_amp_rel", "rel_amplitude", "coupling", "n"]]


def gate(t: pd.DataFrame, seed: int, n_perm: int = 200, min_donors: int = MIN_UNITS) -> dict:
    g = dataset(t, min_donors).gate(n_perm=n_perm, seed=seed)
    return dict(clock_share=g.clock_share, p_pairing=g.p_pairing, p_acrophase=g.p_acrophase, n_pairs=g.n_pairs,
                n_donors=g.n_donors, gate_pass=bool(g.p_pairing < 0.05 and g.p_acrophase < 0.05))


def swapped(t: pd.DataFrame, rng) -> pd.DataFrame:
    """Statuses swapped per unit with probability 1/2 (A58 null)."""
    u = t.donor.unique()
    flip = set(u[rng.random(len(u)) < 0.5])
    s = t.copy()
    f = s.donor.isin(flip)
    c = s.loc[f, "compartment"]
    s.loc[f, "compartment"] = np.where(c.str.startswith("L:"), "N:" + c.str[2:], "L:" + c.str[2:])
    return s


def boot(t: pd.DataFrame, rng) -> pd.DataFrame:
    """Unit bootstrap; resampled units relabelled `<unit>#j`."""
    u = t.donor.unique()
    parts = []
    for j, r in enumerate(rng.choice(u, len(u))):
        x = t[t.donor == r].copy()
        x["donor"] = f"{r}#{j}"
        parts.append(x)
    return pd.concat(parts, ignore_index=True)


def safe(f, *a):
    try:
        return f(*a)
    except Exception:
        return None


def interval(bv, v: float, circular: bool = True) -> tuple[float, float]:
    """95 % percentile interval of bootstrap values around the estimate (offsets wrapped around it)."""
    bv = np.asarray(bv, float)
    bv = bv[np.isfinite(bv)]
    if not len(bv):
        return (np.nan, np.nan)
    d = np.asarray(wrap(bv - v)) if circular else bv - v
    return float(v + np.quantile(d, 0.025)), float(v + np.quantile(d, 0.975))


def p_lower(null, obs: float) -> float:
    return float(S.perm_p(np.asarray(null, float)[np.isfinite(null)], obs, "less", tie=0.0))


def sizes(args, run: dict, n_perm: int, n_boot: int) -> dict:
    q = run["quick"]
    return dict(n_perm=args.n_perm or (30 if q else n_perm), n_boot=args.n_boot or (20 if q else n_boot))
