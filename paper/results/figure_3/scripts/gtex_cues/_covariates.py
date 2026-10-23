"""Helpers of the covariate and disturbance scripts of figure_S07 (A46, A47, A92, A93, A27 part 2): cohorts, protected
phenotypes (in memory only), module statistics and permutations. Method calls are clocksync's; everything marked
"local: not in clocksync" is a small helper for a statistic the package has no public call for."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _shared as S  # noqa: E402
from common import gtex  # noqa: E402

FOLDER = "gtex_cues"
EXTRA = ("Liver", "Pancreas")
ALL_T = tuple(S.TISSUES) + EXTRA


def cohort(hardy, tissues=ALL_T, subsample: int | None = None, seed: int = 0, require_tod: bool = False):
    """(pseudobulk, donors) of GTEx donors of the given Hardy classes (depth = library size, as the legacy A46-A47
    runs); `subsample` keeps that many random donors (checks)."""
    pb, d = gtex.load(tissues=list(tissues), hardy=hardy, genes="clock", depth="library_size", require_tod=require_tod)
    if subsample:
        keep = np.random.default_rng(seed).permutation(sorted(pb.donor.unique()))[:subsample]
        pb, d = pb[pb.donor.isin(keep)].reset_index(drop=True), d.loc[d.index.isin(keep)]
    return pb, d


def cause_group(t):
    """Coarse cause-of-death group from the free text DTHCOD (keyword rules of the legacy `feasibility.cause_group`)."""
    if not isinstance(t, str):
        return np.nan
    t = t.lower()
    for g, keys in (("cerebrovascular", ("cva", "stroke", "intracranial", "ich", "hemorrhage", "aneurysm", "sah")),
                    ("cardiac", ("cardiac", "heart", "mi", "myocard", "arrhythm", "cardio")),
                    ("trauma", ("trauma", "gunshot", "gsw", "mvc", "accident", "fall", "crush", "stab")),
                    ("anoxia_resp", ("anoxia", "overdose", "asphyx", "drown", "hang", "respir", "pulmonary", "copd",
                                     "pneumon"))):
        if any(k in t.replace("/", " ").split() or (len(k) > 3 and k in t) for k in keys):
            return g
    return "other"


def protected() -> pd.DataFrame:
    """Protected phenotypes by donor (controlled access, in memory only): exact age, BMI, ventilation (flag, hours),
    T2D, long-term steroids, sepsis, liver disease, season, cause-of-death group; 99 (unknown) -> NaN."""
    p = gtex.phenotypes(["AGE", "BMI", "DTHVNT", "DTHVNTD", "MHT2D", "MHSTRDLT", "MHSEPSIS", "MHLVRDIS", "DTHSEASON",
                         "DTHCOD"])
    for c in ("DTHVNT", "MHT2D", "MHSTRDLT", "MHSEPSIS", "MHLVRDIS"):
        p[c] = p[c].where(p[c].isin([0, 1]))
    out = pd.DataFrame({"age": pd.to_numeric(p.AGE, errors="coerce"), "bmi": pd.to_numeric(p.BMI, errors="coerce"),
                        "vent": p.DTHVNT, "hours": pd.to_numeric(p.DTHVNTD, errors="coerce"), "t2d": p.MHT2D,
                        "steroid": p.MHSTRDLT, "sepsis": p.MHSEPSIS, "liver_dis": p.MHLVRDIS, "season": p.DTHSEASON})
    out["cause"] = p.DTHCOD.map(cause_group).fillna("unknown")
    return out


def module(vals: dict, plus, minus) -> float:
    """Mean over `plus` minus mean over `minus` of finite values (tissue -> value); NaN if either side is empty."""
    p = [vals[t] for t in plus if np.isfinite(vals.get(t, np.nan))]
    m = [vals[t] for t in minus if np.isfinite(vals.get(t, np.nan))]
    return float(np.mean(p) - np.mean(m)) if p and m else np.nan


def perm_within(groups, rng) -> np.ndarray:
    """Index permutation that shuffles within each group label (local: not in clocksync)."""
    groups = np.asarray(groups)
    out = np.arange(len(groups))
    for g in np.unique(groups):
        i = np.flatnonzero(groups == g)
        out[i] = rng.permutation(i)
    return out


def p_lower(obs: float, null) -> float:
    """One-sided add-one permutation p for statistic <= obs."""
    null = np.asarray(null, float)
    null = null[np.isfinite(null)]
    return float((1 + np.sum(null <= obs + 1e-12)) / (1 + len(null)))


def p_upper(obs: float, null) -> float:
    null = np.asarray(null, float)
    null = null[np.isfinite(null)]
    return float((1 + np.sum(null >= obs - 1e-12)) / (1 + len(null)))


def tertiles(x: pd.Series) -> pd.Series:
    """Tertile labels 0-2 by rank (ties broken by order), as `pd.qcut` on ranks."""
    return pd.qcut(x.rank(method="first"), 3, labels=False)
