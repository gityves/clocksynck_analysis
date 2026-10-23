"""Shared by the scripts of spread_multiplicity: the A20 contrast table and a donor bootstrap of the
contrast estimands (local: not in clocksync)."""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import clocksync as cs  # noqa: E402
from results.figure_4.scripts.resilience import _shared as F5  # noqa: E402

GROUPS = F5.GROUPS
HOURS = 24 / (2 * np.pi)
SE_MAX_H = 4.0  # A23: rows with SE > 4 h stay out of the pool
CLASSES = {  # A20 table: contrast -> disease class (contrasts sharing donors already removed)
    "ild_natri": "fibrosis", "hlca_fibrotic": "fibrosis", "dcm": "fibrosis", "ckd_sc": "fibrosis",
    "skin_ssc": "fibrosis",
    "skin_ad": "inflammation", "myocarditis": "inflammation", "ia_crohn": "inflammation",
    "endometriosis": "inflammation",
    "AD_mssm": "neurodegeneration", "AD_radc": "neurodegeneration", "AD_seaad_dlpfc": "neurodegeneration",
    "DLBD_mssm": "neurodegeneration", "PD_mssm": "neurodegeneration", "VaD_mssm": "neurodegeneration",
    "pd": "neurodegeneration",
    "SCZ_hbcc": "psychiatric", "BD_hbcc": "psychiatric", "SCZ_mssm": "psychiatric",
    "islets_t2d": "metabolic", "DM_radc": "metabolic", "ASCVD_hbcc": "metabolic",
}
COHORT = {"AD_mssm": "MSSM", "DLBD_mssm": "MSSM", "PD_mssm": "MSSM", "VaD_mssm": "MSSM", "SCZ_mssm": "MSSM",
          "SCZ_hbcc": "HBCC", "BD_hbcc": "HBCC", "ASCVD_hbcc": "HBCC", "AD_radc": "RADC", "DM_radc": "RADC"}
ESTIMANDS = {"y1": "d_star", "y2": "log_z2_ratio", "y3": "log_psi_ratio", "y4": "offset_diff_h"}


def wrap(h):
    """Wrap hours to (-12, 12]."""
    return (np.asarray(h, float) + 12) % 24 - 12


def estimands(ds: cs.Dataset, lab: pd.Series, blocks=None, targets=None) -> pd.DataFrame:
    """y1 = d*, y2 = log z2 ratio, y3 = log psi ratio, y4 = relative phase (h) per compartment of one contrast
    (`Dataset.contrast` without resampling; y2, y4 undefined when z2 <= 0 in either group, as legacy)."""
    t = ds.contrast(lab, levels=GROUPS, targets=targets, blocks=blocks, n_perm=0, n_boot=0).table
    return t[list(ESTIMANDS.values())].astype(float)


def boot_se(pb: pd.DataFrame, anchors, lab: pd.Series, targets, n_boot: int, seed: int, point: pd.DataFrame,
            blocks=None) -> pd.DataFrame:
    """local: not in clocksync (the package bootstraps d, not d*, the log ratios or the offset). Donors resampled
    within each group (repeats renamed `<id>.<j>`), the Dataset rebuilt in each replicate (legacy: host and V0 fixed, so
    these SEs include the host refit). SD of each estimand over replicates; the offset is wrapped around the point
    estimate. Needs more than min(20, n_boot / 2) finite replicates."""
    rng = np.random.default_rng(seed)
    by = {d: g for d, g in pb.groupby("donor")}
    reps = []
    for _ in range(n_boot):
        pick = [(g, d) for g in GROUPS for d in rng.choice(lab.index[lab == g].to_numpy(), int((lab == g).sum()))]
        pbb = pd.concat([by[d].assign(donor=f"{d}.{j}") for j, (_, d) in enumerate(pick)], ignore_index=True)
        lb = pd.Series({f"{d}.{j}": g for j, (g, d) in enumerate(pick)})
        try:
            r = estimands(cs.Dataset(pbb, anchors=anchors), lb, "none", targets)
            reps.append(r.reindex(targets).to_numpy(float))
        except (ValueError, KeyError, np.linalg.LinAlgError):
            continue
    cols = list(ESTIMANDS.values())
    a = np.stack(reps) if reps else np.empty((0, len(targets), len(cols)))
    j = cols.index("offset_diff_h")
    p4 = point["offset_diff_h"].reindex(targets).to_numpy(float)
    if len(a):
        a[:, :, j] = p4 + wrap(a[:, :, j] - p4)
    n_ok = np.isfinite(a).sum(0)
    sd = np.where(n_ok > min(20, n_boot // 2), np.nanstd(a, axis=0, ddof=1), np.nan) if len(a) else np.full((len(targets), len(cols)), np.nan)
    out = pd.DataFrame(sd, index=targets, columns=[f"se{i + 1}" for i in range(len(cols))])
    out["n_boot"] = len(a)
    return out
