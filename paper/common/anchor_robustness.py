"""Anchor robustness of a host-referenced contrast (WP0 item 8).

When anchors are chosen by rule (or the dataset holds only disease donors), the host depends on that choice. Two
readouts, both built on the public `clocksync` API (no change to the package is needed):

* `block_agreement`: r* between the host fitted on the anchor block and the host fitted on the remaining compartments
  (`Dataset.agreement`, as the brain-periphery test of A103). r* near 1: the anchors and the rest carry one host, so
  the choice of anchors does not matter; r* clearly < 1: the host is anchor-dependent.
* `host_contrast`: the contrast of the same labels on (a) the anchored host and (b) the host of all compartments
  except the target (`Dataset(anchors=all compartments)`: every compartment is read on its leave-one-out host), on
  one set of permutations per host. `sign_change` flags a target whose case-minus-control difference changes sign.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

import clocksync as cs

KEEP = ["donors_control", "donors_case", "c_control", "c_case", "d", "se_d", "p", "d_star", "p_star",
        "log_z2_ratio", "se_log_z2_ratio", "offset_diff_h", "status"]


def block_agreement(ds: cs.Dataset, anchors=None, n_boot: int = 100, n_split: int = 30, seed: int = 0,
                    min_common: int = 10) -> pd.Series:
    """r* of the anchor block against all other compartments of `ds` (default anchors: `ds.anchors`).

    Needs >= 2 compartments per block (otherwise NaN with a `note`). The split-half reliability needs >= 4 per block:
    with fewer, r* is undefined (`note`) and only the raw correlation `r` (not disattenuated, so a lower bound) is
    reported."""
    anc = list(ds.anchors if anchors is None else anchors)
    rest = [k for k in ds.compartments if k not in anc]
    base = dict(n_anchors=len(anc), n_rest=len(rest))
    if len(anc) < 2 or len(rest) < 2:
        return pd.Series({**base, "r_star": np.nan, "note": "block with < 2 compartments"})
    try:
        r = ds.agreement({"anchors": anc, "rest": rest}, n_boot=n_boot, n_split=n_split, seed=seed,
                         min_common=min_common)
    except (ValueError, np.linalg.LinAlgError, KeyError) as e:
        return pd.Series({**base, "r_star": np.nan, "note": f"agreement failed: {e}"})
    return pd.concat([pd.Series(base), r, pd.Series({"note": "" if np.isfinite(r["r_star"]) else "r* undefined (block < 4 compartments or no valid split); raw r only"})])


def loo_dataset(ds: cs.Dataset) -> cs.Dataset:
    """Same data, every compartment an anchor: each is read on the host of all the others (leave-one-out)."""
    return ds.subset(anchors=ds.compartments)


def sign_change(a, b) -> pd.Series:
    """True where the two estimates have opposite non-zero signs (NaN where either is missing)."""
    a, b = pd.Series(a, dtype=float), pd.Series(b, dtype=float)
    out = pd.Series(np.sign(a) * np.sign(b) < 0, index=a.index, dtype=object)
    out[a.isna() | b.isna()] = np.nan
    return out


def host_contrast(ds: cs.Dataset, groups, targets, levels=("control", "case"), blocks=None, n_perm: int = 200,
                  n_boot: int = 100, seed: int = 0) -> pd.DataFrame:
    """Anchored and all-compartment leave-one-out contrast of `targets`, one row per target.

    Columns suffixed `_anchored` / `_loo` (d, se_d, p, d_star, p_star, log_z2_ratio, ...), `sign_change` (on d),
    `sign_change_z2`; `anchors` lists the anchored host."""
    kw = dict(levels=levels, blocks=blocks, n_perm=n_perm, n_boot=n_boot, seed=seed)
    tabs = {}
    for host, d in (("anchored", ds), ("loo", loo_dataset(ds))):
        t = d.contrast(groups, targets=targets, **kw).table.reindex(targets)
        tabs[host] = t.rename(columns={f"donors_{levels[0]}": "donors_control", f"donors_{levels[1]}": "donors_case",
                                       f"c_{levels[0]}": "c_control", f"c_{levels[1]}": "c_case"})
    out = pd.concat([tabs[h].reindex(columns=[c for c in KEEP if c in tabs[h]]).add_suffix(f"_{h}")
                     for h in ("anchored", "loo")], axis=1)
    out.index.name = "target"
    out["sign_change"] = sign_change(out.d_anchored, out.d_loo)
    out["sign_change_z2"] = sign_change(out.log_z2_ratio_anchored, out.log_z2_ratio_loo)
    out["anchors"] = ";".join(ds.anchors)
    return out.reset_index()
