"""Anchor rule of amendment A17 (step 2): class families and the leave-one-study-out selection of the compartments that
define the host of a dataset. local: not in clocksync (the package ships `suggest_anchors`, a name rule, and takes the
learned anchors as input; see the README, API gaps).

- `family`: cell family of a compartment from its class name (tissue and region prefixes ignored).
- `select`: pool each family's coherence excess over the gate-passing healthy datasets of the other studies (same tissue
  group if any has the family, else all tissues) by DerSimonian-Laird random effects (`clocksync.stats.meta`); anchors =
  the target's qualifying families (pooled lower 95 % limit > 0) in the top half of its families, plus lower qualifying
  families until >= 3 compartments; the target's own pair coherences can then only drop an anchor (an anchor whose pair
  excess is <= 0 with every other anchor); fewer than 3 left: "no anchor".
- `fixed`: the sensitivity rule (stromal-vascular families; neurons in brain and retina).
"""

from __future__ import annotations

import math
import re

import numpy as np
import pandas as pd

from clocksync import stats as S

_ALT = lambda *w: r"(?:^|_)(?:" + "|".join(w) + r")$"
FAMILIES = (  # first match wins
    ("other", _ALT("glia_neuron", "progenitor", "erythroid", "melanocyte")),
    ("mural", _ALT("mural", "smooth_muscle", "mural_smooth_muscle", "pericyte")),
    ("endothelial", _ALT("endothelial", r"EC(?:_\w+)?", "lymphatic_EC", "vascular")),
    ("fibroblast_stromal", _ALT("fibroblast", "stromal", "fibroblast_stromal", "interstitial", "stellate")),
    ("neuron", _ALT("excitatory", r"inhibitory(?:_\w+)?", "glutamatergic", "GABAergic", "neuron", "neuron_other",
                    "amacrine", "bipolar", "horizontal", "ganglion", "rod", "cone")),
    ("glia", _ALT("astrocyte", "oligodendrocyte", "OPC", "Mueller", "ependymal")),
    ("myeloid", _ALT("macrophage", "macrophage_DC", r"\w*monocyte", "cMono", "ncMono", "DC", "cDC", "microglia", "mast",
                     "granulocyte", "myeloid")),
    ("lymphocyte", _ALT("T", r"T_\w+", "CD4_T", "CD8_T", r"CD[48]_\w+", "Treg", "gdT", "NK", "NK_ILC", "ILC", "B",
                        r"B_\w+", "plasma", "lymphocyte", "T_NK", "B_plasma")),
    ("parenchyma", _ALT(r"\w*cardiomyocyte", "proximal_tubule", "loop_of_Henle", "distal", "collecting_duct",
                        "podocyte", "adipocyte", "muscle", "alpha", "beta", "delta")),
    ("epithelial", _ALT("epithelial", "AT1", "AT2", "basal", "secretory", "multiciliated", "ciliated",
                        "keratinocyte", "enterocyte", "colonocyte", "stem", "ductal", "acinar")),
)
STROMAL_VASCULAR = ("fibroblast_stromal", "endothelial", "mural")
CNS = ("brain", "retina")


def family(compartment: str) -> str:
    cls = str(compartment).split(":", 1)[-1]
    return next((f for f, pat in FAMILIES if re.search(pat, cls)), "other")


def pool(y, se) -> tuple[float, float, float]:
    """Random-effects pooled estimate, SE and tau^2 (`clocksync.stats.meta`; one study: the study itself)."""
    m = S.meta(y, se)
    return m["random"], m["se_random"], m["tau2"]


def family_table(classes: pd.DataFrame) -> pd.DataFrame:
    """Per dataset x family: mean excess of its classes, SE = mean of their SEs (conservative)."""
    c = classes.assign(family=classes.compartment.map(family))
    return c.groupby(["dataset", "family"], as_index=False).agg(excess=("excess", "mean"), se=("se", "mean"),
                                                                 n_classes=("compartment", "size"))


def select(target: dict, comps, pairs_t: pd.DataFrame | None, learn: pd.DataFrame, datasets: pd.DataFrame):
    """Anchors of one target. target: dict(study, tissue). comps: its compartments. pairs_t: its own pair coherences
    (k, l, excess) for the confirmation step (None: none). learn: `family_table` per study. datasets: study, tissue,
    gate_pass. Returns (anchors list or [], status, table of the target's families with pooled scores)."""
    meta = datasets.set_index("dataset")
    ok = meta.index[meta.gate_pass & (meta.study != target["study"])]
    L = learn[learn.dataset.isin(ok)].merge(meta[["tissue"]], left_on="dataset", right_index=True)
    fam = pd.Series({k: family(k) for k in comps})
    rows = []
    for f in sorted(set(fam) - {"other"}):
        same = L[(L.family == f) & (L.tissue == target["tissue"])]
        use, scope = (same, "tissue") if len(same) else (L[L.family == f], "all tissues")
        est, se, tau2 = pool(use.excess, use.se)
        rows.append(dict(family=f, scope=scope, n_datasets=len(use), pooled=est, se=se, tau2=tau2,
                         lower=est - 1.96 * se if np.isfinite(est) else np.nan,
                         n_compartments=int((fam == f).sum())))
    tab = pd.DataFrame(rows, columns=["family", "scope", "n_datasets", "pooled", "se", "tau2", "lower",
                                      "n_compartments"])
    if tab.empty:
        return [], "no anchor (no family)", tab
    tab = tab.sort_values("pooled", ascending=False, na_position="last").reset_index(drop=True)
    tab["qualifies"] = tab.lower > 0
    top = math.ceil(len(tab) / 2)
    chosen = [f for f in tab.family[:top] if tab.set_index("family").qualifies[f]]
    for f in tab.family[top:]:
        if sum(int((fam == g).sum()) for g in chosen) >= 3:
            break
        if tab.set_index("family").qualifies[f]:
            chosen.append(f)
    tab["anchor_family"] = tab.family.isin(chosen)
    anchors = [k for k in comps if fam[k] in chosen]
    if pairs_t is not None and len(anchors) >= 2:
        p = pairs_t[pairs_t.k.isin(anchors) & pairs_t.l.isin(anchors)]
        keep = [a for a in anchors if (p[(p.k == a) | (p.l == a)].excess > 0).any()]
        dropped = sorted(set(anchors) - set(keep))
        anchors = keep
    else:
        dropped = []
    if len(anchors) < 3:
        return [], f"no anchor ({len(anchors)} compartments)" + (f"; dropped {dropped}" if dropped else ""), tab
    return anchors, "ok" + (f"; dropped {dropped}" if dropped else ""), tab


def fixed(target: dict, comps) -> list:
    """Sensitivity rule: stromal-vascular families; neurons in brain and retina."""
    fams = ("neuron",) if target["tissue"] in CNS else STROMAL_VASCULAR
    a = [k for k in comps if family(k) in fams]
    return a if len(a) >= 3 else []
