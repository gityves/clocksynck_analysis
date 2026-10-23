"""Frozen anchors of amendment A17 (learned leave one study out from the model-free coherence, or the fixed
stromal-vascular rule), one set per target dataset: `anchors_a17.csv`, copied from the legacy run."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

TABLE = Path(__file__).with_name("anchors_a17.csv")


def anchors(target: str, rule: str = "learned") -> list:
    """Anchor compartments of a target dataset (rule: "learned" or "fixed")."""
    t = pd.read_csv(TABLE)
    t = t[t.target == target]
    if t.empty:
        raise KeyError(f"no A17 anchors for {target!r}; targets: {sorted(pd.read_csv(TABLE).target.unique())[:20]}...")
    col = "anchor" if rule == "learned" else "anchor_fixed"
    return t.loc[t[col].astype(bool), "compartment"].tolist()
