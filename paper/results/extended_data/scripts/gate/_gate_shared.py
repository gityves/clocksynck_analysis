"""Loaders and small helpers shared by the scripts of gate.

Method calls are clocksync's. Everything marked "local: not in clocksync" is for something the package has no public
call for (listed in the folder README). The A25 internal-time estimator and the baboon loaders are those of
windows (imported from there: `window_tools`).
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import clocksync as cs  # noqa: E402
from common import anchors as AN  # noqa: E402
from common import data, gtex  # noqa: E402


sys.path.insert(0, str(Path(__file__).resolve().parents[4]))  # paper/
from common.paths import DATASETS, group_dir  # noqa: E402
TARGETS = DATASETS / "_atlas" / "a17_targets.csv"
W = 2 * np.pi / 24


def _load(path: Path, name: str):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


S02 = group_dir("windows")
sys.path.insert(0, str(S02))
import _internal_time as T  # noqa: E402,F401  (A25 estimator, local to windows)
s02 = _load(S02 / "_shared.py", "s02_shared")  # baboon loaders, timed tissue gains (same-named module of another folder)


def targets() -> pd.DataFrame:
    return pd.read_csv(TARGETS).set_index("name", drop=False)


def anchor_set(name: str, rule: str):
    t = pd.read_csv(AN.TABLE)
    t = t[t.target == name]
    if t.empty:
        return None
    anc = AN.anchors(name, rule)
    if len(anc) < 3 or (rule == "learned" and str(t.status.iloc[0]).startswith("no anchor")):
        return None
    return anc


def family_map() -> dict:
    t = pd.read_csv(AN.TABLE)
    return t.drop_duplicates("compartment").set_index("compartment").family.to_dict()
