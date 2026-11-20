"""Shared pieces of the A336 part-1 wrappers of this agent batch (GTEx 12-tissue machinery, gates, bootstraps, status)."""
import copy
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
EXP = HERE.parents[2]
PAPER = HERE.parents[3]
sys.path.insert(0, str(PAPER))
for d in ("A109_coupling_nulls", "A132_brain_death", "A315_braindeath_gain", "A317_braindeath_confirm",
          "A231_age_coupling", "A335_review_batch"):
    sys.path.insert(0, str(EXP / d / "scripts"))
sys.path.insert(0, str(HERE.parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import a335_gtex as G  # noqa: E402
import _helpers as H  # noqa: E402
import clocksync as cs  # noqa: E402
from a336_common import TRIO, WF, core, trio  # noqa: E402

SEED_BOOT = [336, 1]
N_BOOT_UNIT = 200


def build(n_sets=200, subsample=None):
    """A109 tissues with 14-gene states (clock, matched, responsive) and the core ones, as A335: (T, Tc)."""
    T = H.Tissues(n_sets, 0, subsample)
    T.build_states(np.random.default_rng(1), 0, real=True)
    Zc, Vc = G.build_core(T, None)
    G.add_resp(T, Zc, Vc, n_sets=n_sets)
    Tc = copy.copy(T)
    Tc.Z, Tc.V = Zc, Vc
    return T, Tc


def add_states(T, kind, waveforms):
    """Clock-plane states of `waveforms` on the real clock-gene rows: T.Z[kind][t] (K, rows)."""
    T.Z[kind], T.V[kind] = {}, {}
    for t in H.BODY:
        T.Z[kind][t], T.V[kind][t] = H._state_set(T.rows[t], waveforms, T.rows[t].depth.to_numpy(float))


def gate(T, donors, wf, n_perm=1000):
    """cs.Dataset(...).gate on the 12 body tissues' rows of `donors` for waveform `wf`: (passes, p_pairing, p_acro)."""
    rows = []
    for t in H.BODY:
        r = np.unique(T.R[t][np.asarray(donors)])
        r = r[r >= 0]
        x = T.rows[t].iloc[r].copy()
        x.insert(0, "compartment", t)
        x.insert(0, "donor", np.asarray(T.row_donor[t])[r])
        rows.append(x)
    ds = cs.Dataset(pd.concat(rows, ignore_index=True), anchors=None, waveform=wf, name="trio")
    g = ds.gate(n_perm=n_perm, seed=336)
    return bool(g.p_pairing < 0.05 and g.p_acrophase < 0.05), float(g.p_pairing), float(g.p_acrophase)


def ci90(x):
    x = np.asarray(x, float)
    return float(np.nanquantile(x, 0.05)), float(np.nanquantile(x, 0.95))


def status(met, total=None):
    """holds / at the line / falls from a list of booleans (conditions of the registered criterion)."""
    met = list(met)
    return "holds" if all(met) else "at the line" if any(met) else "falls"
