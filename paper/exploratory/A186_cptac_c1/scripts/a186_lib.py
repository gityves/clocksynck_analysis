"""A186 library: CPTAC pairs in the form A175 / A185 expect (reuses run_a175, q3_lib, h138_8.Cancer by import)."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE.parent / "A175_c1_pairing_strata" / "scripts"))
import run_a175 as r75  # noqa: E402  (puts the A138 scripts on the path)
import q3_lib as q  # noqa: E402
from h138_8 import Cancer  # noqa: E402
from clocksync.genes import gene_columns  # noqa: E402

TYPES = ["LUAD", "CCRCC", "LSCC", "SKCM", "HNSCC", "PDAC", "UCEC", "STAD"]
SEED = r75.SEED


def load_type(tp: str):
    """(cases, tumour counts, normal counts) as DataFrames (cases x gene symbols), from the prep cache."""
    z = np.load(HERE / "output" / "cache" / f"{tp}.npz", allow_pickle=False)
    g = [str(x) for x in z["genes"]]
    cases = [str(x) for x in z["cases"]]
    return (cases, pd.DataFrame(z["ct"].astype(np.float64), index=cases, columns=g),
            pd.DataFrame(z["cn"].astype(np.float64), index=cases, columns=g))


def clock_frame(ct: pd.DataFrame) -> pd.DataFrame:
    """Depth + the 14 clock genes under the panel symbols (as the TCGA table)."""
    cols = q.sh.GENES
    idx = gene_columns(list(ct.columns), cols)
    assert (idx >= 0).all(), [c for c, i in zip(cols, idx) if i < 0]
    d = pd.DataFrame(ct.to_numpy()[:, idx], columns=cols, index=ct.index)
    d["depth"] = ct.sum(axis=1).to_numpy()
    return d


def states(tp: str, i: int, n_gate: int = 200, seed: int = 0) -> dict:
    """As `q3_lib.tcga_states` for one type: gate of the normals (A52 rule, p < 0.05), states of both on one projector."""
    cases, ct, cn = load_type(tp)
    case, ref = clock_frame(ct), clock_frame(cn)
    yr = np.log(np.maximum(ref[q.GENES].astype(float), 1e-300)).sub(np.log(ref.depth.to_numpy(float)), axis=0).to_numpy()
    g = q.reference_gate(yr, q.WF, n_gate, np.random.default_rng([52, seed, i]))
    bt, vt, bn, vn = q.cohort_states(case, ref, True, q.WF)
    return dict(za=q.centred_complex(bt), vA=vt, zb=q.centred_complex(bn), vB=vn, n=len(cases), gate_p=g["p"],
                gated=bool(g["p"] < 0.05), raw=(bt, vt, bn, vn), cases=cases,
                zero_clock_counts_normal=int((ref[q.GENES].to_numpy() == 0).sum()),
                zero_clock_counts_tumour=int((case[q.GENES].to_numpy() == 0).sum()))


def labels(tp: str):
    """G2M score (A97, `run_a175.g2m_scores`) of the tumours of the paired cases and rank tertiles."""
    _, ct, _ = load_type(tp)
    g = r75.g2m_scores(ct)
    return g, r75.tertiles(g)


def companion_type(tp: str, i: int, n_sets: int, lab: np.ndarray):
    """As `run_a175.companion_cancer` (same RNG streams by index i) for a CPTAC type; stats per set on labels `lab` (0 low / 2 high)."""
    from clocksync.config import limit_threads
    from clocksync.paired import cohort_states, pseudo_waveforms
    limit_threads(1)
    cases, ct, cn = load_type(tp)
    lm = np.log(cn.div(cn.sum(axis=1), axis=0).mean() + 1e-12)[cn.mean() >= 1]
    dt, dn = ct.sum(axis=1), cn.sum(axis=1)
    out = {}
    for kind, tag, kw in (("matched", 44, {}), ("responsive", 45, {"kind": "responsive"})):
        ws = pseudo_waveforms(lm, n_sets, np.random.default_rng([SEED, tag, 0, i]), q.WF, **kw)
        res = []
        for w in ws:
            g = list(w.genes)
            bt, vt, bn, vn = cohort_states(ct[g].assign(depth=dt), cn[g].assign(depth=dn), True, w)
            C = Cancer(tp, dict(n=len(cases), raw=(bt, vt, bn, vn)))
            res.append(C)
        out[kind] = res
    return out
