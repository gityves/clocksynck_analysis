"""Shared pieces of the A336 part-1 wrappers of the anchor-free cell-type statistics (A272, A274, A280, A282, A283, A288).

All of these read the same object: a matrix S of clock-plane states (gene sets x donor x compartment rows) per atlas and
class medians of |r| over compartment pairs. One stacked matrix holds, in this row order,

    [0, 401)            14 genes: clock, 200 matched, 200 responsive companions (drawn as originally)
    [401, 802)          core loop: core clock, the same companions restricted to the core positions (C.restrict)
    [802, 802 + nd)     clock states of the three-gene drops (C.drops(n=100): 100 random plus the trio drop)
    802 + nd            the trio alone

so that one placebo / pairing pass serves every readout. Fork pools with module globals (no pickling of file-loaded
modules). Aggregates only: atlases are public (the GTEx rule does not apply).
"""
from __future__ import annotations

import sys
import time
import warnings
from multiprocessing import Pool
from pathlib import Path

HERE = Path(__file__).resolve()
EXP = HERE.parents[2]
PAPER = HERE.parents[3]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(HERE.parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import a336_common as C  # noqa: E402
import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402

N_COMP = 200
I14 = slice(0, 1 + 2 * N_COMP)
IC = slice(1 + 2 * N_COMP, 2 * (1 + 2 * N_COMP))
MATCHED = slice(1, 1 + N_COMP)
RESPONSIVE = slice(1 + N_COMP, 1 + 2 * N_COMP)
FAMILIES = {"matched": MATCHED, "responsive": RESPONSIVE}
T0 = time.time()
G = {}  # module globals shared with forked pool workers


def log(*a):
    print(f"[{time.time() - T0:7.0f} s]", *a, flush=True)


def layout(nd):
    """Row index of each block of the stacked matrix for nd drops."""
    return dict(drops=slice(IC.stop, IC.stop + nd), trio=IC.stop + nd, n=IC.stop + nd + 1)


def extra_waveforms(wfs14, dl):
    """Waveforms appended after the 14-gene block: core block, drops, trio (waveform objects, in stacked order)."""
    core = C.core()
    return [core] + [C.restrict(w, C.WF, core) for w in wfs14[1:]] + [w for _, w in dl] + [C.trio()]


def stack_states(pb, wfs14, dl, statefn):
    """Core + drop + trio state rows (added to the 14-gene block by the caller)."""
    return np.array([statefn(pb, w) for w in extra_waveforms(wfs14, dl)])


def resample(s, rng):
    """Positions resampled within stratum codes `s` (donor bootstrap within stratum, as the A272 offset bootstrap)."""
    return np.concatenate([rng.choice(np.where(s == c)[0], (s == c).sum()) for c in range(s.max() + 1)])


def fast_abs(Z, pair):
    """|r| of a pair for several gene sets Z (sets x rows), states centred within stratum (the A272 statistic)."""
    s = pair["s"]
    a, b = Z[:, pair["ia"]], Z[:, pair["ib"]]
    for x in (a, b):
        for c in range(s.max() + 1):
            m = s == c
            x[:, m] -= x[:, m].mean(1, keepdims=True)
    den = np.sqrt(np.sum(np.abs(a) ** 2, 1) * np.sum(np.abs(b) ** 2, 1))
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.abs(np.sum(a * np.conj(b), 1)) / den


def boot_abs(Z2, pairs, rng):
    """One pair-level donor bootstrap: |r| (2 x pairs) of the two state rows of Z2 on shared resamples."""
    out = np.full((len(Z2), len(pairs)), np.nan)
    for k, pr in enumerate(pairs):
        idx = resample(pr["s"], rng)
        out[:, k] = fast_abs(Z2.copy(), dict(ia=pr["ia"][idx], ib=pr["ib"][idx], s=pr["s"][idx]))
    return out


def gate_atlas(pb, trio, n_perm, name, strata=None):
    """Model-free gate on the trio for one donor x compartment table: (passes, p_pairing, p_acrophase)."""
    warnings.filterwarnings("ignore")
    x = pb[[c for c in ("donor", "compartment", "stratum", "n_cells", "depth") if c in pb] + list(trio.genes)]
    g = cs.Dataset(x, anchors=None, waveform=trio, name=name).gate(n_perm=n_perm, seed=336)
    return bool(g.p_pairing < 0.05 and g.p_acrophase < 0.05), float(g.p_pairing), float(g.p_acrophase)


def run_pool(fn, tasks, jobs=4):
    """Map fn over tasks with a fork pool (module globals G are inherited); progress every ~10 minutes."""
    out, last = [], time.time()
    with Pool(jobs, initializer=limit_threads, initargs=(1,)) as p:
        for i, r in enumerate(p.imap(fn, tasks, chunksize=1)):
            out.append(r)
            if time.time() - last > 600:
                log("tasks", i + 1, "of", len(tasks))
                last = time.time()
    return out


def status(met):
    met = list(met)
    return "holds" if all(met) else "at the line" if any(met) else "falls"


def ci90(x):
    x = np.asarray(x, float)
    return float(np.nanquantile(x, 0.05)), float(np.nanquantile(x, 0.95))


def read_stored(path, **kw):
    return pd.read_csv(path, **kw)
