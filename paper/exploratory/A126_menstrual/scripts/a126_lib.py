"""A126 shared helpers: coherence of the 14 clock genes of a group of pseudobulks (A87 `decompose.components`), expression-matched
pseudo-clock sets (variance-matched draw of A127), Holm."""
import sys
from pathlib import Path
import numpy as np, pandas as pd

ROOT = Path(__file__).resolve().parents[4]
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "paper"))
import clocksync as cs  # noqa: E402
from clocksync.decompose import components, logfrac  # noqa: E402
from clocksync.waveform import Waveform  # noqa: E402
from clocksync._core.controls import CLOCK_LIKE  # noqa: E402
from clocksync.genes import canonical  # noqa: E402
from clocksync._core.controls import CLOCK_LIKE  # noqa: E402,F811

OUT = Path(__file__).resolve().parents[1] / "output"
WF = cs.gtex()
GENES = list(WF.genes)
LUTEAL = ["PGR", "TNFSF11", "WNT4", "AREG"]
DATA = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian")


def wf_for(genes) -> Waveform:
    return Waveform(tuple(genes), WF.amp, WF.psi, None, None, "set")


def coh_counts(C, depth, wf=WF) -> float:
    """Coherence of a group of pseudobulks: C (n x 14 counts), depth (n). Log fractions centred over the group, Poisson term removed."""
    n = np.asarray(C, float)
    y = np.log((n + 0.5) / np.asarray(depth, float)[:, None])
    y = y - y.mean(0)
    return components(y, (1 / (n + 0.5)).mean(0), 1, wf)["coherence"]


def holm(p):
    p = np.asarray(p, float); o = np.argsort(p); out = np.empty(len(p)); run = 0
    for r, k in enumerate(o):
        run = max(run, (len(p) - r) * p[k]); out[k] = min(1, run)
    return out


def pseudo_sets(mean, sd, pool, n_sets, rng):
    """A127 variance-matched draw: each clock gene gets a pseudo gene from the 40 nearest (mean, sd) neighbours, no reuse in a set.
    mean, sd: Series over `pool` plus the clock genes."""
    pm, ps = mean[pool], sd[pool]
    zs = lambda a, ref: (a - ref.mean()) / ref.std()
    P = np.column_stack([zs(pm.values, pm), zs(ps.values, ps)])
    T = np.column_stack([zs(mean[GENES].values, pm), zs(sd[GENES].values, ps)])
    sets = []
    for _ in range(n_sets):
        used, s = set(), []
        for t in rng.permutation(len(GENES)):
            near = [j for j in np.argsort(((P - T[t]) ** 2).sum(1))[:40] if j not in used][:20]
            j = int(rng.choice(near)); used.add(j); s.append((t, pool[j]))
        sets.append([g for _, g in sorted(s)])
    return sets
