"""A336 part 1 shared harness: waveforms, the companion restriction rule, the gene-drop null and the registered reading.

Every part-1 wrapper (`p1_<entry>.py`) imports this module and writes one row per statistic with `write_rows` to
`../output/part1/<entry>.csv`; `merge_part1.py` concatenates them into `../output/two_readouts.csv`.

Registered rules (A336): 14 genes reproduced first; core = cs.gtex().drop(TRIO) with companions drawn as originally and
restricted to the core positions; trio contribution on shared bootstraps (90 %); gene-drop null = the clock statistic on
the 364 three-gene drops (or 100 random + trio, seed 336), magnitude in the direction of the original effect; a fall
reads "carried by PER1, PER2 and NR1D1" only if the trio drop is below the 10th percentile of the drops.
"""
import itertools
from pathlib import Path

import numpy as np
import pandas as pd

import clocksync as cs
from clocksync.waveform import Waveform

OUT = Path(__file__).resolve().parents[1] / "output" / "part1"
TRIO = ("PER1", "PER2", "NR1D1")
WF = cs.gtex()


def core(wf=WF):
    """The core loop of a waveform (its genes minus the trio; amplitudes and phases kept)."""
    return wf.drop([g for g in TRIO if g in wf.genes])


def trio(wf=WF):
    return wf.subset([g for g in wf.genes if g in TRIO])


def restrict(w, base=WF, keep=None, name=None):
    """Companion rule (A329 / A335): a companion waveform `w` drawn against `base` (same gene order) restricted to the
    positions of `keep` (default: base's core genes), with `keep`'s amplitudes and phases."""
    keep = core(base) if keep is None else keep
    idx = [k for k, g in enumerate(base.genes) if g in keep.genes]
    assert len(idx) == len(keep.genes)
    return Waveform(tuple(w.genes[k] for k in idx), keep.amp, keep.psi, None, None, f"{getattr(w, 'name', 'set')}_r")


def drops(wf=WF, n=None, seed=336):
    """[(label, waveform)] of the three-gene drops of wf; n: random subset of n drops plus the trio drop."""
    combos = list(itertools.combinations(wf.genes, 3))
    t = tuple(g for g in wf.genes if g in TRIO)
    if n is not None and n < len(combos):
        rng = np.random.default_rng(seed)
        pick = sorted(set(rng.choice(len(combos), n, replace=False).tolist()) | {combos.index(t)})
        combos = [combos[i] for i in pick]
    return [("+".join(c), wf.drop(list(c))) for c in combos], "+".join(t)


def sv(w):
    r = np.asarray(w.amp) / np.asarray(w.sd)
    return np.linalg.svd(np.c_[r * np.cos(w.psi), r * np.sin(w.psi)], compute_uv=False)


def drop_percentiles(values: dict, trio_label: str, waveforms: dict, direction: int = 1) -> dict:
    """values: drop label -> statistic (NaN = undefined); direction +1 if the original effect is positive (larger =
    stronger), -1 if negative. Percentile = share of valid drops weaker than the trio drop (x 100)."""
    lab = list(values)
    v = np.array([values[k] for k in lab], float) * direction
    ok = np.isfinite(v)
    t = v[lab.index(trio_label)]
    S = np.array([sv(waveforms[k]) for k in lab])
    st = sv(waveforms[trio_label])
    sim = ok & np.all(np.abs(S / st - 1) <= 0.2, axis=1)
    pct = float(np.mean(v[ok] < t) * 100)
    return dict(n_drops=int(ok.sum()), n_undefined=int((~ok).sum()), drop_pct_all=pct, n_similar=int(sim.sum()),
                drop_pct_similar=float(np.mean(v[sim] < t) * 100) if sim.sum() >= 20 else np.nan,
                drop_q10=float(np.quantile(v[ok], 0.10) * direction), drop_median=float(np.median(v[ok]) * direction),
                trio_drop_below_q10=bool(pct < 10))


def reading(core_status: str, below_q10) -> str:
    """Registered reading from the core status ("holds", "falls", "at the line") and the drop percentile flag."""
    if core_status == "holds":
        return "holds"
    if core_status == "at the line":
        return "at the line"
    if below_q10 is True:
        return "falls: carried by PER1, PER2 and NR1D1"
    return "falls: not decidable (precision)"


COLUMNS = ["id", "figure_panel", "statistic", "direction", "value_14_stored", "value_14", "reproduced", "value_core",
           "core_criterion", "core_status", "diff_14_minus_core", "diff_lo90", "diff_hi90", "diff_draws", "n_drops",
           "n_undefined", "drop_pct_all", "n_similar", "drop_pct_similar", "drop_q10", "drop_median",
           "trio_drop_below_q10", "trio_gate_passes", "trio_alone_value", "reading", "blind", "notes", "script"]


def write_rows(entry: str, rows: list):
    """Write the part-1 rows of one entry (dicts with COLUMNS keys; missing keys left empty)."""
    OUT.mkdir(parents=True, exist_ok=True)
    df = pd.DataFrame(rows)
    for c in COLUMNS:
        if c not in df:
            df[c] = np.nan
    df[COLUMNS].to_csv(OUT / f"{entry}.csv", index=False)
    return df[COLUMNS]
