"""Shared helpers for the A115 / A116 feasibility scripts (core-vs-cue offset, one projector per series).

Within-group estimator: per series and group, c_core and c_cue are the clock-plane states of the 14-gene GTEx waveform and
of a cue set (one projector per series and set, built on all samples handed in); the group's core-minus-cue offset is
arg E[c_core c_cue*] (states centred within the group, Poisson-corrected normalisation). Pooling: arg of the sum over
series of w_s * complex correlation. Blinding: callers never pass healthy-control samples.
"""
from __future__ import annotations
import sys
from pathlib import Path
import numpy as np, pandas as pd, pyarrow.parquet as pq

PAPER = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PAPER.parent)); sys.path.insert(0, str(PAPER))
import clocksync as cs
from clocksync import stats as S
from clocksync import paired as PD
from common.paths import DATA

HOURS = S.HOURS
ROOT = DATA / "skin_lesional"
OUT54 = PAPER / "results/figure_5/output/psoriasis"
WF = cs.gtex()
CORE = list(WF.genes)
META = ["patient", "status", "sample", "disease", "platform", "value_type", "depth", "visit", "weeks", "treatment", "arm", "pasi"]


def cue_sets() -> dict:
    ref = pd.read_csv(PAPER / "results/figure_3/output/gtex_cues/a81b_cue_reference.csv", index_col=0)
    sysm = pd.read_csv(OUT54 / "a78_gene_sets.csv", index_col=0)
    sysm = sysm[sysm.set == "systemic"]
    out = {}
    for k in ("cue", "locked"):
        x = ref[ref.set == k]
        out[k] = cs.Waveform(tuple(x.index), x.amp.to_numpy(), x.acro_h.to_numpy() * 2 * np.pi / 24, None, None, k)
    out["systemic"] = cs.Waveform(tuple(sysm.index), sysm.amp_e.to_numpy() * np.log(2), sysm.acro_e.to_numpy() * 2 * np.pi / 24, None, None, "systemic")
    return out


def restrict(wf, genes):
    idx = [wf.genes.index(g) for g in genes]
    return cs.Waveform(tuple(genes), np.asarray(wf.amp)[idx], np.asarray(wf.psi)[idx], None, None, wf.name)


def all_set_genes(sets) -> list:
    return sorted(set(CORE) | {g for w in sets.values() for g in w.genes})


def read_series(path, genes, extra=()):
    """Metadata + the requested gene columns of one processed parquet (BMAL1 -> ARNTL)."""
    names = pq.ParquetFile(path).schema.names
    ren = {"BMAL1": "ARNTL"} if "ARNTL" not in names and "BMAL1" in names else {}
    inv = {v: k for k, v in ren.items()}
    cols = [c for c in list(META) + list(extra) if c in names]
    g = [inv.get(x, x) for x in genes if inv.get(x, x) in names]
    d = pd.read_parquet(path, columns=cols + g).rename(columns=ren)
    return d


def states_for(x: pd.DataFrame, wf, counts: bool):
    """(beta, v) of all rows of x for waveform wf (genes absent from x dropped), one projector for x."""
    g = [y for y in wf.genes if y in x.columns]
    w = wf.subset(g) if wf.name == "core" or set(wf.genes) == set(CORE) else restrict(wf, g)
    L = PD.projector(x[g].astype(float), w, counts=counts)
    return PD.states(x[g].astype(float), L, w, depth=x.depth.to_numpy(float) if counts else None), len(g)


def ccorr_sub(bc, vc, bu, vu, m) -> complex:
    return S.ccorr(bc[m], bu[m], vc[m], vu[m])


def pooled(cs_list, w) -> float:
    return float(np.angle(np.sum(np.asarray(w) * np.asarray(cs_list))) * HOURS)


def wrap(h) -> float:
    return float(S.wrap(h))


def weights54() -> dict:
    """series -> 1 / se^2 of the gated A54 series (discovery + held-out); cohort key 'GSE|psoriasis'."""
    w = {}
    for which in ("discovery", "heldout"):
        ps = pd.read_csv(OUT54 / f"a54_{which}_per_series.csv")
        ps = ps[ps.gated.astype(bool)]
        for s, se in zip(ps.series, ps.se):
            w[s] = 1 / se ** 2
    return w
