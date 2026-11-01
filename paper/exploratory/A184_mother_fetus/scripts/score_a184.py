"""A184 scoring (registered_experiments.md, commit f2445c1): does the term fetus follow its mother's internal time?

H184.1 (one-sided): coherence of (maternal blood, cord blood) pairs > its null from shuffling cord-blood states across
deliveries (5,000). Coherence = |complex correlation| of the paired clock-plane states over deliveries
(`paired.ccorr`, as the power step of prep_d2: `abs(ccorr(bA[j], vA[j], bB, vB))`; there is no separate per-delivery
unit-vector estimator). `clocksync.paired` (log2 arrays, GTEx waveform without PER3, one cohort) supplies the states,
the reference gate (maternal blood), offset (cord minus mother) and the pseudo-clock and responsive companions.
Companions: 200 expression-matched pseudo-clock and 200 responsive sets (means: mean log2 over the maternal-blood
samples), each scored on the same pairs, against the same 5,000 shuffles (the same permutation indices); clock-specific
only above the q95 of both companion coherences. Descriptives: offset with bootstrap interval, without PER1 and NR1D1.

    python score_a184.py [--n_perm 5000] [--n_pseudo 200] [--n_boot 1000]
"""
from __future__ import annotations

import argparse
import gzip
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

from clocksync import paired as P
from clocksync import stats as S
from clocksync.config import limit_threads

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[2] / "scouts_2026-11-01/prep_d2_gse27272/scripts"))
import prep_d2_gse27272 as PR  # noqa: E402  (probe_symbols, SCCC, OUT, SCALE)

OUT = HERE.parents[1] / "output"
SEED = 20261101


def say(*a):
    print(time.strftime("%F %T"), *a, flush=True)


def load():
    """Samples table (one row per delivery x {maternal, cord} blood, all symbols as columns, log2) and per-gene mean."""
    m = pd.read_csv(PR.OUT / "gse27272_samples.csv", dtype={"individual": str})
    with gzip.open(PR.SCCC / "GSE27272_series_matrix.txt.gz", "rt", errors="replace") as f:
        t = f.readlines()
    i0 = [i for i, l in enumerate(t) if l.startswith("!series_matrix_table_begin")][0]
    i1 = [i for i, l in enumerate(t) if l.startswith("!series_matrix_table_end")][0]
    ex = pd.read_csv(pd.io.common.StringIO("".join(t[i0 + 1:i1])), sep="\t", index_col=0).apply(pd.to_numeric, errors="coerce")
    sym = PR.probe_symbols()
    ok = ex.index.isin(sym.index)
    g = ex[ok].groupby(sym.reindex(ex.index[ok]).to_numpy()).mean()  # symbol x gsm, log2
    ct = pd.read_csv(PR.OUT / "clock_table.csv", dtype={"donor": str})  # the prep's clock genes: consistency check
    gc = [c for c in ct.columns if c not in ("donor", "compartment", "depth", "smoking")]
    mm = m.set_index(["individual", "compartment"]).gsm
    chk = np.log2(ct.set_index(["donor", "compartment"])[gc] / PR.SCALE)
    ref = g.loc[gc, mm.reindex(chk.index).to_numpy()].T.set_axis(chk.index)
    assert np.allclose(chk.to_numpy(), ref.to_numpy(), atol=1e-6), "expression differs from the prep table"
    m = m[m.compartment.isin(["maternal_blood", "cord_blood"])]
    both = m.groupby("individual").compartment.nunique() == 2
    m = m[m.individual.isin(both[both].index)]
    x = g[m.gsm].T
    s = pd.concat([pd.DataFrame({"unit": m.individual.to_numpy(), "condition": m.compartment.to_numpy()}, index=m.gsm), x], axis=1)
    means = g[m[m.compartment == "maternal_blood"].gsm].mean(axis=1)
    say("pairs", s.unit.nunique(), "symbols", x.shape[1])
    return s, means, gc


def states(s, wf):
    case = s[s.condition == "cord_blood"].set_index("unit")
    ref = s[s.condition == "maternal_blood"].set_index("unit")
    u = sorted(set(case.index) & set(ref.index))
    return P.cohort_states(case.loc[u], ref.loc[u], False, wf)


def coh(b, j=None):
    bA, vA, bB, vB = b
    if j is not None:
        bA, vA = bA[j], vA[j]
    return float(abs(P.ccorr(bA, vA, bB, vB)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n_perm", type=int, default=5000)
    ap.add_argument("--n_pseudo", type=int, default=200)
    ap.add_argument("--n_boot", type=int, default=1000)
    a = ap.parse_args()
    limit_threads(2)
    s, means, gc = load()
    wf = P.gtex().drop(["PER3"])
    assert list(wf.genes) == [g for g in wf.genes if g in gc] and len(wf.genes) == 13
    res = P.paired(s, unit="unit", condition="condition", levels=("maternal_blood", "cord_blood"), cohort=None, counts=False,
                   waveform=wf, pseudo_means={"all": means}, n_pseudo=a.n_pseudo, n_perm=200, seed=SEED, pseudo_responsive=True)
    res.cohorts.to_csv(OUT / "paired_cohort.csv", index=False)
    say(res.cohorts.round(3).T.to_string())
    b = states(s, wf)
    n = len(b[0])
    obs = coh(b)
    assert abs(obs - float(res.cohorts.coherence[0])) < 1e-9
    rng = np.random.default_rng(SEED)
    perms = np.array([rng.permutation(n) for _ in range(a.n_perm)])
    null = np.array([coh(b, j) for j in perms])
    p = S.perm_p(null, obs, "greater", tie=0.0)
    say(f"coherence {obs:.4f} null mean {null.mean():.4f} q95 {np.quantile(null, .95):.4f} p {p:.4g}")

    # companions on the same pairs and the same shuffles
    rows = []
    for kind in ("matched", "responsive"):
        crng = np.random.default_rng([SEED, 7 if kind == "matched" else 8])
        for j, w in enumerate(P.pseudo_waveforms(means, a.n_pseudo, crng, wf, kind)):
            bj = states(s, w)
            c = coh(bj)
            nj = np.array([coh(bj, q) for q in perms])
            rows.append(dict(kind=kind, set=j, coherence=c, null_q95=float(np.quantile(nj, .95)),
                             p=S.perm_p(nj, c, "greater", tie=0.0), offset_h=float(np.angle(P.ccorr(*bj)) * P.HOURS)))
        say("companion done", kind)
    comp = pd.DataFrame(rows)
    comp.to_csv(OUT / "companions_pseudo_sets.csv", index=False)
    t = dict(hypothesis="H184.1", n_pairs=n, coherence=obs, null_mean=float(null.mean()), null_q95=float(np.quantile(null, .95)),
             p=p, n_perm=a.n_perm, gate_p=float(res.cohorts.gate_p[0]))
    for kind in ("matched", "responsive"):
        v = comp[comp.kind == kind]
        t[f"{kind}_coherence_q95"] = float(np.quantile(v.coherence, .95))
        t[f"{kind}_coherence_median"] = float(np.median(v.coherence))
        t[f"{kind}_pctile_of_clock"] = float(np.mean(v.coherence <= obs))
        t[f"{kind}_frac_sets_p_lt_0.05"] = float(np.mean(v.p < 0.05))
        t[f"above_{kind}_q95"] = bool(obs > t[f"{kind}_coherence_q95"])
    t["passes_0.05"] = bool(p < 0.05)
    t["clock_specific"] = bool(t["passes_0.05"] and t["above_matched_q95"] and t["above_responsive_q95"])
    pd.DataFrame([t]).to_csv(OUT / "tests.csv", index=False)
    say("\n" + pd.Series(t).to_string())

    # descriptives: offset (cord minus mother) with bootstrap interval; without PER1 and NR1D1
    brng = np.random.default_rng([SEED, 9])
    rows = []
    for nme, drop in (("all_genes", []), ("without_PER1_NR1D1", ["PER1", "NR1D1"])):
        bb = states(s, wf.drop(drop) if drop else wf)
        z = P.ccorr(*bb)
        off = float(np.angle(z) * P.HOURS)
        bo = []
        for _ in range(a.n_boot):
            i = brng.integers(0, n, n)
            bo.append(np.angle(P.ccorr(*(x[i] for x in bb))) * P.HOURS)
        d = (np.asarray(bo) - off + 12) % 24 - 12
        lo, hi = np.quantile(d, [0.025, 0.975]) + off
        rows.append(dict(variant=nme, offset_h=off, boot_lo=float(lo), boot_hi=float(hi), coherence=float(abs(z)), n_boot=a.n_boot))
    pd.DataFrame(rows).to_csv(OUT / "descriptives.csv", index=False)
    say("\n" + pd.DataFrame(rows).round(3).to_string())
    say("done")


if __name__ == "__main__":
    main()
