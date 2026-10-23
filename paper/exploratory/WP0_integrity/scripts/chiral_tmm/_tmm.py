"""TMM normalisation factors (edgeR `calcNormFactors(method = "TMM")`) and edgeR-style log2 CPM, in numpy.

Robinson and Oshlack 2010 as in edgeR defaults: logratioTrim 0.3, sumTrim 0.05, doWeighting TRUE, Acutoff -1e10 (zero
counts out). Reference sample = the one whose upper-quartile factor (75 % quantile of count / library size) is closest to
the mean of those factors (first on ties). Factors are rescaled to a geometric mean of 1.

    nf = tmm_factors(counts)                 # counts: samples x genes (as edgeR's columns x rows, transposed)
    L = log2_cpm(counts, nf, prior=0.25)     # edgeR cpm(log = TRUE, prior.count = 0.25) on DGEList(counts) + calcNormFactors

Test: `python _tmm.py` (hand example and invariances). R / edgeR is not installed in the sccc env or system R, so there is
no edgeR cross-check; the hand example below is written out independently of the vectorised code.
"""

from __future__ import annotations

import numpy as np
from scipy.stats import rankdata


def _factor(obs, ref, lib_o, lib_r, trim=0.3, sum_trim=0.05):
    """TMM factor of one sample against the reference (counts over genes)."""
    with np.errstate(divide="ignore", invalid="ignore"):
        log_r = np.log2((obs / lib_o) / (ref / lib_r))
        abs_e = (np.log2(obs / lib_o) + np.log2(ref / lib_r)) / 2
        v = (lib_o - obs) / lib_o / obs + (lib_r - ref) / lib_r / ref
    fin = np.isfinite(log_r) & np.isfinite(abs_e) & (abs_e > -1e10)
    log_r, abs_e, v = log_r[fin], abs_e[fin], v[fin]
    n = len(log_r)
    if n == 0 or np.max(np.abs(log_r)) < 1e-6:
        return 1.0
    lo_l, lo_s = int(np.floor(n * trim)) + 1, int(np.floor(n * sum_trim)) + 1
    hi_l, hi_s = n + 1 - lo_l, n + 1 - lo_s
    r_l, r_s = rankdata(log_r), rankdata(abs_e)
    keep = (r_l >= lo_l) & (r_l <= hi_l) & (r_s >= lo_s) & (r_s <= hi_s)
    f = np.sum(log_r[keep] / v[keep]) / np.sum(1 / v[keep])
    return 2.0 ** (f if np.isfinite(f) else 0.0)


def tmm_factors(counts, trim=0.3, sum_trim=0.05):
    """TMM normalisation factors (geometric mean 1) of the rows (samples) of `counts` (samples x genes)."""
    x = np.asarray(counts, float)
    lib = x.sum(1)
    q75 = np.array([np.quantile(r / l, 0.75) for r, l in zip(x, lib)])
    r = int(np.argmin(np.abs(q75 - q75.mean())))
    f = np.array([1.0 if i == r else _factor(x[i], x[r], lib[i], lib[r], trim, sum_trim) for i in range(len(x))])
    return f / np.exp(np.mean(np.log(f)))


def log2_cpm(counts, nf, prior=0.25):
    """edgeR log2 CPM: effective library = library x nf; the prior is scaled to each library and added twice to it."""
    x = np.asarray(counts, float)
    eff = x.sum(1) * nf
    ps = prior * eff / eff.mean()
    return np.log2((x + ps[:, None]) / (eff + 2 * ps)[:, None] * 1e6)


def _selftest():
    import math
    # 1. scale invariance and identical composition: every factor 1
    a = np.array([100., 200, 300, 400, 500, 600, 700, 800, 900, 1000])
    f = tmm_factors(np.vstack([a, 2 * a, 5 * a]))
    assert np.allclose(f, 1), f
    # 2. hand example: 10 genes, reference = sample 0 (see below), sample 1 has one gene 3x inflated
    ref = np.array([10., 20, 30, 40, 50, 60, 70, 80, 90, 100])
    obs = ref.copy(); obs[9] = 300.
    mid = ref * 1.0
    counts = np.vstack([ref, mid, obs])
    x = counts
    lib = x.sum(1)
    q = [np.quantile(r / l, .75) for r, l in zip(x, lib)]
    r = int(np.argmin([abs(v - np.mean(q)) for v in q]))
    assert r in (0, 1)  # ref and mid share the composition; the first is taken on ties of the distance
    # independent TMM of `obs` vs ref, plain python
    nO, nR = obs.sum(), x[r].sum()
    rows = []
    for o, c in zip(obs, x[r]):
        lr = math.log2((o / nO) / (c / nR)); ae = (math.log2(o / nO) + math.log2(c / nR)) / 2
        v = (nO - o) / nO / o + (nR - c) / nR / c
        rows.append((lr, ae, v))
    n = len(rows)                       # 10: loL = floor(3) + 1 = 4, hiL = 7; loS = floor(.5) + 1 = 1, hiS = 10
    by_lr = sorted(range(n), key=lambda i: rows[i][0])  # no ties except genes 0-8 (same logR): average ranks
    rk = {}
    i = 0
    while i < n:
        j = i
        while j + 1 < n and abs(rows[by_lr[j + 1]][0] - rows[by_lr[i]][0]) < 1e-12:
            j += 1
        for k in range(i, j + 1):
            rk[by_lr[k]] = (i + j) / 2 + 1
        i = j + 1
    keep = [g for g in range(n) if 4 <= rk[g] <= 7]  # absE ranks all inside 1..10
    f_hand = sum(rows[g][0] / rows[g][2] for g in keep) / sum(1 / rows[g][2] for g in keep)
    expect = 2 ** f_hand
    got = _factor(obs, x[r], nO, nR)
    assert abs(got - expect) < 1e-12, (got, expect)
    nf = tmm_factors(counts)
    assert abs(np.exp(np.mean(np.log(nf))) - 1) < 1e-12
    assert nf[2] < nf[0], nf  # an inflated gene raises the library: its factor is below the reference's
    # 3. log2 CPM: with nf = 1 and a single prior it equals the textbook formula
    L = log2_cpm(counts, np.ones(3))
    ps = 0.25 * lib / lib.mean()
    assert np.allclose(L, np.log2((counts + ps[:, None]) / (lib + 2 * ps)[:, None] * 1e6))
    print("tmm self-test passed; factors", np.round(nf, 4))


if __name__ == "__main__":
    _selftest()
