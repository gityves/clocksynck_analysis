"""A352 (registered in `paper/docs/registered_experiments.md`): does the inflamed-skin delay scale with IL-17 / TNF
activity, and less with Th2 activity?  17 gated series of A347, unit = patient, lesional minus non-lesional.  Per patient
the signature score is the lesional minus non-lesional mean z-score (log expression z-scored within series|disease over
the paired samples) of the signature genes present on the platform (series dropped for a signature below 60 % of the
list).  Patients split into tertiles within series (T1 low); pooled paired offset per tertile (series as strata,
weights 1 / se^2 of the stored A53 table, as A347); trend = weighted slope of the pooled offset on the tertile index
(weights 1 / swap-null SD^2 of each tertile).  Readouts: 14-gene waveform and core loop.

  --stage power   swap nulls only (label swaps within patient, never the identity): SD of every tertile offset, of
                  T3 - T1 and of the trend, MDE = 2.49 SD, decidable = MDE(T3 - T1) <= 0.75 h; coverage per signature.
                  Signature scores and tertiles are computed (needed); no offset is computed on the real labels.
  --stage score   needs --go (user's go-ahead): tertile offsets, trend, one-sided permutation p (tertile labels permuted
                  within series, 2,000), matched / responsive companions (200 sets each, same permutations), A308-style q.
  --smoke         score path on lesional labels swapped at random within patient, reduced sizes (output: smoke/).

Seeds: swaps [352, 99, scope], permutations [352, scope, signature], companions [352, scope, family]
(scope 0 psoriasis, 1 AD, 2 both diseases).
D352 (GSE154200 replication row, Tsoi pre-mRNA by tertile) is not part of the power stage (see output notes).
"""
import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq
from scipy.stats import norm

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
sys.path.insert(0, str(PAPER / "exploratory/A347_inflamed_skin/scripts"))
import a347 as M  # noqa: E402  (loader, offsets, complex_pairs, wrap, wsd, swap helpers, A53 weights)
import clocksync as cs  # noqa: E402
from clocksync import paired as P  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402

T, C = M.T, M.C
OUT = HERE.parents[1] / "output"
SIGNATURES = {
    "IL17_TNF": ["IL17A", "IL17F", "IL17C", "IL36G", "IL36A", "DEFB4A", "S100A7A", "LCN2", "CXCL1", "CXCL8", "CCL20",
                 "PI3", "TNF", "TNFAIP6"],
    "Th2": ["IL4", "IL5", "IL13", "IL31", "CCL17", "CCL18", "CCL22", "CCL26", "POSTN", "CISH"],
}
SCOPES = {"psoriasis": 0, "atopic_dermatitis": 1, "both": 2}
MIN_FRAC, MIN_CELL, MDE_K, MDE_MAX = 0.6, 4, 2.49, 0.75
IDX = np.array([1.0, 2.0, 3.0])


def say(*a):
    print(time.strftime("%F %T"), *a, flush=True)


# ---- data ------------------------------------------------------------------------------------------------------

def load_series(dis, keys, b, all_genes):
    """Series -> (case, ref, counts); all_genes False reads only metadata, clock and signature genes."""
    sig = sorted({g for v in SIGNATURES.values() for g in v})
    out = {}
    for k in keys:
        path = b.SRC / f"{k}.parquet"
        if all_genes:
            d = T.read_all_genes(path, T.META)
        else:
            names = pq.ParquetFile(path).schema.names
            gc = T.clock_columns(names)
            keep = [c for c in T.META if c in names] + gc + [g for g in sig if g in names and g not in gc]
            d = pd.read_parquet(path, columns=list(dict.fromkeys(keep)))
        d = d[d.disease == dis]
        gc = T.clock_columns(list(d.columns))
        d = d.rename(columns=dict(zip(gc, T.GENES)))
        d = d.rename(columns={"patient": "unit", "status": "condition"}).assign(unit=lambda z: z["unit"].astype(str))
        case, ref = M.align(d)
        out[k] = (case.dropna(axis=1, how="all"), ref.dropna(axis=1, how="all"), bool((d.value_type == "counts").all()))
    return out


def log_sig(x, genes, counts):
    v = x[genes].to_numpy(float)
    return np.log2(v / x.depth.to_numpy(float)[:, None] * 1e6 + 1) if counts else v


def signature_scores(case, ref, counts, genes):
    """Per patient lesional minus non-lesional mean z of the genes present (None: < 60 % of the list)."""
    have = [g for g in genes if g in case.columns and g in ref.columns]
    if len(have) < MIN_FRAC * len(genes):
        return None, len(have)
    a, r = log_sig(case, have, counts), log_sig(ref, have, counts)
    both = np.vstack([a, r])
    sd = both.std(0)
    ok = sd > 0
    mu = both.mean(0)
    za, zr = (a[:, ok] - mu[ok]) / sd[ok], (r[:, ok] - mu[ok]) / sd[ok]
    return (za - zr).mean(1), len(have)


def tertiles(score):
    """0/1/2 tertile label within a series (ranks, ties by order)."""
    lab = np.empty(len(score), int)
    for t, ix in enumerate(np.array_split(np.argsort(score, kind="stable"), 3)):
        lab[ix] = t
    return lab


def design(scope, data, keys_by_dis):
    """(signature -> series -> tertile labels, coverage rows) for the series of the scope."""
    labs, cov = {}, []
    for sg, genes in SIGNATURES.items():
        labs[sg] = {}
        for dis, keys in keys_by_dis.items():
            for k in keys:
                case, ref, c = data[dis][k]
                s, n_have = signature_scores(case, ref, c, genes)
                n = len(case)
                kept = s is not None and n >= 3 * MIN_CELL
                cov.append(dict(scope_disease=dis, signature=sg, series=k, n_genes_listed=len(genes), n_genes_present=n_have,
                                n_patients=n, kept=kept,
                                why="" if kept else ("<60% genes" if s is None else f"n<{3 * MIN_CELL}")))
                if kept:
                    labs[sg][(dis, k)] = tertiles(s)
    return labs, cov


# ---- statistics ------------------------------------------------------------------------------------------------

def cell_offsets_swaps(pr, w, labs, F, n):
    """(n, 3) pooled swap-null offsets per tertile [h]; pr/labs/F keyed by (dis, series)."""
    out = np.empty((n, 3))
    for t in range(3):
        z = 0
        for key, lab in labs.items():
            i = np.flatnonzero(lab == t)
            a, va, b, vb = pr[key]
            z = z + w[key] * M.ccorr_draws(a[i], va[i], b[i], vb[i], F[key][:, i])
        out[:, t] = np.angle(z) * P.HOURS
    return out


def ccorr_masked(a, va, b, vb, mask):
    """`stats.ccorr` of (a, b) over the rows selected by mask (draws x persons, 0/1 float), vectorised."""
    cnt = mask.sum(1)
    A_ = a[None] - (mask @ a / cnt)[:, None]
    B_ = b[None] - (mask @ b / cnt)[:, None]
    k = cnt - 1
    ea = np.sum(mask * abs(A_) ** 2, 1) / k - mask @ va / cnt
    eb = np.sum(mask * abs(B_) ** 2, 1) / k - mask @ vb / cnt
    return np.sum(mask * A_ * np.conj(B_), 1) / k / np.sqrt(np.maximum(ea, 1e-12) * np.maximum(eb, 1e-12))


def cell_offsets_labels(pr, w, L, chunk=1000):
    """(draws, 3) pooled offsets per tertile [h] for tertile-label draws L[key] (draws x persons, values 0..2)."""
    n = len(next(iter(L.values())))
    out = np.empty((n, 3))
    for s in range(0, n, chunk):
        for t in range(3):
            z = 0
            for key in L:
                z = z + w[key] * ccorr_masked(*pr[key], (L[key][s:s + chunk] == t).astype(float))
            out[s:s + chunk, t] = np.angle(z) * P.HOURS
    return out


def slope(y, wt):
    """Weighted slope of y (draws x 3) on the tertile index, weights wt (3,)."""
    xb = (wt * IDX).sum() / wt.sum()
    yb = (y * wt).sum(1) / wt.sum()
    return ((wt * (IDX - xb)) * (y - yb[:, None])).sum(1) / (wt * (IDX - xb) ** 2).sum()


def rms(x):
    return float(np.sqrt(np.mean(np.asarray(x) ** 2)))


# ---- scopes ----------------------------------------------------------------------------------------------------

def scopes(data, keys_by_dis, wts):
    """scope -> (dict (dis, series) -> data key weights) for psoriasis, AD and both."""
    for sc, diseases in (("psoriasis", ["psoriasis"]), ("atopic_dermatitis", ["atopic_dermatitis"]),
                         ("both", ["psoriasis", "atopic_dermatitis"])):
        yield sc, diseases


def prep(scope_dis, data, keys_by_dis, wts, labs_by_sg, wf):
    pr = {(d, k): M.complex_pairs(*data[d][k], wf) for d in scope_dis for k in keys_by_dis[d]}
    w = {(d, k): wts[d][k] for d in scope_dis for k in keys_by_dis[d]}
    return pr, w


def power(n_swap=10000, smoke=False):
    t0 = time.time()
    b = M.builder()
    keys_by_dis, wts, data = {}, {}, {}
    for dis in M.DISEASES:
        keys_by_dis[dis], wts[dis] = M.gated(dis)
        data[dis] = load_series(dis, keys_by_dis[dis], b, all_genes=False)
    labs_all, cov = design(None, data, keys_by_dis)
    cov = pd.DataFrame(cov)
    rows = []
    for sc, diseases in scopes(data, keys_by_dis, wts):
        for sg in SIGNATURES:
            labs = {key: v for key, v in labs_all[sg].items() if key[0] in diseases}
            sizes = {key: len(v) for key, v in labs.items()}
            rng = np.random.default_rng([352, 99, SCOPES[sc], list(SIGNATURES).index(sg)])
            F = {key: rng.random((n_swap, m)) < 0.5 for key, m in sizes.items()}
            for ro, wf in M.READOUTS.items():
                pr = {key: M.complex_pairs(*data[key[0]][key[1]], wf) for key in labs}
                w = {key: wts[key[0]][key[1]] for key in labs}
                o = M.wrap(cell_offsets_swaps(pr, w, labs, F, n_swap))
                sd = np.array([rms(o[:, t] - 0.0) for t in range(3)])
                diff = M.wrap(o[:, 2] - o[:, 0])
                tw = 1 / sd ** 2
                tr = slope(o, tw)
                cnt = np.array([sum(int((v == t).sum()) for v in labs.values()) for t in range(3)])
                rows.append(dict(scope=sc, signature=sg, readout=ro, n_series=len(labs), n_patients=sum(sizes.values()),
                                 n_per_tertile=int(cnt[0]), n_swaps=n_swap,
                                 sd_T1_h=sd[0], sd_T2_h=sd[1], sd_T3_h=sd[2], sd_T3minusT1_h=rms(diff),
                                 sd_trend_h_per_tertile=rms(tr), mde_T3minusT1_h=MDE_K * rms(diff),
                                 mde_trend_h_per_tertile=MDE_K * rms(tr), mde_trend_span_h=MDE_K * 2 * rms(tr),
                                 decidable=bool(MDE_K * rms(diff) <= MDE_MAX)))
                say(sc, sg, ro, f"series {len(labs)} patients {sum(sizes.values())}  MDE(T3-T1) {rows[-1]['mde_T3minusT1_h']:.2f} h")
    out = OUT / ("power_smoke" if smoke else "power")
    out.mkdir(parents=True, exist_ok=True)
    pw = pd.DataFrame(rows)
    pw.to_csv(out / "power.csv", index=False)
    cov.to_csv(out / "coverage.csv", index=False)
    print(pw.round(3).to_string(index=False))
    print(cov.groupby(["scope_disease", "signature"]).agg(series=("kept", "sum"), listed=("series", "size"),
                                                          patients=("n_patients", lambda s: int(s[cov.loc[s.index, "kept"]].sum()))))
    say(f"runtime {time.time() - t0:.0f} s")


# ---- score -----------------------------------------------------------------------------------------------------

def score(smoke=False):
    t0 = time.time()
    n_perm, n_sets, n_swap = (200, 20, 500) if smoke else (2000, 200, 10000)
    out = OUT / ("smoke" if smoke else "score")
    out.mkdir(parents=True, exist_ok=True)
    b = M.builder()
    means = b.gene_means()
    keys_by_dis, wts, data = {}, {}, {}
    for d, dis in enumerate(M.DISEASES):
        keys_by_dis[dis], wts[dis] = M.gated(dis)
        data[dis] = load_series(dis, keys_by_dis[dis], b, all_genes=True)
        if smoke:
            sw = M.smoke_swap({(dis, k): v for k, v in data[dis].items()}, seed=d)
            data[dis] = {k: v for (_, k), v in sw.items()}
    labs_all, cov = design(None, data, keys_by_dis)
    pd.DataFrame(cov).to_csv(out / "coverage.csv", index=False)
    rows, hyp = [], []
    for sc, diseases in scopes(data, keys_by_dis, wts):
        for s_i, sg in enumerate(SIGNATURES):
            labs = {key: v for key, v in labs_all[sg].items() if key[0] in diseases}
            sizes = {key: len(v) for key, v in labs.items()}
            prng = np.random.default_rng([352, SCOPES[sc], s_i])
            L = {key: np.vstack([lab[None]] + [prng.permutation(lab)[None] for _ in range(n_perm)]) for key, lab in labs.items()}
            srng = np.random.default_rng([352, 99, SCOPES[sc], s_i])
            F = {key: srng.random((n_swap, m)) < 0.5 for key, m in sizes.items()}
            for ro, wf in M.READOUTS.items():
                pr = {key: M.complex_pairs(*data[key[0]][key[1]], wf) for key in labs}
                w = {key: wts[key[0]][key[1]] for key in labs}
                sd = np.array([rms(M.wrap(x)) for x in cell_offsets_swaps(pr, w, labs, F, n_swap).T])
                tw = 1 / sd ** 2
                clock = cell_offsets_labels(pr, w, L)
                tc = slope(M.wrap(clock), tw)
                p_neg = float((1 + np.sum(tc[1:] <= tc[0])) / (1 + n_perm))
                res = dict(scope=sc, signature=sg, readout=ro, n_series=len(labs), n_patients=sum(sizes.values()),
                           off_T1_h=clock[0, 0], off_T2_h=clock[0, 1], off_T3_h=clock[0, 2], trend_h=tc[0],
                           perm_sd_trend=float(tc[1:].std(ddof=1)), p_neg=p_neg, n_perm=n_perm)
                for f, fam in enumerate(["matched", "responsive"]):
                    crng = np.random.default_rng([352, SCOPES[sc], f])
                    per = {key: P.pseudo_waveforms(means[f"{key[1]}|{key[0]}"], n_sets, crng, M.WF14, fam) for key in labs}
                    Mx = np.empty((n_sets, n_perm + 1))
                    for j in range(n_sets):
                        wj = {key: C.restrict(per[key][j], base=M.WF14) if ro == "core" else per[key][j] for key in labs}
                        pj = {key: M.complex_pairs(*data[key[0]][key[1]], wj[key]) for key in labs}
                        Mx[j] = slope(M.wrap(cell_offsets_labels(pj, w, L)), tw)
                    sj = Mx[:, 1:].std(1, ddof=1)
                    sc_ = tc[1:].std(ddof=1)
                    tau2 = np.maximum(0.0, Mx.var(0, ddof=1) - np.mean(sj ** 2))
                    q = norm.cdf((tc - Mx.mean(0)) / np.sqrt(tau2 + sc_ ** 2))
                    res[f"q_{fam}"] = float(q[0])
                    res[f"sets_mean_trend_{fam}"] = float(Mx[:, 0].mean())
                    res[f"calibration_pass_rate_{fam}"] = float(np.mean(q[1:] < 0.05))
                rows.append(res)
                say(sc, sg, ro, "done")
    T_ = pd.DataFrame(rows)
    T_.to_csv(out / "trends.csv", index=False)
    print(T_.round(3).to_string(index=False), f"\nruntime {time.time() - t0:.0f} s")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=["power", "score"], required=True)
    ap.add_argument("--smoke", action="store_true", help="score on randomly swapped labels, reduced sizes")
    ap.add_argument("--go", action="store_true", help="real-label scoring: only with the user's go-ahead")
    a = ap.parse_args()
    if a.stage == "score" and not (a.smoke or a.go):
        sys.exit("score on the real labels needs --go (user's go-ahead); use --smoke to test the code path")
    limit_threads(4)
    if a.stage == "power":
        power(500 if a.smoke else 10000, a.smoke)
    else:
        score(a.smoke)


if __name__ == "__main__":
    main()
