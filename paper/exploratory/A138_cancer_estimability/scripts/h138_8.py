"""A138 note 6, H138.8: does the tumour-normal clock match (H138.4) come from the non-malignant fraction?

Registered in paper/docs/registered_experiments.md (A138 note 6, 2026-10-28). Reuses the data loading, projector, Poisson
correction, cancers, pairs and weights of H138.4 (`q3_lib.tcga_states`, `score.py`). Within each H138.4 cancer the pairs
are split at the median ABSOLUTE purity of the tumour sample; the H138.4 statistic is computed in each half (states
re-centred within the half, one projector per cancer as in H138.4); D = pooled(low) - pooled(high), weights = n pairs of
the half, cancers with >= 10 pairs per half. Null: purity labels permuted within cancer (10,000), p = (1 + #{D_null >= D})
/ 10,001.

    python h138_8.py [--subsample N] [--seed S]
Outputs (output/h138_8/): tests.csv, per_cancer.csv, pair_level.csv, estimate_sensitivity.csv
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
import q3_lib as q  # noqa: E402
from common import cli  # noqa: E402
from power import SEED  # noqa: E402

HERE = Path(__file__).resolve().parents[1]
OUT = HERE / "output" / "h138_8"
PUR_DIR = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian/tcga_purity_absolute")
ABS = PUR_DIR / "mastercalls.txt"
EST_GENES = PUR_DIR / "estimate_pkg" / "x" / "ncomms3612-s2.xlsx"
N_PERM = 10_000
MIN_HALF = 10
T0 = time.time()


def log(*a):
    print(f"[{(time.time() - T0) / 60:7.1f} min]", *a, flush=True)


def load_absolute() -> pd.Series:
    """ABSOLUTE purity by sample barcode (first 15 characters); only 'called' rows with a purity value."""
    t = pd.read_csv(ABS, sep="\t")
    t = t[(t["call status"] == "called") & t.purity.notna()]
    t = t.assign(key=t["sample"].str[:15])
    t = t.sort_values("sample").drop_duplicates("key")
    return t.set_index("key").purity


# ---- ESTIMATE (Yoshihara 2013) -------------------------------------------------------------------------------------

def ssgsea_raw(expr: np.ndarray, idx: np.ndarray, alpha: float = 0.25) -> np.ndarray:
    """Un-normalised ssGSEA enrichment (as the `estimate` package: rank-weighted running sum, summed) of one gene set for
    each row (sample) of `expr` (samples x genes); `idx` = column indices of the set."""
    g = expr.shape[1]
    in_set = np.zeros(g, bool)
    in_set[idx] = True
    out = np.empty(len(expr))
    for i, row in enumerate(expr):
        r = rankdata(row)
        o = np.argsort(-r, kind="stable")
        ins = in_set[o]
        w = np.where(ins, r[o] ** alpha, 0.0)
        p_gw = np.cumsum(w) / w.sum()
        p_ng = np.cumsum(~ins) / (~ins).sum()
        out[i] = np.sum(p_gw - p_ng)
    return out


N_UNIVERSE = 10_412  # genes common to the platforms of the `estimate` package; here the most highly expressed genes


def estimate_purity(counts: pd.DataFrame, sets: dict) -> pd.DataFrame:
    """Stromal / immune ssGSEA scores on the N_UNIVERSE genes with the highest mean count (the published purity formula
    was calibrated on a 10,412-gene universe; the raw score scales with it), purity = cos(0.6049872018 + 0.0001467884 x)."""
    counts = counts[counts.mean().sort_values(ascending=False).index[:N_UNIVERSE]]
    x = np.log2(counts.to_numpy(float) / counts.to_numpy(float).sum(1, keepdims=True) * 1e6 + 1)
    cols = {c: i for i, c in enumerate(counts.columns)}
    sc = {}
    for k, genes in sets.items():
        sc[k] = ssgsea_raw(x, np.array([cols[g] for g in genes if g in cols]))
    est = sc["stromal"] + sc["immune"]
    return pd.DataFrame({"stromal": sc["stromal"], "immune": sc["immune"], "estimate": est,
                         "purity": np.cos(0.6049872018 + 0.0001467884 * est)}, index=counts.index)


# ---- split statistic -----------------------------------------------------------------------------------------------

class Cancer:
    """States of one cancer (raw, H138.4 projector) and the split statistic for arbitrary low-masks."""

    def __init__(self, name: str, u: dict):
        from clocksync import stats as S
        bt, vt, bn, vn = u["raw"]
        self.name, self.n, self.u = name, u["n"], u
        self.a, self.b = S.complex_state(bt), S.complex_state(bn)
        self.va, self.vb = np.asarray(vt, float), np.asarray(vn, float)

    def stat(self, M: np.ndarray) -> np.ndarray:
        """H138.4 statistic on the rows selected by each mask of M (k x n, bool), states re-centred within the subset."""
        M = M.astype(float)
        n = M.sum(1)
        Sa, Sb = M @ self.a, M @ self.b
        cross = M @ (self.a * np.conj(self.b)) - Sa * np.conj(Sb) / n
        ea = (M @ np.abs(self.a) ** 2 - np.abs(Sa) ** 2 / n) / (n - 1) - (M @ self.va) / n
        eb = (M @ np.abs(self.b) ** 2 - np.abs(Sb) ** 2 / n) / (n - 1) - (M @ self.vb) / n
        ok = (ea > 0) & (eb > 0)
        with np.errstate(invalid="ignore", divide="ignore"):
            s = np.abs(cross / (n - 1)) / np.sqrt(ea * eb)
        return np.where(ok, s, np.nan)


def pooled(stats: dict, halves: dict, cancers: list, which: int) -> np.ndarray:
    """Pooled statistic of half `which` (0 low, 1 high) over `cancers`; arrays over permutations."""
    s = np.array([stats[c][which] for c in cancers])
    return q.pool(s, [halves[c][which] for c in cancers])


def split_stats(C: Cancer, purity: np.ndarray, perms: np.ndarray | None):
    """(stat_low, stat_high) for the observed labels (perms None) or per permutation of the purity labels."""
    med = np.median(purity)
    if perms is None:
        low = (purity <= med)[None, :]
    else:
        low = (purity <= med)[perms]
    return C.stat(low), C.stat(~low)


def run_split(cs: dict, pur: dict, label: str, n_perm: int, seed_tag: int, cancers_idx: dict):
    """D, null and per-cancer table for a purity vector per cancer (nan-free subsets already applied in `pur`)."""
    sub, obs, nul, halves = {}, {}, {}, {}
    for c, C in cs.items():
        p = pur[c]["purity"]
        keep = pur[c]["keep"]
        Cs = Cancer.__new__(Cancer)
        Cs.name, Cs.n = c, int(keep.sum())
        Cs.a, Cs.b, Cs.va, Cs.vb = C.a[keep], C.b[keep], C.va[keep], C.vb[keep]
        sub[c] = Cs
        med = np.median(p)
        nl = int((p <= med).sum())
        halves[c] = (nl, len(p) - nl)
    use = [c for c in cs if min(halves[c]) >= MIN_HALF]
    for c in use:
        p = pur[c]["purity"]
        rng = np.random.default_rng([SEED, 8, seed_tag, cancers_idx[c]])
        perms = np.argsort(rng.random((n_perm, len(p))), axis=1)
        obs[c] = tuple(float(x[0]) for x in split_stats(sub[c], p, None))
        lo, hi = split_stats(sub[c], p, perms)
        nul[c] = (lo, hi)
    st_obs = {c: obs[c] for c in use}
    po_l = float(pooled({c: (np.array([obs[c][0]]), np.array([obs[c][1]])) for c in use}, halves, use, 0)[0])
    po_h = float(pooled({c: (np.array([obs[c][0]]), np.array([obs[c][1]])) for c in use}, halves, use, 1)[0])
    D = po_l - po_h
    Dn = pooled(nul, halves, use, 0) - pooled(nul, halves, use, 1)
    return dict(use=use, halves=halves, obs=st_obs, D=D, pooled_low=po_l, pooled_high=po_h, null=Dn, sub=sub)


def tests_row(label: str, r: dict) -> dict:
    s = q.summarise(r["null"])
    x = r["null"][np.isfinite(r["null"])]
    p = float((1 + np.sum(x >= r["D"])) / (len(x) + 1))
    return dict(hypothesis=label, D=r["D"], pooled_low=r["pooled_low"], pooled_high=r["pooled_high"], p=p, null_mean=s["null_mean"],
                null_sd=s["null_sd"], null_q95=s["null_q95"], mde=s["mde"], n_perm=s["n_perm"], n_cancers=len(r["use"]),
                n_pairs_low=sum(r["halves"][c][0] for c in r["use"]), n_pairs_high=sum(r["halves"][c][1] for c in r["use"]),
                cancers=";".join(r["use"]),
                verdict=("pass" if (p < 0.05 and r["D"] > 0) else ("fail, MDE >= observed-scale: not estimable at this design" if False else "fail")))


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.set_defaults(seed=SEED)
    args = ap.parse_args(argv)
    run = cli.setup(args)
    quick = run["quick"]
    out = OUT / "subsample" if quick else OUT
    out.mkdir(parents=True, exist_ok=True)
    n_perm = 50 if quick else N_PERM

    # ---- reproduce H138.4 (pooled 0.3647, 517 pairs, 10 cancers) before any split
    st = q.tcga_states(args.subsample, 0)
    keep_c = [c for c in st if st[c]["gated"]]
    obs4 = {c: q.corr_stat(st[c]["za"], st[c]["vA"], st[c]["zb"], st[c]["vB"]) for c in keep_c}
    pooled4 = float(q.pool(np.array([[obs4[c]] for c in keep_c]), [st[c]["n"] for c in keep_c])[0])
    n4 = sum(st[c]["n"] for c in keep_c)
    log(f"H138.4 reproduced: pooled {pooled4:.4f}, {len(keep_c)} cancers, {n4} pairs")
    if not quick:
        assert abs(pooled4 - 0.3646870339954388) < 5e-4 and n4 == 517 and len(keep_c) == 10, "H138.4 not reproduced"

    samples = q.sh.load_tcga(args.subsample, 0)
    absp = load_absolute()
    sets = pd.read_excel(EST_GENES, sheet_name="gene", header=1)
    sets = {"stromal": sets[sets.Set == "Stromal141_UP"].Gene.tolist(), "immune": sets[sets.Set == "Immune141_UP"].Gene.tolist()}
    cidx = {c: i for i, c in enumerate(keep_c)}
    cs, pa, pe, rows, pair_rows, est_rows = {}, {}, {}, [], [], []
    for c in keep_c:
        case, _ = q.sh.pair_frames(samples, c)
        key = case.barcode.str[:15]
        pv = key.map(absp)
        ct, _ = q.sh.tcga_full_counts(c, case.index)
        ct = ct.loc[case.index]
        ep = estimate_purity(ct, sets)
        C = Cancer(c, st[c])
        assert abs(C.stat(np.ones((1, len(case)), bool))[0] - obs4[c]) < 1e-9, "split statistic differs from H138.4"
        cs[c] = C
        ok = pv.notna().to_numpy()
        pa[c] = dict(purity=pv.to_numpy(float)[ok], keep=ok)
        ok2 = np.isfinite(ep.purity.to_numpy())
        pe[c] = dict(purity=ep.purity.to_numpy(float), keep=ok2)
        sp = spearmanr(pv[ok], ep.purity[ok]) if ok.sum() > 3 else None
        est_rows.append(dict(cancer=c, n_pairs=len(case), n_with_absolute=int(ok.sum()), spearman_absolute_vs_estimate=float(sp.statistic) if sp else np.nan,
                             median_estimate_purity=float(ep.purity.median()), median_absolute_purity=float(pv[ok].median()) if ok.any() else np.nan,
                             spearman_stromal_score_vs_absolute=float(spearmanr(ep.stromal[ok], pv[ok]).statistic) if ok.sum() > 3 else np.nan,
                             spearman_immune_score_vs_absolute=float(spearmanr(ep.immune[ok], pv[ok]).statistic) if ok.sum() > 3 else np.nan))
        # pair-level contribution, whole cancer (states of H138.4: centred over all pairs of the cancer)
        u = st[c]
        za, zb = u["za"], u["zb"]
        theta = np.angle(np.sum(za * np.conj(zb)))
        ea, eb = q.tr_corrected(za, u["vA"]), q.tr_corrected(zb, u["vB"])
        contrib = np.real(za * np.conj(zb) * np.exp(-1j * theta)) / (np.sqrt(ea) * np.sqrt(eb)) if ea > 0 and eb > 0 else np.full(len(za), np.nan)
        s_abs = spearmanr(pv.to_numpy(float)[ok], contrib[ok]) if ok.sum() > 3 else None
        s_est = spearmanr(ep.purity.to_numpy(float), contrib)
        pair_rows.append(dict(cancer=c, n_pairs=len(case), n_with_absolute=int(ok.sum()), spearman_absolute_purity_contribution=float(s_abs.statistic) if s_abs else np.nan,
                              p_two_sided=float(s_abs.pvalue) if s_abs else np.nan, spearman_estimate_purity_contribution=float(s_est.statistic),
                              p_two_sided_estimate=float(s_est.pvalue), theta_c_h=float(theta * q.sh.HOURS)))
        log(c, "absolute match", int(ok.sum()), "/", len(case))

    # primary: ABSOLUTE
    r = run_split(cs, pa, "absolute", n_perm, 0, cidx)
    t = tests_row("H138.8 (ABSOLUTE purity)", r)
    log("H138.8:", {k: t[k] for k in ("D", "p", "null_mean", "null_sd", "mde")})
    # ESTIMATE sensitivity (descriptive)
    re = run_split(cs, pe, "estimate", n_perm, 1, cidx)
    te = tests_row("H138.8 sensitivity (ESTIMATE purity; descriptive)", re)
    te["verdict"] = "descriptive"
    log("ESTIMATE sensitivity:", {k: te[k] for k in ("D", "p", "null_mean", "null_sd", "mde")})
    tt = pd.DataFrame([t, te])
    # tie-aware MDE reading
    tt.loc[0, "verdict"] = ("pass: D > 0, p < 0.05" if (t["p"] < 0.05 and t["D"] > 0) else
                            ("fail, MDE large (>= H138.4 pooled statistic scale?): see mde" if t["mde"] > pooled4 else "fail, MDE small relative to H138.4 pooled statistic"))
    tt.to_csv(out / "tests.csv", index=False)

    pc = []
    for c in keep_c:
        p = pa[c]["purity"]
        row = dict(cancer=c, n_pairs=len(cs[c].a), n_with_purity=int(len(p)), match_rate=len(p) / len(cs[c].a),
                   median_purity=float(np.median(p)) if len(p) else np.nan, used_in_pool=c in r["use"], statistic_all_pairs=obs4[c])
        if c in r["halves"]:
            row |= dict(n_low=r["halves"][c][0], n_high=r["halves"][c][1])
        if c in r["obs"]:
            row |= dict(stat_low=r["obs"][c][0], stat_high=r["obs"][c][1], D=r["obs"][c][0] - r["obs"][c][1])
        if c in re["obs"]:
            row |= dict(estimate_stat_low=re["obs"][c][0], estimate_stat_high=re["obs"][c][1], estimate_D=re["obs"][c][0] - re["obs"][c][1])
        pc.append(row)
    pd.DataFrame(pc).to_csv(out / "per_cancer.csv", index=False)
    pd.DataFrame(pair_rows).to_csv(out / "pair_level.csv", index=False)
    es = pd.DataFrame(est_rows)
    es = pd.concat([es, pd.DataFrame([te | dict(cancer="pooled_split_test")])], ignore_index=True)
    es.to_csv(out / "estimate_sensitivity.csv", index=False)
    pd.DataFrame([dict(h138_4_pooled_reproduced=pooled4, n_pairs=n4, n_cancers=len(keep_c))]).to_csv(out / "h138_4_reproduction.csv", index=False)
    log("done")


if __name__ == "__main__":
    main()
