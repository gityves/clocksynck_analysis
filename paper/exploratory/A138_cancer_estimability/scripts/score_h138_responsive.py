"""A138 note (2026-10-28, roadmap 5.2): H138.4-R, PER1/PER2-free pairing, immediate-early control and the pairing null
permuted within tissue source site / plate (TCGA tumour-normal pairs, 10 cancers).

Registered in paper/docs/registered_experiments.md, section "A109 / A138 note (2026-10-28 ...)". Secondary readouts: they
do not change the registered H138.4 verdict (`score.py`, `q3_lib.py`).
  H138.4-R  pooled (weights n_pairs, gated cancers) |complex correlation| of tumour and own-normal states vs the same
            statistic for 200 responsive-gene sets per cancer (glucocorticoid targets and immediate-early genes,
            expression-matched, carrying the waveform loadings; `pseudo_waveforms(kind="responsive")`), pooled per set
            index; one-sided add-one p < 0.05.
  (i)       PER1/PER2-free pairing (`cohort_states(drop_genes=)`) vs the registered matched pseudo null (200 sets of 14
            genes, as stored) and vs 200 matched sets of 12 genes carrying the reduced waveform; plus its permutation p.
  (ii)      immediate-early positive control (FOS, FOSB, JUN, EGR1, ATF3, IER2, DUSP1, NR4A1; the genes present):
            the pairing statistic of that set, whose 8 genes take the loadings of 8 randomly chosen waveform genes
            (200 random assignments; median and 5-95 % range reported).
  (iii)     pairing null permuted within cancer x tissue source site (TSS, field 2 of the TCGA barcode), within cancer x
            plate of the tumour sample and of the normal sample (field 6), 10,000 permutations, same streams for all.
  Reported  mean Jaccard overlap of the responsive sets (per cancer and overall).
The registered statistic is reproduced first (pooled statistic, matched pseudo median and q95, permutation null summaries
against output/power/h138_4_null.csv and output/score/tests.csv).

    python score_h138_responsive.py [--subsample N] [--jobs J] [--seed S]
Outputs: output/responsive/ (subsample/ for checks): reproduction_check.csv, h138_4_responsive.csv, per_cancer.csv,
stratified_null.csv, jaccard_responsive_sets.csv.  Holm over H109.1-R and H138.4-R: `holm_responsive.py`.
"""

from __future__ import annotations

import itertools
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

sys.path.insert(0, str(Path(__file__).resolve().parent))
import score as sc  # noqa: E402  (registered scoring: SEED, check_null, perm_p, paths)
import q3_lib as q  # noqa: E402
from clocksync.paired import cohort_states, pseudo_waveforms  # noqa: E402
from clocksync.waveform import Waveform  # noqa: E402
from clocksync._core.controls import RESPONSIVE_GENES  # noqa: E402
from common import cli  # noqa: E402

SEED = sc.SEED
HERE = sc.HERE
OUT = HERE / "output" / "responsive"
N_SETS = 200
DROP = ("PER1", "PER2")
IEG = ("FOS", "FOSB", "JUN", "EGR1", "ATF3", "IER2", "DUSP1", "NR4A1")
WF, WF12 = q.WF, q.WF.drop(DROP)
T0 = time.time()


def log(*a):
    print(f"[{(time.time() - T0) / 60:7.1f} min]", *a, flush=True)


def stat_of(ct, cn, dt, dn, w) -> float:
    g = list(w.genes)
    bt, vt, bn, vn = cohort_states(ct[g].assign(depth=dt), cn[g].assign(depth=dn), True, w)
    return q.corr_stat(q.centred_complex(bt), vt, q.centred_complex(bn), vn)


def one_cancer(c: str, i: int, seed: int, subsample, n_sets: int) -> dict:
    """All set-based statistics of one gated cancer (i: its index among the gated cancers, as `score.pseudo_pooled`)."""
    from clocksync.config import limit_threads
    limit_threads(1)
    samples = q.sh.load_tcga(subsample, seed)
    case, ref = q.sh.pair_frames(samples, c)
    ct, cn = q.sh.tcga_full_counts(c, case.index)
    lm = np.log(cn.div(cn.sum(axis=1), axis=0).mean() + 1e-12)[cn.mean() >= 1]
    dt, dn = ct.sum(axis=1), cn.sum(axis=1)
    out = dict(cancer=c)
    # registered matched pseudo null (same stream as score.pseudo_pooled) and its 12-gene companion
    ws = pseudo_waveforms(lm, n_sets, np.random.default_rng([SEED, 44, seed, i]), WF)
    out["matched14"] = [stat_of(ct, cn, dt, dn, w) for w in ws]
    ws = pseudo_waveforms(lm, n_sets, np.random.default_rng([SEED, 46, seed, i]), WF12)
    out["matched12"] = [stat_of(ct, cn, dt, dn, w) for w in ws]
    # responsive sets
    ws = pseudo_waveforms(lm, n_sets, np.random.default_rng([SEED, 45, seed, i]), WF, kind="responsive")
    out["responsive"] = [stat_of(ct, cn, dt, dn, w) for w in ws]
    gs = [set(w.genes) for w in ws]
    jac = [len(a & b) / len(a | b) for a, b in itertools.combinations(gs, 2)]
    out["jaccard"] = dict(cancer=c, n_pool_present=int(sum(g in lm.index for g in RESPONSIVE_GENES)), set_size=len(gs[0]),
                          n_sets=len(gs), mean_jaccard=float(np.mean(jac)), min_jaccard=float(np.min(jac)),
                          max_jaccard=float(np.max(jac)), genes_used=len(set().union(*gs)))
    # PER1/PER2-free registered statistic (the observed pairing frames, API option drop_genes)
    bt, vt, bn, vn = cohort_states(case, ref, True, WF, drop_genes=DROP)
    out["per_free"] = q.corr_stat(q.centred_complex(bt), vt, q.centred_complex(bn), vn)
    # immediate-early control: the IEG present take the loadings of random waveform genes
    ieg = [g for g in IEG if g in lm.index]
    rng = np.random.default_rng([SEED, 47, seed, i])
    wl = list(WF.genes)
    stats = []
    for _ in range(n_sets):
        idx = rng.choice(len(wl), len(ieg), replace=False)
        w = Waveform(tuple(ieg), WF.amp[idx], WF.psi[idx], None, None, "ieg")
        stats.append(stat_of(ct, cn, dt, dn, w))
    out["ieg"] = stats
    out["ieg_genes"] = ";".join(ieg)
    return out


def strat_labels(samples: pd.DataFrame, c: str, n: int) -> dict:
    case, ref = q.sh.pair_frames(samples, c)
    assert len(case) == n and (case.index == ref.index).all()
    bt, bn = case.barcode.astype(str), ref.barcode.astype(str)
    tss = bt.str.split("-").str[1]
    assert (tss == bn.str.split("-").str[1]).all()
    return dict(tss=tss.to_numpy(), plate_tumour=bt.str.split("-").str[5].to_numpy(),
                plate_normal=bn.str.split("-").str[5].to_numpy())


def strat_null(u: dict, labels: np.ndarray, n_perm: int, rng) -> np.ndarray:
    """Statistic with the rows of B permuted within the strata `labels` (each unit gets a normal of its own stratum;
    singleton strata keep their pairing). Same statistic and variances as `q.perm_null`."""
    za, vA, zb, vB = u["za"], u["vA"], u["zb"], u["vB"]
    n = len(za)
    ea, eb = q.tr_corrected(za, vA), q.tr_corrected(zb, vB)
    if ea <= 0 or eb <= 0:
        return np.full(n_perm, np.nan)
    code = pd.factorize(labels)[0].astype(float)
    base = np.argsort(code, kind="stable")
    perms = np.empty((n_perm, n), int)
    perms[:, base] = np.argsort(code[None, :] + rng.random((n_perm, n)) * 0.5, axis=1)
    return np.abs((za[None, :] * np.conj(zb[perms])).sum(1) / (n - 1)) / np.sqrt(ea * eb)


def main(argv=None):
    ap = cli.parser(__doc__)
    args = ap.parse_args(argv)
    run = cli.setup(args)
    quick = run["quick"]
    out = OUT / "subsample" if quick else OUT
    out.mkdir(parents=True, exist_ok=True)
    n_perm = 50 if quick else q.N_PERM
    n_sets = 20 if quick else N_SETS
    samples = q.sh.load_tcga(args.subsample, args.seed)
    st = q.tcga_states(args.subsample, args.seed)
    ref = pd.read_csv(sc.POWER / ("subsample/" if quick else "") / "h138_4_null.csv").set_index("scope")
    # registered statistic and permutation null (streams of score.py)
    obs, nulls = {}, {}
    for i, (c, u) in enumerate(st.items()):
        obs[c] = q.corr_stat(u["za"], u["vA"], u["zb"], u["vB"])
        nulls[c] = q.perm_null(u["za"], u["vA"], u["zb"], u["vB"], n_perm, np.random.default_rng([SEED, 4, args.seed, i]))
        sc.check_null(nulls[c], ref.loc[c], f"H138.4 {c}", quick)
    keep = [c for c in st if np.isfinite(nulls[c]).any() and st[c]["gated"]]
    n = [st[c]["n"] for c in keep]
    pooled_obs = float(q.pool(np.array([[obs[c]] for c in keep]), n)[0])
    pooled_null = q.pool(np.array([nulls[c] for c in keep]), n)
    sc.check_null(pooled_null, ref.loc["pooled_gated"], "H138.4 pooled", quick)
    log("registered pooled statistic", pooled_obs, "cancers", keep)

    res = Parallel(n_jobs=args.jobs)(delayed(one_cancer)(c, i, args.seed, args.subsample, n_sets) for i, c in enumerate(keep))
    pool = lambda key: q.pool(np.array([r[key] for r in res]), n)
    m14, m12, resp, ieg = pool("matched14"), pool("matched12"), pool("responsive"), pool("ieg")
    per_free = float(q.pool(np.array([[r["per_free"]] for r in res]), n)[0])

    # reproduction of the registered values
    stored = pd.read_csv(sc.HERE / "output" / "score" / ("subsample/" if quick else "") / "tests.csv")
    stored = stored[stored.hypothesis == "H138.4"].iloc[0]
    rep = pd.DataFrame([
        dict(quantity="pooled statistic", reproduced=pooled_obs, stored=stored.statistic),
        dict(quantity="permutation p", reproduced=sc.perm_p(pooled_null, pooled_obs), stored=stored.p),
        dict(quantity="matched pseudo median (200 sets)", reproduced=float(np.nanmedian(m14)), stored=stored.companion_median),
        dict(quantity="matched pseudo q95", reproduced=float(np.nanquantile(m14, 0.95)), stored=stored.companion_q95)])
    rep["match"] = [bool(np.isclose(a, b, rtol=1e-7, atol=1e-10)) for a, b in zip(rep.reproduced, rep.stored)]
    rep.to_csv(out / "reproduction_check.csv", index=False)
    print(rep.to_string(index=False), flush=True)
    if not quick:
        assert rep.match.all(), "registered H138.4 not reproduced"

    # permutation null of the PER-free statistic (same streams as the registered null)
    pf_null = []
    for i, c in enumerate(keep):
        bt, vt, bn, vn = cohort_states(*q.sh.pair_frames(samples, c), True, WF, drop_genes=DROP)
        pf_null.append(q.perm_null(q.centred_complex(bt), vt, q.centred_complex(bn), vn, n_perm,
                                   np.random.default_rng([SEED, 4, args.seed, list(st).index(c)])))
    pf_null = q.pool(np.array(pf_null), n)

    def line(test, stat, nul, note=""):
        nul = np.asarray(nul)[np.isfinite(nul)]
        return dict(test=test, statistic=stat, null_median=float(np.median(nul)), null_q95=float(np.quantile(nul, 0.95)),
                    excess=float(stat - np.median(nul)), p_one_sided=sc.perm_p(nul, stat), n_null=len(nul),
                    n_cancers=len(keep), n_pairs=int(sum(n)), note=note)

    rows = [
        line("H138.4 registered (permutation null)", pooled_obs, pooled_null, "10,000 permutations within cancer"),
        line("H138.4 registered vs matched pseudo-clock (14 genes)", pooled_obs, m14, "200 sets, A52 rule"),
        line("H138.4-R: pairing vs responsive sets", pooled_obs, resp, "200 sets per cancer, pooled per set index"),
        line("PER-free pairing vs registered matched null (14 genes)", per_free, m14, "200 sets"),
        line("PER-free pairing vs matched 12-gene sets", per_free, m12, "200 sets"),
        line("PER-free pairing vs permutation null", per_free, pf_null, "10,000 permutations within cancer"),
        line("registered pairing vs matched 12-gene sets (companion)", pooled_obs, m12, "200 sets")]
    ie = dict(test="immediate-early control pairing (positive control)", statistic=float(np.nanmedian(ieg)),
              null_median=float(np.nanmedian(m14)), null_q95=float(np.nanquantile(m14, 0.95)),
              excess=float(np.nanmedian(ieg) - np.nanmedian(m14)),
              p_one_sided=sc.perm_p(m14, float(np.nanmedian(ieg))), n_null=len(m14), n_cancers=len(keep), n_pairs=int(sum(n)),
              note=f"median of {len(ieg)} random loading assignments; 5-95% range {np.nanquantile(ieg, 0.05):.4f}-"
                   f"{np.nanquantile(ieg, 0.95):.4f}; null columns = matched pseudo-clock (14 genes); genes per cancer: "
                   + "|".join(f"{r['cancer']}:{r['ieg_genes']}" for r in res))
    rows.append(ie)
    pd.DataFrame(rows).to_csv(out / "h138_4_responsive.csv", index=False)

    pc = pd.DataFrame([dict(cancer=r["cancer"], n_pairs=n_, registered=obs[r["cancer"]], per_free=r["per_free"],
                            responsive_median=float(np.nanmedian(r["responsive"])),
                            responsive_q95=float(np.nanquantile(r["responsive"], 0.95)),
                            matched14_median=float(np.nanmedian(r["matched14"])),
                            matched12_median=float(np.nanmedian(r["matched12"])), ieg_median=float(np.nanmedian(r["ieg"])),
                            ieg_genes=r["ieg_genes"]) for r, n_ in zip(res, n)])
    pc.to_csv(out / "per_cancer.csv", index=False)
    jac = pd.DataFrame([r["jaccard"] for r in res])
    jac.loc[len(jac)] = dict(cancer="mean over cancers", **jac.drop(columns="cancer").mean().to_dict())
    jac.to_csv(out / "jaccard_responsive_sets.csv", index=False)
    log("responsive/PER-free/IEG done")

    # stratified pairing nulls
    lab = {c: strat_labels(samples, c, st[c]["n"]) for c in keep}
    srows, desc = [], []
    schemes = [("tss", "within cancer x TSS"), ("plate_tumour", "within cancer x tumour plate"),
               ("plate_normal", "within cancer x normal plate")]
    for si, (key, name) in enumerate(schemes):
        nl = {}
        for c in keep:
            i = list(st).index(c)
            lb = lab[c][key]
            nl[c] = strat_null(st[c], lb, n_perm, np.random.default_rng([SEED, 90 + si, args.seed, i]))
            vc = pd.Series(lb).value_counts()
            desc.append(dict(scheme=name, cancer=c, n_pairs=len(lb), n_strata=len(vc), n_singleton_strata=int((vc == 1).sum()),
                             share_pairs_in_multi_strata=float((vc[vc > 1]).sum() / len(lb)), max_stratum=int(vc.max())))
            srows.append(dict(scheme=name, scope=c, statistic=obs[c], null_mean=float(np.nanmean(nl[c])),
                              null_q95=float(np.nanquantile(nl[c], 0.95)),
                              registered_null_mean=float(np.nanmean(nulls[c])), p=sc.perm_p(nl[c], obs[c])))
        pn = q.pool(np.array([nl[c] for c in keep]), n)
        d = pd.DataFrame([x for x in desc if x["scheme"] == name])
        srows.append(dict(scheme=name, scope="pooled_gated", statistic=pooled_obs, null_mean=float(pn.mean()),
                          null_sd=float(pn.std(ddof=1)), null_q95=float(np.quantile(pn, 0.95)),
                          registered_null_mean=float(pooled_null.mean()), registered_null_q95=float(np.quantile(pooled_null, 0.95)),
                          p=sc.perm_p(pn, pooled_obs), n_perm=len(pn), n_strata_total=int(d.n_strata.sum()),
                          share_pairs_in_multi_strata=float((d.share_pairs_in_multi_strata * d.n_pairs).sum() / d.n_pairs.sum())))
        log(name, srows[-1])
    srows.insert(0, dict(scheme="within cancer (registered)", scope="pooled_gated", statistic=pooled_obs,
                         null_mean=float(pooled_null.mean()), null_sd=float(pooled_null.std(ddof=1)),
                         null_q95=float(np.quantile(pooled_null, 0.95)), registered_null_mean=float(pooled_null.mean()),
                         registered_null_q95=float(np.quantile(pooled_null, 0.95)), p=sc.perm_p(pooled_null, pooled_obs),
                         n_perm=len(pooled_null)))
    pd.DataFrame(srows).to_csv(out / "stratified_null.csv", index=False)
    pd.DataFrame(desc).to_csv(out / "stratified_strata.csv", index=False)
    print(pd.DataFrame(rows).drop(columns="note").round(4).to_string(index=False), flush=True)
    print(pd.DataFrame(srows)[lambda d: d.scope == "pooled_gated"].round(4).to_string(index=False), flush=True)
    log("done")


if __name__ == "__main__":
    main()
