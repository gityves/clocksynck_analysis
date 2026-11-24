"""A367 (registered): is the person clock volume a product of the agonal phase? GTEx fast deaths split by Hardy class.

Hardy 1 (violent and fast, terminal phase < 10 min) against Hardy 2 (sudden natural, < 1 h).
- Leverage: A336 part 2 core-loop person leverage (27 peripheral tissues, 20 halvings, median), within a Hardy class;
  companions: A336's core matched / responsive sets (seeds [336, 2, 1] / [336, 2, 2]) within the same donors.
- GC slope: A361's model (A360 terms) on the core-loop volume (A360 rows, 25 tissues), fitted within a class.
Hardy 2 values, differences (donor bootstraps within class), age-matched Hardy 2 subsets and the GC-state level per
class are reported. Aggregates only: nothing per donor is written or printed.

    python a367_hardy.py --stage power [--jobs 8]   # counts, companions and covariates in Hardy 1; no clock value
    python a367_hardy.py --stage score [--jobs 8]   # registered tests (after the power note)
Outputs (`../output/`): power.csv, companions.csv | tests.csv, gc.csv.
"""
import argparse
import sys
import time
import warnings
from pathlib import Path

HERE = Path(__file__).resolve()
EXP = HERE.parents[2]
sys.path.insert(0, str(EXP / "A336_two_readouts" / "scripts"))
sys.path.insert(0, str(EXP / "A361_A363_volume_mechanism" / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

import a336_person_gtex as P  # noqa: E402
import vm_common as V  # noqa: E402
from a361_a363_systemic import A360_MODEL  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.paired import pseudo_waveforms  # noqa: E402
from clocksync.waveform import Waveform  # noqa: E402

OUT = HERE.parents[1] / "output"
LEV_WHOLE, LEV_MATCHED_Q95 = 0.714, 0.621          # A336 part 2 core (registered inputs, not blind)
GC_ALT = {"A361 core": -0.184, "A360 14 genes": -0.235}
N_PERM, N_SUB, N_BOOT, N_PERM_GC = 1000, 200, 2000, 10_000
# A367 addendum (2026-11-24): the user accepted H367.1 at power 0.775 (registered gate 0.8) after the power stage
H367_1_ACCEPTED = True


def leverage_rows():
    """A336 rows and per-donor Hardy class (kept in memory only)."""
    pb, donors = P.load()
    comps = sorted(pb.compartment.unique())
    donors.index = donors.index.astype(str)
    return pb, comps, donors


def companions(pb, comps, jobs):
    genes = [c for c in pb.columns if c not in P.META]
    x = pb[genes].to_numpy(np.float32)
    lm = pd.Series(np.log((x / pb.depth.to_numpy(float)[:, None]).mean(0) + 1e-12), index=genes)[x.mean(0) >= 1]
    lm = lm[~lm.index.duplicated()]
    del x
    wfs = []
    for kind, tag in (("matched", 1), ("responsive", 2)):
        for j, w in enumerate(pseudo_waveforms(lm, P.N_SETS, np.random.default_rng([336, 2, tag]), P.WF, kind=kind)):
            wfs.append(((kind, j), Waveform(tuple(w.genes[k] for k in P.CORE_IDX), P.CORE.amp, P.CORE.psi, None, None,
                                            f"core{j}")))
    return wfs


def setup(pb, comps, wfs):
    keep = sorted({g for _, w in wfs for g in w.genes} | set(P.WF.genes))
    P.S.update(pb=pb[P.META + keep].copy(), comps=comps)


def as_str(st):
    return [(np.asarray(j).astype(str), a, b) for j, a, b in st]


def sub(st, idx_of):
    """States restricted (with repeats allowed) to the donors whose per-halving row indices `idx_of` gives."""
    return [(j[i], a[i], b[i]) for (j, a, b), i in zip(st, idx_of)]


def rows_for(st, draw):
    pos = [{d: k for k, d in enumerate(j)} for j, _, _ in st]
    return [np.array([p[d] for d in draw if d in p], int) for p in pos]


def gc_power(per, h1):
    df = per.loc[per.index.intersection(h1)].dropna(subset=A360_MODEL)
    n, p = len(df), len(A360_MODEL)
    Z = (df[A360_MODEL] - df[A360_MODEL].mean()) / df[A360_MODEL].std(ddof=1)
    others = [t for t in A360_MODEL if t != "gc"]
    X = np.column_stack([np.ones(n), Z[others].to_numpy()])
    y = Z.gc.to_numpy()
    r2x = 1 - float(((y - X @ np.linalg.lstsq(X, y, rcond=None)[0]) ** 2).sum()) / float((y ** 2).sum())
    r2y = float(pd.read_csv(EXP / "A358_A360_gain_sources/output/a360_tests.csv").set_index("test").loc["full model R2", "partial_r2"])
    df_res = n - p - 1
    se = float(np.sqrt((1 - r2y) / (df_res * (1 - r2x))))
    tc = stats.t.ppf(0.95, df_res)
    out = dict(gc_donors_h1=n, gc_vif_h1=1 / (1 - r2x), gc_se_h1=se)
    for lab, b in GC_ALT.items():
        out[f"gc_power_h1|{lab}"] = float(stats.nct.sf(tc, df_res, abs(b) / se))
    return out


def gc_coef(per, v, ids):
    return float(V.coefs(per.loc[per.index.intersection(ids)], v, A360_MODEL).gc)


def main():
    warnings.simplefilter("ignore")
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=["power", "score"], required=True)
    ap.add_argument("--jobs", type=int, default=8)
    a = ap.parse_args()
    limit_threads(1)
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    pb, comps, donors = leverage_rows()
    in_rows = sorted(set(pb.donor.astype(str)))
    hardy = donors.hardy.reindex(in_rows)
    h1, h2 = list(hardy.index[hardy == 1]), list(hardy.index[hardy == 2])
    # GC rows (A360)
    pg, _, gc, ieg = V.a359.load()
    per = V.persons(pg, sorted(pg.compartment.unique()), gc, ieg)
    print(f"leverage rows: {len(h1)} Hardy 1, {len(h2)} Hardy 2; GC rows {len(per)} ({time.time() - t0:.0f}s)", flush=True)

    if a.stage == "power":
        wfs = companions(pb, comps, a.jobs)
        setup(pb, comps, wfs)
        del pb
        st = P.run(wfs, a.jobs)
        comp = pd.DataFrame([dict(kind=lab[0], set=lab[1], leverage_h1=P.lev(as_str(s), keep=h1)[0])
                             for (lab, _), s in zip(wfs, st)])
        comp.to_csv(OUT / "companions.csv", index=False)
        n1 = len(h1)
        n_eff = n1 - 9 - 3
        mq95 = float(comp[comp.kind == "matched"].leverage_h1.quantile(0.95))
        thr = max(mq95, 1.645 / np.sqrt(n_eff))
        zthr = np.arctanh(thr)
        res = dict(lev_donors_h1=n1, lev_donors_h2=len(h2), matched_median_h1=float(comp[comp.kind == "matched"].leverage_h1.median()),
                   matched_q95_h1=mq95, responsive_median_h1=float(comp[comp.kind == "responsive"].leverage_h1.median()),
                   responsive_q95_h1=float(comp[comp.kind == "responsive"].leverage_h1.quantile(0.95)),
                   lev_threshold_h1=thr, n_eff=n_eff,
                   lev_power_h1_whole=float(stats.norm.sf((zthr - np.arctanh(LEV_WHOLE)) * np.sqrt(n_eff - 3))),
                   lev_power_h1_half=float(stats.norm.sf((zthr - np.arctanh(LEV_WHOLE / 2)) * np.sqrt(n_eff - 3))))
        res.update(gc_power(per, h1))
        res["gate_H367.1"] = bool(res["lev_power_h1_whole"] >= 0.8)
        res["gate_H367.2"] = bool(res["gc_power_h1|A361 core"] >= 0.8)
        pd.Series(res).to_csv(OUT / "power.csv", header=["value"])
        print(pd.Series(res).round(4).to_string(), f"\ndone {time.time() - t0:.0f}s", flush=True)
        return

    # ---- score ---------------------------------------------------------------------------------------------------------
    pw = pd.read_csv(OUT / "power.csv", index_col=0)["value"]
    comp = pd.read_csv(OUT / "companions.csv")
    setup(pb, comps, [("core", P.CORE), ("14", P.WF)])
    del pb
    P.S["wfs"] = [("core", P.CORE), ("14 genes", P.WF)]
    sts = {lab: as_str(P.states(i)[1]) for i, (lab, _) in enumerate(P.S["wfs"])}
    rng = np.random.default_rng([367, 1])
    ages = donors.age.reindex(in_rows)
    rows = []
    for gs, st in sts.items():
        L1 = P.lev(st, keep=h1)[0]
        s1 = sub(st, rows_for(st, h1))
        null = np.array([P.lev([(j, a_, b_[rng.permutation(len(b_))]) for j, a_, b_ in s1])[0] for _ in range(N_PERM)])
        r = dict(gene_set=gs, lev_h1=L1, p_perm_h1=float((1 + np.sum(null >= L1)) / (1 + N_PERM)),
                 null_q95_h1=float(np.quantile(null, 0.95)), lev_h2=P.lev(st, keep=h2)[0])
        rs = np.random.default_rng([367, 3])
        r["lev_h2_at_n1"] = float(np.median([P.lev(st, keep=list(rs.choice(h2, len(h1), replace=False)))[0]
                                             for _ in range(N_SUB)]))
        # age-matched Hardy 2 subsets (bracket counts of Hardy 1; with replacement within a bracket if short)
        need = ages.reindex(h1).value_counts()
        pool = {b: [d for d in h2 if ages[d] == b] for b in need.index}
        am = []
        for _ in range(N_SUB):
            draw = [d for b, k in need.items() if pool[b] for d in rs.choice(pool[b], k, replace=len(pool[b]) < k)]
            am.append(P.lev(sub(st, rows_for(st, draw)))[0])
        r["lev_h2_age_matched"] = float(np.median(am))
        rb = np.random.default_rng([367, 4])
        diff = [P.lev(sub(st, rows_for(st, rb.choice(h2, len(h2)))))[0] - P.lev(sub(st, rows_for(st, rb.choice(h1, len(h1)))))[0]
                for _ in range(N_BOOT)]
        r.update(diff_h2_minus_h1=r["lev_h2"] - L1, diff_lo90=float(np.quantile(diff, 0.05)), diff_hi90=float(np.quantile(diff, 0.95)))
        if gs == "core":
            mt, rp = comp[comp.kind == "matched"].leverage_h1, comp[comp.kind == "responsive"].leverage_h1
            r.update(matched_q95_h1=float(mt.quantile(.95)), pct_matched_h1=float((mt < L1).mean() * 100),
                     responsive_q95_h1=float(rp.quantile(.95)), pct_responsive_h1=float((rp < L1).mean() * 100),
                     powered=str(pw["gate_H367.1"]) == "True" or H367_1_ACCEPTED)
            r["passes"] = bool(L1 > r["matched_q95_h1"] and r["p_perm_h1"] < 0.05) if r["powered"] else None
        rows.append(r)
    t = pd.DataFrame(rows)
    t.to_csv(OUT / "tests.csv", index=False)
    print(t.round(4).T.to_string(), flush=True)

    # GC slope (A360 rows)
    comps_g = sorted(pg.compartment.unique())
    hg = donors.hardy.reindex(per.index)
    g1, g2 = list(hg.index[hg == 1]), list(hg.index[hg == 2])
    grows = []
    for gs, w in (("core", V.CORE), ("14 genes", V.WF)):
        v = V.volume(pg, comps_g, w)
        c1, c2 = gc_coef(per, v, g1), gc_coef(per, v, g2)
        d1 = per.loc[per.index.intersection(g1)].assign(log_r=v.reindex(per.index)).dropna(subset=["log_r"] + A360_MODEL)
        rp = np.random.default_rng([367, 2])
        null = []
        for _ in range(N_PERM_GC):
            dp = d1.copy()
            dp["gc"] = rp.permutation(dp.gc.to_numpy())
            null.append(V.coefs(dp, dp.log_r, A360_MODEL).gc)
        rb = np.random.default_rng([367, 4])
        bd = []
        for _ in range(N_BOOT):
            out = []
            for g in (g2, g1):
                ids = rb.choice(g, len(g))
                d = per.loc[ids].assign(log_r=v.reindex(ids).to_numpy()).reset_index(drop=True).dropna(subset=["log_r"] + A360_MODEL)
                out.append(float(V.a360.std_coefs(d, A360_MODEL)[A360_MODEL.index("gc")]))
            bd.append(out[0] - out[1])
        gr = dict(gene_set=gs, gc_donors_h1=len(d1), gc_h1=c1, p_perm_gc_h1=float((1 + np.sum(np.array(null) <= c1)) / (1 + N_PERM_GC)),
                  gc_h2=c2, diff_h2_minus_h1=c2 - c1, diff_lo90=float(np.quantile(bd, 0.05)), diff_hi90=float(np.quantile(bd, 0.95)),
                  gc_state_mean_h1=float(per.gc.reindex(g1).mean()), gc_state_mean_h2=float(per.gc.reindex(g2).mean()),
                  gc_state_sd_h1=float(per.gc.reindex(g1).std()), gc_state_sd_h2=float(per.gc.reindex(g2).std()))
        if gs == "core":
            gr["powered"] = str(pw["gate_H367.2"]) == "True"
            gr["passes"] = bool(c1 < 0 and gr["p_perm_gc_h1"] < 0.05) if gr["powered"] else None
        grows.append(gr)
    g = pd.DataFrame(grows)
    g.to_csv(OUT / "gc.csv", index=False)
    print(g.round(4).T.to_string(), f"\ndone {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
