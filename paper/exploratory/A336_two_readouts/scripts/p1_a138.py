"""A336 part 1, A138 H138.4 (tumour-normal pairing, Figure 6) and H138.8 (low minus high purity D of the pairing).

Original: `A138_cancer_estimability/scripts/{score.py, score_h138_responsive.py, h138_8.py, q3_lib.py}` (imported, not edited).
Same data, seeds (permutation nulls [138, 4, 0, i] with i the index among the cancers of the states; matched [138, 44, 0, i] and
responsive [138, 45, 0, i] sets with i the index among the 14-gene gated cancers; purity-label permutations [138, 8, 0, j]),
10,000 permutations, 200 + 200 sets. Only the waveform changes. Companions: the 14-gene draws restricted to the core
positions (A329 / A335 rule); the cancers pooled are the waveform's own gated, estimable ones (original rule; the 14-gene
index of each cancer is kept for the seeds). Gate = `reference_gate` of the normals on the waveform (rng [52, 0, i], 200).

H138.4 criterion (three conditions): one-sided permutation p < 0.05, above matched q95, above responsive q95 (holds: all;
at the line: some; falls: none). H138.8 ("D > 0, p < 0.05" not met; registered MDE 0.1624 = null q95 + 0.84 SD, 14 genes):
a null; core "holds" if still p >= 0.05 and the 90 % pair-bootstrap upper of D < 0.1624, else "falls". Bootstrap: 200
within-cancer pair resamples (rng [336, 1], shared by both waveforms; a stream per statistic: H138.4 first, then H138.8).
Trio-alone gate (paired design): the H138.4 pairing of the trio against its own permutation null, p < 0.05, on the gated
estimable trio cancers.

    python p1_a138.py [--check] [--jobs 4]
Output: ../output/part1/a138.csv
"""
import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
sys.path.insert(0, str(PAPER))
SC = PAPER / "exploratory/A138_cancer_estimability/scripts"
sys.path.insert(0, str(SC))
sys.path.insert(0, str(HERE.parent))
import q3_lib as q  # noqa: E402
import score as SCO  # noqa: E402
import h138_8 as H8  # noqa: E402
import a336_common as C  # noqa: E402
from clocksync import stats as S  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.paired import cohort_states, pseudo_waveforms, reference_gate  # noqa: E402

ENTRY = "a138"
BASE = q.WF
SEED = SCO.SEED
OUTA = PAPER / "exploratory/A138_cancer_estimability/output"
T4 = pd.read_csv(OUTA / "score/tests.csv").query("hypothesis == 'H138.4'").iloc[0]
RESP = pd.read_csv(OUTA / "responsive/h138_4_responsive.csv").set_index("test")
T8 = pd.read_csv(OUTA / "h138_8/tests.csv").iloc[0]
N_BOOT = 200
SCRIPT = "exploratory/A336_two_readouts/scripts/p1_a138.py"


def log(*a):
    print(time.strftime("%F %T"), *a, flush=True)


def states_wf(samples, wf):
    """As `q.tcga_states` with the waveform as argument: per cancer centred states, Poisson variances, n and gate."""
    g = list(wf.genes)
    out = {}
    for i, c in enumerate(q.sh.CANCERS):
        case, ref = q.sh.pair_frames(samples, c)
        if len(case) < 10:
            continue
        yr = np.log(np.maximum(ref[g].astype(float), 1e-300)).sub(np.log(ref.depth.to_numpy(float)), axis=0).to_numpy()
        gt = reference_gate(yr, wf, 200, np.random.default_rng([52, 0, i]))
        bt, vt, bn, vn = cohort_states(case, ref, True, wf)
        out[c] = dict(za=q.centred_complex(bt), vA=vt, zb=q.centred_complex(bn), vB=vn, n=len(case), gate_p=gt["p"],
                      gated=bool(gt["p"] < 0.05), raw=(bt, vt, bn, vn))
    return out


def obs_stat(st):
    """(per-cancer statistics, keep list) of the gated cancers with positive corrected variances (original keep rule)."""
    obs = {c: q.corr_stat(u["za"], u["vA"], u["zb"], u["vB"]) for c, u in st.items()}
    keep = [c for c in st if np.isfinite(obs[c]) and st[c]["gated"]]
    return obs, keep


def pooled_obs(st):
    obs, keep = obs_stat(st)
    if not keep:
        return np.nan, keep
    return float(q.pool(np.array([[obs[c]] for c in keep]), [st[c]["n"] for c in keep])[0]), keep


def perm_pass(st, n_perm):
    """Permutation nulls (original streams) and pooled null / p over the keep cancers."""
    obs, keep = obs_stat(st)
    nul = {}
    for i, (c, u) in enumerate(st.items()):
        if c in keep:
            nul[c] = q.perm_null(u["za"], u["vA"], u["zb"], u["vB"], n_perm, np.random.default_rng([SEED, 4, 0, i]))
    if not keep:
        return np.nan, np.array([np.nan]), np.nan, keep
    n = [st[c]["n"] for c in keep]
    po = float(q.pool(np.array([[obs[c]] for c in keep]), n)[0])
    pn = q.pool(np.array([nul[c] for c in keep]), n)
    return po, pn, SCO.perm_p(pn, po), keep


def one_cancer_sets(c, i, samples, transform, n_sets):
    """Pooled-statistic inputs of the matched and responsive sets of cancer c (seeds of the originals; i = 14-gene index)."""
    limit_threads(1)
    case, _ = q.sh.pair_frames(samples, c)
    ct, cn = q.sh.tcga_full_counts(c, case.index)
    lm = np.log(cn.div(cn.sum(axis=1), axis=0).mean() + 1e-12)[cn.mean() >= 1]
    dt, dn = ct.sum(axis=1), cn.sum(axis=1)
    res = {}
    for kind, stream in (("matched", 44), ("responsive", 45)):
        ws = pseudo_waveforms(lm, n_sets, np.random.default_rng([SEED, stream, 0, i]), BASE, kind=kind)
        v = []
        for w in ws:
            w = transform(w)
            g = list(w.genes)
            bt, vt, bn, vn = cohort_states(ct[g].assign(depth=dt), cn[g].assign(depth=dn), True, w)
            v.append(q.corr_stat(q.centred_complex(bt), vt, q.centred_complex(bn), vn))
        res[kind] = v
    return c, res


def sets_pass(samples, keep14, keep, n_by, transform, n_sets, jobs):
    """{kind: pooled statistic per set index} over `keep` (sets drawn for the 14-gene keep list with its indices)."""
    out = Parallel(n_jobs=jobs)(delayed(one_cancer_sets)(c, i, samples, transform, n_sets) for i, c in enumerate(keep14))
    d = dict(out)
    return {k: q.pool(np.array([d[c][k] for c in keep]), [n_by[c] for c in keep]) for k in ("matched", "responsive")}


def boot_idx(st, n_boot, rng):
    return {c: [rng.integers(0, u["n"], u["n"]) for _ in range(n_boot)] for c, u in st.items()}


def boot_h4(st, idx, keep):
    """Pooled statistic per bootstrap draw (pairs resampled within cancer, states re-centred)."""
    res = []
    for b in range(len(next(iter(idx.values())))):
        s, w = [], []
        for c in keep:
            bt, vt, bn, vn = st[c]["raw"]
            i = idx[c][b]
            s.append(q.corr_stat(q.centred_complex(bt[i]), vt[i], q.centred_complex(bn[i]), vn[i]))
            w.append(st[c]["n"])
        res.append(float(q.pool(np.array(s)[:, None], w)[0]))
    return np.array(res)


# ---- H138.8 ------------------------------------------------------------------------------------------------------------

def purity_inputs(samples, cancers):
    absp = H8.load_absolute()
    pa = {}
    for c in cancers:
        case, _ = q.sh.pair_frames(samples, c)
        pv = case.barcode.str[:15].map(absp)
        ok = pv.notna().to_numpy()
        pa[c] = dict(purity=pv.to_numpy(float)[ok], keep=ok)
    return pa


def split_D(sub, pur, use, i_by=None):
    """D = pooled(low) - pooled(high) for sub-cancers `sub` (a, b, va, vb arrays, purity) restricted to `use`."""
    lo, hi, hv = [], [], []
    for c in use:
        p = pur[c]
        low = (p <= np.median(p))[None, :]
        lo.append(sub[c].stat(low)[0]); hi.append(sub[c].stat(~low)[0])
        hv.append((int(low.sum()), int((~low).sum())))
    pl = float(q.pool(np.array(lo)[:, None], [h[0] for h in hv])[0])
    ph = float(q.pool(np.array(hi)[:, None], [h[1] for h in hv])[0])
    return pl - ph


def sub_of(st, pa, cancers):
    """Cancer objects restricted to the pairs with ABSOLUTE purity; cancers with >= 10 per half."""
    sub, use = {}, []
    for c in cancers:
        C0 = H8.Cancer(c, st[c])
        k = pa[c]["keep"]
        s = H8.Cancer.__new__(H8.Cancer)
        s.name, s.n = c, int(k.sum())
        s.a, s.b, s.va, s.vb = C0.a[k], C0.b[k], C0.va[k], C0.vb[k]
        sub[c] = s
        p = pa[c]["purity"]
        nl = int((p <= np.median(p)).sum())
        if min(nl, len(p) - nl) >= H8.MIN_HALF:
            use.append(c)
    return sub, use


def boot_D(sub, pa, use, idx):
    res = []
    for b in range(len(idx[use[0]])):
        s2 = {}
        for c in use:
            i = idx[c][b]
            o = H8.Cancer.__new__(H8.Cancer)
            o.name, o.n = c, len(i)
            o.a, o.b, o.va, o.vb = sub[c].a[i], sub[c].b[i], sub[c].va[i], sub[c].vb[i]
            s2[c] = o
        res.append(split_D(s2, {c: pa[c][idx[c][b]] for c in use}, use))
    return np.array(res)


def d_null_p(st, pa, keep_c, n_perm):
    """Registered H138.8: D, one-sided p (labels permuted within cancer, streams [138, 8, 0, j], j = index in keep_c)."""
    cidx = {c: i for i, c in enumerate(keep_c)}
    cs = {c: H8.Cancer(c, st[c]) for c in keep_c}
    pur = {c: dict(purity=pa[c]["purity"], keep=pa[c]["keep"]) for c in keep_c}
    r = H8.run_split(cs, pur, "absolute", n_perm, 0, cidx)
    x = r["null"][np.isfinite(r["null"])]
    return r["D"], float((1 + np.sum(x >= r["D"])) / (len(x) + 1)), r["use"], q.summarise(r["null"])


def trio_gate(st_trio, n_perm):
    po, pn, p, keep = perm_pass(st_trio, n_perm)
    return po, p, keep


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--jobs", type=int, default=4)
    a = ap.parse_args()
    chk = a.check
    limit_threads(1)
    n_perm, n_sets, n_boot = (50, 3, 5) if chk else (q.N_PERM, 200, N_BOOT)
    samples = q.sh.load_tcga()
    core, trio = C.core(BASE), C.trio(BASE)
    # ---- stage A: 14 genes, reproduce ------------------------------------------------------------------------------
    s14 = states_wf(samples, BASE)
    if chk:
        ref = q.tcga_states()
        for c in s14:
            assert np.allclose(s14[c]["za"], ref[c]["za"]) and s14[c]["gated"] == ref[c]["gated"], c
    po14, pn14, p14, keep14 = perm_pass(s14, n_perm)
    n14 = {c: s14[c]["n"] for c in s14}
    sets14 = sets_pass(samples, keep14, keep14, n14, lambda w: w, n_sets, a.jobs)
    pa = purity_inputs(samples, keep14)
    D14, pD14, use14, sm14 = d_null_p(s14, pa, keep14, n_perm)
    if not chk:
        rep = {"stat": np.isclose(po14, T4.statistic, atol=1e-6), "p": np.isclose(p14, T4.p, atol=1e-6),
               "matched_q95": np.isclose(np.quantile(sets14["matched"], 0.95), T4.companion_q95, atol=1e-6),
               "responsive_q95": np.isclose(np.quantile(sets14["responsive"], 0.95), RESP.loc["H138.4-R: pairing vs responsive sets", "null_q95"], atol=1e-6),
               "D": np.isclose(D14, T8.D, atol=1e-6), "pD": np.isclose(pD14, T8.p, atol=1e-6)}
        log("reproduction", {k: bool(v) for k, v in rep.items()})
        r4, r8 = rep["stat"] and rep["p"] and rep["matched_q95"] and rep["responsive_q95"], rep["D"] and rep["pD"]
        if not (r4 and r8):
            rows = []
            for i, stat, stor, ok, v in (("A138_H138.4", "pooled tumour-normal pairing (517 pairs)", T4.statistic, r4, po14),
                                         ("A138_H138.8", "D low minus high purity", T8.D, r8, D14)):
                rows.append(dict(id=i, figure_panel="Fig 6", statistic=stat, direction=1, value_14_stored=stor, value_14=v,
                                 reproduced=bool(ok), notes=None if ok else f"not reproduced: {rep}", script=SCRIPT))
            C.write_rows(ENTRY, rows)
            raise SystemExit(f"A138 not reproduced: {rep}")
    log("14-gene stage done")
    # ---- stage B: core ---------------------------------------------------------------------------------------------
    sc = states_wf(samples, core)
    poc, pnc, pc_, keepc = perm_pass(sc, n_perm)
    setsc = sets_pass(samples, keep14, keepc, {c: sc[c]["n"] for c in sc}, lambda w: C.restrict(w, BASE, core), n_sets, a.jobs)
    mq, rq = float(np.quantile(setsc["matched"], 0.95)), float(np.quantile(setsc["responsive"], 0.95))
    pac = purity_inputs(samples, keepc)
    Dc, pDc, usec, smc = d_null_p(sc, pac, keepc, n_perm)
    log("core done")
    # bootstrap, shared draws (H138.4 stream first, H138.8 second)
    rb = np.random.default_rng([336, 1])
    idx4 = boot_idx(s14, n_boot, rb)
    keep_b = [c for c in keep14 if c in keepc]
    # H138.4 interval of each waveform over its own keep list; the difference over the common cancers
    b14_own, bc_own = boot_h4(s14, idx4, keep14), boot_h4(sc, idx4, keepc)
    b14_com, bc_com = boot_h4(s14, idx4, keep_b), boot_h4(sc, idx4, keep_b)
    sub14, u14 = sub_of(s14, pa, keep14)
    subc, uc = sub_of(sc, pac, keepc)
    comm = [c for c in uc if c in u14]
    idx8 = {c: [rb.integers(0, sub14[c].n, sub14[c].n) for _ in range(n_boot)] for c in comm}
    bD14 = boot_D(sub14, {c: pa[c]["purity"] for c in comm}, comm, idx8)
    bDc = boot_D(subc, {c: pac[c]["purity"] for c in comm}, comm, idx8)
    # ---- gene drops (100 random + trio) ----------------------------------------------------------------------------
    dl, tl = C.drops(BASE, n=100)
    if chk:
        dl = [x for x in dl if x[0] == tl] + dl[:2]
    dv4, dv8, t0 = {}, {}, time.time()
    for n, (lab, w) in enumerate(dl):
        sd = states_wf(samples, w)
        dv4[lab] = pooled_obs(sd)[0]
        _, kd = obs_stat(sd)
        pad = purity_inputs(samples, kd) if kd else {}
        if kd:
            subd, ud = sub_of(sd, pad, kd)
            dv8[lab] = split_D(subd, {c: pad[c]["purity"] for c in ud}, ud) if ud else np.nan
        else:
            dv8[lab] = np.nan
        if n % 10 == 0:
            log("drops", n, len(dl), f"{(time.time() - t0) / 60:.1f} min")
    wfs = dict(dl)
    assert np.isclose(dv4[tl], poc, equal_nan=True) and np.isclose(dv8[tl], Dc, equal_nan=True)
    # ---- trio alone ------------------------------------------------------------------------------------------------
    st_t = states_wf(samples, trio)
    tpo, tp, tkeep = trio_gate(st_t, n_perm)
    gate = bool(len(tkeep) > 0 and tp < 0.05)
    tD = np.nan
    if gate:
        pat = purity_inputs(samples, tkeep)
        subt, ut = sub_of(st_t, pat, tkeep)
        tD = split_D(subt, {c: pat[c]["purity"] for c in ut}, ut) if ut else np.nan
    if chk:
        print("check ok")
        return
    d4 = b14_com - bc_com
    d8 = bD14 - bDc
    lo4, hi4 = np.quantile(bc_own, [0.05, 0.95])
    loD, hiD = np.quantile(bDc, [0.05, 0.95])
    rows = []
    conds = [bool(pc_ < 0.05), bool(poc > mq), bool(poc > rq)]
    st4 = "holds" if all(conds) else "at the line" if any(conds) else "falls"
    pr4 = C.drop_percentiles(dv4, tl, wfs, direction=1)
    rows.append(dict(
        id="A138_H138.4", figure_panel="Fig 6 (tumours follow own normal)", statistic=T4.hypothesis + " pooled tumour-normal |correlation| (gated cancers)",
        direction=1, value_14_stored=T4.statistic, value_14=po14, reproduced=True, value_core=poc,
        core_criterion="permutation p < 0.05 (10,000, within cancer), above matched q95 and above responsive q95 (companions restricted to core)",
        core_status=st4, diff_14_minus_core=po14 - poc, diff_lo90=float(np.quantile(d4, .05)), diff_hi90=float(np.quantile(d4, .95)),
        diff_draws=f"200 pair bootstraps within cancer (rng [336, 1]), shared, common cancers ({len(keep_b)})", **pr4,
        trio_gate_passes=gate, trio_alone_value=tpo if gate else np.nan, reading=C.reading(st4, pr4["trio_drop_below_q10"]),
        blind="core blind",
        notes=(f"core p {pc_:.4g}, null q95 {np.quantile(pnc, .95):.4f}; matched q95 {mq:.4f}, responsive q95 {rq:.4f}; cancers pooled core {len(keepc)} "
               f"({','.join(keepc)}) vs 14 genes {len(keep14)}; core 90 % bootstrap {lo4:.3f} to {hi4:.3f}; trio pairing p {tp:.4f} ({len(tkeep)} cancers gated)"),
        script=SCRIPT))
    mde = float(T8.mde)
    not_det = bool(pDc >= 0.05 or Dc <= 0)
    st8 = "holds" if (not_det and hiD < mde) else "falls"
    pr8 = C.drop_percentiles(dv8, tl, wfs, direction=1)
    rows.append(dict(
        id="A138_H138.8", figure_panel="Fig 6 (tumour-normal match vs purity)", statistic="D low minus high purity of the H138.4 pairing (ABSOLUTE purity)",
        direction=1, value_14_stored=T8.D, value_14=D14, reproduced=True, value_core=Dc,
        core_criterion=f"null: D > 0 not detected (p >= 0.05) and 90 % bootstrap upper of D < registered MDE {mde:.4f}",
        core_status=st8, diff_14_minus_core=D14 - Dc, diff_lo90=float(np.quantile(d8, .05)), diff_hi90=float(np.quantile(d8, .95)),
        diff_draws=f"200 pair bootstraps within cancer (rng [336, 1], second stream), shared, common cancers ({len(comm)})", **pr8,
        trio_gate_passes=gate, trio_alone_value=tD if gate else np.nan, reading=C.reading(st8, pr8["trio_drop_below_q10"]),
        blind="core blind",
        notes=(f"core p {pDc:.4f} (14 genes {pD14:.4f}); core null sd {smc['null_sd']:.4f}, own MDE {smc['mde']:.4f}; cancers used core {len(usec)}, 14 genes {len(use14)}; "
               f"core 90 % bootstrap {loD:.3f} to {hiD:.3f}"), script=SCRIPT))
    df = C.write_rows(ENTRY, rows)
    pd.set_option("display.width", 250)
    print(df.drop(columns=["script", "notes", "core_criterion", "diff_draws"]).round(4).T.to_string())
    for r in rows:
        print(r["id"], r["notes"])


if __name__ == "__main__":
    main()
