"""A336 part 1, Figure 6: A254 H254.1 (pooled lesion offset, time kept) and H254.3 (pooled clock-plane variance ratio r).

Original: `A254_paired_lesions/scripts/compute_a254.py` (imported, not edited). Same data, seeds ([254, 1, key] gate, [254, 2, key]
swap nulls, [254, 6] bootstraps, [254, 3] pooled offset null, [254, tag, key] companions with tag 4 matched / 5 responsive,
200 sets, 10,000 swaps); only the waveform changes. Companions: the original 14-gene draws restricted to the core positions
of each dataset (A329 rule). The gate (`P.reference_gate`) is re-run on the waveform read (14 genes, core, drop, trio); a
dataset enters the pooled statistics if its gate p < 0.05 and it has >= 20 pairs, as originally.

    python p1_a254.py [--check]     # --check: tiny run, assertions only, nothing printed about values nor written
Output: ../output/part1/a254.csv.
"""
import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(PAPER / "exploratory/A254_paired_lesions/scripts"))
sys.path.insert(0, str(HERE.parent))
import compute_a254 as A  # noqa: E402
import a336_common as H  # noqa: E402
from clocksync import paired as P  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402

ENTRY = "a254"
SCRIPT = "paper/exploratory/A336_two_readouts/scripts/p1_a254.py"
ST = pd.read_csv(PAPER / "exploratory/A254_paired_lesions/output/score_pooled.csv").iloc[0]
MARGIN_H = 1.0   # registered detectability limit of H254.1 (not scored if the pooled MDE exceeds 1.0 h)
KINDS = (("matched", 4), ("responsive", 5))


def log(*a):
    print(time.strftime("%F %T"), *a, flush=True)


def gate_for(d, wf):
    yr = P.log_expr(d.ref[list(wf.genes)], d.ref.depth.to_numpy(float) if d.counts else None).to_numpy()
    return P.reference_gate(yr, wf, 200, np.random.default_rng([A.SEED, 1, sum(map(ord, d.key))]))["p"]


def present(d, wf):
    g = [x for x in wf.genes if x in d.case.columns]
    return wf.subset(g)


def clock_stat(d, wf):
    """(gated, states) for waveform wf on dataset d (genes absent from the dataset dropped)."""
    w = present(d, wf)
    return bool(gate_for(d, w) < 0.05 and d.n >= 20), d.states(w), w


def full(D, wfs, n_perm, n_sets, n_boot, comps, with_null=True):
    """All H254.1 / H254.3 quantities for one waveform spec. wfs: dataset -> waveform; comps: kind -> dataset -> list of waveforms."""
    gated, B, W = {}, {}, {}
    for k, d in D.items():
        gated[k], B[k], W[k] = clock_stat(d, wfs[k])
    use = [k for k in D if gated[k]]
    out = dict(use=use, B=B, gated=gated)
    if not use:
        return out
    rng6 = np.random.default_rng([A.SEED, 6])
    nul_r, r, se = {}, {}, {}
    for k, d in D.items():
        n = d.n
        if with_null:
            r2 = np.random.default_rng([A.SEED, 2, sum(map(ord, k))])
            P.swap_offset_null([B[k]], [1.0], n_perm=n_perm, rng=r2)   # drawn first, as in the original (stream kept aligned)
            nul_r[k] = P.swap_null([B[k]], [1.0], n_perm=n_perm, rng=r2)
        for _ in range(n_boot):          # the original's offset bootstrap draws (stream kept aligned)
            rng6.integers(0, n, n)
        r[k] = P.log_var_ratio(*B[k])
        se[k] = P.boot_se(*B[k], n_boot=n_boot, rng=rng6)
    w = np.array([D[k].n for k in use], float)
    obs = P.pooled_offset([P.ccorr(*B[k]) for k in use], w)
    pooled_null = P.swap_offset_null([B[k] for k in use], w, n_perm=n_perm, rng=np.random.default_rng([A.SEED, 3])) if with_null else np.array([np.nan])
    out.update(w=w, pooled_offset_h=obs, null_sd_h=float(np.std(pooled_null)), p_delay=A.p_low(pooled_null, obs), mde_h=A.MDE_K * float(np.std(pooled_null)))
    t = pd.DataFrame(dict(r=[r[k] for k in use], se=[se[k] for k in use]), index=use)
    ok = np.isfinite(t.r) & (t.se > 0)
    iv = 1 / t.se.to_numpy() ** 2
    iv = np.where(ok, iv, 0.0)
    out.update(iv=iv, r=r, se=se, nul_r=nul_r)
    pr = float(np.sum(iv * t.r.to_numpy()) / iv.sum())
    out["pooled_r"] = pr
    out["pooled_r_se"] = float(1 / np.sqrt(iv.sum()))
    if with_null:
        nr = sum(wi * (nul_r[k] - np.nanmean(nul_r[k])) for wi, k in zip(iv / iv.sum(), use))
        out["p_r"] = A.p_two(nr, pr - float(np.sum(iv / iv.sum() * [np.nanmean(nul_r[k]) for k in use])))
    if comps is not None:
        for kind, sets in comps.items():
            po = np.array([P.pooled_offset([P.ccorr(*D[k].states(sets[k][j])) for k in use], w) for j in range(n_sets)])
            v = np.array([np.nansum(iv / iv.sum() * [P.log_var_ratio(*D[k].states(sets[k][j])) for k in use]) for j in range(n_sets)])
            out[f"{kind}_q05_h"] = float(np.nanquantile(po, .05))
            out[f"{kind}_r_q05"], out[f"{kind}_r_q95"] = float(np.quantile(v, .05)), float(np.quantile(v, .95))
    return out


def drop_stats(D, wf):
    """pooled offset and pooled r of a waveform, no companions, no nulls (the gate and the original's bootstrap SE are kept)."""
    gated, B = {}, {}
    for k, d in D.items():
        gated[k], B[k], _ = clock_stat(d, wf)
    use = [k for k in D if gated[k]]
    if not use:
        return np.nan, np.nan, 0
    w = np.array([D[k].n for k in use], float)
    off = P.pooled_offset([P.ccorr(*B[k]) for k in use], w)
    rng6 = np.random.default_rng([A.SEED, 6])
    se_all = {}
    for k in D:                       # the original's stream: offset-bootstrap draws, then the SE draws, dataset by dataset
        for _ in range(A.N_BOOT):
            rng6.integers(0, D[k].n, D[k].n)
        se_all[k] = P.boot_se(*B[k], n_boot=A.N_BOOT, rng=rng6)
    se = np.array([se_all[k] for k in use])
    r = np.array([P.log_var_ratio(*B[k]) for k in use])
    iv = np.where(np.isfinite(r) & (se > 0), 1 / se ** 2, 0.0)
    return off, float(np.sum(iv * np.nan_to_num(r)) / iv.sum()) if iv.sum() > 0 else np.nan, len(use)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    chk = ap.parse_args().check
    limit_threads(4)
    t0 = time.time()
    n_perm, n_sets, n_boot = (30, 2, 20) if chk else (A.N_PERM, A.N_SETS, A.N_BOOT)
    D = {k: A.Pair(k) for k in A.TESTED}
    WF = H.WF
    # ---- stage A: 14 genes
    draws = {kind: {k: D[k].pseudo(kind, n_sets, tag) for k in D} for kind, tag in KINDS}
    wf14 = {k: D[k].wf for k in D}
    comp14 = draws
    T14 = full(D, wf14, n_perm, n_sets, n_boot, comp14)
    log("stage A done")
    exp = dict(pooled_offset_h=ST.pooled_offset_h, null_sd_h=ST.null_sd_h, p_delay=ST.p_delay, matched_q05_h=ST.matched_q05_h,
               responsive_q05_h=ST.responsive_q05_h, pooled_r=ST.pooled_r, pooled_r_se=ST.pooled_r_se, p_r=ST.p_r_two_sided,
               matched_r_q05=ST.matched_r_q05, responsive_r_q05=ST.responsive_r_q05, mde_h=ST.mde_h)
    bad = {k: (T14.get(k), v) for k, v in exp.items() if k not in T14 or not np.isclose(T14[k], v, rtol=0, atol=1e-6)}
    ok = not bad and ";".join(T14["use"]) == ST.datasets
    if chk:
        assert len(T14["use"]) >= 1
    else:
        log("reproduced" if ok else f"NOT reproduced {bad}")
    if not chk and not ok:
        rows = [dict(id=i, figure_panel="Fig 6", statistic=s, direction=dr, value_14_stored=v, value_14=T14.get(c), reproduced=False,
                     notes=f"14 genes did not reproduce score_pooled.csv; no core value computed. {bad}; datasets used {T14['use']}", script=SCRIPT)
                for i, s, dr, v, c in (("A254_H254.1", "pooled lesion offset", "-1", ST.pooled_offset_h, "pooled_offset_h"),
                                       ("A254_H254.3", "pooled log variance ratio r", "-1", ST.pooled_r, "pooled_r"))]
        H.write_rows(ENTRY, rows)
        return
    # ---- stage B
    core, trio = H.core(WF), H.trio(WF)
    wfc = {k: H.core(D[k].wf) for k in D}
    compc = {kind: {k: [H.restrict(w, D[k].wf, H.core(D[k].wf)) for w in draws[kind][k]] for k in D} for kind in draws}
    Tc = full(D, wfc, n_perm, n_sets, n_boot, compc)
    # trio gate (the designs' reference gate on the trio)
    tg = {}
    for k in D:
        try:
            tg[k] = gate_for(D[k], present(D[k], trio))
        except Exception:      # a trio gene absent: the plane of two genes is degenerate
            tg[k] = np.nan
    n_tg = sum(1 for k in D if tg[k] < 0.05 and D[k].n >= 20)
    # unit bootstrap (pairs within dataset), shared by 14 genes and core
    rng = np.random.default_rng([336, 1])
    d_off, d_r, off_c, pr_c = [], [], [], []
    for _ in range(n_boot):
        pick = {k: rng.integers(0, D[k].n, D[k].n) for k in D}

        def one(T):
            if not T["use"]:
                return np.nan, np.nan
            cs_ = [P.ccorr(*[x[pick[k]] for x in T["B"][k]]) for k in T["use"]]
            rb = np.array([P.log_var_ratio(*[x[pick[k]] for x in T["B"][k]]) for k in T["use"]])
            return P.pooled_offset(cs_, T["w"]), float(np.nansum(T["iv"] / T["iv"].sum() * rb))
        o14, r14 = one(T14)
        oc, rc = one(Tc)
        d_off.append(float((o14 - oc + 12) % 24 - 12))
        d_r.append(r14 - rc)
        off_c.append(oc)
        pr_c.append(rc)
    # gene-drop null
    dl, tl = H.drops(WF)
    wfs = dict(dl)
    v_off, v_r = {}, {}
    for lab, w in (dl[:3] if chk else dl):
        try:
            o, r, nuse = drop_stats(D, w)
            v_off[lab] = -abs(o) if nuse >= 3 else np.nan
            v_r[lab] = r
        except Exception:
            v_off[lab] = v_r[lab] = np.nan
    if chk:
        assert tl in wfs and H.drops(WF, n=3)[1] == tl
        print("check ok")
        return
    assert np.isclose(v_r[tl], Tc["pooled_r"]) and np.isclose(v_off[tl], -abs(Tc["pooled_offset_h"])), (v_r[tl], Tc["pooled_r"], v_off[tl], Tc["pooled_offset_h"])
    p_off = H.drop_percentiles(v_off, tl, wfs, +1)
    p_r = H.drop_percentiles(v_r, tl, wfs, -1)
    # H254.1 (equivalence: the registered margin is the 1.0 h detectability limit)
    oc = np.array(off_c)
    dev = (oc - Tc["pooled_offset_h"] + 12) % 24 - 12 + Tc["pooled_offset_h"]
    lo, hi = float(np.nanquantile(dev, .05)), float(np.nanquantile(dev, .95))
    scored = len(Tc["use"]) >= 3 and Tc["mde_h"] <= MARGIN_H
    inside = bool(lo > -MARGIN_H and hi < MARGIN_H)
    spec_delay = bool(Tc["pooled_offset_h"] < Tc["matched_q05_h"] and Tc["pooled_offset_h"] < Tc["responsive_q05_h"] and Tc["p_delay"] < 0.05)
    st1 = "holds" if (scored and inside and not spec_delay) else "falls"
    drd = np.array(d_off)
    n1 = (f"core: pooled offset {Tc['pooled_offset_h']:.4f} h (90% bootstrap interval {lo:.3f} to {hi:.3f}, 200 pair bootstraps), null SD {Tc['null_sd_h']:.4f}, p_delay {Tc['p_delay']:.4g}, "
          f"MDE {Tc['mde_h']:.3f} h, matched q05 {Tc['matched_q05_h']:.4f}, responsive q05 {Tc['responsive_q05_h']:.4f}; datasets gated {Tc['use']}; conditions [scored (>=3 datasets, MDE<=1.0 h), "
          f"interval inside +-{MARGIN_H} h, no clock-specific delay] = [{scored}, {inside}, {not spec_delay}]; drop null on -|pooled offset| (larger = closer to zero), undefined when < 3 gated datasets; "
          f"trio gate: reference gate p per dataset {({k: round(v, 3) for k, v in tg.items()})} (3 genes allow 6 acrophase arrangements: p cannot fall below about 1/6)")
    # H254.3
    c3 = [Tc["p_r"] < 0.05, all((Tc["pooled_r"] < Tc[f"{kd}_r_q05"]) or (Tc["pooled_r"] > Tc[f"{kd}_r_q95"]) for kd in ("matched", "responsive"))]
    st3 = "holds" if all(c3) else "falls" if not any(c3) else "at the line"
    dr3 = np.array(d_r)
    n3 = (f"core: pooled r {Tc['pooled_r']:.4f} (SE {Tc['pooled_r_se']:.4f}), swap p {Tc['p_r']:.4g}, matched q05/q95 {Tc['matched_r_q05']:.4f}/{Tc['matched_r_q95']:.4f}, responsive q05/q95 "
          f"{Tc['responsive_r_q05']:.4f}/{Tc['responsive_r_q95']:.4f}; datasets gated {Tc['use']}; conditions [pooled swap p<0.05, outside both companions' 90% range] = {c3}; "
          "drops: pooled r, direction -1 (more negative = stronger); trio gate as H254.1")
    gate_note = bool(n_tg >= 3) if np.isfinite(list(tg.values())).any() else "not computable"
    rows = [dict(id="A254_H254.1", figure_panel="Fig 6 (chronic lesions: time kept)", statistic="pooled lesion-minus-non-lesion clock offset, 5 gated datasets (h)", direction="-1 (delay); drop null on -|offset|",
                 value_14_stored=ST.pooled_offset_h, value_14=T14["pooled_offset_h"], reproduced=True, value_core=Tc["pooled_offset_h"],
                 core_criterion="equivalence: scored (>=3 gated datasets, pooled MDE <= 1.0 h), core 90% bootstrap interval inside +-1.0 h (registered detectability limit) and no clock-specific delay",
                 core_status=st1, diff_14_minus_core=T14["pooled_offset_h"] - Tc["pooled_offset_h"], diff_lo90=float(np.nanquantile(drd, .05)), diff_hi90=float(np.nanquantile(drd, .95)),
                 diff_draws="own: 200 pair bootstraps (within dataset, shared by both waveforms), seed [336, 1]", **p_off, trio_gate_passes=gate_note,
                 trio_alone_value=np.nan, reading=H.reading(st1, p_off["trio_drop_below_q10"]), blind="core blind", notes=n1, script=SCRIPT),
            dict(id="A254_H254.3", figure_panel="Fig 6 (chronic lesions: loudness)", statistic="pooled log clock-plane variance ratio r (lesion / non-lesion)", direction="-1",
                 value_14_stored=ST.pooled_r, value_14=T14["pooled_r"], reproduced=True, value_core=Tc["pooled_r"],
                 core_criterion="pooled swap p < 0.05 and r outside the central 90% of both companion families (clock-specific)", core_status=st3,
                 diff_14_minus_core=T14["pooled_r"] - Tc["pooled_r"], diff_lo90=float(np.nanquantile(dr3, .05)), diff_hi90=float(np.nanquantile(dr3, .95)),
                 diff_draws="own: 200 pair bootstraps (within dataset, shared by both waveforms), seed [336, 1]", **p_r, trio_gate_passes=gate_note,
                 trio_alone_value=np.nan, reading=H.reading(st3, p_r["trio_drop_below_q10"]), blind="core blind", notes=n3, script=SCRIPT)]
    H.write_rows(ENTRY, rows)
    log("done", f"{(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    main()
