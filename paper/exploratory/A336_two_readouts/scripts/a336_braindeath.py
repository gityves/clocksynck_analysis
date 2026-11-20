"""A336 part 4 (registered): the brain-death order (A317 H317.1, A318 H318.1) on the core loop, with the gene-drop null.

D_rho = rho(post-mortem, fast) - median over 20 organ-donor draws (to 79, seed 317) of rho(organ donors, fast), Fisher-Lee
of pair phases (A317). Waveforms: 14 genes (reproduced first), core loop (`cs.gtex().drop(PER1, PER2, NR1D1)`, A317's
matched and responsive sets restricted to the core positions: A329 / A335 rule), every 11-gene drop of three genes
(364; clock only), and the trio alone (descriptive, only if its gate passes). A318: the same on states residualised on
ischaemic time and RIN. Paired difference 14 genes minus core on 200 shared donor bootstraps (seed [336, 4]).
Aggregates only (GTEx controlled access).

    python a336_braindeath.py --stage power    # core relabelling null only (no observed core D_rho)
    python a336_braindeath.py --stage score
    python a336_braindeath.py --stage check    # 3 drops, 3 sets, 20 permutations; nothing saved or printed
Outputs (`../output/`): braindeath_core_power.csv | braindeath_core.csv, braindeath_drops.csv.
"""
import argparse
import copy
import itertools
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve()
EXP = HERE.parents[2]
for d in ("A109_coupling_nulls", "A132_brain_death", "A317_braindeath_confirm", "A335_review_batch"):
    sys.path.insert(0, str(EXP / d / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import a335_gtex as G  # noqa: E402
import score_a317 as S  # noqa: E402
from score_a317 import B, H, gtex  # noqa: E402
import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402

OUT = HERE.parents[1] / "output"
WF = H.WF
TRIO = ["PER1", "PER2", "NR1D1"]
DROPS = list(itertools.combinations(WF.genes, 3))
TRIO_DROP = DROPS.index(tuple(g for g in WF.genes if g in TRIO))
N_BOOT = 200


def sv(w):
    r = np.asarray(w.amp) / np.asarray(w.sd)
    return np.linalg.svd(np.c_[r * np.cos(w.psi), r * np.sin(w.psi)], compute_uv=False)


def drop_states(T, waveforms, kind):
    """Clock states for a list of waveforms on the real clock-gene rows: T.Z[kind][t] (K, rows)."""
    T.Z[kind], T.V[kind] = {}, {}
    for t in H.BODY:
        T.Z[kind][t], T.V[kind][t] = H._state_set(T.rows[t], waveforms, T.rows[t].depth.to_numpy(float))


def adjust_all(T, h0, attrs):
    """A318 `adjust` applied to every kind of states (clock, sets, drops, trio); returns R without unattributed rows."""
    R = {}
    for t in H.BODY:
        r = T.R[t].copy()
        a = attrs[t]
        isch, rin = a.ischemic_time.to_numpy(float), a.rin.to_numpy(float)
        rows = r[h0][r[h0] >= 0]
        ok = rows[np.isfinite(isch[rows]) & np.isfinite(rin[rows])]
        r[np.isin(r, np.setdiff1d(rows, ok)) & (r >= 0)] = -1
        R[t] = r
        X = np.column_stack([np.ones(len(ok)), (isch[ok] - isch[ok].mean()) / isch[ok].std(),
                             (rin[ok] - rin[ok].mean()) / rin[ok].std()])
        Pm = X @ np.linalg.pinv(X)
        for kind in T.Z:
            if t not in T.Z[kind] or np.ndim(T.Z[kind][t]) < 2 or T.Z[kind][t].shape[0] == 0:
                continue
            sub = T.Z[kind][t][:, ok]
            T.Z[kind][t][:, ok] = sub - (Pm @ sub.T).T + sub.mean(1, keepdims=True)
    return R


def trio_gate(T, fast):
    """clocksync model-free gate of the trio waveform on the fast deaths' 12 body tissues."""
    rows = []
    for t in H.BODY:
        r = T.R[t][fast]
        x = T.rows[t].iloc[r[r >= 0]].copy()
        x.insert(0, "compartment", t)
        x.insert(0, "donor", T.donors[fast][r >= 0])
        rows.append(x)
    pb = pd.concat(rows, ignore_index=True)
    ds = cs.Dataset(pb, anchors=None, waveform=WF.subset(TRIO), name="trio")
    g = ds.gate(n_perm=1000, seed=336)
    return dict(trio_gate_p_pairing=g.p_pairing, trio_gate_p_acrophase=g.p_acrophase,
                trio_gate_passes=bool(g.p_pairing < 0.05 and g.p_acrophase < 0.05))


def design(T, R, fast, od, pm, seed):
    """Organ-donor draws (to the post-mortem n), pair mask from the 14-gene clock presence, relabelling sequence."""
    rng = np.random.default_rng(seed)
    od_draws = [rng.choice(od, len(pm), replace=False) for _ in range(S.N_DRAWS)]
    Zc, Mc = S.raw(T, "clock", 0, R)
    mask = B.pair_mask(Mc, fast, pm, *od_draws[:1]) & B.pair_mask(Mc, od)
    return rng, od_draws, mask


def perm_null(T, kind, k, R, fast, h0, npm, mask, rng, n):
    Z, M = S.raw(T, kind, k, R)
    out = []
    for _ in range(n):
        p = rng.permutation(h0)
        out.append(S.rho(Z, M, fast, p[:npm], mask) - S.rho(Z, M, fast, p[npm:2 * npm], mask))
    return np.array(out)


def readout(T, R, fast, od, pm, h0, seed, n_perm, n_sets, kinds):
    """D_rho, relabelling p and companion q95 for each (kind, set-kind) of `kinds`: {label: (clock kind, matched kind,
    responsive kind)}; relabellings replayed from the same seed for every label (same draws)."""
    out = {}
    for lab, (ck, mk, rk) in kinds.items():
        rng, od_draws, mask = design(T, R, fast, od, pm, seed)
        Z, M = S.raw(T, ck, 0, R)
        D = S.d_rho(Z, M, fast, pm, od_draws, mask)
        rn = perm_null(T, ck, 0, R, fast, h0, len(pm), mask, rng, n_perm)
        res = dict(D_rho=D, p_labels=float((1 + np.sum(rn >= D)) / (1 + n_perm)), null_sd=float(np.nanstd(rn, ddof=1)),
                   pairs=int(mask.sum()))
        for fam, kk in (("matched", mk), ("responsive", rk)):
            v = np.array([S.d_rho(*S.raw(T, kk, j, R), fast, pm, od_draws, mask) for j in range(n_sets)])
            res[f"{fam}_q95"] = float(np.nanquantile(v, 0.95))
            res[f"pct_{fam}"] = float(np.nanmean(v < D))
        res["passes"] = bool(D > 0 and res["p_labels"] < 0.05 and D > res["matched_q95"] and D > res["responsive_q95"])
        out[lab] = res
    return out


def drop_null(T, R, fast, od, pm, seed, kind="drop"):
    _, od_draws, mask = design(T, R, fast, od, pm, seed)
    K = T.Z[kind][H.BODY[0]].shape[0]
    return np.array([S.d_rho(*S.raw(T, kind, k, R), fast, pm, od_draws, mask) for k in range(K)])


def paired_boot(T, R, fast, od, pm, seed, n):
    """D_rho(14) - D_rho(core) on shared donor bootstraps (groups resampled, organ-donor draws from the resample)."""
    _, _, mask = design(T, R, fast, od, pm, seed)
    rng = np.random.default_rng([336, 4])
    Z14, M14 = S.raw(T, "clock", 0, R)
    Zc, Mc = S.raw(T, "core", 0, R)
    out = []
    for _ in range(n):
        f, p, o = rng.choice(fast, len(fast)), rng.choice(pm, len(pm)), rng.choice(od, len(od))
        dr = [rng.choice(o, len(p), replace=False) for _ in range(S.N_DRAWS)]
        out.append((S.d_rho(Z14, M14, f, p, dr, mask), S.d_rho(Zc, Mc, f, p, dr, mask)))
    return np.array(out)


def pct_rows(label, null, D_trio, svs, sv_trio):
    """Percentile of the trio drop (the core loop) among all drops and among drops of similar singular values."""
    sim = np.all(np.abs(svs / sv_trio - 1) <= 0.2, axis=1)
    return dict(readout=label, n_drops=int(np.isfinite(null).sum()), pct_all=float(np.nanmean(null < D_trio) * 100),
                n_similar=int(sim.sum()), pct_similar=float(np.nanmean(null[sim] < D_trio) * 100) if sim.sum() >= 20 else np.nan,
                drop_q10=float(np.nanquantile(null, 0.10)), drop_median=float(np.nanmedian(null)),
                below_q10=bool(np.nanmean(null < D_trio) < 0.10))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=["power", "score", "check"], required=True)
    stage = ap.parse_args().stage
    chk = stage == "check"
    limit_threads(4)
    t0 = time.time()
    log = (lambda *a, **k: None) if chk else (lambda *a, **k: print(*a, **k, flush=True))
    T = H.Tissues(200, 0)
    T.build_states(np.random.default_rng(1), 0, real=True)
    Zc, Vc = G.build_core(T, None)
    G.add_resp(T, Zc, Vc)
    T.Z.update(core=Zc["clock"], core_pseudo=Zc["pseudo"], core_resp=Zc["resp"])
    T.V.update(core=Vc["clock"], core_pseudo=Vc["pseudo"], core_resp=Vc["resp"])
    drops = [WF.drop(list(c)) for c in DROPS]
    if chk:
        drops = drops[:2] + [drops[TRIO_DROP]]
    drop_states(T, drops, "drop")
    drop_states(T, [WF.subset(TRIO)], "trio")
    assert np.allclose(T.Z["drop"][H.BODY[0]][-1 if chk else TRIO_DROP], T.Z["core"][H.BODY[0]][0])
    ph = gtex.phenotypes(B.PH_COLS).reindex(T.donors)
    coh = ph.COHORT.to_numpy()
    fast, h0 = T.group_donors("fast"), T.group_donors("vent")
    od, pm = h0[coh[h0] == B.OPO], h0[coh[h0] == B.PM]
    print(f"states built {time.time() - t0:.0f}s: od {len(od)}, pm {len(pm)}, fast {len(fast)}, drops {len(drops)}", flush=True)
    OUT.mkdir(parents=True, exist_ok=True)
    n_perm, n_sets = (20, 3) if chk else (S.N_PERM, 200)

    if stage == "power":
        rng, _, mask = design(T, T.R, fast, od, pm, S.SEED)
        rn = perm_null(T, "core", 0, T.R, fast, h0, len(pm), mask, rng, n_perm)
        res = dict(pairs=int(mask.sum()), sd_D_rho_core=float(np.nanstd(rn, ddof=1)),
                   mde_D_rho_core=float(2.49 * np.nanstd(rn, ddof=1)),
                   a317_mde_D_rho=float(pd.read_csv(EXP / "A317_braindeath_confirm/output/power.csv", index_col=0).loc["mde_D_rho", "value"]))
        pd.Series(res).to_csv(OUT / "braindeath_core_power.csv", header=["value"])
        print(pd.Series(res).to_string(), flush=True)
        return

    # ---- 14 genes reproduced (A317 H317.1), before any core value -------------------------------------------------
    kinds14 = {"14 genes": ("clock", "pseudo", "resp")}
    r14 = readout(T, T.R, fast, od, pm, h0, S.SEED, n_perm, n_sets, kinds14)["14 genes"]
    if not chk:
        ref = pd.read_csv(EXP / "A317_braindeath_confirm/output/tests.csv").iloc[0]
        rep = {k: (r14[k], ref[k]) for k in ("D_rho", "p_labels", "matched_q95", "responsive_q95")}
        bad = {k: v for k, v in rep.items() if not np.isclose(*v, atol=1e-9)}
        log("A317 reproduction:", {k: (round(a, 6), round(b, 6)) for k, (a, b) in rep.items()})
        if bad:
            raise SystemExit(f"A317 not reproduced: stop before the core loop: {bad}")
    print(f"A317 reproduced {time.time() - t0:.0f}s", flush=True)

    # ---- core, drops, trio, paired difference: unadjusted ---------------------------------------------------------
    rc = readout(T, T.R, fast, od, pm, h0, S.SEED, n_perm, n_sets, {"core": ("core", "core_pseudo", "core_resp")})["core"]
    nd = drop_null(T, T.R, fast, od, pm, S.SEED)
    svs = np.array([sv(w) for w in drops])
    trio_idx = len(drops) - 1 if chk else TRIO_DROP
    assert np.isclose(nd[trio_idx], rc["D_rho"])
    _, od_draws, mask = design(T, T.R, fast, od, pm, S.SEED)
    D_trio_alone = S.d_rho(*S.raw(T, "trio", 0, T.R), fast, pm, od_draws, mask)
    gate = trio_gate(T, fast) if not chk else dict(trio_gate_passes=False)
    pb = paired_boot(T, T.R, fast, od, pm, S.SEED, 5 if chk else N_BOOT)
    diff = pb[:, 0] - pb[:, 1]
    print(f"unadjusted done {time.time() - t0:.0f}s", flush=True)

    # ---- A318 adjustment: reproduce H318.1 on 14 genes, then core and drops -----------------------------------------
    Ta = copy.copy(T)
    Ta.Z = {k: {t: v.copy() for t, v in d.items()} for k, d in T.Z.items()}
    attrs = S.sample_attrs(T)
    Ra = adjust_all(Ta, h0, attrs)
    od_a = od[np.array([any(Ra[t][d] >= 0 for t in H.BODY) for d in od])]
    pm_a = pm[np.array([any(Ra[t][d] >= 0 for t in H.BODY) for d in pm])]
    h0_a = np.concatenate([od_a, pm_a])
    ra = readout(Ta, Ra, fast, od_a, pm_a, h0_a, S.SEED, n_perm, n_sets,
                 {"14 genes": ("clock", "pseudo", "resp"), "core": ("core", "core_pseudo", "core_resp")})
    if not chk:
        ref = pd.read_csv(EXP / "A317_braindeath_confirm/output/a318_tests.csv").iloc[0]
        log("A318 D_rho reproduction:", round(ra["14 genes"]["D_rho"], 6), round(ref.D_rho, 6))
        if not np.isclose(ra["14 genes"]["D_rho"], ref.D_rho, atol=1e-9):
            raise SystemExit("A318 D_rho not reproduced: adjusted core values not written")
    nda = drop_null(Ta, Ra, fast, od_a, pm_a, S.SEED)
    print(f"adjusted done {time.time() - t0:.0f}s", flush=True)
    if chk:
        print("check ok")
        return

    rows = []
    for adj, res, null in (("unadjusted (H317.1)", {"14 genes": r14, "core": rc}, nd),
                           ("ischaemic time and RIN adjusted (H318.1)", ra, nda)):
        p = pct_rows(adj, null, res["core"]["D_rho"], svs, sv(WF.drop(TRIO)))
        for gs in ("14 genes", "core"):
            rows.append(dict(adjustment=adj, gene_set=gs, **res[gs], **(p if gs == "core" else {})))
    tab = pd.DataFrame(rows)
    tab.loc[0, ["diff_14_minus_core", "diff_lo90", "diff_hi90"]] = [r14["D_rho"] - rc["D_rho"], *np.quantile(diff, [0.05, 0.95])]
    tab.loc[0, ["trio_alone_D_rho"]] = D_trio_alone if gate["trio_gate_passes"] else np.nan
    for k, v in gate.items():
        tab.loc[0, k] = v
    core = tab[(tab.gene_set == "core")].iloc[0]
    tab["reading"] = ""
    tab.loc[tab.gene_set == "core", "reading"] = [
        "holds" if r.passes else "falls: carried by PER1, PER2 and NR1D1" if r.below_q10 else "falls: not decidable (precision)"
        for _, r in tab[tab.gene_set == "core"].iterrows()]
    tab.to_csv(OUT / "braindeath_core.csv", index=False)
    pd.DataFrame(dict(drop=["+".join(c) for c in DROPS], sv1=svs[:, 0], sv2=svs[:, 1], D_rho=nd, D_rho_adjusted=nda)).to_csv(
        OUT / "braindeath_drops.csv", index=False)
    pd.set_option("display.width", 220)
    print(tab.round(4).T.to_string())
    print(f"reading (registered): {core.reading}; done {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
