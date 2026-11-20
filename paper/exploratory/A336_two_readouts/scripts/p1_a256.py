"""A336 part 1, Figure 1b: A256 H256.2 (rho, all fresh units) and H256.3 (cells within 2 h of the protocol) on the core loop.

Original: `A256_elapsed_replication/scripts/compute_a256.py` (imported, not edited). Same data, seeds ([256, 1] label
permutations, [256, tag, sum(ord(acc))] companions, 200 sets), same statistics; only the waveform changes. Companions: the
original 14-gene draws restricted to the core positions of each dataset (A329 rule). H256.1 enters only through the Holm
p of H256.2 (as in the original).

    python p1_a256.py [--check]     # --check: tiny run, assertions only, nothing printed about values nor written
Output: ../output/part1/a256.csv.
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
sys.path.insert(0, str(PAPER / "exploratory/A252_elapsed_time/scripts"))
sys.path.insert(0, str(PAPER / "exploratory/A256_elapsed_replication/scripts"))
sys.path.insert(0, str(HERE.parent))
import compute_a256 as A  # noqa: E402
import a336_common as H  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402

ENTRY = "a256"
ST = pd.read_csv(PAPER / "exploratory/A256_elapsed_replication/output/score_tests.csv").set_index("family")
SD = pd.read_csv(PAPER / "exploratory/A256_elapsed_replication/output/score_descriptive.csv")
S23 = ST.loc["H256.2 all fresh units"]
S33 = SD.iloc[0]
N_SETS, N_PERM, N_BOOT = 200, 10_000, 200
SCRIPT = "paper/exploratory/A336_two_readouts/scripts/p1_a256.py"


def log(*a):
    print(time.strftime("%F %T"), *a, flush=True)


def holm2(p1, p2):
    """The original's Holm over H256.1-2 (its formula verbatim); returns the p_holm of H256.2."""
    p = np.array([p1, p2])
    return float(np.minimum(1, [max(2 * p.min(), p.max()) if q == p.max() and len(set(p)) > 1 else 2 * p.min() for q in p])[1])


def make(chk):
    D = {a: A.Data(a) for a in A.DSETS}
    nd = {a: D[a].cls == "nothing_done" for a in A.DSETS}
    return D, nd


def parts_all(D, Z, el=None):
    return A.rho([D[a].parts(Z[a], None if el is None else el[a]) for a in Z])


def count_all(D, Z, el=None):
    return sum(A.within2(D[a].cells(Z[a], None if el is None else el[a])) for a in Z)


def null_pass(D, nd, Zs, n_perm):
    """label permutations with the original stream: per waveform key, H256.1 / H256.2 rho nulls and the count null."""
    rng = np.random.default_rng([A.SEED, 1])
    out = {k: dict(r1=[], r2=[], n=[]) for k in Zs}
    for _ in range(n_perm):
        els = {a: D[a].perm_el(rng) for a in A.DSETS}
        for k, Z in Zs.items():
            out[k]["r1"].append(A.rho([D[a].parts(Z[a], els[a], nd[a]) if nd[a].any() else (0.0, 0.0) for a in Z]))
            out[k]["r2"].append(parts_all(D, Z, els))
            out[k]["n"].append(count_all(D, Z, els))
    return {k: {m: np.array(v) for m, v in d.items()} for k, d in out.items()}


def tests(D, Z, null, comp, n_sets):
    obs = parts_all(D, Z)
    r = dict(rho=obs, p=A.p_up(null["r2"], obs), count=count_all(D, Z))
    r["p_count"] = A.p_up(null["n"], r["count"])
    r["null_mean_count"] = float(null["n"].mean())
    for kind in comp:
        s = np.array([A.rho([D[a].parts(comp[kind][a][j]) for a in A.DSETS]) for j in range(n_sets)])
        r[f"{kind}_q95"] = float(np.quantile(s, .95))
    p1 = A.p_up(null["r1"], A.rho([D[a].parts(Z[a], None, nd_[a]) if nd_[a].any() else (0.0, 0.0) for a in A.DSETS]))
    r["p_holm"] = holm2(p1, r["p"])
    return r


def unit_contrib(d, z):
    zb = np.empty(len(z), complex)
    for i in d.idx:
        zb[i] = z[i][d.el[i] == 0][0]
    c = np.exp(1j * A.W * (d.t0 + d.el)) - np.exp(1j * A.W * d.t0)
    m = d.el > 0
    num = np.bincount(d.unit, weights=np.where(m, ((z - zb) * np.conj(c)).real, 0.0), minlength=len(d.idx))
    den = np.bincount(d.unit, weights=np.where(m, abs(c) ** 2, 0.0), minlength=len(d.idx))
    return num, den, z - zb


def count_cells(d, dz, sel):
    """number of cells within 2 h of the protocol among the samples `sel` (the original's `cells` on a resample)."""
    arm, el, t0 = d.sh.arm.to_numpy()[sel], d.el[sel], d.t0[sel]
    dd, n = dz[sel], 0
    for a in pd.unique(d.sh.arm):
        for e in sorted(set(d.el[d.el > 0])):
            m = (arm == a) & (el == e)
            if m.sum() < 3:
                continue
            c = np.exp(1j * A.W * (t0[m][0] + e)) - np.exp(1j * A.W * t0[m][0])  # noqa: F841
            imp = ((np.angle(dd[m].mean()) - np.pi / 2) / A.W) % 24
            n += abs((imp - (t0[m][0] + e / 2) % 24 + 12) % 24 - 12) <= 2
    return int(n)


nd_ = {}


def main():
    global nd_
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    chk = ap.parse_args().check
    limit_threads(4)
    t_start = time.time()
    n_sets, n_perm, n_boot = (2, 20, 5) if chk else (N_SETS, N_PERM, N_BOOT)
    D, nd_ = make(chk)
    WF = H.WF
    # ---- stage A: 14 genes
    Z14 = {a: D[a].z() for a in A.DSETS}
    draws = {kind: {a: D[a].pseudo(kind, n_sets, tag) for a in A.DSETS} for kind, tag in (("matched", 3), ("responsive", 4))}
    comp14 = {kind: {a: [D[a].z(w) for w in draws[kind][a]] for a in A.DSETS} for kind in draws}
    N14 = null_pass(D, nd_, {"14": Z14}, n_perm)["14"]
    T14 = tests(D, Z14, N14, comp14, n_sets)
    log("stage A done")
    exp2 = dict(rho=S23.rho, p=S23.p, p_holm=S23.p_holm, matched_q95=S23.matched_q95, responsive_q95=S23.responsive_q95)
    ok2 = all(np.isclose(T14[k], v, rtol=0, atol=1e-6) for k, v in exp2.items())
    ok3 = bool(np.isclose(T14["count"], S33.value, atol=1e-6) and np.isclose(T14["p_count"], S33.p, atol=1e-6)
               and np.isclose(T14["null_mean_count"], S33.null_mean, atol=1e-6))
    if chk:
        for a in A.DSETS:
            nu, de, dz = unit_contrib(D[a], Z14[a])
            assert np.isclose(nu.sum(), D[a].parts(Z14[a])[0]) and np.isclose(de.sum(), D[a].parts(Z14[a])[1])
            assert count_cells(D[a], dz, np.arange(len(dz))) == A.within2(D[a].cells(Z14[a]))
    else:
        log("H256.2", "reproduced" if ok2 else "NOT reproduced", "H256.3", "reproduced" if ok3 else "NOT reproduced")
    if not chk and not (ok2 and ok3):
        why = "; ".join(f"{k}: recomputed {T14[k]} vs stored {v}" for k, v in exp2.items()) if not ok2 else ""
        why3 = f"count {T14['count']} (p {T14['p_count']}, null mean {T14['null_mean_count']}) vs stored {S33.value} (p {S33.p}, null mean {S33.null_mean})" if not ok3 else ""
        rows = [dict(id="A256_H256.2", figure_panel="Fig 1b", statistic="rho, all fresh units pooled", direction="+1", value_14_stored=S23.rho,
                     value_14=T14["rho"], reproduced=bool(ok2), notes="" if ok2 else "14 genes did not reproduce; no core value computed. " + why, script=SCRIPT),
                dict(id="A256_H256.3", figure_panel="Fig 1b", statistic="cells within 2 h of protocol midpoint", direction="+1", value_14_stored=S33.value,
                     value_14=T14["count"], reproduced=bool(ok3), notes="" if ok3 else "14 genes did not reproduce; no core value computed. " + why3, script=SCRIPT)]
        H.write_rows(ENTRY, rows)
        return
    # ---- stage B
    core, trio = H.core(WF), H.trio(WF)
    Zc = {a: D[a].z(core) for a in A.DSETS}
    Zt, skipped = {}, []
    for a in A.DSETS:   # the trio is estimable only with >= 3 of its genes present on the platform (GSE9105 lacks PER1)
        try:
            Zt[a] = D[a].z(trio)
        except Exception:
            skipped.append(a)
    compc = {kind: {a: [D[a].z(H.restrict(w, D[a].wf, H.core(D[a].wf))) for w in draws[kind][a]] for a in A.DSETS} for kind in draws}
    NB = null_pass(D, nd_, {"core": Zc, "trio": Zt}, n_perm)
    Tc = tests(D, Zc, NB["core"], compc, n_sets)
    # trio alone
    t_rho, t_cnt = parts_all(D, Zt), count_all(D, Zt)
    p_tr, p_tc = A.p_up(NB["trio"]["r2"], t_rho), A.p_up(NB["trio"]["n"], t_cnt)
    # bootstrap
    U14 = {a: unit_contrib(D[a], Z14[a]) for a in A.DSETS}
    Uc = {a: unit_contrib(D[a], Zc[a]) for a in A.DSETS}
    rng = np.random.default_rng([336, 1])
    d_rho, d_cnt = [], []
    for _ in range(n_boot):
        pick = {a: rng.integers(0, len(D[a].idx), len(D[a].idx)) for a in A.DSETS}
        sel = {a: np.concatenate([D[a].idx[u] for u in pick[a]]) for a in A.DSETS}
        f = lambda U: sum(U[a][0][pick[a]].sum() for a in A.DSETS) / sum(U[a][1][pick[a]].sum() for a in A.DSETS)  # noqa: E731
        d_rho.append(f(U14) - f(Uc))
        d_cnt.append(sum(count_cells(D[a], U14[a][2], sel[a]) - count_cells(D[a], Uc[a][2], sel[a]) for a in A.DSETS))
    # gene-drop null
    dl, tl = H.drops(WF)
    wfs = dict(dl)
    v_rho, v_cnt = {}, {}
    for lab, w in (dl[:3] if chk else dl):
        try:
            Z = {a: D[a].z(w) for a in A.DSETS}
            v_rho[lab], v_cnt[lab] = parts_all(D, Z), float(count_all(D, Z))
        except Exception:
            v_rho[lab] = v_cnt[lab] = np.nan
    if chk:
        assert tl in wfs and H.drops(WF, n=3)[1] == tl
        print("check ok")
        return
    assert np.isclose(v_rho[tl], Tc["rho"]) and v_cnt[tl] == Tc["count"]
    pr, pcn = H.drop_percentiles(v_rho, tl, wfs, +1), H.drop_percentiles(v_cnt, tl, wfs, +1)
    conds = [Tc["rho"] > 0, Tc["p_holm"] < 0.05, Tc["rho"] > Tc["matched_q95"] and Tc["rho"] > Tc["responsive_q95"]]
    st2 = "holds" if all(conds) else "falls" if not any(conds) else "at the line"
    st3 = "holds" if Tc["p_count"] < 0.05 else "falls"
    tg2, tg3 = p_tr < 0.05, p_tc < 0.05
    dr, dc = np.array(d_rho), np.array(d_cnt)
    rows = [dict(id="A256_H256.2", figure_panel="Fig 1b", statistic="rho, all fresh units pooled (7 studies, 166 units)", direction="+1",
                 value_14_stored=S23.rho, value_14=T14["rho"], reproduced=True, value_core=Tc["rho"],
                 core_criterion="rho > 0, p_holm (over H256.1-2, as stored) < 0.05 and rho above matched and responsive q95", core_status=st2,
                 diff_14_minus_core=T14["rho"] - Tc["rho"], diff_lo90=float(np.quantile(dr, .05)), diff_hi90=float(np.quantile(dr, .95)),
                 diff_draws="own: 200 unit bootstraps (units within dataset), seed [336, 1]", **pr, trio_gate_passes=bool(tg2),
                 trio_alone_value=t_rho if tg2 else np.nan, reading=H.reading(st2, pr["trio_drop_below_q10"]), blind="core blind",
                 notes=(f"core: rho {Tc['rho']:.4f}, p {Tc['p']:.4g}, p_holm {Tc['p_holm']:.4g}, matched q95 {Tc['matched_q95']:.4f}, responsive q95 "
                        f"{Tc['responsive_q95']:.4f}; conditions [rho>0, p_holm<0.05, clock-specific] = {conds}; trio alone rho {t_rho:.4f}, permutation p {p_tr:.4g} (datasets without an estimable trio state excluded: {skipped})"),
                 script=SCRIPT),
            dict(id="A256_H256.3", figure_panel="Fig 1b", statistic="cells (dataset x arm x time point) with implied midpoint within 2 h of the protocol's, count",
                 direction="+1", value_14_stored=S33.value, value_14=T14["count"], reproduced=True, value_core=Tc["count"],
                 core_criterion="count within 2 h larger than under permuted time labels (permutation p < 0.05, as registered)", core_status=st3,
                 diff_14_minus_core=T14["count"] - Tc["count"], diff_lo90=float(np.quantile(dc, .05)), diff_hi90=float(np.quantile(dc, .95)),
                 diff_draws="own: 200 unit bootstraps (units within dataset; cells with >= 3 samples recounted), seed [336, 1]", **pcn,
                 trio_gate_passes=bool(tg3), trio_alone_value=t_cnt if tg3 else np.nan, reading=H.reading(st3, pcn["trio_drop_below_q10"]),
                 blind="core blind",
                 notes=(f"core: count {Tc['count']} of 20 cells, permuted-label mean {Tc['null_mean_count']:.4f}, p {Tc['p_count']:.4g}; "
                        f"trio alone count {t_cnt} (of the cells of the other datasets; excluded: {skipped}), permutation p {p_tc:.4g}; magnitude in the drop null = the count"),
                 script=SCRIPT)]
    H.write_rows(ENTRY, rows)
    log("done", f"{(time.time() - t_start) / 60:.1f} min")


if __name__ == "__main__":
    main()
