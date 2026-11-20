"""A336 part 1, A304 (tumour gain and departure on the own normal, TCGA paired): H304.1 gain, H304.1 departure, H304.3, H304.4.

Original: `A304_cancer_gain_departure/scripts/score_a304.py` + `lib_a304.py` (imported, not edited). Same data, seeds (label
permutations [304, contrast index], 2000 draws; companions [304, 7 / 8, cancer index], 200 matched + 200 responsive per
cancer), same statistics (Q_G, Q_R pooled over cancers, Holm families as original); only the waveform changes.
Companions: the original 14-gene pseudo-waveforms restricted to the core positions (A329 / A335 rule).

Nulls (departure rise "not detected"): core "holds" if still not detected (Holm p >= 0.05 or not above matched q95) AND the
upper end of the 90 % pair-bootstrap interval of Q_R (200 draws, rng [336, 1], within-cancer, observed labels) is below the
registered margin = mean Q_R of A304's power simulation at its registered MDE (power/mde.csv mde_k_D; power.csv
mean_Q_R; 14-gene simulation, power stage not rerun). Not detected but interval crossing the margin, or detected: "falls".
Bootstrap keeps V_N^all of each cell fixed. Directions: gain -1 (loss), departure +1 (rise).
Trio-alone gate (paired design): H138.4-style pairing of trio tumour/normal states (`q3_lib.corr_stat`, `perm_null`,
`pool`; 2000 within-cancer permutations, rng [336, 2]) p < 0.05.

    python p1_a304.py [--check] [--jobs 4]
Output: ../output/part1/a304.csv
"""
import argparse
import importlib.util
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
sys.path.insert(0, str(PAPER))
SC = PAPER / "exploratory/A304_cancer_gain_departure/scripts"
sys.path.insert(0, str(SC))
sys.path.insert(0, str(HERE.parent))
import score_a304 as A  # noqa: E402
import lib_a304 as L  # noqa: E402
import a336_common as C  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.paired import cohort_states, pseudo_waveforms  # noqa: E402

ENTRY = "a304"
BASE = L.sh.WF
OUTA = PAPER / "exploratory/A304_cancer_gain_departure/output"
STORED = pd.read_csv(OUTA / "score/tests.csv").set_index(["contrast", "readout"])
MDE = pd.read_csv(OUTA / "power/mde.csv").set_index("contrast")
POW = pd.read_csv(OUTA / "power/power.csv")
ROWS = [("A304_H304.1_gain", "H304.1", "gain", "Q_G", -1, "Q_G: log ratio of gain, high vs low proliferation"),
        ("A304_H304.1_departure", "H304.1", "departure", "Q_R", 1, "Q_R: log ratio of departure, high vs low proliferation"),
        ("A304_H304.3", "H304.3", "departure", "Q_R", 1, "Q_R TP53 altered vs intact (5 cancers)"),
        ("A304_H304.4", "H304.4", "departure", "Q_R", 1, "Q_R luminal A vs other (BRCA)")]
N_BOOT = 200


def log(*a):
    print(time.strftime("%F %T"), *a, flush=True)


def margin(contrast):
    k = MDE.loc[contrast, "mde_k_D"]
    r = POW[(POW.contrast == contrast) & np.isclose(POW.k_D, k) & np.isclose(POW.k_G, 1.0)]
    assert len(r) == 1
    return float(r.mean_Q_R.iloc[0])


def raw_states(samples, wf):
    """Per cancer raw (b_case, v_case, b_ref, v_ref) on the waveform (as q3_lib.tcga_states, without the gate)."""
    out = {}
    for c in L.sh.CANCERS:
        case, ref = L.sh.pair_frames(samples, c)
        if len(case) < 10:
            continue
        out[c] = dict(raw=cohort_states(case, ref, True, wf))
    return out


def observed(states, ann, contrast):
    cells = L.design(states, ann, contrast)
    vals = [L.cell_values(c, c["lab"][None, :]) for c in cells]
    Q = L.pooled(cells, vals)
    return {k: float(v[0]) for k, v in Q.items()}


def set_values_w(c, ann_c, n_sets, samples, transform):
    """A.set_values with the pseudo-waveforms (drawn on the 14-gene base, original seeds) transformed to `transform(w)`."""
    limit_threads(1)
    case, ref = L.sh.pair_frames(samples, c)
    ct, cn = L.sh.tcga_full_counts(c, case.index)
    assert (ct.index == case.index).all() and (cn.index == case.index).all()
    lm = np.log(cn.div(cn.sum(axis=1), axis=0).mean() + 1e-12)[cn.mean() >= 1]
    dt, dn = ct.sum(axis=1), cn.sum(axis=1)
    ci = L.CANCERS.index(c)
    out = {}
    for kind, stream in (("matched", 7), ("responsive", 8)):
        ws = pseudo_waveforms(lm, n_sets, np.random.default_rng([304, stream, ci]), BASE, kind=kind)
        rows = []
        for w in ws:
            w = transform(w)
            g = list(w.genes)
            raw = cohort_states(ct[g].assign(depth=dt), cn[g].assign(depth=dn), True, w)
            r = {}
            for h in L.CONTRASTS:
                cell = L.cell_from_raw(c, raw, ann_c, h)
                if cell is not None:
                    r[h] = (cell["n"], L.cell_values(cell, cell["lab"][None, :]))
            rows.append(r)
        out[kind] = rows
    return out


def readouts(states, ann, samples, n_sets, perms, jobs, transform):
    """Full A304 readout table (as score_a304.main) for the given states and companion transform."""
    qobs, qnull = {}, {}
    for hi, h in enumerate(L.CONTRASTS):
        cells = L.design(states, ann, h)
        rng = np.random.default_rng([304, hi])
        labs = [np.vstack([c["lab"], L.perm_labels(c["lab"], perms, rng)]) for c in cells]
        Q = L.pooled(cells, [L.cell_values(c, P) for c, P in zip(cells, labs)])
        qobs[h] = {k: float(v[0]) for k, v in Q.items()}
        qnull[h] = {k: v[1:] for k, v in Q.items()}
    per = dict(zip(L.CANCERS, Parallel(n_jobs=jobs)(
        delayed(set_values_w)(c, ann[c], n_sets, samples, transform) for c in L.CANCERS)))
    rows = []
    for h in L.CONTRASTS:
        comp = {k: A.pool_sets(per, k, h, n_sets) for k in ("matched", "responsive")}
        for ro, key, side in (("gain", "Q_G", "less"), ("departure", "Q_R", "greater")):
            nl, o = qnull[h][key], qobs[h][key]
            p = L.p_less(nl, o) if side == "less" else L.p_greater(nl, o)
            row = dict(contrast=h, readout=ro, statistic=o, p=p, role=A.TESTS.get((h, ro), "descriptive"))
            for k, d in comp.items():
                s = d[key].dropna()
                row[f"{k}_q05"], row[f"{k}_q95"] = float(s.quantile(0.05)), float(s.quantile(0.95))
                row[f"{k}_pct"] = float((s < o).mean() * 100)
            row["clock_specific"] = bool(o < row["matched_q05"]) if ro == "gain" else bool(o > row["matched_q95"])
            rows.append(row)
    rd = pd.DataFrame(rows)
    for fam in ("primary", "secondary"):
        m = rd.role == fam
        hp = A.holm(dict(zip(rd.index[m], rd.p[m])))
        rd.loc[list(hp), "p_holm"] = list(hp.values())
    rd["detected"] = (rd.p_holm < 0.05) & rd.clock_specific
    return rd.set_index(["contrast", "readout"])


def boot_Q(states, ann, idx):
    """Pooled Q_G, Q_R per bootstrap draw, per contrast: idx[h] = list (draws) of list (cells) of index arrays."""
    out = {}
    for h in L.CONTRASTS:
        cells = L.design(states, ann, h)
        res = []
        for draw in idx[h]:
            cb, vb = [], []
            for c, i in zip(cells, draw):
                c2 = dict(c, zt=c["zt"][i], vt=c["vt"][i], zn=c["zn"][i], vn=c["vn"][i], lab=c["lab"][i])
                cb.append(c2)
                vb.append(L.cell_values(c2, c2["lab"][None, :]))
            Q = L.pooled(cb, vb)
            res.append((float(Q["Q_G"][0]), float(Q["Q_R"][0])))
        out[h] = np.array(res)
    return out


def boot_idx(states, ann, n_boot):
    rb = np.random.default_rng([336, 1])
    idx = {}
    for h in L.CONTRASTS:
        cells = L.design(states, ann, h)
        idx[h] = [[rb.integers(0, c["n"], c["n"]) for c in cells] for _ in range(n_boot)]
    return idx


def trio_gate(samples, trio_wf, n_perm=2000):
    """Pooled H138.4-style pairing of trio tumour / normal states and its within-cancer permutation p."""
    rng = np.random.default_rng([336, 2])
    st, nulls, w = [], [], []
    for c in L.sh.CANCERS:
        case, ref = L.sh.pair_frames(samples, c)
        if len(case) < 10:
            continue
        bt, vt, bn, vn = cohort_states(case, ref, True, trio_wf)
        za, zb = L.q.centred_complex(bt), L.q.centred_complex(bn)
        nl = L.q.perm_null(za, vt, zb, vn, n_perm, rng)
        if not L.q.estimable(za, vt, zb, vn, nl):
            continue
        st.append(L.q.corr_stat(za, vt, zb, vn)); nulls.append(nl); w.append(len(za))
    obs = float(L.q.pool(np.array(st)[:, None], w)[0])
    nul = L.q.pool(np.array(nulls), w)
    return obs, float((1 + np.sum(nul >= obs)) / (1 + len(nul)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--jobs", type=int, default=4)
    a = ap.parse_args()
    chk = a.check
    limit_threads(1)
    n_sets, perms, n_boot = (3, 20, 5) if chk else (200, 2000, N_BOOT)
    samples = L.sh.load_tcga()
    ann = L.annotations(samples)
    core, trio = C.core(BASE), C.trio(BASE)
    # ---- stage A: 14 genes, reproduce -----------------------------------------------------------------------------
    s14 = raw_states(samples, BASE)
    if chk:
        ref14 = L.q.tcga_states()
        for c in s14:
            assert all(np.allclose(x, y) for x, y in zip(s14[c]["raw"], ref14[c]["raw"])), c
    ident = lambda w: w  # noqa: E731
    r14 = readouts(s14, ann, samples, n_sets, perms, a.jobs, ident)
    if not chk:
        rep = {}
        for (h, ro), st in STORED.iterrows():
            if (h, ro) not in r14.index or st.role == "descriptive" and h != "H304.1":
                continue
            new = r14.loc[(h, ro)]
            ok = all(np.isclose(new[k], st[k], atol=1e-6) for k in ("statistic", "p", "matched_q05", "matched_q95", "responsive_pct"))
            ok &= bool(new.detected == st.detected)
            rep[(h, ro)] = ok
        log("reproduction", rep)
        if not all(rep.values()):
            rows = [dict(id=i, figure_panel="Fig 5/6", statistic=s, direction=d, value_14_stored=STORED.loc[(h, ro), "statistic"],
                         value_14=r14.loc[(h, ro), "statistic"], reproduced=False,
                         notes=f"14-gene rerun does not reproduce tests.csv: {rep}", script="exploratory/A336_two_readouts/scripts/p1_a304.py")
                    for i, h, ro, k, d, s in ROWS]
            C.write_rows(ENTRY, rows)
            raise SystemExit(f"A304 not reproduced: {rep}")
    log("14-gene stage done")
    # ---- stage B: core ---------------------------------------------------------------------------------------------
    sc = raw_states(samples, core)
    rc = readouts(sc, ann, samples, n_sets, perms, a.jobs, lambda w: C.restrict(w, BASE, core))
    log("core readouts done")
    # bootstrap (shared indices)
    idx = boot_idx(s14, ann, n_boot)
    b14, bc = boot_Q(s14, ann, idx), boot_Q(sc, ann, idx)
    # gene drops
    dl, tl = C.drops(BASE)
    if chk:
        dl = [x for x in dl if x[0] == tl] + dl[:2]
    dv = {h: {"Q_G": {}, "Q_R": {}} for h in L.CONTRASTS}
    t0 = time.time()
    for n, (lab, w) in enumerate(dl):
        sd = raw_states(samples, w)
        for h in L.CONTRASTS:
            q = observed(sd, ann, h)
            dv[h]["Q_G"][lab], dv[h]["Q_R"][lab] = q["Q_G"], q["Q_R"]
        if n % 40 == 0:
            log("drops", n, len(dl), f"{(time.time() - t0) / 60:.1f} min")
    wfs = dict(dl)
    for h in L.CONTRASTS:
        assert np.isclose(dv[h]["Q_G"][tl], observed(sc, ann, h)["Q_G"], equal_nan=True)
        assert np.isclose(dv[h]["Q_R"][tl], observed(sc, ann, h)["Q_R"], equal_nan=True)
    # trio alone
    tobs, tp = trio_gate(samples, trio, 50 if chk else 2000)
    gate = bool(tp < 0.05)
    st_t = raw_states(samples, trio) if gate else None
    if chk:
        print("check ok")
        return
    rows = []
    for i, h, ro, key, d, stat in ROWS:
        o14, oc = r14.loc[(h, ro)], rc.loc[(h, ro)]
        v14, vc = o14.statistic, oc.statistic
        k = 0 if key == "Q_G" else 1
        diff = b14[h][:, k] - bc[h][:, k]
        lo_b, hi_b = np.nanquantile(bc[h][:, k], [0.05, 0.95])
        pr = C.drop_percentiles(dv[h][key], tl, wfs, direction=d)
        if ro == "gain":
            conds = [bool(oc.statistic < 0 and oc.p_holm < 0.05), bool(oc.clock_specific)]
            status = "holds" if all(conds) else "at the line" if any(conds) else "falls"
            crit = "Q_G < 0, Holm p < 0.05 (family of 2), and below matched q05 (companions restricted to core)"
            note = (f"core p {oc.p:.4f}, Holm {oc.p_holm:.4f}, matched q05 {oc.matched_q05:.4f}, responsive pct {oc.responsive_pct:.1f}; "
                    f"core 90 % pair-bootstrap {lo_b:.3f} to {hi_b:.3f}")
        else:
            mg = margin(h)
            not_det = not bool(oc.detected)
            status = "holds" if (not_det and hi_b < mg) else "falls"
            crit = (f"departure rise not detected (Holm p >= 0.05 or not above matched q95) and 90 % bootstrap upper of Q_R < margin "
                    f"{mg:.3f} (Q_R at registered MDE k_D {MDE.loc[h, 'mde_k_D']})")
            note = (f"core detected {bool(oc.detected)}, p {oc.p:.4f}, Holm {oc.p_holm:.4f}, matched q95 {oc.matched_q95:.4f}; core 90 % "
                    f"pair-bootstrap {lo_b:.3f} to {hi_b:.3f}; margin from 14-gene power simulation (not rerun)")
        if gate:
            tv = observed(st_t, ann, h)[key]
        rows.append(dict(id=i, figure_panel="Fig 5/6", statistic=stat, direction=d, value_14_stored=STORED.loc[(h, ro), "statistic"],
                         value_14=v14, reproduced=True, value_core=vc, core_criterion=crit, core_status=status,
                         diff_14_minus_core=v14 - vc, diff_lo90=float(np.nanquantile(diff, 0.05)),
                         diff_hi90=float(np.nanquantile(diff, 0.95)), diff_draws="200 pair bootstraps within cancer (rng [336, 1]), shared",
                         **pr, trio_gate_passes=gate, trio_alone_value=tv if gate else np.nan,
                         reading=C.reading(status, pr["trio_drop_below_q10"]),
                         blind="core blind", notes=note + f"; trio pairing gate p {tp:.4f} (pooled pairing {tobs:.3f})",
                         script="exploratory/A336_two_readouts/scripts/p1_a304.py"))
    df = C.write_rows(ENTRY, rows)
    pd.set_option("display.width", 250)
    print(df.drop(columns=["script", "notes", "core_criterion", "diff_draws"]).round(4).T.to_string())
    for r in rows:
        print(r["id"], r["notes"])


if __name__ == "__main__":
    main()
