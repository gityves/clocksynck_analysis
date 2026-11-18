"""A304 scoring (A1, A3): gain and departure of the tumour on its own normal sample by tumour state; tumour coupling to
the body's time per cancer.

Registered in paper/docs/registered_experiments.md (A304 and its design addendum). Tests: H304.1 gain and departure
(primary, Holm over 2); H304.3 and H304.4 departure (secondary, Holm over 2); every other readout descriptive.
Companions: 200 expression-matched (decisive, A285) and 200 responsive (reported) pseudo-clock sets per cancer on the
genome-wide counts of the pairs, pooled per set index, observed labels.

    python score_a304.py [--perms 2000] [--sets 200] [--boot 1000] [--jobs 4]
Outputs (output/score/): tests.csv, readouts.csv (every contrast x readout), per_cancer.csv, companions.csv (per set),
sensitivity_purity.csv, a3_per_cancer.csv.
"""

from __future__ import annotations

import argparse
import time

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

import lib_a304 as L
from clocksync.paired import cohort_states, pseudo_waveforms

HOURS = 24 / (2 * np.pi)
T0 = time.time()
TESTS = {("H304.1", "gain"): "primary", ("H304.1", "departure"): "primary",
         ("H304.3", "departure"): "secondary", ("H304.4", "departure"): "secondary"}


def log(*a):
    print(f"[{(time.time() - T0) / 60:6.1f} min]", *a, flush=True)


def holm(p: dict) -> dict:
    k = sorted(p, key=p.get)
    out, run = {}, 0.0
    for i, key in enumerate(k):
        run = max(run, min(1.0, (len(k) - i) * p[key]))
        out[key] = run
    return out


def set_values(c: str, ann_c: pd.DataFrame, n_sets: int, samples) -> dict:
    """Per contrast, per set index: the cell values (case / ref G, R, VN) of cancer c for matched and responsive sets."""
    from clocksync.config import limit_threads
    limit_threads(1)
    case, ref = L.sh.pair_frames(samples, c)
    ct, cn = L.sh.tcga_full_counts(c, case.index)
    assert (ct.index == case.index).all() and (cn.index == case.index).all()
    lm = np.log(cn.div(cn.sum(axis=1), axis=0).mean() + 1e-12)[cn.mean() >= 1]
    dt, dn = ct.sum(axis=1), cn.sum(axis=1)
    ci = L.CANCERS.index(c)
    out = {}
    for kind, stream in (("matched", 7), ("responsive", 8)):
        ws = pseudo_waveforms(lm, n_sets, np.random.default_rng([304, stream, ci]), L.sh.WF, kind=kind)
        rows = []
        for w in ws:
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


def pool_sets(per_cancer: dict, kind: str, contrast: str, n_sets: int) -> pd.DataFrame:
    """Q_G, Q_R, Q_N per set index pooled over the cancers in the contrast."""
    res = []
    for j in range(n_sets):
        cells, vals = [], []
        for c, v in per_cancer.items():
            r = v[kind][j].get(contrast)
            if r is not None:
                cells.append(dict(n=r[0])); vals.append(r[1])
        Q = L.pooled(cells, vals)
        res.append({k: float(x[0]) for k, x in Q.items()})
    return pd.DataFrame(res)


def offset_h(Ccase, Cref) -> float:
    return float(np.angle(Ccase * np.conj(Cref)) * HOURS)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--perms", type=int, default=2000)
    ap.add_argument("--sets", type=int, default=200)
    ap.add_argument("--boot", type=int, default=1000)
    ap.add_argument("--jobs", type=int, default=4)
    a = ap.parse_args()
    from clocksync.config import limit_threads
    limit_threads(1)
    out = L.OUT / "score"
    out.mkdir(parents=True, exist_ok=True)
    samples = L.sh.load_tcga()
    states = L.q.tcga_states()
    ann = L.annotations(samples)
    rel = L.reliability()

    # ---- observed, permutations, per cancer
    readouts, percan, Qobs, Qnull = [], [], {}, {}
    for hi, h in enumerate(L.CONTRASTS):
        cells = L.design(states, ann, h)
        rng = np.random.default_rng([304, hi])
        labs = [np.vstack([c["lab"], L.perm_labels(c["lab"], a.perms, rng)]) for c in cells]
        vals = [L.cell_values(c, P) for c, P in zip(cells, labs)]
        Q = L.pooled(cells, vals)
        Qobs[h] = {k: float(v[0]) for k, v in Q.items()}
        Qnull[h] = {k: v[1:] for k, v in Q.items()}
        for c, v in zip(cells, vals):
            percan.append(dict(contrast=h, cancer=c["cancer"], n=c["n"], n_case=int((c["lab"] == 1).sum()),
                               n_ref=int((c["lab"] == 0).sum()), G_case=v[1]["G"][0], G_ref=v[0]["G"][0],
                               R_case=v[1]["R"][0], R_ref=v[0]["R"][0], VN_case=v[1]["VN"][0], VN_ref=v[0]["VN"][0],
                               offset_h=offset_h(v[1]["C"][0], v[0]["C"][0])))
        num = sum(c["n"] * v[1]["C"][0] / abs(v[1]["C"][0]) * np.conj(v[0]["C"][0] / abs(v[0]["C"][0]))
                  for c, v in zip(cells, vals))
        Qobs[h]["offset_h"] = float(np.angle(num) * HOURS)
        log(h, "observed computed")

    # ---- companions
    per_cancer = dict(zip(L.CANCERS, Parallel(n_jobs=a.jobs)(
        delayed(set_values)(c, ann[c], a.sets, samples) for c in L.CANCERS)))
    log("companion sets computed")
    comp = []
    for h in L.CONTRASTS:
        for kind in ("matched", "responsive"):
            d = pool_sets(per_cancer, kind, h, a.sets)
            d.insert(0, "kind", kind); d.insert(0, "contrast", h); d.insert(2, "set", range(len(d)))
            comp.append(d)
    comp = pd.concat(comp, ignore_index=True)
    comp.to_csv(out / "companions.csv", index=False)

    # ---- readouts and tests
    for h in L.CONTRASTS:
        o, nl = Qobs[h], Qnull[h]
        for ro, key, side in (("gain", "Q_G", "less"), ("departure", "Q_R", "greater"), ("host", "Q_N", "two")):
            pl, pg = L.p_less(nl[key], o[key]), L.p_greater(nl[key], o[key])
            p = pl if side == "less" else pg if side == "greater" else min(1.0, 2 * min(pl, pg))
            row = dict(contrast=h, readout=ro, statistic=o[key], p=p, null_mean=float(np.nanmean(nl[key])),
                       null_sd=float(np.nanstd(nl[key])), n_perm=int(np.isfinite(nl[key]).sum()),
                       role=TESTS.get((h, ro), "descriptive"), pooled_offset_h=o["offset_h"])
            for kind in ("matched", "responsive"):
                s = comp[(comp.contrast == h) & (comp.kind == kind)][key].dropna()
                row[f"{kind}_median"] = float(s.median())
                row[f"{kind}_q05"], row[f"{kind}_q95"] = float(s.quantile(0.05)), float(s.quantile(0.95))
                row[f"{kind}_pct"] = float((s < o[key]).mean() * 100)
                row[f"{kind}_n"] = len(s)
            if ro == "gain":
                row["clock_specific"] = o[key] < row["matched_q05"]
            elif ro == "departure":
                row["clock_specific"] = o[key] > row["matched_q95"]
            readouts.append(row)
    rd = pd.DataFrame(readouts)
    for fam in ("primary", "secondary"):
        m = rd.role == fam
        hp = holm(dict(zip(rd.index[m], rd.p[m])))
        rd.loc[list(hp), "p_holm"] = list(hp.values())
    rd["detected"] = (rd.p_holm < 0.05) & rd.clock_specific.fillna(False).astype(bool)
    rd.to_csv(out / "readouts.csv", index=False)
    rd[rd.role != "descriptive"].to_csv(out / "tests.csv", index=False)
    pd.DataFrame(percan).to_csv(out / "per_cancer.csv", index=False)
    log("readouts\n", rd[["contrast", "readout", "statistic", "p", "p_holm", "matched_q05", "matched_q95",
                          "responsive_pct", "clock_specific", "role", "detected"]].round(4).to_string())

    # ---- sensitivity: H304.1 with cells = cancer x purity half (arms >= 5 per cell)
    cells = []
    for c in L.design(states, ann, "H304.1"):
        ph = L.purity_half(c["ann"])
        for k in (0, 1):
            m = ph == k
            if min((c["lab"][m] == 1).sum(), (c["lab"][m] == 0).sum()) < 5:
                continue
            cells.append(dict(c, zt=c["zt"][m], vt=c["vt"][m], zn=c["zn"][m], vn=c["vn"][m], lab=c["lab"][m],
                              n=int(m.sum()), cancer=f"{c['cancer']}|{k}"))
    rng = np.random.default_rng([304, 11])
    labs = [np.vstack([c["lab"], L.perm_labels(c["lab"], a.perms, rng)]) for c in cells]
    Q = L.pooled(cells, [L.cell_values(c, P) for c, P in zip(cells, labs)])
    pd.DataFrame([dict(cells=len(cells), pairs=sum(c["n"] for c in cells), Q_G=Q["Q_G"][0],
                       p_gain=L.p_less(Q["Q_G"][1:], Q["Q_G"][0]), Q_R=Q["Q_R"][0],
                       p_dep=L.p_greater(Q["Q_R"][1:], Q["Q_R"][0]), Q_N=Q["Q_N"][0])]).to_csv(
        out / "sensitivity_purity.csv", index=False)
    imb = pd.DataFrame([dict(cancer=c["cancer"], purity_case=float(c["ann"].purity[c["lab"] == 1].mean()),
                             purity_ref=float(c["ann"].purity[c["lab"] == 0].mean()))
                        for c in L.design(states, ann, "H304.1")])
    imb.to_csv(out / "purity_imbalance.csv", index=False)

    # ---- A3: per cancer, all pairs
    a3 = []
    for ci, c in enumerate(L.CANCERS):
        u = states[c]
        zt, vt, zn, vn = u["za"], u["vA"], u["zb"], u["vB"]
        n = len(zt)

        def stats(i):
            vN, vT, C = L.arm_moments(zt[i], vt[i], zn[i], vn[i], np.ones((1, len(i)), bool))
            return float(abs(C[0]) / np.sqrt(vN[0] * vT[0])) if vN[0] > 0 and vT[0] > 0 else np.nan, \
                float(abs(C[0]) / vN[0]) if vN[0] > 0 else np.nan
        corr, G = stats(np.arange(n))
        rng = np.random.default_rng([304, 21, ci])
        bs = np.array([stats(rng.integers(0, n, n)) for _ in range(a.boot)])
        r = rel[c] if c in L.TISSUE else np.nan
        a3.append(dict(cancer=c, n=n, rel=r, tissue=L.TISSUE.get(c, ""), corr=corr,
                       corr_lo=np.nanquantile(bs[:, 0], 0.05), corr_hi=np.nanquantile(bs[:, 0], 0.95), G=G,
                       G_lo=np.nanquantile(bs[:, 1], 0.05), G_hi=np.nanquantile(bs[:, 1], 0.95),
                       c_T=corr ** 2 / r, c_T_lo=np.nanquantile(bs[:, 0], 0.05) ** 2 / r,
                       c_T_hi=np.nanquantile(bs[:, 0], 0.95) ** 2 / r, gain_body=G / r,
                       gain_body_lo=np.nanquantile(bs[:, 1], 0.05) / r, gain_body_hi=np.nanquantile(bs[:, 1], 0.95) / r))
    a3 = pd.DataFrame(a3).sort_values("c_T", ascending=False)
    a3.to_csv(out / "a3_per_cancer.csv", index=False)
    log("A3\n", a3.round(3).to_string())


if __name__ == "__main__":
    main()
