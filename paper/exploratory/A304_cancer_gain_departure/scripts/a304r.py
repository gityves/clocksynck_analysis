"""A304r: replication of A304 H304.1 (proliferation gain loss, no departure rise) in CPTAC-3 tumour / normal pairs.

Registered in paper/docs/registered_experiments.md ("A304r"). Data: the A186 cache (8 types, 713 pairs); states, gates and
G2M scores by `a186_lib`; gain / departure statistics by `lib_a304` (identical code path to A304).

  power  A304 simulator at the CPTAC design (no observed contrast): rejection at 0.05 for k_G 0.85 / 0.7 / 0.63 / 0.5 and
         k_D 1.5 / 2                                                     -> output/a304r/{calibration,design,power}.csv
  score  refuses to run without power.csv; observed Q_G / Q_R / Q_N, permutation p, matched / responsive companions,
         per type, core loop (gain), purity proxy                       -> output/a304r/{tests,per_type,companions,core,purity}.csv
    python -I a304r.py --stage power|score [--jobs 6]
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "A186_cptac_c1" / "scripts"))
sys.path.insert(0, str(HERE.parents[1] / "A336_two_readouts" / "scripts"))
import a186_lib as K  # noqa: E402
import lib_a304 as L  # noqa: E402
import power_a304 as PW  # noqa: E402

OUT = L.OUT / "a304r"
SEED = 3042
TISSUE = {"LUAD": "Lung", "LSCC": "Lung", "CCRCC": "Kidney - Cortex", "SKCM": "Skin - Not Sun Exposed (Suprapubic)",
          "STAD": "Stomach"}
SCEN = [("null", 1.0, 1.0)] + [(f"gain x{k}", k, 1.0) for k in (0.85, 0.7, 0.63, 0.5)] + [(f"departure x{k}", 1.0, k) for k in (1.5, 2.0)]
STROMA, IMMUNE = ["COL1A1", "COL1A2", "DCN", "LUM"], ["PTPRC"]
T0 = time.time()


def log(*a):
    print(f"[{(time.time() - T0) / 60:6.1f} min]", *a, flush=True)


def reliability() -> dict:
    a = pd.read_csv(L.A104).set_index("tissue").coupling
    rel = {t: float(a[s]) for t, s in TISSUE.items()}
    med = float(np.median(list(rel.values())))
    return {t: rel.get(t, med) for t in K.TYPES}


def data():
    """Per type: states (A186), annotation frame (prolif = G2M of the tumour), counts."""
    st, ann = {}, {}
    for i, tp in enumerate(K.TYPES):
        s = K.states(tp, i)
        assert s["gated"], tp
        g, _ = K.labels(tp)
        st[tp] = s
        ann[tp] = pd.DataFrame({"prolif": np.asarray(g, float)}, index=s["cases"])
    return st, ann


def cells_of(st, ann, contrast="H304.1"):
    cells = [L.cell_from_raw(tp, st[tp]["raw"], ann[tp], contrast) for tp in K.TYPES]
    return [c for c in cells if c is not None]


def stage_power(a):
    st, ann = data()
    cal = PW.calibrate(st, reliability())
    cal.to_csv(OUT / "calibration.csv")
    cells = cells_of(st, ann)
    pd.DataFrame([dict(type=c["cancer"], n=c["n"], n_high=int((c["lab"] == 1).sum()), n_low=int((c["lab"] == 0).sum()))
                  for c in cells]).to_csv(OUT / "design.csv", index=False)
    rows = []
    for si, (name, kG, kD) in enumerate(SCEN):
        res = Parallel(n_jobs=a.jobs)(delayed(PW.one_rep)(cells, cal, kG, kD, a.perms, [SEED, 2, si, r]) for r in range(a.reps))
        r = pd.DataFrame(res)
        rows.append(dict(scenario=name, k_G=kG, k_D=kD, reps=len(r), power_gain_05=float((r.p_gain < 0.05).mean()),
                         power_dep_05=float((r.p_dep < 0.05).mean()), host_flag_05=float((r.p_host < 0.05).mean()),
                         mean_Q_G=float(r.Q_G.mean()), sd_Q_G=float(r.Q_G.std()), mean_Q_R=float(r.Q_R.mean())))
        log(name, {k: round(v, 3) for k, v in rows[-1].items() if k.startswith(("power", "mean"))})
    P = pd.DataFrame(rows)
    P["scored_H304r_1"] = bool(P.loc[P.k_G == 0.63, "power_gain_05"].iloc[0] >= 0.8)
    P.to_csv(OUT / "power.csv", index=False)
    log("power\n", P.round(3).to_string())


def set_values(tp, i, ann_c, n_sets):
    from clocksync.config import limit_threads
    from clocksync.paired import cohort_states, pseudo_waveforms
    limit_threads(1)
    cases, ct, cn = K.load_type(tp)
    lm = np.log(cn.div(cn.sum(axis=1), axis=0).mean() + 1e-12)[cn.mean() >= 1]
    dt, dn = ct.sum(axis=1), cn.sum(axis=1)
    out = {}
    for kind, stream in (("matched", 7), ("responsive", 8)):
        ws = pseudo_waveforms(lm, n_sets, np.random.default_rng([SEED, stream, i]), L.sh.WF, kind=kind)
        vals = []
        for w in ws:
            g = list(w.genes)
            raw = cohort_states(ct[g].assign(depth=dt), cn[g].assign(depth=dn), True, w)
            cell = L.cell_from_raw(tp, raw, ann_c, "H304.1")
            vals.append(None if cell is None else (cell["n"], L.cell_values(cell, cell["lab"][None, :])))
        out[kind] = vals
    return out


def pooled_sets(per_type, kind, n_sets):
    res = []
    for j in range(n_sets):
        cells, vals = [], []
        for v in per_type.values():
            r = v[kind][j]
            if r is not None:
                cells.append(dict(n=r[0]))
                vals.append(r[1])
        Q = L.pooled(cells, vals)
        res.append({k: float(x[0]) for k, x in Q.items()})
    return pd.DataFrame(res)


def score_from(cells, n_perm, seed):
    rng = np.random.default_rng(seed)
    labs = [np.vstack([c["lab"], L.perm_labels(c["lab"], n_perm, rng)]) for c in cells]
    vals = [L.cell_values(c, P) for c, P in zip(cells, labs)]
    return L.pooled(cells, vals), vals


def gene_score(ct, genes):
    d = ct.sum(axis=1)
    g = [x for x in genes if x in ct.columns]
    lc = np.log2(ct[g].div(d, axis=0) * 1e6 + 1)
    return ((lc - lc.mean()) / lc.std(ddof=0)).mean(axis=1)


def stage_score(a):
    pw = pd.read_csv(OUT / "power.csv")
    scored = bool(pw.scored_H304r_1.iloc[0])
    st, ann = data()
    cells = cells_of(st, ann)
    Q, vals = score_from(cells, a.perms, [SEED, 1])
    obs = {k: float(v[0]) for k, v in Q.items()}
    nul = {k: v[1:] for k, v in Q.items()}
    per_type = dict(zip(K.TYPES, Parallel(n_jobs=a.jobs)(delayed(set_values)(tp, i, ann[tp], a.sets) for i, tp in enumerate(K.TYPES))))
    comp = pd.concat([pooled_sets(per_type, kind, a.sets).assign(kind=kind, set=range(a.sets)) for kind in ("matched", "responsive")],
                     ignore_index=True)
    comp.to_csv(OUT / "companions.csv", index=False)
    rows = []
    for ro, key, side in (("gain", "Q_G", "less"), ("departure", "Q_R", "greater"), ("host", "Q_N", "two")):
        pl, pg = L.p_less(nul[key], obs[key]), L.p_greater(nul[key], obs[key])
        p = pl if side == "less" else pg if side == "greater" else min(1.0, 2 * min(pl, pg))
        r = dict(readout=ro, statistic=obs[key], ratio=float(np.exp(obs[key])), p=p, null_mean=float(np.nanmean(nul[key])),
                 null_sd=float(np.nanstd(nul[key])),
                 hypothesis={"gain": "H304r.1" if scored else "descriptive (power)", "departure": "H304r.2"}.get(ro, "descriptive"))
        for kind in ("matched", "responsive"):
            s = comp[comp.kind == kind][key].dropna()
            r.update({f"{kind}_median": float(s.median()), f"{kind}_q05": float(s.quantile(.05)), f"{kind}_q95": float(s.quantile(.95)),
                      f"{kind}_pct": float((s < obs[key]).mean() * 100), f"{kind}_n": len(s)})
        if ro == "gain":
            r["pass"] = bool(scored and p < 0.05 and obs[key] < r["matched_q05"])
        elif ro == "departure":
            r["pass"] = bool(p < 0.05 and obs[key] > r["matched_q95"])
        rows.append(r)
    T = pd.DataFrame(rows)
    T.to_csv(OUT / "tests.csv", index=False)
    log("tests\n", T.round(4).to_string())
    hours = 24 / (2 * np.pi)
    pt = [dict(type=c["cancer"], n=c["n"], n_high=int((c["lab"] == 1).sum()), n_low=int((c["lab"] == 0).sum()), G_high=v[1]["G"][0],
               G_low=v[0]["G"][0], R_high=v[1]["R"][0], R_low=v[0]["R"][0], VN_high=v[1]["VN"][0], VN_low=v[0]["VN"][0],
               offset_h=float(np.angle(v[1]["C"][0] * np.conj(v[0]["C"][0])) * hours)) for c, v in zip(cells, vals)]
    pd.DataFrame(pt).to_csv(OUT / "per_type.csv", index=False)
    log("per type\n", pd.DataFrame(pt).round(3).to_string())
    # core loop (gain), descriptive
    import a336_common as C
    core = C.core(K.q.WF)
    cc = []
    for i, tp in enumerate(K.TYPES):
        cases, ct, cn = K.load_type(tp)
        raw = K.q.cohort_states(K.clock_frame(ct), K.clock_frame(cn), True, core)
        c = L.cell_from_raw(tp, raw, ann[tp], "H304.1")
        if c is not None:
            cc.append(c)
    Qc, _ = score_from(cc, a.perms, [SEED, 3])
    pd.DataFrame([dict(readout="gain, core loop", statistic=float(Qc["Q_G"][0]), ratio=float(np.exp(Qc["Q_G"][0])),
                       p=L.p_less(Qc["Q_G"][1:], Qc["Q_G"][0]), departure=float(Qc["Q_R"][0]),
                       p_departure=L.p_greater(Qc["Q_R"][1:], Qc["Q_R"][0]))]).to_csv(OUT / "core.csv", index=False)
    log("core", float(Qc["Q_G"][0]), L.p_less(Qc["Q_G"][1:], Qc["Q_G"][0]))
    # purity proxy between arms
    pur = []
    for c in cells:
        cases, ct, _ = K.load_type(c["cancer"])
        ix = [cases.index(x) for x in c["ann"].index]
        s, im = gene_score(ct, STROMA).to_numpy()[ix], gene_score(ct, IMMUNE).to_numpy()[ix]
        pur.append(dict(type=c["cancer"], stromal_high=float(s[c["lab"] == 1].mean()), stromal_low=float(s[c["lab"] == 0].mean()),
                        immune_high=float(im[c["lab"] == 1].mean()), immune_low=float(im[c["lab"] == 0].mean())))
    pd.DataFrame(pur).to_csv(OUT / "purity.csv", index=False)
    log("purity proxy\n", pd.DataFrame(pur).round(3).to_string())


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stage", choices=("power", "score"), required=True)
    ap.add_argument("--jobs", type=int, default=6)
    ap.add_argument("--reps", type=int, default=500)
    ap.add_argument("--perms", type=int, default=2000)
    ap.add_argument("--sets", type=int, default=200)
    a = ap.parse_args()
    from clocksync.config import limit_threads
    limit_threads(1)
    OUT.mkdir(parents=True, exist_ok=True)
    if a.stage == "power":
        a.perms = 500
        stage_power(a)
    else:
        stage_score(a)
    log("done")


if __name__ == "__main__":
    main()
