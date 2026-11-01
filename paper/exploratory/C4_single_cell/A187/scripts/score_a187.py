"""A187 scoring (registered_experiments.md, A187): is the clock of malignant cells phase-locked to the cell cycle?
Pseudobulks from prep_a187.py. `clocksync.paired` (unit = patient, cohort = study, reference = NC), GTEx waveform.

Order (stop rules binding): prep tables, gates -> POWER (MDE = 2.49 x SD of the label-swap null of the pooled G2M-NC offset;
STOP if > 1.0 h) -> H187.1 + companions -> artefact checks (a), (c) -> H187.2 -> descriptives.
Run: python -I score_a187.py [--n_perm 10000] [--min_pairs 10]
A187b (post hoc): --gate all_nc applies the study gate to the NC pseudobulks of every patient with >= 50 NC cells (column
gate_p_all_NC_patients_descriptive of gates.csv); the contrast then runs on the gated studies only, with paired(gated_only=False)
on samples restricted to those studies (paired's own gate on the contrast patients is then informational only). Outputs *_b.csv
(*_b_mp<k>.csv when min_pairs differs from 10).
"""
from __future__ import annotations
import argparse, sys, time
from pathlib import Path
import numpy as np, pandas as pd

from clocksync import paired as P
from clocksync import stats as S
from clocksync.config import limit_threads

HERE = Path(__file__).resolve()
C4 = HERE.parents[2]
OUT = HERE.parents[1] / "output"
STUDIES = ["GSE178341", "GSE132465", "GSE131907", "GSE176078", "GSE161529", "SRP327758", "GSE181919", "SRP416180"]
PILOT = {"GSE178341", "GSE132465", "GSE131907", "GSE176078", "SRP327758"}
MIN_CELLS = 50
SEED = 20261101
SUF = ""          # output suffix ("_b" for A187b)
GATED = True      # paired(gated_only=...): False when the gate was applied upstream (A187b)
CLOCK = list(P.gtex().genes)
META = ["study", "patient", "bin", "n_cells", "depth"]


def say(*a):
    print(time.strftime("%F %T"), *a, flush=True)


def load():
    """All studies: one row per patient x bin, clock + all genes as columns (NaN where a study lacks a gene)."""
    parts = []
    for s in STUDIES:
        g = pd.read_csv(OUT / f"genome_{s}.csv.gz", index_col=[0, 1])
        pb = pd.read_csv(OUT / f"pseudobulk_clock_{s}.csv.gz").set_index(["patient", "bin"])
        pb.index = pb.index.set_levels([pb.index.levels[0].astype(str), pb.index.levels[1]])
        g.index = g.index.set_levels([g.index.levels[0].astype(str), g.index.levels[1]])
        x = g.astype(float).join(pb[["n_cells", "depth"]]).reset_index()
        x.insert(0, "study", s)
        parts.append(x)
    return pd.concat(parts, ignore_index=True)


def cycling_fraction():
    out = []
    for s in STUDIES:
        d = C4 / ("pilot/output" if s in PILOT else "fits/output")
        t = pd.read_csv(d / f"cells_{s}.csv.gz", usecols=["patient", "cycling_prob"])
        f = t.groupby("patient").cycling_prob.apply(lambda x: float((x >= 0.5).mean()))
        out.append(pd.DataFrame({"study": s, "patient": f.index.astype(str), "frac_cycling": f.to_numpy()}))
    return pd.concat(out)


def contrast(sam, case, ref="NC"):
    """Samples of patients with >= MIN_CELLS cells in both bins, columns condition = bin."""
    w = sam[sam.bin.isin([ref, case])].pivot_table(index=["study", "patient"], columns="bin", values="n_cells", fill_value=0)
    ok = w[(w[ref] >= MIN_CELLS) & (w[case] >= MIN_CELLS)].index
    k = pd.MultiIndex.from_frame(sam[["study", "patient"]])
    x = sam[sam.bin.isin([ref, case]) & k.isin(ok)].copy()
    x["condition"] = x["bin"]
    return x


def pseudo_means(x, ref="NC"):
    """study -> mean log fraction of all genes in the reference (NC) pseudobulks of the contrast patients (mean count >= 1)."""
    out = {}
    for s, g in x[x.condition == ref].groupby("study"):
        c = g.drop(columns=META + ["condition"]).dropna(axis=1, how="all")
        c = c.loc[:, c.sum() > 0]
        out[s] = np.log(c.div(g.depth.to_numpy(float), axis=0).mean() + 1e-12)[c.mean() >= 1]
    return out


def run(x, case, drop=None, min_pairs=10, n_perm=10000, pm=None, resp=False, seed=SEED, gated_only=None):
    return P.paired(x, unit="patient", condition="condition", levels=("NC", case), cohort="study", counts=True,
                    n_perm=n_perm, pseudo_means=pm, pseudo_responsive=resp, drop_genes=drop, min_pairs=min_pairs,
                    seed=seed, gated_only=GATED if gated_only is None else gated_only)


def weights(res):
    keys = str(res.pooled.cohorts).split(";")
    c = res.cohorts.set_index("cohort").loc[keys]
    return keys, 1 / c.se.to_numpy() ** 2


def null_sd(res, n_perm, seed):
    keys, w = weights(res)
    null = P.swap_offset_null([res.data[k] for k in keys], w, n_perm=n_perm, rng=seed)
    return null, float(np.std(null, ddof=1))


def pseudo_offsets(x, case, res, pm, n_pseudo, kind, seed):
    """Pooled offset (weights of the clock fit) of each pseudo set; one set drawn per cohort for each j."""
    keys, w = weights(res)
    wf = P.gtex()
    rng = np.random.default_rng([seed, 7 if kind == "matched" else 8])
    pw, fr = {}, {}
    for k in keys:
        g = x[x.study == k]
        c = g[g.condition == case].set_index("patient")
        r = g[g.condition == "NC"].set_index("patient")
        u = sorted(set(c.index) & set(r.index))
        fr[k] = (c.loc[u].dropna(axis=1, how="all"), r.loc[u].dropna(axis=1, how="all"))
        pw[k] = P.pseudo_waveforms(pm[k], n_pseudo, rng, wf, kind)
    out = []
    for j in range(n_pseudo):
        cs_ = [P.ccorr(*P.cohort_states(fr[k][0], fr[k][1], True, pw[k][j])) for k in keys]
        out.append(P.pooled_offset(cs_, w))
    return np.array(out)


def thin(x, seed):
    """G2M bin of each patient binomially thinned (clock genes) to the NC bin's depth; depth set to the NC depth."""
    rng = np.random.default_rng(seed)
    y = x.copy()
    nc = y[y.condition == "NC"].set_index(["study", "patient"]).depth
    for i in y.index[y.condition != "NC"]:
        d0 = nc.loc[(y.at[i, "study"], y.at[i, "patient"])]
        p = d0 / y.at[i, "depth"]
        if p < 1:
            y.loc[i, CLOCK] = rng.binomial(y.loc[i, CLOCK].to_numpy(int), p)
            y.at[i, "depth"] = d0
    return y


def tertile_stat(states, labels, w_fixed=None, minn=3):
    """|pooled offset| high tertile minus low tertile; per study one projector (all contrast patients), weights = group size."""
    def pooled(g):
        cs_, ws = [], []
        for k, b in states.items():
            i = np.flatnonzero(labels[k] == g)
            if len(i) < minn:
                continue
            cs_.append(P.ccorr(*tuple(v[i] for v in b)))
            ws.append(len(i))
        return abs(P.pooled_offset(cs_, ws)) if cs_ else np.nan
    return pooled("high") - pooled("low")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n_perm", type=int, default=10000)
    ap.add_argument("--n_pseudo", type=int, default=200)
    ap.add_argument("--min_pairs", type=int, default=10)
    ap.add_argument("--gate", choices=["contrast", "all_nc"], default="contrast")
    a = ap.parse_args()
    global SUF, GATED
    if a.gate == "all_nc":
        SUF, GATED = "_b" + ("" if a.min_pairs == 10 else f"_mp{a.min_pairs}"), False
    limit_threads(4)
    sam = load()
    say("samples", sam.shape)

    # ---- 1. prep, gates ------------------------------------------------------------------------------------------
    rows = []
    for s in STUDIES:
        w = sam[sam.study == s].pivot_table(index="patient", columns="bin", values="n_cells", fill_value=0)
        m = pd.read_csv(OUT / f"match_{s}.csv").iloc[0]
        rows.append(dict(study=s, patients=len(w), cells_matched=int(m.matched), match_rate=m.match_rate, cells_NC=int(m.cells_NC),
                         cells_S=int(m.cells_S), cells_G2M=int(m.cells_G2M), S_peak_theta=m.S_peak, G2M_peak_theta=m.G2M_peak,
                         patients_NC_G2M=int(((w.NC >= MIN_CELLS) & (w.G2M >= MIN_CELLS)).sum()),
                         patients_NC_S=int(((w.NC >= MIN_CELLS) & (w.S >= MIN_CELLS)).sum())))
    prep = pd.DataFrame(rows)
    prep.to_csv(OUT / f"prep{SUF}.csv", index=False)
    say("\n" + prep.round(3).to_string())

    xg = contrast(sam, "G2M")
    pm = pseudo_means(xg)
    primary = run(xg, "G2M", min_pairs=a.min_pairs, n_perm=100, pm=None)  # gate + states only here; offsets are not read
    wide = run(xg, "G2M", min_pairs=3, n_perm=100, gated_only=False)
    gates = wide.cohorts[["cohort", "n_pairs", "gate_p", "gated"]].rename(columns={"cohort": "study"}).copy()
    gates["gate_p_primary_call"] = gates.study.map(primary.cohorts.set_index("cohort").gate_p) if len(primary.cohorts) else np.nan
    gates["in_primary"] = gates.study.isin(primary.cohorts.cohort) if len(primary.cohorts) else False
    gates["enters"] = gates.in_primary & gates.study.map(primary.cohorts.set_index("cohort").gated).fillna(False).astype(bool)
    # descriptive only (not used for entry): same gate on the NC pseudobulks of every patient with >= 50 NC cells
    allnc = {}
    for st_, g in sam[(sam.bin == "NC") & (sam.n_cells >= MIN_CELLS)].groupby("study"):
        yr = P.log_expr(g[CLOCK].astype(float), g.depth.to_numpy(float)).to_numpy()
        allnc[st_] = (len(g), P.reference_gate(yr, P.gtex(), 200, np.random.default_rng(SEED))["p"])
    gates["n_patients_NC_ge50"] = gates.study.map({k: v[0] for k, v in allnc.items()})
    gates["gate_p_all_NC_patients_descriptive"] = gates.study.map({k: v[1] for k, v in allnc.items()})
    if a.gate == "all_nc":
        gs = set(gates.study[gates.gate_p_all_NC_patients_descriptive < 0.05])
        gates["gated_all_NC"] = gates.study.isin(gs)
        xg = xg[xg.study.isin(gs)]
        pm = {k: v for k, v in pm.items() if k in gs}
        primary = run(xg, "G2M", min_pairs=a.min_pairs, n_perm=100, pm=None)
        say("A187b gated studies:", sorted(gs))
    gates.to_csv(OUT / f"gates{SUF}.csv", index=False)
    say("\n" + gates.round(4).to_string())

    # ---- 2. POWER ------------------------------------------------------------------------------------------------
    prow = []
    for label, mp in (("primary", a.min_pairs), ("sensitivity_min_pairs_4", 4)):
        r0 = primary if mp == a.min_pairs else run(xg, "G2M", min_pairs=mp, n_perm=100)
        if int(r0.pooled.get("n_cohorts", 0)) == 0:
            prow.append(dict(design=label, min_pairs=mp, n_cohorts=0, mde_h=np.nan, proceed=False, note="no gated study: not estimable")); continue
        null, sd = null_sd(r0, a.n_perm, [SEED, 11])
        prow.append(dict(design=label, min_pairs=mp, n_cohorts=int(r0.pooled.n_cohorts), cohorts=r0.pooled.cohorts,
                         n_pairs=int(r0.pooled.n_pairs), null_sd_h=sd, mde_h=2.49 * sd, n_swaps=a.n_perm,
                         proceed=bool(2.49 * sd <= 1.0)))
    power = pd.DataFrame(prow)
    power.to_csv(OUT / f"power{SUF}.csv", index=False)
    say("\n" + power.round(3).to_string())
    if not power.loc[power.design == "primary", "proceed"].fillna(False).iloc[0]:
        say("STOP: MDE > 1.0 h (or no gated cohort): not estimable. Nothing further is scored.")
        return

    # ---- 3. H187.1 + companions ----------------------------------------------------------------------------------
    res = run(xg, "G2M", min_pairs=a.min_pairs, n_perm=a.n_perm, pm=pm, resp=True)
    keys, w = weights(res)
    res.cohorts.to_csv(OUT / f"per_study{SUF}.csv", index=False)
    obs = float(res.pooled.offset_h)
    comp = {k: pseudo_offsets(xg, "G2M", res, pm, a.n_pseudo, k, SEED) for k in ("matched", "responsive")}
    crow = []
    for k, v in comp.items():
        crow.append(dict(kind=k, n_sets=len(v), pseudo_median_abs=float(np.median(abs(v))), pseudo_abs_q95=float(np.quantile(abs(v), 0.95)),
                         abs_offset=abs(obs), pctile_of_abs_offset=float(np.mean(abs(v) <= abs(obs))), beyond_q95=bool(abs(obs) > np.quantile(abs(v), 0.95))))
    crow.append(dict(kind="r_pseudo_clock_paired", pseudo_q95=res.pooled.get("pseudo_q95"), pseudo_median=res.pooled.get("pseudo_median"),
                     responsive_q95=res.pooled.get("responsive_q95"), responsive_median=res.pooled.get("responsive_median"), r_clock=res.pooled.fixed))
    pd.DataFrame(crow).to_csv(OUT / f"companions{SUF}.csv", index=False)
    pd.DataFrame({"kind": np.repeat(list(comp), a.n_pseudo), "offset_h": np.concatenate(list(comp.values()))}).to_csv(OUT / f"companions_pseudo_offsets{SUF}.csv", index=False)
    t1 = dict(hypothesis="H187.1", offset_h=obs, p_two_sided=float(res.pooled.p_offset), n_swaps=a.n_perm, cohorts=res.pooled.cohorts,
              n_pairs=int(res.pooled.n_pairs), r_fixed=float(res.pooled.fixed), r_random=float(res.pooled.get("random", np.nan)),
              pseudo_pctile_matched=crow[0]["pctile_of_abs_offset"], pseudo_pctile_responsive=crow[1]["pctile_of_abs_offset"],
              beyond_q95_matched=crow[0]["beyond_q95"], beyond_q95_responsive=crow[1]["beyond_q95"], mde_h=float(power.mde_h.iloc[0]))
    t1["passes_p05"] = t1["p_two_sided"] < 0.05
    t1["clock_specific"] = bool(t1["passes_p05"] and t1["beyond_q95_matched"])
    say("H187.1", {k: v for k, v in t1.items() if k != "cohorts"})

    # ---- 4. artefacts (descriptive): (a) drop NR1D1, PER1; (c) depth-matched G2M ----------------------------------
    arows = [dict(check="primary", offset_h=obs, p=float(res.pooled.p_offset), n_cohorts=len(keys), n_pairs=int(res.pooled.n_pairs))]
    ra = run(xg, "G2M", drop=("NR1D1", "PER1"), min_pairs=a.min_pairs, n_perm=a.n_perm, seed=SEED + 1)
    arows.append(dict(check="a_drop_NR1D1_PER1", offset_h=float(ra.pooled.offset_h), p=float(ra.pooled.p_offset), n_cohorts=int(ra.pooled.n_cohorts),
                      n_pairs=int(ra.pooled.n_pairs)))
    xt = thin(xg, SEED + 2)
    rc = run(xt, "G2M", min_pairs=a.min_pairs, n_perm=a.n_perm, seed=SEED + 2)
    arows.append(dict(check="c_G2M_thinned_to_NC_depth", offset_h=float(rc.pooled.offset_h), p=float(rc.pooled.p_offset), n_cohorts=int(rc.pooled.n_cohorts),
                      n_pairs=int(rc.pooled.n_pairs)))
    pd.DataFrame(arows).to_csv(OUT / f"artefacts{SUF}.csv", index=False)
    say("\n" + pd.DataFrame(arows).round(3).to_string())

    # ---- 5. H187.2: dose, tertiles of patient cycling fraction ------------------------------------------------------
    cf = cycling_fraction()
    states, labs = {}, {}
    for k in keys:
        g = xg[xg.study == k]
        c = g[g.condition == "G2M"].set_index("patient"); r = g[g.condition == "NC"].set_index("patient")
        u = sorted(set(c.index) & set(r.index))
        states[k] = tuple(np.asarray(v) for v in P.cohort_states(c.loc[u], r.loc[u], True, P.gtex()))
        f = cf[cf.study == k].set_index("patient").frac_cycling.loc[u].to_numpy()
        rk = np.argsort(np.argsort(f, kind="stable"), kind="stable"); n = len(u)
        labs[k] = np.where(rk < n / 3, "low", np.where(rk >= n - n / 3, "high", "mid"))
    obs2 = tertile_stat(states, labs)
    rng = np.random.default_rng([SEED, 12])
    null2 = np.array([tertile_stat(states, {k: rng.permutation(v) for k, v in labs.items()}) for _ in range(a.n_perm)])
    null2 = null2[np.isfinite(null2)]
    t2 = dict(hypothesis="H187.2", statistic_abs_offset_high_minus_low_h=obs2, p_one_sided_greater=S.perm_p(null2, obs2, "greater", tie=0.0),
              n_swaps=len(null2), null_sd_h=float(null2.std(ddof=1)), mde_h=2.49 * float(null2.std(ddof=1)), cohorts=";".join(keys),
              n_patients=int(sum(len(v) for v in labs.values())), n_high=int(sum((v == "high").sum() for v in labs.values())))
    pd.DataFrame([t1, t2]).to_csv(OUT / f"tests{SUF}.csv", index=False)
    say("H187.2", t2)

    # ---- 6. descriptives -----------------------------------------------------------------------------------------
    drows = []
    xs = contrast(sam, "S")
    if a.gate == "all_nc":
        xs = xs[xs.study.isin(set(xg.study))]
    rs_ = run(xs, "S", min_pairs=a.min_pairs, n_perm=a.n_perm, seed=SEED + 3)
    drows.append(dict(quantity="S_minus_NC_pooled", offset_h=float(rs_.pooled.offset_h), p=float(rs_.pooled.p_offset), r=float(rs_.pooled.fixed),
                      n_cohorts=int(rs_.pooled.n_cohorts), n_pairs=int(rs_.pooled.n_pairs)))
    drows.append(dict(quantity="G2M_minus_NC_pooled", offset_h=obs, p=float(res.pooled.p_offset), r=float(res.pooled.fixed),
                      n_cohorts=len(keys), n_pairs=int(res.pooled.n_pairs)))
    ps_all = run(xs, "S", min_pairs=3, n_perm=100, gated_only=False)
    pg_all = run(xg, "G2M", min_pairs=3, n_perm=100, gated_only=False)
    ps = pd.concat([pg_all.cohorts.assign(contrast="G2M-NC"), ps_all.cohorts.assign(contrast="S-NC")])
    ps.to_csv(OUT / f"per_study_all_min3{SUF}.csv", index=False)
    for k in keys:  # leave one study out
        rr = run(xg[xg.study != k], "G2M", min_pairs=a.min_pairs, n_perm=a.n_perm, seed=SEED + 4)
        if int(rr.pooled.n_cohorts) == 0:
            drows.append(dict(quantity=f"LOSO_without_{k}", n_cohorts=0, note="no cohort left")); continue
        drows.append(dict(quantity=f"LOSO_without_{k}", offset_h=float(rr.pooled.offset_h), p=float(rr.pooled.p_offset), r=float(rr.pooled.fixed),
                          n_cohorts=int(rr.pooled.n_cohorts), n_pairs=int(rr.pooled.n_pairs)))
    er = pd.concat([
        pd.read_csv(C4.parent / "scouts_2026-10-31/sc_breast_er/wu_reconciled.csv").rename(columns={"orig.ident": "patient", "final_ER": "er"})[["patient", "er"]].assign(study="GSE176078"),
        pd.read_csv(C4.parent / "scouts_2026-11-01/c4_single_cell/pal_per_patient_malignant.csv")[["patient", "er"]].assign(study="GSE161529")])
    er["patient"] = er.patient.astype(str)
    for lab in ("ER+", "ER-"):
        ids = set(zip(er[er.er == lab].study, er[er.er == lab].patient))
        sub = xg[[(s, p) in ids for s, p in zip(xg.study, xg.patient)]]
        for mp, gated in ((3, False),):
            rr = run(sub, "G2M", min_pairs=mp, n_perm=a.n_perm, seed=SEED + 5, gated_only=gated)
            drows.append(dict(quantity=f"Wu+Pal_{lab}", offset_h=float(rr.pooled.offset_h) if len(rr.cohorts) else np.nan,
                              p=float(rr.pooled.get("p_offset", np.nan)), r=float(rr.pooled.get("fixed", np.nan)), n_cohorts=int(rr.pooled.get("n_cohorts", 0)),
                              n_pairs=int(rr.pooled.get("n_pairs", 0)), note="min_pairs 3, ungated"))
    pd.DataFrame(drows).to_csv(OUT / f"descriptive{SUF}.csv", index=False)
    say("\n" + pd.DataFrame(drows).round(3).to_string())
    say("done")


if __name__ == "__main__":
    main()
