"""A265 (registered in `paper/docs/registered_experiments.md`): COFE against clocksync on the TCGA tumour-normal pairs.

Parts
- `real`: (i) COFE as published on each cancer's tumours (`lib_a265`, `common.benchmark.cofe.run_cofe`), the normals given
  time labels by COFE's rule for samples not used in training (projection on the trained loadings); per cancer: Rayleigh
  test of the tumour-minus-normal label differences, largest empty arc of the tumour labels, CV minimum interior.
  (ii) pairing of COFE's labels: R_c = |mean e^{i (phi_T - phi_N)}| per cancer, pooled (weights n pairs), against 10,000
  permutations of the normals within cancer. clocksync (A138 statistic, `q3_lib.corr_stat`, normals gated by
  `reference_gate`, pooled over gated cancers) on the same pairs, same permutation scheme.
- `sim`: (iii) on the real design: the normals' patient labels permuted within cancer (real pairing broken), every
  patient given a time t (uniform in 08-18 h, or 0-24 h); normal phase t; tumour phase t + N(0, 1 h) (coupled) or uniform
  (uncoupled); counts multiplied by exp(A_g cos(w (phase - phi_g))) for the program genes (`program()`). Both methods as
  in `real`, 1,000 permutations per replicate. Arms: coupled_surgical, uncoupled_surgical, coupled_wholeday.

    python score_a265.py --part real [--jobs 12]
    python score_a265.py --part sim [--reps 10] [--jobs 12]
Outputs (`../output/`): real_per_cancer.csv, real_tests.csv; sim_reps.csv (one row per arm x replicate x cancer, and a
pooled row), sim_summary.csv.
"""
import argparse
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[3]))
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parents[2] / "A138_cancer_estimability" / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy import stats  # noqa: E402

import lib_a265 as L  # noqa: E402
import q3_lib as q  # noqa: E402
from clocksync.paired import cohort_states, reference_gate  # noqa: E402
from common.benchmark.cofe import run_cofe  # noqa: E402

OUT = HERE.parents[1] / "output"
SEED, W = 265, 2 * np.pi / 24
CLOCK = list(q.WF.genes)
COSINOR = HERE.parents[3] / "results" / "figure_2" / "output" / "gene_rhythms" / "gtex_cosinor.csv.gz"
ARMS = {"coupled_surgical": (True, (8.0, 18.0)), "uncoupled_surgical": (False, (8.0, 18.0)),
        "coupled_wholeday": (True, (0.0, 24.0))}


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def program() -> pd.DataFrame:
    """Rhythm program of the simulation, fixed from GTEx (A27 timed cosinor of log CPM on the hour of death, 13 tissues):
    every gene with BH q < 0.05 in >= 1 tissue, plus the 14 clock genes; amplitude = median over its significant tissues
    (all 13 for a clock gene with none), acrophase = circular mean over the same tissues."""
    t = pd.read_csv(COSINOR)
    sig = t[t.q < 0.05]
    rows = []
    for g in sorted(set(sig.gene) | set(CLOCK)):
        x = sig[sig.gene == g]
        x = x if len(x) else t[t.gene == g]
        rows.append(dict(gene=g, amp=float(x.amp.median()),
                         acro_h=float(np.mod(np.angle(np.exp(1j * W * x.acro_h).mean()) / W, 24)), tissues=len(x)))
    return pd.DataFrame(rows).set_index("gene")


def rayleigh(d):
    n, R = len(d), abs(np.exp(1j * d).mean())
    z = n * R ** 2
    return R, float(np.exp(-z) * (1 + (2 * z - z ** 2) / (4 * n) - (24 * z - 132 * z ** 2 + 76 * z ** 3 - 9 * z ** 4) / (288 * n ** 2)))


def largest_gap(h):
    s = np.sort(np.mod(h, 24))
    return float(np.max(np.diff(np.r_[s, s[0] + 24])))


def cofe_one(lt, ln, pairs, seed, jobs):
    r = run_cofe(lt, seed=seed, ncores=jobs, test=ln.loc[pairs])
    pt, pn = r["phase"].loc[pairs].to_numpy() * W, r["test_phase"].loc[pairs].to_numpy() * W
    return pt, pn, r


def clock_frame(counts, depth, pairs):
    f = counts.loc[pairs, L.B.T.clock_columns(counts.columns)].copy()
    f.columns = CLOCK
    f["depth"] = depth.loc[pairs].to_numpy()
    return f


def clocksync_one(ct, dt, cn, dn, pairs, rng):
    case, ref = clock_frame(ct, dt, pairs), clock_frame(cn, dn, pairs)
    yr = np.log(np.maximum(ref[CLOCK].to_numpy(float), 1e-300)) - np.log(ref.depth.to_numpy(float))[:, None]
    g = reference_gate(yr, q.WF, 200, rng)
    bt, vt, bn, vn = cohort_states(case, ref, True, q.WF)
    za, zb = q.centred_complex(bt), q.centred_complex(bn)
    return dict(gate_p=g["p"], gated=bool(g["p"] < 0.05), za=za, vA=vt, zb=zb, vB=vn)


def pooled_tests(per, n_perm, rng):
    """Pooled COFE R_c and clocksync statistic with their within-cancer permutation nulls (same permutations)."""
    keys = list(per)
    w = np.array([len(per[k]["pt"]) for k in keys], float)
    obs_c, null_c, obs_s, null_s, ws = [], [], [], [], []
    for k in keys:
        d = per[k]
        n = len(d["pt"])
        perms = np.argsort(rng.random((n_perm, n)), axis=1)
        obs_c.append(abs(np.exp(1j * (d["pt"] - d["pn"])).mean()))
        null_c.append(np.abs(np.exp(1j * (d["pt"][None, :] - d["pn"][perms])).mean(1)))
        s = d["cs"]
        if s["gated"]:
            ea, eb = q.tr_corrected(s["za"], s["vA"]), q.tr_corrected(s["zb"], s["vB"])
            if ea > 0 and eb > 0:
                obs_s.append(q.corr_stat(s["za"], s["vA"], s["zb"], s["vB"]))
                null_s.append(np.abs((s["za"][None, :] * np.conj(s["zb"][perms])).sum(1) / (n - 1)) / np.sqrt(ea * eb))
                ws.append(n)
    pc, nc = float(np.average(obs_c, weights=w)), np.average(np.array(null_c), axis=0, weights=w)
    out = dict(cofe_R=pc, cofe_null_mean=float(nc.mean()), cofe_p=float((1 + np.sum(nc >= pc)) / (1 + n_perm)),
               cofe_uniform_cancers=int(sum(per[k]["rayleigh_p"] >= 0.05 for k in keys)), cancers=len(keys))
    if ws:
        ps, ns = float(np.average(obs_s, weights=ws)), np.average(np.array(null_s), axis=0, weights=ws)
        out |= dict(clocksync_stat=ps, clocksync_null_mean=float(ns.mean()),
                    clocksync_p=float((1 + np.sum(ns >= ps)) / (1 + n_perm)), clocksync_cancers=len(ws))
    else:
        out |= dict(clocksync_stat=np.nan, clocksync_null_mean=np.nan, clocksync_p=np.nan, clocksync_cancers=0)
    return out


def analyse(data, seed, jobs, n_perm, rng):
    per, rows = {}, []
    for i, (c, d) in enumerate(data.items()):
        lt, ln = L.cofe_matrices(d["counts_t"], d["counts_n"])
        pt, pn, r = cofe_one(lt, ln, d["pairs"], seed + i, jobs)
        R, p = rayleigh(pt - pn)
        cs_ = clocksync_one(d["counts_t"], d["depth_t"], d["counts_n"], d["depth_n"], d["pairs"], rng)
        per[c] = dict(pt=pt, pn=pn, rayleigh_p=p, cs=cs_)
        rows.append(dict(cancer=c, tumours=len(lt), pairs=len(d["pairs"]), genes=lt.shape[1], best_s=r["best_s"],
                         interior_min=r["interior_min"], tumour_largest_gap_h=largest_gap(r["phase"].to_numpy()),
                         diff_R=R, rayleigh_p=p, cofe_uniform=bool(p >= 0.05), clock_gate_p=cs_["gate_p"],
                         clocksync_stat=q.corr_stat(cs_["za"], cs_["vA"], cs_["zb"], cs_["vB"])))
    return pd.DataFrame(rows), pooled_tests(per, n_perm, rng)


def simulate(d, prog, coupled, window, rng):
    """One simulated cancer on the real design (see the module docstring)."""
    pairs = list(d["pairs"])
    cn = d["counts_n"].loc[pairs].copy()
    cn.index = rng.permutation(pairs)  # real pairing broken
    dn = d["depth_n"].loc[cn.index].set_axis(cn.index)
    cn, dn = cn.loc[pairs], dn.loc[pairs]
    ct, dt = d["counts_t"].copy(), d["depth_t"].copy()
    t = pd.Series(rng.uniform(*window, len(ct)), index=ct.index)
    th_n = t.loc[pairs].to_numpy()
    th_t = t.to_numpy() + rng.normal(0, 1, len(t)) if coupled else rng.uniform(0, 24, len(t))
    g = [x for x in prog.index if x in ct.columns]
    A, P = prog.loc[g, "amp"].to_numpy(), prog.loc[g, "acro_h"].to_numpy()
    ft = np.exp(A[None, :] * np.cos(W * (th_t[:, None] - P[None, :])))
    fn = np.exp(A[None, :] * np.cos(W * (th_n[:, None] - P[None, :])))
    ct.loc[:, g] = np.round(ct[g].to_numpy(float) * ft)
    cn.loc[:, g] = np.round(cn[g].to_numpy(float) * fn)
    dt = dt + (ct.sum(1) - d["counts_t"].sum(1))
    dn = dn + (cn.sum(1) - d["counts_n"].loc[cn.index].sum(1).to_numpy())
    return dict(counts_t=ct, counts_n=cn, depth_t=dt, depth_n=dn, pairs=pairs)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--part", choices=["real", "sim"], required=True)
    ap.add_argument("--jobs", type=int, default=12)
    ap.add_argument("--reps", type=int, default=10)
    ap.add_argument("--cancers", nargs="*", default=L.CANCERS)
    a = ap.parse_args()
    OUT.mkdir(exist_ok=True)
    data = {}
    for c in a.cancers:
        data[c] = L.cancer(c)
        log("loaded", c)
    if a.part == "real":
        per, t = analyse(data, SEED, a.jobs, 10_000, np.random.default_rng([SEED, 0]))
        per.to_csv(OUT / "real_per_cancer.csv", index=False)
        pd.DataFrame([t]).to_csv(OUT / "real_tests.csv", index=False)
        pd.set_option("display.width", 250)
        print(per.round(4).to_string(index=False), pd.Series(t).to_string(), sep="\n", flush=True)
        return
    prog = program()
    prog.to_csv(OUT / "sim_program.csv")
    log("program genes", len(prog))
    rows = []
    for rep in range(a.reps):
        for j, (arm, (coupled, window)) in enumerate(ARMS.items()):
            rng = np.random.default_rng([SEED, 1, rep, j])
            sim = {c: simulate(d, prog, coupled, window, rng) for c, d in data.items()}
            per, t = analyse(sim, SEED + 1000 * rep + 100 * j, a.jobs, 1000, rng)
            rows += [dict(arm=arm, rep=rep, **r) for r in per.to_dict("records")]
            rows.append(dict(arm=arm, rep=rep, cancer="pooled", **t))
            pd.DataFrame(rows).to_csv(OUT / "sim_reps.csv", index=False)
            log(arm, rep, {k: round(v, 4) if isinstance(v, float) else v for k, v in t.items()})
    s = pd.DataFrame(rows)
    s = s[s.cancer == "pooled"]
    summ = s.groupby("arm").agg(reps=("rep", "size"), clocksync_pass=("clocksync_p", lambda p: float((p < 0.05).mean())),
                                cofe_pairing_pass=("cofe_p", lambda p: float((p < 0.05).mean())),
                                cofe_published_uncoupled=("cofe_uniform_cancers", lambda u: float((u >= 6).mean())),
                                cofe_uniform_cancers_median=("cofe_uniform_cancers", "median"),
                                clocksync_stat_median=("clocksync_stat", "median"), cofe_R_median=("cofe_R", "median"))
    summ.to_csv(OUT / "sim_summary.csv")
    print(summ.round(3).to_string(), flush=True)


if __name__ == "__main__":
    main()
