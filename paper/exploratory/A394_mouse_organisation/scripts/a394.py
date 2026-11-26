"""A394 (mouse added to the amplitude-coupling rule, A383 part 1, and to the organ gain ranks, A384), registered
2026-11-26.

Deota 2023 (GSE190389; `common/datasets/mouse_ageing_deota`, genes="all"): 48 mice, 22 tissues, feeding arm = stratum.
A383 estimator: `Dataset.fit(gate=False)`, the 16 peripheral tissues as anchors read leave-one-out, brain regions as
targets only, clock = the 14 genes on the GTEx waveform. Sampling times are not used.

Stages:
  power  mouse clock fit with 200 animal bootstraps (`mouse_fit.csv`; read for standard errors and counts only); every
         companion fit (200 matched and 200 responsive mouse sets, the A384 matched sets present in mouse) and their
         statistics; analytic power of H394.1 and simulated power of H394.2. No mouse clock coupling, amplitude, gain or
         correlation is printed or written.
  check  2 sets per family, fits only, nothing written.
  score  H394.1 and H394.2 (each only if its gate passed in `power_a394.csv`) and the registered descriptives.

    RUN_CAP_GB=20 paper/common/run.sh paper/exploratory/A394_mouse_organisation/scripts/a394.py --stage power --jobs 4
    RUN_CAP_GB=20 paper/common/run.sh paper/exploratory/A394_mouse_organisation/scripts/a394.py --stage score
Outputs (`../output/`): power: mouse_fit.csv, mouse_companion_fits.csv, sets_mouse.csv, companions_a394.csv,
power_a394.csv; score: tests_a394.csv, pairs_a394.csv, descriptive_a394.csv. Tissue- and organ-level aggregates only.
"""
import argparse
import sys
import time
from multiprocessing import Pool
from pathlib import Path

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(PAPER / "exploratory" / "A383_species_organisation" / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

import a383_a384 as B  # noqa: E402
import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.paired import pseudo_waveforms  # noqa: E402
from common.datasets.mouse_ageing_deota.build import BRAIN, build, tissue_qc  # noqa: E402

OUT = HERE.parents[1] / "output"
A383_OUT = PAPER / "exploratory" / "A383_species_organisation" / "output"
SEED = 394
WF = B.A.WF
KEYS = ["donor", "compartment", "depth", "stratum"]
FIT_COLS = list(dict.fromkeys(B.FIT_COLS + ["log_amp_rel"]))  # no duplicate column
N_SETS, N_PERM, N_BOOT_B, MIN_A384_SETS = 200, 10000, 2000, 150
OTHERS = ("human", "baboon", "sheep")
# organ universe: the A384 rows plus two baboon-mouse rows; mouse tissue keyed by the baboon code (fixed in A394)
ORGANS = B.ORGANS + [(None, "PAN", None), (None, "KIM", None)]
MOUSE = {"STF": "Stomach", "OMF": "eWAT", "WAS": "iWAT", "ADC": "Adrenal", "HEA": "Heart", "KIC": "Kidney cortex",
         "LIV": "Liver", "LUN": "Lung", "MUG": "Muscle", "DUO": "Duodenum", "ILE": "Ileum", "SPL": "Spleen",
         "PAN": "Pancreas", "KIM": "Kidney medulla"}
MOUSE_KEYS = [MOUSE.get(o[1]) for o in ORGANS]
W = {}


# ---- data and fits -------------------------------------------------------------------------------------------

def mouse_table():
    f = OUT / "mouse_genomewide.parquet"
    if not f.exists():
        pb, _ = build(genes="all")
        pb = pb.loc[:, ~pb.columns.duplicated()]
        pb.to_parquet(f)
    pb = pd.read_parquet(f)
    q = tissue_qc(pb[KEYS + list(WF.genes)])
    kept = list(q.index[q.kept])
    per = [t for t in kept if t not in BRAIN]
    brain = [t for t in kept if t in BRAIN]
    return pb[pb.compartment.isin(kept)].reset_index(drop=True), per, brain


def fit(pb, wf, per, brain, n_boot, seed=SEED, name="mouse"):
    ds = cs.Dataset(pb[KEYS + list(wf.genes)], anchors=per, waveform=wf, name=name, min_donors=5)
    return ds.fit(n_boot=n_boot, gate=False, targets=brain + per, seed=seed).table


def task(job):
    family, i, genes = job
    pb, per, brain = W["pb"], W["per"], W["brain"]
    wf = B.waveform(genes, f"{family}{i}")
    try:
        out = fit(pb, wf, per, brain, 0, name=f"{family}{i}")[FIT_COLS].rename_axis("tissue").reset_index()
    except Exception as e:  # a failed set is counted, not dropped silently
        out = pd.DataFrame([dict(tissue=None, error=f"{type(e).__name__}: {e}"[:200])])
    return out.assign(family=family, set=i)


def log_mean(pb):
    genes = [c for c in pb.columns if c not in KEYS]
    x = pb[genes].to_numpy(np.float32)
    lm = pd.Series(np.log((x / pb.depth.to_numpy(float)[:, None]).mean(0) + 1e-12), index=genes)[x.mean(0) >= 1]
    return lm[~lm.index.duplicated()]


def companion_jobs(pb, n):
    lm = log_mean(pb)
    jobs = []
    for kind, tag in (("matched", 2), ("responsive", 3)):
        ws = pseudo_waveforms(lm, n, np.random.default_rng([SEED, tag]), WF, kind=kind)
        jobs += [(kind, j, list(w.genes)) for j, w in enumerate(ws)]
    s = pd.read_csv(A383_OUT / "sets_a384.csv")
    ok = [(int(i), g.split(";")) for i, g in zip(s.set, s.genes) if all(x in pb.columns for x in g.split(";"))]
    jobs += [("a384", i, g) for i, g in ok[:n]]
    return jobs


def run_fits(jobs, n_jobs, t0):
    parts = []
    with Pool(n_jobs, initializer=limit_threads, initargs=(1,)) as p:
        for k, r in enumerate(p.imap_unordered(task, jobs, chunksize=1)):
            parts.append(r)
            if (k + 1) % 50 == 0 or k + 1 == len(jobs):
                print(f"{k + 1}/{len(jobs)} ({time.time() - t0:.0f}s)", flush=True)
    return pd.concat(parts, ignore_index=True)


# ---- statistics ----------------------------------------------------------------------------------------------

def other_gains(tabs):
    """species -> log gain over the A394 organ universe for human, baboon, sheep (tabs: species -> fit table)."""
    return {sp: np.array([float(tabs[sp].log_amp_rel.get(o[j], np.nan)) if o[j] is not None else np.nan
                          for o in ORGANS]) for j, sp in enumerate(OTHERS)}


def mouse_gains(tab):
    return np.array([float(tab.log_amp_rel.get(k, np.nan)) if k is not None else np.nan for k in MOUSE_KEYS])


def mouse_rhos(m, G):
    out = {}
    for sp in OTHERS:
        k = ~np.isnan(m) & ~np.isnan(G[sp])
        out[sp] = spearmanr(m[k], G[sp][k]).statistic if k.sum() >= 4 else np.nan
    return out


def t_m(m, G):
    return float(np.mean(list(mouse_rhos(m, G).values())))


def perm_mouse(m, G, rng, n):
    """Null of T_m and of each pair: mouse gains permuted over its mapped organs."""
    k = ~np.isnan(m)
    T, P = np.empty(n), {sp: np.empty(n) for sp in OTHERS}
    for i in range(n):
        y = m.copy()
        y[k] = rng.permutation(y[k])
        r = mouse_rhos(y, G)
        T[i] = np.mean(list(r.values()))
        for sp in OTHERS:
            P[sp][i] = r[sp]
    return T, P


def companion_stats(C, per, F_others):
    """Per mouse companion set: rho (own amplitude, coupling), median peripheral coupling, b; per A384 set: T_m."""
    stored = pd.read_csv(A383_OUT / "companion_fits.csv")
    stored = {(sp, i): g.dropna(subset=["tissue"]).set_index("tissue")
              for (sp, f, i), g in stored.groupby(["species", "family", "set"]) if f == "a384"}
    rows = []
    for (fam, i), g in C.dropna(subset=["tissue"]).groupby(["family", "set"]):
        tab = g.set_index("tissue")
        r = dict(family=fam, set=i, rho=B.rho(tab, per)[0], b=B.b_of(tab, per)[0],
                 median_coupling=float(tab.coupling.reindex(per).median()))
        if fam == "a384":
            if all((sp, i) in stored for sp in OTHERS):
                G = other_gains({sp: stored[(sp, i)] for sp in OTHERS})
                rr = mouse_rhos(mouse_gains(tab), G)
                r.update(T_m=float(np.mean(list(rr.values()))), **{f"rho_mouse_{sp}": rr[sp] for sp in OTHERS})
        rows.append(r)
    return pd.DataFrame(rows)


def power_h394_2(se_m, F_others, sd_true, rng, n_sim=2000, n_perm=500, rho=0.6):
    se = {sp: np.array([((F_others[sp].log_amp_hi[o[j]] - F_others[sp].log_amp_lo[o[j]]) / 3.92)
                        if o[j] is not None and o[j] in F_others[sp].index else np.nan for o in ORGANS])
          for j, sp in enumerate(OTHERS)}
    present = {sp: ~np.isnan(se[sp]) for sp in OTHERS}
    present["mouse"], se["mouse"] = ~np.isnan(se_m), se_m
    hit = 0
    for _ in range(n_sim):
        f = rng.normal(0, 1, len(ORGANS))
        v = {}
        for sp in ("mouse",) + OTHERS:
            x = sd_true * (np.sqrt(rho) * f + np.sqrt(1 - rho) * rng.normal(0, 1, len(ORGANS)))
            v[sp] = np.where(present[sp], x + rng.normal(0, np.nan_to_num(se[sp])), np.nan)
        G = {sp: v[sp] for sp in OTHERS}
        T = t_m(v["mouse"], G)
        nT, _ = perm_mouse(v["mouse"], G, rng, n_perm)
        hit += (1 + np.sum(nT >= T)) / (n_perm + 1) < 0.05
    return hit / n_sim


def others():
    F = B.fits()
    return {sp: F[sp][0] for sp in OTHERS}


# ---- stages --------------------------------------------------------------------------------------------------

def power(n_jobs):
    limit_threads(4)
    t0 = time.time()
    pb, per, brain = mouse_table()
    W.update(pb=pb, per=per, brain=brain)
    f = OUT / "mouse_fit.csv"
    if not f.exists():
        fit(pb, WF, per, brain, 200).rename_axis("tissue").to_csv(f)
    tab = pd.read_csv(f, index_col="tissue")  # SEs and counts only
    n1 = int(tab.reindex(per)[["own_amplitude", "coupling"]].notna().all(axis=1).sum())
    se_m = np.array([((tab.log_amp_hi[k] - tab.log_amp_lo[k]) / 3.92) if k is not None and k in tab.index else np.nan
                     for k in MOUSE_KEYS])
    print(f"{pb.donor.nunique()} mice, {len(per)} peripheral, {len(brain)} brain; clock fit ({time.time() - t0:.0f}s)",
          flush=True)
    jobs = companion_jobs(pb, N_SETS)
    pd.DataFrame([dict(family=fm, set=i, genes=";".join(g)) for fm, i, g in jobs]).to_csv(OUT / "sets_mouse.csv",
                                                                                         index=False)
    C = run_fits(jobs, n_jobs, t0)
    C.to_csv(OUT / "mouse_companion_fits.csv", index=False)
    F_others = others()
    S = companion_stats(C, per, F_others)
    S.to_csv(OUT / "companions_a394.csv", index=False)
    q95 = {k: float(np.nanquantile(S[S.family == k].rho, 0.95)) for k in ("matched", "responsive")}
    thr = max(q95.values())
    a384 = S[S.family == "a384"]
    n_a384 = int(a384.T_m.notna().sum()) if "T_m" in a384 else 0
    h = F_others["human"]
    hk = [o[0] for o in ORGANS if o[0] is not None]
    sd_true = float(np.std(h.log_amp_rel.reindex(hk).to_numpy(float), ddof=1))
    p2 = power_h394_2(se_m, F_others, sd_true, np.random.default_rng([SEED, 5]))
    n_pair = {sp: int(sum(1 for j, o in enumerate(ORGANS) if MOUSE_KEYS[j] is not None and MOUSE_KEYS[j] in tab.index
                          and o[OTHERS.index(sp)] is not None and o[OTHERS.index(sp)] in F_others[sp].index))
              for sp in OTHERS}
    P = dict(mice=int(pb.donor.nunique()), peripheral=len(per), brain=len(brain), tissues_h394_1=n1,
             rho_matched_median=float(S[S.family == "matched"].rho.median()), rho_matched_q95=q95["matched"],
             rho_responsive_median=float(S[S.family == "responsive"].rho.median()), rho_responsive_q95=q95["responsive"],
             companion_median_coupling_matched=float(S[S.family == "matched"].median_coupling.median()),
             companion_median_coupling_responsive=float(S[S.family == "responsive"].median_coupling.median()),
             power_h394_1_rho06=B.spearman_power(n1, 0.6, thr), threshold_h394_1=thr,
             **{f"organs_mouse_{sp}": n_pair[sp] for sp in OTHERS},
             median_se_mouse=float(np.nanmedian(se_m)), sd_true_gtex=sd_true, power_h394_2_rho06=p2,
             a384_sets_usable=n_a384,
             T_m_matched_median=float(a384.T_m.median()) if n_a384 else np.nan,
             T_m_matched_q95=float(a384.T_m.quantile(0.95)) if n_a384 else np.nan,
             failed_fits=int(C.tissue.isna().sum()))
    P["gate_h394_1"] = bool(P["power_h394_1_rho06"] >= 0.8)
    P["gate_h394_2"] = bool(p2 >= 0.8 and n_a384 >= MIN_A384_SETS)
    pd.Series(P).to_csv(OUT / "power_a394.csv", header=["value"])
    print(pd.Series(P).to_string(), f"\n({time.time() - t0:.0f}s)")


def check(n_jobs):
    t0 = time.time()
    pb, per, brain = mouse_table()
    W.update(pb=pb, per=per, brain=brain)
    jobs = companion_jobs(pb, 2)
    C = run_fits(jobs, n_jobs, t0)
    print(C.groupby("family").agg(fits=("set", "nunique"), tissues=("tissue", "count")).to_string(),
          f"failed: {int(C.tissue.isna().sum())}; peripheral {len(per)}; genes {pb.shape[1] - 4}; "
          f"({time.time() - t0:.0f}s)", sep="\n")


def score():
    limit_threads(4)
    t0 = time.time()
    P = pd.read_csv(OUT / "power_a394.csv", index_col=0).value
    gate1, gate2 = P["gate_h394_1"] == "True", P["gate_h394_2"] == "True"
    pb, per, brain = mouse_table()
    tab = pd.read_csv(OUT / "mouse_fit.csv", index_col="tissue")
    S = pd.read_csv(OUT / "companions_a394.csv")
    tests, pairs, desc = [], [], []
    # H394.1
    r, t = B.rho(tab, per)
    rng = np.random.default_rng([SEED, 1])
    null = np.array([spearmanr(t.own_amplitude, rng.permutation(t.coupling.to_numpy())).statistic
                     for _ in range(N_PERM)])
    row = dict(test="H394.1", scored=bool(gate1), tissues=len(t), stat=r,
               p_perm=float((1 + np.sum(null >= r)) / (1 + N_PERM)))
    for k in ("matched", "responsive"):
        v = S[S.family == k].rho.to_numpy(float)
        row[f"{k}_q95"], row[f"{k}_pct_below"] = float(np.nanquantile(v, 0.95)), float(100 * np.mean(v[np.isfinite(v)] < r))
    row["passes"] = bool(gate1 and r > 0 and row["p_perm"] < 0.05 and r > row["matched_q95"] and r > row["responsive_q95"])
    tests.append(row)
    # H394.2: not computed when its gate failed (registered: not scored)
    G = other_gains(others())
    m = mouse_gains(tab)
    if not gate2:
        tests.append(dict(test="H394.2", scored=False, tissues=int((~np.isnan(m)).sum())))
    rr = mouse_rhos(m, G) if gate2 else {}
    if gate2:
        T = float(np.mean(list(rr.values())))
        nT, nP = perm_mouse(m, G, np.random.default_rng([SEED, 4]), N_PERM)
        a = S[S.family == "a384"]
        row = dict(test="H394.2", scored=True, tissues=int((~np.isnan(m)).sum()), stat=T,
                   p_perm=float((1 + np.sum(nT >= T)) / (1 + N_PERM)), matched_q95=float(a.T_m.quantile(0.95)),
                   matched_pct_below=float(100 * np.mean(a.T_m.dropna() < T)), matched_valid=int(a.T_m.notna().sum()))
        row["passes"] = bool(T > 0 and row["p_perm"] < 0.05 and T > row["matched_q95"])
        tests.append(row)
    for sp in (OTHERS if gate2 else ()):
        k = ~np.isnan(m) & ~np.isnan(G[sp])
        pairs.append(dict(pair=f"mouse-{sp}", organs=int(k.sum()), rho=rr[sp],
                          p_perm=float((1 + np.sum(nP[sp] >= rr[sp])) / (1 + N_PERM)),
                          matched_q95=float(a[f"rho_mouse_{sp}"].quantile(0.95))))
    # descriptives
    for arm in sorted(pb.stratum.unique()):
        ta = fit(pb[pb.stratum == arm], WF, per, brain, 0, name=f"mouse_{arm}")
        desc.append(dict(item=f"rho H394.1, arm {arm}", value=B.rho(ta, per)[0],
                         n=int(pb[pb.stratum == arm].donor.nunique())))
    b, x, y = B.b_of(tab, per)
    bs = B.boot_slopes(x, y, np.random.default_rng([SEED, 6]), N_BOOT_B)
    desc += [dict(item="b (part 2)", value=b, lo90=float(np.quantile(bs, 0.05)), hi90=float(np.quantile(bs, 0.95)),
                  n=len(x))]
    for k in ("matched", "responsive"):
        v = S[S.family == k].b
        desc.append(dict(item=f"b {k} median", value=float(v.median()), lo90=float(v.quantile(0.05)),
                         hi90=float(v.quantile(0.95))))
    C = pd.read_csv(OUT / "mouse_companion_fits.csv")
    mt = C[C.family == "matched"].dropna(subset=["tissue"])
    med2psi = (2 * mt.psi).groupby(mt.tissue).median().reindex(per)
    lr = np.log((2 * tab.psi.reindex(per)) / med2psi)
    desc += [dict(item="own-term log ratio clock / matched, median", value=float(lr.median()),
                  n=int(lr.notna().sum()), lo90=float((lr > 0).sum())),
             dict(item="median peripheral coupling, clock", value=float(tab.coupling.reindex(per).median()))]
    for k in ("matched", "responsive"):
        desc.append(dict(item=f"median peripheral coupling, {k} (median over sets)",
                         value=float(S[S.family == k].median_coupling.median())))
    gt = pd.DataFrame({"human": [o[0] for o in ORGANS], "baboon": [o[1] for o in ORGANS],
                       "sheep": [o[2] for o in ORGANS], "mouse": MOUSE_KEYS, "log_gain_mouse": m,
                       **{f"log_gain_{sp}": G[sp] for sp in OTHERS}})
    T_, Pp, D = pd.DataFrame(tests), pd.DataFrame(pairs), pd.DataFrame(desc)
    T_.to_csv(OUT / "tests_a394.csv", index=False)
    Pp.to_csv(OUT / "pairs_a394.csv", index=False)
    D.to_csv(OUT / "descriptive_a394.csv", index=False)
    gt.to_csv(OUT / "gains_a394.csv", index=False)
    pd.set_option("display.width", 250)
    print(T_.round(4).to_string(index=False), Pp.round(4).to_string(index=False), D.round(4).to_string(index=False),
          gt.dropna(subset=["mouse"]).round(2).to_string(index=False), f"({time.time() - t0:.0f}s)", sep="\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=["power", "check", "score"], required=True)
    ap.add_argument("--jobs", type=int, default=4)
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    {"power": lambda: power(a.jobs), "check": lambda: check(a.jobs), "score": score}[a.stage]()


if __name__ == "__main__":
    main()
