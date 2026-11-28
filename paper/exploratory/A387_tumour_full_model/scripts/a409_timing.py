"""A409 addendum 1 (POST HOC, DESCRIPTIVE; registered in paper/docs/registered_experiments.md): is the A409 tumour-normal
stroma loudness leverage what shared biopsy time on a hospital-hours arc produces (A366 / A368 for MuTHER)?

  G1 elongation of the pooled within-stratum clouds (A366 `axes`, `combined`) and per stratum (>= 6 patients)
  G2 split-half reliability of each state (10 gene-half splits; host refitted per half)
  G3 stratified timing-only ring simulation through the A409 estimator; GTEx-size and half-size loudness alternatives
  G4 A368 Q on the pooled within-stratum states, its simulated distributions and pseudo-waveform companions
                                                     -> output/a409/timing_{geometry,reliability,sim,tests,companions}.csv
    python a409_timing.py [--jobs 8] [--genome-jobs 2]                  (paper/common/run.sh)
"""
from __future__ import annotations

import argparse
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "A366_muther_timing" / "scripts"))
sys.path.insert(0, str(HERE.parents[1] / "A368_muther_radial" / "scripts"))
import a409 as A  # noqa: E402
import a387_score as SC  # noqa: E402
import a366_muther_timing as M6  # noqa: E402
import a368_muther_radial as M8  # noqa: E402
from a387 import cs, limit_threads  # noqa: E402
from clocksync import deviation, pairwise  # noqa: E402
from clocksync.elasticity import radius_covariance  # noqa: E402

OUT = A.OUT
N_SIM, SIGMA_A, N_SET = 1000, 0.352, 200
PAIRS = {"tumour fibroblast ~ normal-stroma host": (A.TF, A.NORMAL_ANCH),
         "ref: normal fibroblast ~ normal endothelium + mural host": (A.NF, [A.NE, A.NM])}
T0 = time.time()


def log(*a):
    print(f"[{(time.time() - T0) / 60:6.1f} min]", *a, flush=True)


def reliability(table, target, anchors, n_split=10):
    """Split-half reliabilities (target state, host state) over gene-half splits."""
    limit_threads(1)
    warnings.filterwarnings("ignore")
    ra, rb = [], []
    for ga, gb in deviation.split_genes(A.BASE, n_split, np.random.default_rng([4091])):
        halves = []
        for g in (ga, gb):
            ids, a, b, st = A.pair_states(table, target, anchors, A.BASE.subset(list(g)))
            halves.append(pd.DataFrame({"a": a, "b": b}, index=ids))
        j = halves[0].index.intersection(halves[1].index)
        for col, acc in (("a", ra), ("b", rb)):
            x, y = halves[0].loc[j, col].to_numpy(), halves[1].loc[j, col].to_numpy()
            rho = abs(np.sum(x * np.conj(y))) / np.sqrt(np.sum(abs(x) ** 2) * np.sum(abs(y) ** 2))
            acc.append(2 * rho / (1 + rho))
    return float(np.median(ra)), float(np.median(rb))


def sim(st, sd, sigma, fa, fb, rng):
    """Stratified ring: stratum mean phase uniform, patient phases N(mean, sd), log amplitude SD sigma; noise at the
    clock fractions fa / fb; states centred within stratum (as the A409 states)."""
    n = len(st)
    mu = {s: rng.uniform(0, 2 * np.pi) for s in np.unique(st)}
    ph = np.array([mu[s] for s in st]) + rng.normal(0, sd, n)
    T = np.exp(rng.normal(0, sigma, n) - sigma ** 2 / 2) * np.exp(1j * ph)
    T = A.centred_by(T, st)
    v = np.mean(np.abs(T) ** 2)
    out = []
    for f in (fa, fb):
        s = np.sqrt(v * (1 - f) / f / 2)
        out.append(A.centred_by(T + s * (rng.normal(size=n) + 1j * rng.normal(size=n)), st))
    return out


def lev(a, b, st):
    return radius_covariance(a, b, groups=st)["leverage"]


def comp_q(table, wfs):
    out = []
    for fam, j, w in wfs:
        try:
            _, a, b, st = A.pair_states(table, A.TF, A.NORMAL_ANCH, w)
            out.append(dict(family=fam, set=j, Q=M8.Q(a, b)))
        except Exception:  # noqa: BLE001
            out.append(dict(family=fam, set=j, Q=np.nan))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--genome-jobs", type=int, default=2)
    args = ap.parse_args()
    limit_threads(1)
    table = A.R.pooled_table(args.jobs, False)
    geo, rel, sims, tests = [], [], [], []
    for lab, (tgt, anch) in PAIRS.items():
        ids, a, b, st = A.pair_states(table, tgt, anch, A.BASE)
        assert abs(lev(a, b, st) - {"tumour": 0.3612, "ref: n": 0.4853}[lab[:6]]) < 5e-4, "A409 leverage not reproduced"
        c = M6.combined(a, b)
        geo.append(dict(pair=lab, cloud="target", elongation=M6.elongation(a)))
        geo.append(dict(pair=lab, cloud="host", elongation=M6.elongation(b)))
        geo.append(dict(pair=lab, cloud="combined", elongation=M6.elongation(c)))
        for s in pd.unique(st):
            m = st == s
            if m.sum() >= 6:
                geo.append(dict(pair=lab, cloud=f"combined, stratum {s}", n=int(m.sum()), elongation=M6.elongation(c[m])))
        fa, fb = reliability(table, tgt, anch)
        rel.append(dict(pair=lab, rel_target=fa, rel_host=fb))
        target_el = M6.elongation(c)
        rng = np.random.default_rng([4092, len(lab)])
        grid = np.array([np.median([M6.elongation(M6.combined(*sim(st, sd, 0.0, fa, fb, rng))) for _ in range(20)]) for sd in M6.SD_GRID])
        sd_hat = float(M6.SD_GRID[np.argmin(np.abs(grid - target_el))])
        cols = {}
        for name, sg in (("timing_only", 0.0), ("half_gtex", SIGMA_A / 2), ("gtex", SIGMA_A)):
            L, Qv = [], []
            for _ in range(N_SIM):
                x, y = sim(st, sd_hat, sg, fa, fb, rng)
                L.append(lev(x, y, st))
                Qv.append(M8.Q(x, y))
            cols[f"lev_{name}"], cols[f"Q_{name}"] = L, Qv
        S_ = pd.DataFrame(cols).assign(pair=lab)
        sims.append(S_)
        lo, qo = lev(a, b, st), M8.Q(a, b)
        row = dict(pair=lab, patients=len(ids), strata=len(set(st)), elongation_combined=target_el, phase_sd_h=sd_hat * 24 / (2 * np.pi),
                   rel_target=fa, rel_host=fb, leverage=lo, Q=qo)
        for k in ("timing_only", "half_gtex", "gtex"):
            row.update({f"lev_{k}_median": float(np.median(S_[f"lev_{k}"])), f"lev_{k}_q95": float(np.quantile(S_[f"lev_{k}"], .95)),
                        f"lev_pct_{k}": float(100 * np.mean(S_[f"lev_{k}"] < lo)),
                        f"Q_{k}_median": float(np.median(S_[f"Q_{k}"])), f"Q_{k}_q95": float(np.quantile(S_[f"Q_{k}"], .95)),
                        f"Q_pct_{k}": float(100 * np.mean(S_[f"Q_{k}"] < qo))})
        tests.append(row)
        log(lab, {k: round(v, 3) for k, v in row.items() if isinstance(v, float)})
    pd.DataFrame(geo).to_csv(OUT / "timing_geometry.csv", index=False)
    pd.DataFrame(rel).to_csv(OUT / "timing_reliability.csv", index=False)
    pd.concat(sims).to_csv(OUT / "timing_sim.csv", index=False)
    # G4 companions of Q (tumour fibroblast pair), A409 pseudo-waveform construction
    keep = set(zip(table.donor, table.compartment))
    studies = sorted({d.split("|", 1)[0] for d in table.donor})
    Gs = Parallel(n_jobs=args.genome_jobs)(delayed(SC.genome_cohort)(s, {x for x in keep if x[0].startswith(s + "|")}) for s in studies)
    G = pd.concat([g for g in Gs if g is not None and len(g)], ignore_index=True)
    G = G.drop(columns=["stratum"]).merge(table[["donor", "compartment", "stratum"]], on=["donor", "compartment"])
    genes = [c for c in G.columns if c not in ("donor", "compartment", "stratum", "depth", "n_cells")]
    genes = [g for g in genes if G[g].notna().all()]
    ref = G[G.compartment.isin(A.NORMAL_ANCH)]
    frac = ref[genes].div(ref.depth, axis=0)
    lm = np.log(frac.mean() + 1e-12)[(ref[genes].mean() >= 1) | frac.columns.isin(list(A.BASE.genes))]
    from clocksync.paired import pseudo_waveforms
    wfs = [(f, j, w) for code, f in enumerate(("matched", "responsive"))
           for j, w in enumerate(pseudo_waveforms(lm, N_SET, np.random.default_rng([409, code]), A.BASE, f))]
    cols = [c for c in ("donor", "compartment", "stratum", "depth", "n_cells") if c in G.columns]
    chunks = [wfs[i::args.jobs] for i in range(args.jobs)]
    res = Parallel(n_jobs=args.jobs)(delayed(comp_q)(G[cols + sorted({g for _, _, w in ch for g in w.genes} | set(A.BASE.genes))], ch)
                                     for ch in chunks)
    C = pd.DataFrame([r for x in res for r in x])
    C.to_csv(OUT / "timing_companions.csv", index=False)
    T = pd.DataFrame(tests)
    for fam in ("matched", "responsive"):
        v = C[C.family == fam].Q.dropna()
        T.loc[0, f"Q_{fam}_median"], T.loc[0, f"Q_{fam}_q95"] = float(v.median()), float(v.quantile(.95))
        T.loc[0, f"Q_pct_{fam}"] = float(100 * np.mean(v < T.loc[0, "Q"]))
    T.to_csv(OUT / "timing_tests.csv", index=False)
    pd.set_option("display.width", 250, "display.max_columns", 60)
    log("geometry\n" + pd.DataFrame(geo).round(3).to_string(index=False))
    log("tests\n" + T.round(4).T.to_string())
    log("done")


if __name__ == "__main__":
    main()
