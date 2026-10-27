"""A138 power step: null distributions of H138.4 and H138.5 and simulated power of H138.1 (no observed statistic).

Registered in paper/docs/registered_experiments.md (A138, "note 3", power step). NOTHING here computes the observed
statistic of H138.1, H138.4 or H138.5: the script draws the permutation nulls on the real states (the pairing is
broken, the unpermuted pairing is never evaluated) and simulates H138.1 at the design of each Q1 cohort.
- H138.4: per cancer and pooled over the gated cancers (weights n_pairs): null mean, SD, q95, smallest p, MDE.
- H138.5: per cohort (>= 6 paired patients with a stromal host) and pooled (all; gate-passing), plus the fibroblast
  companion; same summary.
- H138.1: `clocksync.simulate` at the design (patients, malignant and anchor depth) of each Q1 cohort (A17 gate passed,
  learned anchors >= 3, >= 5 malignant / tumour-epithelium patients), true coupling in {0, .1, .2, .3, .5}, 100
  replicates, 200 bootstraps; pooled lower bound (DerSimonian-Laird over cohorts) > 0.
MDE: null q95 + 0.84 SD (normal approximation; one-sided 0.05, power 0.8); H138.1: smallest c with the lower bound
> 0 in >= 80 % of replicates ("underpowered at this design" if none).

    python power.py [--subsample N] [--jobs J] [--seed S] [--only h138_1 h138_4 h138_5]
Outputs (paper/exploratory/A138_cancer_estimability/output/power/, `subsample/` for checks): h138_4_null.csv,
h138_5_null.csv, h138_1_sim.csv, h138_1_sim_cohorts.csv, h138_1_designs.csv, power_summary.csv
"""

from __future__ import annotations

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

sys.path.insert(0, str(Path(__file__).resolve().parent))
import q3_lib as q  # noqa: E402
from common import cli  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "output" / "power"
SEED = 138
T0 = time.time()


def log(*a):
    print(f"[{(time.time() - T0) / 60:7.1f} min]", *a, flush=True)


def null_rows(scope: str, kind: str, n_units: int, null, **extra) -> dict:
    return dict(scope=scope, kind=kind, n_units=n_units, **extra, **q.summarise(null))


# ---- H138.4 ------------------------------------------------------------------------------------------------------

def run_h138_4(args, n_perm):
    st = q.tcga_states(args.subsample, args.seed)
    rows, nulls = [], {}
    for i, (c, u) in enumerate(st.items()):
        nulls[c] = q.perm_null(u["za"], u["vA"], u["zb"], u["vB"], n_perm, np.random.default_rng([SEED, 4, args.seed, i]))
        rows.append(null_rows(c, "cancer", u["n"], nulls[c], gate_p=u["gate_p"], gated=u["gated"]))
    defined = [c for c in st if np.isfinite(nulls[c]).any()]  # a cancer with a non-positive corrected variance is left out
    for scope, keep in (("pooled_gated", [c for c in defined if st[c]["gated"]]), ("pooled_all", defined)):
        if keep:
            pooled = q.pool(np.array([nulls[c] for c in keep]), [st[c]["n"] for c in keep])
            rows.append(null_rows(scope, "pooled", sum(st[c]["n"] for c in keep), pooled, n_cohorts=len(keep),
                                  cohorts=";".join(keep)))
    d = pd.DataFrame(rows)
    d.to_csv(OUT / "h138_4_null.csv", index=False)
    log("H138.4 null done:", d[d.kind == "pooled"][["scope", "n_units", "null_mean", "null_sd", "null_q95", "mde"]].round(4).to_string(index=False))
    return d


# ---- H138.5 ------------------------------------------------------------------------------------------------------

def one_study(study: str, n_perm: int, i: int, seed: int):
    from clocksync.config import limit_threads
    limit_threads(1)
    pb, mal = q.study_table(study)
    u = q.sc_units(study, pb, mal)
    res = {}
    for part, x in u.items():
        if x is None:
            continue
        j = 0 if part == "primary" else 1
        null = q.perm_null(x["za"], x["vA"], x["zb"], x["vB"], n_perm, np.random.default_rng([SEED, 5, seed, i, j]))
        res[part] = dict(est=q.estimable(x["za"], x["vA"], x["zb"], x["vB"], null), n=x["n"], null=null, host_source=x["host_source"], n_host=x["n_host"], malignant=mal)
    return study, res


def run_h138_5(args, n_perm):
    studies = list(q.sh.sc_builder().STUDIES)
    if args.subsample:
        studies = studies[:6]
    gate = q.a17_gate()
    out = Parallel(n_jobs=args.jobs)(delayed(one_study)(s, n_perm, i, args.seed) for i, s in enumerate(studies))
    rows, keep = [], {"primary": [], "fibroblast": []}
    for study, res in out:
        for part, r in res.items():
            keep[part].append((study, r))
            rows.append(null_rows(study, part, r["n"], r["null"], host_source=r["host_source"], n_host=r["n_host"],
                                  estimable=r["est"], malignant=r["malignant"] if part == "primary" else "", gate_a17=gate.get(study, np.nan)))
    for part, lst in keep.items():
        lst = [x for x in lst if x[1]["est"]]  # note 4 (ii)
        for scope, sel in (("pooled_all", lst), ("pooled_gate_passing", [x for x in lst if gate.get(x[0], False)])):
            if sel:
                pooled = q.pool(np.array([r["null"] for _, r in sel]), [r["n"] for _, r in sel])
                rows.append(null_rows(scope, part + "_pooled", sum(r["n"] for _, r in sel), pooled, n_cohorts=len(sel),
                                      cohorts=";".join(s for s, _ in sel)))
    d = pd.DataFrame(rows)
    d.to_csv(OUT / "h138_5_null.csv", index=False)
    log("H138.5 null done:", d[d.kind.str.endswith("pooled")][["scope", "kind", "n_units", "n_cohorts", "null_mean", "null_sd", "null_q95", "mde"]].round(4).to_string(index=False))
    return d


# ---- H138.1 ------------------------------------------------------------------------------------------------------

def sim_task(d: dict, ci: int, c: float, rep: int, n_boot: int, seed: int):
    from clocksync.config import limit_threads
    limit_threads(1)
    ss = int(np.random.SeedSequence([SEED, 1, seed, ci, int(round(c * 100)), rep]).generate_state(1)[0])
    t = q.simulate_design(d, c, ss)
    try:
        est, se = q.fit_coupling(t, d, n_boot, ss)
    except Exception as e:  # noqa: BLE001
        est, se = np.nan, np.nan
        print("fit failed", d["study"], c, rep, repr(e), flush=True)
    return dict(study=d["study"], c=c, rep=rep, est=est, se=se)


def run_h138_1(args, n_rep, n_boot):
    gate = q.a17_gate()
    designs = [dd for s in sorted(s for s, g in gate.items() if g) if (dd := q.q1_design(s))]
    if args.subsample:
        designs = designs[:3]
    pd.DataFrame([dict(study=d["study"], malignant=d["mal"], n_patients=d["n"], n_anchors=len(d["anchors"]),
                       depth_malignant=d["depth"][d["mal"]], depth_anchor_median=np.median([d["depth"][a] for a in d["anchors"]]),
                       depth_anchor_min=min(d["depth"][a] for a in d["anchors"]), depth_anchor_max=max(d["depth"][a] for a in d["anchors"]))
                  for d in designs]).to_csv(OUT / "h138_1_designs.csv", index=False)
    tasks = [(d, ci, c, r) for ci, d in enumerate(designs) for c in q.C_GRID for r in range(n_rep)]
    log(f"H138.1: {len(designs)} cohorts, {sum(d['n'] for d in designs)} patients, {len(tasks)} fits")
    res, last = [], time.time()
    gen = Parallel(n_jobs=args.jobs, return_as="generator_unordered", batch_size=1)(
        delayed(sim_task)(d, ci, c, r, n_boot, args.seed) for d, ci, c, r in tasks)
    for k, r in enumerate(gen, 1):
        res.append(r)
        if time.time() - last > 300 or k == len(tasks):
            log(f"H138.1 fits {k}/{len(tasks)}")
            last = time.time()
    coh = pd.DataFrame(res).sort_values(["study", "c", "rep"])
    coh.to_csv(OUT / "h138_1_sim_cohorts.csv", index=False)
    rows = []
    for c, g in coh.groupby("c"):
        pooled = [q.S.meta(x.est, x.se) for _, x in g.groupby("rep")]
        lo = np.array([p["lo"] for p in pooled])
        mr = np.array([p["random"] for p in pooled])
        ok = np.isfinite(lo)
        rows.append(dict(c=c, replicates=int(ok.sum()), share_lo_gt0=float(np.mean(lo[ok] > 0)), pooled_mean=float(mr[ok].mean()),
                         pooled_sd=float(mr[ok].std(ddof=1)), pooled_q95=float(np.quantile(mr[ok], 0.95)),
                         n_cohorts=int(np.mean([p["k"] for p in pooled])), n_boot=n_boot,
                         cohort_est_mean=float(g.est.mean()), cohort_fit_failed=int(g.est.isna().sum())))
    d = pd.DataFrame(rows)
    d.to_csv(OUT / "h138_1_sim.csv", index=False)
    log("H138.1 sim done:\n" + d.round(3).to_string(index=False))
    return d, designs


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--n-rep", type=int, default=100, help="H138.1 replicates per c (registered 100; note 4: 30)")
    ap.add_argument("--n-boot", type=int, default=200, help="H138.1 bootstraps per fit (registered 200; note 4: 100)")
    ap.add_argument("--only", nargs="*", default=["h138_1", "h138_4", "h138_5"], choices=["h138_1", "h138_4", "h138_5"])
    args = ap.parse_args(argv)
    run = cli.setup(args)
    global OUT
    if run["quick"]:
        OUT = OUT / "subsample"
    OUT.mkdir(parents=True, exist_ok=True)
    n_perm, n_rep, n_boot = (50, 2, 10) if run["quick"] else (q.N_PERM, args.n_rep, args.n_boot)
    summ = []
    if "h138_4" in args.only:
        d = run_h138_4(args, n_perm)
        p = d[d.scope == "pooled_gated"].iloc[0]
        summ.append(dict(hypothesis="H138.4", n_units=p.n_units, n_cohorts=p.n_cohorts, **{k: p[k] for k in ("null_mean", "null_sd", "null_q95", "min_p", "mde")},
                         note=f"pooled over gated cancers ({p.cohorts}); |corr| Poisson-corrected, weights n_pairs; {int(p.n_perm)} perms"))
    if "h138_5" in args.only:
        d = run_h138_5(args, n_perm)
        for kind, name in (("primary_pooled", "H138.5"), ("fibroblast_pooled", "H138.5 fibroblast companion")):
            s = d[(d.kind == kind) & (d.scope == "pooled_all")]
            if len(s):
                p = s.iloc[0]
                summ.append(dict(hypothesis=name, n_units=p.n_units, n_cohorts=p.n_cohorts, **{k: p[k] for k in ("null_mean", "null_sd", "null_q95", "min_p", "mde")},
                                 note=f"cohorts with >= {q.MIN_PAIRED} paired patients and a stromal host ({p.cohorts}); weights n"))
    if "h138_1" in args.only:
        d, designs = run_h138_1(args, n_rep, n_boot)
        z = d[d.c == 0].iloc[0]
        ok = d[d.share_lo_gt0 >= 0.8]
        mde = float(ok.c.min()) if len(ok) else np.nan
        summ.append(dict(hypothesis="H138.1", n_units=sum(x["n"] for x in designs), n_cohorts=len(designs), null_mean=z.pooled_mean,
                         null_sd=z.pooled_sd, null_q95=z.pooled_q95, min_p=np.nan, mde=mde,
                         note=(f"simulated pooled DL coupling at c=0; false-positive rate at c=0 = {z.share_lo_gt0:.3f}; "
                               + ("MDE = smallest c on the grid with share(lower bound > 0) >= 0.8" if len(ok) else
                                  f"underpowered at this design (no c <= {max(q.C_GRID)} reaches 0.8)"))))
    if summ:
        f = OUT / "power_summary.csv"
        new = pd.DataFrame(summ)
        if f.exists():
            old = pd.read_csv(f)
            new = pd.concat([old[~old.hypothesis.isin(new.hypothesis)], new], ignore_index=True)
        new.to_csv(f, index=False)
        log("power_summary:\n" + new.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
