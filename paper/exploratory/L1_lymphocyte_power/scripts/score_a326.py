"""A326: could the human lymphocyte designs (A281, A287, A321) have read a lymphocyte clock of mouse size?

Registered in paper/docs/registered_experiments.md (A326). In each design's atlas and compartment pair, the partner
rows' clock-plane states are scaled by r (`clocksync.simulate.scale_clock`, counts: binomial thinning, in-plane noise
scaled with the clock), then thinned to the lymphocyte compartment's median clock-gene count exactly as the design's own
depth control, and read with the design's own statistic (A281: C, 200 permutations; A287: U, pi, p_perm; A321: core-loop
C, 200 permutations). 20 replicates per atlas and r, seed 3260; replicate k uses the generator [3260, design, atlas, r, k]
(for A321 the replicate k is the same in every compartment: the joint replicate).

  validate  r = 1 only: pass rates against the registered thinned controls (A281 pairs.csv, A287 thinned.csv, A321 compartments.csv)
  score     the full r grid: pass_rates.csv, design_power.csv (and replicates.csv)

    python score_a326.py --stage validate | score [--jobs 8]
Outputs (../output/): validation.csv | pass_rates.csv, design_power.csv, replicates.csv
"""
from __future__ import annotations

import importlib.util
import sys
import time
import warnings
from pathlib import Path

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
sys.path.insert(0, str(PAPER))
import cloudpickle  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402
from scipy.stats import norm  # noqa: E402

from clocksync.simulate import scale_clock  # noqa: E402
from common import cli  # noqa: E402


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


EXP = PAPER / "exploratory"
A281 = _load("a326_a281", EXP / "A281_lymphocyte_clock/scripts/score_a281.py")
A287 = _load("a326_a287", EXP / "A287_lymphocyte_covariance/scripts/score_a287.py")
A321 = _load("a326_a321", EXP / "A321_lymphocyte_core_blind/scripts/score_a321.py")
for _m in (A281, A281.A280, A281.A274, A281.A272, A281.A150, A287, A287.A281, A321, A321.A294, A321.A281):
    cloudpickle.register_pickle_by_value(_m)

OUT = HERE.parents[1] / "output"
SEED, N_REP, CHUNK = 3260, 20, 4
RS = [0.10, 0.25, 0.40, 0.55, 0.70, 1.0]
CLOCK, WF, CORE = A281.CLOCK, A281.WF, A321.CORE
DESIGNS = ["A281", "A287", "A321"]
N_PERM_THIN = A281.N_PERM_THIN
WORK_LOG = time.time()


def log(*a):
    print(time.strftime("%H:%M:%S"), *a, flush=True)


def paired(tl, ts):
    d = sorted(set(tl.donor) & set(ts.donor))
    tl = tl[tl.donor.isin(d)].sort_values("donor").reset_index(drop=True)
    ts = ts[ts.donor.isin(d)].sort_values("donor").reset_index(drop=True)
    return tl, ts


def scaled(ts, r, rng):
    """Partner rows with the clock-plane states scaled by r (counts; the design's depth column)."""
    return scale_clock(ts, WF, r, depth=ts.depth.to_numpy(float), rng=rng)


def rep_task(design, key, r, reps, di, ai, ri, ts, target, extra):
    """Replicates `reps` of one (design, atlas / compartment, r): list of dict(rep, p[, structure])."""
    warnings.filterwarnings("ignore")
    out = []
    for k in reps:
        rng = np.random.default_rng([SEED, di, ai, ri, k])
        tt = scaled(ts, r, rng)
        if design == "A281":
            tt = A281.thin(tt, target, rng)
            p = A281.perm_test(A281.y_of(tt, CLOCK), A281.groups_of(tt), N_PERM_THIN, rng)[1]
            out.append(dict(rep=k, p=p, passed=p < 0.05))
        elif design == "A321":
            tt = A281.thin(tt, target, rng)
            p = A281.perm_test(A281.y_of(tt, list(CORE.genes)), A281.groups_of(tt), N_PERM_THIN, rng, wf=CORE)[1]
            out.append(dict(rep=k, p=p, passed=p < 0.05))
        else:  # A287: thin every column of the matched sets, as A287 does
            genes = extra["genes"]
            s = tt[CLOCK].sum(1).to_numpy(float)
            pr = np.minimum(1.0, target / np.maximum(s, 1))
            tt[genes] = rng.binomial(tt[genes].to_numpy(np.int64), pr[:, None])
            tt["depth"] = tt.depth.to_numpy(float) * pr
            res, _ = A287.readout(tt, extra["matched"], N_PERM_THIN, rng)
            out.append(dict(rep=k, p=res["p_perm"], passed=bool(res["structure"]), pi=res["pi"]))
    return dict(design=design, key=key, r=r, rows=out)


def build_inputs(designs=DESIGNS):
    """(tasks specs, registered controls): one spec per (design, atlas / compartment), data loaded once."""
    A281.A150.ATLASES.update(A281.ATLASES)
    keys = list(A281.ATLASES)
    prs = pd.read_csv(EXP / "A281_lymphocyte_clock/output/pairs.csv")
    specs = []
    for ai, row in enumerate(prs.itertuples()):
        k = keys.index(row.atlas)
        A = A281.A150.load_atlas(row.atlas, None, 0)
        pb = A["pb"].reset_index(drop=True)
        log("loaded", row.atlas)
        matched = [list(w.genes) for w in A281.A272.pseudo_sets(A, k, A287.N_SETS, "matched")]
        base = ["donor", "compartment", "stratum", "n_cells", "depth"]
        tl, ts = paired(pb[pb.compartment == row.lymphocyte][base + CLOCK], pb[pb.compartment == row.stromal][base + CLOCK])
        if "A281" in designs:
            specs.append(dict(design="A281", key=row.atlas, ai=ai, ts=ts, target=float(tl[CLOCK].sum(1).median()),
                              extra={}, registered=float(row.thinned_pass_rate)))
        # A287: same pairs; same matched sets as A287 (A281's), columns of the sets kept
        if "A287" not in designs:
            continue
        genes = sorted({g for s in matched for g in s} | set(CLOCK))
        tl2, ts2 = paired(pb[pb.compartment == row.lymphocyte][base + genes], pb[pb.compartment == row.stromal][base + genes])
        specs.append(dict(design="A287", key=row.atlas, ai=ai, ts=ts2, target=float(tl2[CLOCK].sum(1).median()),
                          extra=dict(matched=matched, genes=genes), registered=np.nan))
    th = pd.read_csv(EXP / "A287_lymphocyte_covariance/output/thinned.csv").groupby("atlas").structure.mean()
    for s in specs:
        if s["design"] == "A287":
            s["registered"] = float(th[s["key"]])
    comp = pd.read_csv(EXP / "A321_lymphocyte_core_blind/output/compartments.csv")
    for k, name in enumerate(A321.ATLASES if "A321" in designs else []):
        d = comp[comp.atlas == name]
        if not len(d):
            continue
        pb = A321.load(name)
        log("loaded", name)
        matched = A321.matched_sets(pb, k)
        cols = sorted({g for w in matched for g in w.genes} | set(CLOCK))
        base = ["donor", "compartment", "stratum", "n_cells", "depth"]
        for r in d.itertuples():
            tl, ts = paired(pb[pb.compartment == r.lymphocyte][base + CLOCK], pb[pb.compartment == r.partner][base + CLOCK])
            specs.append(dict(design="A321", key=f"{name}|{r.lymphocyte}", ai=len(specs), ts=ts,
                              target=float(tl[CLOCK].sum(1).median()), extra=dict(informative=bool(r.primary and r.core_power >= 0.8)),
                              registered=float(r.core_power)))
    return specs


def run(specs, rs, jobs):
    tasks = []
    for s in specs:
        di = DESIGNS.index(s["design"])
        for r in rs:
            ri = RS.index(r)
            for c in range(0, N_REP, CHUNK):
                tasks.append(delayed(rep_task)(s["design"], s["key"], r, range(c, c + CHUNK), di, s["ai"], ri, s["ts"],
                                               s["target"], s["extra"]))
    log(f"{len(tasks)} tasks, {jobs} workers")
    res, t0 = [], time.time()
    for i, x in enumerate(Parallel(n_jobs=jobs, return_as="generator_unordered")(tasks), 1):
        res.append(x)
        if i % 5 == 0 or i == len(tasks):
            log(f"{i}/{len(tasks)} tasks, {(time.time() - t0) / 60:.1f} min")
    rows = [dict(design=x["design"], key=x["key"], r=x["r"], **q) for x in res for q in x["rows"]]
    return pd.DataFrame(rows).sort_values(["design", "key", "r", "rep"]).reset_index(drop=True)


def poibin_gt(p, k):
    """P(number of successes > k) of independent Bernoulli(p)."""
    pm = np.array([1.0])
    for q in p:
        pm = np.convolve(pm, [1 - q, q])
    return float(pm[k + 1:].sum())


def pass_rates(reps):
    return reps.groupby(["design", "key", "r"]).passed.mean().rename("pass_rate").reset_index()


def design_power(reps, pr, specs):
    info = {s["key"]: s["registered"] for s in specs}
    out = []
    for r in sorted(pr.r.unique()):
        a = pr[(pr.design == "A281") & (pr.r == r)].set_index("key").pass_rate
        pw = a[[k for k in a.index if info[k] >= 0.8]]  # A281 powered atlases: registered thinned pass rate >= 0.8
        out.append(dict(design="A281", r=r, power=poibin_gt(pw.to_numpy(), 2), n_powered=len(pw)))
        a = pr[(pr.design == "A287") & (pr.r == r)].set_index("key").pass_rate
        pw = a[[k for k in a.index if info[k] >= 0.8]]
        out.append(dict(design="A287", r=r, power=float(1 - np.prod(1 - pw.to_numpy())), n_powered=len(pw)))
        keys = [s["key"] for s in specs if s["design"] == "A321" and s["extra"]["informative"]]
        d = reps[(reps.design == "A321") & (reps.r == r) & reps.key.isin(keys)]
        z = d.groupby("rep").p.apply(lambda p: norm.isf(p.to_numpy()).sum() / np.sqrt(len(p)))
        out.append(dict(design="A321", r=r, power=float(np.mean(norm.sf(z) < 0.05)), n_powered=len(keys)))
    return pd.DataFrame(out).rename(columns={"power": "design_power"})


def main():
    global N_REP
    ap = cli.parser(__doc__)
    ap.add_argument("--stage", choices=["validate", "score"], required=True)
    ap.add_argument("--reps", type=int, default=N_REP, help="replicates (registered: 20; other values are diagnostics)")
    ap.add_argument("--designs", default=",".join(DESIGNS))
    args = ap.parse_args()
    cli.setup(args)
    N_REP = args.reps
    OUT.mkdir(parents=True, exist_ok=True)
    specs = build_inputs(args.designs.split(","))
    rs = [1.0] if args.stage == "validate" else RS
    reps = run(specs, rs, args.jobs)
    pr = pass_rates(reps)
    reg = pd.DataFrame([dict(design=s["design"], key=s["key"], registered=s["registered"]) for s in specs])
    pd.set_option("display.width", 250)
    if args.stage == "validate":
        v = pr.merge(reg, on=["design", "key"]).rename(columns={"pass_rate": "ours"})
        v["difference"] = v.ours - v.registered
        v.drop(columns="r").to_csv(OUT / ("validation.csv" if N_REP == 20 else f"validation_diag_{N_REP}reps.csv"), index=False)
        log("\n" + v.round(3).to_string())
        log("max |difference|", float(v.difference.abs().max()))
        return
    pr.to_csv(OUT / "pass_rates.csv", index=False)
    reps.to_csv(OUT / "replicates.csv", index=False)
    dp = design_power(reps, pr, specs)
    dp.to_csv(OUT / "design_power.csv", index=False)
    log("\n" + pr.round(3).to_string() + "\n" + dp.round(3).to_string())


if __name__ == "__main__":
    main()
