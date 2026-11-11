"""A240 (registered in `paper/docs/registered_experiments.md`): how small, or how synchronised, would a fetal clock
have to be to stay undetected in A237-A239?

The adult feasibility of A237-A239 (adult material cut to each fetal design, `scouts_2026-11-10/feasibility_fetal`)
is rerun, on the same subsamples, with the adult clock made smaller: in every tissue or cell class of the adult pool
the clock-plane state of each sample, centred over the pool, is multiplied by a scale a
(`clocksync.simulate.scale_clock`; counts by binomial thinning, arrays in log space). The pass rate of a readout as a
function of a gives the largest clock it would have missed: a80 is the smallest a of the grid with a pass rate >= 0.8
at it and at every larger a. No fetal expression value is read; the fetal statistics already scored (A237-A239
`l1.csv`) are only compared with the adult ones (`below_obs`: share of adult subsamples whose C is at or below the
fetal C).

Arms (one process each, so that jobs keep the indices, hence the seeds, of the registered feasibility runs)
- cm:     GTEx donors who died at 08-18 h at the Cardoso-Moreira designs (H237.1 L2, H237.2 and H237.3 L1).
- gonad:  adult breast nuclei at the Garcia-Alonso design, thinned to the fetal depths (H238.1 L1, H238.2 L2).
- kang:   GTEx brains (any hour, 12 genes) at the Kang prenatal design (H239.2 L1).
- array:  the postnatal brains of Kang 2011 themselves (exon arrays, 32 brains, the passed control of A239), 20 at a
          time with the registered adjustment: the same platform as the prenatal negative.
- window: the cm L1 readouts at full amplitude with the adult donors restricted to deaths at 08-14 h and 08-12 h,
          where the pool holds >= 1.5 x the fetal n: what a narrow collection window would cost.
At a = 1 the first three arms must reproduce the registered pass rates exactly (`check_<arm>.csv`; the run stops
otherwise).

Controlled access: GTEx time of death is used in memory to select donors; only aggregates are written.

    python amplitude_bound.py --arm cm|gonad|kang|array|window|all [--n-sub 200] [--jobs 8] [--check]

Outputs (`../output/`; `_check` with --check: 3 subsamples, scales 1, 0.5, 0): subsamples_<arm>.csv, curve_<arm>.csv
(arm, label, scale or window, pass rate with exact 95 % CI, kept share of the clock-gene counts, below_obs),
check_<arm>.csv, and with --arm all: curve.csv, bound.csv.
"""
import argparse
import subprocess
import sys
import time
from multiprocessing import Pool
from pathlib import Path

HERE = Path(__file__).resolve()
EXPL = HERE.parents[2]
sys.path.insert(0, str(HERE.parents[3]))
for p in ("scouts_2026-11-10/feasibility_fetal/scripts", "A239_fetal_brain/scripts"):
    sys.path.insert(0, str(EXPL / p))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import feasibility_fetal as ff  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.simulate import scale_clock  # noqa: E402
from common import development as dev  # noqa: E402
from common import gtex  # noqa: E402
from common.benchmark.data import in_window  # noqa: E402

OUT = HERE.parents[1] / "output"
FEAS = EXPL / "scouts_2026-11-10/feasibility_fetal/output"
SEED = 240
SCALES = (1.0, 0.75, 0.5, 0.35, 0.25, 0.15, 0.0)
WINDOWS = {"08-14": (8.0, 14.0), "08-12": (8.0, 12.0)}  # h, GTEx time of death
ARMS = ("cm", "gonad", "kang", "array", "window")
CM_L1 = [f"{o}:4-20wpc" for o in ("heart", "liver", "forebrain", "hindbrain", "kidney", "ovary")] + [
    f"{o}:4-8wpc" for o in ("heart", "liver", "hindbrain", "kidney")] + [
    f"{o}:13-20wpc" for o in ("heart", "liver", "kidney")]  # H237.2, H237.3
LABELS = {"cm": CM_L1 + ["organs:9-20wpc"],
          "gonad": [f"breast:{c}:all" for c in ("mesenchymal_stromal", "supporting", "endothelial", "epithelial")]
          + ["breast:classes4:all"],
          "kang": ["brain_mean:prenatal"], "array": ["brain_mean:postnatal@20"], "window": CM_L1}
REGISTERED = {"cm": "passrates_cm_day_sex.csv", "gonad": "passrates_gonad_sextrend.csv",
              "kang": "passrates_kang_all.csv"}  # the feasibility runs of A237, A238, A239
OBSERVED = {"cm": "A237_fetal_organs", "gonad": "A238_fetal_gonad", "kang": "A239_fetal_brain",
            "array": "A239_fetal_brain"}  # fetal C of the readout (l1.csv), for `below_obs`
N_ARRAY, N_PERM = 20, 500  # brains per array subsample (the prenatal n); permutations, as in the feasibility
POOLS = {}  # the array arm: scale -> (intensities, sheet) of the postnatal brain means
WF_ARRAY = COV_ARRAY = None  # its waveform (12 genes) and covariates, set by `arm_array`


def scaled(pool, a, wf, rng):
    """Pool {donor: table} with the clock of every compartment scaled to `a` over the donors of the pool, and the
    share of the clock-gene counts kept by the thinning."""
    if a == 1:
        return pool, 1.0
    tab, g = pd.concat(pool.values(), ignore_index=True), list(wf.genes)
    out = pd.concat([scale_clock(t, wf, a, depth=t.depth, rng=rng) for _, t in
                     tab.groupby("compartment", sort=False)]).sort_index()
    kept = out[g].to_numpy().sum() / tab[g].to_numpy().sum()
    return {d: t.reset_index(drop=True) for d, t in out.groupby("donor", sort=False)}, kept


def sweep(one, tasks, settings, set_pool, jobs, t0):
    """`one` over the tasks at every setting (the workers are forked after `set_pool(setting)`); adjusted rows."""
    res = []
    for s in settings:
        info = set_pool(s)
        with Pool(jobs, initializer=limit_threads, initargs=(1,)) as p:
            for k, r in enumerate(p.imap_unordered(one, tasks, chunksize=2)):
                res.extend(dict(x, setting=s, **info) for x in r if x.get("adjust", "adjusted") == "adjusted")
                if (k + 1) % 200 == 0:
                    print(f"setting {s}: {k + 1}/{len(tasks)} {time.time() - t0:.0f}s", flush=True)
        print(f"setting {s} done {time.time() - t0:.0f}s {info}", flush=True)
    return pd.DataFrame(res)


def observed(arm):
    """Fetal C of the L1 readouts of the arm, by feasibility label."""
    if arm not in OBSERVED:
        return {}
    c = pd.read_csv(EXPL / OBSERVED[arm] / "output" / "l1.csv").set_index("test").C
    name = {"gonad": lambda x: x.split(":", 1)[1], "array": lambda x: "brain_mean:prenatal"}.get(arm, lambda x: x)
    return {x: c[name(x)] for x in LABELS[arm] if name(x) in c.index}


def curve(res, arm):
    obs, rows = observed(arm), []
    for (label, s), r in res.groupby(["label", "setting"], sort=False):
        k, m = int(r.passed.sum()), len(r)
        lo, hi = ff.ci(k, m)
        rows.append(dict(arm=arm, kind=r.kind.iloc[0], label=label, setting=s, n=int(r.n.median()), n_subsamples=m,
                         pass_rate=k / m, ci_lo=lo, ci_hi=hi, median_stat=r.stat.median(),
                         estimable=r.p.notna().mean(), kept=r.kept.iloc[0], pool=r.pool.iloc[0],
                         fetal_C=obs.get(label, np.nan),
                         below_obs=(r.stat <= obs[label]).mean() if label in obs else np.nan))
    return pd.DataFrame(rows)


def reproduce(res, arm, sfx):
    """Pass rates at a = 1 against the registered feasibility run, on the subsamples both have; stops on a
    difference."""
    reg = pd.read_csv(FEAS / REGISTERED[arm].replace("passrates", "subsamples"))
    reg = reg[reg.adjust == "adjusted"].set_index(["label", "sub"]).passed
    new = res[res.setting == 1].set_index(["label", "sub"]).passed
    both = new.index.intersection(reg.index)
    chk = pd.DataFrame(dict(registered=reg[both].astype(bool), rerun=new[both].astype(bool))).groupby("label").agg(
        n=("rerun", "size"), registered=("registered", "mean"), rerun=("rerun", "mean"),
        differ=("rerun", lambda x: int((x != reg[x.index].astype(bool)).sum())))
    chk.to_csv(OUT / f"check_{arm}{sfx}.csv")
    print(chk.to_string(), flush=True)
    if chk.differ.sum():
        raise SystemExit(f"{arm}: the rerun at a = 1 differs from the registered feasibility")


def tasks_of(jobs, labels, n):
    missing = set(labels) - {j[1] for j in jobs}
    if missing:
        raise SystemExit(f"readouts without a calibration job: {sorted(missing)}")
    return [(j, i) for j in range(len(jobs)) if jobs[j][1] in labels for i in range(n)]


def arm_gtex(arm, scales, n, jobs, t0):
    """cm and kang: GTEx pools of `feasibility_fetal`, as in the registered runs."""
    if arm == "cm":
        ff.SEX = True
        ff.setup_cm("day")
    else:
        ff.WF = ff.WF.drop(["PER3", "NR1D1"])
        ff.GENES = list(ff.WF.genes)
        ff.setup_kang("all")
    base = dict(ff.POOL)

    def set_pool(a):
        new, kept = scaled(base, a, ff.WF, np.random.default_rng([SEED, ARMS.index(arm), int(round(100 * a))]))
        ff.POOL.clear()
        ff.POOL.update(new)
        return dict(kept=kept, pool=len(new))
    return sweep(ff.one, tasks_of(ff.JOBS, LABELS[arm], n), scales, set_pool, jobs, t0)


def arm_gonad(scales, n, jobs, t0):
    import feasibility_gonad as fg
    fg.SEX_TREND = True
    fg.setup()
    base = fg.POOLS["breast"]

    def set_pool(a):
        fg.POOLS["breast"], kept = scaled(base, a, ff.WF, np.random.default_rng([SEED, 1, int(round(100 * a))]))
        return dict(kept=kept, pool=len(base))
    return sweep(fg.one, tasks_of(fg.JOBS, LABELS["gonad"], n), scales, set_pool, jobs, t0)


def arm_window(n, jobs, t0):
    """cm L1 readouts at full amplitude, adult donors of a narrower death-time window (pool >= 1.5 x fetal n)."""
    ff.SEX = True
    ff.setup_cm("day")
    base = dict(ff.POOL)
    _, d = gtex.load(list(ff.ORGAN.values()), hardy=(0, 1, 2, 3, 4), depth="library_size")
    res = []
    for tag, w in WINDOWS.items():
        keep = set(d.index[in_window(d.tod.to_numpy(), w)])
        pool = {k: t for k, t in base.items() if k in keep}
        size = pd.concat(pool.values()).groupby("compartment").donor.nunique()
        labels = [x for x in LABELS["window"] if size.get(x.split(":")[0], 0) >= ff.POOL_FACTOR * next(
            j[2] for j in ff.JOBS if j[1] == x)]
        print(f"window {tag}: {len(pool)} donors, readouts {labels}", flush=True)

        def set_pool(s, pool=pool):
            ff.POOL.clear()
            ff.POOL.update(pool)
            return dict(kept=1.0, pool=len(pool))
        res.append(sweep(ff.one, tasks_of(ff.JOBS, labels, n), [tag], set_pool, jobs, t0))
    return pd.concat(res, ignore_index=True)


def one_array(args):
    a, i = args
    x, obs = POOLS[a]
    rng = np.random.default_rng([SEED, 3, i])  # the same brains and permutations at every scale
    pick = obs.loc[rng.choice(obs.index.to_numpy(), size=N_ARRAY, replace=False)]
    tab, _ = dev.adjust(dev.table(x, pick, WF_ARRAY.genes), WF_ARRAY.genes, COV_ARRAY)
    c, p, _ = dev.coherence(tab, WF_ARRAY, N_PERM, rng)
    return [dict(kind="L1", label=LABELS["array"][0], sub=i, n=N_ARRAY, stat=c, p=p, passed=p < .05)]


def arm_array(scales, n, jobs, t0):
    """Postnatal brain means of Kang 2011 (the control of A239), 20 of the 32 brains per subsample, registered
    adjustment (quadratic of log age, sex); the clock is scaled in log2 intensity over the 32 brains."""
    import score_a239 as a239
    global WF_ARRAY, COV_ARRAY
    WF_ARRAY, COV_ARRAY = a239.WF, a239.cov(2)
    x, obs = a239.load()
    lo, hi = a239.GROUPS["postnatal"]
    obs = obs[(obs.compartment == "brain_mean") & obs.stage.between(lo, hi)]
    logx = np.log2(x.loc[obs.index, list(WF_ARRAY.genes)])
    res = []
    for a in scales:
        POOLS[a] = (2.0 ** scale_clock(logx, WF_ARRAY, a), obs)
        res.append(sweep(one_array, [(a, i) for i in range(n)], [a], lambda s: dict(kept=1.0, pool=len(obs)), jobs,
                         t0))
    return pd.concat(res, ignore_index=True)


def bound(cur):
    """Per readout: pass rate at full amplitude; a80 (a50): the smallest scale with a pass rate >= 0.8 (0.5) at it
    and at every larger scale; a_excl: the same for below_obs <= 0.05 (the fetal C lies in the lowest 5 % of the
    adult subsamples)."""
    def smallest(c, ok):
        a = np.nan
        for s, v in zip(c.setting, ok):
            if not v:
                break
            a = s
        return a
    rows = []
    for (arm, label), c in cur[cur.arm != "window"].groupby(["arm", "label"], sort=False):
        c = c.astype({"setting": float}).sort_values("setting", ascending=False)
        rows.append(dict(arm=arm, kind=c.kind.iloc[0], label=label, n=c.n.iloc[0], pass_rate_full=c.pass_rate.iloc[0],
                         a80=smallest(c, c.pass_rate >= 0.8), a50=smallest(c, c.pass_rate >= 0.5),
                         a_excl=smallest(c, c.below_obs <= 0.05) if c.kind.iloc[0] == "L1" else np.nan,
                         pass_rate_zero=c.pass_rate.iloc[-1], fetal_C=c.fetal_C.iloc[0]))
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--arm", choices=(*ARMS, "all"), required=True)
    ap.add_argument("--n-sub", type=int, default=200)
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    sfx, t0 = "_check" if a.check else "", time.time()
    if a.arm == "all":
        for arm in ARMS:
            subprocess.run([sys.executable, "-u", str(HERE), "--arm", arm, "--n-sub", str(a.n_sub), "--jobs",
                            str(a.jobs)] + ["--check"] * a.check, check=True)
            print(f"arm {arm} done {time.time() - t0:.0f}s", flush=True)
        cur = pd.concat([pd.read_csv(OUT / f"curve_{arm}{sfx}.csv") for arm in ARMS], ignore_index=True)
        cur.to_csv(OUT / f"curve{sfx}.csv", index=False)
        b = bound(cur)
        b.to_csv(OUT / f"bound{sfx}.csv", index=False)
        print(b.round(3).to_string(index=False), flush=True)
        return
    limit_threads(1)
    n, scales = (3, (1.0, 0.5, 0.0)) if a.check else (a.n_sub, SCALES)
    if a.arm in ("cm", "kang"):
        res = arm_gtex(a.arm, scales, n, a.jobs, t0)
    elif a.arm == "gonad":
        res = arm_gonad(scales, n, a.jobs, t0)
    elif a.arm == "array":
        res = arm_array(scales, n, a.jobs, t0)
    else:
        res = arm_window(n, a.jobs, t0)
    res.to_csv(OUT / f"subsamples_{a.arm}{sfx}.csv", index=False)
    if a.arm in REGISTERED:
        reproduce(res, a.arm, sfx)
    cur = curve(res, a.arm)
    cur.to_csv(OUT / f"curve_{a.arm}{sfx}.csv", index=False)
    print(cur.round(3).to_string(index=False), flush=True)
    print(f"runtime {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
