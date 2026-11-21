"""A328 part a power stage (registered in `paper/docs/registered_experiments.md`, A328 and note 1), before any del Valle
clock statistic: the A242 feasibility on the del Valle design.

Read from the fetal data: the sheet (stage, sex, library size) and the cohort median count of each clock gene
(`prep_a328.py`). Each of the 32 fetuses is stood for by one distinct GTEx adrenal gland (death at 08-18 h, any Hardy
class) that receives its stage, sex and library size; the adult clock-gene counts are thinned binomially gene by gene
to target = cohort median count x fetal library size / median library size (lowered by one common factor where the
adult libraries hold fewer reads, `cap`, as A242). Adjustment and test as registered: log fractions residualised on
the stage indicator and sex, signed coherence C, per-gene permutations within stage.

Settings: `a1`, `a0.75`, `a0.5`, the adult clock at relative amplitude a (`clocksync.simulate.scale_clock` with the
14-gene waveform); `shuffled`, every clock gene permuted across the drawn donors (false-positive rate). Tests: `dv:all`
(H328.1, 14 genes) and `dv:core` (H328.2, the 11-gene core waveform). A readout with pass rate < 0.8 at a1 is "not
estimable" and is not computed. Controlled access: GTEx time of death in memory only; outputs are aggregates.

    python power_a328a.py [--n-sub 200] [--jobs 8] [--check]
Outputs (`../output/power_a/`): design.csv, subsamples.csv, passrates.csv (`_check` with --check).
"""
import argparse
import sys
import time
from multiprocessing import Pool
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[3]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import binomtest  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.simulate import scale_clock  # noqa: E402
from common import development as dev  # noqa: E402
from common import gtex  # noqa: E402
from common.benchmark.data import in_window  # noqa: E402

PREP = HERE.parents[1] / "output" / "prep"
OUT = HERE.parents[1] / "output" / "power_a"
SEED, N_PERM = 328, 500
WF = cs.gtex()
GENES = list(WF.genes)
TESTS = {"dv:all": WF, "dv:core": WF.drop(["PER1", "PER2", "NR1D1"])}
SCALES = (1.0, 0.75, 0.5)
S = {}  # state inherited by the workers


def cov(rows) -> np.ndarray:
    """Stage indicator (saturated age model of the four stages) and sex."""
    return dev.design(rows, age=None, factors=("stage", "sex"))


def setup() -> pd.DataFrame:
    obs = pd.read_csv(PREP / "dv_sheet.csv")
    obs["stage"] = obs.stage.astype(str)
    det = pd.read_csv(PREP / "dv_clock_detection.csv", index_col=0).loc["delvalle", GENES].to_numpy(float)
    pb, d = gtex.load(["Adrenal Gland"], hardy=(0, 1, 2, 3, 4), genes="all", depth="library_size")
    pb = pb[pb.donor.isin(set(d.index[in_window(d.tod.to_numpy())]))].reset_index(drop=True)
    depth, rng = pb.depth.to_numpy(float), np.random.default_rng([SEED, 0])
    clock = pb[["donor", "compartment", "depth", *GENES]]
    fmed = np.median(clock[GENES].to_numpy(float) / depth[:, None], axis=0)
    tgt = det[None, :] * (obs.depth.to_numpy(float) / obs.depth.median())[:, None]
    capf = np.minimum(1.0, fmed * depth.min() / tgt.max(0))
    S.update(obs=obs, depth=depth, target=tgt * capf, fmed=fmed,
             x={a: scale_clock(clock, WF, a, depth=clock.depth, rng=rng)[GENES].to_numpy(np.int64) for a in SCALES})
    return pd.DataFrame([dict(test=k, n_fetuses=len(obs), n_covariates=cov(obs).shape[1], pool_donors=len(pb),
                              median_cap=float(np.median(capf)), share_capped=float(np.mean(capf < 1)))
                         for k in TESTS])


def stand_in(x, depth, rng) -> pd.DataFrame:
    p = S["target"] / (S["fmed"] * depth[:, None])
    t = S["obs"][["donor", "compartment", "stage", "sex", "depth"]].copy()
    t[GENES] = rng.binomial(x, np.minimum(p, 1.0)).astype(float)
    return t


def one(i) -> list:
    rng, out = np.random.default_rng([SEED, 1, i]), []
    pick = rng.choice(len(S["depth"]), len(S["obs"]), replace=False)
    depth = S["depth"][pick]
    null = np.column_stack([np.round(rng.permutation(c / depth) * depth) for c in S["x"][1.0][pick].T])
    for setting, x in [*((f"a{a:g}", S["x"][a][pick]) for a in SCALES), ("shuffled", null.astype(np.int64))]:
        tab = dev.adjust(stand_in(x, depth, rng), GENES, cov)[0]
        for j, (name, wf) in enumerate(TESTS.items()):
            c, p, _ = dev.coherence(tab, wf, N_PERM, np.random.default_rng([SEED, 2, i, j]), "stage")
            out.append(dict(test=name, setting=setting, sub=i, C=c, p=p))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--n-sub", type=int, default=200)
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    t0, sfx = time.time(), "_check" if a.check else ""
    limit_threads(1)
    OUT.mkdir(parents=True, exist_ok=True)
    pd.set_option("display.width", 200)
    design = setup()
    design.to_csv(OUT / f"design{sfx}.csv", index=False)
    print(design.round(3).to_string(index=False), flush=True)
    n, res = 4 if a.check else a.n_sub, []
    with Pool(a.jobs, initializer=limit_threads, initargs=(1,)) as pool:
        for k, r in enumerate(pool.imap_unordered(one, range(n))):
            res.extend(r)
            if (k + 1) % 20 == 0 or k + 1 == n:
                print(f"{k + 1}/{n} {time.time() - t0:.0f}s", flush=True)
    res = pd.DataFrame(res).sort_values(["test", "setting", "sub"])
    res.to_csv(OUT / f"subsamples{sfx}.csv", index=False)
    rows = []
    for (test, setting), r in res.groupby(["test", "setting"], sort=False):
        k, m = int((r.p < .05).sum()), len(r)
        ci = binomtest(k, m).proportion_ci(0.95, method="exact")
        rows.append(dict(test=test, setting=setting, n_subsamples=m, pass_rate=k / m, ci_lo=ci.low, ci_hi=ci.high,
                         median_C=r.C.median(), C_q05=r.C.quantile(.05)))
    rates = pd.DataFrame(rows)
    rates.to_csv(OUT / f"passrates{sfx}.csv", index=False)
    print(rates.round(3).to_string(index=False), f"runtime {time.time() - t0:.0f}s", sep="\n", flush=True)


if __name__ == "__main__":
    main()
