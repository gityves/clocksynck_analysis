"""A332 D4 (registered in `paper/docs/registered_experiments.md`, section A332): the developmental loudness curve.

Clock-plane variance V = mean |z|^2 of the centred clock-plane states, per organ and age bin, beside signed coherence
C, in Cardoso-Moreira 2019 (7 organs, bulk RNA-seq; loaders of `A237_fetal_organs/scripts/score_a237.py`) and Kang
2011 brain (exon arrays, mean over the regions of a brain; loaders of `A239_fetal_brain/scripts/score_a239.py`).

Bins (fixed): prenatal < 12 weeks post conception (wpc); prenatal >= 12 wpc; infancy and childhood (birth to 12 y);
adolescence and adulthood (>= 13 y). A bin enters with >= 8 individuals.
Ages: Cardoso-Moreira prenatal = `age_wpc` (weeks post conception); postnatal = parsed from `stage_detail`
("N days | months | year(s) post birth", ages >= 13 y are years; 1 y = 365.25 d). Kang = `age_value` with `age_unit`
PCW (weeks post conception), PMonth (months after birth), Y (years after birth); 12 PCW falls in the second bin.

Per bin: log fractions are residualised within the bin, gene by gene, on the covariates of A237 / A239 restricted to
the bin (constant terms dropped): sex, and a polynomial of age (Kang: log age as in A239, df 2 prenatal, df 1
postnatal; Cardoso-Moreira prenatal: wpc, df 2, df 3 when the bin spans > 8 weeks as in A237; Cardoso-Moreira
postnatal: log2(age in years + 0.1), df 2 for birth-12 y, df 1 for >= 13 y). V uses one projector per table (GTEx
waveform; Kang 12 genes as in A239: PER3, NR1D1 absent), counts path (`paired.projector`, `paired.states`), without
the Poisson correction. Pseudo sets: 200 expression-matched sets (`pseudo_waveforms(kind="matched")`, seed [332, 4]),
each through the same adjustment and V. Reported: log(V_clock / median V_matched), matched q95 of log(V_set / median),
clock percentile among the sets; C (`common.development.coherence`, 2,000 permutations, seed 332) and p.
`pass_rate`: adult pass rate of the same bin in the A237 / A239 feasibility tables (none has these bins: blank);
`pass_rate_nearest` / `pass_label`: for Cardoso-Moreira prenatal bins the nearest design (4-8wpc, 13-20wpc; other
n, other age span), for context only.

    python d4_development.py [--check]    (--check: 200 permutations, 10 sets, `_check` file)

Output: `output/d4_development.csv`, one row per dataset x organ x bin (n libraries or brains, n_individuals).
"""
import argparse
import importlib.util
import re
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve()
EXPL = HERE.parents[1]
sys.path.insert(0, str(HERE.parents[2]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from clocksync.config import limit_threads  # noqa: E402
from clocksync.paired import projector, pseudo_waveforms, states  # noqa: E402
from common import development as dev  # noqa: E402
from common.paths import DATA  # noqa: E402


def _load_module(name, rel):
    spec = importlib.util.spec_from_file_location(name, EXPL / rel)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


A237 = _load_module("score_a237", "A237_fetal_organs/scripts/score_a237.py")
A239 = _load_module("score_a239", "A239_fetal_brain/scripts/score_a239.py")
OUT = HERE.parent / "output"
BINS = ["prenatal<12wpc", "prenatal>=12wpc", "birth-12y", ">=13y"]
FEAS = EXPL / "scouts_2026-11-10/feasibility_fetal/output/passrates_cm_day_sex.csv"
NEAREST = {"prenatal<12wpc": "4-8wpc", "prenatal>=12wpc": "13-20wpc"}
SEED, N_PERM, N_PSEUDO, MIN_N = 332, 2000, 200, 8


def assign(age_y_postnatal, wpc):
    """Bin of each row: `wpc` (weeks post conception, nan if postnatal), `age_y_postnatal` (years, nan if prenatal)."""
    return np.where(wpc.notna(), np.where(wpc < 12, BINS[0], BINS[1]),
                    np.where(age_y_postnatal < 13, BINS[2], BINS[3]))


def cm_ages(obs):
    """Years after birth of the postnatal Cardoso-Moreira libraries, from `stage_detail`."""
    s = pd.read_csv(A237.SRC / "individuals.csv").set_index("sample").loc[obs.index]
    unit = {"day": 1 / 365.25, "month": 1 / 12, "year": 1.0}
    def years(t):
        m = re.match(r"(\d+) (day|month|year)s? post birth", str(t))
        return int(m[1]) * unit[m[2]] if m else np.nan
    return s.stage_detail.map(years).where(~s.prenatal)


def kang_ages(obs):
    """(years after birth, weeks post conception) of the Kang rows, from `age_value` and `age_unit`."""
    s = pd.read_csv(A239.SRC / "samples.csv").drop_duplicates("brain").set_index("brain")
    s = s.loc[obs.donor]
    wpc = pd.Series(np.where(s.age_unit == "PCW", s.age_value, np.nan), index=obs.index)
    y = pd.Series(np.where(s.age_unit == "Y", s.age_value, np.where(s.age_unit == "PMonth", s.age_value / 12, np.nan)),
                  index=obs.index)
    return y, wpc


def v_stat(counts, obs, cov, wf):
    """Mean |z|^2 of the centred clock-plane states of the rows `obs` (adjusted, one projector)."""
    genes = list(wf.genes)
    tab, _ = dev.adjust(dev.table(counts, obs, genes), genes, cov)
    L = projector(tab[genes], wf, counts=True)
    b = states(tab, L, wf, depth=tab.depth.to_numpy(float))[0]
    return float(np.mean(np.sum((b - b.mean(0)) ** 2, axis=1))), tab


TASKS = []  # set before the pool forks (closures and tables are not pickled)


def one(i):
    """Row of task `i` (one dataset x organ x bin; a worker)."""
    name, organ, bin_, counts, obs, cov, wf, n_perm, n_pseudo = TASKS[i]
    t0 = time.time()
    V, tab = v_stat(counts, obs, cov, wf)
    C, p, q95 = dev.coherence(tab, wf, n_perm, np.random.default_rng(SEED))
    lm = dev.log_mean(counts.loc[obs.index], obs.depth)
    vs = np.array([v_stat(counts, obs, cov, w)[0] for w in
                   pseudo_waveforms(lm, n_pseudo, np.random.default_rng([SEED, 4]), wf, "matched")])
    med = np.median(vs)
    return dict(dataset=name, organ=organ, bin=bin_, n=len(obs), n_individuals=obs.donor.nunique(),
                n_covariates=cov(obs).shape[1], V=V, median_V_matched=med, log_V_ratio=np.log(V / med),
                matched_q95_log=np.quantile(np.log(vs / med), .95), clock_percentile=100 * np.mean(vs < V),
                C=C, p=p, seconds=time.time() - t0)


def cm_tasks(n_perm, n_pseudo):
    counts, obs = A237.load()
    yrs = cm_ages(obs)
    obs = obs.assign(bin=assign(yrs, obs.age.where(obs.prenatal)), years=yrs)
    out = []
    for o in A237.ORGANS:
        for b in BINS:
            r = obs[(obs.compartment == o) & (obs.bin == b)]
            if r.donor.nunique() < MIN_N:
                continue
            if b.startswith("prenatal"):
                df = 3 if r.age.max() - r.age.min() > 8 else 2
                r = r.assign(a=r.age)
            else:
                df = 2 if b == "birth-12y" else 1
                r = r.assign(a=np.log2(r.years + .1))
            out.append(("CardosoMoreira2019", o, b, counts, r, _cov(df, "a", "sex"), A237.WF, n_perm, n_pseudo))
    return out


def kang_tasks(n_perm, n_pseudo):
    counts, obs = A239.load()
    obs = obs[obs.compartment == "brain_mean"]
    y, wpc = kang_ages(obs)
    obs = obs.assign(bin=assign(y, wpc))
    out = []
    for b in BINS:
        r = obs[obs.bin == b]
        if r.donor.nunique() >= MIN_N:
            out.append(("Kang2011", "brain", b, counts, r, _cov(2 if b == BINS[1] else 1, "age", "sex"), A239.WF,
                        n_perm, n_pseudo))
    return out


def _cov(df, age, factor):
    def f(g):
        return dev.design(g, age, df, (factor,))
    return f


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--jobs", type=int, default=4)
    a = ap.parse_args()
    limit_threads(1)
    n_perm, n_pseudo = (200, 10) if a.check else (N_PERM, N_PSEUDO)
    TASKS.extend(kang_tasks(n_perm, n_pseudo) + cm_tasks(n_perm, n_pseudo))
    tasks = range(len(TASKS))
    print(len(tasks), "organ x bin tests", flush=True)
    from multiprocessing import Pool
    rows, t0 = [], time.time()
    with Pool(a.jobs, initializer=limit_threads, initargs=(1,)) as pool:
        for i, r in enumerate(pool.imap(one, tasks), 1):
            rows.append(r)
            print(i, r["dataset"], r["organ"], r["bin"], f"n={r['n']} C={r['C']:.3f} "
                  f"logV={r['log_V_ratio']:.2f} {time.time() - t0:.0f}s", flush=True)
    d = pd.DataFrame(rows)
    f = pd.read_csv(FEAS)
    f = f[(f.adjust == "adjusted")].set_index("label").pass_rate
    d["pass_rate"] = np.nan
    d["pass_label"] = [f"{o}:{NEAREST[b]}" if r == "CardosoMoreira2019" and b in NEAREST else ""
                       for r, o, b in zip(d.dataset, d.organ, d.bin)]
    d["pass_rate_nearest"] = d.pass_label.map(f)
    OUT.mkdir(exist_ok=True)
    d.to_csv(OUT / f"d4_development{'_check' if a.check else ''}.csv", index=False)
    pd.set_option("display.width", 250)
    print(d.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
