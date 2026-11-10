"""Feasibility of the bulk human fetal designs (before any registration): adult material cut to the fetal design.

No fetal expression value is read: only the sample sheets of Cardoso-Moreira 2019 (E-MTAB-6814) and Kang 2011
(GSE25219) (individual, organ or region, age). GTEx donors are cut to each design and the share of subsamples where
the planned readout passes is its sensitivity for an adult-like clock (feasible if >= 0.8, the A181 / D1b rule).

Readouts, as they would be registered:
1. L1 per organ: signed clock-gene coherence C (`clocksync.decompose.signed_coherence`: a shift common to the clock
   genes of a sample is fitted out and the sign kept) across individuals, each gene residualised on a polynomial of
   age, against per-gene permutations of the residuals; pass = one-sided p < 0.05. `--stat unsigned`: the coherence R
   of `decompose.components` (the A189 statistic; the first run of this script, `passrates_cm_day.csv`).
2. L2 across organs (or brain regions) of one individual: `Dataset.gate` on the age-residualised log fractions, free
   host, age bins as strata (covariances centred and pairing permuted within a bin); pass = p_pairing < 0.05 and
   p_acrophase < 0.05.
The adult has no developmental age: each adult donor receives the age (and age bin) of the fetal individual it stands
for, so the adjustment costs the same degrees of freedom as in the fetal analysis (`unadjusted` rows: without it).
`--sex` (output suffix `_sex`; the registered A237 design): the adult donor also receives the sex of that individual
as a covariate, the postnatal L1 is adjusted on a quadratic of the stage rank and sex, and a linked library whose sex
differs from the majority sex of its individual counts as unlinked (one liver library).

Designs
- cm (Cardoso-Moreira): L1 per organ at the prenatal n of each age range and at the postnatal n; L2 at the organ
  sets of the individuals with >= 2 linked organs (each fetal individual is stood for by one GTEx donor having its
  organs; `organs_all` counts the unlinked columns as if their guess were right: an upper bound).
- kang: L2 across the brain regions of a donor (GTEx donors with >= 6 of the 13 regions) at the number of brains of
  each age group; L1 on the mean over regions of a brain.
Windows: `day` = GTEx deaths at 08-18 h (terminations are done in hospital hours), `all` = any hour (post-mortem
brains; fetal clocks not entrained to clock time).

Controlled access: GTEx time of death is used in memory to select donors; only aggregates are written.

Kang 2011 is an exon array whose gene-level table lacks PER3 and NR1D1 (`--drop PER3,NR1D1`); array noise is not
modelled, so the postnatal brains of the same platform are the power check of that dataset, not this calibration.

    python feasibility_fetal.py --design cm|kang [--window day|all] [--n-sub 200] [--jobs 4] [--drop G1,G2]
                                [--stat signed|unsigned] [--only L1|L2] [--sex] [--check]
"""
import argparse
import sys
import time
from multiprocessing import Pool
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[4]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import binomtest  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.decompose import components, signed_coherence  # noqa: E402
from common import gtex  # noqa: E402
from common.benchmark.data import in_window  # noqa: E402

OUT = HERE.parents[1] / "output"
NATEXP = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian/natexp")
SEED, N_PERM_L1, N_PERM_GATE, DEPTH = 1110, 500, 500, 1e9
STAT = "signed"  # L1 statistic (`--stat`)
WF = cs.gtex()
GENES = list(WF.genes)
ORGAN = {"heart": "Heart - Left Ventricle", "liver": "Liver", "forebrain": "Brain - Frontal Cortex (BA9)",
         "hindbrain": "Brain - Cerebellum", "kidney": "Kidney - Cortex", "testis": "Testis", "ovary": "Ovary"}
CM_RANGES = {"4-20": (4, 20), "9-20": (9, 20), "13-20": (13, 20), "4-8": (4, 8)}  # weeks post conception
KANG_GROUPS = {"prenatal": (0, 7), "pcw13-22": (4, 6), "pcw<=12": (1, 3), "infant+child": (8, 11), "adult": (12, 15),
               "postnatal": (8, 15)}  # Kang stages (1-7 prenatal, 8-9 first year, 10-11 to 11 y, 12-15 from 13 y)
MIN_REGIONS, POOL_FACTOR = 6, 1.5
KANG_REGIONS = ["A1C", "AMY", "CBC", "DFC", "HIP", "IPC", "ITC", "M1C", "MD", "MFC", "OFC", "S1C", "STC", "STR", "V1C",
                "VFC"]  # the regions sampled from stage 3 onwards (earlier brains have other regions, <= 2 brains each)
CM_BINS, KANG_BINS = [8.5, 12.5], [3.5, 6.5, 7.5, 9.5, 11.5]  # strata: wpc 4-8 | 9-12 | 13-20; Kang stages 1-3 | 4-6 | 7 | 8-9 | 10-11 | 12-15
POSTNATAL = ["newborn", "infant", "toddler", "school", "youngTeenager", "teenager", "oldTeenager", "youngAdult",
             "youngMidAge", "olderMidAge", "senior"]  # Cardoso-Moreira stages after birth, in order
JOBS, POOL = [], {}  # (kind, label, n, ages, sets or None, df, strata, sex or None); donor -> table
SEX = False  # `--sex`


def basis(age, df):
    """Centred polynomial of age (df columns); empty when df = 0 or age is constant."""
    a = np.asarray(age, float)
    if df == 0 or np.ptp(a) == 0:
        return np.empty((len(a), 0))
    a = (a - a.mean()) / a.std()
    return np.column_stack([a ** k for k in range(1, df + 1)])


def covariates(age, df, sex=None):
    """Design of the adjustment: polynomial of age and, when given and not constant, the sex indicator; empty when
    df = 0 (unadjusted)."""
    Z = basis(age, df)
    if df and sex is not None and np.ptp(sex) > 0:
        Z = np.column_stack([Z, np.asarray(sex, float)])
    return Z


def residual(y, Z):
    """y minus its regression on Z (intercept included; the mean is kept)."""
    M = np.column_stack([np.ones(len(y)), Z])
    return y - M @ np.linalg.lstsq(M, y, rcond=None)[0] + y.mean(0)


def logfrac(tab):
    return np.log((tab[GENES].to_numpy(float) + 0.5) / tab.depth.to_numpy(float)[:, None])


def l1(y, rng, n_perm=N_PERM_L1):
    """L1 statistic of log values y (signed coherence C, or the unsigned R with `--stat unsigned`) and its per-gene
    permutation p."""
    y = y - y.mean(0)
    stat = (lambda x: signed_coherence(x, WF)) if STAT == "signed" else (
        lambda x: components(x, np.zeros(x.shape[1]), 1, WF)["coherence"])
    obs = stat(y)
    null = np.array([stat(np.column_stack([rng.permutation(c) for c in y.T])) for _ in range(n_perm)])
    return obs, (1 + np.sum(null >= obs)) / (1 + n_perm)


def adjusted(tab, age, df, sex=None):
    """Table whose log fractions are residualised on age (and sex) within each compartment (counts at a common
    depth)."""
    out = []
    for _, g in tab.groupby("compartment", sort=False):
        g = g.copy()
        r = residual(logfrac(g), covariates(g.donor.map(age).to_numpy(float), df,
                                            None if sex is None else g.donor.map(sex).to_numpy(float)))
        g[GENES] = np.exp(r) * DEPTH
        g["depth"] = DEPTH
        out.append(g)
    return pd.concat(out, ignore_index=True)


def gate(tab, seed):
    keep = tab.groupby("compartment").donor.transform("nunique") >= 5
    g = cs.Dataset(tab[keep], anchors=None, waveform=WF, min_donors=5, name="fetal").gate(
        n_perm=N_PERM_GATE, seed=seed)
    return g.clock_share, g.p_pairing, g.p_acrophase


def draw(sets, rng, pool=None):
    """One adult donor of the pool (default: GTEx) per fetal individual, having all its compartments (rarest sets
    first; an individual whose set no remaining donor has is dropped). Returns [(individual index, donor)]."""
    have = {d: set(t.compartment) for d, t in (POOL if pool is None else pool).items()}
    free, out = set(have), []
    for i in sorted(range(len(sets)), key=lambda i: sum(sets[i] <= have[d] for d in have)):
        ok = sorted(d for d in free if sets[i] <= have[d])
        if ok:
            d = ok[rng.integers(len(ok))]
            free.discard(d)
            out.append((i, d))
    return out


def one(args):
    j, i = args
    kind, label, n, ages, sets, df, strata, sex = JOBS[j]
    rng = np.random.default_rng([SEED, j, i])
    rows = []
    if kind == "L1":  # one compartment (or the mean over a donor's compartments), n donors
        donors = rng.choice(sorted(d for d, t in POOL.items() if t.compartment.isin(sets).any()), size=n, replace=False)
        y = np.vstack([logfrac(POOL[d][POOL[d].compartment.isin(sets)]).mean(0) for d in donors])
        for tag, d in (("adjusted", df), ("unadjusted", 0)):
            r, p = l1(residual(y, covariates(ages, d, sex)), rng)
            rows.append(dict(kind=kind, label=label, sub=i, adjust=tag, n=n, stat=r, p=p, passed=p < .05))
        return rows
    if sets is None:  # all compartments of n donors
        pick = [(k, d) for k, d in enumerate(rng.choice(sorted(POOL), size=n, replace=False))]
        tab = pd.concat([POOL[d] for _, d in pick], ignore_index=True)
    else:
        pick = draw(sets, rng)
        tab = pd.concat([POOL[d][POOL[d].compartment.isin(sets[k])] for k, d in pick], ignore_index=True)
    age, strat = {d: ages[k] for k, d in pick}, {d: f"s{strata[k]}" for k, d in pick}
    sx = None if sex is None else {d: sex[k] for k, d in pick}
    for tag, d in (("adjusted", df), ("unadjusted", 0)):
        try:
            t = adjusted(tab, age, d, sx)
            t["stratum"] = t.donor.map(strat) if d else "one"
            share, pp, pa = gate(t, SEED + i)
        except ValueError:  # no compartment pair shares >= 5 individuals: the gate is not estimable at this design
            share, pp, pa = np.nan, np.nan, np.nan
        rows.append(dict(kind=kind, label=label, sub=i, adjust=tag, n=len(pick), stat=share, p=pp, p_acrophase=pa,
                         passed=bool(pp < .05 and pa < .05)))
    return rows


def pool(tissues, window):
    pb, d = gtex.load(tissues, hardy=(0, 1, 2, 3, 4), depth="library_size")
    if window == "day":
        pb = pb[pb.donor.isin(set(d.index[in_window(d.tod.to_numpy())]))]
    return pb


def add(kind, label, n, ages, sets, df, design, pool_n, strata=None, sex=None):
    ok = n >= 5 and pool_n >= POOL_FACTOR * n
    if ok:
        JOBS.append((kind, label, n, np.asarray(ages, float), sets, df,
                     np.zeros(n, int) if strata is None else np.asarray(strata),
                     np.asarray(sex == "male", float) if SEX and sex is not None else None))
    design.append(dict(kind=kind, label=label, fetal_n=n, age_df=df, adult_pool=pool_n, calibrated=ok))


def setup_cm(window):
    """Columns as linked by `prep_fetal_bulk/scripts/link_cardoso_individuals.py`: organ and age are exact for every
    column (L1); the individual is known for the linked columns (L2 `organs`), the `open` ones stand at a
    guess within their organ, stage and sex (L2 `organs_all`: an upper bound, the guesses are counted as right)."""
    s = pd.read_csv(NATEXP / "cardoso_moreira2019/individuals.csv")
    if SEX:  # an individual has one sex: a linked library of the other sex is not linked
        major = s[s.prenatal & (s.link != "open")].groupby("individual").sex.agg(lambda x: x.mode().iloc[0])
        s.loc[s.prenatal & (s.link != "open") & (s.sex != s.individual.map(major)), "link"] = "open"
    pb = pool(list(ORGAN.values()), window)
    back = {v: k for k, v in ORGAN.items()}
    pb["compartment"] = pb.compartment.map(back)
    POOL.update({d: g.reset_index(drop=True) for d, g in pb.groupby("donor")})
    n_pool = pb.groupby("compartment").donor.nunique()
    design, pre = [], s[s.prenatal]
    for o in ORGAN:
        for tag, (lo, hi) in CM_RANGES.items():
            f = pre[(pre.organ == o) & pre.age_wpc.between(lo, hi)]
            add("L1", f"{o}:{tag}wpc", len(f), f.age_wpc, [o], 3 if hi - lo > 8 else 2, design, n_pool.get(o, 0),
                sex=f.sex)
        f = s[~s.prenatal & (s.organ == o)]
        rank = f.stage.str.lower().map({k.lower(): i for i, k in enumerate(POSTNATAL)})
        add("L1", f"{o}:postnatal", len(f), rank if SEX else np.zeros(len(f)), [o], 2 if SEX else 0, design,
            n_pool.get(o, 0), sex=f.sex)
    have = [set(t.compartment) for t in POOL.values()]
    for name, cols in (("organs", pre[pre.link != "open"]), ("organs_all", pre)):
        ind = cols.groupby("individual").agg(organs=("organ", lambda x: frozenset(x)), age=("age_wpc", "first"),
                                             sex=("sex", lambda x: x.mode().iloc[0]))
        ind = ind[ind.organs.map(len) >= 2]
        for tag, (lo, hi) in CM_RANGES.items():
            f = ind[ind.age.between(lo, hi)]
            sets = [set(x) for x in f.organs]
            cover = sum(any(x <= h for h in have) for x in sets)
            add("L2", f"{name}:{tag}wpc", len(f), f.age, sets, 3 if hi - lo > 8 else 2, design,
                len(POOL) if cover == len(sets) else 0, np.digitize(f.age, CM_BINS), sex=f.sex)
            design[-1]["sets_with_a_donor"] = cover
    return pd.DataFrame(design)


def setup_kang(window):
    s = pd.read_csv(NATEXP / "kang2011/samples.csv")
    s = s[s.region.isin(KANG_REGIONS)]
    b = s.groupby("brain").agg(stage=("stage", "first"), n_regions=("region", "nunique"),
                               age=("age_log", "first"))
    pb = pool(None, window)
    pb = pb[pb.compartment.str.startswith("Brain")]
    pb = pb[pb.groupby("donor").compartment.transform("nunique") >= MIN_REGIONS]
    POOL.update({d: g.reset_index(drop=True) for d, g in pb.groupby("donor")})
    regions, design = sorted(pb.compartment.unique()), []
    for tag, (lo, hi) in KANG_GROUPS.items():
        f = b[b.stage.between(lo, hi) & (b.n_regions >= MIN_REGIONS)]
        df = 2 if tag in ("prenatal", "postnatal") else 1
        add("L2", f"regions:{tag}", len(f), f.age, None, df, design, len(POOL), np.digitize(f.stage, KANG_BINS))
        add("L1", f"brain_mean:{tag}", len(f), f.age, regions, df, design, len(POOL))
    return pd.DataFrame(design)


def ci(k, n):
    r = binomtest(int(k), int(n)).proportion_ci(0.95, method="exact")
    return r.low, r.high


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--design", choices=("cm", "kang"), required=True)
    ap.add_argument("--window", choices=("day", "all"), default="all")
    ap.add_argument("--n-sub", type=int, default=200)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--drop", default="", help="clock genes absent from the fetal platform, comma-separated")
    ap.add_argument("--stat", choices=("signed", "unsigned"), default="signed", help="L1 statistic")
    ap.add_argument("--only", choices=("L1", "L2"), help="run one kind of readout (output suffix _L1 / _L2)")
    ap.add_argument("--sex", action="store_true", help="sex as covariate, adjusted postnatal L1, sex-consistent links")
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    global WF, GENES, STAT, SEX
    STAT, SEX = a.stat, a.sex
    if a.drop:
        WF = WF.drop(a.drop.split(","))
        GENES = list(WF.genes)
    OUT.mkdir(exist_ok=True)
    t0 = time.time()
    limit_threads(1)
    sfx = f"_{a.design}_{a.window}" + (f"_{a.only}" if a.only else "") + ("_unsigned" if a.stat == "unsigned" else "") + (
        "_sex" if a.sex else "") + ("_check" if a.check else "")
    print(f"waveform genes ({len(GENES)}): {', '.join(GENES)}", flush=True)
    design = (setup_cm if a.design == "cm" else setup_kang)(a.window)
    design.to_csv(OUT / f"design{sfx}.csv", index=False)
    print(design.to_string(index=False), flush=True)
    n = 3 if a.check else a.n_sub
    tasks = [(j, i) for j in range(len(JOBS)) for i in range(n) if a.only in (None, JOBS[j][0])]
    print(f"{len(tasks) // n} calibration jobs x {n} subsamples, L1 statistic: {STAT}", flush=True)
    run(tasks, one, a.jobs, sfx, t0)


def run(tasks, one, jobs, sfx, t0):
    """Subsamples in parallel; tables `subsamples{sfx}.csv` and `passrates{sfx}.csv`."""
    res = []
    with Pool(jobs, initializer=limit_threads, initargs=(1,)) as p:
        for k, r in enumerate(p.imap_unordered(one, tasks, chunksize=2)):
            res.extend(r)
            if (k + 1) % 100 == 0:
                print(f"{k + 1}/{len(tasks)} {time.time() - t0:.0f}s", flush=True)
    res = pd.DataFrame(res)
    res.to_csv(OUT / f"subsamples{sfx}.csv", index=False)
    rows = []
    for (kind, label, adj), r in res.groupby(["kind", "label", "adjust"], sort=False):
        k, m = int(r.passed.sum()), len(r)
        lo, hi = ci(k, m)
        rows.append(dict(kind=kind, label=label, adjust=adj, n=int(r.n.median()), n_subsamples=m, pass_rate=k / m,
                         ci_lo=lo, ci_hi=hi, median_stat=r.stat.median(), median_p=r.p.median(),
                         estimable=r.p.notna().mean(), feasible=k / m >= 0.8))
    rates = pd.DataFrame(rows)
    rates.to_csv(OUT / f"passrates{sfx}.csv", index=False)
    print(rates.round(3).to_string(index=False))
    print(f"runtime {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
