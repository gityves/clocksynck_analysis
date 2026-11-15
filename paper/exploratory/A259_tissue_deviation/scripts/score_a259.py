"""A259: how far does a tissue's clock depart from the clock state its person's other tissues share, in hours, net of
measurement noise? (GTEx fast deaths, 11 body tissues; MuTHER fat and skin)

Registered in paper/docs/registered_experiments.md (A259; power note in its addendum). Estimator
`clocksync.deviation.incoherence` (20 gene splits, 20 host splits, 1,000 donor bootstraps, seed 259); companions: 100
matched and 100 stress-responsive pseudo-clock waveforms on the same gene splits. GTEx: aggregates only.

    python score_a259.py [--jobs 4] [--placebo [--smoke]]
Outputs (output/score/): targets.csv (GTEx, one row per target), living.csv (MuTHER and the GTEx fat-skin pairs),
companions.csv (per target and family: kappa_tech q50 / q95, cross-tissue share), tests.csv.

Choices not fixed by the registration
- GTEx table: `common.gtex.load` (fast deaths with time of death, depth = library size), all genes, so the pseudo sets
  can be read from the same rows; the log means for the pseudo sets are over all rows of the 11 tissues.
- kappa_tech of a pseudo set = mean over the 20 gene splits (same positions as the clock genes) of its target's
  within-sample moment |m(beta^A, beta^B)|, divided by the clock's mean within-sample moment over 1 + kappa_k.
- Cross-tissue share of a pseudo set: median over host tissues l of |m(beta_k, beta_l)| (pseudo) / (clock).
- MuTHER arrays: the deposited log2 values (`A254_paired_lesions/output/prep/`), counts=False; pseudo sets from the mean
  over the women of the mean of fat and skin.
"""

from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[2]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(PAPER / "results/figure_3/scripts/gtex_cues"))
import clocksync as cs  # noqa: E402
from clocksync import deviation as DV  # noqa: E402
from clocksync import paired as P  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.genes import canonical  # noqa: E402
from common import gtex  # noqa: E402
import _shared as S  # noqa: E402

OUT = HERE.parent / "output" / "score"
SEED, NS, NH, NB, NP = 259, 20, 20, 1000, 100
WF = cs.gtex()
TARGETS = [t for t in S.BODY if t != "Whole Blood"]
PREP = PAPER / "exploratory/A254_paired_lesions/output/prep"


def _m(x, y):
    return np.mean((x - x.mean()) * np.conj(y - y.mean()))


def within(table, target, wf, splits, counts):
    """Mean over splits (gene positions of the clock) of |m(beta^A, beta^B)| of the target."""
    t = table[table.compartment == target]
    v = []
    for a, b in splits:
        ia, ib = [WF.genes.index(g) for g in a], [WF.genes.index(g) for g in b]
        A = DV.states(t, wf.subset([wf.genes[i] for i in ia]), counts).to_numpy()
        B = DV.states(t, wf.subset([wf.genes[i] for i in ib]), counts).to_numpy()
        v.append(abs(np.mean(A * np.conj(B))))
    return float(np.mean(v))


def cross(table, target, hosts, wf, counts):
    z = DV.states(table, wf, counts)
    Z = z.map(np.real).astype(float).unstack("compartment") + 1j * z.map(np.imag).astype(float).unstack("compartment")
    out = []
    for h in hosts:
        ok = Z[target].notna() & Z[h].notna()
        if ok.sum() >= 20:
            out.append(abs(_m(Z.loc[ok, target].to_numpy(), Z.loc[ok, h].to_numpy())))
    return np.array(out)


def one_target(table, target, hosts, wf_sets, counts):
    limit_threads(1)
    warnings.simplefilter("ignore")
    r = DV.incoherence(table, target, hosts, waveform=WF, counts=counts, n_splits=NS, n_host_splits=NH, n_boot=NB,
                       seed=SEED)
    splits = DV.split_genes(WF, NS, np.random.default_rng(SEED))
    wc = within(table, target, WF, splits, counts)
    clock_var = wc / (1 + r.kappa)
    cc = cross(table, target, hosts, WF, counts)
    rows = []
    for fam, sets in wf_sets.items():
        kt, sh = [], []
        for wf in sets:
            kt.append(within(table, target, wf, splits, counts) / clock_var)
            cp = cross(table, target, hosts, wf, counts)
            sh.append(float(np.median(cp / cc)) if len(cp) == len(cc) and len(cc) else np.nan)
        rows.append(dict(target=target, family=fam, kappa_tech_q50=float(np.median(kt)),
                         kappa_tech_q95=float(np.quantile(kt, 0.95)), cross_share_q50=float(np.nanmedian(sh)),
                         cross_share_q95=float(np.nanquantile(sh, 0.95)), sets=len(sets)))
    return r, pd.DataFrame(rows)


def gtex_table():
    pb, _ = gtex.load(tissues=TARGETS, hardy=(1, 2), genes="all", depth="library_size")
    meta = ["donor", "compartment", "depth"]
    genes = [c for c in pb.columns if c not in meta]
    can = pd.Index([canonical(g) for g in genes])
    keep = ~can.duplicated()
    pb = pd.concat([pb[meta], pb[genes].loc[:, keep].set_axis(can[keep], axis=1)], axis=1)
    x = pb[can[keep]].to_numpy(float)
    lm = pd.Series(np.log((x + 0.5) / pb.depth.to_numpy(float)[:, None]).mean(0), index=can[keep])
    return pb, lm


def muther_table():
    t = {k: pd.read_parquet(PREP / f"muther_{k}_expr.parquet").astype(float) for k in ("fat", "skin")}
    common = sorted(set(t["fat"].columns) & set(t["skin"].columns))
    rows = []
    for k, x in t.items():
        y = x[common].copy()
        y.insert(0, "compartment", k)
        y.insert(0, "donor", x.index.astype(str))
        rows.append(y.reset_index(drop=True))
    tab = pd.concat(rows, ignore_index=True)
    lm = (t["fat"][common].mean() + t["skin"][common].mean()) / 2 * np.log(2.0)
    lm.index = [canonical(g) for g in lm.index]
    tab.columns = [c if c in ("donor", "compartment") else canonical(c) for c in tab.columns]
    tab = tab.loc[:, ~tab.columns.duplicated()]
    return tab, lm[~lm.index.duplicated()]


def families(lm, wf, rng):
    lm = lm.dropna()
    return {"matched": P.pseudo_waveforms(lm, NP, rng, wf, "matched"),
            "responsive": P.pseudo_waveforms(lm, NP, rng, wf, "responsive")}


def shuffle(table, rng):
    t = table.copy()
    for c in t.compartment.unique():
        i = np.flatnonzero((t.compartment == c).to_numpy())
        t.iloc[i, t.columns.get_loc("donor")] = t.donor.to_numpy()[rng.permutation(i)]
    return t


def restrict(table, fam):
    """Only the clock genes and the genes of the pseudo sets (the states read nothing else)."""
    need = list(dict.fromkeys([*WF.genes, *(g for sets in fam.values() for w in sets for g in w.genes)]))
    return table[[c for c in table.columns if c in ("donor", "compartment", "depth")] + need].copy()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--placebo", action="store_true", help="donors shuffled within tissue (code test; output/placebo/)")
    ap.add_argument("--smoke", action="store_true", help="few splits, bootstraps and sets (with --placebo only)")
    a = ap.parse_args()
    global OUT, NS, NH, NB, NP
    if a.placebo:
        OUT = HERE.parent / "output" / "placebo"
    if a.smoke:
        assert a.placebo, "--smoke only with --placebo"
        NS, NH, NB, NP = 4, 4, 20, 5
    limit_threads(1)
    warnings.simplefilter("ignore")
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    pb, lm = gtex_table()
    if a.placebo:
        pb = shuffle(pb, rng)
    n = pb.groupby("compartment").donor.nunique()
    print("GTEx rows", len(pb), n.to_dict(), flush=True)
    fam = families(lm, WF, rng)
    pb = restrict(pb, fam)
    res = Parallel(n_jobs=a.jobs)(delayed(one_target)(pb, t, [h for h in TARGETS if h != t], fam, True)
                                  for t in TARGETS)
    T = pd.DataFrame([r for r, _ in res])
    C = pd.concat([c for _, c in res], ignore_index=True)
    T.to_csv(OUT / "targets.csv", index=False)
    print(T.round(3).to_string(index=False), flush=True)

    # living pair and the GTEx fat-skin pairs (pair mode)
    mt, mlm = muther_table()
    if a.placebo:
        mt = shuffle(mt, rng)
    mfam = families(mlm, WF, rng)
    mt = restrict(mt, mfam)
    rm, cm = one_target(mt, "fat", ["skin"], mfam, False)
    live = [dict(cohort="MuTHER", pair="fat-skin", **rm.to_dict())]
    for sk in ("Skin - Sun Exposed (Lower leg)", "Skin - Not Sun Exposed (Suprapubic)"):
        sub = pb[pb.compartment.isin(["Adipose - Subcutaneous", sk])]
        r = DV.incoherence(sub, "Adipose - Subcutaneous", [sk], waveform=WF, n_splits=NS, n_boot=NB, seed=SEED)
        live.append(dict(cohort="GTEx fast deaths", pair=f"Adipose - Subcutaneous vs {sk}", **r.to_dict()))
    L = pd.DataFrame(live)
    L.to_csv(OUT / "living.csv", index=False)
    C = pd.concat([C, cm.assign(target="MuTHER fat (pair)")], ignore_index=True)
    C.to_csv(OUT / "companions.csv", index=False)

    # tests
    tests = []
    elig = T[T.donors >= 80]
    med_hi = float(np.median(elig.hours_hi))
    tests.append(dict(test="H259.1", readout="median over targets of the upper 90 % bound of sigma (h)", value=med_hi,
                      targets=len(elig), below_1h=int((elig.hours_hi < 1.0).sum()), passed=bool(med_hi < 1.0)))
    holm = cs.stats.holm(elig.p_le_0.to_numpy())
    for (_, r), ph in zip(elig.iterrows(), holm):
        cq = C[C.target == r.target].set_index("family")
        spec = bool(r.kappa > cq.loc["matched", "kappa_tech_q95"] and r.kappa > cq.loc["responsive", "kappa_tech_q95"])
        tests.append(dict(test="H259.2", readout=r.target, value=r.kappa, hours=r.hours, hours_hi=r.hours_hi,
                          p=r.p_le_0, p_holm=ph, above_companions=spec, passed=bool(ph < 0.05 and spec)))
    tests.append(dict(test="H259.3", readout="MuTHER fat-skin upper 90 % bound of sigma_pair (h)",
                      value=float(rm.hours_hi), kappa=rm.kappa, hours=rm.hours, passed=bool(rm.hours_hi < 1.41)))
    out = pd.DataFrame(tests)
    out.to_csv(OUT / "tests.csv", index=False)
    pd.set_option("display.width", 250)
    print(out.round(3).to_string(index=False), "\n", L.round(3).to_string(index=False), "\n", C.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
