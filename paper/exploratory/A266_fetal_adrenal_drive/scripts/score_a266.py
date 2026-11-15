"""A266 (registered in `paper/docs/registered_experiments.md`; NOT BLIND): is the clock-shaped covariation of the fetal adrenal
cortex an entrained clock loop or a hormone-driven response?

Stages
- `power` (before any fetal statistic): adult power of H266.2. From the GTEx adrenal glands (all Hardy classes, any hour),
  draws of N_FETAL donors; the phase relation of the responsive genes to the core clock state in the draw is compared
  with that of the other donors (`relation`, `agreement`), p against 200 permutations of the drawn donors' responsive
  rows; pass rate at alpha 0.05 over 200 draws. Also the false-positive rate (the drawn donors' responsive rows permuted
  once before the test). Writes power_h266_2.csv, adult_relation.csv (per gene, adult reference).
- `score`: H266.1 and descriptives (L1 tests of `design_a266`, `common.development.l1_test`, as A264; scored only if the
  adult pass rate at the full clock is >= 0.8); H266.2 on the fetal cortex; H266.3 receptor table. Writes l1.csv,
  verdicts.csv, relation_genes.csv, receptors.csv.

H266.2 statistic. Per sample, the clock state z (complex) from the CORE waveform (`clocksync.paired.projector` / `states`
on the adjusted counts, centred). Per responsive gene g, the centred log fraction regressed on (Re z, Im z): beta_g =
b1 + i b2, the gene's loading on the clock state. Agreement between two sets of samples: A = Re sum_g beta_g conj(beta'_g)
/ sqrt(sum |beta_g|^2 sum |beta'_g|^2). Fetal: cortex pseudobulks adjusted as A244 (cohort, quadratic age); p against
10,000 permutations of the fetuses' responsive rows (all genes together) within cohort.

    python score_a266.py --stage power|score [--jobs 6]
"""
import argparse
import sys
import time
from multiprocessing import Pool
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[3]))
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parents[2] / "A264_fetal_adrenal_covariates" / "scripts"))
sys.path.insert(0, str(HERE.parents[2] / "A244_fetal_adrenal_classes" / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import design_a266 as d  # noqa: E402
import score_a264 as S264  # noqa: E402
from clocksync._core.controls import RESPONSIVE_GENES  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.paired import projector, states  # noqa: E402
from common import development as dev  # noqa: E402
from common import gtex  # noqa: E402

OUT = d.OUT
N_FETAL, N_DRAW, N_PERM_ADULT, N_PERM = 34, 200, 200, 10_000
CORE = list(d.CORE.genes)
S = {}


def relation(counts: pd.DataFrame, depth, genes) -> np.ndarray:
    """beta_g (complex) of each responsive gene on the centred CORE clock state of the rows."""
    x = counts[CORE].assign(depth=np.asarray(depth, float))
    b, _ = states(x, projector(x, d.CORE), d.CORE, depth=x.depth.to_numpy())
    z = np.c_[b - b.mean(0)]
    y = np.log((counts[list(genes)].to_numpy(float) + 0.5) / np.asarray(depth, float)[:, None])
    y = y - y.mean(0)
    B = np.linalg.lstsq(z, y, rcond=None)[0]
    return B[0] + 1j * B[1]


def agreement(a, b) -> float:
    return float(np.real(np.sum(a * np.conj(b))) / np.sqrt(np.sum(np.abs(a) ** 2) * np.sum(np.abs(b) ** 2)))


def permuted(counts, genes, groups, rng):
    """Responsive rows permuted together within the row groups (the link to the clock state broken)."""
    c = counts.copy()
    v = c[list(genes)].to_numpy(dtype=float, copy=True)
    for i in groups:
        v[i] = v[rng.permutation(i)]
    c.loc[:, list(genes)] = v
    return c


def adult():
    import anndata as ad
    a = ad.read_h5ad(gtex.LABELLED, backed="r")
    have = set(a.var_names)
    a.file.close()
    want = [g for g in dict.fromkeys([*d.WF.genes, *RESPONSIVE_GENES_ON(), *d.RECEPTORS]) if g in have]
    pb, _ = gtex.load(["Adrenal Gland"], hardy=(0, 1, 2, 3, 4), genes=want, depth="library_size", require_tod=False)
    return pb.reset_index(drop=True)


def RESPONSIVE_GENES_ON():
    return [g for g in RESPONSIVE_GENES if g not in ("CYR61",)]  # CCN1 is the current symbol of CYR61


def _draw(i):
    pb, genes = S["pb"], S["genes"]
    rng = np.random.default_rng([d.SEED, 7, i])
    pick = rng.choice(len(pb), N_FETAL, replace=False)
    rest = np.setdiff1d(np.arange(len(pb)), pick)
    ref = relation(pb.iloc[rest], pb.depth.iloc[rest], genes)
    out = []
    for setting in ("adult", "shuffled"):
        x = pb.iloc[pick].reset_index(drop=True)
        if setting == "shuffled":
            x = permuted(x, genes, [np.arange(len(x))], rng)
        obs = agreement(relation(x, x.depth, genes), ref)
        null = [agreement(relation(permuted(x, genes, [np.arange(len(x))], rng), x.depth, genes), ref)
                for _ in range(N_PERM_ADULT)]
        out.append(dict(draw=i, setting=setting, A=obs, p=(1 + np.sum(np.array(null) >= obs)) / (1 + N_PERM_ADULT)))
    return out


def power(jobs):
    pb = adult()
    genes = [g for g in RESPONSIVE_GENES_ON() if g in pb and pb[g].mean() >= d.MIN_COUNT]
    S.update(pb=pb, genes=genes)
    ref = relation(pb, pb.depth, genes)
    pd.DataFrame(dict(gene=genes, beta_abs=np.abs(ref), beta_arg_h=np.angle(ref) * 24 / (2 * np.pi))).to_csv(
        OUT / "adult_relation.csv", index=False)
    with Pool(jobs, initializer=limit_threads, initargs=(1,)) as p:
        rows = [r for rr in p.imap_unordered(_draw, range(N_DRAW)) for r in rr]
    r = pd.DataFrame(rows)
    s = r.groupby("setting").agg(draws=("p", "size"), pass_rate=("p", lambda p: float((p < 0.05).mean())),
                                 median_A=("A", "median")).reset_index()
    s.insert(0, "adult_donors", len(pb))
    s.insert(1, "genes", len(genes))
    s.to_csv(OUT / "power_h266_2.csv", index=False)
    print(s.to_string(index=False), flush=True)


def task(name):
    rows, cov, strata, wf = S["tests"][name]
    return name, dev.l1_test(S["counts"], S["obs"][rows], cov, wf, N_PERM, 200, [d.SEED, 1, list(S["tests"]).index(name)],
                             strata)


def score(jobs):
    t0 = time.time()
    counts, obs = S264.load()
    tests = d.tests(obs)
    rates = pd.read_csv(OUT / "feasibility_passrates.csv")
    a1 = rates[rates.setting == "a1"].set_index("test").pass_rate
    run = [k for k in tests if a1.get(k, 0.0) >= 0.8]
    S.update(counts=counts, obs=obs, tests=tests)
    res = {}
    with Pool(jobs, initializer=limit_threads, initargs=(1,)) as p:
        for name, r in p.imap_unordered(task, run):
            res[name] = r
    l1 = pd.DataFrame([dict(test=k, adult_pass_a1=a1.get(k, np.nan), **res.get(k, {})) for k in tests])
    l1.to_csv(OUT / "l1.csv", index=False)
    # H266.2: fetal cortex, computed only if the adult power is >= 0.8 (registered rule)
    rows = d.d264.d244.pooled(obs, obs.compartment == "cortex")
    f = obs[rows]
    pw_rate = float(pd.read_csv(OUT / "power_h266_2.csv").set_index("setting").loc["adult", "pass_rate"])
    obs_A = p2 = np.nan
    if pw_rate >= 0.8:
        ref_tab = pd.read_csv(OUT / "adult_relation.csv").set_index("gene")
        genes = [g for g in ref_tab.index if counts.loc[f.index, g].mean() >= d.MIN_COUNT]
        ref = (ref_tab.beta_abs * np.exp(1j * ref_tab.beta_arg_h * 2 * np.pi / 24)).loc[genes].to_numpy()
        tab, _ = dev.adjust(dev.table(counts, f, [*CORE, *genes]), [*CORE, *genes], d.d264.cov())
        beta = relation(tab, tab.depth, genes)
        obs_A = agreement(beta, ref)
        rng = np.random.default_rng([d.SEED, 2])
        groups = list(tab.groupby("block", sort=False).indices.values())
        null = np.array([agreement(relation(permuted(tab, genes, groups, rng), tab.depth, genes), ref) for _ in range(N_PERM)])
        p2 = float((1 + np.sum(null >= obs_A)) / (1 + N_PERM))
        pd.DataFrame(dict(gene=genes, fetal_abs=np.abs(beta), fetal_arg_h=np.angle(beta) * 12 / np.pi,
                          adult_abs=np.abs(ref), adult_arg_h=np.angle(ref) * 12 / np.pi)).to_csv(OUT / "relation_genes.csv", index=False)
    # H266.3: receptors (log2 cpm medians by class; adult adrenal for reference)
    rec = []
    for k in ("cortex", "medulla", "stroma"):
        o = obs[d.d264.d244.pooled(obs, obs.compartment == k)]
        for g in d.RECEPTORS:
            rec.append(dict(group=f"fetal {k}", gene=g, n=len(o),
                            median_log2cpm=float(np.median(np.log2(counts.loc[o.index, g] / o.depth * 1e6 + 1)))))
    pb = adult()
    for g in d.RECEPTORS:
        rec.append(dict(group="adult adrenal (GTEx)", gene=g, n=len(pb), median_log2cpm=float(
            np.median(np.log2(pb[g] / pb.depth * 1e6 + 1))) if g in pb else np.nan))  # MTNR1A/B absent from the matrix
    pd.DataFrame(rec).to_csv(OUT / "receptors.csv", index=False)
    v = []
    for k in tests:
        r = res.get(k)
        verdict = ("not estimable" if r is None else "PASS" if r["p"] < 0.05 and r["clock_specific"] else
                   "not clock-specific" if r["p"] < 0.05 else "fail")
        v.append(dict(test=k, registered=k == "cortex:core", n=r and r["n"], stat=r and r["C"], p=r and r["p"],
                      matched_q95=r and r["matched_q95"], responsive_q95=r and r["responsive_q95"],
                      adult_power=a1.get(k, np.nan), verdict=verdict))
    v.append(dict(test="cortex:relation", registered=True, n=len(f), stat=obs_A, p=p2, matched_q95=np.nan,
                  responsive_q95=np.nan, adult_power=pw_rate,
                  verdict="not estimable" if pw_rate < 0.8 else "PASS" if p2 < 0.05 else "fail"))
    v = pd.DataFrame(v)
    v.to_csv(OUT / "verdicts.csv", index=False)
    pd.set_option("display.width", 250)
    print(l1.round(4).to_string(index=False), v.round(4).to_string(index=False), pd.DataFrame(rec).round(2).to_string(index=False),
          f"runtime {time.time() - t0:.0f}s", sep="\n", flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=["power", "score"], required=True)
    ap.add_argument("--jobs", type=int, default=6)
    a = ap.parse_args()
    limit_threads(1)
    power(a.jobs) if a.stage == "power" else score(a.jobs)
