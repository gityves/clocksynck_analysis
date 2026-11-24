"""A359 (registered): are louder organs the ones that read the body's systemic glucocorticoid state more strongly?

GTEx fast deaths, 25 peripheral tissues (A357 rows), 14-gene waveform, every tissue an anchor. Per tissue t: the
reception slope beta_t of the tissue's GC score on the person GC factor of the donor's other tissues, both
residualised on the leave-one-out host; the same for immediate-early genes. H359.1: Spearman(log G_t, beta_GC) > 0;
H359.2: rho_GC - rho_IEG > 0 (tissue bootstrap). Descriptive: adult NR3C1 level per tissue. Aggregates only.

    python a359.py
Output: ../output/a359_tissues.csv, a359_tests.csv.
"""
import sys
import warnings
from pathlib import Path

HERE = Path(__file__).resolve()
EXP = HERE.parents[2]
sys.path.insert(0, str(EXP.parent))
sys.path.insert(0, str(EXP.parent / "results" / "figure_2" / "scripts" / "backbone"))
sys.path.insert(0, str(EXP / "A329_core_departure" / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

import clocksync as cs  # noqa: E402
import score_a329 as A329  # noqa: E402
from clocksync._core import host as HC  # noqa: E402
from clocksync._core.controls import GC_TARGETS, IMMEDIATE_EARLY  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.genes import canonical  # noqa: E402

OUT = HERE.parents[1] / "output"
META = ["donor", "compartment", "depth"]
EXCLUDE = ["Whole Blood", "Testis"]
CLOCK_OUTPUTS = {"NFIL3", "BHLHE40"}
MIN_OTHER = 3


def load():
    pb = A329.load()
    pb = pb[~pb.compartment.isin(EXCLUDE)].reset_index(drop=True)
    wf = cs.gtex()
    clock = {canonical(g) for g in wf.genes}
    gc = [g for g in dict.fromkeys(GC_TARGETS) if g in pb.columns and g not in CLOCK_OUTPUTS and canonical(g) not in clock]
    ieg = [g for g in dict.fromkeys(IMMEDIATE_EARLY) if g in pb.columns and canonical(g) not in clock]
    return pb, wf, gc, ieg


def sample_scores(pb, d, genes):
    """Per Dataset row: mean log fraction over `genes`, centred within tissue x intercept group; and its within-tissue
    standardised version."""
    key = pd.Series(np.arange(len(pb)), index=pd.MultiIndex.from_arrays([pb.donor.astype(str), pb.compartment]))
    rows = key.loc[list(zip(d.resection.astype(str), d.compartment))].to_numpy()
    y = np.log((pb[genes].to_numpy(float)[rows] + 0.5) / pb.depth.to_numpy(float)[rows, None]).mean(1)
    grp = pd.Series(d.compartment.astype(str)) + "|" + pd.Series(d.group.astype(str))
    y = y - pd.Series(y).groupby(grp.to_numpy()).transform("mean").to_numpy()
    z = pd.Series(y).groupby(d.compartment).transform(lambda v: (v - v.mean()) / v.std(ddof=1)).to_numpy()
    idx = pd.MultiIndex.from_arrays([d.resection.astype(str), d.compartment])
    return pd.Series(y, index=idx), pd.Series(z, index=idx)


def reception(d, comps, y, z, hosts):
    out = {}
    for t in comps:
        h = hosts[t]
        H = pd.DataFrame(h.f, index=[str(r) for r in h.resections], columns=["hx", "hy"])
        yt = y[y.index.get_level_values(1) == t].droplevel(1)
        zo = z[z.index.get_level_values(1) != t]
        cnt = zo.groupby(level=0).size()
        x = zo.groupby(level=0).mean()[cnt >= MIN_OTHER]
        sh = yt.index.intersection(x.index).intersection(H.index)
        X = np.column_stack([np.ones(len(sh)), H.loc[sh].to_numpy()])
        ry = yt[sh].to_numpy() - X @ np.linalg.lstsq(X, yt[sh].to_numpy(), rcond=None)[0]
        rx = x[sh].to_numpy() - X @ np.linalg.lstsq(X, x[sh].to_numpy(), rcond=None)[0]
        out[t] = dict(beta=float(ry @ rx / (rx @ rx)), r=float(np.corrcoef(ry, rx)[0, 1]), n=len(sh))
    return pd.DataFrame(out).T


def perm_p(x, y, rng, n=10_000):
    r0 = spearmanr(x, y)[0]
    rs = np.array([spearmanr(x, rng.permutation(y))[0] for _ in range(n)])
    return float(r0), float((1 + np.sum(rs >= r0)) / (n + 1))


def main():
    warnings.filterwarnings("ignore")
    limit_threads(4)
    OUT.mkdir(parents=True, exist_ok=True)
    pb, wf, gc, ieg = load()
    print(f"GC genes {gc}\nIEG genes {ieg}", flush=True)
    comps = sorted(pb.compartment.unique())
    ds = cs.Dataset(pb[META + list(wf.genes)], anchors=comps, waveform=wf, name="a359")
    d = ds.data
    hosts = {t: HC.host(d, [c for c in comps if c != t]) for t in comps}
    tab = {}
    for lab, genes in (("gc", gc), ("ieg", ieg)):
        y, z = sample_scores(pb, d, genes)
        tab[lab] = reception(d, comps, y, z, hosts).add_prefix(f"{lab}_")
    t = tab["gc"].join(tab["ieg"])
    a = pd.read_csv(EXP / "A294_floor/output/tissues.csv", index_col="tissue")
    t["log_gain"] = np.log(np.sqrt(a.S.reindex(t.index)))
    t["D"] = a.own.reindex(t.index)
    cpm = np.log1p(pb["NR3C1"].to_numpy(float) / pb.depth.to_numpy(float) * 1e6)
    t["nr3c1_log_cpm"] = pd.Series(cpm).groupby(pb.compartment.to_numpy()).median().reindex(t.index)
    t = t.astype(float)
    t.to_csv(OUT / "a359_tissues.csv")
    rng = np.random.default_rng(359)
    r_gc, p_gc = perm_p(t.log_gain.to_numpy(), t.gc_beta.to_numpy(), rng)
    r_ieg, p_ieg = perm_p(t.log_gain.to_numpy(), t.ieg_beta.to_numpy(), rng)
    boot = []
    for _ in range(2000):
        i = rng.integers(0, len(t), len(t))
        s = t.iloc[i]
        boot.append(spearmanr(s.log_gain, s.gc_beta)[0] - spearmanr(s.log_gain, s.ieg_beta)[0])
    lo, hi = np.nanquantile(boot, [0.05, 0.95])
    tests = pd.DataFrame([
        dict(test="H359.1", rho=r_gc, p_one_sided=p_gc, passes=bool(p_gc < 0.05 and r_gc > 0), tissues=len(t)),
        dict(test="H359.2", rho=r_ieg, p_one_sided=p_ieg, diff=r_gc - r_ieg, diff_lo90=lo, diff_hi90=hi,
             specific=bool(lo > 0)),
        dict(test="descriptive NR3C1 vs gain", rho=spearmanr(t.log_gain, t.nr3c1_log_cpm)[0]),
        dict(test="descriptive gc_beta vs D", rho=spearmanr(t.D, t.gc_beta)[0]),
        dict(test="descriptive ieg_beta vs D", rho=spearmanr(t.D, t.ieg_beta)[0])])
    tests.to_csv(OUT / "a359_tests.csv", index=False)
    print(t.round(3).sort_values("log_gain").to_string())
    print(tests.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
