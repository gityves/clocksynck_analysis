"""A329 (registered, check N2): is the personal departure the signal layer read through PER1, PER2 and NR1D1?

A294b rows and fit (GTEx fast deaths, 27 peripheral tissues, every tissue an anchor, leave-one-out host), on the
14-gene waveform (`full`) and on the 11-gene core without PER1, PER2, NR1D1 (`core`), each against the 100 matched sets
of A294b (seed [294, 1]); a set's 11-gene companion is its genes at the core positions. Per tissue and set: coupling,
own variance 2 psi, and the R^2 of the per-sample departure on the glucocorticoid-target and immediate-early scores.
Aggregates only (GTEx controlled access): no per-donor value is written.

    python score_a329.py [--jobs 4] [--check]
Outputs (`../output/`): sets.csv (per set and tissue), tissues.csv, tests.csv.
"""
import argparse
import sys
import time
import warnings
from multiprocessing import Pool
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[3]))
sys.path.insert(0, str(HERE.parents[3] / "results" / "figure_2" / "scripts" / "backbone"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import wilcoxon  # noqa: E402

import clocksync as cs  # noqa: E402
import compute_oscillator_gtex as A104  # noqa: E402
from clocksync._core.controls import GC_TARGETS, IMMEDIATE_EARLY  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.genes import canonical  # noqa: E402
from clocksync.paired import pseudo_waveforms  # noqa: E402
from clocksync.waveform import Waveform  # noqa: E402
from common import gtex  # noqa: E402

OUT = HERE.parents[1] / "output"
META = ["donor", "compartment", "depth"]
INDUCED = ("PER1", "PER2", "NR1D1")
A294B_MEDIAN = 0.969
S = {}


def load():
    import anndata as ad
    a = ad.read_h5ad(gtex.LABELLED, backed="r")
    genes = list(a.var_names)
    a.file.close()
    pb, _ = gtex.load(min_donors=A104.MIN_DONORS, genes=genes)
    return pb[~pb.compartment.str.startswith("Brain")].reset_index(drop=True)


def responsive_scores(pb):
    """Per sample: mean over the class genes of the log fraction standardised within tissue (14 clock genes excluded)."""
    clock = {canonical(g) for g in cs.gtex().genes}
    out = pd.DataFrame(index=pb.index)
    for name, genes in (("gc", GC_TARGETS), ("ieg", IMMEDIATE_EARLY)):
        g = [x for x in dict.fromkeys(genes) if x in pb.columns and canonical(x) not in clock]
        y = np.log((pb[g].to_numpy(float) + 0.5) / pb.depth.to_numpy(float)[:, None])
        y = pd.DataFrame(y, index=pb.index).groupby(pb.compartment.to_numpy()).transform(lambda v: (v - v.mean()) / v.std(ddof=1))
        out[name] = y.mean(axis=1)
    return out


def departure_r2(ds, comps, scores, pb):
    """Per tissue: R^2 of the per-sample departure (beta minus its 2 x 2 linear map of the leave-one-out host) on the
    two responsive scores."""
    d = ds.data
    key = pd.Series(np.arange(len(pb)), index=pd.MultiIndex.from_arrays([pb.donor.astype(str), pb.compartment]))
    res = {}
    for t in comps:
        i = np.flatnonzero(d.compartment == t)
        b = d.beta[i].copy()
        g = d.group[i]
        for u in np.unique(g):
            b[g == u] -= b[g == u].mean(0)
        donors = pd.Index(d.resection[i]).astype(str)
        h = ds.host([c for c in comps if c != t])
        h.index = h.index.astype(str)
        ok = donors.isin(h.index)
        b, donors = b[ok], donors[ok]
        H = h.loc[donors, ["x", "y"]].to_numpy()
        H = H - H.mean(0)
        M = np.linalg.lstsq(H, b, rcond=None)[0]
        dep = b - H @ M
        dep = dep - dep.mean(0)
        X = scores.iloc[key.loc[list(zip(donors, [t] * len(donors)))].to_numpy()].to_numpy()
        X = np.column_stack([np.ones(len(X)), X - X.mean(0)])
        fit = X @ np.linalg.lstsq(X, dep, rcond=None)[0]
        res[t] = 1 - ((dep - fit) ** 2).sum() / (dep ** 2).sum()
    return pd.Series(res, name="r2")


def fit(i):
    warnings.filterwarnings("ignore")
    limit_threads(1)
    gene_set, kind, j, w = S["sets"][i]
    pb = S["pb"]
    comps = sorted(pb.compartment.unique())
    try:
        ds = cs.Dataset(pb[META + list(w.genes)], anchors=comps, waveform=w, name=f"{gene_set}{i}")
        pt = A104.per_tissue(ds)
        pt = pt.assign(gene_set=gene_set, kind=kind, set=j, S=pt.z2 * pt.V, own=2 * pt.psi)
        pt = pt.join(departure_r2(ds, comps, S["scores"], pb))
        return pt.rename_axis("tissue").reset_index()
    except Exception as e:  # noqa: BLE001
        return pd.DataFrame([dict(gene_set=gene_set, kind=kind, set=j, error=f"{type(e).__name__}: {e}"[:200])])


def boot_median(x, rng, n=2000):
    b = np.array([np.median(rng.choice(x, len(x))) for _ in range(n)])
    return np.quantile(b, [0.05, 0.95])


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--check", action="store_true", help="3 matched sets; validation only, no core value printed")
    a = ap.parse_args()
    t0 = time.time()
    pb = load()
    limit_threads(1)
    OUT.mkdir(parents=True, exist_ok=True)
    genes = [c for c in pb.columns if c not in META]
    x = pb[genes].to_numpy(np.float32)
    lm = pd.Series(np.log((x / pb.depth.to_numpy(float)[:, None]).mean(0) + 1e-12), index=genes)[x.mean(0) >= 1]
    del x
    wf = cs.gtex()
    core_idx = [k for k, g in enumerate(wf.genes) if canonical(g) not in INDUCED]
    core = wf.subset([wf.genes[k] for k in core_idx])
    n = 3 if a.check else 100
    matched = pseudo_waveforms(lm, n, np.random.default_rng([294, 1]), wf, kind="matched")  # A294b's sets
    sets = [("full", "clock", -1, wf), ("core", "clock", -1, core)]
    for j, w in enumerate(matched):
        sets.append(("full", "matched", j, w))
        sets.append(("core", "matched", j, Waveform(tuple(w.genes[k] for k in core_idx), core.amp, core.psi, None, None, f"core{j}")))
    if a.check:
        sets = [s for s in sets if s[0] == "full"]
    keep = sorted(set(g for *_, w in sets for g in w.genes))
    resp = [g for g in dict.fromkeys(GC_TARGETS + IMMEDIATE_EARLY) if g in pb.columns]
    S.update(pb=pb[META + sorted(set(keep) | set(resp))].copy(), sets=sets)
    S["scores"] = responsive_scores(S["pb"])
    del pb
    print(f"{len(sets)} fits, {len(keep)} genes, responsive genes {len(resp)} {time.time() - t0:.0f}s", flush=True)
    rows = []
    with Pool(a.jobs) as pool:
        for k, r in enumerate(pool.imap_unordered(fit, range(len(sets)))):
            rows.append(r)
            if k % 10 == 0:
                print(f"{k + 1}/{len(sets)} fits {time.time() - t0:.0f}s", flush=True)
    res = pd.concat(rows, ignore_index=True)
    err = res.get("error", pd.Series(dtype=str)).notna().sum()
    full = res[(res.gene_set == "full") & (res.kind == "clock")].set_index("tissue")
    rep = A104.rho(full[["own_amplitude", "coupling"]])
    fm = res[(res.gene_set == "full") & (res.kind == "matched")].dropna(subset=["own"])
    lr_full = np.log(full.own / fm.groupby("tissue").own.median())
    print(f"validation: rho {rep:.3f} (A104 0.734); full median log ratio {np.median(lr_full):.3f} "
          f"(A294b {A294B_MEDIAN}; exact only with 100 sets); errors {err}", flush=True)
    if round(rep, 3) != 0.734 or (not a.check and round(float(np.median(lr_full)), 3) != A294B_MEDIAN):
        raise SystemExit("A294b not reproduced: stop before reading the core loop")
    if a.check:
        print(f"check done {time.time() - t0:.0f}s", flush=True)
        return
    res.to_csv(OUT / "sets.csv", index=False)
    rng = np.random.default_rng(329)
    tab, tests = {}, []
    for gs in ("full", "core"):
        c = res[(res.gene_set == gs) & (res.kind == "clock")].set_index("tissue")
        m = res[(res.gene_set == gs) & (res.kind == "matched")].dropna(subset=["own"])
        t = c[["S", "own", "coupling", "r2"]].add_suffix(f"_{gs}")
        t[f"own_matched_median_{gs}"] = m.groupby("tissue").own.median()
        t[f"log_ratio_{gs}"] = np.log(c.own / t[f"own_matched_median_{gs}"])
        t[f"r2_matched_median_{gs}"] = m.groupby("tissue").r2.median()
        t[f"r2_matched_q95_{gs}"] = m.groupby("tissue").r2.quantile(0.95)
        t[f"r2_excess_{gs}"] = c.r2 - t[f"r2_matched_median_{gs}"]
        t[f"coupling_matched_median_{gs}"] = m.groupby("tissue").coupling.median()
        tab[gs] = t
    t = tab["full"].join(tab["core"])
    # matched-set coupling change at the same positions (11 minus 14 genes), median over sets per tissue
    mf = res[(res.gene_set == "full") & (res.kind == "matched")].set_index(["tissue", "set"]).coupling
    mc = res[(res.gene_set == "core") & (res.kind == "matched")].set_index(["tissue", "set"]).coupling
    t["dc_matched_median"] = (mc - mf).groupby(level="tissue").median()
    t["delta_log_ratio"] = t.log_ratio_core - t.log_ratio_full
    t["delta_coupling"] = t.coupling_core - t.coupling_full
    t.to_csv(OUT / "tissues.csv")
    lo, hi = boot_median(t.log_ratio_core.to_numpy(), rng)
    p1 = wilcoxon(t.delta_log_ratio, alternative="less").pvalue
    med_core = float(t.log_ratio_core.median())
    h1 = dict(test="H329.1", tissues=len(t), median_log_ratio_full=float(t.log_ratio_full.median()),
              median_log_ratio_core=med_core, core_lo90=lo, core_hi90=hi,
              median_delta=float(t.delta_log_ratio.median()), p_wilcoxon_less=float(p1))
    h1["falls"] = bool(h1["median_delta"] < 0 and p1 < 0.05 and med_core <= A294B_MEDIAN / 2)
    p2 = wilcoxon(t.delta_coupling, alternative="greater").pvalue
    rel = float(t.delta_coupling.median() - t.dc_matched_median.median())
    h2 = dict(test="H329.2", tissues=len(t), median_delta_coupling=float(t.delta_coupling.median()),
              matched_median_delta_coupling=float(t.dc_matched_median.median()), clock_minus_matched=rel,
              p_wilcoxon_greater=float(p2))
    h2["rises"] = bool(p2 < 0.05 and rel > 0)
    tests += [h1, h2]
    for gs in ("full", "core"):
        ex = t[f"r2_excess_{gs}"].to_numpy()
        lo, hi = boot_median(ex, rng)
        tests.append(dict(test="H329.3" if gs == "full" else "H329.3 core (descriptive)", tissues=len(t),
                          median_r2=float(t[f"r2_{gs}"].median()), median_r2_matched=float(t[f"r2_matched_median_{gs}"].median()),
                          median_excess=float(np.median(ex)), excess_lo90=lo, excess_hi90=hi,
                          tissues_above_matched_q95=int((t[f"r2_{gs}"] > t[f"r2_matched_q95_{gs}"]).sum()),
                          explained=bool(np.median(ex) > 0 and lo > 0)))
    h3 = tests[2]["explained"]
    reading = ("leakage" if h1["falls"] and h2["rises"] and h3 else
               "clock-specific departure" if not h1["falls"] and not h3 else "partly leakage")
    tests.append(dict(test="reading", reading=reading, errors=int(err)))
    pd.DataFrame(tests).to_csv(OUT / "tests.csv", index=False)
    print(t[["S_full", "log_ratio_full", "log_ratio_core", "coupling_full", "coupling_core", "dc_matched_median",
             "r2_full", "r2_matched_median_full"]].round(3).sort_values("S_full").to_string(), flush=True)
    print(pd.DataFrame(tests).round(4).T.to_string(), flush=True)
    print(f"done {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
