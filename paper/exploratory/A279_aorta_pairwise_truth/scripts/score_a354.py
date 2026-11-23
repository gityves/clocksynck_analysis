"""A354 (registered in `paper/docs/registered_experiments.md`): companion gene sets for the aorta cell-type pairs, clock
intact (aligned arm) against the Bmal1 knockout.

As A279 (library x cell type pseudobulks of SMC, fibroblast, EC, macrophage, >= 30 cells, compartments in >= 10
libraries of the arm), but summed over all genes. Per arm and pair, the A272 pairwise |r| for the 14 clock genes
(must reproduce A279 `pairs.csv`, else the run stops) and for 200 expression-matched and 200 stress-responsive gene
sets (`clocksync.paired.pseudo_waveforms`, the arm's mean log fractions, mouse symbols upper-cased; own projector per
set; seeds [354, arm, family]).
H354.1: knockout stromal pairs (EC-fibroblast, EC-SMC, fibroblast-SMC): clock |r| <= matched q95 in all three.
H354.2: aligned: clock |r| > matched q95 in >= 4 of 6 pairs. D354.3: responsive family, misaligned arm, macrophage pairs.

    python score_a354.py [--n-pseudo 200] [--check]
Outputs: ../output/a354/{pairs,design}.csv, verdicts.csv.
"""
from __future__ import annotations

import argparse
import importlib.util
import itertools
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
sys.path.insert(0, str(PAPER))
import anndata as ad  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scipy.sparse as sp  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.paired import pseudo_waveforms  # noqa: E402
from clocksync._core.controls import RESPONSIVE_GENES  # noqa: E402


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


A279 = _load("a354_a279", HERE.parent / "score_a279.py")
A272, A201 = A279.A272, A279.A201
OUT = HERE.parents[1] / "output" / "a354"
SEED, FAMILIES = 354, ("matched", "responsive")
SV_PAIRS = {("EC", "Fibroblast"), ("EC", "SMC"), ("Fibroblast", "SMC")}
CHUNK = 20000


def genome_pseudobulk():
    """lib (as A201.pseudobulk) and genome-wide counts per (library, cell type), columns = upper-cased symbols."""
    o, Xc, _ = A201.load()
    keep = o.cell_type.astype(str).isin(A279.CLASSES).to_numpy()
    o = o[keep].reset_index(drop=True)
    lib, _ = A201.pseudobulk(o, Xc[keep])
    key = (o.sample_id.astype(str) + "|" + o.cell_type.astype(str)).to_numpy()
    codes = pd.Index(lib.index).get_indexer(key)
    assert (codes >= 0).all()
    a = ad.read_h5ad(A201.H5, backed="r")
    rows = np.flatnonzero(keep)
    G = np.zeros((len(lib), a.n_vars))
    for i in range(0, len(rows), CHUNK):
        X = a.X[rows[i:i + CHUNK]]
        X = X if sp.issparse(X) else sp.csr_matrix(X)
        M = sp.csr_matrix((np.ones(X.shape[0]), (codes[i:i + CHUNK], np.arange(X.shape[0]))), shape=(len(lib), X.shape[0]))
        G += (M @ X).toarray()
    sym = pd.Index(a.var_names.astype(str)).str.upper()
    cnt = pd.DataFrame(G, columns=sym).T.groupby(level=0).sum().T
    return lib, cnt, float(np.abs(G - np.round(G)).max())


def arm_tables(lib, cnt):
    out = {}
    for arm, cond in A279.ARMS.items():
        m = (lib.cond == cond).to_numpy() & (lib.cells >= A279.MIN_CELLS).to_numpy()
        key = lib[m][["sample_id", "cell_type"]].astype(str).reset_index(drop=True)
        pb = cs.from_bulk(cnt[m].round().reset_index(drop=True), key, donor="sample_id", compartment="cell_type",
                          depth=pd.Series(lib.total.to_numpy()[m]), genes="all")
        pb["stratum"] = arm
        need = min(A279.MIN_LIBS, lib[lib.cond == cond].sample_id.nunique())
        n = pb.groupby("compartment").donor.nunique()
        out[arm] = (pb[pb.compartment.isin(n.index[n >= need])].reset_index(drop=True), need)
    return out


def pairs_r(t, need, wf):
    z, v = A272.atlas_states(t, wf)
    rows = {c: g.set_index("donor")["_i"] for c, g in t.assign(_i=np.arange(len(t))).groupby("compartment")}
    out = {}
    for a, b in itertools.combinations(sorted(rows), 2):
        d = rows[a].index.intersection(rows[b].index)
        if len(d) >= need:
            out[(a, b)] = abs(A272.pair_r(z, v, dict(ia=rows[a].loc[d].to_numpy(), ib=rows[b].loc[d].to_numpy(),
                                                     s=np.zeros(len(d), int))))
    return out


def companions(t, need, k, fam, n):
    genes = [g for g in t.columns if g not in ("donor", "compartment", "stratum", "depth")]
    X = t[genes].to_numpy(float)
    frac = (X / t.depth.to_numpy(float)[:, None]).mean(0)
    clock = set(cs.gtex().genes)
    keep = (X.mean(0) >= 1) | np.isin(genes, list(clock))
    lm = pd.Series(np.log(frac + 1e-12), index=genes)[keep]
    wfs = pseudo_waveforms(lm, n, np.random.default_rng([SEED, k, FAMILIES.index(fam)]), cs.gtex(), fam)
    meta = ["donor", "compartment", "stratum", "depth"]
    R = [pairs_r(t[meta + list(w.genes)], need, w) for w in wfs]
    return pd.DataFrame(R), int(lm.index.isin(RESPONSIVE_GENES).sum())


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--n-pseudo", type=int, default=200)
    ap.add_argument("--check", action="store_true", help="pseudobulks, design and the A279 reproduction only")
    a = ap.parse_args()
    limit_threads(4)
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    lib, cnt, frac_err = genome_pseudobulk()
    print(f"pseudobulks {cnt.shape}; max non-integer part {frac_err:.3g}; {time.time() - t0:.0f}s", flush=True)
    tabs = arm_tables(lib, cnt)
    ref = pd.read_csv(HERE.parents[1] / "output" / "pairs.csv").set_index(["arm", "a", "b"]).r_abs
    rows, design = [], []
    for k, (arm, (t, need)) in enumerate(tabs.items()):
        clock = pairs_r(t, need, cs.gtex())
        for (x, y), r in clock.items():  # reproduction of A279 before any companion
            assert np.isclose(r, ref[arm, x, y], atol=1e-6), (arm, x, y, r, ref[arm, x, y])
        print(arm, "clock |r| reproduces A279", flush=True)
        if a.check:
            continue
        comp = {}
        for fam in FAMILIES:
            comp[fam], n_resp = companions(t, need, k, fam, a.n_pseudo)
            print(f"  {arm} {fam}: {len(comp[fam])} sets; {time.time() - t0:.0f}s", flush=True)
        design.append(dict(arm=arm, libraries=t.donor.nunique(), compartments=t.compartment.nunique(),
                           responsive_genes_expressed=n_resp))
        for (x, y), r in clock.items():
            row = dict(arm=arm, a=x, b=y, stromal=(x, y) in SV_PAIRS, clock_r=r)
            for fam in FAMILIES:
                c = comp[fam][(x, y)].to_numpy(float)
                row.update({f"{fam}_median": np.median(c), f"{fam}_q05": np.quantile(c, .05),
                            f"{fam}_q95": np.quantile(c, .95), f"{fam}_pct": 100 * np.mean(c < r)})
            rows.append(row)
    if a.check:
        return
    p = pd.DataFrame(rows)
    p.to_csv(OUT / "pairs.csv", index=False)
    pd.DataFrame(design).to_csv(OUT / "design.csv", index=False)
    ko = p[(p.arm == "bmal1_ko") & p.stromal]
    al = p[p.arm == "aligned"]
    v = pd.DataFrame([
        dict(test="H354.1", arm="bmal1_ko", pairs=len(ko), count=int((ko.clock_r <= ko.matched_q95).sum()),
             rule="clock |r| <= matched q95 in all 3 stromal pairs",
             verdict="PASS" if len(ko) == 3 and (ko.clock_r <= ko.matched_q95).all() else "FAIL"),
        dict(test="H354.2", arm="aligned", pairs=len(al), count=int((al.clock_r > al.matched_q95).sum()),
             rule="clock |r| > matched q95 in >= 4 of 6 pairs",
             verdict="PASS" if (al.clock_r > al.matched_q95).sum() >= 4 else "FAIL")])
    v.to_csv(OUT / "verdicts.csv", index=False)
    pd.set_option("display.width", 250)
    print(pd.DataFrame(design).to_string(index=False), p.round(3).to_string(index=False), v.to_string(index=False),
          sep="\n")
    print(f"runtime {time.time() - t0:.0f}s")


if __name__ == "__main__":
    main()
