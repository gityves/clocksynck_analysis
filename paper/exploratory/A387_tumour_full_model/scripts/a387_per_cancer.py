"""A387 note 1 (POST HOC implementation of D387.1): per-cancer full-model readouts on the pooled host.

Registered in paper/docs/registered_experiments.md (A387 note 1). One host on all patients of the A387 Part 2 table
(anchors: normal fibroblasts, endothelium, mural); each tumour compartment relabelled `tumour:<kind>|<cancer>` for cancers
with >= 8 patients of that compartment with a host, the rest `tumour:<kind>|other`. Per readout (14 genes, core loop):
three fits (reference = each normal counterpart; same bootstrap seed; n_boot 500, 90 % intervals; coupling permutation
p in the first), and companion fits of the fibroblast and endothelial targets on the A387 matched and narrow responsive
sets (seeds [387, 0] and [387, 1], genome tables relabelled the same way). Cancer cells: no companions (as A387).
Output: output/score/per_cancer_{tests,companions}.csv. Descriptive.

    python a387_per_cancer.py [--jobs 12] [--genome-jobs 3]     (paper/common/run.sh, RUN_CAP_GB 40)
"""
from __future__ import annotations

import argparse
import time
import warnings

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

import a387 as R
import a387_score as S
from a387 import ANCHORS, BASE, C, NORMAL, OUT, TARGETS, log
from clocksync.config import limit_threads
from clocksync.paired import pseudo_waveforms

MIN_N, N_BOOT, N_PERM, SEED = 8, 500, 1000, 387


def relabel(table):
    """Tumour compartments relabelled by cancer; returns the table and the new target names per kind."""
    t = table.copy()
    hosted = set(t.donor[t.compartment.isin(ANCHORS)])
    cancer = {d: R.M.cancer(d.split("|", 1)[0]) for d in t.donor.unique()}
    names = {}
    for kind, v in TARGETS.items():
        m = t.compartment == v
        n = t[m & t.donor.isin(hosted)].groupby(t.donor.map(cancer)).donor.nunique()
        keep = set(n[n >= MIN_N].index)
        lab = [f"{v}|{cancer[d] if cancer[d] in keep else 'other'}" for d in t.donor[m]]
        t.loc[m, "compartment"] = lab
        names[kind] = sorted(set(lab))
    return t, names


def fit(table, rd, ref, n_perm, targets):
    limit_threads(1)
    warnings.filterwarnings("ignore")
    f = S.dataset(table, R.READOUTS[rd]).fit(n_boot=N_BOOT, seed=SEED, n_perm=n_perm, gate=False, reference=ref,
                                             targets=targets, level=0.90)
    return rd, ref, f.table


def comp_fits(G, wfs, rd, targets):
    limit_threads(1)
    warnings.filterwarnings("ignore")
    out = []
    for i, w in wfs:
        try:
            t = S.dataset(G, w).fit(n_boot=0, gate=False, targets=targets).table
            out.append(dict(readout=rd, set=i, **{k: float(t.loc[k, "coupling"]) for k in targets}))
        except Exception as e:  # noqa: BLE001
            out.append(dict(readout=rd, set=i, error=f"{type(e).__name__}: {e}"[:200]))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--jobs", type=int, default=12)
    ap.add_argument("--genome-jobs", type=int, default=3)
    args = ap.parse_args()
    limit_threads(1)
    t0 = time.time()
    out = OUT / "score"
    base = R.pooled_table(args.jobs, False)
    table, names = relabel(base)
    targets = [x for v in names.values() for x in v]
    allc = targets + [NORMAL[k] for k in TARGETS]
    log("targets: " + ", ".join(targets))
    jobs = [(rd, NORMAL[k], N_PERM if k == "fibroblast" else 0) for rd in R.READOUTS for k in TARGETS]
    res = Parallel(n_jobs=min(args.jobs, len(jobs)))(delayed(fit)(table, rd, ref, p, allc) for rd, ref, p in jobs)
    fits = {(rd, ref): t for rd, ref, t in res}
    log("fits done")
    # companions (fibroblast and endothelial targets)
    keep = set(zip(base.donor, base.compartment))
    studies = sorted({d.split("|", 1)[0] for d in base.donor})
    Gs = Parallel(n_jobs=args.genome_jobs)(delayed(S.genome_cohort)(s, {x for x in keep if x[0].startswith(s + "|")})
                                           for s in studies)
    G = pd.concat([g for g in Gs if g is not None and len(g)], ignore_index=True)
    G = G.drop(columns=["stratum"]).merge(base[["donor", "compartment", "stratum"]], on=["donor", "compartment"])
    genes = [c for c in G.columns if c not in ("donor", "compartment", "stratum", "depth", "n_cells")]
    genes = [g for g in genes if G[g].notna().all()]
    ref = G[G.compartment.isin(ANCHORS)]
    frac = ref[genes].div(ref.depth, axis=0)
    lm = np.log(frac.mean() + 1e-12)[(ref[genes].mean() >= 1) | frac.columns.isin(list(BASE.genes))]
    sets = [(f, w) for code, f in enumerate(("matched", "responsive"))
            for w in pseudo_waveforms(lm, R.N_SET, np.random.default_rng([SEED, code]), BASE, f)]
    used = sorted({g for _, w in sets for g in w.genes} | set(BASE.genes))
    Gr, _ = relabel(G[["donor", "compartment", "stratum", "depth"] + used])
    ctarg = names["fibroblast"] + names["endothelial"]
    tasks = [(rd, [(i, w if rd == "14" else C.restrict(w, BASE, C.core(BASE))) for i, (_, w) in enumerate(sets)])
             for rd in R.READOUTS]
    comp = Parallel(n_jobs=args.jobs)(delayed(comp_fits)(Gr, ch, rd, ctarg) for rd, ws in tasks
                                      for ch in S.chunks(ws, args.jobs))
    comp = pd.DataFrame([r for x in comp for r in x])
    comp["family"] = [sets[i][0] for i in comp.set]
    comp.to_csv(out / "per_cancer_companions.csv", index=False)
    log("companions done")
    rows = []
    for rd in R.READOUTS:
        t0_ = fits[(rd, NORMAL["fibroblast"])]
        for kind, ns in names.items():
            tr = fits[(rd, NORMAL[kind])]
            for v in ns:
                r = dict(readout=rd, kind=kind, cancer=v.split("|")[1], n=int(t0_.loc[v, "n"]),
                         coupling=t0_.loc[v, "coupling"], coupling_lo90=t0_.loc[v, "coupling_lo"],
                         coupling_hi90=t0_.loc[v, "coupling_hi"], p_coupling=t0_.loc[v, "p_coupling"],
                         boundary_type=t0_.loc[v, "boundary_type"], peak_h=-tr.loc[v, "offset_h"],
                         peak_lo90=-tr.loc[v, "offset_hi"], peak_hi90=-tr.loc[v, "offset_lo"],
                         loudness_ratio=float(np.exp(tr.loc[v, "log_amp_rel"])),
                         loudness_lo90=float(np.exp(tr.loc[v, "log_amp_lo"])),
                         loudness_hi90=float(np.exp(tr.loc[v, "log_amp_hi"])))
                if v in comp.columns:
                    for f in ("matched", "responsive"):
                        c = comp[(comp.readout == rd) & (comp.family == f)][v]
                        r[f"{f}_q95"], r[f"{f}_share_ge"] = float(np.nanquantile(c, .95)), float(np.nanmean(c >= r["coupling"]))
                rows.append(r)
    T = pd.DataFrame(rows)
    T.to_csv(out / "per_cancer_tests.csv", index=False)
    pd.set_option("display.width", 250)
    log("\n" + T.round(3).to_string(index=False))
    log(f"done {(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    main()
