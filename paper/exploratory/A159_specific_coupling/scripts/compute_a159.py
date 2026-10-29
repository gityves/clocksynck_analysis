"""A159: clock-specific coupling (delta c) for every healthy atlas with a genome-wide twin on disk.

Registered in paper/docs/registered_experiments.md (A159). Reuses the A150 code (`compute_a150.py`: genome-wide twin,
reproduction, matched pseudo-clock sets, chunked fits). Differences: relearned anchors (2026-10-28), `anchor_readout="loo"`,
200 matched plus 200 responsive-gene sets per atlas, and the registered exact row check.

  Row check   per atlas, every A17 row (donor x compartment) must be present in the genome-wide table with identical
              clock-gene counts, n_cells and depth. Atlases that fail, cannot be built or do not pass the gate are listed in
              `row_check.csv` and are not covered. (--also-failed fits them too, flagged `covered=False`: descriptive.)
  Fits        per covered atlas: the anchored fit (`fit(gate=False)`, point estimates) under the clock waveform and under
              each pseudo set.
  H159.1-2    over the covered atlases not in A150 (lung_hlca, heart, organ_breast, organ_muscle, organ_colon); Holm.
  D159.1      per atlas x compartment c_clock, matched and responsive median / q95, delta c = c_clock - matched median,
              `clock_specific` = c_clock > matched q95 and > responsive q95; per family median delta c and share clock-specific.

    python compute_a159.py [--stage check|fit] [--n-pseudo 200] [--jobs 4] [--only NAME ...] [--also-failed] [--subsample N]
Outputs: paper/exploratory/A159_specific_coupling/output/{row_check,d159_compartments,d159_families,h159,pseudo_matched,
pseudo_responsive}.csv (aggregates; public atlases).
"""

from __future__ import annotations

import importlib.util
import sys
import time
import warnings
from pathlib import Path

PAPER = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PAPER))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.paired import pseudo_waveforms  # noqa: E402
from common import cli  # noqa: E402
from common.paths import exploratory_dir  # noqa: E402

spec = importlib.util.spec_from_file_location("a150", PAPER / "exploratory/A150_atlas_specificity/scripts/compute_a150.py")
A150 = importlib.util.module_from_spec(spec)
spec.loader.exec_module(A150)

AMENDMENT, SEED = "A159_specific_coupling", 159
ATLASES = dict(A150.ATLASES) | {"organ_adipose": "organ_adipose_genome"}
IN_A150 = set(A150.ATLASES)
META, CLOCK = A150.META, A150.CLOCK
STROMAL_VASCULAR = ("fibroblast_stromal", "endothelial")
T0 = time.time()


def log(*a):
    print(f"[{time.time() - T0:7.0f} s]", *a, flush=True)


def check_rows(name: str, subsample, seed) -> tuple[dict, dict | None]:
    """Registered row check. Returns (row, atlas dict or None). Pass = all rows present and identical counts."""
    row = dict(atlas=name, genome_dataset=ATLASES[name], in_a150=name in IN_A150)
    try:
        T = A150.Target(name, subsample, seed)
        g, _ = A150.GRS.genome_table(ATLASES[name])
    except Exception as e:  # noqa: BLE001
        return row | dict(built=False, error=f"{type(e).__name__}: {str(e)[:150]}", row_check_pass=False), None
    a = T.pb.set_index(["donor", "compartment"])[CLOCK + ["n_cells", "depth"]].astype(float)
    b = g.assign(donor=g.donor.astype(str)).set_index(["donor", "compartment"])
    present = a.index.isin(b.index)
    bb = b.loc[a.index[present], CLOCK + ["n_cells", "depth"]].astype(float)
    diff = np.abs(a[present].to_numpy() - bb.to_numpy())
    bad = (diff.max(1) > 0) if len(diff) else np.array([], bool)
    row |= dict(built=True, rows_a17=len(a), rows_missing_in_genome=int((~present).sum()), rows_with_count_diff=int(bad.sum()),
                max_abs_count_diff=float(diff.max()) if diff.size else np.nan, donors_a17=T.pb.donor.nunique(),
                compartments_a17=T.pb.compartment.nunique(), compartments_missing=int(
                    len(set(T.pb.compartment) - set(b.reset_index().compartment))),
                row_check_pass=bool(present.all() and not bad.any()), error="")
    if not present.all():
        return row, None
    try:
        return row, A150.load_atlas(name, subsample, seed)
    except Exception as e:  # noqa: BLE001
        return row | dict(error=f"{type(e).__name__}: {str(e)[:150]}", row_check_pass=False), None


def fit_table(pb, anchors, donors, wf=None, name=""):
    kw = {} if wf is None else dict(waveform=wf)
    return cs.Dataset(pb, anchors=anchors, donors=donors, name=name, **kw).fit(n_boot=0, gate=False,
                                                                              anchor_readout="loo").table


def sets_for(A: dict, k: int, n: int, kind: str) -> list:
    pb, genes = A["pb"], A["genes"]
    X = pb[genes].to_numpy(float)
    frac = (X / pb.depth.to_numpy(float)[:, None]).mean(0)
    keep = (X.mean(0) >= 1) | np.isin(genes, CLOCK)
    lm = pd.Series(np.log(frac + 1e-12), index=genes)[keep]
    return pseudo_waveforms(lm, n, np.random.default_rng([SEED, k, kind == "responsive"]), cs.gtex(), kind=kind)


def chunk(pb, anchors, donors, wfs, comps):
    limit_threads(1)
    warnings.filterwarnings("ignore")
    out, err = np.full((len(wfs), len(comps)), np.nan), []
    for i, w in enumerate(wfs):
        try:
            out[i] = fit_table(pb[META + list(w.genes)], anchors, donors, w).coupling.reindex(comps).to_numpy(float)
            err.append("")
        except Exception as e:  # noqa: BLE001
            err.append(f"{type(e).__name__}: {str(e)[:120]}")
    return out, err


def p_upper(pseudo: np.ndarray, clock: float) -> float:
    ok = np.isfinite(pseudo)
    return float((1 + np.sum(pseudo[ok] >= clock)) / (1 + ok.sum()))


def holm(p: list) -> list:
    o = np.argsort(p)
    adj, run = np.empty(len(p)), 0.0
    for r, i in enumerate(o):
        run = max(run, (len(p) - r) * p[i])
        adj[i] = min(run, 1.0)
    return adj.tolist()


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--stage", choices=("check", "fit"), default="fit")
    ap.add_argument("--only", nargs="*", default=None)
    ap.add_argument("--n-pseudo", type=int, default=None)
    ap.add_argument("--chunk", type=int, default=5)
    ap.add_argument("--also-failed", action="store_true", help="also fit atlases that fail the row check (descriptive)")
    args = ap.parse_args(argv)
    run = cli.setup(args)
    n_sets = args.n_pseudo or (6 if run["quick"] else 200)
    out = exploratory_dir(AMENDMENT, run["quick"])
    names = [n for n in ATLASES if args.only is None or n in args.only]
    limit_threads(1)
    warnings.filterwarnings("ignore")
    gates = pd.read_csv(PAPER / "results/extended_data/output/host_identity/a17_gates.csv").set_index("name").passed

    rows, atlas = [], {}
    for n in names:
        r, A = check_rows(n, args.subsample, args.seed)
        r["gate_passed"] = bool(gates.get(n, False))
        r["covered"] = bool(r["row_check_pass"] and r["gate_passed"])
        rows.append(r)
        if A is not None:
            atlas[n] = A
        log("row check", {k: r[k] for k in ("atlas", "row_check_pass", "gate_passed", "covered")}, r.get("error", ""))
    pd.DataFrame(rows).to_csv(out / "row_check.csv", index=False)
    if args.stage == "check":
        print(pd.DataFrame(rows).to_string())
        return
    fit_names = [n for n in atlas if (n in atlas and (rows[names.index(n)]["covered"] or
                                                      (args.also_failed and rows[names.index(n)]["gate_passed"])))]
    covered = {r["atlas"] for r in rows if r["covered"]}
    log("fitting", fit_names, "covered:", sorted(covered))

    clock_c, tasks, sets = {}, [], {}
    for k, n in enumerate(fit_names):
        A = atlas[n]
        clock_c[n] = fit_table(A["pb"][META + CLOCK], A["anchors"], A["T"].donors, name=n).coupling
        comps = list(clock_c[n].index)
        for kind in ("matched", "responsive"):
            wfs = sets[(n, kind)] = sets_for(A, k, n_sets, kind)
            for i in range(0, n_sets, args.chunk):
                w = wfs[i:i + args.chunk]
                cols = sorted({g for x in w for g in x.genes})
                tasks.append((n, kind, i, A["pb"][META + cols], A["anchors"], A["T"].donors, w, comps))
    log(f"{len(tasks)} pseudo tasks, {n_sets} sets x 2 kinds x {len(fit_names)} atlases")
    res = Parallel(n_jobs=args.jobs, verbose=5)(delayed(chunk)(*t[3:]) for t in tasks)
    mat = {(n, kind): np.full((n_sets, len(clock_c[n])), np.nan) for n in fit_names for kind in ("matched", "responsive")}
    errs = []
    for t, (m, e) in zip(tasks, res):
        mat[(t[0], t[1])][t[2]:t[2] + len(m)] = m
        errs += [dict(atlas=t[0], kind=t[1], set=t[2] + i, error=x) for i, x in enumerate(e) if x]
    pd.DataFrame(errs, columns=["atlas", "kind", "set", "error"]).to_csv(out / "pseudo_failures.csv", index=False)
    log("pseudo fits done; failed fits:", len(errs))
    for kind in ("matched", "responsive"):
        pd.concat([pd.DataFrame(mat[(n, kind)], columns=list(clock_c[n].index)).assign(atlas=n, set=range(n_sets))
                   for n in fit_names]).to_csv(out / f"pseudo_{kind}.csv", index=False)

    comp = []
    for n in fit_names:
        for j, c in enumerate(clock_c[n].index):
            m, r_ = mat[(n, "matched")][:, j], mat[(n, "responsive")][:, j]
            cc = float(clock_c[n][c])
            comp.append(dict(atlas=n, covered=n in covered, in_a150=n in IN_A150, compartment=c, family=A150.FAM.family(c),
                             c_clock=cc, matched_median=np.nanmedian(m), matched_q95=np.nanquantile(m, 0.95),
                             responsive_median=np.nanmedian(r_), responsive_q95=np.nanquantile(r_, 0.95),
                             delta_c=cc - np.nanmedian(m), delta_c_responsive=cc - np.nanmedian(r_),
                             clock_specific=bool(cc > np.nanquantile(m, 0.95) and cc > np.nanquantile(r_, 0.95)),
                             above_matched_q95=bool(cc > np.nanquantile(m, 0.95)),
                             above_responsive_q95=bool(cc > np.nanquantile(r_, 0.95))))
    pc = pd.DataFrame(comp)
    pc.to_csv(out / "d159_compartments.csv", index=False)
    fam = []
    for scope, g_ in (("covered", pc[pc.covered]), ("all_fitted", pc)):
        for f, g in g_.groupby("family"):
            fam.append(dict(scope=scope, family=f, n_compartments=len(g), n_atlases=g.atlas.nunique(),
                            median_c_clock=g.c_clock.median(), median_delta_c=g.delta_c.median(),
                            median_delta_c_responsive=g.delta_c_responsive.median(),
                            share_clock_specific=g.clock_specific.mean(),
                            share_above_matched_q95=g.above_matched_q95.mean()))
    pd.DataFrame(fam).to_csv(out / "d159_families.csv", index=False)

    # H159.1-2: covered atlases not in A150
    new = [n for n in fit_names if n in covered and n not in IN_A150]
    if not new:
        h = pd.DataFrame([dict(hypothesis=h_, status="not scorable: no covered atlas outside A150", atlases="", n_compartments=0)
                          for h_ in ("H159.1", "H159.2")])
    else:
        idx = [(n, c) for n in new for c in clock_c[n].index]
        col = {n: list(clock_c[n].index) for n in new}
        pooled = np.array([np.nanmean([mat[(n, "matched")][j, col[n].index(c)] for n, c in idx]) for j in range(n_sets)])
        c1 = float(np.mean([clock_c[n][c] for n, c in idx]))
        sv = [(n, c) for n, c in idx if A150.FAM.family(c) in STROMAL_VASCULAR]
        rows_h = [dict(hypothesis="H159.1", n_atlases=len(new), n_compartments=len(idx), clock=c1, pseudo_median=np.nanmedian(pooled),
                       pseudo_q95=np.nanquantile(pooled, 0.95), p=p_upper(pooled, c1))]
        if sv:
            ps = np.array([np.nanmean([mat[(n, "matched")][j, col[n].index(c)] for n, c in sv]) for j in range(n_sets)])
            c2 = float(np.mean([clock_c[n][c] for n, c in sv]))
            rows_h.append(dict(hypothesis="H159.2", n_atlases=len(new), n_compartments=len(sv), clock=c2,
                               pseudo_median=np.nanmedian(ps), pseudo_q95=np.nanquantile(ps, 0.95), p=p_upper(ps, c2)))
        h = pd.DataFrame(rows_h)
        h["p_holm"] = holm(h.p.tolist())
        h["passes"] = h.p_holm < 0.05
        h["atlases"] = ";".join(new)
    h.to_csv(out / "h159.csv", index=False)
    log("done")
    print(pd.DataFrame(rows).to_string()); print(pd.DataFrame(fam).to_string()); print(h.to_string())


if __name__ == "__main__":
    main()
