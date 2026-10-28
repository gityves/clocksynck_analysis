"""A150: is the healthy single-cell coupling of Figure 2 clock-specific, and is it inflated by ambient RNA?

Registered in paper/docs/registered_experiments.md (A150). Covered atlases: the gate-passing A17 healthy atlases with a
genome-wide twin (lung_hlca, organ_breast, organ_muscle, organ_colon, heart if its rows reproduce the A17 `heart`).
Rows, donors, strata and anchors are those of the A17 learned rule (`Target`, `anchor_sets`); counts come from the
genome-wide table restricted to the A17 (donor, compartment) rows (for lung_hlca the genome build differs from the
clock-gene table in part of the rows, see the reproduction table).

  Part 1  200 expression-matched pseudo-clock waveforms per atlas (`clocksync.paired.pseudo_waveforms`, GTEx human
          waveform amplitudes and phases; matched on the log of the mean gene fraction over all rows of the atlas, genes
          with mean count >= 1 plus the clock genes); the A17 anchored fit (point estimates, no gate) per set.
          H150.1: pooled mean coupling, clock vs pseudo sets (set j pooled across atlases), p = (1 + #{pseudo >= clock})/(1 + n).
  Part 2  marker-based ambient over-correction (SoupX-like): soup of a donor = its counts summed over compartments; per
          compartment the 30 non-clock genes with the largest soup / compartment fraction ratio are exclusion genes; rho_dk
          from them, capped at 0.5; clock counts x' = max(x - rho depth s, 0), depth' = depth (1 - rho); refit.

    python compute_a150.py [--subsample N] [--n-pseudo 200] [--jobs 4] [--only NAME ...]
Outputs: paper/exploratory/A150_atlas_specificity/output/{reproduction,per_compartment,per_family,h150_1,d150_2,
pseudo_pooled,pseudo_failures}.csv (aggregates only; healthy atlases are public data).
"""

from __future__ import annotations

import importlib.util
import sys
import time
import warnings
from pathlib import Path

PAPER = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(PAPER / "results/figure_2/scripts/atlas_fits"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync._core.controls import CLOCK_LIKE  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.paired import pseudo_waveforms  # noqa: E402
from common import cli  # noqa: E402
from common.paths import exploratory_dir  # noqa: E402

from _shared import Target, frozen  # noqa: E402  (atlas_fits)


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


GRS = _load("a150_gene_rhythms_shared", PAPER / "results/figure_2/scripts/gene_rhythms/_shared.py")
FAM = _load("a150_anchor_rule", PAPER / "results/extended_data/scripts/host_identity/anchor_rule.py")

AMENDMENT = "A150_atlas_specificity"
SEED = 150
ATLASES = {"lung_hlca": "lung_hlca_genome", "heart": "heart_global_genome", "organ_breast": "organ_breast_genome",
           "organ_muscle": "organ_muscle_genome", "organ_colon": "organ_colon_genome"}
META = ["donor", "compartment", "stratum", "n_cells", "depth"]
CLOCK = list(cs.gtex().genes)
N_EXCL, SOUP_MIN, RHO_CAP, EPS = 30, 1e-5, 0.5, 1e-8
T0 = time.time()


def log(*a):
    print(f"[{time.time() - T0:7.0f} s]", *a, flush=True)


def learned_anchors(name: str, comps) -> list:
    a = frozen()
    a = a[a.target == name]
    anc = a.compartment[a["anchor"].astype(bool)].tolist()
    anc = [c for c in anc if c in set(comps)]
    if len(anc) < 3 or str(a.status.iloc[0]).startswith("no anchor"):
        raise ValueError(f"{name}: no learned anchors")
    return anc


def load_atlas(name: str, subsample, seed) -> dict:
    """A17 target rows with the genome-wide counts (A17 strata)."""
    T = Target(name, subsample, seed)
    g, _ = GRS.genome_table(ATLASES[name])
    genes = [c for c in g.columns if c not in META]
    key = T.pb[["donor", "compartment", "stratum"]]
    G = key.merge(g.drop(columns="stratum"), on=["donor", "compartment"], how="left", validate="1:1")
    if G[genes[0]].isna().any():
        raise ValueError(f"{name}: A17 rows missing from the genome table")
    ok = G.copy()
    # reproduction of the A17 clock-gene table (rows, counts)
    a = T.pb.set_index(["donor", "compartment"])[CLOCK + ["n_cells", "depth"]].astype(float)
    b = G.set_index(["donor", "compartment"])[CLOCK + ["n_cells", "depth"]].astype(float)
    rep = dict(rows_a17=len(a), rows_genome_matched=len(b), rows_with_count_diff=int((np.abs(a.to_numpy() - b.loc[a.index].to_numpy()).max(1) > 0).sum()),
               max_abs_count_diff=float(np.abs(a.to_numpy() - b.loc[a.index].to_numpy()).max()))
    return dict(name=name, T=T, pb=ok, genes=genes, anchors=learned_anchors(name, T.pb.compartment.unique()), rep=rep)


def fit_table(pb: pd.DataFrame, anchors, donors, wf=None, name="") -> pd.DataFrame:
    kw = {} if wf is None else dict(waveform=wf)
    return cs.Dataset(pb, anchors=anchors, donors=donors, name=name, **kw).fit(n_boot=0, gate=False).table


def _pseudo_chunk(pb: pd.DataFrame, anchors, donors, wfs, comps):
    """Couplings (n_sets x n_comp) of a chunk of pseudo-clock waveforms; a failed fit is a row of NaN with its error."""
    limit_threads(1)
    warnings.filterwarnings("ignore")
    out, err = np.full((len(wfs), len(comps)), np.nan), []
    for i, w in enumerate(wfs):
        try:
            t = fit_table(pb[META + list(w.genes)], anchors, donors, w)
            out[i] = t.coupling.reindex(comps).to_numpy(float)
            err.append("")
        except Exception as e:  # noqa: BLE001  (recorded as NaN and counted)
            err.append(f"{type(e).__name__}: {str(e)[:120]}")
    return out, err


def pseudo_sets(A: dict, k: int, n: int) -> list:
    pb, genes = A["pb"], A["genes"]
    X = pb[genes].to_numpy(float)
    mean_count = X.mean(0)
    frac = (X / pb.depth.to_numpy(float)[:, None]).mean(0)
    keep = (mean_count >= 1) | np.isin(genes, CLOCK)
    lm = pd.Series(np.log(frac + 1e-12), index=genes)[keep]
    return pseudo_waveforms(lm, n, np.random.default_rng([SEED, k]), cs.gtex())


# ---- part 2: ambient RNA ---------------------------------------------------------------------------------------

def ambient_correct(A: dict):
    """(corrected table, per-row rho, exclusion genes per compartment) of the marker-based over-correction."""
    pb, genes = A["pb"], A["genes"]
    X = pb[genes].to_numpy(float)
    depth = pb.depth.to_numpy(float)
    total = X.sum()
    soup_atlas = X.sum(0) / total
    drop = np.isin(genes, [g for g in genes if g in CLOCK_LIKE])
    donor_codes, donors = pd.factorize(pb.donor)
    S = np.zeros((len(donors), X.shape[1]))
    np.add.at(S, donor_codes, X)
    S /= S.sum(1, keepdims=True)  # soup fractions per donor
    rho = np.zeros(len(pb))
    excl = {}
    cidx = [genes.index(g) for g in CLOCK]
    for k in pb.compartment.unique():
        m = (pb.compartment == k).to_numpy()
        fk = X[m].sum(0) / X[m].sum()
        ratio = soup_atlas / (fk + EPS)
        ratio[drop | (soup_atlas < SOUP_MIN)] = -np.inf
        N = np.argsort(-ratio)[:N_EXCL]
        excl[k] = [genes[i] for i in N]
        num = X[m][:, N].sum(1)
        den = depth[m] * S[donor_codes[m]][:, N].sum(1)
        r = np.divide(num, den, out=np.zeros_like(num), where=den > 0)
        rho[m] = np.clip(r, 0, RHO_CAP)
    Xc = X[:, cidx]
    soup_clock = S[donor_codes][:, cidx]
    corr = np.maximum(Xc - rho[:, None] * depth[:, None] * soup_clock, 0.0)
    out = pb[META].copy()
    out[CLOCK] = corr
    out["depth"] = depth * (1 - rho)
    return out, rho, excl


# ---- summaries -------------------------------------------------------------------------------------------------

def per_family_tables(pc: pd.DataFrame, P: dict, n_sets: int):
    """Family summaries of D150.1 (family medians of every pseudo set) and D150.2."""
    rows = []
    for fam, g in pc.groupby("family"):
        fm = np.array([np.nanmedian([P[(a, c)][j] for a, c in zip(g.atlas, g.compartment)]) for j in range(n_sets)])
        clock = float(np.median(g.c_clock))
        rows.append(dict(family=fam, n_compartments=len(g), n_atlases=g.atlas.nunique(), c_clock_median=clock,
                         c_pseudo_median=float(np.nanmedian(fm)), c_pseudo_q95=float(np.nanquantile(fm, 0.95)),
                         excess=clock - float(np.nanmedian(fm)), above_q95=bool(clock > np.nanquantile(fm, 0.95)),
                         frac_compartments_above_own_q95=float((g.c_clock > g.c_pseudo_q95).mean()),
                         rho_median=float(g.rho_median.median()), c_corrected_median=float(g.c_corrected.median()),
                         change_median=float((g.c_corrected - g.c_clock).median()),
                         abs_change_median=float((g.c_corrected - g.c_clock).abs().median())))
    return pd.DataFrame(rows).sort_values("c_clock_median", ascending=False)


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--only", nargs="*", default=None)
    ap.add_argument("--n-pseudo", type=int, default=None)
    ap.add_argument("--chunk", type=int, default=5, help="pseudo sets per task")
    args = ap.parse_args(argv)
    run = cli.setup(args)
    n_sets = args.n_pseudo or (6 if run["quick"] else 200)
    out = exploratory_dir(AMENDMENT, run["quick"])
    names = [n for n in ATLASES if args.only is None or n in args.only]
    limit_threads(1)
    warnings.filterwarnings("ignore")

    atlas, rep_rows, comp_rows, clock_c, corr_c, rho_c, excl_rows = {}, [], [], {}, {}, {}, []
    for k, n in enumerate(names):
        A = atlas[n] = load_atlas(n, args.subsample, args.seed)
        T, pb = A["T"], A["pb"]
        t = fit_table(pb[META + CLOCK], A["anchors"], T.donors, name=n)
        clock_c[n] = t.coupling
        ref = pd.read_csv(PAPER / "results/figure_2/output/atlas_fits/a17_healthy_fits.csv")
        ref = ref[(ref.target == n) & (ref.rule == "learned") & (ref.status == "ok")].set_index("compartment").coupling
        both = t.coupling.index.intersection(ref.index)
        d = (t.coupling.loc[both] - ref.loc[both]).abs()
        rep_rows.append(dict(atlas=n, **A["rep"], n_anchors=len(A["anchors"]), compartments_a17=len(ref), compartments_here=len(t),
                             compartments_compared=len(both), max_abs_coupling_diff=float(d.max()),
                             median_abs_coupling_diff=float(d.median()), spearman=float(spearmanr(t.coupling.loc[both], ref.loc[both])[0]),
                             subsample=args.subsample if args.subsample else ""))
        log("clock fit", n, "max |dc| vs A17", round(float(d.max()), 4), A["rep"])
        corr, rho, excl = ambient_correct(A)
        corr_c[n] = fit_table(corr, A["anchors"], T.donors, name=n + "_corrected").coupling
        pbm = pb[["donor", "compartment"]].assign(rho=rho)
        rho_c[n] = pbm.groupby("compartment").rho.median()
        excl_rows += [dict(atlas=n, compartment=c, genes=";".join(g)) for c, g in excl.items()]
        log("ambient", n, "median rho", round(float(np.median(rho)), 4))
    pd.DataFrame(rep_rows).to_csv(out / "reproduction.csv", index=False)
    pd.DataFrame(excl_rows).to_csv(out / "exclusion_genes.csv", index=False)

    # part 1
    tasks, sets = [], {}
    for k, n in enumerate(names):
        A = atlas[n]
        wfs = pseudo_sets(A, k, n_sets)
        comps = list(clock_c[n].index)
        for i in range(0, n_sets, args.chunk):
            w = wfs[i:i + args.chunk]
            cols = sorted({g for x in w for g in x.genes})
            tasks.append((n, i, A["pb"][META + cols], A["anchors"], A["T"].donors, w, comps))
        sets[n] = wfs
    log(f"{len(tasks)} pseudo tasks, {n_sets} sets x {len(names)} atlases")
    res = Parallel(n_jobs=args.jobs, verbose=5)(delayed(_pseudo_chunk)(*t[2:]) for t in tasks)
    pc_mat = {n: np.full((n_sets, len(clock_c[n])), np.nan) for n in names}
    errs = []
    for t, (m, e) in zip(tasks, res):
        pc_mat[t[0]][t[1]:t[1] + len(m)] = m
        errs += [dict(atlas=t[0], set=t[1] + i, error=x) for i, x in enumerate(e) if x]
    pd.DataFrame(errs, columns=["atlas", "set", "error"]).to_csv(out / "pseudo_failures.csv", index=False)
    log("pseudo fits done; failed fits:", len(errs))

    # per compartment
    P = {}
    for n in names:
        for j, c in enumerate(clock_c[n].index):
            P[(n, c)] = pc_mat[n][:, j]
            v = P[(n, c)]
            comp_rows.append(dict(atlas=n, compartment=c, family=FAM.family(c), c_clock=clock_c[n][c],
                                  c_pseudo_median=np.nanmedian(v), c_pseudo_q95=np.nanquantile(v, 0.95),
                                  excess=clock_c[n][c] - np.nanmedian(v), n_pseudo_valid=int(np.isfinite(v).sum()),
                                  rho_median=rho_c[n].get(c, np.nan), c_corrected=corr_c[n].get(c, np.nan)))
    pc = pd.DataFrame(comp_rows)
    pc.to_csv(out / "per_compartment.csv", index=False)

    # H150.1: set j pooled across atlases
    allc = pc[["atlas", "compartment"]].to_records(index=False)
    pooled = np.array([np.nanmean([P[(a, c)][j] for a, c in allc]) for j in range(n_sets)])
    clock = float(pc.c_clock.mean())
    ok = np.isfinite(pooled)
    p = max((1 + np.sum(pooled[ok] >= clock)) / (1 + ok.sum()), 0.005)
    pd.DataFrame(dict(set=range(n_sets), pooled_mean=pooled)).to_csv(out / "pseudo_pooled.csv", index=False)
    pd.DataFrame([dict(atlases=";".join(names), n_atlases=len(names), n_compartments=len(pc), clock_pooled_mean=clock,
                       pseudo_mean=float(np.nanmean(pooled)), pseudo_median=float(np.nanmedian(pooled)),
                       pseudo_q05=float(np.nanquantile(pooled, 0.05)), pseudo_q95=float(np.nanquantile(pooled, 0.95)),
                       pseudo_max=float(np.nanmax(pooled)), n_sets=n_sets, n_valid_sets=int(ok.sum()),
                       n_failed_fits=len(errs), p=p, passes=bool(p < 0.05))]).to_csv(out / "h150_1.csv", index=False)
    fam = per_family_tables(pc, P, n_sets)
    fam.to_csv(out / "per_family.csv", index=False)

    # D150.2
    ok2 = pc.c_corrected.notna()
    rows = []
    for n, g in [("all", pc[ok2])] + [(n, pc[ok2 & (pc.atlas == n)]) for n in names]:
        rows.append(dict(scope=n, n=len(g), spearman_before_after=float(spearmanr(g.c_clock, g.c_corrected)[0]) if len(g) > 2 else np.nan,
                         median_abs_change=float((g.c_corrected - g.c_clock).abs().median()),
                         median_change=float((g.c_corrected - g.c_clock).median()), median_rho=float(g.rho_median.median())))
    d2 = pd.DataFrame(rows)
    fo = fam.set_index("family")
    d2["family_order_before"] = ""
    d2.loc[0, "family_order_before"] = ">".join(fam.sort_values("c_clock_median", ascending=False).family)
    d2["family_order_after"] = ""
    d2.loc[0, "family_order_after"] = ">".join(fam.sort_values("c_corrected_median", ascending=False).family)
    d2.to_csv(out / "d150_2.csv", index=False)
    log("done")
    print(pd.DataFrame(rep_rows).to_string()); print(pd.read_csv(out / "h150_1.csv").T.to_string())
    print(fam.to_string()); print(d2.to_string())


if __name__ == "__main__":
    main()
