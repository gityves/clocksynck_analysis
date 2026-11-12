"""A246: do tumour-microenvironment programs follow the clock state of the tumour, at any phase?

Registered in paper/docs/registered_experiments.md (A246). The statistic, null, pseudo-clock construction and pooling are those of A229
(`compute_a229`), generalised to eight readouts and two states. The 34 A227 cohorts and tables (`compute_a227`).

Readouts y (A146 part 2 programs, A146 gene lists and score; `part2pc_extract.PROGRAMS`): P1 proliferation, P2 MHC-I, P3 interferon response,
P9 PD-L1 in malignant cells; P4 T-cell activation and P10 PD-1 in T_NK, residual on T abundance; P5 myeloid antigen presentation and P6
phagocytosis in myeloid cells, residual on myeloid abundance; P8 log myeloid share of all cells of the sample. Gene score = mean over the program
genes of the within-cohort z of log2(CPM + 1) (genome-wide pseudobulk of the patient x compartment), re-z-scored; residual = score minus its
least-squares fit on the within-cohort z of log2((n_cells + .5) / (N + 1)) (N = cells of the sample in the pseudobulk table, >= 100).
States: (a) `endo`: z of the tumour endothelial pseudobulk (A229); (b) `own`: z of the compartment in which the program is scored (P8: myeloid),
same construction fitted on that compartment's rows of the cohort. Both centred within cohort x chemistry over the patients of the test.
Statistic: r = |sum y conj(z)| / sqrt(sum y^2 sum |z|^2) per cohort (y real, centred), pooled with patient-count weights.
Null: y permuted across patients within cohort, with ONE ordering of all patients of the cohort per iteration shared by all programs (each program
gets the permutation this ordering induces on its patients). Omnibus per state T = mean over the scored programs of (r - null mean) / null SD, its
null from the same permutations. Clock specificity: z from 200 expression-matched pseudo-clock sets (A227 construction, reference = the genome
pseudobulks of the state's compartment), same y; per program and for the omnibus (T with pseudo r standardised by the permutation null mean/SD).

  power  eligibility and the shuffled-y null only; the identity pairing of y and z is never evaluated, no observed statistic is computed, stored or
         printed. MDE = 2.49 x pooled null SD; scored = >= 40 patients and MDE <= 0.25; omnibus scored with >= 5 scored programs; its MDE is in T
         units.                                                           -> output/{eligibility,power}.csv
  tests  (IMPLEMENTED, NOT RUN on real data before the power note is committed; requires --registered-run, or --placebo) observed pooled r,
         permutation p, pseudo q95, Holm (H246.1-2 omnibus over the two states; H246.3 programs within state), pass; descriptives: percentile among
         200 matched gene sets, implied peak and its circular correlation with the A146 hours, cancer types with >= 12 patients, adjacent-normal
         side, breast cohorts alone.      -> output/{tests,pooled,per_cohort,pseudo,matched_sets}.csv
  --placebo shuffles every y once within cohort (fixed seed) before anything else and writes to output/placebo/ (code test only).

Choices not fixed by the registration:
 1. Eligibility: a program needs >= 50 cells per pseudobulk (table rule), median summed count over its genes >= 10 across the compartment's
    patients, and >= half of its genes in the cohort; no program gene is a clock (waveform) gene (asserted). The abundance of T_NK / myeloid and
    P8 use N = sum of n_cells over the side's pseudobulks of the A227 table (as A229 H229.5), not all cells of the sample; patients with N < 100
    are dropped for P4, P10, P5, P6, P8 (A146 rule); log2 as A146 (P8 log2 share, A229 used natural log; r is scale-free).
 2. Malignant compartment = `tumour:malignant` (copy-number calls) where the cohort table has it, else `tumour:epithelial`, as A146 / A149;
    its genome pseudobulks are rebuilt from counts with the same labels. Adjacent-normal side: `adjacent_normal:epithelial` for P1-3, P9.
 3. Gene z and abundance z are taken over the patients having the compartment (and N >= 100), before intersecting with the patients having the state.
 4. A cohort enters a program x state if >= 4 patients have y and the state (A229 MIN_N), no spread excluded; states of a compartment need >= 4 rows.
 5. Shared permutation: random keys (n_perm x patients of the cohort); the ordering of the keys within a program's patients is its permutation.
    Seed 246 (+ cohort, side). Pseudo set i of every compartment is set i of the pseudo-clock draws of that compartment.
 6. Omnibus T of a cancer type / breast subset: programs with >= 12 patients in it, >= 3 programs, null standardised on the subset's null;
    descriptive. Circular correlation (Jammalamadaka) of implied peak with A146 hours over the scored programs; the mirrored value (peak
    reflected) is also reported because the A146 hours are mouse-transferred clock hours and the sense of the axes is not fixed.
 7. Matched sets: per program gene one gene drawn (without replacement) from the 400 nearest in patient-weighted mean log2(CPM + 1) of the compartment
    (genes present in all cohorts; program and waveform genes excluded), scored as the program (same residual); P8 has no gene set; tumour side only.

    python compute_a246.py --stage power|tests [--jobs J] [--n-perm 10000] [--n-pseudo 200] [--smoke] [--placebo] [--out DIR]
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[2]
for p in ("A227_tme_compartment_pairing", "A149_paired_tumour_normal", "A151_fibroblast_pairing", "A226_stroma_pairing_confirm",
          "A229_endothelial_gate", "A146_tme_time"):
    sys.path.insert(0, str(PAPER / "exploratory" / p / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import a149_lib as L  # noqa: E402
import compute_a226 as K  # noqa: E402
import compute_a227 as A  # noqa: E402
import compute_a229 as A29  # noqa: E402
import part2pc_prep as P146  # noqa: E402
from clocksync import stats as S  # noqa: E402
from common import cli  # noqa: E402

OUT = HERE.parent / "output"
SEED = 246
SIDES = ["tumour", "adjacent_normal"]
STATES = {"endo": ("H246.1", "tumour endothelial clock (body time)"), "own": ("H246.2", "own compartment clock")}
# program -> (compartment key, abundance for the residual: compartment key, "share" for P8, or None)
PROGS = {"P1": ("mal", None), "P2": ("mal", None), "P3": ("mal", None), "P9": ("mal", None), "P4": ("T_NK", "T_NK"), "P10": ("T_NK", "T_NK"),
         "P5": ("myeloid", "myeloid"), "P6": ("myeloid", "myeloid"), "P8": ("myeloid", "share")}
COMPS = ["endo", "mal", "T_NK", "myeloid"]
GENES = {p: [L.canonical(g) for g in P146.PROGRAMS[p].split()] for p in PROGS}
assert not {g for v in GENES.values() for g in v} & {L.canonical(g) for g in L.q.WF.genes}, "a program gene is a clock gene"
MIN_PATIENTS, MDE_K, MAX_MDE, MIN_CANCER, MIN_N = A.MIN_PATIENTS, A.MDE_K, A.MAX_MDE, A.MIN_CANCER, L.MIN_PAIRED
MIN_COUNT, MIN_CELLS_ALL, MIN_OMNIBUS = 10, 100, 5
N_MATCH, N_POOL = 200, 400
MAL = L.q.MAL


def seed_of(*tag, study=""):
    return [SEED, *tag, sum(map(ord, study))]


# ---- per-cohort data ---------------------------------------------------------------------------------------------

def comp_names(pb, side):
    mal = MAL if (side == "tumour" and MAL in set(pb.compartment)) else f"{side}:epithelial" if side != "tumour" else L.EPI
    return {"endo": f"{side}:endothelial", "mal": mal, "T_NK": f"{side}:T_NK", "myeloid": f"{side}:myeloid"}


def genome_comps(study, side, comps, pts, calls):
    """Genome-wide pseudobulks (donor x compartment) of `comps` for `pts`; labels of `cancer_sc` (copy-number `calls` for the malignant split)."""
    K.install()
    us, _ = L.ST.units(study)
    us = [u for u in L.ST.majority_chemistry(us) if u[1] in set(pts)]
    sums, depth, strat, genes = {}, {}, {}, None
    for sample, patient, kind, load, chem, _ in us:
        X, _g, umi, comp, barcodes, sym = load()
        lab = L.B._labels(comp, kind, sample, barcodes, calls)
        genes = pd.Index(sym) if genes is None else genes
        for c in comps:
            i = np.flatnonzero(lab == c)
            if len(i):
                k = (patient, c)
                sums[k] = sums.get(k, 0) + np.asarray(X[i].sum(0)).ravel()
                depth[k] = depth.get(k, 0) + float(umi[i].sum())
                strat[k] = chem
    keys = sorted(sums)
    G = pd.DataFrame(np.vstack([sums[k] for k in keys]), columns=genes)
    G.columns = [L.canonical(str(g)) for g in G.columns]
    G = G.T.groupby(level=0).sum().T
    meta = pd.DataFrame({"donor": [k[0] for k in keys], "compartment": [k[1] for k in keys], "stratum": [strat[k] for k in keys],
                         "depth": [depth[k] for k in keys]})
    return pd.concat([meta, G.reset_index(drop=True)], axis=1)


def cohort_data(study, side, n_pseudo=0, full=False):
    """States, program-gene counts and cell counts of one cohort and side. Keys: study, side, pts (union of the patients of the compartments),
    stratum, cells (N), comps {key: dict(has, z, Zp, n_cells, cnt, lc)}, absent {key: reason}. `lc`: log2(CPM + 1) of the program genes
    (all genes if `full`); z / Zp / n_cells padded with NaN over `pts`."""
    from clocksync.config import limit_threads
    from clocksync.paired import projector, pseudo_waveforms, states
    limit_threads(1)
    pb = A.table(study)
    names = comp_names(pb, side)
    rows = {k: L._rows(pb, names[k]) for k in COMPS}
    absent = {k: f"{len(r)} patients with the pseudobulk" for k, r in rows.items() if len(r) < MIN_N}
    use = [k for k in COMPS if k not in absent]
    if not use:
        return dict(study=study, side=side, pts=[], comps={}, absent=absent)
    calls = None
    if names["mal"] == MAL:
        calls = L.new_calls(study) if study in L.NEW else L.ST.read_calls(study)
    G = genome_comps(study, side, [names[k] for k in use], sorted(set().union(*[set(rows[k].index) for k in use])), calls)
    pbk = pb[pb.compartment.str.startswith(side + ":")].drop_duplicates(["donor", "compartment"])
    parts, pts = {}, set()
    for k in use:
        g = G[G.compartment == names[k]].drop_duplicates("donor").set_index("donor")
        p = sorted(set(rows[k].index) & set(g.index))
        if len(p) < MIN_N:
            absent[k] = f"{len(p)} patients in the genome table"
        else:
            parts[k] = (g.loc[p], rows[k].loc[p], p)
            pts |= set(p)
    pts = sorted(pts)
    ix = {p: i for i, p in enumerate(pts)}
    stratum = pd.concat([r.stratum.astype(str) for k, (_, r, _) in parts.items()]).groupby(level=0).first().reindex(pts).to_numpy()
    comps = {}
    for k, (g, r, p) in parts.items():
        z = S.complex_state(states(r, projector(r, L.q.WF, True), L.q.WF, depth=r.depth.to_numpy(float))[0])
        genes = [c for c in g.columns if c not in ("compartment", "stratum", "depth")]
        X = g[genes].to_numpy(float)
        lc = np.log2(X / g.depth.to_numpy(float)[:, None] * 1e6 + 1)
        pg = sorted({x for v in GENES.values() for x in v})
        keep = np.ones(len(genes), bool) if full else np.isin(genes, pg)
        Zp = np.full((n_pseudo, len(pts)), np.nan + 0j)
        if n_pseudo:
            frac = g[genes].div(g.depth, axis=0)
            lm = np.log(frac.mean() + 1e-12)[(g[genes].mean() >= 1) | frac.columns.isin(list(L.q.WF.genes))]
            rng = np.random.default_rng(seed_of(COMPS.index(k), SIDES.index(side), 44, study=study))
            Zp[:, [ix[x] for x in p]] = np.array([S.complex_state(states(g, projector(g, w, True), w, depth=g.depth.to_numpy(float))[0])
                                                  for w in pseudo_waveforms(lm, n_pseudo, rng, L.q.WF)])
        zz = np.full(len(pts), np.nan + 0j)
        zz[[ix[x] for x in p]] = z
        nc = pd.Series(pbk[pbk.compartment == names[k]].set_index("donor").n_cells).reindex(pts).to_numpy(float)
        comps[k] = dict(z=zz, Zp=Zp, n_cells=nc, cnt=g.loc[p, [c for c in genes if c in pg]], lc=pd.DataFrame(lc[:, keep].astype(np.float32), index=p, columns=np.asarray(genes)[keep]))
    return dict(study=study, side=side, pts=pts, stratum=stratum, cells=pbk.groupby("donor").n_cells.sum().reindex(pts).to_numpy(float),
                comps=comps, absent=absent)


# ---- readouts --------------------------------------------------------------------------------------------------

def _z(x):
    return (x - x.mean()) / x.std(ddof=0)


def _y_of(d, prog, genes):
    """y of `prog` over d['pts'] (NaN where undefined) from `genes`; the score of A146 part 2. None if not definable."""
    ck, ab = PROGS[prog]
    c = d["comps"].get(ck)
    if c is None:
        return None
    pts = list(d["pts"])
    n, N = pd.Series(c["n_cells"], index=pts), pd.Series(d["cells"], index=pts)
    ok = c["lc"].index
    if ab:
        ok = ok[(N.reindex(ok) >= MIN_CELLS_ALL).to_numpy()]
    if len(ok) < 2:
        return None
    ab_z = _z(np.log2((n.reindex(ok) + .5) / (N.reindex(ok) + 1))) if ab else None
    if ab == "share":
        s = ab_z
    else:
        lc = c["lc"].loc[ok, genes].astype(float)
        z = (lc - lc.mean()) / lc.std(ddof=0).replace(0, np.nan)
        s = _z(z.mean(axis=1))
        if ab:
            X = np.column_stack([np.ones(len(ok)), ab_z.to_numpy()])
            s = pd.Series(s.to_numpy() - X @ np.linalg.lstsq(X, s.to_numpy(), rcond=None)[0], index=ok)
        s = _z(s)
    return s.reindex(pts).to_numpy(float)


def shuffle_finite(y, rng):
    y = y.copy()
    m = np.flatnonzero(np.isfinite(y))
    y[m] = y[rng.permutation(m)]
    return y


def readout(d, prog, placebo=False):
    """(y over d['pts'], reason). Eligibility: compartment present, >= half of the genes, median summed count >= 10."""
    ck, _ = PROGS[prog]
    c = d["comps"].get(ck)
    n = len(d["pts"])
    if c is None:
        return np.full(n, np.nan), f"{ck} pseudobulk: {d['absent'].get(ck, 'cohort without data')}"
    genes = [g for g in GENES[prog] if g in c["cnt"].columns and c["cnt"][g].notna().any()]
    if prog != "P8":
        if len(genes) < len(GENES[prog]) / 2:
            return np.full(n, np.nan), f"{len(genes)} of {len(GENES[prog])} genes in the cohort"
        med = float(c["cnt"][genes].sum(axis=1).median())
        if med < MIN_COUNT:
            return np.full(n, np.nan), f"median count {med:.1f} < {MIN_COUNT}"
    y = _y_of(d, prog, genes)
    if y is None:
        return np.full(n, np.nan), "fewer than 2 patients with >= 100 cells"
    if placebo:
        y = shuffle_finite(y, np.random.default_rng(seed_of(99, SIDES.index(d["side"]), list(PROGS).index(prog), study=d["study"])))
    return y, ""


def state_z(d, prog, state):
    c = d["comps"].get("endo" if state == "endo" else PROGS[prog][0])
    return (c["z"], c["Zp"], "") if c else (np.full(len(d["pts"]), np.nan + 0j), np.zeros((0, len(d["pts"]))), f"{state} state: pseudobulk absent")


def _centre(v, stratum):
    return v - pd.Series(v).groupby(stratum).transform("mean").to_numpy()


def make_unit(d, y, z):
    """Centred (within stratum, over patients with y and the state) y and z; (None, reason) if < MIN_N patients or no spread."""
    m = np.isfinite(y) & np.isfinite(z)
    if m.sum() < MIN_N:
        return None, f"{int(m.sum())} patients with the readout and the state"
    st = d["stratum"][m]
    yc, zc = _centre(y[m], st), _centre(z[m], st)
    if not (np.sum(yc ** 2) > 0 and np.sum(abs(zc) ** 2) > 0):
        return None, "no spread in y or z"
    return dict(n=int(m.sum()), mask=m, stratum=st, y=yc, z=zc, den=float(np.sqrt(np.sum(yc ** 2) * np.sum(abs(zc) ** 2)))), ""


def units(d, placebo=False):
    """{(prog, state): (unit or None, reason)} of one cohort."""
    out = {}
    for prog in PROGS:
        y, why = readout(d, prog, placebo)
        for st in STATES:
            if why:
                out[(prog, st)] = (None, why)
                continue
            z, _, why2 = state_z(d, prog, st)
            out[(prog, st)] = make_unit(d, y, z) if not why2 else (None, why2)
    return out


stat = A29.stat
cross = A29.cross
peak_h = A29.peak_h


def shared_nulls(d, un, n_perm):
    """Permutation null of every unit of the cohort from one shared ordering of all its patients per iteration."""
    R = np.random.default_rng(seed_of(SIDES.index(d["side"]), study=d["study"])).random((n_perm, len(d["pts"])))
    out = {}
    for k, (u, _) in un.items():
        if u:
            perms = np.argsort(R[:, np.flatnonzero(u["mask"])], axis=1)
            out[k] = np.abs((u["y"][perms] * np.conj(u["z"])[None, :]).sum(1)) / u["den"]
    return out


def pseudo_stats(d, u, Zp):
    return np.array([stat(u, z=_centre(z[u["mask"]], u["stratum"])) for z in Zp])


# ---- power -----------------------------------------------------------------------------------------------------

def power_one(study, side, n_perm, placebo):
    d = cohort_data(study, side)
    base = dict(study=study, cancer=A.cancer(study), side=side)
    if not d["pts"]:
        return [dict(base, program=p, state=s, patients=0, eligible=False, reason="; ".join(d["absent"].values()) or "no compartment")
                for p in PROGS for s in STATES], {}
    un = units(d, placebo)
    elig = [dict(base, program=p, state=s, patients=u["n"] if u else 0, eligible=u is not None, reason=why) for (p, s), (u, why) in un.items()]
    nulls = {k: (un[k][0]["n"], v) for k, v in shared_nulls(d, un, n_perm).items()}
    return elig, nulls


def power_table(nulls, n_perm):
    """nulls: {(side, cohort): {(prog, state): (n, null)}} -> per-unit pooled nulls, power rows, omnibus rows."""
    rows, om, pooled = [], [], {}
    for side in SIDES:
        for st in STATES:
            for prog in PROGS:
                cs = [c for (sd, c), v in nulls.items() if sd == side and (prog, st) in v]
                ns = [nulls[(side, c)][(prog, st)][0] for c in cs]
                nul = A.pool_null(ns, [nulls[(side, c)][(prog, st)][1] for c in cs]) if cs else np.array([np.nan])
                pooled[(side, st, prog)] = nul
                sd_, n = float(nul.std()), int(sum(ns))
                rows.append(dict(side=side, state=st, program=prog, hypothesis="H246.3" if side == "tumour" else "descriptive", patients=n,
                                 cohorts=len(cs), cohort_list=";".join(cs), null_mean=float(nul.mean()), null_sd=sd_, mde=MDE_K * sd_,
                                 n_perm=n_perm, scored=bool(n >= MIN_PATIENTS and MDE_K * sd_ <= MAX_MDE)))
    P = pd.DataFrame(rows)
    for side in SIDES:
        for st in STATES:
            sc = P[(P.side == side) & (P.state == st) & P.scored].program.tolist()
            T = np.mean([(pooled[(side, st, p)] - pooled[(side, st, p)].mean()) / pooled[(side, st, p)].std() for p in sc], 0) if sc else np.array([np.nan])
            om.append(dict(side=side, state=st, program="omnibus", hypothesis=STATES[st][0] if side == "tumour" else "descriptive",
                           patients=np.nan, cohorts=np.nan, cohort_list=";".join(sc),
                           null_mean=float(T.mean()), null_sd=float(T.std()), mde=MDE_K * float(T.std()), n_perm=n_perm, n_scored_programs=len(sc),
                           scored=len(sc) >= MIN_OMNIBUS))
    return pd.concat([P, pd.DataFrame(om)], ignore_index=True)


def stage_power(args):
    cohorts = A.cohort_list(args.smoke)
    res = Parallel(n_jobs=args.jobs)(delayed(power_one)(s, sd, args.n_perm, args.placebo) for s in cohorts for sd in SIDES)
    pd.DataFrame([r for e, _ in res for r in e]).to_csv(args.out / "eligibility.csv", index=False)
    keys = [(sd, s) for s in cohorts for sd in SIDES]
    P = power_table({k: n for k, (_, n) in zip(keys, res)}, args.n_perm)
    P.to_csv(args.out / "power.csv", index=False)
    print(P.drop(columns="cohort_list").round(4).to_string(), flush=True)


# ---- tests (registered statistic; not run on real data before the power note) ----------------------------------

def pool_pseudo(P, w):
    """Weighted pooled pseudo r per set; sets undefined in a cohort are re-weighted over the cohorts that have them (NaN kept for alignment)."""
    Wm = np.where(np.isfinite(P), w[:, None], 0.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.nansum(P * Wm, 0) / Wm.sum(0)


def pool(cells):
    """cells: [dict(n, r, c, null, pseudo)] -> pooled observed r, null, pseudo, complex cross-correlation."""
    w = np.array([x["n"] for x in cells], float)
    return dict(n=int(w.sum()), cohorts=len(cells), obs=float(np.sum(w * [x["r"] for x in cells]) / w.sum()), c=np.sum(w * [x["c"] for x in cells]) / w.sum(),
                nul=A.pool_null(w, [x["null"] for x in cells]), pse=pool_pseudo(np.array([x["pseudo"] for x in cells]), w))


def omnibus(pools):
    """T observed, its permutation null and pseudo values from the pools of the programs {prog: pool}."""
    z = lambda x, p: (x - p["nul"].mean()) / p["nul"].std()  # noqa: E731
    obs = float(np.mean([z(p["obs"], p) for p in pools.values()]))
    nul = np.mean([z(p["nul"], p) for p in pools.values()], 0)
    pse = np.mean([z(p["pse"], p) for p in pools.values()], 0)
    return obs, nul, pse[np.isfinite(pse)]


def circ_corr(a, b):
    a, b = np.asarray(a) * 2 * np.pi / 24, np.asarray(b) * 2 * np.pi / 24
    ma, mb = np.angle(np.exp(1j * a).mean()), np.angle(np.exp(1j * b).mean())
    sa, sb = np.sin(a - ma), np.sin(b - mb)
    return float(np.sum(sa * sb) / np.sqrt(np.sum(sa ** 2) * np.sum(sb ** 2)))


def matched_sets(data, units_by, tumour_pools, args):
    """Percentile of each program (tumour side, both states) among N_MATCH expression-matched gene sets of its compartment used as y."""
    rows, pct = [], {}
    excl = {g for v in GENES.values() for g in v} | {L.canonical(g) for g in L.q.WF.genes}
    for ck in ("mal", "T_NK", "myeloid"):
        ds = [d for (sd, _), d in data.items() if sd == "tumour" and ck in d["comps"]]
        if not ds:
            continue
        names = sorted(set.intersection(*[set(d["comps"][ck]["lc"].columns) for d in ds]))
        w = np.array([len(d["comps"][ck]["lc"]) for d in ds], float)
        mean = pd.Series({g: float(np.sum(w * [d["comps"][ck]["lc"][g].to_numpy(float).mean() for d in ds]) / w.sum()) for g in names})
        for prog in [p for p, (c, _) in PROGS.items() if c == ck and p != "P8"]:
            pg = [g for g in GENES[prog] if g in mean.index]
            cand = {g: (mean.drop([x for x in excl if x in mean.index]) - mean[g]).abs().sort_values().index[:N_POOL].tolist() for g in pg}
            rng = np.random.default_rng([SEED, 7, list(PROGS).index(prog)])
            sets = []
            for _ in range(N_MATCH):
                chosen = []
                for g in pg:
                    chosen.append(rng.choice([x for x in cand[g] if x not in chosen]))
                sets.append(chosen)
            for st in STATES:
                key = ("tumour", st, prog)
                if key not in tumour_pools:
                    continue
                cs = [(d, units_by[("tumour", d["study"])][(prog, st)][0]) for d in ds]
                cs = [(d, u) for d, u in cs if u]
                rs = []
                for j, sel in enumerate(sets):
                    rr, ww = [], []
                    for d, u in cs:
                        y = _y_of(d, prog, sel)
                        if y is None:
                            continue
                        if args.placebo:
                            y = shuffle_finite(y, np.random.default_rng([SEED, 98, j, sum(map(ord, d["study"]))]))
                        yc = _centre(y[u["mask"]], u["stratum"])
                        if np.isfinite(yc).all() and np.sum(yc ** 2) > 0:
                            rr.append(stat(u, y=yc))
                            ww.append(u["n"])
                    rs.append(float(np.average(rr, weights=ww)) if rr else np.nan)
                    rows.append(dict(program=prog, state=st, set=j, genes=";".join(sel), statistic=rs[-1]))
                rs = np.array(rs)
                pct[(prog, st)] = float(np.mean(rs[np.isfinite(rs)] <= tumour_pools[key]["obs"]) * 100)
    return rows, pct


def stage_tests(args):
    assert args.placebo or args.registered_run, "tests on real data only after the power note is committed (--registered-run)"
    pf = next((f for f in (args.out / "power.csv", OUT / "power.csv") if f.exists()), None)
    pw = pd.read_csv(pf) if pf else None
    cohorts = A.cohort_list(args.smoke)
    data = Parallel(n_jobs=args.jobs)(delayed(cohort_data)(s, sd, args.n_pseudo, True) for s in cohorts for sd in SIDES)
    data = {(d["side"], d["study"]): d for d in data if d["pts"]}
    units_by, nulls, cell, per, pseudo_rows = {}, {}, {}, [], []
    for (side, s), d in data.items():
        units_by[(side, s)] = un = units(d, args.placebo)
        nulls[(side, s)] = shared_nulls(d, un, args.n_perm)
        for (prog, st), (u, _) in un.items():
            if u:
                _, Zp, _ = state_z(d, prog, st)
                cell.setdefault((side, st, prog), {})[s] = dict(n=u["n"], r=stat(u), c=cross(u["y"], u["z"]) / u["den"], null=nulls[(side, s)][(prog, st)],
                                                                pseudo=pseudo_stats(d, u, Zp), cancer=A.cancer(s))
                per.append(dict(side=side, state=st, program=prog, cohort=s, cancer=A.cancer(s), n=u["n"], statistic=cell[(side, st, prog)][s]["r"]))

    def scored(side, st, prog):
        if pw is None or args.force:
            return True
        r = pw[(pw.side == side) & (pw.state == st) & (pw.program == prog)]
        return bool(len(r) and r.scored.iloc[0])

    tests, pooled_rows, tpools = [], [], {}
    for side in SIDES:
        for st in STATES:
            progs = [p for p in PROGS if (side, st, p) in cell and scored(side, st, p)]
            pools = {p: pool(list(cell[(side, st, p)].values())) for p in progs}
            for p, pl in pools.items():
                tpools[(side, st, p)] = pl
                pse = pl["pse"][np.isfinite(pl["pse"])]
                tests.append(dict(side=side, state=st, level="program", program=p, hypothesis="H246.3" if side == "tumour" else "descriptive",
                                  patients=pl["n"], cohorts=pl["cohorts"], statistic=pl["obs"], null_mean=float(pl["nul"].mean()), null_sd=float(pl["nul"].std()),
                                  p_pairing=float((1 + np.sum(pl["nul"] >= pl["obs"])) / (1 + len(pl["nul"]))), pseudo_median=float(np.median(pse)),
                                  pseudo_q95=float(np.quantile(pse, .95)), n_sets=len(pse), p_clock=float((1 + np.sum(pse >= pl["obs"])) / (1 + len(pse))),
                                  implied_peak_h=peak_h(pl["c"]), a146_hour=P146.HUMAN_H[p]))
                pseudo_rows += [dict(side=side, state=st, program=p, set=i, value=v) for i, v in enumerate(pl["pse"])]
            if len(progs) >= MIN_OMNIBUS:
                T, nul, pse = omnibus(pools)
                tests.append(dict(side=side, state=st, level="omnibus", program="omnibus", hypothesis=STATES[st][0] if side == "tumour" else "descriptive",
                                  patients=sum(pl["n"] for pl in pools.values()), cohorts=np.nan, n_programs=len(progs), statistic=T, null_mean=float(nul.mean()),
                                  null_sd=float(nul.std()), p_pairing=float((1 + np.sum(nul >= T)) / (1 + len(nul))), pseudo_median=float(np.median(pse)),
                                  pseudo_q95=float(np.quantile(pse, .95)), n_sets=len(pse), p_clock=float((1 + np.sum(pse >= T)) / (1 + len(pse)))))
                pseudo_rows += [dict(side=side, state=st, program="omnibus", set=i, value=v) for i, v in enumerate(pse)]
            if len(progs) >= 3:
                pk = [peak_h(pools[p]["c"]) for p in progs]
                h = [P146.HUMAN_H[p] for p in progs]
                tests.append(dict(side=side, state=st, level="peak_agreement", program="circ_corr_peak_vs_A146_hour", hypothesis="descriptive",
                                  n_programs=len(progs), statistic=circ_corr(pk, h), null_mean=circ_corr([-x for x in pk], h)))  # null_mean: mirrored
            # descriptives: cancer types with >= 12 patients, breast cohorts alone
            subsets = {ca: [s for (sd, s) in data if sd == side and A.cancer(s) == ca] for ca in sorted({A.cancer(s) for (sd, s) in data if sd == side})}
            subsets["breast_cohorts"] = subsets.get("breast", [])
            for name, ss in subsets.items():
                sp = {}
                for p in progs:
                    cs = [x for s, x in cell[(side, st, p)].items() if s in ss]
                    if cs and sum(x["n"] for x in cs) >= (1 if name == "breast_cohorts" else MIN_CANCER):
                        sp[p] = pool(cs)
                        pse = sp[p]["pse"][np.isfinite(sp[p]["pse"])]
                        pooled_rows.append(dict(side=side, state=st, subset=name, program=p, patients=sp[p]["n"], cohorts=sp[p]["cohorts"], statistic=sp[p]["obs"],
                                                null_mean=float(sp[p]["nul"].mean()), null_sd=float(sp[p]["nul"].std()),
                                                pseudo_median=float(np.median(pse)) if len(pse) else np.nan, implied_peak_h=peak_h(sp[p]["c"])))
                if len(sp) >= 3:
                    T, nul, pse = omnibus(sp)
                    pooled_rows.append(dict(side=side, state=st, subset=name, program="omnibus", patients=sum(x["n"] for x in sp.values()), statistic=T, null_mean=float(nul.mean()),
                                            null_sd=float(nul.std()), pseudo_median=float(np.median(pse)) if len(pse) else np.nan))
    Tt = pd.DataFrame(tests)
    Tt["p_holm"] = np.nan
    for lv, grp in (("omnibus", ["side"]), ("program", ["side", "state"])):
        for _, g in Tt[(Tt.level == lv) & (Tt.side == "tumour")].groupby(grp):
            Tt.loc[g.index, "p_holm"] = A.holm(g.p_pairing)
    m = (Tt.side == "tumour") & Tt.level.isin(["omnibus", "program"])
    Tt["pass"] = np.where(m, (Tt.p_holm < 0.05) & (Tt.statistic > Tt.pseudo_q95), np.nan)
    mrows, pct = matched_sets(data, units_by, tpools, args)
    Tt["percentile_among_matched_sets"] = [pct.get((p, s), np.nan) if (sd == "tumour" and lv == "program") else np.nan
                                           for sd, s, p, lv in zip(Tt.side, Tt.state, Tt.program, Tt.level)]
    Tt.to_csv(args.out / "tests.csv", index=False)
    pd.DataFrame(pooled_rows).to_csv(args.out / "pooled.csv", index=False)
    pd.DataFrame(per).to_csv(args.out / "per_cohort.csv", index=False)
    pd.DataFrame(pseudo_rows).to_csv(args.out / "pseudo.csv", index=False)
    pd.DataFrame(mrows).to_csv(args.out / "matched_sets.csv", index=False)
    print(Tt.drop(columns=["a146_hour"]).round(4).to_string(), flush=True)


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--stage", choices=["power", "tests"], required=True)
    ap.add_argument("--n-perm", type=int, default=10_000)
    ap.add_argument("--n-pseudo", type=int, default=200)
    ap.add_argument("--force", action="store_true", help="tests: also run programs flagged not scored")
    ap.add_argument("--smoke", action="store_true", help="3 cohorts, 200 permutations, 5 pseudo sets")
    ap.add_argument("--placebo", action="store_true", help="shuffle every y once within cohort (code test); writes to output/placebo/")
    ap.add_argument("--registered-run", action="store_true", help="tests on real data (only after the power note is committed)")
    ap.add_argument("--out", type=Path, default=None)
    args = ap.parse_args(argv)
    assert args.jobs <= 3, "at most 3 workers for A246"
    cli.setup(args)
    args.out = args.out or OUT / ("placebo" if args.placebo else "smoke" if args.smoke else "")
    args.out.mkdir(parents=True, exist_ok=True)
    if args.smoke:
        args.n_perm, args.n_pseudo = 200, 5
    {"power": stage_power, "tests": stage_tests}[args.stage](args)


if __name__ == "__main__":
    main()
