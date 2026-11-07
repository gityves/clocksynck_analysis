"""A232: do organs with more proliferation run at a different circadian phase? GTEx, organ-map offsets vs a G2M score.

Registered in paper/docs/registered_experiments.md (A232). Two stages:
  power  (no association between an offset and a proliferation score is computed, printed or written)
         tissues per cohort, per-tissue offset bootstrap SE (donors resampled jointly over tissues, host refitted),
         MDE_80 of Spearman rho at n tissues, H232.3 permutation null SD (pooled) and MDE_80 = 2.49 x null SD
         -> output/power.csv, tissue_offsets.csv, power_tissues.csv, draws/.  The H232.3 null is label-independent: the
         tertile labels are shuffled once before `Dataset.contrast`, so no observed top-vs-bottom offset is computed.
  score  H232.1 / H232.2 (Spearman rho offset vs tissue G2M score, 10,000 score permutations over tissues, weighted
         slope), H232.3 (pooled within-tissue top-vs-bottom tertile offset), pseudo-clock and random-gene-set
         companions (200 each), Holm over H232.1-3 -> output/score_*.csv.  NOT to be run before the power note is dated.

    python compute_a232.py --stage power [--subsample N] [--jobs 8] [--n-boot 200] [--n-perm 1000]
Controlled access: only the public sex / Hardy class and the protected time-of-death *filter* (fast cohort keeps donors
with a time of death, as the organ map) are used; no per-donor value is written. Per-tissue offsets, SEs and scores are
aggregates.

Choices not fixed by the registration
- Counts: library-size depth (`depth="library_size"`), as A164 / A208 and the organ map; cohorts fast = Hardy 1-2 with
  time of death, vent = Hardy 0; cultured cells out; first sample of a donor x tissue kept. Data read directly from
  the labelled h5ad (`common.gtex.LABELLED`, same filters as `common.gtex.load`) because genome-wide columns are needed.
- Tissues: >= 50 donors in the cohort (after the filters above). Anchors = the 15 GTEx organ-map tissues present in the
  cohort (also when < 50 donors; they are then anchors, not targets). Offsets: `Dataset.fit(gate=False, n_boot=0,
  targets=eligible tissues)` (circular-mean centred over the targets, anchors on their leave-one-out host, as A164).
- G2M set: Hallmark G2M checkpoint (MSigDB v2024.1 in $SCCC_DATA/reference; the A97 run used v2026.1, not available),
  14 clock genes removed (`figure_5/scripts/cancer/_shared.hallmark_g2m`), genes absent from GTEx dropped. Score =
  mean over genes of the z of log2(CPM + 1) (CPM = counts / library size x 1e6), z over all samples of the cohort
  (all non-cultured tissues, ddof 0). Tissue score = median over its donors. MKI67 / MYC: median log2(CPM + 1).
- Bootstrap SE of an offset: 200 donor resamplings (jointly over tissues, donors copied as id@draw), `Dataset` and host
  refitted, offsets re-centred in each replicate; SE = SD of the offset deviation from the point estimate wrapped to
  [-12, 12) h.
- MDE of rho: (a) Bonett-Wright approximation atanh(rho) = 2.8016 sqrt((1 + rho^2 / 2) / (n - 3)) (80 % power, two-sided
  0.05) and (b) simulation (bivariate normal, Spearman t-approximation test, 4,000 draws per rho on a 0.005 grid).
  Both for exactly measured offsets; `reliability` = 1 - mean(SE^2) / var(offsets) and the MDE divided by
  sqrt(reliability) are reported as the attenuated version.
- H232.3: tertiles within cohort x tissue (top: score >= 2/3 quantile = case; bottom: <= 1/3 quantile = control), host fitted
  on the tertile donors (`Dataset.contrast`, blocks="none", n_boot=0), 1,000 label permutations per tissue (the
  `Contrast.null` draws, host fixed), independent across tissues. Pooling as A208: tissue weights 1 / null SD^2 within
  a cohort, cohort weights 1 / SD^2 of the cohort pooled null, fixed by the null (power stage) and reused in score.
- Random gene sets (score): candidates = genes with mean count >= 1 in the cohort, minus clock and G2M genes; each G2M
  gene is replaced by an unused candidate of the same decile of mean log2(CPM + 1). Pseudo-clock sets: A208 rule
  (`pseudo_waveforms`, log mean fraction of the genes with mean count >= 1, all cohort samples), rng seed 232 per cohort.
  Pseudo-clock sets are run through the whole pipeline: offsets per tissue (fit), the same scores, H232.3 pooled offset.
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[2]
A164 = PAPER / "exploratory/A164_gtex_metabolic_matched/scripts"
sys.path.insert(0, str(A164))
sys.path.insert(0, str(PAPER / "results/figure_5/scripts/cancer"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402
from scipy import stats  # noqa: E402

import compute_a164 as A  # noqa: E402  (say, seed_of, p2, run_chunks, MDE_K)
import _shared as S  # noqa: E402  (gtex_cues: TISSUES = organ-map anchors)
import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.genes import gene_columns  # noqa: E402
from clocksync.paired import pseudo_waveforms  # noqa: E402
from common import cli, gtex  # noqa: E402

import importlib.util  # noqa: E402

_spec = importlib.util.spec_from_file_location("cancer_shared", PAPER / "results/figure_5/scripts/cancer/_shared.py")
CANCER = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(CANCER)

OUT = HERE.parent / "output"
COH = {"fast": ((1, 2), True), "vent": ((0,), False)}  # hardy, require time of death (the filter only)
LEVELS = ("control", "case")  # bottom, top tertile
MIN_DONORS, N_PERM_SCORE, SEED = 50, 10000, 232
Z_MDE = 2.8016  # z(0.975) + z(0.80)
say, MDE_K = A.say, A.MDE_K
wrap = lambda x: (np.asarray(x, float) + 12.0) % 24.0 - 12.0


# ---- data ----------------------------------------------------------------------------------------------------------------
class Coh:
    """One cohort: sample table `meta` (donor, tissue, depth), counts `Xc` and log2(CPM + 1) `Xl` (CSC, rows = meta,
    columns = pool genes and needed genes), gene names, positions of the clock / G2M / MKI67 / MYC columns."""


def load_cohort(name: str, subsample: int | None, seed: int) -> Coh:
    import anndata as ad

    hardy, need_tod = COH[name]
    a = ad.read_h5ad(gtex.LABELLED, backed="r")
    try:
        o = a.obs
        keep = o.DTHHRDY.isin(list(hardy)).to_numpy() & ~o.tissue.astype(str).str.startswith("Cells -").to_numpy()
        if need_tod:
            keep &= o.tod.notna().to_numpy()
        idx = np.flatnonzero(keep)
        df = pd.DataFrame({"donor": o.subject.astype(str).to_numpy()[idx], "tissue": o.tissue.astype(str).to_numpy()[idx]})
        first = ~df.duplicated(["donor", "tissue"]).to_numpy()
        if subsample:
            donors = np.random.default_rng(seed).permutation(sorted(df.donor.unique()))[:subsample]
            first &= df.donor.isin(donors).to_numpy()
        idx, df = idx[first], df[first].reset_index(drop=True)
        lib = o.library_size.to_numpy(float)[idx]
        X = a.X[idx].tocsr()
        var = np.asarray(a.var_names)
    finally:
        a.file.close()
    c = Coh()
    c.name = name
    c.meta = df.assign(depth=lib)
    cnt = np.asarray(X.mean(0)).ravel()
    clock = list(cs.gtex().genes)
    g2m = CANCER.hallmark_g2m()
    named = clock + g2m + ["MKI67", "MYC"]
    full = gene_columns(var, named, missing="ignore")
    c.dropped_g2m = [g for g, k in zip(named, full) if k < 0 and g in g2m]
    pool = np.flatnonzero(cnt >= 1)
    cols = np.unique(np.r_[pool, full[full >= 0]])
    where = {int(k): j for j, k in enumerate(cols)}
    c.names = var[cols]
    c.Xc = X[:, cols].tocsc().astype(np.float32)
    del X
    Xl = c.Xc.multiply(1e6 / lib[:, None]).tocsc()
    Xl.data = np.log2(Xl.data + 1).astype(np.float32)
    c.Xl = Xl.tocsc()
    pos = lambda genes: [where[int(k)] for k in gene_columns(var, genes, missing="ignore") if k >= 0]
    c.clock_names = clock
    c.clock = [where[int(k)] for k in gene_columns(var, clock, missing="raise")]
    c.g2m = pos(g2m)
    c.mki67, c.myc = pos(["MKI67"])[0], pos(["MYC"])[0]
    c.pool = np.flatnonzero(cnt[cols] >= 1)
    return c


def zmean(c: Coh, cols) -> np.ndarray:
    """Per-sample mean over genes of the z (over all samples of the cohort) of log2(CPM + 1)."""
    M = c.Xl[:, list(cols)].toarray().astype(float)
    sd = M.std(0)
    ok = sd > 0
    return ((M[:, ok] - M[:, ok].mean(0)) / sd[ok]).mean(1)


def pseudobulk(c: Coh, cols, names, rows=None) -> pd.DataFrame:
    m = c.meta if rows is None else c.meta[rows]
    X = c.Xc[:, list(cols)]
    X = (X if rows is None else X[np.flatnonzero(rows)]).toarray().astype(float)
    pb = pd.DataFrame(X, columns=list(names))
    pb.insert(0, "donor", m.donor.to_numpy())
    pb.insert(1, "compartment", m.tissue.to_numpy())
    pb.insert(2, "depth", m.depth.to_numpy())
    return pb.reset_index(drop=True)


def log_mean_of(c: Coh) -> pd.Series:
    """Log of the mean count fraction per gene (A208 rule), genes with mean count >= 1; clock genes under their names."""
    frac = np.asarray(c.Xc.multiply(1 / c.meta.depth.to_numpy()[:, None]).mean(0)).ravel()
    lm = pd.Series(np.log(frac + 1e-12), index=c.names)
    lm = lm.iloc[c.pool]
    lm = lm[~lm.index.duplicated()]
    for g, j in zip(c.clock_names, c.clock):
        lm[g] = np.log(frac[j] + 1e-12)
    return lm


def setup(name, subsample, seed, scramble=False):
    """(Coh, clock pseudobulk of eligible tissues + anchors, eligible tissues, anchors, tertile labels, tissue scores)."""
    c = load_cohort(name, subsample, seed)
    n = c.meta.groupby("tissue").donor.nunique()
    min_d = 15 if subsample else MIN_DONORS
    targets = sorted(n.index[n >= min_d])
    anchors = [t for t in S.TISSUES if t in set(c.meta.tissue)]
    rows = c.meta.tissue.isin(set(targets) | set(anchors)).to_numpy()
    pb = pseudobulk(c, c.clock, c.clock_names, rows)
    score = zmean(c, c.g2m)
    if scramble:  # code test only: the G2M score is replaced by noise, so no offset-proliferation association is formed
        score = np.random.default_rng(seed + 7).standard_normal(len(score))
    labs, tscore = {}, {}
    for t in targets:
        i = np.flatnonzero((c.meta.tissue == t).to_numpy())
        s = pd.Series(score[i], index=c.meta.donor.to_numpy()[i])
        lo, hi = s.quantile([1 / 3, 2 / 3])
        labs[t] = pd.concat([pd.Series("control", index=s.index[s <= lo]), pd.Series("case", index=s.index[s >= hi])])
        tscore[t] = dict(g2m=float(s.median()),
                         mki67=float(np.median(c.Xl[:, c.mki67].toarray().ravel()[i])),
                         myc=float(np.median(c.Xl[:, c.myc].toarray().ravel()[i])))
    c.score = score
    return c, pb, targets, anchors, labs, pd.DataFrame(tscore).T


# ---- offsets, bootstrap -----------------------------------------------------------------------------------------------
def offsets(pb, anchors, targets, waveform="human") -> pd.Series:
    ds = cs.Dataset(pb, anchors=[a for a in anchors if a in set(pb.compartment)], waveform=waveform)
    return ds.fit(n_boot=0, gate=False, targets=targets).table.offset_h.reindex(targets)


def boot_chunk(pb, anchors, targets, seeds):
    """Offsets of donor-bootstrap replicates (donors drawn with replacement, copied as id@draw), host refitted."""
    limit_threads(1)
    warnings.filterwarnings("ignore")
    donors = np.array(sorted(pb.donor.unique()))
    out = []
    for sd in seeds:
        draw = np.random.default_rng(sd).choice(len(donors), len(donors))
        m = pd.DataFrame({"orig": donors[draw], "copy": [f"{d}@{j}" for j, d in enumerate(donors[draw])]})
        rep = pb.merge(m, left_on="donor", right_on="orig")
        rep["donor"] = rep["copy"]
        rep = rep.drop(columns=["orig", "copy"]).reset_index(drop=True)
        try:
            out.append(offsets(rep, anchors, targets).to_numpy(float))
        except Exception:
            out.append(np.full(len(targets), np.nan))
    return np.array(out)


# ---- H232.3 contrasts --------------------------------------------------------------------------------------------------
def contrast_task(ds, t, lab, n_perm, seed, want_obs):
    limit_threads(1)
    warnings.filterwarnings("ignore")
    try:
        r = ds.contrast(lab, levels=LEVELS, targets=[t], n_perm=n_perm, n_boot=0, blocks="none", seed=seed)
        obs = float(r.table.loc[t, "offset_diff_h"]) if want_obs and t in r.table.index else np.nan
        nul = r.null["offset_diff_h"][t].to_numpy(float) if r.null is not None else np.array([])
    except Exception:
        obs, nul = np.nan, np.full(n_perm, np.nan)
    return obs, nul


def contrasts(ds, targets, labs, n_perm, jobs, tag, preshuffle):
    """(observed offset_diff_h (targets,) or None, null draws (n_perm, targets)) per tissue contrast. `preshuffle`:
    labels shuffled before the call, so the null is the same distribution and no observed top-vs-bottom value exists."""
    tasks = []
    for t in targets:
        sd = A.seed_of("A232", tag, t)
        lab = labs[t]
        if preshuffle:
            lab = pd.Series(np.random.default_rng(sd).permutation(lab.to_numpy()), index=lab.index)
        tasks.append(delayed(contrast_task)(ds, t, lab, n_perm, sd, not preshuffle))
    res = Parallel(n_jobs=jobs)(tasks)
    obs = None if preshuffle else np.array([r[0] for r in res])
    return obs, np.column_stack([r[1] for r in res]) if n_perm else None


# ---- pooling (A208) --------------------------------------------------------------------------------------------------------
def wmean(x, w):
    ok = np.isfinite(x)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(ok, x, 0) @ w / (ok.astype(float) @ w)


def tissue_weights(null):
    sd = np.nanstd(null, axis=0, ddof=1)
    w = np.where(np.isfinite(sd) & (sd > 0), 1 / sd ** 2, 0.0)
    return w / w.sum()


def pool_nulls(nulls: dict):
    """(pooled null, tissue weights {cohort}, cohort weights, cohort SDs)."""
    tw = {c: tissue_weights(a) for c, a in nulls.items()}
    coh = {c: wmean(a, tw[c]) for c, a in nulls.items()}
    sd = {c: float(np.nanstd(v, ddof=1)) for c, v in coh.items()}
    cw = {c: 1 / s ** 2 for c, s in sd.items()}
    tot = sum(cw.values())
    cw = {c: v / tot for c, v in cw.items()}
    return sum(cw[c] * coh[c] for c in cw), tw, cw, sd


# ---- rho MDE -----------------------------------------------------------------------------------------------------------------
def mde_rho_approx(n: int) -> float:
    from scipy.optimize import brentq
    return float(brentq(lambda r: np.arctanh(r) - Z_MDE * np.sqrt((1 + r ** 2 / 2) / (n - 3)), 1e-6, 0.9999))


def mde_rho_sim(n: int, n_sim: int = 4000, seed: int = SEED) -> float:
    rng = np.random.default_rng(seed)
    crit = stats.t.ppf(0.975, n - 2)
    grid = np.arange(0.005, 1.0, 0.005)
    for r in grid:
        x = rng.standard_normal((n_sim, n))
        y = r * x + np.sqrt(1 - r ** 2) * rng.standard_normal((n_sim, n))
        rk = lambda m: m.argsort(1).argsort(1).astype(float)
        rho = np.array([np.corrcoef(a, b)[0, 1] for a, b in zip(rk(x), rk(y))])
        if np.mean(np.abs(rho) * np.sqrt((n - 2) / (1 - rho ** 2)) > crit) >= 0.8:
            return float(r)
    return np.nan


# ---- stages ------------------------------------------------------------------------------------------------------------------
def bootstrap_offsets(args, out, quick, name, pb, anchors, targets):
    n_boot = args.n_boot or (10 if quick else 200)
    seed0 = A.seed_of("A232boot", name, args.seed)
    f = out / "draws" / f"boot_offsets_{name}.npy"
    if not quick and f.exists():
        return np.load(f)
    arr = A.run_chunks(boot_chunk, n_boot, args.jobs, seed0, f"boot {name}", pb, anchors, targets)
    if not quick:
        np.save(f, arr)
    return arr


def boot_se(point: pd.Series, arr: np.ndarray):
    dev = wrap(arr - point.to_numpy(float)[None, :])
    return np.nanstd(dev, axis=0, ddof=1), np.isfinite(dev).sum(0)


def stage_power(args, out, quick):
    n_perm = args.n_perm or (20 if quick else 1000)
    nulls, trows, prow, orow = {}, [], [], []
    for name in COH:
        c, pb, targets, anchors, labs, _ = setup(name, args.subsample, args.seed)
        say(f"{name}: {c.meta.donor.nunique()} donors, {len(c.meta)} samples, {len(targets)} tissues, "
            f"{len(anchors)} anchors, G2M genes {len(c.g2m)} (dropped: {c.dropped_g2m})")
        off = offsets(pb, anchors, targets)
        arr = bootstrap_offsets(args, out, quick, name, pb, anchors, targets)
        se, nok = boot_se(off, arr)
        nd = pb[pb.compartment.isin(targets)].groupby("compartment").donor.nunique()
        for t, o, s, k in zip(targets, off, se, nok):
            orow.append(dict(cohort=name, tissue=t, n_donors=int(nd[t]), n_top=int((labs[t] == "case").sum()),
                             n_bottom=int((labs[t] == "control").sum()), offset_h=o, boot_se_h=s, n_boot_ok=int(k)))
        ds = cs.Dataset(pb, anchors=anchors, name=f"gtex_a232_{name}")
        _, nul = contrasts(ds, targets, labs, n_perm, args.jobs, f"power{name}", preshuffle=True)
        nulls[name] = nul
        if not quick:
            np.save(out / "draws" / f"null_h232_3_{name}.npy", nul)
        say(f"{name}: null done")
    O = pd.DataFrame(orow)
    O.to_csv(out / "tissue_offsets.csv", index=False)
    pooled, tw, cw, sd = pool_nulls(nulls)
    for name in COH:
        O_c = O[O.cohort == name]
        for t, w, s in zip(O_c.tissue, tw[name], np.nanstd(nulls[name], axis=0, ddof=1)):
            trows.append(dict(cohort=name, tissue=t, null_sd=float(s), tissue_weight=float(w)))
    pd.DataFrame(trows).to_csv(out / "power_tissues.csv", index=False)
    for name in COH:
        O_c = O[O.cohort == name]
        n = len(O_c)
        var_off = float(np.var(wrap(O_c.offset_h - O_c.offset_h.mean()), ddof=1))
        rel = float(np.clip(1 - np.mean(O_c.boot_se_h ** 2) / var_off, 0.05, 1))
        m_a, m_s = mde_rho_approx(n), mde_rho_sim(n)
        prow.append(dict(cohort=name, n_tissues=n, n_donors_median=float(O_c.n_donors.median()),
                         median_offset_se_h=float(O_c.boot_se_h.median()), between_tissue_sd_h=float(np.sqrt(var_off)),
                         reliability=rel, mde_rho_approx=m_a, mde_rho_sim=m_s, mde_rho_attenuated=min(1.0, m_a / np.sqrt(rel)),
                         h3_null_sd_h=sd[name], h3_mde_80_h=MDE_K * sd[name], h3_cohort_weight=cw[name]))
    nt = len(O)
    prow.append(dict(cohort="pooled", n_tissues=nt, median_offset_se_h=float(O.boot_se_h.median()),
                     h3_null_sd_h=float(np.nanstd(pooled, ddof=1)), h3_mde_80_h=MDE_K * float(np.nanstd(pooled, ddof=1))))
    P = pd.DataFrame(prow)
    P.to_csv(out / "power.csv", index=False)
    say("power\n" + P.round(4).to_string(index=False))


def holm(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    adj = np.empty(len(p))
    adj[o] = np.minimum(1, np.maximum.accumulate((len(p) - np.arange(len(p))) * p[o]))
    return adj


def rho_p(off, sc, rng, n_perm=N_PERM_SCORE):
    ok = np.isfinite(off)
    off, sc = np.asarray(off, float)[ok], np.asarray(sc, float)[ok]
    r = stats.spearmanr(off, sc).statistic
    null = np.array([stats.spearmanr(off, rng.permutation(sc)).statistic for _ in range(n_perm)])
    return float(r), float((1 + np.sum(np.abs(null) >= abs(r) - 1e-12)) / (1 + n_perm))


def pseudo_task(i, data):
    """Pseudo-clock set i: per cohort offsets (fit, centred over targets) and per-tissue top-vs-bottom offset difference."""
    limit_threads(1)
    warnings.filterwarnings("ignore")
    res = {}
    for name, (pb, anchors, targets, labs, w) in data.items():
        off = offsets(pb, anchors, targets, waveform=w)
        ds = cs.Dataset(pb, anchors=anchors, waveform=w)
        obs = np.array([contrast_task(ds, t, labs[t], 0, 0, True)[0] for t in targets])
        res[name] = (off.to_numpy(float), obs)
    return i, res


def random_sets(c: Coh, n_sets, rng) -> list:
    """Column sets of the size of the G2M set, matched on deciles of mean log2(CPM + 1), clock and G2M genes excluded."""
    mean = np.asarray(c.Xl.mean(0)).ravel()
    cand = np.setdiff1d(c.pool, np.r_[c.clock, c.g2m])
    edges = np.quantile(mean[cand], np.linspace(0, 1, 11))
    dec = lambda x: np.clip(np.searchsorted(edges, x, side="right") - 1, 0, 9)
    bins = {b: cand[dec(mean[cand]) == b] for b in range(10)}
    want = dec(mean[c.g2m])
    out = []
    for _ in range(n_sets):
        pick = {b: iter(rng.choice(bins[b], int(np.sum(want == b)), replace=False)) for b in np.unique(want)}
        out.append([int(next(pick[b])) for b in want])
    return out


def stage_score(args, src, quick):
    out = src / "scramble" if args.scramble else src
    out.mkdir(parents=True, exist_ok=True)
    n_perm = args.n_perm or (20 if quick else 1000)
    n_sets = 5 if quick else 200
    n_rank = 200 if quick else N_PERM_SCORE
    pw = pd.read_csv(src / "power_tissues.csv")
    pp = pd.read_csv(src / "power.csv").set_index("cohort")
    boot = pd.read_csv(src / "tissue_offsets.csv")
    C, tis, gates = {}, [], []
    for name in COH:
        c, pb, targets, anchors, labs, ts = setup(name, args.subsample, args.seed, args.scramble)
        ds = cs.Dataset(pb, anchors=anchors, name=f"gtex_a232_{name}")
        gates.append(dict(cohort=name, **ds.gate(n_perm=50 if quick else 1000).to_series().to_dict()))
        off = offsets(pb, anchors, targets)
        obs, nul = contrasts(ds, targets, labs, n_perm, args.jobs, f"score{name}", preshuffle=False)
        wt = pw[pw.cohort == name].set_index("tissue").tissue_weight.reindex(targets).to_numpy()
        C[name] = dict(c=c, pb=pb, targets=targets, anchors=anchors, labs=labs, ts=ts, off=off, obs=obs, nul=nul, wt=wt)
        b = boot[boot.cohort == name].set_index("tissue")
        tis.append(pd.DataFrame(dict(cohort=name, tissue=targets, offset_h=off.to_numpy(), boot_se_h=b.boot_se_h.reindex(targets).to_numpy(),
                                     g2m=ts.g2m.to_numpy(), mki67_median=ts.mki67.to_numpy(), myc_median=ts.myc.to_numpy(),
                                     top_minus_bottom_h=obs, null_sd_h=np.nanstd(nul, axis=0, ddof=1))))
        say(f"{name}: observed offsets and contrasts done")
    pd.concat(tis).to_csv(out / "score_tissues.csv", index=False)
    pd.DataFrame(gates).to_csv(out / "score_gate.csv", index=False)
    cw = {n: float(pp.loc[n, "h3_cohort_weight"]) for n in COH}
    pooled = lambda per: sum(cw[n] * float(wmean(per[n], C[n]["wt"])) for n in COH)
    # companions: pseudo-clock sets through the whole pipeline
    rng = np.random.default_rng(SEED)
    pseudo = {}
    for name in COH:
        c = C[name]["c"]
        pseudo[name] = pseudo_waveforms(log_mean_of(c), n_sets, np.random.default_rng(SEED), cs.gtex())
    rows_in = {n: C[n]["c"].meta.tissue.isin(set(C[n]["targets"]) | set(C[n]["anchors"])).to_numpy() for n in COH}
    pos = {n: {g: j for j, g in reversed(list(enumerate(C[n]["c"].names)))} for n in COH}

    def gen():
        for i in range(n_sets):
            d = {}
            for n in COH:
                w = pseudo[n][i]
                pb = pseudobulk(C[n]["c"], [pos[n][g] for g in w.genes], list(w.genes), rows_in[n])
                d[n] = (pb, C[n]["anchors"], C[n]["targets"], C[n]["labs"], w)
            yield delayed(pseudo_task)(i, d)
    ps = []
    for k, (i, res) in enumerate(Parallel(n_jobs=args.jobs, return_as="generator")(gen()), 1):
        row = dict(set=i)
        for n in COH:
            off, ob = res[n]
            row[f"rho_{n}"] = stats.spearmanr(off, C[n]["ts"].g2m.to_numpy(), nan_policy="omit").statistic
        row["pooled_top_minus_bottom_h"] = pooled({n: res[n][1] for n in COH})
        ps.append(row)
        if k % 10 == 0:
            say(f"pseudo-clock set {k}/{n_sets}")
    PS = pd.DataFrame(ps).sort_values("set")
    PS.to_csv(out / "score_pseudo_clock.csv", index=False)
    # companions: random gene sets of G2M size
    rs = []
    for n in COH:
        c = C[n]["c"]
        sets = random_sets(c, n_sets, np.random.default_rng(SEED + 1))
        for j, cols in enumerate(sets):
            sc = zmean(c, cols)
            med = pd.Series(sc).groupby(c.meta.tissue.to_numpy()).median().reindex(C[n]["targets"]).to_numpy()
            rs.append(dict(cohort=n, set=j, rho=stats.spearmanr(C[n]["off"].to_numpy(float), med, nan_policy="omit").statistic))
    RS = pd.DataFrame(rs)
    RS.to_csv(out / "score_random_sets.csv", index=False)
    # tests
    rng = np.random.default_rng(A.seed_of("A232rank", args.seed))
    tests = []
    for h, n in (("H232.1", "fast"), ("H232.2", "vent")):
        t = C[n]
        r, p = rho_p(t["off"].to_numpy(float), t["ts"].g2m.to_numpy(), rng, n_rank)
        se = boot[boot.cohort == n].set_index("tissue").boot_se_h.reindex(t["targets"]).to_numpy()
        ok = np.isfinite(t["off"].to_numpy(float)) & np.isfinite(se) & (se > 0)
        slope = np.polyfit(t["ts"].g2m.to_numpy()[ok], t["off"].to_numpy(float)[ok], 1, w=1 / se[ok])[0]
        q_c = float(np.nanquantile(np.abs(PS[f"rho_{n}"]), 0.95))
        q_p = float(np.nanquantile(np.abs(RS[RS.cohort == n].rho), 0.95))
        tests.append(dict(hypothesis=h, cohort=n, n_tissues=len(t["targets"]), estimate=r, p=p, weighted_slope_h_per_unit=slope,
                          pseudo_clock_q95=q_c, random_set_q95=q_p, clock_pass=bool(abs(r) > q_c), prolif_pass=bool(abs(r) > q_p)))
    per = {n: C[n]["obs"] for n in COH}
    D = pooled(per)
    # same weights as the power stage (cohort weights from power.csv; tissue weights from power_tissues.csv)
    pn = sum(cw[n] * wmean(C[n]["nul"], C[n]["wt"]) for n in COH)
    p3 = float((1 + np.sum(np.abs(pn[np.isfinite(pn)]) >= abs(D) - 1e-12)) / (1 + np.isfinite(pn).sum()))
    q3 = float(np.nanquantile(np.abs(PS.pooled_top_minus_bottom_h), 0.95))
    tests.append(dict(hypothesis="H232.3", cohort="pooled", n_tissues=sum(len(C[n]["targets"]) for n in COH), estimate=D, p=p3,
                      null_sd=float(np.nanstd(pn, ddof=1)), mde_80=MDE_K * float(np.nanstd(pn, ddof=1)),
                      cohort_estimates=";".join(f"{n}:{float(wmean(per[n], C[n]['wt'])):.4f}" for n in COH),
                      pseudo_clock_q95=q3, clock_pass=bool(abs(D) > q3), prolif_pass=np.nan))
    T = pd.DataFrame(tests)
    T["p_holm"] = holm(T.p)
    T["specific"] = (T.p_holm < 0.05) & T.clock_pass & (T.prolif_pass.fillna(True).astype(bool))
    r1 = T.loc[T.hypothesis == "H232.1", "estimate"].iloc[0]
    T["same_sign_as_H232.1"] = np.sign(T.estimate) == np.sign(r1)
    T.to_csv(out / "score_tests.csv", index=False)
    say("tests\n" + T.round(4).to_string(index=False))


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--stage", choices=("power", "score"), required=True)
    ap.add_argument("--n-perm", type=int, default=None)
    ap.add_argument("--scramble", action="store_true", help="score stage code test on noise scores (needs --subsample)")
    ap.add_argument("--n-boot", type=int, default=None)
    args = ap.parse_args(argv)
    quick = cli.setup(args)["quick"]
    out = OUT / "subsample" if quick else OUT
    (out / "draws").mkdir(parents=True, exist_ok=True)
    warnings.filterwarnings("ignore", message="below the paper settings")
    {"power": stage_power, "score": stage_score}[args.stage](args, out, quick)


if __name__ == "__main__":
    main()
