"""A426: blind replication of the colder-louder lead (registered in paper/docs/registered_experiments.md, "A426").

Claim: within cohorts, tumours whose stromal clock is louder are immune-colder. Replication patients R1-R4 (`a426_data`: on-disk
patients A403 did not score, 48 CRC-atlas donors, Navin TNBC, GSE312983) with the A403 / A403c readout:
  host      per cancer, tumour-stroma host (anchors fibroblast / endothelial / mural, GTEx waveform, cohort x chemistry strata)
            fitted on the union of the cancer's A403 patients and its replication patients; a cancer without A403 patients
            enters only if its anchors-only gate (1,000 permutations, p < 0.05 on pairing and acrophase) passes with >= 8
            replication patients
  loudness  log host radius centred within stratum, made phase-free by the A403c C1 residual (two harmonics of the host phase,
            per cancer over the union patients, stratum-centred)
  readout   A403 immune-hot score (CD8A, GZMB, IFNG, CXCL9, CXCL10, PRF1; log2 CPM + 1 of the patient's all-compartment sum;
            mean within-stratum z)
  covariates  number of stromal anchors and log median stromal depth, residualised within stratum (A403 note 1)
  slope     pooled (patient weights) within-stratum standardised slope on the REPLICATION patients only; strata with < 3
            replication patients dropped
  H426.1    slope < 0, one-sided permutation p < 0.05 (readout shuffled within stratum, 10,000, seed [426, 1]) and below the q05 of
            the same slope on 200 matched pseudo-hosts per cancer (phase-free likewise; seed [426, 0], the A401 draw rule)

  --stage power  BLIND. Counts per source / cancer / stratum, gates of the new cancers (anchors only), null SD of the slope and
                 power at -0.10 and -0.126 (normal approximation), companion q05, all from READOUTS SHUFFLED WITHIN STRATUM
                 BEFORE ANY SLOPE; the real host is never fitted on the replication patients (companion pseudo-hosts only).
                 No loudness, host radius or phase, immune-hot value or slope of the replication patients is printed or stored.
                                                        -> output/a426/power/{patients,counts,strata,gates,cancers,companions,power}.csv
  --stage score  NOT TO BE RUN before the power note is written (one run). H426.1 and the descriptive items: A403 raw-loudness
                 slope, per cancer with I^2, R1 held-out patients of A403 cohorts vs independent cohorts, scRNA vs snRNA (R4),
                 T/NK and myeloid shares, fixed-effect combination with A403c.        -> output/a426/score/
  --smoke        code check on a few small sources, 6 pseudo-hosts, 200 permutations; the readout is shuffled within stratum in both
                 stages and nothing observed is written.

    python a426.py --stage power [--jobs 4] [--genome-jobs 2] [--smoke]       (long: paper/common/run.sh)
"""
from __future__ import annotations

import argparse
import sys
import time
import warnings
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402
from scipy.stats import norm  # noqa: E402

import a401 as U  # noqa: E402
import a402_a405 as W  # noqa: E402
import a403_a405_score as S  # noqa: E402
import a403b as B3  # noqa: E402
import a403c as C3  # noqa: E402
import a426_data as D  # noqa: E402
import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402

OUT = U.OUT / "a426"
CACHE = D.DATA / "a426_cache"
BASE, META, ANCH, MALT = U.BASE, U.META, U.ANCH, U.MALT
IMMUNE_HOT = W.IMMUNE_HOT
COVS = ["n_anchor", "log_depth"]
SEED_SET, SEED_PERM, GATE_SEED = [426, 0], [426, 1], 426
N_PERM, N_SET, GATE_PERM = 10_000, 200, 1_000
MIN_STRATUM, MIN_NEW, B_TRUE = 3, 8, (-0.10, -0.126)
SMOKE_SOURCES = ("SRP414286", "bassez2021", "navin_tnbc_cxg", "GSE125449")
log = W.log


# ---- A403 side: anchors of the 12 P401 cancers ----------------------------------------------------------------------------

def a403_pool(args, out):
    """(anchor rows of the A403 pool, A403 donors), the A401 pool of the 12 P401 cancers on 14 genes (as A403c); cached."""
    f = CACHE / ("a403_anchor_pool_smoke.pkl" if args.smoke else "a403_anchor_pool.pkl")
    if f.exists():
        return pd.read_pickle(f)
    frames, donors, _ = U.pool_sources(args, out)
    from joblib.externals.loky import get_reusable_executor
    get_reusable_executor().shutdown(wait=True)
    p401 = pd.read_csv(U.OUT / "power" / "p401.csv")
    cancers = [c for c in p401[(p401.readout.astype(str) == "14") & p401.p401_pass].cancer if c in set(donors.cancer) and c != "unspecified"]
    donors = donors[donors.cancer.isin(cancers)].reset_index(drop=True)
    ids = set(donors.donor)
    anch = {s: g[g.donor.isin(ids) & g.compartment.isin(ANCH)].reset_index(drop=True) for s, g in frames.items()}
    res = ({s: g for s, g in anch.items() if len(g)}, donors[["donor", "source", "cancer"]])
    frames.clear()
    W.trim()
    CACHE.mkdir(parents=True, exist_ok=True)
    pd.to_pickle(res, f)
    return res


# ---- replication side -------------------------------------------------------------------------------------------------------

def replication(args):
    """Eligible replication patients: unit rule (>= 2 stromal classes, >= 50 immune cells), no duplicate id, source measures the 14
    clock genes. Returns (all-compartment frames of eligible patients, patients table, sources table)."""
    frames, pat, srcs = D.load_replication(only=SMOKE_SOURCES if args.smoke else None, verbose=log)
    pat = D.overlap_with_a403(pat)
    pat["eligible"] = pat.unit & ~pat.id_overlap
    ok = set(pat.donor[pat.eligible])
    frames = {s: f[f.donor.isin(ok)].reset_index(drop=True) for s, f in frames.items()}
    return {s: f for s, f in frames.items() if len(f)}, pat, srcs


def readout_inputs(frames, pat):
    """Per eligible patient: counts of the six immune-hot genes and the depth of the all-compartment sum, and the cell counts behind
    the T/NK and myeloid shares. Counts, not scores; used only by `immune_hot` / `shuffled_readout`."""
    rows = []
    for s, f in frames.items():
        have = [g for g in IMMUNE_HOT if g in f.columns]
        t = f.groupby("donor")[have + ["depth", "n_cells"]].sum()
        t["tnk_cells"] = f[f.compartment == "tumour:T_NK"].groupby("donor").n_cells.sum().reindex(t.index).fillna(0)
        t["myeloid_cells"] = f[f.compartment == "tumour:myeloid"].groupby("donor").n_cells.sum().reindex(t.index).fillna(0)
        t["stratum"] = f.drop_duplicates("donor").set_index("donor").stratum.reindex(t.index)
        rows.append(t)
    return pd.concat(rows)


def immune_hot(ro):
    """A403 immune-hot score: mean over the available genes (>= 3) of the within-stratum z of log2 CPM + 1."""
    genes = [g for g in IMMUNE_HOT if g in ro.columns and ro[g].notna().any()]
    if len(genes) < 3:
        return pd.Series(np.nan, index=ro.index)
    lc = np.log2(ro[genes].div(ro.depth, axis=0) * 1e6 + 1)
    g = lc.groupby(ro.stratum)
    z = (lc - g.transform("mean")) / g.transform(lambda x: x.std(ddof=0)).replace(0, np.nan)
    return z.mean(axis=1)


def shuffled_readout(ro, rng, col="immune_hot"):
    """The readout shuffled within stratum, in one step: the unshuffled values never leave this function."""
    y = immune_hot(ro) if col == "immune_hot" else ro[col]
    out = y.copy()
    for _, i in y.groupby(ro.stratum).indices.items():
        out.iloc[i] = y.iloc[i].to_numpy()[rng.permutation(len(i))]
    return out


# ---- tables, gates, companions ----------------------------------------------------------------------------------------------

def cancer_tables(frames, donors, cancers):
    """As `a402_a405.cancer_hosts` without fitting the host: per cancer the stroma-only table of the classes with >= 5 patients."""
    out = {}
    for c in cancers:
        ids = set(donors[donors.cancer == c].donor)
        parts = [f[f.donor.isin(ids) & f.compartment.isin(ANCH)] for f in frames.values()]
        parts = [p for p in parts if len(p)]
        Gc = pd.concat(parts, join="inner", ignore_index=True)
        classes = [a for a in ANCH if Gc[Gc.compartment == a].donor.nunique() >= 5]
        if len(classes) < 2:
            continue
        Gc = Gc[Gc.compartment.isin(classes)].reset_index(drop=True)
        tab = Gc[META + list(BASE.genes)]
        out[c] = dict(classes=classes, Gc=Gc, tab=tab, stratum=tab.drop_duplicates("donor").set_index("donor").stratum)
    return out


def gate_new(c, frames, donors, n_perm):
    """Anchors-only gate of a cancer on its replication patients (14 genes, mixed strata)."""
    ids = set(donors[donors.cancer == c].donor)
    parts = [f[f.donor.isin(ids) & f.compartment.isin(ANCH)] for f in frames.values()]
    parts = [p for p in parts if len(p)]
    Gc = pd.concat(parts, join="inner", ignore_index=True)
    classes = [a for a in ANCH if Gc[Gc.compartment == a].donor.nunique() >= 5]
    base = dict(cancer=c, patients=len(ids), classes="+".join(k.split(":")[1] for k in classes))
    if len(classes) < 2:
        return dict(base, status="not estimable: < 2 stromal classes with >= 5 patients", enters=False)
    tab = Gc[Gc.compartment.isin(classes)][META + list(BASE.genes)]
    limit_threads(1)
    g = cs.Dataset(tab, anchors=classes, waveform=BASE, mixed_strata=True).gate(n_perm=n_perm, seed=GATE_SEED)
    return dict(base, status="ok", p_pairing=g.p_pairing, p_acrophase=g.p_acrophase, clock_share=g.clock_share,
                share_null_mean=g.share_null_mean, gate_passes=bool(g.passes), enters=bool(g.passes and len(ids) >= MIN_NEW))


def covariates(tabs):
    return S.covariates({c: dict(tab=t["tab"]) for c, t in tabs.items()})


def companion_loudness(tabs, args, n_set):
    """200 matched pseudo-hosts per cancer on the union stroma (seed [426, 0], the A401 draw rule). Returns [(raw, phase-free)] per set,
    each a Series donor -> loudness (log radius centred within stratum; C1 residual)."""
    from clocksync.paired import pseudo_waveforms
    per_set = [[] for _ in range(n_set)]
    for c, t in tabs.items():
        Gc = t["Gc"]
        genes = [g for g in Gc.columns if g not in META]
        frac = Gc[genes].div(Gc.depth, axis=0)
        lm = np.log(frac.mean() + 1e-12)[(Gc[genes].mean() >= 1) | frac.columns.isin(list(BASE.genes))]
        ws = pseudo_waveforms(lm, n_set, np.random.default_rng(SEED_SET), BASE, "matched")
        res = Parallel(n_jobs=args.jobs)(delayed(C3.pseudo_states)(Gc[META + sorted(set(w.genes) | set(BASE.genes))], t["classes"], w) for w in ws)
        for j, r in enumerate(res):
            if len(r):
                per_set[j].append(r.assign(cancer=c))
        log(f"companion pseudo-hosts {c}")
    stratum = pd.concat([t["stratum"] for t in tabs.values()])
    out = []
    for parts in per_set:
        if not parts:
            out.append((pd.Series(dtype=float), pd.Series(dtype=float)))
            continue
        d = pd.concat(parts)
        d = d[np.isfinite(d.radius.astype(float)) & (d.radius.astype(float) > 0)]
        st = stratum.reindex(d.index)
        lr = pd.Series(W.centred(np.log(d.radius.astype(float).to_numpy()), st), index=d.index)
        pf, _ = C3.phase_free(lr, np.angle(d.z.to_numpy().astype(complex)), st.to_numpy(), d.cancer.to_numpy())
        out.append((lr, pf.dropna()))
    return out


# ---- the slope --------------------------------------------------------------------------------------------------------------

def units(df, y_col, loud_col, covs=COVS, min_n=MIN_STRATUM):
    """Per stratum (>= `min_n` patients): stratum-centred readout and loudness, residualised on the covariates when n > #covariates + 2,
    z-scored (A403 `units`; strata play the role of cohorts)."""
    out = []
    for s, g in df.dropna(subset=[y_col, loud_col] + covs).groupby("stratum"):
        if len(g) < min_n:
            continue
        y = W.centred(g[y_col].to_numpy(float), g.stratum)
        l_ = W.centred(g[loud_col].to_numpy(float), g.stratum)
        X = np.column_stack([W.centred(g[c].to_numpy(float), g.stratum) for c in covs]) if covs else np.empty((len(g), 0))
        if len(g) > X.shape[1] + 2:
            y, l_ = S.residual(y, X), S.residual(l_, X)
        if y.std() < 1e-12 or l_.std() < 1e-12:
            continue
        out.append(dict(stratum=s, n=len(g), y=y / y.std(), l=l_ / l_.std()))
    return out


def null_slopes(un, rng, n_perm):
    perms = [np.argsort(rng.random((n_perm, u["n"])), axis=1) for u in un]
    return W.pooled_slope(un, perms)


def analytic_null_sd(un):
    """SD of the pooled slope under within-unit permutation of standardised variables: sqrt(sum n_u^2 / (n_u - 1)) / sum n_u."""
    n = np.array([u["n"] for u in un], float)
    return float(np.sqrt(np.sum(n ** 2 / (n - 1))) / n.sum())


# ---- power stage --------------------------------------------------------------------------------------------------------------

def setup(args, out):
    """Union tables of the cancers that enter, the replication table of their patients and the gates of the new cancers."""
    anch403, d403 = a403_pool(args, out)
    frames, pat, srcs = replication(args)
    srcs.to_csv(out / "sources.csv", index=False)
    rep = pat[pat.eligible]
    a403_cancers = set(d403.cancer)
    gates = []
    cancers = sorted(set(rep.cancer))
    entering = []
    for c in cancers:
        n = int((rep.cancer == c).sum())
        if c in a403_cancers:
            entering.append(c)
            gates.append(dict(cancer=c, patients=n, status="A403 cancer: host on the union, no gate", enters=True))
        elif n >= MIN_NEW:
            gates.append(gate_new(c, frames, rep, 100 if args.smoke else GATE_PERM))
            if gates[-1]["enters"]:
                entering.append(c)
            log(f"gate {c}: {gates[-1].get('status')} passes {gates[-1].get('gate_passes')}")
        else:
            gates.append(dict(cancer=c, patients=n, status=f"< {MIN_NEW} replication patients: not gated, does not enter", enters=False))
    gates = pd.DataFrame(gates)
    gates.to_csv(out / "gates.csv", index=False)
    # union tables: A403 patients are replaced by their replication rows when both exist (Wu / Pal)
    rep_ids = set(rep.donor)
    un_frames = {f"A403:{s}": f[~f.donor.isin(rep_ids)] for s, f in anch403.items()}
    un_frames.update({f"R:{s}": f for s, f in frames.items()})
    un_donors = pd.concat([d403[~d403.donor.isin(rep_ids)][["donor", "source", "cancer"]],
                           rep[rep.cancer.isin(entering)][["donor", "source", "cancer"]]], ignore_index=True)
    un_donors = un_donors[un_donors.cancer.isin(entering)]
    tabs = cancer_tables(un_frames, un_donors, entering)
    return dict(frames=frames, pat=pat, rep=rep, gates=gates, tabs=tabs, entering=entering, d403=d403, ro=readout_inputs(frames, rep))


def rep_frame(st):
    """Replication patients with a stroma table (>= 1 retained class): cancer, stratum, covariates. No clock quantity."""
    cov = covariates(st["tabs"])
    strat = pd.concat([t["stratum"] for t in st["tabs"].values()])
    rep = st["rep"][st["rep"].donor.isin(strat.index) & st["rep"].donor.isin(cov.donor)]
    df = rep[["donor", "cancer", "source", "tier"]].assign(stratum=lambda d: strat.reindex(d.donor).to_numpy()).merge(cov, on="donor")
    n = df.groupby("stratum").donor.transform("size")
    return df.assign(stratum_n=n, in_stratum_ge3=n >= MIN_STRATUM)


def stage_power(args):
    out = OUT / ("power_smoke" if args.smoke else "power")
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    n_set, n_perm = (6, 200) if args.smoke else (N_SET, N_PERM)
    st = setup(args, out)
    df = rep_frame(st)
    st["pat"].to_csv(out / "patients.csv", index=False)
    # counts per source, cancer and stratum (patients with a stroma table; strata of < 3 are dropped from the slope)
    cnt = df.groupby(["source", "cancer", "stratum"]).agg(patients=("donor", "size"), stratum_patients=("stratum_n", "first"),
                                                         used=("in_stratum_ge3", "first")).reset_index()
    cnt.to_csv(out / "counts.csv", index=False)
    sa = df.groupby("stratum").agg(patients=("donor", "size"), cancers=("cancer", lambda x: "+".join(sorted(set(x)))), source=("source", "first")).reset_index()
    sa["kept"] = sa.patients >= MIN_STRATUM
    sa.to_csv(out / "strata.csv", index=False)
    ca = df[df.in_stratum_ge3].groupby("cancer").agg(patients=("donor", "size"), strata=("stratum", "nunique"), sources=("source", "nunique")).reset_index()
    ca["a403_patients"] = ca.cancer.map(st["d403"].groupby("cancer").size()).fillna(0).astype(int)
    ca["classes"] = ca.cancer.map(lambda c: "+".join(k.split(":")[1] for k in st["tabs"][c]["classes"]) if c in st["tabs"] else "")
    ca.to_csv(out / "cancers.csv", index=False)
    log("counts per cancer\n" + ca.to_string(index=False))
    # companion loudness (pseudo-hosts only) and slopes of SHUFFLED readouts
    comp = companion_loudness(st["tabs"], args, n_set)
    W.trim()
    rows = []
    for j, (lr, pf) in enumerate(comp):
        y = shuffled_readout(st["ro"], np.random.default_rng([426, 3, j])).rename("y")
        for kind, loud in (("raw", lr), ("phase_free", pf)):
            d = df.merge(y, left_on="donor", right_index=True).merge(loud.rename("loud"), left_on="donor", right_index=True)
            un = units(d, "y", "loud")
            if not un:
                continue
            nul = null_slopes(un, np.random.default_rng([426, 4, j, int(kind == "raw")]), 1000 if not args.smoke else 200)
            rows.append(dict(set=j, loudness=kind, slope_shuffled_readout=W.pooled_slope(un), null_sd=float(nul.std()), analytic_null_sd=analytic_null_sd(un),
                             patients=int(sum(u["n"] for u in un)), strata=len(un)))
    C = pd.DataFrame(rows)
    C.to_csv(out / "companions.csv", index=False)
    summ = []
    for kind, g in C.groupby("loudness"):
        sd = float(g.null_sd.median())
        q05, q95, med = (float(g.slope_shuffled_readout.quantile(p)) for p in (.05, .95, .5))
        for b in B_TRUE:
            thr_perm = norm.ppf(0.05) * sd
            thr = min(thr_perm, q05)
            summ.append(dict(loudness=kind, true_slope=b, patients=int(g.patients.median()), strata=int(g.strata.median()), sets=len(g), null_sd_median=sd,
                             null_sd_min=float(g.null_sd.min()), null_sd_max=float(g.null_sd.max()), analytic_null_sd=float(g.analytic_null_sd.median()),
                             companion_q05=q05, companion_median=med, companion_q95=q95, companion_sd=float(g.slope_shuffled_readout.std()),
                             power_permutation_only=float(norm.cdf((thr_perm - b) / sd)), threshold_both=thr,
                             power_permutation_and_companion=float(norm.cdf((thr - b) / sd)),
                             mde_80=float(-(norm.ppf(0.8) - norm.ppf(0.05)) * sd)))
    P = pd.DataFrame(summ)
    P.to_csv(out / "power.csv", index=False)
    pd.set_option("display.width", 250, "display.max_columns", 40)
    log("power\n" + P.round(4).to_string(index=False))
    log(f"power stage done {(time.time() - t0) / 60:.1f} min")


# ---- score stage (NOT TO BE RUN before the power note) -------------------------------------------------------------------------

def stage_score(args):
    out = OUT / ("score_smoke" if args.smoke else "score")
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    n_set, n_perm = (6, 200) if args.smoke else (N_SET, N_PERM)
    st = setup(args, out)
    # real hosts on the union (this is the first fit of the host on the replication patients)
    anch403, d403 = a403_pool(args, out)
    rep = st["rep"]
    rep_ids = set(rep.donor)
    frames = {f"A403:{s}": f[~f.donor.isin(rep_ids)] for s, f in anch403.items()}
    frames.update({f"R:{s}": f for s, f in st["frames"].items()})
    donors = pd.concat([d403[~d403.donor.isin(rep_ids)], rep[["donor", "source", "cancer"]]], ignore_index=True)
    hosts = W.cancer_hosts(frames, donors, st["entering"])
    allloud = W.loudness_frame(hosts)
    allloud["loud_pf"], _ = C3.phase_free(allloud.loud, np.angle(allloud.z.to_numpy()), allloud.stratum.to_numpy(), allloud.cancer.to_numpy())
    df = allloud[allloud.donor.isin(rep_ids)].drop(columns=["z", "cohort"]).merge(S.covariates(hosts), on="donor", how="left")
    ro = st["ro"]
    df = df.merge(immune_hot(ro).rename("immune_hot"), left_on="donor", right_index=True)
    for c in ("tnk_cells", "myeloid_cells", "n_cells"):
        df[c] = df.donor.map(ro[c])
    for nm, c in (("tnk_share", "tnk_cells"), ("myeloid_share", "myeloid_cells")):
        df[nm] = (df[c] / df.n_cells).where(df.n_cells >= 100)
    df = df.merge(rep[["donor", "source", "tier"]], on="donor")
    df["cohort"] = df.stratum
    if args.smoke:  # code check: readout shuffled within stratum, nothing observed is written
        df["immune_hot"] = shuffled_readout(ro, np.random.default_rng([426, 9])).reindex(df.donor).to_numpy()
    df = df[df.groupby("stratum").donor.transform("size") >= MIN_STRATUM].reset_index(drop=True)
    comp = companion_loudness(hosts, args, n_set)
    # H426.1 (primary): phase-free loudness
    un = units(df, "immune_hot", "loud_pf")
    obs = W.pooled_slope(un)
    nul = null_slopes(un, np.random.default_rng(SEED_PERM), n_perm)
    p_one = float((np.sum(nul <= obs + 1e-12) + 1) / (len(nul) + 1))
    sl = []
    for lr, pf in comp:
        d2 = df.drop(columns="loud_pf").merge(pf.rename("loud_pf"), left_on="donor", right_index=True)
        u2 = units(d2, "immune_hot", "loud_pf") if len(d2) else []
        sl.append(W.pooled_slope(u2) if u2 else np.nan)
    sl = np.array(sl)[np.isfinite(sl)]
    q05 = float(np.quantile(sl, .05))
    pd.DataFrame([dict(hypothesis="H426.1", loudness="phase_free", patients=int(sum(u["n"] for u in un)), strata=len(un), slope=obs, null_sd=float(nul.std()),
                       p_one_sided_negative=p_one, companion_n=len(sl), companion_median=float(np.median(sl)), companion_q05=q05, companion_q95=float(np.quantile(sl, .95)),
                       passes=bool(p_one < 0.05 and obs < q05))]).to_csv(out / "h426.csv", index=False)
    pd.DataFrame({"set": range(len(sl)), "slope": sl}).to_csv(out / "companions.csv", index=False)
    # descriptive
    desc = []
    for lab, col, sub in (("raw loudness (A403 quantity)", "loud", df), ("phase-free, R1 held-out patients of A403 cohorts", "loud_pf", df[df.tier == "R1"]),
                          ("phase-free, R2-R4 independent cohorts", "loud_pf", df[df.tier != "R1"]),
                          ("phase-free, R4 (snRNA)", "loud_pf", df[df.tier == "R4"]), ("phase-free, scRNA (R1-R3)", "loud_pf", df[df.tier != "R4"])):
        u2 = units(sub, "immune_hot", col)
        if u2:
            o2 = W.pooled_slope(u2)
            n2 = null_slopes(u2, np.random.default_rng(SEED_PERM), n_perm)
            desc.append(dict(item=lab, patients=int(sum(u["n"] for u in u2)), strata=len(u2), slope=o2, null_sd=float(n2.std()),
                             p_one_sided_negative=float((np.sum(n2 <= o2 + 1e-12) + 1) / (len(n2) + 1))))
    pd.DataFrame(desc).to_csv(out / "descriptive.csv", index=False)
    rows = []
    for rd in ("immune_hot", "tnk_share", "myeloid_share"):
        rows += B3.per_cancer(df, rd, COVS, {}, f"a426{rd}", loud_col="loud_pf")
    per = pd.DataFrame(rows)
    per.to_csv(out / "per_cancer.csv", index=False)
    c2 = pd.read_csv(U.OUT / "a403c" / "c2.csv")
    a = c2[(c2.cancer == "pooled") & (c2.loudness == "phase_free") & (c2.readout == "immune_hot")].iloc[0]
    r = per[(per.cancer == "pooled") & (per.readout == "immune_hot")].iloc[0]
    w = np.array([1 / a.boot_sd ** 2, 1 / r.boot_sd ** 2])
    s_ = np.array([a.slope, r.slope])
    pd.DataFrame([dict(item="fixed-effect combination of A403c (post hoc discovery) and A426 (blind replication)", slope_a403c=a.slope, slope_a426=r.slope,
                       combined=float((w * s_).sum() / w.sum()), combined_se=float(1 / np.sqrt(w.sum())))]).to_csv(out / "combined.csv", index=False)
    log(f"score stage done {(time.time() - t0) / 60:.1f} min")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stage", choices=("power", "score"), required=True)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--genome-jobs", type=int, default=2)
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()
    warnings.filterwarnings("ignore")
    limit_threads(1)
    stage_power(args) if args.stage == "power" else stage_score(args)


if __name__ == "__main__":
    main()
