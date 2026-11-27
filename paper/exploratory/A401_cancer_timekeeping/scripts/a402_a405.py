"""A402-A405 power stages (registered in paper/docs/registered_experiments.md, "A402-A405"). POWER ONLY.

NO TARGET STATISTIC is computed: no observed A246 omnibus (A402), no observed slope of a readout on loudness (A403), no
observed immune-compartment coupling (A404), no observed group difference (A405). Computed: loudness and readouts
themselves, nulls from shuffled pairings, companion distributions, counts, the blind simulator
(`clocksync.simulate.like(..., loading=z)`).

Data: the A401 pooled tumour table (`a401.pool_sources`: A388 cohorts, GEO builds, Census, BIOKEY; the same deduplication),
cancers passing P401 on 14 genes (`output/power/p401.csv`), host = the tumour's own stroma. The per-source slim tables of
all tumour compartments (clock, A246 program, immune-hot and G2M genes) are cached in `output/cache_a402_a405/`.

  A402  A246 omnibus null with state = tumour-stroma host (posterior mean of the cancer's pooled host), on the patients of
        cohorts new relative to A246 (A246 used the A227 cohorts): programs by `compute_a246` (A146 genes and scores), r of
        A229 per cohort pooled with patient weights, states shuffled within cohort (one ordering per iteration shared by
        all programs, 2,000), omnibus T = mean over scored programs of (r - null mean) / null SD; MDE = 2.49 SD; decision
        MDE(T) <= 0.43.                                                        -> a402_power.csv
  A403  per-patient loudness (log radius of the host, centred within cohort x chemistry) and readouts (immune-hot score,
        T/NK share, myeloid share, G2M proliferation index); slope (z-scored readout on z-scored loudness, within cohort x
        chemistry, pooled with patient weights); null = readout permuted within cohort (10,000); MDE = 2.49 SD; slope of
        shuffled readouts on the loudness of 200 matched pseudo-hosts (companion loudness) -> a403_power.csv, a403_counts.csv
  A404  patients per immune compartment on (a) the A387 normal-stroma table, (b) the A401 tumour-stroma table; blind-
        simulator power of the coupling (|z| 0.5, c 0.1 / 0.2 / 0.3, null threshold from c 0.02); companion q95 of the
        immune couplings (200 matched, 200 narrow responsive) on the real data           -> a404_power.csv
  A405  BIOKEY group sizes (pre; both time points; expansion x cohort) and the MDE of a pre-treatment E vs NE difference of
        loudness from 2,000 label shuffles (SD of the difference x 2.49)                 -> a405_power.csv

    python a402_a405.py --stage power [--jobs 4] [--genome-jobs 2] [--smoke]       (long: paper/common/run.sh)
"""
from __future__ import annotations

import argparse
import itertools
import sys
import time
import warnings
from pathlib import Path

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
EXP = PAPER / "exploratory"
sys.path.insert(0, str(HERE.parent))
for p in ("A246_tme_clock_state/scripts", "A229_endothelial_gate/scripts", "A146_tme_time/scripts", "A387_tumour_full_model/scripts"):
    sys.path.insert(0, str(EXP / p))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import a401 as U  # noqa: E402
import a388 as B  # noqa: E402
import compute_a246 as P246  # noqa: E402
import compute_a227 as M  # noqa: E402
import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.genes import canonical  # noqa: E402

OUT = U.OUT / "power_a402_a405"
import os  # noqa: E402
CACHE = U.OUT / os.environ.get("A402_CACHE", "cache_a402_a405")  # a402b: own cache with extra genes (env reaches the loky workers)
MALT, META, BASE, CORE, ANCH = U.MALT, U.META, U.BASE, U.CORE, U.ANCH
IMMUNE_HOT = ["CD8A", "GZMB", "IFNG", "CXCL9", "CXCL10", "PRF1"]
IMMUNE = ["tumour:myeloid", "tumour:T_NK", "tumour:B_plasma"]
SEED = 402
MDE_K = 2.49


def trim():
    """Return freed memory to the OS (the genome-wide companion builds of A404 need it)."""
    import ctypes
    import gc
    gc.collect()
    try:
        ctypes.CDLL("libc.so.6").malloc_trim(0)
    except OSError:
        pass


def log(*a):
    print(time.strftime("%F %T"), *a, flush=True)


def g2m_genes():
    import importlib.util
    spec = importlib.util.spec_from_file_location("fig5_shared", PAPER / "results/figure_5/scripts/cancer/_shared.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return [canonical(g) for g in m.hallmark_g2m()]


G2M = g2m_genes()
NEEDED = set(canonical(g) for g in BASE.genes) | {g for v in P246.GENES.values() for g in v} | set(IMMUNE_HOT) | set(G2M)
NEEDED |= {canonical(g) for g in os.environ.get("A402_EXTRA_GENES", "").split(",") if g}  # a402b: T/NK subtype markers


# ---- slim tables of all tumour compartments ------------------------------------------------------------------------------

def slim_a388(study, patients, mal):
    """Donor x tumour-compartment pseudobulks of the needed genes of one A388 cohort rebuilt from counts (labels of
    `cancer_sc`; copy-number calls for the malignant split where the cohort has them); compartments with >= 50 cells."""
    import a149_lib as L
    limit_threads(1)
    warnings.filterwarnings("ignore")
    import compute_a226 as K226
    K226.install()
    f = CACHE / f"slim_{study}.parquet"
    if f.exists():
        return study, pd.read_parquet(f), ""
    try:
        us, _ = L.ST.units(study)
        us = [u for u in L.ST.majority_chemistry(us) if u[1] in set(patients)]
        calls = None
        if mal == L.MAL:
            calls = L.new_calls(study) if study in L.NEW else L.ST.read_calls(study)
        sums, ncell, depth, strat, cols = {}, {}, {}, {}, None
        for sample, patient, kind, load, chem, _ in us:
            if kind != "tumour":
                continue
            X, _g, umi, comp, barcodes, sym = load()
            lab = L.B._labels(comp, kind, sample, barcodes, calls)
            if cols is None:
                canon = np.array([canonical(str(g)) for g in sym])
                cols = np.flatnonzero(np.isin(canon, list(NEEDED)))
                names = canon[cols]
            for c in set(lab):
                if c.endswith(":unassigned"):
                    continue
                i = np.flatnonzero(lab == c)
                k = (patient, c)
                sums[k] = sums.get(k, 0) + np.asarray(X[i][:, cols].sum(0)).ravel()
                ncell[k] = ncell.get(k, 0) + len(i)
                depth[k] = depth.get(k, 0) + float(umi[i].sum())
                strat[k] = chem
        keys = sorted(sums)
        G = pd.DataFrame(np.vstack([sums[k] for k in keys]).astype(np.float32), columns=names)
        G = G.T.groupby(level=0).sum().T
        meta = pd.DataFrame({"donor": [f"{study}|{k[0]}" for k in keys], "compartment": [k[1] for k in keys],
                             "stratum": [f"{study}|{strat[k]}" for k in keys], "n_cells": [ncell[k] for k in keys],
                             "depth": [depth[k] for k in keys]})
        G = pd.concat([meta, G.reset_index(drop=True)], axis=1)
        G["compartment"] = G.compartment.where(G.compartment != "tumour:epithelial", MALT)
        G = G[G.n_cells >= U.MIN_CELLS].reset_index(drop=True)
        CACHE.mkdir(parents=True, exist_ok=True)
        G.to_parquet(f)
        return study, G, ""
    except Exception as e:  # noqa: BLE001
        return study, None, f"{type(e).__name__}: {e}"[:300]


def slim_parquet(s, donors):
    f = CACHE / f"slim_{s}.parquet"
    if f.exists():
        return pd.read_parquet(f)
    if s in U.GEO_CANCER:
        G, _ = U.build_geo(s, True, NEEDED)
    elif s == "bassez2021":
        G, _ = U.build_bassez(True, NEEDED)
    else:
        G, _ = U.build_census(s, True, NEEDED)
    G = G[G.donor.isin(donors)].reset_index(drop=True)
    CACHE.mkdir(parents=True, exist_ok=True)
    G.to_parquet(f)
    return G


def load_slim(sources, donors, f388, args):
    """Slim all-compartment tables of every source with donors of the P401 cancers (patients of `donors`)."""
    out, errs = {}, []
    a388 = [s for s in f388 if s in set(donors.source) and s not in ("GSE176078", "GSE161529")]
    others = [s for s in sources if s in set(donors.source) and s not in f388]
    inv = pd.read_csv(B.FEAS).set_index("study")
    mals = {}
    for s in a388:
        r = B.cohort_table(s, None, 1)
        mals[s] = r[2] if not isinstance(r, str) else None
    res = Parallel(n_jobs=args.genome_jobs)(
        delayed(slim_a388)(s, [d.split("|", 1)[1] for d in donors[donors.source == s].donor], mals[s]) for s in a388 if mals[s])
    for s, g, e in res:
        if g is None:
            errs.append(dict(source=s, error=e))
        else:
            out[s] = g[g.donor.isin(set(donors.donor))]
    for s in others:
        try:
            out[s] = slim_parquet(s, set(donors[donors.source == s].donor))
        except Exception as e:  # noqa: BLE001
            errs.append(dict(source=s, error=f"{type(e).__name__}: {e}"[:300]))
    return out, pd.DataFrame(errs, columns=["source", "error"])


# ---- hosts, loudness --------------------------------------------------------------------------------------------------

def cancer_hosts(frames, donors, cancers, ws_n=0):
    """Per cancer: stroma-only table, classes, host (14-gene readout) of the pooled tumour stroma."""
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
        h = cs.Dataset(tab, anchors=classes, waveform=BASE, mixed_strata=True).host()
        strat = tab.drop_duplicates("donor").set_index("donor").stratum
        out[c] = dict(classes=classes, Gc=Gc, tab=tab, host=h, stratum=strat)
    return out


def centred(v, stratum):
    v = pd.Series(v)
    return (v - v.groupby(stratum.to_numpy()).transform("mean")).to_numpy()


def loudness_frame(hosts):
    """Donor table: cancer, cohort, stratum, loudness (log radius centred within cohort x chemistry), host state z."""
    rows = []
    for c, h in hosts.items():
        hh = h["host"]
        st = h["stratum"].reindex(hh.index)
        L_ = centred(np.log(hh.radius.to_numpy()), st)
        rows.append(pd.DataFrame({"donor": hh.index, "cancer": c, "stratum": st.to_numpy(), "loud": L_,
                                  "z": (hh.x + 1j * hh.y).to_numpy(), "cohort": [d.split("|")[0] for d in hh.index]}))
    return pd.concat(rows, ignore_index=True)


# ---- A402 --------------------------------------------------------------------------------------------------------------

def cohort_dict(G, ids):
    """A246-format cohort data (`compute_a246.cohort_data`) of one source from its slim table, patients `ids`."""
    G = G[G.donor.isin(ids)]
    pts = sorted(set(G.donor))
    ix = {p: i for i, p in enumerate(pts)}
    strat = G.drop_duplicates("donor").set_index("donor").stratum.reindex(pts).to_numpy()
    cells = G.groupby("donor").n_cells.sum().reindex(pts).to_numpy(float)
    pg = sorted({g for v in P246.GENES.values() for g in v})
    comps, absent = {}, {}
    for key, name in (("mal", MALT), ("T_NK", "tumour:T_NK"), ("myeloid", "tumour:myeloid")):
        g = G[G.compartment == name].drop_duplicates("donor").set_index("donor")
        if len(g) < P246.MIN_N:
            absent[key] = f"{len(g)} patients with the pseudobulk"
            continue
        cols = [c for c in pg if c in g.columns]
        lc = np.log2(g[cols].to_numpy(float) / g.depth.to_numpy(float)[:, None] * 1e6 + 1)
        nc = pd.Series(g.n_cells).reindex(pts).to_numpy(float)
        comps[key] = dict(z=None, Zp=None, n_cells=nc, cnt=g[cols], lc=pd.DataFrame(lc.astype(np.float32), index=g.index, columns=cols))
    return dict(study=G.donor.iloc[0].split("|")[0], side="tumour", pts=pts, stratum=strat, cells=cells, comps=comps, absent=absent)


def a402_units(d, zhost):
    """{prog: unit} with y of A246 and the host state z over d['pts']."""
    z = pd.Series(zhost).reindex(d["pts"]).to_numpy()
    out = {}
    for prog in P246.PROGS:
        y, why = P246.readout(d, prog)
        out[prog] = (P246.make_unit(d, y, z) if not why else (None, why))
    return out


def a402_null(d, un, n_perm, seed):
    """States shuffled within cohort, one ordering per iteration shared by all programs."""
    R = np.random.default_rng(seed).random((n_perm, len(d["pts"])))
    out = {}
    for k, (u, _) in un.items():
        if u:
            perms = np.argsort(R[:, np.flatnonzero(u["mask"])], axis=1)
            out[k] = (u["n"], np.abs((u["y"][None, :] * np.conj(u["z"])[perms]).sum(1)) / u["den"])
    return out


def a402_table(per_cohort, label, n_perm):
    """per_cohort: {cohort: {prog: (n, null)}} -> program rows and the omnibus row."""
    rows, pooled = [], {}
    for prog in P246.PROGS:
        cs_ = [c for c, v in per_cohort.items() if prog in v]
        ns = [per_cohort[c][prog][0] for c in cs_]
        nul = M.pool_null(ns, [per_cohort[c][prog][1] for c in cs_]) if cs_ else np.array([np.nan])
        pooled[prog] = nul
        sd, n = float(np.std(nul)), int(sum(ns))
        rows.append(dict(set=label, program=prog, patients=n, cohorts=len(cs_), cohort_list=";".join(cs_), null_mean=float(np.mean(nul)),
                         null_sd=sd, mde=MDE_K * sd, n_perm=n_perm, scored=bool(n >= 40 and MDE_K * sd <= 0.25)))
    P = pd.DataFrame(rows)
    sc = P[P.scored].program.tolist()
    if sc:
        T = np.mean([(pooled[p] - pooled[p].mean()) / pooled[p].std() for p in sc], 0)
        om = dict(set=label, program="omnibus", patients=np.nan, cohorts=np.nan, cohort_list=";".join(sc), null_mean=float(T.mean()),
                  null_sd=float(T.std()), mde=MDE_K * float(T.std()), n_perm=n_perm, scored=len(sc) >= 5, n_scored_programs=len(sc))
    else:
        om = dict(set=label, program="omnibus", scored=False, n_scored_programs=0)
    om["decision_mde_max"] = 0.43
    om["decision_meets"] = bool(om.get("scored", False) and om.get("mde", np.inf) <= 0.43)
    return pd.concat([P, pd.DataFrame([om])], ignore_index=True)


def stage_a402(slim, donors, loud, a227, args, out):
    n_perm = 100 if args.smoke else 2000
    z_of = loud.set_index("donor").z
    new, allp = {}, {}
    for s, G in slim.items():
        ids = set(G.donor) & set(z_of.index)
        if len(ids) < 4:
            continue
        d = cohort_dict(G, ids)
        un = a402_units(d, z_of)
        nl = a402_null(d, un, n_perm, [SEED, sum(map(ord, s))])
        allp[s] = nl
        if s not in a227:
            new[s] = nl
    T = pd.concat([a402_table(new, "new relative to A246 (cohorts outside the A227 list)", n_perm),
                   a402_table(allp, "all patients of the P401 cancers (descriptive)", n_perm)], ignore_index=True)
    T.to_csv(out / "a402_power.csv", index=False)
    log("A402\n" + T[["set", "program", "patients", "cohorts", "null_sd", "mde", "scored", "n_scored_programs", "decision_meets"]].round(4).to_string(index=False))
    return T


# ---- A403 --------------------------------------------------------------------------------------------------------------

def readouts(slim, loud):
    """Per donor: immune-hot score, T/NK share, myeloid share, G2M index (A391: malignant, natural log1p CPM)."""
    rows = []
    for s, G in slim.items():
        G = G[G.donor.isin(set(loud.donor))]
        if G.empty:
            continue
        ih = [g for g in IMMUNE_HOT if g in G.columns]
        tot = G.groupby("donor")[ih + ["depth"]].sum() if ih else None
        N = G.groupby("donor").n_cells.sum()
        r = pd.DataFrame(index=sorted(set(G.donor)))
        if tot is not None and len(ih) >= 3:
            lc = np.log2(tot[ih].div(tot.depth, axis=0) * 1e6 + 1)
            z = (lc - lc.mean()) / lc.std(ddof=0).replace(0, np.nan)
            r["immune_hot"] = z.mean(axis=1)
        for nm, comp in (("tnk_share", "tumour:T_NK"), ("myeloid_share", "tumour:myeloid")):
            n = G[G.compartment == comp].groupby("donor").n_cells.sum().reindex(r.index).fillna(0)
            r[nm] = (n / N.reindex(r.index)).where(N.reindex(r.index) >= 100)
        m = G[G.compartment == MALT].drop_duplicates("donor").set_index("donor")
        gg = [g for g in G2M if g in m.columns]
        if len(m) >= 4 and len(gg) >= 10:
            lc = np.log1p(m[gg].div(m.depth, axis=0) * 1e6)
            z = (lc - lc.mean()) / lc.std(ddof=0).replace(0, np.nan)
            r["prolif"] = z.mean(axis=1)
        r["cohort"] = s
        rows.append(r.rename_axis("donor").reset_index())
    return pd.concat(rows, ignore_index=True)


def slope_units(df, col, loud_col="loud"):
    """Per cohort: centred (within stratum) z-scored readout and loudness; units with >= 4 patients and spread."""
    units = []
    for s, g in df.dropna(subset=[col, loud_col]).groupby("cohort"):
        if len(g) < 4:
            continue
        y = centred(g[col].to_numpy(float), g.stratum)
        l = centred(g[loud_col].to_numpy(float), g.stratum)
        if y.std() == 0 or l.std() == 0:
            continue
        units.append(dict(cohort=s, n=len(g), y=y / y.std(), l=l / l.std(), donor=g.donor.to_numpy()))
    return units


def pooled_slope(units, perms=None):
    n = np.array([u["n"] for u in units], float)
    if perms is None:
        r = np.array([np.mean(u["y"] * u["l"]) for u in units])
        return float((n * r).sum() / n.sum())
    r = np.array([np.mean(u["y"][p] * u["l"][None, :], axis=1) for u, p in zip(units, perms)])
    return (n[:, None] * r).sum(0) / n.sum()


def stage_a403(slim, hosts, loud, ro, frames, donors, args, out):
    from clocksync.paired import pseudo_waveforms
    n_perm = 200 if args.smoke else 10_000
    n_set = 6 if args.smoke else 200
    df = loud.merge(ro.drop(columns="cohort"), on="donor", how="inner").assign(cohort=lambda d: d.donor.str.split("|").str[0])
    cnt = df.groupby(["cancer", "cohort"]).agg(patients=("donor", "size"), **{c: (c, lambda x: int(x.notna().sum()))
                                                                          for c in ("immune_hot", "tnk_share", "myeloid_share", "prolif")}).reset_index()
    cnt.to_csv(out / "a403_counts.csv", index=False)
    # companion loudness of 200 matched pseudo-hosts per cancer (same draw as P401: the cancer's pooled anchors)
    comp_loud = {}
    for c, h in hosts.items():
        Gc = h["Gc"]
        genes = [g for g in Gc.columns if g not in META]
        frac = Gc[genes].div(Gc.depth, axis=0)
        lm = np.log(frac.mean() + 1e-12)[(Gc[genes].mean() >= 1) | frac.columns.isin(list(BASE.genes))]
        ws = pseudo_waveforms(lm, n_set, np.random.default_rng([U.SEED, 0]), BASE, "matched")
        used = sorted({g for w in ws for g in w.genes} | set(BASE.genes))
        Gu = Gc[META + used]
        res = Parallel(n_jobs=args.jobs)(delayed(pseudo_host)(Gu, h["classes"], w, h["stratum"]) for w in ws)
        comp_loud[c] = res
        log(f"A403 companion loudness {c}")
    rows = []
    for label, cancers in [("pooled over the P401 cancers", list(hosts))] + [(f"cancer {c}", [c]) for c in hosts]:
        sub = df[df.cancer.isin(cancers)]
        for col, hyp in (("immune_hot", "H403.1"), ("prolif", "H403.2"), ("tnk_share", "descriptive"), ("myeloid_share", "descriptive")):
            un = slope_units(sub, col)
            if not un:
                continue
            rng = np.random.default_rng([403, sum(map(ord, label)), sum(map(ord, col))])
            perms = [np.argsort(rng.random((n_perm, u["n"])), axis=1) for u in un]
            nul = pooled_slope(un, perms)
            sd = float(nul.std())
            row = dict(scope=label, readout=col, hypothesis=hyp, patients=int(sum(u["n"] for u in un)), cohorts=len(un), null_sd=sd,
                       mde_slope_std_units=MDE_K * sd, n_perm=n_perm)
            if label.startswith("pooled") and col in ("immune_hot", "prolif"):
                sl = []
                for j in range(n_set):
                    lj = pd.concat([pd.Series(comp_loud[c][j]) for c in hosts])  # donor -> companion loudness
                    d2 = sub.drop(columns="loud").merge(lj.rename("loud"), left_on="donor", right_index=True)
                    u2 = slope_units(d2, col)
                    r2 = np.random.default_rng([403, 7, j])
                    sh = [r2.permutation(u["n"]) for u in u2]  # one shuffle of the readout per companion set
                    sl.append(float(pooled_slope(u2, [p[None, :] for p in sh])[0]))
                sl = np.array(sl)
                row.update(companion_n_sets=n_set, companion_slope_sd=float(sl.std()), companion_q05=float(np.quantile(sl, .05)),
                           companion_q95=float(np.quantile(sl, .95)), companion_median=float(np.median(sl)))
            rows.append(row)
    T = pd.DataFrame(rows)
    T.to_csv(out / "a403_power.csv", index=False)
    log("A403\n" + T.round(4).to_string(index=False))
    return T


def pseudo_host(Gu, classes, w, stratum):
    limit_threads(1)
    warnings.filterwarnings("ignore")
    h = cs.Dataset(Gu, anchors=classes, waveform=w, mixed_strata=True).host()
    return pd.Series(centred(np.log(h.radius.to_numpy()), stratum.reindex(h.index)), index=h.index)


# ---- A404 --------------------------------------------------------------------------------------------------------------

def sim_power(table, anchors, targets, n_sim, tag, rows, readouts_=("14", "core")):
    for rd in readouts_:
        wf = BASE if rd == "14" else CORE
        ds = cs.Dataset(table, anchors=anchors, waveform=wf, mixed_strata=True)
        hosted = set(ds.host().index.astype(str))
        for t in targets:
            donors_t = set(ds.table.resection[ds.table.compartment == t].astype(str))
            base = dict(host=tag, readout=rd, target=t, patients=len(donors_t), patients_with_host=len(donors_t & hosted),
                        cohorts=len({x.split("|")[0] for x in donors_t}))
            if base["patients_with_host"] < 5:
                rows.append(dict(base, note="< 5 patients with a host"))
                continue
            try:
                sim = cs.simulate.like(ds, t, loading=0.5)
            except Exception as e:  # noqa: BLE001
                rows.append(dict(base, note=f"{type(e).__name__}: {e}"[:150]))
                continue
            rng = np.random.default_rng([404, rd == "core", sum(map(ord, tag + t))])
            reps = {c: np.array([sim.draw(c, rng)["coupling"] for _ in range(n_sim)], float) for c in (0.02, 0.1, 0.2, 0.3)}
            thr = float(np.nanquantile(reps[0.02], .95))
            pw = {c: float(np.nanmean(reps[c] > thr)) for c in (0.1, 0.2, 0.3)}
            ok = [c for c in (0.1, 0.2, 0.3) if pw[c] >= .8]
            rows.append(dict(base, z_abs=0.5, null_threshold=thr, power_c01=pw[0.1], power_c02=pw[0.2], power_c03=pw[0.3],
                             mde_c=min(ok) if ok else np.nan, n_sim=n_sim, note=""))


def comp_chunk(G, anchors, targets, ws, rd, tag):
    limit_threads(1)
    warnings.filterwarnings("ignore")
    out = []
    for i, w in ws:
        try:
            t = cs.Dataset(G, anchors=anchors, waveform=w, mixed_strata=True).fit(n_boot=0, gate=False, targets=targets).table
            out += [dict(host=tag, readout=rd, set=i, target=k, coupling=float(t.loc[k, "coupling"])) for k in targets]
        except Exception as e:  # noqa: BLE001
            out.append(dict(host=tag, readout=rd, set=i, error=f"{type(e).__name__}: {e}"[:150]))
    return out


def companion_q95(G, anchors, targets, tag, args, n_set):
    """200 matched and 200 narrow responsive sets (seeds [404, 0] / [404, 1]) drawn on the pooled anchors of `G`."""
    from clocksync.paired import pseudo_waveforms
    genes = [c for c in G.columns if c not in META and c != "median_umi"]
    genes = [g for g in genes if G[g].notna().all()]
    ref = G[G.compartment.isin(anchors)]
    frac = ref[genes].div(ref.depth, axis=0)
    lm = np.log(frac.mean() + 1e-12)[(ref[genes].mean() >= 1) | frac.columns.isin(list(BASE.genes))]
    sets = [(f, w) for code, f in enumerate(("matched", "responsive"))
            for w in pseudo_waveforms(lm, n_set, np.random.default_rng([404, code]), BASE, f)]
    used = sorted({g for _, w in sets for g in w.genes} | set(BASE.genes))
    Gu = G[META + used]
    tasks = []
    for rd in ("14", "core"):
        ws = [(i, w if rd == "14" else U.C.restrict(w, BASE, CORE)) for i, (_, w) in enumerate(sets)]
        for k in range(args.jobs * 3):
            if ws[k::args.jobs * 3]:
                tasks.append((Gu, anchors, targets, ws[k::args.jobs * 3], rd, tag))
    comp = pd.DataFrame([r for x in Parallel(n_jobs=args.jobs)(delayed(comp_chunk)(*t) for t in tasks) for r in x])
    comp["family"] = [sets[i][0] for i in comp.set]
    return comp


def a387_immune_table(smoke):
    """A387 normal-stroma design with the tumour immune compartments as targets (A227 cohorts and tables)."""
    ANCH_N = ["adjacent_normal:fibroblast", "adjacent_normal:endothelial", "adjacent_normal:mural"]
    T = []
    for s in M.cohort_list(smoke):
        try:
            pb = M.table(s)
        except Exception:  # noqa: BLE001
            continue
        pb = pb[pb.compartment.isin(ANCH_N + IMMUNE)].copy()
        if pb.empty:
            continue
        pb["donor"] = s + "|" + pb.donor.astype(str)
        pb["stratum"] = s + "|" + (pb.stratum.astype(str) if "stratum" in pb else "0")
        T.append(pb[META + list(BASE.genes)] if all(g in pb.columns for g in BASE.genes) else pb[META + [g for g in BASE.genes if g in pb.columns]])
    return pd.concat(T, ignore_index=True), ANCH_N


def genome_a387(study, ANCH_N, pts):
    limit_threads(1)
    warnings.filterwarnings("ignore")
    try:
        G = M.genome_pair(study, tuple(ANCH_N + IMMUNE), pts)
        G["donor"] = study + "|" + G.donor.astype(str)
        return study, G
    except Exception as e:  # noqa: BLE001
        return study, f"{type(e).__name__}: {e}"[:200]


def stage_a404(slim, frames, donors, hosts, args, out):
    n_sim = 20 if args.smoke else 200
    n_set = 6 if args.smoke else 200
    rows = []
    # (a) A387 normal-stroma host
    ta, ANCH_N = a387_immune_table(args.smoke)
    sim_power(ta, ANCH_N, IMMUNE, n_sim, "A387 normal-stroma host", rows)
    log("A404 (a) simulated")
    # (b) tumour-stroma host: the stroma of the P401 cancers (clock genes) plus the immune compartments of the slim tables
    parts = []
    for c, h in hosts.items():
        parts.append(h["tab"])
        h.pop("Gc", None)
    stroma = pd.concat(parts, ignore_index=True)
    imm = pd.concat([G[G.compartment.isin(IMMUNE) & G.donor.isin(set(stroma.donor))][META + list(BASE.genes)] for G in slim.values()
                     if all(g in G.columns for g in BASE.genes)], ignore_index=True)
    tb = pd.concat([stroma, imm], ignore_index=True)
    census = {p.name[:-8] for p in U.CEN.glob("*.parquet")}
    sources_b = [s for s in slim if s in U.GEO_CANCER or s == "bassez2021" or s in census]
    log(f"A404 (b) sources for the genome-wide companions: {sources_b} (slim: {sorted(slim)})")
    slim.clear()
    sim_power(tb, ANCH, IMMUNE, n_sim, "A401 tumour-stroma host (pooled P401 cancers)", rows)
    log("A404 (b) simulated")
    P = pd.DataFrame(rows)
    # companions on the real data (comparators): (a) A227 cohorts via genome_pair; (b) GEO / Census / BIOKEY genome-wide tables
    comp = []
    donors_a = ta.groupby("donor").size().index
    studies = sorted({d.split("|")[0] for d in donors_a})
    from joblib.externals.loky import get_reusable_executor
    get_reusable_executor().shutdown(wait=True)
    res = Parallel(n_jobs=args.genome_jobs, max_nbytes=None)(delayed(genome_a387)(s, ANCH_N, [d.split("|", 1)[1] for d in donors_a if d.startswith(s + "|")])
                                             for s in studies)
    Ga = pd.concat([g for s, g in res if not isinstance(g, str)], ignore_index=True)
    gcols = [c for c in Ga.columns if c not in META]
    Ga[gcols] = Ga[gcols].astype(np.float32)
    Ga = Ga.drop(columns=[c for c in ("stratum",) if c in Ga]).merge(ta[["donor", "compartment", "stratum", "n_cells"]].drop_duplicates(["donor", "compartment"]),
                                                                     on=["donor", "compartment"])
    from joblib.externals.loky import get_reusable_executor
    get_reusable_executor().shutdown(wait=True)
    fa = out / "a404_companions_a.csv"
    if fa.exists():  # part (a) of an earlier run of this stage
        comp.append(pd.read_csv(fa))
    else:
        comp.append(companion_q95(Ga, ANCH_N, [t for t in IMMUNE if t in set(Ga.compartment)], "A387 normal-stroma host", args, n_set))
        comp[-1].to_csv(fa, index=False)
    del Ga, res
    trim()
    log("A404 (a) companions done")
    gen, common = [], None
    for s in sources_b:
        try:
            if s in U.GEO_CANCER:
                G, _ = U.build_geo(s, True, None)
            elif s == "bassez2021":
                G, _ = U.build_bassez(True, None)
            else:
                G, _ = U.build_census(s, True, None)
            G = G[G.donor.isin(set(donors.donor)) & G.compartment.isin(ANCH + IMMUNE)].reset_index(drop=True)
            log(f"A404 genome-wide {s}: {G.shape}, immune rows {int(G.compartment.isin(IMMUNE).sum())}")
            if G.compartment.isin(IMMUNE).any():
                common = list(G.columns) if common is None else [c for c in common if c in set(G.columns)]
                gen = [g[common] for g in gen] + [G[common]]  # inner join kept incrementally (bounded memory)
            del G
            trim()
        except Exception as e:  # noqa: BLE001
            log(f"A404 genome-wide {s}: {type(e).__name__}: {e}"[:200])
    log(f"A404 (b) genome-wide sources: {len(gen)}")
    Gb = pd.concat(gen, ignore_index=True)
    del gen
    trim()
    comp.append(companion_q95(Gb, ANCH, [t for t in IMMUNE if t in set(Gb.compartment)], "A401 tumour-stroma host (GEO / Census / BIOKEY sources only)", args, n_set))
    C_ = pd.concat(comp, ignore_index=True)
    C_.to_csv(out / "a404_companions.csv", index=False)
    S = C_.dropna(subset=["coupling"]).groupby(["host", "readout", "target", "family"]).coupling.agg(
        n_ok="size", median="median", q95=lambda x: float(np.quantile(x, .95))).reset_index()
    S = S.pivot_table(index=["host", "readout", "target"], columns="family", values=["q95", "median", "n_ok"]).reset_index()
    S.columns = ["_".join(c).strip("_") for c in S.columns]
    host_key = {"A387 normal-stroma host": "A387 normal-stroma host"}
    P["_k"] = P.host.str.slice(0, 4)
    S["_k"] = S.host.str.slice(0, 4)
    P = P.merge(S.drop(columns="host"), on=["_k", "readout", "target"], how="left").drop(columns="_k")
    P.to_csv(out / "a404_power.csv", index=False)
    log("A404\n" + P.round(3).to_string(index=False))
    return P


# ---- A405 --------------------------------------------------------------------------------------------------------------

def stage_a405(loud, args, out):
    n_shuf = 200 if args.smoke else 2000
    G = pd.read_parquet(U.BASSEZ, columns=["donor", "sample", "timepoint", "compartment", "n_cells", "expansion", "cohort"])
    G = G[G.n_cells >= U.MIN_CELLS]
    rows = []
    for tp in ("Pre", "On"):
        g = G[G.timepoint == tp]
        mal = set(g[g.compartment == "tumour:malignant"].donor)
        st = set(g[g.compartment.isin(ANCH)].donor)
        pts = mal & st
        rows.append((tp, pts))
    pre, on = rows[0][1], rows[1][1]
    meta = G.drop_duplicates("donor").set_index("donor")[["expansion", "cohort"]]
    out_rows = []
    for label, pts in (("pre-treatment", pre), ("both time points", pre & on), ("on-treatment", on)):
        m = meta.loc[sorted(pts)]
        for (e, c), n in m.groupby(["expansion", "cohort"]).size().items():
            out_rows.append(dict(group=label, expansion=e, cohort=c, patients=int(n)))
        out_rows.append(dict(group=label, expansion="all", cohort="all", patients=len(pts)))
    # MDE of the pre-treatment E vs NE difference in loudness (loudness of the A401 breast host; label shuffles)
    lb = loud[loud.cancer == "BRCA"].set_index("donor")
    pre_ids = [d for d in pre if d in lb.index and meta.expansion.get(d) in ("E", "NE")]
    y = lb.loc[pre_ids, "loud"].to_numpy()
    lab = np.array([meta.expansion[d] for d in pre_ids]) == "E"
    n_e, n_ne = int(lab.sum()), int((~lab).sum())
    res = dict(group="pre-treatment E vs NE (loudness)", patients=len(pre_ids), n_E=n_e, n_NE=n_ne)
    if n_e >= 2 and n_ne >= 2:
        rng = np.random.default_rng([405, 1])
        d = np.array([y[(p := rng.permutation(lab))].mean() - y[~p].mean() for _ in range(n_shuf)])
        sd = float(d.std())
        res.update(null_sd_difference=sd, mde_difference=MDE_K * sd, loudness_sd=float(y.std(ddof=1)), mde_in_loudness_sd=MDE_K * sd / float(y.std(ddof=1)),
                   n_shuffles=n_shuf)
    out_rows.append(res)
    T = pd.DataFrame(out_rows)
    T.to_csv(out / "a405_power.csv", index=False)
    log("A405\n" + T.round(4).to_string(index=False))
    return T


# ---- main --------------------------------------------------------------------------------------------------------------

def stage_power(args):
    OUT_ = OUT / ("smoke" if args.smoke else "")
    OUT_.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    frames, donors, f388 = U.pool_sources(args, OUT_)
    from joblib.externals.loky import get_reusable_executor
    get_reusable_executor().shutdown(wait=True)
    p401 = pd.read_csv(U.OUT / "power" / "p401.csv")
    cancers = p401[(p401.readout.astype(str) == "14") & p401.p401_pass].cancer.tolist()
    # registered set: the matched P401 cancers plus OV, lung_mixed, PAAD, CESC ('unspecified' is not a cancer type)
    cancers = [c for c in cancers if c in set(donors.cancer) and c != "unspecified"]
    log(f"P401 cancers: {cancers}")
    donors = donors[donors.cancer.isin(cancers)].reset_index(drop=True)
    hosts = cancer_hosts(frames, donors, cancers)
    sources, f388 = list(frames), set(f388)  # names only: the frames themselves are released
    frames.clear()  # the genome-wide pooled frames are not needed any more (source names kept for the slim tables)
    trim()
    loud = loudness_frame(hosts)
    log(f"hosts: {len(hosts)} cancers, {len(loud)} patients with a host")
    slim, errs = load_slim(sources, donors, f388, args)
    errs.to_csv(OUT_ / "slim_errors.csv", index=False)
    log(f"slim tables: {len(slim)} sources, errors: {len(errs)}")
    a227 = set(M.cohort_list())
    st = set(args.stages.split(","))
    if "a402" in st:
        stage_a402(slim, donors, loud, a227, args, OUT_)
    if "a403" in st:
        ro = readouts(slim, loud)
        stage_a403(slim, hosts, loud, ro, frames, donors, args, OUT_)
    if "a405" in st:
        stage_a405(loud, args, OUT_)
    get_reusable_executor().shutdown(wait=True)
    if "a404" not in st:
        return
    trim()
    stage_a404(slim, frames, donors, hosts, args, OUT_)
    log(f"power done {(time.time() - t0) / 60:.1f} min")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stage", choices=("power",), required=True)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--genome-jobs", type=int, default=2)
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--stages", default="a402,a403,a405,a404", help="comma list (a rerun of one stage after a failure)")
    args = ap.parse_args()
    limit_threads(1)
    stage_power(args)


if __name__ == "__main__":
    main()
