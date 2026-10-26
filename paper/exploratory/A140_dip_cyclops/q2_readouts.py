"""A140 Q2 parts (i) and (iii) on GSE266919, as registered in paper/docs/registered_experiments.md (A140 Q2 readouts,
dated note 2026-10-26). Part (ii), the CYCLOPS pipeline, is not run here.

Inputs: output/q2/ from `build_q2.py` (pseudobulks, cell counts, sample table). Outputs: output/q2/*.csv, aggregates only
(no per-sample host phases or per-sample rows).

Order:
 (a) orientation check on a simulated timed cohort: is a larger `Dataset.host()` phase_h later in time? (orientation.csv;
     the sign is applied to the host phase in (iii))
 (b) part (i) on the primary set (44 Pre-treatment biopsies, strata = breast vs other) and on the all-78 sensitivity set
     (strata = Group x site): gate (1,000 permutations), A17 learned anchors (frozen learning tables, target study
     zhang2025, tissue breast, pairs from `Dataset.coherence()` of the primary set) and the fixed rule, `Dataset.fit`
     (200 bootstraps), `stability` (200) and `free_agreement`
 (c) part (iii) only if the gate passes and stability >= 0.8 on the set: P140.1-3 (Holm) and the descriptive readouts,
     OLS slope per h of the readout on the unwrapped host phase with one intercept per stratum, 95 % CI from 1,000 sample
     bootstraps with host refit, p from 1,000 within-stratum permutations of the host phase, one-sided in the predicted
     direction. Otherwise "not estimable" rows (never "decoupled").

    python q2_readouts.py [--subsample N] [--jobs J]     (--subsample: check run, 50 permutations, 20 bootstraps)
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[1]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(PAPER / "results" / "extended_data" / "scripts" / "host_identity"))

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

import clocksync as cs
from clocksync import stats as S
from clocksync.anchors import family
from clocksync.simulate import Compartment
from clocksync.waveform import gtex
from common import cli
from common.paths import DATASETS, RESULTS

import anchor_rule as AR  # noqa: E402

COMPS = ["CD4_T", "CD8_T", "myeloid", "B", "NK"]
HOST = RESULTS / "extended_data" / "output" / "host_identity"
MIN_CELLS = 20
APC = ["HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "CD74", "CIITA", "HLA-A", "HLA-B", "HLA-C",
       "B2M", "TAP1", "TAP2", "TAPBP", "PSMB8", "PSMB9"]
CLASS2, CLASS1 = APC[:8], APC[8:]
RECOG = ["CD3D", "CD3E", "LCK", "TRBC1"]
ISG = ["ISG15", "IFIT1", "IFIT2", "IFIT3", "IFI6", "IFI44L", "MX1", "OAS1", "STAT1", "IRF7"]  # choice, not registered
OUT = HERE / "output" / "q2"


# ---- data ------------------------------------------------------------------------------------------------------------

def load(out: Path):
    sm = pd.read_csv(out / "samples.csv", index_col=0)
    cells = pd.read_csv(out / "cells_per_sample.csv", index_col=0)[COMPS]
    sm["site"] = np.where(sm.Tissue == "breast", "breast", "other")
    sm["stratum_all"] = sm.Group + ":" + sm.site
    pbs = {c: pd.read_csv(out / f"pb_{c}.csv.gz", index_col=0) for c in COMPS}
    return sm, cells, pbs


def clock_table(pbs, samples, stratum: pd.Series) -> pd.DataFrame:
    """Long pseudobulk table (clocksync format) of the clock genes; depth = all-gene UMI sum; >= 20 cells."""
    genes = gtex().genes
    rows = []
    for c, pb in pbs.items():
        p = pb.loc[pb.index.intersection(samples)]
        p = p[p.n_cells >= MIN_CELLS]
        cnt = p.drop(columns="n_cells")
        t = pd.DataFrame({"donor": p.index, "compartment": c, "stratum": stratum.reindex(p.index).to_numpy(),
                          "n_cells": p.n_cells.to_numpy(), "depth": cnt.sum(1).to_numpy()})
        sub = pd.DataFrame({g: (cnt[g].to_numpy() if g in cnt else np.zeros(len(cnt))) for g in genes})
        rows.append(pd.concat([t, sub], axis=1))
    return pd.concat(rows, ignore_index=True)


def missing_clock_genes(pbs) -> dict:
    genes = gtex().genes
    return {c: [g for g in genes if g not in pb.columns] for c, pb in pbs.items()}


# ---- (a) orientation -------------------------------------------------------------------------------------------------

def orientation(seeds=(0, 1, 2, 3, 4), n=60) -> pd.DataFrame:
    """Simulated timed cohort (5 compartments, known hours): resultant length of (phase - t) and (phase + t); the larger
    one gives the orientation (phase - t concentrated: a larger host phase_h is later)."""
    wf = gtex()
    rows = []
    for seed in seeds:
        comps = {c: Compartment(z=0.8 + 0.1 * i, psi=0.05) for i, c in enumerate(COMPS)}
        tab, t = cs.simulate(wf, comps, n=n, window_h=(0.0, 24.0), jitter_h=0.0, seed=seed)
        ds = cs.Dataset(tab.rename(columns={"resection": "donor"}), anchors=COMPS)
        h = ds.host()
        tt = pd.Series(t, index=[f"r{i:03d}" for i in range(n)]).reindex(h.index)
        ph, th = h.phase_h.to_numpy() * 2 * np.pi / 24, tt.to_numpy() * 2 * np.pi / 24
        r_same, r_rev = abs(np.exp(1j * (ph - th)).mean()), abs(np.exp(1j * (ph + th)).mean())
        rows.append(dict(seed=seed, n=len(h), R_phase_minus_t=r_same, R_phase_plus_t=r_rev,
                         orientation="later = larger" if r_same > r_rev else "reversed"))
    return pd.DataFrame(rows)


# ---- (b) part (i) ----------------------------------------------------------------------------------------------------

def learning_tables():
    """Frozen learning tables of A17 (as `compute_anchors.learn`): family table of the gate-passing non-target studies."""
    gates = pd.read_csv(HOST / "a17_gates.csv").set_index("name")
    classes = pd.read_csv(HOST / "a17_coherence_classes.csv")
    reg = pd.read_csv(DATASETS / "_atlas" / "a17_targets.csv").set_index("name", drop=False)
    names = [n for n in gates.index if gates.passed[n] and reg.tissue[n] != "multi"]
    ft = AR.family_table(classes[classes.name.isin(names)].rename(columns={"name": "dataset"}))
    ft["study"] = ft.dataset.map(reg.study)
    ft["tissue"] = ft.dataset.map(reg.tissue)
    per = ft.groupby(["study", "tissue", "family"], as_index=False).agg(excess=("excess", "mean"), se=("se", "mean"))
    learning = per.rename(columns={"study": "dataset"})[["dataset", "family", "excess", "se"]]
    datasets = per[["study", "tissue"]].drop_duplicates().assign(dataset=lambda d: d.study, gate_pass=True)
    return learning, datasets


def part_i(tag, table, donors, rc, learning, datasets, pairs_primary=None):
    """Gate, anchors, fit, stability, agreement of one set. Returns (summary dict, tables dict, Dataset or None)."""
    ds0 = cs.Dataset(table, anchors=None, donors=donors, name=tag)
    comps = ds0.compartments
    pairs = pairs_primary if pairs_primary is not None else ds0.coherence(n_perm=rc["n_perm_coh"], n_boot=200)[0]
    anchors, status, fam_tab = AR.select(dict(study="zhang2025", tissue="breast"), comps, pairs, learning, datasets)
    fixed = AR.fixed(dict(tissue="breast"), comps)
    g = ds0.gate(n_perm=rc["n_perm"], seed=0)
    summ = dict(set=tag, n_samples=len(ds0.donor_ids), n_rows=len(ds0.table), compartments="|".join(comps),
                gate_clock_cov=g.clock_cov, gate_clock_share=g.clock_share, gate_share_null_mean=g.share_null_mean,
                gate_p_pairing=g.p_pairing, gate_p_acrophase=g.p_acrophase, gate_n_pairs=g.n_pairs,
                gate_passes=bool(g.passes), anchor_status=status, anchors="|".join(anchors),
                fixed_rule_anchors="|".join(fixed) if fixed else "none")
    tabs = {f"{tag}_families.csv": fam_tab, f"{tag}_pairs.csv": pairs.assign(set=tag)}
    ds = None
    if anchors:
        ds = cs.Dataset(table, anchors=anchors, donors=donors, name=tag)
        st = ds.stability(n_boot=200, seed=0)
        summ.update(stability=st["stability"], median_similarity=st["median_similarity"], stability_n=st["n"],
                    free_agreement=ds.free_agreement())
        if g.passes:
            fit = ds.fit(n_boot=rc["n_boot"], seed=0, gate=g)
            tabs[f"{tag}_fit.csv"] = fit.table.reset_index().assign(set=tag)
    else:
        summ.update(stability=np.nan, median_similarity=np.nan, stability_n=0, free_agreement=np.nan)
    return summ, tabs, ds


# ---- (c) part (iii) --------------------------------------------------------------------------------------------------

def zscore_logcpm(pb, samples, genes):
    """Mean over genes of the z-scored (over samples) log2(CPM + 1); genes absent from the object are skipped."""
    p = pb.loc[samples]
    cnt = p.drop(columns="n_cells")
    have = [g for g in genes if g in cnt.columns]
    lc = np.log2(cnt[have].div(cnt.sum(1), axis=0) * 1e6 + 1)
    z = (lc - lc.mean()) / lc.std(ddof=0).replace(0, np.nan)
    return z.mean(1), have


def readouts(sm, cells, pbs, samples):
    """Per-readout vectors over `samples` (Series, NaN where not defined) and metadata (predicted sign, gene lists)."""
    tot = cells.loc[samples].sum(1)
    logit = lambda c: np.log((cells.loc[samples, c] + 0.5) / (tot - cells.loc[samples, c] + 0.5))
    R, meta = {}, {}
    R["P140.1_CD8_fraction_logit"], meta["P140.1_CD8_fraction_logit"] = logit("CD8_T"), (1, "")
    R["P140.2_myeloid_fraction_logit"], meta["P140.2_myeloid_fraction_logit"] = logit("myeloid"), (-1, "")
    ok_m = cells.loc[samples, "myeloid"] >= MIN_CELLS
    ok = lambda c: cells.loc[samples, c] >= MIN_CELLS

    def score(name, comp, genes, sign):
        s = pd.Series(np.nan, index=samples)
        idx = [x for x in samples if cells.loc[x, comp] >= MIN_CELLS]
        z, have = zscore_logcpm(pbs[comp], idx, genes)
        s.loc[idx] = z
        R[name], meta[name] = s, (sign, "|".join(have))

    score("P140.3_myeloid_antigen_presentation", "myeloid", APC, -1)
    R["D_CD4_fraction_logit"], meta["D_CD4_fraction_logit"] = logit("CD4_T"), (1, "")
    score("D_myeloid_APC_class_II", "myeloid", CLASS2, -1)
    score("D_myeloid_APC_class_I", "myeloid", CLASS1, -1)
    for c in ("CD4_T", "CD8_T"):
        score(f"D_{c}_interferon_response", c, ISG, 1)
        score(f"D_{c}_recognition", c, RECOG, 1)
    return pd.DataFrame(R), meta


def unwrap(ph_h, centre=None):
    ph = np.asarray(ph_h, float) * 2 * np.pi / 24
    mu = np.angle(np.exp(1j * ph).mean()) if centre is None else centre
    return (np.angle(np.exp(1j * (ph - mu))) * 24 / (2 * np.pi))


def slope(u, y, strata):
    """OLS slope of y on u with one intercept per stratum (within-stratum centring)."""
    u, y, s = np.asarray(u, float), np.asarray(y, float), np.asarray(strata)
    ok = np.isfinite(u) & np.isfinite(y)
    u, y, s = u[ok], y[ok], s[ok]
    uc, yc = u.copy(), y.copy()
    for q in np.unique(s):
        m = s == q
        uc[m] -= u[m].mean()
        yc[m] -= y[m].mean()
    den = (uc ** 2).sum()
    return float((uc * yc).sum() / den) if den > 0 else np.nan


def boot_chunk(ds, Y, strat, sign, ids, seeds):
    """Bootstrap slopes (host refit on a donor resample) for every readout; Y: samples x readouts."""
    out = []
    for sd in seeds:
        rng = np.random.default_rng([0, sd])
        try:
            rep = ds.resample(rng)
            h = rep.host()
        except (ValueError, np.linalg.LinAlgError):
            out.append(np.full(Y.shape[1], np.nan))
            continue
        base = h.index.str.split("#").str[0]
        u = sign * unwrap(h.phase_h.to_numpy())
        row = []
        for j in range(Y.shape[1]):
            y = Y.iloc[:, j].reindex(base).to_numpy()
            row.append(slope(u, y, strat.reindex(base).to_numpy()))
        out.append(np.array(row))
    return np.array(out)


def part_iii(tag, ds, sm, cells, pbs, strat_col, orient_sign, rc, jobs):
    h = ds.host()
    keep = [d for d in h.index if (cells.loc[d] >= MIN_CELLS).sum() >= 2]
    h = h.loc[keep]
    strat = sm[strat_col].reindex(h.index)
    Y, meta = readouts(sm, cells, pbs, list(h.index))
    sign = orient_sign
    u = sign * unwrap(h.phase_h.to_numpy())
    ph = np.exp(1j * h.phase_h.to_numpy() * 2 * np.pi / 24)
    hs = dict(set=tag, n_hosted=len(h), host_resultant_R=float(abs(ph.mean())),
              host_circular_sd_h=float(np.sqrt(-2 * np.log(abs(ph.mean()))) * 24 / (2 * np.pi)),
              unwrapped_sd_h=float(np.std(u)), orientation_sign=sign)
    n_boot, n_perm = rc["n_boot_iii"], rc["n_perm_iii"]
    chunks = np.array_split(np.arange(n_boot), max(1, jobs))
    res = Parallel(n_jobs=jobs)(delayed(boot_chunk)(ds, Y, sm[strat_col], sign, None, c) for c in chunks)
    B = np.vstack(res)
    rng = np.random.default_rng(1)
    rows = []
    for j, name in enumerate(Y.columns):
        y = Y[name].to_numpy()
        obs = slope(u, y, strat.to_numpy())
        ok = np.isfinite(y)
        pred = meta[name][0]
        null = []
        for _ in range(n_perm):
            up = u.copy()
            for q in strat.unique():
                m = np.flatnonzero((strat == q).to_numpy())
                up[m] = u[rng.permutation(m)]
            null.append(slope(up, y, strat.to_numpy()))
        null = np.array(null)
        b = B[:, j][np.isfinite(B[:, j])]
        rows.append(dict(set=tag, readout=name, n=int(ok.sum()), predicted_sign=pred, slope_per_h=obs,
                         ci_lo=np.percentile(b, 2.5), ci_hi=np.percentile(b, 97.5), n_boot_ok=len(b),
                         same_direction=bool(np.sign(obs) == pred),
                         p_one_sided=(1 + np.sum(pred * null >= pred * obs)) / (1 + len(null)),
                         p_two_sided=(1 + np.sum(np.abs(null) >= abs(obs))) / (1 + len(null)),
                         genes=meta[name][1]))
    t = pd.DataFrame(rows)
    m = t.readout.str.startswith("P140.")
    p = t.loc[m, "p_one_sided"].to_numpy()
    order = np.argsort(p)
    adj = np.empty(len(p))
    run = 0.0
    for r, i in enumerate(order):
        run = max(run, min(1.0, (len(p) - r) * p[i]))
        adj[i] = run
    t.loc[m, "p_holm"] = adj
    t["replicated"] = np.where(m, t.p_holm < 0.05, np.nan)
    return t, hs


def main(argv=None):
    ap = cli.parser(__doc__)
    args = ap.parse_args(argv)
    run = cli.setup(args)
    data = HERE / "output" / "q2"
    out = data / "subsample" if args.subsample else data
    out.mkdir(exist_ok=True)
    rc = dict(n_perm=run["n_perm"], n_boot=run["n_boot"], n_perm_coh=200 if not run["quick"] else 20,
              n_boot_iii=20 if run["quick"] else 1000, n_perm_iii=50 if run["quick"] else 1000)
    warnings.filterwarnings("ignore")
    sm, cells, pbs = load(data)
    miss = missing_clock_genes(pbs)
    print("clock genes absent from objects (counts set to 0):", {k: v for k, v in miss.items() if v}, flush=True)

    # (a)
    ori = orientation()
    ori.to_csv(out / "orientation.csv", index=False)
    flip = ori.orientation.eq("reversed").sum() > len(ori) / 2
    sign = -1 if flip else 1
    print(ori.to_string(), "\nsign applied:", sign, flush=True)

    # (b)
    learning, datasets = learning_tables()
    pre = sm.index[sm.Group == "Pre-treatment"]
    sets = {"primary": (pre, "site"), "all78": (sm.index, "stratum_all")}
    summaries, hosts, iii = [], [], []
    pairs_primary = None
    for tag, (samples, scol) in sets.items():
        if args.subsample:
            samples = pd.Index(np.random.default_rng(args.seed).permutation(samples)[:args.subsample])
        table = clock_table(pbs, list(samples), sm[scol])
        donors = sm.loc[sm.index.isin(table.donor), :]
        summ, tabs, ds = part_i(tag, table, donors, rc, learning, datasets, None)
        summaries.append(summ)
        for k, v in tabs.items():
            v.to_csv(out / k, index=False)
        print(tag, summ, flush=True)
        ok = summ["gate_passes"] and ds is not None and summ["stability"] >= 0.8
        if ok:
            t, hs = part_iii(tag, ds, sm, cells, pbs, scol, sign, rc, args.jobs)
            iii.append(t)
            hosts.append(hs)
            print(t.drop(columns="genes").to_string(), flush=True)
        else:
            why = ("gate failed" if not summ["gate_passes"] else "no anchors" if ds is None
                   else f"host unstable (stability {summ['stability']:.3f})")
            names = ["P140.1_CD8_fraction_logit", "P140.2_myeloid_fraction_logit", "P140.3_myeloid_antigen_presentation"]
            iii.append(pd.DataFrame(dict(set=tag, readout=names, status=f"not estimable: {why}")))
    pd.DataFrame(summaries).to_csv(out / "part_i_summary.csv", index=False)
    pd.concat(iii, ignore_index=True).to_csv(out / "part_iii_readouts.csv", index=False)
    if hosts:
        pd.DataFrame(hosts).to_csv(out / "part_iii_host_summary.csv", index=False)
    print("done", flush=True)


if __name__ == "__main__":
    main()
