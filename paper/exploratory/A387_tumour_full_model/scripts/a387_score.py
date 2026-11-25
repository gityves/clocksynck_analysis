"""A387 score stage: complex-pooled pairing (Part 1) and full model (Part 2) of tumour compartments against the
patient's normal tissue. Registered in paper/docs/registered_experiments.md (A387; power note of 2026-11-25).

Shared definitions (units, R_C, nulls, placebo seeds, pooled table) come from `a387.py`. One scoring run.

Part 1  (A346 states; placebo draws cached by the power stage)
  tests       R_C and pooled phase of the clock row, 14 genes and core loop; within-cohort permutation p (same seeds as the
              power null); companion percentile pi and p_cal per family (`clocksync.pairwise`); Holm over the two stromal
              compartments; phase 90 % interval from 2,000 half-samples without replacement within cohort (shared by both
              readouts; the half-sample variance equals the full-sample one, so the rescaling factor is 1); the A346
              weighted mean of moduli recomputed beside R_C.
  per cohort / per cancer   R_C and phase of the clock row (D387.1).
  trio        14 minus core on the shared half-samples; the 364 three-gene drops from the 14 genes (states recomputed per
              unit); the reliability-predicted core value R_14 x sqrt(rel_T,core rel_N,core / (rel_T,14 rel_N,14)),
              rel = 1 - mean Poisson trace / mean |z|^2, pooled over units with patient weights (D387.2).
Part 2  (one Dataset over all A227 cohorts; anchors = normal fibroblasts, endothelium, mural)
  gate        on the anchors (G387), as the power stage.
  fit         n_boot 500 (90 % intervals), n_perm 1000, anchor_readout "loo", attenuation "df"; offsets and log
              amplitudes against each target's normal counterpart (one fit per reference, same bootstrap seed).
  companions  200 matched and 200 narrow responsive waveforms (`clocksync.paired.pseudo_waveforms`, seeds [387, 0] and
              [387, 1]) drawn on the pooled normal anchors' expression, genome-wide pseudobulks of the anchors and the
              tumour fibroblasts and endothelium (`compute_a227.genome_pair`), fitted without bootstrap; core readout =
              the same sets on the core positions (`a336_common.restrict`). The malignant compartment needs copy-number
              labels that `genome_pair` does not apply: it gets no companions (descriptive, D387.1).
  trio        coupling 14 minus core on the shared bootstrap replicates; the 364 three-gene drops (point estimates);
              host reliability on both readouts.
  pairs       coupling differences between tumour compartments on the patients having both (n_boot 200).
  S387        each stromal target on a host without its normal counterpart (two anchors).
  per cancer  fits on the donors of one cancer type with >= 15 patients of the target (D387.1).

  --smoke  3 cohorts and the pairing scrambled (tumour rows permuted within stratum): a code check that never reads a
           real target value.

    python a387_score.py [--jobs 12] [--genome-jobs 3] [--smoke]     (paper/common/run.sh, RUN_CAP_GB 40; the first run
    with 12 genome workers was killed for memory at 42.6 GB, before any Part 2 output)
"""
from __future__ import annotations

import argparse
import itertools
import pickle
import time
import warnings

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

import a387 as R
from a387 import BASE, BLOCK, C, HOURS, KCODE, M, N_SET, NB, OUT, PRIMARY, TARGETS, NORMAL, ANCHORS, BOUND_H, cs, log
from clocksync.config import limit_threads
from clocksync.pairwise import calibrated_p, percentile

import a149_lib as L  # noqa: E402  (path set by a387)

CORE = C.core(BASE)
FAM_SL = {f: slice(1 + i * N_SET, 1 + (i + 1) * N_SET) for i, f in enumerate(("matched", "responsive", "responsive_broad"))}
DROPS = list(itertools.combinations(BASE.genes, 3))
TRIO_DROP = tuple(g for g in BASE.genes if g in C.TRIO)
N_HALF, N_BOOT, N_PERM2, N_BOOT_SMALL, MIN_CANCER, N_COMP = 2000, 500, 1000, 200, 15, N_SET
SEED = 387


def holm(p):
    return M.holm(np.asarray(p, float))


def wrap_h(x):
    return (np.asarray(x, float) + 12.0) % 24.0 - 12.0


# ---- Part 1 -------------------------------------------------------------------------------------------------------

def load_units(smoke):
    with open(R.A346_CACHE, "rb") as f:
        U = pickle.load(f)
    if smoke:
        rng = np.random.default_rng(0)
        U = {k: v[:3] for k, v in U.items()}
        for units in U.values():
            for u in units:
                u["zb"] = u["zb"][:, rng.permutation(u["n"])]
    return U


def half_samples(units, rows, seed):
    """R_C and phase of `rows` on N_HALF half-samples (floor(n/2) patients per cohort, without replacement)."""
    rng = np.random.default_rng(seed)
    Rs, Ps = [], []
    for _ in range(N_HALF):
        num, ea, eb = 0j, 0.0, 0.0
        for u in units:
            i = rng.choice(u["n"], max(u["n"] // 2, 2), replace=False)
            c, a, b = R._sums(u["za"][rows][:, i], u["zb"][rows][:, i])
            num, ea, eb = num + c, ea + a, eb + b
        Rs.append(np.abs(num) / np.sqrt(ea * eb))
        Ps.append(np.angle(num))
    return np.array(Rs), np.array(Ps)


def mean_of_moduli(units, row):
    """The A346 statistic (patient-weighted mean of per-cohort uncorrected |r|) of one row."""
    v = [(M.raw_stat(u["za"][row], u["zb"][row]), u["n"]) for u in units]
    v = [(s, n) for s, n in v if np.isfinite(s)]
    return float(np.average([s for s, _ in v], weights=[n for _, n in v]))


def unit_drops(study, kind, pts, smoke):
    """States of one unit for the 14 genes, the core loop and the 364 three-gene drops; reliabilities."""
    limit_threads(1)
    warnings.filterwarnings("ignore")
    if kind == "epithelium":
        pb, a = L.table(study)
        b = L.NORM
    else:
        pb = M.table(study)
        a, b = M.PAIRS[kind]
    perm = np.random.default_rng([0, sum(map(ord, study))]).permutation(len(pts)) if smoke else None
    out = {}
    for name, w in [("14", BASE), ("core", CORE)] + [(d, BASE.drop(list(d))) for d in DROPS]:
        p = L.pair_states(pb, a, b, wf=w, pts=pts)
        if p is None:
            out[name] = None
            continue
        zb, vb = (p["zb"][perm], p["vB"][perm]) if smoke else (p["zb"], p["vB"])
        out[name] = (p["za"], zb, p["vA"], vb)
    return study, kind, out


def rel(z, v):
    z = z - z.mean()
    return 1.0 - float(np.mean(v) / np.mean(np.abs(z) ** 2))


def part1(args, cache, out):
    U = load_units(args.smoke)
    tests, per_cohort, per_cancer, trio, drops = [], [], [], [], []
    jobs = [(u["study"], k) for k, units in U.items() for u in units]
    log(f"part 1: states of the {len(DROPS)} drops in {len(jobs)} units")
    res = Parallel(n_jobs=args.jobs)(delayed(_unit_job)(s, k, args.smoke) for s, k in jobs)
    S = {(s, k): o for s, k, o in res}
    for kind, units in U.items():
        z = np.load(cache / f"part1_draws_{kind}.npz")
        D = z["draws"]
        n = int(sum(u["n"] for u in units))
        hs_R, hs_P = half_samples(units, np.array([BLOCK["14"], BLOCK["core"]]), [SEED, KCODE[kind], 2])
        obs = {}
        for rd in BLOCK:
            rows = np.arange(BLOCK[rd], BLOCK[rd] + NB)
            Rr, Pr = R.cpool(units, rows)
            if not args.smoke:
                assert np.allclose(Rr[1:], z[f"comp{rd}"], equal_nan=True), "companions differ from the power stage"
            nul = R.perm_null(units, BLOCK[rd], kind)
            j = 0 if rd == "14" else 1
            dev = wrap_h((hs_P[:, j] - Pr[0]) * HOURS)
            lo, hi = np.quantile(dev, [.05, .95])
            off = float(Pr[0] * HOURS)
            r = dict(kind=kind, readout=rd, cohorts=len(units), patients=n, R_C=float(Rr[0]),
                     p_pairing=float((1 + np.sum(nul >= Rr[0])) / (1 + len(nul))), null_mean=float(nul.mean()),
                     offset_h=off, offset_lo90=off + lo, offset_hi90=off + hi,
                     peak_h=-off, peak_lo90=-(off + hi), peak_hi90=-(off + lo),
                     mean_of_moduli_A346=mean_of_moduli(units, BLOCK[rd]))
            Db = D[:, rows]
            for f, sl in FAM_SL.items():
                q95 = np.nanquantile(Db[:, sl], .95, axis=1)
                r.update({f"{f}_q95": float(np.nanquantile(Rr[sl], .95)), f"{f}_pi": float(percentile(Rr, sl)[0]),
                          f"{f}_pcal": float(calibrated_p(percentile(Rr, sl)[0], percentile(Db, sl))[0]),
                          f"{f}_placebo_rate": float(np.mean(Db[:, 0] > q95))})
            tests.append(r)
            obs[rd] = r
            for u in units:
                Ru, Pu = R.cpool([u], np.array([BLOCK[rd]]))
                per_cohort.append(dict(kind=kind, readout=rd, study=u["study"], cancer=M.cancer(u["study"]), n=u["n"],
                                       R_C=float(Ru[0]), peak_h=float(-Pu[0] * HOURS)))
            for cn in sorted({M.cancer(u["study"]) for u in units}):
                us = [u for u in units if M.cancer(u["study"]) == cn]
                Rc, Pc = R.cpool(us, np.array([BLOCK[rd]]))
                per_cancer.append(dict(kind=kind, readout=rd, cancer=cn, cohorts=len(us), n=int(sum(u["n"] for u in us)),
                                       R_C=float(Rc[0]), peak_h=float(-Pc[0] * HOURS)))
        # trio (D387.2): shared half-samples, 364 drops, reliability prediction
        d_obs = obs["14"]["R_C"] - obs["core"]["R_C"]
        d_b = hs_R[:, 0] - hs_R[:, 1]
        lo, hi = np.quantile(d_b - d_b.mean(), [.05, .95])
        st = [S[(u["study"], kind)] for u in units]
        dv = {}
        for d in DROPS:
            num, ea, eb = 0j, 0.0, 0.0
            for o in st:
                if o.get(d) is None:
                    continue
                c, a, b = R._sums(o[d][0][None], o[d][1][None])
                num, ea, eb = num + c[0], ea + a[0], eb + b[0]
            dv[d] = float(np.abs(num) / np.sqrt(ea * eb))
            drops.append(dict(part=1, kind=kind, drop="+".join(d), value=dv[d]))
        vals = np.array(list(dv.values()))
        tr = dv[TRIO_DROP]
        if not args.smoke:
            assert np.isclose(tr, obs["core"]["R_C"], atol=1e-6), "trio drop differs from the core readout"
        w = np.array([u["n"] for u in units], float)
        rl = {rd: [(rel(o[rd][0], o[rd][2]), rel(o[rd][1], o[rd][3])) for o in st] for rd in ("14", "core")}
        fac = np.sqrt(np.average([a * b for a, b in rl["core"]], weights=w) / np.average([a * b for a, b in rl["14"]], weights=w))
        trio.append(dict(part=1, kind=kind, value_14=obs["14"]["R_C"], value_core=obs["core"]["R_C"], diff=d_obs,
                         diff_lo90=d_obs + lo, diff_hi90=d_obs + hi, trio_drop_percentile=float(np.mean(vals <= tr)),
                         drops_median=float(np.median(vals)), reliability_factor=float(fac),
                         predicted_core=float(obs["14"]["R_C"] * fac)))
        log(f"part 1 {kind} scored")
    T = pd.DataFrame(tests)
    for col in ("p_pairing", "matched_pcal", "responsive_pcal", "responsive_broad_pcal"):
        T[f"{col}_holm"] = np.nan
        for rd in BLOCK:
            sel = T.index[(T.readout == rd) & T.kind.isin(PRIMARY)]
            T.loc[sel, f"{col}_holm"] = holm(T.loc[sel, col])
    prim = T.kind.isin(PRIMARY)
    T["H387_1"] = prim & (T.p_pairing_holm < 0.05) & (T.matched_pcal_holm < 0.05)
    T["H387_2"] = T.H387_1 & (T.peak_lo90 > -BOUND_H) & (T.peak_hi90 < BOUND_H)
    T["H387_2_straddles"] = T.H387_1 & ~T.H387_2 & (T.peak_lo90 < BOUND_H) & (T.peak_hi90 > -BOUND_H)
    T.to_csv(out / "part1_tests.csv", index=False)
    pd.DataFrame(per_cohort).to_csv(out / "part1_per_cohort.csv", index=False)
    pd.DataFrame(per_cancer).to_csv(out / "part1_per_cancer.csv", index=False)
    return T, pd.DataFrame(trio), pd.DataFrame(drops)


def _unit_job(study, kind, smoke):
    """Patients of a unit as A346 builds it, then its drop states."""
    limit_threads(1)
    warnings.filterwarnings("ignore")
    if kind == "epithelium":
        pb, a = L.table(study)
        b = L.NORM
    else:
        pb = M.table(study)
        a, b = M.PAIRS[kind]
    u = L.pair_states(pb, a, b)
    return unit_drops(study, kind, u["pts"], smoke)


# ---- Part 2 -------------------------------------------------------------------------------------------------------

def scramble(table, smoke):
    if not smoke:
        return table
    t = table.copy()
    rng = np.random.default_rng(1)
    for k in TARGETS.values():
        for st, g in t[t.compartment == k].groupby("stratum"):
            t.loc[g.index, "donor"] = g.donor.to_numpy()[rng.permutation(len(g))]
    return t


def dataset(table, wf, anchors=ANCHORS):
    return cs.Dataset(table, anchors=anchors, waveform=wf, mixed_strata=True)


def main_fit(table, rd, ref, n_perm):
    limit_threads(1)
    warnings.filterwarnings("ignore")
    f = dataset(table, R.READOUTS[rd]).fit(n_boot=N_BOOT, seed=SEED, n_perm=n_perm, gate=False, reference=ref,
                                           attenuation="df", level=0.90)
    return rd, ref, f.table, f.replicates


def genome_cohort(study, keep):
    """Genome-wide pseudobulks of the stromal compartments of one cohort, rows restricted to `keep` (donor, comp)."""
    limit_threads(1)
    warnings.filterwarnings("ignore")
    comps = tuple(ANCHORS + [TARGETS["fibroblast"], TARGETS["endothelial"]])
    pts = sorted({d.split("|", 1)[1] for d, _ in keep})
    if not pts:
        return None
    G = M.genome_pair(study, comps, pts)
    G["donor"] = study + "|" + G.donor.astype(str)
    G = G[[(d, c) in keep for d, c in zip(G.donor, G.compartment)]]
    return G


def companion_fits(G, wfs, rd):
    limit_threads(1)
    warnings.filterwarnings("ignore")
    out = []
    for i, w in wfs:
        try:
            t = dataset(G, w).fit(n_boot=0, gate=False, targets=[TARGETS["fibroblast"], TARGETS["endothelial"]]).table
            out.append(dict(readout=rd, set=i, **{k: float(t.loc[TARGETS[k], "coupling"]) for k in PRIMARY}))
        except Exception as e:  # noqa: BLE001
            out.append(dict(readout=rd, set=i, error=f"{type(e).__name__}: {e}"[:200]))
    return out


def drop_fits(table, drops):
    limit_threads(1)
    warnings.filterwarnings("ignore")
    out = []
    for d in drops:
        t = dataset(table, BASE.drop(list(d))).fit(n_boot=0, gate=False, targets=list(TARGETS.values())).table
        out += [dict(part=2, kind=k, drop="+".join(d), value=float(t.loc[v, "coupling"])) for k, v in TARGETS.items()]
    return out


def chunks(x, n):
    return [x[i::n] for i in range(n) if x[i::n]]


def part2(args, out):
    table = scramble(R.pooled_table(args.jobs, args.smoke), args.smoke)
    gates = []
    for rd, wf in R.READOUTS.items():
        g = dataset(table, wf).subset(compartments=ANCHORS).gate(n_perm=1000, seed=SEED)
        gates.append(dict(readout=rd, p_pairing=g.p_pairing, p_acrophase=g.p_acrophase, passes=bool(g.passes)))
    gates = pd.DataFrame(gates).set_index("readout")
    log("part 2 gate\n" + gates.to_string())
    # main fits: three references x two readouts, same bootstrap seed
    jobs = [(rd, NORMAL[k], N_PERM2 if k == "fibroblast" else 0) for rd in R.READOUTS for k in TARGETS]
    res = Parallel(n_jobs=min(args.jobs, len(jobs)))(delayed(main_fit)(table, *j) for j in jobs)
    fits = {(rd, ref): (t, rep) for rd, ref, t, rep in res}
    log("part 2 main fits done")
    # companions
    keep = set(zip(table.donor, table.compartment))
    studies = sorted({d.split("|", 1)[0] for d in table.donor})
    Gs = Parallel(n_jobs=args.genome_jobs)(delayed(genome_cohort)(s, {x for x in keep if x[0].startswith(s + "|")}) for s in studies)
    G = pd.concat([g for g in Gs if g is not None and len(g)], ignore_index=True)
    G = G.drop(columns=["stratum"]).merge(table[["donor", "compartment", "stratum"]], on=["donor", "compartment"])
    G = scramble(G, args.smoke)
    genes = [c for c in G.columns if c not in ("donor", "compartment", "stratum", "depth", "n_cells")]
    genes = [g for g in genes if G[g].notna().all()]   # genes measured in every cohort
    ref = G[G.compartment.isin(ANCHORS)]
    frac = ref[genes].div(ref.depth, axis=0)
    lm = np.log(frac.mean() + 1e-12)[(ref[genes].mean() >= 1) | frac.columns.isin(list(BASE.genes))]
    from clocksync.paired import pseudo_waveforms
    sets = [(f, w) for code, f in enumerate(("matched", "responsive"))
            for w in pseudo_waveforms(lm, N_COMP, np.random.default_rng([SEED, code]), BASE, f)]
    used = sorted({g for _, w in sets for g in w.genes} | set(BASE.genes))
    Gu = G[["donor", "compartment", "stratum", "depth"] + used]
    log(f"part 2 companions: {len(sets)} sets, {len(used)} genes, {len(Gu)} rows")
    tasks = [(rd, [(i, w if rd == "14" else C.restrict(w, BASE, CORE)) for i, (_, w) in enumerate(sets)]) for rd in R.READOUTS]
    comp = Parallel(n_jobs=args.jobs)(delayed(companion_fits)(Gu, ch, rd) for rd, ws in tasks for ch in chunks(ws, args.jobs))
    comp = pd.DataFrame([r for x in comp for r in x])
    comp["family"] = [sets[i][0] for i in comp.set]
    comp.to_csv(out / "part2_companions.csv", index=False)
    log("part 2 companions done")
    # tests
    rows = []
    for rd in R.READOUTS:
        t0, _ = fits[(rd, NORMAL["fibroblast"])]
        for k, v in TARGETS.items():
            tr, _ = fits[(rd, NORMAL[k])]
            r = dict(kind=k, readout=rd, gate_passes=bool(gates.loc[rd, "passes"]), n=int(t0.loc[v, "n"]),
                     coupling=t0.loc[v, "coupling"], coupling_lo90=t0.loc[v, "coupling_lo"],
                     coupling_hi90=t0.loc[v, "coupling_hi"], p_coupling=t0.loc[v, "p_coupling"],
                     boundary_type=t0.loc[v, "boundary_type"], host_reliability=t0.loc[v, "host_reliability"],
                     coupling_corrected=t0.loc[v, "coupling_corrected"],
                     normal_coupling=t0.loc[NORMAL[k], "coupling"],
                     offset_h=tr.loc[v, "offset_h"], offset_lo90=tr.loc[v, "offset_lo"], offset_hi90=tr.loc[v, "offset_hi"],
                     log_amp_ratio=tr.loc[v, "log_amp_rel"], log_amp_lo90=tr.loc[v, "log_amp_lo"],
                     log_amp_hi90=tr.loc[v, "log_amp_hi"], n_boot_failed=int(tr.loc[v, "n_boot_failed"]))
            r.update(peak_h=-r["offset_h"], peak_lo90=-r["offset_hi90"], peak_hi90=-r["offset_lo90"])
            if k in PRIMARY:
                for f in ("matched", "responsive"):
                    c = comp[(comp.readout == rd) & (comp.family == f)][k]
                    r[f"{f}_q95"], r[f"{f}_median"] = float(np.nanquantile(c, .95)), float(np.nanmedian(c))
                    r[f"{f}_share_ge"] = float(np.nanmean(c >= r["coupling"]))
                    r[f"{f}_n_ok"] = int(c.notna().sum())
            rows.append(r)
    T = pd.DataFrame(rows)
    T["p_coupling_holm"] = np.nan
    for rd in R.READOUTS:
        sel = T.index[(T.readout == rd) & T.kind.isin(PRIMARY)]
        T.loc[sel, "p_coupling_holm"] = holm(T.loc[sel, "p_coupling"])
    prim = T.kind.isin(PRIMARY)
    T["H387_3"] = prim & T.gate_passes & (T.p_coupling_holm < 0.05) & (T.coupling > T.matched_q95)
    T["H387_4"] = T.H387_3 & (T.peak_lo90 > -BOUND_H) & (T.peak_hi90 < BOUND_H)
    T.to_csv(out / "part2_tests.csv", index=False)
    # trio (D387.2)
    trio, drops = [], []
    r14, rc = fits[("14", NORMAL["fibroblast"])][1], fits[("core", NORMAL["fibroblast"])][1]
    dres = Parallel(n_jobs=args.jobs)(delayed(drop_fits)(table, ch) for ch in chunks(DROPS, args.jobs))
    drops = pd.DataFrame([r for x in dres for r in x])
    for k, v in TARGETS.items():
        a = T[(T.kind == k) & (T.readout == "14")].iloc[0]
        b = T[(T.kind == k) & (T.readout == "core")].iloc[0]
        dv = drops[drops.kind == k].set_index("drop").value
        tr = dv["+".join(TRIO_DROP)]
        aligned = len(r14) == len(rc) and a.n_boot_failed == 0 and b.n_boot_failed == 0
        d = (r14[v] - rc[v]).to_numpy() if aligned else np.full(1, np.nan)
        trio.append(dict(part=2, kind=k, value_14=a.coupling, value_core=b.coupling, diff=a.coupling - b.coupling,
                         diff_lo90=float(np.nanquantile(d, .05)), diff_hi90=float(np.nanquantile(d, .95)),
                         replicates_aligned=aligned, trio_drop_percentile=float(np.mean(dv <= tr)),
                         drops_median=float(dv.median()), host_reliability_14=a.host_reliability,
                         host_reliability_core=b.host_reliability))
    log("part 2 drops done")
    # pairs on common patients
    pairs = []
    ks = list(TARGETS)
    for rd in R.READOUTS:
        ds = dataset(table, R.READOUTS[rd])
        for x, y in itertools.combinations(ks, 2):
            a, b = TARGETS[x], TARGETS[y]
            both = sorted(set(table.donor[table.compartment == a]) & set(table.donor[table.compartment == b]))
            try:
                f = ds.subset(donors=both).fit(n_boot=N_BOOT_SMALL, seed=SEED, gate=False, targets=[a, b], level=0.90)
                c = f.compare({x: [a], y: [b]})
                pairs.append(dict(readout=rd, a=x, b=y, n_both=len(both), c_a=f.table.loc[a, "coupling"],
                                  c_b=f.table.loc[b, "coupling"], **{k: c.iloc[-1][k] for k in c.columns if k != "set"}))
            except Exception as e:  # noqa: BLE001
                pairs.append(dict(readout=rd, a=x, b=y, n_both=len(both), error=str(e)[:200]))
    # S387: host without the normal counterpart
    sens = []
    for rd in R.READOUTS:
        for k in PRIMARY:
            anc = [x for x in ANCHORS if x != NORMAL[k]]
            try:
                f = dataset(table, R.READOUTS[rd], anc).fit(n_boot=N_BOOT_SMALL, seed=SEED, gate=False,
                                                           targets=[TARGETS[k], NORMAL[k]], reference=NORMAL[k], level=0.90)
                t = f.table.loc[TARGETS[k]]
                sens.append(dict(readout=rd, kind=k, anchors="+".join(anc), n=int(t.n), coupling=t.coupling,
                                 coupling_lo90=t.coupling_lo, coupling_hi90=t.coupling_hi, peak_h=-t.offset_h,
                                 peak_lo90=-t.offset_hi, peak_hi90=-t.offset_lo, log_amp_ratio=t.log_amp_rel))
            except Exception as e:  # noqa: BLE001
                sens.append(dict(readout=rd, kind=k, error=str(e)[:200]))
    # per cancer (D387.1)
    canc = []
    cancer_of = {d: M.cancer(d.split("|", 1)[0]) for d in table.donor.unique()}
    for rd in R.READOUTS:
        for k, v in TARGETS.items():
            dk = table.donor[table.compartment == v].unique()
            for cn in sorted({cancer_of[d] for d in dk}):
                don = [d for d in table.donor.unique() if cancer_of[d] == cn]
                n_t = sum(cancer_of[d] == cn for d in dk)
                if n_t < MIN_CANCER:
                    continue
                try:
                    f = dataset(table[table.donor.isin(don)], R.READOUTS[rd]).fit(
                        n_boot=N_BOOT_SMALL, seed=SEED, gate=False, targets=[v, NORMAL[k]], reference=NORMAL[k], level=0.90)
                    t = f.table.loc[v]
                    canc.append(dict(readout=rd, kind=k, cancer=cn, n=int(t.n), coupling=t.coupling,
                                     coupling_lo90=t.coupling_lo, coupling_hi90=t.coupling_hi, peak_h=-t.offset_h,
                                     peak_lo90=-t.offset_hi, peak_hi90=-t.offset_lo))
                except Exception as e:  # noqa: BLE001
                    canc.append(dict(readout=rd, kind=k, cancer=cn, n=n_t, error=str(e)[:200]))
    pd.DataFrame(pairs).to_csv(out / "part2_pairs.csv", index=False)
    pd.DataFrame(sens).to_csv(out / "part2_sensitivity.csv", index=False)
    pd.DataFrame(canc).to_csv(out / "part2_per_cancer.csv", index=False)
    return T, pd.DataFrame(trio), drops


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--jobs", type=int, default=12)
    ap.add_argument("--genome-jobs", type=int, default=3, help="workers for the genome-wide pseudobulks (memory-bound)")
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()
    limit_threads(1)
    global N_HALF, N_BOOT, N_PERM2, N_BOOT_SMALL, N_COMP, DROPS
    if args.smoke:
        R.N_PERM = 200
        N_HALF, N_BOOT, N_PERM2, N_BOOT_SMALL, N_COMP = 50, 20, 50, 10, 6
        DROPS = DROPS[:6] + [TRIO_DROP]
    cache = OUT / ("cache_smoke" if args.smoke else "cache")
    out = OUT / ("score_smoke" if args.smoke else "score")
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    T1, trio1, drops1 = part1(args, cache, out)
    log("part 1 done")
    T2, trio2, drops2 = part2(args, out)
    pd.concat([trio1, trio2]).to_csv(out / "trio.csv", index=False)
    pd.concat([drops1, drops2]).to_csv(out / "drops.csv", index=False)
    pd.set_option("display.width", 250)
    log("part 1\n" + T1.round(4).to_string(index=False))
    log("part 2\n" + T2.round(4).to_string(index=False))
    log("trio\n" + pd.concat([trio1, trio2]).round(4).to_string(index=False))
    log(f"done {(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    main()
