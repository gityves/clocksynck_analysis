"""A412 (registered in `paper/docs/registered_experiments.md`): under intensive care, does the organ's own clock noise psi
stay when the shared rhythm S = z2 V falls (uptake, kappa = 0) or scale with S (kappa = 1)?

GTEx fast deaths vs the 79 ventilated post-mortem donors with a living brain (A349 part 2 design, matched primary,
pooled descriptive), 12 body tissues, host anchored on the 15 A164 tissues, 14 genes (core loop descriptive).
Per tissue dpsi_k = log psi_fast - log psi_vent, dS_k = log S_fast - log S_vent (S = z2 V per group on the pooled-host
contrast, `clocksync._core.host.on_host` as `Dataset.contrast`); companion excess dpsi*_k = dpsi_k - median over the 200
A263 expression-matched GTEx sets of their dpsi_k (same donors, same contrast, set read at the clock positions);
kappa = median_k dpsi*_k / median_k dS_k.

  --stage power  (before any clock psi contrast is looked at; prints no per-tissue psi)
      G412.1  spread-matched null: matched fast controls restricted to the circular time-of-death windows whose host
              spread V (median over tissues, pooled host of the fast controls) is closest to the ventilated V (chosen on
              V alone, within 10 % of it, >= 40 donors, midpoints >= 4 h apart, best 5); pseudo-contrast full fast vs window (nested groups,
              own contrast on the same `on_host`), kappa with 50 A263 sets as companions; bias must lie within +-0.15.
      G412.2  simulation at the observed design (`clocksync._core.host.simulator`, host and donors of each group fixed,
              the group's own S): truth A psi_vent = psi_fast, truth B psi_vent = psi_fast S_vent / S_fast (the
              simulation psi is the observed fast-group psi, never written out), 100 replicates each; the companion
              baseline is a property of the data and cancels, only its 50-of-200 subsampling noise is added; power =
              share of A replicates with kappa < 0.5 - 1.645 SD_A; median kappa under B must be >= 0.7.
      companions: dpsi and dS of the 200 sets (matched and pooled designs) -> output/companions.csv
  --stage score  (written now, NOT run before the power note): kappa and its 90 % interval (500 donor bootstraps within
      groups, host fixed, companion median recomputed in a 50-set subsample per draw) for 14 genes and core loop,
      matched and pooled; descriptive per-tissue table, Spearman of dpsi* against dS, coupling implied at constant psi.
      The core loop uses the 14-gene companions' medians (labelled: "not decided: precision or the removed genes").
Aggregates only (GTEx controlled access): no donor-level output.

    python a412.py --stage power|score [--jobs 4] [--check]
"""
import argparse
import sys
import time
from multiprocessing import Pool

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

HERE_DIR = __import__("pathlib").Path(__file__).resolve().parent
EXP = HERE_DIR.parents[1]
sys.path.insert(0, str(EXP / "A349_icu_limits" / "scripts"))
sys.path.insert(0, str(EXP / "A231_age_coupling" / "scripts"))
import a349_common as C  # noqa: E402
from a349_common import B, H, S, gtex  # noqa: E402
import clocksync as cs  # noqa: E402
from clocksync._core import host as A  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.waveform import Waveform  # noqa: E402
import matched_label as ML  # noqa: E402
import _shared as SH  # noqa: E402

OUT = HERE_DIR.parent / "output"
SETS = EXP / "A263_companions" / "output" / "sets_human.csv"
TRIO = ["PER1", "PER2", "NR1D1"]
WF = cs.gtex()
GROUPS = ("fast", "ventilated")
G = {}  # worker globals: pb, anchors
N_SUB = 50


# ---- data and design --------------------------------------------------------------------------------------------

def setup():
    """Pseudobulk with the clock genes and all companion genes, donor table, matched / pooled designs, companion sets."""
    T, fast, pm = C.tissues()
    ids_f, ids_p = T.donors[fast].astype(str), T.donors[pm].astype(str)
    ph = gtex.phenotypes(B.PH_COLS + ["AGE", "SEX"]).reindex(T.donors)
    attrs = S.sample_attrs(T)
    med = lambda col: np.array([np.nanmedian([attrs[t][col].to_numpy(float)[T.R[t][d]] for t in H.BODY if T.R[t][d] >= 0]  # noqa: E731
                                             or [np.nan]) for d in range(len(T.donors))])
    q = pd.DataFrame(dict(age=pd.to_numeric(ph.AGE, errors="coerce").to_numpy(),
                          male=(pd.to_numeric(ph.SEX, errors="coerce") == 1).astype(int).to_numpy(), hardy=0,
                          ischaemic_h=med("ischemic_time") / 60), index=np.arange(len(T.donors)))
    sets = ML.match_ratio(q, pd.Series(np.isin(q.index, pm), index=q.index), pd.Series(np.isin(q.index, fast), index=q.index))
    cd = T.donors[np.array([s[0] for s in sets])].astype(str)
    ctl = T.donors[np.array(sorted({x for s in sets for x in s[1:]}))].astype(str)
    designs = {"matched": dict(fast=np.asarray(ctl), ventilated=np.asarray(cd)),
               "pooled": dict(fast=np.asarray(ids_f), ventilated=np.asarray(ids_p))}
    cset = pd.read_csv(SETS)
    cset = cset[cset.kind == "matched"].reset_index(drop=True)
    cgenes = [g.split(";") for g in cset.genes]
    allg = sorted(set(WF.genes) | {g for x in cgenes for g in x})
    pb, donors = gtex.load(tissues=list(SH.TISSUES), hardy=(0, 1, 2), genes=allg, depth="library_size", require_tod=True)
    keep = set(ids_f) | set(ids_p)
    pb = pb[pb.donor.isin(keep)].reset_index(drop=True)
    anc = [a for a in SH.TISSUES if a in set(pb.compartment)]
    return pb, donors, anc, designs, cgenes


def dataset(pb, anc, genes=None, drop=None, name="a412"):
    wf = WF if genes is None else Waveform(tuple(genes), WF.amp, WF.psi, None, None, name)
    cols = ["donor", "compartment", "depth"] + list(wf.genes)
    return cs.Dataset(pb[cols], anchors=anc, waveform=wf, name=name, drop_genes=drop)


def restrict(ds, anc, host_ids):
    """Data on the host donors and the hosts of the 12 body tissues (anchors leave-one-out), as `Dataset.contrast`."""
    d = ds.data
    d = d.take(np.flatnonzero(np.isin(d.resection, host_ids)))
    Hs = A._hosts(d, [a for a in anc if a in set(d.compartment)], H.BODY)
    return d, Hs


def tissue_donors(d, Hs, k):
    return pd.Index(np.intersect1d(d.donors(k), Hs[k].resections))


def estimates(d, Hs, groups):
    """{tissue: {group: on_host dict}} for donor groups (may overlap) on the hosts."""
    out = {}
    for k in H.BODY:
        if Hs.get(k) is None or k not in set(d.compartment):
            continue
        idx = tissue_donors(d, Hs, k)
        out[k] = {g: A.on_host(d, k, Hs[k], np.intersect1d(ids, idx)) for g, ids in groups.items()}
    return out


def lratio(a, b):
    return float(np.log(a / b)) if np.isfinite(a) and np.isfinite(b) and a > 0 and b > 0 else np.nan


def dtable(est, a="fast", b="ventilated"):
    """Per tissue: dpsi = log psi_a - log psi_b, dS = log S_a - log S_b (S = z2 V of the group's own estimate)."""
    rows = []
    for k, e in est.items():
        s = [e[g]["z2"] * e[g]["V"] for g in (a, b)]
        rows.append(dict(tissue=k, dpsi=lratio(e[a]["psi"], e[b]["psi"]), dS=lratio(s[0], s[1])))
    return pd.DataFrame(rows).set_index("tissue")


def kappa(dpsi, dS, comp):
    """median(dpsi - comp) / median(dS) over tissues (arrays aligned; comp = companion median per tissue)."""
    num, den = np.nanmedian(np.asarray(dpsi) - np.asarray(comp)), np.nanmedian(dS)
    return float(num / den) if np.isfinite(num) and np.isfinite(den) and den != 0 else np.nan


# ---- companions -----------------------------------------------------------------------------------------------

def comp_task(job):
    """One A263 set: dpsi, dS of every design in `job['designs']` ({name: (host ids, groups)})."""
    i, genes, designs = job
    try:
        ds = dataset(G["pb"], G["anc"], genes, name=f"c{i}")
        rows, cache = [], {}
        for name, (host, groups) in designs.items():
            key = hash(tuple(sorted(host)))
            if key not in cache:
                cache[key] = restrict(ds, G["anc"], host)
            d, Hs = cache[key]
            t = dtable(estimates(d, Hs, groups)).reset_index()
            rows.append(t.assign(set=i, design=name))
        return pd.concat(rows), None
    except Exception as e:  # counted, not dropped silently
        return None, (i, f"{type(e).__name__}: {e}"[:200])


def companions(designs, cgenes, jobs, n_sets, t0, tag):
    todo = [(i, cgenes[i], designs) for i in range(n_sets)]
    parts, failed = [], []
    with Pool(jobs, initializer=limit_threads, initargs=(1,)) as p:
        for k, (r, f) in enumerate(p.imap_unordered(comp_task, todo)):
            parts.append(r) if r is not None else failed.append(f)
            if (k + 1) % 10 == 0 or k + 1 == len(todo):
                print(f"companions {tag}: {k + 1}/{len(todo)} ({time.time() - t0:.0f}s)", flush=True)
    return (pd.concat(parts, ignore_index=True) if parts else pd.DataFrame()), failed


def comp_matrix(ct, design, tissues):
    """sets x tissues matrix of companion dpsi."""
    x = ct[ct.design == design].pivot(index="set", columns="tissue", values="dpsi")
    return x.reindex(columns=tissues)


# ---- bootstrap -------------------------------------------------------------------------------------------------

def boot_kappa(d, Hs, groups, cm, n, seed, n_sub=N_SUB):
    """Donor bootstrap (within groups, host fixed): draws of (kappa, per-tissue dpsi, dS); the companion median is
    recomputed in a `n_sub`-set subsample per draw. `cm`: sets x tissues companion dpsi."""
    rng = np.random.default_rng(seed)
    ks = [k for k in H.BODY if Hs.get(k) is not None and k in set(d.compartment) and k in cm.columns]
    don = {}
    for k in ks:
        idx = tissue_donors(d, Hs, k)
        don[k] = {g: np.intersect1d(ids, idx) for g, ids in groups.items()}
    rows = {k: {r: np.flatnonzero((d.compartment == k) & (d.resection == r)) for g in GROUPS for r in don[k][g]} for k in ks}
    M = cm[ks].to_numpy(float)
    out = []
    for _ in range(n):
        dp, dS = [], []
        for k in ks:
            e = {g: A.on_host_drawn(d, k, Hs[k], rng.choice(don[k][g], len(don[k][g])), rows[k]) for g in GROUPS}
            dp.append(lratio(e["fast"]["psi"], e["ventilated"]["psi"]))
            dS.append(lratio(e["fast"]["z2"] * e["fast"]["V"], e["ventilated"]["z2"] * e["ventilated"]["V"]))
        sub = rng.choice(M.shape[0], min(n_sub, M.shape[0]), replace=False)
        out.append(kappa(dp, dS, np.nanmedian(M[sub], axis=0)))
    return np.array(out)


# ---- G412.1 window null ----------------------------------------------------------------------------------------

def windows(d, Hs, ctl, tod, v_target, min_n=40, n_best=5):
    """Circular time-of-death windows of the fast controls ranked by |log V_window / V_target| (median V over tissues on
    the controls' host), chosen on V alone; greedy with midpoints >= 4 h apart."""
    cand = []
    for w in range(2, 21):
        for s in range(24):
            ids = ctl[((tod.reindex(ctl).to_numpy() - s) % 24) < w]
            if len(ids) < min_n:
                continue
            v = np.nanmedian([A.on_host(d, k, Hs[k], np.intersect1d(ids, tissue_donors(d, Hs, k)))["V"]
                              for k in H.BODY if Hs.get(k) is not None])
            cand.append((abs(np.log(v / v_target)), s, w, v, ids))
    cand = sorted((c for c in cand if c[0] <= np.log(1.1)), key=lambda x: x[0])  # V within 10 % of the target
    pick = []
    for c in cand:
        mid = (c[1] + c[2] / 2) % 24
        if all(min(abs(mid - m), 24 - abs(mid - m)) >= 4 for m in [(p[1] + p[2] / 2) % 24 for p in pick]):
            pick.append(c)
        if len(pick) == n_best:
            break
    return pick


# ---- G412.2 simulation -----------------------------------------------------------------------------------------

def simulate_power(d, Hs, groups, comp_noise, n_rep, seed):
    """Kappa replicates under truth A (psi equal) and B (psi scaled with S) at the observed design."""
    sims, est = {}, {}
    for k in H.BODY:
        if Hs.get(k) is None or k not in set(d.compartment):
            continue
        idx = tissue_donors(d, Hs, k)
        ok = True
        for g in GROUPS:
            s, e = A.simulator(d, k, Hs[k], np.intersect1d(groups[g], idx))
            sims[k, g], est[k, g] = s, e
            ok &= np.isfinite(e["psi"]) and e["z2"] > 0
        if not ok:
            for g in GROUPS:
                sims.pop((k, g)), est.pop((k, g))
    ks = sorted({k for k, _ in sims})
    S0 = {k: {g: est[k, g]["z2"] * est[k, g]["V"] for g in GROUPS} for k in ks}
    psi_f = {k: est[k, "fast"]["psi"] for k in ks}  # simulation input only, never reported
    out = []
    for truth in ("A", "B"):
        rng = np.random.default_rng([seed, ord(truth)])
        psi_v = {k: psi_f[k] * (S0[k]["ventilated"] / S0[k]["fast"] if truth == "B" else 1.0) for k in ks}
        for r in range(n_rep):
            dp, dS = [], []
            for k in ks:
                cf = S0[k]["fast"] / (S0[k]["fast"] + 2 * psi_f[k])
                cv = S0[k]["ventilated"] / (S0[k]["ventilated"] + 2 * psi_v[k])
                a, b = sims[k, "fast"](cf, rng), sims[k, "ventilated"](cv, rng)
                dp.append(lratio(a["psi"], b["psi"]))
                dS.append(lratio(a["z2"] * a["V"], b["z2"] * b["V"]))
            noise = comp_noise(rng)  # per simulated tissue (aligned to ks)
            out.append(dict(truth=truth, rep=r, kappa=kappa(dp, dS, noise),
                            median_dS=float(np.nanmedian(dS))))
    return pd.DataFrame(out), ks, {k: S0[k] for k in ks}


# ---- stages ---------------------------------------------------------------------------------------------------

def power(a):
    limit_threads(4)
    t0 = time.time()
    OUT.mkdir(exist_ok=True)
    pb, donors, anc, designs, cgenes = setup()
    n_sets = 5 if a.check else len(cgenes)
    print(f"data ({time.time() - t0:.0f}s): matched {len(designs['matched']['fast'])} fast / {len(designs['matched']['ventilated'])}"
          f" ventilated; {len(cgenes)} companion sets", flush=True)
    G.update(pb=pb, anc=anc)
    ds = dataset(pb, anc, name="a412_14")
    m = designs["matched"]
    host_m = np.union1d(m["fast"], m["ventilated"])
    # windows (host = matched fast controls); V of the ventilated group on the matched host (V is reported in A349)
    d_f, H_f = restrict(ds, anc, m["fast"])
    d_m, H_m = restrict(ds, anc, host_m)
    est_m = estimates(d_m, H_m, m)
    v_vent = float(np.nanmedian([e["ventilated"]["V"] for e in est_m.values()]))
    v_fast = float(np.nanmedian([e["fast"]["V"] for e in est_m.values()]))
    pick = windows(d_f, H_f, m["fast"], donors.tod, v_vent)
    wdesigns = {f"window{j}": (m["fast"], dict(fast=m["fast"], ventilated=p[4])) for j, p in enumerate(pick)}
    print(f"matched V fast {v_fast:.2f}, ventilated {v_vent:.2f}; windows: "
          + "; ".join(f"[{p[1]}h,+{p[2]}h) n={len(p[4])} V={p[3]:.2f}" for p in pick), flush=True)
    # companions: real designs (all sets) and window designs (first 50 sets)
    real = {n: (np.union1d(x["fast"], x["ventilated"]), x) for n, x in designs.items()}
    ct, failed = companions(real, cgenes, a.jobs, n_sets, t0, "real")
    ct.to_csv(OUT / ("companions_check.csv" if a.check else "companions.csv"), index=False)
    cw, failed_w = companions(wdesigns, cgenes, a.jobs, min(N_SUB, n_sets), t0, "windows")
    cw.to_csv(OUT / ("companions_windows_check.csv" if a.check else "companions_windows.csv"), index=False)
    # G412.1
    rows = []
    for j, p in enumerate(pick):
        grp = dict(fast=m["fast"], ventilated=p[4])
        est = estimates(d_f, H_f, grp)
        t = dtable(est)
        cm = comp_matrix(cw, f"window{j}", list(t.index))
        k_obs = kappa(t.dpsi, t.dS, np.nanmedian(cm.to_numpy(float), axis=0))
        bk = boot_kappa(d_f, H_f, grp, cm, 20 if a.check else 100, [412, 1, j])
        rows.append(dict(window=j, start_h=p[1], width_h=p[2], n_window=len(p[4]), n_full=len(m["fast"]), V_window=p[3],
                         V_target=v_vent, median_dS=float(t.dS.median()), kappa=k_obs, boot_sd=float(np.nanstd(bk, ddof=1)),
                         n_boot_ok=int(np.isfinite(bk).sum())))
    W = pd.DataFrame(rows)
    W.to_csv(OUT / ("power_null_windows_check.csv" if a.check else "power_null_windows.csv"), index=False)
    # G412.2
    cmat = comp_matrix(ct, "matched", H.BODY)
    n_rep = 5 if a.check else 100
    sim, ks, S0 = simulate_power_aligned(d_m, H_m, m, cmat, n_rep, 412)
    sim.to_csv(OUT / ("power_sim_check.csv" if a.check else "power_sim.csv"), index=False)
    sd_a = float(sim[sim.truth == "A"].kappa.std(ddof=1))
    ka, kb = sim[sim.truth == "A"].kappa, sim[sim.truth == "B"].kappa
    summ = dict(n_tissues_sim=len(ks), V_fast_median=v_fast, V_vent_median=v_vent,
                S_ratio_median=float(np.nanmedian([np.log(S0[k]["fast"] / S0[k]["ventilated"]) for k in ks])),
                G412_1_bias=float(W.kappa.median()), G412_1_sd_windows=float(W.kappa.std(ddof=1)) if len(W) > 1 else np.nan,
                G412_1_median_boot_sd=float(W.boot_sd.median()), G412_1_pass=bool(abs(W.kappa.median()) <= 0.15),
                kappa_A_median=float(ka.median()), kappa_A_sd=sd_a,
                power_A=float((ka < 0.5 - 1.645 * sd_a).mean()),
                kappa_B_median=float(kb.median()), kappa_B_sd=float(kb.std(ddof=1)),
                G412_2_power_pass=bool((ka < 0.5 - 1.645 * sd_a).mean() >= 0.8), G412_2_B_pass=bool(kb.median() >= 0.7),
                companions_matched=int(ct[ct.design == "matched"].set.nunique()), companions_pooled=int(ct[ct.design == "pooled"].set.nunique()),
                companions_failed=len(failed), companion_windows_failed=len(failed_w), n_sets_requested=n_sets)
    pd.Series(summ).to_csv(OUT / ("power_summary_check.csv" if a.check else "power_summary.csv"), header=["value"])
    print(pd.Series(summ).to_string(), f"\ndone {time.time() - t0:.0f}s", flush=True)
    if failed or failed_w:
        print("failed draws:", failed[:5], failed_w[:5])


def simulate_power_aligned(d, Hs, groups, cmat, n_rep, seed):
    """`simulate_power` with the companion subsampling noise aligned to the simulated tissues."""
    tis = list(cmat.columns)
    M = cmat.to_numpy(float)
    full = np.nanmedian(M, axis=0)
    holder = {}

    def noise(rng):
        sub = rng.choice(len(M), min(N_SUB, len(M)), replace=False)
        v = pd.Series(np.nanmedian(M[sub], axis=0) - full, index=tis)
        return v.reindex(holder["ks"]).to_numpy()

    class _N:  # the simulated tissues are those with finite S and psi (`_sim_tissues`)
        def __call__(self, rng):
            return noise(rng)

    sims_ks = _sim_tissues(d, Hs, groups)
    holder["ks"] = sims_ks
    return simulate_power(d, Hs, groups, _N(), n_rep, seed)


def _sim_tissues(d, Hs, groups):
    ks = []
    for k in H.BODY:
        if Hs.get(k) is None or k not in set(d.compartment):
            continue
        idx = tissue_donors(d, Hs, k)
        e = [A.on_host(d, k, Hs[k], np.intersect1d(groups[g], idx)) for g in GROUPS]
        if all(np.isfinite(x["psi"]) and x["z2"] > 0 for x in e):
            ks.append(k)
    return ks


def score(a):
    limit_threads(4)
    t0 = time.time()
    pb, donors, anc, designs, cgenes = setup()
    ct = pd.read_csv(OUT / "companions.csv")
    rows, tabs = [], []
    for gs, (genes, drop) in {"14 genes": (None, None), "core": (None, TRIO)}.items():
        ds = dataset(pb, anc, None, drop, name=f"a412_{gs}")
        for design, grp in designs.items():
            d, Hs = restrict(ds, anc, np.union1d(grp["fast"], grp["ventilated"]))
            est = estimates(d, Hs, grp)
            t = dtable(est)
            cm = comp_matrix(ct, design, list(t.index))
            cmed = np.nanmedian(cm.to_numpy(float), axis=0)
            t["comp_dpsi"], t["dpsi_star"] = cmed, t.dpsi - cmed
            bk = boot_kappa(d, Hs, grp, cm, 500, [412, 3, gs == "core", design == "pooled"])
            k = kappa(t.dpsi, t.dS, cmed)
            c_vent = {kk: e["ventilated"]["z2"] * e["ventilated"]["V"] for kk, e in est.items()}
            implied = {kk: c_vent[kk] / (c_vent[kk] + 2 * e["fast"]["psi"]) if kk in c_vent else np.nan for kk, e in est.items()}
            observed = {kk: c_vent[kk] / (c_vent[kk] + 2 * e["ventilated"]["psi"]) for kk, e in est.items()}
            t["c_vent_implied_const_psi"] = pd.Series(implied)
            t["c_vent_observed"] = pd.Series(observed)
            rho = spearmanr(t.dpsi_star, t.dS, nan_policy="omit")
            rows.append(dict(gene_set=gs, design=design, kappa=k, lo90=np.nanquantile(bk, 0.05), hi90=np.nanquantile(bk, 0.95),
                             n_boot_ok=int(np.isfinite(bk).sum()), median_dpsi=float(t.dpsi.median()),
                             median_dpsi_star=float(t.dpsi_star.median()), median_dS=float(t.dS.median()),
                             spearman_dpsi_star_dS=float(rho[0]), spearman_p=float(rho[1]), n_tissues=len(t)))
            tabs.append(t.reset_index().assign(gene_set=gs, design=design))
            print(f"{gs} / {design} ({time.time() - t0:.0f}s)", flush=True)
    pd.DataFrame(rows).to_csv(OUT / "score_kappa.csv", index=False)
    pd.concat(tabs).to_csv(OUT / "score_tissues.csv", index=False)
    print(pd.DataFrame(rows).round(3).to_string(index=False))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=["power", "score"], required=True)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    (power if a.stage == "power" else score)(a)


if __name__ == "__main__":
    main()
