"""A357 (registered): is the signal-gene part of the personal departure a person-level systemic signal that tissues
read unequally?

GTEx fast deaths, A104 periphery rows without Whole Blood and Testis (25 tissues), A294b 14-gene waveform, every tissue
an anchor. Per gene set (the clock, or one of the 100 matched sets of A294b, seed [294, 1]):
- per-gene residuals of y = log((count + 0.5) / depth) on the body's clock state (posterior host mean, strata intercepts),
  host from all tissues except the pair (k, l); corrected signed correlation of the residuals of k and l over the shared
  donors (host-error leakage b_k' C b_l removed, C the mean posterior host covariance); rho_g = median over pairs;
- H357.1: Delta = mean rho over the sig positions (PER1, PER2, NR1D1) minus mean over the 11 core positions;
- H357.2: per tissue k, R^2 of the clock-plane departure d_k (beta minus its 2 x 2 map of the host without k) on the
  person signal factor s(-k) (mean standardised sig residual of the donor's other tissues l, each on the host without
  k and l; donors with >= 3 other tissues).
Stage `power`: the 100 matched sets only (no clock-gene residual is computed). Stage `score`: the clock set, read
against the stored matched sets. Aggregates only (GTEx controlled access): per pair, per gene, per tissue; no per-donor
value is written.

    python score_a357.py --stage power [--jobs 12] [--check]
    python score_a357.py --stage score
Outputs (`../output/`): power_sets.csv, power_tissues.csv, power_genes.csv, power.csv (stage power);
score_genes.csv, score_pairs.csv, score_tissues.csv, score_tests.csv (stage score).
"""
import argparse
import sys
import time
import warnings
from multiprocessing import Pool
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[3]))
sys.path.insert(0, str(HERE.parents[3] / "results" / "figure_2" / "scripts" / "backbone"))
sys.path.insert(0, str(HERE.parents[2] / "A329_core_departure" / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import norm  # noqa: E402

import clocksync as cs  # noqa: E402
import score_a329 as A329  # noqa: E402
from clocksync._core import host as HC  # noqa: E402
from clocksync._core.controls import GC_TARGETS, IMMEDIATE_EARLY  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.genes import canonical  # noqa: E402
from clocksync.paired import pseudo_waveforms  # noqa: E402

OUT = HERE.parents[1] / "output"
META = ["donor", "compartment", "depth"]
EXCLUDE = ["Whole Blood", "Testis"]
SIG = ("PER1", "PER2", "NR1D1")
MIN_SHARED, MIN_OTHER = 30, 3
EFFECT_H1, EFFECT_H2 = 0.10, 0.03
S = {}


def log_y(pb, genes):
    """log((count + 0.5) / depth) per sample, indexed by (donor, compartment)."""
    y = np.log((pb[genes].to_numpy(float) + 0.5) / pb.depth.to_numpy(float)[:, None])
    return pd.DataFrame(y, columns=genes, index=pd.MultiIndex.from_arrays([pb.donor.astype(str), pb.compartment]))


def residuals(d, Y, k, h):
    """Per-gene residuals of tissue k on host h (strata intercepts): DataFrame donors x genes, coefficients (2 x genes),
    per-donor posterior host covariance (n x 2 x 2)."""
    pos = pd.Series(np.arange(len(h.resections)), index=[str(r) for r in h.resections])
    i = np.flatnonzero((d.compartment == k) & np.isin(d.resection.astype(str), pos.index))
    donors = d.resection[i].astype(str)
    j = pos[donors].to_numpy()
    y = Y.loc[list(zip(donors, [k] * len(i)))].to_numpy()
    Z = np.column_stack([pd.get_dummies(pd.Series(d.group[i])).to_numpy(float), h.f[j]])
    coef = np.linalg.lstsq(Z, y, rcond=None)[0]
    e = y - Z @ coef
    return pd.DataFrame(e, index=donors, columns=Y.columns), coef[-2:], h.f_cov[j]


def corrected_r(ea, ba, eb, bb, C):
    """Signed per-gene correlation of two residual matrices over shared donors with the host-error leakage removed."""
    a, b = ea.to_numpy(), eb.to_numpy()
    m = (a * b).mean(0) - np.einsum("ig,ij,jg->g", ba, C, bb)
    va = (a ** 2).mean(0) - np.einsum("ig,ij,jg->g", ba, C, ba)
    vb = (b ** 2).mean(0) - np.einsum("ig,ij,jg->g", bb, C, bb)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where((va > 0) & (vb > 0), m / np.sqrt(va * vb), np.nan)


def departures(d, comps):
    """Per tissue: clock-plane departure d_k (beta minus its 2 x 2 map of the host without k), indexed by donor."""
    out = {}
    for k in comps:
        h = HC.host(d, [c for c in comps if c != k])
        pos = pd.Series(np.arange(len(h.resections)), index=[str(r) for r in h.resections])
        i = np.flatnonzero((d.compartment == k) & np.isin(d.resection.astype(str), pos.index))
        b = d.beta[i].copy()
        g = d.group[i]
        for u in np.unique(g):
            b[g == u] -= b[g == u].mean(0)
        H = h.f[pos[d.resection[i].astype(str)].to_numpy()]
        H = H - H.mean(0)
        dep = b - H @ np.linalg.lstsq(H, b, rcond=None)[0]
        out[k] = pd.DataFrame(dep - dep.mean(0), index=d.resection[i].astype(str), columns=["x", "y"])
    return out


def r2(dep, s):
    """R^2 of the 2-column departure on one factor with an intercept, over the donors present in both."""
    sh = dep.index.intersection(s.index)
    if len(sh) < MIN_SHARED:
        return np.nan, len(sh)
    D = dep.loc[sh].to_numpy()
    D = D - D.mean(0)
    X = np.column_stack([np.ones(len(sh)), s[sh].to_numpy() - s[sh].mean()])
    fit = X @ np.linalg.lstsq(X, D, rcond=None)[0]
    return float(1 - ((D - fit) ** 2).sum() / (D ** 2).sum()), len(sh)


def analyse(name, wf, pb, extra=(), core_wf=None):
    """Pair correlations per gene, person signal factors and R^2 per tissue for one gene set."""
    comps = sorted(pb.compartment.unique())
    genes = list(wf.genes)
    ds = cs.Dataset(pb[META + genes], anchors=comps, waveform=wf, name=name)
    d = ds.data
    Y = log_y(pb, genes + [g for g in extra if g not in genes])
    keyc = pd.MultiIndex.from_arrays([d.resection.astype(str), d.compartment])
    Yd = Y.loc[keyc]
    grp = pd.Series(d.compartment.astype(str)) + "|" + pd.Series(d.group.astype(str))
    Y = pd.DataFrame(Yd.to_numpy() - Yd.groupby(grp.to_numpy()).transform("mean").to_numpy(), index=keyc,
                     columns=Y.columns)
    sig = [g for g in genes if canonical(g) in SIG] if name == "clock" else [genes[k] for k in S["sig_idx"]]
    core = [g for g in genes if g not in sig]
    pair_rows, sig_res, core_res = [], {}, {}
    for a_i, a in enumerate(comps):
        for b in comps[a_i + 1:]:
            h = HC.host(d, [c for c in comps if c not in (a, b)])
            ea, ba, Ca = residuals(d, Y, a, h)
            eb, bb, _ = residuals(d, Y, b, h)
            for k, l, e in ((a, b, eb), (b, a, ea)):  # residuals of l on the host without k and l, standardised
                z = (e - e.mean()) / e.std(ddof=1)
                sig_res[(k, l)] = z[sig].mean(axis=1)
                core_res[(k, l)] = z[core].mean(axis=1)
            sh = ea.index.intersection(eb.index)
            if len(sh) < MIN_SHARED:
                continue
            pa = pd.Series(np.arange(len(ea)), index=ea.index)[sh].to_numpy()
            C = Ca[pa].mean(0)
            r = corrected_r(ea.loc[sh], ba, eb.loc[sh], bb, C)
            pair_rows.append(pd.DataFrame(dict(a=a, b=b, n=len(sh), gene=list(Y.columns), r=r)))
    pairs = pd.concat(pair_rows, ignore_index=True)
    pairs["role"] = np.where(pairs.gene.isin(sig), "sig", np.where(pairs.gene.isin(core), "core", "responsive"))
    pairs["position"] = pairs.gene.map({g: p for p, g in enumerate(genes)})
    deps = departures(d, comps)
    if core_wf is not None:
        dsc = cs.Dataset(pb[META + list(core_wf.genes)], anchors=comps, waveform=core_wf, name=f"{name}_core")
        deps_core = departures(dsc.data, comps)
    tis = []
    for k in comps:
        row = dict(tissue=k)
        for lab, res in (("sig", sig_res), ("core", core_res)):
            parts = [res[(k, l)] for l in comps if l != k]
            cat = pd.concat(parts, axis=1)
            s = cat.mean(axis=1)[cat.notna().sum(axis=1) >= MIN_OTHER]
            row[f"r2_{lab}"], row["n"] = r2(deps[k], s)
            if lab == "sig":
                S.setdefault("factors", {})[k] = s
            if core_wf is not None:
                row[f"r2_{lab}_coredep"], _ = r2(deps_core[k], s)
        tis.append(row)
    return pairs, pd.DataFrame(tis)


def summarise(pairs):
    g = pairs.groupby(["position", "gene", "role"]).r.median().rename("rho").reset_index()
    rs = g[g.role == "sig"].rho.mean()
    rc = g[g.role == "core"].rho.mean()
    return g, dict(rho_sig=float(rs), rho_core=float(rc), delta=float(rs - rc))


def set_task(j):
    warnings.filterwarnings("ignore")
    limit_threads(1)
    try:
        pairs, tis = analyse(f"m{j}", S["sets"][j], S["pb"])
        g, row = summarise(pairs)
        row.update(set=j, median_r2_sig=float(tis.r2_sig.median()), median_r2_core=float(tis.r2_core.median()),
                   boot_sd_r2=float(np.std([np.nanmedian(np.random.default_rng([357, j, b]).choice(
                       tis.r2_sig.to_numpy(), len(tis))) for b in range(500)])), pairs=int(pairs[["a", "b"]].drop_duplicates().shape[0]))
        return row, tis.assign(set=j), g.assign(set=j)
    except Exception as e:  # noqa: BLE001
        return dict(set=j, error=f"{type(e).__name__}: {e}"[:200]), None, None


def boot_median(x, rng, n=2000):
    b = np.array([np.median(rng.choice(x, len(x))) for _ in range(n)])
    return np.quantile(b, [0.05, 0.95])


def prepare(check):
    pb = A329.load()
    pb = pb[~pb.compartment.isin(EXCLUDE)].reset_index(drop=True)
    genes = [c for c in pb.columns if c not in META]
    x = pb[genes].to_numpy(np.float32)
    lm = pd.Series(np.log((x / pb.depth.to_numpy(float)[:, None]).mean(0) + 1e-12), index=genes)[x.mean(0) >= 1]
    del x
    wf = cs.gtex()
    S["sig_idx"] = [k for k, g in enumerate(wf.genes) if canonical(g) in SIG]
    matched = pseudo_waveforms(lm, 100, np.random.default_rng([294, 1]), wf, kind="matched")  # A294b's sets
    S["sets"] = matched[:2] if check else matched
    S["wf"] = wf
    resp = [g for g in dict.fromkeys(GC_TARGETS + IMMEDIATE_EARLY) if g in pb.columns]
    keep = sorted(set(g for w in matched for g in w.genes) | set(wf.genes) | set(resp))
    S["pb"] = pb[META + keep].copy()
    S["resp"] = resp
    return wf


def power_stage(a, t0):
    rows, tis, gen = [], [], []
    with Pool(a.jobs) as pool:
        for k, (row, t, g) in enumerate(pool.imap_unordered(set_task, range(len(S["sets"])))):
            rows.append(row)
            if t is not None:
                tis.append(t)
                gen.append(g)
            print(f"{k + 1}/{len(S['sets'])} sets {time.time() - t0:.0f}s", flush=True)
    sets = pd.DataFrame(rows).sort_values("set")
    err = int(sets.get("error", pd.Series(dtype=str)).notna().sum())
    ok = sets.dropna(subset=["delta"])
    sd1 = max(ok.delta.std(ddof=1), ok.rho_sig.std(ddof=1))
    sd2 = max(ok.median_r2_sig.std(ddof=1), ok.boot_sd_r2.median())
    pw = pd.DataFrame([
        dict(test="H357.1", effect=EFFECT_H1, sd_delta=ok.delta.std(ddof=1), sd_rho_sig=ok.rho_sig.std(ddof=1), sd=sd1,
             power=float(norm.cdf(EFFECT_H1 / sd1 - 1.645)), q95_delta=ok.delta.quantile(0.95),
             q95_rho_sig=ok.rho_sig.quantile(0.95), sets=len(ok), errors=err),
        dict(test="H357.2", effect=EFFECT_H2, sd_median_r2=ok.median_r2_sig.std(ddof=1),
             boot_sd_median_r2=ok.boot_sd_r2.median(), sd=sd2, power=float(norm.cdf(EFFECT_H2 / sd2 - 1.645)),
             q95_median_r2=ok.median_r2_sig.quantile(0.95), matched_median_r2=ok.median_r2_sig.median(), sets=len(ok),
             errors=err)])
    pw["decidable"] = pw.power >= 0.8
    suffix = "_check" if a.check else ""
    sets.to_csv(OUT / f"power_sets{suffix}.csv", index=False)
    if tis:
        pd.concat(tis).to_csv(OUT / f"power_tissues{suffix}.csv", index=False)
        pd.concat(gen).to_csv(OUT / f"power_genes{suffix}.csv", index=False)
    pw.to_csv(OUT / f"power{suffix}.csv", index=False)
    print(pw.round(4).T.to_string(), flush=True)


def score_stage(t0):
    pw = pd.read_csv(OUT / "power.csv").set_index("test")
    ms = pd.read_csv(OUT / "power_sets.csv").dropna(subset=["delta"])
    mt = pd.read_csv(OUT / "power_tissues.csv")
    core_wf = S["wf"].subset([g for g in S["wf"].genes if canonical(g) not in SIG])
    pairs, tis = analyse("clock", S["wf"], S["pb"], extra=S["resp"], core_wf=core_wf)
    g, row = summarise(pairs)
    pairs.to_csv(OUT / "score_pairs.csv", index=False)
    resp = pairs[pairs.role == "responsive"].groupby("gene").r.median().rename("rho")
    g = pd.concat([g, resp.reset_index().assign(
        role=lambda x: np.where(x.gene.isin(GC_TARGETS), "glucocorticoid", "immediate-early"))], ignore_index=True)
    g.to_csv(OUT / "score_genes.csv", index=False)
    # descriptive: person signal factor against the glucocorticoid / immediate-early scores of the other tissues
    sc = A329.responsive_scores(S["pb"])
    sc.index = pd.MultiIndex.from_arrays([S["pb"].donor.astype(str), S["pb"].compartment])
    cors = {}
    for k, s in S["factors"].items():
        oth = sc[sc.index.get_level_values(1) != k].groupby(level=0).mean()
        sh = s.index.intersection(oth.index)
        cors[k] = {f"cor_{c}": float(np.corrcoef(s[sh], oth.loc[sh, c])[0, 1]) for c in ("gc", "ieg")}
    mmed = mt.groupby("tissue").r2_sig.median().rename("r2_matched_median")
    mq95 = mt.groupby("tissue").r2_sig.quantile(0.95).rename("r2_matched_q95")
    tis = tis.set_index("tissue").join(mmed).join(mq95).join(pd.DataFrame(cors).T)
    tis["excess"] = tis.r2_sig - tis.r2_matched_median
    tis.to_csv(OUT / "score_tissues.csv")
    tests = []
    p1 = pw.loc["H357.1"]
    h1 = dict(test="H357.1", decidable=bool(p1.decidable), **row, q95_delta=float(ms.delta.quantile(0.95)),
              q95_rho_sig=float(ms.rho_sig.quantile(0.95)),
              pct_delta=float((ms.delta < row["delta"]).mean()), pct_rho_sig=float((ms.rho_sig < row["rho_sig"]).mean()))
    h1["passes"] = bool(h1["delta"] > h1["q95_delta"] and h1["rho_sig"] > h1["q95_rho_sig"])
    ex = tis.excess.dropna().to_numpy()
    lo, hi = boot_median(ex, np.random.default_rng(357))
    med_r2 = float(tis.r2_sig.median())
    p2 = pw.loc["H357.2"]
    h2 = dict(test="H357.2", decidable=bool(p2.decidable), tissues=len(ex), median_r2=med_r2,
              median_r2_matched=float(ms.median_r2_sig.median()), q95_median_r2=float(ms.median_r2_sig.quantile(0.95)),
              median_excess=float(np.median(ex)), excess_lo90=lo, excess_hi90=hi,
              pct_median_r2=float((ms.median_r2_sig < med_r2).mean()),
              tissues_above_q95=int((tis.r2_sig > tis.r2_matched_q95).sum()))
    h2["passes"] = bool(np.median(ex) > 0 and lo > 0 and med_r2 > h2["q95_median_r2"])
    desc = dict(test="descriptive", median_r2_core_factor=float(tis.r2_core.median()),
                median_r2_sig_on_core_departure=float(tis.r2_sig_coredep.median()),
                median_r2_core_on_core_departure=float(tis.r2_core_coredep.median()),
                median_cor_gc=float(tis.cor_gc.median()), median_cor_ieg=float(tis.cor_ieg.median()),
                rho_gc_median=float(g[g.role == "glucocorticoid"].rho.median()),
                rho_ieg_median=float(g[g.role == "immediate-early"].rho.median()))
    tests = pd.DataFrame([h1, h2, desc])
    tests.to_csv(OUT / "score_tests.csv", index=False)
    print(g.round(3).to_string(index=False), flush=True)
    print(tis.round(3).to_string(), flush=True)
    print(tests.round(4).T.to_string(), flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", choices=["power", "score"], required=True)
    ap.add_argument("--jobs", type=int, default=12)
    ap.add_argument("--check", action="store_true", help="power stage on 2 matched sets (timing only)")
    a = ap.parse_args()
    t0 = time.time()
    prepare(a.check)
    limit_threads(1)
    OUT.mkdir(parents=True, exist_ok=True)
    print(f"{S['pb'].compartment.nunique()} tissues, {len(S['sets'])} matched sets {time.time() - t0:.0f}s", flush=True)
    if a.stage == "power":
        power_stage(a, t0)
    else:
        score_stage(t0)
    print(f"done {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
