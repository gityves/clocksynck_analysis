"""A336c (descriptive check of A336 part 2, registered 2026-11-26): is the GTEx person loudness the time of death?

A336 part 2 leverage (`clocksync.elasticity.radius_covariance` on the 20 halvings of the 27 peripheral tissues, fast
deaths, median over halvings), with in addition: V1 harmonics 1-4 of the recorded time of death; V2 2 h time-of-death
bins (primary); V3 1 h bins. Clock (14 genes, core loop) with a donors-permuted null within 2 h bins; the first 100
matched and 100 responsive A336 sets adjusted the same way. Descriptive: R^2 of the half-host radii on the 2 h bins
(raw, and after A336's phase harmonics), against the chance level (k - 1) / (n - 1) and the companions' R^2.
Aggregates only (time of death is controlled access).

    RUN_CAP_GB=20 paper/common/run.sh paper/exploratory/A336_two_readouts/scripts/a336c_time_of_death.py --jobs 6
    python a336c_time_of_death.py --check        # 2 sets, 20 permutations; nothing written
Outputs (`../output/`): c_leverages.csv, c_summary.csv.
"""
import argparse
import sys
import time
import warnings
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import a336_person_gtex as G  # noqa: E402
from clocksync import elasticity as E  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.paired import pseudo_waveforms  # noqa: E402

OUT = HERE.parents[1] / "output"
VARIANTS = ("V0 A336", "V1 time-of-death harmonics 1-4", "V2 2 h bins", "V3 1 h bins")
N_SETS, N_PERM = 100, 1000
A336 = {"14 genes": 0.765144, "core": 0.713847}  # person_gtex.csv, 6 decimals


def adjust(j, tod, variant):
    """(covariates, groups) of `radius_covariance` for donors j."""
    t = tod.reindex(j).to_numpy(float)
    if variant == "V1 time-of-death harmonics 1-4":
        w = 2 * np.pi * t / 24
        return np.column_stack([g(k * w) for k in range(1, 5) for g in (np.cos, np.sin)]), None
    if variant == "V2 2 h bins":
        return None, np.floor(t / 2).astype(int)
    if variant == "V3 1 h bins":
        return None, np.floor(t).astype(int)
    return None, None


def leverage(st, tod, variant, perm=None):
    L = []
    for h, (j, a, b) in enumerate(st):
        if perm is not None:
            b = b[perm[h]]
        X, g = adjust(j, tod, variant)
        L.append(E.radius_covariance(a, b, covariates=X, groups=g)["leverage"])
    return float(np.median(L))


def within_bin_perms(st, tod, rng, n):
    """Per halving, a permutation of B among donors of the same 2 h bin."""
    out = []
    for _ in range(n):
        ps = []
        for j, _, _ in st:
            b = np.floor(tod.reindex(j).to_numpy(float) / 2)
            p = np.arange(len(j))
            for v in np.unique(b):
                m = np.flatnonzero(b == v)
                p[m] = m[rng.permutation(len(m))]
            ps.append(p)
        out.append(ps)
    return out


def radius_r2(st, tod):
    """Median over halvings and halves of the R^2 of the radius on 2 h bins: raw, and after A336's phase harmonics."""
    raw, res = [], []
    for j, a, b in st:
        B = pd.get_dummies(pd.Series(np.floor(tod.reindex(j).to_numpy(float) / 2))).to_numpy(float)
        H = np.column_stack([np.ones(len(a))] + [g(k * p) for p in (np.angle(a), np.angle(b)) for k in (1, 2)
                                                 for g in (np.cos, np.sin)])
        for z in (a, b):
            r = np.abs(z)
            e = r - H @ np.linalg.lstsq(H, r, rcond=None)[0]
            for v, acc in ((r, raw), (e, res)):
                f = v - B @ np.linalg.lstsq(B, v, rcond=None)[0]
                acc.append(1 - f.var() / v.var())
    return float(np.median(raw)), float(np.median(res))


def main():
    warnings.simplefilter("ignore")
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--jobs", type=int, default=6)
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    limit_threads(1)
    t0 = time.time()
    pb, donors = G.load()
    comps = sorted(pb.compartment.unique())
    genes = [c for c in pb.columns if c not in G.META]
    x = pb[genes].to_numpy(np.float32)
    lm = pd.Series(np.log((x / pb.depth.to_numpy(float)[:, None]).mean(0) + 1e-12), index=genes)[x.mean(0) >= 1]
    lm = lm[~lm.index.duplicated()]
    del x
    # the A336 sets (same seeds, 200 drawn), first N_SETS of each kind
    sets = {kind: pseudo_waveforms(lm, G.N_SETS, np.random.default_rng([336, 2, tag]), G.WF, kind=kind)[:N_SETS]
            for kind, tag in (("matched", 1), ("responsive", 2))}
    n = 2 if a.check else N_SETS
    core_of = lambda w, j: G.Waveform(tuple(w.genes[k] for k in G.CORE_IDX), G.CORE.amp, G.CORE.psi, None, None,  # noqa: E731
                                      f"core{j}")
    wfs = [(("clock", 0, "14 genes"), G.WF), (("clock", 0, "core"), G.CORE)]
    wfs += [((kind, j, rd), w if rd == "14 genes" else core_of(w, j))
            for kind, ws in sets.items() for j, w in enumerate(ws[:n]) for rd in ("14 genes", "core")]
    keep = sorted({g for _, w in wfs for g in w.genes})
    G.S.update(pb=pb[G.META + keep].copy(), comps=comps)
    del pb
    tod = donors.tod.astype(float)
    tod.index = tod.index.astype(str)
    print(f"{len(comps)} tissues, {G.S['pb'].donor.nunique()} donors, {len(wfs)} waveforms ({time.time() - t0:.0f}s)",
          flush=True)
    st = G.run(wfs, a.jobs)
    st = [[(np.asarray(j).astype(str), za, zb) for j, za, zb in s] for s in st]
    don = np.unique(np.concatenate([j for j, _, _ in st[0]]))
    bins = np.floor(tod.reindex(don).to_numpy(float) / 2)
    if np.isnan(bins).any():
        raise SystemExit("donor without a time of death in the A336 set")
    nb = pd.Series(bins).value_counts()
    chance = (len(nb) - 1) / (len(don) - 1)
    for lab, s in zip((w[0] for w in wfs[:2]), st[:2]):
        v0 = leverage(s, tod, VARIANTS[0])
        print(f"A336 reproduction {lab[2]}: {v0:.6f} vs {A336[lab[2]]:.6f}", flush=True)
        if not a.check and not np.isclose(v0, A336[lab[2]], atol=5e-7):
            raise SystemExit("A336 part 2 not reproduced: stop before any adjusted value")

    rows = []
    for (kind, j, rd), s in zip((w[0] for w in wfs), st):
        r2_raw, r2_res = radius_r2(s, tod)
        for v in VARIANTS:
            rows.append(dict(kind=kind, set=j, readout=rd, variant=v, leverage=leverage(s, tod, v),
                             r2_bins_raw=r2_raw, r2_bins_resid=r2_res))
    C = pd.DataFrame(rows)
    rng = np.random.default_rng([336, 3, 1])
    perms = within_bin_perms(st[0], tod, rng, 20 if a.check else N_PERM)
    summ = []
    for k, rd in enumerate(("14 genes", "core")):
        for v in VARIANTS:
            c = C[(C.readout == rd) & (C.variant == v)]
            clk = c[c.kind == "clock"].iloc[0]
            null = np.array([leverage(st[k], tod, v, perm=p) for p in perms])
            q = {kk: float(np.quantile(c[c.kind == kk].leverage, 0.95)) for kk in ("matched", "responsive")}
            comp = c[c.kind != "clock"]
            p_perm = float((1 + np.sum(null >= clk.leverage)) / (1 + len(null)))
            summ.append(dict(readout=rd, variant=v, clock=clk.leverage, p_perm=p_perm, null_q95=float(np.quantile(null, 0.95)),
                             matched_q95=q["matched"], responsive_q95=q["responsive"],
                             pct_matched_below=float(100 * np.mean(comp[comp.kind == "matched"].leverage < clk.leverage)),
                             pct_responsive_below=float(100 * np.mean(comp[comp.kind == "responsive"].leverage < clk.leverage)),
                             holds=bool(clk.leverage > q["matched"] and p_perm < 0.05),
                             above_both=bool(clk.leverage > max(q.values())),
                             clock_r2_bins_raw=clk.r2_bins_raw, clock_r2_bins_resid=clk.r2_bins_resid,
                             companion_r2_bins_raw_median=float(comp.r2_bins_raw.median()),
                             companion_r2_bins_resid_median=float(comp.r2_bins_resid.median()),
                             chance_r2=chance, donors=len(don), bins=len(nb), min_bin=int(nb.min())))
        print(f"{rd} done ({time.time() - t0:.0f}s)", flush=True)
    Sm = pd.DataFrame(summ)
    if a.check:
        print("check ok", flush=True)
        return
    OUT.mkdir(parents=True, exist_ok=True)
    C.to_csv(OUT / "c_leverages.csv", index=False)
    Sm.to_csv(OUT / "c_summary.csv", index=False)
    pd.set_option("display.width", 250)
    print(Sm.round(4).to_string(index=False), f"\n({time.time() - t0:.0f}s)", flush=True)


if __name__ == "__main__":
    main()
