"""A402b (registered in paper/docs/registered_experiments.md, "A402b"): tumour T/NK and myeloid shares against the phase of
the tumour-stroma clock (per readout, r units; replaces A402's omnibus rule).

Pool, P401 cancers, hosts and slim tables as A402-A405 (`a402_a405.py`); the slim tables are rebuilt in their own cache
(`output/cache_a402b/`) with the T/NK subtype markers added. Statistic r = |sum y conj z| / sqrt(sum y^2 sum |z|^2) per cohort
(A229 / A246), y and z centred within cohort x chemistry, pooled with patient weights; null: the pairing shuffled within
cohort, one ordering per cohort per iteration shared by all readouts. Effect = r - null mean.

  power  nulls only (no identity pairing evaluated): MDE = 2.49 x pooled null SD of the two shares on the primary set
         (cohorts outside the A227 list) and on all patients; gates 0.076 (T/NK) and 0.053 (myeloid).
                                                                          -> output/a402b/power.csv
  score  H402b.1-2 (scored readouts: permutation p, Holm, matched / responsive pseudo-host q95), D1 per cancer and phase
         agreement, D2 glucocorticoid score, D3 T/NK subtype markers, D4 all patients.
                                                                          -> output/a402b/{tests,d1_cancer,d1_agreement,d2,d3}.csv

    python a402b.py --stage power|score [--jobs 16] [--genome-jobs 4]      (long: paper/common/run.sh)
"""
from __future__ import annotations

import os

SUBTYPES = {"CD8": ["CD8A", "CD8B"], "CD4": ["CD4", "CD40LG"], "Treg": ["FOXP3", "IL2RA"], "NK": ["KLRF1", "NCR1", "NCAM1"]}
os.environ["A402_CACHE"] = "cache_a402b"  # before importing a402_a405 (module constants; the env reaches the loky workers)
os.environ["A402_EXTRA_GENES"] = ",".join(g for v in SUBTYPES.values() for g in v)

import argparse  # noqa: E402
import time  # noqa: E402
import warnings  # noqa: E402

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import a402_a405 as P  # noqa: E402
from a402_a405 import BASE, META, MDE_K, U, cs, log, trim, centred, cancer_hosts, loudness_frame, load_slim  # noqa: E402
from clocksync._core.controls import GC_TARGETS  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402

OUT = U.OUT / "a402b"
N_PERM, N_BOOT, N_SET = 10_000, 2000, 200
GATE = {"tnk_share": 0.076, "myeloid_share": 0.053}
HYP = {"tnk_share": "H402b.1", "myeloid_share": "H402b.2"}
COMP = {"tnk_share": "tumour:T_NK", "myeloid_share": "tumour:myeloid"}
MIN_N, MIN_ALL = 4, 100
SEED = 4022


# ---- table ---------------------------------------------------------------------------------------------------------------

def zcoh(df, cols):
    """Within-cohort z of each column (cohort = donor prefix)."""
    coh = np.array([d.split("|")[0] for d in df.index])
    g = df.groupby(coh)
    return (df - g.transform("mean")) / g.transform(lambda x: x.std(ddof=0)).replace(0, np.nan)


def readouts(slim, z_of):
    """Per donor: log shares (A246 P8 construction, also for T/NK), T/NK subtype marker scores."""
    rows = []
    for s, G in slim.items():
        G = G[G.donor.isin(set(z_of.index))]
        if G.empty:
            continue
        N = G.groupby("donor").n_cells.sum()
        r = pd.DataFrame(index=pd.Index(sorted(set(G.donor)), name="donor"))
        for nm, comp in COMP.items():
            n = G[G.compartment == comp].drop_duplicates("donor").set_index("donor").n_cells
            ok = n.index[(N.reindex(n.index) >= MIN_ALL).to_numpy()]
            r[nm] = np.log2((n.reindex(ok) + .5) / (N.reindex(ok) + 1))
        t = G[G.compartment == COMP["tnk_share"]].drop_duplicates("donor").set_index("donor")
        if len(t) >= MIN_N:
            for nm, genes in SUBTYPES.items():
                gg = [g for g in genes if g in t.columns]
                if len(gg) >= len(genes) / 2:
                    lc = np.log2(t[gg].div(t.depth, axis=0) * 1e6 + 1)
                    r[nm] = zcoh(lc, gg).mean(axis=1)
        r["cohort"] = s
        rows.append(r.reset_index())
    return pd.concat(rows, ignore_index=True)


def gc_score(hosts):
    """Per donor: glucocorticoid-target score in the summed tumour stroma (mean within-cohort z of log2 CPM)."""
    rows = []
    for h in hosts.values():
        Gc = h["Gc"]
        genes = [g for g in GC_TARGETS if g in Gc.columns]
        t = Gc.groupby("donor")[genes + ["depth"]].sum()
        z = zcoh(np.log2(t[genes].div(t.depth, axis=0) * 1e6 + 1), genes)
        z = z.loc[:, z.notna().mean() > 0.5]
        rows.append(pd.DataFrame({"donor": z.index, "gc": z.mean(axis=1).to_numpy()}))
    return pd.concat(rows, ignore_index=True)


# ---- statistic -----------------------------------------------------------------------------------------------------------

def cohort_units(df, col, zcol="z"):
    """{cohort: unit} with y, z centred within stratum over the patients with both (>= 4, spread)."""
    out = {}
    for s, g in df.dropna(subset=[col]).groupby("cohort"):
        g = g[np.isfinite(g[zcol].to_numpy())]
        if len(g) < MIN_N:
            continue
        y = centred(g[col].to_numpy(float), g.stratum)
        z = centred(g[zcol].to_numpy(), g.stratum)
        den = float(np.sqrt(np.sum(y ** 2) * np.sum(abs(z) ** 2)))
        if den > 0:
            out[s] = dict(n=len(g), donor=g.donor.to_numpy(), cancer=g.cancer.iloc[0], y=y, z=z, den=den)
    return out


def cvals(u, y=None, z=None):
    y, z = (u["y"] if y is None else y), (u["z"] if z is None else z)
    return np.sum(y * np.conj(z)) / u["den"]


def perm_c(u, order):
    """Complex cross-correlations under the shared cohort ordering `order` (n_perm x patients of the cohort, keyed by donor)."""
    keys = order.loc[u["donor"]].to_numpy().T  # n_perm x patients
    perms = np.argsort(keys, axis=1)
    return (u["y"][perms] * np.conj(u["z"])[None, :]).sum(1) / u["den"]


def orders(df, n_perm):
    """One random key per (iteration, patient) per cohort, shared by every readout."""
    out = {}
    for s, g in df.groupby("cohort"):
        rng = np.random.default_rng([SEED, sum(map(ord, s))])
        out[s] = pd.DataFrame(rng.random((len(g), n_perm)).astype(np.float32), index=g.donor.to_numpy())
    return out


def pooled(units, ords):
    """Pooled observed r, complex c and the permutation null (|.| per cohort pooled, and complex pooled)."""
    w = np.array([u["n"] for u in units.values()], float)
    c = np.array([cvals(u) for u in units.values()])
    nc = np.array([perm_c(u, ords[s]) for s, u in units.items()])
    return dict(n=int(w.sum()), cohorts=len(units), r=float(np.sum(w * abs(c)) / w.sum()), c=complex(np.sum(w * c) / w.sum()),
                null=(w[:, None] * abs(nc)).sum(0) / w.sum(), nullc=(w[:, None] * nc).sum(0) / w.sum(), w=w)


def peak_h(c):
    return float((-np.angle(c)) % (2 * np.pi) * 24 / (2 * np.pi))


def companions(units, Z, w):
    """Pooled r for each pseudo-host set: Z = {fam: [Series donor -> z]}."""
    out = {}
    for fam, sets in Z.items():
        vals = []
        for zs in sets:
            num, den = 0.0, 0.0
            for u, wi in zip(units.values(), w):
                zz = zs.reindex(u["donor"]).to_numpy()
                if not np.all(np.isfinite(zz)):
                    continue
                zc = centred(zz, u["strat"])
                d = float(np.sqrt(np.sum(u["y"] ** 2) * np.sum(abs(zc) ** 2)))
                if d > 0:
                    num += wi * abs(np.sum(u["y"] * np.conj(zc)) / d)
                    den += wi
            vals.append(num / den if den else np.nan)
        out[fam] = np.array(vals)
    return out


def pseudo_states(tab, classes, w, stratum):
    limit_threads(1)
    warnings.filterwarnings("ignore")
    try:
        h = cs.Dataset(tab, anchors=classes, waveform=w, mixed_strata=True).host()
        return pd.Series((h.x + 1j * h.y).to_numpy(), index=h.index)
    except Exception:  # noqa: BLE001  (pseudo-host not fitted: the set is undefined for this cancer)
        return pd.Series(dtype=complex)


def companion_states(hosts, args):
    """{family: [Series donor -> pseudo-host z]} on each cancer's pooled stroma; the A403 sets (seeds [401, 0] / [401, 1])."""
    from clocksync.paired import pseudo_waveforms
    out = {"matched": [[] for _ in range(N_SET)], "responsive": [[] for _ in range(N_SET)]}
    for c, h in hosts.items():
        Gc = h["Gc"]
        genes = [g for g in Gc.columns if g not in META]
        frac = Gc[genes].div(Gc.depth, axis=0)
        lm = np.log(frac.mean() + 1e-12)[(Gc[genes].mean() >= 1) | frac.columns.isin(list(BASE.genes))]
        for code, fam in enumerate(("matched", "responsive")):
            ws = pseudo_waveforms(lm, N_SET, np.random.default_rng([U.SEED, code]), BASE, fam)
            res = Parallel(n_jobs=args.jobs)(delayed(pseudo_states)(Gc[META + sorted(set(w.genes) | set(BASE.genes))], h["classes"], w,
                                                                    h["stratum"]) for w in ws)
            for j, r in enumerate(res):
                out[fam][j].append(r)
        log(f"companion states {c}")
    return {f: [pd.concat(v) if v else pd.Series(dtype=complex) for v in sets] for f, sets in out.items()}


# ---- stages --------------------------------------------------------------------------------------------------------------

def build(args):
    frames, donors, f388 = U.pool_sources(args, OUT)
    from joblib.externals.loky import get_reusable_executor
    get_reusable_executor().shutdown(wait=True)
    p401 = pd.read_csv(U.OUT / "power" / "p401.csv")
    cancers = [c for c in p401[(p401.readout.astype(str) == "14") & p401.p401_pass].cancer if c in set(donors.cancer) and c != "unspecified"]
    donors = donors[donors.cancer.isin(cancers)].reset_index(drop=True)
    hosts = cancer_hosts(frames, donors, cancers)
    sources, f388 = list(frames), set(f388)
    frames.clear()
    trim()
    loud = loudness_frame(hosts)
    assert len(hosts) == 12 and len(loud) == 574, "pool differs from A402-A405"
    slim, errs = load_slim(sources, donors, f388, args)
    assert len(slim) == 51 and errs.empty, errs
    ro = readouts(slim, loud.set_index("donor").z)
    df = loud.merge(ro.drop(columns="cohort"), on="donor", how="inner")
    df["cohort"] = df.donor.str.split("|").str[0]
    a227 = set(P.M.cohort_list())
    df["new"] = ~df.cohort.isin(a227)
    log(f"table: {len(df)} patients ({int(df.new.sum())} new); T/NK share {df.tnk_share.notna().sum()}, myeloid {df.myeloid_share.notna().sum()}")
    return df, hosts


def stage_power(args):
    df, _ = build(args)
    ords = orders(df, 2000)
    rows = []
    for lab, sub in (("primary (new relative to A246)", df[df.new]), ("all patients (descriptive)", df)):
        for col in COMP:
            un = cohort_units(sub, col)
            w = np.array([u["n"] for u in un.values()], float)
            nul = (w[:, None] * np.array([abs(perm_c(u, ords[s])) for s, u in un.items()])).sum(0) / w.sum()
            rows.append(dict(set=lab, readout=col, patients=int(w.sum()), cohorts=len(un), null_mean=float(nul.mean()), null_sd=float(nul.std()),
                             mde=MDE_K * float(nul.std()), gate=GATE[col], scored=bool(MDE_K * nul.std() <= GATE[col]) if lab.startswith("primary") else False))
    T = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    T.to_csv(OUT / "power.csv", index=False)
    log("power\n" + T.round(4).to_string(index=False))


def test_row(units, ords, Z, label, col, extra=None):
    pl = pooled(units, ords)
    comp = companions(units, Z, pl["w"]) if Z else {}
    p = float((np.sum(pl["null"] >= pl["r"] - 1e-12) + 1) / (len(pl["null"]) + 1))
    r = dict(set=label, readout=col, patients=pl["n"], cohorts=pl["cohorts"], r=pl["r"], null_mean=float(pl["null"].mean()),
             null_sd=float(pl["null"].std()), excess=pl["r"] - float(pl["null"].mean()), p=p, implied_peak_h=peak_h(pl["c"]))
    for fam, v in comp.items():
        v = v[np.isfinite(v)]
        r.update({f"{fam}_n": len(v), f"{fam}_median": float(np.median(v)), f"{fam}_q95": float(np.quantile(v, .95))})
    return r | (extra or {}), pl


def with_strata(units, df):
    st = df.set_index("donor").stratum
    for u in units.values():
        u["strat"] = pd.Series(st.reindex(u["donor"]).to_numpy())
    return units


def boot_peak(units, rng):
    """90 % arc of the implied peak from within-cohort patient bootstraps (re-centred)."""
    pk = []
    for _ in range(N_BOOT):
        num, wsum = 0j, 0.0
        for u in units.values():
            i = rng.integers(0, u["n"], u["n"])
            y, z = u["y"][i] - u["y"][i].mean(), u["z"][i] - u["z"][i].mean()
            d = np.sqrt(np.sum(y ** 2) * np.sum(abs(z) ** 2))
            if d > 0:
                num += u["n"] * np.sum(y * np.conj(z)) / d
                wsum += u["n"]
        pk.append(peak_h(num / wsum) if wsum else np.nan)
    pk = np.array(pk)
    centre = np.angle(np.exp(1j * pk * 2 * np.pi / 24).mean()) * 24 / (2 * np.pi)
    dev = (pk - centre + 12) % 24 - 12
    return float((centre + np.quantile(dev, .05)) % 24), float((centre + np.quantile(dev, .95)) % 24)


def stage_score(args):
    pw = pd.read_csv(OUT / "power.csv")
    scored = set(pw[pw.scored].readout)
    df, hosts = build(args)
    df = df.merge(gc_score(hosts), on="donor", how="left")
    ords = orders(df, N_PERM)
    Z = companion_states(hosts, args)
    from joblib.externals.loky import get_reusable_executor
    get_reusable_executor().shutdown(wait=True)
    new, rows, d1, agree = df[df.new], [], [], []
    # H402b.1-2 (and descriptive when not scored), D4 all patients
    for lab, sub in (("primary", new), ("all", df)):
        for col in COMP:
            un = with_strata(cohort_units(sub, col), sub)
            r, pl = test_row(un, ords, Z, lab, col, dict(hypothesis=HYP[col] if (lab == "primary" and col in scored) else "descriptive"))
            r["peak_lo90"], r["peak_hi90"] = boot_peak(un, np.random.default_rng([SEED, 1, sum(map(ord, col + lab))]))
            rows.append(r)
            # D1: per cancer and phase agreement across cancers
            cancers = sorted({u["cancer"] for u in un.values()})
            Cc, Nc, Wc = [], [], []
            for c in cancers:
                uc = {s: u for s, u in un.items() if u["cancer"] == c}
                rc, plc = test_row(uc, ords, Z, lab, col, dict(cancer=c))
                rc["peak_lo90"], rc["peak_hi90"] = boot_peak(uc, np.random.default_rng([SEED, 2, sum(map(ord, col + lab + c))]))
                d1.append(rc)
                Cc.append(plc["c"] * plc["n"])
                Nc.append(plc["nullc"] * plc["n"])
            Cc, Nc = np.array(Cc), np.array(Nc)
            A_obs = float(abs(Cc.sum()) / abs(Cc).sum())
            A_null = abs(Nc.sum(0)) / abs(Nc).sum(0)
            agree.append(dict(set=lab, readout=col, cancers=len(cancers), agreement=A_obs, null_mean=float(A_null.mean()),
                              null_q95=float(np.quantile(A_null, .95)), p=float((np.sum(A_null >= A_obs - 1e-12) + 1) / (len(A_null) + 1))))
            log(f"{lab} {col} done")
    T = pd.DataFrame(rows)
    prim = T[T.hypothesis.str.startswith("H")]
    if len(prim):
        order = prim.p.sort_values().index
        holm, run = {}, 0.0
        for k, i in enumerate(order):
            run = max(run, min(1.0, (len(order) - k) * T.loc[i, "p"]))
            holm[i] = run
        T["holm"] = pd.Series(holm)
        T["pass"] = ((T.holm < 0.05) & (T.r > T.matched_q95)).where(T.hypothesis.str.startswith("H"))
    T.to_csv(OUT / "tests.csv", index=False)
    pd.DataFrame(d1).to_csv(OUT / "d1_cancer.csv", index=False)
    pd.DataFrame(agree).to_csv(OUT / "d1_agreement.csv", index=False)
    log("tests\n" + T.round(4).to_string(index=False))
    log("agreement\n" + pd.DataFrame(agree).round(4).to_string(index=False))
    # D2 glucocorticoids
    d2 = []
    for lab, sub in (("primary", new), ("all", df)):
        un = with_strata(cohort_units(sub, "gc"), sub)
        r, _ = test_row(un, ords, Z, lab, "gc_stroma", dict(kind="gc score on host state"))
        d2.append(r)
        for col in COMP:
            s2 = sub.dropna(subset=[col, "gc"]).copy()
            res, cors = [], []
            for s, g in s2.groupby("cohort"):
                a, b = centred(g[col].to_numpy(float), g.stratum), centred(g.gc.to_numpy(float), g.stratum)
                if len(g) >= MIN_N and b.std() > 0 and a.std() > 0:
                    beta = np.sum(a * b) / np.sum(b * b)
                    res.append(pd.Series(a - beta * b, index=g.index))
                    cors.append((len(g), float(np.corrcoef(a, b)[0, 1])))
            s2[col + "_gcres"] = pd.concat(res) if res else np.nan
            un = with_strata(cohort_units(s2, col + "_gcres"), s2)
            r, _ = test_row(un, ords, Z, lab, col, dict(kind="share residualised on gc"))
            n_, c_ = np.array(cors).T
            r["r_share_gc_within_cohort"] = float((n_ * c_).sum() / n_.sum())
            d2.append(r)
    pd.DataFrame(d2).to_csv(OUT / "d2.csv", index=False)
    log("D2\n" + pd.DataFrame(d2).round(4).to_string(index=False))
    # D3 T/NK subtype markers
    d3 = []
    for lab, sub in (("primary", new), ("all", df)):
        for nm in SUBTYPES:
            if nm not in sub or sub[nm].notna().sum() < 20:
                continue
            un = with_strata(cohort_units(sub, nm), sub)
            r, _ = test_row(un, ords, Z, lab, nm)
            d3.append(r)
    pd.DataFrame(d3).to_csv(OUT / "d3.csv", index=False)
    log("D3\n" + pd.DataFrame(d3).round(4).to_string(index=False))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stage", choices=("power", "score"), required=True)
    ap.add_argument("--jobs", type=int, default=16)
    ap.add_argument("--genome-jobs", type=int, default=4)
    args = ap.parse_args()
    args.smoke = False
    limit_threads(1)
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    (stage_power if args.stage == "power" else stage_score)(args)
    log(f"done {(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    main()
