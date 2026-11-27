"""A410 (registered in paper/docs/registered_experiments.md, "A410"): cancer-cell coupling contribution against the retention of
the cancer's lineage program, within cancer types.

Pool, P401 cancers and tumour-stroma hosts as A403; the genome-wide malignant pseudobulks of the 7 cancers with a registered
lineage program are kept from the pooled frames. y = per-patient contribution Re(m conj(b h)) / mean |m|^2 (states centred
within stratum; b per cancer); x = lineage score (mean within-cohort z of log2 CPM of the program genes).
  power  null of the pooled slope from within-cohort shuffles of x (the identity pairing is not evaluated); MDE
                                                                               -> output/a410/power.csv
  score  H410.1, matched / responsive companions (pseudo-host and malignant states on the same pseudo waveform),
         per cancer with bootstrap intervals, tumour-fibroblast check, core loop  -> output/a410/{tests,per_cancer,companions}.csv
    python a410.py --stage power|score [--jobs 12] [--genome-jobs 4]          (long: paper/common/run.sh)
"""
from __future__ import annotations

import argparse
import time
import warnings

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

import a403_a405_score as S
from a402_a405 import BASE, CORE, G2M, META, U, cs, log, trim, centred, cancer_hosts, loudness_frame
from clocksync.config import limit_threads
from clocksync.deviation import states as dev_states
from clocksync.genes import canonical

OUT = U.OUT / "a410"  # A411 (--x prolif): output/a411
COUPLED = ("LUAD", "THCA", "KIRC")  # cancers whose cancer cells keep time (A401)
SEED = 410
MALT, FIB = U.MALT, "tumour:fibroblast"
PROGRAMS = {"LUAD": "SFTPC SFTPB SFTPA1 SFTPA2 NAPSA LAMP3 SLC34A2 ABCA3",
            "THCA": "TG TPO TSHR PAX8 SLC26A4 DIO1 DIO2 FOXE1 NKX2-1 SLC5A5",
            "BRCA": "ESR1 PGR GATA3 FOXA1 XBP1 TFF1 AGR2 KRT18",
            "PRAD": "AR KLK3 KLK2 NKX3-1 TMPRSS2 FOLH1 ACP3",
            "COAD": "CDX1 CDX2 VIL1 KRT20 CEACAM5 CA2 SLC26A3 GUCA2A",
            "STAD": "MUC5AC TFF1 TFF2 GKN1 GKN2 PGC CLDN18 LIPF",
            "KIRC": "LRP2 CUBN SLC34A1 SLC22A6 SLC3A1 SLC5A12 ALDOB"}
PROGRAMS = {c: [canonical(g) for g in v.split()] for c, v in PROGRAMS.items()}
COVS = ["log_mal_depth", "log_mal_cells", "n_anchor", "log_depth"]
N_PERM, N_BOOT, N_SET, GATE = 10_000, 2000, 200, 0.20


def zcoh(df):
    coh = np.array([d.split("|")[0] for d in df.index])
    g = df.groupby(coh)
    return (df - g.transform("mean")) / g.transform(lambda x: x.std(ddof=0)).replace(0, np.nan)


def lineage(Mc, cancer):
    genes = [g for g in PROGRAMS[cancer] if g in Mc.columns]
    t = Mc.drop_duplicates("donor").set_index("donor")
    if len(genes) < len(PROGRAMS[cancer]) / 2:
        return pd.Series(dtype=float), genes
    lc = np.log2(t[genes].div(t.depth, axis=0) * 1e6 + 1)
    return zcoh(lc).mean(axis=1), genes


def prolif(Mc):
    """Malignant G2M index (A391 / A403): mean within-cohort z of natural log1p CPM of the Hallmark G2M genes."""
    t = Mc.drop_duplicates("donor").set_index("donor")
    genes = [g for g in G2M if g in t.columns]
    if len(genes) < 10:
        return pd.Series(dtype=float), genes
    return zcoh(np.log1p(t[genes].div(t.depth, axis=0) * 1e6)).mean(axis=1), genes


def contrib(rows, host, stratum, wf):
    """Per-patient contribution of the compartment `rows` (one per donor) on the host state (both centred within stratum)."""
    s = dev_states(rows[META + list(wf.genes)], wf)
    m = pd.Series(s.to_numpy(), index=s.index.get_level_values("donor"))
    hz = pd.Series((host.x + 1j * host.y).to_numpy(), index=host.index.astype(str))
    ids = m.index.intersection(hz.index)
    st = stratum.reindex(ids)
    m, h = centred(m[ids].to_numpy(), st), centred(hz[ids].to_numpy(), st)
    b = (m * np.conj(h)).sum() / (np.abs(h) ** 2).sum()
    return pd.Series(np.real(m * np.conj(b * h)) / np.mean(np.abs(m) ** 2), index=ids)


def units(df, ycol, xcol, covs):
    out = []
    for s, g in df.dropna(subset=[ycol, xcol] + covs).groupby("cohort"):
        if len(g) < 4:
            continue
        y, x = centred(g[ycol].to_numpy(float), g.stratum), centred(g[xcol].to_numpy(float), g.stratum)
        X = np.column_stack([centred(g[c].to_numpy(float), g.stratum) for c in covs])
        if len(g) > X.shape[1] + 2:
            y, x = S.residual(y, X), S.residual(x, X)
        if y.std() < 1e-12 or x.std() < 1e-12:
            continue
        out.append(dict(cohort=s, n=len(g), y=y / y.std(), x=x / x.std()))
    return out


def pooled_slope(un):
    w = np.array([u["n"] for u in un], float)
    return float(sum(u["n"] * np.mean(u["y"] * u["x"]) for u in un) / w.sum())


def perm_null(un, n_perm, seed):
    rng = np.random.default_rng(seed)
    w = np.array([u["n"] for u in un], float)
    nul = np.zeros(n_perm)
    for u in un:
        P = np.argsort(rng.random((n_perm, u["n"])), axis=1)
        nul += u["n"] * (u["y"][None, :] * u["x"][P]).mean(1)
    return nul / w.sum()


def build(args):
    frames, donors, f388 = U.pool_sources(args, OUT)
    from joblib.externals.loky import get_reusable_executor
    get_reusable_executor().shutdown(wait=True)
    p401 = pd.read_csv(U.OUT / "power" / "p401.csv")
    cancers = [c for c in p401[(p401.readout.astype(str) == "14") & p401.p401_pass].cancer if c in set(donors.cancer) and c != "unspecified"]
    donors = donors[donors.cancer.isin(cancers)].reset_index(drop=True)
    hosts = cancer_hosts(frames, donors, cancers)
    mal = {}
    for c in (PROGRAMS if args.x == "lineage" else list(hosts)):
        if c not in hosts:
            continue
        ids = set(donors[donors.cancer == c].donor)
        parts = [f[f.donor.isin(ids) & (f.compartment == MALT)] for f in frames.values()]
        parts = [p for p in parts if len(p)]
        Mc = pd.concat(parts, join="inner", ignore_index=True)
        keep = sorted((set(hosts[c]["Gc"].columns) | set(PROGRAMS.get(c, [])) | set(G2M)) & set(Mc.columns))
        mal[c] = Mc[keep].reset_index(drop=True)
    frames.clear()
    trim()
    loud = loudness_frame(hosts)
    assert len(hosts) == 12 and len(loud) == 574, "pool differs from A403"
    cov = S.covariates(hosts)
    rows, info = [], []
    for c, Mc in mal.items():
        h = hosts[c]
        x, genes = lineage(Mc, c) if args.x == "lineage" else prolif(Mc)
        y = contrib(Mc.drop_duplicates("donor"), h["host"], h["stratum"], BASE)
        Gc = h["Gc"]
        yf = contrib(Gc[Gc.compartment == FIB].drop_duplicates("donor"), h["host"], h["stratum"], BASE) if (Gc.compartment == FIB).any() else pd.Series(dtype=float)
        m1 = Mc.drop_duplicates("donor").set_index("donor")
        d = pd.DataFrame({"y": y, "x": x, "y_fib": yf}).reindex(y.index)
        d["log_mal_depth"] = np.log(m1.depth.reindex(d.index).astype(float))
        d["log_mal_cells"] = np.log(m1.n_cells.reindex(d.index).astype(float))
        d["cancer"] = c
        rows.append(d.rename_axis("donor").reset_index())
        info.append(dict(cancer=c, program_genes=" ".join(genes), n_genes=len(genes), n_registered=len(PROGRAMS.get(c, G2M)), patients=len(d)))
    df = pd.concat(rows, ignore_index=True).merge(cov, on="donor", how="left")
    df["cohort"] = df.donor.str.split("|").str[0]
    df["stratum"] = df.donor.map(loud.set_index("donor").stratum)
    pd.DataFrame(info).to_csv(OUT / "programs.csv", index=False)
    log(f"table: {len(df)} patients; x in {df.x.notna().sum()}, y in {df.y.notna().sum()}")
    return df, hosts, mal


def stage_power(args):
    df, _, _ = build(args)
    un = units(df, "y", "x", COVS)
    nul = perm_null(un, 2000, [SEED, 9])
    T = pd.DataFrame([dict(patients=int(sum(u["n"] for u in un)), cohorts=len(un), null_mean=float(nul.mean()), null_sd=float(nul.std()),
                           mde=2.49 * float(nul.std()), gate=GATE, scored=bool(2.49 * nul.std() <= GATE))])
    T.to_csv(OUT / "power.csv", index=False)
    log("power\n" + T.round(4).to_string(index=False))


def pseudo_one(c, Gc, Mc, classes, stratum, w):
    limit_threads(1)
    warnings.filterwarnings("ignore")
    try:
        genes = sorted(set(w.genes) | set(BASE.genes))
        h = cs.Dataset(Gc[META + genes], anchors=classes, waveform=w, mixed_strata=True).host()
        return contrib(Mc.drop_duplicates("donor"), h, stratum, w)
    except Exception:  # noqa: BLE001  (pseudo waveform without a fit: no value)
        return pd.Series(dtype=float)


def stage_score(args):
    pw = pd.read_csv(OUT / "power.csv")
    scored = bool(pw.scored.iloc[0])
    df, hosts, mal = build(args)
    un = units(df, "y", "x", COVS)
    obs = pooled_slope(un)
    nul = perm_null(un, N_PERM, [SEED, 1])
    lower = args.x == "prolif"  # A411 predicts a negative slope
    p = float((np.sum((nul <= obs + 1e-12) if lower else (nul >= obs - 1e-12)) + 1) / (N_PERM + 1))
    # companions
    from clocksync.paired import pseudo_waveforms
    comp = {"matched": [[] for _ in range(N_SET)], "responsive": [[] for _ in range(N_SET)]}
    for c, Mc in mal.items():
        h = hosts[c]
        Gc = h["Gc"]
        genes = [g for g in Gc.columns if g not in META]
        frac = Gc[genes].div(Gc.depth, axis=0)
        lm = np.log(frac.mean() + 1e-12)[(Gc[genes].mean() >= 1) | frac.columns.isin(list(BASE.genes))]
        lm = lm[lm.index.isin(Mc.columns)]
        for code, fam in enumerate(("matched", "responsive")):
            ws = pseudo_waveforms(lm, N_SET, np.random.default_rng([SEED, code]), BASE, fam)
            res = Parallel(n_jobs=args.jobs)(delayed(pseudo_one)(c, Gc, Mc, h["classes"], h["stratum"], w) for w in ws)
            for j, r in enumerate(res):
                comp[fam][j].append(r)
        log(f"companions {c}")
    crow = []
    for fam, sets in comp.items():
        for j, parts in enumerate(sets):
            yj = pd.concat(parts) if parts else pd.Series(dtype=float)
            d2 = df.drop(columns="y").merge(yj.rename("y"), left_on="donor", right_index=True)
            u2 = units(d2, "y", "x", COVS)
            crow.append(dict(family=fam, set=j, slope=pooled_slope(u2) if u2 else np.nan))
    C = pd.DataFrame(crow)
    C.to_csv(OUT / "companions.csv", index=False)
    row = dict(hypothesis=(f"H{SEED}.1" if scored else "descriptive (power)"), patients=int(sum(u["n"] for u in un)), cohorts=len(un),
               slope=obs, null_sd=float(nul.std()), p_one_sided=p)
    for fam in ("matched", "responsive"):
        v = C[C.family == fam].slope.dropna().to_numpy()
        row.update({f"{fam}_n": len(v), f"{fam}_median": float(np.median(v)), f"{fam}_q05": float(np.quantile(v, .05)),
                    f"{fam}_q95": float(np.quantile(v, .95)), f"{fam}_pct": float(100 * np.mean(v < obs))})
    row["pass"] = bool(scored and p < 0.05 and (obs < row["matched_q05"] if lower else obs > row["matched_q95"]))
    rows = [row]
    # tumour-fibroblast check and core loop (descriptive)
    for lab, ycol in (("tumour fibroblast contribution (check)", "y_fib"),):
        u3 = units(df, ycol, "x", COVS)
        n3 = perm_null(u3, N_PERM, [SEED, 2])
        o3 = pooled_slope(u3)
        rows.append(dict(hypothesis=lab, patients=int(sum(u["n"] for u in u3)), cohorts=len(u3), slope=o3, null_sd=float(n3.std()),
                         p_one_sided=float((np.sum((n3 <= o3 + 1e-12) if lower else (n3 >= o3 - 1e-12)) + 1) / (N_PERM + 1))))
    if lower:
        u5 = units(df[df.cancer.isin(COUPLED)], "y", "x", COVS)
        n5 = perm_null(u5, N_PERM, [SEED, 7])
        o5 = pooled_slope(u5)
        rows.append(dict(hypothesis="coupled cancers LUAD, THCA, KIRC (descriptive)", patients=int(sum(u["n"] for u in u5)), cohorts=len(u5),
                         slope=o5, null_sd=float(n5.std()), p_one_sided=float((np.sum(n5 <= o5 + 1e-12) + 1) / (N_PERM + 1))))
    core_parts = []
    for c, Mc in mal.items():
        h = hosts[c]
        hc = cs.Dataset(h["tab"][META + list(CORE.genes)], anchors=h["classes"], waveform=CORE, mixed_strata=True).host()
        core_parts.append(contrib(Mc.drop_duplicates("donor"), hc, h["stratum"], CORE))
    dcore = df.drop(columns="y").merge(pd.concat(core_parts).rename("y"), left_on="donor", right_index=True)
    u4 = units(dcore, "y", "x", COVS)
    n4 = perm_null(u4, N_PERM, [SEED, 3])
    o4 = pooled_slope(u4)
    rows.append(dict(hypothesis="core loop (descriptive)", patients=int(sum(u["n"] for u in u4)), cohorts=len(u4), slope=o4,
                     null_sd=float(n4.std()), p_one_sided=float((np.sum((n4 <= o4 + 1e-12) if lower else (n4 >= o4 - 1e-12)) + 1) / (N_PERM + 1))))
    T = pd.DataFrame(rows)
    T.to_csv(OUT / "tests.csv", index=False)
    pd.set_option("display.width", 250, "display.max_columns", 30)
    log("tests\n" + T.round(4).to_string(index=False))
    # per cancer, bootstrap within cohort
    pc = []
    for c, d in df.groupby("cancer"):
        uc = units(d, "y", "x", COVS)
        if not uc:
            continue
        rng = np.random.default_rng([SEED, 5, sum(map(ord, c))])
        bs = []
        for _ in range(N_BOOT):
            idx = np.concatenate([rng.choice(g.index.to_numpy(), len(g)) for _, g in d.groupby("cohort")])
            ub = units(d.loc[idx].reset_index(drop=True), "y", "x", COVS)
            bs.append(pooled_slope(ub) if ub else np.nan)
        nc = perm_null(uc, N_PERM, [SEED, 6, sum(map(ord, c))])
        oc = pooled_slope(uc)
        pc.append(dict(cancer=c, patients=int(sum(u["n"] for u in uc)), cohorts=len(uc), slope=oc, lo90=float(np.nanquantile(bs, .05)),
                       hi90=float(np.nanquantile(bs, .95)), p_two_sided=float((np.sum(np.abs(nc) >= abs(oc) - 1e-12) + 1) / (N_PERM + 1)),
                       mean_contribution=float(d.y.mean())))
    P = pd.DataFrame(pc)
    P.to_csv(OUT / "per_cancer.csv", index=False)
    log("per cancer\n" + P.round(4).to_string(index=False))


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stage", choices=("power", "score"), required=True)
    ap.add_argument("--jobs", type=int, default=12)
    ap.add_argument("--genome-jobs", type=int, default=4)
    ap.add_argument("--x", choices=("lineage", "prolif"), default="lineage", help="lineage (A410) or proliferation (A411)")
    args = ap.parse_args()
    args.smoke = False
    global OUT, SEED
    if args.x == "prolif":
        OUT, SEED = U.OUT / "a411", 411
    limit_threads(1)
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    (stage_power if args.stage == "power" else stage_score)(args)
    log(f"done {(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    main()
