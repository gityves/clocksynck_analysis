"""A419 (POST HOC, DESCRIPTIVE; registered in paper/docs/registered_experiments.md, "A419"): the endothelial colder-louder
slope of A416: matched / narrow responsive / broad responsive companions with placebo-calibrated p (S1), effect size as an
amplitude ratio between immune-hot tertiles (S2), IFN-gamma and adhesion readouts in the endothelium (S3).
                                                      -> output/a419/{s1,s1_sets,s2,s3}.csv
    python a419.py [--jobs 4] [--genome-jobs 2]                          (long: paper/common/run.sh)
"""
from __future__ import annotations

import argparse
import time

import numpy as np
import pandas as pd

import a403_a405_score as S
import a403b as B
import a416 as K
from a402_a405 import BASE, META, U, log, trim, centred, cancer_hosts, loudness_frame, load_slim, readouts
from clocksync import deviation, pairwise
from clocksync.config import limit_threads
from clocksync.genes import canonical

OUT = U.OUT / "a419"
ENDO = "tumour:endothelial"
IFN = [canonical(g) for g in ("STAT1", "IRF1", "CXCL9", "CXCL10")]
ADH = [canonical(g) for g in ("VCAM1", "ICAM1", "SELE")]
ISG = [canonical(g) for g in ("STAT1", "IRF1", "GBP1", "GBP4", "PSMB9", "TAP1")]  # A421: endothelium-intrinsic
CONTAM = [canonical(g) for g in ("PTPRC", "CD3E", "LYZ")]
N_SET, N_BOOT = 200, 2000
FAMS = {"matched": ("matched", 0), "responsive": ("responsive", 1), "responsive_broad": ("responsive_broad", 2), "placebo": ("matched", 7)}


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
    slim, errs = load_slim(sources, donors, f388, args)
    assert len(slim) == 51 and errs.empty
    ro = readouts(slim, loud)
    df = loud.merge(ro.drop(columns="cohort"), on="donor").merge(S.covariates(hosts), on="donor", how="left")
    df["cohort"] = df.donor.str.split("|").str[0]
    df["loud_pf"], _ = K.C3.phase_free(df.loud, np.angle(df.z.to_numpy()), df.stratum.to_numpy(), df.cancer.to_numpy())
    return df, hosts


def tertile_ratio(d, col):
    """Hot minus cold (within-cohort immune-hot tertiles) mean of `col`, cohort-weighted; exp = amplitude ratio."""
    num, w = 0.0, 0
    for _, g in d.dropna(subset=["immune_hot", col]).groupby("cohort"):
        if len(g) < 6:
            continue
        q = g.immune_hot.rank(pct=True)
        hot, cold = g[col][q > 2 / 3], g[col][q <= 1 / 3]
        if len(hot) and len(cold):
            num += len(g) * (hot.mean() - cold.mean())
            w += len(g)
    return num / w if w else np.nan


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--genome-jobs", type=int, default=2)
    ap.add_argument("--variant", choices=("main", "intrinsic"), default="main", help="intrinsic: A421 (S3 only, intrinsic ISG score)")
    args = ap.parse_args()
    global OUT
    if args.variant == "intrinsic":
        OUT = U.OUT / "a421"
    args.smoke = False
    limit_threads(1)
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    df, hosts = build(args)
    from clocksync.paired import pseudo_waveforms
    obs, sets, sc, rel = [], {f: [[] for _ in range(N_SET)] for f in FAMS}, [], []
    for c, h in hosts.items():
        Gc = h["Gc"]
        rows = Gc[Gc.compartment == ENDO].drop_duplicates("donor")
        if rows.donor.nunique() < 6:
            continue
        obs.append(K.class_loud(rows[META + list(BASE.genes)], BASE, c))
        genes = [g for g in Gc.columns if g not in META]
        frac = Gc[genes].div(Gc.depth, axis=0)
        lm = np.log(frac.mean() + 1e-12)[(Gc[genes].mean() >= 1) | frac.columns.isin(list(BASE.genes))]
        for f, (kind, code) in (FAMS.items() if args.variant == "main" else []):
            for j, w in enumerate(pseudo_waveforms(lm, N_SET, np.random.default_rng([U.SEED, code]), BASE, kind)):
                try:
                    sets[f][j].append(K.class_loud(rows[META + list(w.genes)], w, c))
                except Exception:  # noqa: BLE001
                    pass
        t = rows.set_index("donor")
        sc.append(pd.DataFrame({"ifn": K.score(t, IFN), "adh": K.score(t, ADH), "isg": K.score(t, ISG),
                                "contam": K.score(t, CONTAM)}).rename_axis("donor").reset_index())
        # split-half reliability of the endothelial state (S2 disattenuation)
        st = rows.stratum.astype(str)
        for ga, gb in (deviation.split_genes(BASE, 5, np.random.default_rng([4191, sum(map(ord, c))])) if args.variant == "main" else []):
            za = centred(pairwise.states(rows, BASE.subset(list(ga))), st)
            zb = centred(pairwise.states(rows, BASE.subset(list(gb))), st)
            rho = abs(np.sum(za * np.conj(zb))) / np.sqrt(np.sum(abs(za) ** 2) * np.sum(abs(zb) ** 2))
            rel.append((len(rows), 2 * rho / (1 + rho)))
        log(f"endothelium {c}")
    df = df.merge(pd.concat(obs).rename("loud_endo"), left_on="donor", right_index=True, how="left")
    df = df.merge(pd.concat(sc, ignore_index=True), on="donor", how="left")
    if args.variant == "intrinsic":
        out = []
        for a, b in (("isg", "immune_hot"), ("isg", "loud_endo"), ("isg", "ifn"), ("contam", "immune_hot"), ("contam", "loud_endo")):
            rr = []
            for _, g in df.dropna(subset=[a, b]).groupby("cohort"):
                if len(g) >= 4:
                    x, y = centred(g[a].to_numpy(float), g.stratum), centred(g[b].to_numpy(float), g.stratum)
                    if x.std() > 0 and y.std() > 0:
                        rr.append((len(g), float(np.corrcoef(x, y)[0, 1])))
            n3, r3 = np.array(rr).T
            out.append(dict(pair=f"{a} ~ {b}", patients=int(n3.sum()), r_within_cohort=float((n3 * r3).sum() / n3.sum())))
        for lab, covs in (("endothelial slope", B.BASECOV), ("+ intrinsic ISG", B.BASECOV + ["isg"]),
                          ("+ contamination", B.BASECOV + ["contam"]), ("+ both", B.BASECOV + ["isg", "contam"]),
                          ("+ IFN (A419) + contamination", B.BASECOV + ["ifn", "contam"])):
            r = B.describe(df, "immune_hot", covs, None, [4211, len(covs)], loud_col="loud_endo")
            out.append(dict(pair=lab, patients=r["patients"], slope=r["slope"], lo90=r["lo90"], hi90=r["hi90"], p_two_sided=r["p_two_sided"]))
        T = pd.DataFrame(out)
        T.to_csv(OUT / "s3_intrinsic.csv", index=False)
        pd.set_option("display.width", 250)
        log("A421\n" + T.round(4).to_string(index=False))
        return
    o = B.slope(df, "immune_hot", B.BASECOV, "loud_endo")[0]
    sl = {}
    for f in FAMS:
        v = []
        for parts in sets[f]:
            if not parts:
                v.append(np.nan)
                continue
            d2 = df.drop(columns="loud_endo").merge(pd.concat(parts).rename("loud_endo"), left_on="donor", right_index=True, how="left")
            v.append(B.slope(d2, "immune_hot", B.BASECOV, "loud_endo")[0])
        sl[f] = np.array(v)
    pd.DataFrame(sl).to_csv(OUT / "s1_sets.csv", index=False)
    rows_ = []
    plac = sl["placebo"][np.isfinite(sl["placebo"])]
    for f in ("matched", "responsive", "responsive_broad"):
        v = sl[f][np.isfinite(sl[f])]
        pct = float(np.mean(v <= o))
        pl_pct = np.array([np.mean(v <= x) for x in plac])
        rows_.append(dict(family=f, n_sets=len(v), observed=o, median=float(np.median(v)), q05=float(np.quantile(v, .05)),
                          percentile_low=100 * pct, p_cal=float((np.sum(pl_pct <= pct) + 1) / (len(pl_pct) + 1))))
    S1 = pd.DataFrame(rows_)
    S1.to_csv(OUT / "s1.csv", index=False)
    pd.set_option("display.width", 250, "display.max_columns", 30)
    log("S1\n" + S1.round(4).to_string(index=False))
    # S2 effect size
    n_, r_ = np.array(rel).T
    rel_endo = float((n_ * r_).sum() / n_.sum())
    rng = np.random.default_rng([4192])
    out2 = []
    for lab, col in (("endothelial state", "loud_endo"), ("composite host", "loud_pf")):
        d = df.dropna(subset=[col, "immune_hot"]).reset_index(drop=True)
        est = tertile_ratio(d, col)
        bs = np.array([tertile_ratio(d.loc[B.boot_index(d, rng)].reset_index(drop=True), col) for _ in range(N_BOOT)])
        out2.append(dict(loudness=lab, log_ratio_hot_minus_cold=est, amplitude_ratio=float(np.exp(est)),
                         ratio_lo90=float(np.exp(np.nanquantile(bs, .05))), ratio_hi90=float(np.exp(np.nanquantile(bs, .95))),
                         state_reliability=rel_endo if col == "loud_endo" else np.nan))
    S2 = pd.DataFrame(out2)
    S2.to_csv(OUT / "s2.csv", index=False)
    log("S2\n" + S2.round(4).to_string(index=False))
    # S3 mechanism readouts
    out3 = []
    for a, b in (("ifn", "loud_endo"), ("ifn", "immune_hot"), ("adh", "loud_endo"), ("adh", "immune_hot")):
        rr = []
        for _, g in df.dropna(subset=[a, b]).groupby("cohort"):
            if len(g) >= 4:
                x, y = centred(g[a].to_numpy(float), g.stratum), centred(g[b].to_numpy(float), g.stratum)
                if x.std() > 0 and y.std() > 0:
                    rr.append((len(g), float(np.corrcoef(x, y)[0, 1])))
        n3, r3 = np.array(rr).T
        out3.append(dict(pair=f"{a} ~ {b}", patients=int(n3.sum()), r_within_cohort=float((n3 * r3).sum() / n3.sum())))
    for lab, covs in (("endothelial slope", B.BASECOV), ("+ IFN-gamma score", B.BASECOV + ["ifn"]),
                      ("+ adhesion score", B.BASECOV + ["adh"]), ("+ both", B.BASECOV + ["ifn", "adh"])):
        r = B.describe(df, "immune_hot", covs, None, [4193, len(covs)], loud_col="loud_endo")
        out3.append(dict(pair=lab, patients=r["patients"], slope=r["slope"], lo90=r["lo90"], hi90=r["hi90"], p_two_sided=r["p_two_sided"]))
    S3 = pd.DataFrame(out3)
    S3.to_csv(OUT / "s3.csv", index=False)
    log("S3\n" + S3.round(4).to_string(index=False))
    log(f"done {(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    main()
