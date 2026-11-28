"""A416 (POST HOC, DESCRIPTIVE; registered in paper/docs/registered_experiments.md, "A416"): is the colder-louder lead present
within stromal cell types, or a composition effect of the stroma composite?

Pool, hosts, readouts and covariates as A403c.
  K1  per-class phase-free loudness (fibroblast, endothelium, mural; class state on its own rows) against immune-hot, with
      matched / responsive pseudo-waveform class loudness as companions; per cancer          -> output/a416/k1.csv
  K2  A403c phase-free composite slope with composition covariates                           -> output/a416/k2.csv
  K3  fibroblast subtype scores (myCAF-like, iCAF-like): correlations and adjusted slope      -> output/a416/k3.csv
    python a416.py [--jobs 12] [--genome-jobs 4]                        (long: paper/common/run.sh)
"""
from __future__ import annotations

import argparse
import time

import numpy as np
import pandas as pd

import a403_a405_score as S
import a403b as B
import a403c as C3
from a402_a405 import BASE, META, U, log, trim, centred, cancer_hosts, loudness_frame, load_slim, readouts
from clocksync import pairwise
from clocksync.config import limit_threads
from clocksync.genes import canonical

OUT = U.OUT / "a416"
CLASSES = {"fibroblast": "tumour:fibroblast", "endothelial": "tumour:endothelial", "mural": "tumour:mural"}
IMM = [canonical(g) for g in ("PTPRC", "LYZ", "CD3E", "CD74")]
MYCAF = [canonical(g) for g in ("ACTA2", "TAGLN", "POSTN", "COL11A1", "MMP11")]
ICAF = [canonical(g) for g in ("IL6", "CXCL12", "PDGFRA", "CFD", "DPT")]
N_SET = 200


def zcoh(df):
    coh = np.array([d.split("|")[0] for d in df.index])
    g = df.groupby(coh)
    return (df - g.transform("mean")) / g.transform(lambda x: x.std(ddof=0)).replace(0, np.nan)


def score(t, genes):
    g = [x for x in genes if x in t.columns]
    if len(g) < len(genes) / 2:
        return pd.Series(np.nan, index=t.index)
    lc = np.log2(t[g].div(t.depth, axis=0) * 1e6 + 1)
    return zcoh(lc).mean(axis=1)


def class_loud(rows, wf, cancer):
    """Phase-free log radius of a class's state (rows: one per donor), centred within stratum; strata of one dropped."""
    st = rows.stratum.astype(str).to_numpy()
    keep = pd.Series(st).map(pd.Series(st).value_counts()).to_numpy() >= 2
    rows, st = rows[keep], st[keep]
    if len(rows) < 6:
        return pd.Series(dtype=float)
    z = pairwise.states(rows, wf)
    z = centred(z, pd.Series(st))
    ok = np.abs(z) > 0
    lr = pd.Series(np.log(np.abs(z[ok])), index=rows.donor.to_numpy()[ok])
    lr = pd.Series(centred(lr.to_numpy(), pd.Series(st[ok])), index=lr.index)
    pf, _ = C3.phase_free(lr, np.angle(z[ok]), st[ok], np.array([cancer] * ok.sum()))
    return pf.dropna()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--jobs", type=int, default=12)
    ap.add_argument("--genome-jobs", type=int, default=4)
    args = ap.parse_args()
    args.smoke = False
    limit_threads(1)
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
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
    assert len(hosts) == 12 and len(loud) == 574
    slim, errs = load_slim(sources, donors, f388, args)
    assert len(slim) == 51 and errs.empty
    ro = readouts(slim, loud)
    df = loud.merge(ro.drop(columns="cohort"), on="donor").merge(S.covariates(hosts), on="donor", how="left")
    df["cohort"] = df.donor.str.split("|").str[0]
    assert abs(B.slope(df, "immune_hot", B.BASECOV)[0] - (-0.10228)) < 5e-4, "A403 not reproduced"
    df["loud_pf"], _ = C3.phase_free(df.loud, np.angle(df.z.to_numpy()), df.stratum.to_numpy(), df.cancer.to_numpy())
    log(f"table {len(df)}; A403 reproduced")

    # per-class loudness (observed and companions), composition covariates, fibroblast subtype scores
    from clocksync.paired import pseudo_waveforms
    obs = {k: [] for k in CLASSES}
    comp = {k: {"matched": [[] for _ in range(N_SET)], "responsive": [[] for _ in range(N_SET)]} for k in CLASSES}
    extra = []
    for c, h in hosts.items():
        Gc = h["Gc"]
        genes = [g for g in Gc.columns if g not in META]
        frac = Gc[genes].div(Gc.depth, axis=0)
        lm = np.log(frac.mean() + 1e-12)[(Gc[genes].mean() >= 1) | frac.columns.isin(list(BASE.genes))]
        ws = {fam: pseudo_waveforms(lm, N_SET, np.random.default_rng([U.SEED, code]), BASE, fam)
              for code, fam in enumerate(("matched", "responsive"))}
        for k, comp_name in CLASSES.items():
            rows = Gc[Gc.compartment == comp_name].drop_duplicates("donor")
            if rows.donor.nunique() < 6:
                continue
            obs[k].append(class_loud(rows[META + list(BASE.genes)], BASE, c))
            for fam, wl in ws.items():
                for j, w in enumerate(wl):
                    try:
                        comp[k][fam][j].append(class_loud(rows[META + list(w.genes)], w, c))
                    except Exception:  # noqa: BLE001
                        pass
        n = Gc.pivot_table(index="donor", columns="compartment", values="n_cells", aggfunc="sum").fillna(0)
        tot = Gc.groupby("donor")[[g for g in IMM if g in Gc.columns] + ["depth"]].sum()
        fib = Gc[Gc.compartment == CLASSES["fibroblast"]].drop_duplicates("donor").set_index("donor")
        e = pd.DataFrame(index=n.index)
        e["log_fib_endo"] = np.log((n.get(CLASSES["fibroblast"], 0) + .5) / (n.get(CLASSES["endothelial"], 0) + .5))
        e["log_stroma_cells"] = np.log(n.sum(axis=1))
        e["immune_contam"] = score(tot, IMM)
        e["mycaf"] = score(fib, MYCAF).reindex(e.index)
        e["icaf"] = score(fib, ICAF).reindex(e.index)
        extra.append(e.rename_axis("donor").reset_index())
        log(f"classes {c}")
    for k in CLASSES:
        df = df.merge((pd.concat(obs[k]) if obs[k] else pd.Series(dtype=float)).rename(f"loud_{k}"), left_on="donor", right_index=True, how="left")
    df = df.merge(pd.concat(extra, ignore_index=True), on="donor", how="left")
    pd.set_option("display.width", 300, "display.max_columns", 40)
    # K1
    K1 = []
    for k in CLASSES:
        cm = {fam: [pd.concat(p) if p else pd.Series(dtype=float) for p in sets] for fam, sets in comp[k].items()}
        for r in B.per_cancer(df, "immune_hot", B.BASECOV, cm, f"k1{k}", loud_col=f"loud_{k}"):
            K1.append(dict(r, cls=k))
        log(f"K1 {k} done")
    K1 = pd.DataFrame(K1)
    K1.to_csv(OUT / "k1.csv", index=False)
    log("K1 pooled\n" + K1[K1.cancer == "pooled"].round(4).to_string(index=False))
    # K2
    K2 = []
    for lab, covs in (("A403c phase-free (reference)", B.BASECOV),
                      ("+ composition", B.BASECOV + ["log_fib_endo", "log_stroma_cells", "immune_contam"])):
        r = B.describe(df, "immune_hot", covs, None, [4161, len(covs)], loud_col="loud_pf")
        K2.append(dict(r, model=lab))
    K2 = pd.DataFrame(K2)
    K2.to_csv(OUT / "k2.csv", index=False)
    log("K2\n" + K2.round(4).to_string(index=False))
    # K3
    K3 = []
    for a, b in (("mycaf", "immune_hot"), ("icaf", "immune_hot"), ("mycaf", "loud_fibroblast"), ("icaf", "loud_fibroblast"),
                 ("immune_contam", "immune_hot"), ("log_fib_endo", "immune_hot")):
        rr = []
        for s, g in df.dropna(subset=[a, b]).groupby("cohort"):
            if len(g) >= 4:
                x, y = centred(g[a].to_numpy(float), g.stratum), centred(g[b].to_numpy(float), g.stratum)
                if x.std() > 0 and y.std() > 0:
                    rr.append((len(g), float(np.corrcoef(x, y)[0, 1])))
        if rr:
            n_, r_ = np.array(rr).T
            K3.append(dict(pair=f"{a} ~ {b}", patients=int(n_.sum()), r_within_cohort=float((n_ * r_).sum() / n_.sum())))
    r = B.describe(df, "immune_hot", B.BASECOV + ["mycaf", "icaf"], None, [4162], loud_col="loud_fibroblast")
    if r:
        K3.append(dict(pair="fibroblast slope adjusted for mycaf + icaf", patients=r["patients"], slope=r["slope"],
                       lo90=r["lo90"], hi90=r["hi90"], p_two_sided=r["p_two_sided"]))
    K3 = pd.DataFrame(K3)
    K3.to_csv(OUT / "k3.csv", index=False)
    log("K3\n" + K3.round(4).to_string(index=False))
    log(f"done {(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    main()
