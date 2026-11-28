"""A422 (POST HOC, DESCRIPTIVE; registered in paper/docs/registered_experiments.md, "A422"): the colder-louder slope in low- vs
high-immune-contamination halves of the stromal (composite) and endothelial pseudobulks.
                                                     -> output/a422/d1.csv
    python a422.py [--jobs 4] [--genome-jobs 2]                          (long: paper/common/run.sh)
"""
from __future__ import annotations

import argparse
import time

import numpy as np
import pandas as pd

import a403b as B
import a416 as K
import a419 as E
from a402_a405 import BASE, META, U, log
from clocksync.config import limit_threads

OUT = U.OUT / "a422"
CONTAM = K.IMM  # PTPRC, LYZ, CD3E, CD74


def boot_diff(d, col, loud, rng, n=2000):
    """Bootstrap (within cohort) of slope(low half) - slope(high half)."""
    out = []
    for _ in range(n):
        db = d.loc[B.boot_index(d, rng)].reset_index(drop=True)
        lo = B.slope(db[db[col] == "low"], "immune_hot", B.BASECOV, loud)[0]
        hi = B.slope(db[db[col] == "high"], "immune_hot", B.BASECOV, loud)[0]
        out.append(lo - hi)
    return np.array(out)


def split(df, score, lab):
    q = df.groupby("cohort")[score].rank(pct=True)
    df[lab] = np.where(q <= 0.5, "low", np.where(q > 0.5, "high", None))
    df[lab + "_t"] = np.where(q > 2 / 3, "top", "rest")
    return df


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--genome-jobs", type=int, default=2)
    args = ap.parse_args()
    args.smoke = False
    limit_threads(1)
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    df, hosts = E.build(args)
    obs, sc = [], []
    for c, h in hosts.items():
        Gc = h["Gc"]
        tot = Gc.groupby("donor")[[g for g in CONTAM if g in Gc.columns] + ["depth"]].sum()
        s = pd.DataFrame({"contam_stroma": K.score(tot, CONTAM)})
        rows = Gc[Gc.compartment == E.ENDO].drop_duplicates("donor")
        if rows.donor.nunique() >= 6:
            obs.append(K.class_loud(rows[META + list(BASE.genes)], BASE, c))
            s = s.join(pd.DataFrame({"contam_endo": K.score(rows.set_index("donor"), CONTAM)}), how="left")
        sc.append(s.rename_axis("donor").reset_index())
    df = df.merge(pd.concat(obs).rename("loud_endo"), left_on="donor", right_index=True, how="left")
    df = df.merge(pd.concat(sc, ignore_index=True), on="donor", how="left")
    rng = np.random.default_rng([4221])
    rows = []
    for lab, loud, score in (("composite", "loud_pf", "contam_stroma"), ("endothelium", "loud_endo", "contam_endo")):
        d = split(df.dropna(subset=[loud, score, "immune_hot"]).copy(), score, "half")
        for part, sub in (("all", d), ("low contamination half", d[d.half == "low"]), ("high contamination half", d[d.half == "high"]),
                          ("without top contamination tertile", d[d.half_t == "rest"])):
            r = B.describe(sub.reset_index(drop=True), "immune_hot", B.BASECOV, None, [4222, len(part), len(lab)], loud_col=loud)
            rows.append(dict(loudness=lab, subset=part, patients=r["patients"], cohorts=r["cohorts"], slope=r["slope"], lo90=r["lo90"],
                             hi90=r["hi90"], p_two_sided=r["p_two_sided"]))
        bd = boot_diff(d.reset_index(drop=True), "half", loud, rng)
        rows.append(dict(loudness=lab, subset="low minus high", slope=float(np.nanmedian(bd)), lo90=float(np.nanquantile(bd, .05)),
                         hi90=float(np.nanquantile(bd, .95))))
        log(f"{lab} done")
    T = pd.DataFrame(rows)
    T.to_csv(OUT / "d1.csv", index=False)
    pd.set_option("display.width", 250)
    log("D1\n" + T.round(4).to_string(index=False))
    log(f"done {(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    main()
