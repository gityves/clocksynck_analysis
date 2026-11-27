"""A414 (registered in paper/docs/registered_experiments.md, "A414"): the immune programs that follow the surgery hour in A146
(P2 MHC-I, P4 T-cell activation, P10 PD-1) against the tumour-stroma clock phase.

Pool, hosts, statistic, shared permutations and pseudo-host companions as A402b (`a402b.py`); program scores by
`compute_a246.readout` on the A402 slim tables (`a402_a405.cohort_dict`).
                                                       -> output/a414/{tests,per_cancer}.csv
    python a414.py [--jobs 12] [--genome-jobs 4]                         (long: paper/common/run.sh)
"""
from __future__ import annotations

import argparse
import time

import numpy as np
import pandas as pd

import a402b as A  # sets the a402b cache before a402_a405 is imported
import a402_a405 as P
from a402_a405 import U, log, trim, cancer_hosts, loudness_frame, load_slim

OUT = U.OUT / "a414"
PROGS = {"P2": "H414.1", "P4": "H414.2", "P10": "H414.3"}
A146_H = {"P2": 17.0, "P4": 13.0, "P10": 23.0}


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
    z_of = loud.set_index("donor").z
    rows = []
    for s, G in slim.items():
        ids = set(G.donor) & set(z_of.index)
        if len(ids) < 4:
            continue
        d = P.cohort_dict(G, ids)
        r = pd.DataFrame(index=pd.Index(d["pts"], name="donor"))
        for prog in PROGS:
            y, why = P.P246.readout(d, prog)
            r[prog] = y if not why else np.nan
        rows.append(r.reset_index())
    ro = pd.concat(rows, ignore_index=True)
    df = loud.merge(ro, on="donor", how="inner")
    df["cohort"] = df.donor.str.split("|").str[0]
    df["new"] = ~df.cohort.isin(set(P.M.cohort_list()))
    log(f"table: {len(df)} patients ({int(df.new.sum())} new); " + ", ".join(f"{p} {df[p].notna().sum()}" for p in PROGS))
    return df, hosts


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--jobs", type=int, default=12)
    ap.add_argument("--genome-jobs", type=int, default=4)
    args = ap.parse_args()
    args.smoke = False
    from clocksync.config import limit_threads
    limit_threads(1)
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    df, hosts = build(args)
    ords = A.orders(df, A.N_PERM)
    Z = A.companion_states(hosts, args)
    from joblib.externals.loky import get_reusable_executor
    get_reusable_executor().shutdown(wait=True)
    rows, pc = [], []
    for lab, sub in (("primary", df[df.new]), ("all", df)):
        for prog, hyp in PROGS.items():
            un = A.with_strata(A.cohort_units(sub, prog), sub)
            if not un:
                continue
            r, _ = A.test_row(un, ords, Z, lab, prog, dict(hypothesis=hyp if lab == "primary" else "descriptive", a146_hour=A146_H[prog]))
            r["peak_lo90"], r["peak_hi90"] = A.boot_peak(un, np.random.default_rng([414, 1, sum(map(ord, prog + lab))]))
            rows.append(r)
            if lab == "all":
                for c in sorted({u["cancer"] for u in un.values()}):
                    uc = {s: u for s, u in un.items() if u["cancer"] == c}
                    rc, _ = A.test_row(uc, ords, None, lab, prog, dict(cancer=c))
                    pc.append(rc)
            log(f"{lab} {prog} done")
    T = pd.DataFrame(rows)
    prim = T[T.hypothesis.str.startswith("H")]
    order = prim.p.sort_values().index
    holm, run = {}, 0.0
    for k, i in enumerate(order):
        run = max(run, min(1.0, (len(order) - k) * T.loc[i, "p"]))
        holm[i] = run
    T["holm"] = pd.Series(holm)
    T["pass"] = ((T.holm < 0.05) & (T.r > T.matched_q95)).where(T.hypothesis.str.startswith("H"))
    T.to_csv(OUT / "tests.csv", index=False)
    pd.DataFrame(pc).to_csv(OUT / "per_cancer.csv", index=False)
    pd.set_option("display.width", 300, "display.max_columns", 40)
    log("tests\n" + T.round(4).to_string(index=False))
    log(f"done {(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    main()
