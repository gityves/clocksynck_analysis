"""A402c (POST HOC, DESCRIPTIVE; registered in paper/docs/registered_experiments.md, "A402c"): precision of the single-cell
tumour-stroma time reference and the composition rhythm A402b could detect.

  R1  per cancer, 20 gene-half splits: half-host states, split-half complex correlation, Spearman-Brown reliability of the
      full host, approximate phase error SD in hours                        -> output/a402c/r1.csv (+ r1_splits.csv)
  S1  power of the A402b statistic for a rhythm A cos(phi_true - phi_0) added to the within-cohort-permuted share, at the
      observed host phase (ceiling) and with the R1 phase error (measured error) -> output/a402c/s1.csv
    python a402c.py [--jobs 12] [--genome-jobs 4]                        (long: paper/common/run.sh)
"""
from __future__ import annotations

import argparse
import time
import warnings

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

import a402b as A  # sets the a402b cache before a402_a405 is imported
from a402_a405 import BASE, META, U, cs, log, centred
from clocksync import deviation
from clocksync.config import limit_threads

OUT = U.OUT / "a402c"
N_SPLIT, N_REP, N_PERM = 20, 200, 1000
AMPS = (0.125, 0.25, 0.5, 0.75)


def half_state(tab, classes, genes):
    limit_threads(1)
    warnings.filterwarnings("ignore")
    try:
        h = cs.Dataset(tab[META + list(genes)], anchors=classes, waveform=BASE.subset(list(genes)), mixed_strata=True).host()
        return pd.Series((h.x + 1j * h.y).to_numpy(), index=h.index)
    except Exception:  # noqa: BLE001  (a half that cannot be fitted: no value for this split)
        return pd.Series(dtype=complex)


def r1(hosts, args):
    tasks = []
    for c, h in hosts.items():
        for j, (ga, gb) in enumerate(deviation.split_genes(BASE, N_SPLIT, np.random.default_rng([4023, sum(map(ord, c))]))):
            tasks += [(c, j, "A", ga), (c, j, "B", gb)]
    res = Parallel(n_jobs=args.jobs)(delayed(half_state)(hosts[c]["tab"], hosts[c]["classes"], g) for c, j, s, g in tasks)
    got = {(c, j, s): r for (c, j, s, _), r in zip(tasks, res)}
    rows = []
    for c, h in hosts.items():
        st = h["stratum"]
        for j in range(N_SPLIT):
            a, b = got[(c, j, "A")], got[(c, j, "B")]
            ids = a.index.intersection(b.index)
            if len(ids) < 10:
                continue
            s = st.reindex(ids)
            za, zb = centred(a.reindex(ids).to_numpy(), s), centred(b.reindex(ids).to_numpy(), s)
            rho = float(abs(np.sum(za * np.conj(zb))) / np.sqrt(np.sum(abs(za) ** 2) * np.sum(abs(zb) ** 2)))
            rows.append(dict(cancer=c, split=j, patients=len(ids), rho_half=rho))
    S = pd.DataFrame(rows)
    S.to_csv(OUT / "r1_splits.csv", index=False)
    R = S.groupby("cancer").agg(patients=("patients", "max"), splits=("split", "size"), rho_half=("rho_half", "median"),
                                 rho_q10=("rho_half", lambda x: x.quantile(.1)), rho_q90=("rho_half", lambda x: x.quantile(.9))).reset_index()
    R["rel_full"] = 2 * R.rho_half / (1 + R.rho_half)
    R["phase_sd_h"] = 12 / np.pi * np.sqrt((1 - R.rel_full) / (2 * R.rel_full))
    w = R.patients.to_numpy(float)
    R = pd.concat([R, pd.DataFrame([dict(cancer="pooled", patients=int(w.sum()), splits=int(R.splits.sum()),
                                         rho_half=float(np.sum(w * R.rho_half) / w.sum()), rel_full=float(np.sum(w * R.rel_full) / w.sum()),
                                         phase_sd_h=float(np.sum(w * R.phase_sd_h) / w.sum()))])], ignore_index=True)
    R.to_csv(OUT / "r1.csv", index=False)
    return R


def s1(df, R):
    sd_rad = {r.cancer: r.phase_sd_h * np.pi / 12 for r in R.itertuples() if r.cancer != "pooled"}
    new = df[df.new]
    ords = A.orders(new, N_PERM)
    rows = []
    for col in A.COMP:
        un = A.with_strata(A.cohort_units(new, col), new)
        prep = []
        for s, u in un.items():
            perms = np.argsort(ords[s].loc[u["donor"]].to_numpy().T, axis=1)
            prep.append(dict(n=u["n"], y=u["y"], z=u["z"], strat=u["strat"], perms=perms, sd=sd_rad[u["cancer"]],
                             phi=np.angle(u["z"])))
        w = np.array([p["n"] for p in prep], float)
        for ref in ("ceiling", "measured error"):
            for A_ in AMPS:
                rng = np.random.default_rng([4024, sum(map(ord, col + ref)), int(A_ * 1000)])
                hits = []
                for _ in range(N_REP):
                    phi0 = rng.uniform(0, 2 * np.pi)
                    obs, nul = 0.0, np.zeros(N_PERM)
                    for p in prep:
                        phit = p["phi"] + (rng.normal(0, p["sd"], p["n"]) if ref == "measured error" else 0.0)
                        y = p["y"][rng.permutation(p["n"])] + A_ * np.cos(phit - phi0)
                        y = centred(y, p["strat"])
                        den = np.sqrt(np.sum(y ** 2) * np.sum(abs(p["z"]) ** 2))
                        obs += p["n"] * abs(np.sum(y * np.conj(p["z"]))) / den
                        nul += p["n"] * np.abs((y[p["perms"]] * np.conj(p["z"])[None, :]).sum(1)) / den
                    obs, nul = obs / w.sum(), nul / w.sum()
                    hits.append((np.sum(nul >= obs) + 1) / (N_PERM + 1) < 0.05)
                rows.append(dict(readout=col, reference=ref, amplitude_log2=A_, fold_peak_trough=float(2 ** (2 * A_)),
                                 patients=int(w.sum()), cohorts=len(prep), power=float(np.mean(hits)), reps=N_REP))
                log(rows[-1])
    T = pd.DataFrame(rows)
    T.to_csv(OUT / "s1.csv", index=False)
    return T


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--jobs", type=int, default=12)
    ap.add_argument("--genome-jobs", type=int, default=4)
    args = ap.parse_args()
    args.smoke = False
    limit_threads(1)
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    df, hosts = A.build(args)
    R = r1(hosts, args)
    from joblib.externals.loky import get_reusable_executor
    get_reusable_executor().shutdown(wait=True)
    pd.set_option("display.width", 250)
    log("R1\n" + R.round(4).to_string(index=False))
    T = s1(df, R)
    log("S1\n" + T.round(4).to_string(index=False))
    log(f"done {(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    main()
