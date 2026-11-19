"""A335 part 4 (registered): does the GTEx clock still read time of death when its waveform is fitted on other donors?

The shipped waveform was fitted on the same 244 fast deaths that A270 tests on. The donors are split once into halves
(seed [335, 4]); per half, the waveform is refitted on the other half (`fit_waveform`, tissue as group, as
`compute_gtex_benchmark.cross_fit_waveform`) and the clock's R with time of death on this half is read with
`Dataset.timed` (1,000 permutations), beside the shipped waveform and A270's 200 responsive and 200 matched sets
(same seeds), all on the same half. Aggregates only (GTEx time of death is controlled access).

    python a335_4_a270_crossfit.py [--jobs 4]
Outputs (`../output/`): part4_sets.csv, part4_tests.csv, part4_hours.csv.
"""
import argparse
import sys
import time
from multiprocessing import Pool
from pathlib import Path

HERE = Path(__file__).resolve()
EXP = HERE.parents[2]
sys.path.insert(0, str(EXP / "A270_responsive_time" / "scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import score_a270 as A  # noqa: E402  (sets sys.path for paper/ and the organs scripts)
from clocksync.config import limit_threads  # noqa: E402
from clocksync.paired import pseudo_waveforms  # noqa: E402
from clocksync.waveform import fit_waveform  # noqa: E402

OUT = HERE.parents[1] / "output"
CLOCK = list(A.WF.genes)


def cross_fit(pb, hours, donors):
    """Waveform refitted on the donors not in `donors` (tissue as group)."""
    t = pb[~pb.donor.isin(donors)]
    y = np.log(np.maximum(t[CLOCK].to_numpy(float), 0.5) / t.depth.to_numpy(float)[:, None])
    return fit_waveform(y, hours.reindex(t.donor).to_numpy(float), CLOCK, t.compartment.to_numpy(), name="cross-fitted")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--jobs", type=int, default=4)
    a = ap.parse_args()
    t0 = time.time()
    pb, hours = A.load("human")
    limit_threads(1)
    OUT.mkdir(parents=True, exist_ok=True)
    lm = A.log_mean(pb)
    donors = np.array(sorted(pb.donor.unique()))
    perm = np.random.default_rng([335, 4]).permutation(len(donors))
    halves = {"half1": set(donors[perm[: len(donors) // 2]]), "half2": set(donors[perm[len(donors) // 2:]])}
    hrs = []
    for h, d in halves.items():
        b = np.floor(hours.reindex(sorted(d)).to_numpy(float) / 4).astype(int) * 4
        hrs.append({"half": h, "donors": len(d), **{f"{k:02d}-{k + 4:02d}h": int((b == k).sum()) for k in range(0, 24, 4)}})
    pd.DataFrame(hrs).to_csv(OUT / "part4_hours.csv", index=False)
    print(pd.DataFrame(hrs).to_string(index=False), flush=True)
    comp = []
    for j, kind in enumerate(("responsive", "matched"), 1):
        comp += [(kind, w) for w in pseudo_waveforms(lm, 200, np.random.default_rng([A.SEED, j]), A.WF, kind)]
    genes = sorted({g for _, w in comp for g in w.genes} | set(CLOCK))
    rows, tests = [], []
    for h, d in halves.items():
        wcf = cross_fit(pb, hours, d)
        sets = [("clock_crossfit", wcf), ("clock_shipped", A.WF)] + comp
        sub = pb[pb.donor.isin(d)]
        A.S.update(sets=sets, pb=sub[A.META + genes].copy(), hours=hours.reindex(sorted(d)), n_perm=1000)
        with Pool(a.jobs, initializer=limit_threads, initargs=(1,)) as p:
            res = list(p.imap_unordered(A.task, range(len(sets))))
        R = pd.DataFrame(res).sort_values("set").assign(half=h)
        rows.append(R)
        cf, sh = R[R.kind == "clock_crossfit"].iloc[0], R[R.kind == "clock_shipped"].iloc[0]
        c = R[R.kind.isin(["responsive", "matched"]) & R.R.notna()]
        t = dict(half=h, donors=len(d), R_crossfit=cf.R, p_R_crossfit=cf.p_R, R_shipped=sh.R, p_R_shipped=sh.p_R,
                 companions_valid=len(c), companion_R_max=float(c.R.max()), companion_R_q95=float(c.R.quantile(0.95)),
                 responsive_R_q95=float(c[c.kind == "responsive"].R.quantile(0.95)),
                 matched_R_q95=float(c[c.kind == "matched"].R.quantile(0.95)),
                 pct_below_crossfit=float(100 * (c.R <= cf.R).mean()),
                 waveform_phase_corr=float(np.abs(np.mean(np.exp(1j * (wcf.psi - A.WF.psi))))))
        t["above_all"] = bool(cf.p_R < 0.05 and cf.R > t["companion_R_max"])
        t["above_q95"] = bool(cf.p_R < 0.05 and cf.R > t["companion_R_q95"])
        tests.append(t)
        print(pd.Series(t).to_string(), f"{time.time() - t0:.0f}s", sep="\n", flush=True)
    pd.concat(rows).to_csv(OUT / "part4_sets.csv", index=False)
    T = pd.DataFrame(tests)
    verdict = "holds" if T.above_all.all() else "holds at q95" if T.above_q95.all() else "falls"
    T["verdict"] = verdict
    T.to_csv(OUT / "part4_tests.csv", index=False)
    print("verdict:", verdict, f"runtime {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
