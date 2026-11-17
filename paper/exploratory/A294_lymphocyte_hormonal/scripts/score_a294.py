"""A294: is the lymphocyte clock-gene coherence of breast and heart (A281's 2 of 9) carried by the hormone-induced genes?

Registered in paper/docs/registered_experiments.md (A294). The nine A281 pairs (lymphocyte, fibroblast-stromal on shared
donors), A281's statistic and loaders (signed coherence C, per-gene permutations within strata):
  core loop  C on the 11-gene waveform without PER1, PER2, NR1D1 (`Waveform.subset`), 1000 permutations; power control:
             fibroblast-stromal rows thinned to the lymphocytes' clock-gene count (20 replicates, 200 permutations each).
  H294.1     breast CD4 T and heart lymphocytes, where the power control passes (>= 80 %): core-loop p < 0.05.
  H294.2     the 14-gene C above the 95th percentile of 200 stress-responsive sets.
Descriptive: 5th / 50th / 95th percentile of A281's 200 matched sets (reproduces A281's median and percentile).

    python score_a294.py [--jobs 6]
Outputs: paper/exploratory/A294_lymphocyte_hormonal/output/{pairs,tests}.csv
"""

from __future__ import annotations

import importlib.util
import sys
import warnings
from pathlib import Path

PAPER = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(PAPER))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

from clocksync.decompose import signed_coherence  # noqa: E402
from common import cli  # noqa: E402
from common.paths import exploratory_dir  # noqa: E402

spec = importlib.util.spec_from_file_location("a294_a281", PAPER / "exploratory/A281_lymphocyte_clock/scripts/score_a281.py")
A281 = importlib.util.module_from_spec(spec)
sys.modules["a294_a281"] = A281
spec.loader.exec_module(A281)
import cloudpickle  # noqa: E402

for _m in (A281, A281.A280, A281.A274, A281.A272, A281.A150):
    cloudpickle.register_pickle_by_value(_m)

AMENDMENT = "A294_lymphocyte_hormonal"
SEED = 294
INDUCED = ("PER1", "PER2", "NR1D1")
CORE = A281.WF.subset([g for g in A281.CLOCK if g not in INDUCED])
DECIDE = ("organ_breast", "heart")


def pair_task(n, ly, sv, tl, ts, matched, responsive, seed):
    warnings.filterwarnings("ignore")
    rng = np.random.default_rng(seed)
    d = sorted(set(tl.donor) & set(ts.donor))
    tl = tl[tl.donor.isin(d)].sort_values("donor").reset_index(drop=True)
    ts = ts[ts.donor.isin(d)].sort_values("donor").reset_index(drop=True)
    target = float(tl[A281.CLOCK].sum(1).median())
    y14 = A281.y_of(tl, A281.CLOCK)
    c14 = signed_coherence(y14, A281.WF)
    core = list(CORE.genes)
    c11, p11 = A281.perm_test(A281.y_of(tl, core), A281.groups_of(tl), A281.N_PERM, rng, wf=CORE)
    thin_p = []
    for _ in range(A281.N_THIN):
        tt = A281.thin(ts, target, rng)
        thin_p.append(A281.perm_test(A281.y_of(tt, core), A281.groups_of(tt), A281.N_PERM_THIN, rng, wf=CORE)[1])
    pm = np.array([signed_coherence(A281.y_of(tl, list(w.genes)), w) for w in matched])
    pr = np.array([signed_coherence(A281.y_of(tl, list(w.genes)), w) for w in responsive])
    pm, pr = pm[np.isfinite(pm)], pr[np.isfinite(pr)]
    return dict(atlas=n, lymphocyte=ly, stromal=sv, shared_donors=len(d), C14=c14, C11=c11, p11=p11,
                core_power=float(np.mean(np.array(thin_p) < 0.05)),
                matched_q05=float(np.quantile(pm, 0.05)), matched_median=float(np.median(pm)),
                matched_q95=float(np.quantile(pm, 0.95)), matched_pct=float(np.mean(pm < c14)),
                responsive_q95=float(np.quantile(pr, 0.95)), responsive_median=float(np.median(pr)),
                responsive_pct=float(np.mean(pr < c14)), responsive_valid=len(pr))


def main():
    ap = cli.parser(__doc__)
    args = ap.parse_args()
    cli.setup(args)
    out = exploratory_dir(AMENDMENT)
    out.mkdir(parents=True, exist_ok=True)
    A281.A150.ATLASES.update(A281.ATLASES)
    a281 = pd.read_csv(exploratory_dir("A281_lymphocyte_clock") / "pairs.csv")
    tasks = []
    for k, n in enumerate(A281.ATLASES):
        if n not in set(a281.atlas):
            continue
        r = a281.set_index("atlas").loc[n]
        A = A281.A150.load_atlas(n, None, 0)
        pb = A["pb"].reset_index(drop=True)
        matched = A281.A272.pseudo_sets(A, k, A281.N_PSEUDO, "matched")
        responsive = A281.A272.pseudo_sets(A, k, A281.N_PSEUDO, "responsive")
        cols = sorted({g for w in matched + responsive for g in w.genes} | set(A281.CLOCK))
        sub = pb[["donor", "compartment", "stratum", "n_cells", "depth"] + cols]
        tasks.append(delayed(pair_task)(n, r.lymphocyte, r.stromal, sub[sub.compartment == r.lymphocyte],
                                        sub[sub.compartment == r.stromal], matched, responsive, [SEED, k]))
        print(n, "queued", flush=True)
    pairs = pd.DataFrame(Parallel(n_jobs=args.jobs)(tasks))
    rep = pairs.merge(a281[["atlas", "C_lymph", "matched_C_median", "lymph_pct_among_matched"]], on="atlas")
    pairs["repro_C"] = np.isclose(rep.C14, rep.C_lymph)
    pairs["repro_matched"] = np.isclose(rep.matched_median, rep.matched_C_median) & np.isclose(rep.matched_pct,
                                                                                                 rep.lymph_pct_among_matched)
    pairs.to_csv(out / "pairs.csv", index=False)
    rows = []
    for n in DECIDE:
        r = pairs.set_index("atlas").loc[n]
        informative = bool(r.core_power >= 0.8)
        h1 = bool(r.p11 < 0.05) if informative else None
        h2 = bool(r.C14 > r.responsive_q95)
        verdict = ("undecided (power control fails)" if not informative else
                   "hormone-driven, not a lymphocyte clock" if (not h1 and not h2) else
                   "a weak lymphocyte clock signal" if (h1 and h2) else "undecided (mixed)")
        rows.append(dict(atlas=n, lymphocyte=r.lymphocyte, core_power=r.core_power, informative=informative, C11=r.C11,
                         p11=r.p11, H294_1=h1, C14=r.C14, responsive_q95=r.responsive_q95, responsive_pct=r.responsive_pct,
                         H294_2=h2, verdict=verdict))
    tests = pd.DataFrame(rows)
    tests.to_csv(out / "tests.csv", index=False)
    pd.set_option("display.width", 250)
    print(pairs.round(3).to_string())
    print(tests.round(3).to_string())


if __name__ == "__main__":
    main()
