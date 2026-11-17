"""A295 (registered in `paper/docs/registered_experiments.md`): within a tissue, does a sample's composition change the
gain of its clock state but not its phase?

GTEx fast deaths, A104 periphery rows. Per tissue (A292 design: states x, leave-one-out host h, strata intercepts) and
score s (lymphocyte or stromal-vascular marker score, standardised within tissue): complex least squares
x = mu + c s + a h + b (s h). Gain change per SD Re(b/a); phase change per SD Im(b/a) x 24 / (2 pi) h. Donor bootstrap
SE (200); random-effects pooling over tissues. H295.1: pooled lymphocyte gain change < 0 (one-sided p < 0.05).
H295.2: pooled phase change within +-0.5 h (90 % interval) for both scores. Aggregates only.

    python score_a295.py [--n-boot 200]
Outputs (`../output/`): tissues.csv, tests.csv.
"""
import argparse
import importlib.util
import sys
import warnings
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[3]))
sys.path.insert(0, str(HERE.parents[3] / "results" / "figure_2" / "scripts" / "backbone"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import norm  # noqa: E402

import clocksync as cs  # noqa: E402
import compute_oscillator_gtex as A104  # noqa: E402
from clocksync import stats as ST  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import gtex  # noqa: E402

spec = importlib.util.spec_from_file_location("a295_a292", HERE.parents[2] / "A292_external_tracking/scripts/score_a292.py")
A292 = importlib.util.module_from_spec(spec)
sys.modules["a295_a292"] = A292
spec.loader.exec_module(A292)

OUT = HERE.parents[1] / "output"
W = 2 * np.pi / 24
SCORES = {"lymphocyte": ["CD3D", "CD3E", "CD2", "CD247", "MS4A1", "CD79A", "CD19"],
          "stromal_vascular": ["DCN", "LUM", "COL1A2", "PECAM1", "VWF", "CDH5"]}
BOUND_H = 0.5


def fit(x, h, s, G):
    Z = np.column_stack([G, s, h, s * h])
    coef, *_ = np.linalg.lstsq(Z, x, rcond=None)
    a, b = coef[-2], coef[-1]
    r = b / a
    return float(np.real(r)), float(np.imag(r) / W)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--n-boot", type=int, default=200)
    a = ap.parse_args()
    warnings.filterwarnings("ignore")
    limit_threads(1)
    OUT.mkdir(parents=True, exist_ok=True)
    markers = sorted({g for v in SCORES.values() for g in v})
    clock = list(cs.gtex().genes)
    pbm, donors = gtex.load(min_donors=A104.MIN_DONORS, genes=clock + markers)
    pbm = pbm[~pbm.compartment.str.startswith("Brain")].reset_index(drop=True)
    sc = {}
    for name, g in SCORES.items():
        v = np.log((pbm[g].sum(1) + 0.5) / pbm.depth)
        sc[name] = v.groupby(pbm.compartment).transform(lambda y: (y - y.mean()) / y.std())
    sc = pd.DataFrame(sc).set_index([pbm.donor.astype(str), pbm.compartment])
    pb = pbm[["donor", "compartment", "depth"] + clock]
    ds, pt, comps, D = A292.prepared(pb, donors.tod, "periphery")
    rows = []
    for k in comps:
        e = D[k]
        for name in SCORES:
            s = sc[name].reindex(list(zip(e["donor"], [k] * len(e["donor"])))).to_numpy(float)
            ok = np.isfinite(s)
            x, h, G, s_ = e["x"][ok], e["h"][ok], e["G"][ok], s[ok]
            g0, p0 = fit(x, h, s_, G)
            rng = np.random.default_rng([295, comps.index(k), list(SCORES).index(name)])
            bs = []
            for _ in range(a.n_boot):
                i = rng.integers(0, len(x), len(x))
                try:
                    bs.append(fit(x[i], h[i], s_[i], G[i]))
                except np.linalg.LinAlgError:
                    continue
            bs = np.array(bs)
            rows.append(dict(tissue=k, score=name, n=int(ok.sum()), gain_change=g0, gain_se=float(np.nanstd(bs[:, 0], ddof=1)),
                             phase_change_h=p0, phase_se_h=float(np.nanstd(bs[:, 1], ddof=1))))
    t = pd.DataFrame(rows)
    t.to_csv(OUT / "tissues.csv", index=False)
    tests = []
    for name in SCORES:
        g = t[t.score == name]
        mg = ST.meta(g.gain_change.to_numpy(), g.gain_se.to_numpy())
        mp = ST.meta(g.phase_change_h.to_numpy(), g.phase_se_h.to_numpy())
        z = mg["random"] / mg["se_random"]
        p_lo, p_hi = mp["random"] - 1.645 * mp["se_random"], mp["random"] + 1.645 * mp["se_random"]
        tests.append(dict(score=name, tissues=len(g), pooled_gain_change=mg["random"], gain_se=mg["se_random"],
                          gain_p_one_sided_negative=float(norm.cdf(z)), gain_tau=float(np.sqrt(mg["tau2"])),
                          pooled_phase_change_h=mp["random"], phase_lo90=p_lo, phase_hi90=p_hi,
                          phase_tau_h=float(np.sqrt(mp["tau2"])), phase_within_bound=bool(-BOUND_H < p_lo and p_hi < BOUND_H)))
    tests = pd.DataFrame(tests)
    tests["H295_1"] = np.where(tests.score == "lymphocyte",
                               np.where((tests.pooled_gain_change < 0) & (tests.gain_p_one_sided_negative < 0.05), "PASS", "FAIL"), "")
    tests.to_csv(OUT / "tests.csv", index=False)
    print(t.round(3).to_string(index=False))
    print(tests.round(3).T.to_string())
    print("H295.2", "PASS" if tests.phase_within_bound.all() else "FAIL")


if __name__ == "__main__":
    main()
