"""A299 (registered in `paper/docs/registered_experiments.md`): does the medulla content of a GTEx adrenal sample move
its clock phase towards the medulla's?

A295 model on the adrenal only: x = mu + c s + a h + b (s h), the adrenal's leave-one-out host (every peripheral tissue
an anchor, A104 rows), s a marker score standardised within the adrenal. Phase change per SD Im(b/a) x 24 / (2 pi) h
(positive = advanced), gain change Re(b/a); donor bootstrap 2000. Primary donors Hardy 1-4; fast deaths descriptive.
H299.1: medulla phase change > 0, one-sided bootstrap p < 0.05. Aggregates only.

    python score_a299.py
Outputs (`../output/`): tests.csv.
"""
import importlib.util
import sys
import warnings
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[3]))
sys.path.insert(0, str(HERE.parents[3] / "results" / "figure_2" / "scripts" / "backbone"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
import compute_oscillator_gtex as A104  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import gtex  # noqa: E402


def _load(name, rel):
    spec = importlib.util.spec_from_file_location(name, HERE.parents[2] / rel)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


A292 = _load("a299_a292", "A292_external_tracking/scripts/score_a292.py")
A295 = _load("a299_a295", "A295_composition/scripts/score_a295.py")
OUT = HERE.parents[1] / "output"
TISSUE = "Adrenal Gland"
SCORES = {"medulla": ["TH", "DBH", "PNMT", "CHGA", "CHGB"], "cortex": ["STAR", "CYP11A1", "CYP11B1", "CYP17A1", "CYP21A2"]}
N_BOOT = 2000


def run(hardy, label):
    clock = list(cs.gtex().genes)
    markers = sorted({g for v in SCORES.values() for g in v})
    pbm, donors = gtex.load(min_donors=A104.MIN_DONORS, genes=clock + markers, hardy=hardy)
    pbm = pbm[~pbm.compartment.str.startswith("Brain")].reset_index(drop=True)
    ad = pbm[pbm.compartment == TISSUE]
    _, _, comps, D = A292.prepared(pbm[["donor", "compartment", "depth"] + clock], donors.tod, "periphery")
    e = D[TISSUE]
    rows = []
    for name, g in SCORES.items():
        v = np.log((ad[g].sum(1) + 0.5) / ad.depth)
        v = pd.Series(((v - v.mean()) / v.std()).to_numpy(), index=ad.donor.astype(str).to_numpy())
        s = v.reindex(e["donor"]).to_numpy(float)
        ok = np.isfinite(s)
        x, h, G, s_ = e["x"][ok], e["h"][ok], e["G"][ok], s[ok]
        g0, p0 = A295.fit(x, h, s_, G)
        rng = np.random.default_rng([299, list(SCORES).index(name), hardy[0]])
        bs = np.array([A295.fit(x[i], h[i], s_[i], G[i]) for i in (rng.integers(0, len(x), len(x)) for _ in range(N_BOOT))])
        rows.append(dict(donors=label, score=name, n=int(ok.sum()), tissues_in_host=len(comps) - 1,
                         phase_change_h=p0, phase_lo90=float(np.quantile(bs[:, 1], 0.05)), phase_hi90=float(np.quantile(bs[:, 1], 0.95)),
                         p_phase_positive=float((1 + np.sum(bs[:, 1] <= 0)) / (1 + N_BOOT)),
                         gain_change=g0, gain_lo90=float(np.quantile(bs[:, 0], 0.05)), gain_hi90=float(np.quantile(bs[:, 0], 0.95))))
    return rows


def main():
    warnings.filterwarnings("ignore")
    limit_threads(1)
    OUT.mkdir(parents=True, exist_ok=True)
    rows = run((1, 2, 3, 4), "Hardy 1-4 (primary)") + run((1, 2), "fast deaths (descriptive)")
    t = pd.DataFrame(rows)
    t["H299_1"] = ""
    m = (t.donors == "Hardy 1-4 (primary)") & (t.score == "medulla")
    t.loc[m, "H299_1"] = np.where((t.loc[m, "phase_change_h"] > 0) & (t.loc[m, "p_phase_positive"] < 0.05), "PASS", "FAIL")
    t.to_csv(OUT / "tests.csv", index=False)
    print(t.round(3).T.to_string())


if __name__ == "__main__":
    main()
