"""A335 part 3 (registered): equivalence bounds for H293.3 (cell-class phase order) and H295.1 (composition and loudness),
from the saved A293 Part 2 atlas offsets and the saved A295 pooled estimate. (H315.1 is in `a335_gtex.py`.)

    python a335_3_equivalence.py
Output: `../output/part3_equivalence.csv`.
"""
from pathlib import Path

import numpy as np
import pandas as pd

EXP = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parents[1] / "output"
W = 2 * np.pi / 24


def circ_mean(h):
    return float(np.angle(np.mean(np.exp(1j * W * np.asarray(h)))) / W)


def wrap(h):
    return (np.asarray(h) + 12) % 24 - 12


def main():
    rows = []
    # H293.3: the A293 Part 2 bootstrap, reproduced from the six atlas offsets (seed 293, first class pair, 2,000 draws)
    at = pd.read_csv(EXP / "A293_organisation_levels/output/part2/atlases.csv")
    off = at[at.class_pair == "endothelial-fibroblast_stromal"].offset_h.to_numpy()
    rng = np.random.default_rng(293)
    boot = wrap([circ_mean(off[rng.integers(0, len(off), len(off))]) for _ in range(2000)])
    lo, hi = np.quantile(boot, [0.05, 0.95])
    orig = pd.read_csv(EXP / "A293_organisation_levels/output/part2/tests.csv").iloc[0]
    for margin in (1.0, 0.6):
        rows.append(dict(test="H293.3", statistic="pooled endothelial minus fibroblast offset (h)", atlases=len(off),
                         estimate=circ_mean(off), lo90=lo, hi90=hi, margin=margin,
                         reproduces_A293=bool(np.isclose(lo, orig.boot_lo90) and np.isclose(hi, orig.boot_hi90)),
                         equivalent=bool(-margin < lo and hi < margin), primary=margin == 1.0))
    # H295.1: random-effects pooled relative-gain change per SD of lymphocyte score (margin post hoc)
    t = pd.read_csv(EXP / "A295_composition/output/tests.csv").set_index("score").loc["lymphocyte"]
    est, se = float(t.pooled_gain_change), float(t.gain_se)
    rows.append(dict(test="H295.1", statistic="pooled relative-gain change per SD lymphocyte score", atlases=int(t.tissues),
                     estimate=est, lo90=est - 1.645 * se, hi90=est + 1.645 * se, margin=0.10,
                     equivalent=bool(-0.10 < est - 1.645 * se and est + 1.645 * se < 0.10), primary=True,
                     note="margin post hoc"))
    o = pd.DataFrame(rows)
    OUT.mkdir(parents=True, exist_ok=True)
    o.to_csv(OUT / "part3_equivalence.csv", index=False)
    print(o.round(4).to_string())


if __name__ == "__main__":
    main()
