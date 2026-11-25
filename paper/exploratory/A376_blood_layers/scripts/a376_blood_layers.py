"""A376 (registered 2026-11-25): does the time-of-death signal of GTEx whole blood sit in the glucocorticoid /
stress-responsive layer rather than in the clock genes?

GTEx v10 fast deaths with recorded time of death, the 20 tissues of `blood_composition_posthoc` (whole blood, testis,
13 brain regions, five solid tissues). Per tissue, on log2-normalised counts of expressed genes (mean log2 >= 3):
per-gene cosinor R^2 on time of death; a gene set's rhythmicity = median over its genes of the gene's percentile among
all expressed genes. Sets: the clock genes of the GTEx waveform; 200 expression-matched and 200 responsive sets
(clocksync controls, the rules behind `paired.pseudo_waveforms`; seed 376); descriptive: the expressed glucocorticoid
targets as one set, and every statistic after regressing out the blood-cell composition scores.
H376.1 (blood): clock <= responsive q05 and responsive median > matched q95. H376.2: clock > responsive q95 in >= 15 of
the 19 other tissues. Controlled access: only per-tissue aggregates are written.

    python a376_blood_layers.py      # about a minute
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
sys.path[:0] = [str(PAPER), str(PAPER / "results/figure_1/scripts/method"),
                str(PAPER / "exploratory/blood_composition_posthoc/scripts")]
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from blood_composition import CLOCK, LOW_GAIN, MARKERS, TISSUES, cosinor, z  # noqa: E402
from _shared import load  # noqa: E402
from clocksync._core.controls import GC_TARGETS, pseudo_clock_genes, responsive_pseudo_genes  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common.benchmark.data import log2_normalised  # noqa: E402

OUT = HERE.parents[1] / "output"
N_SETS, SEED = 200, 376


def rhythmicity(pct: pd.Series, sets) -> np.ndarray:
    return np.array([np.median(pct[list(s)]) for s in sets])


def tissue(t, counts, tod, rng):
    tod = tod.dropna()
    L = log2_normalised(counts.loc[tod.index])
    L = L.loc[:, L.mean() >= 3]
    h = tod.to_numpy(float)
    comp = np.column_stack([z(L[[g for g in v if g in L.columns]]).mean(1) for v in MARKERS.values()
                            if sum(g in L.columns for g in v) >= 2])
    cg = [g for g in CLOCK if g in L.columns]
    gc = [g for g in GC_TARGETS if g in L.columns]
    lm = L.mean()
    matched = pseudo_clock_genes(lm, N_SETS, rng, genes=cg)
    responsive = responsive_pseudo_genes(lm, N_SETS, rng, genes=cg)
    pool = sorted({g for s in responsive for g in s})
    rows = []
    for adj, C in (("none", None), ("composition", comp)):
        r2, _, _ = cosinor(L.to_numpy(float), h, C)
        pct = pd.Series(pd.Series(r2).rank(pct=True).to_numpy() * 100, index=L.columns)
        m, r = rhythmicity(pct, matched), rhythmicity(pct, responsive)
        rows.append(dict(tissue=t, adjustment=adj, donors=len(L), clock_genes=len(cg), responsive_pool=len(pool),
                         gc_genes=len(gc), clock=float(np.median(pct[cg])), gc_targets=float(np.median(pct[gc])),
                         matched_median=np.median(m), matched_q95=np.quantile(m, 0.95),
                         responsive_q05=np.quantile(r, 0.05), responsive_median=np.median(r),
                         responsive_q95=np.quantile(r, 0.95), clock_share_below_responsive=float(np.mean(r > np.median(pct[cg])))))
    return rows


def main():
    limit_threads()
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    tissues = list(TISSUES) + list(LOW_GAIN)
    data = load(tissues)
    res = pd.DataFrame([r for t in tissues if t in data for r in tissue(t, *data[t], rng)])
    res.to_csv(OUT / "tissues.csv", index=False)
    a = res[res.adjustment == "none"].set_index("tissue")
    b = a.loc["Whole Blood"]
    h1a, h1b = b.clock <= b.responsive_q05, b.responsive_median > b.matched_q95
    other = a.drop("Whole Blood")
    n2 = int((other.clock > other.responsive_q95).sum())
    tests = pd.DataFrame([dict(hypothesis="H376.1", part_a_clock_le_responsive_q05=bool(h1a),
                               part_b_responsive_gt_matched_q95=bool(h1b), verdict="PASS" if h1a and h1b else "FAIL"),
                          dict(hypothesis="H376.2", tissues_clock_gt_responsive_q95=n2, of=len(other),
                               verdict="PASS" if n2 >= 15 else "FAIL")])
    tests.to_csv(OUT / "tests.csv", index=False)
    pd.set_option("display.width", 250)
    print(res.round(3).to_string())
    print(tests.to_string())


if __name__ == "__main__":
    main()
