"""Post hoc check of H189.1 (2026-11-10, not registered): is the coherence of A189 a clock or a shift common to the
clock genes of a sample?

The registered statistic, the unsigned coherence R of `decompose.components`, is near 1 when the clock genes of a
sample move together (the leave-one-gene-out state then points away from the left-out gene: coherence in antiphase).
`decompose.signed_coherence` C fits a per-sample offset on the other genes and keeps the sign. Same compartments,
permutation scheme, pseudo sets and seeds as `score_a189.py`.

    python posthoc_signed.py [--n-perm 10000] [--n-pseudo 200]
"""
import argparse
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import score_a189 as a  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.decompose import logfrac, signed_coherence  # noqa: E402
from clocksync.paired import pseudo_waveforms  # noqa: E402


def signed(tab, wf):
    y = logfrac(tab, list(wf.genes))[0]
    return signed_coherence(y, wf), y


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--n-perm", type=int, default=10000)
    ap.add_argument("--n-pseudo", type=int, default=200)
    args = ap.parse_args()
    limit_threads(1)
    ref = pd.read_csv(a.OUT / "h189_1_coherence.csv").set_index("compartment")
    rows = []
    for i, (comp, level) in enumerate(a.L1.items()):
        tab = a.table(level, [comp], a.GENES)
        c, y = signed(tab, a.WF)
        rng = np.random.default_rng([a.SEED, 1, i])
        null = np.array([signed_coherence(np.column_stack([rng.permutation(g) for g in y.T]), a.WF)
                         for _ in range(args.n_perm)])
        row = dict(compartment=comp, n_fetuses=tab.donor.nunique(), R_registered=ref.R[comp], C=c,
                   p=(1 + np.sum(null >= c)) / (1 + args.n_perm), null_q95=np.quantile(null, .95))
        lm = a.log_mean(level, [comp])
        for kind, tag in (("matched", 2), ("responsive", 3)):
            ws = pseudo_waveforms(lm, args.n_pseudo, np.random.default_rng([a.SEED, tag, i]), a.WF, kind)
            cs_ = np.array([signed(a.table(level, [comp], list(w.genes)), w)[0] for w in ws])
            row[f"{kind}_q95"], row[f"{kind}_pct"] = np.quantile(cs_, .95), 100 * np.mean(cs_ < c)
        logo = {g: signed(tab, a.WF.drop([g]))[0] for g in a.GENES}
        row["logo_min"], row["logo_min_gene"] = min(logo.values()), min(logo, key=logo.get)
        rows.append(row)
        print(comp, round(c, 3), row["p"], flush=True)
    out = pd.DataFrame(rows)
    out["p_holm"] = a.holm(out.p)
    out.to_csv(a.OUT / "h189_1_signed_posthoc.csv", index=False)
    print(out.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
