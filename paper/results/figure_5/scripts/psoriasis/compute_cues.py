"""A78: is the psoriatic lesion delay in the local clock or in the systemic timing cues?

Pre-registered in docs/PLAN_internal_desync.md (archive) (A78 H78.1). Reference planes from the timed healthy skin (GSE205155,
A67): a systemic gene set (rhythmic in epidermis and dermis, in phase within 0.5 h, top 40 epidermal amplitude) and an
epidermal-output set (rhythmic in the epidermis only, top 40), each with its own epidermal acrophases and amplitudes;
the core clock is the 14-gene GTEx waveform. The A54 estimator is applied per gated psoriasis series and set (one
projector per series and set), pooled with the A54 weights (1 / se^2), with identical per-pair label swaps across sets
(n_perm). Estimand: core minus systemic pooled offset (h, wrapped). H78.1: < 0 and one-sided swap p < 0.05
(delay stronger in the core clock than in the cue plane).

    python compute_cues.py [--subsample N]      (needs compute_replication.py run with the same --subsample)
Outputs: paper/results/figure_5/output/psoriasis/a78_gene_sets.csv, a78_psoriasis_sets.csv, a78_psoriasis_verdict.csv
Legacy: legacy/desync/validation/paired_damping/scripts/a78_local_vs_cues.py (+ a67_timed_layers.py for the sets)
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import pandas as pd  # noqa: E402
import _shared as H  # noqa: E402
from common import cli  # noqa: E402

SEED = 78


def main(argv=None):
    args = H.parser(__doc__).parse_args(argv)
    run = cli.setup(args)
    out = H.out_dir(H.FOLDER, run["quick"])
    n_perm = H.sizes(args, run)["n_perm"]
    sets = H.reference_sets(out)
    samples, w = H.gated_samples(out, args.subsample, args.seed)
    print(len(w), "gated series,", samples.groupby("cohort").unit.nunique().sum(), "units", flush=True)
    # reference = systemic: differences of core and of epidermal output relative to it
    s1, d1 = H.compare_cue_sets(samples, w, sets, "systemic", n_perm, SEED)
    d1 = d1.set_index("set")
    core = d1.loc["core"]
    other = d1.loc["epidermal_output"]
    v = dict(core_minus_systemic_h=core.difference_h, p_lower=core.p_lower, p_two=core.p_two,
             null_q025=core.null_q025, null_q975=core.null_q975)
    # core minus epidermal output (two-sided): compare_sets differences are set minus reference, so recompute with
    # the epidermal-output set as the reference on the same seed (identical swaps)
    s2, d2 = H.compare_cue_sets(samples, w, sets, "epidermal_output", n_perm, SEED)
    c2 = d2.set_index("set").loc["core"]
    v.update(core_minus_other_h=c2.difference_h, p_two_other=c2.p_two, series=int(core.cohorts), pairs=int(core.pairs))
    v["H78_1"] = "PASS" if (v["core_minus_systemic_h"] < 0 and v["p_lower"] < 0.05) else "FAIL"
    s1.to_csv(out / "a78_psoriasis_sets.csv", index=False)
    pd.Series(v).to_csv(out / "a78_psoriasis_verdict.csv", header=False)
    print(s1.round(3).to_string(index=False), "\n", pd.Series(v).to_string())


if __name__ == "__main__":
    main()
