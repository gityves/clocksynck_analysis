"""A81b H81.2 and H81.3: the GTEx cue set in healthy timed skin and in psoriasis.

Pre-registered in docs/PLAN_internal_desync.md (archive) (A81b; written after the A81 selection and evaluation, before any A81b skin,
psoriasis or other-disease scoring). The cue reference is `paper/results/figure_3/output/gtex_cues/a81b_cue_reference.csv`
(cue = low lock to the host in ventilated donors, rhythmic across tissues in fast deaths; locked = high lock; 40 genes each;
frame check failed at 0.48 and was overridden post hoc in A81b, as in the legacy).
  H81.2 (GSE205155, timed): per subject, epidermis-minus-dermis offset of the core clock and of the cue state (A76 complex
    regression on clock time); mean core minus mean cue > 0, subject bootstrap (10,000) one-sided p < 0.05.
  H81.3 (psoriasis): the A78 estimator (gated A54 series, weights 1 / se^2, identical per-pair label swaps across sets, 10,000)
    with the cue plane as reference: core minus cue < 0 at one-sided p < 0.05. Sets: core, cue, locked, A78 systemic.
Needs compute_replication.py run with the same --subsample (gated series and weights).

    python compute_cues_gtex.py [--subsample N] [--n-perm 10000] [--n-boot 10000]
Outputs (paper/results/figure_5/output/psoriasis/): a81b_h2_verdict.csv, a81b_timed_skin_subjects_summary.csv, a81b_h3_sets.csv,
a81b_h3_verdict.csv
Legacy: legacy/desync/validation/paired_damping/scripts/a82_cues_diseases.py (`timed_skin`, `main`), a81b_cues.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import pandas as pd  # noqa: E402

import _cues as Q  # noqa: E402
import _shared as H  # noqa: E402
from common import cli  # noqa: E402

SEED = 78  # the legacy swaps used the A78 stream for every A78-type comparison


def main(argv=None):
    args = H.parser(__doc__).parse_args(argv)
    run = cli.setup(args)
    out = H.out_dir(H.FOLDER, run["quick"])
    n = H.sizes(args, run, 10000, 10000)
    wfs = Q.cue_waveforms(H.reference_sets()["systemic"])
    v2, d = Q.timed_core_cue(wfs, n["n_boot"])
    pd.Series(v2).to_csv(out / "a81b_h2_verdict.csv", header=False)
    d.to_csv(out / "a81b_timed_skin_subjects_summary.csv")  # public GEO data, per subject
    print(pd.Series(v2).to_string(), flush=True)
    samples, w = H.gated_samples(out, args.subsample, args.seed)
    sets, diff = H.compare_cue_sets(samples, w, wfs, "cue", n["n_perm"], SEED)
    c = diff.set_index("set").loc["core"]
    _, d2 = H.compare_cue_sets(samples, w, wfs, "locked", n["n_perm"], SEED)
    o = d2.set_index("set").loc["core"]
    v = dict(core_minus_cue_h=c.difference_h, p_lower=c.p_lower, p_two=c.p_two, null_q025=c.null_q025,
             null_q975=c.null_q975, core_minus_locked_h=o.difference_h, p_two_locked=o.p_two, series=int(c.cohorts),
             pairs=int(c.pairs))
    v["H81_3"] = "PASS" if (v["core_minus_cue_h"] < 0 and v["p_lower"] < 0.05) else "FAIL"
    sets.to_csv(out / "a81b_h3_sets.csv", index=False)
    pd.Series(v).to_csv(out / "a81b_h3_verdict.csv", header=False)
    print(sets.round(3).to_string(index=False), "\n", pd.Series(v).to_string())


if __name__ == "__main__":
    main()
