"""A79: specificity companion of A78 (other skin diseases and atopic dermatitis against the timing-cue planes).

Pre-registered in docs/PLAN_internal_desync.md (archive) (A79 H79.1). The A78 comparison (core clock vs the systemic timed-skin cue
plane, identical per-pair label swaps across sets; see psoriasis/compute_cues.py) applied to the gated series of
the other skin conditions (A59, pooled) and of atopic dermatitis (A53, secondary). H79.1 (specificity): the other
diseases do NOT show core minus systemic < 0 at one-sided p < 0.05 (PASS = no such contrast).

    python compute_cues_companion.py [--subsample N]   (needs compute_other_diseases.py run with the same --subsample)
Outputs: paper/results/figure_5/output/psoriasis_details/a79_other_diseases_{sets,verdict}.csv,
  a79_atopic_dermatitis_{sets,verdict}.csv, a79_verdict.csv
Legacy: legacy/desync/validation/paired_damping/scripts/a79_companion.py (+ a78_local_vs_cues.compare)
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import pandas as pd  # noqa: E402
import _details as D  # noqa: E402
import _shared as H  # noqa: E402
from common import cli  # noqa: E402

SEED = 78


def compare(samples, w, sets, n_perm, tag, out):
    s1, d1 = H.compare_cue_sets(samples, w, sets, "systemic", n_perm, SEED)
    _, d2 = H.compare_cue_sets(samples, w, sets, "epidermal_output", n_perm, SEED)
    c, o = d1.set_index("set").loc["core"], d2.set_index("set").loc["core"]
    v = dict(core_minus_systemic_h=c.difference_h, p_lower=c.p_lower, p_two=c.p_two, null_q025=c.null_q025,
             null_q975=c.null_q975, core_minus_other_h=o.difference_h, p_two_other=o.p_two, series=int(c.cohorts),
             pairs=int(c.pairs))
    s1.to_csv(out / f"a79_{tag}_sets.csv", index=False)
    pd.Series(v).to_csv(out / f"a79_{tag}_verdict.csv", header=False)
    print(tag, "\n", s1.round(3).to_string(index=False), "\n", pd.Series(v).to_string(), flush=True)
    return v


def gated(samples, per_series: pd.DataFrame):
    ps = per_series[per_series.gated.astype(bool)].set_index("cohort")
    s = samples[samples.cohort.isin(ps.index)]
    return s, {k: 1 / ps.loc[k, "se"] ** 2 for k in ps.index if k in set(s.cohort)}


def main(argv=None):
    args = H.parser(__doc__).parse_args(argv)
    run = cli.setup(args)
    out = H.out_dir(D.FOLDER, run["quick"])
    n_perm = H.sizes(args, run)["n_perm"]
    sets = H.reference_sets()
    other, _ = D.other_skin(args.subsample, args.seed)
    other = other[other.source == "standard"]
    ps = pd.read_csv(out / "a59_all_per_series.csv")
    v = compare(*gated(other, ps), sets, n_perm, "other_diseases", out)
    v["H79_1"] = "PASS" if not (v["core_minus_systemic_h"] < 0 and v["p_lower"] < 0.05) else "FAIL"
    pd.Series(v).to_csv(out / "a79_verdict.csv", header=False)
    disc, _ = D.discovery_all(args.subsample, args.seed)
    ad = disc[disc.disease == "atopic_dermatitis"]
    pa = pd.read_csv(out / "a53_skin_per_series.csv")
    pa = pa[pa.disease == "atopic_dermatitis"]
    compare(*gated(ad, pa), sets, n_perm, "atopic_dermatitis", out)


if __name__ == "__main__":
    main()
