"""A82: core clock vs timing cues in other diseases (GTEx cue reference).

Pre-registered in docs/PLAN_internal_desync.md (archive) (A82; sets from A81b). The A78 estimator (identical per-pair label swaps across
gene sets, pooled offsets with the A54 weights 1 / se^2 of the gated series) with the cue plane as reference, in four units:
atopic dermatitis (A53, 8 gated series), other skin diseases pooled (A59, standard-scale series), inflamed gut (A55, gated
series) and TCGA tumour vs adjacent normal (10 cancers, A52 gating from cancer). Sets: core, cue, locked (A81b
reference `paper/results/figure_3/output/gtex_cues/a81b_cue_reference.csv`), systemic (A78). Estimand: core minus cue
pooled offset (h, wrapped), two-sided swap p (n_perm 10,000), Holm over the four units (H82: offset if p_holm < 0.05); core
minus locked as the control. Secondary, descriptive: A59 conditions with >= 30 gated pairs (2,000 swaps), IBD CD and UC
(2,000), each TCGA cancer (1,000).
Needs compute_other_diseases.py (and cancer/compute_tcga_paired.py for TCGA) run with the same --subsample.

    python compute_cues_diseases.py [--subsample N] [--n-perm 10000]
Outputs (paper/results/figure_5/output/psoriasis_details/): a82_<unit>_{sets,verdict}.csv, a82_summary.csv
Legacy: legacy/desync/validation/paired_damping/scripts/a82_cues_diseases.py (+ a79_companion.py loaders, a78 `compare`)
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))  # paper/
from common.paths import group_dir, group_out  # noqa: E402
sys.path.insert(0, str(group_dir("psoriasis")))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import _cues as Q  # noqa: E402
import _details as D  # noqa: E402
import _shared as H  # noqa: E402
from clocksync import stats as S  # noqa: E402
from common import cli  # noqa: E402


SEED = 78


def compare(samples, w, wfs, tag, n_perm, out, prefix="a82", ref="cue", other="locked"):
    """Pooled offsets per set with identical swaps; core minus `ref` (and minus `other`); writes <prefix>_<tag>_*."""
    sets, d1 = H.compare_cue_sets(samples, w, wfs, ref, n_perm, SEED)
    _, d2 = H.compare_cue_sets(samples, w, wfs, other, n_perm, SEED)
    c, o = d1.set_index("set").loc["core"], d2.set_index("set").loc["core"]
    v = dict(core_minus_cue_h=c.difference_h, p_lower=c.p_lower, p_two=c.p_two, null_q025=c.null_q025, null_q975=c.null_q975,
             core_minus_locked_h=o.difference_h, p_two_locked=o.p_two, series=int(c.cohorts), pairs=int(c.pairs))
    sets.to_csv(out / f"{prefix}_{tag}_sets.csv", index=False)
    pd.Series(v).to_csv(out / f"{prefix}_{tag}_verdict.csv", header=False)
    print(tag, "\n", sets.round(3).to_string(index=False), "\n", pd.Series(v).to_string(), flush=True)
    return v


def holm_summary(prim: dict, out, prefix, label):
    p = pd.Series({k: x["p_two"] for k, x in prim.items()}).sort_values()
    holm = np.minimum(1, np.maximum.accumulate(p.to_numpy() * (len(p) - np.arange(len(p)))))
    summ = pd.DataFrame(prim).T.assign(p_holm=pd.Series(holm, index=p.index))
    summ[label] = np.where(summ.p_holm < 0.05, "offset", "no offset")
    summ.to_csv(out / f"{prefix}_summary.csv")
    print(summ.round(3).to_string())
    return summ


def cancer_table(quick: bool):
    """A52 per-cancer gating of cancer (documented output): the subsample table in quick mode."""
    for f in ((group_out("cancer") / "subsample" / "a52_per_cancer.csv",) if quick else ()) + (
            group_out("cancer") / "a52_per_cancer.csv",):
        if f.exists():
            return pd.read_csv(f)
    return None


def units(args, run, out, wfs):
    """The four primary units and the descriptive splits: {tag: (samples, weights)}."""
    sets_genes = Q.all_set_genes(wfs)
    prim, sec = {}, {}
    disc, _ = D.discovery_all(args.subsample, args.seed)
    ad = disc[disc.disease == "atopic_dermatitis"]
    prim["atopic_dermatitis"] = Q.gated_cohorts(ad, pd.read_csv(out / "a53_skin_per_series.csv").query(
        "disease == 'atopic_dermatitis'"))
    other, _ = D.other_skin(args.subsample, args.seed)
    other = other[other.source == "standard"]
    s, w = Q.gated_cohorts(other, pd.read_csv(out / "a59_all_per_series.csv"))
    prim["other_skin_pooled"] = (s, w)
    cond = s.drop_duplicates("cohort").set_index("cohort").disease.astype(str)
    for c in cond.unique():
        cs_ = [k for k in cond.index[cond == c] if k in w]
        sub = s[s.cohort.isin(cs_)]
        if sub.groupby("cohort").unit.nunique().sum() >= (10 if run["quick"] else 30):
            sec[f"skin_{c}"] = (sub, {k: w[k] for k in cs_}, 2000)
    ibd, _ = D.ibd(args.subsample, args.seed)
    ps = pd.read_csv(out / "a55_ibd_per_series.csv")
    prim["ibd"] = Q.gated_cohorts(ibd, ps)
    for dis in ("CD", "UC"):
        p2 = pd.read_csv(out / f"a55_ibd_{dis}_per_series.csv")
        sec[f"ibd_{dis}"] = (*Q.gated_cohorts(ibd[ibd.disease == dis], p2), 2000)
    pc = cancer_table(run["quick"])
    if pc is not None:
        t, w = Q.tcga_table(pc, sets_genes, args.subsample, args.seed)
        prim["tcga"] = (t, w)
        for c in w:
            sec[f"tcga_{c}"] = (t[t.cohort == c], {c: w[c]}, 1000)
    else:
        print("TCGA arm skipped: cancer/compute_tcga_paired.py outputs (a52_per_cancer.csv) not found", flush=True)
    return prim, sec


def main(argv=None):
    args = H.parser(__doc__).parse_args(argv)
    run = cli.setup(args)
    out = H.out_dir(D.FOLDER, run["quick"])
    n_perm = H.sizes(args, run, 10000, 1)["n_perm"]
    wfs = Q.cue_waveforms(H.reference_sets()["systemic"])
    prim, sec = units(args, run, out, wfs)
    res = {k: compare(*u, wfs, k, n_perm, out) for k, u in prim.items()}
    holm_summary(res, out, "a82", "H82")
    for k, (s, w, n) in sec.items():
        compare(s, w, wfs, k, min(n, n_perm), out)


if __name__ == "__main__":
    main()
