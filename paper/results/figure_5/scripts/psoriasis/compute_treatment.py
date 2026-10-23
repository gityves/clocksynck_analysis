"""A56 H56.2: does the psoriatic lesion delay resolve after treatment? (and the post hoc treated-vs-placebo contrast)

Pre-registered in docs/PLAN_internal_desync.md (archive) (A56 H56.2). Patients with lesional and uninvolved biopsies at baseline
and after treatment (GSE117468 week 12, GSE201827 week 12; treated arms). Per series the lesional-vs-uninvolved offset
(A54 estimator, one projector per visit pair) at baseline and post; pooled over series with weights 1 / se^2 of the
baseline variance ratio (patient bootstrap, n_boot). Estimand: pooled baseline offset minus pooled post offset (h,
wrapped; < 0: the delay resolves). Null: baseline/post labels swapped per patient (n_perm), one-sided p (lower).
Secondary: placebo arms, GSE201827 week 52, GSE117468 by drug. Post hoc (not pre-specified, labelled): treated-minus-
placebo change of the offset, arm labels permuted within series.

    python compute_treatment.py [--subsample N]
Outputs: paper/results/figure_5/output/psoriasis/a56_treatment.csv, a56_treatment_series.csv, a56_interaction_posthoc.csv
Legacy: legacy/desync/validation/paired_damping/scripts/a56_treatment.py, a56_interaction_posthoc.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import _shared as H  # noqa: E402
from common import cli  # noqa: E402

PD = H.PD
RUNS = {
    "treated_primary": [("GSE117468.parquet", True, None), ("GSE201827.parquet", True, None)],
    "placebo": [("GSE117468.parquet", False, None), ("GSE201827.parquet", False, None)],
    "GSE201827_w52_treated": [("GSE201827_w52.parquet", True, None)],
    "brodalumab": [("GSE117468.parquet", True, "brodalumab")],
    "ustekinumab": [("GSE117468.parquet", True, "ustekinumab")],
}


def thin_visits(t: dict, n, seed):
    """--subsample: keep n random patients of the aligned visit tables."""
    if n is None:
        return t
    idx = t[("baseline", "lesional")].index
    keep = np.random.default_rng(seed).permutation(sorted(idx))[:n]
    m = idx.isin(keep)
    return {k: v[m] for k, v in t.items()}


def swap_visits(b, p, sw):
    f = lambda x: sw[:, None] if x.ndim == 2 else sw
    return (tuple(np.where(f(x), y, x) for x, y in zip(b, p)), tuple(np.where(f(x), x, y) for x, y in zip(b, p)))


def analysis(specs, label, rng, n_perm, n_boot, n, seed, min_n):
    sets = []
    for file, treated, drug in specs:
        t, counts = H.visits(file, treated, drug)
        t = thin_visits(t, n, seed)
        k = len(t[("baseline", "lesional")])
        if k < min_n:
            continue
        b, p = H.visit_states(t, counts)
        se = PD.boot_se(*b, n_boot=n_boot, rng=rng)
        sets.append((file, b, p, 1 / se ** 2, k))
    w = np.array([s[3] for s in sets])
    ob = PD.pooled_offset([PD.ccorr(*s[1]) for s in sets], w)
    op = PD.pooled_offset([PD.ccorr(*s[2]) for s in sets], w)
    stat = H.wrap1(ob - op)
    null = []
    for _ in range(n_perm):
        cb, cp = [], []
        for _, b, p, _, k in sets:
            bb, pp = swap_visits(b, p, rng.random(k) < 0.5)
            cb.append(PD.ccorr(*bb))
            cp.append(PD.ccorr(*pp))
        null.append(H.wrap1(PD.pooled_offset(cb, w) - PD.pooled_offset(cp, w)))
    null = np.asarray(null)
    per = [dict(set=label, series=f, n=k, offset_baseline_h=H.offset_h(b), offset_post_h=H.offset_h(p),
                r_baseline=PD.log_var_ratio(*b), r_post=PD.log_var_ratio(*p)) for f, b, p, _, k in sets]
    v = dict(set=label, n_patients=int(sum(s[4] for s in sets)), offset_baseline_h=ob, offset_post_h=op,
             baseline_minus_post_h=stat, p_lower=float((1 + np.sum(null <= stat)) / (1 + n_perm)))
    return v, per


def interaction(rng, n_perm, n_boot, n, seed):
    """Post hoc: treated minus placebo change (baseline - post) of the pooled offset; arm labels permuted in series."""
    data = []
    for f in ("GSE117468.parquet", "GSE201827.parquet"):
        t, counts = H.visits(f, None)
        t = thin_visits(t, n, seed)
        arm = t[("baseline", "lesional")].arm.astype(str).str.lower().str.contains("placebo").to_numpy()
        b, p = H.visit_states(t, counts)
        data.append((b, p, arm, 1 / PD.boot_se(*b, n_boot=n_boot, rng=rng) ** 2))
    w = np.array([d[3] for d in data])

    def change(sel):
        ob = PD.pooled_offset([PD.ccorr(*H.sub(d[0], m)) for d, m in zip(data, sel)], w)
        op = PD.pooled_offset([PD.ccorr(*H.sub(d[1], m)) for d, m in zip(data, sel)], w)
        return H.wrap1(ob - op)

    def stat(arms):
        tr, pl = change([~a for a in arms]), change(list(arms))
        return H.wrap1(tr - pl), tr, pl

    obs, tr, pl = stat([d[2] for d in data])
    null = np.array([stat([rng.permutation(d[2]) for d in data])[0] for _ in range(n_perm)])
    return dict(change_treated_h=tr, change_placebo_h=pl, treated_minus_placebo_h=obs,
                p_lower=float((1 + np.sum(null <= obs)) / (1 + n_perm)))


def main(argv=None):
    args = H.parser(__doc__).parse_args(argv)
    run = cli.setup(args)
    out = H.out_dir(H.FOLDER, run["quick"])
    z = H.sizes(args, run, 2000, 500)
    n_perm, n_boot = z["n_perm"], z["n_boot"]
    min_n = 3 if run["quick"] else 5
    rng = np.random.default_rng(562)
    rows, per = [], []
    for label, specs in RUNS.items():
        v, p = analysis(specs, label, rng, n_perm, n_boot, args.subsample, args.seed, min_n)
        if label == "treated_primary":
            v["H56_2"] = "PASS" if (v["baseline_minus_post_h"] < 0 and v["p_lower"] < 0.05) else "FAIL"
        rows.append(v)
        per += p
        print(v, flush=True)
    pd.DataFrame(rows).to_csv(out / "a56_treatment.csv", index=False)
    pd.DataFrame(per).to_csv(out / "a56_treatment_series.csv", index=False)
    v = interaction(np.random.default_rng(563), n_perm, n_boot, args.subsample, args.seed)
    v["label"] = "post hoc"
    pd.Series(v).to_csv(out / "a56_interaction_posthoc.csv", header=False)
    print(v)


if __name__ == "__main__":
    main()
