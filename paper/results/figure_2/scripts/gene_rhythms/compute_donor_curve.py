"""A39: is donor count why only the lung atlas predicts gene rhythms? HLCA lung donors subsampled through the A27 pipeline.

Pre-registered (A39, written before running). HLCA lung (200 donors) is the only atlas whose non-clock rhythmic genes match
GTEx (R 0.56 in A31b). Draws of n = 25, 45, 100 donors (10 each, seed 39), the full 200, and (n = 25 only) the same draws with
the class rule lowered to >= 10 donors; the A27 pipeline unchanged (host from the A17 anchors present, `Dataset.gene_rhythms`,
bulk recombination, 20 host permutations per draw), scored against GTEx Lung v1 (set b: p < 0.01 non-clock-plane genes; set c:
clock plane). Primary: at n = 45 the median R over draws exceeds the median null 95 % quantile (if yes, donor count does not
explain the other atlases' failures). Secondary: R vs n, offset stability. A draw with < 2 anchors or classes is infeasible.

    python compute_donor_curve.py [--subsample N]     (N: only checks; 2 draws per n and 10 permutations)
Outputs: paper/results/figure_2/output/gene_rhythms/a39_{draws_scores,summary}.csv
Legacy: gonogo_new/scripts/a39_subsample.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import _shared as S  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import cli  # noqa: E402
from common.paths import out_dir  # noqa: E402

FOLDER = "gene_rhythms"
SEED, NS, NDRAW, NPERM = 39, (25, 45, 100), 10, 20
_CACHE = {}


def draws(donors, ndraw):
    rng = np.random.default_rng(SEED)
    rows = [dict(draw=f"n{n}_d{d}", n=n, min_donors=20, donors=list(rng.choice(donors, n, replace=False)))
            for n in NS for d in range(ndraw)]
    rows.append(dict(draw="n200_full", n=len(donors), min_donors=20, donors=list(donors)))
    rows += [dict(r, draw=r["draw"] + "_md10", min_donors=10) for r in rows if r["n"] == 25]
    return rows


def one(r, genes, truth, n_perm):
    limit_threads(1)
    if "pb" not in _CACHE:
        _CACHE["pb"], _ = S.genome_table("lung_hlca_genome")
    pb = _CACHE["pb"]
    pb = S.class_rule(pb[pb.donor.isin(r["donors"])], r["min_donors"])
    comps = sorted(pb.compartment.unique())
    anchors = S.atlas_anchors("lung_hlca", comps)
    info = dict(draw=r["draw"], n=r["n"], min_donors=r["min_donors"], donors_table=pb.donor.nunique(), classes=len(comps),
                n_anchors=len(anchors), status="ok")
    if len(anchors) < 2 or len(comps) < 2:
        return dict(info, status=f"infeasible ({len(anchors)} anchors, {len(comps)} classes)")
    try:
        pred, nulls, *_ = S.fit_atlas(pb, anchors, genes, n_perm, seed=0)
    except (ValueError, np.linalg.LinAlgError) as e:  # e.g. a leave-one-out host with < 2 anchors
        return dict(info, status=f"failed: {type(e).__name__}: {e}")
    for key, which in (("a", "p<0.01 non-frame"), ("frame", "frame genes")):
        s = S.score_vs_null(pred, nulls, truth, S.gene_set(truth, which))
        info.update({f"{key}_n": s["n"]})
        if "R_err" in s:
            info.update({f"{key}_R": s["R_err"], f"{key}_offset_h": s["R_rot_h"], f"{key}_null95": s["null_R_95"],
                         f"{key}_p": s["p_R"], f"{key}_pass": s["R_err"] > s["null_R_95"]})
    return info


def circmean_h(x):
    x = np.asarray(x, float)
    x = x[~np.isnan(x)]
    return float(np.angle(np.exp(1j * S.W * x).mean()) / S.W % 24) if len(x) else np.nan


def main(argv=None):
    args = cli.parser(__doc__).parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    tr = pd.read_csv(out / "gtex_cosinor.csv.gz")
    genes = sorted(set(tr.gene))
    truth = S.truth_c(tr[tr.tissue == "Lung"])
    pb, _ = S.genome_table("lung_hlca_genome")
    donors = sorted(pb.donor.unique())
    n_perm = 10 if run["quick"] else NPERM
    rows = draws(donors, 2 if run["quick"] else NDRAW)
    per = pd.DataFrame(Parallel(n_jobs=args.jobs)(delayed(one)(r, genes, truth, n_perm) for r in rows))
    per["tag"] = np.where(per.draw.str.endswith("_md10"), "n" + per.n.astype(str) + "_md10", "n" + per.n.astype(str))
    per["feasible"] = per.status == "ok"
    per.to_csv(out / "a39_draws_scores.csv", index=False)
    summ = []
    for tag, q in per.groupby("tag", sort=False):
        f = q[q.feasible]
        g = lambda c: f[c] if c in f else pd.Series(dtype=float)  # noqa: E731
        summ.append(dict(tag=tag, n=int(q.n.iloc[0]), draws=len(q), feasible=len(f), classes_med=g("classes").astype(float).median(),
                         anchors_med=g("n_anchors").astype(float).median(), a_genes_med=g("a_n").median(), a_R_med=g("a_R").median(),
                         a_null95_med=g("a_null95").median(), a_frac_pass=float(g("a_pass").fillna(False).astype(bool).sum() / len(q)),
                         a_offset_circmean_h=circmean_h(g("a_offset_h")), a_offset_median_h=g("a_offset_h").median(),
                         frame_R_med=g("frame_R").median(), frame_null95_med=g("frame_null95").median(),
                         frame_offset_circmean_h=circmean_h(g("frame_offset_h")),
                         frame_frac_pass=float(g("frame_pass").fillna(False).astype(bool).sum() / len(q))))
    summ = pd.DataFrame(summ)
    p = summ.set_index("tag").loc["n45"]
    summ["primary_n45_pass"] = bool(p.feasible >= (2 if run["quick"] else 5) and p.a_R_med > p.a_null95_med)
    summ.to_csv(out / "a39_summary.csv", index=False)
    pd.set_option("display.width", 250)
    print(summ.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
