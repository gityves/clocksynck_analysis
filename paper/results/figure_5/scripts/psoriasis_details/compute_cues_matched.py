"""A84: tissue-matched GTEx cue sets for gut, tumours and skin.

Pre-registered in docs/PLAN_internal_desync.md (archive) (A84). The A81b cue / locked pools (eligible genes of the 11-tissue selection
table, split at the median lock) are matched to the target tissue T from the 26-tissue selection table
(`paper/results/figure_3/output/gtex_cues/selection_genes_26.csv.gz`; rhythm from all fast deaths): candidates expressed
in T with |acrophase_T - global| < 3 h and p_fast,T < 0.2, top 40 by amplitude in T, waveform from T (colon: complex mean of
the transverse and sigmoid segments; fewer than 15 genes: the global A81b set restricted to genes expressed in T). Tissues:
colon (IBD, COAD), liver (LIHC), breast (BRCA), prostate (PRAD), kidney cortex (KIRC, KIRP), uterus (UCEC), lung (LUAD),
thyroid (THCA), skin (psoriasis, atopic dermatitis); BLCA keeps the global sets. H84: IBD (colon sets) and TCGA (per-cancer
sets, pooled with the A52 weights) core minus matched cue, two-sided swap p (n_perm 10,000), Holm over 2. Secondary: each
TCGA cancer (1,000 swaps), psoriasis and atopic dermatitis with the skin sets (2,000).
Needs compute_other_diseases.py, compute_replication.py (psoriasis) and cancer/compute_tcga_paired.py
outputs with the same --subsample.

    python compute_cues_matched.py [--subsample N] [--n-perm 10000]
Outputs (paper/results/figure_5/output/psoriasis_details/): a84_set_<tissue>_{cue,locked}.csv, a84_sets_info.csv,
a84_<unit>_{sets,verdict}.csv, a84_summary.csv
Legacy: legacy/desync/validation/paired_damping/scripts/a84_matched_cues.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))  # paper/
from common.paths import group_dir  # noqa: E402
sys.path.insert(0, str(group_dir("psoriasis")))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import pandas as pd  # noqa: E402

import _cues as Q  # noqa: E402
import _details as D  # noqa: E402
import _shared as H  # noqa: E402
from common import cli  # noqa: E402
from compute_cues_diseases import cancer_table, compare, holm_summary  # noqa: E402

FOLDER6 = H.FOLDER


def matched_sets(out):
    g81b, g26 = Q.selection_tables()
    pool = Q.pools(g81b)
    ref = Q.reference()
    S, info = {}, []
    for k, tis in Q.TISSUE.items():
        S[k], r = Q.matched(g26, pool, tis, ref)
        info.append(r.assign(tissue=k))
        for name in ("cue", "locked"):
            S[k][name].table().assign(tissue=k, set=name).to_csv(out / f"a84_set_{k}_{name}.csv")
    info = pd.concat(info)
    info.to_csv(out / "a84_sets_info.csv", index=False)
    print(info.to_string(), flush=True)
    return S


def main(argv=None):
    args = H.parser(__doc__).parse_args(argv)
    run = cli.setup(args)
    out = H.out_dir(D.FOLDER, run["quick"])
    out6 = H.out_dir(FOLDER6, run["quick"])
    n_perm = H.sizes(args, run, 10000)["n_perm"]
    S = matched_sets(out)
    glob = {k: v for k, v in Q.cue_waveforms().items() if k in ("core", "cue", "locked")}
    genes = Q.all_set_genes({f"{t}_{k}": w for t, d in S.items() for k, w in d.items()}) + Q.all_set_genes(glob)
    prim = {}
    ibd, _ = D.ibd(args.subsample, args.seed)
    prim["ibd"] = (*Q.gated_cohorts(ibd, pd.read_csv(out / "a55_ibd_per_series.csv")), S["colon"])
    pc = cancer_table(run["quick"])
    if pc is not None:
        t, w = Q.tcga_table(pc, genes, args.subsample, args.seed)
        sets_of = lambda c: S[Q.CANCER[c]] if Q.CANCER.get(c) else glob
        prim["tcga"] = (t, w, sets_of)
    res = {k: compare(s, w, wf, k, n_perm, out, prefix="a84") for k, (s, w, wf) in prim.items()}
    holm_summary(res, out, "a84", "H84")
    if pc is not None:
        for c in w:
            compare(t[t.cohort == c], {c: w[c]}, sets_of(c), f"tcga_{c}", min(1000, n_perm), out, prefix="a84")
    samples, w6 = H.gated_samples(out6, args.subsample, args.seed)
    compare(samples, w6, S["skin"], "psoriasis", min(2000, n_perm), out, prefix="a84")
    disc, _ = D.discovery_all(args.subsample, args.seed)
    ad = disc[disc.disease == "atopic_dermatitis"]
    compare(*Q.gated_cohorts(ad, pd.read_csv(out / "a53_skin_per_series.csv").query("disease == 'atopic_dermatitis'")),
            S["skin"], "atopic_dermatitis", min(2000, n_perm), out, prefix="a84")


if __name__ == "__main__":
    main()
