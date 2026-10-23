"""A57, A74 (and the descriptive treatment part of A73): does the change of the psoriatic lesion delay follow clinical
response, and proliferation beyond response?

Pre-registered in docs/PLAN_internal_desync.md (archive) (A57 H57.1, A73 descriptive, A74 H74.1). Estimand: the change of the
lesional-vs-uninvolved offset (baseline minus post, h, wrapped; < 0: the delay resolves) in one patient group minus the
other. Null: group labels permuted across patients (n_perm), one-sided p (lower).
A57: GSE117468, all arms, responders (PASI improvement >= 75 %) vs non-responders (primary); median split of the PASI
improvement; placebo arm (>= 50 %, descriptive); GSE201827 split at the median reduction of the lesional IL-17 score.
A73 (treatment, descriptive): above vs below the median reduction of the lesional proliferation score.
A74: within PASI75 responders (and all patients), above vs below the median reduction of proliferation (H74.1: high -
low < 0, p < 0.05), IL-17, IFN-gamma, and differentiation restoration; Spearman of the reduction with PASI improvement.

    python compute_response.py [--subsample N]
Outputs: paper/results/figure_5/output/psoriasis/a57_response.csv, a73_treatment.csv, a74_response_prolif.csv
Legacy: legacy/desync/validation/paired_damping/scripts/a57_response.py, a74_response_prolif.py, a73_proliferation.py
        (treatment part)
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from scipy.stats import spearmanr  # noqa: E402
import _shared as H  # noqa: E402
from common import cli  # noqa: E402
from compute_treatment import thin_visits  # noqa: E402


def change(b, p, m) -> float:
    return H.wrap1(H.offset_h(H.sub(b, m)) - H.offset_h(H.sub(p, m)))


def test(b, p, grp, rng, label, n_perm):
    def stat(g):
        return H.wrap1(change(b, p, g) - change(b, p, ~g))

    obs = stat(grp)
    null = np.array([stat(rng.permutation(grp)) for _ in range(n_perm)])
    return dict(set=label, n_resp=int(grp.sum()), n_nonresp=int((~grp).sum()), change_resp_h=change(b, p, grp),
                change_nonresp_h=change(b, p, ~grp), resp_minus_nonresp_h=obs,
                p_lower=float((1 + np.sum(null <= obs)) / (1 + n_perm)))


def mean_score(x, genes, counts):
    g = [v for v in genes if v in x.columns]
    return H.logexpr(x, g, counts).mean(1)


def get(file, treated, n, seed):
    t, counts = H.visits(file, treated, all_genes=True)
    return thin_visits(t, n, seed), counts


def main(argv=None):
    args = H.parser(__doc__).parse_args(argv)
    run = cli.setup(args)
    out = H.out_dir(H.FOLDER, run["quick"])
    n_perm = H.sizes(args, run, 2000)["n_perm"]
    n, seed = args.subsample, args.seed

    # A57
    rng = np.random.default_rng(57)
    rows = []
    t, counts = get("GSE117468.parquet", None, n, seed)
    lb, lp = t[("baseline", "lesional")], t[("post", "lesional")]
    imp = 1 - lp.pasi.to_numpy(float) / lb.pasi.to_numpy(float)
    b, p = H.visit_states(t, counts)
    r = test(b, p, imp >= 0.75, rng, "GSE117468_PASI75_primary", n_perm)
    r["H57_1"] = "PASS" if (r["resp_minus_nonresp_h"] < 0 and r["p_lower"] < 0.05) else "FAIL"
    rows.append(r)
    rows.append(test(b, p, imp > np.median(imp), rng, "GSE117468_median_split", n_perm))
    pl = lb.arm.astype(str).str.contains("Placebo").to_numpy()
    bp, pp = H.visit_states({k: v[pl] for k, v in t.items()}, counts)
    if (imp[pl] >= 0.5).sum() >= 3:
        rows.append(test(bp, pp, imp[pl] >= 0.5, rng, "GSE117468_placebo_PASI50", n_perm))
    t2, c2 = get("GSE201827.parquet", None, n, seed)
    red = (mean_score(t2[("baseline", "lesional")], H.IL17, c2) - mean_score(t2[("post", "lesional")], H.IL17, c2))
    b2, p2 = H.visit_states(t2, c2)
    rows.append(test(b2, p2, red > np.median(red), rng, "GSE201827_IL17_reduction", n_perm))
    d = pd.DataFrame(rows)
    d.to_csv(out / "a57_response.csv", index=False)
    print(d.round(3).T.to_string(), flush=True)

    # A73 treatment (descriptive): a fresh stream, as the legacy script
    rng = np.random.default_rng(73)
    rows = []
    for file, treated, label in (("GSE117468.parquet", None, "GSE117468_all"),
                                 ("GSE117468.parquet", True, "GSE117468_treated"),
                                 ("GSE201827.parquet", None, "GSE201827_all")):
        t, counts = get(file, treated, n, seed)
        a, z = t[("baseline", "lesional")], t[("post", "lesional")]
        red = mean_score(a, H.PROLIF, counts) - mean_score(z, H.PROLIF, counts)
        b, p = H.visit_states(t, counts)
        r = test(b, p, red > np.median(red), rng, label + "_prolif_reduction", n_perm)
        rows.append(dict(score=label, n_pairs=r["n_resp"] + r["n_nonresp"], high_minus_low_h=r["resp_minus_nonresp_h"],
                         offset_high_h=r["change_resp_h"], offset_low_h=r["change_nonresp_h"], p_lower=r["p_lower"]))
    pd.DataFrame(rows).to_csv(out / "a73_treatment.csv", index=False)
    print(pd.DataFrame(rows).round(3).T.to_string(), flush=True)

    # A74
    rng = np.random.default_rng(74)
    t, counts = get("GSE117468.parquet", None, n, seed)
    lb, lp = t[("baseline", "lesional")], t[("post", "lesional")]
    imp = 1 - lp.pasi.to_numpy(float) / lb.pasi.to_numpy(float)
    red = {name: mean_score(lb, g, counts) - mean_score(lp, g, counts)
           for name, g in (("proliferation", H.PROLIF), ("IL17", H.IL17), ("IFNG", H.IFNG))}
    red["differentiation_restoration"] = ((mean_score(lp, H.DIFF, counts) - mean_score(lp, H.BASAL, counts)) -
                                          (mean_score(lb, H.DIFF, counts) - mean_score(lb, H.BASAL, counts)))
    resp = imp >= 0.75
    rows = []
    for name, rr_all in red.items():
        for label, m in (("responders", resp), ("all", np.ones_like(resp))):
            bs, ps = H.visit_states({k: v[m] for k, v in t.items()}, counts)
            rr = rr_all[m]
            o = test(bs, ps, rr > np.median(rr), rng, f"{name}_{label}", n_perm)
            rows.append(dict(score=name, patients=label, n=int(m.sum()), change_high_h=o["change_resp_h"],
                             change_low_h=o["change_nonresp_h"], high_minus_low_h=o["resp_minus_nonresp_h"],
                             p_lower=o["p_lower"], spearman_with_pasi_improvement=spearmanr(rr_all, imp).statistic))
    d = pd.DataFrame(rows)
    pr = d[(d.score == "proliferation") & (d.patients == "responders")].iloc[0]
    d["H74_1"] = np.where((d.score == "proliferation") & (d.patients == "responders"),
                          "PASS" if (pr.high_minus_low_h < 0 and pr.p_lower < 0.05) else "FAIL", "")
    d.to_csv(out / "a74_response_prolif.csv", index=False)
    print(d.round(3).to_string(), flush=True)


if __name__ == "__main__":
    main()
