"""A52: is the clock-plane variance of TCGA tumours lower than that of the adjacent normal of the same patient?

Pre-registered in docs/PLAN_internal_desync.md (archive) (A52, secondaries). Units: patients of ten cancers with one primary
tumour and one solid-tissue normal. Estimand: r = log tr Cov(beta_tumour) - log tr Cov(beta_normal) (one projector
for both, Poisson-corrected), per cancer and pooled over the cancers whose normals pass the clock gate (inverse-variance
weights fixed from the data; random-effects interval reported).
Null: per-patient swaps of the tumour / normal labels (states centred within sample type first), same weights;
specificity: 200 expression-matched pseudo-clock gene sets carrying the waveform (H52.2: pooled clock r below the
5th percentile of the pooled pseudo-clock r).
Interval: patient bootstrap (1,000) per cancer.
Secondaries: pooled r without LIHC and KIRC; r by pathological stage (early I-II, late III-IV; >= 15 pairs per cancer,
all cancers pooled, descriptive); pooled tumour-minus-normal offset (circular mean, weights n x coherence).

    python compute_tcga_paired.py [--subsample N] [--seed S]
Outputs (paper/results/figure_5/output/cancer/): a52_per_cancer.csv, a52_pooled.csv, a52_pseudo.csv, a52_secondary.csv,
a52_stage.csv
Legacy: legacy/desync/validation/paired_damping/scripts/a52_tcga.py and a52_secondary.py.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.paired import cohort_states, log_var_ratio, pseudo_waveforms, swap_null  # noqa: E402
from common import cli  # noqa: E402
from common.paths import out_dir  # noqa: E402

import _shared as sh  # noqa: E402

FOLDER = "cancer"
SEED = 52


def stage_group(s: pd.Series) -> pd.Series:
    """Pathological stage -> early (I, II) / late (III, IV), as the legacy `a52_secondary.stage`."""
    roman = s.astype(str).str.extract(r"Stage (I{1,3}V?|IV)")[0]
    return roman.map({"I": "early", "II": "early", "III": "late", "IV": "late"})


def pseudo_r(cancer: str, units, n: int, rng) -> np.ndarray:
    """Pseudo-clock r of one cancer: expression-matched non-clock gene sets carrying the waveform's amplitudes and
    phases (local: `paired` takes the pseudo genes from the sample table, which holds the 14 clock genes only, so the
    genome-wide counts of the paired patients are read here; the match is on the mean log fraction of the normals)."""
    ct, cn = sh.tcga_full_counts(cancer, units)
    lm = np.log(cn.div(cn.sum(axis=1), axis=0).mean() + 1e-12)[cn.mean() >= 1]
    dt, dn = ct.sum(axis=1), cn.sum(axis=1)
    out = []
    for w in pseudo_waveforms(lm, n, rng, sh.WF):
        g = list(w.genes)
        out.append(log_var_ratio(*cohort_states(ct[g].assign(depth=dt), cn[g].assign(depth=dn), True, w)))
    return np.array(out)


def main(argv=None):
    args = cli.parser(__doc__).parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    n_perm, n_boot, n_pseudo = (50, 20, 20) if run["quick"] else (10_000, 1000, 200)
    samples = sh.load_tcga(args.subsample, args.seed)
    rng = np.random.default_rng(SEED + args.seed)

    res = cs.paired(samples, unit="unit", condition="condition", levels=sh.LEVELS, cohort="cohort",
                    counts="value_type", n_perm=n_perm, n_boot=n_boot, n_gate=200, seed=SEED + args.seed)
    d = res.cohorts.rename(columns={"cohort": "cancer", "gate_p": "gate_p", "x_corr": "corr",
                                    "x_offset_h": "offset_h_cross"})
    pseudo = {c: pseudo_r(c, sh.pair_frames(samples, c)[0].index, n_pseudo, rng) for c in d.cancer}
    d["p_pseudo"] = [cs.stats.perm_p(pseudo[c], r, "less", tie=0.0) for c, r in zip(d.cancer, d.r)]
    d["pseudo_median"] = [float(np.nanmedian(pseudo[c])) for c in d.cancer]
    d.to_csv(out / "a52_per_cancer.csv", index=False)
    pd.DataFrame(pseudo).to_csv(out / "a52_pseudo.csv", index=False)

    keys = res.pooled.get("cohorts", "").split(";") if res.pooled.get("n_cohorts", 0) else []
    g = d.set_index("cancer").loc[keys]
    w = 1 / g.se.to_numpy() ** 2
    null = swap_null([res.data[c] for c in keys], w, n_perm=n_perm, rng=rng)
    ps_pool = np.nansum(pd.DataFrame(pseudo)[keys].to_numpy() * (w / w.sum()), axis=1)
    r_pool = float(res.pooled["fixed"])
    q05 = float(np.nanquantile(ps_pool, 0.05))
    p_swap = cs.stats.perm_p(null, r_pool, "less", tie=0.0)
    pooled = pd.Series(dict(
        cancers=";".join(keys), k=len(keys), r_fixed=r_pool, se_fixed=res.pooled["se_fixed"],
        r_random=res.pooled["random"], lo_random=res.pooled["lo"], hi_random=res.pooled["hi"],
        tau2=res.pooled["tau2"], p_swap_lower=p_swap, p_swap_upper=res.pooled["p_swap_upper"],
        pseudo_pooled_q05=q05, pseudo_pooled_median=float(np.nanmedian(ps_pool)),
        p_pseudo=cs.stats.perm_p(ps_pool, r_pool, "less", tie=0.0),
        H52_1="PASS" if (r_pool < 0 and p_swap < 0.05) else "FAIL", H52_2="PASS" if r_pool < q05 else "FAIL",
        offset_h_pooled_swap=res.pooled["offset_h"], p_offset=res.pooled["p_offset"], n_swap=n_perm))
    pooled.to_csv(out / "a52_pooled.csv", header=False)
    print(d[["cancer", "n_pairs", "gate_p", "r", "se", "p_pseudo"]].round(3).to_string(index=False))
    print(pooled.to_string(), flush=True)

    # ---- secondaries ------------------------------------------------------------------------------------------
    rows = []
    sub = d[~d.cancer.isin(["LIHC", "KIRC"])]
    rows.append(dict(analysis="without_LIHC_KIRC", **cs.stats.meta(sub.r, sub.se)))
    st = samples.assign(stage_group=stage_group(samples.stage))
    per = []
    for lev in ("early", "late"):
        x = st[st.stage_group == lev]
        r = cs.paired(x, unit="unit", condition="condition", levels=sh.LEVELS, cohort="cohort", counts="value_type",
                      n_perm=1, n_boot=min(n_boot, 500), n_gate=1, gated_only=False, min_pairs=15 if not run["quick"] else 6,
                      seed=SEED + 1 + args.seed)
        if len(r.cohorts):
            per.append(r.cohorts.assign(stage=lev)[["cohort", "stage", "n_pairs", "r", "se"]].rename(
                columns={"cohort": "cancer", "n_pairs": "n"}))
            rows.append(dict(analysis=f"stage_{lev}", **cs.stats.meta(r.cohorts.r, r.cohorts.se)))
    pd.concat(per).to_csv(out / "a52_stage.csv", index=False)
    wt = d.n_pairs * d["corr"]
    m = np.sum(wt * np.exp(1j * d.offset_h_cross / sh.HOURS)) / np.sum(wt)
    rows.append(dict(analysis="offset_h_weighted_mean", fixed=float(np.angle(m) * sh.HOURS), k=len(d)))
    sec = pd.DataFrame(rows)
    sec.to_csv(out / "a52_secondary.csv", index=False)
    print(sec.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
