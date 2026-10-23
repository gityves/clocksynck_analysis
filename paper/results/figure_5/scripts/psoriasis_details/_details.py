"""Helpers of psoriasis_details: sample tables (all gene columns) of the other paired skin and gut series,
and the A53 per-cohort damping pipeline built from the public building blocks of `clocksync.paired`.

The A54 pipeline and the table helpers are `psoriasis/_shared.py` (imported as `_shared`).
"""

from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(Path(__file__).resolve().parents[4]))  # paper/
from common.paths import group_dir  # noqa: E402
sys.path.insert(0, str(group_dir("psoriasis")))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import _shared as H  # noqa: E402

FOLDER = "psoriasis_details"
SRC_SKIN = H.DATA / "skin_lesional"
SRC_IBD = H.DATA / "ibd_bulk" / "processed"
PD, S, GENES, WF = H.PD, H.S, H.GENES, H.WF


def thin_cohorts(samples, n, seed):
    return H.thin(samples, n, seed)


def discovery_all(subsample=None, seed=0):
    """Every sample of the 17 discovery series (psoriasis and AD), all gene columns, cohort `series|disease`, and the
    expression-matched means per cohort (builder `skin_paired_discovery`)."""
    b = H.builder("skin_paired_discovery")
    samples = H.lesional_all_genes(b.SRC)
    return H.thin(samples, subsample, seed), b.gene_means()


def other_skin(subsample=None, seed=0):
    """Series of `processed_other` (standard scale) and `processed_other_linear`, all genes; column `source`."""
    b = H.builder("skin_paired_other")
    parts = []
    for tag, folder in b.FOLDERS.items():
        t = H.lesional_all_genes(b.ROOT / folder)
        t.insert(3, "source", tag)
        parts.append(t)
    return H.thin(pd.concat(parts, ignore_index=True), subsample, seed), b.gene_means()


def ibd(subsample=None, seed=0):
    b = H.builder("gut_bulk_ibd")
    samples = H.lesional_all_genes(b.SRC, by_disease=False)
    return H.thin(samples, subsample, seed), b.gene_means()


def pairs_of(samples: pd.DataFrame, cohort) -> tuple:
    case, ref = H.split_pairs(samples[samples.cohort == cohort])
    return case, ref, H.is_counts(samples[samples.cohort == cohort])


def p_lower(null, obs) -> float:
    return float(S.perm_p(null, obs, "less", tie=0.0))


def damping(samples: pd.DataFrame, means: dict, seed: int, n_perm: int, n_boot: int, n_pseudo: int = 200,
            min_pairs: int = 5, diseases=None):
    """A53 (H53.2): per cohort the reference gate, r with bootstrap se, cross-moment bounds and the expression-matched
    pseudo-clock r; pooled by disease over gated cohorts with the lower-tail swap p and pseudo q05 (the A53 hypothesis:
    r < 0). One random stream over the cohorts in file x disease order (the legacy order), as `clocksync.paired` does.

    local: not in clocksync (`paired` pools all cohorts together and tests r > 0 only): the per-cohort loop is the body
    of `paired` on the public blocks; the lower-tail tests and the pseudo q05 are local.
    """
    rng = np.random.default_rng(seed)
    rows, data, pseudo = [], {}, {}
    for series in sorted(samples.series.unique()):
        s = samples[samples.series == series]
        for dis in sorted(s.disease.astype(str).unique()):
            key = f"{series}|{dis}"
            if key not in set(s.cohort):
                continue
            case, ref, counts = pairs_of(s, key)
            if len(case) < min_pairs:
                continue
            yn = H.logexpr(ref, GENES, counts)
            g = PD.reference_gate(yn, WF, 200, rng)
            b = PD.cohort_states(case, ref, counts, WF)
            r, se = PD.log_var_ratio(*b), PD.boot_se(*b, n_boot=n_boot, rng=rng)
            ps = np.array([PD.log_var_ratio(*PD.cohort_states(case, ref, counts, w))
                           for w in PD.pseudo_waveforms(means[key], n_pseudo, rng, WF)])
            data[key], pseudo[key] = b, ps
            rows.append(dict(cohort=key, series=series, disease=dis, n_pairs=len(case), counts=counts,
                             gate_share=g["share"], gate_null=g["share_null_mean"], gate_p=g["p"],
                             gated=g["p"] < 0.05, r=r, se=se, p_pseudo=p_lower(ps, r),
                             pseudo_median=float(np.nanmedian(ps)), **{f"x_{k}": v for k, v in PD.cross(*b).items()}))
    d = pd.DataFrame(rows)
    out = []
    for dis in (diseases or sorted(d.disease.unique())):
        g = d[(d.disease == dis) & d.gated & np.isfinite(d.r) & (d.se > 0)]
        if g.empty:
            continue
        keys, w = list(g.cohort), 1 / g.se.to_numpy() ** 2
        pool = S.meta(g.r, g.se)
        null = PD.swap_null([data[k] for k in keys], w, n_perm=n_perm, rng=rng)
        ps = pd.DataFrame(pseudo)[keys].to_numpy()
        ps_pool = np.nansum(ps * (w / w.sum()), axis=1)
        q05 = float(np.nanquantile(ps_pool, 0.05))
        p_swap = p_lower(null, pool["fixed"])
        out.append(dict(disease=dis, series=";".join(g.series), n_pairs=int(g.n_pairs.sum()), r=pool["fixed"],
                        se_fixed=pool["se_fixed"], r_re=pool["random"], lo_re=pool["lo"], hi_re=pool["hi"],
                        tau2=pool["tau2"], k=pool["k"], p_swap_lower=p_swap, pseudo_q05=q05,
                        pseudo_median=float(np.nanmedian(ps_pool)), p_pseudo=p_lower(ps_pool, pool["fixed"]),
                        H53_2="PASS" if (pool["fixed"] < 0 and p_swap < 0.05 and pool["fixed"] < q05) else "FAIL"))
    return d, pd.DataFrame(out)
