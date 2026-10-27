"""A141 steps 1 and 2: coverage of the cycle by the brain hosts and absolute amplitude against GTEx brain.

Step 1 (needs the output of step3_simulate.py for the uniform-sampling reference). Per cohort (dataset; controls = the
control set of the contrast with the most controls, A135; "all donors" = the labelled donors of that contrast):
Rayleigh R of the host phases, host spread V and eccentricity, for the A17-anchor host and the all-compartment host.
Reference: (i) GTEx fast-death brain tissues (Hardy 1-2, timed, any hour): host on all brain tissues with >= 40 donors,
aggregates only; (ii) the uniform-sampling simulations of step 3 (same donors, fitted amplitudes and noise): 2.5 / 97.5
percentiles of R. "Full cycle": R inside that 95 % range.
Step 2: per compartment, absolute coupled amplitude sqrt(z2 V) (host-referenced amplitude times host amplitude, in the
clock-plane units of the timed GTEx fit) on the control donors with the A17 host, against the range of `abs_amp` of the
GTEx brain regions (tissues.csv, cohort fast). Reported, not tested.
GTEx time of death is controlled access: only aggregates are written, no per-donor phase.

    python step12_coverage_amplitude.py [--subsample N]
Output: output/step1_coverage.csv, step1_gtex.csv, step2_amplitude.csv, step2_summary.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _a141 as H  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402

OUT = H.OUT


def gtex_brain() -> pd.DataFrame:
    org = H._mod("organs_shared", H.PAPER / "results/figure_2/scripts/organs/_shared.py")
    pb, _ = org.cohort("fast")
    brain = sorted(k for k in pb.compartment.unique() if k.startswith("Brain -"))
    pb = pb[pb.compartment.isin(brain)]
    ds = org.dataset(pb, anchors=brain)
    rows = []
    for label, tis in (("all brain tissues", brain),
                       ("cortical (Cortex, Frontal Cortex BA9, ACC BA24)",
                        [k for k in brain if k in ("Brain - Cortex", "Brain - Frontal Cortex (BA9)",
                                                   "Brain - Anterior cingulate cortex (BA24)")])):
        st = H.host_stats(ds.data, tis)
        rows.append(dict(reference=f"GTEx fast, {label}", n_tissues=len(tis), **st))
    return pd.DataFrame(rows)


def main():
    limit_threads(1)
    obs = pd.read_csv(OUT / "step3" / "observed.csv").set_index("contrast")
    sims = pd.read_csv(OUT / "step3" / "sims.csv")
    g = gtex_brain()
    g.to_csv(OUT / "step1_gtex.csv", index=False)
    print(g.round(3).to_string(), flush=True)
    rows = []
    for coh, (_, ctl) in H.COHORTS.items():
        for hn in ("a17", "allcomp"):
            for sn in ("control", "all"):
                col = lambda s: f"{hn}_{sn}_{s}"  # noqa: E731
                sm = sims[sims.contrast == ctl]
                lo, med, hi = np.quantile(sm[col("R")], [0.025, 0.5, 0.975])
                R = obs.loc[ctl, col("R")]
                rows.append(dict(cohort=coh, contrast_set=ctl, host=hn, donors=sn, n=int(obs.loc[ctl, col("n")]), R=R,
                                 p_rayleigh=obs.loc[ctl, col("p_rayleigh")], V=obs.loc[ctl, col("V")],
                                 ecc=obs.loc[ctl, col("ecc")], sim_R_lo=lo, sim_R_median=med, sim_R_hi=hi,
                                 sim_V_median=sm[col("V")].median(), sim_ecc_median=sm[col("ecc")].median(),
                                 n_sim=len(sm), full_cycle=bool(lo <= R <= hi), R_above_range=bool(R > hi)))
    c = pd.DataFrame(rows)
    c.to_csv(OUT / "step1_coverage.csv", index=False)
    print(c.round(3).to_string(), flush=True)

    # step 2
    t = pd.read_csv(H.TISSUES)
    t = t[(t.cohort == "fast") & t.tissue.str.startswith("Brain")]
    lo, med, hi = t.abs_amp.min(), t.abs_amp.median(), t.abs_amp.max()
    rows = []
    for coh, (_, ctl) in H.COHORTS.items():
        C = H.Contrast(ctl)
        ctrl = C.lab.index[C.lab == "control"]
        par = H.fitted_params(C.ds.subset(donors=ctrl))
        par = par.reindex([k for k in C.targets if k in par.index])
        for k, r in par.iterrows():
            rows.append(dict(cohort=coh, compartment=k, n=r.n, anchor=r.anchor, log_amp_rel=r.log_amp_rel, host_amplitude=np.sqrt(r.V),
                             log_amp_abs=np.log(r.amp) if r.amp > 0 else np.nan, abs_amp=r.amp, coupling=r.coupling, psi=r.psi,
                             gtex_min=lo, gtex_median=med, gtex_max=hi, ratio_to_gtex_median=r.amp / med,
                             within_gtex_range=bool(lo <= r.amp <= hi)))
    a = pd.DataFrame(rows)
    a.to_csv(OUT / "step2_amplitude.csv", index=False)
    s = a.groupby("cohort").agg(n_comp=("abs_amp", "size"), median_abs_amp=("abs_amp", "median"),
                                min_abs_amp=("abs_amp", "min"), max_abs_amp=("abs_amp", "max"),
                                n_within=("within_gtex_range", "sum"), n_above=("abs_amp", lambda x: int((x > hi).sum())),
                                n_below=("abs_amp", lambda x: int((x < lo).sum())))
    s["gtex_min"], s["gtex_median"], s["gtex_max"], s["n_gtex_regions"] = lo, med, hi, len(t)
    s["median_ratio_to_gtex_median"] = s.median_abs_amp / med
    allv = a.abs_amp
    s.loc["all cohorts"] = [len(a), allv.median(), allv.min(), allv.max(), int(a.within_gtex_range.sum()),
                            int((allv > hi).sum()), int((allv < lo).sum()), lo, med, hi, len(t), allv.median() / med]
    s.to_csv(OUT / "step2_summary.csv")
    print(s.round(3).to_string(), flush=True)


if __name__ == "__main__":
    main()
