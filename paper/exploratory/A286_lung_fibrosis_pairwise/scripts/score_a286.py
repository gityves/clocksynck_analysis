"""A286: does fibrosis change how much clock state the lung's cell classes share? HLCA healthy vs fibrotic donors.

Registered in paper/docs/registered_experiments.md (A286, staged). Anchor-free pairwise readout (`clocksync.pairwise`):
per arm, pairs re-derived on the arm's rows; S1 = T_SV (median |r| of SV_X pairs), S2 = D = T_SV - T_LY, S3 = F =
median |r| of fibroblast-X minus endothelial-X pairs (X = non-SV classes); contrast = case minus control.

  power  arm labels permuted across donors within stratum (2000 draws, seeds [286, d]): null SD and MDE (2.49 SD) of
         each contrast, all controls, living controls only, and within study (donors of the strata holding both arms;
         power-stage amendment). No observed contrast is computed.

    python score_a286.py --stage power [--cohort hlca|natri] [--draws 2000] [--jobs 6]
Outputs (../output/, natri: ../output/natri/): design.csv, power.csv
"""

from __future__ import annotations

import importlib
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
sys.path.insert(0, str(PAPER))
sys.path.insert(0, str(PAPER / "exploratory/A276_pairwise_calibration/scripts"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

from clocksync import pairwise as pw  # noqa: E402
from common import cli  # noqa: E402

AMENDMENT = HERE.parents[1]
SEED, MDE_K = 286, 2.49
LIVING = ("surgical_resection", "biopsy", "biopsy_cryopreserved", "lung_explant")
SV = ("fibroblast_stromal", "endothelial", "mural")


def family_of():
    import calibrate_a276 as C
    return C.A150.FAM.family


COHORTS = {"hlca": ("hlca_fibrotic_genome", "hlca_fibrotic"), "natri": ("lung_natri_genome", "ild_natri")}


def load(cohort):
    ds, col = COHORTS[cohort]
    pb, d = importlib.import_module(f"common.datasets.{ds}.build").build()
    d = d.set_index("donor") if "donor" in d else d
    return pb.reset_index(drop=True), d.rename(columns={col: "arm"})


def arm_stats(z, key, fam, members) -> dict:
    """S1, T_LY, S2, S3 of one arm (rows of `members` donors)."""
    m = key.donor.isin(members).to_numpy()
    k = key[m].reset_index(drop=True)
    pairs = pw.pair_index(k)
    r = pw.abs_r(z[m], pairs)[0]
    fa = np.array([fam[p["a"]] for p in pairs])
    fb = np.array([fam[p["b"]] for p in pairs])
    sv_a, sv_b = np.isin(fa, SV), np.isin(fb, SV)
    ly_a, ly_b = fa == "lymphocyte", fb == "lymphocyte"
    sv_x = (sv_a & ~sv_b & ~ly_b) | (sv_b & ~sv_a & ~ly_a)
    ly_x = (ly_a & ~sv_b & ~ly_b) | (ly_b & ~sv_a & ~ly_a)
    fib_x = ((fa == "fibroblast_stromal") & ~sv_b) | ((fb == "fibroblast_stromal") & ~sv_a)
    end_x = ((fa == "endothelial") & ~sv_b) | ((fb == "endothelial") & ~sv_a)
    med = lambda s: float(np.nanmedian(r[s])) if s.any() else np.nan  # noqa: E731
    t_sv, t_ly = med(sv_x), med(ly_x)
    return dict(S1=t_sv, T_LY=t_ly, S2=t_sv - t_ly, S3=med(fib_x) - med(end_x),
                n_pairs=len(pairs), n_SV_X=int(sv_x.sum()), n_LY_X=int(ly_x.sum()), n_fib_X=int(fib_x.sum()),
                n_end_X=int(end_x.sum()), donors=int(k.donor.nunique()))


def permute(lab: pd.Series, strata: pd.Series, rng) -> pd.Series:
    """Arm labels permuted across donors within stratum (strata with one arm keep their labels)."""
    out = lab.copy()
    for _, idx in strata.groupby(strata).groups.items():
        v = lab[idx]
        if v.nunique() > 1:
            out[idx] = rng.permutation(v.to_numpy())
    return out


def contrast(z, key, fam, lab, d):
    c = arm_stats(z, key, fam, lab.index[lab == "case"])
    k = arm_stats(z, key, fam, lab.index[lab == "control"])
    return {s: c[s] - k[s] for s in ("S1", "S2", "S3")}


def main():
    ap = cli.parser(__doc__)
    ap.add_argument("--stage", choices=("power",), required=True)
    ap.add_argument("--draws", type=int, default=2000)
    ap.add_argument("--cohort", choices=sorted(COHORTS), default="hlca")
    args = ap.parse_args()
    cli.setup(args)
    out = AMENDMENT / "output" / ("" if args.cohort == "hlca" else args.cohort)
    out.mkdir(parents=True, exist_ok=True)
    pb, d = load(args.cohort)
    family = family_of()
    fam = {c: family(c) for c in pb.compartment.unique()}
    z = pw.states(pb)
    key = pb[["donor", "compartment", "stratum"]]
    strata = pb.drop_duplicates("donor").set_index("donor").stratum
    lab0 = d.arm.dropna()
    lab0 = lab0[lab0.index.isin(strata.index)]
    both = lab0.groupby(strata.reindex(lab0.index)).nunique()
    within = lab0[strata.reindex(lab0.index).isin(both.index[both > 1]).to_numpy()]
    designs = {"all controls": lab0, "within study (strata with both arms)": within}
    if args.cohort == "hlca":
        designs["living controls"] = d.arm[(d.arm == "case") | d.tissue_sampling_method_major.isin(LIVING)].dropna()  # power-stage amendment: the only unconfounded design
    rows, power = [], []
    for name, lab in designs.items():
        lab = lab[lab.index.isin(strata.index)]
        st = strata.reindex(lab.index)
        for arm in ("control", "case"):
            rows.append(dict(design=name, arm=arm, **{k: v for k, v in arm_stats(z, key, fam, lab.index[lab == arm]).items()
                                                       if k.startswith("n_") or k == "donors"}))
        mixed = st.groupby(st).apply(lambda s: lab[s.index].nunique() > 1)
        res = Parallel(n_jobs=args.jobs)(delayed(contrast)(z, key, fam, permute(lab, st, np.random.default_rng([SEED, d_])), d)
                                         for d_ in range(args.draws))
        for s in ("S1", "S2", "S3"):
            v = np.array([r[s] for r in res], float)
            power.append(dict(design=name, statistic=s, draws=int(np.isfinite(v).sum()), null_mean=float(np.nanmean(v)),
                              null_sd=float(np.nanstd(v)), mde=MDE_K * float(np.nanstd(v)),
                              strata_with_both_arms=int(mixed.sum()), donors_in_those=int(st.isin(mixed.index[mixed]).sum())))
        print(name, "done", flush=True)
    D, P = pd.DataFrame(rows), pd.DataFrame(power)
    D.to_csv(out / "design.csv", index=False)
    P.to_csv(out / "power.csv", index=False)
    pd.set_option("display.width", 200)
    print(D.to_string(), P.round(4).to_string(), sep="\n")


if __name__ == "__main__":
    main()
