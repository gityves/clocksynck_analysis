"""A17 step 6: waveform sensitivity of the healthy fits (descriptive).

Pre-registered in docs/PLAN_internal_desync.md (archive) (A17 step 6: "Gate and anchored couplings are rerun with the tissue's own GTEx
waveform where GTEx has one"; descriptive, no hypothesis id). Each healthy fit of `atlas_fits` (learned anchors, frozen table) is
run with the global human waveform and with the GTEx waveform of its tissue (`Dataset(waveform="<GTEx tissue>")`): gate (500
permutations, free host as in the fits) and the anchored couplings (point estimates, no bootstrap). Estimand: change of the gate verdict, the
correlation of the couplings (all compartments; per fit) and their median absolute difference. GTEx has no ileum entry: ileum fits are not rerun.

    python compute_waveform_sensitivity.py [--subsample N] [--jobs 4] [--only NAME ...] [--n-perm P]
Outputs: paper/results/extended_data/output/gate/a17_waveform_sens_gates.csv, a17_waveform_sens.csv, a17_waveform_sens_summary.csv
Legacy: legacy/desync/common/waveform_sens.py -> common_output/a17_waveform_sens{,_gates}.csv
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))  # paper/
from common.paths import group_dir  # noqa: E402
sys.path.insert(0, str(group_dir("atlas_fits")))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

from _shared import Target  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import cli  # noqa: E402
from common.paths import out_dir  # noqa: E402
from compute_atlas_fits import anchor_sets  # noqa: E402

FOLDER = "gate"
WF = {"psychad_hbcc": "Brain - Frontal Cortex (BA9)", "psychad_mssm": "Brain - Frontal Cortex (BA9)",
      "seaad_dlpfc": "Brain - Frontal Cortex (BA9)", "seaad_mtg": "Brain - Cortex", "pd_control": "Brain",
      "lung_hlca": "Lung", "lung_natri_control": "Lung", "ia_lung": "Lung", "heart": "Heart - Left Ventricle",
      "kidney_sc": "Kidney - Cortex", "organ_breast": "Breast - Mammary Tissue",
      "organ_adipose": "Adipose - Subcutaneous", "organ_muscle": "Muscle - Skeletal",
      "organ_colon": "Colon - Transverse", "kong_colon": "Colon - Transverse", "ibd_rectum": "Colon - Sigmoid",
      "skin_ad_control": "Skin - Not Sun Exposed (Suprapubic)", "skin_ssc_control": "Skin - Not Sun Exposed (Suprapubic)",
      "islets": "Pancreas", "endometrium": "Uterus"}


def one(name, subsample, seed, n_perm, i):
    limit_threads(1)
    warnings.filterwarnings("ignore")
    try:
        T = Target(name, subsample, seed)
        anc = [a for a in anchor_sets(name)["learned"] or [] if a in set(T.pb.compartment)]
        gate, rows = {}, []
        for label, wf in (("tissue", WF[name]), ("global", "human")):
            g = T.dataset(anchors=None, waveform=wf).gate(n_perm=n_perm, seed=i)
            gate.update({f"p_pairing_{label}": g.p_pairing, f"p_acrophase_{label}": g.p_acrophase, f"passed_{label}": bool(g.passes)})
            if len(anc) >= 3 and g.passes:
                t = T.dataset(anchors=anc, waveform=wf).fit(n_boot=0, gate=False).table
                rows.append(t[["coupling", "own_amplitude", "offset_h", "psi", "V", "n", "anchor", "boundary"]].reset_index()
                            .assign(waveform=label))
        gate = dict(target=name, waveform=WF[name], **gate)
        if len(rows) == 2:
            a, b = (r.set_index("compartment") for r in rows)
            cmp_ = a[["coupling"]].join(b[["coupling"]], lsuffix="_tissue", rsuffix="_global", how="inner").reset_index().assign(target=name)
            return gate, cmp_
        return gate, None
    except Exception as e:
        return dict(target=name, error=f"{type(e).__name__}: {e}"), None


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--only", nargs="*", default=None)
    ap.add_argument("--n-perm", type=int, default=None)
    args = ap.parse_args(argv)
    run = cli.setup(args)
    n_perm = args.n_perm or (100 if run["quick"] else 500)
    out = out_dir(FOLDER, run["quick"])
    names = [n for n in WF if args.only is None or n in args.only]
    res = Parallel(n_jobs=args.jobs)(delayed(one)(n, args.subsample, args.seed, n_perm, i) for i, n in enumerate(names))
    gates = pd.DataFrame([r[0] for r in res])
    gates.to_csv(out / "a17_waveform_sens_gates.csv", index=False)
    cp = pd.concat([r[1] for r in res if r[1] is not None], ignore_index=True)
    cp.to_csv(out / "a17_waveform_sens.csv", index=False)
    d = cp.coupling_tissue - cp.coupling_global
    per = cp.groupby("target").apply(lambda x: x.coupling_tissue.corr(x.coupling_global), include_groups=False)
    ok = gates.get("passed_tissue", pd.Series(dtype=bool))
    summ = pd.DataFrame([dict(fits=len(gates), verdict_unchanged=int((gates.passed_tissue == gates.passed_global).sum()),
                              passes_tissue=int(gates.passed_tissue.sum()), passes_global=int(gates.passed_global.sum()),
                              compartments=len(cp), r=cp.coupling_tissue.corr(cp.coupling_global), median_abs_diff=float(d.abs().median()),
                              min_r_per_fit=float(per.min()), fit_min_r=per.idxmin())])
    summ.to_csv(out / "a17_waveform_sens_summary.csv", index=False)
    print(gates.round(3).to_string(index=False)); print(summ.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
