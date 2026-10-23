"""A58: which cells carry the psoriatic lesion delay? Status-split host on single-cell pseudobulks.

Pre-registered in docs/PLAN_internal_desync.md (archive) (A58 H58.1, H58.2). Whole-skin psoriasis single-cell and single-nucleus
data (GSE278330, GSE228421, GSE173706; builder `psoriasis_sc`): unit = study | patient | visit, compartments
`L:<class>` (lesional) and `N:<class>` (uninvolved), pseudobulks >= 20 cells, classes present in both statuses in >= 5
units, one stratum per study; the host is fitted on all compartments (free host, `Dataset(anchors=None)`).
Estimand (H58.1): the circular mean over classes of offset(L:k) - offset(N:k) (h) < 0. Null: statuses swapped per unit
(n_perm, one-sided lower p < 0.05); the gate must pass (a failed gate is uninformative). H58.2 (descriptive): keratinocyte
offset minus the stromal-vascular mean (fibroblast, endothelial, mural) and minus the immune mean (T, myeloid), per-class
offsets and log amplitude ratios with unit-bootstrap percentile intervals (offsets wrapped around the estimate).

    python compute_sc_status.py [--subsample N] [--n-perm 2000] [--n-boot 500] [--jobs 4]
Outputs: paper/results/figure_5/output/psoriasis/a58_{counts,gate,classes,verdict}.csv
Legacy: legacy/desync/validation/paired_damping/scripts/a58_sc.py (+ build_psoriasis_sc.py for the table)
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import _shared as H  # noqa: E402
import _singlecell as C  # noqa: E402
from common import cli  # noqa: E402

SEED = 58
SV, IMM = ("fibroblast", "endothelial", "mural"), ("T", "myeloid")


def stat(t, ks) -> dict:
    tab = C.fit_states(t[t.compartment.str[2:].isin(ks)].reset_index(drop=True))
    off = {k: float(C.wrap(tab.offset_h[f"L:{k}"] - tab.offset_h[f"N:{k}"])) for k in ks}
    amp = {k: float(tab.log_amp_rel[f"L:{k}"] - tab.log_amp_rel[f"N:{k}"]) for k in ks}
    out = dict(stat=C.cmean(list(off.values())), off=off, amp=amp)
    ker = off.get("keratinocyte")
    if ker is not None:
        sv = [off[k] for k in SV if k in off]
        im = [off[k] for k in IMM if k in off]
        out["ker_minus_sv"] = float(C.wrap(ker - C.cmean(sv))) if sv else np.nan
        out["ker_minus_imm"] = float(C.wrap(ker - C.cmean(im))) if im else np.nan
    return out


def main(argv=None):
    args = H.parser(__doc__).parse_args(argv)
    run = cli.setup(args)
    out = H.out_dir(H.FOLDER, run["quick"])
    n = C.sizes(args, run, 2000, 500)
    pb, _ = C.load("psoriasis_sc", args.subsample, args.seed)
    t = C.status_split(pb)
    ks = C.classes(t)
    t = t[t.compartment.str[2:].isin(ks)].reset_index(drop=True)
    t.groupby("compartment").donor.nunique().to_csv(out / "a58_counts.csv")
    gate = C.gate(t, SEED)
    pd.Series(gate).to_csv(out / "a58_gate.csv", header=False)
    obs = stat(t, ks)
    print("units", t.donor.nunique(), "classes", ks, "\n", gate, "\n", {k: obs[k] for k in ("stat", "ker_minus_sv")},
          flush=True)
    seeds = np.random.default_rng(SEED).integers(0, 2**31, n["n_perm"] + n["n_boot"])
    par = Parallel(n_jobs=args.jobs)
    null = par(delayed(C.safe)(lambda s: stat(C.swapped(t, np.random.default_rng(s)), ks)["stat"], s)
               for s in seeds[:n["n_perm"]])
    null = np.array([x for x in null if x is not None], float)
    bts = [b for b in par(delayed(C.safe)(lambda s: stat(C.boot(t, np.random.default_rng(s)), ks), s)
                          for s in seeds[n["n_perm"]:]) if b]
    rows = []
    for k in ks:
        lo, hi = C.interval([b["off"][k] for b in bts], obs["off"][k])
        alo, ahi = C.interval([b["amp"][k] for b in bts], obs["amp"][k], circular=False)
        rows.append(dict(cls=k, offset_h=obs["off"][k], off_lo=lo, off_hi=hi, log_amp=obs["amp"][k], amp_lo=alo,
                         amp_hi=ahi))
    cls = pd.DataFrame(rows)
    cls.to_csv(out / "a58_classes.csv", index=False)
    p = C.p_lower(null, obs["stat"])
    v = dict(units=t.donor.nunique(), studies=t.stratum.nunique(), classes=";".join(ks), gate_pass=gate["gate_pass"],
             stat_h=obs["stat"], p_swap=p, n_null=len(null), n_boot=len(bts),
             ker_minus_sv_h=obs.get("ker_minus_sv"), ker_minus_imm_h=obs.get("ker_minus_imm"),
             H58_1=("UNINFORMATIVE (gate)" if not gate["gate_pass"] else
                    "PASS" if (obs["stat"] < 0 and p < 0.05) else "FAIL"))
    for key in ("ker_minus_sv", "ker_minus_imm"):
        if obs.get(key) is not None:
            v[f"{key}_lo"], v[f"{key}_hi"] = C.interval([b[key] for b in bts if b.get(key) is not None], obs[key])
    pd.Series(v).to_csv(out / "a58_verdict.csv", header=False)
    print(cls.round(2).to_string(), "\n", pd.Series(v).to_string())


if __name__ == "__main__":
    main()
