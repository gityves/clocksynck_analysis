"""A69: how much of the psoriatic delay does the keratinocyte state shift explain, with healthy-skin phases? (descriptive)

Pre-registered in docs/PLAN_internal_desync.md (archive) (A69). Bin shares of keratinocyte mRNA and bin clock-gene CPM per status
from the psoriasis units with bins in both statuses (builder `psoriasis_sc_bins`, A65 data; kc_cycling, kc_d1..kc_d4); bin
states z (complex, in the free-host frame) from the A63 healthy fit (`healthy_skin_bins`; sensitivity: healthy plus
uninvolved psoriatic donors). Composition-only prediction (A60 mixing, `_attribution.py`): arg(sum_k g_k^L z_k / sum_k
g_k^N z_k) with the same z in both statuses, gains from the bulk projector (A54 bulk series at 30 M reads). Donor bootstrap of
the A63 fit (500); weights from all units or without GSE278330.

    python compute_bin_mix.py [--subsample N] [--n-boot 500] [--jobs 4]
Outputs (paper/results/figure_5/output/psoriasis_details/): a69_bin_mix.csv, a69_bin_shares.csv
Legacy: legacy/desync/validation/paired_damping/scripts/a69_bin_mix.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import _attribution as M  # noqa: E402
import _details as D  # noqa: E402
import _shared as H  # noqa: E402
import _singlecell as C  # noqa: E402
from compute_differentiation import KC, healthy_table  # noqa: E402
from common import cli  # noqa: E402

SEED = 69


def weights(pb, exclude=()):
    d = pb[pb.compartment.isin(KC) & ~pb.stratum.isin(exclude)]
    out, shares = {}, {}
    for st, s in C.STATUS.items():
        x = d[d.status == st].groupby("compartment")
        dep = x.depth.sum().reindex(KC)
        cpm = x[H.GENES].sum().reindex(KC).div(dep, axis=0).T * 1e6
        shares[s] = dep / dep.sum()
        out[s] = (shares[s], cpm)
    return out, pd.DataFrame(shares)


def zs(t) -> dict:
    tab = C.fit_states(t)
    return {k: tab.rel_amplitude[k] * np.exp(1j * tab.offset_h[k] / H.HOURS) for k in KC if k in tab.index}


def predict(z, g) -> float:
    ks = [k for k in KC if k in z]
    num = sum(g["L"][k] * z[k] for k in ks)
    den = sum(g["N"][k] * z[k] for k in ks)
    return float(np.angle(num / den) * H.HOURS)


def main(argv=None):
    args = H.parser(__doc__).parse_args(argv)
    run = cli.setup(args)
    out = H.out_dir(D.FOLDER, run["quick"])
    n_boot = C.sizes(args, run, 0, 500)["n_boot"]
    bins, _ = C.load("psoriasis_sc_bins", args.subsample, args.seed)
    healthy, _ = C.load("healthy_skin_bins", args.subsample, args.seed)
    L = M.bulk_projector(M.bulk())
    rows, share_tab = [], []
    for wlabel, excl in (("all_units", ()), ("no_GSE278330", ("GSE278330",))):
        w, sh = weights(bins, excl)
        share_tab.append(sh.assign(weights=wlabel))
        g = {s: M.gains(M.shares(w[s][0], w[s][1]), L) for s in "LN"}
        for flabel, status in (("healthy", ["healthy"]), ("all_donors", ["healthy", "nonlesional"])):
            t = healthy_table(healthy, status)
            v = predict(zs(t), g)
            seeds = np.random.default_rng(SEED).integers(0, 2**31, n_boot)
            bv = np.array([x for x in Parallel(n_jobs=args.jobs)(
                delayed(C.safe)(lambda s: predict(zs(C.boot(t, np.random.default_rng(s))), g), s) for s in seeds)
                if x is not None])
            lo, hi = C.interval(bv, v)
            rows.append(dict(weights=wlabel, phases=flabel, composition_only_h=v, lo=lo, hi=hi,
                             frac_negative=float(np.mean(bv < 0)), n_boot=len(bv)))
            print(rows[-1], flush=True)
    pd.DataFrame(rows).to_csv(out / "a69_bin_mix.csv", index=False)
    pd.concat(share_tab).to_csv(out / "a69_bin_shares.csv")
    print(pd.concat(share_tab).round(3).to_string())


if __name__ == "__main__":
    main()
