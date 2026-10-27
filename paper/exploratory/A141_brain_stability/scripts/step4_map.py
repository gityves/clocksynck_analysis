"""A141 step 4 (only for contrasts that are not noise-driven in step 3): which compartments agree across the two hosts.

Expectations fixed in the registration: agree = the neuron classes among themselves, and astrocytes with
oligodendrocytes; noise-driven = compartments with the lowest coupling or < 20 donors per group.
Choices not fixed there (set before scoring): agree = same sign of d on the A17 and the all-LOO host and
|d_A17 - d_LOO| <= the pooled bootstrap SE sqrt((se_A17^2 + se_LOO^2) / 2) of the WP0 table; "lowest coupling" = the
lowest third of the targets of the contrast by control coupling on the A17 host (`Dataset.contrast` c_control);
neuron classes = excitatory / glutamatergic / inhibitory_* / GABAergic / neuron_other (regions of `pd` kept
separate); astrocyte / oligodendrocyte likewise per region. Expected-agree set = neurons and astro/oligo; the
neuron and glia sets are reported separately. Every compartment is labelled by expectation (agree, noise, both, none);
the 2 x 2 table is expected-agree vs expected-noise against observed agree.

    python step4_map.py
Output: output/step4_map.csv, step4_summary.csv (only when step3_summary.csv has a contrast with noise_driven False)
"""

from __future__ import annotations

import re
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _a141 as H  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402

NEURON = re.compile(r"^brain:(?:(?:DMV|GPi|M1|PFC|V1)_)?(excitatory|glutamatergic|inhibitory.*|GABAergic|neuron_other)$")
GLIA = re.compile(r"^brain:(?:(?:DMV|GPi|M1|PFC|V1)_)?(astrocyte|oligodendrocyte)$")


def main():
    limit_threads(1)
    s3 = pd.read_csv(H.OUT / "step3" / "step3_summary.csv")
    todo = s3[~s3.noise_driven].contrast.tolist()
    print("not noise-driven:", todo, flush=True)
    if not todo:
        print("nothing to map: every contrast is noise-driven")
        return
    rows = []
    for name in todo:
        C = H.Contrast(name)
        o = pd.read_csv(H.OUT / "step3" / f"observed_d_{name}.csv").set_index("target")
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            t = C.ds.contrast("lab", targets=C.targets, levels=H.GROUPS, blocks=C.blocks, n_perm=0, n_boot=0).table
        o["c_control"] = t.c_control.reindex(o.index)
        o["donors_control"], o["donors_case"] = t.donors_control.reindex(o.index), t.donors_case.reindex(o.index)
        cut = o.c_control.quantile(1 / 3)
        se = np.sqrt((o.se_anchored_wp0 ** 2 + o.se_loo_wp0 ** 2) / 2)
        o["same_sign"] = np.sign(o.d_anchored) * np.sign(o.d_loo) > 0
        o["abs_diff"], o["pooled_se"] = (o.d_anchored - o.d_loo).abs(), se
        o["agree"] = o.same_sign & (o.abs_diff <= se)
        o["neuron"] = [bool(NEURON.match(k)) for k in o.index]
        o["glia_astro_oligo"] = [bool(GLIA.match(k)) for k in o.index]
        o["expect_agree"] = o.neuron | o.glia_astro_oligo
        o["expect_noise"] = (o.c_control <= cut) | (o.donors_control < 20) | (o.donors_case < 20)
        o["contrast"] = name
        rows.append(o.reset_index())
    m = pd.concat(rows, ignore_index=True)
    m["expectation"] = np.select([m.expect_agree & ~m.expect_noise, m.expect_noise & ~m.expect_agree,
                                  m.expect_agree & m.expect_noise], ["agree", "noise", "both"], "none")
    m.to_csv(H.OUT / "step4_map.csv", index=False)
    s = []
    for name, g in m.groupby("contrast"):
        for lab, sel in (("all targets", g.index == g.index), ("expected agree (neurons, astro/oligo)", g.expect_agree),
                         ("neurons", g.neuron), ("astrocyte/oligodendrocyte", g.glia_astro_oligo),
                         ("expected noise (low coupling or < 20 donors per group)", g.expect_noise),
                         ("expected agree and not expected noise", g.expect_agree & ~g.expect_noise),
                         ("expected noise and not expected agree", g.expect_noise & ~g.expect_agree)):
            x = g[sel]
            s.append(dict(contrast=name, set=lab, n=len(x), n_agree=int(x.agree.sum()), n_same_sign=int(x.same_sign.sum()),
                          share_agree=x.agree.mean() if len(x) else np.nan))
    pd.DataFrame(s).to_csv(H.OUT / "step4_summary.csv", index=False)
    print(pd.DataFrame(s).round(3).to_string(), flush=True)


if __name__ == "__main__":
    main()
