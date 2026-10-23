"""A41 readout 4c: the A34 primary (sheep tissue offsets against the human timed tissue phase) with the baboon waveform as the reference.

Pre-registered as A41 (2026-10-16; docs/PLAN_internal_desync.md (archive)): "Readout 4 (secondary): sheep (A34) refitted with the full-day baboon reference:
gate; tissue couplings (Spearman with the human-reference couplings); sheep-vs-human tissue-offset r (A34's primary). If these barely change, A34 does not
depend on which primate reference is used." Readouts 4a, 4b are in `compute_baboon_reference.py`. Here (4c): sheep (`sheep` builder, 343 animals, 45
tissues as anchors, BioProject strata) fitted with the human GTEx waveform and with the baboon waveform fitted on the 12 full-day animals (same inputs as the
A28 waveforms: `fit_waveform` of log((count + 0.5) / depth), ZT, group = tissue); the 21 matched tissues of A34 (`species/compute_sheep.MATCH`); the
human truth d_t = -arg g_t of the GTEx fast-death gains (table `a34_human_gain.csv` of `species`, run first; aggregates only). Estimand:
circular correlation (Jammalamadaka-SenGupta) of the sheep offsets with d_t, one-sided permutation p (10,000), and the Fisher-Lee r (descriptive), for each reference.

    python compute_sheep_reference_primary.py [--subsample N] [--n-perm P]
Outputs: paper/results/extended_data/output/gate/a41_sheep_primary.csv, a41_sheep_offsets.csv (per tissue, both references)
Legacy: legacy/desync/validation/baboon/scripts/a41_sheep.py (4c) -> output/a41_sheep_readouts.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
sys.path.insert(0, str(HERE))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import cli  # noqa: E402
from common.paths import out_dir  # noqa: E402

import _gate_shared as sh  # noqa: E402
from compute_baboon_reference import G, baboon_waveform  # noqa: E402
sys.path.insert(0, str(Path(__file__).resolve().parents[4]))  # paper/
from common.paths import group_dir, group_out  # noqa: E402
sys.path.insert(0, str(group_dir("species")))  # after the S02 helpers: both folders have an `_shared`
sys.modules.pop("_shared", None)
import compute_sheep as CS  # noqa: E402

FOLDER = "gate"


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--n-perm", type=int, default=None)
    args = ap.parse_args(argv)
    run = cli.setup(args)
    limit_threads(4)
    out = out_dir(FOLDER, run["quick"])
    n_perm = args.n_perm or (1000 if run["quick"] else 10_000)
    pb, donors = sh.s02.baboon()
    wf_b = baboon_waveform(pb, donors, sorted(pb.donor.unique()), "baboon_full")
    spb, _ = sh.s02.built("sheep")
    if args.subsample:
        keep = np.random.default_rng(args.seed).permutation(sorted(spb.donor.unique()))[: max(args.subsample, 30)]
        spb = spb[spb.donor.isin(keep)]
    cols = [c for c in ("donor", "compartment", "depth", "stratum") if c in spb] + G
    hum = pd.read_csv(group_out("species") / ("subsample" if run["quick"] else "") / "a34_human_gain.csv").set_index("tissue")
    rows, offs = [], {}
    for lab, wf in (("human_ref", cs.gtex()), ("baboon_ref", wf_b)):
        ds = cs.Dataset(spb[cols], anchors=sorted(spb.compartment.unique()), waveform=wf, mixed_strata=True, min_donors=5, name="sheep")
        off = ds.fit(n_boot=0, gate=False).table
        offs[lab] = off
        r, m = CS.primary(off, hum, CS.MATCH, n_perm)
        a = np.array([off.offset_h[k] for k, _ in m]) / CS.HOURS
        b = np.array([hum.d_h[v] for _, v in m]) / CS.HOURS
        rows.append(dict(reference=lab, n=r["n"], r_js=r["r"], p_one_sided=r["p"], r_fisher_lee=cs.stats.fisher_lee(a, b),
                         jackknife_min=r.get("jackknife_min"), jackknife_max=r.get("jackknife_max"), tissues=";".join(k for k, _ in m)))
    pd.DataFrame({f"{k}_{c}": v[c] for k, v in offs.items() for c in ("offset_h", "coupling")}).to_csv(out / "a41_sheep_offsets.csv")
    R = pd.DataFrame(rows)
    R.to_csv(out / "a41_sheep_primary.csv", index=False)
    print(R.drop(columns="tissues").round(4).to_string(index=False))


if __name__ == "__main__":
    main()
