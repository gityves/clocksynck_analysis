"""A42 readout 1 and A103.3 at full size: the brain-vs-periphery r* of sheep and baboon with the registered bootstraps, in parallel.

Pre-registered in docs/PLAN_internal_desync.md (archive) (A42 readout 1; A103 H103.3). Same estimand, null and decision rule as
`compute_brain_periphery_species.py` (r* = |complex correlation| of the brain host and the peripheral host over animals,
disattenuated by the split-half reliabilities; `Dataset.agreement`, 50 splits; animal bootstrap, 200 replicates for the primary
analyses and 100 for the sensitivities; pass if r* < 1 and the upper bound < 1). That script runs them serially (hours); this one
splits the independent analyses into tasks that run in parallel and checkpoint (a finished task is skipped on restart):
  sheep_primary, sheep_identical_studies, sheep_pineal_pituitary, a103_3_sheep_no_per12, baboon_primary, baboon_sens.
THE FULL RUN TAKES HOURS (about one core-hour per 10 bootstrap replicates of a primary analysis; expect 3-6 h with 4 jobs); it has
not been run here (the serial script's subsample checks completed). `--subsample N` and `--n-boot` give a check in minutes.
The baboon timed readout and readouts 2-3 are cheap and stay in the serial script (`--skip-r1`).

    python compute_species_rstar_full.py [--jobs 4] [--n-boot 200] [--n-boot-sens 100] [--subsample N] [--only TASK ...]
Outputs (paper/results/figure_2/output/organs/[subsample/]): r1_full/<task>.csv (one row per task), r1_sheep_readout_full.csv,
a103_sheep_readout_full.csv, r1_baboon_readout_full.csv
Legacy: legacy/desync/validation/cross_species/scripts/r1.py, brain_artefact/scripts/a103_sheep.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

from clocksync.config import limit_threads  # noqa: E402
from common import cli, data  # noqa: E402
from common.paths import out_dir  # noqa: E402

import _shared as S  # noqa: E402
import compute_brain_periphery_species as B  # noqa: E402

FOLDER = "organs"
TASKS = ("sheep_primary", "sheep_identical_studies", "sheep_pineal_pituitary", "a103_3_sheep_no_per12", "baboon_primary",
         "baboon_sens")


def task(name, spb, bpb, nb, nbs, n_split, ck):
    """One independent analysis; returns its row(s) as a DataFrame and writes the checkpoint."""
    limit_threads(1)
    f = ck / f"{name}.csv"
    if f.exists():
        return name, pd.read_csv(f)
    rows = []
    if name.startswith("baboon"):
        ds = B.dataset(bpb, name="baboon")
        sens = name == "baboon_sens"
        brain, per, _ = B.groups("baboon", sorted(ds.compartments), sens)
        el = B.eligible(ds.table, brain, per)
        B.readout1("sens: pineal/pituitary/ocular as brain" if sens else "PRIMARY", ds, brain, per, nbs if sens else nb,
                   n_split, len(el), rows)
    else:
        g12 = [g for g in S.CLOCK if g not in ("PER1", "PER2")] if name.startswith("a103") else None
        ds = B.dataset(spb, genes=g12, mixed=True, name="sheep")
        comps = sorted(ds.compartments)
        if name == "sheep_pineal_pituitary":
            brain, per, _ = B.groups("sheep", comps, True)
            B.readout1("sens: pineal/pituitary as brain", ds, brain, per, nbs, n_split, 10, rows)
        else:
            brain, per, _ = B.groups("sheep", comps)
            el = B.eligible(ds.table, brain, per)
            if name == "sheep_identical_studies":
                t = ds.table[ds.table.resection.isin(el) & ds.table.compartment.isin(brain + per)]
                tt = pd.DataFrame(dict(a=t.resection, s=t.stratum, b=t.compartment.isin(brain)))
                sb, sp = tt[tt.b].groupby("a").s.apply(set), tt[~tt.b].groupby("a").s.apply(set)
                same = [a for a in el if sb[a] == sp[a]]
                if len(same) < 12:
                    rows.append(dict(analysis="sens: identical brain/periphery study sets", status=f"skipped ({len(same)} < 12 animals)"))
                else:
                    B.readout1(f"sens: identical brain/periphery study sets ({len(same)} animals)", ds.subset(donors=same),
                               brain, per, nbs, n_split, 10, rows)
            elif len(el) < 20 and nb > 3:
                rows.append(dict(analysis="PRIMARY", status=f"infeasible ({len(el)} < 20 animals)"))
            else:
                B.readout1("PRIMARY H103.3: no PER1/PER2" if g12 else "PRIMARY", ds, brain, per, nb, n_split, 10, rows)
    df = pd.DataFrame(rows)
    df.to_csv(f, index=False)
    return name, df


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--n-boot", type=int, default=200)
    ap.add_argument("--n-boot-sens", type=int, default=100)
    ap.add_argument("--only", nargs="*", choices=TASKS, default=None)
    args = ap.parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    ck = out / "r1_full"
    ck.mkdir(exist_ok=True)
    nb, nbs, n_split = (3, 2, 6) if run["quick"] and args.n_boot == 200 else (args.n_boot, args.n_boot_sens, 50)
    spb, _ = data._builder("sheep").build()
    if args.subsample:
        keep = np.random.default_rng(args.seed).permutation(sorted(spb.donor.unique()))[:args.subsample]
        spb = spb[spb.donor.isin(keep)]
    bpb, _ = data._builder("baboon").build()
    names = list(args.only or TASKS)
    res = dict(Parallel(n_jobs=min(args.jobs, 4))(delayed(task)(n, spb, bpb, nb, nbs, n_split, ck) for n in names))
    for fn, keys in (("r1_sheep_readout_full.csv", ("sheep_primary", "sheep_identical_studies", "sheep_pineal_pituitary")),
                     ("a103_sheep_readout_full.csv", ("a103_3_sheep_no_per12",)),
                     ("r1_baboon_readout_full.csv", ("baboon_primary", "baboon_sens"))):
        parts = [res[k] for k in keys if k in res]
        if parts:
            pd.concat(parts, ignore_index=True).to_csv(out / fn, index=False)
    pd.set_option("display.width", 220)
    print(pd.concat(res.values(), ignore_index=True).round(4).to_string(index=False))


if __name__ == "__main__":
    main()
