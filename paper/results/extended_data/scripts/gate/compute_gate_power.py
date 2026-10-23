"""Gate power and what a failed gate means: for every healthy atlas fit whose gate fails, the probability that the gate
would pass if a shared clock of the size seen in lung were present (A17 step 5; reading rules of A17).

Pre-registered in docs/PLAN_internal_desync.md (archive) (A17 step 5, gate power; A10 for the reading of failed gates). A failed
gate is uninformative, never "decoupled". Per healthy target (the 34 A17 healthy fits): the model-free gate of the pooled
table (`Dataset.gate`); for every failing target, `power.gate_power`: a shared clock of host covariance Sigma_f of the
anchored HLCA lung host and loading |z| = the median anchored amplitude of the HLCA immune classes (lymphocyte, myeloid)
is injected into every compartment of the real table by binomial thinning (per donor f ~ N(0, Sigma_f), clock-plane change
|z| f), and the gate is rerun (n_rep injections, 200 permutations); power = share of injections where both gate p values
are < 0.05. Reading: power < 0.5 "low power"; power >= 0.8 with p_pairing >= 0.05 "no shared signal"; power >= 0.8 with
p_acrophase >= 0.05 "shared, not clock-shaped"; else "indeterminate". Caveat of the pre-registration: the injected spread is
that of lung (a wide window), so power is overstated for hospital-hours cohorts.

    python compute_gate_power.py [--subsample N] [--seed S] [--jobs 4] [--n-rep N]
Outputs (paper/results/extended_data/output/gate/): gate_power_gates.csv (every healthy target), gate_power.csv (failing targets),
gate_power_lung_host.csv (Sigma_f and |z| used)
Legacy: legacy/desync/common/run_anchored.py (`gatepower` step) -> common_output/a17_gatepower.csv, a17_gates.csv
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402

import clocksync as cs  # noqa: E402
from clocksync import power  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import cli, data  # noqa: E402
from common.paths import out_dir  # noqa: E402

import _gate_shared as sh  # noqa: E402

FOLDER = "gate"
IMMUNE = ("lymphocyte", "myeloid")


def load(row, subsample, seed):
    pb, donors = data.load(row.dataset, subsample=subsample, seed=seed)
    if row.name == "lung_hlca" and "healthy" in donors:  # the A4 healthy table
        pb = pb[pb.donor.isin(donors.index[donors.healthy.astype(bool)])]
    return pb.reset_index(drop=True), donors


def lung_host(args, run):
    """Sigma_f of the anchored HLCA lung host and the median anchored amplitude |z| of its immune classes."""
    row = sh.targets().loc["lung_hlca"]
    pb, donors = load(row, None, args.seed)  # the lung host is always fitted on all donors
    anc = [a for a in sh.anchor_set("lung_hlca", "learned") if a in set(pb.compartment)]
    ds = cs.Dataset(pb, anchors=anc, name="lung_hlca")
    e, ft = sh.T.host_estimate(ds, anc)
    fam = pd.Series(sh.family_map())
    # `Fit.table` of all compartments (anchors on leave-one-out hosts, the others on the full host)
    full = ds.fit(n_boot=0, gate=False).table
    imm = [k for k in full.index if fam.get(k) in IMMUNE]
    z_abs = float(np.sqrt(full.loc[imm, "z2"]).median())
    return np.asarray(e["sigma_f"], float), z_abs, imm


def one(name, sigma_f, z_abs, n_rep, run, args):
    limit_threads(1)
    row = sh.targets().loc[name]
    pb, donors = load(row, args.subsample, args.seed)
    ds = cs.Dataset(pb, anchors=None, name=name)
    g = ds.gate(n_perm=run["n_perm"], seed=args.seed)
    head = dict(target=name, donors=len(ds.donor_ids), classes=len(ds.compartments), gate_p_pairing=g.p_pairing,
                gate_p_acrophase=g.p_acrophase, gate_clock_share=g.clock_share, gate_passed=g.passes,
                clock_umi_median=float(ds.table[list(ds.waveform.genes)].sum(axis=1).median()))
    if g.passes:
        return head
    gp = power.gate_power(ds.study, sigma_f, z_abs, n_rep=n_rep, n_perm=200, seed=args.seed)
    reading = ("low power" if gp < 0.5 else "no shared signal" if gp >= 0.8 and g.p_pairing >= 0.05 else
               "shared, not clock-shaped" if gp >= 0.8 and g.p_acrophase >= 0.05 else "indeterminate")
    return head | dict(gate_power=gp, reading=reading, z_abs=z_abs)


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--n-rep", type=int, default=None, help="injections per failing target (default 100; 5 with --subsample)")
    args = ap.parse_args(argv)
    run = cli.setup(args)
    out = out_dir(FOLDER, run["quick"])
    n_rep = args.n_rep if args.n_rep is not None else (5 if run["quick"] else 100)
    sigma_f, z_abs, imm = lung_host(args, run)
    pd.DataFrame([dict(sigma_xx=sigma_f[0, 0], sigma_xy=sigma_f[0, 1], sigma_yy=sigma_f[1, 1], z_abs=z_abs,
                       immune_classes=";".join(imm))]).to_csv(out / "gate_power_lung_host.csv", index=False)
    names = [n for n in sh.targets().index if sh.targets().loc[n, "kind"] == "healthy"]
    if run["quick"]:
        names = ["blood_onek1k", "organ_adipose", "ia_spleen", "lung_hlca", "organ_breast"]
    rows = Parallel(n_jobs=args.jobs)(delayed(one)(n, sigma_f, z_abs, n_rep, run, args) for n in names)
    d = pd.DataFrame(rows)
    d.to_csv(out / "gate_power_gates.csv", index=False)
    fail = d[~d.gate_passed.astype(bool)]
    fail.to_csv(out / "gate_power.csv", index=False)
    pd.set_option("display.width", 220)
    print(f"lung host: |z| {z_abs:.3f}, Sigma_f {np.round(sigma_f, 3).tolist()}")
    print(f"{int((~d.gate_passed.astype(bool)).sum())} of {len(d)} healthy gates fail")
    print(fail[["target", "donors", "classes", "gate_p_pairing", "gate_p_acrophase", "gate_power", "reading"]].round(3).to_string(index=False))


if __name__ == "__main__":
    main()
