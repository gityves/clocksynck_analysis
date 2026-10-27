"""WP0 / PI Q3 (2026-10-27): recompute two figure panels on the anchored A17+ estimator.

Re-estimation of existing results with the current estimator (anchors learned by A17, `cs.Dataset(..., anchors=...)`,
`Dataset.fit`; targets read on their host, an anchor read on the leave-one-out host of the other anchors). Not a new
hypothesis. Three panels were named by the PI; panel 1 (GTEx departures) is not run: it is a timed
analysis against recorded time of death with no host, so it has no anchored version (PI, 2026-10-27: kept, readout named).

  panel2  Figure 2, "Neurons carry the brain host, glia less; retina". Every healthy brain atlas (PsychAD HBCC and MSSM,
          SEA-AD DLPFC and MTG, PD control) and the retina, under both frozen rules of `common/anchors_a17.csv`
          (learned, fixed neuronal/stromal-vascular), `Dataset.fit` with `n_boot` donor bootstraps; neuron-minus-glia
          mean coupling with a percentile interval from the same replicates (`Fit.compare`). The family of a
          compartment is the `family` column of the frozen table. The gate is read from
          figure_2/output/atlas_fits/a17_healthy_gates.csv (all six pass); it is not rerun.
  panel3  Figure 5, "Fibrotic-lung T cells lose relative amplitude, two cohorts". The A98 decomposition
          (`figure_5/scripts/lung_breaks/compute_decomposition.py`, imported unchanged) rerun for Natri ILD and HLCA fibrotic
          on the A22-confirmed anchors: learned anchors minus the ones flagged `drop` in
          figure_4/output/spread_multiplicity/a22_flags.csv (Natri: lung:EC_capillary, which the M3 A98 and A17 tables
          still carry). Quantities: relative amplitude |z| (`amp`), coupling, own noise psi, with donor-bootstrap
          intervals and the fibroblast specificity companion on the same draws. The anchored d* contrast of the A22
          confirmed set is read from a23_inputs.csv, not recomputed.

    python q3_anchored.py [--subsample N] [--seed S] [--jobs 6] [--only panel2 panel3] [--n-boot B]
Outputs (paper/exploratory/WP0_integrity/output/q3_anchored/, subsample/ under it for the check):
  panel2_fits.csv, panel2_neuron_glia.csv, panel3_groups.csv, panel3_differences.csv, panel3_anchors.csv
"""

from __future__ import annotations

import importlib.util
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

ROOT = Path(__file__).resolve().parents[4]
PAPER = ROOT / "paper"
sys.path.insert(0, str(PAPER))
import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import cli  # noqa: E402

OUT = Path(__file__).resolve().parents[1] / "output" / "q3_anchored"
RES = PAPER / "results"
TARGETS2 = ["psychad_hbcc", "psychad_mssm", "seaad_dlpfc", "seaad_mtg", "pd_control", "retina"]
RULES = (("learned", "anchor"), ("fixed", "anchor_fixed"))
CONTRASTS3 = ["ild_natri", "hlca_fibrotic"]
# old numbers (archive, free fit / pre-A17) quoted next to the new ones in the README
OLD = {"neuron_coupling_brain": "about 0.9 (archive desync/README.md, four cohorts, free fit)",
       "retina_neuron_minus_glia": "+0.17 [0.04, 0.35] (A2 primary, 95 donors, free fit)"}


def _load(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def _atlas():
    """`Target` of the atlas_fits group (loaded under its own name: `_shared` exists in several groups)."""
    d = RES / "figure_2/scripts/atlas_fits"
    sys.path.insert(0, str(d))
    try:
        return _load("q3_atlas_shared", d / "_shared.py")
    finally:
        sys.path.remove(str(d))


def _decomposition():
    d = RES / "figure_5/scripts/lung_breaks"
    sys.path.insert(0, str(d))
    sys.modules.pop("_shared", None)
    return _load("q3_decomposition", d / "compute_decomposition.py")


# panel 2 ---------------------------------------------------------------------------------------------------------

def fit_one(name: str, rule: str, col: str, n_boot: int, subsample, seed: int):
    limit_threads(1)
    warnings.filterwarnings("ignore")
    A = _atlas()
    T = A.Target(name, subsample, seed)
    fro = A.frozen()
    fro = fro[fro.target == name]
    anc = [a for a in fro.compartment[fro[col].astype(bool)] if a in set(T.pb.compartment)]
    if len(anc) < 3:
        return None
    fit = T.dataset(anchors=anc).fit(n_boot=n_boot, seed=seed, gate=False)
    fam = fro.drop_duplicates("compartment").set_index("compartment").family
    tab = fit.table.reset_index().rename(columns={"index": "compartment"})
    tab = tab.assign(target=name, rule=rule, anchor=tab.compartment.isin(anc), family=tab.compartment.map(fam),
                     n_anchors=len(anc), n_donors=T.pb.donor.nunique())
    sets = {f: [c for c in tab.compartment[tab.family == f] if c in fit.replicates.columns] for f in ("neuron", "glia")}
    cmp = fit.compare({k: v for k, v in sets.items() if v})
    cmp = cmp.assign(target=name, rule=rule, n_neuron=len(sets["neuron"]), n_glia=len(sets["glia"]),
                     anchors_neuron=int(tab[tab.family == "neuron"].anchor.sum()),
                     anchors_glia=int(tab[tab.family == "glia"].anchor.sum()))
    print(f"panel2 {name} {rule}: {len(anc)} anchors", flush=True)
    return tab, cmp


def panel2(args, run):
    jobs = [(n, r, c) for n in TARGETS2 for r, c in RULES]
    res = Parallel(n_jobs=args.jobs)(delayed(fit_one)(n, r, c, args.n_boot or (50 if run["quick"] else 500),
                                                      args.subsample, args.seed) for n, r, c in jobs)
    res = [x for x in res if x is not None]
    fits = pd.concat([a for a, _ in res], ignore_index=True)
    cmp = pd.concat([b for _, b in res], ignore_index=True)
    gates = pd.read_csv(RES / "figure_2/output/atlas_fits/a17_healthy_gates.csv").set_index("target")
    cmp["gate_passed"] = cmp.target.map(gates.passed)
    fits.to_csv(OUT / "panel2_fits.csv", index=False)
    cmp.to_csv(OUT / "panel2_neuron_glia.csv", index=False)
    print(cmp[cmp.set == "neuron - glia"][["target", "rule", "estimate", "lo", "hi", "anchors_neuron", "anchors_glia"]]
          .round(3).to_string(index=False))


# panel 3 ---------------------------------------------------------------------------------------------------------

def confirmed_anchors() -> pd.DataFrame:
    """Learned anchors of the two lung contrasts with the A22 drops removed (a22_flags.csv `drop`)."""
    fl = pd.read_csv(RES / "figure_4/output/spread_multiplicity/a22_flags.csv")
    fl = fl[fl.contrast.isin(CONTRASTS3)]
    return fl.assign(kept=~fl["drop"].astype(bool))[["contrast", "anchor", "d_star", "p_star", "drop", "kept"]]


def lung_one(name: str, n_boot: int, subsample, seed: int):
    D = _decomposition()
    fl = confirmed_anchors()
    keep = fl[(fl.contrast == name) & fl.kept].anchor.tolist()
    D.sh.anchor_set = lambda n, rule="learned": keep  # the decomposition reads its anchors through this call
    return D.run_contrast(name, n_boot, subsample, seed)


def panel3(args, run):
    nb = args.n_boot or (20 if run["quick"] else 200)
    res = Parallel(n_jobs=min(args.jobs, len(CONTRASTS3)))(
        delayed(lung_one)(n, nb, args.subsample, args.seed) for n in CONTRASTS3)
    est = pd.concat([e for e, _ in res], ignore_index=True)
    dif = pd.concat([d for _, d in res], ignore_index=True)
    est.to_csv(OUT / "panel3_groups.csv", index=False)
    dif.to_csv(OUT / "panel3_differences.csv", index=False)
    confirmed_anchors().to_csv(OUT / "panel3_anchors.csv", index=False)
    print(dif[(dif.quantity.isin(["amp", "coupling", "psi"]))][["contrast", "role", "compartment", "quantity", "control",
                                                              "case", "diff", "lo", "hi", "reading"]]
          .round(3).to_string(index=False))


def main(argv=None):
    ap = cli.parser(__doc__)
    ap.add_argument("--only", nargs="+", default=["panel2", "panel3"], choices=["panel2", "panel3"])
    ap.add_argument("--n-boot", type=int, default=None,
                    help="bootstrap draws (panel2 500, panel3 200; 50 and 20 with --subsample)")
    args = ap.parse_args(argv)
    run = cli.setup(args)
    global OUT
    OUT = OUT / "subsample" if run["quick"] else OUT
    OUT.mkdir(parents=True, exist_ok=True)
    for p in args.only:
        {"panel2": panel2, "panel3": panel3}[p](args, run)


if __name__ == "__main__":
    main()
