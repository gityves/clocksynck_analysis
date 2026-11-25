"""Slide data (results, drives, 4.2c): the reference map of slide 3.1 for fast deaths, slow deaths, ventilated donors and
brain death.

GTEx v10 donors with a recorded time of death in four groups (as A293 part 1, brain death with all its donors): fast deaths (Hardy 1-2), slow deaths
(Hardy 4), ventilated post-mortem donors with a living brain (Hardy 0, post-mortem cohort) and brain-dead organ donors
(Hardy 0, organ-donor cohort). Each group is fitted on
its own (`Dataset.fit(gate=False)`, 14 clock genes, peripheral tissues as anchors read leave-one-out, brain regions as
targets; A263 R2 readout as `organs_compute.py`), tissues with >= 30 donors in the group, donor bootstrap for per-tissue
intervals. Offsets and gains are relative to each group's own host: the faces compare arrangements, not loudness
between groups (that is slide 4.2b). Descriptive, post hoc display. Aggregates only (one row per group and tissue).

    RUN_CAP_GB=30 paper/common/run.sh paper/presentation/scripts/organs_groups_compute.py [--n-boot 200] [--groups ...]
Output: paper/presentation/output/organs/couplings_<group>.csv (one file per group)
"""
import argparse
import sys
import time
from pathlib import Path

PAPER = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PAPER))
import clocksync as cs  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import gtex  # noqa: E402

OUT = PAPER / "presentation" / "output" / "organs"
MIN_DONORS = 30
SEED = 263
GROUPS = {"fast": (1, 2), "slow": (4,), "ventilated": (0,), "brain_death": (0,)}
COHORT = {"ventilated": "Postmortem", "brain_death": "Organ Donor (OPO)"}  # A132 PM and OPO cohorts


def load(group):
    pb, _ = gtex.load(hardy=GROUPS[group], genes="clock", depth="library_size", require_tod=True)
    if group in COHORT:
        coh = gtex.phenotypes(["COHORT"]).COHORT
        pb = pb[pb.donor.map(coh).eq(COHORT[group]).to_numpy()]
    n = pb.groupby("compartment").donor.nunique()
    return pb[pb.compartment.isin(n.index[n >= MIN_DONORS])].reset_index(drop=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-boot", type=int, default=200)
    ap.add_argument("--groups", nargs="+", choices=list(GROUPS), default=list(GROUPS))
    args = ap.parse_args()
    limit_threads(4)
    OUT.mkdir(parents=True, exist_ok=True)
    wf = cs.gtex()
    for g in args.groups:
        t0 = time.time()
        pb = load(g)
        comps = sorted(pb.compartment.unique())
        brain = [c for c in comps if c.startswith("Brain -")]
        per = [c for c in comps if c not in brain and c != "Pituitary"]
        ds = cs.Dataset(pb[["donor", "compartment", "depth", *wf.genes]], anchors=per, waveform=wf, name=f"groups_{g}",
                        min_donors=5)
        tab = ds.fit(n_boot=args.n_boot, gate=False, targets=comps, seed=SEED).table
        tab = tab.assign(death=g, donors_group=pb.donor.nunique(),
                         group=["brain" if t in brain else "excluded" if t not in per else "periphery" for t in tab.index])
        tab.rename_axis("tissue").to_csv(OUT / f"couplings_{g}.csv")
        print(f"{g}: {pb.donor.nunique()} donors, {len(comps)} tissues, {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
