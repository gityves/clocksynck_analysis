"""Slide data (results, R2.1): coupling of every GTEx tissue to the body's clock state (peripheral-anchor host).

The A263 R2 fit (registered readout; A263 saved only the brain and periphery medians): GTEx v10 fast deaths, tissues with
>= 40 donors, 14 clock genes, `Dataset.fit(gate=False)` with the peripheral tissues as anchors (peripheral tissues read
leave-one-out; brain = the "Brain -" regions, pituitary excluded from both groups but read as a target), here with a donor
bootstrap for per-tissue intervals. Aggregates only (one row per tissue; GTEx controlled access).

    RUN_CAP_GB=30 paper/common/run.sh paper/presentation/scripts/organs_compute.py [--n-boot 200]
Output: paper/presentation/output/organs/couplings.csv
"""
import argparse
import sys
from pathlib import Path

PAPER = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(PAPER / "exploratory" / "A263_companions" / "scripts"))
import clocksync as cs  # noqa: E402
import score_a263 as A  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402

OUT = PAPER / "presentation" / "output" / "organs"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-boot", type=int, default=200)
    args = ap.parse_args()
    limit_threads(4)
    OUT.mkdir(parents=True, exist_ok=True)
    pb = A.table("human")
    brain, per, excl = A.groups("human", sorted(pb.compartment.unique()))
    cols = ["donor", "compartment", "depth", *A.WF.genes]
    ds = cs.Dataset(pb[cols], anchors=per, waveform=A.WF, name="human", min_donors=5)
    tab = ds.fit(n_boot=args.n_boot, gate=False, targets=brain + per + excl, seed=A.SEED).table
    tab = tab.assign(group=["brain" if t in brain else "excluded" if t in excl else "periphery" for t in tab.index])
    tab.rename_axis("tissue").to_csv(OUT / "couplings.csv")
    m = tab.groupby("group").coupling.median()
    print(m.round(3).to_dict(), f"{len(tab)} tissues", flush=True)


if __name__ == "__main__":
    main()
