"""Adult power of the A266 L1 readouts (core-loop waveform), before any fetal statistic: the A264 feasibility
(`A264/scripts/feasibility_a264.py`: distinct GTEx adrenal glands thinned to the cortex counts of each cohort, adult clock at
1 / 0.75 / 0.5 and shuffled, 200 subsamples, 500 permutations) on the tests of `design_a266`.
Each class is thinned to its own counts (`clock_detection_classes.csv`).

    python feasibility_a266.py [--n-sub 200] [--jobs 8]
Outputs (`../output/`): feasibility_design.csv, feasibility_subsamples.csv, feasibility_passrates.csv, adult_scores.csv.
"""
import sys
from pathlib import Path

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parents[2] / "A264_fetal_adrenal_covariates" / "scripts"))
import design_a266 as d  # noqa: E402
import feasibility_a264 as F  # noqa: E402

if __name__ == "__main__":
    F.d.tests = d.tests
    F.OUT, F.SEED = d.OUT, d.SEED
    F.main()
