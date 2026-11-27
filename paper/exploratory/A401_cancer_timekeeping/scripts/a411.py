"""A411 (registered in paper/docs/registered_experiments.md, "A411"): `a410.py` with x = malignant G2M index (own log name).

    python a411.py --stage power|score [--jobs 6] [--genome-jobs 2]      (long: paper/common/run.sh)
"""
import sys

from a410 import main

if __name__ == "__main__":
    sys.argv += ["--x", "prolif"]
    main()
