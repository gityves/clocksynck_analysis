"""A421 (registered in paper/docs/registered_experiments.md, "A421"): `a419.py --variant intrinsic` under its own log name.

    python a421.py [--jobs 4] [--genome-jobs 2]                          (long: paper/common/run.sh)
"""
import sys

from a419 import main

if __name__ == "__main__":
    sys.argv += ["--variant", "intrinsic"]
    main()
