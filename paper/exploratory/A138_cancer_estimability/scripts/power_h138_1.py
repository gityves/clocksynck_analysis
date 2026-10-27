"""A138 note 4: H138.1 power simulation, 30 replicates per c, 100 bootstraps per fit (wrapper of power.py, own log name).

    python power_h138_1.py [--subsample N] [--jobs J]
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import power  # noqa: E402

if __name__ == "__main__":
    power.main(["--only", "h138_1", "--n-rep", "30", "--n-boot", "100", *sys.argv[1:]])
