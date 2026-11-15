"""A265 part (iii), the simulation: `score_a265.py --part sim` under its own name (own run.sh log, so that it can run beside
the real-data part).

    python simulate_a265.py [--reps 10] [--jobs 12]
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import score_a265  # noqa: E402

if __name__ == "__main__":
    sys.argv = [sys.argv[0], "--part", "sim", *sys.argv[1:]]
    score_a265.main()
