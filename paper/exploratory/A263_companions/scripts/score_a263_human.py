"""A263 for one species under its own run.sh log (`score_a263.py --species human`)."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import score_a263  # noqa: E402

if __name__ == "__main__":
    sys.argv[1:1] = ["--species", "human"]
    score_a263.main()
