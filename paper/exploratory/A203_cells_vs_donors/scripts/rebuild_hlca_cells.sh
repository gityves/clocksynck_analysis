#!/usr/bin/env bash
# A203 fix: rebuild the HLCA per-cell clock-gene counts with ARNTL (Census symbol BMAL1) into a fresh folder, then check
# that cell ids and the 13 other genes equal the stored file. The stored file is replaced by hand after the check.
set -euo pipefail
repo=$(cd "$(dirname "$0")/../../../.." && pwd)
cd "$repo"
data=${SCCC_DATA:-/mnt/max8tb/data/CoPhaser_data/cancer_circadian}
out="$data/cells/hlca_full_rebuild_2026-11-04"
PYTHONPATH="$repo/clocksync/src" "$HOME/miniconda3/envs/census_env/bin/python" -u "$repo/paper/common/datasets/cells_census/build.py" hlca_full --threads 8 --out "$out"
"$(conda info --base)/envs/sccc/bin/python" -u - "$out" <<'PY'
import sys; from pathlib import Path
sys.path.insert(0, "paper/common/datasets/cells_census"); sys.path.insert(0, "paper")
import build as B, pandas as pd
out = Path(sys.argv[1])
print("check_columns (ids and 13 genes equal):", B.check_columns("hlca_full", out))
print("ARNTL total:", int(pd.read_parquet(out / "cells.parquet", columns=["ARNTL"]).ARNTL.sum()))
PY
