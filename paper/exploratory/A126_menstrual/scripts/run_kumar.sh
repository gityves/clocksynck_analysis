#!/usr/bin/env bash
# A126 Kumar pseudobulk extraction (census_env); launched through paper/common/run.sh
exec /home/maxine/miniconda3/envs/census_env/bin/python -u "$(dirname "$0")/extract_kumar_pb.py" --threads 4
