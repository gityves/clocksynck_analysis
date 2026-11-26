#!/usr/bin/env bash
# A401 Census pseudobulk build (census_env, 4 workers); launch through paper/common/run.sh
exec /home/maxine/miniconda3/envs/census_env/bin/python -u "$(dirname "$0")/build_census_tumour.py" --workers 4 --threads 2 "$@"
