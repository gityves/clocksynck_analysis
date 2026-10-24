"""Common command line of the compute scripts: `--subsample N` (check on N random donors, fewer permutations,
outputs under .../subsample/), `--seed`, `--jobs`."""

from __future__ import annotations

import argparse

from clocksync.config import limit_threads


def parser(doc: str) -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(description=(doc or "").split("\n")[0])
    ap.add_argument("--subsample", type=int, default=None, help="check on N random donors per dataset")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--jobs", type=int, default=1, help="parallel jobs (at most 8 workers in total on the shared machine)")
    return ap


def setup(args) -> dict:
    """Thread cap; resampling sizes for a full run or a quick subsample check."""
    limit_threads(1 if args.jobs > 1 else 4)
    quick = args.subsample is not None
    return dict(n_perm=50 if quick else 1000, n_boot=20 if quick else 200, n_null=10 if quick else 200,
                quick=quick)
