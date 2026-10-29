"""A161 shared helpers: GEO series-matrix / SOFT metadata parsing (python -I; downloads parsed, never executed)."""
import gzip, re
from pathlib import Path
import numpy as np, pandas as pd

DATA = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian")
EN = DATA / "endometrium_gse51981"
SC = DATA / "scout_endometrium"
OUT = Path(__file__).resolve().parents[1] / "output"
HERE = Path(__file__).resolve().parent


def header(path) -> pd.DataFrame:
    """Sample table (one row per sample) from the '!Sample_' header lines of a series matrix; characteristics split key: value."""
    rows = {}
    with gzip.open(path, "rt", errors="replace") as f:
        for ln in f:
            if ln.startswith("!series_matrix_table_begin"):
                break
            if not ln.startswith("!Sample_"):
                continue
            k, *v = ln.rstrip("\n").split("\t")
            k = k[len("!Sample_"):]
            v = [x.strip('"') for x in v]
            if k == "characteristics_ch1":
                for j, x in enumerate(v):
                    kk, _, vv = x.partition(": ") if ": " in x else ("char", "", x)
                    rows.setdefault((j, kk.strip()), vv)
                continue
            for j, x in enumerate(v):
                rows.setdefault((j, k), x)
    df = pd.Series(rows).unstack()
    return df


def soft_samples(path) -> pd.DataFrame:
    """Sample table from a family SOFT (only ^SAMPLE / !Sample_ lines; expression tables not read)."""
    recs, cur = [], None
    with gzip.open(path, "rt", errors="replace") as f:
        for ln in f:
            if ln.startswith("^SAMPLE"):
                cur = {"gsm": ln.split("=")[1].strip()}; recs.append(cur)
            elif cur is not None and ln.startswith("!Sample_"):
                k, _, v = ln.rstrip("\n").partition(" = ")
                k = k[len("!Sample_"):]
                if k == "characteristics_ch1":
                    kk, _, vv = v.partition(": ")
                    k, v = (kk, vv) if vv else ("char", v)
                if k not in cur:
                    cur[k] = v
                else:
                    cur[k] += " | " + v
            elif ln.startswith("^PLATFORM") or ln.startswith("^SERIES"):
                cur = None
    return pd.DataFrame(recs).set_index("gsm")


def matrix(path) -> pd.DataFrame:
    """Probes x samples values of a series matrix."""
    with gzip.open(path, "rt", errors="replace") as f:
        for ln in f:
            if ln.startswith("!series_matrix_table_begin"):
                break
        return pd.read_csv(f, sep="\t", index_col=0, skipfooter=1, engine="python").astype(float)
