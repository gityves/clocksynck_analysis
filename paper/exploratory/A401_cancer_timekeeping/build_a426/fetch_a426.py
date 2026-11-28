"""A426 data acquisition: GSE312983 (untreated primary CRC snRNA, CellBender-filtered h5) and the Navin TNBC h5ad.

Only metadata and the files named below are fetched, each into its own new directory under $SCCC_DATA, `raw/`; every file
gets a sha256 line in `raw/SHA256SUMS`. Nothing downloaded is executed or unpacked by code (the tar members of GSE312983 are
read by HTTP range at offsets computed from the GEO file list, header and HDF5 signature verified). Generic client agent.

    fetch_a426.py gse312983 [--workers 3]
    fetch_a426.py navin
"""
from __future__ import annotations

import argparse
import hashlib
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pandas as pd
import requests

DATA = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian")
HERE = Path(__file__).resolve().parent
BUILD = HERE.parent / "output" / "a426" / "build"
UA = {"User-Agent": "python-requests/2.31"}
TAR = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE312nnn/GSE312983/suppl/GSE312983_RAW.tar"
NAVIN = "https://datasets.cellxgene.cziscience.com/e94bd3cc-6271-424a-baac-12f8eb320a0e.h5ad"


def log(*a):
    print(time.strftime("%F %T"), *a, flush=True)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""):
            h.update(b)
    return h.hexdigest()


def get(headers, tries=5):
    for k in range(tries):
        try:
            r = requests.get(TAR, headers={**UA, **headers}, timeout=300)
            r.raise_for_status()
            return r.content
        except Exception as e:  # noqa: BLE001
            log("retry", k, repr(e)[:100])
            time.sleep(10 * (k + 1))
    raise RuntimeError("fetch failed")


def member(row, raw: Path):
    out = raw / row.Name
    if out.exists() and out.stat().st_size == row.Size:
        return row.Name, sha256(out), row.Size
    hdr = get({"Range": f"bytes={row.offset}-{row.offset + 511}"})
    assert hdr[:100].split(b"\0")[0].decode() == row.Name, row.Name
    d = get({"Range": f"bytes={row.offset + 512}-{row.offset + 511 + row.Size}"})
    assert len(d) == row.Size and d[:8] == b"\x89HDF\r\n\x1a\n", (row.Name, len(d))
    tmp = out.with_suffix(".part")
    tmp.write_bytes(d)
    tmp.rename(out)
    return row.Name, hashlib.sha256(d).hexdigest(), row.Size


def gse312983(workers):
    raw = DATA / "gse312983" / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    ix = pd.read_csv(BUILD / "GSE312983_tar_index.csv")
    m = pd.read_csv(BUILD / "GSE312983_samples.csv", index_col=0)
    untreated = set(m.index[m["treatment status"] == "No"])
    sel = ix[ix.gsm.isin(untreated) & (ix.kind == "cellbender_filtered")]
    log(f"{len(sel)} untreated CellBender-filtered files, {sel.Size.sum() / 1e9:.1f} GB")
    res = []
    with ThreadPoolExecutor(workers) as ex:
        for k, r in enumerate(ex.map(lambda t: member(t, raw), [r for r in sel.itertuples()])):
            res.append(r)
            if k % 10 == 0:
                log(f"{k + 1}/{len(sel)} {r[0]}")
    (raw / "SHA256SUMS").write_text("".join(f"{h}  {n}\n" for n, h, _ in sorted(res)))
    m.loc[sorted(untreated)].to_csv(raw.parent / "samples_untreated.csv")
    log("done")


def navin():
    raw = DATA / "navin_tnbc_cxg" / "raw"
    raw.mkdir(parents=True, exist_ok=True)
    out = raw / "e94bd3cc-6271-424a-baac-12f8eb320a0e.h5ad"
    if not out.exists():
        tmp = out.with_suffix(".part")
        with requests.get(NAVIN, headers=UA, stream=True, timeout=300) as r:
            r.raise_for_status()
            with open(tmp, "wb") as f:
                for b in r.iter_content(1 << 22):
                    f.write(b)
        tmp.rename(out)
    (raw / "SHA256SUMS").write_text(f"{sha256(out)}  {out.name}\n")
    log("done", out.stat().st_size)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("what", choices=("gse312983", "navin"))
    ap.add_argument("--workers", type=int, default=3)
    a = ap.parse_args()
    gse312983(a.workers) if a.what == "gse312983" else navin()
