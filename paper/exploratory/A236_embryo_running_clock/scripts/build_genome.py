"""A236: genome-wide cell type x embryo pseudobulks of the late window (E15.0-E18.75) of the Qiu et al. 2024 atlas.

Embryos embryo_57 ... embryo_72 (16) lie in runs 20, 21, 22, 23, 26, 27 of GSE228590. Each run matrix (MatrixMarket,
genes x cells, sorted by cell, ~2-5 GB compressed) is streamed `curl | gzip -dc | python -I build_genome.py stream`
and never written to disk. The parser reads stdin only, runs isolated (`-I`) from a directory without downloaded
files, and sums the UMI counts of all 49,585 genes per (embryo, cell type) key, plus the number of cells and the
all-gene depth (same definitions as the stored panel table). Stream checks as in the archive builder: line count ==
nnz, 49,585 genes, cell count == annotation, no non-integer token, nothing out of range; a failure aborts the build.
Per-run partials (sparse int) are merged into one table. No clock quantity is computed.

Inputs  $SCCC_DATA/natexp/qiu2024/prep/{keys.csv, codes.run_<r>.npy, panel.csv, embryo_runs.csv} (archive prep),
        gene annotation of GSE228590.
Output  $SCCC_DATA/natexp/qiu2024/genome_late/pseudobulk_celltype.h5ad (obs: donor, compartment, trajectory, n_cells,
        depth; var: gene id, symbol; X = counts, CSR float32), PROVENANCE.md, partial/ (per-run partials and checks).

    RUN_CAP_GB=30 paper/common/run.sh paper/exploratory/A236_embryo_running_clock/scripts/build_genome.py [--jobs 4]
    build_genome.py --test-lines 3000000 --runs 21     # parser test on the head of one stream (relaxed check)
    build_genome.py --check                            # panel comparison only (after a build)
Runtime: ~1-2 h (network bound, 4 streams); RAM ~1 GB per stream.
"""
import argparse
import datetime
import hashlib
import json
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian/natexp/qiu2024")
ANNOT = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian/scout_natexp/mouse_timelapse/GSE228590_gene_annotation.csv.gz")
BASE = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE228nnn/GSE228590/suppl"
RUNS = ["20", "21", "22", "23", "26", "27"]
EMBRYOS = [f"embryo_{i}" for i in range(57, 73)]
NG = 49585
LOG_EVERY = 50_000_000  # matrix lines per progress line


def late_keys():
    import pandas as pd
    keys = pd.read_csv(ROOT / "prep/keys.csv")
    return keys, keys.embryo_id.isin(EMBRYOS).to_numpy()


# ---- stream parser (child process, python -I, stdin only) ---------------------------------------------------------

def _make_parser():
    import numba
    import numpy as np

    @numba.njit(cache=False)
    def parse(buf, codes, S, depth, st):
        n = len(buf)
        i = 0
        v3 = np.zeros(3, np.int64)
        while i < n:
            for k in range(3):
                while i < n and (buf[i] == 32 or buf[i] == 10 or buf[i] == 13):
                    i += 1
                x = 0
                while i < n and buf[i] >= 48 and buf[i] <= 57:
                    x = x * 10 + (buf[i] - 48)
                    i += 1
                if i < n and buf[i] != 32 and buf[i] != 10 and buf[i] != 13:
                    st[4] += 1  # non-integer token
                    while i < n and buf[i] != 32 and buf[i] != 10 and buf[i] != 13:
                        i += 1
                v3[k] = x
            while i < n and (buf[i] == 32 or buf[i] == 10 or buf[i] == 13):
                i += 1
            g, c, v = v3[0], v3[1], v3[2]
            st[0] += 1
            if g > st[1]:
                st[1] = g
            if c > st[2]:
                st[2] = c
            if c < 1 or c > len(codes) or g < 1 or g > S.shape[1]:
                st[5] += 1
                continue
            k = codes[c - 1]
            if k < 0:
                st[3] += 1  # cell of another embryo (or not annotated)
                continue
            depth[k] += v
            S[k, g - 1] += v
    return parse


def stream(r, url, n_keys, test_lines=0):
    """Child: parse the decompressed mtx on stdin, write partial/<r>.npz and .json."""
    import numpy as np
    import scipy.sparse as sp
    codes = np.load(ROOT / f"genome_late/partial/codes.run_{r}.npy")
    parse = _make_parser()
    S = np.zeros((n_keys, NG), np.int64)
    depth = np.zeros(n_keys, np.int64)
    st = np.zeros(6, np.int64)  # lines, max gene, max cell, lines of other cells, non-integer, out of range
    inp = sys.stdin.buffer
    line = inp.readline()
    while line.startswith(b"%"):
        line = inp.readline()
    ng, nc, nz = map(int, line.split())
    rest, nxt, t0 = b"", LOG_EVERY, time.time()
    while True:
        blk = inp.read(1 << 26)
        if not blk:
            break
        blk = rest + blk
        cut = blk.rfind(b"\n") + 1
        rest = blk[cut:]
        parse(np.frombuffer(blk[:cut], np.uint8), codes, S, depth, st)
        if st[0] >= nxt:
            print(f"run_{r}: {st[0]:,} / {nz:,} lines ({time.time() - t0:.0f}s)", flush=True)
            nxt += LOG_EVERY
    if rest.strip():
        parse(np.frombuffer(rest + b"\n", np.uint8), codes, S, depth, st)
    complete = st[0] == nz and ng == NG and nc == len(codes) and st[1] <= ng and st[2] <= nc
    ok = bool(st[4] == 0 and st[5] == 0 and (complete or test_lines))
    meta = dict(run=r, url=url, header=[ng, nc, nz], lines=int(st[0]), max_gene=int(st[1]), max_cell=int(st[2]),
                lines_other_cells=int(st[3]), noninteger_tokens=int(st[4]), out_of_range=int(st[5]),
                cells_annotation=len(codes), cells_in_window=int((codes >= 0).sum()), test_lines=test_lines,
                ok=ok, date=str(datetime.datetime.now(datetime.timezone.utc).date()))
    print(json.dumps(meta), flush=True)
    if not ok:
        sys.exit(f"run_{r}: stream check failed")
    if S.max() >= 2**31:
        sys.exit("count overflow int32")
    out = ROOT / "genome_late/partial"
    sp.save_npz(out / f"{r}.S.npz", sp.csr_matrix(S.astype(np.int32)))
    np.savez(out / f"{r}.vec.npz", depth=depth, ncell=np.bincount(codes[codes >= 0], minlength=n_keys))
    (out / f"{r}.json").write_text(json.dumps(meta))


# ---- parent ---------------------------------------------------------------------------------------------------------

def prepare(runs):
    """Late-window key table and per-run codes remapped to the late keys (-1 = other embryo)."""
    import numpy as np
    keys, late = late_keys()
    lk = keys[late].reset_index(drop=True)
    remap = np.full(len(keys), -1, np.int32)
    remap[np.flatnonzero(late)] = np.arange(len(lk))
    out = ROOT / "genome_late/partial"
    out.mkdir(parents=True, exist_ok=True)
    runs_all = sorted(p.name[10:-4] for p in (ROOT / "prep").glob("codes.run_*.npy"))
    for r in runs_all:  # every other run must hold no cell of the window
        c = np.load(ROOT / f"prep/codes.run_{r}.npy")
        n = int((c >= 0).sum() and (remap[c[c >= 0]] >= 0).sum())
        if r in RUNS:
            np.save(out / f"codes.run_{r}.npy", np.where(c >= 0, remap[np.maximum(c, 0)], -1).astype(np.int32))
            print(f"run_{r}: {n:,} cells of the window", flush=True)
        elif n:
            raise SystemExit(f"run_{r} holds {n} cells of the window; add it to RUNS")
    return lk


def launch(r, n_keys, test_lines):
    url = f"{BASE}/GSE228590_gene_count.run_{r}.mtx.gz"
    head = f" | head -n {test_lines}" if test_lines else ""
    cmd = (f"set -o pipefail; curl -s --retry 5 --retry-all-errors {url} | gzip -dc{head} | "
           f"{sys.executable} -I {Path(__file__).resolve()} stream {r} {url} {n_keys} {test_lines}")
    for attempt in range(1, 4):
        p = subprocess.run(["bash", "-c", cmd], cwd="/", stdout=sys.stdout, stderr=sys.stderr)
        # with head, curl/gzip get SIGPIPE: only the parser's status counts
        if p.returncode == 0 or (test_lines and (ROOT / f"genome_late/partial/{r}.json").exists()):
            return r
        print(f"run_{r}: attempt {attempt} failed (exit {p.returncode})", flush=True)
    raise SystemExit(f"run_{r}: failed after 3 attempts")


def merge(lk, runs):
    import anndata as ad
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp
    p = ROOT / "genome_late/partial"
    metas = [json.loads((p / f"{r}.json").read_text()) for r in runs]
    if not all(m["ok"] for m in metas):
        raise SystemExit("a stream check failed")
    S = sum(sp.load_npz(p / f"{r}.S.npz").astype(np.float64) for r in runs).tocsr()
    depth = sum(np.load(p / f"{r}.vec.npz")["depth"] for r in runs)
    ncell = sum(np.load(p / f"{r}.vec.npz")["ncell"] for r in runs)
    keep = ncell > 0
    lk, S, depth, ncell = lk[keep].reset_index(drop=True), S[keep], depth[keep], ncell[keep]
    ga = pd.read_csv(ANNOT)
    obs = pd.DataFrame({"donor": lk.embryo_id.to_numpy(), "compartment": lk.cell_type.to_numpy(),
                        "trajectory": lk.trajectory.to_numpy(), "n_cells": ncell.astype(int),
                        "depth": depth.astype(float)}, index=lk.key.to_numpy())
    var = pd.DataFrame({"symbol": ga.gene_short_name.astype(str).to_numpy()}, index=ga.gene_ID.astype(str).to_numpy())
    A = ad.AnnData(X=S.astype(np.float32), obs=obs, var=var)
    out = ROOT / "genome_late/pseudobulk_celltype.h5ad"
    A.write_h5ad(out)
    mean = np.asarray(A.X.mean(0)).ravel()
    print(f"obs {A.n_obs}, genes {A.n_vars}, genes with mean count >= 1: {int((mean >= 1).sum())}, "
          f"cells {int(ncell.sum()):,}, size {out.stat().st_size / 1e6:.0f} MB", flush=True)
    sha = hashlib.sha256(out.read_bytes()).hexdigest()
    lines = ["# Genome-wide late-window pseudobulks of GSE228590 (Qiu et al. 2024), A236",
             f"- Built {datetime.date.today()} by `paper/exploratory/A236_embryo_running_clock/scripts/build_genome.py`: "
             "the run matrices were streamed (`curl | gzip -dc | python -I`), parsed only, and not stored.",
             f"- Window: embryo_57 ... embryo_72 ({len(set(obs.donor))} embryos), {A.n_obs} embryo x cell type rows, "
             f"{int(ncell.sum()):,} cells, {A.n_vars} genes (gene annotation GSE228590, `symbol` = gene_short_name; not unique).",
             f"- sha256 of `pseudobulk_celltype.h5ad`: `{sha}`", "- Per run: URL, mtx header (genes, cells, nonzeros), "
             "stream check (lines parsed == nnz, 49,585 genes, cells == annotation, no non-integer token, none out of range):",
             "", "| run | URL | header | lines parsed | cells in annotation | cells in window | check | date |",
             "|---|---|---|---|---|---|---|---|"]
    for m in metas:
        lines.append(f"| {m['run']} | {m['url']} | {' x '.join(map(str, m['header']))} | {m['lines']} | "
                     f"{m['cells_annotation']} | {m['cells_in_window']} | {'ok' if m['ok'] else 'FAILED'} | {m['date']} |")
    (ROOT / "genome_late/PROVENANCE.md").write_text("\n".join(lines) + "\n")


def check():
    """Stored panel table (late window) vs the new table: n_cells, depth and the 67 panel-gene counts, exactly."""
    import anndata as ad
    import numpy as np
    import pandas as pd
    A = ad.read_h5ad(ROOT / "genome_late/pseudobulk_celltype.h5ad")
    old = ad.read_h5ad(ROOT / "pseudobulk_celltype/pseudobulk.h5ad")
    pan = pd.read_csv(ROOT / "prep/panel.csv")
    old = old[old.obs.donor_id.isin(EMBRYOS)]
    new_idx = pd.Series(np.arange(A.n_obs), index=A.obs.donor.astype(str) + "|" + A.obs.compartment.astype(str))
    oi = (old.obs.donor_id.astype(str) + "|" + old.obs.cell_type.astype(str)).to_numpy()
    missing = [k for k in oi if k not in new_idx.index]
    extra = set(new_idx.index) - set(oi)
    ix = new_idx.reindex(oi).dropna().astype(int).to_numpy()
    Xn = A.X[ix][:, pan.row.to_numpy()].toarray()
    Xo = old.X.toarray() if hasattr(old.X, "toarray") else np.asarray(old.X)
    sel = ~pd.Series(oi).isin(missing).to_numpy()
    bad_x = int((Xn != Xo[sel]).any(1).sum())
    bad_n = int((A.obs.n_cells.to_numpy()[ix] != old.obs.n_cells.to_numpy()[sel]).sum())
    bad_d = int((A.obs.depth.to_numpy()[ix] != old.obs.depth.to_numpy()[sel]).sum())
    print(f"panel check: rows compared {len(oi)}, missing in new {len(missing)}, extra in new {len(extra)}; "
          f"mismatches: panel counts {bad_x}, n_cells {bad_n}, depth {bad_d}", flush=True)
    return len(missing) + len(extra) + bad_x + bad_n + bad_d


def main():
    if len(sys.argv) > 1 and sys.argv[1] == "stream":
        r, url, nk, tl = sys.argv[2:6]
        return stream(r, url, int(nk), int(tl))
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--runs", nargs="+", default=RUNS)
    ap.add_argument("--test-lines", type=int, default=0, help="parse only the first N matrix lines, relaxed check, no merge")
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args()
    if a.check:
        sys.exit(1 if check() else 0)
    lk = prepare(a.runs)
    print(f"{len(lk)} late-window keys", flush=True)
    with ThreadPoolExecutor(min(a.jobs, 4)) as ex:
        for r in ex.map(lambda r: launch(r, len(lk), a.test_lines), a.runs):
            print(f"run_{r}: done", flush=True)
    if a.test_lines:
        return
    merge(lk, RUNS)
    if check():
        sys.exit("panel check found mismatches")


if __name__ == "__main__":
    main()
