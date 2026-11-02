"""Qiu et al. 2024 (Nature): mouse embryo time-lapse atlas E8 to P0, sci-RNA-seq3, 74 embryos, 11.4 M cells.

Source:   GEO GSE228590: 17 per-run MatrixMarket count matrices (genes x cells, 49,585 genes; `.mtx.gz`, about 600 GB of
          text in total, NOT stored: streamed once, `curl | gzip -dc`), per-run cell annotation (`ca.run_*.csv.gz`),
          gene annotation, sample sheet and `df_cell.csv` (cell type, trajectory, embryo, E-day of every cell), in
          `$SCCC_DATA/natexp/qiu2024/` and `$SCCC_DATA/scout_natexp/mouse_timelapse/`.
Replaces: legacy/desync/validation/natural_experiments/scripts/build_qiu.py and fix_qiu_var.py, build_qiu.md, run_qiu.sh.
Used by:  A108, A108b, A108c (embryonic clock emergence, time-of-day following, lineage gates).

Pipeline (three steps; the first and the last are cheap)
1. `prepare()`: cell -> (embryo, cell type) key codes per run from `df_cell.csv` and the run annotations, and the
   matrix rows of the 14 clock genes (mouse symbols matched case-insensitively; aliases BMAL1 to Arntl, CIART to
   Ciart / Gm129 / Chrono);
2. `stream_run(mtx_stream, codes, rows, n_keys)` per run: sums the clock-gene counts, the genome-wide depth and the
   number of non-zero entries of every key, from a decompressed MatrixMarket text stream (numba parser; refuses to
   return unless lines, header, indices and integer tokens check); `save_partial` writes the small result;
3. `build(level)`: adds the partial results and returns the standard table. Without `partials=` it reads the
   partials of the legacy run (`$SCCC_DATA/natexp/qiu2024/partial/*.npz`, 67 panel columns, selected to the 14 clock genes
   through `prep/panel.csv`), so no matrix is streamed again.
No feature table is written to disk, so the empty `feature_name` that `fix_qiu_var.py` repaired cannot occur.

Choices (as the legacy builder)
- unit (donor) = embryo (`embryo_id`); key = embryo x `celltype_update`; level "celltype" (default) or "trajectory" (sum of
  the cell types of each of the 26 `major_trajectory`s; every cell type belongs to one);
- runs 4, 13-16, 17_sub1, 17_sub2, 18-27 (runs 28 and 29, a birth-mode experiment, are absent from `df_cell`);
- stratum: `batch`, the sci-RNA-seq3 run(s) of the embryo (one value per embryo; analyses use it as stratum only in
  windows that hold several batches); every pseudobulk with at least one cell is kept (analyses use >= 50 cells and
  compartments present in >= 75 % of the embryos of a window); extra column `nnz` (non-zero entries).

Donor table (embryo): `day` (E-day; P0 = 19.5, a label only), `day_label`, `batch`, `t_clock` = 24 x frac(E-day) h
(the clock time of the sample of an embryo at E15.0-E18.75, one sample per embryo; NaN for P0; not interpretable before
E15), sample-sheet columns (RT_group, somite_count, Plug_Date, Date_Harvested, litter_ID, shipment_ID), `n_cells`.

Clock time checked against the Qiu 2024 Methods (PMC10901739, 2026-11-02): plugs checked in the morning, noon of the plug
day = E0.5, and collection at E15.0-E18.75 "was performed precisely at 00:00, 06:00, 12:00 and 18:00", so frac 0, .25,
.5, .75 = 00, 06, 12, 18 h local time (not ZT: the light schedule is not given). The sample sheet agrees: E(n).0 is
harvested on the date before E(n).25-.75 of the same plug-to-harvest span, i.e. at midnight. E10.25-E14.75 are staged
morphologically (eMOSS), so their fraction is not a clock time.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numba
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve()
sys.path[:0] = [str(HERE.parents[1]), str(Path(__file__).resolve().parents[3])]
from clocksync.genes import canonical  # noqa: E402
from common.paths import DATA  # noqa: E402
from natexp_common import GENES, compare, read_legacy, standard  # noqa: E402

ROOT = DATA / "natexp" / "qiu2024"
SCOUT = DATA / "scout_natexp" / "mouse_timelapse"
LEVELS = ("celltype", "trajectory")
N_GENES = 49585
ALIAS = {"ARNTL": ("Arntl",), "CIART": ("Ciart", "Gm129", "Chrono")}
SHEET_COLS = ["RT_group", "somite_count", "Plug_Date", "Date_Harvested", "litter_ID", "shipment_ID"]


def day_num(label: str) -> float:
    return 19.5 if label == "P0" else float(label[1:])  # P0 = noon of the day of birth (label only)


# --- 1. keys, codes, gene rows ---------------------------------------------------------------------------------

def cell_table() -> pd.DataFrame:
    cols = ["cell_id", "day", "embryo_id", "experimental_batch", "celltype_update", "major_trajectory"]
    c = pd.read_csv(SCOUT / "df_cell.csv", usecols=cols).rename(
        columns={"celltype_update": "cell_type", "major_trajectory": "trajectory"})
    c["key"] = c.embryo_id + "|" + c.cell_type
    return c


def key_table(cell: pd.DataFrame) -> pd.DataFrame:
    keys = pd.DataFrame({"key": sorted(cell.key.unique())})
    keys["embryo_id"] = keys.key.str.split("|").str[0]
    keys["cell_type"] = keys.key.str.split("|", n=1).str[1]
    keys["trajectory"] = keys.cell_type.map(cell.drop_duplicates("cell_type").set_index("cell_type").trajectory)
    return keys


def run_codes(cell: pd.DataFrame, keys: pd.DataFrame, annot_dir: Path = ROOT / "raw_annot") -> dict:
    """Key code of every cell (matrix column) of every run, -1 for cells absent from `df_cell`."""
    code = pd.Series(np.arange(len(keys)), index=keys.key)
    cid = pd.Series(code.reindex(cell.key).to_numpy(), index=cell.cell_id.to_numpy())
    out = {}
    for f in sorted(annot_dir.glob("ca.run_*.csv.gz")):
        ids = pd.read_csv(f, usecols=["cell_id"]).cell_id
        c = cid.reindex(ids).to_numpy()
        out[f.name[3:-7]] = np.where(np.isnan(c), -1, c).astype(np.int32)
    return out


def gene_rows(annotation: Path = SCOUT / "GSE228590_gene_annotation.csv.gz") -> np.ndarray:
    """Matrix row of each of the 14 clock genes (first match of the mouse symbol, case-insensitive, then aliases)."""
    sym = pd.read_csv(annotation).gene_short_name.astype(str)
    low = {}
    for i, s in enumerate(sym):
        low.setdefault(s.lower(), i)
    rows = []
    for g in GENES:
        hit = next((low[s.lower()] for s in (g.capitalize(),) + ALIAS.get(g, ()) if s.lower() in low), None)
        if hit is None:
            raise KeyError(f"clock gene {g} not found in the gene annotation")
        rows.append(hit)
    return np.array(rows)


# --- 2. one run, streamed ------------------------------------------------------------------------------------------

@numba.njit(cache=True)
def _parse(buf, codes, rowmap, S, depth, nnz, st):
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
        if c < 1 or c > len(codes) or g < 1 or g > len(rowmap):
            st[5] += 1
            continue
        k = codes[c - 1]
        if k < 0:
            st[3] += 1
            continue
        depth[k] += v
        nnz[k] += 1
        j = rowmap[g - 1]
        if j >= 0:
            S[k, j] += v


def stream_run(stream, codes: np.ndarray, rows: np.ndarray, n_keys: int, n_genes: int = N_GENES) -> dict:
    """Sum one run from a binary text stream of a MatrixMarket file (genes x cells, coordinate format; e.g.
    `gzip.open(...)`, or `curl | gzip -dc` on stdin). Returns S (keys x 14), depth, nnz and ncell per key plus the
    stream check; raises if lines parsed differ from the header, an index is out of range or a token is not an integer."""
    rowmap = np.full(n_genes, -1, np.int32)
    rowmap[rows] = np.arange(len(rows))
    S = np.zeros((n_keys, len(rows)), np.int64)
    depth, nnz = np.zeros(n_keys, np.int64), np.zeros(n_keys, np.int64)
    st = np.zeros(6, np.int64)  # lines, max gene, max cell, lines of unmapped cells, non-integer tokens, out of range
    line = stream.readline()
    while line.startswith(b"%"):
        line = stream.readline()
    ng, nc, nz = map(int, line.split())
    rest = b""
    while True:
        blk = stream.read(1 << 26)
        if not blk:
            break
        blk = rest + blk
        cut = blk.rfind(b"\n") + 1
        rest = blk[cut:]
        _parse(np.frombuffer(blk[:cut], np.uint8), codes, rowmap, S, depth, nnz, st)
    if rest.strip():
        _parse(np.frombuffer(rest + b"\n", np.uint8), codes, rowmap, S, depth, nnz, st)
    check = dict(header=[ng, nc, nz], lines=int(st[0]), noninteger_tokens=int(st[4]), out_of_range=int(st[5]),
                 cells_annotation=len(codes), cells_mapped=int((codes >= 0).sum()))
    if not (st[0] == nz and ng == n_genes and nc == len(codes) and st[4] == 0 and st[5] == 0):
        raise ValueError(f"stream check failed: {check}")
    return dict(S=S, depth=depth, nnz=nnz, ncell=np.bincount(codes[codes >= 0], minlength=n_keys), check=check)


def save_partial(path: Path, part: dict) -> None:
    np.savez(path.with_suffix(".npz"), S=part["S"], depth=part["depth"], nnz=part["nnz"], ncell=part["ncell"],
             genes=np.array(GENES))
    path.with_suffix(".json").write_text(json.dumps(part["check"]))


# --- 3. merge -----------------------------------------------------------------------------------------------------

def read_partials(directory: Path = ROOT / "partial", panel: Path = ROOT / "prep" / "panel.csv") -> dict:
    """Sum the partial results of a directory; columns are selected / ordered to the 14 clock genes (partials written
    by `save_partial` carry `genes`; those of the legacy run have the 67 columns of `panel`)."""
    tot = None
    for p in sorted(Path(directory).glob("*.npz")):
        z = np.load(p)
        if "genes" in z:
            cols = [list(z["genes"]).index(g) for g in GENES]
        else:
            names = [canonical(g) for g in pd.read_csv(panel).feature_name]
            cols = [names.index(g) for g in GENES]
        part = dict(S=z["S"][:, cols], depth=z["depth"], nnz=z["nnz"], ncell=z["ncell"])
        tot = part if tot is None else {k: tot[k] + part[k] for k in part}
    if tot is None:
        raise FileNotFoundError(f"no partial results in {directory}")
    return tot


def embryo_table(cell: pd.DataFrame, n_cells: pd.Series) -> pd.DataFrame:
    first = cell.drop_duplicates("embryo_id").set_index("embryo_id")
    runs = cell.groupby("embryo_id").experimental_batch.agg(lambda s: ";".join(sorted(s.unique())))
    ss = pd.read_csv(SCOUT / "sample_sheet.txt", sep="\t").set_index("embryo_id")
    emb = pd.DataFrame(index=sorted(cell.embryo_id.unique()))
    emb["day_label"] = first.day.reindex(emb.index)
    emb["day"] = emb.day_label.map(day_num)
    emb["batch"] = runs.reindex(emb.index)
    emb["t_clock"] = np.where(emb.day_label == "P0", np.nan, 24 * (emb.day % 1))
    for c in SHEET_COLS:
        emb[c] = ss[c].reindex(emb.index)
    emb["n_cells"] = n_cells.reindex(emb.index)
    return emb.rename_axis("donor")


def build(level: str = "celltype", partials: Path | None = None) -> tuple[pd.DataFrame, pd.DataFrame]:
    """(pseudobulk, donors) from the partial results (see the module docstring)."""
    cell = cell_table()
    keys = key_table(cell)
    part = read_partials(partials or ROOT / "partial")
    if len(part["ncell"]) != len(keys):
        raise ValueError("partial results and key table differ in length")
    keep = part["ncell"] > 0
    t = keys[keep].reset_index(drop=True)
    t["n_cells"] = part["ncell"][keep]
    t["depth"] = part["depth"][keep].astype(float)
    t["nnz"] = part["nnz"][keep]
    t[GENES] = part["S"][keep].astype(float)
    if level == "trajectory":
        t = t.groupby(["embryo_id", "trajectory"], as_index=False)[["n_cells", "depth", "nnz"] + GENES].sum()
        t = t.rename(columns={"trajectory": "compartment"})
    else:
        t = t.rename(columns={"cell_type": "compartment"})
    emb = embryo_table(cell, pd.Series(part["ncell"][keep], index=keys.embryo_id[keep]).groupby(level=0).sum())
    pb = t.rename(columns={"embryo_id": "donor"})
    pb["stratum"] = pb.donor.map(emb.batch)
    return standard(pb, extra=("nnz",)), emb


def from_legacy(level: str = "celltype") -> tuple[pd.DataFrame, pd.DataFrame]:
    """Standard table from the legacy `pseudobulk_<level>/pseudobulk.h5ad` and `samples.csv`. The legacy var table
    has an empty `feature_name` (repaired once by `fix_qiu_var.py`); only the var index is used here."""
    d = ROOT / f"pseudobulk_{level}"
    pb = read_legacy(d / "pseudobulk.h5ad", extra=("batch", "nnz"))
    pb["stratum"] = pb.pop("batch")
    emb = pd.read_csv(d / "samples.csv").set_index("embryo_id").rename_axis("donor")
    emb["t_clock"] = np.where(emb.day_label == "P0", np.nan, 24 * (emb.day % 1))
    return standard(pb, extra=("nnz",)), emb


if __name__ == "__main__":
    import io
    import tempfile

    for lev in LEVELS:
        a, da = build(lev)
        b, db = from_legacy(lev)
        print(lev, "full table vs legacy:", compare(a, b), "| nnz equal:", bool((a.merge(b, on=["donor", "compartment"]).pipe(
            lambda m: (m.nnz_x == m.nnz_y).all()))), "| donors:", len(da), "| embryo n_cells equal:",
            bool((da.n_cells.astype(int) == db.n_cells.reindex(da.index).astype(int)).all()))
    leg = pd.read_csv(ROOT / "prep" / "panel.csv").set_index("feature_name")
    print("gene rows equal legacy panel:", bool((gene_rows() == leg.loc[["BMAL1" if g == "ARNTL" else g for g in GENES], "row"].to_numpy()).all()),
          "| keys equal legacy:", bool(key_table(cell_table()).key.equals(pd.read_csv(ROOT / "prep" / "keys.csv").key)))
    # stream_run on a small synthetic MatrixMarket stream (the real runs are not stored)
    rng = np.random.default_rng(1)
    ng, nc, nk = 60, 500, 7
    rows = np.sort(rng.choice(ng, len(GENES), replace=False))
    M = rng.poisson(0.7, (ng, nc))
    codes = rng.integers(-1, nk, nc).astype(np.int32)
    g, c = np.nonzero(M)
    text = f"%%MatrixMarket matrix coordinate integer general\n{ng} {nc} {len(g)}\n" + \
        "".join(f"{i + 1} {j + 1} {M[i, j]}\n" for i, j in zip(g, c))
    part = stream_run(io.BytesIO(text.encode()), codes, rows, nk, n_genes=ng)
    ref_S = np.zeros((nk, len(rows)), np.int64)
    for k in range(nk):
        ref_S[k] = M[rows][:, codes == k].sum(1)
    ok = (part["S"] == ref_S).all() and (part["depth"] == [M[:, codes == k].sum() for k in range(nk)]).all()
    with tempfile.TemporaryDirectory() as t:
        save_partial(Path(t) / "x", part)
        back = read_partials(Path(t))
    print("stream_run on synthetic mtx:", "ok" if ok and (back["S"] == ref_S).all() else "FAILED")
