"""Stream the count matrix and cell metadata out of a large Seurat `.rds` without R and in bounded memory
(GSE266265 endometrial atlas, 42 GB gzip; amendment A11 part 3). Self-contained (run with `python -I`; the input is
downloaded data). The file is parsed as R serialization format (XDR, version 2 or 3); nothing in it is executed.

An `.rds` has no index, so it cannot be opened by patient: `readRDS` needs the whole object in memory. This reader
walks the stream once:
- the `i` (row index) and `x` (value) vectors of the counts matrix of assay `--assay-index` (default the first) are
  written to disk as they stream past (`counts_i.int32`, `counts_x.float32`), with a running check of integrality;
- every other large vector (> `BIG` elements: normalised data, scaled data, other assays) is read and discarded;
- small objects (dimensions, gene and cell names, column pointers, `meta.data`) are kept and written out.

Output (`OUT_DIR`): `counts_i.int32`, `counts_x.float32`, `counts_p.npy` (column pointers, cells are columns),
`genes.txt`, `cells.txt`, `meta.csv` (Seurat `meta.data`), `assays.txt`, `structure.txt` (slot tree with lengths),
`x_check.json` (non-integer values, max, sum).

    gzip -dc FILE.rds | python -I stream_rds.py OUT_DIR [--assay-index 0]
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

BIG = 20_000_000          # elements; larger vectors are spilled (counts) or discarded
CHUNK = 1 << 24           # elements per read of a large vector
NIL, GLOBALENV, UNBOUND, MISSING, BASENS, NSSXP, PKGSXP, PERSIST = 254, 253, 252, 251, 250, 249, 248, 247
EMPTYENV, BASEENV, ATTRLANG, ATTRLIST, ALTREP, REF = 242, 241, 240, 239, 238, 255
DTYPE = {10: ">i4", 13: ">i4", 14: ">f8", 15: ">c16", 24: "u1"}


class R:
    """A parsed R object: SEXP type, value (list, numpy array, str, list of (tag, value) for pairlists), attributes."""
    __slots__ = ("type", "val", "attr")

    def __init__(self, type_, val=None, attr=None):
        self.type, self.val, self.attr = type_, val, attr

    def a(self, name):
        return None if self.attr is None else self.attr.get(name)


class Reader:
    def __init__(self, stream, out: Path, assay_index: int):
        self.f, self.out, self.assay_index = stream, out, assay_index
        self.refs, self.path, self.spilled, self.check = [], [], {}, dict(nonint=0, max=0.0, sum=0.0)

    # -- primitives
    def read(self, n: int) -> bytes:
        b = self.f.read(n)
        while len(b) < n:  # pipes can return short reads
            more = self.f.read(n - len(b))
            if not more:
                raise EOFError("stream ended early")
            b += more
        return b

    def int(self) -> int:
        return int.from_bytes(self.read(4), "big", signed=True)

    def length(self) -> int:
        n = self.int()
        if n == -1:
            hi, lo = self.int(), self.int()
            n = (hi << 32) + lo
        return n

    def string(self) -> str | None:
        n = self.int()
        return None if n == -1 else self.read(n).decode("utf-8", "replace")

    # -- large vectors
    def spill_target(self, t: int) -> str | None:
        """File name if this vector is the i or x slot of the chosen assay's counts matrix (Seurat v3/v4 `counts`
        slot or v5 `layers$counts`)."""
        p = self.path
        if len(p) < 3 or p[0] != "assays" or p[1] != f"[{self.assay_index}]" or "counts" not in p:
            return None
        return {("i", 13): "counts_i.int32", ("x", 14): "counts_x.float32"}.get((p[-1], t))

    def big_vector(self, t: int, n: int) -> R:
        dt = np.dtype(DTYPE[t])
        name = self.spill_target(t)
        fh = open(self.out / name, "wb") if name else None
        left = n
        while left:
            k = min(left, CHUNK)
            v = np.frombuffer(self.read(k * dt.itemsize), dt)
            if fh is not None:
                if t == 14:
                    self.check["nonint"] += int(np.count_nonzero(v != np.round(v)))
                    self.check["max"] = max(self.check["max"], float(v.max()))
                    self.check["sum"] += float(v.sum())
                    fh.write(v.astype("<f4").tobytes())
                else:
                    fh.write(v.astype("<i4").tobytes())
            left -= k
        if fh is not None:
            fh.close()
            self.spilled["/".join(self.path)] = (name, n)
        print(("spilled " if name else "skipped ") + "/".join(self.path), n, flush=True)
        return R(t, ("<large>", n, name))

    # -- items
    def attributes(self, has_attr: bool):
        if not has_attr:
            return None
        pl = self.item()
        return {tag: v for tag, v in pl.val} if pl is not None and pl.type == 2 else None

    def pairlist(self, t: int, flags: int) -> R:
        items, attr = [], None
        while True:
            has_attr, has_tag = flags & (1 << 9), flags & (1 << 10)
            a = self.attributes(has_attr)
            attr = attr or a
            tag = None
            if has_tag:
                s = self.item()
                tag = s.val if isinstance(s, R) else None
            self.path.append(tag or f"[{len(items)}]")
            car = self.item()
            self.path.pop()
            items.append((tag, car))
            flags = self.int()
            nt = flags & 0xFF
            if nt == NIL:
                return R(2 if t in (2, ATTRLIST) else t, items, attr)
            if nt not in (2, ATTRLIST):  # improper list tail
                self.dispatch(flags)
                return R(t, items, attr)

    def item(self):
        return self.dispatch(self.int())

    def dispatch(self, flags: int):
        t = flags & 0xFF
        has_attr = bool(flags & (1 << 9))
        if t == NIL:
            return None
        if t in (GLOBALENV, UNBOUND, MISSING, BASENS, EMPTYENV, BASEENV):
            return R(t)
        if t == REF:
            k = (flags >> 8) or self.int()
            return self.refs[k - 1]
        if t in (NSSXP, PKGSXP, PERSIST):
            assert self.int() == 0
            names = [self.item().val for _ in range(self.int())]
            r = R(t, names)
            self.refs.append(r)
            return r
        if t == 1:  # SYMSXP
            r = R(1, self.item().val)
            self.refs.append(r)
            return r
        if t == 4:  # ENVSXP
            locked = self.int()
            r = R(4, None)
            self.refs.append(r)
            enclos, frame, hashtab, attr = self.item(), self.item(), self.item(), self.item()
            r.val = dict(locked=locked, frame=frame, hashtab=hashtab)
            return r
        if t in (2, 3, 5, 6, 17, ATTRLANG, ATTRLIST):
            return self.pairlist(t, flags)
        if t == 9:  # CHARSXP
            return R(9, self.string())
        if t in (7, 8):
            return R(t, self.read(self.int()).decode())
        if t in DTYPE:
            n = self.length()
            if n > BIG:
                r = self.big_vector(t, n)
            else:
                r = R(t, np.frombuffer(self.read(n * np.dtype(DTYPE[t]).itemsize), DTYPE[t]).copy())
            r.attr = self.attributes(has_attr)
            return r
        if t == 16:  # STRSXP
            n = self.length()
            v = [self.item().val for _ in range(n)]
            return R(16, v, self.attributes(has_attr))
        if t in (19, 20):  # VECSXP, EXPRSXP
            n = self.length()
            v = []
            for k in range(n):
                self.path.append(f"[{k}]")
                v.append(self.item())
                self.path.pop()
            return R(t, v, self.attributes(has_attr))
        if t == 25:  # S4SXP
            return R(25, None, self.attributes(has_attr))
        if t == 22:  # EXTPTRSXP
            r = R(22)
            self.refs.append(r)
            self.item(); self.item()
            r.attr = self.attributes(has_attr)
            return r
        if t == 23:  # WEAKREFSXP
            r = R(23)
            self.refs.append(r)
            r.attr = self.attributes(has_attr)
            return r
        if t == ALTREP:
            info, state, attr = self.item(), self.item(), self.item()
            return R(ALTREP, (info, state), {tag: v for tag, v in attr.val} if attr is not None else None)
        raise ValueError(f"unsupported SEXP type {t} at {'/'.join(map(str, self.path))}")

    def header(self):
        fmt = self.read(2)
        if fmt != b"X\n":
            raise ValueError(f"not XDR serialization: {fmt!r}")
        version = self.int()
        self.int(); self.int()
        if version == 3:
            self.read(self.int())  # native encoding


def names_of(r: R) -> list:
    n = r.a("names")
    return list(n.val) if n is not None else []


def compact_or_vector(r: R) -> np.ndarray | list:
    """Values of a vector, expanding ALTREP compact sequences and wrappers."""
    if r.type == ALTREP:
        cls = r.val[0].val[0][1].val  # info: pairlist (class sym, package sym, type)
        state = r.val[1]
        if cls in ("compact_intseq", "compact_realseq"):
            n, start, step = state.val[:3]
            return start + step * np.arange(int(n))
        if cls.startswith("wrap_"):
            return compact_or_vector(state.val[0][1] if state.type == 2 else state.val[0])
        if cls == "deferred_string":
            return list(compact_or_vector(state.val[0][1] if state.type == 2 else state.val[0]).astype(str))
        raise ValueError(f"ALTREP class {cls}")
    return r.val


def column(r: R):
    """A meta.data column: factors become their labels."""
    v = compact_or_vector(r)
    lev = r.a("levels")
    if lev is not None:
        lv = np.array(lev.val, dtype=object)
        v = np.asarray(v)
        return np.where(v > 0, lv[np.clip(v - 1, 0, len(lv) - 1)], None)
    return v


def tree(r, depth=0, name="", lines=None, max_depth=6):
    lines = [] if lines is None else lines
    if depth > max_depth or not isinstance(r, R):
        return lines
    n = len(r.val) if isinstance(r.val, (list, np.ndarray)) else (r.val[1] if isinstance(r.val, tuple) and
                                                                  r.val[0] == "<large>" else "")
    cls = r.a("class")
    lines.append(f"{'  ' * depth}{name} type={r.type} len={n} class={cls.val if cls is not None else ''}")
    if r.attr:
        for k, v in r.attr.items():
            if k not in ("names", "row.names", "class"):
                tree(v, depth + 1, f"@{k}", lines, max_depth)
    if r.type == 19:
        nm = names_of(r)
        for k, v in enumerate(r.val[:60]):
            tree(v, depth + 1, f"${nm[k] if k < len(nm) else k}", lines, max_depth)
    return lines


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("out_dir")
    ap.add_argument("--assay-index", type=int, default=0)
    args = ap.parse_args(argv)
    out = Path(args.out_dir)
    out.mkdir(parents=True, exist_ok=True)
    rd = Reader(sys.stdin.buffer, out, args.assay_index)
    rd.header()
    obj = rd.item()
    (out / "structure.txt").write_text("\n".join(tree(obj)) + "\n")
    assays = obj.a("assays")
    (out / "assays.txt").write_text("\n".join(map(str, names_of(assays))) + "\n")
    a = assays.val[args.assay_index]
    counts = a.a("counts")
    if counts is None:  # Seurat v5 Assay5: layers$counts
        layers = a.a("layers")
        counts = layers.val[names_of(layers).index("counts")]
    dim = counts.a("Dim").val
    dn = counts.a("Dimnames").val
    genes, cells = (list(compact_or_vector(d)) if d is not None else [] for d in dn)
    if not cells:  # v5: names held by the assay's cells / features maps
        raise ValueError("counts matrix has no column names; see structure.txt")
    np.save(out / "counts_p.npy", np.asarray(counts.a("p").val, dtype=np.int64))
    (out / "genes.txt").write_text("\n".join(genes) + "\n")
    (out / "cells.txt").write_text("\n".join(cells) + "\n")
    md = obj.a("meta.data")
    cols = {k: column(v) for k, v in zip(names_of(md), md.val)}
    rn = compact_or_vector(md.a("row.names"))
    meta = pd.DataFrame(cols, index=pd.Index(rn if len(rn) == len(next(iter(cols.values()))) else None, name="cell"))
    meta.to_csv(out / "meta.csv")
    check = dict(rd.check, dim=[int(x) for x in dim], nnz=rd.spilled.get("assays/[%d]/counts/x" % args.assay_index,
                                                                         (None, None))[1], spilled=rd.spilled)
    (out / "x_check.json").write_text(json.dumps(check, indent=1))
    print(json.dumps(check, indent=1), meta.shape, flush=True)


if __name__ == "__main__":
    main()
