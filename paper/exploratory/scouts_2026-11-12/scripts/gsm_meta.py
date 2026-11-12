"""Parse GEO series-matrix sample metadata (and optionally the expression table)."""
import gzip, io, re
import pandas as pd
def read_meta(path):
    meta, tbl = {}, []
    inside = False
    with gzip.open(path, "rt", errors="replace") as fh:
        for line in fh:
            if line.startswith("!series_matrix_table_begin"): inside = True; continue
            if line.startswith("!series_matrix_table_end"): inside = False; continue
            if inside: tbl.append(line); continue
            if line.startswith("!Sample_"):
                k, *v = line.rstrip("\n").split("\t")
                v = [x.strip('"') for x in v]
                k = k[len("!Sample_"):]
                if k == "characteristics_ch1":
                    for i, x in enumerate(v):
                        m = re.match(r"([^:]+):\s*(.*)", x)
                        key, val = (m.group(1).strip(), m.group(2)) if m else (f"char{i}", x)
                        meta.setdefault("c_" + key, [None] * len(v))[i] = val
                else:
                    meta.setdefault(k, v)
    df = pd.DataFrame(meta)
    expr = None
    if tbl:
        expr = pd.read_csv(io.StringIO("".join(tbl)), sep="\t", index_col=0, quotechar='"')
    return df, expr
