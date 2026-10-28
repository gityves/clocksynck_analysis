"""A156 stage 1 (blind): parse the GSE48113 series-matrix header and the GPL15331 platform table (text only, python -I).

Writes output/samples.csv, output/probe_map.csv.gz and output/probe_missing.csv.gz (per-probe fraction of non-missing values;
no clock statistic and no gene-level expression is computed).
"""
import csv, gzip
from pathlib import Path

GEO = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian/geo/GSE48113")
OUT = Path(__file__).resolve().parents[1] / "output"; OUT.mkdir(exist_ok=True)
chars, rows, n_probe, keep = [], {}, 0, []
with gzip.open(GEO / "GSE48113_series_matrix.txt.gz", "rt") as f, \
        gzip.open(OUT / "probe_missing.csv.gz", "wt", newline="") as o:
    w = csv.writer(o); w.writerow(["probe", "n_present"])
    on = False
    for line in f:
        if line.startswith("!series_matrix_table_begin"):
            on = True; hdr = next(f).rstrip("\n").split("\t"); continue
        if line.startswith("!series_matrix_table_end"):
            break
        if not on:
            if line.startswith("!Sample_"):
                k, *v = line.rstrip("\n").split("\t"); v = [x.strip('"') for x in v]
                if k == "!Sample_characteristics_ch1":
                    chars.append(v)
                else:
                    rows[k] = v
            continue
        p = line.rstrip("\n").split("\t")
        w.writerow([p[0].strip('"'), sum(1 for x in p[1:] if x not in ("", "null"))]); n_probe += 1
gsm = rows["!Sample_geo_accession"]; meta = {"gsm": gsm, "title": rows["!Sample_title"]}
for c in chars:
    meta[c[0].split(":")[0].strip()] = [x.split(":", 1)[1].strip() for x in c]
with open(OUT / "samples.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(list(meta))
    for i in range(len(gsm)):
        w.writerow([meta[k][i] for k in meta])
with open(GEO / "GPL15331_full.txt", encoding="ascii", errors="replace", newline="") as f, \
        gzip.open(OUT / "probe_map.csv.gz", "wt", newline="") as o:
    w = csv.writer(o); w.writerow(["probe", "control_type", "symbol"]); on = False; h = None
    for line in f:
        line = line.rstrip("\r\n")
        if line.startswith("!platform_table_begin"): on = True; continue
        if line.startswith("!platform_table_end"): break
        if not on: continue
        p = line.split("\t")
        if h is None: h = {k: i for i, k in enumerate(p)}; continue
        p += [""] * (len(h) - len(p))
        w.writerow([p[0], p[h["CONTROL_TYPE"]], p[h["GENE_SYMBOL"]]])
print("samples", len(gsm), "probes in matrix", n_probe)
