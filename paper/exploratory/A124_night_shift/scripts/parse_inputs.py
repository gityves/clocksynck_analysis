"""A124 stage 1 (blind): parse the GSE106800 series matrix header and the GPL17692 platform table (text only, run with python -I).

Writes to output/: samples.csv (metadata), probe_map.csv.gz (probe -> gene symbols). No expression value is read here.
"""
import csv, gzip, re, sys
from pathlib import Path

GEO = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian/geo/GSE106800")
OUT = Path(__file__).resolve().parents[1] / "output"
OUT.mkdir(exist_ok=True)

rows, chars = {}, []
with gzip.open(GEO / "GSE106800_series_matrix.txt.gz", "rt") as f:
    for line in f:
        if line.startswith("!series_matrix_table_begin"):
            break
        if line.startswith("!Sample_"):
            k, *v = line.rstrip("\n").split("\t")
            v = [x.strip('"') for x in v]
            if k == "!Sample_characteristics_ch1":
                chars.append(v)
            else:
                rows[k] = v
gsm, title = rows["!Sample_geo_accession"], rows["!Sample_title"]
meta = {"gsm": gsm, "title": title, "supp": rows["!Sample_supplementary_file"]}
for c in chars:
    key = c[0].split(":")[0].strip()
    meta[key] = [x.split(":", 1)[1].strip() if ":" in x else x for x in c]
with open(OUT / "samples.csv", "w", newline="") as f:
    w = csv.writer(f); keys = list(meta); w.writerow(keys)
    for i in range(len(gsm)):
        w.writerow([meta[k][i] for k in keys])

# platform table: ID, gene_assignment, category
pat = re.compile(r"^(\S+) // (\S+) // ")
n = 0
with open(GEO / "GPL17692_full.txt", encoding="ascii", errors="replace", newline="") as f, \
        gzip.open(OUT / "probe_map.csv.gz", "wt", newline="") as o:
    w = csv.writer(o); w.writerow(["probe", "category", "symbols_all", "n_total_probes"])
    on = False; hdr = None
    for line in f:
        line = line.rstrip("\r\n")
        if line.startswith("!platform_table_begin"):
            on = True; hdr = None; continue
        if line.startswith("!platform_table_end"):
            break
        if not on:
            continue
        p = line.split("\t")
        if hdr is None:
            hdr = {k: i for i, k in enumerate(p)}; continue
        ga = p[hdr["gene_assignment"]]
        syms = []
        for part in ga.split(" /// "):
            q = part.split(" // ")
            if len(q) > 1 and q[1] != "---":
                syms.append(q[1])
        w.writerow([p[0], p[hdr["category"]], ";".join(dict.fromkeys(syms)), p[hdr["total_probes"]]])
        n += 1
print("samples", len(gsm), "probes", n)
