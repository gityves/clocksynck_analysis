"""Download MuTHER (E-TABM-1140) processed zips + sample sheet from the EBI ArrayExpress FTP; PROVENANCE with sha256."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from download import ROOT, DATE, fetch, sha
B = "https://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/TABM/E-TABM-1140/"
d = f"{ROOT}/living_multitissue/E-TABM-1140"; os.makedirs(d, exist_ok=True)
names = ["E-TABM-1140.README.txt", "E-TABM-1140.idf.txt", "E-TABM-1140.sdrf.txt", "E-TABM-1140.processed.2.zip", "E-TABM-1140.processed.3.zip", "E-TABM-1140.processed.4.zip"]
rows = []
for n in names:
    p = f"{d}/{n}"
    if n.endswith(".txt") and os.path.exists(p) and os.path.getsize(p) == 0: os.remove(p)
    fetch(B + n, p); rows.append((n, B + n, os.path.getsize(p), sha(p)))
with open(f"{d}/PROVENANCE.md", "w") as f:
    f.write(f"# E-TABM-1140 (MuTHER)\n\nSource: ArrayExpress E-TABM-1140. Downloaded {DATE}. Processed (normalised) Illumina HT-12 matrices as zip, left zipped; raw intensity files not downloaded.\n\n| file | URL | bytes | sha256 |\n|---|---|---|---|\n")
    for r in rows: f.write(f"| {r[0]} | {r[1]} | {r[2]} | {r[3]} |\n")
print("done", [(r[0], r[2]) for r in rows])
