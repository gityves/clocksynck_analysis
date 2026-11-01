"""Extract the 21 ER-labelled Pal 2021 tumour samples (10x triplets) from HumanBreast10X.zip into a new empty directory.
Parse only (zipfile, no execution of archive contents).   python -I extract_pal.py"""
import pandas as pd, zipfile
from pathlib import Path
D = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian/sc_breast_er/GSE161529")
DEST = D / "hb10x_tumours"
HERE = Path(__file__).resolve().parents[2]
units = pd.read_csv(HERE / "../scouts_2026-11-01/c4_single_cell/pal_per_patient_malignant.csv")
samples = [p.replace("_", "-") for p in units.patient]
assert len(samples) == 21
DEST.mkdir()  # fails if it exists: must be new
z = zipfile.ZipFile(D / "figshare17058077/HumanBreast10X.zip")
for s in samples:
    for f in ("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz"):
        name = f"Data/{s}/{f}"
        (DEST / s).mkdir(exist_ok=True)
        (DEST / s / f).write_bytes(z.read(name))
print("extracted", len(samples), "samples to", DEST)
