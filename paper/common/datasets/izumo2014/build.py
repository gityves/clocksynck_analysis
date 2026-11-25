"""Izumo et al. 2014 eLife (forebrain Bmal1 knockout, SCN clockless): per-mouse spread of the six peripheral organ peak
phases (pituitary, liver, kidney, heart, lung, spleen), parsed from the published source-data tables (Figure 5 source
data b and c: LD, DD; Figure 8 source data: DD with restricted feeding).

The archive is an untrusted download: run with `python -I`; only the docx XML text is read (zipfile), nothing executed.

    python -I paper/common/datasets/izumo2014/build.py
Input: $SCCC_DATA/izumo2014/elife04617_supplementary.zip (PROVENANCE.md). Output: $SCCC_DATA/izumo2014/per_mouse_spread.csv
(condition, genotype control / knockout, mouse, mean_vector_deg, circular_sd_deg, circular_variance).
"""
import csv
import io
import os
import re
import zipfile
from pathlib import Path

DIR = Path(os.environ.get("SCCC_DATA", "/mnt/max8tb/data/CoPhaser_data/cancer_circadian")) / "izumo2014"
TABLES = {"elife04617s001.docx": ("LD", "DD"), "elife04617s002.docx": ("FR",)}
TEXT = re.compile(r"<w:t[^>]*>([^<]*)</w:t>")


def rows(doc: bytes):
    x = zipfile.ZipFile(io.BytesIO(doc)).read("word/document.xml").decode("utf8")
    for r in re.findall(r"<w:tr[ >].*?</w:tr>", x, flags=re.S):
        yield [("".join(TEXT.findall(c))).strip() for c in re.findall(r"<w:tc>.*?</w:tc>", r, flags=re.S)]


def per_mouse(doc: bytes, conditions):
    """Rows of the per-mouse blocks: a header 'LD | Mean Vector ...' then 'WT<i> | ... | BKO<i> | ...'."""
    out, cond = [], None
    for c in rows(doc):
        if len(c) == 8 and c[1].startswith("Mean Vector"):
            cond = c[0]
            continue
        if cond in conditions and len(c) == 8 and c[0].startswith("WT"):
            for geno, cells in (("control", c[:4]), ("knockout", c[4:])):
                if cells[1]:
                    out.append(dict(condition=cond, genotype=geno, mouse=cells[0], mean_vector_deg=float(cells[1]),
                                    circular_sd_deg=float(cells[2]), circular_variance=float(cells[3])))
    return out


def main():
    z = zipfile.ZipFile(DIR / "elife04617_supplementary.zip")
    out = [r for name, conds in TABLES.items() for r in per_mouse(z.read(name), conds)]
    with open(DIR / "per_mouse_spread.csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(out[0]))
        w.writeheader()
        w.writerows(out)
    n = {}
    for r in out:
        n[(r["condition"], r["genotype"])] = n.get((r["condition"], r["genotype"]), 0) + 1
    print(n)


if __name__ == "__main__":
    main()
