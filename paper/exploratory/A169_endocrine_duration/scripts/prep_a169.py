"""A169 data preparation (blind to the clock statistic): GSE111563, GSE59515 and GSE55374 series matrices parsed with the A168 parser
(`A168_endocrine/scripts/prep_a168.py`: log2 check, probes averaged per symbol, A97 / oestrogen scores) into ../output/.
Time points are coded 0 (pre), 1 (early: 2 weeks or tp2, 0-120 d) and 2 (late: 3 months or tp4, > 120 d).
Run: python -I prep_a169.py   (reads downloaded text files only; nothing downloaded is executed)"""
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1] / "A168_endocrine" / "scripts"))
import pandas as pd  # noqa: E402

import a168_lib as L  # noqa: E402
import prep_a168 as P  # noqa: E402

OUT = HERE.parent / "output"
OUT.mkdir(exist_ok=True)
L.OUT = OUT  # prep_a168.finish writes to L.OUT
P.PLATFORM.update({"GSE111563": "GPL10558", "GSE59515": "GPL10558", "GSE55374": "GPL10558"})


def char(h, key):
    """Values of the characteristics line starting with `key:` (one per sample)."""
    for r in h["!Sample_characteristics_ch1"]:
        if r[0].startswith(key + ":"):
            return [x.split(":", 1)[1].strip() for x in r]
    raise KeyError(key)


def table(gse, rows):
    return pd.DataFrame(rows, columns=["gsm", "patient", "arm", "tp", "response", "title"]).set_index("gsm")


def gse111563():
    gse = "GSE111563"; h = P.header(P.D / gse / f"{gse}_series_matrix.txt.gz")
    rows = []
    for g, t in zip(h["!Sample_geo_accession"][0], h["!Sample_title"][0]):
        p, k = re.fullmatch(r"Breasttumour_Patient(\w+)_Timepoint(\d)", t).groups()
        rows.append((g, p, "LET", {"1": 0, "2": 1, "4": 2}[k], "", t))
    P.finish(gse, "gse111563", table(gse, rows), P.expression(gse))


def gse59515():
    gse = "GSE59515"; h = P.header(P.D / gse / f"{gse}_series_matrix.txt.gz")
    tp = {"pre-treatment": 0, "2 wks of letrozole treatment": 1, "3 months of letrozole treatment": 2}
    rows = list(zip(h["!Sample_geo_accession"][0], char(h, "patient id"), ["LET"] * 75, [tp[x] for x in char(h, "time point")],
                    char(h, "clinical response"), h["!Sample_title"][0]))
    P.finish(gse, "gse59515", table(gse, rows), P.expression(gse))


def gse55374():
    gse = "GSE55374"; h = P.header(P.D / gse / f"{gse}_series_matrix.txt.gz")
    tp = {"pre-treatment": 0, "2 wk endocrine therapy": 1, "3 mo endocrine therapy": 2}
    rows = list(zip(h["!Sample_geo_accession"][0], [s.replace("patient ", "") for s in char(h, "subject")], ["LET"] * 36,
                    [tp[x] for x in char(h, "treatment")], char(h, "clinical response"), h["!Sample_title"][0]))
    P.finish(gse, "gse55374", table(gse, rows), P.expression(gse))


if __name__ == "__main__":
    gse111563(); gse59515(); gse55374()
    pd.DataFrame(P.LOG, columns=["series", "item", "value"]).to_csv(OUT / "prep_log.csv", index=False)
    cnt = pd.concat([pd.read_parquet(OUT / f"{k}_meta.parquet").assign(series=k) for k in ("gse111563", "gse59515", "gse55374")])
    cnt = cnt.groupby(["series", "tp"]).size().rename("n_samples")
    cnt.to_csv(OUT / "sample_counts.csv"); print(cnt.to_string())
