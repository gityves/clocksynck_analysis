"""Convert the GSE141549 normalised expression xlsx (no values in the GEO series matrices) to parquet, parsing only.
python -I s00_convert_xlsx.py   (reads the first sheet; columns: ID_REF, Transcript, ILMN_Gene, PROBE_SEQUENCE, SAMPLE ...)"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
import openpyxl, pandas as pd
import a161_lib as L
(L.EN / "derived").mkdir(exist_ok=True)
for tag, f in (("HT12", "GSE141549_Illumina_HumanHT-12_V4.0_expression_beadchip_normalized.xlsx"),
               ("WG6", "GSE141549_Illumina_HumanWG-6_v2.0_expression_beadchip_normalized.xlsx")):
    wb = openpyxl.load_workbook(L.EN / f, read_only=True)
    ws = wb.worksheets[0]
    it = ws.iter_rows(values_only=True)
    hdr = next(it)
    rows = [r for r in it if r[0] is not None]
    df = pd.DataFrame(rows, columns=[str(h) for h in hdr])
    ann = df[["ID_REF", "ILMN_Gene"]]
    num = df[[c for c in df.columns if c.startswith("SAMPLE")]].apply(pd.to_numeric, errors="coerce")
    num.insert(0, "ILMN_Gene", df["ILMN_Gene"].values); num.index = df["ID_REF"].values
    num.to_parquet(L.EN / "derived" / f"GSE141549_{tag}.parquet")
    print(tag, num.shape, flush=True)
