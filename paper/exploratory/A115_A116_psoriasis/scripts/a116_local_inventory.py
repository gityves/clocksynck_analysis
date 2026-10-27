"""A116 scouting of what is already on disk: patients per arm and visit with lesional and uninvolved biopsies (counts only; no
offset, no outcome-linked statistic). Reads metadata columns of the processed parquets.

    python -I a116_local_inventory.py
Output: output/a116_local_inventory.csv
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import pandas as pd, pyarrow.parquet as pq
import _lib as L
rows = []
for folder in ("processed", "processed_heldout", "processed_treatment"):
    for f in sorted((L.ROOT / folder).glob("GSE*.parquet")):
        names = pq.ParquetFile(f).schema.names
        cols = [c for c in ("patient", "status", "disease", "visit", "weeks", "treatment", "arm", "pasi") if c in names]
        d = pd.read_parquet(f, columns=cols)
        d = d[d.disease == "psoriasis"] if "disease" in d else d
        if d.empty: continue
        has = {c: c in names for c in ("visit", "treatment", "arm", "pasi")}
        g = [c for c in ("visit", "arm") if c in d]
        if not g: d = d.assign(none="-"); g = ["none"]
        t = d.groupby(g + ["status"]).patient.nunique().unstack("status").fillna(0).astype(int).reset_index()
        t.insert(0, "file", f"{folder}/{f.name}")
        for c, v in has.items(): t[f"has_{c}"] = v
        if "pasi" in d and g: t["pasi_nonnull"] = d.groupby(g).pasi.apply(lambda x: x.notna().sum()).values
        rows.append(t)
r = pd.concat(rows, ignore_index=True); r.to_csv(L.PAPER / "exploratory/A115_A116_psoriasis/output/a116_local_inventory.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_rows", 500); print(r.to_string())
