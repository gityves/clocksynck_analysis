"""Parse series-matrix headers (metadata only) into per-sample tables; print platform, n and characteristic key/value summaries."""
import sys, re, glob, collections
from pathlib import Path
import pandas as pd
H = Path(sys.argv[1]) if len(sys.argv) > 1 else None

def parse(p):
    h = {}
    for line in open(p, errors="replace"):
        if line.startswith("!"):
            k, *v = line.rstrip("\n").split("\t"); h.setdefault(k, []).append([x.strip('"') for x in v])
    return h

def table(h):
    n = len(h["!Sample_geo_accession"][0])
    d = {"gsm": h["!Sample_geo_accession"][0], "title": h["!Sample_title"][0]}
    ch = h.get("!Sample_characteristics_ch1", [])
    for i, row in enumerate(ch):
        row = row + [""] * (n - len(row))
        keys = {r.split(":", 1)[0].strip() for r in row if ":" in r}
        k = collections.Counter(r.split(":", 1)[0].strip() for r in row if ":" in r).most_common(1)[0][0] if keys else f"char{i}"
        col = k
        while col in d: col += "_"
        d[col] = [r.split(":", 1)[1].strip() if ":" in r else r for r in row]
    d["platform"] = h["!Sample_platform_id"][0]
    d["source"] = h.get("!Sample_source_name_ch1", [[""] * n])[0]
    return pd.DataFrame(d)

if __name__ == "__main__":
    for p in sorted(H.glob("*.hdr")):
        h = parse(p); t = table(h)
        print("=" * 100); print(p.stem, "n=", len(t), "platforms", t.platform.value_counts().to_dict())
        print("title:", h.get("!Series_title", [[""]])[0][0][:150])
        print("summary:", " ".join(h.get("!Series_summary", [[""]])[0])[:300])
        print("source:", t.source.value_counts().head(6).to_dict())
        for c in t.columns[2:]:
            if c in ("platform", "source"): continue
            vc = t[c].value_counts()
            print(f"  [{c}] nuniq={len(vc)}", dict(vc.head(7)) if len(vc) <= 40 or True else "")
