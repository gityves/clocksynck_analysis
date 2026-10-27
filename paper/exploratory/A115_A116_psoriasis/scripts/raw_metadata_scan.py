"""Metadata-only scan of every raw GEO series matrix of skin_lesional: per-sample status counts (no expression read)."""
import gzip, re, sys, collections
from pathlib import Path
import pandas as pd
ROOT = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian/skin_lesional")
rows = []
for d in sorted(ROOT.glob("GSE*")):
    f = next(iter(d.glob("*series_matrix.txt.gz")), None)
    if f is None: continue
    fields = collections.defaultdict(list)
    with gzip.open(f, "rt", errors="replace") as h:
        for line in h:
            if line.startswith("!series_matrix_table_begin"): break
            if line.startswith("!Sample_"):
                p = line.rstrip("\n").split("\t"); fields[p[0]].append([x.strip('"') for x in p[1:]])
    n = len(fields["!Sample_title"][0])
    txt = [" | ".join(col[i] for k in ("!Sample_title","!Sample_source_name_ch1","!Sample_characteristics_ch1") for col in fields[k]) for i in range(n)]
    low = [t.lower() for t in txt]
    healthy = [bool(re.search(r"healthy|normal skin|ctrl|control|disease: c\b|nn_sample|\bnn\b|normal \(c\)", t)) and not re.search(r"non-?lesion|\bnl\b|nonlesion", t) for t in low]
    les = [bool(re.search(r"\bls\b|lesional skin|_lesional|lesion: true|disease: lp|psoriasis lesion|lesional", t)) and not re.search(r"non-?lesion|\bnl\b", t) for t in low]
    nl = [bool(re.search(r"non-?lesion|\bnl\b|nonlesion|disease: nlp|pn_sample|uninvolved", t)) for t in low]
    title = fields["!Sample_title"][0]
    ser_title = ""
    with gzip.open(f, "rt", errors="replace") as h:
        for line in h:
            if line.startswith("!Series_title"): ser_title = line.split("\t",1)[1].strip().strip('"'); break
    rows.append(dict(series=d.name, n_samples=n, n_flag_lesional=sum(les), n_flag_nonlesional=sum(nl), n_flag_healthy=sum(healthy),
                     example_healthy=" ; ".join(t[:70] for t,h in zip(title,healthy) if h)[:140], series_title=ser_title[:110]))
df = pd.DataFrame(rows); df.to_csv("output/raw_metadata_scan.csv", index=False)
pd.set_option("display.width",250); pd.set_option("display.max_colwidth",60)
print(df.drop(columns=["example_healthy"]).to_string())
print(df[df.n_flag_healthy>0][["series","n_flag_healthy","example_healthy"]].to_string())
