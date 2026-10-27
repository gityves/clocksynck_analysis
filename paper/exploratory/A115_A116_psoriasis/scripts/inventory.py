"""A115 inventory: sample counts per status for every processed skin series (metadata only; no clock statistic)."""
import sys
from pathlib import Path
import pandas as pd, pyarrow.parquet as pq
ROOT = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian/skin_lesional")
META = ["patient","status","sample","dataset","disease","platform","value_type","depth","visit","treatment","arm","pasi"]
rows = []
for folder in ["processed","processed_heldout","processed_other","processed_other_linear","processed_treatment"]:
    for f in sorted((ROOT/folder).glob("GSE*.parquet")):
        names = pq.ParquetFile(f).schema.names
        cols = [c for c in META if c in names]
        d = pd.read_parquet(f, columns=cols)
        up = {n.upper() for n in names}
        for dis, g in d.groupby("disease"):
            vc = g.status.value_counts().to_dict()
            rows.append(dict(folder=folder, file=f.name, series=f.stem.split("_")[0], disease=dis,
                platform=g.platform.iloc[0] if "platform" in g else "", value_type=g.value_type.iloc[0],
                n_samples=len(g), n_patients=g.patient.nunique(), statuses=str(vc),
                n_lesional=vc.get("lesional",0), n_nonlesional=vc.get("nonlesional",0),
                n_healthy=sum(v for k,v in vc.items() if k not in ("lesional","nonlesional")),
                n_genes=len(names), has_clock14=all(x in up or (x=="ARNTL" and "BMAL1" in up) for x in
                  "ARNTL CLOCK CRY1 CRY2 DBP NPAS2 NR1D1 NR1D2 PER1 PER2 PER3 TEF CIART HLF".split()) ,
                columns=",".join(cols)))
pd.DataFrame(rows).to_csv("output/inventory_processed.csv", index=False)
pd.set_option("display.width",250); pd.set_option("display.max_colwidth",70)
print(pd.DataFrame(rows).drop(columns=["columns","n_genes"]).to_string())
