"""A115 inventory table: per psoriasis / atopic dermatitis series (processed matrices) the lesional, uninvolved (pairs) and
healthy-control counts, platform, and the presence of the 14 GTEx clock genes, the A81b cue set (40), locked set (40) and the A78
systemic set (40). Healthy counts come from the raw GEO sample metadata (output/raw_metadata_scan.csv); healthy samples are not in
the processed parquets. No expression statistic is computed.   Output: output/a115_inventory.csv"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import pandas as pd, pyarrow.parquet as pq
import _lib as L
OUT = Path(__file__).resolve().parents[1] / "output"
sets = L.cue_sets(); raw = pd.read_csv(OUT / "raw_metadata_scan.csv").set_index("series")
H = {"GSE13355": (64, "series_matrix (GPL570), same table"), "GSE121212": (38, "GSE121212_readcount.txt.gz, CTRL_*_healthy columns"),
     "GSE201827": (17, "series_matrix (GPL570), 'Untreated.Control.Normal.*'"), "GSE268958": (19, "GSE268958_raw_counts.tsv.gz, 'condition: Healthy' (NRML)"),
     "GSE66511": (12, "GSE66511_Psoriasis_counts.txt.gz, disease: C")}
rows = []
for folder in ("processed", "processed_heldout"):
    for f in sorted((L.ROOT / folder).glob("GSE*.parquet")):
        names = set(pq.ParquetFile(f).schema.names); names |= {"ARNTL"} if "BMAL1" in names else set()
        d = pd.read_parquet(f, columns=["patient", "status", "disease", "platform", "value_type"])
        for dis, g in d.groupby("disease"):
            if dis not in ("psoriasis", "atopic_dermatitis"): continue
            pat = set(g[g.status == "lesional"].patient) & set(g[g.status == "nonlesional"].patient)
            s = f.stem; h = H.get(s, (0, "")) if True else None
            rows.append(dict(series=s, disease=dis, folder=folder, platform=g.platform.iloc[0], value_type=g.value_type.iloc[0],
                             n_lesional=int((g.status == "lesional").sum()), n_uninvolved=int((g.status == "nonlesional").sum()), n_pairs=len(pat),
                             n_healthy_raw=h[0], healthy_source=h[1], healthy_in_processed=0,
                             core14=sum(x in names for x in L.CORE), cue40=sum(x in names for x in sets["cue"].genes),
                             locked40=sum(x in names for x in sets["locked"].genes), systemic40=sum(x in names for x in sets["systemic"].genes)))
r = pd.DataFrame(rows); r.to_csv(OUT / "a115_inventory.csv", index=False)
pd.set_option("display.width", 250); print(r.drop(columns=["healthy_source", "folder"]).to_string())
