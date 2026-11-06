"""A226 scout step 1a: screen the user's GEO database table (sample titles) for studies with tumour and adjacent/normal samples.
Metadata only. Writes scout/db_screen.csv."""
import re, sys
from pathlib import Path
import pandas as pd
ROOT = Path(__file__).resolve().parents[4]
DB = Path("/home/maxine/Documents/paychere/scRNAseq_database/output/datasets")
S = pd.read_csv(DB / "build_private_human_atlas_samples.csv")
St = pd.read_csv(DB / "build_private_human_atlas_studies.csv").set_index("accession")
used = set(pd.read_csv(ROOT / "paper/exploratory/A151_fibroblast_pairing/output/units.csv").study) | \
       set(pd.read_csv(ROOT / "paper/exploratory/A149_paired_tumour_normal/output/units.csv").study)
T = re.compile(r"tumou?r|cancer|carcinoma|malignan|\bPT\b|primary|\bTT\b|\bLC\b|\bCA\b|[-_ ]T\d*\b|[-_]T$|HCC|ICC|PDAC|NSCLC|LUAD|LUSC", re.I)
N = re.compile(r"normal|adjacent|para-?(tumou?r|cancer)|non-?(tumou?r|malignant|cancer)|healthy|benign|[-_ ]N\d*\b|[-_]N$|\bNAT\b|\bANT\b|\bNT\b|control|margin|peritumou?r", re.I)
rows = []
for st, g in S.groupby("study"):
    t = g.title.fillna("") + " " + g.description.fillna("")
    nt = t.map(lambda x: bool(T.search(x))).sum(); nn = t.map(lambda x: bool(N.search(x))).sum()
    both = t.map(lambda x: bool(T.search(x)) and bool(N.search(x))).sum()
    if nt >= 3 and nn >= 3:
        rows.append(dict(study=st, n_gsm=len(g), n_tumour_like=nt, n_normal_like=nn, ambiguous_both=both,
                         used_A149_A151=st in used, title=St.title.get(st, "")[:110], date=St.submission_date.get(st, ""),
                         example=" | ".join(g.title.head(4).astype(str))[:150]))
D = pd.DataFrame(rows).sort_values(["used_A149_A151", "n_gsm"], ascending=[True, False])
D.to_csv(ROOT / "paper/exploratory/A226_stroma_pairing_confirm/scout/db_screen.csv", index=False)
print(len(S.study.unique()), "studies;", len(D), "tumour+normal-like;", (~D.used_A149_A151).sum(), "not used")
