"""A226 scout 1c: GEO (NCBI E-utilities, metadata only) series 2022-2026, human, single-cell/nucleus, tumour with adjacent/paired
normal. Writes scout/geo_series.csv (title, n samples, sample-title tumour/normal counts, supplementary types)."""
import json, re, sys, time, urllib.parse, urllib.request
from pathlib import Path
import pandas as pd
ROOT = Path(__file__).resolve().parents[4]
OUT = Path(__file__).resolve().parents[1] / "scout"
E = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
def get(u, **p):
    for i in range(5):
        try:
            return urllib.request.urlopen(E + u + "?" + urllib.parse.urlencode(p), timeout=90).read()
        except Exception:
            time.sleep(3)
    raise RuntimeError(u)
term = ('(single cell OR single-cell OR single nucleus OR single-nucleus OR snRNA-seq OR scRNA-seq) AND (tumor OR tumour OR cancer OR carcinoma) '
        'AND (adjacent normal OR paired normal OR matched normal OR para-tumor OR paracancerous OR adjacent tissue OR normal tissue OR non-tumor OR '
        'nontumor OR peritumor) AND "Homo sapiens"[Organism] AND gse[Entry Type] AND 2018/01/01:2026/12/31[PDAT]')
ids = []
for st in range(0, 20000, 500):
    r = json.loads(get("esearch.fcgi", db="gds", retmax=500, retstart=st, retmode="json", term=term))["esearchresult"]
    ids += r["idlist"]
    if st + 500 >= int(r["count"]): break
print(len(ids), "records", flush=True)
used = set(pd.read_csv(ROOT / "paper/exploratory/A151_fibroblast_pairing/output/units.csv").study) | set(pd.read_csv(ROOT / "paper/exploratory/A149_paired_tumour_normal/output/units.csv").study)
T = re.compile(r"tumou?r|cancer|carcinoma|malignan|\bPT\b|HCC|NSCLC|LUAD|LUSC|[-_ .]T\d*$|[-_ .]T[-_ .\d]|CRC|PDAC|ccRCC", re.I)
N = re.compile(r"normal|adjacent|para-?(tumou?r|cancer|carcinoma)|non-?(tumou?r|malignant|cancer)|benign|[-_ .]N\d*$|[-_ .]N[-_ .\d]|\bNAT\b|\bANT\b|\bNT\b|margin|peritumou?r|\bnon.?cancer", re.I)
rows = []
for i in range(0, len(ids), 150):
    d = json.loads(get("esummary.fcgi", db="gds", id=",".join(ids[i:i+150]), retmode="json"))["result"]
    for k in d.get("uids", []):
        x = d[k]
        if not x.get("accession", "").startswith("GSE"): continue
        smp = [s.get("title", "") for s in x.get("samples", [])]
        nt = sum(bool(T.search(s)) and not N.search(s) for s in smp); nn = sum(bool(N.search(s)) for s in smp)
        rows.append(dict(study=x["accession"], date=x.get("pdat"), n_samples=len(smp), n_t=nt, n_n=nn, used=x["accession"] in used,
                         suppl=x.get("suppfile", ""), gpl=x.get("gpl", ""), title=x.get("title", "")[:150], etype=x.get("gdstype", ""),
                         summary=x.get("summary", "")[:600], sample_titles=" | ".join(smp[:60])))
    print(i, flush=True)
D = pd.DataFrame(rows).drop_duplicates("study")
D.to_csv(OUT / "geo_series.csv", index=False)
print(len(D), (D.n_t.ge(3) & D.n_n.ge(3)).sum())
