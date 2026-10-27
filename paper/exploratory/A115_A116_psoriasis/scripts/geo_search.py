"""List GEO series matching psoriasis + placebo/treatment queries (titles and sample counts only). Output output/geo_search_hits.csv"""
import json, time, urllib.request, urllib.parse
from pathlib import Path
import pandas as pd
Q = ['psoriasis AND placebo AND (biopsy OR skin) AND "Homo sapiens"[Organism] AND gse[Entry Type]',
     'psoriasis AND (randomized OR clinical trial) AND (lesional OR non-lesional) AND gse[Entry Type]',
     'psoriasis AND (secukinumab OR ixekizumab OR brodalumab OR guselkumab OR risankizumab OR tildrakizumab OR bimekizumab OR ustekinumab) AND gse[Entry Type]',
     'psoriasis AND (apremilast OR tofacitinib OR baricitinib OR deucravacitinib OR methotrexate OR etanercept OR adalimumab OR infliximab OR efalizumab) AND gse[Entry Type]',
     'psoriasis AND (time course OR "week 12" OR "week 16" OR "week 4") AND gse[Entry Type]']
def get(u):
    for _ in range(3):
        try: return json.load(urllib.request.urlopen(u, timeout=40))
        except Exception as e: time.sleep(2)
ids = []
for q in Q:
    r = get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&retmax=200&retmode=json&term=" + urllib.parse.quote(q))
    ids += r["esearchresult"]["idlist"]; time.sleep(0.4)
ids = sorted(set(ids))
rows = []
for i in range(0, len(ids), 100):
    r = get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gds&retmode=json&id=" + ",".join(ids[i:i + 100]))
    for k, v in r["result"].items():
        if k == "uids": continue
        rows.append(dict(accession=v.get("accession"), gdstype=v.get("gdstype"), n_samples=v.get("n_samples"), taxon=v.get("taxon"), gpl=v.get("gpl"), title=v.get("title")[:140], summary=v.get("summary", "")[:300]))
    time.sleep(0.4)
d = pd.DataFrame(rows); d = d[d.accession.str.startswith("GSE", na=False) & (d.taxon == "Homo sapiens")]
d.to_csv(Path(__file__).resolve().parents[1] / "output/geo_search_hits.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 110); pd.set_option("display.max_rows", 500)
print(len(d)); print(d[["accession", "n_samples", "gdstype", "title"]].to_string())
