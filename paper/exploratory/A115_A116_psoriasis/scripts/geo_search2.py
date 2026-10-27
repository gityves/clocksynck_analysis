"""Second GEO scan: psoriasis series whose title or summary mentions placebo / vehicle / double-blind (titles only).
Output output/geo_search_placebo_hits.csv"""
import json, time, urllib.request, urllib.parse
from pathlib import Path
import pandas as pd
Q = ['psoriasis AND placebo AND gse[Entry Type]', 'psoriasis AND "double-blind" AND gse[Entry Type]', 'psoriasis AND "randomized" AND gse[Entry Type]',
     'psoriasis AND vehicle AND gse[Entry Type]', '"plaque psoriasis" AND (PASI) AND gse[Entry Type]', 'psoriasis AND biologic AND skin AND (week OR day) AND gse[Entry Type]',
     'psoriasis AND (risankizumab OR tildrakizumab OR bimekizumab OR certolizumab OR spesolimab OR "PF-06700841" OR "BI 730357" OR "JNJ-77242113") AND gse[Entry Type]',
     'psoriasis AND (LAS41008 OR "ixekizumab" OR "secukinumab") AND placebo AND gse[Entry Type]']
def get(u):
    for _ in range(3):
        try: return json.load(urllib.request.urlopen(u, timeout=40))
        except Exception: time.sleep(2)
ids = []
for q in Q:
    r = get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=gds&retmax=300&retmode=json&term=" + urllib.parse.quote(q)); ids += r["esearchresult"]["idlist"]; time.sleep(0.4)
ids = sorted(set(ids)); rows = []
for i in range(0, len(ids), 100):
    r = get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=gds&retmode=json&id=" + ",".join(ids[i:i + 100]))
    for k, v in r["result"].items():
        if k != "uids": rows.append(dict(accession=v.get("accession"), n_samples=v.get("n_samples"), taxon=v.get("taxon"), gpl=v.get("gpl"), type=v.get("gdstype"), title=v.get("title"), summary=v.get("summary", "")))
    time.sleep(0.4)
d = pd.DataFrame(rows); d = d[d.accession.str.startswith("GSE", na=False) & (d.taxon == "Homo sapiens")]
t = (d.title + " " + d.summary).str.lower()
d["placebo_flag"] = t.str.contains("placebo|vehicle|double-blind|double blind|randomi")
d.to_csv(Path(__file__).resolve().parents[1] / "output/geo_search_placebo_hits.csv", index=False)
pd.set_option("display.width", 250); pd.set_option("display.max_colwidth", 120); pd.set_option("display.max_rows", 500)
x = d[d.placebo_flag]; print(len(d), len(x)); print(x[["accession", "n_samples", "type", "title"]].to_string())
