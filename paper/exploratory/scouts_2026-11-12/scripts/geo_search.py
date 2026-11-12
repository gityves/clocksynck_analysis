"""List GEO series (accession, n samples, taxon, title) for an Entrez gds query. Usage: geo_search.py "<query>" [retmax]"""
import sys, json, urllib.request, urllib.parse
S = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
import time
def get(u):
    time.sleep(1.2)
    return json.load(urllib.request.urlopen(u, timeout=60))
q = sys.argv[1]; n = sys.argv[2] if len(sys.argv) > 2 else "100"
ids = get(S + "esearch.fcgi?db=gds&retmode=json&retmax=%s&term=%s" % (n, urllib.parse.quote(q)))["esearchresult"]["idlist"]
for i in range(0, len(ids), 50):
    d = get(S + "esummary.fcgi?db=gds&retmode=json&id=" + ",".join(ids[i:i+50]))["result"]
    for k in d["uids"]:
        r = d[k]; print(r["accession"], r.get("n_samples"), r.get("taxon"), "|", r["title"][:110])
