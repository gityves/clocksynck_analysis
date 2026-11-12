"""NCBI E-utilities GDS search: print accession, n samples, title for each query (metadata only)."""
import sys, json, urllib.request, urllib.parse, time
B = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
def get(u):
    for i in range(4):
        try: return json.loads(urllib.request.urlopen(u, timeout=60).read())
        except Exception as e: time.sleep(2)
    return {}
for q in sys.argv[1:]:
    term = f'({q}) AND gse[Entry Type] AND "Homo sapiens"[Organism]'
    r = get(B + "esearch.fcgi?db=gds&retmax=8&retmode=json&term=" + urllib.parse.quote(term))
    ids = r.get("esearchresult", {}).get("idlist", [])
    print("##", q, len(ids))
    if not ids: continue
    s = get(B + "esummary.fcgi?db=gds&retmode=json&id=" + ",".join(ids))
    for i in ids:
        d = s["result"].get(i, {})
        print(f'{d.get("accession")}\tn={d.get("n_samples")}\t{d.get("gdstype","")[:30]}\t{d.get("title","")[:95]}')
    time.sleep(0.4)
