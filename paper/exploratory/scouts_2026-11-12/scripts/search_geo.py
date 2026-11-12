"""GEO search via E-utilities: list human GSE series matching repeat-biopsy design queries (metadata only)."""
import json, sys, time, urllib.parse, urllib.request
B = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/"
def get(u):
    for _ in range(3):
        try:
            return urllib.request.urlopen(u, timeout=40).read().decode()
        except Exception:
            time.sleep(2)
    return "{}"
QUERIES = [
 "skeletal muscle biopsy exercise before after human[orgn] RNA-seq",
 "skeletal muscle acute exercise vastus lateralis pre post biopsy",
 "muscle biopsy resistance exercise time course 1h 3h human",
 "muscle biopsy insulin clamp before after human",
 "human adipose tissue biopsy insulin clamp before after",
 "subcutaneous adipose biopsy fasting fed paired human",
 "cardiac surgery right atrial appendage biopsy before after cardiopulmonary bypass",
 "human heart biopsy ischemia reperfusion before after cardioplegia",
 "human liver biopsy ischemia reperfusion transplant pre post reperfusion",
 "kidney transplant biopsy time zero post reperfusion human",
 "surgery start end biopsy same patient muscle human",
 "tourniquet skeletal muscle biopsy human",
 "bed rest muscle biopsy human",
 "human muscle biopsy placebo control repeated biopsies same day",
 "skeletal muscle biopsy anesthesia human",
 "human skin biopsy repeated same subject hours",
 "human adipose tissue biopsy time course after meal",
 "human muscle biopsy around the clock time of day",
 "endurance exercise muscle biopsy 3 h post human transcriptome",
 "high intensity interval exercise muscle biopsy transcriptome human",
 "muscle biopsy before and after exercise training acute bout",
 "human intestinal mucosa biopsy before after human",
 "human endometrium biopsy time course hours",
 "human adipose biopsy before after exercise",
 "skeletal muscle biopsy glucose ingestion insulin infusion human microarray",
 "muscle biopsy amino acid ingestion time course human",
 "lipid infusion intralipid muscle biopsy human",
 "cold exposure human adipose muscle biopsy",
 "human muscle sleep deprivation biopsy",
 "human adipose sleep restriction biopsy",
]
seen = {}
for q in QUERIES:
    t = urllib.parse.quote(f"({q}) AND Homo sapiens[orgn] AND gse[ETYP]")
    r = json.loads(get(f"{B}esearch.fcgi?db=gds&term={t}&retmax=40&retmode=json") or "{}")
    ids = r.get("esearchresult", {}).get("idlist", [])
    if not ids: continue
    s = json.loads(get(f"{B}esummary.fcgi?db=gds&id={','.join(ids)}&retmode=json"))
    for i in ids:
        d = s["result"].get(i)
        if not d or d["accession"] in seen: continue
        seen[d["accession"]] = d
        print(f"{d['accession']}\tn={d['n_samples']}\t{d['gdstype'][:30]}\t{d['title'][:110]}\t[{q[:30]}]")
        sys.stdout.flush()
    time.sleep(0.4)
json.dump(seen, open(sys.argv[1], "w"))
