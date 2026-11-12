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
 "right atrial appendage cardiac surgery gene expression",
 "myocardium biopsy cardioplegia before after ischemia human",
 "left ventricular biopsy aortic cross-clamp reperfusion human",
 "ischemic preconditioning human heart biopsy",
 "heart transplant donor biopsy ischemic time reperfusion",
 "kidney transplant biopsy pre-implantation post-reperfusion gene expression",
 "liver biopsy time zero post-reperfusion transplantation microarray",
 "liver biopsy beginning end of surgery human",
 "adipose tissue biopsy beginning and end of surgery human",
 "omental subcutaneous adipose paired biopsy surgery duration",
 "muscle biopsy during surgery orthopedic hip knee arthroplasty paired",
 "sevoflurane propofol anesthesia tissue gene expression human biopsy",
 "remote ischemic preconditioning human biopsy",
 "muscle biopsy 4 hours after meal human RNA-seq",
 "muscle biopsy fasted fed state human transcriptome",
 "adipose tissue biopsy fasted postprandial human transcriptome",
 "skeletal muscle biopsies 0 h 3 h 6 h 24 h exercise human",
 "muscle biopsy time course after exercise 1 h 4 h 8 h human",
 "muscle biopsy control leg rest no exercise paired same day",
 "metaMEx exercise skeletal muscle human transcriptome",
 "human skeletal muscle circadian time of day biopsy",
 "human adipose tissue circadian biopsy time of day",
 "skin biopsy time course after UV irradiation human",
 "liver biopsy intraoperative repeated human gene expression",
 "adipose tissue biopsy 2 h insulin infusion euglycemic clamp human RNA-seq",
 "growth hormone or glucocorticoid infusion muscle biopsy human",
 "human muscle biopsy lipopolysaccharide endotoxin infusion",
 "high fat meal muscle adipose biopsy time course human",
 "stomach duodenum biopsy before after meal human",
 "human colon mucosa biopsy repeated same day",
 "heat stress water immersion muscle biopsy human",
 "hypoxia exposure muscle biopsy human",
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
