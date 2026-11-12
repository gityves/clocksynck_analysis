"""Fetch GEO brief text for series+samples (metadata only) and print a compact digest."""
import sys, urllib.request, re, os, time
out = sys.argv[1]; os.makedirs(out, exist_ok=True)
for acc in sys.argv[2:]:
    u = f"https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={acc}&targ=all&form=text&view=brief"
    p = f"{out}/{acc}.txt"
    if not os.path.exists(p):
        for _ in range(3):
            try:
                open(p, "w").write(urllib.request.urlopen(u, timeout=90).read().decode(errors="replace")); break
            except Exception as e:
                time.sleep(3)
        time.sleep(0.4)
    t = open(p).read().splitlines()
    print("=====", acc)
    for l in t:
        if l.startswith(("!Series_title", "!Series_overall_design", "!Series_summary", "!Series_platform_id", "!Series_pubmed_id")):
            print(l[:700])
    sm = [i for i, l in enumerate(t) if l.startswith("^SAMPLE")]
    print("n_samples", len(sm))
    for i in sm[:int(os.environ.get("NS", 6))]:
        blk = []
        for l in t[i:i+60]:
            if l.startswith("^") and blk: break
            if l.startswith(("^SAMPLE", "!Sample_title", "!Sample_characteristics", "!Sample_source_name")): blk.append(l[:160])
        print(" | ".join(blk))
