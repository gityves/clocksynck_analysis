"""Counts for TALLY.md from tally.csv (pandas-free)."""
import csv, collections as C
rows = list(csv.DictReader(open("tally.csv")))
N = len(rows)
def pct(a,b): return f"{a}/{b} ({100*a/b:.0f}%)" if b else "0/0"
def tab(title, keyf, rs=rows):
    print("\n##", title)
    g = C.defaultdict(list)
    for r in rs: g[keyf(r)].append(r)
    print("group | n | pos | mixed | null | unresolved | headline-positive of resolved")
    for k,v in sorted(g.items(), key=lambda x:-len(x[1])):
        c=C.Counter(r["claim"] for r in v); res=len(v)-c["unresolved"]
        print(f"{k} | {len(v)} | {c['pos']} | {c['mixed']} | {c['null']} | {c['unresolved']} | {pct(c['pos']+c['mixed'],res)}")
tab("by cancer group", lambda r:r["group"])
tab("breast vs other", lambda r:"breast" if r["group"]=="Breast" else "non-breast")
tab("claim basis", lambda r:r["read_level"])
for fld in ["validation","adj_prolif","adj_subtype","random_ctrl","self_null"]:
    c=C.Counter(r[fld].split(" ")[0] for r in rows); print("\n##",fld,dict(c))
# funnel
print("\n## funnel (n known)")
bins=[(0,100,"<100"),(100,250,"100-249"),(250,600,"250-599"),(600,10**9,">=600")]
for lo,hi,l in bins:
    v=[r for r in rows if r["n_num"] and lo<=float(r["n_num"])<hi]
    c=C.Counter(r["claim"] for r in v); res=len(v)-c["unresolved"]
    print(l,len(v),dict(c),pct(c["pos"]+c["mixed"],res), "median n"," ")
print("n unknown:",sum(1 for r in rows if not r["n_num"]))
# restricted to claim_basis not inferred
for lab,f in [("read (F/A/meta_text)",lambda r:r["claim_basis"] in("F","A","meta_text")),("inferred only",lambda r:r["claim_basis"]=="I")]:
    v=[r for r in rows if f(r)]; c=C.Counter(r["claim"] for r in v); print("\n",lab,len(v),dict(c))
# type of design
tab("design", lambda r:"signature/LASSO/score" if any(k in r["readout"].lower() for k in["signature","score","model","index","pca"]) else "single gene/protein", rows)
tab("single-gene series by n<250", lambda r:"n<250" if r["n_num"] and float(r["n_num"])<250 else "n>=250/unk", [r for r in rows if r["read_level"]=="M"])
print("\n## breast rows: adjustment")
b=[r for r in rows if r["group"]=="Breast"]
for f in ["adj_prolif","adj_subtype","validation"]: print(f,dict(C.Counter(r[f].split(" ")[0] for r in b)))
print("breast rows with read basis F:",[r["first_author"] for r in b if r["read_level"]=="F"])
print("distinct years", min(int(r["year"]) for r in rows), max(int(r["year"]) for r in rows))
y=C.Counter((int(r["year"])//5*5) for r in rows); print(sorted(y.items()))
