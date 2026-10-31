"""Per-patient cell-type counts from GSE206785_metadata.txt.gz (cell labels only; no expression). python -I."""
import gzip, csv, collections as C
f = "/mnt/max8tb/data/CoPhaser_data/cancer_circadian/scout_sc_gastric/GSE206785/GSE206785_metadata.txt.gz"
rd = csv.DictReader(gzip.open(f, "rt"), delimiter="\t")
cnt = C.Counter(); types = C.Counter(); ann = C.Counter()
for r in rd:
    cnt[(r["Patient"], r["Tissue"], r["Type"])] += 1; types[(r["Tissue"], r["Type"])] += 1
    if r["Type"] in ("Epithelial",) or "pithel" in r["Type"] or "umor" in r["Annotation"]: ann[(r["Tissue"], r["Type"], r["Annotation"])] += 1
print(sorted(types.items()))
print(len({k[0] for k in cnt}), "patients")
pats = sorted({k[0] for k in cnt}); tys = sorted({k[2] for k in cnt})
with open("GSE206785_patient_x_type.csv", "w", newline="") as o:
    w = csv.writer(o); w.writerow(["patient","tissue"]+tys)
    for p in pats:
        for t in ("Normal","Tumor"):
            w.writerow([p,t]+[cnt.get((p,t,y),0) for y in tys])
for k,v in sorted(ann.items()): print(k,v)
