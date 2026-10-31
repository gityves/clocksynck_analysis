"""Tabulate tumour sample counts and matrix sizes from soft/filelists (metadata only). python -I."""
import csv, re, collections as C
R="/mnt/max8tb/data/CoPhaser_data/cancer_circadian/scout_sc_gastric"
rows=list(csv.DictReader(open("samples_all.csv")))
# GSE270680
t=C.Counter(); donors=set(); st=C.Counter(); sub=C.Counter()
for r in rows:
    if r["gse"]!="GSE270680": continue
    m=re.search(r"tissue: ([^;]+)",r["characteristics"]).group(1); t[m]+=1
    if m=="Gastric cancer tumor":
        donors.add(r["title"].split(",")[0]); 
        st[re.search(r"Stage: ([^;]+)",r["characteristics"]).group(1)]+=1; sub[re.search(r"subtype: (.+)$",r["characteristics"]).group(1)]+=1
print("GSE270680",t,len(donors),"tumour donors",st,sub)
size=C.Counter()
for ln in open(f"{R}/GSE270680/filelist.txt"):
    p=ln.rstrip().split("\t")
    if p[0]=="File" and "matrix" in p[1]:
        tag=re.search(r"_D\d+([A-Z])_",p[1]).group(1); size[tag]+=int(p[3])
print({k:round(v/1e9,2) for k,v in size.items()})
for g,key in (("GSE183904","(Tumor)"),("GSE150290","cancer stomach"),("GSE268238","gastric cancer"),("GSE201347","gastric cancer"),("GSE231539","gastric cancer"),("GSE254188","gastric cancer"),("GSE249874","adenocarcinoma")):
    s=[r for r in rows if r["gse"]==g]; n=sum(key in (r["title"]+r["characteristics"]) for r in s); print(g,len(s),"tumour-like",n)
c=C.Counter(re.search(r"lauren classification: (\w+)",r["characteristics"]).group(1) for r in rows if r["gse"]=="GSE150290" and "lauren" in r["characteristics"]); print(c)
c=C.Counter(re.search(r"diagnosis: (\w+)",r["characteristics"]).group(1) for r in rows if r["gse"]=="GSE150290" and "diagnosis" in r["characteristics"]); print(c)
