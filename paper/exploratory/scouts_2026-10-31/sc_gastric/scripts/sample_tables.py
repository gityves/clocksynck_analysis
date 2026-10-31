"""Parse GEO family soft files (metadata only) into one per-sample table. Run with python -I."""
import gzip, csv, sys, re
root = "/mnt/max8tb/data/CoPhaser_data/cancer_circadian/scout_sc_gastric"
gses = ["GSE183904","GSE134520","GSE163558","GSE167297","GSE150290","GSE206785","GSE246662","GSE270680","GSE268238","GSE249874","GSE201347","GSE231539","GSE254188","GSE281274","GSE234209","GSE158631","GSE112302","GSE140182"]
rows = []
for g in gses:
    cur = None
    for ln in gzip.open(f"{root}/{g}/{g}_family.soft.gz", "rt", errors="replace"):
        ln = ln.rstrip("\n")
        if ln.startswith("^SAMPLE"):
            cur = {"gse": g, "gsm": ln.split("= ")[1], "chars": [], "supp": []}; rows.append(cur)
        elif cur is not None and ln.startswith("!Sample_title"): cur["title"] = ln.split("= ",1)[1]
        elif cur is not None and ln.startswith("!Sample_characteristics_ch1"): cur["chars"].append(ln.split("= ",1)[1])
        elif cur is not None and ln.startswith("!Sample_instrument_model"): cur["inst"] = ln.split("= ",1)[1]
        elif cur is not None and ln.startswith("!Sample_supplementary_file"): cur["supp"].append(ln.split("= ",1)[1].split("/")[-1])
        elif ln.startswith("^PLATFORM"): cur = None
with open("samples_all.csv", "w", newline="") as f:
    w = csv.writer(f); w.writerow(["gse","gsm","title","characteristics","instrument","supp"])
    for r in rows: w.writerow([r["gse"], r["gsm"], r.get("title",""), "; ".join(r["chars"]), r.get("inst",""), "|".join(r["supp"])])
print(len(rows))
