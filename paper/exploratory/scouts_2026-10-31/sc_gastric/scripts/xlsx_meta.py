"""Dump an xlsx (metadata table) to csv with the stdlib only. python -I script.py in.xlsx out.csv"""
import sys, zipfile, re, csv
import xml.etree.ElementTree as ET
ns = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
z = zipfile.ZipFile(sys.argv[1])
ss = []
if "xl/sharedStrings.xml" in z.namelist():
    for si in ET.fromstring(z.read("xl/sharedStrings.xml")).findall("m:si", ns):
        ss.append("".join(t.text or "" for t in si.iter("{%s}t" % ns["m"])))
names = sorted(n for n in z.namelist() if n.startswith("xl/worksheets/sheet"))
def col(ref): 
    s = re.match(r"[A-Z]+", ref).group(); n = 0
    for c in s: n = n*26 + ord(c)-64
    return n-1
w = csv.writer(open(sys.argv[2], "w", newline=""))
for nm in names[:1]:
    for row in ET.fromstring(z.read(nm)).iter("{%s}row" % ns["m"]):
        out = []
        for c in row.findall("m:c", ns):
            i = col(c.get("r")); out += [""]*(i-len(out))
            v = c.find("m:v", ns); t = c.get("t")
            out.append(ss[int(v.text)] if (t == "s" and v is not None) else (v.text if v is not None else "".join(x.text or "" for x in c.iter("{%s}t" % ns["m"]))))
        w.writerow(out)
print(names)
