"""Fetch GEO series-matrix headers only (metadata) and print a compact sample-sheet summary. Parse only."""
import sys, gzip, io, urllib.request, re, json, os
SCR = "/tmp/claude-1000/-home-maxine-Documents-paychere-scCircadianOverview/a7711f70-ea70-4615-a0ec-0b08050a8353/scratchpad/sm"
def urls(gse):
    n = gse[3:]; stub = "GSE" + (n[:-3] if len(n) > 3 else "") + "nnn"
    base = f"https://ftp.ncbi.nlm.nih.gov/geo/series/{stub}/{gse}/matrix/"
    return base
def list_matrix(gse):
    html = urllib.request.urlopen(urls(gse), timeout=60).read().decode()
    return re.findall(r'href="([^"]+series_matrix\.txt\.gz)"', html)
def header(gse, fname):
    p = f"{SCR}/{fname}.hdr"
    if os.path.exists(p): return open(p).read()
    r = urllib.request.urlopen(urls(gse) + fname, timeout=120)
    g = gzip.GzipFile(fileobj=r); out = []
    for line in io.TextIOWrapper(g, errors="replace"):
        if line.startswith("!series_matrix_table_begin"): break
        out.append(line)
    open(p, "w").write("".join(out)); return "".join(out)
def parse(h):
    d = {}
    for l in h.splitlines():
        k, _, v = l.partition("\t")
        d.setdefault(k, []).append(v.split("\t") if v else [])
    return d
if __name__ == "__main__":
    for gse in sys.argv[1:]:
        try:
            fs = list_matrix(gse)
        except Exception as e:
            print(gse, "ERR", e); continue
        for f in fs:
            h = parse(header(gse, f))
            t = h.get("!Series_title", [[""]])[0]
            ids = h.get("!Sample_geo_accession", [[]])[0]
            plat = h.get("!Sample_platform_id", [[]])[0]
            print("="*80); print(f, "|", t, "| n=", len(ids), "| plat", set(plat))
            tit = h.get("!Sample_title", [[]])[0]
            print(" titles:", [x.strip('"') for x in tit[:6]])
            ch = h.get("!Sample_characteristics_ch1", [])
            for row in ch:
                vals = [x.strip('"') for x in row]
                u = sorted(set(vals))
                print("  CH:", (vals[0][:60]), "| nuniq", len(u), "|", u[:6] if len(u) <= 8 else u[:3])
