"""List GEO FTP matrix/suppl directories (name, size) for series accessions."""
import sys, re, urllib.request
for a in sys.argv[1:]:
    g = a[:-3] + "nnn"
    print("==", a)
    for d in ("matrix", "suppl"):
        try:
            h = urllib.request.urlopen(f"https://ftp.ncbi.nlm.nih.gov/geo/series/{g}/{a}/{d}/", timeout=40).read().decode()
        except Exception as e:
            print("  ", d, "ERR", e); continue
        for m in re.finditer(r'<a href="([^"?/][^"]*)">[^<]*</a>\s*([\d-]+ [\d:]+)\s+(\S+)', h):
            print("  ", d, m.group(1), m.group(3))
