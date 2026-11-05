"""Map the Zotero export to references.bib (DOI, then title), extract full texts, write index.tsv.

Usage: python -I lit_index.py <export dir> <references.bib> <out index.tsv>
"""
import html
import re
import subprocess
import sys
from html.parser import HTMLParser
from pathlib import Path


def entries(path):
    txt = Path(path).read_text(encoding="utf-8", errors="replace")
    out = []
    for m in re.finditer(r"@(\w+)\{([^,]+),(.*?)\n\}", txt, re.S):
        body = m.group(3)
        f = {k.lower(): v for k, v in re.findall(r"^\s*(\w+)\s*=\s*\{(.*?)\},?\s*$", body, re.M)}
        out.append(dict(type=m.group(1), key=m.group(2).strip(), **f))
    return out


def norm_title(t):
    return re.sub(r"[^a-z0-9]", "", html.unescape(re.sub(r"[{}\\]", "", t or "")).lower())[:60]


def norm_doi(d):
    return (d or "").lower().strip().removeprefix("https://doi.org/")


class Strip(HTMLParser):
    def __init__(self):
        super().__init__()
        self.parts, self.skip = [], 0

    def handle_starttag(self, tag, a):
        self.skip += tag in ("script", "style")

    def handle_endtag(self, tag):
        self.skip -= tag in ("script", "style") and self.skip > 0
        if tag in ("p", "div", "h1", "h2", "h3", "li", "tr", "br"):
            self.parts.append("\n")

    def handle_data(self, d):
        if not self.skip:
            self.parts.append(d)


def extract(src, dst):
    if src.suffix.lower() == ".pdf":
        subprocess.run(["pdftotext", str(src), str(dst)], check=False, capture_output=True)
    else:
        p = Strip()
        p.feed(src.read_text(encoding="utf-8", errors="replace"))
        dst.write_text(re.sub(r"\n\s*\n+", "\n\n", "".join(p.parts)))
    return len(dst.read_text(errors="replace")) if dst.exists() else 0


MANUAL = {"hammarlund2025cyclops": "_cyclops_"}  # thesis record saved as an HTML snapshot


def main(exp, ours, out):
    exp = Path(exp)
    zot = entries(exp / "clocksync.bib")
    ref = entries(ours)
    by_doi = {norm_doi(e.get("doi")): e for e in zot if e.get("doi")}
    by_title = {norm_title(e.get("title")): e for e in zot}
    by_key = {e["key"]: e for e in zot}
    (exp / "text").mkdir(exist_ok=True)
    used, rows = set(), []
    for r in ref:
        z = by_doi.get(norm_doi(r.get("doi"))) or by_title.get(norm_title(r.get("title"))) or by_key.get(MANUAL.get(r["key"]))
        if z is None:
            rows.append((r["key"], r.get("doi", ""), "", "", 0, "missing"))
            continue
        used.add(z["key"])
        files = [f.replace("\\;", ";") for f in re.split(r"(?<!\\);", z.get("file", "")) if f.strip()]
        best = sorted(files, key=lambda f: not f.lower().endswith(".pdf"))
        if not best:
            rows.append((r["key"], r.get("doi", ""), z["key"], "", 0, "no file"))
            continue
        src = exp / best[0].strip()
        n = extract(src, exp / "text" / f"{r['key']}.txt")
        kind = "pdf" if src.suffix.lower() == ".pdf" else "html"
        status = "ok" if n > 15000 else ("short" if n > 0 else "empty")
        rows.append((r["key"], r.get("doi", ""), z["key"], best[0].strip(), n, f"{status} ({kind})"))
    with open(out, "w") as fh:
        fh.write("citekey\tdoi\tzotero_key\tfile\ttext_chars\tstatus\n")
        for row in rows:
            fh.write("\t".join(map(str, row)) + "\n")
        for z in zot:
            if z["key"] not in used:
                fh.write(f"-\t{z.get('doi', '')}\t{z['key']}\t{z.get('file', '')}\t0\tzotero only: {z.get('title', '')[:80]}\n")


if __name__ == "__main__":
    main(*sys.argv[1:4])
