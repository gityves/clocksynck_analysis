"""Which of the clock genes (GENES), Hallmark G2M genes and purity genes each platform annotates (annotation tables only, no expression)."""
import sys, gzip, io, re
from pathlib import Path
import pandas as pd
ROOT = Path("/home/maxine/Documents/paychere/scCircadianOverview")
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "paper")); sys.path.insert(0, str(ROOT / "paper/exploratory/A127_scanb/scripts"))
import a127_lib as A
from clocksync.genes import canonical
GENES, g2m, PUR = A.GENES, set(A.sh.hallmark_g2m()), list(A.sh.PURITY)
P = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian/scout_luad_lihc/platforms")
SYM = {"GPL10687": "GeneSymbol", "GPL13667": "Gene Symbol", "GPL15048": "GeneSymbol", "GPL5474": "Symbol", "GPL6244": "gene_assignment",
       "GPL7015": "Gene symbol", "GPL96": "Gene Symbol"}
def canon(s):
    c = canonical(s); return c if c in GENES else s
def load(g):
    if (P / f"{g}.annot.gz").exists():
        lines = [l for l in gzip.open(P / f"{g}.annot.gz", "rt", errors="replace") if not l.startswith(("^", "!", "#"))]
        col = "Gene symbol"
    else:
        txt = open(P / f"{g}.txt", errors="replace").read()
        lines = txt.split("!platform_table_begin\n")[1].split("!platform_table_end")[0].splitlines(True); col = SYM[g]
    d = pd.read_csv(io.StringIO("".join(lines)), sep="\t", usecols=["ID", col], dtype=str, quoting=3).dropna()
    if g == "GPL6244":
        d["s"] = d[col].map(lambda x: [e.split("//")[1].strip() for e in x.split("///") if e.count("//") >= 2])
    else:
        d["s"] = d[col].str.split(r"\s*///\s*|\s*//\s*|;|,", regex=True)
    s = set(canon(x.strip()) for l in d.s for x in l if x.strip() and x.strip() != "---")
    return s, d.ID.nunique()
rows = []
for g in ["GPL96", "GPL570", "GPL6884", "GPL7015", "GPL15048", "GPL3921", "GPL571", "GPL5474", "GPL10558", "GPL6480", "GPL10687", "GPL13667", "GPL6244"]:
    s, n = load(g)
    rows.append(dict(platform=g, n_probes=n, n_symbols=len(s), n_clock=sum(x in s for x in GENES), n_clock_total=len(GENES),
                     clock_missing=";".join(x for x in GENES if x not in s), n_g2m=sum(x in s for x in g2m), n_g2m_total=len(g2m),
                     n_purity=sum(x in s for x in PUR), n_purity_total=len(PUR), purity_missing=";".join(x for x in PUR if x not in s)))
df = pd.DataFrame(rows); df.to_csv(Path(__file__).resolve().parents[1] / "platform_coverage.csv", index=False)
print(df.drop(columns=["purity_missing"]).to_string())
