"""Download processed GEO files (series matrix, supplementary count tables, GPL annotation) into the data tree and write
PROVENANCE.md with URL, date and sha256. Download only; nothing is executed."""
import sys, os, hashlib, urllib.request, shutil
ROOT = "/mnt/max8tb/data/CoPhaser_data/cancer_circadian"
DATE = "2026-11-12"
G = "https://ftp.ncbi.nlm.nih.gov/geo/"
def stub(g): return "GSE" + g[3:-3] + "nnn"
def sha(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 22), b""): h.update(b)
    return h.hexdigest()
def fetch(url, dest):
    if os.path.exists(dest) and os.path.getsize(dest) > 0: return
    tmp = dest + ".part"
    with urllib.request.urlopen(url, timeout=300) as r, open(tmp, "wb") as f: shutil.copyfileobj(r, f, 1 << 22)
    os.rename(tmp, dest)
def run(group, gse, files, note):
    d = f"{ROOT}/{group}/{gse}"; os.makedirs(d, exist_ok=True)
    rows = []
    for rel, name in files:
        url = G + rel; fetch(url, f"{d}/{name}")
        rows.append((name, url, os.path.getsize(f"{d}/{name}"), sha(f"{d}/{name}")))
    with open(f"{d}/PROVENANCE.md", "w") as f:
        f.write(f"# {gse}\n\nSource: NCBI GEO (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc={gse}). Downloaded {DATE}. Processed matrices only (no FASTQ, no CEL).\n\n{note}\n\n| file | URL | bytes | sha256 |\n|---|---|---|---|\n")
        for r in rows: f.write(f"| {r[0]} | {r[1]} | {r[2]} | {r[3]} |\n")
    print(gse, "ok", [(r[0], r[2]) for r in rows])
def sm(g, suffix=""): return (f"series/{stub(g)}/{g}/matrix/{g}{suffix}_series_matrix.txt.gz", f"{g}{suffix}_series_matrix.txt.gz")
def su(g, n): return (f"series/{stub(g)}/{g}/suppl/{n}", n)
def gpl(p, kind="annot"):
    st = f"GPL{p[3:-3]}nnn"
    return (f"platforms/{st}/{p}/annot/{p}.annot.gz", f"{p}.annot.gz") if kind == "annot" else (f"platforms/{st}/{p}/soft/{p}_family.soft.gz", f"{p}_family.soft.gz")
if __name__ == "__main__":
    which = sys.argv[1:]
    J = {
     "GSE43292": ("paired_lesion", [sm("GSE43292")], "Carotid atheroma vs macroscopically intact tissue, 32 patients (Affymetrix HuGene 1.0 ST)."),
     "GSE57218": ("paired_lesion", [sm("GSE57218"), su("GSE57218","GSE57218_Non-normalized_data.txt.gz")], "RAAK knee/hip OA cartilage, OA vs preserved, illumina HT-12."),
     "GSE151243": ("paired_lesion", [sm("GSE151243"), su("GSE151243","GSE151243_Raw_gene_counts_matrix.txt.gz")], "Hidradenitis suppurativa lesion vs perilesion skin, RNA-seq counts."),
     "GSE192352": ("paired_lesion", [sm("GSE192352"), su("GSE192352","GSE192352_RAW_counts.csv.gz")], "Uterine leiomyoma vs own myometrium (benign tumour), RNA-seq counts."),
     "GSE10334": ("paired_lesion", [sm("GSE10334")], "Gingiva, periodontitis-affected vs unaffected sites, Affymetrix U133 Plus 2."),
     "GSE16134": ("paired_lesion", [sm("GSE16134")], "Gingiva, periodontitis-affected vs unaffected sites, Affymetrix U133 Plus 2."),
     "GSE26051": ("paired_lesion", [sm("GSE26051")], "Tendinopathy lesional vs non-lesional tendon, Affymetrix U133 Plus 2."),
     "GSE39491": ("paired_lesion", [sm("GSE39491")], "Barrett metaplasia vs matched squamous oesophagus and gastric cardia, same BE patients, Affymetrix U133A 2.0 (premalignant lesion)."),
     "GSE119717": ("paired_lesion", [sm("GSE119717"), su("GSE119717","GSE119717_Non-normalized_data.txt.gz")], "AAA: dilated vs non-dilated perivascular adipose tissue of the same subject, plus omental and subcutaneous fat (cross-listed in living_multitissue)."),
     "GSE24293": ("living_multitissue", [sm("GSE24293")], "Greenawalt/MGH RYGB liver, Agilent two-colour (log ratio vs pool)."),
     "GSE24294": ("living_multitissue", [sm("GSE24294")], "Greenawalt/MGH RYGB omental fat, Agilent two-colour (log ratio vs pool)."),
     "GSE24297": ("living_multitissue", [sm("GSE24297")], "Greenawalt/MGH RYGB subcutaneous fat, Agilent two-colour (log ratio vs pool)."),
     "GSE24335": ("living_multitissue", [sm("GSE24335")], "Greenawalt/MGH RYGB (see sample sheet for tissue)."),
     "GSE73439": ("living_multitissue", [sm("GSE73439","-GPL5175")], "Pregnancy/labour, subcutaneous vs visceral adipose of the same woman, gene-level Affymetrix HuEx (exon-level series skipped)."),
     "GSE24883": ("living_multitissue", [sm("GSE24883")], "Paired SAT/VAT, 32 women, lean/overweight/obese/obese-MetS, Agilent two-colour."),
     "GSE320593": ("living_multitissue", [sm("GSE320593"), su("GSE320593","GSE320593_AT.counts.csv.gz"), su("GSE320593","GSE320593_liver.counts.csv.gz")], "Liver and subcutaneous adipose biopsies with/without perioperative dexamethasone, surgical vs non-surgical, RNA-seq counts."),
     "GSE133979": ("living_multitissue", [sm("GSE133979"), su("GSE133979","GSE133979_Rectus_Abdominis_Muscle_variance_stabilized_counts.txt.gz"), su("GSE133979","GSE133979_Subcutaneous_adipose_variance_stabilized_counts.txt.gz")], "Rectus abdominis muscle and SAT, PDAC cachexia vs controls (cancer patients)."),
    }
    ADD = {'GSE39491': 'GPL571', 'GSE43292': 'GPL6244', 'GSE57218': 'GPL6947', 'GSE10334': 'GPL570', 'GSE16134': 'GPL570', 'GSE26051': 'GPL570', 'GSE119717': 'GPL10558', 'GSE24883': 'GPL4133', 'GSE24293': 'GPL4372', 'GSE24294': 'GPL4372', 'GSE24297': 'GPL4372', 'GSE24335': 'GPL4372'}
    for w in which:
        g, f, n = J[w]
        if w in ADD: f = f + [gpl(ADD[w])]
        if w == 'GSE73439': f = f + [gpl('GPL5175', 'soft')]
        run(g, w, f, n)
