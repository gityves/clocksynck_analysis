"""A215 / A216: checkpoint-deregulated proliferation and tumour clock coherence (registered in paper/docs/registered_experiments.md,
A215 and A216; mirrors A213).
  --stage power : mutation / copy-number status tables, cell eligibility and counts, permutation nulls (labels only) of every statistic,
                  MDE = 2.49 x null SD (status_*.csv, status_fractions.csv, eligibility.csv, counts.csv, power.csv, null_draws.csv).
                  No coherence of any real group is computed. The A216 precondition (P0 below its own pseudo 5th percentile) needs the
                  observed P0 and is evaluated in the tests stage only.
  --stage tests : observed statistics, one-sided permutation p, Holm (H215.1-3; H216.1-2), pseudo-clock companion (200 sets), S215.5,
                  attenuation A = P1 / P0 with a tumour bootstrap, descriptives, A176 components.
    python -I compute_a215.py --stage power|tests [--nperm 2000] [--npseudo 200] [--nboot 1000] [--jobs 4]
H215.1 TCGA TP53, H215.2 METABRIC TP53, H215.3 TCGA RB pathway: cells = stratum (TCGA project; METABRIC PAM50) x proliferation tertile,
T = mean over cells of coherence(altered) - coherence(intact), labels permuted within cell.
H215.4: within project, I = [coh(top) - coh(bottom prolif tertile)]_altered - [same]_intact, mean over projects, TP53 status permuted within
project x pooled proliferation tertile. S215.5: H215.1 on proliferation-matched pairs, swap null.
H216.1-2: strata = project (PAM50) x FGA tertile x TP53 status; P1 = mean over strata of coh(top) - coh(bottom proliferation tertile);
P0 = same with strata = project (PAM50) on the tumours of the eligible strata. Coherence, expression, proliferation, pseudo sets as A210 / A213."""
import argparse
import csv
import gzip
import os
import sys
import zlib
from pathlib import Path

for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "1")
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(HERE.parents[1] / "A168_endocrine" / "scripts"))
sys.path.insert(0, str(ROOT / "paper"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402
from scipy.stats import norm  # noqa: E402

import a168_lib as L  # noqa: E402
import clocksync as cs  # noqa: E402
from clocksync.decompose import components  # noqa: E402
from common import data  # noqa: E402
from common.paths import DATA  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--stage", required=True, choices=["power", "tests"])
ap.add_argument("--nperm", type=int, default=2000); ap.add_argument("--npseudo", type=int, default=200)
ap.add_argument("--nboot", type=int, default=1000); ap.add_argument("--jobs", type=int, default=4)
args = ap.parse_args()
OUT = HERE.parent / "output"; OUT.mkdir(exist_ok=True)
A210 = HERE.parents[1] / "A210_myc_coherence" / "output"
A213 = HERE.parents[1] / "A213_heterogeneity_coherence" / "output"
CACHE = A210 / "tcga_cache"  # per-project clock counts and proliferation (A210; open-access TCGA)
K, MIN_ARM, MDE_MAX = norm.ppf(0.95) + norm.ppf(0.80), 10, 0.25
SEED = {"A215": zlib.crc32(b"A215") % (2 ** 31), "A216": zlib.crc32(b"A216") % (2 ** 31)}
NCHUNK = 20
PROJ_A97 = ["BLCA", "BRCA", "COAD", "KIRC", "KIRP", "LIHC", "LUAD", "OV", "PRAD", "THCA", "UCEC"]
PROJ_A172 = ["CESC", "ESCA", "GBM", "HNSC", "LGG", "LUSC", "PAAD", "PCPG", "READ", "SARC", "STAD", "TGCT"]
PROJECTS = PROJ_A97 + PROJ_A172
GENES, tertile, LN2 = list(L.GENES), L.A166.tertile, L.LN2
pl = lambda nl, o: cs.stats.perm_p(nl, o, "less", tie=0.0)
say = lambda *a: print(pd.Timestamp.now().strftime("%F %T"), *a, flush=True)
XENA = DATA / "tcga" / "xena_gistic" / "Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz"
MC3 = DATA / "tcga" / "mc3" / "mc3.v0.2.8.PUBLIC.maf.gz"
MBD = DATA / "metabric" / "brca_metabric"
ABS = DATA / "tcga" / "absolute_pancanatlas" / "TCGA_mastercalls.abs_tables_JSedit.fixed.txt"
NONSILENT = {"Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Splice_Site",
             "Translation_Start_Site", "Nonstop_Mutation"}
TRUNC = {"Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site"}
MISS = {"Missense_Mutation", "In_Frame_Del", "In_Frame_Ins"}
MC3_FILTER_OK = {"PASS", "wga", "native_wga_mix"}  # PanCanAtlas MC3 recommendation (choice left open by the registration)
GENE3 = ("TP53", "RB1", "CDKN2A")


# ---- status tables (text parse only) ----------------------------------------------------------------------------------------------
def tp53_class(vcs):
    return "trunc" if vcs & TRUNC else "missense" if vcs & MISS else "other"


def maf_status(path, opener, sample_col, sample_fn, filt=False):
    """Per-sample non-silent TP53 / RB1 / CDKN2A mutation flags and the TP53 variant class; sequenced = every barcode present in the MAF."""
    seq, hits, n_all, n_flt = set(), {}, 0, 0
    cols = ["Hugo_Symbol", "Variant_Classification", sample_col] + (["FILTER"] if filt else [])
    with opener(path, "rt") as f:
        for ch in pd.read_csv(f, sep="\t", usecols=cols, dtype=str, chunksize=500_000, comment="#", low_memory=False):
            ch["s"] = ch[sample_col].map(sample_fn); seq.update(ch.s.unique())
            h = ch[ch.Hugo_Symbol.isin(GENE3) & ch.Variant_Classification.isin(NONSILENT)]
            n_all += len(h)
            if filt:
                h = h[h.FILTER.fillna("").map(lambda x: set(x.split(";")) <= MC3_FILTER_OK)]
            n_flt += len(h)
            for s, g, v in zip(h.s, h.Hugo_Symbol, h.Variant_Classification):
                hits.setdefault(s, {}).setdefault(g, set()).add(v)
    rows = [dict(sample=s, sequenced=1, tp53_mut=int("TP53" in hits.get(s, {})), rb1_mut=int("RB1" in hits.get(s, {})),
                 cdkn2a_mut=int("CDKN2A" in hits.get(s, {})), tp53_cls=tp53_class(hits.get(s, {}).get("TP53", set())) if "TP53" in hits.get(s, {}) else "")
            for s in sorted(seq)]
    say(path.name, "sequenced", len(seq), "non-silent hits (before / after FILTER)", n_all, n_flt)
    return pd.DataFrame(rows).set_index("sample")


def gistic_rows(path, opener, skip):
    """GISTIC state of TP53, RB1, CDKN2A per sample (text parse; first column = gene symbol)."""
    out = {}
    with opener(path, "rt") as f:
        r = csv.reader(f, delimiter="\t"); head = next(r)
        for row in r:
            if row[0] in GENE3:
                out[row[0]] = pd.Series(pd.to_numeric(row[skip:], errors="coerce"), index=head[skip:])
    return pd.DataFrame(out)


def status_table(name):
    f = OUT / f"status_{name}.csv"
    if not f.exists():
        if name == "tcga":
            m = maf_status(MC3, gzip.open, "Tumor_Sample_Barcode", lambda x: x[:15], filt=True)
            g = gistic_rows(XENA, gzip.open, 1)
        else:
            m = maf_status(MBD / "data_mutations.txt", open, "Tumor_Sample_Barcode", lambda x: x)
            ids = [l for l in open(MBD / "case_lists" / "cases_sequenced.txt") if l.startswith("case_list_ids:")][0].split(":", 1)[1].split()
            m = m.reindex(sorted(set(ids) | set(m.index))); m["sequenced"] = m.sequenced.fillna(0).astype(int) | m.index.isin(ids).astype(int)
            m = m[m.sequenced == 1].fillna({"tp53_mut": 0, "rb1_mut": 0, "cdkn2a_mut": 0, "tp53_cls": ""})
            g = gistic_rows(MBD / "data_cna.txt", open, 2)
        g.columns = [c.lower() + "_gistic" for c in g.columns]
        d = m.join(g, how="left"); d.index.name = "sample"
        # altered: non-silent mutation or GISTIC -2; intact: none of these and a copy-number call; otherwise not defined
        for gene, pre in (("tp53", ["tp53"]), ("rb", ["rb1", "cdkn2a"])):
            mut = np.column_stack([d[f"{x}_mut"] == 1 for x in pre]).any(1)
            dele = np.column_stack([d[f"{x}_gistic"] == -2 for x in pre]).any(1)
            call = np.column_stack([d[f"{x}_gistic"].notna() for x in pre]).all(1)
            d[gene] = np.where(mut | dele, 1.0, np.where(call, 0.0, np.nan))
        d["tp53_class"] = np.where(d.tp53 != 1, "", np.where(d.tp53_mut == 1, d.tp53_cls, "deletion"))  # trunc / missense / other / deletion only
        d.to_csv(f)
    return pd.read_csv(f, index_col=0, keep_default_na=False, na_values=[""]).drop(columns="tp53_cls")


# ---- cohort tables ----------------------------------------------------------------------------------------------------------------
def primaries(B, p):
    s = B.samples_of(p)
    return s[s.sample_type == "Primary Tumor"].sort_values("barcode").drop_duplicates("patient")


def load_lc(B, p):
    """log2 CPM of genes with mean CPM >= 1 for the primaries of p (A213 `load_lc`)."""
    ct = B.counts_of(p, primaries(B, p)); depth = ct.sum(axis=1)
    lc = np.log2(ct.div(depth, axis=0) * 1e6 + 1)
    return lc.loc[:, (2 ** lc - 1).mean() >= 1].astype(np.float32)


def tcga_table(p, st, ft, ab):
    d = pd.read_parquet(CACHE / f"{p}.parquet").reset_index(drop=True)
    d["fga"] = ft.reindex(d["sample"]).to_numpy(); d["purity"] = ab["purity"].reindex(d["sample"]).to_numpy()
    for c in ("tp53", "rb", "tp53_class"):
        d[c] = st[c].reindex(d["sample"]).to_numpy()
    d["tp53_class"] = d.tp53_class.where(d.tp53.notna(), None)
    return d


def metabric(st, fm):
    m, X = L.A166.load("metabric"); m = m[m.pam50.notna()].copy(); X = X.loc[m.index]
    m = m.reset_index(drop=True); X = X.reset_index(drop=True); m["stratum"] = m.pam50.astype(str)
    m["fga"] = fm.reindex(m.patient).to_numpy()
    for c in ("tp53", "rb", "tp53_class"):
        m[c] = st[c].reindex(m.patient).to_numpy()
    return m, X


# ---- designs: each has draw(rng) -> state and values(Y, wf, state=None) -> one value per unit -------------------------------------
def diffs(Y, wf, cells, perms=None):
    out = []
    for k, c in enumerate(cells):
        l = c["lab"] if perms is None else c["lab"][perms[k]]
        out.append(L.coh(Y, c["idx"][l == 1], wf) - L.coh(Y, c["idx"][l == 0], wf))
    return np.array(out)


class CellD:
    """Labelled arms within cells (label 1 minus label 0); labels permuted within cell."""
    def __init__(s, cells): s.cells = cells
    def draw(s, rng): return [rng.permutation(len(c["lab"])) for c in s.cells]
    def values(s, Y, wf, st=None): return diffs(Y, wf, s.cells, st)


class SwapD:
    """Matched pairs (a = altered, z = intact); the two members swapped with probability 1/2."""
    def __init__(s, cells): s.cells = cells
    def draw(s, rng): return [rng.random(len(c["a"])) < 0.5 for c in s.cells]
    def values(s, Y, wf, st=None):
        out = []
        for k, c in enumerate(s.cells):
            f = np.zeros(len(c["a"]), bool) if st is None else st[k]
            out.append(L.coh(Y, np.where(f, c["z"], c["a"]), wf) - L.coh(Y, np.where(f, c["a"], c["z"]), wf))
        return np.array(out)


class InterD:
    """H215.4 for one project: pos (rows of the project table), status (1 altered / 0 intact), prolif; the status is permuted within the pooled
    proliferation tertile `grp`; tertiles are recomputed within project x status."""
    def __init__(s, pos, status, prolif):
        s.pos, s.status, s.prolif = pos, status.astype(int), prolif; s.grp = tertile(prolif)
    def draw(s, rng):
        o = s.status.copy()
        for g in np.unique(s.grp):
            i = np.flatnonzero(s.grp == g); o[i] = rng.permutation(o[i])
        return o
    def arms(s, status):
        out = {}
        for k in (1, 0):
            i = np.flatnonzero(status == k); t = tertile(s.prolif[i]) if len(i) >= 3 else np.full(len(i), 1)
            out[k] = (s.pos[i[t == 2]], s.pos[i[t == 0]])
        return out
    def values(s, Y, wf, st=None):
        a = s.arms(s.status if st is None else st)
        d = {k: L.coh(Y, a[k][0], wf) - L.coh(Y, a[k][1], wf) for k in (1, 0)}
        return np.array([d[1] - d[0]])


def cells_status(m, col, sub=None):
    """Stratum x proliferation-tertile cells of the tumours of m with a status `col` (positions of m kept); arms altered (label 1) / intact (0).
    `sub` (boolean mask over m) restricts the altered arm. Returns (eligibility rows, eligible cells)."""
    m = m[m[col].notna()].copy(); m["t"] = m.groupby("stratum").prolif.transform(tertile); rows, cells = [], []
    for (s, t), g in m.groupby(["stratum", "t"]):
        if sub is not None:
            g = g[(g[col] == 0) | sub.reindex(g.index).fillna(False).to_numpy()]
        n1, n0 = int((g[col] == 1).sum()), int((g[col] == 0).sum()); ok = min(n1, n0) >= MIN_ARM
        rows.append(dict(stratum=s, cell=int(t), n=len(g), n_arm1=n1, n_arm0=n0, eligible=ok))
        if ok:
            cells.append(dict(key=f"{s}|{int(t)}", stratum=s, idx=g.index.to_numpy(), lab=g[col].to_numpy().astype(int)))
    return pd.DataFrame(rows), cells


def cells_inter(m, col="tp53"):
    """H215.4: one InterD per project with >= MIN_ARM tumours in each of the four arms (top / bottom tertile within project x status)."""
    m = m[m[col].notna()]; rows, des = [], {}
    for s, g in m.groupby("stratum"):
        d = InterD(g.index.to_numpy(), g[col].to_numpy(), g.prolif.to_numpy(float)); a = d.arms(d.status)
        n = {k: (len(a[k][0]), len(a[k][1])) for k in (1, 0)}; ok = min(min(n[1]), min(n[0])) >= MIN_ARM
        rows.append(dict(stratum=s, cell=-1, n=len(g), n_arm1=min(n[1]), n_arm0=min(n[0]), eligible=ok))
        if ok:
            des[s] = d
    return pd.DataFrame(rows), des


def match(pa, pz, caliper):
    """1:1 nearest-neighbour matching without replacement on the continuous score: all pairs within the caliper taken in increasing
    distance (ties by position), each tumour used once (A213)."""
    D = np.abs(pa[:, None] - pz[None, :]); i, j = np.nonzero(D <= caliper)
    o = np.lexsort((j, i, D[i, j])); ua, uz, A, Z = set(), set(), [], []
    for k in o:
        if i[k] not in ua and j[k] not in uz:
            ua.add(i[k]); uz.add(j[k]); A.append(i[k]); Z.append(j[k])
    return np.array(A, int), np.array(Z, int)


def cells_matched(m, cells, col, sd):
    """S215.5: within each eligible H215.1 cell, altered tumours matched 1:1 to intact ones on prolif (caliper 0.1 x `sd`); >= MIN_ARM pairs."""
    rows, out = [], []
    for c in cells:
        g = m.loc[c["idx"]]; ga, gz = g[g[col] == 1], g[g[col] == 0]
        a, z = match(ga.prolif.to_numpy(float), gz.prolif.to_numpy(float), 0.1 * sd)
        rows.append(dict(stratum=c["stratum"], cell=c["key"], n_arm1=len(ga), n_arm0=len(gz), n_pairs=len(a), eligible=len(a) >= MIN_ARM))
        if len(a) >= MIN_ARM:
            out.append(dict(key=c["key"], stratum=c["stratum"], a=ga.index.to_numpy()[a], z=gz.index.to_numpy()[z]))
    return pd.DataFrame(rows), out


def cells_state(m):
    """A216: strata = stratum x FGA tertile (within stratum, over tumours with TP53 status and FGA) x TP53 status; within each stratum top
    (label 1) vs bottom (label 0) proliferation tertile, >= MIN_ARM per arm. Returns (eligibility rows, P1 cells, P0 cells, strata with all
    tumours for the bootstrap). P0: strata = stratum only, on the tumours of the eligible P1 strata."""
    m = m[m.tp53.notna() & m.fga.notna()].copy(); m["f"] = m.groupby("stratum").fga.transform(tertile); rows, cells, strata = [], [], []
    for (s, f, t), g in m.groupby(["stratum", "f", "tp53"]):
        pos = g.index.to_numpy(); e = tertile(g.prolif.to_numpy(float)) if len(g) >= 3 else np.full(len(g), 1)
        top, bot = pos[e == 2], pos[e == 0]; ok = min(len(top), len(bot)) >= MIN_ARM
        key = f"{s}|f{int(f)}|tp53_{int(t)}"
        rows.append(dict(stratum=s, cell=key, n=len(g), n_arm1=len(top), n_arm0=len(bot), eligible=ok))
        if ok:
            cells.append(dict(key=key, stratum=s, idx=np.r_[top, bot], lab=np.r_[np.ones(len(top), int), np.zeros(len(bot), int)]))
            strata.append(dict(key=key, stratum=s, pos=pos))
    p0 = []
    for s in sorted({c["stratum"] for c in strata}):
        pos = np.concatenate([c["pos"] for c in strata if c["stratum"] == s]); e = tertile(m.prolif.loc[pos].to_numpy(float))
        top, bot = pos[e == 2], pos[e == 0]
        p0.append(dict(key=s, stratum=s, idx=np.r_[top, bot], lab=np.r_[np.ones(len(top), int), np.zeros(len(bot), int)]))
    return pd.DataFrame(rows), cells, p0, strata


# ---- nulls ------------------------------------------------------------------------------------------------------------------------
def null_chunk(blocks, seed, n):
    """blocks: (Y, wf, design); the statistic of a draw is the mean over all units of all blocks."""
    rng = np.random.default_rng(seed); out = np.empty(n)
    for i in range(n):
        out[i] = np.mean(np.concatenate([d.values(Y, wf, d.draw(rng)) for Y, wf, d in blocks]))
    return out


def draws(fn, blocks, seed, tag, n):
    sizes = [n // NCHUNK + (j < n % NCHUNK) for j in range(NCHUNK)]
    return np.concatenate(Parallel(n_jobs=args.jobs)(delayed(fn)(blocks, np.random.SeedSequence([seed, tag, j]), s) for j, s in enumerate(sizes)))


# ---- build everything (no coherence of a real group) ------------------------------------------------------------------------------
st_tc, st_mb = status_table("tcga"), status_table("metabric")
ft = pd.read_csv(A213 / "fga_tcga.csv", index_col=0).fga; fm = pd.read_csv(A213 / "fga_metabric.csv", index_col=0).fga  # A213 FGA (GISTIC != 0)
ab = pd.read_csv(ABS, sep="\t"); ab = ab[ab["call status"] == "called"].drop_duplicates("array").set_index("array")
mb, Xmb = metabric(st_mb, fm)
tc = {p: tcga_table(p, st_tc, ft, ab) for p in PROJECTS}
Ymb = LN2 * Xmb[GENES].to_numpy(float)
Yof = lambda counts, genes: np.log((counts[list(genes)].to_numpy(float) + 0.5) / counts.depth.to_numpy(float)[:, None])
Ytc = {p: Yof(tc[p], GENES) for p in PROJECTS}
TAB = {"TCGA": tc, "METABRIC": {"METABRIC": mb}}
sd_prolif = {p: float(tc[p].loc[tc[p].tp53.notna(), "prolif"].std()) for p in PROJECTS}; sd_prolif["METABRIC"] = float(mb.loc[mb.tp53.notna(), "prolif"].std())

# H215.1-3: cells per stratum
HYP215 = {"H215.1": ("TCGA", "tp53"), "H215.2": ("METABRIC", "tp53"), "H215.3": ("TCGA", "rb")}
el, cells_h = [], {}
for h, (coh_, col) in HYP215.items():
    cells_h[h] = {}
    for k, tab in TAB[coh_].items():
        tab = tab.assign(stratum=k) if coh_ == "TCGA" else tab
        e, c = cells_status(tab, col); e.insert(0, "hypothesis", h); el.append(e); cells_h[h][k] = c
# H215.4 and S215.5
inter, inter_el = {}, []
for p in PROJECTS:
    e, d = cells_inter(tc[p].assign(stratum=p)); e.insert(0, "hypothesis", "H215.4"); inter_el.append(e); inter.update({p: d[p]} if p in d else {})
el += inter_el
s5 = {}
for p in PROJECTS:
    e, c = cells_matched(tc[p], cells_h["H215.1"][p], "tp53", sd_prolif[p]); e.insert(0, "hypothesis", "S215.5"); el.append(e); s5[p] = c
# H216
HYP216 = {"H216.1": "TCGA", "H216.2": "METABRIC"}
st216, p1_cells, p0_cells, strata216 = {}, {}, {}, {}
for h, coh_ in HYP216.items():
    p1_cells[h], p0_cells[h], strata216[h] = {}, {}, {}
    for k, tab in TAB[coh_].items():
        tab = tab.assign(stratum=k) if coh_ == "TCGA" else tab
        e, c, p0, sr = cells_state(tab); e.insert(0, "hypothesis", h); el.append(e)
        p1_cells[h][k], p0_cells[h][k], strata216[h][k] = c, p0, sr
elig = pd.concat(el, ignore_index=True); elig.to_csv(OUT / "eligibility.csv", index=False)

ncells = lambda cc: int(sum(len(c["lab"]) for c in cc))
cnt = []
for h, (coh_, col) in HYP215.items():
    for k, c in cells_h[h].items():
        tab = TAB[coh_][k]; cnt.append(dict(hypothesis=h, stratum=k, n_with_status=int(tab[col].notna().sum()), n_altered=int((tab[col] == 1).sum()),
                                            n_cells=len(c), n_tumours=ncells(c)))
for p in PROJECTS:
    t = tc[p][tc[p].tp53.notna()]
    cnt.append(dict(hypothesis="H215.4", stratum=p, n_with_status=len(t), n_altered=int((t.tp53 == 1).sum()), n_cells=int(p in inter), n_tumours=len(t) if p in inter else 0))
    cnt.append(dict(hypothesis="S215.5", stratum=p, n_with_status=len(t), n_altered=int((t.tp53 == 1).sum()), n_cells=len(s5[p]), n_tumours=2 * int(sum(len(c["a"]) for c in s5[p]))))
for h, coh_ in HYP216.items():
    for k, tab in TAB[coh_].items():
        t = tab[tab.tp53.notna() & tab.fga.notna()]
        cnt.append(dict(hypothesis=h, stratum=k, n_with_status=len(t), n_altered=int((t.tp53 == 1).sum()), n_cells=len(p1_cells[h][k]), n_tumours=ncells(p1_cells[h][k]),
                        n_tumours_P0=sum(len(s["pos"]) for s in strata216[h][k])))
pd.DataFrame(cnt).to_csv(OUT / "counts.csv", index=False)
fr = []
for coh_, tabs in TAB.items():
    for k, t in tabs.items():
        fr.append(dict(cohort=coh_, stratum=k, n_rows=len(t), n_sequenced_with_tp53_status=int(t.tp53.notna().sum()), tp53_altered=int((t.tp53 == 1).sum()),
                       n_with_rb_status=int(t.rb.notna().sum()), rb_altered=int((t.rb == 1).sum())))
fr = pd.DataFrame(fr); fr["tp53_fraction"] = fr.tp53_altered / fr.n_sequenced_with_tp53_status; fr["rb_fraction"] = fr.rb_altered / fr.n_with_rb_status
fr.to_csv(OUT / "status_fractions.csv", index=False)

# blocks of (Y, wf, design) per hypothesis
blocks = {}
for h, (coh_, col) in HYP215.items():
    blocks[h] = [(Ymb, L.WF, CellD(cells_h[h]["METABRIC"]))] if coh_ == "METABRIC" else [(Ytc[p], L.WF, CellD(cells_h[h][p])) for p in PROJECTS if cells_h[h][p]]
blocks["H215.4"] = [(Ytc[p], L.WF, inter[p]) for p in PROJECTS if p in inter]
blocks["S215.5"] = [(Ytc[p], L.WF, SwapD(s5[p])) for p in PROJECTS if s5[p]]
blocks["H216.1"] = [(Ytc[p], L.WF, CellD(p1_cells["H216.1"][p])) for p in PROJECTS if p1_cells["H216.1"][p]]
blocks["H216.2"] = [(Ymb, L.WF, CellD(p1_cells["H216.2"]["METABRIC"]))] if p1_cells["H216.2"]["METABRIC"] else []
HYPS = ["H215.1", "H215.2", "H215.3", "H215.4", "S215.5", "H216.1", "H216.2"]
amend = lambda h: "A216" if h.startswith("H216") else "A215"
COHORT = {"H215.1": "TCGA", "H215.2": "METABRIC", "H215.3": "TCGA", "H215.4": "TCGA", "S215.5": "TCGA", "H216.1": "TCGA", "H216.2": "METABRIC"}
MEAS = {"H215.1": "tp53", "H215.2": "tp53", "H215.3": "rb", "H215.4": "tp53 x prolif", "S215.5": "tp53 matched", "H216.1": "P1", "H216.2": "P1"}
nunits = lambda h: int(sum(len(d.cells) if hasattr(d, "cells") else 1 for _, _, d in blocks[h]))
ntum = lambda h: int(sum(ncells(d.cells) if isinstance(d, CellD) else 2 * sum(len(c["a"]) for c in d.cells) if isinstance(d, SwapD) else len(d.pos)
                         for _, _, d in blocks[h]))
say({h: (len(blocks[h]), nunits(h), ntum(h)) for h in HYPS})

if args.stage == "power":
    nl = {}
    for i, h in enumerate(HYPS):
        if blocks[h]:
            nl[h] = draws(null_chunk, blocks[h], SEED[amend(h)], i, args.nperm); say(h, "null done", nl[h].mean().round(5), nl[h].std(ddof=1).round(5))
    pd.DataFrame(nl).to_csv(OUT / "null_draws.csv", index=False)
    pw = pd.DataFrame([dict(hypothesis=h, cohort=COHORT[h], measure=MEAS[h], n_projects=len(blocks[h]), n_units=nunits(h), n_tumours=ntum(h),
                            n_perm=args.nperm, null_mean=nl[h].mean() if h in nl else np.nan, null_sd=nl[h].std(ddof=1) if h in nl else np.nan,
                            MDE_80pct=K * nl[h].std(ddof=1) if h in nl else np.nan) for h in HYPS])
    pw["scored"] = pw.MDE_80pct <= MDE_MAX
    pw["precondition"] = np.where(pw.hypothesis.str.startswith("H216"), "P0 below its pseudo 5th percentile: evaluated in the tests stage", "")
    pw.to_csv(OUT / "power.csv", index=False); print(pw.round(4).to_string(), flush=True)
    say("done power"); sys.exit(0)

# ---- tests (not run before the dated power note) ----------------------------------------------------------------------------------
pw = pd.read_csv(OUT / "power.csv").set_index("hypothesis"); nl = pd.read_csv(OUT / "null_draws.csv"); assert len(nl) == args.nperm
obs_u = {h: np.concatenate([d.values(Y, wf) for Y, wf, d in blocks[h]]) if blocks[h] else np.array([np.nan]) for h in HYPS}
obs = {h: float(obs_u[h].mean()) for h in HYPS}
# P0 on the tumours of the eligible P1 strata
p0_blocks = {h: ([(Ytc[p], L.WF, CellD(p0_cells[h][p])) for p in PROJECTS if p0_cells[h][p]] if h == "H216.1" else [(Ymb, L.WF, CellD(p0_cells[h]["METABRIC"]))])
             for h in HYP216}
obs_p0 = {h: float(np.mean(np.concatenate([d.values(Y, wf) for Y, wf, d in p0_blocks[h]]))) for h in HYP216}

# pseudo-clock companion: 200 sets per project (A213 construction), one family per amendment, shared by the hypotheses of the amendment
B = data._builder("tcga_paired")
sets_mb = pd.read_csv(L.A166.OUT / "pseudo_sets_metabric.csv").to_numpy()[: args.npseudo].tolist()
pseudo_mb = lambda design: np.array([np.mean(design.values(LN2 * Xmb[s].to_numpy(float), L.A166.wf_for(s))) for s in sets_mb])
pseudo = {"H215.2": pseudo_mb(CellD(cells_h["H215.2"]["METABRIC"])), "H216.2": pseudo_mb(CellD(p1_cells["H216.2"]["METABRIC"])),
          "H216.2_P0": pseudo_mb(CellD(p0_cells["H216.2"]["METABRIC"]))}
pc = {h: [] for h in ("H215.1", "H215.3", "H215.4", "H216.1", "H216.1_P0")}
for j, p in enumerate(PROJECTS):
    lc = load_lc(B, p); t = primaries(B, p)
    assert (t.patient.to_numpy() == tc[p].patient.to_numpy()).all()
    lc = lc.set_axis(t.patient.to_numpy()); dep = tc[p].depth.to_numpy(float)[:, None]
    des = {"H215.1": CellD(cells_h["H215.1"][p]), "H215.3": CellD(cells_h["H215.3"][p]), "H215.4": inter.get(p), "H216.1": CellD(p1_cells["H216.1"][p]),
           "H216.1_P0": CellD(p0_cells["H216.1"][p])}
    for amd, hs in (("A215", ("H215.1", "H215.3", "H215.4")), ("A216", ("H216.1", "H216.1_P0"))):
        if not any(des[h] is not None and getattr(des[h], "cells", [1]) for h in hs):
            continue
        sets = L.pseudo_sets(lc, L.A166.pool_genes(lc), GENES, args.npseudo, np.random.default_rng(np.random.SeedSequence([SEED[amd], 1000 + j])))
        for h in hs:
            if des[h] is not None and getattr(des[h], "cells", [1]):
                pc[h].append(np.array([des[h].values(np.log((2 ** lc[list(g)].to_numpy(float) - 1) / 1e6 + 0.5 / dep), L.wf_pseudo(L.WF, g)) for g in sets]))
    say(p, "pseudo done")
for h in pc:
    pseudo[h] = np.concatenate(pc[h], axis=1).mean(axis=1)
pd.DataFrame(pseudo).to_csv(OUT / "pseudo_T.csv", index=False)

# S215.5 pseudo is not registered (companion is for T and I only)
q05 = lambda h: float(np.quantile(pseudo[h], 0.05))
rows = []
for h in HYPS:
    sc = bool(pw.loc[h, "scored"]); r = dict(hypothesis=h, cohort=COHORT[h], measure=MEAS[h], scored_power=sc, stat=obs[h], null_mean=pw.loc[h, "null_mean"],
                                              null_sd=pw.loc[h, "null_sd"], MDE_80pct=pw.loc[h, "MDE_80pct"], excess=obs[h] - pw.loc[h, "null_mean"])
    if h.startswith("H216"):
        r["P0"] = obs_p0[h]; r["P0_pseudo_q05"] = q05(h + "_P0"); r["P0_below_pseudo_q05"] = bool(obs_p0[h] < q05(h + "_P0"))
        sc = sc and r["P0_below_pseudo_q05"]
    r["scored"] = sc; r["p_one_lower"] = pl(nl[h].to_numpy(), obs[h]) if sc or h in ("H215.4", "S215.5") else np.nan
    if h in pseudo:
        r.update(pseudo_percentile=100 * np.mean(pseudo[h] <= obs[h]), pseudo_q05=q05(h), below_pseudo_q05=bool(obs[h] < q05(h)))
    rows.append(r)
res = pd.DataFrame(rows)


def holm(df, hs):
    s = df.hypothesis.isin(hs) & df.scored; ps = df.p_one_lower[s].to_numpy(); o = np.argsort(ps)
    adj = np.empty(len(ps)); adj[o] = np.minimum(1, np.maximum.accumulate(np.sort(ps) * (len(ps) - np.arange(len(ps)))))
    df.loc[s, "p_holm"] = adj


res["p_holm"] = np.nan; holm(res, ["H215.1", "H215.2", "H215.3"]); holm(res, ["H216.1", "H216.2"])
fam = res.hypothesis.isin(["H215.1", "H215.2", "H215.3", "H216.1", "H216.2"])
res["verdict"] = np.where(fam & ~res.scored, "not scored", np.where(fam, np.where(res.p_holm < 0.05, "pass", "fail"), ""))
i4 = res.hypothesis == "H215.4"
res.loc[i4, "verdict"] = np.where(~res.scored[i4], "not scored", np.where(res.p_one_lower[i4] < 0.05, "pass", "fail"))
res["specificity"] = np.where((res.verdict == "pass") & res.below_pseudo_q05.notna(), np.where(res.below_pseudo_q05.fillna(False), "clock-specific", "not specific to the clock"), "")
res.to_csv(OUT / "tests.csv", index=False); print(res.round(4).to_string(), flush=True)
say("tests done")


# ---- A216 attenuation A = P1 / P0 (tumour bootstrap within strata) ---------------------------------------------------------------
def boot_chunk(bl, seed, n):
    rng = np.random.default_rng(seed); out = np.empty((n, 2))
    for i in range(n):
        v1, v0 = [], []
        for Y, wf, prolif, strata in bl:
            allp = {}
            for s in strata:
                pos = rng.choice(s["pos"], len(s["pos"])); t = tertile(prolif[pos]); v1.append(L.coh(Y, pos[t == 2], wf) - L.coh(Y, pos[t == 0], wf))
                allp.setdefault(s["stratum"], []).append(pos)
            for s, ps in allp.items():
                pos = np.concatenate(ps); t = tertile(prolif[pos]); v0.append(L.coh(Y, pos[t == 2], wf) - L.coh(Y, pos[t == 0], wf))
        out[i] = np.mean(v1), np.mean(v0)
    return out


att = []
for h, coh_ in HYP216.items():
    if coh_ == "TCGA":
        bl = [(Ytc[p], L.WF, tc[p].prolif.to_numpy(float), strata216[h][p]) for p in PROJECTS if strata216[h][p]]
    else:
        bl = [(Ymb, L.WF, mb.prolif.to_numpy(float), strata216[h]["METABRIC"])]
    sizes = [args.nboot // NCHUNK + (j < args.nboot % NCHUNK) for j in range(NCHUNK)]
    bs = np.vstack(Parallel(n_jobs=args.jobs)(delayed(boot_chunk)(bl, np.random.SeedSequence([SEED["A216"], 99, j, 1 if coh_ == "TCGA" else 2]), s)
                                              for j, s in enumerate(sizes)))
    A = bs[:, 0] / bs[:, 1]; a = obs[h] / obs_p0[h]
    att.append(dict(hypothesis=h, P1=obs[h], P0=obs_p0[h], A=a, A_lo=np.percentile(A, 2.5), A_hi=np.percentile(A, 97.5), n_boot=len(A),
                    P1_boot_sd=bs[:, 0].std(ddof=1), P0_boot_sd=bs[:, 1].std(ddof=1), frac_boot_P0_nonneg=float(np.mean(bs[:, 1] >= 0))))
att = pd.DataFrame(att)
res216 = res.set_index("hypothesis"); att["reading"] = [
    "not scored" if not res216.loc[h, "scored"] else
    "persists" if (res216.loc[h, "verdict"] == "pass" and res216.loc[h, "specificity"] == "clock-specific" and a.A >= 0.5) else
    "carried by state variables" if (res216.loc[h, "verdict"] == "fail" and a.A_hi < 0.5) else "inconclusive" for h, a in att.set_index("hypothesis", drop=False).iterrows()]
att.to_csv(OUT / "attenuation.csv", index=False); print(att.round(4).to_string(), flush=True)
say("attenuation done")

# S215.5 pairs report (the statistic and its swap-null p are in tests.csv; pair counts in eligibility.csv)

# ---- descriptives -----------------------------------------------------------------------------------------------------------------
def parts(Y, wf, top, bot):
    """(dR, dlog_plane, dlog_off) of top vs bottom, each group centred (A176 Contrast.parts)."""
    c = []
    for i in (top, bot):
        y = Y[i]; y = y - y.mean(0); c.append(components(y, np.zeros(y.shape[1]), 1, wf))
    return [c[0]["coherence"] - c[1]["coherence"], c[0]["log_plane"] - c[1]["log_plane"], c[0]["log_off"] - c[1]["log_off"]]


def desc_cells(h, cells, tab, Y, k, extra):
    out = []
    for c in cells:
        top, bot = c["idx"][c["lab"] == 1], c["idx"][c["lab"] == 0]; pr = parts(Y, L.WF, top, bot); pro = tab.prolif.to_numpy()
        r = dict(hypothesis=h, stratum=k, cell=c["key"], n_top=len(top), n_bottom=len(bot), T_cell=L.coh(Y, top, L.WF) - L.coh(Y, bot, L.WF),
                 dR_components=pr[0], dlog_plane=pr[1], dlog_off=pr[2], ds=pr[1] - pr[2], prolif_top_minus_bottom=pro[top].mean() - pro[bot].mean())
        for name in extra:
            if name in tab:
                r[name + "_top_minus_bottom"] = np.nanmean(tab[name].to_numpy(float)[top]) - np.nanmean(tab[name].to_numpy(float)[bot])
        out.append(r)
    return out


dc = []
for h, (coh_, col) in HYP215.items():
    for k, c in cells_h[h].items():
        tab = tc[k] if coh_ == "TCGA" else mb
        dc += desc_cells(h, c, tab, Ytc[k] if coh_ == "TCGA" else Ymb, k, ["purity"])
for h, coh_ in HYP216.items():
    for k, c in p1_cells[h].items():
        tab = tc[k] if coh_ == "TCGA" else mb
        dc += desc_cells(h, c, tab, Ytc[k] if coh_ == "TCGA" else Ymb, k, ["fga", "tp53", "purity"])
dc = pd.DataFrame(dc); dc.to_csv(OUT / "descriptive_cells.csv", index=False)
num = [c for c in dc.columns if c in ("T_cell", "dR_components", "dlog_plane", "dlog_off", "ds") or c.endswith("_top_minus_bottom")]
dsum = dc.groupby("hypothesis")[num].mean(); dsum.insert(0, "n_cells", dc.groupby("hypothesis").size()); dsum.to_csv(OUT / "descriptive.csv"); print(dsum.round(4).to_string(), flush=True)
pp = dc.groupby(["hypothesis", "stratum"]).agg(n_cells=("T_cell", "size"), T=("T_cell", "mean")).reset_index(); pp.to_csv(OUT / "descriptive_per_project.csv", index=False)

# truncating vs missense TP53 (descriptive; same cells, altered arm restricted to the class, intact arm unchanged)
cls = []
for coh_ in ("TCGA", "METABRIC"):
    for kind in ("trunc", "missense"):
        v = []
        for k, tab in TAB[coh_].items():
            tab = tab.assign(stratum=k) if coh_ == "TCGA" else tab
            _, c = cells_status(tab, "tp53", sub=(tab.tp53_class == kind))
            v += list(diffs(Ytc[k] if coh_ == "TCGA" else Ymb, L.WF, c))
        cls.append(dict(cohort=coh_, altered_class=kind, n_cells=len(v), T=float(np.mean(v)) if v else np.nan))
pd.DataFrame(cls).to_csv(OUT / "descriptive_tp53_class.csv", index=False); print(pd.DataFrame(cls).round(4).to_string(), flush=True)
say("done tests")
