"""A224: does cis copy number of the clock genes carry the TP53 and subtype coherence contrasts? (registered in paper/docs/registered_experiments.md, A224;
rebuilds A219 H219.2 and A215 H215.1 on dosage-adjusted expression.)
  --stage power : dosage fits (dosage.csv: beta_g with SE, CN x phase interaction, share of between-tumour variance carried by CN), eligibility / counts, and
                  permutation nulls (labels only, on adjusted expression) of D_adj (H224.1), T_adj (H224.2) and the phase-adjusted PER1 / PER2 mutant - wild-type
                  difference (H224.3); MDE = 2.49 x null SD (power.csv, null_draws.csv). No D_adj, T_adj, mutant-vs-wild-type difference or attenuation is computed.
  --stage tests : observed D_adj / T_adj / H224.3 with permutation p, pseudo-clock companions (200 sets adjusted alike), Holm over H224.1-3, attenuation
                  (D - D_adj) / D with a leave-one-cell-out jackknife, descriptives (naive mutant vs wild-type, 2 x 2 mutation x 17p loss).
  --smoke --out DIR : both stages on random noise expression (real labels, cells and copy numbers) in DIR, to test the code path.
    python -I compute_a224.py --stage power|tests [--nperm 2000] [--npseudo 200] [--jobs 8]
Dosage: CN_g = GISTIC thresholded state (-2..2); per cohort (TCGA: the 23 projects of A215 = A219 TCGA projects + BRCA; METABRIC), per gene g, on natural-log expression Y
(A219 / A215 definition), OLS  Y_g ~ project + beta_g CN_g + cos(phi) + sin(phi), phi = per-tumour phase of the leave-gene-out projection of the other 13 clock genes
(project-centred) on the GTEx waveform plane (complex weights amp x exp(-1j psi)). Adjusted Y_g = Y_g - beta_g CN_g (CN missing: no adjustment). Pseudo genes: own CN, own beta,
phase of the 14 real clock genes (nothing left out). Cells, labels, permutations, seeds and pseudo sets are those of A219 (H219.2) and A215 (H215.1)."""
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
EXP = HERE.parents[1]
sys.path.insert(0, str(EXP / "A168_endocrine" / "scripts"))
sys.path.insert(0, str(ROOT / "paper"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402
from scipy.stats import f as fdist, norm, t as tdist  # noqa: E402

import a168_lib as L  # noqa: E402
import clocksync as cs  # noqa: E402
from common import data  # noqa: E402
from common.paths import DATA  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--stage", required=True, choices=["power", "tests"])
ap.add_argument("--nperm", type=int, default=2000); ap.add_argument("--npseudo", type=int, default=200)
ap.add_argument("--jobs", type=int, default=8)
ap.add_argument("--smoke", action="store_true"); ap.add_argument("--out", default=None)
args = ap.parse_args()
assert args.jobs <= 8, "budget: at most 8 workers"
if args.smoke:
    assert args.out, "--smoke needs --out DIR"; args.nperm, args.npseudo = min(args.nperm, 40), min(args.npseudo, 3)
OUT = Path(args.out) if args.out else HERE.parent / "output"; OUT.mkdir(parents=True, exist_ok=True)
A215O, A217, A219O = EXP / "A215_checkpoint_coherence" / "output", EXP / "A217_subtype_matched" / "output", EXP / "A219_subtype_coherence" / "output"
CACHE = EXP / "A210_myc_coherence" / "output" / "tcga_cache"
XENA = DATA / "tcga" / "xena_gistic" / "Gistic2_CopyNumber_Gistic2_all_thresholded.by_genes.gz"
MBD = DATA / "metabric" / "brca_metabric"
K, MIN_SUB, MIN_ARM, MDE_COH, MDE_LEVEL = norm.ppf(0.95) + norm.ppf(0.80), 15, 10, 0.25, 0.3
SEED = zlib.crc32(b"A224") % (2 ** 31)                       # new draws (H224.3 labels, jackknife is deterministic)
SEED219, SEED215 = zlib.crc32(b"A219") % (2 ** 31), zlib.crc32(b"A215") % (2 ** 31)   # reused: same permutations / pseudo sets as A219 / A215
NCHUNK = 20
GENES, LN2 = list(L.GENES), L.LN2
CW = L.WF.amp * np.exp(-1j * L.WF.psi)
ALIAS = {"CIART": ["CIART", "C1orf51"]}
pg = lambda nl, o: cs.stats.perm_p(nl, o, "greater", tie=0.0)
pl = lambda nl, o: cs.stats.perm_p(nl, o, "less", tie=0.0)
say = lambda *a: print(pd.Timestamp.now().strftime("%F %T"), *a, flush=True)
tertile = lambda x: L.A166.tertile(x) if len(x) >= 3 else np.ones(len(x), int)
rng_noise = np.random.default_rng(0)
noise = lambda n, m=len(GENES): rng_noise.normal(size=(n, m))

# ---- cohorts ----------------------------------------------------------------------------------------------------------------------
pj = pd.read_csv(A217 / "projects.csv"); TP = list(pj.project[pj.in_tests]); TCGA219 = TP + ["BRCA"]
PROJ_A97 = ["BLCA", "BRCA", "COAD", "KIRC", "KIRP", "LIHC", "LUAD", "OV", "PRAD", "THCA", "UCEC"]
PROJ_A172 = ["CESC", "ESCA", "GBM", "HNSC", "LGG", "LUSC", "PAAD", "PCPG", "READ", "SARC", "STAD", "TGCT"]
PROJECTS = PROJ_A97 + PROJ_A172
assert set(TCGA219) <= set(PROJECTS)
sub_tab = pd.read_csv(A217 / "subtypes_tcga.csv", keep_default_na=False)
st_tc = pd.read_csv(A215O / "status_tcga.csv", index_col=0, keep_default_na=False, na_values=[""])
tc, Y0 = {}, {}   # project -> table (patient, sample, depth, prolif, sub, tp53, tp53_mut, tp53_gistic), Y (n x 14, natural log)
for p in PROJECTS:
    d = pd.read_parquet(CACHE / f"{p}.parquet").reset_index(drop=True)
    s = sub_tab[sub_tab.project == p].set_index("patient").subtype.replace("", np.nan)
    d["sub"] = d.patient.map(s)
    if p == "BRCA":
        d["sub"] = d["sub"].str.replace("BRCA_", "", regex=False)
    for c in ("tp53", "tp53_mut", "tp53_gistic"):
        d[c] = st_tc[c].reindex(d["sample"]).to_numpy()
    tc[p] = d
    Y0[p] = noise(len(d)) if args.smoke else np.log((d[GENES].to_numpy(float) + 0.5) / d.depth.to_numpy(float)[:, None])
mb, Xmb = L.A166.load("metabric"); mb = mb.reset_index(drop=True); Xmb = Xmb.reset_index(drop=True)
mb["sub"] = mb.pam50.where(mb.pam50.notna())
Ymb0 = noise(len(mb)) if args.smoke else LN2 * Xmb[GENES].to_numpy(float)


# ---- copy number ------------------------------------------------------------------------------------------------------------------
def read_xena(genes, samples):
    """GISTIC thresholded state (genes x samples, NaN where absent) of `genes` (by Xena name) for the 15-char `samples`; text parse only."""
    want, out = set(genes), {}
    with gzip.open(XENA, "rt") as f:
        r = csv.reader(f, delimiter="\t"); head = next(r); pos = {s: i + 1 for i, s in enumerate(head[1:])}
        idx = [pos.get(s, -1) for s in samples]
        for row in r:
            if row[0] in want and row[0] not in out:
                v = np.array([float(row[i]) if i > 0 and row[i] != "" else np.nan for i in idx]); out[row[0]] = v
    return pd.DataFrame(out, index=list(samples)).T


def cn_tcga(genes):
    """genes x all TCGA samples of PROJECTS; clock-gene aliases resolved (CIART = C1orf51); genes with no row are dropped."""
    names = {g: next((a for a in ALIAS.get(g, [g])), g) for g in genes}
    ext = {g: ALIAS.get(g, [g]) for g in genes}
    samples = [s for p in PROJECTS for s in tc[p]["sample"]]
    raw = read_xena(sorted({a for v in ext.values() for a in v}), samples)
    rows = {}
    for g in genes:
        a = next((x for x in ext[g] if x in raw.index), None)
        if a is not None:
            rows[g] = raw.loc[a]
    return pd.DataFrame(rows).T


def cn_metabric(genes):
    c = pd.read_csv(MBD / "data_cna.txt", sep="\t").drop_duplicates("Hugo_Symbol").set_index("Hugo_Symbol").drop(columns="Entrez_Gene_Id")
    rows = {}
    for g in genes:
        a = next((x for x in ALIAS.get(g, [g]) if x in c.index), None)
        if a is not None:
            rows[g] = c.loc[a].reindex(mb.patient).to_numpy(float)
    return pd.DataFrame(rows, index=mb.patient).T


# ---- phase and dosage fit ---------------------------------------------------------------------------------------------------------
def centre(Y, grp):
    """Subtract the group mean (columns of Y; grp integer codes)."""
    Y = np.asarray(Y, float); n = np.bincount(grp)
    if Y.ndim == 1:
        return Y - (np.bincount(grp, Y) / n)[grp]
    return Y - np.vstack([np.bincount(grp, Y[:, j]) / n for j in range(Y.shape[1])]).T[grp]


def phases(Yc, leave_out=None):
    """Phase (rad) of the projection of centred Y (n x 14) on the complex waveform weights; `leave_out` = gene index left out."""
    w = CW.copy()
    if leave_out is not None:
        w[leave_out] = 0
    return np.angle(Yc @ w)


def ols(y, X, grp):
    """Within-group OLS of y on X (group dummies absorbed). Returns coef, se, dof, residual SS."""
    yc, Xc = centre(y, grp), centre(X, grp)
    beta, *_ = np.linalg.lstsq(Xc, yc, rcond=None); res = yc - Xc @ beta
    dof = len(y) - len(np.unique(grp)) - X.shape[1]; s2 = (res ** 2).sum() / dof
    se = np.sqrt(s2 * np.diag(np.linalg.pinv(Xc.T @ Xc)))
    return beta, se, dof, float((res ** 2).sum())


def fit_gene(y, cn, ph, grp, interaction=False):
    """Fit Y ~ group + beta CN + cos + sin (rows with finite y, CN, phase). Returns dict, or None when CN has no within-group variance."""
    ok = np.isfinite(y) & np.isfinite(cn) & np.isfinite(ph); y, cn, ph, g = y[ok], cn[ok], ph[ok], np.unique(grp[ok], return_inverse=True)[1]
    if ok.sum() < 30 or centre(cn, g).std() == 0:
        return None
    X0 = np.column_stack([np.cos(ph), np.sin(ph)]); X1 = np.column_stack([cn, X0])
    b, se, dof, ss1 = ols(y, X1, g); _, _, _, ss0 = ols(y, X0, g); _, _, _, ssc = ols(y, cn[:, None], g)
    yc = centre(y, g)
    r = dict(n=int(ok.sum()), n_groups=int(g.max() + 1), beta_ln=b[0], se_ln=se[0], beta=b[0] / LN2, se=se[0] / LN2, phase_cos=b[1] / LN2, phase_sin=b[2] / LN2,
             cn_var_share=float(centre(b[0] * cn, g).var() / yc.var()), cn_partial_r2=float(1 - ss1 / ss0), cn_alone_r2=float(1 - ssc / (yc ** 2).sum()),
             p_beta=float(2 * tdist.sf(abs(b[0] / se[0]), dof)), cn_sd=float(centre(cn, g).std()))
    if interaction:
        X2 = np.column_stack([X1, cn * X0[:, 0], cn * X0[:, 1]]); b2, se2, dof2, ss2 = ols(y, X2, g)
        F = ((ss1 - ss2) / 2) / (ss2 / dof2)
        r.update(beta_inter=b2[0] / LN2, int_cos=b2[3] / LN2, int_cos_se=se2[3] / LN2, int_sin=b2[4] / LN2, int_sin_se=se2[4] / LN2, int_joint_p=float(fdist.sf(F, 2, dof2)),
                 int_amp=float(np.hypot(b2[3], b2[4]) / LN2))
    return r


# ---- designs (A219 H219.2 and A215 H215.1; copied, Y replaced) -------------------------------------------------------------------------
def cell_value(Y, wf, c, lab):
    return L.coh(Y, c["idx"][lab == 1], wf) - L.coh(Y, c["idx"][lab == 0], wf)


class CellD:
    """Labelled arms within cells (label 1 minus label 0); labels permuted within cell."""
    def __init__(s, cells): s.cells = cells
    def draw(s, rng): return [rng.permutation(len(c["lab"])) for c in s.cells]
    def values(s, Y, wf, st=None):
        return np.array([cell_value(Y, wf, c, c["lab"] if st is None else c["lab"][st[k]]) for k, c in enumerate(s.cells)])


def make_cells_219(tab):
    """A219 H219.2 cells of one cohort: proliferation tertile over tumours with a subtype; LumA vs the pooled entering others (LumB / Her2 / Basal), >= 15 per arm-subtype."""
    t = tab[tab["sub"].notna()].copy(); t["pos"] = t.index
    t["g"] = t.groupby(lambda i: 0).prolif.transform(tertile)
    rows, cells = [], []
    for g, c in t.groupby("g"):
        n = c["sub"].value_counts(); enter = [s for s in n.index if n[s] >= MIN_SUB and s in ("LumA", "LumB", "Her2", "Basal")]
        ok = len(enter) >= 2 and "LumA" in enter
        for s in sorted(n.index):
            rows.append(dict(cell=int(g), subtype=s, n=int(n[s]), enters=s in enter, cell_eligible=ok))
        if ok:
            k = c[c["sub"].isin(enter)]
            cells.append(dict(key=f"{int(g)}", idx=k.pos.to_numpy(), lab=(k["sub"] == "LumA").astype(int).to_numpy(), sub=k["sub"].to_numpy()))
    return pd.DataFrame(rows), cells


def cells_status(m, col):
    """A215: stratum x proliferation-tertile cells of the tumours of m with status `col`; arms altered (1) / intact (0), >= MIN_ARM each."""
    m = m[m[col].notna()].copy(); m["t"] = m.groupby("stratum").prolif.transform(tertile); rows, cells = [], []
    for (s, t), g in m.groupby(["stratum", "t"]):
        n1, n0 = int((g[col] == 1).sum()), int((g[col] == 0).sum()); ok = min(n1, n0) >= MIN_ARM
        rows.append(dict(stratum=s, cell=int(t), n=len(g), n_arm1=n1, n_arm0=n0, eligible=ok))
        if ok:
            cells.append(dict(key=f"{s}|{int(t)}", stratum=s, idx=g.index.to_numpy(), lab=g[col].to_numpy().astype(int)))
    return pd.DataFrame(rows), cells


mb_tab = mb[["sub", "prolif"]].copy()
e_mb, c_mb = make_cells_219(mb_tab); e_br, c_br = make_cells_219(tc["BRCA"][["sub", "prolif"]])
cells_p = {}; el215 = []
for p in PROJECTS:
    e, c = cells_status(tc[p].assign(stratum=p), "tp53"); e.insert(0, "project", p); el215.append(e); cells_p[p] = c
pd.concat([e_mb.assign(cohort="METABRIC", hypothesis="H224.1"), e_br.assign(cohort="BRCA", hypothesis="H224.1")] + [e.assign(hypothesis="H224.2") for e in el215],
          ignore_index=True).to_csv(OUT / "eligibility.csv", index=False)

# ---- H224.3 eligibility: TP53 non-silent mutant vs wild type, no PER1 / TP53 loss and no PER2 change -------------------------------
IPER1, IPER2 = GENES.index("PER1"), GENES.index("PER2")
need = sorted(set(GENES))
CNt = None


def load_cn():
    global CNt, CNm
    CNt = cn_tcga(need); CNm = cn_metabric(need)
    say("CN genes found: TCGA", len(CNt), "METABRIC", len(CNm), "missing TCGA:", sorted(set(GENES) - set(CNt.index)), "METABRIC:", sorted(set(GENES) - set(CNm.index)))


load_cn()
cn_of = lambda CN, g, ids: CN.loc[g, ids].to_numpy(float) if g in CN.index else np.full(len(ids), np.nan)


def dosage_clock():
    """Fit the 14 clock genes in TCGA (pooled over PROJECTS), TCGA BRCA alone (descriptive) and METABRIC.
    Returns dosage rows, adjusted Y dicts, per-tumour leave-gene-out phases and the pooled phase coefficients."""
    rows, adj, phs, coefs = [], {}, {}, {}
    sets = {"TCGA": (PROJECTS, [Y0[p] for p in PROJECTS], [tc[p]["sample"].to_numpy() for p in PROJECTS], CNt),
            "METABRIC": (["METABRIC"], [Ymb0], [mb.patient.to_numpy()], CNm)}
    for coh, (names, Ys, ids, CN) in sets.items():
        grp = np.concatenate([np.full(len(y), i) for i, y in enumerate(Ys)])
        Yall = np.vstack(Ys); Yc = np.vstack([y - y.mean(0) for y in Ys]); idsall = np.concatenate(ids)
        ph = np.column_stack([phases(Yc, j) for j in range(len(GENES))]); phs[coh] = ph
        for j, g in enumerate(GENES):
            cn = cn_of(CN, g, idsall)
            for tag, m in ((coh, np.ones(len(cn), bool)),) + ((("TCGA_BRCA", grp == names.index("BRCA")),) if coh == "TCGA" else ()):
                r = fit_gene(Yall[m, j], cn[m], ph[m, j], grp[m], interaction=True)
                if r is None:
                    rows.append(dict(gene=g, cohort=tag, n=int(np.isfinite(cn[m]).sum()), note="no CN call or no CN variance: not adjusted"))
                else:
                    rows.append(dict(gene=g, cohort=tag, n_cn_nonzero=int((np.nan_to_num(cn[m]) != 0).sum()), note="", **r))
                if tag == coh:
                    coefs[(coh, g)] = r; adj.setdefault(coh, np.empty_like(Yall))[:, j] = Yall[:, j] - (0 if r is None else r["beta_ln"]) * np.nan_to_num(cn)
    return pd.DataFrame(rows), adj, phs, coefs


dos, ADJ, PH, COEF = dosage_clock()
dos.to_csv(OUT / "dosage.csv", index=False); say("dosage fitted")
print(dos[dos.cohort != "TCGA_BRCA"][["gene", "cohort", "n", "beta", "se", "cn_var_share", "cn_partial_r2", "int_joint_p"]].round(4).to_string(), flush=True)
ofs = np.cumsum([0] + [len(tc[p]) for p in PROJECTS])
Yt = {p: ADJ["TCGA"][ofs[i]:ofs[i + 1]] for i, p in enumerate(PROJECTS)}
Ymb = ADJ["METABRIC"]
PHt = {p: PH["TCGA"][ofs[i]:ofs[i + 1]] for i, p in enumerate(PROJECTS)}

# H224.3 data: phase-adjusted, dosage-adjusted log2 PER1 and PER2 (mean of the two) among tumours with GISTIC PER1 >= 0, TP53 >= 0, PER2 == 0 and a TP53 mutation call
h3, h3rows = {}, []
for i, p in enumerate(PROJECTS):
    d = tc[p]; c1, c2 = cn_of(CNt, "PER1", d["sample"].to_numpy()), cn_of(CNt, "PER2", d["sample"].to_numpy())
    lvl = []
    for j in (IPER1, IPER2):
        r = COEF[("TCGA", GENES[j])]; ph = PHt[p][:, j]
        lvl.append((Yt[p][:, j] - r["phase_cos"] * LN2 * np.cos(ph) - r["phase_sin"] * LN2 * np.sin(ph)) / LN2)
    y = np.mean(lvl, axis=0)
    ok = (d.tp53_mut.notna() & (d.tp53_gistic >= 0) & (c1 >= 0) & (c2 == 0)).to_numpy()
    n1, n0 = int((ok & (d.tp53_mut == 1)).sum()), int((ok & (d.tp53_mut == 0)).sum())
    h3rows.append(dict(project=p, n_tp53_mut=n1, n_tp53_wt=n0, eligible=min(n1, n0) >= MIN_ARM))
    if min(n1, n0) >= MIN_ARM:
        h3[p] = dict(y=y[ok], lab=d.tp53_mut.to_numpy()[ok].astype(int))
pd.DataFrame(h3rows).to_csv(OUT / "eligibility_h224_3.csv", index=False)


def fe_diff(y, lab):
    """Mutant minus wild-type difference with project fixed effects (within-project OLS slope on the label), pooled over projects."""
    num = den = 0.0
    for yy, ll in zip(y, lab):
        lc = ll - ll.mean(); num += (lc * (yy - yy.mean())).sum(); den += (lc ** 2).sum()
    return num / den


def h3_chunk(seed, n):
    rng = np.random.default_rng(seed); ys = [h3[p]["y"] for p in h3]; ls = [h3[p]["lab"] for p in h3]
    return np.array([fe_diff(ys, [rng.permutation(l) for l in ls]) for _ in range(n)])


# ---- blocks -----------------------------------------------------------------------------------------------------------------------
ncells = lambda cc: int(sum(len(c["lab"]) for c in cc))
blocks = {"H224.1": [(Ymb, L.WF, CellD(c_mb))] * bool(c_mb) + [(Yt["BRCA"], L.WF, CellD(c_br))] * bool(c_br),
          "H224.2": [(Yt[p], L.WF, CellD(cells_p[p])) for p in PROJECTS if cells_p[p]]}
nunits = lambda h: int(sum(len(d.cells) for _, _, d in blocks[h]))
ntum = lambda h: int(sum(ncells(d.cells) for _, _, d in blocks[h]))
pd.DataFrame([dict(hypothesis=h, n_cohorts=len(blocks[h]), n_cells=nunits(h), n_tumours=ntum(h)) for h in blocks]
             + [dict(hypothesis="H224.3", n_cohorts=len(h3), n_cells=len(h3), n_tumours=sum(len(v["y"]) for v in h3.values()))]).to_csv(OUT / "counts.csv", index=False)
say({h: (len(blocks[h]), nunits(h), ntum(h)) for h in blocks}, "H224.3 projects", len(h3))


def null_chunk(bl, seed, n):
    rng = np.random.default_rng(seed); out = np.empty(n)
    for i in range(n):
        out[i] = np.mean(np.concatenate([d.values(Y, wf, d.draw(rng)) for Y, wf, d in bl]))
    return out


def draws(bl, seed, tag, n):
    sizes = [n // NCHUNK + (j < n % NCHUNK) for j in range(NCHUNK)]
    return np.concatenate(Parallel(n_jobs=args.jobs)(delayed(null_chunk)(bl, np.random.SeedSequence([seed, tag, j]), s) for j, s in enumerate(sizes)))


HYPS = ["H224.1", "H224.2", "H224.3"]
MEAS = {"H224.1": "D_adj: LumA - others at matched proliferation, METABRIC + TCGA BRCA (A219 H219.2 on adjusted expression)",
        "H224.2": "T_adj: TCGA TP53-altered - intact at matched proliferation (A215 H215.1 on adjusted expression)",
        "H224.3": "TP53 mutant - wild-type, phase- and dosage-adjusted mean log2 PER1 / PER2, no PER1 / TP53 loss, PER2 CN 0 (log2 units)"}
if args.stage == "power":
    nl = {}
    # same seeds / tags as the original power stages: A219 H219.2 (tag 1), A215 H215.1 (tag 0)
    nl["H224.1"] = draws(blocks["H224.1"], SEED219, 1, args.nperm); say("H224.1 null", nl["H224.1"].mean().round(5), nl["H224.1"].std(ddof=1).round(5))
    nl["H224.2"] = draws(blocks["H224.2"], SEED215, 0, args.nperm); say("H224.2 null", nl["H224.2"].mean().round(5), nl["H224.2"].std(ddof=1).round(5))
    sizes = [args.nperm // NCHUNK + (j < args.nperm % NCHUNK) for j in range(NCHUNK)]
    nl["H224.3"] = np.concatenate([h3_chunk(np.random.SeedSequence([SEED, 3, j]), s) for j, s in enumerate(sizes)]) if h3 else np.array([np.nan])
    say("H224.3 null", np.mean(nl["H224.3"]).round(5), np.std(nl["H224.3"], ddof=1).round(5))
    pd.DataFrame({h: pd.Series(v) for h, v in nl.items()}).to_csv(OUT / "null_draws.csv", index=False)
    pw = pd.DataFrame([dict(hypothesis=h, test=MEAS[h], n_cohorts=len(blocks[h]) if h in blocks else len(h3), n_units=nunits(h) if h in blocks else len(h3),
                            n_tumours=ntum(h) if h in blocks else sum(len(v["y"]) for v in h3.values()), n_perm=args.nperm, null_mean=np.mean(nl[h]), null_sd=np.std(nl[h], ddof=1),
                            MDE_80pct=K * np.std(nl[h], ddof=1), limit=MDE_LEVEL if h == "H224.3" else MDE_COH) for h in HYPS])
    pw["scored"] = pw.MDE_80pct <= pw.limit
    pw.to_csv(OUT / "power.csv", index=False); print(pw.round(4).to_string(), flush=True)
    say("done power"); sys.exit(0)

# ---- tests (NOT run before the dated power note) ---------------------------------------------------------------------------------------
pw = pd.read_csv(OUT / "power.csv").set_index("hypothesis"); nl = pd.read_csv(OUT / "null_draws.csv"); assert len(nl) == args.nperm
Y_raw = {"METABRIC": Ymb0, **Y0}
raw_blocks = {"H224.1": [(Ymb0, L.WF, CellD(c_mb))] * bool(c_mb) + [(Y0["BRCA"], L.WF, CellD(c_br))] * bool(c_br),
              "H224.2": [(Y0[p], L.WF, CellD(cells_p[p])) for p in PROJECTS if cells_p[p]]}
unit = lambda bl: np.concatenate([d.values(Y, wf) for Y, wf, d in bl])
u_adj = {h: unit(blocks[h]) for h in blocks}; u_raw = {h: unit(raw_blocks[h]) for h in blocks}
obs = {h: float(u_adj[h].mean()) for h in blocks}; obs_raw = {h: float(u_raw[h].mean()) for h in blocks}
obs["H224.3"] = fe_diff([h3[p]["y"] for p in h3], [h3[p]["lab"] for p in h3]) if h3 else np.nan


def jackknife(a, b):
    """Attenuation (mean(a) - mean(b)) / mean(a) with a leave-one-cell-out jackknife (a: unadjusted, b: adjusted per-cell values)."""
    n = len(a); f = lambda i: (a[i].mean() - b[i].mean()) / a[i].mean()
    full = f(np.arange(n)); loo = np.array([f(np.delete(np.arange(n), k)) for k in range(n)])
    return full, float(np.sqrt((n - 1) / n * ((loo - loo.mean()) ** 2).sum()))


# pseudo-clock companions: genes with their own CN and beta, adjusted alike
def beta_pseudo(genes, Ylist_of, CNtab, ids_list, phase_list):
    """Pooled beta (ln) per pseudo gene: Ylist_of(g) -> list of per-group vectors (None where the gene is absent)."""
    out = {}
    for g in genes:
        if g not in CNtab.index:
            continue
        ys, cns, phs_, grp = [], [], [], []
        for k, (y, ids, ph) in enumerate(zip(Ylist_of(g), ids_list, phase_list)):
            if y is None:
                continue
            ys.append(y); cns.append(CNtab.loc[g, ids].to_numpy(float)); phs_.append(ph); grp.append(np.full(len(y), k))
        if not ys:
            continue
        r = fit_gene(np.concatenate(ys), np.concatenate(cns), np.concatenate(phs_), np.concatenate(grp))
        if r is not None:
            out[g] = r["beta_ln"]
    return out


def pseudo_values(des, Y, wf):
    return des.values(Y, wf)


def par(it):
    return np.array(Parallel(n_jobs=args.jobs)(delayed(pseudo_values)(d, Y, wf) for d, Y, wf in it))


ph14_t = {p: phases(Y0[p] - Y0[p].mean(0)) for p in PROJECTS}; ph14_m = phases(Ymb0 - Ymb0.mean(0))   # phase of the 14 real clock genes (nothing left out)
sets_mb = pd.read_csv(L.A166.OUT / "pseudo_sets_metabric.csv").to_numpy()[: args.npseudo].tolist()
pseudo, beta_rows = {}, []
B = None if args.smoke else data._builder("tcga_paired")


def primaries(B, p):
    s = B.samples_of(p)
    return s[s.sample_type == "Primary Tumor"].sort_values("barcode").drop_duplicates("patient")


def load_lc(p):
    """log2 CPM of genes with mean CPM >= 1 for the primaries of p (A213 `load_lc`), rows = patients."""
    pr = primaries(B, p); ct = B.counts_of(p, pr); depth = ct.sum(axis=1)
    lc = np.log2(ct.div(depth, axis=0) * 1e6 + 1)
    return lc.loc[:, (2 ** lc - 1).mean() >= 1].astype(np.float32).set_axis(pr.patient.to_numpy())


# METABRIC pseudo sets (stored A166 sets): own CN, own pooled beta
if args.smoke:
    allm = GENES; CNmp = CNm; bmb = {g: (COEF[("METABRIC", g)] or {"beta_ln": 0.0})["beta_ln"] for g in GENES}
    sets_mb = [GENES] * args.npseudo; Xm = lambda s: noise(len(mb))
else:
    allm = sorted({g for s in sets_mb for g in s}); CNmp = cn_metabric(allm)
    bmb = beta_pseudo(allm, lambda g: [LN2 * Xmb[g].to_numpy(float)], CNmp, [mb.patient.to_numpy()], [ph14_m]); Xm = lambda s: LN2 * Xmb[s].to_numpy(float)
beta_rows += [dict(cohort="METABRIC", gene=g, beta_ln=b) for g, b in bmb.items()]


def adj_with(Y, CNtab, ids, s, bet):
    cn = np.column_stack([np.nan_to_num(CNtab.loc[g, ids].to_numpy(float)) if g in CNtab.index else np.zeros(len(ids)) for g in s])
    return Y - np.array([bet.get(g, 0.0) for g in s]) * cn


mb_units = par((CellD(c_mb), adj_with(Xm(s), CNmp, mb.patient.to_numpy(), s, bmb), L.wf_pseudo(L.WF, s)) for s in sets_mb)
say("METABRIC pseudo done")
# TCGA: sets per project (A219 seed for BRCA / H224.1; A215 seeds for all projects / H224.2); beta pooled over all projects
if args.smoke:
    sets219 = [GENES] * args.npseudo; sets215 = {p: [GENES] * args.npseudo for p in PROJECTS}; CNtp = CNt
    btc = {g: (COEF[("TCGA", g)] or {"beta_ln": 0.0})["beta_ln"] for g in GENES}
    Ytc_of = lambda p, s: noise(len(tc[p]))
else:
    LC, sets215 = {}, {}
    for j, p in enumerate(PROJECTS):
        LC[p] = load_lc(p); assert (LC[p].index.to_numpy() == tc[p].patient.to_numpy()).all()
        sets215[p] = L.pseudo_sets(LC[p], L.A166.pool_genes(LC[p]), GENES, args.npseudo, np.random.default_rng(np.random.SeedSequence([SEED215, 1000 + j])))
    sets219 = L.pseudo_sets(LC["BRCA"], L.A166.pool_genes(LC["BRCA"]), GENES, args.npseudo, np.random.default_rng(np.random.SeedSequence([SEED219, 1000 + TCGA219.index("BRCA")])))
    say("TCGA pseudo sets generated")
    allg = sorted({g for s in sets219 for g in s} | {g for p in PROJECTS for s in sets215[p] for g in s}); CNtp = cn_tcga(allg)
    dep = {p: tc[p].depth.to_numpy(float) for p in PROJECTS}
    Yg = lambda g: [np.log((2 ** LC[p][g].to_numpy(float) - 1) / 1e6 + 0.5 / dep[p]) if g in LC[p].columns else None for p in PROJECTS]
    btc = beta_pseudo(allg, Yg, CNtp, [tc[p]["sample"].to_numpy() for p in PROJECTS], [ph14_t[p] for p in PROJECTS])
    Ytc_of = lambda p, s: np.log((2 ** LC[p][list(s)].to_numpy(float) - 1) / 1e6 + 0.5 / dep[p][:, None])
beta_rows += [dict(cohort="TCGA", gene=g, beta_ln=b) for g, b in btc.items()]
pd.DataFrame(beta_rows).to_csv(OUT / "pseudo_beta.csv", index=False); say("pseudo betas fitted")
tcga_adj = lambda p, s: adj_with(Ytc_of(p, s), CNtp, tc[p]["sample"].to_numpy(), s, btc)
br_units = par((CellD(c_br), tcga_adj("BRCA", s), L.wf_pseudo(L.WF, s)) for s in sets219) if c_br else np.empty((args.npseudo, 0))
pseudo["H224.1"] = np.concatenate([mb_units, br_units], axis=1).mean(axis=1); say("H224.1 pseudo done")
pc2 = []
for p in PROJECTS:
    if cells_p[p]:
        pc2.append(par((CellD(cells_p[p]), tcga_adj(p, s), L.wf_pseudo(L.WF, s)) for s in sets215[p]))
    say(p, "H224.2 pseudo done")
pseudo["H224.2"] = np.concatenate(pc2, axis=1).mean(axis=1)
pd.DataFrame(pseudo).to_csv(OUT / "pseudo_stat.csv", index=False)

rows = []
for h in HYPS:
    sc = bool(pw.loc[h, "scored"])
    r = dict(hypothesis=h, measure=pw.loc[h, "test"], scored=sc, stat=obs[h], null_mean=pw.loc[h, "null_mean"], null_sd=pw.loc[h, "null_sd"], MDE_80pct=pw.loc[h, "MDE_80pct"],
             excess=obs[h] - pw.loc[h, "null_mean"])
    if h == "H224.1":
        r["p"] = pg(nl[h].to_numpy(), obs[h]) if sc else np.nan; q = float(np.quantile(pseudo[h], 0.95)); r.update(pseudo_percentile=100 * np.mean(pseudo[h] <= obs[h]), pseudo_q95=q, specific=bool(obs[h] > q))
    elif h == "H224.2":
        r["p"] = pl(nl[h].to_numpy(), obs[h]) if sc else np.nan; q = float(np.quantile(pseudo[h], 0.05)); r.update(pseudo_percentile=100 * np.mean(pseudo[h] <= obs[h]), pseudo_q05=q, specific=bool(obs[h] < q))
    else:
        d = nl[h].to_numpy(); r["p"] = (1 + np.sum(np.abs(d - pw.loc[h, "null_mean"]) >= abs(obs[h] - pw.loc[h, "null_mean"]))) / (len(d) + 1) if sc else np.nan
    if h in blocks:
        r["stat_unadjusted"] = obs_raw[h]
        a, s_ = jackknife(u_raw[h], u_adj[h]); r.update(attenuation=a, attenuation_jackknife_se=s_, attenuation_lo=a - 1.96 * s_, attenuation_hi=a + 1.96 * s_)
    rows.append(r)
res = pd.DataFrame(rows)
s = res.scored.to_numpy(); ps = res.p[s].to_numpy(); o = np.argsort(ps)
adj = np.empty(len(ps)); adj[o] = np.minimum(1, np.maximum.accumulate(np.sort(ps) * (len(ps) - np.arange(len(ps))))); res["p_holm"] = np.nan; res.loc[s, "p_holm"] = adj
res["verdict"] = np.where(~res.scored, "not scored", np.where(res.p_holm < 0.05, "pass", "fail"))
res.to_csv(OUT / "tests.csv", index=False); print(res.round(4).to_string(), flush=True)

# descriptives: naive mutant vs wild-type over all tumours, 2 x 2 (mutation x 17p loss) of the phase-adjusted level (project-centred)
dsc, grid = [], []
ys_all, labs_all, loss_all, grp_all = [], [], [], []
for i, p in enumerate(PROJECTS):
    d = tc[p]; lv = []
    for j in (IPER1, IPER2):
        r = COEF[("TCGA", GENES[j])]; ph = PHt[p][:, j]
        lv.append((Y0[p][:, j] - r["phase_cos"] * LN2 * np.cos(ph) - r["phase_sin"] * LN2 * np.sin(ph)) / LN2)   # unadjusted for dosage (naive)
    ok = d.tp53_mut.notna().to_numpy()
    y = np.mean(lv, axis=0)[ok]; y = y - y.mean()
    ys_all.append(y); labs_all.append(d.tp53_mut.to_numpy()[ok].astype(int)); loss_all.append((d.tp53_gistic.to_numpy()[ok] < 0).astype(int)); grp_all.append(np.full(ok.sum(), i))
y, lab, loss = map(np.concatenate, (ys_all, labs_all, loss_all))
dsc.append(dict(contrast="naive mutant - wild type, all tumours (phase-adjusted, no dosage adjustment, project FE)", value=fe_diff(ys_all, labs_all), n=len(y)))
for m_ in (0, 1):
    for l_ in (0, 1):
        k = (lab == m_) & (loss == l_)
        grid.append(dict(tp53_mutant=m_, tp53_gistic_loss=l_, n=int(k.sum()), mean_log2_per1_per2=float(y[k].mean()) if k.any() else np.nan))
pd.DataFrame(dsc).to_csv(OUT / "descriptive_naive.csv", index=False); pd.DataFrame(grid).to_csv(OUT / "descriptive_2x2.csv", index=False)
print(pd.DataFrame(dsc).to_string(), pd.DataFrame(grid).to_string(), flush=True)
say("done tests")
