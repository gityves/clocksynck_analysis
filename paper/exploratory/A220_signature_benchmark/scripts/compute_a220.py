"""A220: are published-style circadian risk signatures better than size- and expression-matched random gene signatures?
(registered in paper/docs/registered_experiments.md, A220).
  --stage power : TCGA-CDR OS and project eligibility (>= 60 OS events), the expressed GOBP_CIRCADIAN_RHYTHM list per cohort, 200 random lists per
                  cohort (expression-decile matched, outside the circadian list and the 14 clock genes), the pipeline on every random list (10 stratified
                  70/30 splits: univariate Cox p < 0.05 filter, L1 Cox with 5-fold CV, held-out Harrell C, incremental C over age + stage + proliferation
                  [+ subtype]), the null of the statistic P (each random list ranked among the others) and MDE = 2.49 x null SD
                  (eligibility.csv, power.csv, null_percentiles.csv, cache/rand_<cohort>.parquet). The clock list is NOT run.
  --stage tests : the same pipeline on the GO circadian list and on the 14 waveform genes, percentiles among the cached random lists, P for H220.1-3,
                  permutation-style p against the random-list null of P, Holm over the scored hypotheses (power.csv `scored`), descriptive tables
                  (tests.csv, per_cohort.csv, selected_genes.csv).
  --smoke --out DIR : both stages on N(0,1) noise expression (real survival, covariates, splits) with few lists and splits, in DIR.
    python -I compute_a220.py --stage power|tests [--nlists 200] [--nsplit 10] [--jobs 8] [--cohorts BRCA,LUAD,...] [--smoke --out DIR]
Cohorts: TCGA recount3 projects (one primary tumour per patient as A210, OS from TCGA-CDR), METABRIC (cBioPortal OS), SCAN-B GSE96058 (OS)."""
import argparse
import os
import sys
import time
import zlib
from pathlib import Path

for v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(v, "1")
HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(HERE.parents[1] / "A168_endocrine" / "scripts")); sys.path.insert(0, str(ROOT / "paper"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402
from scipy.stats import norm, spearmanr  # noqa: E402

import a168_lib as L  # noqa: E402
import a220_lib as K  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common import data  # noqa: E402
from common.paths import DATA  # noqa: E402

ap = argparse.ArgumentParser()
ap.add_argument("--stage", required=True, choices=["power", "tests"])
ap.add_argument("--nlists", type=int, default=200); ap.add_argument("--nsplit", type=int, default=10)
ap.add_argument("--jobs", type=int, default=8); ap.add_argument("--cohorts", default=None, help="comma list; default all eligible")
ap.add_argument("--smoke", action="store_true"); ap.add_argument("--out", default=None)
args = ap.parse_args()
limit_threads(1)
if args.smoke:
    assert args.out, "--smoke needs --out DIR"; args.nlists, args.nsplit = min(args.nlists, 6), min(args.nsplit, 2)
args.jobs = min(args.jobs, 8)
OUT = Path(args.out) if args.out else HERE.parent / "output"; OUT.mkdir(parents=True, exist_ok=True)
CACHE = OUT / "cache"; CACHE.mkdir(exist_ok=True)
EXP = HERE.parents[1]
SEED = zlib.crc32(b"A220") % (2 ** 31)
MIN_EVENTS, MDE_MAX, KMDE, TEST_FRAC, NFOLD, MIN_COVER = 60, 0.20, norm.ppf(0.95) + norm.ppf(0.80), 0.30, 5, 0.5
GENES = list(L.GENES)
say = lambda *a: print(pd.Timestamp.now().strftime("%F %T"), *a, flush=True)
GMT = DATA / "reference" / "c5.go.bp.v2024.1.Hs.symbols.gmt"
CDR = DATA / "tcga" / "cdr" / "TCGA-CDR-SupplementalTableS1.xlsx"
SUBFILE = EXP / "A217_subtype_matched" / "output" / "subtypes_tcga.csv"
PROJECTS = sorted(p.name for p in (DATA / "tcga" / "recount3").iterdir() if p.is_dir() and p.name.isupper() and p.name != "LAML")
BREAST = ["METABRIC", "SCANB"]


# ---- clinical tables ----------------------------------------------------------------------------------------------------------------
def cdr_table():
    """TCGA-CDR (public Table S1) reduced to patient, OS, OS.time (days), age, stage (1-4 or NaN); cached in output/cache."""
    f = CACHE / "cdr_os.csv"
    if not f.exists():
        d = pd.read_excel(CDR, sheet_name="TCGA-CDR")
        st = d.ajcc_pathologic_tumor_stage.astype(str).str.replace("Stage ", "", regex=False).str.strip()
        stage = st.map(lambda s: 4 if s.startswith("IV") else 3 if s.startswith("III") else 2 if s.startswith("II") else 1 if s in ("I", "IA", "IB", "IC") else np.nan)
        pd.DataFrame(dict(patient=d.bcr_patient_barcode, os=d.OS, os_t=d["OS.time"], age=pd.to_numeric(d.age_at_initial_pathologic_diagnosis, errors="coerce"),
                          stage=stage)).to_csv(f, index=False)
    return pd.read_csv(f).set_index("patient")


def go_list():
    for line in GMT.read_text().splitlines():
        r = line.split("\t")
        if r[0] == "GOBP_CIRCADIAN_RHYTHM":
            return r[2:]
    raise KeyError("GOBP_CIRCADIAN_RHYTHM")


def primaries(B, p):
    s = B.samples_of(p)
    return s[s.sample_type == "Primary Tumor"].sort_values("barcode").drop_duplicates("patient")


def tcga_clinical(B, p, cdr):
    """Primary tumours of p with valid CDR OS (patient index): os_t, os_e, age, stage."""
    t = primaries(B, p).set_index("patient"); c = cdr.reindex(t.index)
    c = c[c.os.notna() & (c.os_t > 0)]
    return t.loc[c.index], c


def onehot(s, min_n=10):
    """Dummies of a categorical series (missing = own level; levels with < min_n pooled into 'other'; first level dropped); empty frame if < 2 levels."""
    if s.dropna().nunique() < 2:
        return pd.DataFrame(index=s.index)
    s = s.fillna("NA").astype(str); vc = s.value_counts(); s = s.where(s.map(vc) >= min_n, "other")
    if s.nunique() < 2:
        return pd.DataFrame(index=s.index)
    return pd.get_dummies(s, drop_first=True, dtype=float).add_prefix("sub_")


def covariates(df, num, cats=(), pre=()):
    """Covariate matrix (rows of df): numeric columns median-imputed (a missing indicator is added when 0 < missing < 100 %; a column observed in
    < 50 % of tumours is dropped), categorical columns one-hot. Constant columns dropped."""
    out = {}
    for c in num:
        x = pd.to_numeric(df[c], errors="coerce"); miss = x.isna().mean()
        if miss > 1 - MIN_COVER:
            continue
        if miss > 0:
            out[c + "_missing"] = x.isna().astype(float)
        out[c] = x.fillna(x.median())
    Z = pd.DataFrame(out, index=df.index)
    for c in cats:
        Z = Z.join(onehot(df[c]))
    return Z.loc[:, Z.std() > 0]


class Cohort:
    """One analysis cohort sorted by time: X (n x G float32 expression), genes, medians, t, d, covariates Z, proliferation z."""
    def __init__(self, name, X, genes, t, d, Z, prolif, source):
        o = np.argsort(t, kind="stable")
        self.name, self.genes, self.source = name, np.array(genes), source
        self.X, self.t, self.d, self.Z, self.prolif = X[o], np.asarray(t, float)[o], np.asarray(d, float)[o], Z.to_numpy(float)[o], np.asarray(prolif, float)[o]
        self.zcols = list(Z.columns); self.index = {g: i for i, g in enumerate(self.genes)}
        self.med = np.median(self.X, axis=0)
        self.expressed = np.where(self.med > np.percentile(self.med, 25))[0]

    @property
    def n(self):
        return len(self.t)

    @property
    def events(self):
        return int(self.d.sum())


def zmean(lc, genes):
    g = [x for x in genes if x in lc.columns]
    z = ((lc[g] - lc[g].mean()) / lc[g].std(ddof=0).replace(0, np.nan)).mean(axis=1)
    return ((z - z.mean()) / z.std()).to_numpy()


def noise(co_args, name):
    """--smoke: N(0,1) expression of the same shape (a fixed seed per cohort)."""
    X, genes, *rest = co_args
    rng = np.random.default_rng([SEED, zlib.crc32(name.encode())])
    return (rng.normal(size=X.shape).astype(np.float32), genes, *rest)


def load_tcga(p, B, cdr, sub):
    t, c = tcga_clinical(B, p, cdr)
    pm = primaries(B, p); ct = B.counts_of(p, pm); ct.index = pm.patient.to_numpy(); depth = ct.sum(axis=1)
    lc = np.log2(ct.div(depth, axis=0) * 1e6 + 1); lc = lc.loc[:, (2 ** lc - 1).mean() >= 1]
    pr = pd.Series(zmean(lc, L.sh.hallmark_g2m()), index=ct.index)  # proliferation metagene of A210 (all primaries, then restricted)
    lc = lc.loc[c.index].astype(np.float32); pr = pr.loc[c.index]
    df = pd.DataFrame(dict(age=c.age, stage=c.stage, prolif=pr), index=c.index)
    s = sub[sub.project == p].set_index("patient").subtype.replace("", np.nan) if len(sub) else pd.Series(dtype=object)
    df["subtype"] = s.reindex(df.index)
    cats = ["subtype"] if df.subtype.notna().mean() >= MIN_COVER else []
    Z = covariates(df, ["age", "stage", "prolif"], cats)
    a = (lc.to_numpy(), list(lc.columns), c.os_t.to_numpy(), c.os.to_numpy(), Z, df.prolif.to_numpy(), "TCGA recount3 + TCGA-CDR")
    if args.smoke:
        a = noise(a, p)
    return Cohort(p, *a)


def load_metabric():
    m, X = L.A166.load("metabric"); m = m[m.pam50.notna()].copy(); X = X.loc[m.index]
    d = DATA / "metabric" / "brca_metabric"
    pt = pd.read_csv(d / "data_clinical_patient.txt", sep="\t", comment="#").set_index("PATIENT_ID")
    smp = pd.read_csv(d / "data_clinical_sample.txt", sep="\t", comment="#").set_index("SAMPLE_ID")
    c = m.join(pt[["AGE_AT_DIAGNOSIS", "LYMPH_NODES_EXAMINED_POSITIVE", "OS_MONTHS", "OS_STATUS"]], on="patient").join(smp[["TUMOR_SIZE"]])
    c["nodes"] = np.log1p(pd.to_numeric(c.LYMPH_NODES_EXAMINED_POSITIVE, errors="coerce"))
    c["os_t"], c["os_e"] = pd.to_numeric(c.OS_MONTHS, errors="coerce"), c.OS_STATUS.astype(str).str.startswith("1").astype(float)
    ok = (c.os_t > 0) & c.OS_STATUS.notna(); c, X = c[ok], X.loc[ok.to_numpy()]
    Z = covariates(c.rename(columns={"AGE_AT_DIAGNOSIS": "age", "TUMOR_SIZE": "size"}), ["age", "size", "nodes", "prolif"], ["pam50"])
    a = (X.to_numpy(np.float32), list(X.columns), c.os_t.to_numpy(), c.os_e.to_numpy(), Z, c.prolif.to_numpy(), "METABRIC cBioPortal")
    return Cohort("METABRIC", *(noise(a, "METABRIC") if args.smoke else a))


def load_scanb():
    m, X = L.A166.A.load()
    c = m.copy(); c["os_t"], c["os_e"] = pd.to_numeric(c["overall survival days"]), pd.to_numeric(c["overall survival event"])
    c["age"], c["size"] = pd.to_numeric(c["age at diagnosis"], errors="coerce"), pd.to_numeric(c["tumor size"], errors="coerce")
    c["node"] = (c["lymph node status"] == "NodePositive").astype(float).where(c["lymph node status"].isin(["NodePositive", "NodeNegative"]))
    ok = (c.os_t > 0) & c.os_e.notna(); c, X = c[ok], X.loc[ok.to_numpy()]
    Z = covariates(c, ["age", "size", "node", "prolif"], ["pam50"])
    a = (X.to_numpy(np.float32), list(X.columns), c.os_t.to_numpy(), c.os_e.to_numpy(), Z, c.prolif.to_numpy(), "SCAN-B GSE96058")
    return Cohort("SCANB", *(noise(a, "SCANB") if args.smoke else a))


# ---- eligibility ----------------------------------------------------------------------------------------------------------------------
def eligibility():
    """eligibility.csv rows for every TCGA project (from CDR + sample table only) and the two breast cohorts (filled once loaded)."""
    cdr, B = cdr_table(), data._builder("tcga_paired"); rows = []
    for p in PROJECTS:
        t, c = tcga_clinical(B, p, cdr)
        rows.append(dict(cohort=p, source="TCGA", n_tumours_expression=len(primaries(B, p)), n_os=len(t), events=int(c.os.sum()),
                         eligible=bool(c.os.sum() >= MIN_EVENTS)))
    return pd.DataFrame(rows)


def get_cohorts(names):
    cdr, B = cdr_table(), data._builder("tcga_paired")
    sub = pd.read_csv(SUBFILE, keep_default_na=False) if SUBFILE.exists() else pd.DataFrame()
    for nm in names:
        t0 = time.time()
        co = load_metabric() if nm == "METABRIC" else load_scanb() if nm == "SCANB" else load_tcga(nm, B, cdr, sub)
        say(f"loaded {nm}: n={co.n} events={co.events} genes={len(co.genes)} covariates={co.zcols} ({time.time() - t0:.0f}s)")
        yield co


# ---- gene lists -----------------------------------------------------------------------------------------------------------------------
def circadian_idx(co):
    """Indices of the GO list genes expressed in the cohort (and of all GO genes present, for the exclusion from the random pool)."""
    go = [g for g in go_list() if g in co.index]
    return np.array([co.index[g] for g in go if co.index[g] in set(co.expressed)]), np.array([co.index[g] for g in go])


def random_lists(co, n):
    """n lists of len(circ) expressed genes, one per circadian gene from its median-expression decile (deciles of the expressed genes), outside the
    GO list and the 14 clock genes, no gene twice within a list. Seeded by crc32('A220') and the cohort."""
    circ, go_all = circadian_idx(co)
    ex = co.expressed; dec = pd.qcut(pd.Series(co.med[ex]).rank(method="first"), 10, labels=False).to_numpy()
    dec_of = dict(zip(ex, dec)); excl = set(go_all) | {co.index[g] for g in GENES if g in co.index}
    pool = {d: np.array([g for g in ex[dec == d] if g not in excl]) for d in range(10)}
    need = np.bincount([dec_of[g] for g in circ], minlength=10)
    assert all(len(pool[d]) >= need[d] for d in range(10)), "decile pool too small"
    rng = np.random.default_rng([SEED, zlib.crc32(co.name.encode())])
    return circ, [np.concatenate([rng.choice(pool[d], need[d], replace=False) for d in range(10)]) for _ in range(n)]


# ---- pipeline -------------------------------------------------------------------------------------------------------------------------
def make_splits(co, nsplit):
    """nsplit stratified 70/30 splits (same for every list of the cohort): (train idx, test idx, CV fold of train rows)."""
    out = []
    for s in range(nsplit):
        rng = np.random.default_rng([SEED, zlib.crc32(co.name.encode()), 1000 + s]); te = np.zeros(co.n, bool); fold = np.zeros(co.n, int)
        for e in (0, 1):
            idx = rng.permutation(np.where(co.d == e)[0]); k = int(round(TEST_FRAC * len(idx))); te[idx[:k]] = True
            fold[idx[k:]] = rng.permutation(np.arange(len(idx) - k) % NFOLD)
        tr = np.where(~te)[0]; out.append((tr, np.where(te)[0], fold[tr]))
    return out


def one_split(X, t, d, Z, tr, te, folds):
    """The published-style pipeline on one split: returns (C, incremental C, log-rank p of the held-out median split, selected columns, beta)."""
    keep, (Xtr, Xte) = K.standardise(X[tr].astype(float), X[te].astype(float))
    ttr, dtr, tte, dte = t[tr], d[tr], t[te], d[te]; f, l = K.ties(ttr)
    p = K.univariate_p(Xtr, dtr, f, l); sel = np.where(p < 0.05)[0]
    if len(sel) < 5:
        sel = np.argsort(p)[:min(5, len(p))]
    sel = np.sort(sel); beta, _, _ = K.lasso_cv(Xtr[:, sel], ttr, dtr, folds)
    s_tr, s_te = Xtr[:, sel] @ beta, Xte[:, sel] @ beta
    c = K.harrell_c(tte, dte, s_te)
    pl = K.logrank_p(tte, dte, s_te > np.median(s_te)) if np.ptp(s_te) > 0 else 1.0
    inc = np.nan
    if Z.shape[1]:
        _, (Ztr, Zte) = K.standardise(Z[tr], Z[te]); bz, _, _ = K.cox_newton(Ztr, dtr, f, l)
        base = K.harrell_c(tte, dte, Zte @ bz)
        if np.ptp(s_tr) > 0:
            sd = s_tr.std(); J, Jte = np.column_stack([Ztr, s_tr / sd]), np.column_stack([Zte, s_te / sd])
            bj, _, _ = K.cox_newton(J, dtr, f, l); inc = K.harrell_c(tte, dte, Jte @ bj) - base
        else:
            inc = 0.0
    return c, inc, pl, sel, keep, beta


def run_chunk(X, t, d, Z, splits, lists, ids):
    """Pipeline for the gene-index lists `lists` (ids = their labels): DataFrame list, split, c, inc, p_lr, n_pass_filter, n_selected."""
    rows = []
    for i, cols in zip(ids, lists):
        Xl = X[:, cols]
        for s, (tr, te, folds) in enumerate(splits):
            c, inc, pl, sel, keep, beta = one_split(Xl, t, d, Z, tr, te, folds)
            rows.append((i, s, c, inc, pl, len(sel), int((beta != 0).sum())))
    return pd.DataFrame(rows, columns=["list", "split", "c", "inc", "p_lr", "n_filter", "n_selected"])


def run_lists(co, lists, ids, splits, tag):
    """Chunks of lists in parallel; logs each completed chunk (timestamped)."""
    per = max(1, int(np.ceil(len(lists) / (4 * args.jobs)))); chunks = [(lists[i:i + per], ids[i:i + per]) for i in range(0, len(lists), per)]
    t0, parts = time.time(), []
    gen = Parallel(n_jobs=args.jobs, return_as="generator")(delayed(run_chunk)(co.X, co.t, co.d, co.Z, splits, a, b) for a, b in chunks)
    for k, df in enumerate(gen):
        parts.append(df); say(f"  {tag}: chunk {k + 1}/{len(chunks)} done ({time.time() - t0:.0f}s)")
    return pd.concat(parts, ignore_index=True)


def summarise(df):
    """Mean over splits per list: c, inc; validation rate; median-split validation (median log-rank p < 0.05)."""
    g = df.groupby("list")
    return pd.DataFrame(dict(c=g.c.mean(), inc=g.inc.mean(), p_lr_median=g.p_lr.median(), valid_rate=g.p_lr.apply(lambda x: (x < 0.05).mean()),
                             n_filter=g.n_filter.mean(), n_selected=g.n_selected.mean()))


def percentile(x, ref):
    """(# ref < x + 0.5 # ref == x) / len(ref)."""
    ref = np.asarray(ref); return ((ref < x).sum() + 0.5 * (ref == x).sum()) / len(ref)


def null_percentiles(vals):
    """Percentile of each list among the others."""
    v = np.asarray(vals); return np.array([percentile(v[i], np.delete(v, i)) for i in range(len(v))])


# ---- stages ---------------------------------------------------------------------------------------------------------------------------
def stage_power():
    el = eligibility(); say("TCGA projects:", el.cohort[el.eligible].tolist())
    rows, nulls = [], {}
    for co in get_cohorts(args.cohorts.split(",") if args.cohorts else list(el[el.eligible].cohort) + BREAST):
        circ, lists = random_lists(co, args.nlists); splits = make_splits(co, args.nsplit)
        say(f"{co.name}: circadian list {len(circ)} expressed genes; {len(lists)} random lists x {len(splits)} splits")
        f = CACHE / f"rand_{co.name}.parquet"
        if f.exists() and (len(pd.read_parquet(f)) == len(lists) * len(splits)):
            df = pd.read_parquet(f); say("  cached")
        else:
            df = run_lists(co, lists, list(range(len(lists))), splits, co.name); df.to_parquet(f)
        S = summarise(df); nulls[co.name] = dict(c=null_percentiles(S.c), inc=null_percentiles(S.inc), n=len(S))
        rows.append(dict(cohort=co.name, source=co.source, n=co.n, events=co.events, n_expressed_genes=len(co.expressed), n_genes=len(co.genes),
                         circ_list_size=len(circ), n_covariates=len(co.zcols), covariates=";".join(co.zcols), n_random_lists=len(S), n_splits=len(splits),
                         random_c_mean=S.c.mean(), random_c_sd=S.c.std(), random_inc_mean=S.inc.mean(), random_inc_sd=S.inc.std(),
                         random_valid_share=(S.p_lr_median < 0.05).mean()))
    E = pd.DataFrame(rows); isb = E.cohort.isin(BREAST)
    E = pd.concat([el.merge(E.drop(columns=["source", "n", "events"]), on="cohort", how="left"),
                   E[isb].rename(columns={"n": "n_os"}).assign(eligible=True)], ignore_index=True)
    E.to_csv(OUT / "eligibility.csv", index=False)
    tc = [c for c in nulls if c not in BREAST]; pr = []
    nl = min(v["n"] for v in nulls.values()); N = lambda key, cs: np.mean([nulls[c][key][:nl] for c in cs], axis=0)
    for h, key, cs in (("H220.1", "c", tc), ("H220.2", "inc", tc), ("H220.3", "inc", [c for c in BREAST if c in nulls])):
        if not cs:
            continue
        P = N(key, cs); sd = P.std(ddof=1)
        pr.append(dict(hypothesis=h, n_cohorts=len(cs), cohorts=";".join(cs), n_null_lists=nl, null_mean=P.mean(), null_sd=sd, mde=KMDE * sd, scored=bool(KMDE * sd <= MDE_MAX)))
    pd.DataFrame(pr).to_csv(OUT / "power.csv", index=False)
    pd.DataFrame({f"{c}_{k}": nulls[c][k][:nl] for c in nulls for k in ("c", "inc")}).to_csv(OUT / "null_percentiles.csv", index=False)
    say("power stage done"); print(pd.DataFrame(pr).to_string())


def holm(p):
    p = np.asarray(p, float); o = np.argsort(p); m = len(p); adj = np.empty(m); run = 0.0
    for r, i in enumerate(o):
        run = max(run, (m - r) * p[i]); adj[i] = min(1.0, run)
    return adj


def stage_tests():
    el = pd.read_csv(OUT / "eligibility.csv"); pw = pd.read_csv(OUT / "power.csv").set_index("hypothesis")
    cohorts = list(el.cohort[el.eligible & el.circ_list_size.notna()])
    pct, desc, selrows = {}, [], []
    for co in get_cohorts(cohorts):
        circ, lists = random_lists(co, args.nlists); splits = make_splits(co, args.nsplit)
        S = summarise(pd.read_parquet(CACHE / f"rand_{co.name}.parquet"))
        clk14 = np.array([co.index[g] for g in GENES if g in co.index])
        res = run_lists(co, [circ, clk14], ["go", "clock14"], splits, co.name + " clock lists")
        R = summarise(res); go, c14 = R.loc["go"], R.loc["clock14"]
        pct[co.name] = dict(c=percentile(go.c, S.c), inc=percentile(go.inc, S.inc))
        desc.append(dict(cohort=co.name, n=co.n, events=co.events, circ_list_size=len(circ), clock_c=go.c, clock_inc=go.inc, pct_c=pct[co.name]["c"],
                         pct_inc=pct[co.name]["inc"], clock_validated=go.p_lr_median < 0.05, clock_valid_rate=go.valid_rate,
                         random_valid_share=(S.p_lr_median < 0.05).mean(), random_c_mean=S.c.mean(), clock14_c=c14.c, clock14_inc=c14.inc,
                         clock14_pct_c=percentile(c14.c, S.c), clock14_pct_inc=percentile(c14.inc, S.inc)))
        # selected genes: frequency of non-zero LASSO coefficient across splits, correlation of the gene with proliferation
        cnt = np.zeros(len(circ))
        for tr, te, folds in splits:
            keep, (Xtr, Xte) = K.standardise(co.X[tr][:, circ].astype(float), co.X[te][:, circ].astype(float)); ttr, dtr = co.t[tr], co.d[tr]
            f, l = K.ties(ttr); p = K.univariate_p(Xtr, dtr, f, l); sel = np.sort(np.where(p < 0.05)[0] if (p < 0.05).sum() >= 5 else np.argsort(p)[:5])
            b, _, _ = K.lasso_cv(Xtr[:, sel], ttr, dtr, folds); cnt[np.where(keep)[0][sel[b != 0]]] += 1
        for j, g in enumerate(circ):
            if cnt[j]:
                selrows.append(dict(cohort=co.name, gene=co.genes[g], share_of_splits=cnt[j] / len(splits), rho_proliferation=spearmanr(co.X[:, g], co.prolif)[0]))
    D = pd.DataFrame(desc); D.to_csv(OUT / "per_cohort.csv", index=False); pd.DataFrame(selrows).to_csv(OUT / "selected_genes.csv", index=False)
    nulls = pd.read_csv(OUT / "null_percentiles.csv")
    tc = [c for c in D.cohort if c not in BREAST]; rows = []
    for h, key, cs in (("H220.1", "c", tc), ("H220.2", "inc", tc), ("H220.3", "inc", [c for c in BREAST if c in set(D.cohort)])):
        if not cs or h not in pw.index:
            continue
        Pobs = np.mean([pct[c][key] for c in cs]); nullP = np.mean([nulls[f"{c}_{key}"].to_numpy() for c in cs], axis=0)
        rows.append(dict(hypothesis=h, n_cohorts=len(cs), P=Pobs, null_mean=nullP.mean(), null_sd=nullP.std(ddof=1), mde=pw.loc[h, "mde"], scored=bool(pw.loc[h, "scored"]),
                         p_one_sided=(1 + (nullP >= Pobs).sum()) / (1 + len(nullP))))
    T = pd.DataFrame(rows); sc = T.scored.to_numpy()
    T["p_holm"] = np.nan; T.loc[sc, "p_holm"] = holm(T.p_one_sided[sc])
    T.to_csv(OUT / "tests.csv", index=False); say("tests stage done"); print(T.to_string())


if __name__ == "__main__":
    {"power": stage_power, "tests": stage_tests}[args.stage]()
