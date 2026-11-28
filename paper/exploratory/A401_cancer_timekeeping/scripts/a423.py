"""A423: is the per-cancer cancer-cell coupling a contamination of loosely labelled cancer cells by normal epithelium?
Registered in paper/docs/registered_experiments.md (A423). Data: the A401 pooled tumour table (`a401.pool_sources`, its
cancer labels and deduplication); host = the cancer's pooled tumour-stroma anchors, fit as A401 (`Dataset.fit`, cohort x
chemistry strata, 14 genes deciding, core loop beside).

Within-patient split. Cells = every tumour-sample epithelial cell of a patient (marker rule of `cancer_common.compartments`,
the class-E definition in every source). Split score, the 14 clock genes removed from every profile and score:
  - THCA: mean within-patient z of log CP10k of LGALS3 KRT19 FN1 CITED1 SERPINA1 TIMP1 MET minus that of TFF3 TPO DIO1 SLC26A7
    (clonality is reported beside it);
  - else copy-number clonality (`cancer_common.cnv`): the stored calls `$SCCC_DATA/cnv{,_a149}/<study>.csv.gz` where they
    exist, otherwise computed against the baseline of the cohort's own non-epithelial tumour cells (stromal and immune;
    classic inferCNV; ranking only).
Cells are split at the within-patient median into HIGH (more aberrant) and LOW; a patient enters with >= 50 cells per half.
Pseudobulks (genome-wide counts, canonical symbols) of HIGH, LOW, ALL (all epithelial cells) and, where stored calls give
>= 50 cells of each, the called-normal / called-malignant cells are cached in `output/a423/cache/split_<source>.parquet`.

  power  (blind: counts, labels, split validation, simulation; NO clock fit of any split compartment)
         patients.csv, counts.csv (cancer x source x split score), label_classes.csv, validation_patients.csv,
         validation.csv (D423.2: share of called-normal cells in LOW vs HIGH), unsplit.csv, power_h1.csv, power_h2.csv,
         power_summary.csv
  score  H423.1 (Delta = c_LOW - c_HIGH, patient-weighted over cancers with >= 5 split patients; label-swap null within
         patient, 1,000; 200 matched sets per cancer; matched q95), H423.2 (Stouffer z of c_HIGH vs matched over LUAD, THCA,
         KIRC), Holm; D423.1 (A401 target by label class), D423.2 (called-normal vs called-malignant), D423.3 (per cancer
         c_ALL / c_HIGH / c_LOW, z, median UMI per half, H401.2 with e_mal from c_HIGH)
         NOT TO BE RUN before the power note is written.
  --smoke  code check on a few small cohorts: the split targets' donor labels are permuted within stratum and nothing
         about a fitted split target is printed or written.

    python a423.py --stage power [--jobs 4] [--genome-jobs 2] [--smoke]       (long: paper/common/run.sh)
"""
from __future__ import annotations

import argparse
import gc
import re
import sys
import time
import warnings
from pathlib import Path

HERE = Path(__file__).resolve()
PAPER = HERE.parents[3]
EXP = PAPER / "exploratory"
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parents[1] / "build_geo"))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import scipy.sparse as sp  # noqa: E402
from joblib import Parallel, delayed  # noqa: E402
from scipy.stats import norm, spearmanr  # noqa: E402

import a401 as U  # noqa: E402
import clocksync as cs  # noqa: E402
from cancer_common import cnv  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from clocksync.genes import canonical  # noqa: E402

B, C = U.B, U.C
OUT = U.OUT / "a423"
BASE, CORE, META, MALT, ANCH, READOUTS = U.BASE, U.CORE, U.META, U.MALT, U.ANCH, U.READOUTS
SEED = 423
HIGH, LOW, ALL = "tumour:epi_high", "tumour:epi_low", "tumour:epi_all"
CALLN, CALLM = "tumour:epi_called_normal", "tumour:epi_called_malignant"
MIN_HALF = 50
CLOCK = {canonical(g) for g in BASE.genes}
BASELINE = ("fibroblast", "endothelial", "mural", "T_NK", "B_plasma", "myeloid", "mast")
THCA_UP = ["LGALS3", "KRT19", "FN1", "CITED1", "SERPINA1", "TIMP1", "MET"]
THCA_DOWN = ["TFF3", "TPO", "DIO1", "SLC26A7"]
H2_CANCERS = ("LUAD", "THCA", "KIRC")
AUTHOR_GEO = {"GSE125449", "GSE202051", "GSE212527"}  # author malignant labels (class A)
NOT_SPLIT = {"GSE212527": "sarcoma: no epithelial compartment", "GSE176078": "Wu atlas: bespoke loader (not split)",
             "GSE161529": "Pal atlas: bespoke loader (not split)"}
SMOKE_SPLIT = ("GSE281736", "GSE222315", "SRP416180", "GSE164789", "GSE163203")
C_NORMAL, SHARES, C_HIGH_GRID, C_NULL, Z_GRID = 0.43, (0.2, 0.4), (0.15, 0.25), 0.02, (0.5, 1.0)
N_SIM = 200
DATA = U.DATA


def log(*a):
    print(time.strftime("%F %T"), *a, flush=True)


# ---- cell-level readers -------------------------------------------------------------------------------------------------

def a388_samples(study, patients):
    """Tumour samples (sample, patient, load) of an A388 cohort for `patients` (the A388 unit rules, majority chemistry)."""
    B.K.install()
    us, _ = B.L.ST.units(study)
    return [(s, p, ld) for s, p, k, ld, _c, _n in B.L.ST.majority_chemistry(us) if p in patients and k == "tumour"]


def _loader(read):
    """GEO loader returning the tuple of `cancer_common.studies.units` loaders: counts, genes, UMI (sum of the cell's
    counts, as the pseudobulk builds), marker-rule compartment, barcodes, symbols."""
    from geobuild import marker_comp

    def load():
        X, sym, bc = read()
        X = sp.csr_matrix(X)
        umi = np.asarray(X.sum(1)).ravel().astype(float)
        return X, sym, umi, marker_comp(X, sym, umi), bc, sym
    return load


def _slice(whole, rows):
    """Loader of a subset of the cells of a matrix read once (memoised)."""
    def read():
        X, sym, bc = whole()
        return X[rows], sym, bc[rows]
    return read


def _memo(f):
    box = []

    def g():
        if not box:
            box.append(f())
        return box[0]
    return g


def geo_samples(g, patients):
    """Tumour samples (sample, patient, load) of a GEO build, from the raw files with the logic of `build_geo/build_<g>.py`."""
    from geobuild import read_10x, read_dense, read_dense_stream, read_h5_10x
    raw = U.GEO / g / "raw"
    x = raw / "x"
    out = []
    if g == "GSE207493":
        for d in sorted(p for p in x.iterdir() if p.is_dir()):
            out.append((d.name, d.name, _loader(lambda d=d: read_10x(d))))
    elif g == "GSE163203":
        for f in sorted(x.glob("*_dge.txt.gz")):
            m = re.match(r"GSM\d+_(Adj_)?PTC(with|without)HT_(\d+)(?:_\d)?_dge", f.name)
            if not m.group(1):
                out.append((f"PTC_{m.group(2)}HT_{m.group(3)}", f"P{m.group(3)}", _loader(lambda f=f: read_dense(f))))
    elif g == "GSE191288":
        pat = {"T1L": "P1", "T1R": "P1", "T2L": "P2", "T2R": "P2", "T3L": "P3", "T3R": "P3"}
        for f in sorted(x.glob("*.h5")):
            s = f.stem.split("_", 1)[1]
            if s in pat:
                out.append((s, pat[s], _loader(lambda f=f: read_h5_10x(f))))
    elif g == "GSE135337":
        for f in sorted(x.glob("*_gene_cell_exprs_table.*gz")):
            s = re.match(r"GSM\d+_(BC\w+?)_gene", f.name).group(1)
            if s != "BCN":
                out.append((s, s, _loader(lambda f=f: read_dense(f, sym_col="Symbol", skip_cols=2))))
    elif g == "GSE148071":
        for f in sorted(x.glob("*_exp.txt.gz")):
            n = re.search(r"_(P\d+)_", f.name).group(1)
            out.append((n, n, _loader(lambda f=f: read_dense(f))))
    elif g == "GSE155698":
        for d in sorted(p for p in x.iterdir() if p.is_dir()):
            if d.name.startswith("AdjNorm"):
                continue
            pat = "P" + re.sub(r"[AB]$", "", d.name.split("TISSUE_")[1])

            def read(d=d):
                if (d / "filtered_feature_bc_matrix").is_dir():
                    return read_10x(d / "filtered_feature_bc_matrix")
                if (d / "filtered_gene_bc_matrices").is_dir():
                    return read_10x(d / "filtered_gene_bc_matrices", gz=False)
                return read_h5_10x(d / "filtered_feature_bc_matrix.h5")
            out.append((d.name, pat, _loader(read)))
    elif g == "GSE197461":
        for s in ("SCC_1", "SCC_2", "SCC_3", "ADC_1", "ADC_2", "ADC_3", "ADC_4", "ADC_5"):
            pre = next(x.glob(f"GSM*_{s}_barcodes.tsv.gz")).name.replace("barcodes.tsv.gz", "")
            out.append((s, s, _loader(lambda pre=pre: read_10x(x / pre))))
    elif g == "GSE141445":
        whole = _memo(lambda: read_dense_stream(x / "GSM4203181_data.raw.matrix.txt.gz"))
        suf = _memo(lambda: pd.Series(whole()[2]).str.extract(r"-(\d+)$", expand=False).astype(int).to_numpy())
        for i in range(1, 14):
            out.append((f"S{i}", f"S{i}", _loader(lambda i=i: _slice(whole, np.flatnonzero(suf() == i))())))
    elif g == "GSE125449":
        sets = {s: _memo(lambda s=s: read_10x(raw / f"GSE125449_{s}_")) for s in ("Set1", "Set2")}
        for s in sets:
            smp = pd.read_csv(raw / f"GSE125449_{s}_samples.txt.gz", sep="\t")
            for name in smp.Sample.unique():
                out.append((name, name.split("_")[1], _loader(lambda s=s, n=name, smp=smp: _slice(
                    sets[s], np.flatnonzero((smp.Sample == n).to_numpy()))())))
    elif g == "GSE202051":
        def whole():
            import anndata as ad
            import h5py
            f = x / "GSE202051_totaldata-final-toshare.h5ad"
            a = ad.read_h5ad(f, backed="r")
            o = a.obs[["sampleid", "pid", "treatment_status"]].copy()
            with h5py.File(f, "r") as h:
                X = ad.io.read_elem(h["layers/counts"]).tocsr()
            keep = np.flatnonzero((o.treatment_status == "Untreated").to_numpy())
            return X[keep], np.array(a.var_names), np.array(o.pid.astype(str))[keep]
        whole = _memo(whole)
        for p in sorted(patients):
            out.append((p, p, _loader(lambda p=p: _by_pid(whole, p))))
    else:
        raise KeyError(g)
    return [o for o in out if o[1] in patients]


def _by_pid(whole, p):
    X, sym, pid = whole()
    i = np.flatnonzero(pid == p)
    return X[i], sym, np.array([f"{p}_{k}" for k in range(len(i))])


# ---- cohort cells ---------------------------------------------------------------------------------------------------------

class Cohort:
    """Tumour-sample epithelial cells of a cohort per patient, on one canonical-symbol gene axis; baseline = mean log CP10k of
    the non-epithelial (stromal / immune) tumour cells of the cohort (None when not requested)."""

    def __init__(self):
        self.ids, self.names, self.pat, self.base, self.n_base = {}, [], {}, None, 0

    def gene_ids(self, sym):
        out = []
        for s in sym:
            c = canonical(str(s))
            if c not in self.ids:
                self.ids[c] = len(self.names)
                self.names.append(c)
            out.append(self.ids[c])
        return np.array(out)

    def read(self, samples, baseline):
        bs = np.zeros(0)
        for sample, patient, load in samples:
            X, _g, umi, comp, bcs, sym = load()
            ids = self.gene_ids(sym)
            P = sp.csr_matrix((np.ones(len(ids)), (np.arange(len(ids)), ids)), shape=(len(ids), len(self.names)))
            epi = np.flatnonzero(comp == "epithelial")
            d = self.pat.setdefault(patient, dict(X=[], umi=[], bc=[], sample=[]))
            d["X"].append((sp.csr_matrix(X[epi]) @ P).tocsr().astype(np.float32))
            d["umi"].append(np.asarray(umi, float)[epi])
            d["bc"].append(np.asarray(bcs)[epi].astype(str))
            d["sample"].append(np.full(len(epi), sample))
            if baseline:
                nb = np.flatnonzero(np.isin(comp, BASELINE))
                if len(bs) < len(self.names):
                    bs = np.pad(bs, (0, len(self.names) - len(bs)))
                for k in range(0, len(nb), 2000):
                    j = nb[k:k + 2000]
                    E = cnv.log_norm(sp.csr_matrix(X[j]) @ P, np.asarray(umi, float)[j], np.arange(len(self.names)))
                    bs[:E.shape[1]] += E.sum(0)
                self.n_base += len(nb)
            del X
            gc.collect()
        n = len(self.names)
        for d in self.pat.values():
            for X in d["X"]:
                X.resize((X.shape[0], n))
        if baseline and self.n_base:
            self.base = np.pad(bs, (0, n - len(bs))) / self.n_base
        for p, d in self.pat.items():
            d["X"], d["umi"] = sp.vstack(d["X"]).tocsr(), np.concatenate(d["umi"])
            d["bc"], d["sample"] = np.concatenate(d["bc"]), np.concatenate(d["sample"])
        return self


def stored_calls(study):
    """Stored copy-number calls of a study (clonality, cnv_call) indexed by (sample, barcode), or None."""
    dirs = [B.L.CNV_A149, B.L.ST.CNV] if study in B.L.NEW else [B.L.ST.CNV, B.L.CNV_A149]
    for d in dirs:
        f = d / f"{study}.csv.gz"
        if f.exists():
            t = pd.read_csv(f, usecols=["sample", "barcode", "clonality", "cnv_call"]).set_index(["sample", "barcode"])
            return t[~t.index.duplicated()]
    return None


def clonality_computed(coh):
    """Clonality (`cnv.clonality`: correlation with the patient's own consensus minus the mean correlation with the other
    patients') of every epithelial cell against the cohort's stromal / immune baseline, the 14 clock genes removed; streamed in
    three passes over the cells (scores, consensus of the upper half by score, correlations). {patient: array} or None."""
    if coh.base is None or len(coh.pat) < 2:
        return None
    pos = cnv.positions()
    cols, chrom = cnv.genome_order(np.array(coh.names), pos.assign(symbol=pos.symbol.str.upper()))
    keep = ~np.isin(np.array(coh.names)[cols], list(CLOCK))
    cols, chrom = cols[keep], chrom[keep]
    s = sum(np.asarray(d["X"][:, cols].sum(0)).ravel() for d in coh.pat.values())
    expressed = s / sum(d["X"].shape[0] for d in coh.pat.values()) >= 0.1
    base = coh.base[cols]

    def chunks(d, size=3000):
        for k in range(0, d["X"].shape[0], size):
            sl = slice(k, k + size)
            E = cnv.log_norm(d["X"][sl], d["umi"][sl], cols)
            yield sl, cnv.profiles(E, base, chrom, expressed)

    score = {p: np.concatenate([cnv.score(P) for _, P in chunks(d)]) for p, d in coh.pat.items()}
    cons = {}
    for p, d in coh.pat.items():
        sel, tot, n = score[p] >= np.median(score[p]), 0.0, 0
        for sl, P in chunks(d):
            m = sel[sl]
            tot, n = tot + P[m].sum(0), n + int(m.sum())
        cons[p] = tot / n
    out = {}
    for p, d in coh.pat.items():
        v = []
        for _, P in chunks(d):
            other = [cnv.corr(P, c) for q, c in cons.items() if q != p]
            v.append(cnv.corr(P, cons[p]) - np.mean(other, 0))
        out[p] = np.concatenate(v)
    return out


def marker_score(d, names):
    """THCA papillary-carcinoma score: mean within-patient z of log CP10k of the up genes minus that of the down genes."""
    ix = {g: i for i, g in enumerate(names)}

    def z(genes):
        j = [ix[canonical(g)] for g in genes if canonical(g) in ix]
        if not j:
            return None
        E = np.log1p(1e4 * d["X"][:, j].toarray() / np.maximum(d["umi"][:, None], 1.0))
        sd = E.std(0)
        return ((E - E.mean(0)) / np.where(sd > 0, sd, 1.0)).mean(1)
    up, dn = z(THCA_UP), z(THCA_DOWN)
    return None if up is None or dn is None else up - dn


def halves(score):
    """Boolean HIGH / LOW masks over cells with a finite score: split at the within-patient median (rank, ties in cell order)."""
    ok = np.flatnonzero(np.isfinite(score))
    rk = pd.Series(score[ok]).rank(method="first").to_numpy()
    high = np.zeros(len(score), bool)
    low = np.zeros(len(score), bool)
    high[ok[rk > len(ok) / 2]] = True
    low[ok[rk <= len(ok) / 2]] = True
    return high, low


def split_source(source, patients, cancer_of, stratum_of, cache, smoke=False):
    """Cell-level split of one source: patient table and cached pseudobulks. `patients`: ids within the source."""
    limit_threads(1)
    warnings.filterwarnings("ignore")
    f, ft = cache / f"split_{source}.parquet", cache / f"split_{source}_patients.csv"
    if f.exists() and ft.exists():
        t = pd.read_csv(ft)
        if set(t.patient.astype(str)) == set(patients):
            return source, t, ""
    try:
        is_geo = source in U.GEO_CANCER
        samples = geo_samples(source, set(patients)) if is_geo else a388_samples(source, set(patients))
        calls = stored_calls(source)
        coh = Cohort().read(samples, baseline=calls is None)
        comp = clonality_computed(coh) if calls is None else None
        rows, pats, genes = [], [], np.array(coh.names)
        for p, d in coh.pat.items():
            donor, cancer = f"{source}|{p}", cancer_of.get(f"{source}|{p}")
            n = d["X"].shape[0]
            r = dict(source=source, patient=p, donor=donor, cancer=cancer, n_epi=n, split_score="", reason="",
                     clonality_source="stored_calls" if calls is not None else "computed")
            clon, call = None, None
            if calls is not None:
                m = calls.reindex(pd.MultiIndex.from_arrays([d["sample"], d["bc"]]))
                clon, call = m.clonality.to_numpy(float), m.cnv_call.to_numpy(object)
                r["n_with_calls"] = int(np.isfinite(clon).sum())
            elif comp is not None:
                clon = comp[p]
            if cancer == "THCA":
                score, r["split_score"] = marker_score(d, coh.names), "marker_score"
                if score is None:
                    r["reason"] = "THCA marker genes missing"
                elif clon is not None and np.isfinite(clon).sum() > 10:
                    ok = np.isfinite(clon)
                    r["spearman_marker_clonality"] = float(spearmanr(score[ok], clon[ok])[0])
            else:
                score, r["split_score"] = clon, "stored_calls_clonality" if calls is not None else "computed_clonality"
                if score is None:
                    r["reason"] = "no clonality (single patient or no stromal baseline)"
            if score is None or n < 2 * MIN_HALF:
                r["reason"] = r["reason"] or f"< {2 * MIN_HALF} epithelial cells"
                pats.append(r)
                continue
            hi, lo = halves(np.asarray(score, float))
            r.update(n_high=int(hi.sum()), n_low=int(lo.sum()))
            if min(hi.sum(), lo.sum()) < MIN_HALF:
                r["reason"] = f"< {MIN_HALF} scored cells in a half"
                pats.append(r)
                continue
            groups = {HIGH: hi, LOW: lo, ALL: np.ones(n, bool)}
            if call is not None:
                cn, cm = call == "normal", call == "malignant"
                r.update(n_called_normal=int(cn.sum()), n_called_malignant=int(cm.sum()),
                         normal_in_low=int((cn & lo).sum()), normal_in_high=int((cn & hi).sum()),
                         malignant_in_low=int((cm & lo).sum()), malignant_in_high=int((cm & hi).sum()),
                         n_low_with_call=int((lo & (call != "nan") & pd.notna(call)).sum()),
                         n_high_with_call=int((hi & (call != "nan") & pd.notna(call)).sum()))
                if min(cn.sum(), cm.sum()) >= MIN_HALF:
                    groups.update({CALLN: cn, CALLM: cm})
            r["median_umi_high"], r["median_umi_low"] = float(np.median(d["umi"][hi])), float(np.median(d["umi"][lo]))
            for name, m in groups.items():
                i = np.flatnonzero(m)
                rows.append((dict(donor=donor, compartment=name, stratum=stratum_of[donor], n_cells=len(i),
                                  depth=float(d["umi"][i].sum()), median_umi=float(np.median(d["umi"][i]))),
                             np.asarray(d["X"][i].sum(0)).ravel().astype(np.float32)))
            pats.append(r)
        T = pd.DataFrame(pats)
        if rows:
            G = pd.concat([pd.DataFrame([m for m, _ in rows]), pd.DataFrame(np.vstack([v for _, v in rows]), columns=genes)], axis=1)
            cache.mkdir(parents=True, exist_ok=True)
            G.to_parquet(f)
        else:
            T["no_rows"] = True
        T.to_csv(ft, index=False)
        return source, T, ""
    except Exception as e:  # noqa: BLE001
        return source, None, f"{type(e).__name__}: {e}"[:300]


# ---- label classes, counts, validation ------------------------------------------------------------------------------------

def label_classes(donors, f388):
    """Class (C copy-number calls, A author malignant labels, E all tumour epithelium) of every source of the A401 table."""
    coh = pd.read_csv(B.FEAS.parent / "cohorts.csv").set_index("study").malignant_label
    cen = pd.read_csv(U.CEN / "celltype_mapping.csv").groupby("dataset").tumour_label.first()
    rows = []
    for s in sorted(set(donors.source)):
        if s in f388:
            c = "A" if s in B.BREAST else "C" if coh.get(s) == "tumour:malignant" else "E"
            why = "Wu / Pal author labels" if s in B.BREAST else f"A388 cohort label {coh.get(s)}"
        elif s in U.GEO_CANCER:
            c, why = ("A", "author malignant labels") if s in AUTHOR_GEO else ("E", "marker-rule epithelium")
        elif s == "bassez2021":
            c, why = "A", "BIOKEY 3CA author labels"
        else:
            c, why = ("A", "CELLxGENE malignant cell") if cen.get(s) == "malignant" else ("E", "CELLxGENE epithelial"), None
            why = why or "celltype_mapping.csv tumour_label"
        rows.append(dict(source=s, label_class=c, basis=why))
    return pd.DataFrame(rows)


def power_tables(P, donors, classes, out):
    """counts, unsplit, validation from the patient table (no clock value)."""
    P = P.copy()
    P["split"] = P.n_high.fillna(0).ge(MIN_HALF) & P.n_low.fillna(0).ge(MIN_HALF)
    P = P.merge(classes, on="source", how="left")
    P.to_csv(out / "patients.csv", index=False)
    tab = donors.groupby(["cancer", "source"]).donor.nunique().rename("patients_in_table").reset_index()
    sp_ = P[P.split].groupby(["cancer", "source", "label_class", "split_score"]).patient.nunique().rename("patients_split").reset_index()
    cnt = tab.merge(classes, on="source", how="left").merge(sp_, on=["cancer", "source", "label_class"], how="left").fillna({"patients_split": 0})
    cnt["patients_split"] = cnt.patients_split.astype(int)
    cnt["can_split"] = ~cnt.source.isin(NOT_SPLIT) & cnt.source.isin(set(P.source))
    cnt["cell_level"] = cnt.source.isin(set(P.source)) | cnt.source.isin(NOT_SPLIT)
    cnt.sort_values(["cancer", "source"]).to_csv(out / "counts.csv", index=False)
    per = cnt.groupby("cancer").agg(patients_in_table=("patients_in_table", "sum"), patients_split=("patients_split", "sum"),
                                    sources_split=("patients_split", lambda x: int((x > 0).sum()))).reset_index()
    per["h423_1_eligible"] = per.patients_split >= 5
    per.to_csv(out / "counts_cancer.csv", index=False)
    P[~P.split][["source", "patient", "cancer", "n_epi", "split_score", "reason"]].to_csv(out / "unsplit.csv", index=False)
    v = P[P.split & P.n_called_normal.notna()].copy()
    for k, col in (("low", "normal_in_low"), ("high", "normal_in_high")):
        v[f"share_normal_{k}"] = v[col] / v[f"n_{k}_with_call"].clip(lower=1)
    v["low_minus_high"] = v.share_normal_low - v.share_normal_high
    v[["source", "patient", "cancer", "label_class", "split_score", "n_low", "n_high", "n_called_normal", "n_called_malignant",
       "share_normal_low", "share_normal_high", "low_minus_high"]].to_csv(out / "validation_patients.csv", index=False)
    g = v.groupby(["cancer", "source", "split_score"]).agg(
        patients=("patient", "nunique"), cells_called_normal=("n_called_normal", "sum"), share_normal_low=("share_normal_low", "mean"),
        share_normal_high=("share_normal_high", "mean"), mean_low_minus_high=("low_minus_high", "mean"),
        patients_low_gt_high=("low_minus_high", lambda x: int((x > 0).sum()))).reset_index()
    g.to_csv(out / "validation.csv", index=False)
    return P, cnt, per, g


# ---- tables per cancer -------------------------------------------------------------------------------------------------------

def load_split_rows(cache, sources):
    return {s: pd.read_parquet(cache / f"split_{s}.parquet") for s in sources if (cache / f"split_{s}.parquet").exists()}


def cancer_frame(c, frames, donors, split, P, targets=(HIGH, LOW, ALL), extra=(), mal_key=None):
    """Table of cancer `c`: pooled stromal anchors of all its patients (A401 join over the cohorts' common genes), the split
    rows of its split patients and, with `mal_key` ({donor: class}), the A401 target rows relabelled `tumour:malignant|<class>`
    (classes with >= 5 patients); classes with >= 5 patients; or None."""
    ids = set(donors[donors.cancer == c].donor)
    keep_c = ANCH + ([MALT] if mal_key else [])
    parts = [f[f.donor.isin(ids) & f.compartment.isin(keep_c)] for f in frames.values()]
    parts = [p for p in parts if len(p)]
    allc = pd.concat(parts, join="inner", ignore_index=True)
    classes = [a for a in ANCH if allc[allc.compartment == a].donor.nunique() >= 5]
    ok = set(P[P.split & (P.cancer == c)].donor)
    keep = list(targets) + list(extra)
    sp_ = [s[s.donor.isin(ok) & s.compartment.isin(keep)] for s in split.values()]
    sp_ = [s for s in sp_ if len(s)]
    if len(classes) < 2 or not sp_:
        return None
    genes = [g for g in allc.columns if g not in META]
    for s in sp_:
        genes = [g for g in genes if g in s.columns]
    if not set(BASE.genes) <= set(genes):
        return None
    anc = allc[allc.compartment.isin(classes)]
    mal = allc[allc.compartment == MALT]
    rows = [anc[META + genes]] + [s[META + genes] for s in sp_]
    if mal_key:
        mal = mal.assign(compartment=[f"{MALT}|{mal_key[d]}" for d in mal.donor])
        n = mal.groupby("compartment").donor.nunique()
        rows.append(mal[mal.compartment.isin(n.index[n >= 5])][META + genes])
    G = pd.concat(rows, ignore_index=True)
    G[genes] = G[genes].astype(np.float32)
    return dict(G=G, anc=anc, classes=classes, genes=genes, n_split=len(ok & set(G[G.compartment == LOW].donor)),
                mode="loo" if len(classes) >= 3 else "full")


# ---- power: simulation at the real design ---------------------------------------------------------------------------------

def c_eff(share, c_n=C_NORMAL, ratio=1.0):
    """Coupling of a mixture whose clock RNA is `share` normal (coupled at c_n; loading share * z) and 1 - share cancer cells
    (uncoupled; own clock-plane variance `ratio` times the normal cells' own noise): the clock-plane states mix linearly,
    state = share (z f + e_n) + (1 - share) e_c, so c = share^2 / (share^2 + q (share^2 + ratio (1 - share)^2)),
    q = (1 - c_n) / c_n."""
    q = (1 - c_n) / c_n
    return share ** 2 / (share ** 2 + q * (share ** 2 + ratio * (1 - share) ** 2))


RATIOS = {"equal_total_variance": (1 + (1 - C_NORMAL) / C_NORMAL) / ((1 - C_NORMAL) / C_NORMAL), "equal_own_noise": 1.0}


def sim_task(tab, classes, rd, z, seed, n_sim, tag):
    """Simulated coupling draws of HIGH / LOW at the real design (`simulate.like(..., loading=z)`: the targets are never
    regressed on the host). Returns {scenario: array}."""
    limit_threads(1)
    warnings.filterwarnings("ignore")
    ds = cs.Dataset(tab, anchors=classes, waveform=READOUTS[rd], mixed_strata=True)
    rng = np.random.default_rng(seed)
    out = {}

    def draws(target, c, loading):
        sim = cs.simulate.like(ds, target, loading=loading)
        return np.array([sim.draw(c, rng)["coupling"] for _ in range(n_sim)], float)
    out["null_high"], out["null_low"] = draws(HIGH, C_NULL, z), draws(LOW, C_NULL, z)
    out["null_high2"] = draws(HIGH, C_NULL, z)
    for name, ratio in RATIOS.items():
        for s in SHARES:
            out[f"low_share{s}_{name}"] = draws(LOW, c_eff(s, ratio=ratio), s * z)
    for c in C_HIGH_GRID:
        out[f"high_c{c}"] = draws(HIGH, c, z)
    return out


def stage_sim(tables, args, out):
    n_sim = 10 if args.smoke else N_SIM
    jobs = [(c, rd, z) for c in tables for rd in READOUTS for z in Z_GRID]
    res = Parallel(n_jobs=args.jobs, return_as="generator")(
        delayed(sim_task)(tables[c]["tab"], tables[c]["classes"], rd, z, [SEED, 5, i], n_sim, c) for i, (c, rd, z) in enumerate(jobs))
    R = {j: r for j, r in zip(jobs, res)}
    w = {c: tables[c]["n_split"] for c in tables}
    rows1, rows2 = [], []
    for rd in READOUTS:
        for z in Z_GRID:
            cs_ = list(tables)
            n = np.array([w[c] for c in cs_], float)

            def pooled(low, high):
                L = np.vstack([R[(c, rd, z)][low] for c in cs_])
                H = np.vstack([R[(c, rd, z)][high] for c in cs_])
                return (n[:, None] * (L - H)).sum(0) / n.sum()
            null = pooled("null_low", "null_high")
            thr = float(np.nanquantile(null, .95))
            rows1.append(dict(hypothesis="H423.1", readout=rd, z_abs=z, scenario="null", share_normal=0.0, assumption="", delta_median=float(np.nanmedian(null)),
                              null_q95=thr, power=float(np.nanmean(null > thr)), cancers="+".join(cs_), patients=int(n.sum()), n_sim=n_sim))
            for name in RATIOS:
                for s in SHARES:
                    d = pooled(f"low_share{s}_{name}", "null_high")
                    rows1.append(dict(hypothesis="H423.1", readout=rd, z_abs=z, scenario=f"normal share {s} coupled at {C_NORMAL}",
                                      share_normal=s, assumption=name, c_eff_low=c_eff(s, ratio=RATIOS[name]), delta_median=float(np.nanmedian(d)),
                                      null_q95=thr, power=float(np.nanmean(d > thr)), cancers="+".join(cs_), patients=int(n.sum()), n_sim=n_sim))
            h2 = [c for c in H2_CANCERS if c in tables]
            if not h2:
                continue
            n2 = np.array([w[c] for c in h2], float)
            mu = {c: np.nanmean(R[(c, rd, z)]["null_high"]) for c in h2}
            sd = {c: np.nanstd(R[(c, rd, z)]["null_high"], ddof=1) for c in h2}

            def stouffer(key, weights):
                Z = np.vstack([(R[(c, rd, z)][key] - mu[c]) / sd[c] for c in h2])
                return (weights[:, None] * Z).sum(0) / np.sqrt((weights ** 2).sum())
            for weights, wname in ((n2, "patient-weighted"), (np.ones(len(h2)), "unweighted")):
                t1 = stouffer("null_high2", weights)
                rows2.append(dict(hypothesis="H423.2", readout=rd, z_abs=z, weighting=wname, c_high=C_NULL, power=np.nan,
                                  type_I=float(np.nanmean(t1 > 1.645)), cancers="+".join(h2), patients=int(n2.sum()), n_sim=n_sim))
                for c_high in C_HIGH_GRID:
                    Z = stouffer(f"high_c{c_high}", weights)
                    rows2.append(dict(hypothesis="H423.2", readout=rd, z_abs=z, weighting=wname, c_high=c_high, power=float(np.nanmean(Z > 1.645)),
                                      type_I=np.nan, cancers="+".join(h2), patients=int(n2.sum()), n_sim=n_sim))
    P1, P2 = pd.DataFrame(rows1), pd.DataFrame(rows2)
    P1.to_csv(out / "power_h1.csv", index=False)
    P2.to_csv(out / "power_h2.csv", index=False)
    per = pd.DataFrame([dict(cancer=c, readout=rd, z_abs=z, patients=w[c],
                             cells_low_median=float(tables[c]["tab"][tables[c]["tab"].compartment == LOW].n_cells.median()),
                             cells_high_median=float(tables[c]["tab"][tables[c]["tab"].compartment == HIGH].n_cells.median()),
                             null_sd_high=float(np.nanstd(R[(c, rd, z)]["null_high"], ddof=1)),
                             null_sd_low=float(np.nanstd(R[(c, rd, z)]["null_low"], ddof=1)),
                             null_median_high=float(np.nanmedian(R[(c, rd, z)]["null_high"]))) for (c, rd, z) in R])
    per.to_csv(out / "power_per_cancer.csv", index=False)
    return P1, P2


# ---- power stage ----------------------------------------------------------------------------------------------------------

def slim(frames):
    """Frames reduced to the meta columns and the 14 clock genes (the power stage needs no other gene)."""
    return {s: f[META + list(BASE.genes)].copy() for s, f in frames.items()}


def stage_power(args):
    out = OUT / ("power_smoke" if args.smoke else "power")
    cache = OUT / ("cache_smoke" if args.smoke else "cache")
    out.mkdir(parents=True, exist_ok=True)
    cache.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    pool = cache / "pool_slim.pkl"  # donors, strata and clock-gene frames of the A401 pool: a small parent for the split workers
    if pool.exists():
        frames, donors, f388, stratum_of = pd.read_pickle(pool)
    else:
        frames, donors, f388 = U.pool_sources(args, out)
        from joblib.externals.loky import get_reusable_executor
        get_reusable_executor().shutdown(wait=True)
        stratum_of = {d: s for f in frames.values() for d, s in zip(f.donor, f.stratum)}
        frames = slim(frames)
        gc.collect()
        pd.to_pickle((frames, donors, list(f388), stratum_of), pool)
    donors = donors[donors.cancer.notna()].reset_index(drop=True)
    classes = label_classes(donors, f388)
    classes.to_csv(out / "label_classes.csv", index=False)
    cancer_of = dict(zip(donors.donor, donors.cancer))
    srcs = [s for s in sorted(set(donors.source)) if (s in f388 or s in U.GEO_CANCER) and s not in NOT_SPLIT]
    if args.smoke:
        srcs = [s for s in srcs if s in SMOKE_SPLIT]
    log("split sources:", ", ".join(srcs))
    todo = [(s, [d.split("|", 1)[1] for d in donors[donors.source == s].donor]) for s in srcs]
    res = Parallel(n_jobs=args.jobs, return_as="generator")(
        delayed(split_source)(s, p, cancer_of, stratum_of, cache) for s, p in todo)
    tabs, errs = [], []
    for s, T, e in res:
        if T is None:
            errs.append(dict(source=s, error=e))
            log(f"{s}: FAILED {e}")
        else:
            tabs.append(T)
            log(f"{s}: {len(T)} patients, {int(T.n_high.fillna(0).ge(MIN_HALF).sum()) if 'n_high' in T else 0} split")
    pd.DataFrame(errs, columns=["source", "error"]).to_csv(out / "split_errors.csv", index=False)
    P = pd.concat(tabs, ignore_index=True)
    P["patient"] = P.patient.astype(str)
    P, cnt, per, val = power_tables(P, donors, classes, out)
    pd.set_option("display.width", 250, "display.max_columns", 40, "display.max_rows", 300)
    log("counts per cancer\n" + per.to_string(index=False))
    log("validation\n" + val.round(4).to_string(index=False))
    # simulation at the real design: cancers with >= 5 split patients and a host
    split = load_split_rows(cache, srcs)
    tables = {}
    for c in per[per.h423_1_eligible].cancer:
        t = cancer_frame(c, frames, donors, split, P, targets=(HIGH, LOW))
        if t is None:
            log(f"{c}: no host or no rows: not simulated")
            continue
        t["tab"] = t["G"][META + list(BASE.genes)]
        tables[c] = t
        log(f"{c}: {t['n_split']} split patients, classes {'+'.join(k.split(':')[1] for k in t['classes'])}")
    if args.smoke:
        tables = {c: t for c, t in tables.items()}
    P1, P2 = stage_sim(tables, args, out)
    log("power H423.1\n" + P1.round(4).to_string(index=False))
    log("power H423.2\n" + P2.round(4).to_string(index=False))
    log(f"power done {(time.time() - t0) / 60:.1f} min")


# ---- score ---------------------------------------------------------------------------------------------------------------

def matched_chunk(Gu, classes, targets, ws, rd):
    limit_threads(1)
    warnings.filterwarnings("ignore")
    out = []
    for i, w in ws:
        try:
            t = cs.Dataset(Gu, anchors=classes, waveform=w, mixed_strata=True).fit(n_boot=0, gate=False, targets=targets).table
            out += [dict(set=i, readout=rd, target=k, coupling=float(t.loc[k, "coupling"])) for k in targets]
        except Exception as e:  # noqa: BLE001
            out.append(dict(set=i, readout=rd, error=f"{type(e).__name__}: {e}"[:200]))
    return out


def swap_chunk(tab, classes, js):
    """Couplings of LOW and HIGH with the two labels of each patient swapped at random (draw j: seed [423, 1, j])."""
    limit_threads(1)
    warnings.filterwarnings("ignore")
    don = np.array(sorted(set(tab[tab.compartment == LOW].donor)))
    comp0 = tab.compartment.to_numpy()
    out = []
    for j in js:
        flip = set(don[np.random.default_rng([SEED, 1, j]).random(len(don)) < 0.5])
        m = tab.donor.isin(flip).to_numpy()
        comp = np.where(m & (comp0 == LOW), HIGH, np.where(m & (comp0 == HIGH), LOW, comp0))
        try:
            t = cs.Dataset(tab.assign(compartment=comp), anchors=classes, waveform=BASE, mixed_strata=True).fit(
                n_boot=0, gate=False, targets=[LOW, HIGH]).table
            out.append(dict(draw=j, c_low=float(t.loc[LOW, "coupling"]), c_high=float(t.loc[HIGH, "coupling"])))
        except Exception as e:  # noqa: BLE001
            out.append(dict(draw=j, error=f"{type(e).__name__}: {e}"[:200]))
    return out


def chunks_of(items, n):
    return [items[k::n] for k in range(n) if items[k::n]]


def score_cancer(c, t, args, n_set, n_swap):
    """Observed fits, matched sets and label-swap draws of one cancer. Returns (fits, matched, swaps)."""
    from clocksync.paired import pseudo_waveforms
    G, classes, genes = t["G"], t["classes"], t["genes"]
    tab = G[META + list(BASE.genes)]
    targets = [k for k in sorted(set(G.compartment) - set(classes)) if G[G.compartment == k].donor.nunique() >= 5]
    anc = t["anc"]
    frac = anc[genes].div(anc.depth, axis=0)
    lm = np.log(frac.mean() + 1e-12)[(anc[genes].mean() >= 1) | frac.columns.isin(list(BASE.genes))]
    ws = list(enumerate(pseudo_waveforms(lm, n_set, np.random.default_rng([SEED, 0]), BASE, "matched")))
    used = sorted({g for _, w in ws for g in w.genes} | set(BASE.genes))
    Gu = G[META + used]
    fits = Parallel(n_jobs=min(args.jobs, len(READOUTS)))(
        delayed(U.fit_one)(tab, classes, rd, targets, None, 0, 0, t["mode"]) for rd in READOUTS)
    mt = []
    for rd in READOUTS:
        w_rd = [(i, w if rd == "14" else C.restrict(w, BASE, CORE)) for i, w in ws]
        mt += [(Gu, classes, targets, ch, rd) for ch in chunks_of(w_rd, args.jobs * 2)]
    matched = pd.DataFrame([r for x in Parallel(n_jobs=args.jobs)(delayed(matched_chunk)(*a) for a in mt) for r in x])
    gt = tab[tab.compartment.isin(classes + [HIGH, LOW])]
    swaps = pd.DataFrame([r for x in Parallel(n_jobs=args.jobs)(
        delayed(swap_chunk)(gt, classes, ch) for ch in chunks_of(list(range(n_swap)), args.jobs * 2)) for r in x])
    return dict(zip(READOUTS, fits)), matched, swaps


def holm_rows(p):
    return B.holm(np.array(p, float))


def stage_score(args):
    """Scoring run (one run, after the power note). `--smoke`: split targets' donor labels permuted within stratum; fitted values
    of the targets are neither printed nor written."""
    out = OUT / ("score_smoke" if args.smoke else "score")
    cache = OUT / ("cache_smoke" if args.smoke else "cache")
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    n_set, n_swap = (4, 8) if args.smoke else (200, 1000)
    frames, donors, f388 = U.pool_sources(args, out)
    from joblib.externals.loky import get_reusable_executor
    get_reusable_executor().shutdown(wait=True)
    donors = donors[donors.cancer.notna()].reset_index(drop=True)
    P = pd.read_csv(OUT / ("power_smoke" if args.smoke else "power") / "patients.csv", dtype={"patient": str})
    P["split"] = P.n_high.fillna(0).ge(MIN_HALF) & P.n_low.fillna(0).ge(MIN_HALF)
    split = load_split_rows(cache, sorted(set(P.source)))
    cls = pd.read_csv(OUT / ("power_smoke" if args.smoke else "power") / "label_classes.csv").set_index("source").label_class
    mal_key = {d: cls[s] for d, s in zip(donors.donor, donors.source) if s in cls.index}
    per = P[P.split].groupby("cancer").donor.nunique()
    cancers = [c for c in per.index if per[c] >= 5]
    rows, matched_all, swaps_all, status = [], {}, {}, []
    for c in cancers:
        t = cancer_frame(c, frames, donors, split, P, extra=(CALLN, CALLM), mal_key=mal_key)
        if t is None:
            status.append(dict(cancer=c, status="no host"))
            log(f"{c}: no host")
            continue
        if args.smoke:  # donor labels of every non-anchor row permuted within stratum: code path only
            rng = np.random.default_rng(1)
            for st, g in t["G"][~t["G"].compartment.isin(t["classes"])].groupby(["compartment", "stratum"]):
                t["G"].loc[g.index, "donor"] = g.donor.to_numpy()[rng.permutation(len(g))]
        fits, matched, swaps = score_cancer(c, t, args, n_set, n_swap)
        ok = {k: bool(len(v)) for k, v in dict(fits=fits["14"], matched=matched, swaps=swaps).items()}
        status.append(dict(cancer=c, status="ok", n_split=t["n_split"], classes="+".join(k.split(":")[1] for k in t["classes"]),
                           n_sets=int(matched.set.nunique()) if "set" in matched else 0, n_swaps=len(swaps), **ok))
        log(f"{c}: scored ({t['n_split']} split patients; {n_set} sets, {n_swap} swaps)")
        if args.smoke:
            continue
        matched_all[c], swaps_all[c] = matched, swaps
        matched.assign(cancer=c).to_csv(out / f"matched_{c}.csv", index=False)
        swaps.assign(cancer=c).to_csv(out / f"swaps_{c}.csv", index=False)
        for rd, f in fits.items():
            m = matched[matched.readout == rd]
            for k, r in f.iterrows():
                v = m[m.target == k].coupling.dropna().to_numpy() if "target" in m else np.array([])
                if k in t["classes"] or len(v) < min(10, n_set):
                    continue
                pp = P[P.split & (P.cancer == c)]
                rows.append(dict(cancer=c, readout=rd, target=k.split(":", 1)[1], patients=int(r.n), coupling=float(r.coupling),
                                 matched_mean=float(v.mean()), matched_sd=float(v.std(ddof=1)), matched_median=float(np.median(v)),
                                 matched_q95=float(np.quantile(v, .95)), z=float((r.coupling - v.mean()) / v.std(ddof=1)),
                                 e=max(float(r.coupling) - float(np.median(v)), 0.0), share_matched_ge=float(np.mean(v >= r.coupling)),
                                 host_reliability=float(r.host_reliability), n_matched_ok=len(v),
                                 median_umi_high=float(pp.median_umi_high.median()), median_umi_low=float(pp.median_umi_low.median())))
    pd.DataFrame(status).to_csv(out / "status.csv", index=False)
    if args.smoke:
        log(f"smoke done {(time.time() - t0) / 60:.1f} min: {len(status)} cancers, no value written")
        return
    R = pd.DataFrame(rows)
    R.to_csv(out / "per_cancer.csv", index=False)
    tests(R, matched_all, swaps_all, P, out)
    log(f"score done {(time.time() - t0) / 60:.1f} min")


def tests(R, matched_all, swaps_all, P, out):
    """H423.1, H423.2 (Holm), D423.1-D423.3 from the per-cancer fits."""
    pd.set_option("display.width", 250, "display.max_columns", 40, "display.max_rows", 300)
    w = P[P.split].groupby("cancer").donor.nunique()
    H, rows = [], []
    for rd in READOUTS:
        r = R[R.readout == rd].set_index(["cancer", "target"])
        cs_ = [c for c in w.index if (c, "epi_high") in r.index and (c, "epi_low") in r.index]
        if not cs_:
            continue
        n = np.array([w[c] for c in cs_], float)
        d_obs = float(sum(n[k] * (r.coupling[(c, "epi_low")] - r.coupling[(c, "epi_high")]) for k, c in enumerate(cs_)) / n.sum())
        sets = pd.concat([matched_all[c][matched_all[c].readout == rd].pivot_table(index="set", columns="target", values="coupling")
                          .assign(cancer=c) for c in cs_])
        dm = pd.DataFrame({c: sets[sets.cancer == c][LOW].sub(sets[sets.cancer == c][HIGH]) for c in cs_}).dropna()
        d_m = (dm * (n / n.sum())).sum(1).to_numpy()
        row = dict(hypothesis="H423.1", readout=rd, cancers="+".join(cs_), patients=int(n.sum()), delta=d_obs,
                   matched_median=float(np.median(d_m)), matched_q95=float(np.quantile(d_m, .95)), n_matched=len(d_m),
                   p_matched=float((1 + np.sum(d_m >= d_obs)) / (len(d_m) + 1)))
        if rd == "14":
            nl = []
            for k in range(1000):
                v = [swaps_all[c].set_index("draw").reindex([k]) for c in cs_]
                nl.append(sum(n[i] * (v[i].c_low.iloc[0] - v[i].c_high.iloc[0]) for i in range(len(cs_))) / n.sum())
            nl = np.array(nl, float)
            row.update(swap_q95=float(np.nanquantile(nl, .95)), p_swap=float((1 + np.sum(nl[np.isfinite(nl)] >= d_obs)) / (np.isfinite(nl).sum() + 1)))
        H.append(row)
        h2 = [c for c in H2_CANCERS if c in cs_]
        if h2:
            z = np.array([r.z[(c, "epi_high")] for c in h2])
            for wname, wt in (("patient-weighted", np.array([w[c] for c in h2], float)), ("unweighted", np.ones(len(h2)))):
                Z = float((wt * z).sum() / np.sqrt((wt ** 2).sum()))
                H.append(dict(hypothesis="H423.2", readout=rd, weighting=wname, cancers="+".join(h2), patients=int(sum(w[c] for c in h2)),
                              z_per_cancer=";".join(f"{c}:{x:.3f}" for c, x in zip(h2, z)), Z=Z, p_one_sided=float(norm.sf(Z)),
                              passes=bool(Z > 1.645)))
    Hd = pd.DataFrame(H)
    p14 = Hd[(Hd.readout == "14") & (Hd.hypothesis == "H423.1")]
    p2 = Hd[(Hd.readout == "14") & (Hd.hypothesis == "H423.2") & (Hd.weighting == "patient-weighted")]
    if len(p14) and len(p2):
        ps = [max(float(p14.p_swap.iloc[0]), float(p14.p_matched.iloc[0])), float(p2.p_one_sided.iloc[0])]
        Hd.loc[p14.index, "p_primary"], Hd.loc[p2.index, "p_primary"] = ps
        Hd.loc[p14.index, "holm_p"], Hd.loc[p2.index, "holm_p"] = holm_rows(ps)
        Hd.loc[p14.index, "passes"] = bool(ps[0] < .05 and Hd.loc[p14.index, "delta"].iloc[0] > Hd.loc[p14.index, "matched_q95"].iloc[0])
        Hd["passes_holm"] = np.nan
        Hd.loc[p14.index, "passes_holm"] = bool(Hd.loc[p14.index, "holm_p"].iloc[0] < .05 and Hd.loc[p14.index, "passes"].iloc[0])
        Hd.loc[p2.index, "passes_holm"] = bool(Hd.loc[p2.index, "holm_p"].iloc[0] < .05)
    Hd.to_csv(out / "h423.csv", index=False)
    log("H423\n" + Hd.round(4).to_string(index=False))
    R[R.target.isin(["malignant|C", "malignant|A", "malignant|E"])].to_csv(out / "d423_1_label_class.csv", index=False)
    R[R.target.isin(["epi_called_normal", "epi_called_malignant"])].to_csv(out / "d423_2_called.csv", index=False)
    d3 = R[R.target.isin(["epi_all", "epi_high", "epi_low"])].pivot_table(
        index=["cancer", "readout"], columns="target", values=["coupling", "z", "patients", "median_umi_high", "median_umi_low"]).round(4)
    d3.columns = [f"{a}_{b}" for a, b in d3.columns]
    d3.reset_index().to_csv(out / "d423_3_per_cancer.csv", index=False)
    # H401.2 with e_mal from c_HIGH (A401 selection: matched cancers passing P401 on the readout)
    a401 = pd.read_csv(U.OUT / "score" / "per_cancer.csv")
    for rd in READOUTS:
        a = a401[(a401.readout.astype(str) == rd) & a401.p401_pass & a401.matched & a401.P.notna()].set_index("cancer")
        e = R[(R.readout == rd) & (R.target == "epi_high")].set_index("cancer").e
        both = [c for c in a.index if c in e.index]
        rho, p = U.spearman_exact(e.reindex(both), a.a83_share.reindex(both)) if len(both) >= 3 else (np.nan, np.nan)
        mix = a.e_mal.copy()
        mix.loc[both] = e.reindex(both)
        rho2, p2_ = U.spearman_exact(mix, a.a83_share)
        rows.append(dict(readout=rd, set="cancers with a split", n=len(both), cancers="+".join(both), rho=rho, p_exact_one_sided=p))
        rows.append(dict(readout=rd, set="all A401 cancers (e_mal replaced where split)", n=len(a), cancers="+".join(a.index), rho=rho2, p_exact_one_sided=p2_))
        rows.append(dict(readout=rd, set="A401 e_mal (reference)", n=len(a), cancers="+".join(a.index), rho=U.spearman_exact(a.e_mal, a.a83_share)[0],
                         p_exact_one_sided=U.spearman_exact(a.e_mal, a.a83_share)[1]))
    pd.DataFrame(rows).to_csv(out / "d423_3_h401_2.csv", index=False)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stage", choices=("power", "score"), required=True)
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--genome-jobs", type=int, default=2)
    ap.add_argument("--smoke", action="store_true")
    args = ap.parse_args()
    limit_threads(1)
    stage_score(args) if args.stage == "score" else stage_power(args)


if __name__ == "__main__":
    main()
