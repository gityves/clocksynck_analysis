"""GTEx cue reference and the disease items of A81b (H81.2-3), A82, A84 and A88.

The cue reference (A81b: cue = rhythm follows timing cues, locked = clock-output-like) and the 26-tissue selection table are
outputs of `gtex_cues/compute_cue_reference.py` (documented files, aggregates over genes and tissues).
Local helpers (not in clocksync): the per-gene aggregation of the selection table and the tissue-matched sets of A84 (copies
of the S07 / legacy rules, since the S07 helper module is not importable next to `_shared`), the per-subject timed-skin
offsets of A81b H81.2, the TCGA all-gene paired table, a GTEx row reader for A88.
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import numpy as np
import pandas as pd
from scipy import stats as st

import clocksync as cs
from clocksync.genes import canonical, gene_columns
from common.paths import DATA

import _shared as H

W = H.W
sys.path.insert(0, str(Path(__file__).resolve().parents[4]))  # paper/
from common.paths import group_out  # noqa: E402
SEL = group_out("gtex_cues")
REF_FILE, SEL26 = SEL / "a81b_cue_reference.csv", SEL / "selection_genes_26.csv.gz"
N_SET = 40
CLOCK = set(H.GENES) | {"BMAL1"}
CLOCK_OUTPUT = {"HLF", "NFIL3", "BHLHE40", "BHLHE41", "RORC", "NPAS2", "ARNTL2", "CRY2", "TSPAN4", "USP2", "LONRF3", "TEF",
                "NAMPT", "KLF10", "DDIT4L", "WEE1", "PPARA", "CIPC"}
EXCLUDE = CLOCK | CLOCK_OUTPUT
TISSUE = {"colon": ["Colon - Transverse", "Colon - Sigmoid"], "liver": ["Liver"], "breast": ["Breast - Mammary Tissue"],
          "prostate": ["Prostate"], "kidney": ["Kidney - Cortex"], "uterus": ["Uterus"], "lung": ["Lung"],
          "thyroid": ["Thyroid"], "skin": ["Skin - Not Sun Exposed (Suprapubic)", "Skin - Sun Exposed (Lower leg)"]}
CANCER = {"COAD": "colon", "LIHC": "liver", "BRCA": "breast", "PRAD": "prostate", "KIRC": "kidney", "KIRP": "kidney",
          "UCEC": "uterus", "LUAD": "lung", "THCA": "thyroid", "BLCA": None}


def _wf(x: pd.DataFrame, name: str) -> cs.Waveform:
    return cs.Waveform(tuple(canonical(g) or g for g in x.index), x.amp.to_numpy(), x.acro_h.to_numpy() * W, None, None, name)


def reference() -> pd.DataFrame:
    return pd.read_csv(REF_FILE, index_col=0)


def cue_waveforms(systemic: cs.Waveform | None = None) -> dict:
    """core (14-gene GTEx waveform), cue and locked (A81b reference) and, when given, the A78 systemic set."""
    ref = reference()
    out = {"core": H.WF, **{k: _wf(ref[ref.set == k], k) for k in ("cue", "locked")}}
    if systemic is not None:
        out["systemic"] = systemic
    return out


# ---- disease items --------------------------------------------------------------------------------------------------

def gated_cohorts(samples: pd.DataFrame, per_series: pd.DataFrame, key="cohort"):
    """(samples of the gated cohorts, {cohort: 1 / se^2}) from a per-series table (`gated`, `se`)."""
    ps = per_series[per_series.gated.astype(bool)].set_index(key)
    s = samples[samples.cohort.isin(ps.index)]
    return s, {k: 1 / ps.loc[k, "se"] ** 2 for k in ps.index if k in set(s.cohort)}


def tcga_table(per_cancer: pd.DataFrame, genes, subsample: int | None = None, seed: int = 0):
    """Gated TCGA cancers as one paired sample table with the columns `genes` (and every clock gene): unit = patient,
    condition normal / tumour, cohort = cancer, depth = all-gene count sum. `per_cancer`: cancer
    `a52_per_cancer.csv` (cohort `cancer`, `gated`, `se`)."""
    b = H.builder("tcga_paired")
    pc = per_cancer[per_cancer.gated.astype(bool)].set_index("cancer")
    parts, w = [], {}
    rng = np.random.default_rng(seed)
    for c in pc.index:
        t, n, ct, cn = b._pairs(c)
        ct.index, cn.index = t.patient.to_numpy(), n.patient.to_numpy()
        u = sorted(set(ct.index) & set(cn.index))
        if subsample:
            u = sorted(rng.permutation(u)[:subsample])
        keep = [g for g in dict.fromkeys(list(genes) + H.GENES) if g in ct.columns]
        for cond, x in (("tumour", ct), ("normal", cn)):
            x = x.loc[u]
            d = x[keep].copy()
            d.insert(0, "depth", x.sum(1).to_numpy())
            d.insert(0, "cohort", c)
            d.insert(0, "condition", cond)
            d.insert(0, "unit", u)
            d["value_type"] = "counts"
            parts.append(d)
        w[c] = 1 / pc.loc[c, "se"] ** 2
    return pd.concat(parts, ignore_index=True), w


def all_set_genes(wfs: dict) -> list[str]:
    return sorted({g for k, wf in wfs.items() for g in wf.genes})


# ---- A84: tissue-matched sets ----------------------------------------------------------------------------------------

def per_gene(g: pd.DataFrame) -> pd.DataFrame:
    """local (copy of S07 `_genes.per_gene`): across tissues number of tissues, Fisher p, coherence Rw, lock (median over
    tissues with p_fast < 0.05), amplitude and acrophase of the summed complex loading."""
    g = g.assign(c=g.amp_fast * np.exp(1j * W * g.acro_h))
    rows = {}
    for gene, d in g.groupby("gene"):
        rows[gene] = dict(nt=len(d), fisher=st.chi2.sf(-2 * np.log(d.p_fast).sum(), 2 * len(d)),
                          Rw=abs(d.c.sum()) / d.c.abs().sum(), lock=d.loc[d.p_fast < 0.05, "lock"].median(),
                          amp=abs(d.c.sum()) / len(d), acro_h=np.angle(d.c.sum()) / W % 24)
    return pd.DataFrame.from_dict(rows, orient="index").rename_axis("gene")


def selection_tables() -> tuple[pd.DataFrame, pd.DataFrame]:
    """(A81b tissues, all 26 tissues) of the gene x tissue selection table; the A81b pool uses the tissues of
    `a81b_selection_genes.csv.gz` (11 tissues, rhythm from all fast deaths, lock from the ventilated selection half)."""
    g26 = pd.read_csv(SEL26)
    t81 = pd.read_csv(SEL / "a81b_selection_genes.csv.gz", usecols=["tissue"]).tissue.unique()
    return g26[g26.tissue.isin(t81)], g26


def pools(g81b: pd.DataFrame) -> dict:
    """A84 candidate pools: the eligible A81b genes (nt >= 8, Fisher p < 1e-4, Rw >= 0.6, not clock / clock output) split
    at their median lock."""
    agg = per_gene(g81b)
    el = agg[(agg.nt >= 8) & (agg.fisher < 1e-4) & (agg.Rw >= 0.6) & ~agg.index.isin(EXCLUDE)]
    med = el.lock.median()
    return {"cue": el[el.lock < med], "locked": el[el.lock >= med]}


def matched(g: pd.DataFrame, pool: dict, tissues, fallback: pd.DataFrame):
    """Tissue-matched cue / locked sets: candidates expressed in the target tissue, |acrophase - global| < 3 h,
    p_fast,T < 0.2, top 40 by amplitude in T (waveform from T); < 15 genes: the global A81b set restricted to genes
    expressed in T (`fallback`: the reference table)."""
    t = g[g.tissue.isin(tissues)].assign(c=lambda d: d.amp_fast * np.exp(1j * W * d.acro_h))
    t = t.groupby("gene").agg(c=("c", "mean"), p=("p_fast", "min"))  # colon: complex mean of the two segments
    out, rows = {"core": H.WF}, []
    for k, cand in pool.items():
        x = t.loc[t.index.intersection(cand.index)].copy()
        x["amp"], x["acro_h"] = x.c.abs(), np.angle(x.c) / W % 24
        dphi = np.abs(((x.acro_h - cand.loc[x.index, "acro_h"]) + 12) % 24 - 12)
        sel = x[(dphi < 3) & (x.p < 0.2)].nlargest(N_SET, "amp")
        fb = len(sel) < 15
        if fb:
            sel = fallback[(fallback.set == k) & fallback.index.isin(x.index)][["amp", "acro_h"]]
        out[k] = _wf(sel, k)
        rows.append(dict(set=k, genes=len(sel), fallback=fb))
    return out, pd.DataFrame(rows)


# ---- A81b H81.2: timed skin, per-subject epidermis-minus-dermis offset of the core and the cue state ------------------

def timed_core_cue(wfs: dict, n_boot: int, seed: int = 82) -> tuple[dict, pd.DataFrame]:
    """Per subject (GSE205155, 11 subjects) the epidermis-minus-dermis offset of the core and of the cue state (A76
    complex regression of the state on clock time); mean core minus mean cue, subject bootstrap, one-sided p (upper)."""
    s, m = H.timed_skin()
    off = {}
    for k in ("core", "cue"):
        g = [x for x in wfs[k].genes if x in m.index]
        wf = wfs[k].subset(g) if k == "core" else H.restrict(wfs[k], g)
        y = m.loc[g, s.gsm].T * np.log(2.0)
        y.index = s.index
        L = H.PD.projector(y, wf, counts=False)
        b = (y - y.mean(1).to_numpy()[:, None]).to_numpy() @ L.T
        c = pd.Series(b[:, 0] + 1j * b[:, 1], index=s.index)
        o = {}
        for subj, d in s.groupby("subject"):
            a = {}
            for layer, x in d.groupby("layer"):
                e = np.exp(1j * W * x.clock_h.to_numpy(float))
                Xd = np.column_stack([np.ones_like(e), e, np.conj(e)])
                a[layer] = np.linalg.lstsq(Xd, c.loc[x.index].to_numpy(), rcond=None)[0][1]
            o[subj] = np.angle(a["epidermis"] / a["dermis"]) / W
        off[k] = (pd.Series(o), len(g))
    d = pd.DataFrame({k: v[0] for k, v in off.items()})
    cmean = lambda x: float(np.angle(np.mean(np.exp(1j * W * np.asarray(x)))) / W)
    diff = float(cs.stats.wrap(cmean(d.core) - cmean(d.cue)))
    rng = np.random.default_rng(seed)
    boot = np.array([cs.stats.wrap(cmean(d.core.iloc[i]) - cmean(d.cue.iloc[i]))
                     for i in (rng.integers(0, len(d), len(d)) for _ in range(n_boot))])
    v = dict(subjects=len(d), core_genes=off["core"][1], cue_genes=off["cue"][1], core_h=cmean(d.core), cue_h=cmean(d.cue),
             core_minus_cue_h=diff, p_upper=float((1 + np.sum(boot <= 0)) / (1 + n_boot)),
             lo=float(np.quantile(boot, 0.025)), hi=float(np.quantile(boot, 0.975)),
             subjects_core_gt_cue=int((cs.stats.wrap(d.core - d.cue) > 0).sum()))
    v["H81_2"] = "PASS" if (v["core_minus_cue_h"] > 0 and v["p_upper"] < 0.05) else "FAIL"
    return v, d


# ---- A88: GTEx rows ---------------------------------------------------------------------------------------------------

def gtex_rows(tissue: str, hardy, subjects, genes, depth: str = "library_size") -> pd.DataFrame:
    """local: counts of `genes` (symbols, aliases matched; absent ones dropped) of the labelled GTEx rows of one tissue and
    Hardy group(s) with a time of death, donors in `subjects` (None: all), plus `depth` (annotated library size, or the
    count sum over all genes), `hardy`, `tod`, `sex`, `age_bracket`, `subject` (first sample of a donor). Controlled access:
    per-donor values stay in memory."""
    import anndata as ad
    from common import gtex
    a = ad.read_h5ad(gtex.LABELLED, backed="r")
    try:
        o = a.obs
        keep = (o.DTHHRDY.isin(list(hardy)) & (o.tissue.astype(str) == tissue) & o.tod.notna()
                & (o.subject.astype(str).isin(subjects) if subjects is not None else True)).to_numpy()
        idx = np.flatnonzero(keep)
        cols = gene_columns(np.asarray(a.var_names), list(genes))
        ok = cols >= 0
        Xa = a.X[idx]
        tot = np.asarray(Xa.sum(1)).ravel()
        X = Xa[:, cols[ok]]
        X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
        ob = o.iloc[idx]
        d = pd.DataFrame(X, columns=np.asarray(list(genes))[ok])
        d["depth"] = tot if depth == "sum" else ob.library_size.to_numpy(float)
        d["hardy"] = ob.DTHHRDY.astype(int).to_numpy()
        d["tod"] = ob.tod.to_numpy(float) % 24
        d["sex"] = ob.sex.astype(str).to_numpy()
        d["age_bracket"] = ob.age_bracket.astype(str).to_numpy()
        d["subject"] = ob.subject.astype(str).to_numpy()
    finally:
        a.file.close()
    return d[~d.subject.duplicated().to_numpy()].reset_index(drop=True)
