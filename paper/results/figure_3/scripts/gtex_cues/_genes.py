"""Genome-wide GTEx helpers of the cue family (A81, A81b, A89-A91, A100-A102b): tissue lists, gene classes, per-tissue
expression with all annotations, the donor halves of A81, and the gene x tissue selection table (rhythm in fast deaths,
lock to the host in ventilated donors). Host clock states come from `clocksync` (`Dataset.fit(...).host` through
`_shared.host_states`); the per-gene regressions are local ("local: not in clocksync"). Controlled access: per-donor values
stay in memory.
"""

from __future__ import annotations

import os
import sys
from pathlib import Path

# joblib workers import numpy afresh: cap their BLAS threads (threadpool_limits in the parent does not reach them)
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")

import numpy as np
import pandas as pd
from joblib import Parallel, delayed

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import clocksync as cs  # noqa: E402
from clocksync.genes import canonical, gene_columns  # noqa: E402
from common import gtex  # noqa: E402
from common.paths import DATA  # noqa: E402

import _shared as S  # noqa: E402
from _shared import TISSUES, W  # noqa: E402

FOLDER = "gtex_cues"
CORE_GENES = tuple(cs.gtex().genes)
CLOCK = set(CORE_GENES) | {"BMAL1"}
# D6 clock outputs (legacy glucocorticoid/diag6_classes.CLASSES["clock output (non-frame)"])
CLOCK_OUTPUT = {"HLF", "NFIL3", "BHLHE40", "BHLHE41", "RORC", "NPAS2", "ARNTL2", "CRY2", "TSPAN4", "USP2", "LONRF3",
                "TEF", "NAMPT", "KLF10", "DDIT4L", "WEE1", "PPARA", "CIPC"}
EXCLUDE = CLOCK | CLOCK_OUTPUT
A81_TISSUES = [t for t in TISSUES if t != "Pituitary" and not t.startswith("Brain")]   # 12 (with whole blood)
A81B_TISSUES = [t for t in A81_TISSUES if t != "Whole Blood"]                          # 11
A84_ADDED = ["Colon - Transverse", "Colon - Sigmoid", "Liver", "Breast - Mammary Tissue", "Prostate",
             "Kidney - Cortex", "Uterus"]
A89_EXTRA = ["Stomach", "Esophagus - Muscularis", "Esophagus - Gastroesophageal Junction", "Adrenal Gland",
             "Heart - Left Ventricle", "Artery - Aorta", "Artery - Coronary", "Testis"]
ALL26 = A81B_TISSUES + A84_ADDED + A89_EXTRA
N_SET, SEED_HALVES = 40, 81


# ---- donors and rows ---------------------------------------------------------------------------------------------

def obs_table() -> pd.DataFrame:
    """Sample annotations of the labelled GTEx table (no counts read)."""
    import anndata as ad
    a = ad.read_h5ad(gtex.LABELLED, backed="r")
    try:
        o = a.obs[["DTHHRDY", "tissue", "tod", "subject"]].copy()
    finally:
        a.file.close()
    o["subject"] = o.subject.astype(str)
    o["tissue"] = o.tissue.astype(str)
    return o


def pick_donors(o: pd.DataFrame, n: int | None, seed: int = 0):
    """Subsample check: n random donors with a time of death in each of fast deaths and ventilated cohorts (a set), or
    None for the full data."""
    if n is None:
        return None
    rng = np.random.default_rng(seed)
    keep = set()
    for hardy in ((1, 2), (0,)):
        s = np.array(sorted(o.loc[o.DTHHRDY.isin(hardy) & o.tod.notna(), "subject"].unique()))
        keep |= set(rng.permutation(s)[:n])
    return keep


def halves(o: pd.DataFrame, donors=None) -> dict:
    """A81 selection / evaluation donors per Hardy group (donors with a time of death), split at random (seed 81:
    fast first, then ventilated, from one generator). `donors`: restrict to a subsample."""
    rng = np.random.default_rng(SEED_HALVES)
    out = {}
    for grp, hardy in (("fast", [1, 2]), ("vent", [0])):
        s = np.array(sorted(o.loc[o.DTHHRDY.isin(hardy) & o.tod.notna(), "subject"].unique()))
        if donors is not None:
            s = np.array([x for x in s if x in donors])
        s = s[rng.permutation(len(s))]
        out[grp] = {"sel": set(s[: len(s) // 2]), "eval": set(s[len(s) // 2:])}
    return out


def load_rows(tissue: str, hardy, subjects=None, need_tod: bool = True, genes=None):
    """CPM (library size; float) of the rows of one tissue and their annotations. `genes`: only these (by symbol,
    aliases matched) in the order given, else every gene (local: not in clocksync)."""
    import anndata as ad
    a = ad.read_h5ad(gtex.LABELLED, backed="r")
    try:
        o = a.obs
        keep = (o.DTHHRDY.isin(list(hardy)) & (o.tissue.astype(str) == tissue)).to_numpy()
        if need_tod:
            keep = keep & o.tod.notna().to_numpy()
        if subjects is not None:
            keep = keep & o.subject.astype(str).isin(subjects).to_numpy()
        idx = np.flatnonzero(keep)
        names = np.asarray(a.var_names)
        if genes is None:
            X = a.X[idx]
            gsel = names
        else:
            cols = gene_columns(names, genes)
            ok = cols >= 0
            X = a.X[idx][:, cols[ok]]
            gsel = np.asarray(genes)[ok]
        X = X.toarray() if hasattr(X, "toarray") else np.asarray(X)
        ob = o.iloc[idx].copy()
    finally:
        a.file.close()
    cpm = X / ob.library_size.to_numpy(float)[:, None] * 1e6
    obs = pd.DataFrame({"subject": ob.subject.astype(str).to_numpy(), "hardy": ob.DTHHRDY.astype(int).to_numpy(),
                        "tod": ob.tod.to_numpy(float) % 24, "sex": ob.sex.astype(str).to_numpy(),
                        "age_bracket": ob.age_bracket.astype(str).to_numpy(),
                        "ischemic_min": ob.ischemic_time.astype(float).to_numpy(),
                        "rin": ob.rin.astype(float).to_numpy(), "center": ob.center.astype(str).to_numpy(),
                        "library_size": ob.library_size.to_numpy(float)}, index=ob.index)
    return cpm, gsel, obs


def hosts_for(tissues, hardy, donors=None, jobs: int = 1) -> dict:
    """tissue -> host DataFrame (x, y by donor) from the clock genes of the other `TISSUES` (leave one out; tissues
    outside `TISSUES` share the full 15-tissue host), clocksync `Dataset.fit(...).host`."""
    pb, _ = S.clock_table(hardy=hardy)
    sets = {t: tuple(x for x in TISSUES if x != t) for t in tissues}
    uniq = sorted(set(sets.values()))
    fits = {s: S.host_states(pb, list(s), donors=donors) for s in uniq}
    return {t: fits[sets[t]] for t in tissues}


# ---- selection table (A81 / A81b / A84 / A89) ---------------------------------------------------------------------

def select_tissue(tis, Yf, obf, Yv, obv, genes, fv) -> pd.DataFrame:
    """Per gene of one tissue: p of the timed cosinor in fast deaths, its amplitude and acrophase (h), and the lock =
    |host loading| in ventilated donors over |timed loading|, relative to the median of that ratio over the tissue's
    rhythmic clock-plane genes (p_fast < 0.01)."""
    tf = obf.tod.to_numpy(float) % 24
    c_fast, p_fast = S.timed_cosinor(Yf, tf, S.d2_covariates(obf))
    subj = obv.subject.astype(str).to_numpy()
    m = np.isin(subj, fv.index)
    f = fv[["x", "y"]].loc[subj[m]].to_numpy()
    f = f - f.mean(0)
    Z = S.d2_covariates(obv[m])
    Z = Z[:, Z.std(0) > 0]
    B, _ = S.lstsq_fit(Yv[m], f, Z)
    c_int = B[-2] + 1j * B[-1]
    ratio = np.abs(c_int) / np.abs(c_fast)
    frame = np.isin(genes, list(CLOCK)) & (p_fast < 0.01)
    return pd.DataFrame(dict(tissue=tis, gene=genes, p_fast=p_fast, amp_fast=np.abs(c_fast),
                             acro_h=np.angle(c_fast) / W % 24, lock=ratio / np.median(ratio[frame])))


def _one_selection(tis, fast_sub, vent_sub, fv):
    Xf, gf, obf = load_rows(tis, (1, 2), fast_sub)
    Xv, gv, obv = load_rows(tis, (0,), vent_sub)
    expr = (np.median(Xf, 0) >= 5) & (np.median(Xv, 0) >= 5)
    return select_tissue(tis, np.log(Xf[:, expr] + 0.5), obf, np.log(Xv[:, expr] + 0.5), obv, gf[expr], fv)


def selection(tissues, fast_sub, vent_sub, vent_hosts_donors, donors=None, jobs: int = 1) -> pd.DataFrame:
    """Selection table of `tissues`: rhythm from the `fast_sub` donors, lock from the `vent_sub` ventilated donors with
    the host fitted on `vent_hosts_donors` (the same donors)."""
    H = hosts_for(tissues, (0,), donors=vent_hosts_donors)
    res = Parallel(n_jobs=jobs)(delayed(_one_selection)(t, fast_sub, vent_sub, H[t]) for t in tissues)
    return pd.concat(res, ignore_index=True)


def selection_tables(out: Path, subsample: int | None = None, seed: int = 0, jobs: int = 1) -> pd.DataFrame:
    """A81b selection on all fast deaths and the ventilated selection half (A81 split) for the 26 tissues of A81b, A84
    and A89, cached as `selection_genes_26.csv.gz` in `out`; also writes `a81b_selection_genes.csv.gz` (the 11 tissues
    of A81b, as the legacy file)."""
    f = out / "selection_genes_26.csv.gz"
    if f.exists():
        return pd.read_csv(f)
    o = obs_table()
    sub = pick_donors(o, subsample, seed)
    H = halves(o, sub)
    fast = H["fast"]["sel"] | H["fast"]["eval"]
    g = selection(ALL26, fast, H["vent"]["sel"], H["vent"]["sel"], jobs=jobs)
    g.to_csv(f, index=False)
    g[g.tissue.isin(A81B_TISSUES)].to_csv(out / "a81b_selection_genes.csv.gz", index=False)
    return g


# ---- per-gene aggregation across tissues (A81b) -----------------------------------------------------------------

def per_gene(g: pd.DataFrame) -> pd.DataFrame:
    """Across tissues: number of tissues, Fisher combined p, weighted coherence Rw, lock (median over tissues with
    p_fast < 0.05), amplitude |sum c| / nt and acrophase (h) of the summed complex loading."""
    from scipy import stats
    g = g.assign(c=g.amp_fast * np.exp(1j * W * g.acro_h))
    rows = {}
    for gene, d in g.groupby("gene"):
        rows[gene] = dict(nt=len(d), fisher=stats.chi2.sf(-2 * np.log(d.p_fast).sum(), 2 * len(d)),
                          Rw=abs(d.c.sum()) / d.c.abs().sum(), lock=d.loc[d.p_fast < 0.05, "lock"].median(),
                          amp=abs(d.c.sum()) / len(d), acro_h=np.angle(d.c.sum()) / W % 24)
    return pd.DataFrame.from_dict(rows, orient="index").rename_axis("gene")


def cue_sets(agg: pd.DataFrame) -> tuple[pd.DataFrame, float, int]:
    """A81b eligible genes (nt >= 8, Fisher p < 1e-4, Rw >= 0.6, not clock plane / D6 clock output) and the cue (lock
    below the median) and locked (at or above) sets, top 40 by amplitude each."""
    el = agg[(agg.nt >= 8) & (agg.fisher < 1e-4) & (agg.Rw >= 0.6) & ~agg.index.isin(EXCLUDE)]
    med = el.lock.median()
    sets = pd.concat([el[el.lock < med].nlargest(N_SET, "amp").assign(set="cue"),
                      el[el.lock >= med].nlargest(N_SET, "amp").assign(set="locked")])
    return sets, float(med), len(el)


def waveform_of(x: pd.DataFrame, name: str) -> cs.waveform.Waveform:
    """Reference waveform of a gene set: amplitude and acrophase (GTEx time-of-death frame; acro_h hours)."""
    return cs.waveform.Waveform(tuple(x.index), x.amp.to_numpy(), x.acro_h.to_numpy() * W, None, None, name)


def circ_corr_js(a, b) -> float:
    """Jammalamadaka-SenGupta circular correlation (radians in); the guard statistic of the A81b frame check."""
    sa = np.sin(a - np.angle(np.exp(1j * a).mean()))
    sb = np.sin(b - np.angle(np.exp(1j * b).mean()))
    return float(np.sum(sa * sb) / np.sqrt(np.sum(sa ** 2) * np.sum(sb ** 2)))


def frame_check(agg: pd.DataFrame) -> dict:
    """A81b guard: the 14 clock genes processed like every other gene must reproduce the fixed GTEx waveform (circular
    correlation of acrophases > 0.8; the legacy run gave 0.48 and was overridden post hoc)."""
    core = cs.gtex().table()
    idx = {canonical(x): x for x in agg.index}
    common = [g for g in core.index if g in idx]
    a = W * core.loc[common, "psi_h"].to_numpy()
    b = W * agg.loc[[idx[g] for g in common], "acro_h"].to_numpy()
    return dict(genes=len(common), circ_corr=circ_corr_js(a, b),
                mean_shift_h=float(np.angle(np.exp(1j * (b - a)).mean()) / W))


def circ_R(h) -> float:
    return float(abs(np.exp(1j * W * np.asarray(h)).mean()))


def hallmarks():
    """Hallmark gene sets (MSigDB gmt under $SCCC_DATA/reference), or None when not present."""
    for f in sorted((DATA / "reference").glob("h.all.*.Hs.symbols.gmt"), reverse=True):
        return {l.split("\t")[0]: set(l.strip().split("\t")[2:]) for l in open(f)}
    return None
