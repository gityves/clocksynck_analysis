"""Loaders and small local helpers shared by cancer and cancer_pipeline.

Method calls are clocksync's. Everything marked "local: not in clocksync" is a short helper for a statistic the
package has no public call for (listed in the folder READMEs).
"""

from __future__ import annotations

import importlib.util
import os
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[4]))
import clocksync as cs  # noqa: E402
from clocksync.decompose import logfrac  # noqa: E402
from common import data  # noqa: E402
from common.paths import DATA, PSEUDOBULK  # noqa: E402

WF = cs.gtex()
GENES = list(WF.genes)
HOURS = 24 / (2 * np.pi)
CANCERS = ["BLCA", "BRCA", "COAD", "KIRC", "KIRP", "LIHC", "LUAD", "PRAD", "THCA", "UCEC"]
LEVELS = ("normal", "tumour")  # reference, case
MAL, NORM = "tumour:malignant", "adjacent_normal:epithelial"
STROMA = ("fibroblast", "endothelial", "mural")


# ---- TCGA ---------------------------------------------------------------------------------------------------------

def tcga_builder():
    """paper/common/datasets/tcga_paired/build.py as a module (for `tumours`, `gene_means`, `samples_of`, `counts_of`)."""
    return data._builder("tcga_paired")


def load_tcga(subsample: int | None = None, seed: int = 0) -> pd.DataFrame:
    """Paired tumour / normal samples of the ten cancers (one row per sample; unit = patient, condition, cohort =
    cancer, depth, clinical fields, 14 clock genes). Standard table `$SCCC_DATA/pseudobulk/tcga_paired.parquet`
    (built from recount3 by the builder when missing). `subsample`: that many random patients per cancer."""
    f = PSEUDOBULK / "tcga_paired.parquet"
    if f.exists():
        s = pd.read_parquet(f)
    else:
        s, _ = tcga_builder().build()
        PSEUDOBULK.mkdir(parents=True, exist_ok=True)
        s.to_parquet(f, index=False)
    if subsample is not None:
        rng = np.random.default_rng(seed)
        keep = []
        for _, g in s.groupby("cohort", sort=False):
            u = sorted(g.unit.unique())
            keep += list(rng.permutation(u)[:subsample])
        s = s[s.unit.isin(keep)]
    return s.reset_index(drop=True)


def tcga_full_counts(cancer: str, units) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Genome-wide scaled counts (patients x gene symbols) of the primary tumours and solid-tissue normals of the
    given patients of one cancer (recount3, as the builder; index = patient)."""
    t, n, ct, cn = tcga_builder()._pairs(cancer)
    ct.index, cn.index = t.patient.to_numpy(), n.patient.to_numpy()
    u = [x for x in units if x in ct.index]
    return ct.loc[u], cn.loc[u]


def pair_frames(samples: pd.DataFrame, cohort: str, levels=LEVELS):
    """(case, reference) rows of one cancer aligned by patient (as `clocksync.paired`)."""
    x = samples[samples.cohort == cohort]
    case = x[x.condition == levels[1]].drop_duplicates("unit").set_index("unit")
    ref = x[x.condition == levels[0]].drop_duplicates("unit").set_index("unit")
    u = sorted(set(case.index) & set(ref.index))
    return case.loc[u], ref.loc[u]


PROJECTS = CANCERS + ["OV"]
GMT = Path(os.environ.get("SCCC_HALLMARK_GMT", DATA / "reference" / "h.all.v2024.1.Hs.symbols.gmt"))  # the A97 run used v2026.1
PURITY = ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRB", "PTPRC", "CD3E", "CD14", "LYZ", "CD74"]
FOLDER = "cancer"


def _zmean(lcpm: pd.DataFrame, genes) -> pd.Series:
    g = [x for x in genes if x in lcpm.columns]
    z = (lcpm[g] - lcpm[g].mean()) / lcpm[g].std(ddof=0).replace(0, np.nan)
    return z.mean(axis=1)


def hallmark_g2m() -> list[str]:
    """Hallmark G2M checkpoint genes without the clock genes (MSigDB Hallmark, `$SCCC_DATA/reference`)."""
    for line in GMT.read_text().splitlines():
        f = line.split("\t")
        if f[0] == "HALLMARK_G2M_CHECKPOINT":
            return [g for g in f[2:] if g not in GENES]
    raise KeyError("HALLMARK_G2M_CHECKPOINT")


sys.path.insert(0, str(Path(__file__).resolve().parents[4]))  # paper/
from common.paths import group_out  # noqa: E402
def load_tumours(project: str, subsample: int | None = None, seed: int = 0):
    """All primary tumours of a TCGA project (one per patient, lowest barcode, any histology; unpaired) and the 14-gene
    counts of its solid-tissue normals. Returns (tumours, normals): tumours has unit, barcode, depth, the 14 clock genes,
    `prolif` (mean z of Hallmark G2M, clock genes removed) and `purity` (mean z of stromal / immune markers; high = low
    purity) of log2 CPM z-scored within the project, the clinical fields and `er`; normals has unit, depth, genes. The
    genome-wide counts are reduced once and cached under paper/results/figure_5/output/cancer/cache/. `subsample`: that many
    random tumours (and no normals reduction)."""
    cache = group_out(FOLDER) / "cache"
    cache.mkdir(parents=True, exist_ok=True)
    f = cache / f"tcga_tumours_{project}.parquet"
    fn = cache / f"tcga_normals_{project}.parquet"
    if not f.exists():
        b = tcga_builder()
        s = b.samples_of(project)
        t = s[s.sample_type == "Primary Tumor"].sort_values("barcode").drop_duplicates("patient")
        n = s[s.sample_type == "Solid Tissue Normal"].sort_values("barcode").drop_duplicates("patient")
        ct = b.counts_of(project, t)
        depth = ct.sum(axis=1)
        gm = hallmark_g2m()
        need = [g for g in sorted(set(gm) | set(PURITY)) if g in ct.columns]
        lc = np.log2(ct[need].div(depth, axis=0) * 1e6 + 1)
        md = b._md(project, "tcga").drop_duplicates("tcga_barcode").set_index("tcga_barcode").reindex(t.barcode)
        er = md["xml_breast_carcinoma_estrogen_receptor_status"].to_numpy() if "xml_breast_carcinoma_estrogen_receptor_status" in md else np.nan
        clin = md["cgc_case_clinical_stage"].to_numpy() if "cgc_case_clinical_stage" in md else np.nan
        tum = pd.DataFrame({"unit": t.patient.to_numpy(), "barcode": t.barcode.to_numpy(), "depth": depth.to_numpy(),
                            **{c: t[c].to_numpy() for c in ("age", "stage", "vital", "t_death", "t_fu")},
                            "clinical_stage": clin, "er": er,
                            "prolif": _zmean(lc, gm).to_numpy(), "purity": _zmean(lc, PURITY).to_numpy()})
        clock = b.T.clock_columns(ct.columns)
        tum = pd.concat([tum, pd.DataFrame(ct[clock].to_numpy(float), columns=GENES)], axis=1)
        tum.to_parquet(f, index=False)
        if len(n):
            cn = b.counts_of(project, n)
            nor = pd.DataFrame({"unit": n.patient.to_numpy(), "depth": cn.sum(axis=1).to_numpy(),
                                **pd.DataFrame(cn[clock].to_numpy(float), columns=GENES).to_dict("series")})
            nor.to_parquet(fn, index=False)
    tum = pd.read_parquet(f)
    nor = pd.read_parquet(fn) if fn.exists() else None
    if subsample is not None:
        tum = tum.sample(min(subsample, len(tum)), random_state=seed).reset_index(drop=True)
    return tum, nor


# ---- single-cell cohorts --------------------------------------------------------------------------------------------

def sc_builder():
    return data._builder("cancer_sc")


def load_sc(subsample: int | None = None, seed: int = 0):
    """Pooled table of the 29 cohorts with copy-number calls (donor = `study|patient`, stratum = `study|stratum`)
    and its donor table. `subsample` N: N random patients with malignant and normal epithelium and a stromal host,
    plus N random others (check only)."""
    pb, don = data.load("cancer_sc")
    if subsample is not None:
        rng = np.random.default_rng(seed)
        paired = sc_builder().paired_epithelium(pb)
        others = sorted(set(pb.donor) - set(paired))
        keep = set(rng.permutation(paired)[:subsample]) | set(rng.permutation(others)[:subsample])
        pb = pb[pb.donor.isin(keep)].reset_index(drop=True)
        don = don.loc[don.index.isin(keep)]
    return pb, don


# ---- local: not in clocksync -----------------------------------------------------------------------------------------

def plane(waveform: cs.Waveform = WF) -> np.ndarray:
    """local: orthogonal projector (genes x genes) on the clock plane of the gene-centred waveform design (the
    package keeps it private; same formula as `clocksync._core.study.plane`)."""
    Xc = waveform.X - waveform.X.mean(0)
    return Xc @ np.linalg.solve(Xc.T @ Xc, Xc.T)


def share_cov(S: np.ndarray, P: np.ndarray) -> float:
    """local: clock-plane share ||P S P||^2 / ||S||^2 of a gene x gene covariance."""
    return float(np.sum((P @ S @ P) ** 2) / np.sum(S ** 2))


def pooled_share(y: np.ndarray, pv: np.ndarray, k: int, P: np.ndarray) -> float:
    """local: share of the pooled-unit covariance with the Poisson diagonal subtracted (A86); y = log fractions
    centred within k groups, pv = mean Poisson variance per gene."""
    return share_cov(y.T @ y / (len(y) - k) - np.diag(pv), P)


def centred_units(d: pd.DataFrame, group: str = "stratum"):
    """local: log fractions of the clock genes of the rows of `d`, centred within `group` (study x stratum), groups of
    a single unit dropped (A86). Returns (y, mean Poisson variance per gene, number of groups)."""
    grp = d[group].astype(str)
    d = d[grp.map(grp.value_counts()).to_numpy() >= 2]
    return logfrac(d, GENES, d[group].astype(str).to_numpy())


def wrap_h(h):
    return (np.asarray(h, float) + 12.0) % 24.0 - 12.0


def mean_ci(x) -> tuple[float, float]:
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    return (float(np.quantile(x, 0.025)), float(np.quantile(x, 0.975))) if len(x) else (np.nan, np.nan)


# ---- local: not in clocksync (the joint factor model of the cohort pipeline) ---------------------------------------------
# `Dataset.study` is the package's core joint-factor `Study`; the public `Dataset.fit` is the anchored estimator. The
# pre-registered M4 pipeline (docs/cancer/ANALYSIS_PLAN.md) uses the joint model with a reference compartment (E1-E3,
# contrasts H1-H3), the gene loadings (E4) and the own rhythm (E5), none of which has a public call yet.

def joint_dataset(pb: pd.DataFrame, ref: str, name: str = "") -> cs.Dataset:
    """Dataset whose joint-model reference compartment is `ref` (anchors = [ref]; the anchored `fit` is not used)."""
    return cs.Dataset(pb, anchors=[ref], name=name)


def joint_estimates(ds: cs.Dataset, n_boot: int, seed: int = 0) -> pd.DataFrame:
    """local: E1-E3 of the joint model (reference compartment amplitude 1, offset 0) with resection-bootstrap intervals."""
    return ds.study.estimates(n_boot=n_boot, seed=seed)


def joint_contrasts(ds: cs.Dataset, pairs, n_boot: int, seed: int = 0) -> pd.DataFrame:
    """local: log amplitude ratio and offset difference (h) of compartment pairs (H2, H3), bootstrap intervals."""
    return ds.study.contrasts(pairs, n_boot=n_boot, seed=seed)


def own_rhythm(ds: cs.Dataset, k: str, n_perm: int, seed: int = 0) -> pd.Series:
    """local: E5, clock-plane variance of compartment k's residual after its host-coupled part (rho > 0: free-running)."""
    return ds.study.own_rhythm(k, n_perm=n_perm, seed=seed)


def gene_loadings(ds: cs.Dataset, k: str, n_perm: int, seed: int = 0):
    """local: E4, per-gene timing deviation (dev_h) and amplitude ratio of compartment k from the conserved waveform."""
    return ds.study.gene_loadings(k, n_perm=n_perm, seed=seed)

