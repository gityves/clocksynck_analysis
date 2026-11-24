"""Shared pieces of A361-A363 (what the person's clock volume is): GTEx fast-death rows, person covariates and the
registered A360 regression, the person volume of any waveform, gene sets, and the design simulator.

- Rows: A357 / A360 (`a359.load`: 25 peripheral tissues, Whole Blood and Testis excluded) for A361 and A363; A312 rows
  (27 peripheral tissues) for A362.
- Person volume of waveform w: log |posterior-mean host| over all tissues (`clocksync._core.host`, every tissue an
  anchor), donors with >= 5 tissues (A360).
- Simulator: A312 `calib_level.py` cohort at the GTEx design (tissue S and departure from the A294b table, person
  volume SD 0.352 in log, uniform clock time, Poisson counts), extended with a per-tissue elasticity exponent and a
  per-gene extra exponent. No GTEx expression enters a simulation.
GTEx phenotypes and time of death are controlled access: nothing here writes a per-donor value.
"""
from __future__ import annotations

import sys
import warnings
from pathlib import Path

HERE = Path(__file__).resolve()
EXP = HERE.parents[2]
PAPER = EXP.parent
for p in (PAPER, EXP / "A358_A360_gain_sources" / "scripts", EXP / "A143_brain_artefact" / "scripts",
          PAPER / "results" / "figure_2" / "scripts" / "backbone"):
    sys.path.insert(0, str(p))
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

import a359  # noqa: E402
import a360  # noqa: E402
import clocksync as cs  # noqa: E402
from clocksync._core import host as HC  # noqa: E402
from clocksync._core.controls import GC_TARGETS, IMMEDIATE_EARLY  # noqa: E402
from clocksync.genes import canonical  # noqa: E402
from common import gtex  # noqa: E402
from common.paths import DATA  # noqa: E402

OUT = HERE.parents[1] / "output"
META = a359.META
WF = cs.gtex()
TRIO = ["PER1", "PER2", "NR1D1"]
CORE = WF.drop(TRIO)
CORE_IDX = [k for k, g in enumerate(WF.genes) if g in CORE.genes]
A294 = EXP / "A294_floor" / "output" / "tissues.csv"
HALLMARK = DATA / "msigdb_hallmark" / "h.all.v2024.1.Hs.symbols.gmt"
SIGMA_A, DEPTH, TAU = 0.352, 2e7, 0.1


def hallmark(name: str) -> list:
    for line in HALLMARK.read_text().splitlines():
        f = line.split("\t")
        if f[0] == name:
            return f[2:]
    raise KeyError(name)


def inflammatory_genes(columns) -> list:
    """Hallmark TNFA_SIGNALING_VIA_NFKB minus immediate-early genes, glucocorticoid targets and clock genes (A363)."""
    clock = {canonical(g) for g in WF.genes}
    drop = set(IMMEDIATE_EARLY) | set(GC_TARGETS)
    cols = set(columns)
    return [g for g in dict.fromkeys(hallmark("HALLMARK_TNFA_SIGNALING_VIA_NFKB"))
            if g in cols and g not in drop and canonical(g) not in clock]


def lm_of(pb) -> pd.Series:
    """Log mean expression fraction per gene (genes with mean count >= 1): the base of the pseudo-clock draws (A336)."""
    genes = [c for c in pb.columns if c not in META]
    x = pb[genes].to_numpy(np.float32)
    lm = pd.Series(np.log((x / pb.depth.to_numpy(float)[:, None]).mean(0) + 1e-12), index=genes)[x.mean(0) >= 1]
    return lm[~lm.index.duplicated()]


def core_of(w, name="core"):
    """Core-loop positions of a 14-gene (pseudo-)waveform (A336 rule)."""
    from clocksync.waveform import Waveform
    return Waveform(tuple(w.genes[k] for k in CORE_IDX), CORE.amp, CORE.psi, None, None, name)


def volume(pb, comps, w) -> pd.Series:
    """log host radius per donor for waveform w (all tissues anchors)."""
    d = cs.Dataset(pb[META + list(w.genes)], anchors=comps, waveform=w, name=w.name or "w").data
    h = HC.host(d, comps)
    return pd.Series(np.log(np.hypot(h.f[:, 0], h.f[:, 1])), index=[str(r) for r in h.resections])


def persons(pb, comps, gc, ieg, inf=None) -> pd.DataFrame:
    """Person covariates exactly as `a360.main` (base, technical, systemic, demographic), plus the inflammatory state
    (mean over the person's tissues of the within-tissue standardised score) when `inf` is given."""
    d = cs.Dataset(pb[META + list(WF.genes)], anchors=comps, waveform=WF, name="a36x").data
    per = pd.DataFrame(index=pd.Index(sorted(set(d.resection.astype(str)))))
    per["n_tissues"] = pd.Series(d.resection.astype(str)).value_counts()
    for lab, genes in (("gc", gc), ("ieg", ieg)) + ((("inf", inf),) if inf else ()):
        _, z = a359.sample_scores(pb, d, genes)
        per[lab] = z.groupby(level=0).mean()
    _, donors = gtex.load(min_donors=None, genes=list(WF.genes))
    att = a360.A143.sample_attributes(pb)
    per["rin"] = pd.Series(att.SMRIN.to_numpy(float)).groupby(pb.donor.astype(str).to_numpy()).median()
    dn = donors.copy()
    dn.index = dn.index.astype(str)
    per["tod_cos"], per["tod_sin"] = np.cos(2 * np.pi * dn.tod / 24), np.sin(2 * np.pi * dn.tod / 24)
    per["ischaemic_h"], per["age"], per["male"] = dn.ischaemic_h, dn.age, dn.male
    per["bmi"] = gtex.phenotypes(["BMI"]).BMI.astype(float)
    return per[per.n_tissues >= a360.MIN_TISSUES]


BASE, TECH, SYS = a360.BLOCKS["base"], a360.BLOCKS["technical"], a360.BLOCKS["systemic"]


def coefs(per, log_r, terms) -> pd.Series:
    """Standardised coefficients of log_r on `terms` (A360 `std_coefs`) over complete rows."""
    df = per.assign(log_r=log_r.reindex(per.index)).dropna(subset=["log_r"] + terms)
    return pd.Series(a360.std_coefs(df, terms), index=terms)


# ---------------------------------------------------------------------------------------------------------------------
# design simulator (A312 calib_level, extended)

def tissue_table(tissues) -> pd.DataFrame:
    return pd.read_csv(A294, index_col="tissue")[["S", "own"]].loc[list(tissues)]


def simulate(table: pd.DataFrame, presence: pd.Series, wf, rng, n: int = 244, expo: dict | None = None,
             gene_extra: dict | None = None, noise: dict | None = None, person_sd=None) -> pd.DataFrame:
    """Pseudobulk counts (donor, compartment, depth, genes) of a cohort at the GTEx design.
    expo: tissue -> elasticity exponent e_t (clock amplitude of tissue t in person p proportional to a_p ** e_t, scaled
    to unit mean square); gene_extra: gene -> extra exponent x_g (that gene's clock term further multiplied by
    a_p ** x_g / mean); default 1 and 0 (every tissue and gene in proportion).
    noise: tissue -> gene x gene covariance of the per-sample log-expression noise (replaces the default independent
    SD `TAU`); person_sd: per-gene SD of a person offset shared by all of the person's tissues (A363 addendum)."""
    expo, gene_extra = expo or {}, gene_extra or {}
    po = None if person_sd is None else rng.normal(0, 1, (n, len(wf.genes))) * np.asarray(person_sd, float)[None, :]
    a = np.exp(rng.normal(0, SIGMA_A, n))
    theta = rng.uniform(0, 2 * np.pi, n)
    base = np.full(len(wf.genes), -10.0) if wf.mesor is None else np.asarray(wf.mesor, float)
    gx = np.array([gene_extra.get(g, 0.0) for g in wf.genes])
    gfac = a[:, None] ** gx[None, :]
    gfac /= np.sqrt(np.mean(gfac ** 2, 0, keepdims=True))
    rows = []
    for k, r in table.iterrows():
        e = expo.get(k, 1.0)
        ak = a ** e
        ak /= np.sqrt(np.mean(ak ** 2))
        z = np.sqrt(r.S) * np.exp(1j * rng.uniform(0, 2 * np.pi))
        off = rng.normal(0, 0.5, len(wf.genes))
        keep = rng.random(n) < presence[k]
        for p in np.flatnonzero(keep):
            c = ak[p] * z * np.exp(1j * theta[p])
            beta = np.array([c.real, c.imag]) + rng.normal(0, np.sqrt(r.own / 2), 2)
            eps = rng.normal(0, TAU, len(wf.genes)) if noise is None else rng.multivariate_normal(
                np.zeros(len(wf.genes)), noise[k], method="eigh")
            lr = base + off + (wf.X @ beta) * gfac[p] + eps + (0.0 if po is None else po[p])
            rows.append([f"p{p:03d}", k, DEPTH, *rng.poisson(DEPTH * np.exp(lr))])
    return pd.DataFrame(rows, columns=["donor", "compartment", "depth", *wf.genes])


def quiet():
    warnings.filterwarnings("ignore")
