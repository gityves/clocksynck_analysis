"""A332 D1: loudness (clock-plane variance ratio) of lesional against non-lesional skin, psoriasis and eczema.

Design: paper/docs/registered_experiments.md, "## A332", part D1 (descriptive, not blind; no verdict). Estimator of A254
H254.3 / A183 (`clocksync.paired`): both samples of a series on one projector, states centred within sample type,
r = log ratio of the clock-plane variance lesion / non-lesion, amplitude ratio exp(r / 2); series pooled by inverse
variance (bootstrap SE); two-sided swap p (10,000 per-person label swaps, null centred on its mean as A254); 200
expression-matched and 200 responsive pseudo-clock sets per series (the clock's weights, same pooling); "beyond
companions" = r outside the central 90 % of both. Seed 332.

Data
  GSE154200 (A183 preparation, `A183 score_a183.load`): psoriasis, nummular eczema, atopic eczema, pairs per batch as A183
    (batches with >= 10 pairs of a diagnosis, reference gate p < 0.05).
  A115 replication: four psoriasis series (GSE13355, GSE121212, GSE201827, GSE268958), lesional against uninvolved skin of
    the same person. The lesional and uninvolved rows are those of skill `skin_lesional/processed/` (heldout for
    GSE13355), whose uninvolved rows `processed_a115/` was verified to equal (A115 acceptance); `processed_a115/` itself
    holds only uninvolved and healthy rows.

Usage    python d1_skin_loudness.py [--smoke]     (--smoke: 300 swaps, 5 pseudo sets, prints only)
Outputs  output/d1_series.csv (per diagnosis / series / pooled), output/d1_cohorts.csv (per batch or series)
"""
from __future__ import annotations

import argparse
import importlib.util
import sys
import time
import zlib
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow.parquet as pq

HERE = Path(__file__).resolve().parent
PAPER = HERE.parents[2]
sys.path.insert(0, str(PAPER))
import clocksync as cs  # noqa: E402
from clocksync import paired as P  # noqa: E402
from clocksync.config import limit_threads  # noqa: E402
from common.paths import DATA  # noqa: E402

OUT = HERE.parent / "output"
SEED = 332
WF = cs.gtex()
ROOT = DATA / "skin_lesional"
A115 = {"GSE13355": ROOT / "processed_heldout/GSE13355.parquet", "GSE121212": ROOT / "processed/GSE121212.parquet",
        "GSE201827": ROOT / "processed/GSE201827.parquet", "GSE268958": ROOT / "processed/GSE268958.parquet"}
NONGENE = {"depth", "unit", "condition", "cohort", "diagnosis", "visit", "weeks", "treatment", "arm", "pasi", "platform", "dataset"}
DIAGS = {"psoriasis": "GSE154200 psoriasis", "Nummular Eczema": "GSE154200 nummular eczema", "Atopic eczema": "GSE154200 atopic eczema"}


def say(*a):
    print(time.strftime("%F %T"), *a, flush=True)


def p_two(nul, obs):
    nul = np.asarray(nul)
    return float((1 + np.sum(np.abs(nul) >= abs(obs))) / (1 + len(nul)))


def a183():
    spec = importlib.util.spec_from_file_location("score_a183", PAPER / "exploratory/A183_skin_gse154200/scripts/score_a183.py")
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def gene_means(x, counts):
    """Mean log expression of all genes in the non-lesional samples (as A254): counts log mean fraction (mean count >= 1)."""
    g = [c for c in x.columns if c not in NONGENE and pd.api.types.is_numeric_dtype(x[c])]
    if counts:
        return np.log(x[g].div(x.depth.to_numpy(float), axis=0).mean() + 1e-12)[x[g].mean() >= 1]
    return (x[g].mean() * np.log(2.0)).dropna()


def a115_samples(series, counts_of):
    path = A115[series]
    names = pq.ParquetFile(path).schema.names
    d = pd.read_parquet(path)
    if "BMAL1" in d.columns and "ARNTL" not in d.columns:
        d = d.rename(columns={"BMAL1": "ARNTL"})
    d = d[d.status.isin(["lesional", "nonlesional"]) & (d.disease == "psoriasis")].rename(columns={"patient": "unit", "status": "condition"})
    d["unit"] = d.unit.astype(str)
    d["cohort"] = series
    counts_of[series] = bool((d.value_type == "counts").all())
    meta = ["sample", "disease", "value_type"]
    return d.drop(columns=[c for c in meta if c in d]).reset_index(drop=True)


class Group:
    """Series of one row group: per cohort the clock's r, SE, gate; companions; pooling and swap null."""

    def __init__(self, name, samples, counts, means, n_perm, n_sets, tag):
        self.name, self.n_perm, self.n_sets = name, n_perm, n_sets
        present = [g for g in WF.genes if g in samples.columns]
        self.wf = WF.subset(present)
        res = P.paired(samples, unit="unit", condition="condition", levels=("nonlesional", "lesional"), cohort="cohort", counts=counts,
                       waveform=self.wf, n_perm=1, n_boot=1000, seed=SEED + tag)
        self.cohorts = res.cohorts
        self.res = res
        self.keys = [k for k in res.cohorts[res.cohorts.gated & np.isfinite(res.cohorts.r) & (res.cohorts.se > 0)].cohort]
        self.data, self.comp = res.data, {}
        if not self.keys:
            return
        self.cd = {}
        for k, x in samples.groupby("cohort"):
            if k in self.keys:
                case, ref = P._split(x, "unit", "condition", ("nonlesional", "lesional"), counts)
                self.cd[k] = (case, ref)
                c = counts if isinstance(counts, bool) else bool((x[counts] == "counts").all())
                for kind, t in (("matched", 4), ("responsive", 5)):
                    rng = np.random.default_rng([SEED, tag, t])
                    self.comp.setdefault(kind, {})[k] = np.array([P.log_var_ratio(*P.cohort_states(case, ref, c, w))
                                                                  for w in P.pseudo_waveforms(means[k], n_sets, rng, self.wf, kind)])
        self.rng = np.random.default_rng([SEED, tag, 6])

    def row(self, label, keys=None):
        keys = keys or self.keys
        t = self.cohorts.set_index("cohort").loc[keys]
        iv = 1 / t.se.to_numpy() ** 2
        wn = iv / iv.sum()
        r, se = float(np.sum(wn * t.r)), float(1 / np.sqrt(iv.sum()))
        nul = P.swap_null([self.data[k] for k in keys], iv, n_perm=self.n_perm, rng=np.random.default_rng([SEED, zlib.crc32(label.encode()) % 1000, 3]))
        mu = float(np.nanmean(nul))
        out = dict(group=label, cohorts=";".join(map(str, keys)), pairs=int(t.n_pairs.sum()), r=r, se=se, amplitude_ratio=float(np.exp(r / 2)),
                   p_swap_two_sided=p_two(nul - mu, r - mu), n_swaps=self.n_perm)
        for kind in ("matched", "responsive"):
            v = np.nansum(np.array([self.comp[kind][k] for k in keys]) * wn[:, None], axis=0)
            out.update({f"{kind}_q05": float(np.quantile(v, .05)), f"{kind}_median": float(np.median(v)), f"{kind}_q95": float(np.quantile(v, .95)),
                        f"{kind}_percentile": float(np.mean(v <= r))})
        out["beyond_companions"] = bool(all((r < out[f"{k}_q05"]) or (r > out[f"{k}_q95"]) for k in ("matched", "responsive")))
        return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--smoke", action="store_true")
    a = ap.parse_args(argv)
    limit_threads(4)
    n_perm, n_sets = (300, 5) if a.smoke else (10_000, 200)
    rows, coh = [], []

    # ---- GSE154200 per diagnosis (A183 preparation) -------------------------------------------------------------
    A = a183()
    s = A.load()
    means = A.gene_means(s)
    for i, (d, label) in enumerate(DIAGS.items()):
        g = Group(label, s[s.diagnosis == d], True, means, n_perm, n_sets, tag=10 + i)
        c = g.cohorts.assign(group=label)
        coh.append(c)
        if not g.keys:
            rows.append(dict(group=label, cohorts="", pairs=0, note="no gated batch with >= 10 pairs"))
            say(label, "not estimable")
            continue
        rows.append(g.row(label))
        say(label, {k: round(v, 3) if isinstance(v, float) else v for k, v in rows[-1].items()})
    # ---- A115 replication ---------------------------------------------------------------------------------------
    counts_of = {}
    frames, mm = [], {}
    series = {}
    for i, k in enumerate(A115):
        x = a115_samples(k, counts_of)
        nl = x[x.condition == "nonlesional"]
        mm[k] = gene_means(nl.assign(depth=nl.depth.to_numpy(float)) if counts_of[k] else nl, counts_of[k])
        series[k] = Group(f"A115 {k}", x, counts_of[k], mm, n_perm, n_sets, tag=20 + i)
        coh.append(series[k].cohorts.assign(group=f"A115 {k}"))
        if series[k].keys:
            rows.append(series[k].row(f"A115 {k}"))
        else:
            rows.append(dict(group=f"A115 {k}", cohorts="", pairs=0, note="gate failed"))
        say(rows[-1]["group"], {q: round(v, 3) if isinstance(v, float) else v for q, v in rows[-1].items()})
    # pooled across the four series (inverse variance, gated series), companions and null on the same weights
    ok = [k for k in A115 if series[k].keys]
    pool = Group.__new__(Group)
    pool.n_perm, pool.n_sets = n_perm, n_sets
    pool.cohorts = pd.concat([series[k].cohorts for k in ok], ignore_index=True)
    pool.data = {k: series[k].data[k] for k in ok}
    pool.comp = {kind: {k: series[k].comp[kind][k] for k in ok} for kind in ("matched", "responsive")}
    pool.keys = ok
    rows.append(pool.row("A115 pooled (inverse variance)", ok))
    say(rows[-1]["group"], {q: round(v, 3) if isinstance(v, float) else v for q, v in rows[-1].items()})
    if not a.smoke:
        pd.DataFrame(rows).to_csv(OUT / "d1_series.csv", index=False)
        pd.concat(coh, ignore_index=True).to_csv(OUT / "d1_cohorts.csv", index=False)
    print(pd.DataFrame(rows).round(3).T.to_string())


if __name__ == "__main__":
    main()
