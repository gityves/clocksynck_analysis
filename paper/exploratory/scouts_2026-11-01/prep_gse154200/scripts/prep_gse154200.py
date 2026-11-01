"""GSE154200 lesional / non-lesional skin: counts -> clocksync two-compartment tables, gate and power. NO contrast.

Stage 1 (`--stage prep`): parse the count matrix (untrusted download: parsed only), map libraries MUC#### -> GSM ->
patient / state / diagnosis (`output/gse154200_samples.csv`), check the clock genes, build the pseudobulk table
(donor = patient, compartment = lesional | nonlesional, stratum = batch), write `output/clock_table.csv` and
`output/clock_gene_coverage.csv`.
Stage 2 (`--stage gate`): `Dataset.gate` per diagnosis and pooled -> `output/gates.csv`.
Stage 3 (`--stage power`): label-permutation MDE on d* (`clocksync.power.mde`, A153 conversion) for the between-disease
contrasts, and the within-person swap null of the paired estimand (A115/A54 sense) -> `output/power.csv`.
No observed between-group statistic is computed or written anywhere.

    python prep_gse154200.py --stage prep|gate|power [--jobs 4]
"""
from __future__ import annotations

import argparse
import gzip
import sys
import time
import zlib
from pathlib import Path

import numpy as np
import pandas as pd

import clocksync as cs
from clocksync.config import limit_threads
from clocksync.genes import CLOCK_GENES, canonical
from clocksync import paired as P

SCCC = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian/skin_gse154200")
OUT = Path(__file__).resolve().parents[1] / "output"
DIAG = {"psoriasis": "psoriasis", "nummular": "Nummular Eczema", "atopic": "Atopic eczema"}
POWER, SD_TO_MDE = 0.8, 2.49  # one-sided 0.05 and 80 %: 1.645 + 0.842


def say(*a):
    print(time.strftime("%F %T"), *a, flush=True)


def titles() -> pd.DataFrame:
    rows, cur = [], None
    for l in open(SCCC / "GSE154200_gsm_brief.txt", encoding="utf-8", errors="replace"):
        if l.startswith("^SAMPLE"):
            cur = {"gsm": l.split("=")[1].strip()}
            rows.append(cur)
        elif l.startswith("!Sample_title") and cur is not None:
            cur["title"] = l.split("=", 1)[1].strip()
    return pd.DataFrame(rows)


def prep():
    m = pd.read_csv(OUT / "gse154200_samples.csv").merge(titles(), on="gsm")
    with gzip.open(SCCC / "GSE154200_countMat.txt.gz", "rt") as f:
        c = pd.read_csv(f, sep="\t", index_col=0)
    say("count matrix", c.shape)
    c = c.groupby(level=0).sum()  # duplicated symbols summed
    miss = sorted(set(m.title) - set(c.columns))
    m = m[m.title.isin(c.columns)].copy()
    say(f"libraries in matrix {c.shape[1]}, with metadata {len(m)} (metadata without a column: {len(miss)})")
    m = m.rename(columns={"patient id": "patient", "batch id": "batch"})
    m = m[m.state.isin(["lesional", "nonlesional"])]
    both = m.groupby("patient").state.nunique() == 2
    m = m[m.patient.isin(both[both].index)]
    m["patient"] = m.patient.astype(str)
    m = m.drop_duplicates(["patient", "state"]).set_index("title")  # one library per person x state
    nb = m.groupby("patient").batch.nunique()  # clocksync needs one stratum per donor
    odd = nb[nb > 1].index
    say(f"persons with different batches for the two biopsies, dropped: {len(odd)} ({m[m.patient.isin(odd)].diagnosis.unique().tolist()})")
    m = m[~m.patient.isin(odd)]
    lib = c.sum(axis=0)
    cl = c.T.loc[m.index]
    tab = cs.from_bulk(cl.assign(**{}), m, donor="patient", compartment="state", stratum="batch",
                       depth=lib.reindex(m.index), genes="clock")
    tab["diagnosis"] = tab.donor.map(m.drop_duplicates("patient").set_index("patient").diagnosis)
    genes = [g for g in tab.columns if g not in ("donor", "compartment", "stratum", "depth", "diagnosis")]
    cov = pd.DataFrame({"gene": genes, "median_count": [float(tab[g].median()) for g in genes],
                        "frac_zero": [float((tab[g] == 0).mean()) for g in genes]})
    cov["present"] = True
    absent = [g for g in map(canonical, CLOCK_GENES) if g not in genes]
    cov = pd.concat([cov, pd.DataFrame({"gene": absent, "present": False})], ignore_index=True)
    cov.to_csv(OUT / "clock_gene_coverage.csv", index=False)
    tab.to_csv(OUT / "clock_table.csv", index=False)
    say(f"clock genes present {len(genes)}/{len(CLOCK_GENES)}: absent {absent}; table {tab.shape}, "
        f"persons {tab.donor.nunique()}, depth median {tab.depth.median():.3g}")
    say(tab.groupby("diagnosis").donor.nunique().to_dict())


def load():
    tab = pd.read_csv(OUT / "clock_table.csv", dtype={"donor": str, "stratum": str})
    donors = tab.drop_duplicates("donor").set_index("donor")[["diagnosis"]]
    return tab.drop(columns="diagnosis"), donors


def dataset(tab, donors, diag):
    keep = donors.index[donors.diagnosis.isin(diag)]
    return cs.Dataset(tab[tab.donor.isin(keep)], anchors=["nonlesional"], donors=donors.loc[keep],
                      name="gse154200_" + "+".join(diag))


def gate(n_perm=1000):
    tab, donors = load()
    rows = []
    for name, diag in {**{k: [v] for k, v in DIAG.items()}, "pooled": list(DIAG.values())}.items():
        ds = dataset(tab, donors, diag)
        g = ds.gate(n_perm=n_perm, seed=zlib.crc32(name.encode()) % 2 ** 31)
        rows.append(dict(scope=name, persons=ds.table.resection.nunique(), n_perm=n_perm, p_pairing=g.p_pairing,
                         p_acrophase=g.p_acrophase, clock_share=g.clock_share, gate_passed=bool(g.passes)))
        say(rows[-1])
    pd.DataFrame(rows).to_csv(OUT / "gates.csv", index=False)


def states_all(tab, donors):
    """Paired states (lesional = case, non-lesional = reference) of all persons on ONE projector (all diagnoses)."""
    wf = P.gtex()
    g = list(wf.genes)
    t = tab.pivot(index="donor", columns="compartment")
    pick = lambda c: pd.DataFrame({x: t[(x, c)] for x in g}).assign(depth=t[("depth", c)])
    bA, vA, bB, vB = P.cohort_states(pick("lesional"), pick("nonlesional"), True, wf)
    return t.index.astype(str), (bA, vA, bB, vB), t[("stratum", "lesional")].astype(str).to_numpy()


def stat(b, i):
    """(offset in h, coherence, r) of the persons `i`; the same estimator as `paired.paired` (A54/A115 paired)."""
    bA, vA, bB, vB = (x[i] for x in b)
    z = P.ccorr(bA, vA, bB, vB)
    return float(np.angle(z) * P.HOURS), float(abs(z)), float(P.log_var_ratio(bA, vA, bB, vB))


def wrap(h):
    return (h + 12.0) % 24.0 - 12.0


def label_perm(idx, other, strata, n_perm, rng):
    """Null of (offset difference in h, coherence difference, r difference) between two person sets, from permuting
    the diagnosis labels within batch (group sizes kept). Only the null is kept; the observed difference is never
    formed."""
    pool = np.concatenate([idx, other])
    isA = np.r_[np.ones(len(idx), bool), np.zeros(len(other), bool)]
    out = np.empty((n_perm, 3))
    for p in range(n_perm):
        lab = isA.copy()
        for s in np.unique(strata[pool]):
            j = np.flatnonzero(strata[pool] == s)
            lab[j] = rng.permutation(isA[j])
        a, c = stat(B, pool[lab]), stat(B, pool[~lab])
        out[p] = (wrap(a[0] - c[0]), a[1] - c[1], a[2] - c[2])
    return out


def swap_null(idx, n_perm, rng):
    """Within-person swap null (lesional <-> non-lesional label) of one diagnosis: offset (h) and r under no
    difference between the two biopsies; only the null is returned."""
    bA, vA, bB, vB = (x[idx] for x in B)
    out = np.empty((n_perm, 2))
    for p in range(n_perm):
        sw = rng.random(len(idx)) < 0.5
        b = P._swap((bA - bA.mean(0), vA, bB - bB.mean(0), vB), sw)
        out[p] = (float(np.angle(P.ccorr(*b)) * P.HOURS), float(P.log_var_ratio(*b)))
    return out


def power_stage(n_perm):
    global B
    tab, donors = load()
    ids, B, strata = states_all(tab, donors)
    diag = donors.diagnosis.reindex(ids).to_numpy()
    sets = {k: np.flatnonzero(diag == v) for k, v in DIAG.items()}
    rng = np.random.default_rng(20261101)
    rows = []
    cons = [("psoriasis", "nummular"), ("psoriasis", "atopic"), ("atopic", "nummular")]
    for a, b in cons:
        nul = label_perm(sets[a], sets[b], strata, n_perm, rng)
        rows.append(dict(contrast=f"{a}_vs_{b}", readout="offset_h", n_a=len(sets[a]), n_b=len(sets[b]), n_perm=n_perm,
                         null_sd=nul[:, 0].std(ddof=1)))
        rows.append(dict(contrast=f"{a}_vs_{b}", readout="coherence", n_a=len(sets[a]), n_b=len(sets[b]), n_perm=n_perm,
                         null_sd=nul[:, 1].std(ddof=1)))
        rows.append(dict(contrast=f"{a}_vs_{b}", readout="r_log_var_ratio", n_a=len(sets[a]), n_b=len(sets[b]),
                         n_perm=n_perm, null_sd=nul[:, 2].std(ddof=1)))
    for a in DIAG:
        rest = np.concatenate([sets[k] for k in DIAG if k != a])
        nul = label_perm(sets[a], rest, strata, n_perm, rng)
        for j, r in enumerate(("offset_h", "coherence", "r_log_var_ratio")):
            rows.append(dict(contrast=f"{a}_vs_other_two", readout=r, n_a=len(sets[a]), n_b=len(rest), n_perm=n_perm,
                             null_sd=nul[:, j].std(ddof=1)))
    for a in DIAG:  # lesional vs non-lesional within one diagnosis (A115 paired estimand, swap null)
        nul = swap_null(sets[a], n_perm, rng)
        for j, r in enumerate(("offset_h", "r_log_var_ratio")):
            rows.append(dict(contrast=f"{a}_lesional_vs_nonlesional", readout=r, n_a=len(sets[a]), n_b=len(sets[a]),
                             n_perm=n_perm, null_sd=nul[:, j].std(ddof=1)))
    d = pd.DataFrame(rows)
    d["mde"] = SD_TO_MDE * d.null_sd
    d["unit"] = d.readout.map({"offset_h": "hours", "coherence": "coherence units", "r_log_var_ratio": "log variance ratio"})
    d["dstar_estimable"] = False
    d.to_csv(OUT / "power.csv", index=False)
    say("\n" + d.round(3).to_string())


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", required=True, choices=["prep", "gate", "power"])
    ap.add_argument("--jobs", type=int, default=4)
    ap.add_argument("--n_rep", type=int, default=100)
    ap.add_argument("--n_perm", type=int, default=2000)
    a = ap.parse_args()
    limit_threads(1 if a.jobs > 1 else 4)
    {"prep": prep, "gate": gate}.get(a.stage, lambda: power_stage(a.n_perm))()
