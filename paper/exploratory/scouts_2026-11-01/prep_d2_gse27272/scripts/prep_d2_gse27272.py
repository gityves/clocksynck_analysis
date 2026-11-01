"""D2 scout: mother and fetus at term (GSE27272, Illumina HumanRef-8 v3, GPL6883): prep, gate and power ONLY. NO contrast value.

Design: delivery (`individual`) x compartment {M maternal blood, D cord blood, PL placenta}; M is the host (the mother's internal time),
D and PL are the fetal targets (L2). Design factor in the metadata: maternal smoking (all deliveries vaginal: no labour / delivery-mode
factor; no time of day or labour duration recorded). Values: deposited lumi-quantile log2 (each compartment normalised separately, only
probes with detection p < 0.01 in all samples kept), probes averaged per symbol on the platform annotation (`GPL6883_full.txt`).
clocksync takes counts: the log2 values are put on the linear scale (1000 x 2^x; the centred log fraction undoes it exactly, depth cancels, and
the Poisson measurement variance of the count model becomes negligible, so all noise is carried by the compartment's own psi).

Stages (parse only, nothing downloaded is executed):
  prep   samples table, counts, clock-gene coverage, `clock_table.csv`
  gate   `Dataset.gate` (model-free): all compartments; fetal targets x maternal host; per smoking arm (gate only)
  power  label-permutation nulls of d* (smoking contrast, per fetal target) and target-donor-shuffle null of d* (one-group readout);
         MDE = 2.49 x null SD; stop rule MDE(d*) > 0.3 -> not estimable. No observed coupling or group difference is ever formed.

    python -I prep_d2_gse27272.py --stage prep|gate|power [--n_perm 1000]
"""
from __future__ import annotations

import argparse
import gzip
import re
import time
import zlib
from pathlib import Path

import numpy as np
import pandas as pd

import clocksync as cs
from clocksync import dataset as DS
from clocksync import stats as S
from clocksync import paired as P
from clocksync.config import limit_threads
from clocksync.genes import CLOCK_GENES, canonical

SCCC = Path("/mnt/max8tb/data/CoPhaser_data/cancer_circadian/fetal_gse27272")
OUT = Path(__file__).resolve().parents[1] / "output"
COMP = {"M": "maternal_blood", "D": "cord_blood", "PL": "placenta"}
TARGETS = ["cord_blood", "placenta"]
SD_TO_MDE, DSTAR_STOP, SCALE = 2.49, 0.3, 1000.0


def say(*a):
    print(time.strftime("%F %T"), *a, flush=True)


def header():
    """Series-matrix header -> per-sample table (GSM, title, source, characteristics)."""
    rows = {}
    with gzip.open(SCCC / "GSE27272_series_matrix.txt.gz", "rt", errors="replace") as f:
        for l in f:
            if l.startswith("!series_matrix_table_begin"):
                break
            if not l.startswith("!Sample_"):
                continue
            k, *v = l.rstrip("\n").split("\t")
            v = [x.strip('"') for x in v]
            k = k[len("!Sample_"):]
            if k == "characteristics_ch1":
                for x in v:
                    key, _, val = x.partition(":")
                    rows.setdefault("char_" + key.strip(), []).append(val.strip())
                # one list entry per sample: characteristics rows are separate lines, keep by key
            elif k in ("title", "geo_accession", "source_name_ch1", "platform_id"):
                rows[k] = v
    # characteristics lines repeat per key: rebuild with equal length
    n = len(rows["geo_accession"])
    m = pd.DataFrame({k: v for k, v in rows.items() if len(v) == n})
    return m


def probe_symbols():
    with open(SCCC / "GPL6883_full.txt", errors="replace") as fh:
        lines = [l for l in fh if not l.startswith(("^", "!", "#"))]
    a = pd.read_csv(pd.io.common.StringIO("".join(lines)), sep="\t", dtype=str, usecols=["ID", "Symbol"]).dropna()
    a = a[a.Symbol.str.strip().ne("") & a.Symbol.ne("---")]
    a["sym"] = a.Symbol.map(canonical)
    return a.set_index("ID")["sym"]


def prep():
    m = header()
    m = m.rename(columns={"geo_accession": "gsm", "source_name_ch1": "source"})
    m["suffix"] = m.title.str.extract(r"(\d+)(M|PL|D)$")[1]
    m["individual"] = m.title.str.extract(r"(\d+)(?:M|PL|D)$")[0]
    assert m.suffix.notna().all() and (m["char_individual"] == m.individual).all()
    m["compartment"] = m.suffix.map(COMP)
    m["smoking"] = m.title.str.extract(r"^(non-smoker|smoker)")[0]
    assert (m.smoking == m.source.str.extract(r"(non-smoker|smoker)$")[0]).all()
    m["delivery_mode"] = m["char_mode of delivery"]
    num = lambda c: pd.to_numeric(m[c], errors="coerce")
    for c, n in (("char_age (years)", "age"), ("char_maternal bmi", "bmi"), ("char_gestational age (weeks)", "gest_weeks"), ("char_parity", "parity")):
        m[n] = num(c)
    keep = ["gsm", "title", "individual", "compartment", "smoking", "delivery_mode", "age", "bmi", "gest_weeks", "parity"]
    m[keep].to_csv(OUT / "gse27272_samples.csv", index=False)
    # expression
    with gzip.open(SCCC / "GSE27272_series_matrix.txt.gz", "rt", errors="replace") as f:
        t = [l for l in f]
    i0 = [i for i, l in enumerate(t) if l.startswith("!series_matrix_table_begin")][0]
    i1 = [i for i, l in enumerate(t) if l.startswith("!series_matrix_table_end")][0]
    ex = pd.read_csv(pd.io.common.StringIO("".join(t[i0 + 1:i1])), sep="\t", index_col=0).apply(pd.to_numeric, errors="coerce")
    say("matrix probes x samples", ex.shape, "max", float(np.nanmax(ex.values)), "NaN frac", float(ex.isna().mean().mean()))
    assert np.nanmax(ex.values) < 30  # log2 scale
    sym = probe_symbols()
    inmat = ex.index.isin(sym.index)
    say("probes with a symbol", int(inmat.sum()), "of", len(ex))
    g = ex[inmat].groupby(sym.reindex(ex.index[inmat]).to_numpy()).mean()  # probes -> symbol, mean of log2
    nprobe = ex[inmat].groupby(sym.reindex(ex.index[inmat]).to_numpy()).size()
    say("symbols", g.shape[0])
    clock = [canonical(c) for c in CLOCK_GENES]
    cov = pd.DataFrame({"gene": clock, "present": [c in g.index for c in clock],
                        "n_probes": [int(nprobe.get(c, 0)) for c in clock]})
    for c in COMP.values():
        gs = m[m.compartment == c].gsm
        cov["mean_log2_" + c] = [float(g.loc[x, gs].mean()) if x in g.index else np.nan for x in clock]
        cov["sd_across_samples_" + c] = [float(g.loc[x, gs].std()) if x in g.index else np.nan for x in clock]
    cov.to_csv(OUT / "clock_gene_coverage.csv", index=False)
    say("clock genes present", int(cov.present.sum()), "/", len(clock), "absent", cov.gene[~cov.present].tolist())
    present = [c for c in clock if c in g.index]
    # pseudo-counts: linear scale of the log2 values (clr removes the depth)
    tab = pd.DataFrame({"donor": m.individual.to_numpy(), "compartment": m.compartment.to_numpy()})
    lin = SCALE * 2.0 ** g.loc[present, m.gsm].T.to_numpy(float)  # scaled so the count-model measurement variance is negligible
    tab["depth"] = 100.0 * lin.sum(1)  # >= panel sum; cancels in the centred log fraction
    tab = pd.concat([tab, pd.DataFrame(lin, columns=present)], axis=1)
    tab["smoking"] = m.smoking.to_numpy()
    tab.to_csv(OUT / "clock_table.csv", index=False)
    # counts
    per = m.groupby("individual").compartment.nunique()
    n_comp = m.groupby("compartment").individual.nunique()
    rows = [("samples", len(m)), ("deliveries", m.individual.nunique())]
    rows += [(f"n_{c}", int(n_comp.get(c, 0))) for c in COMP.values()]
    rows += [(f"deliveries_with_{k}_compartments", int((per == k).sum())) for k in (1, 2, 3)]
    rows += [("deliveries_ge2_compartments", int((per >= 2).sum()))]
    sets = {c: set(m[m.compartment == c].individual) for c in COMP.values()}
    for a, b in (("maternal_blood", "cord_blood"), ("maternal_blood", "placenta"), ("cord_blood", "placenta")):
        rows.append((f"deliveries_{a}_and_{b}", len(sets[a] & sets[b])))
    rows.append(("deliveries_all_three", len(sets["maternal_blood"] & sets["cord_blood"] & sets["placenta"])))
    rows.append(("deliveries_fetal_target_with_maternal_host", len((sets["cord_blood"] | sets["placenta"]) & sets["maternal_blood"])))
    d = m.drop_duplicates("individual")
    sm = m.groupby("individual").smoking.nunique()
    rows.append(("deliveries_inconsistent_smoking_label", int((sm > 1).sum())))
    for s in ("smoker", "non-smoker"):
        ids = set(d[d.smoking == s].individual)
        rows.append((f"deliveries_{s}", len(ids)))
        for c in ("cord_blood", "placenta"):
            rows.append((f"deliveries_{s}_{c}_with_host", len(ids & sets[c] & sets["maternal_blood"])))
    rows.append(("delivery_mode", ";".join(f"{k}:{v}" for k, v in m.drop_duplicates("individual").delivery_mode.value_counts().items())))
    pd.DataFrame(rows, columns=["item", "value"]).to_csv(OUT / "counts.csv", index=False)
    say("\n" + pd.DataFrame(rows, columns=["item", "value"]).to_string())


def load():
    tab = pd.read_csv(OUT / "clock_table.csv", dtype={"donor": str})
    donors = tab.drop_duplicates("donor").set_index("donor")[["smoking"]]
    return tab.drop(columns="smoking"), donors


def dataset(tab, donors, keep=None):
    present = [c for c in tab.columns if c not in ("donor", "compartment", "depth")]
    wf = cs.waveform.gtex()
    drop = [g for g in wf.genes if g not in present]
    ids = donors.index if keep is None else keep
    t = tab[tab.donor.isin(ids)]
    return cs.Dataset(t, anchors=["maternal_blood"], donors=donors.loc[sorted(set(t.donor))], drop_genes=drop or None,
                      name="gse27272")


def gate(n_perm):
    tab, donors = load()
    rows = []

    def run(scope, ds, **kw):
        g = ds.gate(n_perm=n_perm, seed=zlib.crc32(scope.encode()) % 2 ** 31, **kw)
        rows.append(dict(scope=scope, donors=ds.table.resection.nunique(), n_pairs=g.n_pairs, n_perm=n_perm,
                         p_pairing=g.p_pairing, p_acrophase=g.p_acrophase, clock_share=g.clock_share,
                         share_null_mean=g.share_null_mean, gate_passed=bool(g.passes)))
        say(rows[-1])

    ds = dataset(tab, donors)
    run("all_compartments_all_deliveries", ds)
    run("fetal_targets_x_maternal_host", ds, targets=TARGETS, hosts=["maternal_blood"])
    for t in TARGETS:
        run(f"{t}_x_maternal_host", ds, targets=[t], hosts=["maternal_blood"])
    run("cord_blood_x_placenta", ds, targets=["cord_blood"], hosts=["placenta"])
    for s in ("non-smoker", "smoker"):  # arms: gate only (needed before any contrast)
        ks = donors.index[donors.smoking == s]
        d2 = dataset(tab, donors, ks)
        run(f"arm_{s}_all_compartments", d2)
        run(f"arm_{s}_fetal_x_host", d2, targets=TARGETS, hosts=["maternal_blood"])
    pd.DataFrame(rows).to_csv(OUT / "gates.csv", index=False)


def null_dstar(ds, lab, targets, n_perm, seed, shuffle_target=False):
    """Null draws of d* (case - reference group, `Dataset.contrast` estimand) per target from permuting the labels;
    the observed difference is never formed. Hosts: pooled anchors, fixed across permutations (as `contrast`)."""
    d = ds.data
    d = d.take(np.flatnonzero(np.isin(d.resection, lab.index)))
    levels = tuple(sorted(set(lab)))
    H = DS.A._hosts(d, [a for a in ds.anchors if a in set(d.compartment)], list(targets))
    idx = {k: pd.Index(np.intersect1d(d.donors(k), H[k].resections)) for k in targets}
    V0 = {k: DS.A.on_host(d, k, H[k], idx[k].to_numpy())["V"] for k in targets}

    def cs_(k, dn):
        if len(dn) < ds.min_donors:
            return np.nan
        e = DS.A.on_host(d, k, H[k], dn)
        s = e["z2"] * V0[k]
        return s / (s + 2 * e["psi"]) if np.isfinite(s) else np.nan

    rng = np.random.default_rng(seed)
    out = {k: np.empty(n_perm) for k in targets}
    for p in range(n_perm):
        lp = S.permute_within(lab, None, rng)
        for k in targets:
            lk = lp.reindex(idx[k]).dropna()
            out[k][p] = cs_(k, lk.index[lk == levels[1]].to_numpy()) - cs_(k, lk.index[lk == levels[0]].to_numpy())
    return out, {k: {l: int((lab.reindex(idx[k]).dropna() == l).sum()) for l in levels} for k in targets}


def pair_states(tab, donors, target):
    """Paired states of (target, maternal blood) per delivery on one projector (arrays: counts=False; log2 x ln 2)."""
    wf = cs.waveform.gtex()
    wf = wf.drop([g for g in wf.genes if g not in tab.columns])
    g = list(wf.genes)
    log2 = lambda c: tab[tab.compartment == c].set_index("donor")[g].pipe(lambda x: np.log2(x / SCALE))
    A_, B_ = log2(target), log2("maternal_blood")
    ids = sorted(set(A_.index) & set(B_.index))
    A_, B_ = A_.loc[ids].assign(depth=1.0), B_.loc[ids].assign(depth=1.0)
    return ids, P.cohort_states(A_, B_, False, wf)


def stat(b, i):
    """(coherence |z|, offset in h) of the deliveries `i` (fetal minus maternal), the A54/A115 paired estimator."""
    z = P.ccorr(*(x[i] for x in b))
    return float(abs(z)), float(np.angle(z) * P.HOURS)


def wrap(h):
    return (h + 12.0) % 24.0 - 12.0


def pair_nulls(b, smoke, n_perm, rng):
    """Smoking-label permutation null of the coherence and offset differences (group sizes kept), and, for the one-group
    readout, the null of the coherence when the fetal states are shuffled across deliveries. Observed values are never formed."""
    n = len(smoke)
    out, one = np.empty((n_perm, 2)), np.empty(n_perm)
    bA, vA, bB, vB = b
    for p in range(n_perm):
        lab = rng.permutation(smoke)
        a, c = stat(b, np.flatnonzero(lab)), stat(b, np.flatnonzero(~lab))
        out[p] = (a[0] - c[0], wrap(a[1] - c[1]))
        j = rng.permutation(n)
        one[p] = abs(P.ccorr(bA[j], vA[j], bB, vB))
    return out, one


def power(n_perm):
    tab, donors = load()
    rows = []
    rng = np.random.default_rng(20261101)
    for t in TARGETS:  # (A) pairwise estimand: fetal compartment against the mother's blood
        ids, b = pair_states(tab, donors, t)
        smoke = (donors.smoking.reindex(ids) == "smoker").to_numpy()
        nul, one = pair_nulls(b, smoke, n_perm, rng)
        base = dict(target=t, host="maternal_blood", estimand="pairwise_complex_coherence", n_perm=n_perm, n_deliveries=len(ids),
                    n_smoker=int(smoke.sum()), n_nonsmoker=int((~smoke).sum()))
        for r, j, unit in (("coherence_diff_smoker_minus_nonsmoker", 0, "coherence units"), ("offset_diff_h", 1, "hours")):
            sd = float(nul[:, j].std(ddof=1))
            rows.append(dict(base, contrast="smoker_vs_nonsmoker", readout=r, unit=unit, null_sd=sd, mde=SD_TO_MDE * sd))
        rows.append(dict(base, contrast="one_group_vs_shuffled_fetal_states", readout="coherence_null_q95", unit="coherence units",
                         null_sd=float(one.std(ddof=1)), mde=float(np.quantile(one, 0.95))))
        say(rows[-3:])
    # (B) d* (spread-adjusted coupling difference): needs a host fitted on >= 2 anchor compartments; the mother has ONE
    # compartment, so the only available host is the free host of all three (each target on the other two, leave-one-out)
    ds = cs.Dataset(tab[tab.compartment.isin(COMP.values())], anchors=None, donors=donors, drop_genes=[g for g in cs.waveform.gtex().genes if g not in tab.columns] or None,
                    name="gse27272_free")
    nul, n = null_dstar(ds, donors.smoking.rename_axis(None), TARGETS, n_perm, seed=20261102)
    for k in TARGETS:
        sd = float(np.nanstd(nul[k], ddof=1))
        rows.append(dict(target=k, host="free host (maternal blood + the other fetal compartment)", estimand="d_star",
                         n_perm=n_perm, n_deliveries=None, n_smoker=n[k].get("smoker"), n_nonsmoker=n[k].get("non-smoker"),
                         contrast="smoker_vs_nonsmoker", readout="d_star_diff", unit="d* units", null_sd=sd, mde=SD_TO_MDE * sd,
                         dstar_stop_rule=DSTAR_STOP, dstar_estimable=bool(SD_TO_MDE * sd <= DSTAR_STOP)))
        say(rows[-1])
    pd.DataFrame(rows).to_csv(OUT / "power.csv", index=False)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", required=True, choices=["prep", "gate", "power"])
    ap.add_argument("--n_perm", type=int, default=1000)
    a = ap.parse_args()
    limit_threads(4)
    {"prep": prep, "gate": lambda: gate(a.n_perm), "power": lambda: power(a.n_perm)}[a.stage]()
