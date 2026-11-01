"""A183 scoring (registered_experiments.md, commit f2445c1): does the psoriatic lesion delay replicate blind in GSE154200,
and is it specific to psoriasis?

H183.1  psoriasis offset (lesional minus non-lesional, h; negative = lesion delayed) < 0: `clocksync.paired`
        (cohort = batch, GTEx waveform, per-person label-swap null, 10,000).
H183.2  psoriasis offset - mean(nummular, atopic offsets) < 0: diagnosis labels permuted within batch (10,000).
        Offsets per diagnosis: pooled over batches (one projector per batch built from all persons of the batch, so
        diagnoses share a frame), weighted resultant with weights = persons (fixed under within-batch permutation).
Holm over the two. Companions: 200 expression-matched pseudo-clock and 200 responsive sets (`pseudo_waveforms`),
per batch, same statistic; clock-specific only beyond the q05 of both. Descriptives as registered.

    python score_a183.py [--n_perm 10000] [--n_pseudo 200] [--n_boot 1000]
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

import clocksync as cs
from clocksync import paired as P
from clocksync import stats as S
from clocksync.config import limit_threads
from clocksync.genes import canonical

HERE = Path(__file__).resolve()
sys.path.insert(0, str(HERE.parents[2] / "scouts_2026-11-01/prep_gse154200/scripts"))
import prep_gse154200 as PR  # noqa: E402  (titles, SCCC, DIAG; reuse of the prep code)

OUT = HERE.parents[1] / "output"
DIAGS = ["psoriasis", "Nummular Eczema", "Atopic eczema"]
MIN_N = 3  # persons of a diagnosis in a batch for it to contribute to the own-framework pooled offset
SEED = 20261101


def say(*a):
    print(time.strftime("%F %T"), *a, flush=True)


def wrap(h):
    return (np.asarray(h, float) + 12.0) % 24.0 - 12.0


def load():
    """Samples table (one row per person x state, all genes as columns) of the prep's persons."""
    import gzip
    m = pd.read_csv(PR.OUT / "gse154200_samples.csv").merge(PR.titles(), on="gsm")
    with gzip.open(PR.SCCC / "GSE154200_countMat.txt.gz", "rt") as f:
        c = pd.read_csv(f, sep="\t", index_col=0)
    c.index = [canonical(g) for g in c.index]
    c = c.groupby(level=0).sum()
    m = m[m.title.isin(c.columns)].rename(columns={"patient id": "patient", "batch id": "batch"})
    m = m[m.state.isin(["lesional", "nonlesional"])]
    both = m.groupby("patient").state.nunique() == 2
    m = m[m.patient.isin(both[both].index)]
    m["patient"] = m.patient.astype(str)
    m = m.drop_duplicates(["patient", "state"]).set_index("title")
    odd = m.groupby("patient").batch.nunique()
    m = m[~m.patient.isin(odd[odd > 1].index)]
    tab = pd.read_csv(PR.OUT / "clock_table.csv", dtype={"donor": str, "stratum": str})
    assert set(tab.donor) == set(m.patient), "persons differ from the prep table"
    depth = c.sum(axis=0)
    x = c.T.loc[m.index].astype(float)
    s = pd.DataFrame({"unit": m.patient.to_numpy(), "condition": m.state.to_numpy(), "cohort": m.batch.astype(str).to_numpy(),
                      "diagnosis": m.diagnosis.to_numpy(), "depth": depth.loc[m.index].to_numpy(float)}, index=m.index)
    s = pd.concat([s, x], axis=1)
    say("samples", s.shape, s.drop_duplicates("unit").diagnosis.value_counts().to_dict())
    return s


def gene_means(s):
    """batch -> mean log fraction of all genes in the non-lesional samples of the batch (mean count >= 1)."""
    meta = ["unit", "condition", "cohort", "diagnosis", "depth"]
    out = {}
    for b, g in s[s.condition == "nonlesional"].groupby("cohort"):
        x = g.drop(columns=meta)
        out[b] = np.log(x.div(g.depth.to_numpy(float), axis=0).mean() + 1e-12)[x.mean() >= 1]
    return out


class Design:
    """Per batch: persons (all diagnoses), diagnosis, paired states on one projector for a waveform."""

    def __init__(self, s):
        self.b = {}
        for k, g in s.groupby("cohort"):
            case = g[g.condition == "lesional"].set_index("unit")
            ref = g[g.condition == "nonlesional"].set_index("unit")
            u = sorted(set(case.index) & set(ref.index))
            self.b[k] = (case.loc[u], ref.loc[u], case.loc[u].diagnosis.to_numpy())
        self.keys = list(self.b)

    def states(self, wf):
        return {k: P.cohort_states(c, r, True, wf) for k, (c, r, _) in self.b.items()}

    def labels(self):
        return {k: v[2] for k, v in self.b.items()}


def offset_of(st, lab, d, idx=None, boot=None, rng=None):
    """Pooled offset, coherence and r of diagnosis d (batches with >= MIN_N persons; weights = persons)."""
    cs_, w, rs = [], [], []
    for k, b in st.items():
        i = np.flatnonzero(lab[k] == d)
        if len(i) < MIN_N:
            continue
        if boot is not None:
            i = rng.choice(i, len(i))
        bb = tuple(x[i] for x in b)
        cs_.append(P.ccorr(*bb))
        w.append(len(i))
        rs.append(P.log_var_ratio(*bb))
    if not cs_:
        return np.nan, np.nan, np.nan
    w = np.asarray(w, float)
    z = np.sum(w * np.asarray(cs_)) / w.sum()
    r = np.nansum(w * np.asarray(rs)) / w[np.isfinite(rs)].sum() if np.isfinite(rs).any() else np.nan
    return float(np.angle(z) * P.HOURS), float(abs(z)), float(r)


def pooled_diff(st, lab, boot=False, rng=None):
    o = [offset_of(st, lab, d, boot=True if boot else None, rng=rng)[0] for d in DIAGS]
    return float(wrap(o[0] - 0.5 * (o[1] + o[2]))), o


def perm_null(st, lab, n_perm, rng):
    """Diagnosis labels permuted within batch (group sizes kept): null of the psoriasis-minus-eczema offset."""
    keys = list(st)
    out = np.empty(n_perm)
    for p in range(n_perm):
        l = {k: rng.permutation(lab[k]) for k in keys}
        out[p] = pooled_diff(st, l)[0]
    return out


def pseudo_states(D, means, wfs):
    """Per batch states of every pseudo set j: list over j of {batch: states}."""
    return [{k: P.cohort_states(D.b[k][0], D.b[k][1], True, w[k]) for k in D.keys} for w in wfs]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n_perm", type=int, default=10000)
    ap.add_argument("--n_pseudo", type=int, default=200)
    ap.add_argument("--n_boot", type=int, default=1000)
    a = ap.parse_args()
    limit_threads(1)
    s = load()
    means = gene_means(s)
    wf = P.gtex()
    D = Design(s)
    lab = D.labels()
    rng = np.random.default_rng(SEED)

    # ---- H183.1: clocksync.paired on psoriasis, cohort = batch ----------------------------------------------------
    ps = s[s.diagnosis == "psoriasis"]
    res = P.paired(ps, unit="unit", condition="condition", levels=("nonlesional", "lesional"), cohort="cohort", counts=True,
                   pseudo_means=means, n_pseudo=a.n_pseudo, n_perm=a.n_perm, seed=SEED, pseudo_responsive=True)
    res.cohorts.to_csv(OUT / "h183_1_paired_cohorts.csv", index=False)
    res.pooled.to_frame("value").to_csv(OUT / "h183_1_paired_pooled.csv")
    say(res.summary())
    say(res.cohorts.round(3).to_string())
    used = str(res.pooled.cohorts).split(";")
    use = res.cohorts.set_index("cohort").loc[used]
    w = 1 / use.se.to_numpy() ** 2
    obs1 = float(res.pooled.offset_h)

    # ---- own-framework statistics ---------------------------------------------------------------------------------
    st = D.states(wf)
    obs2, offs = pooled_diff(st, lab)
    say("own-framework offsets", dict(zip(DIAGS, np.round(offs, 3))), "difference", round(obs2, 3))
    null2 = perm_null(st, lab, a.n_perm, rng)
    p2 = S.perm_p(null2, obs2, "less", tie=0.0)
    p1 = float(res.pooled.p_offset_lower)

    # ---- companions: same statistics on pseudo-clock and responsive sets (per batch, drawn per batch) ---------------
    comp_rows = []
    for kind in ("matched", "responsive"):
        crng = np.random.default_rng([SEED, 7 if kind == "matched" else 8])
        sets = {k: P.pseudo_waveforms(means[k], a.n_pseudo, crng, wf, kind) for k in D.keys}
        for j in range(a.n_pseudo):
            wj = {k: sets[k][j] for k in D.keys}
            # H183.1 statistic: psoriasis-only cohorts of `paired`, weights of the clock fit
            c1 = []
            for k in used:
                case = ps[(ps.cohort == k) & (ps.condition == "lesional")].set_index("unit")
                ref = ps[(ps.cohort == k) & (ps.condition == "nonlesional")].set_index("unit")
                u = sorted(set(case.index) & set(ref.index))
                c1.append(P.ccorr(*P.cohort_states(case.loc[u], ref.loc[u], True, wj[k])))
            o1 = P.pooled_offset(c1, w)
            # H183.2 statistic: all-diagnosis batch projector, within-batch labels as observed
            stj = {k: P.cohort_states(D.b[k][0], D.b[k][1], True, wj[k]) for k in D.keys}
            o2, _ = pooled_diff(stj, lab)
            comp_rows.append(dict(kind=kind, set=j, h183_1_offset=o1, h183_2_diff=o2))
        say("companion done", kind)
    comp = pd.DataFrame(comp_rows)
    comp.to_csv(OUT / "companions_pseudo_sets.csv", index=False)
    rows = []
    for h, obs, p, col in (("H183.1", obs1, p1, "h183_1_offset"), ("H183.2", obs2, p2, "h183_2_diff")):
        row = dict(hypothesis=h, statistic_h=obs, p=p, n_perm=a.n_perm)
        for kind in ("matched", "responsive"):
            v = comp[comp.kind == kind][col].to_numpy()
            row[f"{kind}_q05"] = float(np.quantile(v, 0.05))
            row[f"{kind}_median"] = float(np.median(v))
            row[f"{kind}_pctile_of_clock"] = float(np.mean(v <= obs))
            row[f"beyond_{kind}_q05"] = bool(obs < row[f"{kind}_q05"])
        rows.append(row)
    t = pd.DataFrame(rows)
    order = np.argsort(t.p.to_numpy())
    holm = np.empty(2)
    run = 0.0
    for r_, i in enumerate(order):
        run = max(run, min(1.0, (2 - r_) * t.p.iloc[i]))
        holm[i] = run
    t["p_holm"] = holm
    t["passes_holm_0.05"] = t.p_holm < 0.05
    t["clock_specific"] = t["passes_holm_0.05"] & t.beyond_matched_q05 & t.beyond_responsive_q05
    t["n_persons_used"] = [int(use.n_pairs.sum()), int(sum(len(v) for v in lab.values()))]
    t.to_csv(OUT / "tests.csv", index=False)
    say("\n" + t.round(4).T.to_string())

    # ---- descriptives ---------------------------------------------------------------------------------------------
    brng = np.random.default_rng([SEED, 9])
    boots = []
    for _ in range(a.n_boot):
        o = [offset_of(st, lab, d, boot=True, rng=brng)[0] for d in DIAGS]
        boots.append(o + [float(wrap(o[0] - 0.5 * (o[1] + o[2])))])
    boots = np.array(boots)
    point = offs + [obs2]
    names = DIAGS + ["psoriasis_minus_mean_eczema"]
    rows = []
    for j, nme in enumerate(names):
        dlt = wrap(boots[:, j] - point[j])
        lo, hi = np.quantile(dlt, [0.025, 0.975]) + point[j]
        row = dict(quantity=nme, offset_h=point[j], boot_lo=float(lo), boot_hi=float(hi), n_boot=a.n_boot)
        if j < 3:
            o, c, r = offset_of(st, lab, DIAGS[j])
            row.update(coherence=c, r=r, n_persons=int((np.concatenate(list(lab.values())) == DIAGS[j]).sum()))
        rows.append(row)
    pd.DataFrame(rows).to_csv(OUT / "descriptives_offsets.csv", index=False)
    say("\n" + pd.DataFrame(rows).round(3).to_string())

    rows = []
    for k, b in st.items():
        for d in DIAGS:
            i = np.flatnonzero(lab[k] == d)
            if len(i) < 2:
                rows.append(dict(batch=k, diagnosis=d, n=len(i)))
                continue
            bb = tuple(x[i] for x in b)
            z = P.ccorr(*bb)
            rows.append(dict(batch=k, diagnosis=d, n=len(i), offset_h=float(np.angle(z) * P.HOURS), coherence=float(abs(z)),
                             r=P.log_var_ratio(*bb), in_pooled=len(i) >= MIN_N))
    pd.DataFrame(rows).to_csv(OUT / "descriptives_by_batch.csv", index=False)

    rows = []
    variants = {"all_genes": []} | {f"without_{g}": [g] for g in wf.genes} | {"without_PER1_NR1D1": ["PER1", "NR1D1"]}
    for nme, drop in variants.items():
        w_ = wf.drop(drop) if drop else wf
        d_, o_ = pooled_diff(D.states(w_), lab)
        rows.append(dict(variant=nme, psoriasis=o_[0], nummular=o_[1], atopic=o_[2], psoriasis_minus_mean_eczema=d_))
    pd.DataFrame(rows).to_csv(OUT / "descriptives_gene_sensitivity.csv", index=False)
    say("\n" + pd.DataFrame(rows).round(3).to_string())
    say("done")


if __name__ == "__main__":
    main()
